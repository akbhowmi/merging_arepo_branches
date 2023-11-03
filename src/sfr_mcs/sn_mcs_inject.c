/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/sfr_mcs/sn_mcs_inject.c
 * \date        08/2018
 * \author      Matthew C Smith
 * \brief
 * \details     Private to M C Smith, but collaboration encouraged.
                Originally developed in 2015, ported into main repo 2018.
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../allvars.h"
#include "../proto.h"

#if defined(SFR_MCS) && defined(SN_MCS)

#if !defined(VORONOI_DYNAMIC_UPDATE)
#error "SN_MCS requires VORONOI_DYNAMIC_UPDATE."
#endif

#define N_NEIGHBOUR_MAX 256 /* This should never really be reached */

static struct snhost_in
{
  int HostTask;
  int HostIndex;
  double mass_deposited;
  double energy_deposited;
  double starvel[3];
} * SnhostIn, *SnhostGet;

static struct iso_weights
{
  double r;
  double x_plus[3];
  double x_minus[3];
  double omega;
  double w[3];
  int Task;
  int Index;
} * Neighbour_weights;

static struct snvar_in
{
  int Task;
  int Index;
  double dm;
  double dp[3];
  double dE;
#ifdef MECHANICAL_FEEDBACK
  double E_51;
  double p_ej_tot;
#endif
} * SnvarIn, *SnvarGet;

/* Functions used with qsort */
static int snhost_compare(const void *a, const void *b);
static int snvar_compare(const void *a, const void *b);

/* This function finds the host cells of SNe. Distributes SNe variables to hosts,
whether local or on other tasks.
int n_sn_on_task: the number of star particles with a SNe event this timestep on
this task. */
void find_sn_host_cells_and_distribute(int n_sn_on_task)
{
  int i, idx, j;
  int n_sn_processed, nexport, nimport, ngrp, sendTask, recvTask;
  mesh_search_data *searchdata;
  int *starindex;
  MPI_Status status;

  mpi_printf("SN_MCS: Finding SN host cells...\n");

  searchdata = mymalloc("searchdata", n_sn_on_task * sizeof(mesh_search_data));
  starindex  = mymalloc("starindex", n_sn_on_task * sizeof(int));

  n_sn_processed = 0;

  /* Load search data for all star particles with a SNe event this timestep,
  then search.
  To Do: loop over star particle struct. */
  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if((P[i].Type == 4) && (P[i].N_SN > 0))
        {
          searchdata[n_sn_processed].Pos[0] = P[i].Pos[0];
          searchdata[n_sn_processed].Pos[1] = P[i].Pos[1];
          searchdata[n_sn_processed].Pos[2] = P[i].Pos[2];
          starindex[n_sn_processed]         = i;
          n_sn_processed++;
        }
    }

  assert(n_sn_processed == n_sn_on_task); /* Check the number is correct */

  find_nearest_meshpoint_global(searchdata, n_sn_processed, 0, 0);

  /* Check how many particles need to be exported */
  for(i = 0; i < NTask; i++)
    Send_count[i] = 0;

  for(i = 0, nexport = 0; i < n_sn_processed; i++)
    {
      if(searchdata[i].Task != ThisTask)
        {
          Send_count[searchdata[i].Task] += 1;
          nexport++;
        }
    }

  /* Send the export counts to other tasks and get receive counts */
  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(i = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; i < NTask; i++)
    {
      nimport += Recv_count[i];

      if(i > 0)
        {
          Send_offset[i] = Send_offset[i - 1] + Send_count[i - 1];
          Recv_offset[i] = Recv_offset[i - 1] + Recv_count[i - 1];
        }
    }

  /* Allocate export/import buffers */
  SnhostGet = (struct snhost_in *)mymalloc("SnhostGet", nimport * sizeof(struct snhost_in));
  SnhostIn  = (struct snhost_in *)mymalloc("SnhostIn", nexport * sizeof(struct snhost_in));

  /* Distribute SN variables either locally or to the export buffer */
  for(i = 0, j = 0; i < n_sn_processed; i++)
    {
      if(searchdata[i].Task == ThisTask)
        {
          if(P[searchdata[i].u.Index].ID != 0)
            {
              SphP[searchdata[i].u.Index].mass_deposited   = P[starindex[i]].M_ej;
              SphP[searchdata[i].u.Index].energy_deposited = P[starindex[i]].N_SN * All.SupernovaEnergy * All.cf_atime * All.cf_atime;
              SphP[searchdata[i].u.Index].starvel[0] += P[starindex[i]].Vel[0];
              SphP[searchdata[i].u.Index].starvel[1] += P[starindex[i]].Vel[1];
              SphP[searchdata[i].u.Index].starvel[2] += P[starindex[i]].Vel[2];
              SphP[searchdata[i].u.Index].N_SN_hosted += 1;
#ifdef SN_MCS_LOG
              sn_add_to_log(searchdata[i].u.Index);
#endif
            }
        }
      else
        {
          SnhostIn[j].HostTask         = searchdata[i].Task;
          SnhostIn[j].HostIndex        = searchdata[i].u.Index;
          SnhostIn[j].mass_deposited   = P[starindex[i]].M_ej;
          SnhostIn[j].energy_deposited = P[starindex[i]].N_SN * All.SupernovaEnergy * All.cf_atime * All.cf_atime;
          SnhostIn[j].starvel[0]       = P[starindex[i]].Vel[0];
          SnhostIn[j].starvel[1]       = P[starindex[i]].Vel[1];
          SnhostIn[j].starvel[2]       = P[starindex[i]].Vel[2];
          j++;
        }

      /* Reset star properties */
      P[starindex[i]].M_ej = 0;
      P[starindex[i]].N_SN = 0;
    }

  assert(j == nexport); /* Another check that buffers have been filled correctly */

  /* Get export buffer into correct order */
  qsort(SnhostIn, nexport, sizeof(struct snhost_in), snhost_compare);

  /* Exchange particle data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&SnhostIn[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct snhost_in), MPI_BYTE, recvTask,
                           TAG_DENS_A, &SnhostGet[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct snhost_in), MPI_BYTE,
                           recvTask, TAG_DENS_A, MPI_COMM_WORLD, &status);
            }
        }
    }

  myfree(SnhostIn);

  /* Now distribute imported SN variables */
  for(i = 0; i < nimport; i++)
    {
      if(P[SnhostGet[i].HostIndex].ID != 0)
        {
          assert(SnhostGet[i].HostTask == ThisTask);   /* Check we have right task */
          assert(P[SnhostGet[i].HostIndex].Type == 0); /* Check we have a gas cell */
          SphP[SnhostGet[i].HostIndex].mass_deposited   = SnhostGet[i].mass_deposited;
          SphP[SnhostGet[i].HostIndex].energy_deposited = SnhostGet[i].energy_deposited;
          SphP[SnhostGet[i].HostIndex].starvel[0] += SnhostGet[i].starvel[0];
          SphP[SnhostGet[i].HostIndex].starvel[1] += SnhostGet[i].starvel[1];
          SphP[SnhostGet[i].HostIndex].starvel[2] += SnhostGet[i].starvel[2];
          SphP[SnhostGet[i].HostIndex].N_SN_hosted += 1;
#ifdef SN_MCS_LOG
          sn_add_to_log(SnhostGet[i].HostIndex);
#endif
        }
    }

  myfree(SnhostGet);
  myfree(starindex);
  myfree(searchdata);
}

/* Sorts by task, then by index */
static int snhost_compare(const void *a, const void *b)
{
  if(((struct snhost_in *)a)->HostTask < (((struct snhost_in *)b)->HostTask))
    return -1;

  if(((struct snhost_in *)a)->HostTask > (((struct snhost_in *)b)->HostTask))
    return +1;

  if(((struct snhost_in *)a)->HostIndex < (((struct snhost_in *)b)->HostIndex))
    return -1;

  if(((struct snhost_in *)a)->HostIndex > (((struct snhost_in *)b)->HostIndex))
    return +1;

  return 0;
}

void stellar_feedback_distribute(void)
{
  int i, idx, j, k;
  int n_sn_host_on_task, nexport, nimport, ngrp, sendTask, recvTask;
  int q, dpi, vf, particle, n_neighbour, index;
  double dx, dy, dz;
  double sum_omega_x_plus[3], sum_omega_x_minus[3], factor_plus[3], factor_minus[3], sum_weights;
  double p_ej_tot, mod_weight, dp[3], dm, dp_sq0, kin_energy, u0, u1, dv[3];
  MPI_Status status;

#ifdef MECHANICAL_FEEDBACK
  double n, p_terminal, E_51, Zsun, boost_factor;
#endif
  mpi_printf("SN_MCS: Distributing feedback...\n");

  /* Find number of SNe hosting cells on this task, also compute average of velocity vector */
  n_sn_host_on_task = 0;
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;
      if(SphP[i].N_SN_hosted > 0)
        {
          SphP[i].starvel[0] /= SphP[i].N_SN_hosted;
          SphP[i].starvel[1] /= SphP[i].N_SN_hosted;
          SphP[i].starvel[2] /= SphP[i].N_SN_hosted;
          n_sn_host_on_task++;
        }
    }

  for(i = 0; i < NTask; i++)
    Send_count[i] = 0;

  SnvarIn = (struct snvar_in *)mymalloc("SnvarIn", N_NEIGHBOUR_MAX * n_sn_host_on_task * sizeof(struct snvar_in));

  Neighbour_weights = (struct iso_weights *)mymalloc("Neighbour_weights", N_NEIGHBOUR_MAX * sizeof(struct iso_weights));

  nexport = 0;

  /* Find cells hosting SNe and compute weights */
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;
      if(SphP[i].N_SN_hosted > 0)
        {
          q           = SphP[i].first_connection;
          n_neighbour = 0;
          while(q >= 0)
            {
              if(n_neighbour == N_NEIGHBOUR_MAX)
                terminate("Number of neighbouring cells too large, increase N_NEIGHBOUR_MAX. n_neighbour %d\n", n_neighbour);

              dpi      = DC[q].dp_index;
              vf       = DC[q].vf_index;
              particle = Mesh.DP[dpi].index;

              if(particle < 0)
                {
                  q = DC[q].next;
                  continue;
                }

              if(Mesh.DP[dpi].task == ThisTask)
                {
                  if(particle >= NumGas && Mesh.DP[dpi].task == ThisTask)
                    particle -= NumGas;

                  Neighbour_weights[n_neighbour].Index = particle;
                }
              else
                Neighbour_weights[n_neighbour].Index = Mesh.DP[dpi].originalindex;

              Neighbour_weights[n_neighbour].Task = Mesh.DP[dpi].task;

              dx = Mesh.DP[dpi].x - P[i].Pos[0];
              dy = Mesh.DP[dpi].y - P[i].Pos[1];
              dz = Mesh.DP[dpi].z - P[i].Pos[2];

              Neighbour_weights[n_neighbour].r = sqrt(dx * dx + dy * dy + dz * dz);

              Neighbour_weights[n_neighbour].x_plus[0] = fmax(dx, 0) / Neighbour_weights[n_neighbour].r;
              Neighbour_weights[n_neighbour].x_plus[1] = fmax(dy, 0) / Neighbour_weights[n_neighbour].r;
              Neighbour_weights[n_neighbour].x_plus[2] = fmax(dz, 0) / Neighbour_weights[n_neighbour].r;

              Neighbour_weights[n_neighbour].x_minus[0] = fmin(dx, 0) / Neighbour_weights[n_neighbour].r;
              Neighbour_weights[n_neighbour].x_minus[1] = fmin(dy, 0) / Neighbour_weights[n_neighbour].r;
              Neighbour_weights[n_neighbour].x_minus[2] = fmin(dz, 0) / Neighbour_weights[n_neighbour].r;

              Neighbour_weights[n_neighbour].omega =
                  0.5 * (1.0 - 1.0 / (sqrt(1.0 + (4.0 * Mesh.VF[vf].area /
                                                  (M_PI * Neighbour_weights[n_neighbour].r * Neighbour_weights[n_neighbour].r)))));

              n_neighbour++;

              if(q == SphP[i].last_connection)
                break;

              q = DC[q].next;
            }

          for(k = 0; k < 3; k++)
            {
              sum_omega_x_plus[k]  = 0;
              sum_omega_x_minus[k] = 0;
            }

          for(j = 0; j < n_neighbour; j++)
            {
              for(k = 0; k < 3; k++)
                {
                  sum_omega_x_plus[k] += Neighbour_weights[j].omega * fabs(Neighbour_weights[j].x_plus[k]);
                  sum_omega_x_minus[k] += Neighbour_weights[j].omega * fabs(Neighbour_weights[j].x_minus[k]);
                }
            }

          for(k = 0; k < 3; k++)
            {
              factor_plus[k] =
                  sqrt(0.5 * (1.0 + sum_omega_x_minus[k] * sum_omega_x_minus[k] / (sum_omega_x_plus[k] * sum_omega_x_plus[k])));
              factor_minus[k] =
                  sqrt(0.5 * (1.0 + sum_omega_x_plus[k] * sum_omega_x_plus[k] / (sum_omega_x_minus[k] * sum_omega_x_minus[k])));
            }

          for(j = 0; j < n_neighbour; j++)
            {
              for(k = 0; k < 3; k++)
                {
                  Neighbour_weights[j].w[k] =
                      Neighbour_weights[j].x_plus[k] * factor_plus[k] + Neighbour_weights[j].x_minus[k] * factor_minus[k];
                  Neighbour_weights[j].w[k] *= Neighbour_weights[j].omega;
                }
            }

          sum_weights = 0;

          for(j = 0; j < n_neighbour; j++)
            {
              sum_weights +=
                  sqrt(Neighbour_weights[j].w[0] * Neighbour_weights[j].w[0] + Neighbour_weights[j].w[1] * Neighbour_weights[j].w[1] +
                       Neighbour_weights[j].w[2] * Neighbour_weights[j].w[2]);
            }

          for(j = 0; j < n_neighbour; j++)
            {
              Neighbour_weights[j].w[0] *= (1.0 - All.HostCellFeedbackFraction) / sum_weights;
              Neighbour_weights[j].w[1] *= (1.0 - All.HostCellFeedbackFraction) / sum_weights;
              Neighbour_weights[j].w[2] *= (1.0 - All.HostCellFeedbackFraction) / sum_weights;
#ifdef SN_MCS_WEIGHTS_VERBOSE
              printf("SN_MCS_WEIGHTS_VERBOSE: %d %d %g %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n", ThisTask,
                     All.NumCurrentTiStep, All.Time, P[i].ID, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2], Neighbour_weights[j].r,
                     factor_plus[0], factor_plus[1], factor_plus[2], factor_minus[0], factor_minus[1], factor_minus[2],
                     Neighbour_weights[j].omega, Neighbour_weights[j].x_plus[0], Neighbour_weights[j].x_plus[1],
                     Neighbour_weights[j].x_plus[2], Neighbour_weights[j].x_minus[0], Neighbour_weights[j].x_minus[1],
                     Neighbour_weights[j].x_minus[2], sum_weights);
#endif
            }

          p_ej_tot = sqrt(2.0 * All.SNKineticRatio * SphP[i].energy_deposited * SphP[i].mass_deposited);
#ifdef MECHANICAL_FEEDBACK
          E_51 = (SphP[i].energy_deposited * All.SNKineticRatio) / (All.SupernovaEnergy * All.cf_atime * All.cf_atime);
#endif

          for(j = 0; j < n_neighbour; j++)
            {
              if(Neighbour_weights[j].Task == ThisTask)
                {
                  mod_weight = sqrt(Neighbour_weights[j].w[0] * Neighbour_weights[j].w[0] +
                                    Neighbour_weights[j].w[1] * Neighbour_weights[j].w[1] +
                                    Neighbour_weights[j].w[2] * Neighbour_weights[j].w[2]);

                  index = Neighbour_weights[j].Index;
                  if(P[index].ID == 0)
                    {
                      printf("ID0 NGB INJECT LOCAL: Task %d Sync-point %d index %d ID %d NumGas %d Mass %g\n", ThisTask,
                             All.NumCurrentTiStep, index, P[index].ID, NumGas, P[index].Mass);
                      continue;
                    }

                  u0 = (SphP[index].Energy -
                        0.5 *
                            (SphP[index].Momentum[0] * SphP[index].Momentum[0] + SphP[index].Momentum[1] * SphP[index].Momentum[1] +
                             SphP[index].Momentum[2] * SphP[index].Momentum[2]) /
                            P[index].Mass) /
                       P[index].Mass;
#ifdef MCS_UTHERM_CATCH
                  /*Debug code */
                  if(u0 < 0)
                    {
                      printf(
                          "UTHERM_CATCH_1: u0 < 0: Sync-point %d Time %g Task %d Energy %g Momentum %g %g %g Mass %g u0 %g index %d "
                          "ID %d\n",
                          All.NumCurrentTiStep, All.Time, ThisTask, SphP[index].Energy, SphP[index].Momentum[0],
                          SphP[index].Momentum[1], SphP[index].Momentum[2], P[index].Mass, u0, index, P[index].ID);
                      SphP[index].Energy += (-u0 + All.MinEgySpec * All.cf_atime * All.cf_atime);
                    }
#endif

                  dm = SphP[i].mass_deposited * mod_weight;

#ifndef SN_NO_ENERGY

                  SphP[index].Energy += mod_weight * SphP[i].energy_deposited;

                  dp[0] = Neighbour_weights[j].w[0] * p_ej_tot;
                  dp[1] = Neighbour_weights[j].w[1] * p_ej_tot;
                  dp[2] = Neighbour_weights[j].w[2] * p_ej_tot;

                  if(dm > 0) /*Guard against machine precision errors, occasionally occurs if omega ~ 0 becuase of distorted cell*/
                    {
                      dp_sq0 = dp[0] * dp[0] + dp[1] * dp[1] + dp[2] * dp[2];

                      /* Now move from star rest frame to coordinate frame */
                      dp[0] += dm * SphP[i].starvel[0];
                      dp[1] += dm * SphP[i].starvel[1];
                      dp[2] += dm * SphP[i].starvel[2];

                      SphP[index].Energy += (dp[0] * dp[0] + dp[1] * dp[1] + dp[2] * dp[2] - dp_sq0) / (2.0 * dm);

#ifdef MECHANICAL_FEEDBACK
                      n = SphP[index].Density * All.cf_a3inv * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;
                      n *= HYDROGEN_MASSFRAC / PROTONMASS;
                      Zsun       = fmax(SphP[index].Metallicity / 0.0127, 0.01);
                      p_terminal = All.SupernovaTerminalMomentum * pow(E_51, 16.0 / 17.0) * pow(n, -2.0 / 17.0) * pow(Zsun, -0.14) *
                                   All.cf_atime;
                      boost_factor = fmin(sqrt(1.0 + P[index].Mass / dm), (p_terminal / p_ej_tot));
                      if(boost_factor < 1.0)
                        {
                          terminate(
                              "SN_MCS MECHANICAL_FEEDBACK: boost_factor < 1, this shouldn't happen! p_terminal %g p_ej_tot %g "
                              "P[index].Mass %g dm %g n %g Zsun %g E_51 %g\n",
                              p_terminal, p_ej_tot / All.cf_atime, P[index].Mass, dm, n, Zsun, E_51);
                        }

                      dp[0] *= boost_factor;
                      dp[1] *= boost_factor;
                      dp[2] *= boost_factor;
#endif
                      SphP[index].Momentum[0] += dp[0];
                      SphP[index].Momentum[1] += dp[1];
                      SphP[index].Momentum[2] += dp[2];
                    }
#endif
                  P[index].Mass += dm;

#ifdef METALS
#ifdef METAL_ERROR_CATCH
                  /*Debug code*/
                  if((P[index].Metallicity < 0.0) || (SphP[index].Metallicity < 0.0) || (SphP[index].MassMetallicity < 0.0))
                    {
                      mpi_printf(
                          "METAL_ERROR_CATCH: sn_evaluate Task %d Time %g ID %d P[index].Mass %g P[index].Metallicity %g "
                          "SphP[index].Metallicity %g SphP[index].MassMetallicity %g Density %g Pos %g %g %g\n",
                          ThisTask, All.Time, P[index].ID, P[index].Mass, P[index].Metallicity, SphP[index].Metallicity,
                          SphP[index].MassMetallicity, SphP[index].Density, P[index].Pos[0], P[index].Pos[1], P[index].Pos[2]);
#ifdef MIN_METALLICITY_ON_STARTUP
                      P[index].Metallicity        = All.MinimumMetallicityOnStartUp;
                      SphP[index].Metallicity     = All.MinimumMetallicityOnStartUp;
                      SphP[index].MassMetallicity = P[index].Mass * All.MinimumMetallicityOnStartUp;
#else
                      P[index].Metallicity        = 0.0;
                      SphP[index].Metallicity     = 0.0;
                      SphP[index].MassMetallicity = 0.0;
#endif
                    }
#endif
                  SphP[index].MassMetallicity += dm * All.SNEjectaMetallicity;
                  SphP[index].Metallicity = SphP[index].MassMetallicity / P[index].Mass;
                  assert((SphP[index].Metallicity >= 0) && (SphP[index].Metallicity <= 1.0));
                  P[index].Metallicity = SphP[index].Metallicity;
#endif

#ifdef SN_MCS_CELL_VERBOSE
                  printf("SN_MCS_CELL_VERBOSE: %d %d %g %d %d %d %d %d %g %g %g %g %g %g %g %g %g %g %g %g\n", ThisTask,
                         All.NumCurrentTiStep, All.Time, P[i].ID, SphP[i].N_SN_hosted, n_neighbour, j + 1, Neighbour_weights[j].Task,
                         P[i].Pos[0], P[i].Pos[1], P[i].Pos[2], SphP[i].starvel[0], SphP[i].starvel[1], SphP[i].starvel[2], dm,
                         mod_weight * SphP[i].energy_deposited * (1.0 - All.SNKineticRatio), dp[0], dp[1], dp[2],
                         Neighbour_weights[j].r);
#endif
                  u1 = (SphP[index].Energy -
                        0.5 *
                            (SphP[index].Momentum[0] * SphP[index].Momentum[0] + SphP[index].Momentum[1] * SphP[index].Momentum[1] +
                             SphP[index].Momentum[2] * SphP[index].Momentum[2]) /
                            P[index].Mass) /
                       P[index].Mass;
                  if(u1 < u0)
                    SphP[index].Energy = u0 * P[index].Mass + 0.5 *
                                                                  (SphP[index].Momentum[0] * SphP[index].Momentum[0] +
                                                                   SphP[index].Momentum[1] * SphP[index].Momentum[1] +
                                                                   SphP[index].Momentum[2] * SphP[index].Momentum[2]) /
                                                                  P[index].Mass;
                }
              else
                {
                  mod_weight = sqrt(Neighbour_weights[j].w[0] * Neighbour_weights[j].w[0] +
                                    Neighbour_weights[j].w[1] * Neighbour_weights[j].w[1] +
                                    Neighbour_weights[j].w[2] * Neighbour_weights[j].w[2]);

                  dm = SphP[i].mass_deposited * mod_weight;

#ifndef SN_NO_ENERGY

                  SnvarIn[nexport].dE = mod_weight * SphP[i].energy_deposited;

                  dp[0] = Neighbour_weights[j].w[0] * p_ej_tot;
                  dp[1] = Neighbour_weights[j].w[1] * p_ej_tot;
                  dp[2] = Neighbour_weights[j].w[2] * p_ej_tot;

                  if(dm > 0) /*Guard against machine precision errors, occasionally occurs if omega ~ 0 becuase of distorted cell*/
                    {
                      dp_sq0 = dp[0] * dp[0] + dp[1] * dp[1] + dp[2] * dp[2];

                      /* Now move from star rest frame to coordinate frame */
                      dp[0] += dm * SphP[i].starvel[0];
                      dp[1] += dm * SphP[i].starvel[1];
                      dp[2] += dm * SphP[i].starvel[2];

                      SnvarIn[nexport].dE += (dp[0] * dp[0] + dp[1] * dp[1] + dp[2] * dp[2] - dp_sq0) / (2.0 * dm);
                    }

                  SnvarIn[nexport].dp[0] = dp[0];
                  SnvarIn[nexport].dp[1] = dp[1];
                  SnvarIn[nexport].dp[2] = dp[2];

#ifdef MECHANICAL_FEEDBACK
                  SnvarIn[nexport].p_ej_tot = p_ej_tot;
                  SnvarIn[nexport].E_51     = E_51;
#endif
#endif

                  SnvarIn[nexport].dm   = dm;
                  SnvarIn[nexport].Task = Neighbour_weights[j].Task;
                  Send_count[Neighbour_weights[j].Task] += 1;
                  SnvarIn[nexport].Index = Neighbour_weights[j].Index;
                  nexport++;
#ifdef SN_MCS_CELL_VERBOSE
                  printf("SN_MCS_CELL_VERBOSE: %d %d %g %d %d %d %d %d %g %g %g %g %g %g %g %g %g %g %g %g\n", ThisTask,
                         All.NumCurrentTiStep, All.Time, P[i].ID, SphP[i].N_SN_hosted, n_neighbour, j + 1, Neighbour_weights[j].Task,
                         P[i].Pos[0], P[i].Pos[1], P[i].Pos[2], SphP[i].starvel[0], SphP[i].starvel[1], SphP[i].starvel[2], dm,
                         mod_weight * SphP[i].energy_deposited * (1.0 - All.SNKineticRatio), dp[0], dp[1], dp[2],
                         Neighbour_weights[j].r);
#endif
                }
            }

          /* Give host cell energy, mass and metals */
          if(P[i].ID != 0)
            {
              P[i].Mass += SphP[i].mass_deposited * All.HostCellFeedbackFraction;
#ifndef SN_NO_ENERGY
              SphP[i].Energy += SphP[i].energy_deposited * All.HostCellFeedbackFraction;
              dv[0] = SphP[i].starvel[0] - P[i].Vel[0];
              dv[1] = SphP[i].starvel[1] - P[i].Vel[1];
              dv[2] = SphP[i].starvel[2] - P[i].Vel[2];
              kin_energy =
                  All.HostCellFeedbackFraction * 0.5 * SphP[i].mass_deposited * (dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]);
              SphP[i].Energy += kin_energy;
#endif

#ifdef METALS
#ifdef METAL_ERROR_CATCH
              /*Debug code*/
              if((P[i].Metallicity < 0.0) || (SphP[i].Metallicity < 0.0) || (SphP[i].MassMetallicity < 0.0))
                {
                  mpi_printf(
                      "METAL_ERROR_CATCH: sn_evaluate Task %d Time %g ID %d P[i].Mass %g P[i].Metallicity %g SphP[i].Metallicity %g "
                      "SphP[i].MassMetallicity %g Density %g Pos %g %g %g\n",
                      ThisTask, All.Time, P[i].ID, P[i].Mass, P[i].Metallicity, SphP[i].Metallicity, SphP[i].MassMetallicity,
                      SphP[i].Density, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
#ifdef MIN_METALLICITY_ON_STARTUP
                  P[i].Metallicity        = All.MinimumMetallicityOnStartUp;
                  SphP[i].Metallicity     = All.MinimumMetallicityOnStartUp;
                  SphP[i].MassMetallicity = P[i].Mass * All.MinimumMetallicityOnStartUp;
#else
                  P[i].Metallicity = 0.0;
                  SphP[i].Metallicity = 0.0;
                  SphP[i].MassMetallicity = 0.0;
#endif
                }
#endif
              SphP[i].MassMetallicity += SphP[i].mass_deposited * All.HostCellFeedbackFraction * All.SNEjectaMetallicity;
              SphP[i].Metallicity = SphP[i].MassMetallicity / P[i].Mass;
              assert((SphP[i].Metallicity >= 0) && (SphP[i].Metallicity <= 1.0));
              P[i].Metallicity = SphP[i].Metallicity;
#endif

#ifdef SN_MCS_CELL_VERBOSE
              printf("SN_MCS_CELL_VERBOSE: %d %d %g %d %d %d %d %d %g %g %g %g %g %g %g %g %g %g %g %g\n", ThisTask,
                     All.NumCurrentTiStep, All.Time, P[i].ID, SphP[i].N_SN_hosted, n_neighbour, 0, ThisTask, P[i].Pos[0], P[i].Pos[1],
                     P[i].Pos[2], SphP[i].starvel[0], SphP[i].starvel[1], SphP[i].starvel[2],
                     SphP[i].mass_deposited * All.HostCellFeedbackFraction, SphP[i].energy_deposited * All.HostCellFeedbackFraction,
                     0.0, 0.0, 0.0, 0.0);
#endif
            }
          else
            {
              /*Debug code*/
              printf("ID0 NGB HOSTINJECT LOCAL: Task %d Sync-point %d index %d ID %d NumGas %d Mass %g\n", ThisTask,
                     All.NumCurrentTiStep, i, P[i].ID, NumGas, P[i].Mass);
            }

          /* Reset host quantities */
          SphP[i].N_SN_hosted      = 0;
          SphP[i].mass_deposited   = 0;
          SphP[i].energy_deposited = 0;
          SphP[i].starvel[0]       = 0;
          SphP[i].starvel[1]       = 0;
          SphP[i].starvel[2]       = 0;
        }
    }

  myfree(Neighbour_weights);

  /* Get export buffer into correct order */
  qsort(SnvarIn, nexport, sizeof(struct snvar_in), snvar_compare);

  /* Send the export counts to other tasks and get receive counts */
  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(i = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; i < NTask; i++)
    {
      nimport += Recv_count[i];

      if(i > 0)
        {
          Send_offset[i] = Send_offset[i - 1] + Send_count[i - 1];
          Recv_offset[i] = Recv_offset[i - 1] + Recv_count[i - 1];
        }
    }

  SnvarGet = (struct snvar_in *)mymalloc("SnvarGet", nimport * sizeof(struct snvar_in));

  /* Exchange particle data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&SnvarIn[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct snvar_in), MPI_BYTE, recvTask,
                           TAG_DENS_A, &SnvarGet[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct snvar_in), MPI_BYTE,
                           recvTask, TAG_DENS_A, MPI_COMM_WORLD, &status);
            }
        }
    }

  /* Distribute imported quantitites */
  for(i = 0; i < nimport; i++)
    {
      assert(SnvarGet[i].Task == ThisTask); /* Check we have right task */
      index = SnvarGet[i].Index;
      if(index > NumGas)
        index -= NumGas;
      assert(P[index].Type == 0); /* Check we have a gas cell */
      if(P[index].ID == 0)
        {
          /*Debug code*/
          printf("ID0 NGB INJECT LOCAL: Task %d Sync-point %d index %d ID %d NumGas %d Mass %g\n", ThisTask, All.NumCurrentTiStep,
                 index, P[index].ID, NumGas, P[index].Mass);
          continue;
        }

      u0 = (SphP[index].Energy -
            0.5 *
                (SphP[index].Momentum[0] * SphP[index].Momentum[0] + SphP[index].Momentum[1] * SphP[index].Momentum[1] +
                 SphP[index].Momentum[2] * SphP[index].Momentum[2]) /
                P[index].Mass) /
           P[index].Mass;
#ifdef MCS_UTHERM_CATCH
      /*Debug code */
      if(u0 < 0)
        {
          printf("UTHERM_CATCH_2: u0 < 0: Sync-point %d Time %g Task %d Energy %g Momentum %g %g %g Mass %g u0 %g index %d ID %d\n",
                 All.NumCurrentTiStep, All.Time, ThisTask, SphP[index].Energy, SphP[index].Momentum[0], SphP[index].Momentum[1],
                 SphP[index].Momentum[2], P[index].Mass, u0, index, P[index].ID);
          SphP[index].Energy += (-u0 + All.MinEgySpec * All.cf_atime * All.cf_atime);
        }
#endif

#ifndef SN_NO_ENERGY
      SphP[index].Energy += SnvarGet[i].dE;

#ifdef MECHANICAL_FEEDBACK
      if(SnvarGet[i].dm > 0) /*Guard against machine precision errors, occasionally occurs if omega ~ 0 becuase of distorted cell*/
        {
          n = SphP[index].Density * All.cf_a3inv * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;
          n *= HYDROGEN_MASSFRAC / PROTONMASS;
          Zsun       = fmax(SphP[index].Metallicity / 0.0127, 0.01);
          p_terminal = All.SupernovaTerminalMomentum * pow(SnvarGet[i].E_51, 16.0 / 17.0) * pow(n, -2.0 / 17.0) * pow(Zsun, -0.14) *
                       All.cf_atime;
          boost_factor = fmin(sqrt(1.0 + P[index].Mass / SnvarGet[i].dm), (p_terminal / SnvarGet[i].p_ej_tot));
          if(boost_factor < 1.0)
            {
              terminate(
                  "SN_MCS MECHANICAL_FEEDBACK: boost_factor < 1, this shouldn't happen! p_terminal %g p_ej_tot %g P[index].Mass %g dm "
                  "%g n %g Zsun %g E_51 %g\n",
                  p_terminal, SnvarGet[i].p_ej_tot / All.cf_atime, P[index].Mass, SnvarGet[i].dm, n, Zsun, SnvarGet[i].E_51);
            }
          SnvarGet[i].dp[0] *= boost_factor;
          SnvarGet[i].dp[1] *= boost_factor;
          SnvarGet[i].dp[2] *= boost_factor;
        }
#endif
      SphP[index].Momentum[0] += SnvarGet[i].dp[0];
      SphP[index].Momentum[1] += SnvarGet[i].dp[1];
      SphP[index].Momentum[2] += SnvarGet[i].dp[2];

#endif

      P[index].Mass += SnvarGet[i].dm;

#ifdef METALS
#ifdef METAL_ERROR_CATCH
      /*Debug code*/
      if((P[index].Metallicity < 0.0) || (SphP[index].Metallicity < 0.0) || (SphP[index].MassMetallicity < 0.0))
        {
          mpi_printf(
              "METAL_ERROR_CATCH: sn_evaluate Task %d Time %g ID %d P[i].Mass %g P[i].Metallicity %g SphP[i].Metallicity %g "
              "SphP[i].MassMetallicity %g Density %g Pos %g %g %g\n",
              ThisTask, All.Time, P[i].ID, P[index].Mass, P[index].Metallicity, SphP[index].Metallicity, SphP[index].MassMetallicity,
              SphP[index].Density, P[index].Pos[0], P[index].Pos[1], P[index].Pos[2]);
#ifdef MIN_METALLICITY_ON_STARTUP
          P[index].Metallicity        = All.MinimumMetallicityOnStartUp;
          SphP[index].Metallicity     = All.MinimumMetallicityOnStartUp;
          SphP[index].MassMetallicity = P[index].Mass * All.MinimumMetallicityOnStartUp;
#else
          P[index].Metallicity = 0.0;
          SphP[index].Metallicity = 0.0;
          SphP[index].MassMetallicity = 0.0;
#endif
        }
#endif
      SphP[index].MassMetallicity += SnvarGet[i].dm * All.SNEjectaMetallicity;
      SphP[index].Metallicity = SphP[index].MassMetallicity / P[index].Mass;
      assert((SphP[index].Metallicity >= 0) && (SphP[index].Metallicity <= 1.0));
      P[index].Metallicity = SphP[index].Metallicity;
#endif

      u1 = (SphP[index].Energy -
            0.5 *
                (SphP[index].Momentum[0] * SphP[index].Momentum[0] + SphP[index].Momentum[1] * SphP[index].Momentum[1] +
                 SphP[index].Momentum[2] * SphP[index].Momentum[2]) /
                P[index].Mass) /
           P[index].Mass;
      if(u1 < u0)
        SphP[index].Energy = u0 * P[index].Mass + 0.5 *
                                                      (SphP[index].Momentum[0] * SphP[index].Momentum[0] +
                                                       SphP[index].Momentum[1] * SphP[index].Momentum[1] +
                                                       SphP[index].Momentum[2] * SphP[index].Momentum[2]) /
                                                      P[index].Mass;
    }

  myfree(SnvarGet);
  myfree(SnvarIn);
  update_primitive_variables();
}

/* Sorts by task, then by index */
static int snvar_compare(const void *a, const void *b)
{
  if(((struct snvar_in *)a)->Task < (((struct snvar_in *)b)->Task))
    return -1;

  if(((struct snvar_in *)a)->Task > (((struct snvar_in *)b)->Task))
    return +1;

  if(((struct snvar_in *)a)->Index < (((struct snvar_in *)b)->Index))
    return -1;

  if(((struct snvar_in *)a)->Index > (((struct snvar_in *)b)->Index))
    return +1;

  return 0;
}

#endif

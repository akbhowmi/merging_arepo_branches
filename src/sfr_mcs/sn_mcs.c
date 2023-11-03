/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/sfr_mcs/sn_mcs.c
 * \date        08/2018
 * \author     	Matthew C Smith
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

#include <gsl/gsl_randist.h>

#if !defined(METALS) && !defined(SB99_FIXED_Z)
#error "If you are not using METALS, then Starburst99 requires a metallicity, use SB99_FIXED_Z to achieve this"
#endif

void check_for_supernovae(int *local_n_sn_event, int *global_n_sn_event)
{
  int i, idx;
  int n_sn_event     = 0;
  int n_sn           = 0;
  int n_sn_tot       = 0;
  int n_sn_event_tot = 0;
  int n_sn_i, n_sn_init;
  double dt, dtime, t_star, raw_sn_rate, norm_sn_rate, n_sn_expect, snr_global;
  double m_ej;

  mpi_printf("SN_MCS: Checking for SNe\n");

  /* To Do: Loop over star particle struct instead of all active particles */
  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Type != 4)
        continue;

      if(P[i].ID == 0 && P[i].Mass == 0)
        continue;

      if(All.ComovingIntegrationOn && (P[i].Mass > All.MaxSNStarMass))
        continue;

      dt    = (P[i].TimeBinGrav ? (((integertime)1) << P[i].TimeBinGrav) : 0) * All.Timebase_interval;
      dtime = dt / All.cf_hubble_a;

      t_star = get_time_difference_in_Gyr(P[i].StellarAge, All.Time) * 1.0e09;
      if((t_star <= sb99.t_min) || (t_star >= sb99.t_max)) /* No SNe go off outside of this time range */
        continue;

#ifndef SB99_FIXED_Z
      raw_sn_rate = get_sn_rate(t_star, P[i].Metallicity); /* returns log10(SN Msun^-1 yr^-1) */
#else
      raw_sn_rate = get_sn_rate(t_star); /* returns log10(SN Msun^-1 yr^-1) */
#endif
      norm_sn_rate = pow(10.0, raw_sn_rate) *
                     (All.UnitTime_in_s * All.UnitMass_in_g / (All.HubbleParam * All.HubbleParam * SEC_PER_YEAR * SOLAR_MASS)) *
                     P[i].InitialMass;

      n_sn_expect = norm_sn_rate * dtime;
      n_sn_i      = (int)gsl_ran_poisson(random_generator, n_sn_expect);
      n_sn_init   = n_sn_i;

      if(n_sn_i > 0)
        {
          m_ej = All.SNMassReturn * n_sn_i;
          if((P[i].Mass - m_ej) < 0.1 * P[i].InitialMass)
            {
              P[i].M_ej = P[i].Mass;
              P[i].Mass = 0.0;
            }
          else
            {
              P[i].M_ej = m_ej;
              P[i].Mass -= m_ej;
            }
          P[i].N_SN = n_sn_i;
          P[i].N_SN_cum += n_sn_i;

          assert(P[i].M_ej >= 0);
          assert(P[i].Mass >= 0);

#ifdef SN_MCS_DEBUG
#ifndef METALS
          printf(
              "Task %d, Time %g, dtime = %g, ID = %d, t_star %g, raw_sn_rate = %g, InitialMass = %g, mass = %g, norm_sn_rate = %g, "
              "n_sn_expect = %g, n_sn_init = %d, n_sn_i = %d, N_SN = %d, mass returned = %g\n",
              ThisTask, All.Time, dtime, P[i].ID, t_star, raw_sn_rate, P[i].InitialMass, P[i].Mass, norm_sn_rate, n_sn_expect,
              n_sn_init, n_sn_i, P[i].N_SN, P[i].M_ej);
#else
          printf("SN_MCS_DEBUG %d %g %g %d %g %g %g %g %g %g %g %d %d %d %d %g\n", ThisTask, All.Time, dtime, P[i].ID, t_star,
                 raw_sn_rate, P[i].InitialMass, P[i].Mass, P[i].Metallicity, norm_sn_rate, n_sn_expect, n_sn_init, n_sn_i, P[i].N_SN,
                 P[i].N_SN_cum, P[i].M_ej);
#endif
#endif
          n_sn_event++;
          P[i].N_SN_event_cum += 1;
          n_sn += n_sn_i;
        }
      else
        {
          P[i].N_SN = 0;
          P[i].M_ej = 0;
        }

    } /* Main loop of active particles */

  MPI_Reduce(&n_sn, &n_sn_tot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD); /*This is used for log, so only Task 0 needs this */
  MPI_Allreduce(&n_sn_event, &n_sn_event_tot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD); /*Every task needs this to trigger SN routines */
  if(ThisTask == 0)
    {
      if(All.TimeStep > 0)
        snr_global = n_sn_tot / (All.TimeStep / All.cf_time_hubble_a);
      else
        snr_global = 0.0;
      fprintf(FdSnr, "%14e %14i %14i %14e\n", All.Time, n_sn_tot, n_sn_event_tot, snr_global);
      myflush(FdSnr);
    }

  *local_n_sn_event  = n_sn_event;
  *global_n_sn_event = n_sn_event_tot;
}

/*Tidy up, reset things for next timestep, largely for safety */
void sn_finish(void)
{
  int i, idx;

  /* To do: loop over star particle struct */
  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;
      if(P[i].Type != 4)
        continue;
      if(P[i].ID == 0 && P[i].Mass == 0)
        continue;

      /* Prepare 'dead' star particles for deletion */
      if(P[i].Mass < (0.1 * P[i].InitialMass))
        {
          P[i].ID     = 0;
          P[i].Mass   = 0;
          P[i].Vel[0] = 0;
          P[i].Vel[1] = 0;
          P[i].Vel[2] = 0;
        }

      P[i].N_SN = 0; /* Redundant, these should be taken care of elsewhere, but kept just in case */
      P[i].M_ej = 0.0;
    }
}

void do_supernovae(void)
{
  int n_sn_event, n_sn_event_global;

  check_for_supernovae(&n_sn_event, &n_sn_event_global);

  if(n_sn_event_global > 0) /*i.e. some task somewhere has at least one SN */
    {
      find_sn_host_cells_and_distribute(n_sn_event);
      stellar_feedback_distribute();
    }

  sn_finish(); /*Tidy up */
}

#endif

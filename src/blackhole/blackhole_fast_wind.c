

#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../allvars.h"
#include "../proto.h"

/*! \file blackhole_feedback.c
 *  \brief distribute feedback energy
 */

#ifdef BLACK_HOLES
#ifdef BH_FAST_WIND

void blackhole_update_wind_affected_cells(void)
{ 
  int idx, i, j;
  
  mpi_printf("BH_FAST_WIND: Updating momenta of affected cells.\n");
  
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    { 
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;
      
      if(P[i].Mass == 0 && P[i].ID == 0)
        continue;
      
      if(P[i].Type == 0)
        { 
          SphP[i].Energy -= 0.5 * P[i].Mass * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);
          P[i].Mass += SphP[i].Injected_BH_Wind_Mass;
          
          for(j = 0; j < 3; j++) 
            { 
              SphP[i].Momentum[j] += SphP[i].Injected_BH_Wind_Momentum[j];
              P[i].Vel[j] = SphP[i].Momentum[j] / P[i].Mass;
              
              SphP[i].Injected_BH_Wind_Momentum[j] = 0;
            }
          
          SphP[i].Energy += 0.5 * P[i].Mass * (P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);
          if(SphP[i].Injected_BH_Wind_Mass > 0)
            { 
              SphP[i].MassMetallicity = SphP[i].Metallicity * P[i].Mass;
              for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
                SphP[i].MassMetals[j] = SphP[i].MetalsFraction[j] * P[i].Mass;
              
              printf("CellInjectionData: Time = %.10g, ID = %llu, Mass = %g, Position = (%g|%g|%g), Velocity = (%g|%g|%g), Mom = (%g|%g|%g), Energy = %g, Utherm = %g, Ne = %g, Density = %g, Volume = %g, Metallicity = %g, MassMetallicity = %g , Inj_mass = %g\n",All.Time, P[i].ID, P[i].Mass, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2], P[i].Vel[0], P[i].Vel[1], P[i].Vel[2], SphP[i].Momentum[0], SphP[i].Momentum[1], SphP[i].Momentum[2], SphP[i].Energy, SphP[i].Utherm, SphP[i].Ne, SphP[i].Density, SphP[i].Volume, SphP[i].Metallicity, SphP[i].MassMetallicity, SphP[i].Injected_BH_Wind_Mass);
              myflush(stdout);
            }
          SphP[i].Injected_BH_Wind_Mass = 0;
        }
    }
}



static int blackhole_evaluate(int target, int mode, int threadid);

/* local data structure for collecting particle/cell data that is sent to other processors if needed */
typedef struct
{
  MyDouble Pos[3];
  MyFloat BH_Hsml;
  MyFloat BH_WindEnergy;
  MyFloat BH_Density;
  unsigned int bh_pid;
#ifdef BH_FAST_WIND_STOCHASTIC
  //int totalngbnum;
  MyFloat totalngbmass;
  //MyFloat totalngbweight;
#else
  int randomtask;
#endif

  int Firstnode;
} data_in;

static data_in *DataGet;

/* routine that fills the relevant particle/cell data into the input structure defined above */
static void particle2in(data_in *in, int i, int firstnode)
{
  for(int k = 0; k < 3; k++)
    in->Pos[k] = P[i].Pos[k];

  in->BH_Hsml       = BPP(i).BH_Hsml;
  in->BH_WindEnergy = BPP(i).BH_WindEnergy;
  in->BH_Density    = BPP(i).BH_Density;
  //in->bh_pid        = BPP(i).PID;
  in->bh_pid        = P[i].ID;
#ifdef BH_FAST_WIND_STOCHASTIC
  //in->totalngbnum   = BPP(i).Total_Ngb_Num;
  in->totalngbmass  = BPP(i).Total_Ngb_Mass;
  //in->totalngbweight = BPP(i).Total_Ngb_Weight;
#else
  in->randomtask    = BPP(i).Random_Task;
#endif

  in->Firstnode = firstnode;
}

/* local data structure that holds results acquired on remote processors */
typedef struct
{
  int dummy;
} data_out;

static data_out *DataResult;

/* routine to store or combine result data */
static void out2particle(data_out *out, int i, int mode) {}

#include "../generic_comm_helpers2.h"

/*static int blackhole_wind_test_for_injection(int idx)
{
#if(defined(GFM_WINDS_VARIABLE) && (GFM_WINDS_VARIABLE == 1)) || defined(GFM_WINDS_LOCAL)
  double egycmp = All.RadioFeedbackReiorientationFactor * 0.5 * pow(BPP(idx).BH_DMVelDisp, 2) * BPP(idx).BH_Density *
                  (4.0 * M_PI / 3.0) * pow(BPP(idx).BH_Hsml, 3);
#else
  double egycmp =
      All.RadioFeedbackReiorientationFactor * BPP(idx).BH_U * BPP(idx).BH_Density * (4.0 * M_PI / 3.0) * pow(BPP(idx).BH_Hsml, 3);
#endif

  if(BPP(idx).BH_WindEnergy >= egycmp)
    return 1;
  else
    return 0;
}*/

static void kernel_local(void)
{
  int i;
#ifdef GENERIC_ASYNC
  int flag = 0;
#endif
#pragma omp parallel private(i)
  {
    int j, threadid = get_thread_num();
#ifdef GENERIC_ASYNC
    int count = 0;
#endif
    for(j = 0; j < NTask; j++)
      Thread[threadid].Exportflag[j] = -1;

    while(1)
      {
        if(Thread[threadid].ExportSpace < MinSpace)
          break;

#ifdef GENERIC_ASYNC
        if(threadid == 0)
          {
            if((count & POLLINGINTERVAL) == 0)
              if(generic_polling_primary(count, TimeBinsBHAccretion.NActiveParticles))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        i = NextParticle++;

        if(i >= TimeBinsBHAccretion.NActiveParticles)
          break;

        int idx = TimeBinsBHAccretion.ActiveParticleList[i];

        if(idx < 0)
          continue;

        if(BPP(idx).BH_WindEnergy == 0)
          continue;

        blackhole_evaluate(idx, MODE_LOCAL_PARTICLES, threadid);
      }
  }
}


static void kernel_imported(void)
{
  /* now do the particles that were sent to us */
  int i, cnt = 0;
#pragma omp parallel private(i)
  {
    int threadid = get_thread_num();
#ifdef GENERIC_ASYNC
    int count = 0;
#endif

    while(1)
      {
#pragma omp atomic capture
        i = cnt++;

        if(i >= Nimport)
          break;

#ifdef GENERIC_ASYNC
        if(threadid == 0)
          {
            if((count & POLLINGINTERVAL) == 0)
              generic_polling_secondary();
          }

        count++;
#endif

        blackhole_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}



/* This function distributes the feedback energy and look for gas cells that we might stochastically swallow
 */


void blackhole_blow_wind(void)
{
  mpi_printf("BH_FAST_WIND: Start assigning BH momentum feedback\n");
  double t0 = second();

  generic_set_MaxNexport();


  generic_comm_pattern(TimeBinsBHAccretion.NActiveParticles, kernel_local, kernel_imported);


  for(int i = 0; i < TimeBinsBHAccretion.NActiveParticles; i++)
  {
    int idx = TimeBinsBHAccretion.ActiveParticleList[i];

    if (idx < 0)
      continue;

    BPP(idx).BH_WindEnergy = 0;
  }


  double t1 = second();
  mpi_printf("BH_FAST_WIND: Done assigning BH momentum feedback\n", timediff(t0, t1));
}

static int blackhole_evaluate(int target, int mode, int threadid)
{
  int numnodes, *firstnode;

  data_in local, *in;
  data_out out;

  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      in = &local;

      numnodes  = 1;
      firstnode = NULL;
    }
  else
    {
      in = &DataGet[target];

      generic_get_numnodes(target, &numnodes, &firstnode);
    }

  MyDouble *pos     = in->Pos;
  double h_i        = in->BH_Hsml;
  double windenergy = in->BH_WindEnergy;
  double dens       = in->BH_Density;
#ifndef BH_FAST_WIND_STOCHASTIC
  int BH_randomtask = in->randomtask;
#endif
  double hinv  = 1.0 / h_i;
  double hinv3 = hinv * hinv * hinv;

  int nfound = ngb_treefind_variable_threads(pos, h_i, target, mode, threadid, numnodes, firstnode);

  int temp_flag=0;
  int count = 0;

  if (windenergy ==0) {
    printf("FastWind Injection1: Time = %.10g, BH_PID = %llu, Task = %d, BH_WindEnergy = 0!!\n", All.Time, in->bh_pid, ThisTask);
    myflush(stdout);
  }
#ifdef BH_FAST_WIND_STOCHASTIC
  //int total_ngb_num   = in->totalngbnum;
  //double total_ngb_weight = in->totalngbweight;
  double total_ngb_mass  = in-> totalngbmass;
  
  for(int n = 0; n < nfound; n++)
    {
      int j = Thread[threadid].Ngblist[n];

      if(P[j].Mass > 0 && P[j].ID != 0)
        {
          double r2 = Thread[threadid].R2list[n];
          double r  = sqrt(r2);
          double u  = r * hinv;
          double dx = P[j].Pos[0] - pos[0];
          double dy = P[j].Pos[1] - pos[1];
          double dz = P[j].Pos[2] - pos[2];
          double nx, ny, nz;

          /*if(u < 0.5)
            wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
          else
            wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);*/

          double wind_vel = All.FastWindVelocity * CLIGHT / All.UnitVelocity_in_cm_per_s;
          double totenergy  = 0.5 * total_ngb_mass * wind_vel * wind_vel;
          double prob = windenergy / totenergy;
          //prob *= (wk / total_ngb_weight) * total_ngb_num; 
          
          if(All.FastWindKickDirection == 0)
	   {
	     double theta = acos(2 * get_random_number() - 1);
	     double phi   = 2 * M_PI * get_random_number();
	     nx    = sin(theta) * cos(phi);
	     ny    = sin(theta) * sin(phi);
	     nz    = cos(theta);
	   }
	   else if(All.FastWindKickDirection == 1)
	   {
	     nx = dx / r;
	     ny = dy / r;
	     nz = dz / r;
	   }
	   if(get_random_number() < prob && u < 1.0)
	   {
             SphP[j].Injected_BH_Wind_Momentum[0] += P[j].Mass * wind_vel * nx;
             SphP[j].Injected_BH_Wind_Momentum[1] += P[j].Mass * wind_vel * nx;
             SphP[j].Injected_BH_Wind_Momentum[2] += P[j].Mass * wind_vel * nx;
             
             printf("FastWindStochasticInjection: Time = %.10g, BH_PID = %llu, NTask = %d, ThisTask = %d, P[j].ID = %lld, P[j].Mass = %g, SphP[j].Density = %g, SphP[j].Energy = %g, SphP[j].Utherm = %g, SphP[j].Mom = (%g|%g|%g), P[j].Vel = (%g|%g|%g), SphP[j].Metallicity = %g, SphP[j].MassMetallicity = %g,  wind_energy = %g, nfound = %d, nx = %g, ny = %g, nz = %g, dx = %g, dy = %g, dz = %g, probanility = %g, \n", All.Time, in->bh_pid, NTask, ThisTask, P[j].ID, P[j].Mass, SphP[j].Density, SphP[j].Energy, SphP[j].Utherm, SphP[j].Momentum[0], SphP[j].Momentum[1], SphP[j].Momentum[2], P[j].Vel[0], P[j].Vel[1], P[j].Vel[2], SphP[j].Metallicity, SphP[j].MassMetallicity, windenergy, nfound, nx, ny, nz, dx, dy, dz, prob);
             myflush(stdout);
             fprintf(FdBlackHolesFWstoch, "BH=%llu %.10g %g %g %g %g %g %g %g %g %g %g\n", (long long)in->bh_pid, All.Time, windenergy, totenergy, prob, P[j].Mass, SphP[j].Density, nx, ny, nz, r, h_i);
           }
           //fprintf(FdBlackHolesFWstoch, "BH=%llu %.10g %g %g %g %g %d %g %g %g %g %g %g %g %g %g\n", (long long)in->bh_pid, All.Time, windenergy, totenergy, wk, total_ngb_weight, total_ngb_num, prob, randnum, P[j].Mass, SphP[j].Density, nx, ny, nz, r, h_i);
        }
    }
#else
  while (ThisTask == BH_randomtask && temp_flag == 0 && windenergy > 0)
  {
    if(nfound == 0)
      terminate("FastWind terminate: Time = %.10g. nfound is zero for BH_randomtask = %d. BH_PID = %llu", All.Time, BH_randomtask, in->bh_pid);

    int random_index = nfound * get_random_number();
    int j = Thread[threadid].Ngblist[random_index];
    if (P[j].Mass > 0 && P[j].ID != 0) {
      double dx = P[j].Pos[0] - pos[0];
      double dy = P[j].Pos[1] - pos[1];
      double dz = P[j].Pos[2] - pos[2];
      double r2 = dx * dx + dy * dy + dz * dz;
      double r  = sqrt(r2);

      double wind_vel = All.FastWindVelocity * CLIGHT / All.UnitVelocity_in_cm_per_s;
      double wind_mom = 2 * windenergy / wind_vel;
      double wind_mass = wind_mom / wind_vel;
      double nx, ny, nz;

      if(All.FastWindKickDirection == 0)
      {
        double theta = acos(2 * get_random_number() - 1);
        double phi   = 2 * M_PI * get_random_number();
        nx    = sin(theta) * cos(phi);
        ny    = sin(theta) * sin(phi);
        nz    = cos(theta);
      }
      else if(All.FastWindKickDirection == 1)
      {
        nx = dx / r;
        ny = dy / r;
        nz = dz / r;
      }
      SphP[j].Injected_BH_Wind_Mass += wind_mass;
      SphP[j].Injected_BH_Wind_Momentum[0] += wind_mom * nx;
      SphP[j].Injected_BH_Wind_Momentum[1] += wind_mom * ny;
      SphP[j].Injected_BH_Wind_Momentum[2] += wind_mom * nz;
      temp_flag++;
      printf("FastWind MomentumInjection: Time = %.10g, BH_PID = %llu, NTask = %d, ThisTask = %d, P[j].ID = %lld, P[j].Mass = %g, SphP[j].Density = %g, SphP[j].Energy = %g, SphP[j].Utherm = %g, SphP[j].Mom = (%g|%g|%g), P[j].Vel = (%g|%g|%g), SphP[j].Metallicity = %g, SphP[j].MassMetallicity = %g, Inj_mass = %g, wind_mom = %g, nfound = %d, random index = %d, nx = %g, ny = %g, nz = %g, r = %g \n", All.Time, in->bh_pid, NTask, ThisTask, P[j].ID, P[j].Mass, SphP[j].Density, SphP[j].Energy, SphP[j].Utherm, SphP[j].Momentum[0], SphP[j].Momentum[1], SphP[j].Momentum[2], P[j].Vel[0], P[j].Vel[1], P[j].Vel[2], SphP[j].Metallicity, SphP[j].MassMetallicity, wind_mass, wind_mom, nfound, random_index, nx, ny, nz, r);
      myflush(stdout);
    }
    else
    {
      count++;
      printf("FastWind Injection failed, All.Time = %.10g, BH_PID = %llu, ThisTask = %d, nfound = %d, random_index = %d, P[j].Mass = %g, P[j].ID = %lld \n", All.Time, in->bh_pid, ThisTask, nfound, random_index, P[j].Mass, P[j].ID);
      myflush(stdout);
      if(count == 1000)
        terminate("FW_terminate: Time = %.10g, BH_PID = %llu, ThisTask = %d, nfound = %d", All.Time, in->bh_pid, ThisTask, nfound);
    }
  }
#endif

  out.dummy = 0;

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}


#endif
#endif

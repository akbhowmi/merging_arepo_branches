

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

static int blackhole_find_tasknum(int target, int mode, int threadid);

/* local data structure for collecting particle/cell data that is sent to other processors if needed */
typedef struct
{
  MyDouble Pos[3];
  MyFloat BH_Hsml;
  unsigned int bh_pid;
  int Firstnode;
} data_in;

static data_in *DataGet;

/* routine that fills the relevant particle/cell data into the input structure defined above */
static void particle2in(data_in *in, int i, int firstnode)
{
  for(int k = 0; k < 3; k++)
    in->Pos[k] = P[i].Pos[k];

  in->BH_Hsml       = BPP(i).BH_Hsml;
  //in->BH_WindEnergy = BPP(i).BH_WindEnergy;
  //in->BH_Density    = BPP(i).BH_Density;
  in->bh_pid        = P[i].ID;
  //for(int k = 0; k < 3; k++)
    //in->WindDir[k] = BPP(i).WindDir[k];

  in->Firstnode = firstnode;
}

/* local data structure that holds results acquired on remote processors */
typedef struct
{
  int tasknumber;
  
#ifdef BH_FAST_WIND_STOCHASTIC
  MyFloat total_ngb_mass;
#endif

} data_out;

static data_out *DataResult;

/* routine to store or combine result data */
static void out2particle(data_out *out, int i, int mode)
{
#ifdef BH_FAST_WIND_STOCHASTIC
  if(mode == MODE_LOCAL_PARTICLES) /* initial store */
  {
    BPP(i).Total_Ngb_Mass = out->total_ngb_mass;
  }
  else
  {
    BPP(i).Total_Ngb_Mass += out->total_ngb_mass;
  }
#else
  if (out->tasknumber != -1)
  {
    BPP(i).BHTaskList[out->tasknumber] = 1;
  }
#endif
}

#include "../generic_comm_helpers2.h"



static void kernel_local_findtasks(void)
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

        blackhole_find_tasknum(idx, MODE_LOCAL_PARTICLES, threadid);
      }
  }

}


static void kernel_imported_findtasks(void)
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

        blackhole_find_tasknum(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}


/* This function distributes the feedback energy and look for gas cells that we might stochastically swallow
 */

void blackhole_broadcast_tasknum(void)
{
  //int *BHTaskNums = (int *)mymalloc("BHTaskNums", NTask * sizeof(int));
  //memset(BHTaskNums, 0, NTask * sizeof(int));

#ifndef BH_FAST_WIND_STOCHASTIC
  for(int i = 0; i < TimeBinsBHAccretion.NActiveParticles; i++)
  {
    int idx = TimeBinsBHAccretion.ActiveParticleList[i];

    if (idx < 0)
      continue;

    BPP(idx).Random_Task = -1;
    for(int j = 0; j < 1024; j++)
      BPP(idx).BHTaskList[j] = 0;
    printf("FastWind Task1: All.Time = %.10g, idx = %d, P[idx].ID = %llu, ThisTask = %d. BHTaskList reset to zero. random_task set to -1. \n", All.Time, idx, P[idx].ID, ThisTask);
    myflush(stdout);
  }
#endif

  generic_set_MaxNexport();

  generic_comm_pattern(TimeBinsBHAccretion.NActiveParticles, kernel_local_findtasks, kernel_imported_findtasks);

#ifndef BH_FAST_WIND_STOCHASTIC
  for(int i = 0; i < TimeBinsBHAccretion.NActiveParticles; i++)
  {
    int idx = TimeBinsBHAccretion.ActiveParticleList[i];

    if (idx < 0 || BPP(idx).BH_WindEnergy == 0)
      continue;

    int sum = 0;
    for(int j = 0; j < NTask; j++)
      sum += BPP(idx).BHTaskList[j];
    if(sum == 0)
    {
      terminate("FastWind terminating All.Time = %.10g, ThisTask = %d, idx = %d, P[idx].ID = %llu  All elements of the array BHTaskList are 0. No tasks with BH neighbors found!! \n", All.Time, ThisTask, idx, P[idx].ID);
      //return;
    }
  }


  for(int i = 0; i < TimeBinsBHAccretion.NActiveParticles; i++)
  {
    int idx = TimeBinsBHAccretion.ActiveParticleList[i];

    if (idx < 0 || BPP(idx).BH_WindEnergy == 0)
      continue;

    printf("FastWind Task2: Time = %.10g, Tasks with Ngbs for P[idx].ID = %lld are : ", All.Time, P[idx].ID);
    for(int i = 0; i < NTask; i++)
    {
      if(BPP(idx).BHTaskList[i] == 1)
        printf("%d ", i);
    }
    printf("\n");
    myflush(stdout);


    int j = NTask * get_random_number();

    while(BPP(idx).BHTaskList[j] == 0)
    {
      j = NTask * get_random_number();
    }
    BPP(idx).Random_Task = j;
    printf("FastWind Task3: All.Time = %.10g, idx = %d, P[idx].ID = %llu, ThisTask = %d, random_task set to %d \n", All.Time, idx, P[idx].ID, ThisTask, BPP(idx).Random_Task);
    myflush(stdout);
  }
#endif

}



static int blackhole_find_tasknum(int target, int mode, int threadid)
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
  double hinv, hinv3, u, wk;

  int nfound = ngb_treefind_variable_threads(pos, h_i, target, mode, threadid, numnodes, firstnode);
  int sum = 0;
#ifdef BH_FAST_WIND_STOCHASTIC
  MyFloat totmass = 0;
#endif
  if(nfound != 0)
  {
    for(int i = 0; i < nfound; i++)
    {
      int j = Thread[threadid].Ngblist[i];
      if(P[j].Mass > 0 && P[j].ID != 0)
      {
        sum++;
#ifdef BH_FAST_WIND_STOCHASTIC
        totmass += P[j].Mass;
#endif       
      }
    }
    if(sum > 0)
    {
      out.tasknumber = ThisTask;
      printf("FastWind Task4 BHTaskListUpdate: Time = %.10g, BH_PID = %llu, ThisTask = %d, nfound = %d \n", All.Time, in->bh_pid, ThisTask, nfound);
      myflush(stdout);
    }
    else
      out.tasknumber = -1;
  }
  else
    out.tasknumber = -1;

#ifdef BH_FAST_WIND_STOCHASTIC
  out.total_ngb_mass = totmass;
#endif

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

#endif
#endif

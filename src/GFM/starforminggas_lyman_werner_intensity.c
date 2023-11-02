/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/GFM/agn_radiation.c
 * \date        MM/YYYY
 * \author
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../allvars.h"
#include "../proto.h"

#if defined(BLACK_HOLES) && defined(CALCULATE_LYMAN_WERNER_INTENSITY_LOCAL_STARFORMINGGAS)

static int assign_starforminggas_lyman_werner_intensity_evaluate(int target, int mode, int threadid);
static int reset_starforminggas_lyman_werner_intensity_evaluate(int target, int mode, int threadid);
static void kernel_local_reset(void);
static void kernel_imported_reset(void);

typedef struct
{
  MyDouble Pos[3];
  MyFloat Mass; 
  int Firstnode;
  MyFloat hsml , metallicity , Sfr;
} data_in;

static data_in *DataGet;

/* routine that fills the relevant particle/cell data into the input structure defined above */
static void particle2in(data_in *in, int i, int firstnode)
{
  in->Pos[0] = P[i].Pos[0];
  in->Pos[1] = P[i].Pos[1];
  in->Pos[2] = P[i].Pos[2];
  in->Mass = P[i].Mass;
  in->hsml = SphP[i].Hsml;
  in->Firstnode = firstnode;
  in->metallicity = SphP[i].Metallicity;
  in->Sfr = SphP[i].Sfr;
}

/* local data structure that holds results acquired on remote processors */
typedef struct
{
  char dummy;
} data_out;

static data_out *DataResult;

/* routine to store or combine result data */
static void out2particle(data_out *out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES) /* initial store */
    {
    }
  else /* merge */
    {
    }
}

#include "../generic_comm_helpers2.h"

static void kernel_local(void)
{
  int idx;
#ifdef GENERIC_ASYNC
  int flag = 0;
#endif

#pragma omp parallel private(i, idx)
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
        idx = NextParticle++;

        if(idx >= NumPart)
          break;

        int i = idx;
        if((i < 0) || ((P[i].Type != 0) || (SphP[i].Sfr == 0)))
          continue;
        assign_starforminggas_lyman_werner_intensity_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
      }
  }
}

static void kernel_local_reset(void)
{
  int idx;
#ifdef GENERIC_ASYNC
  int flag = 0;
#endif

#pragma omp parallel private(i, idx)
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
        idx = NextParticle++;

        if(idx >= NumPart)
          break;

        int i = idx;
        if((i < 0) || ((P[i].Type != 0) || (SphP[i].Sfr == 0)))
          continue;
        reset_starforminggas_lyman_werner_intensity_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
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

        assign_starforminggas_lyman_werner_intensity_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

static void kernel_imported_reset(void)
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

        reset_starforminggas_lyman_werner_intensity_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

void assign_starforminggas_lyman_werner_intensity(void)
{
  int i;
  for(i = 0; i < NumGas; i++)
   {
     SphP[i].StarFormingGasLymanWernerIntensity_type2 = 0;
     SphP[i].StarFormingGasLymanWernerIntensity_type3 = 0;
   }
  generic_set_MaxNexport();
  generic_comm_pattern(NumPart, kernel_local, kernel_imported);
}

void reset_starforminggas_lyman_werner_intensity(void)
{
  generic_set_MaxNexport();
  generic_comm_pattern(NumPart, kernel_local_reset, kernel_imported_reset);
}




int assign_starforminggas_lyman_werner_intensity_evaluate(int target, int mode, int threadid)
{
  int j, n, numnodes, *firstnode, find_neighbors = 0, population_type;
  double dx, dy, dz, r2, h, h2, hinc, hinc2, wp , h_inv, u;
  MyDouble *pos;
  MyFloat mass, metallicity, sfr;
  double starforminggas_lum_lyman_werner, local_intensity = 0;
#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
#endif

  data_in local, *target_data;

  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      target_data = &local;

      numnodes  = 1;
      firstnode = NULL;
    }
  else
    {
      target_data = &DataGet[target];

      generic_get_numnodes(target, &numnodes, &firstnode);
    }

  sfr   = target_data->Sfr;
  pos    = target_data->Pos;
  h      = target_data->hsml;
  h2     = h * h;
  hinc   = All.LymanWernerRadiusOfInclusion * h;
  metallicity = target_data->metallicity;
  hinc2 = hinc * hinc;
  mass = sfr * 5e6 * 1e-10 * All.HubbleParam; 

  if (metallicity < 0.001 * GFM_SOLAR_METALLICITY) /* This is a Pop III star */
      {
        population_type = 3;
        starforminggas_lum_lyman_werner = 0.66 * mass * 1.0e7 * All.HubbleParam;
        find_neighbors = 1;
      }
  else if((metallicity >= 0.001 * GFM_SOLAR_METALLICITY) & (metallicity < 0.1 * GFM_SOLAR_METALLICITY)) /* This is a Pop II star */
      {
        population_type = 2;
        starforminggas_lum_lyman_werner =  0.66 * mass * 1.0e7 * All.HubbleParam;
        find_neighbors = 1;
      }


  if (find_neighbors == 1)
     {
       int nfound = ngb_treefind_variable_threads(pos, hinc, target, mode, threadid, numnodes, firstnode);
       for(n = 0; n < nfound; n++)
         {
           j = Thread[threadid].Ngblist[n];

           dx = NGB_PERIODIC_LONG_X(pos[0] - P[j].Pos[0]);
           dy = NGB_PERIODIC_LONG_Y(pos[1] - P[j].Pos[1]);
           dz = NGB_PERIODIC_LONG_Z(pos[2] - P[j].Pos[2]);

           r2 = dx * dx + dy * dy + dz * dz;

           if(r2 < hinc2 && r2 > 0)
             {
               if (r2 > h2)
                 {
                   local_intensity = starforminggas_lum_lyman_werner / r2;
                 } 
               else
                 {
                   h_inv  = 1.0 / h;
                   u      = sqrt(r2) * h_inv;

                   if(u < 0.5)
                      {
                          wp  = h_inv * (-2.8 + u * u * (5.333333333333 + u * u * (6.4 * u - 9.6)));
                      }
                   else
                     {
                          wp = h_inv * (-3.2 + 0.066666666667 / u + u * u * (10.666666666667 + u * (-16.0 + u * (9.6 - 2.133333333333 * u))));
                     }
                     local_intensity = starforminggas_lum_lyman_werner * wp * wp;
                  }  
               local_intensity /= pow(All.cf_atime , 2);
               if(population_type == 2)
                  SphP[j].StarFormingGasLymanWernerIntensity_type2 += local_intensity;
               if(population_type == 3)
                  SphP[j].StarFormingGasLymanWernerIntensity_type3 += local_intensity;
             }
         }
     }
  return 0;
}

int reset_starforminggas_lyman_werner_intensity_evaluate(int target, int mode, int threadid)
{
  int j, n, numnodes, *firstnode, find_neighbors = 0, population_type;
  double dx, dy, dz, r2, h, h2, hinc, hinc2;
  MyDouble *pos;
  MyFloat mass, metallicity;
#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
#endif

  data_in local, *target_data;

  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      target_data = &local;

      numnodes  = 1;
      firstnode = NULL;
    }
  else
    {
      target_data = &DataGet[target];

      generic_get_numnodes(target, &numnodes, &firstnode);
    }

  pos    = target_data->Pos;
  h      = target_data->hsml;
  hinc   = All.LymanWernerRadiusOfInclusion * h;
  metallicity = target_data->metallicity;
  h2 = h * h;
  hinc2 = hinc * hinc;

  if (metallicity < 0.001 * GFM_SOLAR_METALLICITY) /* This is a Pop III star */
     {
       population_type = 3;
       find_neighbors = 1;
     }
  else if((metallicity >= 0.001 * GFM_SOLAR_METALLICITY) & (metallicity < 0.1 * GFM_SOLAR_METALLICITY)) /* This is a Pop II star */
     {
       population_type = 2;
       find_neighbors = 1;
     }

  if (find_neighbors == 1)
     {
        int nfound = ngb_treefind_variable_threads(pos, hinc, target, mode, threadid, numnodes, firstnode);
        for(n = 0; n < nfound; n++)
          {
            j = Thread[threadid].Ngblist[n];

            dx = NGB_PERIODIC_LONG_X(pos[0] - P[j].Pos[0]);
            dy = NGB_PERIODIC_LONG_Y(pos[1] - P[j].Pos[1]);
            dz = NGB_PERIODIC_LONG_Z(pos[2] - P[j].Pos[2]);

            r2 = dx * dx + dy * dy + dz * dz;

            if(r2 < hinc2 && r2 > 0)
              {
                 if (population_type == 3)   
                    SphP[j].StarFormingGasLymanWernerIntensity_type3 = 0;

                 if (population_type == 2)
                    SphP[j].StarFormingGasLymanWernerIntensity_type2 = 0;
              }
          }
    }
  return 0;
}
#endif

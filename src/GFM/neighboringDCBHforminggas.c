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
#ifdef CHECK_FOR_ENOUGH_GAS_MASS_IN_DCBH_FORMING_POCKETS 

static int neighboringDCBHforminggas_evaluate(int target, int mode, int threadid);

typedef struct
{
  MyDouble Pos[3];
  MyFloat Mass; 
  int Firstnode;
  MyFloat hsml , metallicity , Sfr , lymanwernerintensity;
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
  MyFloat NeighboringDCBHFormingGasMass;
  char dummy;
} data_out;

static data_out *DataResult;

/* routine to store or combine result data */
static void out2particle(data_out *out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES) /* initial store */
    {
      SphP[i].NeighboringDCBHFormingGasMass = out->NeighboringDCBHFormingGasMass;
    }
  else /* merge */
    {
      SphP[i].NeighboringDCBHFormingGasMass += out->NeighboringDCBHFormingGasMass;
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
        MyFloat lymanwernerintensity = SphP[i].StarFormingGasLymanWernerIntensity_type2 + SphP[i].StarFormingGasLymanWernerIntensity_type3;
        if((i < 0) || ((P[i].Type != 0) || (lymanwernerintensity <= All.MinLymanWernerFluxForNewSeed)))
          continue;
        if ((SphP[i].GasIsDense != 1) || (SphP[i].MassMetallicity >= All.MaxMetallicityForAssumingMetalFree * P[i].Mass * GFM_SOLAR_METALLICITY))
          continue;
        neighboringDCBHforminggas_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
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

        neighboringDCBHforminggas_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}


void neighboringDCBHforminggas(void)
{
  int i;
  for(i = 0; i < NumGas; i++)
   {
     SphP[i].NeighboringDCBHFormingGasMass = 0;
   }
  generic_set_MaxNexport();
  generic_comm_pattern(NumPart, kernel_local, kernel_imported);
}




int neighboringDCBHforminggas_evaluate(int target, int mode, int threadid)
{
  int j, n, numnodes, *firstnode, find_neighbors = 0, population_type;
  double dx, dy, dz, r2, h, h2, hinc, hinc2, wp , h_inv, u;
  MyDouble *pos;
  MyFloat mass , neighboringDCBHforminggasmass;
#ifdef PERIODIC
  double xtmp, ytmp, ztmp;
#endif

  data_in local, *target_data;
  data_out out;

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
  h2     = h * h;
  hinc   = All.LymanWernerRadiusOfInclusion * h * 0.1;
  hinc2 = hinc * hinc;
  neighboringDCBHforminggasmass = 0.;
  int nfound = ngb_treefind_variable_threads(pos, hinc, target, mode, threadid, numnodes, firstnode);
  for(n = 0; n < nfound; n++)
         {
           j = Thread[threadid].Ngblist[n];

           dx = NGB_PERIODIC_LONG_X(pos[0] - P[j].Pos[0]);
           dy = NGB_PERIODIC_LONG_Y(pos[1] - P[j].Pos[1]);
           dz = NGB_PERIODIC_LONG_Z(pos[2] - P[j].Pos[2]);

           r2 = dx * dx + dy * dy + dz * dz;

           if(r2 < hinc2)
             {
                  if (SphP[j].StarFormingGasLymanWernerIntensity_type2 + SphP[j].StarFormingGasLymanWernerIntensity_type3 > All.MinLymanWernerFluxForNewSeed)
                    if (SphP[j].MassMetallicity < All.MaxMetallicityForAssumingMetalFree * P[j].Mass * GFM_SOLAR_METALLICITY)
#ifdef SUPPRESS_STARFORMATION_ABOVE_CRITICAL_LYMANWERNERFLUX
                      if(SphP[j].GasIsDense == 1)
#else
                      if(SphP[j].Sfr > 0)
#endif
#ifdef REFINEMENT_HIGH_RES_GAS
                        neighboringDCBHforminggasmass += SphP[j].HighResMass;
#else
                        neighboringDCBHforminggasmass += P[j].Mass; 
#endif
             }
         }
  out.NeighboringDCBHFormingGasMass = neighboringDCBHforminggasmass;
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;
  return 0;
}
#endif
#endif

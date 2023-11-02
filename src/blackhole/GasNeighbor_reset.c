/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/density.c
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
#include "../domain.h"
#include "../proto.h"


#ifdef PREVENT_SPURIOUS_RESEEDING2
/*! Structure for communication during the density computation. Holds data that is sent to other processors.
 */
static int GasNeighbor_evaluate_reset(int target, int mode, int threadid);

static MyFloat *NumNgb, *DhsmlDensityFactor;
#ifdef FIX_SPH_PARTICLES_AT_IDENTICAL_COORDINATES
static MyFloat *MinDist;
#endif

/* local data structure for collecting particle/cell data that is sent to other processors if needed */
typedef struct
{
  MyDouble Pos[3];
  MyFloat Hsml;
#ifdef FIX_SPH_PARTICLES_AT_IDENTICAL_COORDINATES
  MyIDType ID;
#endif

  int Firstnode;
} data_in;

static data_in *DataGet;

/*! \brief Routine that fills the relevant particle/cell data into the input
 *         structure defined above. Needed by generic_comm_helpers2.
 *
 *  \param[out] in Data structure to fill.
 *  \param[in] i Index of particle in P and SphP arrays.
 *  \param[in] firstnode First note of communication.
 *
 *  \return void
 */
static void particle2in(data_in *in, int i, int firstnode)
{
  in->Pos[0] = P[i].Pos[0];
  in->Pos[1] = P[i].Pos[1];
  in->Pos[2] = P[i].Pos[2];
  in->Hsml   = BPP(i).BH_Hsml;
#ifdef FIX_SPH_PARTICLES_AT_IDENTICAL_COORDINATES
  in->ID = P[i].ID;
#endif

  in->Firstnode = firstnode;
}

/* local data structure that holds results acquired on remote processors */
typedef struct
{
 int NeighborOfBlackhole;
} data_out;

static data_out *DataResult;

/*! \brief Routine to store or combine result data. Needed by
 *         generic_comm_helpers2.
 *
 *  \param[in] out Data to be moved to appropriate variables in global
 *  particle and cell data arrays (P, SphP,...)
 *  \param[in] i Index of particle in P and SphP arrays
 *  \param[in] mode Mode of function: local particles or information that was
 *  communicated from other tasks and has to be added locally?
 *
 *  \return void
 */
static void out2particle(data_out *out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES) /* initial store */
    {
      if(P[i].Type == 0)
        {
      	    SphP[i].NeighborOfBlackhole = out->NeighborOfBlackhole;
       }
    }
  else /* combine */
    {
      if(P[i].Type == 0)
        {
           SphP[i].NeighborOfBlackhole = out->NeighborOfBlackhole;
      }
    }
}

#define NO_GENERIC_COMM_PATTERN_FOR_GIVEN_PARTICLES
#include "../generic_comm_helpers2.h"

/*! \brief Routine that defines what to do with local particles.
 *
 *  Calls the *_evaluate function in MODE_LOCAL_PARTICLES.
 *
 *  \return void
 */
static void kernel_local(void)
{
  int idx;
#ifdef GENERIC_ASYNC
  int flag = 0;
#endif

#pragma omp parallel private(idx)
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
              if(generic_polling_primary(count, NumPart))
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
        if(i < 0)
          continue;

        if(GasNeighbor_isactive_reset(i))
          GasNeighbor_evaluate_reset(i, MODE_LOCAL_PARTICLES, threadid);
      }
  }
}

/*! \brief Routine that defines what to do with imported particles.
 *
 *  Calls the *_evaluate function in MODE_IMPORTED_PARTICLES.
 *
 *  \return void
 */
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

        GasNeighbor_evaluate_reset(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

/*! \file density.c
 *  \brief SPH density computation and smoothing length determination
 *
 *  This file contains the "first SPH loop", where the SPH densities and some
 *  auxiliary quantities are computed.  There is also functionality that
 *  corrects the smoothing length if needed.
 */

/*! This function computes the local density for each active SPH particle, the
 * number of neighbours in the current smoothing radius, and the divergence
 * and rotation of the velocity field.  The pressure is updated as well.  If a
 * particle with its smoothing region is fully inside the local domain, it is
 * not exported to the other processors. The function also detects particles
 * that have a number of neighbours outside the allowed tolerance range. For
 * these particles, the smoothing length is adjusted accordingly, and the
 * density() computation is called again.  Note that the smoothing length is
 * not allowed to fall below the lower bound set by MinGasHsml (this may mean
 * that one has to deal with substantially more than normal number of
 * neighbours.)
 */

/*! \brief Main function of SPH density calculation.
 *
 *  This function computes the local density for each active SPH particle and
 *  the number of weighted neighbors in the current smoothing radius. If a
 *  particle with its smoothing region is fully inside the local domain, it is
 *  not exported to the other processors. The function also detects particles
 *  that have a number of neighbors outside the allowed tolerance range. For
 *  these particles, the smoothing length is adjusted accordingly, and the
 *  computation is called again.
 *
 *  \return void
 */

void GasNeighbor_reset(void)
{
  CPU_Step[CPU_MISC] += measure_time();

  generic_set_MaxNexport();

  generic_comm_pattern(NumPart, kernel_local, kernel_imported);
  /* collect some timing information */
  CPU_Step[CPU_INIT] += measure_time();
}


/* This function represents the core of the SPH density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */
static int GasNeighbor_evaluate_reset(int target, int mode, int threadid)
{
  int j, n;
  int numngb, numnodes, *firstnode;
  double h, h2, hinv, hinv3, hinv4;
  MyFloat rho;
  double wk, dwk;
  double dx, dy, dz, r, r2, u, mass_j;
  MyFloat weighted_numngb;
  MyFloat dhsmlrho;
  MyDouble *pos;
#ifdef OTVET_FLUXLIMITER
  int k;
  double mj_dwk_r;
  MyFloat grad_ngamma[3][OT_N_BINS];
#endif

#ifdef FIX_SPH_PARTICLES_AT_IDENTICAL_COORDINATES
  MyFloat mindist = MAX_REAL_NUMBER;
  MyIDType ID;
#endif
#ifdef OTVET_FLUXLIMITER
  for(k = 0; k < OT_N_BINS; k++)
    for(j = 0; j < 3; j++)
      grad_ngamma[j][k] = 0;
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

  pos = target_data->Pos;
  h   = target_data->Hsml;
#ifdef FIX_SPH_PARTICLES_AT_IDENTICAL_COORDINATES
  ID = target_data->ID;
#endif

  h2   = h * h;
  hinv = 1.0 / h;
#ifndef TWODIMS
  hinv3 = hinv * hinv * hinv;
#else
  hinv3 = hinv * hinv / boxSize_Z;
#endif
  hinv4 = hinv3 * hinv;

  numngb = 0;
  rho = weighted_numngb = dhsmlrho = 0;

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, threadid, numnodes, firstnode);
  
  for(int n = 0; n < nfound; n++)
    {
      int j = Thread[threadid].Ngblist[n];
      if(P[j].Mass > 0 && P[j].ID != 0)
        {
          double r2 = Thread[threadid].R2list[n];
          if (r2 < h2)
            SphP[j].NeighborOfBlackhole = 1;
        }
    }
 
// /* Now collect the result at the right place */
//  if(mode == MODE_LOCAL_PARTICLES)
//    out2particle(&out, target, MODE_LOCAL_PARTICLES);
//  else
//    DataResult[target] = out;

  return 0;
}

int GasNeighbor_isactive_reset(int n)
{

  if(P[n].Type == 5)
    return 1;

  return 0;	

}

#endif

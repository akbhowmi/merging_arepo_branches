/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/fof/fof_spinmeasurement.c
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
#include <sys/stat.h>
#include <sys/types.h>

#include "../allvars.h"
#include "../proto.h"

#if (defined(CALCULATE_SPIN_STARFORMINGGAS)) || (defined(FOF) && defined(GFM_BIPOLAR_WINDS) && (GFM_BIPOLAR_WINDS == 3))

#include "fof.h"

#ifdef GFM_BIPOLAR_WINDS
#define SPIN_MEASURE_RADIUS_FRACTION 0.1
#endif

#ifdef CALCULATE_SPIN_STARFORMINGGAS
#define SPIN_MEASURE_RADIUS_FRACTION 1.0
#endif

/* local data structure for collecting particle/cell data that is sent to other processors if needed */
typedef struct
{
  MyDouble Pos[3];
  MyFloat Mass;

  int Firstnode;
} data_in;

static data_in *DataGet;

/* routine that fills the relevant particle/cell data into the input structure defined above */
static void particle2in(data_in *in, int i, int firstnode)
{
  in->Pos[0] = Group[i].Pos_MinPotential[0];
  in->Pos[1] = Group[i].Pos_MinPotential[1];
  in->Pos[2] = Group[i].Pos_MinPotential[2];
  in->Mass   = Group[i].Mass;
  in->Firstnode = firstnode;
}

/* local data structure that holds results acquired on remote processors */
typedef struct
{
  MyFloat DensGasMass;
  MyFloat DensGasCenter[3];
  MyFloat DensGasMomentum[3];
  MyFloat DensGasAngMomentum[3];
#ifdef CALCULATE_SPIN_STARFORMINGGAS
  MyFloat Mean_Temperature;
#endif
} data_out;

static data_out *DataResult;

/* routine to store or combine result data */
static void out2particle(data_out *out, int i, int mode)
{
  int k;

  if(mode == MODE_LOCAL_PARTICLES) /* initial store */
    {
      Group[i].DensGasMass = out->DensGasMass;
      for(k = 0; k < 3; k++)
        {
          Group[i].DensGasCenter[k]      = out->DensGasCenter[k];
          Group[i].DensGasMomentum[k]    = out->DensGasMomentum[k];
          Group[i].DensGasAngMomentum[k] = out->DensGasAngMomentum[k];
        }
#ifdef CALCULATE_SPIN_STARFORMINGGAS
      Group[i].MeanTemperature = out->Mean_Temperature;
#endif    
     }
  else /* combine */
    {
      Group[i].DensGasMass += out->DensGasMass;
      for(k = 0; k < 3; k++)
        {
          Group[i].DensGasCenter[k] += out->DensGasCenter[k];
          Group[i].DensGasMomentum[k] += out->DensGasMomentum[k];
          Group[i].DensGasAngMomentum[k] += out->DensGasAngMomentum[k];
        }
#ifdef CALCULATE_SPIN_STARFORMINGGAS
      Group[i].MeanTemperature += out->Mean_Temperature;
#endif
    }
}

#include "../generic_comm_helpers2.h"

static int fof_spin_evaluate(int target, int mode, int threadid);

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
              if(generic_polling_primary(count, Ngroups))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        i = NextParticle++;

        if(i >= Ngroups)
          break;

        fof_spin_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
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

        fof_spin_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

void fof_spin_measurement(void)
{
  generic_set_MaxNexport();

  generic_comm_pattern(Ngroups, kernel_local, kernel_imported);

  /* do final operations on results */
  for(int i = 0; i < Ngroups; i++)
    {
#ifdef CALCULATE_SPIN_STARFORMINGGAS
      Group[i].DensGasDimensionlessSpin = 0;
      Group[i].DensGasDimensionlessSpin_Max = 0;
      
      double mass	   = Group[i].Mass;
      double vvir	   = pow(10 * All.G * All.cf_H * mass, 1.0 / 3);
      double rvir_phys = vvir / (All.cf_H * 10.);     /* physical */
      double rvir	   =   (rvir_phys / All.cf_atime); /* comoving */
      double mass_unit_conversion = 1e10;
      double sectoyrs = 3.171e-08;
      double Mpctokm = 3.0857e19;
      double Grav_constant = 6.67e-11;
      double kpctoMpc = 0.001;
      double kmtom = 1000;
      double Mass_of_sun = 2e30;
      double rvir_in_m = rvir * kpctoMpc * Mpctokm * kmtom;
      double HaloCircVel = sqrt(Grav_constant * mass * mass_unit_conversion * Mass_of_sun / rvir_in_m) / kmtom / sqrt(All.Time);
      double mu = 0.59;     //Mean molcular weight of primordial gas
      double J_circ = (HaloCircVel * rvir) * sqrt(2); //Specific
      Group[i].VirialTemperature=4e4 * (mu / 1.2) * pow(mass * 1e2, 2.0/3) / All.Time / 10;
      Group[i].RvirEstimate = rvir;
      Group[i].MeanTemperature /= Group[i].DensGasMass;
#endif
      if(Group[i].DensGasMass > 0)
        {
          Group[i].DensGasCenter[0] /= Group[i].DensGasMass;
          Group[i].DensGasCenter[1] /= Group[i].DensGasMass;
          Group[i].DensGasCenter[2] /= Group[i].DensGasMass;
          Group[i].DensGasAngMomentum[0] -=
              Group[i].DensGasCenter[1] * Group[i].DensGasMomentum[2] - Group[i].DensGasCenter[2] * Group[i].DensGasMomentum[1];
          Group[i].DensGasAngMomentum[1] -=
              Group[i].DensGasCenter[2] * Group[i].DensGasMomentum[0] - Group[i].DensGasCenter[0] * Group[i].DensGasMomentum[2];
          Group[i].DensGasAngMomentum[2] -=
              Group[i].DensGasCenter[0] * Group[i].DensGasMomentum[1] - Group[i].DensGasCenter[1] * Group[i].DensGasMomentum[0];
#ifdef CALCULATE_SPIN_STARFORMINGGAS
          Group[i].DensGasAngMomentum[0] /= All.Time;     // Converting to peculiar velocity units
          Group[i].DensGasAngMomentum[1] /= All.Time;
          Group[i].DensGasAngMomentum[2] /= All.Time;
          for(int k = 0; k < 3; k++)
            {
               Group[i].DensGasDimensionlessSpin += Group[i].DensGasAngMomentum[k] * Group[i].DensGasAngMomentum[k];
            }
          Group[i].DensGasDimensionlessSpin = sqrt(Group[i].DensGasDimensionlessSpin);
          Group[i].DensGasDimensionlessSpin /= Group[i].DensGasMass; //Specific 
          Group[i].DensGasDimensionlessSpin /= J_circ;
          double md = 0.05;
          double jd = 0.05;
          double Q_c = 2;
          Group[i].DensGasDimensionlessSpin_Max = md * md * Q_c / 8 / jd * sqrt(Group[i].VirialTemperature / Group[i].MeanTemperature);
#endif
        }
    }
}

/*! This function represents the core of the blackhole density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */
static int fof_spin_evaluate(int target, int mode, int threadid)
{
  int j, k, n;
  int numnodes, *firstnode;
  double h, h2, mass;
  double dx, dy, dz, r2;
  MyDouble *pos;
  MyFloat xyz[3];
  MyFloat DensGasMass;
  MyFloat DensGasCenter[3];
  MyFloat DensGasMomentum[3];
  MyFloat DensGasAngMomentum[3];
#ifdef CALCULATE_SPIN_STARFORMINGGAS
  MyFloat mu, Temperature, Mean_Temperature;
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

  pos  = target_data->Pos;
  mass = target_data->Mass;

  DensGasMass = 0;
  for(j = 0; j < 3; j++)
    {
      DensGasCenter[j]      = 0;
      DensGasMomentum[j]    = 0;
      DensGasAngMomentum[j] = 0;
    }
#ifdef CALCULATE_SPIN_STARFORMINGGAS
  Mean_Temperature = 0; 
#endif


  /* estimate the virial radius based on the mass */
  double vvir      = pow(10 * All.G * All.cf_H * mass, 1.0 / 3);
  double rvir_phys = vvir / (10 * All.cf_H);     /* physical */
  double rvir      = (rvir_phys / All.cf_atime); /* comoving */

  h  = SPIN_MEASURE_RADIUS_FRACTION * rvir;
  h2 = h * h;

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, threadid, numnodes, firstnode);

  for(n = 0; n < nfound; n++)
    {
      j = Thread[threadid].Ngblist[n];

      dx = FOF_NEAREST_LONG_X(pos[0] - P[j].Pos[0]);
      dy = FOF_NEAREST_LONG_Y(pos[1] - P[j].Pos[1]);
      dz = FOF_NEAREST_LONG_Z(pos[2] - P[j].Pos[2]);

      r2 = dx * dx + dy * dy + dz * dz;

      if(r2 < h2 && P[j].Type == 0)
#ifndef CALCULATE_SPIN_STARFORMINGGAS
        if(SphP[j].Sfr > 0)
#endif
          {
            DensGasMass += P[j].Mass;
            for(k = 0; k < 3; k++)
              {
                xyz[k] = fof_periodic_nearest(P[j].Pos[k] - pos[k]);
                DensGasCenter[k] += P[j].Mass * xyz[k];
                DensGasMomentum[k] += P[j].Mass * P[j].Vel[k];
              }

            DensGasAngMomentum[0] += P[j].Mass * (xyz[1] * P[j].Vel[2] - xyz[2] * P[j].Vel[1]);
            DensGasAngMomentum[1] += P[j].Mass * (xyz[2] * P[j].Vel[0] - xyz[0] * P[j].Vel[2]);
            DensGasAngMomentum[2] += P[j].Mass * (xyz[0] * P[j].Vel[1] - xyz[1] * P[j].Vel[0]);
#ifdef CALCULATE_SPIN_STARFORMINGGAS
            mu = 4./(1 + 3 * H_FRACTION + 4 * H_FRACTION * SphP[j].Ne) * MASS_OF_PROTON_CGS;
            Temperature = (5./3 - 1) * SphP[j].Utherm * 1e10  * mu / BOLTZMANN_CONSTANT_CGS;
            Mean_Temperature += Temperature * P[j].Mass;
#endif
          }
    }
#ifdef CALCULATE_SPIN_STARFORMINGGAS
  out.Mean_Temperature = Mean_Temperature; 
#endif

  out.DensGasMass = DensGasMass;
  for(k = 0; k < 3; k++)
    {
      out.DensGasCenter[k]      = DensGasCenter[k];
      out.DensGasMomentum[k]    = DensGasMomentum[k];
      out.DensGasAngMomentum[k] = DensGasAngMomentum[k];
    }

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

#endif

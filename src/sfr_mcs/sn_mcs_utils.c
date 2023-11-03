/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/sfr_mcs/sn_mcs.c
 * \date        04/2018
 * \author      Matthew C Smith
 * \brief
 * \details     Originally developed in 2015, ported into main repo 2018.
                Please contact the author before use.
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

#ifdef SN_NO_ENERGY
#define SN_ENERGY 0.0
#else
#define SN_ENERGY 1e51 /* ergs */
#endif

#ifdef MECHANICAL_FEEDBACK
#define P_TERMINAL 3e5 /* Msun km/s */
#endif

void init_sne(void)
{
  if(All.ComovingIntegrationOn)
    All.MaxSNStarMass = All.TargetGasMass * All.MaxSNStarMassFac;

  All.MaxSNEvalTimestep *= (SEC_PER_MEGAYEAR / All.UnitTime_in_s * All.HubbleParam);

  All.SupernovaEnergy = SN_ENERGY * All.HubbleParam / All.UnitEnergy_in_cgs;
#ifdef MECHANICAL_FEEDBACK
  All.SupernovaTerminalMomentum = P_TERMINAL * SOLAR_MASS * 1e5 * All.HubbleParam / (All.UnitMass_in_g * All.UnitVelocity_in_cm_per_s);
#endif
}

void read_sb99_tables(void)
{
  int i;
  double *tmpdbl;
  hid_t file_id, dataset, datatype;
  char fname[MAXLEN_PATH], setname[MAXLEN_PATH];

#ifndef SB99_FIXED_Z
  int j;
#endif
  file_id = my_H5Fopen(All.SB99TablesPath, H5F_ACC_RDONLY, H5P_DEFAULT);

  /* read number of elements */

  dataset = my_H5Dopen(file_id, "Number_of_timesteps");
  my_H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &sb99.N_t, "Number_of_timesteps");
  my_H5Dclose(dataset, "Number_of_timesteps");

#ifndef SB99_FIXED_Z
  dataset = my_H5Dopen(file_id, "Number_of_metallicities");
  my_H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &sb99.N_Z, "Number_of_metallicities");
  my_H5Dclose(dataset, "Number_of_metallicities");
#endif
  mpi_printf("SN_MCS: Starburst99 header complete...\n");

  /* Allocate arrays */
  sb99.Timesteps = (MyFloat *)mymalloc("sb99.Timesteps", sb99.N_t * sizeof(MyFloat));
#ifndef SB99_FIXED_Z
  sb99.Metallicities = (MyFloat *)mymalloc("sb99.Metallicities", sb99.N_Z * sizeof(MyFloat));
  sb99.Rates         = (MyFloat **)mymalloc("sb99.Rates", sb99.N_Z * sizeof(MyFloat *));
  for(i = 0; i < sb99.N_Z; i++)
    sb99.Rates[i] = (MyFloat *)mymalloc("sb99.Rates", sb99.N_t * sizeof(MyFloat));
#else
  sb99.Rates = (MyFloat *)mymalloc("sb99.Rates", sb99.N_t * sizeof(MyFloat));
#endif

  mpi_printf("SN_MCS: Starburst99 arrays allocated...\n");

  tmpdbl  = (double *)mymalloc("Timesteps", sb99.N_t * sizeof(double));
  dataset = my_H5Dopen(file_id, "Timesteps");
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, "Timesteps");
  my_H5Dclose(dataset, "Timesteps");
  for(i = 0; i < sb99.N_t; i++)
    sb99.Timesteps[i] = tmpdbl[i];
  myfree(tmpdbl);

  sb99.t_min = sb99.Timesteps[0];
  sb99.t_max = sb99.Timesteps[sb99.N_t - 1];

#ifndef SB99_FIXED_Z
  tmpdbl  = (double *)mymalloc("Metallicities", sb99.N_Z * sizeof(double));
  dataset = my_H5Dopen(file_id, "Metallicities");
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, "Metallicities");
  my_H5Dclose(dataset, "Metallicities");
  for(i = 0; i < sb99.N_Z; i++)
    sb99.Metallicities[i] = tmpdbl[i];
  myfree(tmpdbl);

  char *tempname[sb99.N_Z];
  datatype = my_H5Tcopy(H5T_C_S1);
  my_H5Tset_size(datatype, H5T_VARIABLE);
  dataset = my_H5Dopen(file_id, "Metallicity_Names");
  my_H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, tempname, "Metallicity_Names");
  my_H5Dclose(dataset, "Metallicity_Names");
  my_H5Tclose(datatype);

  tmpdbl = (double *)mymalloc("Rates", sb99.N_t * sizeof(double));
  for(i = 0; i < sb99.N_Z; i++)
    {
      sprintf(setname, "/Rates/%s", tempname[i]);
      dataset = my_H5Dopen(file_id, setname);
      my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, setname);
      my_H5Dclose(dataset, setname);

      for(j = 0; j < sb99.N_t; j++)
        sb99.Rates[i][j] = tmpdbl[j];
    }
  myfree(tmpdbl);
#else
  tmpdbl     = (double *)mymalloc("Rates", sb99.N_t * sizeof(double));
  sprintf(setname, "Rates/%s", All.SB99_metallicity_name);
  dataset = my_H5Dopen(file_id, setname);
  my_H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpdbl, setname);
  my_H5Dclose(dataset, setname);
  for(i = 0; i < sb99.N_t; i++)
    sb99.Rates[i] = tmpdbl[i];
  myfree(tmpdbl);
#endif
  my_H5Fclose(file_id, fname);

  mpi_printf("SN_MCS: Starburst99 setup complete...\n");
}

#ifndef SB99_FIXED_Z
int get_sb99_z_index(MyFloat metallicity)
{
  int left;
  MyFloat middle_val;

  left = 0;

  while(left < (sb99.N_Z - 1))
    {
      middle_val = (sb99.Metallicities[left] + sb99.Metallicities[left + 1]) / 2.0;
      if(metallicity < middle_val)
        return left;

      left++;
    }
  return left;
}
#endif

void get_sb99_t_indicies(MyFloat t, MyFloat *timesteps, int N_t, int *it_low, int *it_high, MyFloat *delta)
{
  MyFloat deltat, dt_local;
  int i1, i2;

  if(t > timesteps[0])
    {
      for(i1 = 0; i1 < N_t - 1 && t > timesteps[i1 + 1]; i1++)
        ;

      i2 = i1 + 1;

      if(i2 >= N_t)
        i2 = N_t - 1;

      if(t >= timesteps[0] && t <= timesteps[N_t - 1])
        dt_local = t - timesteps[i1];
      else
        dt_local = 0;

      deltat = timesteps[i2] - timesteps[i1];

      if(delta > 0)
        dt_local = dt_local / deltat;
      else
        dt_local = 0;
    }
  else
    {
      i1       = 0;
      i2       = 0;
      dt_local = 0.0;
    }

  *it_low  = i1;
  *it_high = i2;
  *delta   = dt_local;
}

#ifndef SB99_FIXED_Z
double get_sn_rate(MyFloat t, MyFloat metallicity)
{
  int it_low, it_high, iz;
  double delta, rate;

  get_sb99_t_indicies(t, sb99.Timesteps, sb99.N_t, &it_low, &it_high, &delta);
  iz   = get_sb99_z_index(metallicity);
  rate = sb99.Rates[iz][it_low] + (sb99.Rates[iz][it_high] - sb99.Rates[iz][it_low]) * delta;
  return rate;
}
#else
double get_sn_rate(MyFloat t)
{
  int it_low, it_high;
  double delta, rate;

  get_sb99_t_indicies(t, sb99.Timesteps, sb99.N_t, &it_low, &it_high, &delta);
  rate = sb99.Rates[it_low] + (sb99.Rates[it_high] - sb99.Rates[it_low]) * delta;
  return rate;
}
#endif

#ifdef SN_MCS_LOG
#if !defined(SN_MCS_LOG_MAX) || !defined(SN_MCS_LOG_MIN) || !defined(SN_MCS_LOG_N)
#error "SN_MCS_LOG requires SN_MCS_LOG_MAX, SN_MCS_LOG_MIN and SN_MCS_LOG_N to be defined"
#endif
void setup_sn_log(void)
{
  int i;

  sn_dens_hist = gsl_histogram_alloc(SN_MCS_LOG_N);
  gsl_histogram_set_ranges_uniform(sn_dens_hist, SN_MCS_LOG_MIN, SN_MCS_LOG_MAX);
  /* Write bin edges to file */
  if(ThisTask == 0)
    {
      for(i = 0; i <= sn_dens_hist->n; i++)
        fprintf(FdSNdens, "%g ", sn_dens_hist->range[i]);

      fprintf(FdSNdens, "\n");
      fflush(FdSNdens);
    }
}

void sn_add_to_log(int i)
{
  double n;

  n = SphP[i].Density * All.cf_a3inv * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;
  n *= HYDROGEN_MASSFRAC / PROTONMASS;

  if(gsl_histogram_increment(sn_dens_hist, log10(n)))
    {
      printf("SN_MCS_LOG: OUT OF BOUNDS on Task: %d Time: %g log10(n): %g\n", ThisTask, All.Time, log10(n));
      /* out of bounds, add to top or bottom bin */
      if(log10(n) >= gsl_histogram_max(sn_dens_hist))
        sn_dens_hist->bin[sn_dens_hist->n - 1] += 1;
      else if(log10(n) >= gsl_histogram_min(sn_dens_hist))
        sn_dens_hist->bin[0] += 0;
    }
}

void write_sn_dens_log(void)
{
  int i;
  double *denshistogram;

  mpi_printf("SN_MCS_LOG: Writing density histogram...\n");

  if(ThisTask == 0)
    denshistogram = (double *)mymalloc("denshistogram", sn_dens_hist->n * sizeof(double));

  /* Sum histograms from all tasks and reset */
  MPI_Reduce(sn_dens_hist->bin, denshistogram, sn_dens_hist->n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  gsl_histogram_reset(sn_dens_hist);

  if(ThisTask == 0)
    {
      fprintf(FdSNdens, "%e", All.Time);
      for(i = 0; i < sn_dens_hist->n; i++)
        fprintf(FdSNdens, " %g", denshistogram[i]);

      fprintf(FdSNdens, "\n");
      fflush(FdSNdens);

      myfree(denshistogram);
    }
}
#endif  // SN_MCS_LOG

#endif

/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/sfr_mcs/sfr_mcs.c
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

#ifdef SFR_MCS

#ifndef USE_SFR
#error "SFR_MCS requires USE_SFR"
#endif

#ifndef STELLARAGE
#error "SFR_MCS requires STELLARAGE"
#endif

#if defined(JEANS_PRESSURE_LIMIT_MCS) && (defined(JEANS_PRESSURE_LIMIT) || defined(ENFORCE_JEANS_STABILITY_OF_CELLS))
#error \
    "You have defined JEANS_PRESSURE_LIMIT_MCS along with some other Jeans criteria based EOS,\nwhich you probably don't want to do."
#endif

void init_star_formation(void)
{
  if(RestartFlag != RESTART_RESTART)
    {
      All.DensThreshold *= ((PROTONMASS / HYDROGEN_MASSFRAC) / (All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam));
      if(All.ComovingIntegrationOn)
        All.OverDensThresh = All.CritOverDensity * All.OmegaBaryon * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);
#ifdef SN_MCS
      init_sne();
#endif
    }
}

#ifdef SFR_MCS_LOG
#if !defined(SFR_MCS_LOG_MAX) || !defined(SFR_MCS_LOG_MIN) || !defined(SFR_MCS_LOG_N)
#error "SFR_MCS_LOG requires SFR_MCS_LOG_MAX, SFR_MCS_LOG_MIN and SFR_MCS_LOG_N to be defined"
#endif
void setup_sf_log(void)
{
  int i;

  sf_dens_hist = gsl_histogram_alloc(SFR_MCS_LOG_N);
  if(RestartFlag != RESTART_RESTART)
    {
      gsl_histogram_set_ranges_uniform(sf_dens_hist, SFR_MCS_LOG_MIN, SFR_MCS_LOG_MAX);
      /* Write bin edges to file */
      if(ThisTask == 0)
        {
          for(i = 0; i <= sf_dens_hist->n; i++)
            fprintf(FdSFdens, "%g ", sf_dens_hist->range[i]);

          fprintf(FdSFdens, "\n");
          fflush(FdSFdens);
        }
    }
}

void sf_add_to_log(double rho)
{
  double n;
  n = rho * All.cf_a3inv * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam;
  n *= HYDROGEN_MASSFRAC / PROTONMASS;
  if(gsl_histogram_increment(sf_dens_hist, log10(n)))
    {
      printf("SFR_MCS_LOG: OUT OF BOUNDS on Task: %d Time: %g log10(n): %g\n", ThisTask, All.Time, log10(n));
      /* out of bounds, add to top or bottom bin */
      if(log10(n) >= gsl_histogram_max(sf_dens_hist))
        sf_dens_hist->bin[sf_dens_hist->n - 1] += 1;
      else if(log10(n) >= gsl_histogram_min(sf_dens_hist))
        sf_dens_hist->bin[0] += 0;
    }
}

void write_sf_dens_log(void)
{
  int i;
  double *denshistogram;

  mpi_printf("SFR_MCS_LOG: Writing density histogram...\n");

  if(ThisTask == 0)
    denshistogram = (double *)mymalloc("denshistogram", sf_dens_hist->n * sizeof(double));

  /* Sum histograms from all tasks and reset */
  MPI_Reduce(sf_dens_hist->bin, denshistogram, sf_dens_hist->n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  gsl_histogram_reset(sf_dens_hist);

  if(ThisTask == 0)
    {
      fprintf(FdSFdens, "%e", All.Time);
      for(i = 0; i < sf_dens_hist->n; i++)
        fprintf(FdSFdens, " %g", denshistogram[i]);

      fprintf(FdSFdens, "\n");
      fflush(FdSFdens);

      myfree(denshistogram);
    }
}
#endif

void cooling_and_starformation(void)
{
  CPU_Step[CPU_MISC] += measure_time();

  int idx, i, flag_sf;
  double dens, t_ff;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Mass == 0 && P[i].ID == 0)
        continue; /* skip cells that have been swallowed or eliminated */

      flag_sf = 0;

      dens = SphP[i].Density;

      SphP[i].Sfr = 0.0;

      if(dens * All.cf_a3inv >= All.DensThreshold)
        flag_sf = 1;

      if(All.ComovingIntegrationOn) /* to protect against SF at too high redshift */
        if(dens < All.OverDensThresh)
          flag_sf = 0;

      cool_cell(i);

      if(flag_sf == 1)
        {
          t_ff        = sqrt(3.0 * M_PI / (32.0 * All.G * dens * All.cf_a3inv));
          SphP[i].Sfr = All.SfrEfficiency * P[i].Mass / t_ff;
          SphP[i].Sfr *= (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);
          TimeBinSfr[P[i].TimeBinHydro] += SphP[i].Sfr;
        }
    } /* End of loop over active particles */

  CPU_Step[CPU_COOLINGSFR] += measure_time();
}

double get_starformation_rate(int i)
{
  double t_ff, dens;
  double rateOfSF = 0.0;

  dens = SphP[i].Density;

  if(dens * All.cf_a3inv >= All.DensThreshold)
    {
      t_ff     = sqrt(3.0 * M_PI / (32.0 * All.G * dens * All.cf_a3inv));
      rateOfSF = All.SfrEfficiency * P[i].Mass / t_ff;
      rateOfSF *= (All.UnitMass_in_g / SOLAR_MASS) / (All.UnitTime_in_s / SEC_PER_YEAR);
    }

  return rateOfSF;
}
#endif /*SFR_MCS*/

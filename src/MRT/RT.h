/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/anisotropic_RT/RT.h
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

#include "../allvars.h"

#ifdef MRT

#define fTOLERENCE 1.000000001

#if defined(MRT_RIEMANN_HLLE) || defined(MRT_RIEMANN_HLLE_NEW)
double lambda1[101][101];
double lambda2[101][101];
double lambda3[101][101];
double lambda4[101][101];

void readin_hlle_eingenvalues(void);

#endif

#ifdef MRT_SINGLE_STAR

void read_stellar_table(void);

#endif

void mpi_printf_rt(const int flag, const char *fmt, ...);

void init_RT(void);
void RT_initialize_cell(int i);
void add_source_fluxes(void);
void set_VET_single(int, struct sph_particle_data *);
void mrt_update_chemistry(void);

#ifdef MRT_COMOVING
void cell_do_lorentz_boost(int, struct sph_particle_data *, double, double, double);
void do_comoving_frame_source_terms(void);
#endif

#ifdef MRT_IR
#ifdef MRT_IR_GRAIN_KAPPA
void read_grain_kappa_data(void);
extern int IR_N_pts;
extern double *IR_logT, *IR_logkappaP, *IR_logkappaR;
extern gsl_interp_accel *accIR_kappaP;
extern gsl_interp_accel *accIR_kappaR;
extern gsl_spline *splineIR_kappaP;
extern gsl_spline *splineIR_kappaR;
#endif
double mrt_update_IR_cooling(int, double);
#ifdef MRT_IR_LTE_GSL
int mrt_IR_rate_ODEs(double, const double *, double *, void *);
int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
#endif
#endif

#ifdef MRT_LOCAL_FEEDBACK
double inject_photons_from_star(int, double, double);
#endif

#if defined(MRT_IR) || defined(MRT_UV_ONLY_DUST)
void set_kappa_times_rho_IR(int, struct sph_particle_data *);
#endif

#ifdef MRT_UV_ONLY_DUST
void mrt_update_chemistry_only_dust(void);
#endif

#ifndef MRT_NO_UV
void initialize_ionic_species(void);
#endif

void set_full_ionization_mrt(int);

#ifdef MRT_RADIATION_PRESSURE
void do_radiation_pressure_source_terms(void);
#endif

#ifdef COOLING
int do_cooling_mrt(int i);
#ifdef MRT_EQUIL_CHEM_COOL
void update_radiation_state_mrt(int i, PhotoCurrent *ppc);
void update_chem_mrt(int i, double dtcool, GasState *pgs);
#endif
#endif

#if defined(MRT_COOLING_HEATING)
double mrt_DoHeating(int, double);
double mrt_DoCooling(int, double);
double mrt_get_cooling_rate(int, double);
double mrt_get_heating_rate(int);
double mrt_GetCoolingTime(int, double, double, double *);
#endif

#ifdef MRT_COUPLED_THERMOCHEMISTRY
void mrt_thermochemistry(void);
#endif

#ifdef MRT_CHEMISTRY_PS2009
void mrt_update_chemistry_ps2009(void);
void mrt_write_stats(void);
void mrt_IR_chemistry(void);
#endif

#ifdef MRT_IR_ONLY_CHEMISTRY
void mrt_IR_chemistry(void);
#endif

#ifdef MRT_SETUP_SPECIAL_BOUNDARIES
void mrt_setup(int);
void exchange_vector(void);
#endif

#ifdef MRT_CHEMISTRY_PS2011
int mrt_rate_ODEs(double, const double *, double *, void *);
void mrt_update_chemistry_ps2011(void);
#endif

#ifdef MRT_CHEM_SG
void set_rates_sgchem(void);
#endif

#ifndef MRT_NO_UV
void mrt_get_sigma(void);
void mrt_get_luminosities(void);
#endif

#ifdef MRT_TIME_EXTRAPOLATION
void calculate_div_F(struct state *dl, struct state *st);
#endif

#if defined(MRT_SOURCES) || defined(MRT_LOCAL_FEEDBACK)
extern int Nsource;
extern double PhotRate;
extern int Nfreq, N_age, N_metallicity;
extern double *PhotEnergy, *LumPhotons;
extern double *LogAge, *LogMetallicity;
extern double *lum_tab;

void mrt_load_spectrum_table(void);
void mrt_free_spectrum_table(void);

#ifdef MRT_STARS
void start_stellar_sources(void);
void end_stellar_sources(void);
void do_ionizing_stellar_sources(void);
void add_ionizing_stellar_radiation(void);
#endif

#ifdef MRT_BH
void start_blackhole_sources(void);
void end_blackhole_sources(void);
void do_ionizing_blackhole_sources(void);
void add_ionizing_blackhole_radiation(void);
#ifdef MRT_BH_BIPOLAR
int is_cell_within_cone(MyIDType, MyDouble *);
#endif
#endif

#endif /* MRT_SOURCES */

#endif /* MRT */

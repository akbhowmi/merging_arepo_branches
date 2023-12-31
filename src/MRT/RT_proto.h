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
//#include "../mesh.h"
#include "../proto.h"

#ifdef MRT

void mrt_run(void);
void mrt_run_sources(void);
void exchange_primitive_variables_RT(void);
void exchange_primitive_variables_and_gradients_RT(void);
void update_primitive_variables_RT(void);
void calculate_gradients_RT(void);
void compute_interface_fluxes_RT(tessellation *T);
int face_get_state_RT(tessellation *T, int p, int i, struct state *st);
void state_convert_to_local_frame_RT(struct state *st, double *vel_face, double hubble_a, double atime);
void face_do_time_extrapolation_RT(struct state *delta, struct state *st, double atime);
void face_do_spatial_extrapolation_RT(struct state *delta, struct state *st, struct state *st_other);
void face_add_extrapolations_RT(struct state *st_face, struct state *delta_time, struct state *delta_space);
void face_add_extrapolation_RT(struct state *st_face, struct state *delta);
void face_turn_velocities_RT(struct state *st, struct geometry *geom);
void face_turnback_velocities_RT(struct state_face *st_face, struct geometry *geom);
void face_limit_fluxes_RT(struct state *st_L, struct state *st_R, struct state *st_center_L, struct state *st_center_R,
                          struct fluxes *flux, double dt, double *count, double *count_reduced);
void apply_flux_list_RT(void);
double face_timestep_RT(struct state *state_L, struct state *state_R, double *hubble_a, double *atime);

#ifdef MRT_SINGLE_STAR
extern struct stellarparameters
{
  double *Age;
  double *Nphot;
  double **frac;

  double **sigH2;
  double **EH2;
  double **PH2;

  double **sigH;
  double **EH;
  double **PH;

  double **sigHe;
  double **EHe;
  double **PHe;

} StellarParameters;
#endif

static inline double interpolate_2d(double (*table)[101], int i, int j, double dx, double dy)
{
  return (1 - dx) * (1 - dy) * table[i][j] + (1 - dx) * dy * table[i][j + 1] + dx * (1 - dy) * table[i + 1][j] +
         dx * dy * table[i + 1][j + 1];
}

static inline MyFloat interpolate_4d(MyFloat ****table, int i, int j, int k, int l, double dx, double dy, double dz, double dw)
{
  int il = i, jl = j, kl = k, ll = l;
  int ir = i + 1, jr = j + 1, kr = k + 1, lr = l + 1;

  double dxl = 1 - dx, dyl = 1 - dy, dzl = 1 - dz, dwl = 1 - dw;
  double dxr = dx, dyr = dy, dzr = dz, dwr = dw;
  if(dxr == 0)
    ir = i;

  if(dyr == 0)
    jr = j;

  if(dzr == 0)
    kr = k;

  if(dwr == 0)
    lr = l;

  return dxl * dyl * dzl * dwl * table[il][jl][kl][ll] + dxl * dyl * dzl * dwr * table[il][jl][kl][lr] +
         dxl * dyl * dzr * dwl * table[il][jl][kr][ll] + dxl * dyl * dzr * dwr * table[il][jl][kr][lr] +
         dxl * dyr * dzl * dwl * table[il][jr][kl][ll] + dxl * dyr * dzl * dwr * table[il][jr][kl][lr] +
         dxl * dyr * dzr * dwl * table[il][jr][kr][ll] + dxl * dyr * dzr * dwr * table[il][jr][kr][lr] +
         dxr * dyl * dzl * dwl * table[ir][jl][kl][ll] + dxr * dyl * dzl * dwr * table[ir][jl][kl][lr] +
         dxr * dyl * dzr * dwl * table[ir][jl][kr][ll] + dxr * dyl * dzr * dwr * table[ir][jl][kr][lr] +
         dxr * dyr * dzl * dwl * table[ir][jr][kl][ll] + dxr * dyr * dzl * dwr * table[ir][jr][kl][lr] +
         dxr * dyr * dzr * dwl * table[ir][jr][kr][ll] + dxr * dyr * dzr * dwr * table[ir][jr][kr][lr];
}

double get_spectrum(double e, int i_age, int i_metallicity);

/*Required structures*/

extern struct rt_face_data
{
  struct geometry geom;
  int face_active;
  int face_responsibility;
  double face_timestep;
  double vel_face[3];
  double vel_face_turned[3];
  // double Ldx , Ldy, Ldz, Rdx, Rdy, Rdz ;
} * rtfacedata;

#endif

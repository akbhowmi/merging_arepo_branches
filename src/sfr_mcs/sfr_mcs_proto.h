/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/sfr_mcs/sfr_mcs_proto.h
 * \date        10/2015
 * \author     	Matthew C Smith
 * \brief
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#ifndef SFR_MCS_PROTO_H
#define SFR_MCS_PROTO_H

void init_star_formation(void);

#ifdef SFR_MCS_LOG
void setup_sf_log(void);
void sf_add_to_log(double rho);
void write_sf_dens_log(void);
#endif

#ifdef GFM_COOLING_METAL

#define GFM_MIN_METAL -20 /* This is otherwise defined in stellar_evolution_vars.h */

void metal_cool_mcs_init(double localMetallicityFloor, double localLogMetallicityFloorInSolar);
#ifdef UVB_SELF_SHIELDING
void update_radiation_state(MyFloat rho, MyFloat metallicity, PhotoCurrent *localpc);
#endif
void update_gas_state(MyFloat rho, MyFloat metallicity, GasState *localgs);

/* This function declared in stellar_evolution_proto.h but needed in cooling_metal.c, so reproduced here. */
static inline MyFloat interpol_4d(MyFloat ****table, int i, int j, int k, int l, double dx, double dy, double dz, double dw)
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
#endif

#ifdef SN_MCS
void init_sne(void);
void read_sb99_tables(void);

#ifndef SB99_FIXED_Z
int get_sb99_z_index(MyFloat metallicity);
#endif
void get_sb99_t_indicies(MyFloat t, MyFloat *timesteps, int N_t, int *it_low, int *it_high, MyFloat *delta);

#ifndef SB99_FIXED_Z
double get_sn_rate(MyFloat t, MyFloat metallicity);
#else
double get_sn_rate(MyFloat t);
#endif

void check_for_supernovae(int *local_n_sn_event, int *global_n_sn_event);
void find_sn_host_cells_and_distribute(int n_sn_event);
void stellar_feedback_distribute(void);
void sn_finish(void);
void do_supernovae(void);

#ifdef SN_MCS_LOG
void setup_sn_log(void);
void sn_add_to_log(int i);
void write_sn_dens_log(void);
#endif  // SN_MCS_LOG

#endif  // SN_MCS
#endif  // SFR_MCS
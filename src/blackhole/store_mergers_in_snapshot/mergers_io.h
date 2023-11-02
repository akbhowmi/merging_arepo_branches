/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/fof/fof.h
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

#include "../../allvars.h"

extern int Nmergers, TotNmergers;

extern struct merger_properties
{
  long long ID1;
  long long ID2;
  int TaskID;
  float Time;
  MyFloat BH_Mass1;
  MyFloat BH_Mass2;
#ifdef OUTPUT_HOST_PROPERTIES_FOR_BH_MERGERS 
  MyFloat BH_Mdot1;
  MyFloat BH_Mdot2;
  MyFloat BH_HostHaloMass1;
  MyFloat BH_HostHaloMass2;
  MyFloat BH_HostStellarMass1;
  MyFloat BH_HostStellarMass2; 
  MyFloat BH_HostGasMass1;
  MyFloat BH_HostGasMass2;
  MyFloat BH_HostSFR1;
  MyFloat BH_HostSFR2;
#endif
} * MergerEvents, *DomainMergerBuf;


extern struct merger_header
{
  int Nmergers;
  int TotNmergers;
  int num_files;
  double time;
  double redshift;
  double HubbleParam;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  int flag_doubleprecision;
} merger_catalogue_header;

enum merger_iofields
{
  IO_MERGER_ID1,
  IO_MERGER_ID2,
  IO_MERGER_TASKID,
  IO_MERGER_TIME,
  IO_MERGER_BHMASS1,
  IO_MERGER_BHMASS2,
#ifdef OUTPUT_HOST_PROPERTIES_FOR_BH_MERGERS
  IO_MERGER_BHMDOT1,
  IO_MERGER_BHMDOT2,
  IO_MERGER_BHHOSTHALOMASS1,
  IO_MERGER_BHHOSTHALOMASS2,
  IO_MERGER_BHHOSTSTELLARMASS1,
  IO_MERGER_BHHOSTSTELLARMASS2,
  IO_MERGER_BHHOSTGASMASS1,
  IO_MERGER_BHHOSTGASMASS2,
  IO_MERGER_BHHOSTSFR1,
  IO_MERGER_BHHOSTSFR2,
#endif
  IO_MERGER_LASTENTRY
};

void merger_fill_write_buffer(enum merger_iofields blocknr, int *startindex, int pc);
int merger_blockpresent(enum merger_iofields blocknr);
int merger_get_datatype(enum merger_iofields blocknr);
int merger_get_bytes_per_blockelement(enum merger_iofields blocknr);
int merger_get_particles_in_block(enum merger_iofields blocknr);
void merger_get_dataset_name(enum merger_iofields blocknr, char *label);
void merger_get_Tab_IO_Label(enum merger_iofields blocknr, char *label);
int merger_get_dataset_group(enum merger_iofields blocknr);
//void merger_fill_write_buffer(enum merger_iofields blocknr, int *startindex, int pc);
int merger_get_values_per_blockelement(enum merger_iofields blocknr);
//void add_merger_event(long long Merging_ID1 , long long Merging_ID2);

#ifdef OUTPUT_HOST_PROPERTIES_FOR_BH_MERGERS
void add_merger_event(int ThisTask, float Merger_Time, long long Merging_ID1 , MyFloat Merging_BH_Mass1, long long Merging_ID2, MyFloat Merging_BH_Mass2 , MyFloat Merging_BH_Mdot1 , MyFloat Merging_BH_Mdot2 , MyFloat Merging_BH_HostHaloMass1 , MyFloat Merging_BH_HostHaloMass2 , MyFloat Merging_BH_HostStellarMass1 , MyFloat Merging_BH_HostStellarMass2 , MyFloat Merging_BH_HostGasMass1 , MyFloat Merging_BH_HostGasMass2 , MyFloat Merging_BH_HostSFR1 , MyFloat Merging_BH_HostSFR2);
#else
void add_merger_event(int ThisTask, float Merger_Time, long long Merging_ID1 , MyFloat Merging_BH_Mass1, long long Merging_ID2, MyFloat Merging_BH_Mass2);
#endif

void merger_write_file(char *fname, int writeTask, int lastTask);
void merger_prepare_ID_list(void);
void save_mergers(int num);

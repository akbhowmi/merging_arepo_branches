/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/io.c
 * \date        MM/YYYY
 * \author
 * \brief       Output of a snapshot (or an image file) file to disk.
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include <errno.h>
#include <math.h>
#include <mpi.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#include "allvars.h"
#include "fof/fof.h"
#include "gitversion/version.h"
#include "proto.h"
#include "voronoi.h"

#ifdef HAVE_HDF5
#include <hdf5.h>
static void write_dataset_attributes(hid_t hdf5_dataset, enum iofields blocknr);
#endif

#ifdef TOLERATE_WRITE_ERROR
static char alternative_fname[MAXLEN_PATH];
#endif

#ifdef OUTPUT_XDMF
/* needs to be included after allvars.h */
#include <libgen.h> /* for basename() function */
static void write_xdmf(char *fname);
#endif /* #ifdef OUTPUT_XDMF */

static int n_type[NTYPES]; /**< contains the local (for a single task) number of particles of each type in the snapshot file */
static long long ntot_type_all[NTYPES]; /**< contains the global number of particles of each type in the snapshot file */
static int subbox_dump = 0;

/*! \brief Function for registering an output field.
 *
 *  Don't forget to add the new IO_FLAG to allvars.h.
 *
 *  \param[in] field Specifies the field as an enumeration type iofields
 *             (allvars.h), e.g. IO_POS. Don't forget to insert new fields
 *             also in allvars.h.
 *  \param[in] label The label of the dataset (4 characters).
 *  \param[in] datasetname The name of the hdf5 dataset (maximum 256
 *             characters).
 *  \param[in] type_in_memory The type of the field in the memory (use
 *             MEM_NONE if specifying io_func).
 *  \param[in] type_in_file_output The output type in the hdf5 file.
 *  \param[in] type_in_file_input The input type in the hdf5 file (use
 *             FILE_MY_OUTPUT_TYPE for MyInputType, input is disabled with
 *             FILE_NONE).
 *  \param[in] values_per_block The number of values per field, e.g. 1 for
 *             mass, 3 for velocities.
 *  \param[in] array The array in which the value is stored. For an io_func
 *             this influences the particle index, the default (A_NONE) is an
 *             index into P/SphP, can be changed to indexes into StarP/BHP/TLL
 *             using A_STARP/A_BHP/A_TLL.
 *  \param[in] pointer_to_field A Pointer to the field in one of the global
 *             arrays, e.g. &SphP[0].Density, or &P[0].Vel[0].
 *  \param[in] io_func Alternatively, if the value to output/input is not a
 *             simple field, you can define a function which handles i/o.
 *  \param[in] typelist_bitmask Specifies for which particle type the field is
 *             present, e.g. 1+2+8 => field present for particle types 0,1,3
 *             (or use ALL_TYPES, GAS_ONLY,...).
 *
 *  \return void
 */
void init_field(enum iofields field, const char *label, const char *datasetname, enum types_in_memory type_in_memory,
                enum types_in_file type_in_file_output, enum types_in_file type_in_file_input, int values_per_block, enum arrays array,
                void *pointer_to_field, void (*io_func)(int, int, void *, int), int typelist_bitmask)
{
  int alloc_step = 5;

  if(Max_IO_Fields == 0)
    {
      IO_Fields     = (IO_Field *)mymalloc("IO_Fields", alloc_step * sizeof(IO_Field));
      Max_IO_Fields = alloc_step;
    }
  else if(Max_IO_Fields == N_IO_Fields)
    {
      Max_IO_Fields = ((Max_IO_Fields / alloc_step) + 1) * alloc_step;
      IO_Fields     = (IO_Field *)myrealloc(IO_Fields, Max_IO_Fields * sizeof(IO_Field));
    }

  IO_Fields[N_IO_Fields].field = field;
  memcpy(IO_Fields[N_IO_Fields].label, label, IO_LABEL_SIZE);
  strncpy(IO_Fields[N_IO_Fields].datasetname, datasetname, IO_DATASET_NAME_SIZE);
  IO_Fields[N_IO_Fields].datasetname[IO_DATASET_NAME_SIZE - 1] = '\0';
  IO_Fields[N_IO_Fields].type_in_memory                        = type_in_memory;
  IO_Fields[N_IO_Fields].type_in_file_output                   = type_in_file_output;
  IO_Fields[N_IO_Fields].type_in_file_input                    = type_in_file_input;
  IO_Fields[N_IO_Fields].values_per_block                      = values_per_block;
  IO_Fields[N_IO_Fields].snap_type                             = SN_FULL;
  IO_Fields[N_IO_Fields].typelist                              = typelist_bitmask;

  IO_Fields[N_IO_Fields].array = array;

  if(array == A_NONE)
    {
      IO_Fields[N_IO_Fields].offset = 0;
    }
  else if(array == A_SPHP)
    {
      IO_Fields[N_IO_Fields].offset = (size_t)pointer_to_field - (size_t)SphP;
    }
  else if(array == A_P)
    {
      IO_Fields[N_IO_Fields].offset = (size_t)pointer_to_field - (size_t)P;
    }
#ifdef TRACER_MC
  else if(array == A_TLL)
    {
      IO_Fields[N_IO_Fields].offset = (size_t)pointer_to_field - (size_t)TracerLinkedList;
    }
#endif
#if defined(GFM) || defined(SFR_MCS)
  else if(array == A_STARP)
    {
      IO_Fields[N_IO_Fields].offset = (size_t)pointer_to_field - (size_t)StarP;
    }
#endif
#ifdef BLACK_HOLES
  else if(array == A_BHP)
    {
      IO_Fields[N_IO_Fields].offset = (size_t)pointer_to_field - (size_t)BHP;
    }
#endif
#ifdef DUST_LIVE
  else if(array == A_DUSTP)
    {
      IO_Fields[N_IO_Fields].offset = (size_t)pointer_to_field - (size_t)DustP;
    }
#endif
  else if(array == A_PS)
    {
      IO_Fields[N_IO_Fields].offset = (size_t)pointer_to_field - (size_t)PS;
    }

  IO_Fields[N_IO_Fields].io_func = io_func;

  /* validate types */
  if(type_in_memory == MEM_INT &&
     ((type_in_file_input != FILE_NONE && type_in_file_input != FILE_INT) || type_in_file_output != FILE_INT))
    {
      terminate("combination of datatypes not supported (field %s)", datasetname);
    }

  if(type_in_memory == MEM_MY_ID_TYPE &&
     ((type_in_file_input != FILE_NONE && type_in_file_input != FILE_MY_ID_TYPE) || type_in_file_output != FILE_MY_ID_TYPE))
    {
      terminate("combination of datatypes not supported (field %s)", datasetname);
    }

  if((type_in_memory == MEM_FLOAT || type_in_memory == MEM_MY_SINGLE || type_in_memory == MEM_DOUBLE) &&
     ((type_in_file_input != FILE_NONE && (type_in_file_input == FILE_MY_ID_TYPE || type_in_file_input == FILE_INT)) ||
      type_in_file_output == FILE_INT || type_in_file_output == FILE_MY_ID_TYPE))
    {
      terminate("combination of datatypes not supported (field %s)", datasetname);
    }

  IO_Fields[N_IO_Fields].a       = 0.;
  IO_Fields[N_IO_Fields].h       = 0.;
  IO_Fields[N_IO_Fields].L       = 0.;
  IO_Fields[N_IO_Fields].M       = 0.;
  IO_Fields[N_IO_Fields].V       = 0.;
  IO_Fields[N_IO_Fields].c       = 0.;
  IO_Fields[N_IO_Fields].hasunit = 0;

  N_IO_Fields++;
}

/*! \brief Function for adding units to output field.
 *
 *  This only works for fields registered with init_field.
 *
 *  \param[in] field Specifies the field as an enumeration type iofields
 *             (allvars.h), e.g. IO_POS.
 *  \param[in] a the exponent of the cosmological a factor.
 *  \param[in] h the exponent of the hubble parameter.
 *  \param[in] L the length unit scaling.
 *  \param[in] M the mass unit scaling.
 *  \param[in] V the velocity unit scaling.
 *  \param[in] c conversion factor to cgs units (zero indicates dimensionless
 *             quantity, integer count, etc).
 *
 *  \return void
 */
void init_units(const enum iofields field, const double a, const double h, const double L, const double M, const double V,
                const double c)
{
  for(int i = 0; i < N_IO_Fields; i++)
    if(IO_Fields[i].field == field)
      {
        IO_Fields[i].hasunit = 1;
        IO_Fields[i].a       = a;
        IO_Fields[i].h       = h;
        IO_Fields[i].L       = L;
        IO_Fields[i].M       = M;
        IO_Fields[i].V       = V;
        IO_Fields[i].c       = c;
        break;
      }
}

/*! \brief Function for determining whether a field is dumped in snapshot.
 *
 *  This only works for fields registered with init_field.
 *  The member snap_type is initialized to SN_FULL in init_field.
 *
 *  \param[in] field Specifies the field as an enumeration type iofields
 *             (allvars.h), e.g. IO_POS.
 *  \param[in] type In which snapshot types this field should be present
 *             (e.g. SN_FULL).
 *
 *  \return void
 */
void init_snapshot_type(const enum iofields field, const enum sn_type type)
{
  for(int i = 0; i < N_IO_Fields; i++)
    if(IO_Fields[i].field == field)
      IO_Fields[i].snap_type = type;
}

#ifdef TOLERATE_WRITE_ERROR
/*! \brief Print information about a write error.
 *
 *  If a write error occurs, this function prints some useful debug information
 *  and sets to 1 the variable WriteErrorFlag so that the write operation that
 *  caused the error can be performed again.
 *
 *  \param[in] check Flag that indicates where the function was called [0 and 1
 *             in my_fwrite(), 2 in my_hdf5_error_handler(), 3 in
 *             hdf5_header_error_handler()].
 *  \param[in] nwritten Number of elements actually written.
 *  \param[in] nmemb Number of elements that should be written.
 *
 *  \return void
 */
void write_error(const int check, const size_t nwritten, const size_t nmemb)
{
  if(!WriteErrorFlag)
    {
      int len;
      char hostname[MPI_MAX_PROCESSOR_NAME];
      MPI_Get_processor_name(hostname, &len);

      printf("TOLERATE_WRITE_ERROR: write failed node=%s  nwritten=%lld  nmemb=%lld  errno=%s  task=%d  check=%d\n", hostname,
             (long long)nwritten, (long long)nmemb, strerror(errno), ThisTask, check);
      myflush(stdout);
      WriteErrorFlag = 1;
    }
}
#endif

/*! \brief Checks if a snapshot should be saved.
 *
 *  This function checks whether a snapshot file or other kinds of output
 *  files, such as a projection, should be saved at the current time-step.
 *  If that is the case, the appropriate functions to produce the desired
 *  file are called and the parameter controlling the output are updated
 *  accordingly.
 *
 *  \return void
 */
void create_snapshot_if_desired(void)
{
  WroteSnapThisTimestep = 0;

#ifdef OUTPUT_EVERY_STEP
  All.Ti_nextoutput = All.Ti_Current;
#endif

#ifdef SUBBOX_SNAPSHOTS
  FILE *fd;
  if(All.Time <= All.SubboxMaxTime && All.Time >= All.SubboxMinTime)
    if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin)
      {
        if(All.SubboxSyncCounter % All.SubboxSyncModulo == 0)
          {
            int tmp_1 = All.NumFilesPerSnapshot, tmp_2 = All.NumFilesWrittenInParallel;

            All.NumFilesPerSnapshot       = All.SubboxNumFilesPerSnapshot;
            All.NumFilesWrittenInParallel = All.SubboxNumFilesWrittenInParallel;

            if(ThisTask == 0)
              {
                fd = fopen("subboxes_times.txt", "a");
                fprintf(fd, "%d %g\n", All.SubboxSnapshotFileCount, All.Time);
                fclose(fd);
              }

            mpi_printf("SUBBOX_SNAPSHOTS: Writing subboxes snapshots %03d (%06d)\n", All.SubboxSnapshotFileCount,
                       All.SubboxSyncCounter);

            DumpFlag = DUMP_BOTH; /* make subboxes full snapshots */

            for(int isubbox = 1; isubbox <= Nsubboxes; isubbox++)
              savepositions(All.SubboxSnapshotFileCount, isubbox);

            All.SubboxSnapshotFileCount++;
            All.NumFilesPerSnapshot       = tmp_1;
            All.NumFilesWrittenInParallel = tmp_2;
          }
        All.SubboxSyncCounter++;
      }
#endif

#ifdef HIGH_FREQUENCY_OUTPUT_STARS
  if(All.HighFreqStarsSnapshotCount < All.HighFreqStarsSnapshotNum &&
     All.Time > All.HighFreqStarsOutputTimes[All.HighFreqStarsSnapshotCount])
    {
      for(int p = 0; p < NumPart; p++)
        {
          if(P[p].Type == PTYPE_STARS)
            drift_particle(p, All.Ti_Current);
        }

      mpi_printf("HIGH_FREQUENCY_OUTPUT_STARS: Writing star snapshot %03d at t=%g)\n", All.HighFreqStarsSnapshotCount,
                 All.HighFreqStarsOutputTimes[All.HighFreqStarsSnapshotCount]);

      savepositions(All.HighFreqStarsSnapshotCount, 42);

      All.HighFreqStarsSnapshotCount++;
    }
#endif

#ifdef TRACER_TRAJECTORY
  tracer_write_output_if_needed();
#endif

#ifdef SINK_PARTICLES_OUTPUT_EVERY_NEW_SINK
  if(NSinksLastSnapshot < NSinksAllTasks)
    {
      All.Ti_nextoutput = All.Ti_Current;
      mpi_printf(
          "SINK_PARTICLES_OUTPUT_EVERY_NEW_SINK: A new sink was formed sind the last global timestep, we now have %03d sinks and set "
          "the next output time to now, waiting for the global timestep to write snapshot\n",
          NSinksAllTasks);
    }
  if(NSinksLastSnapshot > NSinksAllTasks)
    {
      NSinksLastSnapshot = NSinksAllTasks;
      mpi_printf("SINK_PARTICLES_OUTPUT_EVERY_NEW_SINK: A new sink was lost, we are reducing NSinksLastSnapshot \n", NSinksAllTasks);
    }

#endif

  if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin) /* allow only top-level synchronization points */
    if(All.Ti_Current >= All.Ti_nextoutput && All.Ti_nextoutput >= 0)
      {
        DumpFlag = DumpFlagNextSnap;
        produce_dump();

        All.Ti_nextoutput = find_next_outputtime(All.Ti_Current + 1);

#ifdef SINK_PARTICLES_OUTPUT_EVERY_NEW_SINK
        NSinksLastSnapshot = NSinksAllTasks;
#endif

        WroteSnapThisTimestep = 1;
      }

#ifdef VORONOI_FREQUENT_IMAGES
  if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin) /* allow only top-level synchronization points */
    {
      if(All.ImageTimeBin == 0) /* on start-up, this is initialized to zero, and so we calculate the value once */
        {
          integertime ti_image = (All.TimeBetweenImages / All.Timebase_interval);

          /* make it a power 2 subdivision */
          integertime ti_min = TIMEBASE;
          while(ti_min > ti_image)
            ti_min >>= 1;
          ti_image = ti_min;

          All.ImageTimeBin = get_timestep_bin(ti_image);
        }

      if(TimeBinSynchronized[All.ImageTimeBin])
        {
          exchange_primitive_variables();
          calculate_gradients();
          exchange_primitive_variables_and_gradients();
          make_voronoi_image(All.ImageCount);
          make_voronoi_image_slice(All.ImageCount);
          All.ImageCount++;
        }
    }
#endif
}

/*! \brief A wrapper function used to create a snapshot.
 *
 *  This function wraps together savepositions(), the function that
 *  saves the snapshot file to the disk, with functions used for
 *  special output needs.
 *
 *  \return void
 */
void produce_dump(void)
{
#ifdef UPDATE_GRADIENTS_FOR_OUTPUT
  exchange_primitive_variables();
  calculate_gradients();
#endif

#ifdef SHOCK_FINDER_BEFORE_OUTPUT
  shock_finder_on_the_fly();
#endif

#ifdef SINK_PARTICLES
  dump_sink_particle_snapshot(All.SnapshotFileCount);
#endif

  savepositions(All.SnapshotFileCount++, 0); /* write snapshot file */

#ifdef BH_NEW_LOGS
  close_bh_logfiles();   // Close the BH logfiles for the previous snapshot
  if(All.Ti_Current < TIMEBASE) // Only create new logfiles if we're not at the final time
    open_bh_logfiles();    // Open new BH logfiles for the next snapshot
#endif

#ifdef POWERSPECTRUM_ON_THE_FLY
  {
    int typeflag[NTYPES];
    for(int type = 0; type < NTYPES; type++)
      typeflag[type] = 1;
    calculate_power_spectra_and_ntot(All.SnapshotFileCount - 1, typeflag);
  }
#endif

#ifdef VORONOI_FIELD_DUMP_PIXELS_X
  exchange_primitive_variables();
  calculate_gradients();
  exchange_primitive_variables_and_gradients();
  do_special_dump(All.SnapshotFileCount - 1, 1);
#ifdef VORONOI_NOGRADS
  do_special_dump(All.SnapshotFileCount - 1, 0);
#endif
#endif

#ifdef VORONOI_IMAGES_FOREACHSNAPSHOT
  exchange_primitive_variables();
  calculate_gradients();
  exchange_primitive_variables_and_gradients();
  make_voronoi_image(All.SnapshotFileCount - 1);
  make_voronoi_image_slice(All.SnapshotFileCount - 1);
#endif

#ifdef SFR_MCS_LOG
  write_sf_dens_log();
#endif
#ifdef SN_MCS_LOG
  write_sn_dens_log();
#endif
}

/*! \brief Saves snapshot to disk.
 *
 *  This function writes a snapshot of the particle distribution to one or
 *  several files. If NumFilesPerSnapshot>1, the snapshot is distributed
 *  into several files, which are written simultaneously. Each file contains
 *  data from a group of processors of size roughly NTask/NumFilesPerSnapshot.
 *
 *  \param[in] num The snapshot number.
 *  \param[in] subbox_flag If greater than 0 instructs the code to output only
 *             a subset of the whole domain.
 *
 *  \return void
 */
void savepositions(const int num, const int subbox_flag)
{
  char buf[MAXLEN_PATH];

  const double t0 = second();
  CPU_Step[CPU_MISC] += measure_time();

  if(DumpFlag)
    {
      subbox_dump = 0;

      if(subbox_flag > 0)
        {
          mpi_printf("\nwriting small subbox #%d snapshot file #%d @ time %g ... \n", subbox_flag - 1, num, All.Time);
          subbox_dump = 1;
        }
      else
        mpi_printf("\nwriting snapshot file #%d @ time %g ... (DumpFlag=%d)\n", num, All.Time, DumpFlag);

#ifdef FOF
      if(RestartFlag != RESTART_FOF_SUBFIND && RestartFlag != RESTART_SHOCK_FINDER && RestartFlag != RESTART_RECALC_POTENTIAL &&
         RestartFlag != RESTART_CALC_VORONOI_DM_DENSITY && subbox_flag == 0 && DumpFlag != DUMP_ONLY_SNAP)
        {
#ifdef TGSET
          if(TGD.NHMax < 1e4)
#endif
            {
              mpi_printf("\nWe shall first compute a group catalogue for this snapshot file\n");
              fof_fof(num);
            }
        }
#endif

      int ngroups = All.NumFilesPerSnapshot / All.NumFilesWrittenInParallel;
      if(All.NumFilesPerSnapshot % All.NumFilesWrittenInParallel)
        ngroups++;

      /* assign processors to output files */
      int filenr, masterTask, lastTask;
      distribute_file(All.NumFilesPerSnapshot, &filenr, &masterTask, &lastTask);

      if(DumpFlag != DUMP_ONLY_HALOS)
        {
          CommBuffer = (void *)mymalloc("CommBuffer", COMMBUFFERSIZE);

#ifndef HAVE_HDF5
          if(All.SnapFormat == SNAP_FORMAT_HDF5)
            mpi_terminate("Code wasn't compiled with HDF5 support enabled!");
#endif
          /* determine global and local particle numbers */
          for(int n = 0; n < NTYPES; n++)
            n_type[n] = 0;

          for(int n = 0; n < NumPart; n++)
            {
#ifdef SUBBOX_SNAPSHOTS
              if(subbox_flag > 0)
                if(P[n].Pos[0] < SubboxXmin[subbox_flag - 1] || P[n].Pos[0] > SubboxXmax[subbox_flag - 1] ||
                   P[n].Pos[1] < SubboxYmin[subbox_flag - 1] || P[n].Pos[1] > SubboxYmax[subbox_flag - 1] ||
                   P[n].Pos[2] < SubboxZmin[subbox_flag - 1] || P[n].Pos[2] > SubboxZmax[subbox_flag - 1])
                  continue;
#endif
#ifdef HCOUTPUT
              if(subbox_flag > 0)
                {
                  double prad = 0.;
                  for(int k = 0; k < 3; k++)
                    prad += (P[n].Pos[k] - All.HCOutput_Halo_Center[k]) * (P[n].Pos[k] - All.HCOutput_Halo_Center[k]);
                  prad = pow(prad, 0.5);
                  if(prad > All.HCOutput_RadialCut)
                    continue;
                  else if(!((1 << P[n].Type) & (HCOUTPUT)))
                    continue;
                }
#endif

#ifdef HIGH_FREQUENCY_OUTPUT_STARS
              if(subbox_flag > 0)
                if(P[n].Type != PTYPE_STARS || StarP[P[n].AuxDataID].BirthTime <= 0)
                  continue;
#endif
              n_type[P[n].Type]++;
            }

#ifdef TRACER_MC
          if(n_type[TRACER_MC])
            terminate("ERROR! Output particle type specified for TRACER_MC is already in use!");

          n_type[TRACER_MC] = N_tracer;

#ifdef SUBBOX_SNAPSHOTS
          if(subbox_flag > 0)
            {
              n_type[TRACER_MC] = 0;

              for(int n = 0; n < NumPart; n++)
                {
                  if(P[n].Pos[0] < SubboxXmin[subbox_flag - 1] || P[n].Pos[0] > SubboxXmax[subbox_flag - 1] ||
                     P[n].Pos[1] < SubboxYmin[subbox_flag - 1] || P[n].Pos[1] > SubboxYmax[subbox_flag - 1] ||
                     P[n].Pos[2] < SubboxZmin[subbox_flag - 1] || P[n].Pos[2] > SubboxZmax[subbox_flag - 1])
                    continue;

                  n_type[TRACER_MC] += P[n].NumberOfTracers;
                }
            }
#endif

#ifdef HIGH_FREQUENCY_OUTPUT_STARS
          if(subbox_flag > 0)
            n_type[TRACER_MC] = 0;
#endif

          tracer_cellids = (MyIDType *)mymalloc("tracer_cellids", sizeof(MyIDType) * All.MaxPartTracer);
          memset(tracer_cellids, 0, sizeof(MyIDType) * All.MaxPartTracer);
          for(int i = 0; i < NumPart; i++)
            {
              int next = P[i].TracerHead;
              while(next >= 0)
                {
                  tracer_cellids[next] = P[i].ID;

#ifdef SUBBOX_SNAPSHOTS
                  if(subbox_flag > 0)
                    if(P[i].Pos[0] < SubboxXmin[subbox_flag - 1] || P[i].Pos[0] > SubboxXmax[subbox_flag - 1] ||
                       P[i].Pos[1] < SubboxYmin[subbox_flag - 1] || P[i].Pos[1] > SubboxYmax[subbox_flag - 1] ||
                       P[i].Pos[2] < SubboxZmin[subbox_flag - 1] || P[i].Pos[2] > SubboxZmax[subbox_flag - 1])
                      tracer_cellids[next] = 0;
#endif

#ifdef HIGH_FREQUENCY_OUTPUT_STARS
                  if(subbox_flag > 0)
                    tracer_cellids[next] = 0;
#endif
                  next = TracerLinkedList[next].Next;
                }
            }

#endif

#ifdef GFM_WINDS_SAVE_PARTTYPE
          /* write wind particles as a different particle type, instead of the usual 4 (stars) */
          if(n_type[GFM_WINDS_SAVE_PARTTYPE])
            terminate("ERROR! Output particle type specified for GFM_WINDS_SAVE_PARTTYPE is already in use!");

          int N_wind = 0;

          for(int n = 0; n < NumPart; n++)
            if(P[n].Type == PTYPE_STARS && STP(n).BirthTime <= 0)
              {
                N_wind++;
                P[n].Type = GFM_WINDS_SAVE_PARTTYPE;
              }

          n_type[GFM_WINDS_SAVE_PARTTYPE] = N_wind;
          n_type[PTYPE_STARS] -= N_wind;
#endif

          sumup_large_ints(NTYPES, n_type, ntot_type_all);

          mpi_printf("All.NumFilesPerSnapshot=%d, subbox_flag=%d\n", All.NumFilesPerSnapshot, subbox_flag);
          if(All.NumFilesPerSnapshot > 1)
            {
              if(ThisTask == 0)
                {
                  file_path_sprintf(buf, "%s/snapdir_%03d", All.OutputDir, num);
#ifdef SUBBOX_SNAPSHOTS
                  if(subbox_flag > 0)
                    {
                      file_path_sprintf(buf, "%s/subbox%d", All.OutputDir, subbox_flag - 1);
                      mkdir(buf, MKDIR_MODE);
                      file_path_sprintf(buf, "%s/subbox%d/snapdir_subbox%d_%03d", All.OutputDir, subbox_flag - 1, subbox_flag - 1,
                                        num);
                    }
#endif
#ifdef HIGH_FREQUENCY_OUTPUT_STARS
                  if(subbox_flag > 0)
                    file_path_sprintf(buf, "%s/snapdir_highfreqstars_%03d", All.OutputDir, num);
#endif
#ifdef HCOUTPUT
                  if(subbox_flag == 0)
                    mkdir(buf, MKDIR_MODE);
#endif
#ifndef HCOUTPUT
                  mkdir(buf, MKDIR_MODE);
#endif

#ifdef HCOUTPUT
                  if(subbox_flag > 0)
                    {
                      mkdir(All.HCOutput_Directory, MKDIR_MODE);
                      file_path_sprintf(buf, "%s/snipdir_%03d", All.HCOutput_Directory, num);
                      mkdir(buf, MKDIR_MODE);
                    }
#endif

#ifdef TOLERATE_WRITE_ERROR
                  file_path_sprintf(alternative_fname, "%s/snapdir_%03d", AlternativeOutputDir, num);
#ifdef SUBBOX_SNAPSHOTS
                  if(subbox_flag > 0)
                    {
                      file_path_sprintf(alternative_fname, "%s/subbox%d", AlternativeOutputDir, subbox_flag - 1);
                      mkdir(alternative_fname, MKDIR_MODE);
                      file_path_sprintf(alternative_fname, "%s/subbox%d/snapdir_subbox%d_%03d", AlternativeOutputDir, subbox_flag - 1,
                                        subbox_flag - 1, num);
                    }
#endif
                  mkdir(alternative_fname, MKDIR_MODE);
#endif
                }

              MPI_Barrier(MPI_COMM_WORLD);
            }

          if(All.NumFilesPerSnapshot > 1)
            file_path_sprintf(buf, "%s/snapdir_%03d/%s_%03d.%d", All.OutputDir, num, All.SnapshotFileBase, num, filenr);
          else
            file_path_sprintf(buf, "%s/%s_%03d", All.OutputDir, All.SnapshotFileBase, num);

#ifdef TOLERATE_WRITE_ERROR
          if(All.NumFilesPerSnapshot > 1)
            file_path_sprintf(alternative_fname, "%s/snapdir_%03d/%s_%03d.%d", AlternativeOutputDir, num, All.SnapshotFileBase, num,
                              filenr);
          else
            file_path_sprintf(alternative_fname, "%s/%s_%03d", AlternativeOutputDir, All.SnapshotFileBase, num);
#endif

#ifdef SUBBOX_SNAPSHOTS
          if(subbox_flag > 0)
            {
              if(All.NumFilesPerSnapshot > 1)
                file_path_sprintf(buf, "%s/subbox%d/snapdir_subbox%d_%03d/%s_subbox%d_%03d.%d", All.OutputDir, subbox_flag - 1,
                                  subbox_flag - 1, num, All.SnapshotFileBase, subbox_flag - 1, num, filenr);
              else
                file_path_sprintf(buf, "%s/%s_subbox%d_%03d", All.OutputDir, All.SnapshotFileBase, subbox_flag - 1, num);
            }
#ifdef TOLERATE_WRITE_ERROR
          if(subbox_flag > 0)
            {
              if(All.NumFilesPerSnapshot > 1)
                file_path_sprintf(alternative_fname, "%s/subbox%d/snapdir_subbox%d_%03d/%s_subbox%d_%03d.%d", AlternativeOutputDir,
                                  subbox_flag - 1, subbox_flag - 1, num, All.SnapshotFileBase, subbox_flag - 1, num, filenr);
              else
                file_path_sprintf(alternative_fname, "%s/%s_subbox%d_%03d", AlternativeOutputDir, All.SnapshotFileBase,
                                  subbox_flag - 1, num);
            }
#endif
#endif

#ifdef HIGH_FREQUENCY_OUTPUT_STARS
          if(subbox_flag > 0)
            {
              if(All.NumFilesPerSnapshot > 1)
                file_path_sprintf(buf, "%s/snapdir_highfreqstars_%03d/%s_highfreqstars_%03d.%d", All.OutputDir, num,
                                  All.SnapshotFileBase, num, filenr);
              else
                file_path_sprintf(buf, "%s/%s_highfreqstars_%03d", All.OutputDir, All.SnapshotFileBase, num);
            }
#endif

#ifdef HCOUTPUT
          if(subbox_flag > 0)
            {
              if(All.NumFilesPerSnapshot > 1)
                file_path_sprintf(buf, "%s/snipdir_%03d/snipshot_%03d.%d", All.HCOutput_Directory, num, num, filenr);
              else
                file_path_sprintf(buf, "%s/snipshot_%03d", All.HCOutput_Directory, num);
            }
#endif

          if(RestartFlag == RESTART_FOF_SUBFIND)
            {
#ifndef FOF_STOREIDS
              if(All.NumFilesPerSnapshot > 1)
                file_path_sprintf(buf, "%s/snapdir_%03d/%s-groupordered_%03d.%d", All.OutputDir, num, All.SnapshotFileBase, num,
                                  filenr);
              else
                file_path_sprintf(buf, "%s/%s-groupordered_%03d", All.OutputDir, All.SnapshotFileBase, num);
#else
              if(All.NumFilesPerSnapshot > 1)
                file_path_sprintf(buf, "%s/snapdir_%03d/%s-storeids_%03d.%d", All.OutputDir, num, All.SnapshotFileBase, num, filenr);
              else
                file_path_sprintf(buf, "%s/%s-storeids_%03d", All.OutputDir, All.SnapshotFileBase, num);
#endif
            }

          if(RestartFlag == RESTART_TRACER_POWER_SPECTRA)
            {
#ifdef TRACER_PARTICLE
              if(All.NumFilesPerSnapshot > 1)
                file_path_sprintf(buf, "%s/snapdir_%03d/%s-TRpart_%03d.%d", All.OutputDir, num, All.SnapshotFileBase, num, filenr);
              else
                file_path_sprintf(buf, "%s/%s-TRpart_%03d", All.OutputDir, All.SnapshotFileBase, num);
#endif
            }

#ifdef ADDBACKGROUNDGRID
          if(All.NumFilesPerSnapshot > 1)
            file_path_sprintf(buf, "%s-with-grid.%d", All.InitCondFile, filenr);
          else
            file_path_sprintf(buf, "%s-with-grid", All.InitCondFile);
#endif

          for(int gr = 0; gr < ngroups; gr++)
            {
              if(filenr / All.NumFilesWrittenInParallel == gr) /* ok, it's this processor's turn */
                {
                  if(ThisTask == masterTask && filenr % All.NumFilesWrittenInParallel == 0)
                    {
                      printf("writing snapshot files group %d out of %d - files %d-%d (total of %d files): '%s'\n", gr + 1, ngroups,
                             filenr, imin(filenr + All.NumFilesWrittenInParallel, All.NumFilesPerSnapshot) - 1,
                             All.NumFilesPerSnapshot, buf);
                      myflush(stdout);
                    }
                  write_file(buf, masterTask, lastTask, subbox_flag);
#ifdef OUTPUT_XDMF
                  if(All.SnapFormat == SNAP_FORMAT_HDF5)
                    write_xdmf(buf);
#endif
                }
              MPI_Barrier(MPI_COMM_WORLD);
            }

#ifdef GFM_WINDS_SAVE_PARTTYPE
          /* restore wind particles to the usual parttype 4 (stars) */
          for(int n = 0; n < NumPart; n++)
            if(P[n].Type == GFM_WINDS_SAVE_PARTTYPE)
              P[n].Type = PTYPE_STARS;
#endif

#ifdef TRACER_MC
          myfree(tracer_cellids);
#endif

          myfree(CommBuffer);
          CommBuffer = NULL;

          const double t1 = second();
          CPU_Step[CPU_SNAPSHOT] += measure_time();

          mpi_printf("done with writing snapshot (took %g sec).\n", timediff(t0, t1));
        }
      else
        {
          mpi_printf("done with writing files: no dump of snapshot (DumpFlag = %d).\n", DumpFlag);
        } /* if(DumpFlag != DUMP_ONLY_HALOS) */

#ifdef FOF
      if(RestartFlag != RESTART_FOF_SUBFIND && RestartFlag != RESTART_SHOCK_FINDER && RestartFlag != RESTART_SNAP_CONVERSION &&
         RestartFlag != RESTART_RECALC_POTENTIAL && subbox_flag == 0 && DumpFlag != DUMP_ONLY_SNAP)
        {
#ifdef TGSET
          if(TGD.NHMax < 1e4)
#endif
            {
#ifndef FOF_STOREIDS
              /* now revert from output order to the original order */
              for(int n = 0; n < NumPart; n++)
                {
                  PS[n].TargetTask  = PS[n].OriginTask;
                  PS[n].TargetIndex = PS[n].OriginIndex;
                }

              fof_subfind_exchange(MPI_COMM_WORLD);

              myfree(PS);

              /* do resize because subfind may have increased these limits */
              if(All.MaxPart != fof_OldMaxPart)
                {
                  All.MaxPart = fof_OldMaxPart;
                  reallocate_memory_maxpart_ignore_timebins();
                }
              if(All.MaxPartSph != fof_OldMaxPartSph)
                {
                  All.MaxPartSph = fof_OldMaxPartSph;
                  reallocate_memory_maxpartsph_ignore_timebins();
                }
#if defined(GFM) || defined(SFR_MCS)
              if(All.MaxPartStar != fof_OldMaxPartStar)
                {
                  All.MaxPartStar = fof_OldMaxPartStar;
                  reallocate_memory_maxpartstar();
                }
#endif
#ifdef BLACK_HOLES
              if(All.MaxPartBHs != fof_OldMaxPartBHs)
                {
                  All.MaxPartBHs = fof_OldMaxPartBHs;
                  reallocate_memory_maxpartBHs();
                }
#endif
#ifdef DUST_LIVE
              if(All.MaxPartDust != fof_OldMaxPartDust)
                {
                  All.MaxPartDust = fof_OldMaxPartDust;
                  reallocate_memory_maxpartdust();
                }
#endif

              CPU_Step[CPU_FOF] += measure_time();
#endif

              /* recreate the mesh that we had free to reduce peak memory usage */
              create_mesh();
              mesh_setup_exchange();
            }
        }
#endif

#ifdef VORONOI_MESHOUTPUT /* after mesh is recreated */
      for(int gr = 0; gr < ngroups; gr++)
        {
          if(All.NumFilesPerSnapshot > 1)
            file_path_sprintf(buf, "%s/snapdir_%03d/voronoi_mesh_%03d.%d", All.OutputDir, num, num, filenr);
          else
            file_path_sprintf(buf, "%s/voronoi_mesh_%03d", All.OutputDir, num);
          if(RestartFlag == RESTART_TRACER_POWER_SPECTRA)
            {
#ifdef TRACER_PARTICLE
              if(All.NumFilesPerSnapshot > 1)
                file_path_sprintf(buf, "%s/snapdir_%03d/voronoi_mesh-TRpart_%03d.%d", All.OutputDir, num, num, filenr);
              else
                file_path_sprintf(buf, "%s/voronoi_mesh-TRpart_%03d", All.OutputDir, num);
#endif
            }
          write_voronoi_mesh(&Mesh, buf, masterTask, lastTask);
        }
#endif

#if(defined(TRACER_MC) || defined(TRACER_PARTICLE)) && !defined(TRACER_NO_RESET_EACH_SNAP)
      /* reset 'maximum' type tracer properties, for full snapshot dumps only */
      if(subbox_flag == 0 && (DumpFlag == DUMP_BOTH || DumpFlag == DUMP_ONLY_SNAP))
        reset_tracer_parent_fluid_properties();
#endif

      All.Ti_lastoutput = All.Ti_Current;

      CPU_Step[CPU_SNAPSHOT] += measure_time();
    }
}

/*! \brief This function fills the write buffer with particle data.
 *
 *  \param[out] buffer Buffer to be filled.
 *  \param[in] blocknr ID of the output block (i.e. position, velocities...).
 *  \param[in, out] startindex Pointer containing the offset in write buffer.
 *  \param[in] pc Number of particle to be put in the buffer.
 *  \param[in] type Particle type.
 *  \param[in] subbox_flag If greater than 0 instructs the code to output
 *             only a subset of the whole domain.
 *
 *  \return void
 */
void fill_write_buffer(void *const buffer, const enum iofields blocknr, int *const startindex, const int pc, const int type,
                       const int subbox_flag)
{
  /* determine which field we are working on */
  int field = -1;
  for(int f = 0; f < N_IO_Fields; f++)
    {
      if(IO_Fields[f].field == blocknr)
        {
          field = f;
          break;
        }
    }
  if(field < 0)
    terminate("IO field=%d not registered with init_field()", (int)blocknr);

  set_cosmo_factors_for_current_time();

  MyOutputFloat *fp = (MyOutputFloat *)buffer;
  MyIDType *ip      = (MyIDType *)buffer;
  int *intp         = (int *)buffer;
  double *doublep   = (double *)buffer;
  float *floatp     = (float *)buffer;

  int pindex = *startindex;

  for(int n = 0; n < pc; pindex++)
    {
      /* SUBBOX_SNAPSHOTS specialized output */
#ifdef TRACER_MC
      if(type == TRACER_MC)
        if(tracer_cellids[pindex] == 0)
          continue;
#endif

#ifdef SUBBOX_SNAPSHOTS
#ifdef TRACER_MC
      if(type != TRACER_MC)
#endif
        if(subbox_flag > 0)
          if(P[pindex].Pos[0] < SubboxXmin[subbox_flag - 1] || P[pindex].Pos[0] > SubboxXmax[subbox_flag - 1] ||
             P[pindex].Pos[1] < SubboxYmin[subbox_flag - 1] || P[pindex].Pos[1] > SubboxYmax[subbox_flag - 1] ||
             P[pindex].Pos[2] < SubboxZmin[subbox_flag - 1] || P[pindex].Pos[2] > SubboxZmax[subbox_flag - 1])
            continue;
#endif

#ifdef HIGH_FREQUENCY_OUTPUT_STARS
      if(subbox_flag > 0)
#ifdef TRACER_MC
        if(type != TRACER_MC)
#endif
          if(P[pindex].Type != 4 || StarP[P[pindex].AuxDataID].BirthTime <= 0)
            continue;
#endif

#ifdef HCOUTPUT
      if(subbox_flag > 0)
        {
          double prad = 0.;
          for(int k = 0; k < 3; k++)
            prad += (P[pindex].Pos[k] - All.HCOutput_Halo_Center[k]) * (P[pindex].Pos[k] - All.HCOutput_Halo_Center[k]);
          prad = pow(prad, 0.5);
          if(prad > All.HCOutput_RadialCut)
            continue;
          else if(!((1 << P[pindex].Type) & (HCOUTPUT)))
            continue;
        }
#endif

        /* normal particle output */
#ifdef TRACER_MC
      if(type == TRACER_MC || P[pindex].Type == type)
#else
      if(P[pindex].Type == type)
#endif
        {
          if(IO_Fields[field].io_func)
            {
              int particle;
              switch(IO_Fields[field].array)
                {
                  case A_NONE:
                  case A_SPHP:
                  case A_P:
#ifdef TRACER_MC
                  case A_TLL:
#endif
                    particle = pindex;
                    break;
#if defined(GFM) || defined(SFR_MCS)
                  case A_STARP:
                    particle = P[pindex].AuxDataID;
                    break;
#endif
#ifdef BLACK_HOLES
                  case A_BHP:
                    particle = P[pindex].AuxDataID;
                    break;
#endif
#ifdef DUST_LIVE
                  case A_DUSTP:
                    particle = P[pindex].AuxDataID;
                    break;
#endif
                  case A_PS:
                    terminate("Not good, trying to read into PS[]?");
                    break;
                  default:
                    terminate("ERROR in fill_write_buffer: Array not found!");
                    break;
                }

              switch(IO_Fields[field].type_in_file_output)
                {
                  case FILE_NONE:
                    terminate("error");
                    break;
                  case FILE_INT:
                    IO_Fields[field].io_func(particle, IO_Fields[field].values_per_block, intp, 0);
                    intp += IO_Fields[field].values_per_block;
                    n++;
                    break;
                  case FILE_MY_ID_TYPE:
                    IO_Fields[field].io_func(particle, IO_Fields[field].values_per_block, ip, 0);
                    ip += IO_Fields[field].values_per_block;
                    n++;
                    break;
                  case FILE_MY_IO_FLOAT:
                    IO_Fields[field].io_func(particle, IO_Fields[field].values_per_block, fp, 0);
                    fp += IO_Fields[field].values_per_block;
                    n++;
                    break;
                  case FILE_DOUBLE:
                    IO_Fields[field].io_func(particle, IO_Fields[field].values_per_block, doublep, 0);
                    doublep += IO_Fields[field].values_per_block;
                    n++;
                    break;
                  case FILE_FLOAT:
                    IO_Fields[field].io_func(particle, IO_Fields[field].values_per_block, floatp, 0);
                    floatp += IO_Fields[field].values_per_block;
                    n++;
                    break;
                }
            }
          else
            {
              void *array_pos;

              switch(IO_Fields[field].array)
                {
                  case A_NONE:
                    array_pos = 0;
                    break;

                  case A_SPHP:
                    array_pos = SphP + pindex;
                    break;

                  case A_P:
                    array_pos = P + pindex;
                    break;
#if defined(GFM) || defined(SFR_MCS)
                  case A_STARP:
                    array_pos = StarP + P[pindex].AuxDataID;
                    break;
#endif
#ifdef BLACK_HOLES
                  case A_BHP:
                    array_pos = BHP + P[pindex].AuxDataID;
                    break;
#endif
#ifdef DUST_LIVE
                  case A_DUSTP:
                    array_pos = DustP + P[pindex].AuxDataID;
                    break;
#endif
                  case A_PS:
                    array_pos = PS + pindex;
                    break;

                  default:
                    terminate("ERROR in fill_write_buffer: Array not found!");
                    break;
                }

              for(int k = 0; k < IO_Fields[field].values_per_block; k++)
                {
                  double value = 0.;

                  switch(IO_Fields[field].type_in_memory)
                    {
                      case MEM_INT:
                        *intp = *((int *)((size_t)array_pos + IO_Fields[field].offset + k * sizeof(int)));
                        intp++;
                        break;

                      case MEM_MY_ID_TYPE:
                        *ip = *((MyIDType *)((size_t)array_pos + IO_Fields[field].offset + k * sizeof(MyIDType)));
                        ip++;
                        break;

                      case MEM_FLOAT:
                        value = *((float *)((size_t)array_pos + IO_Fields[field].offset + k * sizeof(float)));
                        break;

                      case MEM_DOUBLE:
                        value = *((double *)((size_t)array_pos + IO_Fields[field].offset + k * sizeof(double)));
                        break;

                      case MEM_MY_SINGLE:
                        value = *((MySingle *)((size_t)array_pos + IO_Fields[field].offset + k * sizeof(MySingle)));
                        break;

                      case MEM_MY_FLOAT:
                        value = *((MyFloat *)((size_t)array_pos + IO_Fields[field].offset + k * sizeof(MyFloat)));
                        break;

                      case MEM_MY_DOUBLE:
                        value = *((MyDouble *)((size_t)array_pos + IO_Fields[field].offset + k * sizeof(MyDouble)));
                        break;

                      case MEM_NONE:
                        terminate("ERROR in fill_write_buffer: reached MEM_NONE with no io_func specified!");
                        break;

                      default:
                        terminate("ERROR in fill_write_buffer: Type not found!");
                        break;
                    }

                  switch(IO_Fields[field].type_in_file_output)
                    {
                      case FILE_MY_IO_FLOAT:
                        *fp = value;
                        fp++;
                        break;

                      case FILE_DOUBLE:
                        *doublep = value;
                        doublep++;
                        break;

                      case FILE_FLOAT:
                        *floatp = value;
                        floatp++;
                        break;

                      default:
                        break;
                    }
                }

              n++;
            } /* end io_func/not */
        }     /* end type if */
    }         /* end particle loop */

  *startindex = pindex;
}

/*! \brief This function tells the size in bytes of one data entry in each of
 *         the blocks defined for the output file.
 *
 *  \param[in] blocknr ID of the output block (i.e. position, velocities...).
 *  \param[in] mode Used to distinguish whether the function is called in input
 *             mode (mode > 0) or in output mode (mode = 0). The size of one
 *             data entry may vary depending on the mode.
 *
 *  \return Size of the data entry in bytes.
 */
int get_bytes_per_blockelement(const enum iofields blocknr, const int mode)
{
  int bytes_per_blockelement = 0;

  for(int f = 0; f < N_IO_Fields; f++)
    {
      if(IO_Fields[f].field == blocknr)
        {
          if(mode)
            {
              switch(IO_Fields[f].type_in_file_input)
                {
                  case FILE_NONE:
                    terminate("error");
                    break;
                  case FILE_INT:
                    bytes_per_blockelement = IO_Fields[f].values_per_block * sizeof(int);
                    break;
                  case FILE_MY_ID_TYPE:
                    bytes_per_blockelement = IO_Fields[f].values_per_block * sizeof(MyIDType);
                    break;
                  case FILE_MY_IO_FLOAT:
                    bytes_per_blockelement = IO_Fields[f].values_per_block * sizeof(MyInputFloat);
                    break;
                  case FILE_DOUBLE:
                    bytes_per_blockelement = IO_Fields[f].values_per_block * sizeof(double);
                    break;
                  case FILE_FLOAT:
                    bytes_per_blockelement = IO_Fields[f].values_per_block * sizeof(float);
                    break;
                }
            }
          else
            {
              switch(IO_Fields[f].type_in_file_output)
                {
                  case FILE_NONE:
                    terminate("error");
                    break;
                  case FILE_INT:
                    bytes_per_blockelement = IO_Fields[f].values_per_block * sizeof(int);
                    break;
                  case FILE_MY_ID_TYPE:
                    bytes_per_blockelement = IO_Fields[f].values_per_block * sizeof(MyIDType);
                    break;
                  case FILE_MY_IO_FLOAT:
                    bytes_per_blockelement = IO_Fields[f].values_per_block * sizeof(MyOutputFloat);
                    break;
                  case FILE_DOUBLE:
                    bytes_per_blockelement = IO_Fields[f].values_per_block * sizeof(double);
                    break;
                  case FILE_FLOAT:
                    bytes_per_blockelement = IO_Fields[f].values_per_block * sizeof(float);
                    break;
                }
            }
          break;
        }
    }

  return bytes_per_blockelement;
}

/*! \brief This function determines the type of one data entry in each of the
 *         blocks defined for the output file.
 *
 *  Used only if output in HDF5 format is enabled.
 *
 *  \param[in] blocknr ID of the output block (i.e. position, velocities...).
 *  \param[in] mode For input mode > 0, for output mode = 0.
 *
 *  \return typekey, a flag that indicates the type of the data entry.
 */
int get_datatype_in_block(const enum iofields blocknr, const int mode)
{
  for(int f = 0; f < N_IO_Fields; f++)
    if(IO_Fields[f].field == blocknr)
      {
        int typekey;
        if(mode)
          typekey = IO_Fields[f].type_in_file_input;
        else
          typekey = IO_Fields[f].type_in_file_output;
        return typekey;
      }
  terminate("error invalid field");
}

/*! \brief This function determines the number of elements composing one data
 *         entry in each of the blocks defined for the output file.
 *
 *  Used only if output in HDF5 format is enabled.
 *
 *  \param[in] blocknr ID of the output block (i.e. position, velocities...).
 *
 *  \return Number of elements of one data entry.
 */
int get_values_per_blockelement(const enum iofields blocknr)
{
  for(int f = 0; f < N_IO_Fields; f++)
    if(IO_Fields[f].field == blocknr)
      {
        const int values = IO_Fields[f].values_per_block;
        return values;
      }
  terminate("reached last entry in switch - strange.");
}

/*! \brief Gets particle number in an output block.
 *
 *  This function determines how many particles there are in a given block,
 *  based on the information in the header structure. It also flags particle
 *  types that are present in the block in the typelist array.
 *
 *  \param[in] blocknr ID of the output block (i.e. position, velocities...).
 *  \param[out] typelist Array of length NTYPES which indicates whether
 *                       particles of the respective types are present in the
 *                       given block.
 *
 *  \return The total number of particles in the block.
 */
int get_particles_in_block(const enum iofields blocknr, int *const typelist)
{
  int npart = 0;
  switch(blocknr)
    {
      case IO_MASS:
        for(int i = 0; i < NTYPES; i++)
          {
            typelist[i] = 0;
#ifndef SIDM
            /* for SIDM always write out individual particle masses since those can change for the inelastic case */
            if(All.MassTable[i] == 0)
#endif
              if(header.npart[i] > 0)
                {
#ifdef TRACER_MC
                  if(i == TRACER_MC)
                    /* temporary, only needed for old snapshots with
                     * MassTable[TRACER_MC] = 0 */
                    continue;
#endif
                  typelist[i] = 1;
                  npart += header.npart[i];
                }
          }
        return npart; /* with masses */
        break;

      /* TODO: this IC-type-dependent logic has to go */
      case IO_AGE:
      case IO_LOCFBEVENT:
        for(int i = 0; i < NTYPES; i++)
          if((All.ICFormat == SNAP_FORMAT_HDF5 && (i == PTYPE_GAS || i == PTYPE_HALO || i == PTYPE_BNDRY)) ||
             (All.ICFormat != SNAP_FORMAT_HDF5 && i != PTYPE_STARS))
            typelist[i] = 0;
        return header.npart[PTYPE_STARS] +
               (All.ICFormat == SNAP_FORMAT_HDF5 ? (header.npart[PTYPE_DISK] + header.npart[PTYPE_BULGE]) : 0);
        break;

      case IO_Z:
        for(int i = 0; i < NTYPES; i++)
          {
            if((All.ICFormat == SNAP_FORMAT_HDF5 && (i == PTYPE_HALO || i == PTYPE_BNDRY)) ||
               (All.ICFormat != SNAP_FORMAT_HDF5 && (i != PTYPE_GAS && i != PTYPE_STARS)))
              typelist[i] = 0;
            else
              typelist[i] = 1;
          }
        return header.npart[PTYPE_GAS] + header.npart[PTYPE_STARS] +
               (All.ICFormat == SNAP_FORMAT_HDF5 ? (header.npart[PTYPE_DISK] + header.npart[PTYPE_BULGE]) : 0);
        break;

      case IO_LASTENTRY:
        terminate("reached last entry in switch - strange.");
        break;

      default:
        for(int f = 0; f < N_IO_Fields; f++)
          {
            if(IO_Fields[f].field == blocknr)
              {
                for(int i = 0; i < NTYPES; i++)
                  {
                    if((IO_Fields[f].typelist & (1 << i)) && header.npart[i] > 0)
                      {
                        typelist[i] = 1;
                        npart += header.npart[i];
                      }
                    else
                      typelist[i] = 0;
                  }

                return npart;
              }
          }
        break;
    } /* end switch */

  terminate("reached end of function - this should not happen");
}

/*! \brief Checks if a block is expected for file input or output.
 *
 *  This function tells whether a block in the input/output file is requested
 *  or not. Because the blocks processed in the two cases are different, the
 *  mode is indicated with the flag write (1=write, 0=read).
 *
 *  \param[in] blocknr ID of the output block (i.e. position, velocities...).
 *  \param[in] write If 0 the function is in read mode, if 1 the function is
 *             in write mode.
 *
 *  \return 0 if the block is not present, 1 otherwise.
 */
int blockpresent(const enum iofields blocknr, const int write)
{
  /* TODO: get rid of these special conditions */
#ifndef READ_IN_ALL_IC_FIELDS
  if(!write)
    {
#ifdef SGS_TURBULENCE_IN_ICS
/* the modules have to be in the same order as the corresponding fields in the ICs */
#ifdef MHD
      if(RestartFlag == RESTART_IC && blocknr == IO_BFLD)
        return 1;
#endif
#ifdef COSMIC_RAYS_IN_ICS
      if(RestartFlag == RESTART_IC && blocknr == IO_CRENERGY)
        return 1;
#endif
      if(RestartFlag == RESTART_IC && blocknr == IO_SGS_T_SPECIFICENERGY)
        return 1;
#endif

#ifdef PASSIVE_SCALARS
      if(RestartFlag == RESTART_IC && blocknr == IO_PASS)
        return 1;
#endif

#ifdef TGSET
      if((RestartFlag == RESTART_IC || RestartFlag == RESTART_FOF_SUBFIND) && blocknr > IO_MASS)
#else
#if defined(EOS_DEGENERATE) || defined(EOS_OPAL)
#ifdef MHD
      if(RestartFlag == RESTART_IC && (blocknr > IO_U && blocknr != IO_EOSXNUC && blocknr != IO_BFLD))
#else
      if(RestartFlag == RESTART_IC && (blocknr > IO_U && blocknr != IO_EOSXNUC))
#endif
#else
#ifdef CR_IC
      if(RestartFlag == RESTART_IC && ((blocknr > IO_CR_Q0 && blocknr != IO_BFLD) || (blocknr >= IO_RHO && blocknr <= IO_ACCEL)))
#else
      /* STANDARD CODE PATH HERE. We only shortcut the file read
       * if we aren't using HDF5. With HDF5, we simply try to read
       * all the blocks that could be present. */
#if defined(MHD) && !defined(MHD_SEEDFIELD)
      if(All.ICFormat != SNAP_FORMAT_HDF5 && RestartFlag == RESTART_IC && (blocknr > IO_U && blocknr != IO_BFLD))
#else
      if(All.ICFormat != SNAP_FORMAT_HDF5 && RestartFlag == RESTART_IC && blocknr > IO_U)
#endif
#endif /* CR_IC */
#endif /* EOS_DEGENERATE */
#endif /* TGSET */
        /* ignore all other blocks in non-HDF5 initial conditions */
        return 0;
    }
#endif /* #ifndef READ_IN_ALL_IC_FIELDS */

#ifdef HCOUTPUT
  if(DumpFlag == DUMP_HCOUTPUT)
    {
      if(blocknr == IO_POS || blocknr == IO_VEL || blocknr == IO_ID || blocknr == IO_MASS || blocknr == IO_GFM_AGE)
        return 1;
      else
        return 0;
    }
#endif
  for(int f = 0; f < N_IO_Fields; f++)
    {
      if(IO_Fields[f].field == blocknr)
        {
          if(!write)
            {
              if(IO_Fields[f].type_in_file_input != FILE_NONE)
                return 1;
            }
          else
            {
              if(IO_Fields[f].type_in_file_output == FILE_NONE)
                return 0;

              /* subboxes: write all fields except those marked by SN_NO_SUBBOX or SN_MINI_ONLY
                 (must come first to ignore DumpFlag) */
              if(subbox_dump)
                {
                  if(IO_Fields[f].snap_type == SN_NO_SUBBOX || IO_Fields[f].snap_type == SN_MINI_ONLY)
                    return 0;

                  return 1;
                }

              /* normal full snapshot (with or without groupcat): only skip fields marked by SN_MINI_ONLY */
              if(DumpFlag == DUMP_BOTH || DumpFlag == DUMP_ONLY_SNAP)
                {
                  if(IO_Fields[f].snap_type == SN_MINI_ONLY)
                    return 0;

                  return 1;
                }

              /* mini-snaps: write only those fields marked by either SN_MINI or SN_MINI_ONLY */
              if(DumpFlag == DUMP_BOTH_MINI)
                {
                  if(IO_Fields[f].snap_type == SN_MINI || IO_Fields[f].snap_type == SN_MINI_ONLY)
                    return 1;

                  if(IO_Fields[f].typelist == BHS_ONLY)
                    /* temporarily hard-coded that all BH fields are included
                     * in mini-snaps */
                    return 1;

                  /* specifically do not include any other fields in mini-snaps */
                  return 0;
                }
            }
          return 0;
        }
    }

  /* default: not present */
  return 0;
}

/*! \brief This function associates a short 4-character block name with each
 *         block number.
 *
 *   This is stored in front of each block for snapshot FileFormat=2.
 *
 *  \param[in] blocknr ID of the output block (i.e. position, velocities...).
 *  \param[in] label string containing the dataset name.
 *
 *  \return void
 */
void get_Tab_IO_Label(const enum iofields blocknr, char *const label)
{
  for(int f = 0; f < N_IO_Fields; f++)
    {
      if(IO_Fields[f].field == blocknr)
        {
          memcpy(label, IO_Fields[f].label, IO_LABEL_SIZE);
          return;
        }
    }
  terminate("error invalid field");
}

/*! \brief This function associates a dataset name with each block number.
 *
 *   This is needed to name the dataset if the output is written in HDF5
 *   format.
 *
 *  \param[in] blocknr ID of the output block (i.e. position, velocities...).
 *  \param[in] buf String containing the dataset name.
 *
 *  \return void
 */
void get_dataset_name(const enum iofields blocknr, char *const buf)
{
  for(int f = 0; f < N_IO_Fields; f++)
    {
      if(IO_Fields[f].field == blocknr)
        {
          snprintf(buf, IO_DATASET_NAME_SIZE, "%s", IO_Fields[f].datasetname);
          return;
        }
    }
  terminate("error invalid field");
}

/*! \brief Actually write the snapshot file to the disk.
 *
 *  This function writes a snapshot file containing the data from processors
 *  'writeTask' to 'lastTask'. 'writeTask' is the one that actually writes.
 *  Each snapshot file contains a header and cell/particle details. The
 *  output fields for each particle type depend on included physics
 *  and compile-time flags.
 *
 *  \param[in] fname String containing the file name.
 *  \param[in] writeTask The rank of the task in a writing group that which
 *             is responsible for the output operations.
 *  \param[in] lastTask The rank of the last task in a writing group.
 *  \param[in] subbox_flag If greater than 0: instructs the code to output
 *             only a subset of the whole domain.
 *
 *  \return void
 */
void write_file(const char *fname, const int writeTask, const int lastTask, const int subbox_flag)
{
  int nextblock;
  int n_for_this_task, pc;
  int ntot_type[NTYPES], nn[NTYPES];
  char label[IO_LABEL_SIZE];
  int blksize;
  MPI_Status status;
  FILE *fd  = NULL;
  int pcsum = 0;

#ifdef HAVE_HDF5
  hid_t hdf5_file = 0, hdf5_grp[NTYPES], hdf5_headergrp = 0, hdf5_dataspace_memory;
  hid_t hdf5_datatype = 0, hdf5_dataspace_in_file = 0, hdf5_dataset = 0;
  hsize_t dims[2], count[2], start[2];
  int rank = 0;
  char buf[MAXLEN_PATH];
  char dataset_name[IO_DATASET_NAME_SIZE];
#ifdef HDF5_FILTERS
  hid_t hdf5_properties = 0;
#endif
  hid_t hdf5_paramsgrp = 0;
  hid_t hdf5_configgrp = 0;
#endif

#define SKIP                                 \
  {                                          \
    my_fwrite(&blksize, sizeof(int), 1, fd); \
  }

  myassert(writeTask <= lastTask);

#ifdef TOLERATE_WRITE_ERROR
  for(int try_io = 0; try_io < 2; try_io++)
    {
      WriteErrorFlag = 0;
#endif

      /* determine particle numbers of each type in file */

      if(ThisTask == writeTask)
        {
          for(int n = 0; n < NTYPES; n++)
            ntot_type[n] = n_type[n];

          for(int task = writeTask + 1; task <= lastTask; task++)
            {
              MPI_Recv(&nn[0], NTYPES, MPI_INT, task, TAG_LOCALN, MPI_COMM_WORLD, &status);
              for(int n = 0; n < NTYPES; n++)
                ntot_type[n] += nn[n];
            }

          for(int task = writeTask + 1; task <= lastTask; task++)
            MPI_Send(&ntot_type[0], NTYPES, MPI_INT, task, TAG_N, MPI_COMM_WORLD);
        }
      else
        {
          MPI_Send(&n_type[0], NTYPES, MPI_INT, writeTask, TAG_LOCALN, MPI_COMM_WORLD);
          MPI_Recv(&ntot_type[0], NTYPES, MPI_INT, writeTask, TAG_N, MPI_COMM_WORLD, &status);
        }

#ifdef HCOUTPUT
      if(subbox_flag > 0)
        {
          for(int n = 0; n < NTYPES; n++)
            {
              if(!((1 << P[n].Type) & (HCOUTPUT)))
                {
                  n_type[n]    = 0;
                  ntot_type[n] = 0;
                }
            }
        }
#endif

      /* fill file header */

      for(int n = 0; n < NTYPES; n++)
        {
          header.npart[n]              = ntot_type[n];
          header.npartTotal[n]         = (unsigned int)ntot_type_all[n];
          header.npartTotalHighWord[n] = (unsigned int)(ntot_type_all[n] >> 32);
          header.mass[n]               = All.MassTable[n];
        }

      header.time = All.Time;

      if(All.ComovingIntegrationOn)
        header.redshift = 1.0 / All.Time - 1;
      else
        header.redshift = 0;

      header.flag_sfr        = 0;
      header.flag_feedback   = 0;
      header.flag_cooling    = 0;
      header.flag_stellarage = 0;
      header.flag_metals     = 0;

      header.flag_tracer_field = 0;
#ifdef TRACER_FIELD
      header.flag_tracer_field = 1;
#endif

#ifdef COOLING
      header.flag_cooling = 1;
#endif

#ifdef USE_SFR
      header.flag_sfr      = 1;
      header.flag_feedback = 1;
#ifdef STELLARAGE
      header.flag_stellarage = 1;
#endif
#ifdef METALS
      header.flag_metals = 1;
#endif
#endif

      header.num_files   = All.NumFilesPerSnapshot;
      header.BoxSize     = All.BoxSize;
      header.Omega0      = All.Omega0;
      header.OmegaLambda = All.OmegaLambda;
      header.HubbleParam = All.HubbleParam;

#ifdef OUTPUT_IN_DOUBLEPRECISION
      header.flag_doubleprecision = 1;
#else
  header.flag_doubleprecision = 0;
#endif

      /* open file and write header */
      if(ThisTask == writeTask)
        {
          if(All.SnapFormat == SNAP_FORMAT_HDF5)
            {
#ifdef HAVE_HDF5
              file_path_sprintf(buf, "%s.hdf5", fname);
              hdf5_file = my_H5Fcreate(buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
              mpi_printf("filename=%s\n", buf);
              hdf5_headergrp = my_H5Gcreate(hdf5_file, "/Header", 0);

              for(int type = 0; type < NTYPES; type++)
                if(header.npart[type] > 0)
                  {
                    sprintf(buf, "/PartType%d", type);
                    hdf5_grp[type] = my_H5Gcreate(hdf5_file, buf, 0);
                  }

              write_header_attributes_in_hdf5(hdf5_headergrp);

              hdf5_paramsgrp = my_H5Gcreate(hdf5_file, "/Parameters", 0);
              write_parameters_attributes_in_hdf5(hdf5_paramsgrp);

              hdf5_configgrp = my_H5Gcreate(hdf5_file, "/Config", 0);
              write_compile_time_options_in_hdf5(hdf5_configgrp);
#endif
            }
          else
            {
              if(!(fd = fopen(fname, "w")))
                terminate("file open error: can't open file `%s' for writing snapshot.", fname);

              if(All.SnapFormat == SNAP_FORMAT_GADGET_VARIANT)
                {
                  blksize = sizeof(int) + IO_LABEL_SIZE * sizeof(char);
                  SKIP;
                  my_fwrite((void *)"HEAD", sizeof(char), IO_LABEL_SIZE, fd);
                  nextblock = sizeof(header) + 2 * sizeof(int);
                  my_fwrite(&nextblock, sizeof(int), 1, fd);
                  SKIP;
                }

              blksize = sizeof(header);
              SKIP;
              my_fwrite(&header, sizeof(header), 1, fd);
              SKIP;
            }
        }

      /* write blocks with particle data */
      for(int bnr = 0; bnr < IO_LASTENTRY; bnr++)
        {
          const enum iofields blocknr = (enum iofields)bnr;

          if(blockpresent(blocknr, 1))
            {
              const int bytes_per_blockelement = get_bytes_per_blockelement(blocknr, 0);
              const int blockmaxlen            = (int)(COMMBUFFERSIZE / bytes_per_blockelement);

              int typelist[NTYPES];
              const int npart = get_particles_in_block(blocknr, typelist);
              if(npart > 0)
                {
                  if(ThisTask == 0)
                    {
                      char block_name[IO_DATASET_NAME_SIZE];
                      get_dataset_name(blocknr, block_name);
                      if(subbox_flag == 0)
                        mpi_printf("writing block %d (%s)...\n", blocknr, block_name);
                    }

                  if(ThisTask == writeTask && (All.SnapFormat == SNAP_FORMAT_GADGET || All.SnapFormat == SNAP_FORMAT_GADGET_VARIANT))
                    {
                      if(All.SnapFormat == SNAP_FORMAT_GADGET_VARIANT)
                        {
                          blksize = sizeof(int) + IO_LABEL_SIZE * sizeof(char);
                          SKIP;
                          get_Tab_IO_Label(blocknr, label);
                          my_fwrite(label, sizeof(char), IO_LABEL_SIZE, fd);
                          nextblock = npart * bytes_per_blockelement + 2 * sizeof(int);
                          my_fwrite(&nextblock, sizeof(int), 1, fd);
                          SKIP;
                        }

                      blksize = npart * bytes_per_blockelement;
                      SKIP;
                    }

                  for(int type = 0; type < NTYPES; type++)
                    {
                      if(typelist[type])
                        {
#ifdef HAVE_HDF5
                          if(ThisTask == writeTask && All.SnapFormat == SNAP_FORMAT_HDF5 && header.npart[type] > 0)
                            {
                              switch(get_datatype_in_block(blocknr, 0))
                                {
                                  case FILE_INT:
#ifndef SPECIAL_BOUNDARY
                                    hdf5_datatype = my_H5Tcopy(H5T_NATIVE_UINT);
#else
                                    hdf5_datatype = my_H5Tcopy(H5T_NATIVE_INT);
#endif
                                    break;
                                  case FILE_MY_IO_FLOAT:
#ifdef OUTPUT_IN_DOUBLEPRECISION
                                    hdf5_datatype = my_H5Tcopy(H5T_NATIVE_DOUBLE);
#else
                                    hdf5_datatype = my_H5Tcopy(H5T_NATIVE_FLOAT);
#endif
                                    break;
                                  case FILE_MY_ID_TYPE:
#ifdef LONGIDS
#ifndef SPECIAL_BOUNDARY
                                    hdf5_datatype = my_H5Tcopy(H5T_NATIVE_UINT64);
#else
                                    hdf5_datatype = my_H5Tcopy(H5T_NATIVE_INT64);
#endif
#else
#ifndef SPECIAL_BOUNDARY
                                    hdf5_datatype = my_H5Tcopy(H5T_NATIVE_UINT32);
#else
                                    hdf5_datatype = my_H5Tcopy(H5T_NATIVE_INT32);
#endif
#endif
                                    break;
                                  case FILE_DOUBLE:
                                    hdf5_datatype = my_H5Tcopy(H5T_NATIVE_DOUBLE);
                                    break;
                                  case FILE_FLOAT:
                                    hdf5_datatype = my_H5Tcopy(H5T_NATIVE_FLOAT);
                                    break;
                                }

                              dims[0] = header.npart[type];
                              dims[1] = get_values_per_blockelement(blocknr);
                              if(dims[1] == 1)
                                rank = 1;
                              else
                                rank = 2;

                              get_dataset_name(blocknr, dataset_name);

                              hdf5_dataspace_in_file = my_H5Screate_simple(rank, dims, NULL);
#ifdef HDF5_FILTERS
                              hdf5_properties = my_H5Pcreate(H5P_DATASET_CREATE);
                              /* set chunk size */
                              my_H5Pset_chunk(hdf5_properties, rank, dims);
                              /* reshuffle bytes to get better compression ratio */
                              my_H5Pset_shuffle(hdf5_properties);
                              /* gzip compression level 9 */
                              my_H5Pset_deflate(hdf5_properties, 9);
                              /* Fletcher32 checksum on dataset */
                              my_H5Pset_fletcher32(hdf5_properties);

                              if(my_H5Pall_filters_avail(hdf5_properties))
                                hdf5_dataset =
                                    my_H5Dcreate(hdf5_grp[type], dataset_name, hdf5_datatype, hdf5_dataspace_in_file, hdf5_properties);
                              else
                                {
                                  printf("HDF5_FILTERS: Warning selected filters not available! Writing data without filters!\n");
                                  myflush(stdout);
                                  hdf5_dataset =
                                      my_H5Dcreate(hdf5_grp[type], dataset_name, hdf5_datatype, hdf5_dataspace_in_file, H5P_DEFAULT);
                                }
#else
                              hdf5_dataset =
                                  my_H5Dcreate(hdf5_grp[type], dataset_name, hdf5_datatype, hdf5_dataspace_in_file, H5P_DEFAULT);
#endif
                              write_dataset_attributes(hdf5_dataset, blocknr);
                            }
#endif

                          pcsum               = 0;
                          int remaining_space = blockmaxlen;
                          int bufferstart     = 0;

                          int offset = 0;
                          for(int task = writeTask; task <= lastTask; task++)
                            {
                              if(task == ThisTask)
                                {
                                  n_for_this_task = n_type[type];
                                  for(int p = writeTask; p <= lastTask; p++)
                                    {
                                      if(p != ThisTask)
                                        MPI_Send(&n_for_this_task, 1, MPI_INT, p, TAG_NFORTHISTASK, MPI_COMM_WORLD);
                                    }
                                }
                              else
                                {
                                  MPI_Recv(&n_for_this_task, 1, MPI_INT, task, TAG_NFORTHISTASK, MPI_COMM_WORLD, &status);
                                }

                              while(n_for_this_task > 0)
                                {
                                  pc = n_for_this_task;

                                  if(pc > remaining_space)
                                    pc = remaining_space;

                                  void *buffer = (void *)((char *)CommBuffer + bufferstart * bytes_per_blockelement);

                                  if(ThisTask == task)
                                    fill_write_buffer(buffer, blocknr, &offset, pc, type, subbox_flag);

                                  if(ThisTask == writeTask && task != writeTask)
                                    MPI_Recv(buffer, bytes_per_blockelement * pc, MPI_BYTE, task, TAG_PDATA, MPI_COMM_WORLD, &status);

                                  if(ThisTask != writeTask && task == ThisTask)
                                    MPI_Ssend(buffer, bytes_per_blockelement * pc, MPI_BYTE, writeTask, TAG_PDATA, MPI_COMM_WORLD);

                                  remaining_space -= pc;
                                  bufferstart += pc;

                                  if(remaining_space == 0)
                                    {
                                      /* write stuff (number of elements equal to bufferstart) */
                                      if(ThisTask == writeTask)
                                        {
                                          if(All.SnapFormat == SNAP_FORMAT_HDF5)
                                            {
#ifdef HAVE_HDF5
                                              start[0] = pcsum;
                                              start[1] = 0;

                                              count[0] = bufferstart;
                                              count[1] = get_values_per_blockelement(blocknr);

                                              my_H5Sselect_hyperslab(hdf5_dataspace_in_file, H5S_SELECT_SET, start, NULL, count, NULL);

                                              dims[0]               = bufferstart;
                                              dims[1]               = get_values_per_blockelement(blocknr);
                                              hdf5_dataspace_memory = my_H5Screate_simple(rank, dims, NULL);

                                              my_H5Dwrite(hdf5_dataset, hdf5_datatype, hdf5_dataspace_memory, hdf5_dataspace_in_file,
                                                          H5P_DEFAULT, CommBuffer, dataset_name);

                                              my_H5Sclose(hdf5_dataspace_memory, H5S_SIMPLE);
#endif
                                            }
                                          else
                                            {
                                              my_fwrite(CommBuffer, bytes_per_blockelement, bufferstart, fd);
                                            }
                                        }

                                      pcsum += bufferstart;
                                      remaining_space = blockmaxlen;
                                      bufferstart     = 0;
                                    }

                                  n_for_this_task -= pc;
                                }
                            }

                          if(bufferstart > 0)
                            {
                              /* write remaining stuff (number of elements equal to bufferstart) */
                              if(ThisTask == writeTask)
                                {
                                  if(All.SnapFormat == SNAP_FORMAT_HDF5)
                                    {
#ifdef HAVE_HDF5
                                      start[0] = pcsum;
                                      start[1] = 0;

                                      count[0] = bufferstart;
                                      count[1] = get_values_per_blockelement(blocknr);

                                      my_H5Sselect_hyperslab(hdf5_dataspace_in_file, H5S_SELECT_SET, start, NULL, count, NULL);

                                      dims[0]               = bufferstart;
                                      dims[1]               = get_values_per_blockelement(blocknr);
                                      hdf5_dataspace_memory = my_H5Screate_simple(rank, dims, NULL);

                                      my_H5Dwrite(hdf5_dataset, hdf5_datatype, hdf5_dataspace_memory, hdf5_dataspace_in_file,
                                                  H5P_DEFAULT, CommBuffer, dataset_name);

                                      my_H5Sclose(hdf5_dataspace_memory, H5S_SIMPLE);
#endif
                                    }
                                  else
                                    {
                                      my_fwrite(CommBuffer, bytes_per_blockelement, bufferstart, fd);
                                    }
                                }

                              pcsum += bufferstart;
                              remaining_space = blockmaxlen;
                              bufferstart     = 0;
                            }

#ifdef HAVE_HDF5
                          if(ThisTask == writeTask && All.SnapFormat == SNAP_FORMAT_HDF5 && header.npart[type] > 0)
                            {
                              my_H5Dclose(hdf5_dataset, dataset_name);
#ifdef HDF5_FILTERS
                              my_H5Pclose(hdf5_properties);
#endif
                              my_H5Sclose(hdf5_dataspace_in_file, H5S_SIMPLE);
                              my_H5Tclose(hdf5_datatype);
                            }
#endif
                        }
                    }

                  if(ThisTask == writeTask && (All.SnapFormat == SNAP_FORMAT_GADGET || All.SnapFormat == SNAP_FORMAT_GADGET_VARIANT))
                    SKIP;
                }

#ifdef TOLERATE_WRITE_ERROR
              if(ThisTask == writeTask)
                {
                  for(int p = writeTask; p <= lastTask; p++)
                    if(p != ThisTask)
                      MPI_Send(&WriteErrorFlag, 1, MPI_INT, p, TAG_KEY, MPI_COMM_WORLD);
                }
              else
                MPI_Recv(&WriteErrorFlag, 1, MPI_INT, writeTask, TAG_KEY, MPI_COMM_WORLD, &status);
#endif
            }

#ifdef TOLERATE_WRITE_ERROR
          if(WriteErrorFlag) /* don't write further blocks in this case */
            break;
#endif
        }

      if(ThisTask == writeTask)
        {
          if(All.SnapFormat == SNAP_FORMAT_HDF5)
            {
#ifdef HAVE_HDF5
              for(int type = NTYPES - 1; type >= 0; type--)
                if(header.npart[type] > 0)
                  my_H5Gclose(hdf5_grp[type], buf);
              my_H5Gclose(hdf5_headergrp, "/Header");
              my_H5Gclose(hdf5_paramsgrp, "/Parameters");
              my_H5Gclose(hdf5_configgrp, "/Config");

              file_path_sprintf(buf, "%s.hdf5", fname);
              my_H5Fclose(hdf5_file, buf);
#endif
            }
          else
            fclose(fd);
        }

#ifdef TOLERATE_WRITE_ERROR
      if(WriteErrorFlag == 0)
        break;

      if(try_io == 0)
        {
          if(ThisTask == writeTask)
            {
              printf(
                  "TOLERATE_WRITE_ERROR: Try to write to alternative file: masterTask=%d  lastTask=%d  try_io=%d "
                  "alternative-filename='%s'\n",
                  writeTask, lastTask, try_io, alternative_fname);
              myflush(stdout);
            }
          /* try on a different output directory */
          fname = alternative_fname;
        }
      else
        {
          terminate("TOLERATE_WRITE_ERROR: Second try with alternative file failed too.");
        }
    }
#endif
}

#ifdef HAVE_HDF5
/*! \brief Write the fields contained in the header group of the HDF5 snapshot file
 *
 *  This function stores the fields of the structure io_header as attributes belonging
 *  to the header group of the HDF5 file.
 *
 *  \param handle contains a reference to the header group
 */
void write_header_attributes_in_hdf5(hid_t handle)
{
  hsize_t adim[1] = {NTYPES};
  hid_t hdf5_dataspace, hdf5_attribute;

  hdf5_dataspace = my_H5Screate(H5S_SIMPLE);
  my_H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL, "NumPart_ThisFile");
  hdf5_attribute = my_H5Acreate(handle, "NumPart_ThisFile", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_INT, header.npart, "NumPart_ThisFile");
  my_H5Aclose(hdf5_attribute, "NumPart_ThisFile");
  my_H5Sclose(hdf5_dataspace, H5S_SIMPLE);

  hdf5_dataspace = my_H5Screate(H5S_SIMPLE);
  my_H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL, "NumPart_Total");
  hdf5_attribute = my_H5Acreate(handle, "NumPart_Total", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, header.npartTotal, "NumPart_Total");
  my_H5Aclose(hdf5_attribute, "NumPart_Total");
  my_H5Sclose(hdf5_dataspace, H5S_SIMPLE);

  hdf5_dataspace = my_H5Screate(H5S_SIMPLE);
  my_H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL, "NumPart_Total_HighWord");
  hdf5_attribute = my_H5Acreate(handle, "NumPart_Total_HighWord", H5T_NATIVE_UINT, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_UINT, header.npartTotalHighWord, "NumPart_Total_HighWord");
  my_H5Aclose(hdf5_attribute, "NumPart_Total_HighWord");
  my_H5Sclose(hdf5_dataspace, H5S_SIMPLE);

  hdf5_dataspace = my_H5Screate(H5S_SIMPLE);
  my_H5Sset_extent_simple(hdf5_dataspace, 1, adim, NULL, "MassTable");
  hdf5_attribute = my_H5Acreate(handle, "MassTable", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, header.mass, "MassTable");
  my_H5Aclose(hdf5_attribute, "MassTable");
  my_H5Sclose(hdf5_dataspace, H5S_SIMPLE);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Time", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.time, "Time");
  my_H5Aclose(hdf5_attribute, "Time");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Redshift", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.redshift, "Redshift");
  my_H5Aclose(hdf5_attribute, "Redshift");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "BoxSize", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.BoxSize, "BoxSize");
  my_H5Aclose(hdf5_attribute, "BoxSize");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "NumFilesPerSnapshot", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.num_files, "NumFilesPerSnapshot");
  my_H5Aclose(hdf5_attribute, "NumFilesPerSnapshot");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Omega0", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.Omega0, "Omega0");
  my_H5Aclose(hdf5_attribute, "Omega0");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "OmegaLambda", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.OmegaLambda, "OmegaLambda");
  my_H5Aclose(hdf5_attribute, "OmegaLambda");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "OmegaBaryon", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.OmegaBaryon, "OmegaBaryon");
  my_H5Aclose(hdf5_attribute, "OmegaBaryon");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "HubbleParam", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &header.HubbleParam, "HubbleParam");
  my_H5Aclose(hdf5_attribute, "HubbleParam");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Flag_Sfr", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_sfr, "Flag_Sfr");
  my_H5Aclose(hdf5_attribute, "Flag_Sfr");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Flag_Cooling", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_cooling, "Flag_Cooling");
  my_H5Aclose(hdf5_attribute, "Flag_Cooling");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Flag_StellarAge", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_stellarage, "Flag_StellarAge");
  my_H5Aclose(hdf5_attribute, "Flag_StellarAge");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Flag_Metals", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_metals, "Flag_Metals");
  my_H5Aclose(hdf5_attribute, "Flag_Metals");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Flag_Feedback", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_feedback, "Flag_Feedback");
  my_H5Aclose(hdf5_attribute, "Flag_Feedback");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Flag_DoublePrecision", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.flag_doubleprecision, "Flag_DoublePrecision");
  my_H5Aclose(hdf5_attribute, "Flag_DoublePrecision");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Composition_vector_length", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &header.composition_vector_length, "Composition_vector_length");
  my_H5Aclose(hdf5_attribute, "Composition_vector_length");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hid_t atype = my_H5Tcopy(H5T_C_S1);

  my_H5Tset_size(atype, strlen(GIT_COMMIT));
  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Git_commit", atype, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, atype, GIT_COMMIT, "Git_commit");
  my_H5Aclose(hdf5_attribute, "Git_commit");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  my_H5Tset_size(atype, strlen(GIT_DATE));
  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Git_date", atype, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, atype, GIT_DATE, "Git_date");
  my_H5Aclose(hdf5_attribute, "Git_date");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

#if defined(DG) || defined(DG_TEST_PROBLEM)
  int out        = NUMDIMS;
  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Dims", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &out, "Dims");
  my_H5Aclose(hdf5_attribute, "Dims");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  out            = DEGREE_K;
  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Degree_K", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &out, "Degree_K");
  my_H5Aclose(hdf5_attribute, "Degree_K");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  double deltat  = timediff(StartOfRun, second());
  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Run_Time", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &deltat, "Run_Time");
  my_H5Aclose(hdf5_attribute, "Run_Time");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);
#endif

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "UnitLength_in_cm", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.UnitLength_in_cm, "UnitLength_in_cm");
  my_H5Aclose(hdf5_attribute, "UnitLength_in_cm");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "UnitMass_in_g", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.UnitMass_in_g, "UnitMass_in_g");
  my_H5Aclose(hdf5_attribute, "UnitMass_in_g");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "UnitVelocity_in_cm_per_s", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &All.UnitVelocity_in_cm_per_s, "UnitVelocity_in_cm_per_s");
  my_H5Aclose(hdf5_attribute, "UnitVelocity_in_cm_per_s");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);
}

/*! \brief Write the parameters read from the parameter file in the HDF5 snapshot file
 *
 *  This function stores the parameter io_header as attributes belonging
 *  to the parameter group of the HDF5 file.
 *
 *  \param handle contains a reference to the parameter group
 */
void write_parameters_attributes_in_hdf5(const hid_t handle)
{
  hid_t hdf5_dataspace, hdf5_attribute, atype = my_H5Tcopy(H5T_C_S1);

  my_H5Tset_size(atype, MAXLEN_PARAM_VALUE);

  for(int i = 0; i < All.NParameters; i++)
    {
      switch(ParametersType[i])
        {
          case PARAM_REAL:
            hdf5_dataspace = my_H5Screate(H5S_SCALAR);
            hdf5_attribute = my_H5Acreate(handle, Parameters[i], H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
            my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, ParametersValue[i], Parameters[i]);
            my_H5Aclose(hdf5_attribute, Parameters[i]);
            my_H5Sclose(hdf5_dataspace, H5S_SCALAR);
            break;
          case PARAM_STRING:
            hdf5_dataspace = my_H5Screate(H5S_SCALAR);
            hdf5_attribute = my_H5Acreate(handle, Parameters[i], atype, hdf5_dataspace, H5P_DEFAULT);
            my_H5Awrite(hdf5_attribute, atype, ParametersValue[i], Parameters[i]);
            my_H5Aclose(hdf5_attribute, Parameters[i]);
            my_H5Sclose(hdf5_dataspace, H5S_SCALAR);
            break;
          case PARAM_INT:
            hdf5_dataspace = my_H5Screate(H5S_SCALAR);
            hdf5_attribute = my_H5Acreate(handle, Parameters[i], H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
            my_H5Awrite(hdf5_attribute, H5T_NATIVE_INT, ParametersValue[i], Parameters[i]);
            my_H5Aclose(hdf5_attribute, Parameters[i]);
            my_H5Sclose(hdf5_dataspace, H5S_SCALAR);
            break;
        }
    }

  my_H5Tclose(atype);
}

/*! \brief A simple error handler for HDF5.
 *
 *  This function prints the HDF5 error stack and, if write errors are tolerated,
 *  calls the write_error() function to print information about the error and
 *  returns a positive integer to allow the repetition of the write operation
 *  (see also the HDF5 documentation).
 *
 *  \param[in] unused The parameter is not used, but it is necessary for
 *             compatibility with the HDF5 library.
 *
 *  \return 1 if the write error is tolerated, otherwise the run is terminated.
 */
herr_t my_hdf5_error_handler(void *unused)
{
  herr_t result;
#ifdef TOLERATE_WRITE_ERROR
  if(FlagNyt == 0)
    write_error(2, 0, 0);
  result = 1;
#else
  result = 0;
#endif
  H5Eprint(stdout);
  return result;
}

void write_dataset_attributes(const hid_t hdf5_dataset, const enum iofields blocknr)
{
  int ind = -1;
  for(int f = 0; f < N_IO_Fields; f++)
    {
      if(IO_Fields[f].field == blocknr)
        {
          ind = f;
          break;
        }
    }

  if(ind < 0)
    return;

  if(IO_Fields[ind].hasunit == 0)
    return;

  if(All.ComovingIntegrationOn)
    {
      hid_t hdf5_dataspace = my_H5Screate(H5S_SCALAR);
      hid_t hdf5_attribute = my_H5Acreate(hdf5_dataset, "a_scaling", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
      my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &IO_Fields[ind].a, "a_scaling");
      my_H5Aclose(hdf5_attribute, "a_scaling");
      my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

      hdf5_dataspace = my_H5Screate(H5S_SCALAR);
      hdf5_attribute = my_H5Acreate(hdf5_dataset, "h_scaling", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
      my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &IO_Fields[ind].h, "h_scaling");
      my_H5Aclose(hdf5_attribute, "h_scaling");
      my_H5Sclose(hdf5_dataspace, H5S_SCALAR);
    }
  else
    {
      double zero          = 0;
      hid_t hdf5_dataspace = my_H5Screate(H5S_SCALAR);
      hid_t hdf5_attribute = my_H5Acreate(hdf5_dataset, "a_scaling", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
      my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &zero, "a_scaling");
      my_H5Aclose(hdf5_attribute, "a_scaling");
      my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

      hdf5_dataspace = my_H5Screate(H5S_SCALAR);
      hdf5_attribute = my_H5Acreate(hdf5_dataset, "h_scaling", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
      my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &zero, "h_scaling");
      my_H5Aclose(hdf5_attribute, "h_scaling");
      my_H5Sclose(hdf5_dataspace, H5S_SCALAR);
    }

  hid_t hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hid_t hdf5_attribute = my_H5Acreate(hdf5_dataset, "length_scaling", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &IO_Fields[ind].L, "length_scaling");
  my_H5Aclose(hdf5_attribute, "length_scaling");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(hdf5_dataset, "mass_scaling", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &IO_Fields[ind].M, "mass_scaling");
  my_H5Aclose(hdf5_attribute, "mass_scaling");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(hdf5_dataset, "velocity_scaling", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &IO_Fields[ind].V, "velocity_scaling");
  my_H5Aclose(hdf5_attribute, "velocity_scaling");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(hdf5_dataset, "to_cgs", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &IO_Fields[ind].c, "to_cgs");
  my_H5Aclose(hdf5_attribute, "to_cgs");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);
}

#endif

#ifdef OUTPUT_XDMF
/*! \brief Outputs a xdmf file corresponding to this snapshot.
 *
 *  This xdmf file can be used to load the snapshot into programs like visit.
 *  This option only works with output format 3 (hdf5).
 *
 *  \param[in] fname Name of the snapshot.
 *
 *  \return void
 */
static void write_xdmf(char *const fname)
{
  FILE *f;
  char buf[MAXLEN_PATH], buf2[IO_DATASET_NAME_SIZE];
  int npresent[NTYPES];

  for(int i = 0; i < NTYPES; i++)
    npresent[i] = 0;

#ifdef OUTPUT_IN_DOUBLEPRECISION
  int prec = 8;
#else
  int prec = 4;
#endif

  file_path_sprintf(buf, "%s.xmf", fname);
  f = fopen(buf, "w");

  fprintf(f, "<?xml version=\"1.0\" ?>\n");
  fprintf(f, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
  fprintf(f, "<Xdmf Version=\"2.0\">\n");
  fprintf(f, " <Domain>");

  /* hdf5 file path relative to xmf file, uses basename function of libgen.h,
   * i.e. POSIX version of basename() */
  file_path_sprintf(buf, "./%s.hdf5", basename(fname));

  for(int type = 0; type < NTYPES; type++)
    {
      int bnr;

      for(bnr = 0; bnr < 1000; bnr++)
        {
          enum iofields i = (enum iofields)bnr;

          if(i == IO_LASTENTRY)
            break;

          if(blockpresent(i, 1))
            {
              if(header.npart[type] > 0)
                {
                  if(i == IO_POS)
                    {
                      fprintf(f, "  <Grid Name=\"PartType%d\" GridType=\"Uniform\">\n", type);
                      fprintf(f, "   <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"%d\"/>\n", header.npart[type]);
                      fprintf(f, "   <Geometry GeometryType=\"XYZ\">\n");
                      fprintf(f, "    <DataItem Dimensions=\"%d 3\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n",
                              header.npart[type], prec);
                      fprintf(f, "     %s:/PartType0/Coordinates\n", buf);
                      fprintf(f, "    </DataItem>\n");
                      fprintf(f, "   </Geometry>\n");

                      npresent[type] = 1;
                    }
                  else
                    {
                      int dim   = get_values_per_blockelement(i);
                      int dtype = get_datatype_in_block(i, 0);
                      get_dataset_name(i, buf2);

                      if(dim == 1 || dim == 3)
                        {
                          if(dtype == 1)
                            {
                              if(dim == 1)
                                {
                                  fprintf(f, "   <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Node\">\n", buf2);
                                  fprintf(f, "    <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n",
                                          header.npart[type], prec);
                                }
                              else
                                {
                                  fprintf(f, "   <Attribute Name=\"%s\" AttributeType=\"Vector\" Center=\"Node\">\n", buf2);
                                  fprintf(f,
                                          "    <DataItem Dimensions=\"%d\ 3\" NumberType=\"Float\" Precision=\"%d\" Format=\"HDF\">\n",
                                          header.npart[type], prec);
                                }

                              fprintf(f, "     %s:/PartType%d/%s\n", buf, type, buf2);
                              fprintf(f, "    </DataItem>\n");
                              fprintf(f, "   </Attribute>\n");
                            }
                        }
                    }
                }
            }
        }
      if(npresent[type] == 1)
        {
          fprintf(f, "  </Grid>\n");
        }
    }

  fprintf(f, " </Domain>\n");
  fprintf(f, "</Xdmf>");

  fclose(f);
}
#endif

/*! \brief  A wrapper for the fwrite() function.
 *
 *  This catches I/O errors occuring for fwrite(). In this case we
 *  better stop. If stream is null, no attempt at writing is done.
 *
 *  \param[in] ptr Pointer to the beginning of data to write.
 *  \param[in] size Size in bytes of a single data element.
 *  \param[in] nmemb Number of elements to be written.
 *  \param[in] stream Pointer to the output stream.
 *
 *  \return Number of elements written to stream.
 */
size_t my_fwrite_fullinfo(const void *data, const size_t size, const size_t nmemb, FILE *stream, const char *func, const char *file,
                          const int line)
{
#ifdef TOLERATE_WRITE_ERROR
  if(WriteErrorFlag)
    return 0;
#endif

  if(!stream)
    return 0;

  size_t nwritten = 0;
  if(size * nmemb > 0)
    {
      if((nwritten = fwrite(data, size, nmemb, stream)) != nmemb)
        {
#ifdef TOLERATE_WRITE_ERROR
          write_error(0, nwritten, nmemb);
#else
          terminate(
              "I/O error (fwrite) has occured: %s\n"
              "  at %s()/%s/line %d.",
              strerror(errno), func, file, line);
#endif
        }
    }

#ifdef TOLERATE_WRITE_ERROR
  if(ferror(stream))
    write_error(1, nwritten, nmemb);
#endif

  return nwritten;
}

/*! \brief  A wrapper for the fread() function.
 *
 *  This catches I/O errors occuring for fread(). In this case we
 *  better stop. If stream is null, no attempt at readingis done.
 *
 *  \param[out] ptr Pointer to the beginning of memory location where to
 *              store data.
 *  \param[in] size Size in bytes of a single data element.
 *  \param[in] nmemb Number of elements to be read.
 *  \param[in] stream Pointer to the input stream.
 *
 *  \return Number of elements read from stream.
 */
size_t my_fread_fullinfo(void *data, const size_t size, const size_t nmemb, FILE *stream, const char *func, const char *file,
                         const int line)
{
  if(!stream)
    return 0;

  size_t nread = 0;
  if(size * nmemb > 0)
    {
      if((nread = fread(data, size, nmemb, stream)) != nmemb)
        {
          char *err;
          if(feof(stream))
            err = "end of file";
          else
            err = strerror(errno);
          terminate(
              "I/O error (fread) has occured: %s\n"
              "  at %s()/%s/line %d.",
              err, func, file, line);
        }
    }

  return nread;
}

/*! \brief A wrapper for the printf() function.
 *
 *  This function has the same functionalities of the standard printf()
 *  function. However, data is written to the standard output only for
 *  the task with rank 0.
 *
 *  \param[in] fmt String that contains format arguments.
 *
 *  \return void
 */
void mpi_printf(const char *fmt, ...)
{
  if(ThisTask == 0)
    {
      va_list l;
      va_start(l, fmt);
      vprintf(fmt, l);
      myflush(stdout);
      va_end(l);
    }
}

/*! \brief A wrapper for the fprintf() function.
 *
 *  This function has the same functionalities of the standard fprintf()
 *  function. However, data is written to the standard output only for
 *  the task with rank 0.
 *
 *  \param[in,out] stream File stream to write to.
 *  \param[in] fmt String that contains format arguments.
 *  \param[in] ... Data to print according to the format string.
 *
 *  \return void
 */
void mpi_fprintf(FILE *stream, const char *fmt, ...)
{
  if(ThisTask == 0)
    {
      va_list l;
      va_start(l, fmt);
      vfprintf(stream, fmt, l);
      myflush(stream);
      va_end(l);
    }
}

/*! \brief A function for printing debug information in parallel.
 *
 *  This function works like printf(), however it takes care
 *  that the output is contigous in the stdout from task 0 to task NTask - 1.
 *  Run this debug function only in code parts which all tasks reach.
 *
 *
 *  \param[in] fmt String that contains format arguments.
 *  \param[in] ... Data to print according to the format string.
 *
 *  \return void
 */
void mpi_printf_each(const char *fmt, ...)
{
  char buffer[2048];

  va_list l;
  va_start(l, fmt);
  vsnprintf(buffer, sizeof(buffer), fmt, l);
  va_end(l);

  if(ThisTask == 0)
    {
      /* print own message */
      printf("%s", buffer);
      /* print message from other tasks */
      for(int i = 1; i < NTask; i++)
        {
          MPI_Recv(buffer, sizeof(buffer), MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          printf("%s", buffer);
        }
    }
  else
    {
      MPI_Send(buffer, strlen(buffer) + 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    }
}

/*! \brief A wrapper for the system() function.
 *
 *  This function has the same functionalities of the standard system()
 *  function. However, a warning is emitted (using #warn()) if system()
 *  returns a non-zero status code. Moreover, system() will only be
 *  called if the macro NOCALLSOFSYSTEM is unset.
 *
 *  \param[in] command String that contains the command line to be executed.
 *
 *  \return void
 */
int my_system(const char *const command)
{
  int status = 0;
#ifndef NOCALLSOFSYSTEM
  status = system(command);
  if(status != 0)
    warn("calling system(\"%s\") returned exit code %d", command, status);
#endif
  return status;
}

/*! \brief A function for creating file paths which makes sure that the length
 *         of the path does not exceed #MAXLEN_PATH.
 *
 *  This function works like sprintf(), however it checks if the length of the
 *  resulting path exceeds #MAXLEN_PATH and terminates if this is the case.
 *
 *
 *  \param[in,out] buf Buffer which the path is written to.
 *  \param[in] fmt String that contains format arguments.
 *
 *  \return void
 */
void file_path_sprintf(char *const buf, const char *const fmt, ...)
{
  va_list l;
  va_start(l, fmt);
  int status = vsnprintf(buf, MAXLEN_PATH, fmt, l);
  if(status >= MAXLEN_PATH)
    terminate("length of file path exceeds MAXLEN_PATH = %d: %s", MAXLEN_PATH, buf);
  va_end(l);
}

/*! \brief Opens the requested file name and returns the file descriptor.
 *
 *  If opening fails, an error is printed and the file descriptor is
 *  null.
 *
 *  \param[in] fname The file name.
 *
 *  \return A file descriptor to the file.
 */
FILE *open_file(const char *const fname)
{
  FILE *fd;
  if(!(fd = fopen(fname, "w")))
    printf("can't open file `%s' for writing.\n", fname);
  return fd;
}

/*! \brief Open the file descriptors needed for the output image files.
 *
 *   If the pointer to the FILE* is null, that file is not opened
 *
 *   \param id1    type of the image, ie "slice", "proj", or "grid"
 *   \param num    number of the image
 *   \param dens   pointer to the file descriptor for density
 *   \param temp   pointer to the file descriptor for temperature
 *   \param met    pointer to the file descriptor for metallicity
 *   \param vel    pointer to the file descriptor for velocity
 *   \param mag    pointer to the file descriptor for magnetic field
 *   \param vort   pointer to the file descriptor for vorticity
 *   \param phot   pointer to the file descriptor for photon density
 *   \param chem   pointer to the file descriptor for the GFM chemical elements
 *   \param denstr pointer to the file descriptor for the MC tracer density
 *   \param dust pointer to the file descriptor for the dust temperature
 *   \param h2   pointer to the file descriptor for the H2   abundance
 *   \param hp   pointer to the file descriptor for the H+   abundance
 *   \param co   pointer to the file descriptor for the CO   abundance
 *   \param chx  pointer to the file descriptor for the CHx  abundance
 *   \param ohx  pointer to the file descriptor for the OHx  abundance
 *   \param hcop pointer to the file descriptor for the HCO+ abundance
 *   \param cp   pointer to the file descriptor for the C+   abundance
 *   \param mp   pointer to the file descriptor for the M+   abundance
 *   \param hep  pointer to the file descriptor for the He+  abundance
 *   \param rih  pointer to the file descriptor for the H photoionization rate
 *   \param hrih pointer to the file descriptor for the heating rate due to H photoionization
 *   \param sxrates pointer to the file descriptor for the full set of SIMPLEX photorates
 */

void open_image_files(char *id1, int num, FILE **dens, FILE **temp, FILE **met, FILE **vel, FILE **mag, FILE **vort, FILE **phot,
                      FILE **chem, FILE **denstr, FILE **dust, FILE **h2, FILE **hp, FILE **co, FILE **chx, FILE **ohx, FILE **hcop,
                      FILE **cp, FILE **mp, FILE **hep, FILE **rih, FILE **hrih, FILE **sxrates)

{
  if(ThisTask != 0)
    return;

  char buf[MAXLEN_PATH];

  if(dens)
    {
      file_path_sprintf(buf, "%s/density_%s_%03d", All.OutputDir, id1, num);
      *dens = open_file(buf);
    }

  if(temp)
    {
      file_path_sprintf(buf, "%s/temp_%s_%03d", All.OutputDir, id1, num);
      *temp = open_file(buf);
    }

  if(met)
    {
      file_path_sprintf(buf, "%s/metal_%s_%03d", All.OutputDir, id1, num);
      *met = open_file(buf);
    }

  if(vel)
    {
      file_path_sprintf(buf, "%s/velocity_%s_%03d", All.OutputDir, id1, num);
      *vel = open_file(buf);
    }

  if(mag)
    {
      file_path_sprintf(buf, "%s/magnetic_%s_%03d", All.OutputDir, id1, num);
      *mag = open_file(buf);
    }

  if(vort)
    {
      file_path_sprintf(buf, "%s/vorticity_%s_%03d", All.OutputDir, id1, num);
      *vort = open_file(buf);
    }

  if(phot)
    {
      file_path_sprintf(buf, "%s/photondens_%s_%03d", All.OutputDir, id1, num);
      *phot = open_file(buf);
    }

  if(chem)
    {
      file_path_sprintf(buf, "%s/chemelements_%s_%03d", All.OutputDir, id1, num);
      *chem = open_file(buf);
    }

  if(denstr)
    {
      file_path_sprintf(buf, "%s/density_trmc_%s_%03d", All.OutputDir, id1, num);
      *denstr = open_file(buf);
    }

  if(dust)
    {
      file_path_sprintf(buf, "%s/tdust_%s_%03d", All.OutputDir, id1, num);
      *dust = open_file(buf);
    }

  if(h2)
    {
      file_path_sprintf(buf, "%s/xH2_%s_%03d", All.OutputDir, id1, num);
      *h2 = open_file(buf);
    }

  if(hp)
    {
      file_path_sprintf(buf, "%s/xHP_%s_%03d", All.OutputDir, id1, num);
      *hp = open_file(buf);
    }

  if(co)
    {
      file_path_sprintf(buf, "%s/xCO_%s_%03d", All.OutputDir, id1, num);
      *co = open_file(buf);
    }

  if(chx)
    {
      file_path_sprintf(buf, "%s/xCHX_%s_%03d", All.OutputDir, id1, num);
      *chx = open_file(buf);
    }

  if(ohx)
    {
      file_path_sprintf(buf, "%s/xOHX_%s_%03d", All.OutputDir, id1, num);
      *ohx = open_file(buf);
    }

  if(hcop)
    {
      file_path_sprintf(buf, "%s/xHCOP_%s_%03d", All.OutputDir, id1, num);
      *hcop = open_file(buf);
    }

  if(cp)
    {
      file_path_sprintf(buf, "%s/xCP_%s_%03d", All.OutputDir, id1, num);
      *cp = open_file(buf);
    }

  if(mp)
    {
      file_path_sprintf(buf, "%s/xMP_%s_%03d", All.OutputDir, id1, num);
      *mp = open_file(buf);
    }

  if(hep)
    {
      file_path_sprintf(buf, "%s/xHeP_%s_%03d", All.OutputDir, id1, num);
      *hep = open_file(buf);
    }

  if(rih)
    {
      file_path_sprintf(buf, "%s/rih_%s_%03d", All.OutputDir, id1, num);
      *rih = open_file(buf);
    }
  if(hrih)
    {
      file_path_sprintf(buf, "%s/hrih_%s_%03d", All.OutputDir, id1, num);
      *hrih = open_file(buf);
    }
  if(sxrates)
    {
      file_path_sprintf(buf, "%s/sxrates_%s_%03d", All.OutputDir, id1, num);
      *sxrates = open_file(buf);
    }
}

/*! \brief Write the header for the output image files.
 *
 *  The header contains the dimensions in pixel of the image
 *
 *  \param fd pointer to the image file
 *  \param nx number of pixels in the x direction
 *  \param ny number of pixels in the y direction
 *  \param nz number of pixels in the z direction
 */
void write_image_header(FILE *fd, const int nx, const int ny, const int nz)
{
  my_fwrite(&nx, sizeof(int), 1, fd);
  my_fwrite(&ny, sizeof(int), 1, fd);
  if(nz > 0)
    my_fwrite(&nz, sizeof(int), 1, fd);
}

/*! \brief Write the footer for the output image files.
 *
 *  The footer contains the units and limits of the image
 *
 *  \param fd pointer to the image file
 */
void write_image_footer(FILE *fd, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax)
{
  my_fwrite(&xmin, sizeof(double), 1, fd);
  my_fwrite(&xmax, sizeof(double), 1, fd);
  my_fwrite(&ymin, sizeof(double), 1, fd);
  my_fwrite(&ymax, sizeof(double), 1, fd);
  my_fwrite(&zmin, sizeof(double), 1, fd);
  my_fwrite(&zmax, sizeof(double), 1, fd);
  /* units */
  my_fwrite(&All.UnitMass_in_g, sizeof(double), 1, fd);
  my_fwrite(&All.UnitLength_in_cm, sizeof(double), 1, fd);
  my_fwrite(&All.UnitTime_in_s, sizeof(double), 1, fd);
}

#ifdef SUBBOX_SNAPSHOTS
/*! \brief Reads the coordinates from the subbox file.
 *
 *  \param[in] fname Name of the subbox file.
 *
 *  \return void
 */
void read_subbox_coordinates(const char *const fname)
{
  int iter, i;
  FILE *fsubboxes;
  double dummy;

  Nsubboxes = 0;

  for(iter = 0, i = 0; iter < 2; iter++)
    {
      if(!(fsubboxes = fopen(fname, "r")))
        terminate("SUBBOX_SNAPSHOTS: cannot read subboxes coordinates in file `%s'", fname);
      if(iter == 0)
        while(fscanf(fsubboxes, "%lg %lg %lg %lg %lg %lg", &dummy, &dummy, &dummy, &dummy, &dummy, &dummy) != EOF)
          Nsubboxes++;
      if(iter == 1)
        while(fscanf(fsubboxes, "%lg %lg %lg %lg %lg %lg", &(SubboxXmin[i]), &(SubboxXmax[i]), &(SubboxYmin[i]), &(SubboxYmax[i]),
                     &(SubboxZmin[i]), &(SubboxZmax[i])) != EOF)
          i++;
      fclose(fsubboxes);

      if(iter == 0)
        {
          SubboxXmin = (double *)mymalloc("SubboxXmin", Nsubboxes * sizeof(double));
          SubboxXmax = (double *)mymalloc("SubboxXmax", Nsubboxes * sizeof(double));
          SubboxYmin = (double *)mymalloc("SubboxYmin", Nsubboxes * sizeof(double));
          SubboxYmax = (double *)mymalloc("SubboxYmax", Nsubboxes * sizeof(double));
          SubboxZmin = (double *)mymalloc("SubboxZmin", Nsubboxes * sizeof(double));
          SubboxZmax = (double *)mymalloc("SubboxZmax", Nsubboxes * sizeof(double));
          mpi_printf("SUBBOX_SNAPSHOTS: read subboxes coordinates with %d entries in file `%s'.\n", Nsubboxes, fname);
        }
    }

  for(i = 0; i < Nsubboxes; i++)
    mpi_printf("SUBBOX_SNAPSHOTS: subbox #%d coordinates: %g %g %g %g %g %g\n", i, SubboxXmin[i], SubboxXmax[i], SubboxYmin[i],
               SubboxYmax[i], SubboxZmin[i], SubboxZmax[i]);
  mpi_printf("\n");
}
#endif

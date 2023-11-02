/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/fof/fof_io.c
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
#include <inttypes.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "../../allvars.h"
#include "../../domain.h"
#include "../../proto.h"
#include "../../gitversion/version.h"
#include "mergers_io.h"


#ifdef HAVE_HDF5
#include <hdf5.h>
static void merger_write_header_attributes_in_hdf5(hid_t handle);
#endif

static double wrap_position(double pos, int direction);

#ifdef STORE_MERGERS_IN_SNAPSHOT

/*! \brief Main routine for merger output.
 *
 *  \param[in] num Index of group file (snapshot index for this output).
 *
 *  \return void
 */
void save_mergers(int num)
{
  int filenr, mrg, nmrgs, masterTask, lastTask;
  double t0, t1;
  char buf[MAXLEN_PATH];

  t0 = second();

  CommBuffer = mymalloc("CommBuffer", COMMBUFFERSIZE);

  if(NTask < All.NumFilesPerSnapshot)
    {
      warn(
          "Number of processors must be larger or equal than All.NumFilesPerSnapshot! Reducing All.NumFilesPerSnapshot "
          "accordingly.\n");
      All.NumFilesPerSnapshot = NTask;
    }

  if(All.SnapFormat < SNAP_FORMAT_GADGET || All.SnapFormat > SNAP_FORMAT_HDF5)
    mpi_printf("Unsupported File-Format. All.SnapFormat=%d\n", All.SnapFormat);

#ifndef HAVE_HDF5
  if(All.SnapFormat == SNAP_FORMAT_HDF5)
    mpi_terminate("Code wasn't compiled with HDF5 support enabled!\n");
#endif

  /* assign processors to output files */
  distribute_file(All.NumFilesPerSnapshot, &filenr, &masterTask, &lastTask);

  if(All.NumFilesPerSnapshot > 1)
    {
      if(ThisTask == 0)
        {
          file_path_sprintf(buf, "%s/mergers_%03d", All.OutputDir, num);
          mkdir(buf, MKDIR_MODE);
        }
      MPI_Barrier(MPI_COMM_WORLD);
    }

  if(All.NumFilesPerSnapshot > 1)
    file_path_sprintf(buf, "%s/mergers_%03d/mergers_tab_%03d.%d", All.OutputDir, num, num, filenr);
  else
    file_path_sprintf(buf, "%s/mergers_tab_%03d", All.OutputDir, num);

  nmrgs = All.NumFilesPerSnapshot / All.NumFilesWrittenInParallel;
  if((All.NumFilesPerSnapshot % All.NumFilesWrittenInParallel))
    nmrgs++;

  for(mrg = 0; mrg < nmrgs; mrg++)
    {
      if((filenr / All.NumFilesWrittenInParallel) == mrg) /* ok, it's this processor's turn */
        merger_write_file(buf, masterTask, lastTask);

      MPI_Barrier(MPI_COMM_WORLD);
    }

  myfree(CommBuffer);
  CommBuffer = NULL;


  t1 = second();

  mpi_printf("Merger catalogues saved. took = %g sec\n", timediff(t0, t1));
}

/*! \brief Writes a file with name fname containing data from writeTask to
 *         lastTask.
 *
 *  \param[in] fname Filename of the output file.
 *  \param[in] writeTask Task responsible for writing the file.
 *  \param[in] lastTask Last task whose data is still in this file.
 *
 *  \return void
 */
void merger_write_file(char *fname, int writeTask, int lastTask)
{
  mpi_printf("merger_write_file has been called");
  int bytes_per_blockelement, npart, nextblock;
  int n_for_this_task, n, p, pc, offset = 0, task;
  int blockmaxlen, n_type_mergers, ntot_type_mergers, nn;
  enum merger_iofields blocknr;
  char label[IO_LABEL_SIZE];
  int bnr;
  int blksize;
  MPI_Status status;
  FILE *fd = 0;
#ifdef HAVE_HDF5
  hid_t hdf5_file = 0, hdf5_grp, hdf5_headergrp = 0, hdf5_dataspace_memory;
  hid_t hdf5_datatype = 0, hdf5_dataspace_in_file = 0, hdf5_dataset = 0;
  hid_t hdf5_paramsgrp = 0, hdf5_configgrp = 0;
  herr_t hdf5_status;
  hsize_t dims[2], count[2], start[2];
  int rank = 0, pcsum = 0;
  char dataset_name[IO_DATASET_NAME_SIZE];
#endif

#define SKIP                                 \
  {                                          \
    my_fwrite(&blksize, sizeof(int), 1, fd); \
  }

  /* determine group/id numbers of each type in file */

  n_type_mergers = Nmergers;

  if(ThisTask == writeTask)
    {
     // for(n = 0; n < 3; n++)
      ntot_type_mergers = n_type_mergers;
      for(task = writeTask + 1; task <= lastTask; task++)
        {
          MPI_Recv(&nn, 1, MPI_INT, task, TAG_LOCALN, MPI_COMM_WORLD, &status);
          //for(n = 0; n < 3; n+
          ntot_type_mergers += nn;
        }

      for(task = writeTask + 1; task <= lastTask; task++)
        MPI_Send(&ntot_type_mergers, 1, MPI_INT, task, TAG_N, MPI_COMM_WORLD);
    }
  else
    {
      MPI_Send(&n_type_mergers, 1, MPI_INT, writeTask, TAG_LOCALN, MPI_COMM_WORLD);
      MPI_Recv(&ntot_type_mergers, 1, MPI_INT, writeTask, TAG_N, MPI_COMM_WORLD, &status);
    }

  /* fill file header */

  MPI_Allreduce(&Nmergers, &TotNmergers, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  merger_catalogue_header.Nmergers    = ntot_type_mergers;

  merger_catalogue_header.TotNmergers    = TotNmergers;

  merger_catalogue_header.num_files = All.NumFilesPerSnapshot;

  merger_catalogue_header.time = All.Time;
  if(All.ComovingIntegrationOn)
    merger_catalogue_header.redshift = 1.0 / All.Time - 1;
  else
    merger_catalogue_header.redshift = 0;
  merger_catalogue_header.HubbleParam = All.HubbleParam;
  merger_catalogue_header.BoxSize     = All.BoxSize;
  merger_catalogue_header.Omega0      = All.Omega0;
  merger_catalogue_header.OmegaLambda = All.OmegaLambda;

#ifdef OUTPUT_IN_DOUBLEPRECISION
  merger_catalogue_header.flag_doubleprecision = 1;
#else
  merger_catalogue_header.flag_doubleprecision = 0;
#endif

  /* open file and write header */

  if(ThisTask == writeTask)
    {
      if(All.SnapFormat == SNAP_FORMAT_HDF5)
        {
#ifdef HAVE_HDF5
          char buf[MAXLEN_PATH];
          file_path_sprintf(buf, "%s.hdf5", fname);
          hdf5_file = my_H5Fcreate(buf, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
          mpi_printf("MERGER: writing merger catalogue: '%s' (file 1 of %d)\n", fname, All.NumFilesPerSnapshot);
          hdf5_headergrp = my_H5Gcreate(hdf5_file, "/Header", 0);
          hdf5_grp = my_H5Gcreate(hdf5_file, "/Merger", 0);
          merger_write_header_attributes_in_hdf5(hdf5_headergrp);
#endif
        }
      else
        {
          if(!(fd = fopen(fname, "w")))
            {
              printf("can't open file `%s' for writing snapshot.\n", fname);
              terminate("file open error");
            }

          mpi_printf("MERGERS: writing merger catalogue: '%s' (file 1 of %d)\n", fname, All.NumFilesPerSnapshot);

          if(All.SnapFormat == SNAP_FORMAT_GADGET_VARIANT)
            {
              blksize = sizeof(int) + 4 * sizeof(char);
              SKIP;
              my_fwrite((void *)"HEAD", sizeof(char), 4, fd);
              nextblock = sizeof(merger_catalogue_header) + 2 * sizeof(int);
              my_fwrite(&nextblock, sizeof(int), 1, fd);
              SKIP;
            }

          blksize = sizeof(merger_catalogue_header);
          SKIP;
          my_fwrite(&merger_catalogue_header, sizeof(merger_catalogue_header), 1, fd);
          SKIP;
        }
    }

  for(bnr = 0; bnr < 1000; bnr++)
    {
      blocknr = (enum merger_iofields)bnr;

      if(blocknr == IO_MERGER_LASTENTRY)
        break;

      if(merger_blockpresent(blocknr))
        {
          bytes_per_blockelement = merger_get_bytes_per_blockelement(blocknr);

          blockmaxlen = (int)(COMMBUFFERSIZE / bytes_per_blockelement);

          npart   = merger_get_particles_in_block(blocknr);

          if(npart > 0)
            {
              if(ThisTask == 0)
                {
                  char tmp[IO_DATASET_NAME_SIZE];

                  merger_get_dataset_name(blocknr, tmp);
                  printf("MERGER: writing block %d (%s)...\n", (int)blocknr, tmp);
                }

              if(ThisTask == writeTask)
                {
                  if(All.SnapFormat == SNAP_FORMAT_GADGET || All.SnapFormat == SNAP_FORMAT_GADGET_VARIANT)
                    {
                      if(All.SnapFormat == SNAP_FORMAT_GADGET_VARIANT)
                        {
                          blksize = sizeof(int) + IO_LABEL_SIZE * sizeof(char);
                          SKIP;
                          merger_get_Tab_IO_Label(blocknr, label);
                          my_fwrite(label, sizeof(char), IO_LABEL_SIZE, fd);
                          nextblock = npart * bytes_per_blockelement + 2 * sizeof(int);
                          my_fwrite(&nextblock, sizeof(int), 1, fd);
                          SKIP;
                        }

                      blksize = npart * bytes_per_blockelement;
                      SKIP;
                    }
                  else if(All.SnapFormat == SNAP_FORMAT_HDF5)
                    {
#ifdef HAVE_HDF5
                      switch(merger_get_datatype(blocknr))
                        {
                          case 0:
                            hdf5_datatype = my_H5Tcopy(H5T_NATIVE_INT);
                            break;
                          case 1:
#ifdef OUTPUT_IN_DOUBLEPRECISION
                            hdf5_datatype = my_H5Tcopy(H5T_NATIVE_DOUBLE);
#else
                            hdf5_datatype = my_H5Tcopy(H5T_NATIVE_FLOAT);
#endif
                            break;
                          case 2:
                            hdf5_datatype = my_H5Tcopy(H5T_NATIVE_UINT64);
                            break;
                        }

                      dims[0] = ntot_type_mergers;
                      dims[1] = merger_get_values_per_blockelement(blocknr);
                      if(dims[1] == 1)
                        rank = 1;
                      else
                        rank = 2;

                      merger_get_dataset_name(blocknr, dataset_name);

                      hdf5_dataspace_in_file = my_H5Screate_simple(rank, dims, NULL);

                      hdf5_dataset = my_H5Dcreate(hdf5_grp, dataset_name, hdf5_datatype, hdf5_dataspace_in_file, H5P_DEFAULT);

                      pcsum = 0;
#endif
                    }
                }

              for(task = writeTask, offset = 0; task <= lastTask; task++)
                {
                  if(task == ThisTask)
                    {
                      n_for_this_task = n_type_mergers;

                      for(p = writeTask; p <= lastTask; p++)
                        if(p != ThisTask)
                          MPI_Send(&n_for_this_task, 1, MPI_INT, p, TAG_NFORTHISTASK, MPI_COMM_WORLD);
                    }
                  else
                    MPI_Recv(&n_for_this_task, 1, MPI_INT, task, TAG_NFORTHISTASK, MPI_COMM_WORLD, &status);

                  while(n_for_this_task > 0)
                    {
                      pc = n_for_this_task;

                      if(pc > blockmaxlen)
                        pc = blockmaxlen;

                      if(ThisTask == task)
                        merger_fill_write_buffer(blocknr, &offset, pc);

                      if(ThisTask == writeTask && task != writeTask)
                        MPI_Recv(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, task, TAG_PDATA, MPI_COMM_WORLD, &status);

                      if(ThisTask != writeTask && task == ThisTask)
                        MPI_Ssend(CommBuffer, bytes_per_blockelement * pc, MPI_BYTE, writeTask, TAG_PDATA, MPI_COMM_WORLD);

                      if(ThisTask == writeTask)
                        {
                          if(All.SnapFormat == SNAP_FORMAT_HDF5)
                            {
#ifdef HAVE_HDF5
                              start[0] = pcsum;
                              start[1] = 0;

                              count[0] = pc;
                              count[1] = merger_get_values_per_blockelement(blocknr);
                              pcsum += pc;

                              my_H5Sselect_hyperslab(hdf5_dataspace_in_file, H5S_SELECT_SET, start, NULL, count, NULL);

                              dims[0]               = pc;
                              dims[1]               = merger_get_values_per_blockelement(blocknr);
                              hdf5_dataspace_memory = my_H5Screate_simple(rank, dims, NULL);

                              hdf5_status = my_H5Dwrite(hdf5_dataset, hdf5_datatype, hdf5_dataspace_memory, hdf5_dataspace_in_file,
                                                        H5P_DEFAULT, CommBuffer, dataset_name);

                              (void)hdf5_status;

                              my_H5Sclose(hdf5_dataspace_memory, H5S_SIMPLE);
#endif
                            }
                          else
                            {
                              my_fwrite(CommBuffer, bytes_per_blockelement, pc, fd);
                            }
                        }

                      n_for_this_task -= pc;
                    }
                }

              if(ThisTask == writeTask)
                {
                  if(All.SnapFormat == 3)
                    {
#ifdef HAVE_HDF5
                      my_H5Dclose(hdf5_dataset, dataset_name);
                      my_H5Sclose(hdf5_dataspace_in_file, H5S_SIMPLE);
                      my_H5Tclose(hdf5_datatype);
#endif
                    }
                  else
                    SKIP;
                }
            }
        }
    }

  if(ThisTask == writeTask)
    {
      if(All.SnapFormat == SNAP_FORMAT_HDF5)
        {
#ifdef HAVE_HDF5
          my_H5Gclose(hdf5_grp, "/Mergers");
          my_H5Gclose(hdf5_headergrp, "/Header");
          my_H5Fclose(hdf5_file, fname);
#endif
        }
      else
        fclose(fd);
    }
}

/*! \brief Copies data from global group array to appropriate output buffer.
 *
 *  \param[in] blocknr Number (identifier) of the field to be written.
 *  \param[in] startindex First particle index to be included.
 *  \param[in] pc Particle count; number of particles to be written.
 *
 *  \return void
 */
void merger_fill_write_buffer(enum merger_iofields blocknr, int *startindex, int pc)
{
  mpi_printf("merger_fill_write_buffer has been called with pc : %d \n", pc);

  int n, k, pindex, *ip;
  MyOutputFloat *fp;
  MyIDType *idp;

  fp  = (MyOutputFloat *)CommBuffer;
  ip  = (int *)CommBuffer;
  idp = (MyIDType *)CommBuffer;

//#ifdef FOF_FUZZ_SORT_BY_NEAREST_GROUP
  long long *llp = (long long *)CommBuffer;
//#endif

  pindex = *startindex;

  for(n = 0; n < pc; pindex++, n++)
    {
      switch(blocknr)
        {
//          case IO_FOF_LEN:
//            *ip++ = Group[pindex].Len;
//            break;
//          case IO_FOF_MTOT:
//            *fp++ = Group[pindex].Mass;
//            break;
//          case IO_FOF_POS:
//            for(k = 0; k < 3; k++)
//#ifdef SUBFIND
//              *fp++ = wrap_position(Group[pindex].Pos[k] - All.GlobalDisplacementVector[k], k);
//#else
//              *fp++ = wrap_position(Group[pindex].CM[k] - All.GlobalDisplacementVector[k], k);
//#endif
//            break;

          case IO_MERGER_ID1:
   //         mpi_printf("Printing ID1: %llu with pindex \n: %d", pindex, MergerEvents[pindex].ID1);
            *llp++ = MergerEvents[pindex].ID1;
            break;
          case IO_MERGER_ID2:
   //         mpi_printf("Printing ID2: %llu with pindex \n: %d", pindex, MergerEvents[pindex].ID2);  
            *llp++ = MergerEvents[pindex].ID2;
            break;
          case IO_MERGER_TASKID:
            *ip++ = MergerEvents[pindex].TaskID;
            break;
          case IO_MERGER_TIME:
            *fp++ = MergerEvents[pindex].Time;
            break;
          case IO_MERGER_BHMASS1:
            *fp++ = MergerEvents[pindex].BH_Mass1;
            break;
          case IO_MERGER_BHMASS2:
            *fp++ = MergerEvents[pindex].BH_Mass2;
            break;
#ifdef OUTPUT_HOST_PROPERTIES_FOR_BH_MERGERS
          case IO_MERGER_BHMDOT1:
            *fp++ = MergerEvents[pindex].BH_Mdot1;
            break;
          case IO_MERGER_BHMDOT2:
            *fp++ = MergerEvents[pindex].BH_Mdot2;
            break;
          case IO_MERGER_BHHOSTHALOMASS1:
            *fp++ = MergerEvents[pindex].BH_HostHaloMass1;
            break;
          case IO_MERGER_BHHOSTHALOMASS2:
            *fp++ = MergerEvents[pindex].BH_HostHaloMass2;
            break;
          case IO_MERGER_BHHOSTSTELLARMASS1:
            *fp++ = MergerEvents[pindex].BH_HostStellarMass1;
            break;
          case IO_MERGER_BHHOSTSTELLARMASS2:
            *fp++ = MergerEvents[pindex].BH_HostStellarMass2;
            break;
          case IO_MERGER_BHHOSTGASMASS1:
            *fp++ = MergerEvents[pindex].BH_HostGasMass1;
            break;
          case IO_MERGER_BHHOSTGASMASS2:
            *fp++ = MergerEvents[pindex].BH_HostGasMass2;
            break;
          case IO_MERGER_BHHOSTSFR1:
            *fp++ = MergerEvents[pindex].BH_HostSFR1;
            break;
          case IO_MERGER_BHHOSTSFR2:
            *fp++ = MergerEvents[pindex].BH_HostSFR2;
            break;
#endif
        }
   }
}

/*! \brief Associates the output variable blocknumber with its name.
 *
 *  \param[in] blocknr Number (identifier) of the field to be written.
 *  \param[out] label Name of field.
 *
 *  \return void
 */
void merger_get_dataset_name(enum merger_iofields blocknr, char *label)
{
  switch(blocknr)
    {
      case IO_MERGER_ID1:
        strcpy(label, "ID1");
        break;
      case IO_MERGER_ID2:
        strcpy(label, "ID2");
        break;
      case IO_MERGER_TASKID:
        strcpy(label, "TaskID");
        break;
      case IO_MERGER_TIME:
        strcpy(label, "Time");
        break;
      case IO_MERGER_BHMASS1:
        strcpy(label, "BH_Mass1");
        break;
      case IO_MERGER_BHMASS2:
        strcpy(label, "BH_Mass2");
        break;
#ifdef OUTPUT_HOST_PROPERTIES_FOR_BH_MERGERS
      case IO_MERGER_BHMDOT1:
        strcpy(label, "BH_Mdot1");
        break;
      case IO_MERGER_BHMDOT2:
        strcpy(label, "BH_Mdot2");
        break;
      case IO_MERGER_BHHOSTHALOMASS1:
        strcpy(label, "BH_HostHaloMass1");
        break;
      case IO_MERGER_BHHOSTHALOMASS2:
        strcpy(label, "BH_HostHaloMass2");
        break;
      case IO_MERGER_BHHOSTSTELLARMASS1:
        strcpy(label, "BH_HostStellarMass1");
        break;
      case IO_MERGER_BHHOSTSTELLARMASS2:
        strcpy(label, "BH_HostStellarMass2");
        break;
      case IO_MERGER_BHHOSTGASMASS1:
        strcpy(label, "BH_HostGasMass1");
        break;
      case IO_MERGER_BHHOSTGASMASS2:
        strcpy(label, "BH_HostGasMass2");
        break;
      case IO_MERGER_BHHOSTSFR1:
        strcpy(label, "BH_HostSFR1");
        break;
      case IO_MERGER_BHHOSTSFR2:
        strcpy(label, "BH_HostSFR2");
        break;
#endif
      case IO_MERGER_LASTENTRY:
        terminate("should not be reached");
        break;
    }
}

/*! \brief Is this output field a group or subhalo property?
 *
 *  \param[in] blocknr Number (identifier) of the field to be written.
 *
 *  \return 0: merger property;
 */
int merger_get_dataset_group(enum merger_iofields blocknr)
{
  switch(blocknr)
    {
      case IO_MERGER_ID1:
      case IO_MERGER_ID2:
      case IO_MERGER_TASKID:
      case IO_MERGER_TIME:
      case IO_MERGER_BHMASS1:
      case IO_MERGER_BHMASS2:
#ifdef OUTPUT_HOST_PROPERTIES_FOR_BH_MERGERS
      case IO_MERGER_BHMDOT1:
      case IO_MERGER_BHMDOT2:
      case IO_MERGER_BHHOSTHALOMASS1:
      case IO_MERGER_BHHOSTHALOMASS2:
      case IO_MERGER_BHHOSTSTELLARMASS1:
      case IO_MERGER_BHHOSTSTELLARMASS2:
      case IO_MERGER_BHHOSTGASMASS1:
      case IO_MERGER_BHHOSTGASMASS2:
      case IO_MERGER_BHHOSTSFR1:
      case IO_MERGER_BHHOSTSFR2:
#endif
        return 0;
        break;
      case IO_MERGER_LASTENTRY:
        terminate("reached last entry in switch - strange.");
        break;
    }
  terminate("reached end of function - this should not happen");
  return 0;
}

/*! \brief Returns number of particles of specific field.
 *
 *  \param[in] blocknr Number (identifier) of the field to be written.
 *
 *  \return Number of entries of this property.
 */
int merger_get_particles_in_block(enum merger_iofields blocknr)
{
  switch(blocknr)
    {
      case IO_MERGER_ID1:
      case IO_MERGER_ID2:
      case IO_MERGER_TASKID:
      case IO_MERGER_TIME:
      case IO_MERGER_BHMASS1:
      case IO_MERGER_BHMASS2:
#ifdef OUTPUT_HOST_PROPERTIES_FOR_BH_MERGERS
      case IO_MERGER_BHMDOT1:
      case IO_MERGER_BHMDOT2:
      case IO_MERGER_BHHOSTHALOMASS1:
      case IO_MERGER_BHHOSTHALOMASS2:
      case IO_MERGER_BHHOSTSTELLARMASS1:
      case IO_MERGER_BHHOSTSTELLARMASS2:
      case IO_MERGER_BHHOSTGASMASS1:
      case IO_MERGER_BHHOSTGASMASS2:
      case IO_MERGER_BHHOSTSFR1:
      case IO_MERGER_BHHOSTSFR2:
#endif
        return merger_catalogue_header.Nmergers;
        break;
      case IO_MERGER_LASTENTRY:
        terminate("reached last entry in switch - strange.");
        break;
    }
  terminate("reached end of function - this should not happen");
  return 0;
}

/*! \brief Returns the number of elements per entry of a given property.
 *
 *  \param[in] blocknr Number (identifier) of the field to be written.
 *
 *  \return Number of values per element of the specified property.
 */
int merger_get_values_per_blockelement(enum merger_iofields blocknr)
{
  int values = 0;
  switch(blocknr)
    {
      case IO_MERGER_ID1:
      case IO_MERGER_ID2:
      case IO_MERGER_TASKID:
      case IO_MERGER_TIME:
      case IO_MERGER_BHMASS1:
      case IO_MERGER_BHMASS2:
#ifdef OUTPUT_HOST_PROPERTIES_FOR_BH_MERGERS
      case IO_MERGER_BHMDOT1:
      case IO_MERGER_BHMDOT2:
      case IO_MERGER_BHHOSTHALOMASS1:
      case IO_MERGER_BHHOSTHALOMASS2:
      case IO_MERGER_BHHOSTSTELLARMASS1:
      case IO_MERGER_BHHOSTSTELLARMASS2:
      case IO_MERGER_BHHOSTGASMASS1:
      case IO_MERGER_BHHOSTGASMASS2:
      case IO_MERGER_BHHOSTSFR1:
      case IO_MERGER_BHHOSTSFR2:
#endif
        values = 1;
        break;
      case IO_MERGER_LASTENTRY:
        terminate("reached last entry in switch - should not get here");
        break;
    }
  return values;
}

/*! \brief Returns the number of bytes per element of a given property.
 *
 *  \param[in] blocknr Number (identifier) of the field to be written.
 *
 *  \return Number of bytes per element for this property.
 */
int merger_get_bytes_per_blockelement(enum merger_iofields blocknr)
{
  int bytes_per_blockelement = 0;
  switch(blocknr)
    {
      case IO_MERGER_ID1:
      case IO_MERGER_ID2:
        bytes_per_blockelement = sizeof(long long);
        break;
      case IO_MERGER_TASKID:
        bytes_per_blockelement = sizeof(int);
        break;
      case IO_MERGER_TIME:
      case IO_MERGER_BHMASS1:
      case IO_MERGER_BHMASS2:
#ifdef OUTPUT_HOST_PROPERTIES_FOR_BH_MERGERS
      case IO_MERGER_BHMDOT1:
      case IO_MERGER_BHMDOT2:
      case IO_MERGER_BHHOSTHALOMASS1:
      case IO_MERGER_BHHOSTHALOMASS2:
      case IO_MERGER_BHHOSTSTELLARMASS1:
      case IO_MERGER_BHHOSTSTELLARMASS2:
      case IO_MERGER_BHHOSTGASMASS1:
      case IO_MERGER_BHHOSTGASMASS2:
      case IO_MERGER_BHHOSTSFR1:
      case IO_MERGER_BHHOSTSFR2:
#endif
        bytes_per_blockelement = sizeof(MyFloat);
        break;
      case IO_MERGER_LASTENTRY:
        terminate("reached last entry in switch - should not get here");
        break;
    }
  return bytes_per_blockelement;
}

/*! \brief Returns key for datatype of element of a given property.
 *
 *  \param[in] blocknr Number (identifier) of the field to be written.
 *
 *  \return Key for datatype: 0: int, 1: (output)float, 2: long long.
 */
int merger_get_datatype(enum merger_iofields blocknr)
{
  int typekey = 0;

  switch(blocknr)
    {
      case IO_MERGER_ID1:
      case IO_MERGER_ID2:
        typekey = 2; /* native long long */
        break;
      case IO_MERGER_TASKID:
        typekey = 0;
        break;
      case IO_MERGER_TIME:
      case IO_MERGER_BHMASS1:
      case IO_MERGER_BHMASS2:
#ifdef OUTPUT_HOST_PROPERTIES_FOR_BH_MERGERS
      case IO_MERGER_BHMDOT1:
      case IO_MERGER_BHMDOT2:
      case IO_MERGER_BHHOSTHALOMASS1:
      case IO_MERGER_BHHOSTHALOMASS2:
      case IO_MERGER_BHHOSTSTELLARMASS1:
      case IO_MERGER_BHHOSTSTELLARMASS2:
      case IO_MERGER_BHHOSTGASMASS1:
      case IO_MERGER_BHHOSTGASMASS2:
      case IO_MERGER_BHHOSTSFR1:
      case IO_MERGER_BHHOSTSFR2:
#endif
        typekey = 1;
        break;
      case IO_MERGER_LASTENTRY:
        terminate("should not be reached");
        break;
    }

  return typekey;
}

/*! \brief Determines if block is present in the current code configuration.
 *
 *  \param[in] blocknr Number (identifier) of the field to be written.
 *
 *  \return 0: not present; 1: present.
 */
int merger_blockpresent(enum merger_iofields blocknr)
{
  int present = 0;

  switch(blocknr)
    {
      case IO_MERGER_ID1:
      case IO_MERGER_ID2:
      case IO_MERGER_TASKID:
      case IO_MERGER_TIME:
      case IO_MERGER_BHMASS1:
      case IO_MERGER_BHMASS2:
#ifdef OUTPUT_HOST_PROPERTIES_FOR_BH_MERGERS
      case IO_MERGER_BHMDOT1:
      case IO_MERGER_BHMDOT2:
      case IO_MERGER_BHHOSTHALOMASS1:
      case IO_MERGER_BHHOSTHALOMASS2:
      case IO_MERGER_BHHOSTSTELLARMASS1:
      case IO_MERGER_BHHOSTSTELLARMASS2:
      case IO_MERGER_BHHOSTGASMASS1:
      case IO_MERGER_BHHOSTGASMASS2:
      case IO_MERGER_BHHOSTSFR1:
      case IO_MERGER_BHHOSTSFR2:
#endif
        present = 1;
        break;
      case IO_MERGER_LASTENTRY:
        terminate("should not be reached");
        break;
    }
  return present;
}

/*! \brief Get the 4 letter IO label for a given output field.
 *
 *  \param[in] blocknr Number (identifier) of the field to be written.
 *  \param[out] label String with the label.
 *
 *  \return void
 */
void merger_get_Tab_IO_Label(enum merger_iofields blocknr, char *label)
{
  switch(blocknr)
    {
      case IO_MERGER_ID1:
        memcpy(label, "MONE", IO_LABEL_SIZE);
        break;
      case IO_MERGER_ID2:
        memcpy(label, "MTWO", IO_LABEL_SIZE);
        break;
      case IO_MERGER_TASKID:
        memcpy(label, "MTID", IO_LABEL_SIZE);
        break;
      case IO_MERGER_TIME:
        memcpy(label, "MTME", IO_LABEL_SIZE);
        break;
      case IO_MERGER_BHMASS1:
        memcpy(label, "MAS1", IO_LABEL_SIZE);
        break;
      case IO_MERGER_BHMASS2:
        memcpy(label, "MAS2", IO_LABEL_SIZE);
        break;
#ifdef OUTPUT_HOST_PROPERTIES_FOR_BH_MERGERS
      case IO_MERGER_BHMDOT1:
        memcpy(label, "MDT1", IO_LABEL_SIZE);
        break;
      case IO_MERGER_BHMDOT2:
        memcpy(label, "MDT2", IO_LABEL_SIZE);
        break;
      case IO_MERGER_BHHOSTHALOMASS1:
        memcpy(label, "HHM1", IO_LABEL_SIZE);
        break;
      case IO_MERGER_BHHOSTHALOMASS2:
        memcpy(label, "HHM2", IO_LABEL_SIZE);
        break;
      case IO_MERGER_BHHOSTSTELLARMASS1:
        memcpy(label, "HSM1", IO_LABEL_SIZE);
        break;
      case IO_MERGER_BHHOSTSTELLARMASS2:
        memcpy(label, "HSM2", IO_LABEL_SIZE);
        break;
      case IO_MERGER_BHHOSTGASMASS1:
        memcpy(label, "HGM1", IO_LABEL_SIZE);
        break;
      case IO_MERGER_BHHOSTGASMASS2:
        memcpy(label, "HGM2", IO_LABEL_SIZE);
        break;
      case IO_MERGER_BHHOSTSFR1:
        memcpy(label, "HSF1", IO_LABEL_SIZE);
        break;
      case IO_MERGER_BHHOSTSFR2:
        memcpy(label, "HSF2", IO_LABEL_SIZE);
        break;
#endif
      case IO_MERGER_LASTENTRY:
        terminate("should not be reached");
        break;
   }
}

#ifdef HAVE_HDF5
void merger_write_header_attributes_in_hdf5(hid_t handle)
{
  hid_t hdf5_dataspace, hdf5_attribute;
  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Nmergers_ThisFile", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &merger_catalogue_header.Nmergers, "Nmergers_ThisFile");
  my_H5Aclose(hdf5_attribute, "Nmergers_ThisFile");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "NumFiles", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &merger_catalogue_header.num_files, "NumFiles");
  my_H5Aclose(hdf5_attribute, "NumFiles");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Time", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &merger_catalogue_header.time, "Time");
  my_H5Aclose(hdf5_attribute, "Time");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Redshift", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &merger_catalogue_header.redshift, "Redshift");
  my_H5Aclose(hdf5_attribute, "Redshift");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "HubbleParam", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &merger_catalogue_header.HubbleParam, "HubbleParam");
  my_H5Aclose(hdf5_attribute, "HubbleParam");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "BoxSize", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &merger_catalogue_header.BoxSize, "BoxSize");
  my_H5Aclose(hdf5_attribute, "BoxSize");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "Omega0", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &merger_catalogue_header.Omega0, "Omega0");
  my_H5Aclose(hdf5_attribute, "Omega0");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "OmegaLambda", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &merger_catalogue_header.OmegaLambda, "OmegaLambda");
  my_H5Aclose(hdf5_attribute, "OmegaLambda");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

  hdf5_dataspace = my_H5Screate(H5S_SCALAR);
  hdf5_attribute = my_H5Acreate(handle, "FlagDoubleprecision", H5T_NATIVE_INT, hdf5_dataspace, H5P_DEFAULT);
  my_H5Awrite(hdf5_attribute, H5T_NATIVE_INT, &merger_catalogue_header.flag_doubleprecision, "FlagDoubleprecision");
  my_H5Aclose(hdf5_attribute, "FlagDoubleprecision");
  my_H5Sclose(hdf5_dataspace, H5S_SCALAR);
}
#endif

static double wrap_merger_position(double pos, int dim)
{
#if defined(REFLECTIVE_X)
  if(dim == 0)
    return pos;
#endif

#if defined(REFLECTIVE_Y)
  if(dim == 1)
    return pos;
#endif

#if defined(REFLECTIVE_Z)
  if(dim == 2)
    return pos;
#endif

  double boxsize = All.BoxSize;

#ifdef LONG_X
  if(dim == 0)
    boxsize *= LONG_X;
#endif
#ifdef LONG_Y
  if(dim == 1)
    boxsize *= LONG_Y;
#endif
#ifdef LONG_Z
  if(dim == 2)
    boxsize *= LONG_Z;
#endif

  while(pos < 0)
    pos += boxsize;

  while(pos >= boxsize)
    pos -= boxsize;

  return pos;
}

#endif

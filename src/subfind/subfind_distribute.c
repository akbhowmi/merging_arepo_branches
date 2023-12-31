/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/subfind/subfind_distribute.c
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

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../allvars.h"
#include "../fof/fof.h"
#include "../proto.h"
#include "subfind.h"

#ifdef SUBFIND

static struct group_properties *send_Group
#ifdef SEED_HALO_ENVIRONMENT_CRITERION
,*send_Group2
#endif
;

/*! \brief Distributes groups equally on MPI tasks.
 *
 *  \return void
 */
void subfind_distribute_groups(void)
{
  int i, nexport, nimport, target, ngrp, recvTask;

  /* count how many we have of each task */
  for(i = 0; i < NTask; i++)
    Send_count[i] = 0;

  for(i = 0; i < Ngroups; i++)
    {
      target = Group[i].TargetTask;

      if(target < 0 || target >= NTask)
        terminate("target < 0 || target >= NTask");

      if(target != ThisTask)
        Send_count[target]++;
    }

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(i = 0, nexport = 0, nimport = 0, Recv_offset[0] = Send_offset[0] = 0; i < NTask; i++)
    {
      nimport += Recv_count[i];
      nexport += Send_count[i];

      if(i > 0)
        {
          Send_offset[i] = Send_offset[i - 1] + Send_count[i - 1];
          Recv_offset[i] = Recv_offset[i - 1] + Recv_count[i - 1];
        }
    }

  send_Group = (struct group_properties *)mymalloc_movable(&send_Group, "send_Group", nexport * sizeof(struct group_properties));

#ifdef SEED_HALO_ENVIRONMENT_CRITERION
  send_Group2 = (struct group_properties *)mymalloc_movable(&send_Group2, "send_Group2", nexport * sizeof(struct group_properties));
#endif

  for(i = 0; i < NTask; i++)
    Send_count[i] = 0;


  for(i = 0; i < Ngroups; i++)
    {
      target = Group[i].TargetTask;

      if(target != ThisTask)
        {
          send_Group[Send_offset[target] + Send_count[target]] = Group[i];
          Group[i] = Group[Ngroups - 1];
#ifdef SEED_HALO_ENVIRONMENT_CRITERION
          send_Group2[Send_offset[target] + Send_count[target]] = Group2[i];
          Group2[i] = Group2[Ngroups - 1];
#endif
          Send_count[target]++;
          Ngroups--;
          i--;
        }
    }

  if(Ngroups + nimport > MaxNgroups)
    {
#ifdef VERBOSE
      printf("SUBFIND: Task=%d: (Ngroups=%d) + (nimport=%d) > (MaxNgroups=%d). Will increase MaxNgroups.\n", ThisTask, Ngroups,
             nimport, MaxNgroups);
#endif
      MaxNgroups = Ngroups + nimport;
      Group      = (struct group_properties *)myrealloc_movable(Group, sizeof(struct group_properties) * MaxNgroups);
#ifdef SEED_HALO_ENVIRONMENT_CRITERION
      Group2	 = (struct group_properties *)myrealloc_movable(Group2, sizeof(struct group_properties) * MaxNgroups);
#endif
    }

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the group info */
              MPI_Sendrecv(&send_Group[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct group_properties), MPI_BYTE,
                           recvTask, TAG_DENS_A, &Group[Ngroups + Recv_offset[recvTask]],
                           Recv_count[recvTask] * sizeof(struct group_properties), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD,
                           MPI_STATUS_IGNORE);
#ifdef SEED_HALO_ENVIRONMENT_CRITERION
              MPI_Sendrecv(&send_Group2[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct group_properties), MPI_BYTE,
                           recvTask, TAG_DENS_A2, &Group2[Ngroups + Recv_offset[recvTask]],
                           Recv_count[recvTask] * sizeof(struct group_properties), MPI_BYTE, recvTask, TAG_DENS_A2, MPI_COMM_WORLD,
                           MPI_STATUS_IGNORE);
#endif

            }
        }
    }

  Ngroups += nimport;
#ifdef SEED_HALO_ENVIRONMENT_CRITERION
  myfree_movable(send_Group2);
#endif
  myfree_movable(send_Group);


}

static struct particle_data *partBuf;
static struct subfind_data *subBuf;

/* \brief Distributes particles on MPI tasks.
 *
 *  This function redistributes the particles in P[] and PS[] according to what
 *  is stored in PS[].TargetTask, and PS[].TargetIndex. NOTE: The associated
 *  SphP[] is not moved, i.e. the association is broken until the particles are
 *  moved back into the original order!
 *
 *  \param[in] Communicator MPI communicator.
 *
 *  \return void
 */
void subfind_distribute_particles(MPI_Comm Communicator)
{
  int nimport, nexport;
  int i, j, n, ngrp, target;
  int max_load, load;
  int CommThisTask, CommNTask;

  MPI_Comm_size(Communicator, &CommNTask);
  MPI_Comm_rank(Communicator, &CommThisTask);

  for(n = 0; n < CommNTask; n++)
    Send_count[n] = 0;

  for(n = 0; n < NumPart; n++)
    {
      target = PS[n].TargetTask;

      if(target != CommThisTask)
        {
          if(target < 0 || target >= CommNTask)
            terminate("n=%d targettask=%d", n, target);

          Send_count[target]++;
        }
    }

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, Communicator);

  for(j = 0, nimport = 0, nexport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < CommNTask; j++)
    {
      nexport += Send_count[j];
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  /* for resize */
  load = (NumPart + nimport - nexport);
  MPI_Allreduce(&load, &max_load, 1, MPI_INT, MPI_MAX, Communicator);

  partBuf = (struct particle_data *)mymalloc_movable(&partBuf, "partBuf", nexport * sizeof(struct particle_data));
  subBuf  = (struct subfind_data *)mymalloc_movable(&subBuf, "subBuf", nexport * sizeof(struct subfind_data));

  for(i = 0; i < CommNTask; i++)
    Send_count[i] = 0;

  for(n = 0; n < NumPart; n++)
    {
      target = PS[n].TargetTask;

      if(target != CommThisTask)
        {
          partBuf[Send_offset[target] + Send_count[target]] = P[n];
          subBuf[Send_offset[target] + Send_count[target]]  = PS[n];

          P[n]  = P[NumPart - 1];
          PS[n] = PS[NumPart - 1];

          Send_count[target]++;
          NumPart--;
          n--;
        }
    }

  /* do resize */
  if(max_load > (1.0 - ALLOC_TOLERANCE) * All.MaxPart)
    {
      All.MaxPart = max_load / (1.0 - 2 * ALLOC_TOLERANCE);
      reallocate_memory_maxpart_ignore_timebins();
      PS = (struct subfind_data *)myrealloc_movable(PS, All.MaxPart * sizeof(struct subfind_data));
    }

  for(i = 0; i < CommNTask; i++)
    Recv_offset[i] += NumPart;

#ifndef NO_ISEND_IRECV_IN_DOMAIN

  MPI_Request *requests = (MPI_Request *)mymalloc("requests", 8 * CommNTask * sizeof(MPI_Request));
  int n_requests        = 0;

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      target = CommThisTask ^ ngrp;

      if(target < CommNTask)
        {
          if(Recv_count[target] > 0)
            {
              MPI_Irecv(P + Recv_offset[target], Recv_count[target] * sizeof(struct particle_data), MPI_BYTE, target, TAG_PDATA,
                        Communicator, &requests[n_requests++]);
              MPI_Irecv(PS + Recv_offset[target], Recv_count[target] * sizeof(struct subfind_data), MPI_BYTE, target, TAG_KEY,
                        Communicator, &requests[n_requests++]);
            }
        }
    }

  MPI_Barrier(Communicator); /* not really necessary, but this will guarantee that all receives are
                                posted before the sends, which helps the stability of MPI on
                                bluegene, and perhaps some mpich1-clusters */

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      target = CommThisTask ^ ngrp;

      if(target < CommNTask)
        {
          if(Send_count[target] > 0)
            {
              MPI_Isend(partBuf + Send_offset[target], Send_count[target] * sizeof(struct particle_data), MPI_BYTE, target, TAG_PDATA,
                        Communicator, &requests[n_requests++]);
              MPI_Isend(subBuf + Send_offset[target], Send_count[target] * sizeof(struct subfind_data), MPI_BYTE, target, TAG_KEY,
                        Communicator, &requests[n_requests++]);
            }
        }
    }

  MPI_Waitall(n_requests, requests, MPI_STATUSES_IGNORE);
  myfree(requests);

#else
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      target = CommThisTask ^ ngrp;

      if(target < CommNTask)
        {
          if(Send_count[target] > 0 || Recv_count[target] > 0)
            {
              MPI_Sendrecv(partBuf + Send_offset[target], Send_count[target] * sizeof(struct particle_data), MPI_BYTE, target,
                           TAG_PDATA, P + Recv_offset[target], Recv_count[target] * sizeof(struct particle_data), MPI_BYTE, target,
                           TAG_PDATA, Communicator, MPI_STATUS_IGNORE);

              MPI_Sendrecv(subBuf + Send_offset[target], Send_count[target] * sizeof(struct subfind_data), MPI_BYTE, target, TAG_KEY,
                           PS + Recv_offset[target], Recv_count[target] * sizeof(struct subfind_data), MPI_BYTE, target, TAG_KEY,
                           Communicator, MPI_STATUS_IGNORE);
            }
        }
    }
#endif

  NumPart += nimport;
  myfree_movable(subBuf);
  myfree_movable(partBuf);

  /* finally, let's also address the desired local order according to PS[].TargetIndex */

  struct fof_local_sort_data *sd;
  int *Id;

  sd = (struct fof_local_sort_data *)mymalloc("mp", sizeof(struct fof_local_sort_data) * (NumPart));
  Id = (int *)mymalloc("Id", sizeof(int) * (NumPart));

  for(i = 0; i < NumPart; i++)
    {
      sd[i].index       = i;
      sd[i].targetindex = PS[i].TargetIndex;
    }

  qsort(sd, NumPart, sizeof(struct fof_local_sort_data), fof_compare_local_sort_data_targetindex);

  for(i = 0; i < NumPart; i++)
    Id[sd[i].index] = i;

  subfind_reorder_P(Id, 0, NumPart);

  for(i = 0; i < NumPart; i++)
    Id[sd[i].index] = i;

  subfind_reorder_PS(Id, 0, NumPart);

  myfree(Id);
  myfree(sd);
}

/*! \brief Reorders elements in the P array.
 *
 * \param[in] Id Array containing ordering.
 * \param[in] Nstart Start index (in Id and P).
 * \param[in] N Final element index + 1.
 *
 *  \return void
 */
void subfind_reorder_P(int *Id, int Nstart, int N)
{
  int i;
  struct particle_data Psave, Psource;
  int idsource, idsave, dest;

  for(i = Nstart; i < N; i++)
    {
      if(Id[i] != i)
        {
          Psource  = P[i];
          idsource = Id[i];

          dest = Id[i];

          do
            {
              Psave  = P[dest];
              idsave = Id[dest];

              P[dest]  = Psource;
              Id[dest] = idsource;

#ifdef GFM
              if(P[dest].Type == 4)
                StarP[P[dest].AuxDataID].PID = dest;
#endif
#ifdef BLACK_HOLES
              if(P[dest].Type == 5)
                BHP[P[dest].AuxDataID].PID = dest;
#endif
#ifdef DUST_LIVE
              if(P[dest].Type == DUST_LIVE)
                DustP[P[dest].AuxDataID].PID = dest;
#endif

              if(dest == i)
                break;

              Psource  = Psave;
              idsource = idsave;

              dest = idsource;
            }
          while(1);
        }
    }
}

/*! \brief Reorders elements in the PS array.
 *
 * \param[in] Id Array containing ordering.
 * \param[in] Nstart Start index (in Id and P).
 * \param[in] N Final element index + 1.
 *
 *  \return void
 */
void subfind_reorder_PS(int *Id, int Nstart, int N)
{
  int i;
  struct subfind_data PSsave, PSsource;
  int idsource, idsave, dest;

  for(i = Nstart; i < N; i++)
    {
      if(Id[i] != i)
        {
          PSsource = PS[i];

          idsource = Id[i];
          dest     = Id[i];

          do
            {
              PSsave = PS[dest];
              idsave = Id[dest];

              PS[dest] = PSsource;
              Id[dest] = idsource;

              if(dest == i)
                break;

              PSsource = PSsave;
              idsource = idsave;

              dest = idsource;
            }
          while(1);
        }
    }
}

#endif

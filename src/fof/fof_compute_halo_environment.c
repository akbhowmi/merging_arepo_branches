/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/fof/fof.c
 * \date        MM/YYYY
 * \author
 * \brief       parallel FoF group finder
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

#include "../allvars.h"
#include "../domain.h"
#include "../proto.h"
#include "../subfind/subfind.h"
#include "../blackhole/store_mergers_in_snapshot/mergers_io.h"
#include "fof.h"

#ifdef PROBABILISTIC_SEEDING
#include <time.h>
#endif


#ifdef FOF
#ifdef SEED_HALO_ENVIRONMENT_CRITERION

static MyIDType *MinID;
static int *Head, *Len, *Next, *Tail, *MinIDTask;


/*---------------------------------*/
/*! \brief Computes all kinds of properties of groups.
 *
 *  Not complete after calling this. There is still the function
 *  fof_finish_group_properties, which finalizes the calculation
 *  (with normalization, averages, unit conversions and other operations).
 *
 *  \param[in] gr Index in Group array.
 *  \param[in] start Start index in FOF_PList.
 *  \param[in] len Number of particles in this group.
 *
 *  \return void
 */

void fof_compute_halo_environment(int gr, int start, int len)
{
  int j, k, index, type, start_index = FOF_PList[start].Pindex;
  double xyz[3];

#ifdef BLACK_HOLES
#ifdef SEED_HALO_ENVIRONMENT_CRITERION
 double gr_Mass = 0;
 int SelectThisParticle = 0;
#if(PTYPE_USED_FOR_ENVIRONMENT_BASED_SEEDING == 0)
 int gr_IsThisTheDensestCell_sum = 0;
#else
 int gr_IsThisTheMinPotential_sum = 0;  
#endif
 int gr_NumberOfMajorNeighbors  = 0; 
#ifdef CREATE_SUBFOFS
  MyFloat gr_HostHaloMass = 0;
#endif
#endif
#endif
 /* calculate */
  for(k = 0; k < len; k++)
    {
      index = FOF_PList[start + k].Pindex;
#ifdef BLACK_HOLES
#ifdef SEED_HALO_ENVIRONMENT_CRITERION
      gr_Mass += P[index].Mass;
#if(PTYPE_USED_FOR_ENVIRONMENT_BASED_SEEDING == 0)
      if(P[index].Type == 0)
         {
          gr_IsThisTheDensestCell_sum += SphP[index].IsThisTheDensestCell;
          if(SphP[index].IsThisTheDensestCell == 1)
             gr_NumberOfMajorNeighbors = P[index].no_of_BHs_ngb;
         }
#else
      if(P[index].Type == 1)
         {
          gr_IsThisTheMinPotential_sum += P[index].IsThisTheMinPotential;
          if(P[index].IsThisTheMinPotential == 1)
             gr_NumberOfMajorNeighbors = P[index].no_of_BHs_ngb;
         }
#endif

  

#ifdef CREATE_SUBFOFS
#if(PTYPE_USED_FOR_ENVIRONMENT_BASED_SEEDING == 0)
      if (All.SubFOF_mode == 0)
       	    SelectThisParticle = SphP[index].IsThisTheDensestCell;
      else if(All.SubFOF_mode == 1)
            SelectThisParticle = SphP[index].IsThisTheDensestCell_bFOF;
#else
      if (All.SubFOF_mode == 0)
            SelectThisParticle = P[index].IsThisTheMinPotential;
      else if(All.SubFOF_mode == 1)
            SelectThisParticle = P[index].IsThisTheMinPotential_bFOF;
#endif
#if(PTYPE_USED_FOR_ENVIRONMENT_BASED_SEEDING == 0)
      if(P[index].Type == 0) && (SelectThisParticle == 1)
         gr_HostHaloMass += SphP[index].HostHaloMass;
#else
      if((P[index].Type == 1) && (SelectThisParticle == 1))
         gr_HostHaloMass += P[index].HostHaloMass;
#endif
#endif
#endif
#endif
    }
#ifdef BLACK_HOLES
#ifdef SEED_HALO_ENVIRONMENT_CRITERION
#if(PTYPE_USED_FOR_ENVIRONMENT_BASED_SEEDING == 0)
  Group2[gr].IsThisTheDensestCell_sum = gr_IsThisTheDensestCell_sum;
#else
  Group2[gr].IsThisTheMinPotential_sum = gr_IsThisTheMinPotential_sum; 
#endif
  Group2[gr].NumberOfMajorNeighbors = gr_NumberOfMajorNeighbors;
  Group2[gr].Mass = gr_Mass;

#ifdef CREATE_SUBFOFS
  Group2[gr].HostHaloMass = gr_HostHaloMass;
#endif
#endif
#endif
}


/*! \brief Global exchange of identified groups to their appropriate task.
 *
 *  \return void
 */
void fof_exchange_halo_environment_data(void)
{
  struct group_properties *get_Group;
  int i, j, ngrp, recvTask, nimport, start;
  double xyz[3];

  /* sort the groups according to task */
  mysort(Group2, NgroupsExt, sizeof(struct group_properties), fof_compare_Group_MinIDTask);

  /* count how many we have of each task */
  for(i = 0; i < NTask; i++)
    Send_count[i] = 0;
  for(i = 0; i < NgroupsExt; i++)
    Send_count[FOF_GList[i].MinIDTask]++;

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      if(j == ThisTask) /* we will not exchange the ones that are local */
        Recv_count[j] = 0;
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  get_Group = (struct group_properties *)mymalloc("get_Group", sizeof(struct group_properties) * nimport);

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the group data */
              MPI_Sendrecv(&Group2[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct group_properties), MPI_BYTE, recvTask,
                           TAG_ENVIRONMENT, &get_Group[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct group_properties),
                           MPI_BYTE, recvTask, TAG_ENVIRONMENT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  /* sort the groups again according to MinID */
  mysort(Group2, NgroupsExt, sizeof(struct group_properties), fof_compare_Group_MinID);
  mysort(get_Group, nimport, sizeof(struct group_properties), fof_compare_Group_MinID);

  /* now add in the partial imported group data to the main ones */
//  printf("\n Value of Group2[0].MinID %d at task %d",Group2[0].MinID,ThisTask);
  for(i = 0, start = 0; i < nimport; i++)
    {
      while(Group2[start].MinID < get_Group[i].MinID)
        {
          start++;
          if(start >= NgroupsExt)
            terminate("start >= NgroupsExt");
        }
#ifdef BLACK_HOLES
#ifdef SEED_HALO_ENVIRONMENT_CRITERION
#if(PTYPE_USED_FOR_ENVIRONMENT_BASED_SEEDING == 0)
      Group2[start].IsThisTheDensestCell_sum += get_Group[i].IsThisTheDensestCell_sum;
#else
      Group2[start].IsThisTheMinPotential_sum += get_Group[i].IsThisTheMinPotential_sum; 
#endif
      Group2[start].NumberOfMajorNeighbors += get_Group[i].NumberOfMajorNeighbors;
      Group2[start].Mass += get_Group[i].Mass;

#ifdef CREATE_SUBFOFS
      Group2[start].HostHaloMass += get_Group[i].HostHaloMass;
#endif

#endif
#endif
    }

  myfree(get_Group);
}

/*! \brief Finalizes group property calculation.
 *
 *  Called after a loop over all particles of a group is already completed.
 *
 *  \return void
 */
void fof_finish_halo_environment(void)
{
  /* eliminate the non-local groups */
  int ngr = NgroupsExt;
  for(int i = 0; i < ngr; i++)
    {
      if(Group2[i].MinIDTask != (MyIDType)ThisTask)
        {
          Group2[i] = Group2[ngr - 1];
          i--;
          ngr--;
        }
    }

  if(ngr != Ngroups)
    terminate("ngr != Ngroups");

  mysort(Group2, Ngroups, sizeof(struct group_properties), fof_compare_Group_MinID);
}

#endif
#endif /* of FOF */

/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/fof/fof_gfm.c
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

#include "../allvars.h"
#include "../domain.h"
#include "../proto.h"
#include "../subfind/subfind.h"
#include "fof.h"

#ifdef FOF

#ifdef SEED_HALO_ENVIRONMENT_CRITERION
struct group_mass_MinID
{
  MyFloat mass;
  int index_maxdens, task_maxdens;
  int index_MinPot, task_MinPot;
  MyIDType MinID; 
};

int compare_group_mass_ID_densest_gas_cell(const void *a, const void *b)
{
  if(((struct group_mass_MinID *)a)->MinID < (((struct group_mass_MinID *)b)->MinID))
    return -1;

  if(((struct group_mass_MinID *)a)->MinID > (((struct group_mass_MinID *)b)->MinID))
    return +1;

  return 0;
}

void fof_tag_densest_gas_cell(
    void) /* assigns mass of host FoF group to SphP[].w.HostHaloMass for SPH particles and/or to BHP[].HostHaloMass for BH particles */
{
  int i, j, k, start, lenloc, nimport;
  struct group_mass_MinID *required_groups, *groups_to_export;
  for(i = 0; i < NTask; i++)
    Send_count[i] = 0;
  for(i = 0; i < NgroupsExt; i++)         /* loop over all groups for which at least one particle is on this task */
    Send_count[FOF_GList[i].MinIDTask]++; /* its FoF group properties are stored on Task = MinIDTask */

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  qsort(FOF_GList, NgroupsExt, sizeof(struct fof_group_list), fof_compare_FOF_GList_MinIDTask_MinID);

  required_groups = (struct group_mass_MinID *)mymalloc("required_groups", NgroupsExt * sizeof(struct group_mass_MinID));

  for(i = 0; i < NgroupsExt; i++)
    required_groups[i].MinID = FOF_GList[i].MinID;

  groups_to_export = (struct group_mass_MinID *)mymalloc("groups_to_export", nimport * sizeof(struct group_mass_MinID));

  int slen = sizeof(struct group_mass_MinID);
  for(i = 0; i < NTask; i++)
    {
      Send_count[i] *= slen;
      Send_offset[i] *= slen;
      Recv_count[i] *= slen;
      Recv_offset[i] *= slen;
    }

  MPI_Alltoallv(required_groups, Send_count, Send_offset, MPI_BYTE, groups_to_export, Recv_count, Recv_offset, MPI_BYTE,
                MPI_COMM_WORLD);

  for(i = 0; i < NTask; i++)
    {
      Send_count[i] /= slen;
      Send_offset[i] /= slen;
      Recv_count[i] /= slen;
      Recv_offset[i] /= slen;
    }

  for(j = 0, start = 0; j < NTask; j++)
    {
      i = 0;
      k = 0;

      while(i < Recv_count[j] && k < Ngroups)
        {
          if(groups_to_export[start].MinID == Group[k].MinID)
            {
              groups_to_export[start].mass = Group[k].Mass;
              groups_to_export[start].index_maxdens = Group[k].index_maxdens;
              groups_to_export[start].task_maxdens = Group[k].task_maxdens;
              groups_to_export[start].index_MinPot = Group[k].index_MinPot;
              groups_to_export[start].task_MinPot = Group[k].task_MinPot;
              i++;
              k++;
              start++;
            }
          else
            k++;
        }
    }
  if(start != nimport)
    terminate("start != nimport");

  for(i = 0; i < NTask; i++)
    {
      Send_count[i] *= slen;
      Send_offset[i] *= slen;
      Recv_count[i] *= slen;
      Recv_offset[i] *= slen;
    }

  MPI_Alltoallv(groups_to_export, Recv_count, Recv_offset, MPI_BYTE, required_groups, Send_count, Send_offset, MPI_BYTE,
                MPI_COMM_WORLD);

  for(i = 0; i < NTask; i++)
    {
      Send_count[i] /= slen;
      Send_offset[i] /= slen;
      Recv_count[i] /= slen;
      Recv_offset[i] /= slen;
    }

  myfree(groups_to_export);

  qsort(required_groups, NgroupsExt, sizeof(struct group_mass_MinID), compare_group_mass_ID_densest_gas_cell);

#ifdef CREATE_SUBFOFS 
#if(PTYPE_USED_FOR_ENVIRONMENT_BASED_SEEDING == 0)
  if(All.SubFOF_mode == 0)
   {
    for(i = 0; i < NumGas; i++)
     {
        SphP[i].HostHaloMass = 0;
        SphP[i].IsThisTheDensestCell = 0;
     }
   }
  else if(All.SubFOF_mode == 1)
   {
    for(i = 0; i < NumGas; i++)
     {
        SphP[i].HostHaloMass_bFOF = 0;
        SphP[i].IsThisTheDensestCell_bFOF = 0;
     }
   }
#else
  if(All.SubFOF_mode == 0)
   {
    for(i = 0; i < NumPart; i++)
     {
        P[i].HostHaloMass = 0;
        P[i].IsThisTheMinPotential = 0;
     }
   }
  else if(All.SubFOF_mode == 1)
   {
     for(i = 0; i < NumPart; i++)
      {
        P[i].HostHaloMass_bFOF = 0;
        P[i].IsThisTheMinPotential_bFOF = 0;
      }
   }
#endif
#endif

#ifndef CREATE_SUBFOFS
#if(PTYPE_USED_FOR_ENVIRONMENT_BASED_SEEDING == 0)
  for(i = 0; i < NumGas; i++)
    {
      	SphP[i].HostHaloMass = 0;
        SphP[i].IsThisTheDensestCell = 0;
    }
#else
  for(i = 0; i < NumPart; i++)
    {
      	P[i].HostHaloMass = 0;
        P[i].IsThisTheMinPotential = 0;
    }
#endif
#endif

  for(i = 0, start = 0; i < NgroupsExt; i++)
    {
      while(FOF_PList[start].MinID < required_groups[i].MinID)
        {
          start++;
          if(start > NumPart)
            terminate("start > NumPart");
        }

      if(FOF_PList[start].MinID != required_groups[i].MinID)
        terminate("FOF_PList[start].MinID != required_groups[i].MinID");

      for(lenloc = 0; start + lenloc < NumPart;)
        if(FOF_PList[start + lenloc].MinID == required_groups[i].MinID)
          {
#if(PTYPE_USED_FOR_ENVIRONMENT_BASED_SEEDING == 0)
            if(P[FOF_PList[start + lenloc].Pindex].Type == 0)
              {
                if(required_groups[i].mass > All.MinHaloMassForTaggingMinPotential)
                  {
#ifdef CREATE_SUBFOFS
                   if(All.SubFOF_mode == 0)
#endif
                      SphP[FOF_PList[start + lenloc].Pindex].HostHaloMass = required_groups[i].mass;

                   if((FOF_PList[start + lenloc].Pindex == required_groups[i].index_maxdens) && (required_groups[i].task_maxdens == ThisTask))
                     {
#ifdef CREATE_SUBFOFS
       	       	      if(All.SubFOF_mode == 1)  
       	       	       	{
                          SphP[FOF_PList[start + lenloc].Pindex].IsThisTheDensestCell_bFOF = 1;
                          SphP[FOF_PList[start + lenloc].Pindex].HostHaloMass_bFOF = required_groups[i].mass;
       	       	       	}
                      else
                        {
#endif                       
                          SphP[FOF_PList[start + lenloc].Pindex].IsThisTheDensestCell = 1;
#ifdef CREATE_SUBFOFS                         
                        }
#endif
                     }
                   }
              }
#else
            if(P[FOF_PList[start + lenloc].Pindex].Type == 1)
             {
               if(required_groups[i].mass > All.MinHaloMassForTaggingMinPotential)
                 {
#ifdef CREATE_SUBFOFS
                  if(All.SubFOF_mode == 0)
#endif
                      P[FOF_PList[start + lenloc].Pindex].HostHaloMass = required_groups[i].mass;  

                  if((FOF_PList[start + lenloc].Pindex == required_groups[i].index_MinPot) && (required_groups[i].task_MinPot == ThisTask))
                    {
#ifdef CREATE_SUBFOFS
                      if(All.SubFOF_mode == 1)
                        {
                          P[FOF_PList[start + lenloc].Pindex].IsThisTheMinPotential_bFOF = 1;
                          P[FOF_PList[start + lenloc].Pindex].HostHaloMass_bFOF = required_groups[i].mass;
                        }
                      else
                        {
#endif
                          P[FOF_PList[start + lenloc].Pindex].IsThisTheMinPotential = 1;
#ifdef CREATE_SUBFOFS
                        }
#endif
                    }
                 }
             }
#endif
            lenloc++;
          }
        else
          break;

      start += lenloc;
    }

  myfree(required_groups);

  qsort(FOF_GList, NgroupsExt, sizeof(struct fof_group_list), fof_compare_FOF_GList_MinID); /* restore original order */
}
#endif

#endif

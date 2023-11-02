/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/domain_toplevel.c
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
#include <strings.h>

#include "allvars.h"
#include "bsd_tree.h"
#include "domain.h"
#include "proto.h"
#include "voronoi.h"
#include "fof/fof.h"

/* do some preparation work for use of red-black ordered binary tree based on BSD macros */

/* define structure of my tree nodes */
struct mydata_groups
{
  double workload;
  int topnode_index;

  /* this creates the linkage pointers needed by the RB tree, using symbolic name 'linkage' */
  RB_ENTRY(mydata_groups) linkage_groups;
};

/*! \brief Comparison function of tree elements.
 *
 *  Compares elements workload and topnode_index.
 *
 *  \param[in] lhs pointer to left hand side top level tree node.
 *  \param[in] rhs pointer to right hand side top level tree node.
 *
 *  \return -1: left is larger or lower topnode index, 1 opposite, 0 equal.
 */
static int mydata_cmp_groups(struct mydata_groups *lhs, struct mydata_groups *rhs)
{
  if(lhs->workload > rhs->workload)
    return -1;
  else if(lhs->workload < rhs->workload)
    return 1;
  else if(lhs->topnode_index < rhs->topnode_index)
    return -1;
  else if(lhs->topnode_index > rhs->topnode_index)
    return 1;

  return 0;
}

/* the following macro declares 'struct mytree_groups', which is the header element needed as handle for a tree */
RB_HEAD(mytree_groups, mydata_groups);

static struct mydata_groups *nload_groups;
static struct mytree_groups queue_load_groups;

/* the following macros declare appropriate function prototypes and functions needed for this type of tree */
RB_PROTOTYPE_STATIC(mytree_groups, mydata_groups, linkage_groups, mydata_cmp_groups);
RB_GENERATE_STATIC(mytree_groups, mydata_groups, linkage_groups, mydata_cmp_groups);

//static double *list_cost_groups, *list_sphcost_groups;



/*! \brief Free arrays needed in domain decomposition.
 *
 *  This is the counterpart to domain_allocate; need to free arrays in reverse
 *  allocation order.
 *
 * \return void
 */
void domain_free_groups(void)
{
//  if(!DomainStartList)
//    terminate("domain storage not allocated");

  myfree_movable(DomainTask_groups);
  myfree_movable(TopNodes_groups);
  DomainTask_groups             = NULL;
  TopNodes_groups               = NULL;
}



void domain_allocate_groups(void)
{
  MaxTopNodes_groups = (int)(All.TopNodeAllocFactor_groups * MaxNgroups + 1);
  printf("\n MaxNgroups = %d, All.MaxPart = %d, All.TopNodeAllocFactor = %g, All.TopNodeAllocFactor_groups = %g", MaxNgroups,All.MaxPart,All.TopNodeAllocFactor,All.TopNodeAllocFactor_groups);
  TopNodes_groups               = (struct topnode_data *)mymalloc_movable(&TopNodes_groups, "TopNodes_groups", (MaxTopNodes_groups * sizeof(struct topnode_data)));
  DomainTask_groups             = (int *)mymalloc_movable(&DomainTask_groups, "DomainTask_groups", (MaxTopNodes_groups * sizeof(int)));
}

void domain_allocate_lists_groups(void)
{
  DomainLeaveNode_groups = (struct domain_cost_data *)mymalloc_movable(&DomainLeaveNode_groups, "DomainLeaveNode_groups",
                                                                (MaxTopNodes_groups * sizeof(struct domain_cost_data)));
}

void domain_free_lists_groups(void)
{
  myfree(DomainLeaveNode_groups);
}


/*! \brief Construct top-level tree.
 *
 *  This function constructs the global top-level tree node that is used
 *  for the domain decomposition. This is done by considering the string of
 *  Peano-Hilbert keys for all particles, which is recursively chopped off
 *  in pieces of eight segments until each segment holds at most a certain
 *  number of particles.
 *
 *  \return 0
 */
int domain_determineTopTree_groups(void)
{
  double t0 = second();
  int count = 0, message_printed = 0;

  mp_groups           = (struct domain_peano_hilbert_data *)mymalloc_movable(&mp_groups, "mp_groups", sizeof(struct domain_peano_hilbert_data) * Ngroups);
//  list_cost_groups    = (double *)mymalloc_movable(&list_cost_groups, "list_cost_groups", sizeof(double) * Ngroups);
//  list_sphcost_groups = (double *)mymalloc_movable(&list_sphcost_groups, "listsph_cost_groups", sizeof(double) * Ngroups);

  for(int i = 0; i < Ngroups; i++)
    {
      peano1D xb = domain_double_to_int(((Group[i].CM_unwrapped[0] - DomainCorner[0]) * DomainInverseLen) + 1.0);
      peano1D yb = domain_double_to_int(((Group[i].CM_unwrapped[1] - DomainCorner[1]) * DomainInverseLen) + 1.0);
      peano1D zb = domain_double_to_int(((Group[i].CM_unwrapped[2] - DomainCorner[2]) * DomainInverseLen) + 1.0);

      mp_groups[count].key = peano_hilbert_key(xb, yb, zb, BITS_PER_DIMENSION);
      mp_groups[count].index        = i;
      count++;

//      list_cost_groups[i]    = domain_grav_tot_costfactor(i);
//      list_sphcost_groups[i] = domain_hydro_tot_costfactor(i);
    }



  mysort_domain(mp_groups, count, sizeof(struct domain_peano_hilbert_data));

  NTopnodes_groups            = 1;
  NTopleaves_groups           = 1;
  topNodes_groups[0].Daughter = -1;
  topNodes_groups[0].Parent   = -1;
  topNodes_groups[0].Size     = PEANOCELLS;
  topNodes_groups[0].StartKey = 0;
  topNodes_groups[0].PIndex   = 0;
  topNodes_groups[0].Count    = count;
  topNodes_groups[0].Cost     = gravcost;
  topNodes_groups[0].SphCost  = sphcost;

  int limitNTopNodes = 2 * imax(1 + (NTask / 7 + 1) * 8, All.TopNodeFactor * All.MultipleDomains * NTask);

  while(limitNTopNodes > MaxTopNodes_groups)
    {
      mpi_printf("DOMAIN: Increasing TopNodeAllocFactor_groups=%g  ", All.TopNodeAllocFactor_groups);
      All.TopNodeAllocFactor_groups *= 1.3;
      mpi_printf("new value=%g\n", All.TopNodeAllocFactor_groups);
      printf("limitNTopNodes=%d",limitNTopNodes);
      if(All.TopNodeAllocFactor_groups > 1e13)
        terminate("something seems to be going seriously wrong here. Stopping.\n");

      MaxTopNodes_groups = (int)(All.TopNodeAllocFactor_groups * MaxNgroups + 1);

      topNodes_groups        = (struct local_topnode_data *)myrealloc_movable(topNodes_groups, (MaxTopNodes_groups * sizeof(struct local_topnode_data)));
      TopNodes_groups        = (struct topnode_data *)myrealloc_movable(TopNodes_groups, (MaxTopNodes_groups * sizeof(struct topnode_data)));
      DomainTask_groups      = (int *)myrealloc_movable(DomainTask_groups, (MaxTopNodes_groups * sizeof(int)));
      DomainLeaveNode_groups = (struct domain_cost_data *)myrealloc_movable(DomainLeaveNode_groups, (MaxTopNodes_groups * sizeof(struct domain_cost_data)));
    }


  RB_INIT(&queue_load_groups);
  limitNTopNodes = 1;
  nload_groups     = (struct mydata_groups *)mymalloc("nload_groups", limitNTopNodes * sizeof(struct mydata_groups));
  int *list = (int *)mymalloc("list", limitNTopNodes * sizeof(int));
  double limit = 1.0 / (All.TopNodeFactor * All.MultipleDomains * NTask);

  /* insert the root node */

  nload_groups[0].workload      = 1.0;
  nload_groups[0].topnode_index = 0;
  RB_INSERT(mytree_groups, &queue_load_groups, &nload_groups[0]);
  int iter = 0;

  myfree(list);
  myfree(nload_groups);

//  myfree(list_sphcost_groups);
//  myfree(list_cost_groups);
  myfree(mp_groups);

  return 0;

}






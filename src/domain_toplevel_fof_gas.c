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

static double *list_cost_groups, *list_sphcost_groups;
#ifdef AMR
static int *list_level_groups;
#endif


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
  MaxTopNodes_groups = (int)(All.TopNodeAllocFactor_groups * All.MaxPart + 1);
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

  mp_groups           = (struct domain_peano_hilbert_data *)mymalloc_movable(&mp_groups, "mp_groups", sizeof(struct domain_peano_hilbert_data) * NumPart);
  list_cost_groups    = (double *)mymalloc_movable(&list_cost_groups, "list_cost_groups", sizeof(double) * NumPart);
  list_sphcost_groups = (double *)mymalloc_movable(&list_sphcost_groups, "listsph_cost_groups", sizeof(double) * NumPart);
#ifdef AMR
  list_level_groups                = (int *)mymalloc_movable(&list_level_groups, "list_level_groups", sizeof(int) * NumPart);
  int message_printed_level_groups = 0;
#endif

  for(int i = 0; i < NumPart; i++)
    {
      peano1D xb = domain_double_to_int(((P[i].Pos[0] - DomainCorner[0]) * DomainInverseLen) + 1.0);
      peano1D yb = domain_double_to_int(((P[i].Pos[1] - DomainCorner[1]) * DomainInverseLen) + 1.0);
      peano1D zb = domain_double_to_int(((P[i].Pos[2] - DomainCorner[2]) * DomainInverseLen) + 1.0);

      mp_groups[count].key = peano_hilbert_key(xb, yb, zb, BITS_PER_DIMENSION);
      mp_groups[count].index        = i;
      count++;

      list_cost_groups[i]    = domain_grav_tot_costfactor(i);
      list_sphcost_groups[i] = domain_hydro_tot_costfactor(i);

#ifdef AMR
      if(P[i].Type == 0)
        {
          list_level_groups[i] = SphP[i].Level;
        }
      else
	{
          list_level_groups[i] = AMR_MAX_REFLEVEL;
        }
#endif
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
#if defined(AMR) || defined(MODGRAV)
  topNodes_groups[0].Level = 0;
#endif
#ifdef AMR
  topNodes_groups[0].MinCellLevel = All.MinRefLevel;
#endif
  int limitNTopNodes = 2 * imax(1 + (NTask / 7 + 1) * 8, All.TopNodeFactor * All.MultipleDomains * NTask);

#ifdef ADDBACKGROUNDGRID
  limitNTopNodes = imax(limitNTopNodes, 2 * All.GridSize * All.GridSize * All.GridSize);
#endif

  while(limitNTopNodes > MaxTopNodes_groups)
    {
      mpi_printf("DOMAIN: Increasing TopNodeAllocFactor_groups=%g  ", All.TopNodeAllocFactor_groups);
      All.TopNodeAllocFactor_groups *= 1.3;
      mpi_printf("new value=%g\n", All.TopNodeAllocFactor_groups);
      printf("limitNTopNodes=%d",limitNTopNodes);
      if(All.TopNodeAllocFactor_groups > 1e13)
        terminate("something seems to be going seriously wrong here. Stopping.\n");

      MaxTopNodes_groups = (int)(All.TopNodeAllocFactor_groups * All.MaxPart + 1);

      topNodes_groups        = (struct local_topnode_data *)myrealloc_movable(topNodes_groups, (MaxTopNodes_groups * sizeof(struct local_topnode_data)));
      TopNodes_groups        = (struct topnode_data *)myrealloc_movable(TopNodes_groups, (MaxTopNodes_groups * sizeof(struct topnode_data)));
      DomainTask_groups      = (int *)myrealloc_movable(DomainTask_groups, (MaxTopNodes_groups * sizeof(int)));
      DomainLeaveNode_groups = (struct domain_cost_data *)myrealloc_movable(DomainLeaveNode_groups, (MaxTopNodes_groups * sizeof(struct domain_cost_data)));
    }

  RB_INIT(&queue_load_groups);
  nload_groups     = (struct mydata_groups *)mymalloc("nload_groups", limitNTopNodes * sizeof(struct mydata_groups));
  int *list = (int *)mymalloc("list", limitNTopNodes * sizeof(int));
#ifdef ADDBACKGROUNDGRID
  peanokey MaxTopleaveSize = (PEANOCELLS / (All.GridSize * All.GridSize * All.GridSize));
#else
  double limit = 1.0 / (All.TopNodeFactor * All.MultipleDomains * NTask);
#endif

  /* insert the root node */
  nload_groups[0].workload      = 1.0;
  nload_groups[0].topnode_index = 0;
  RB_INSERT(mytree_groups, &queue_load_groups, &nload_groups[0]);
  int iter = 0;

  do
    {
      printf("\n Entering ThisTask = %d ",ThisTask);
      count = 0;

      double first_workload = 0;

      for(struct mydata_groups *nfirst = RB_MIN(mytree_groups, &queue_load_groups); nfirst != NULL; nfirst = RB_NEXT(mytree_groups, &queue_load_groups, nfirst))
        {
          if(topNodes_groups[nfirst->topnode_index].Size >= 8)
            {
              first_workload = nfirst->workload;
              break;
            }
        }

      for(struct mydata_groups *np = RB_MIN(mytree_groups, &queue_load_groups); np != NULL; np = RB_NEXT(mytree_groups, &queue_load_groups, np))
        {
#ifndef ADDBACKGROUNDGRID
          if(np->workload < 0.125 * first_workload)
            break;

          if(NTopnodes_groups + 8 * (count + 1) >= limitNTopNodes)
            break;
#endif

#ifdef ADDBACKGROUNDGRID
          if(topNodes_groups[np->topnode_index].Size > MaxTopleaveSize)
#elif defined(MODGRAV)
          if((np->workload > limit || (NTopleaves_groups < All.MultipleDomains * NTask && count == 0) ||
              topNodes_groups[np->topnode_index].Level >= All.MinLevelTopLeaf) &&
             topNodes_groups[np->topnode_index].Level <= All.MaxLevelTopLeaf)
#else
          if(np->workload > limit || (NTopleaves_groups < All.MultipleDomains * NTask && count == 0))
#endif
            {
              if(topNodes_groups[np->topnode_index].Size < 8)
                {
                  if(!message_printed)
                    {
                      mpi_printf("DOMAIN: Note: we would like to refine top-level tree, but Peano grid is not fine enough\n");

#ifdef ADDBACKGROUNDGRID
                      mpi_printf("DOMAIN: workload %g, Topleaves %d, Domains %d, count %d\n", np->workload, NTopleaves_groups,
                                 All.MultipleDomains * NTask, count);
#else
                      mpi_printf("DOMAIN: limit %g, workload %g, Topleaves %d, Domains %d, count %d\n", limit, np->workload,
                                 NTopleaves_groups, All.MultipleDomains * NTask, count);
#endif

#ifndef OVERRIDE_PEANOGRID_WARNING
                      terminate(
                          "Consider setting BITS_PER_DIMENSION (current value=%d) up to a value of 42 to get a fine enough Peano "
                          "grid, or force a continuation by activating OVERRIDE_PEANOGRID_WARNING",
                          BITS_PER_DIMENSION);
#endif

                      message_printed = 1;
                    }
                }
#ifdef AMR
              else if(topNodes_groups[np->topnode_index].Level + 1 >= topNodes_groups[np->topnode_index].MinCellLevel)
                {
                  if(message_printed_level_groups == 0)
                    {
                      mpi_printf("DOMAIN: Note: we would like to refine top-tree, but cells are too course\n");
                      message_printed_level_groups = 1;
                    }
                }
#endif
              else
                {
                  list[count] = np->topnode_index;
                  count++;
                }
            }
        }

      if(count > 0)
        {
          printf("\n Count = %d,trying to refine ThisTask = %d",count,ThisTask);
//          domain_do_local_refine_groups(count, list);
          iter++;
        }
      else
        printf("\n Count = %d ThisTask = %d, exiting",count,ThisTask);
    }
  while(count > 0);

  myfree(list);
  myfree(nload_groups);
#ifdef AMR
  myfree(list_level_groups);
#endif
  myfree(list_sphcost_groups);
  myfree(list_cost_groups);
  myfree(mp_groups);

  /* count the number of top leaves */
  NTopleaves_groups = 0;
  domain_walktoptree_groups(0);
#ifndef AMR
//  if(NTopleaves_groups < All.MultipleDomains * NTask)
//    terminate("NTopleaves_groups = %d < All.MultipleDomains * NTask = %d * %d = %d", NTopleaves, All.MultipleDomains, NTask,
//              All.MultipleDomains * NTask);
#endif

  double t1 = second();
  mpi_printf("DOMAIN: NTopleaves_groups=%d, determination of top-level fof tree involved %d iterations and took %g sec\n", NTopleaves_groups, iter,
             timediff(t0, t1));

  t0 = second();

//  domain_sumCost();

  t1 = second();
  mpi_printf("DOMAIN group: cost summation for top-level tree took %g sec\n", timediff(t0, t1));

  return 0;
}

/*! \brief Refine top-level tree locally.
 *
 *  \param[in] n Number of nodes that should be refined.
 *  \param[in] list List of node indices that should be refined.
 *
 *  \return void
 */
void domain_do_local_refine_groups(int n, int *list) /* In list[], we store the node indices hat should be refined, N is their number */
{
  double *worktotlist = (double *)mymalloc("worktotlist", 8 * n * sizeof(double));
  double *worklist    = (double *)mymalloc("worklist", 8 * n * sizeof(double));

#ifdef AMR
  int *levellist_min = (int *)mymalloc("levellist_min", 8 * n * sizeof(int));
  int *levellist     = (int *)mymalloc("levellist", 8 * n * sizeof(int));
#endif

  double non_zero = 0, non_zero_tot;

  /* create the new nodes */
  for(int k = 0; k < n; k++)
    {
      int i                = list[k];
      topNodes_groups[i].Daughter = NTopnodes_groups;
      NTopnodes_groups += 8;
      NTopleaves_groups += 7;

      for(int j = 0; j < 8; j++)
        {
          int sub = topNodes_groups[i].Daughter + j;

          topNodes_groups[sub].Daughter = -1;
          topNodes_groups[sub].Parent   = i;
          topNodes_groups[sub].Size     = (topNodes_groups[i].Size >> 3);
          topNodes_groups[sub].StartKey = topNodes_groups[i].StartKey + j * topNodes_groups[sub].Size;
          topNodes_groups[sub].PIndex   = topNodes_groups[i].PIndex;
          topNodes_groups[sub].Count    = 0;
          topNodes_groups[sub].Cost     = 0;
          topNodes_groups[sub].SphCost  = 0;
#ifdef AMR
          topNodes_groups[sub].Level        = topNodes_groups[i].Level + 1;
          topNodes_groups[sub].MinCellLevel = AMR_MAX_REFLEVEL;
#endif
        }

      int sub = topNodes_groups[i].Daughter;

      for(int p = topNodes_groups[i].PIndex, j = 0; p < topNodes_groups[i].PIndex + topNodes_groups[i].Count; p++)
        {

          printf("p = %d %d %lli, ThisTask = %d \n",p,topNodes_groups[i].PIndex, topNodes_groups[i].Count,ThisTask);
 //         printf("Keys: %g,%g",mp_groups[p].key,topNodes_groups[sub + 1].StartKey);
          if(j < 7)
            while(mp_groups[p].key >= topNodes_groups[sub + 1].StartKey)
              {
//                printf("Keys: %d,%d",mp_groups[p].key,topNodes_groups[sub + 1].StartKey);
                j++;
                sub++;
                topNodes_groups[sub].PIndex = p;
                printf("Entered while loop p = %d %d %lli \n",p,topNodes_groups[i].PIndex, topNodes_groups[i].Count);
                if(j >= 7)
                  break;
              }

          printf("Line 2 j= %d p = %d %d %lli, ThisTask = %d \n",j,p,topNodes_groups[i].PIndex, topNodes_groups[i].Count,ThisTask);
          topNodes_groups[sub].Cost += list_cost_groups[mp_groups[p].index];
          printf("Line 3 j= %d p = %d %d %lli, ThisTask = %d \n",j,p,topNodes_groups[i].PIndex, topNodes_groups[i].Count,ThisTask);
          topNodes_groups[sub].SphCost += list_sphcost_groups[mp_groups[p].index];
          printf("Line 4 j= %d p = %d %d %lli, ThisTask = %d \n",j,p,topNodes_groups[i].PIndex, topNodes_groups[i].Count,ThisTask);
#ifdef AMR
          topNodes_groups[sub].MinCellLevel = imin(topNodes_groups[sub].MinCellLevel, list_level_groups[mp_groups[p].index]);
          printf("Line 5 j= %d p = %d %d %lli, ThisTask = %d \n",j,p,topNodes_groups[i].PIndex, topNodes_groups[i].Count,ThisTask);
#endif
          topNodes_groups[sub].Count++;

          printf("Last line j= %d p = %d %d %lli, ThisTask = %d \n",j,p,topNodes_groups[i].PIndex, topNodes_groups[i].Count,ThisTask);
        }

//      printf("Exit loop p = %d %d \n ",p,topNodes_groups[i].PIndex + topNodes_groups[i].Count);


      for(int j = 0; j < 8; j++)
        {
          int sub_ = topNodes_groups[i].Daughter + j;
          worklist[k * 8 + j] =
              fac_work * topNodes_groups[sub_].Cost + fac_worksph * topNodes_groups[sub_].SphCost + fac_load * topNodes_groups[sub_].Count;

          if(worklist[k * 8 + j] != 0)
            non_zero++;

#ifdef AMR
          levellist[k * 8 + j] = topNodes_groups[sub_].MinCellLevel;
#endif
        }
    }

  MPI_Allreduce(&non_zero, &non_zero_tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if(non_zero_tot > 0.05 * (NTask * 8 * n))
    MPI_Allreduce(worklist, worktotlist, 8 * n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  else
    allreduce_sparse_double_sum(worklist, worktotlist, 8 * n);

#ifdef AMR
  MPI_Allreduce(levellist, levellist_min, 8 * n, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

  for(int k = 0; k < n; k++)
    {
      int i = list[k];
      for(int j = 0; j < 8; j++)
        {
          int sub                    = topNodes_groups[i].Daughter + j;
          topNodes_groups[sub].MinCellLevel = levellist_min[k * 8 + j];
        }
    }

  myfree(levellist);
  myfree(levellist_min);
#endif

  for(int k = 0; k < n; k++)
    {
      int i = list[k];
      if(nload_groups[i].topnode_index==1)
             printf("just");
      RB_REMOVE(mytree_groups, &queue_load_groups, &nload_groups[i]);
    }

  for(int k = 0, l = 0; k < n; k++)
    {
      int i = list[k];

      for(int j = 0; j < 8; j++, l++)
        {
          int sub = topNodes_groups[i].Daughter + j;

          /* insert the node */
          nload_groups[sub].workload      = worktotlist[l];
          nload_groups[sub].topnode_index = sub;
          RB_INSERT(mytree_groups, &queue_load_groups, &nload_groups[sub]);
        }
    }

  myfree(worklist);
  myfree(worktotlist);
}

/*! \brief Walks top level tree recursively.
 *
 *  This function walks the global top tree in order to establish the
 *  number of leaves it has, and for assigning the leaf numbers along the
 *  Peano-Hilbert Curve. These leaves are later combined to domain pieces,
 *  which are distributed to different processors.
 *
 *  \param[in] no Present node.
 *
 *  \return void
 */
void domain_walktoptree_groups(int no)
{
  if(topNodes_groups[no].Daughter == -1)
    {
      topNodes_groups[no].Leaf = NTopleaves_groups;
      NTopleaves_groups++;
    }
  else
    {
      for(int i = 0; i < 8; i++)
        domain_walktoptree_groups(topNodes_groups[no].Daughter + i);
    }
}

/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/ngbtree.c
 * \date        MM/YYYY
 * \author
 * \brief       construct neighbor tree
 * \details     This file contains the neighbor tree construction
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
#include <time.h>

#include "allvars.h"
#include "domain.h"
#include "forcetree.h"
#include "proto.h"
#include "fof/fof.h"

static void ngb_record_topnode_siblings_groups(int no, int sib);
static int ngb_treebuild_construct_groups(int ngroups);
static void ngb_update_node_recursive_groups(int no, int sib, int father, int *last, int mode);
static void ngb_exchange_topleafdata_groups(void);
static int ngb_create_empty_nodes_groups(int no, int topnode, int bits, int x, int y, int z);
static void ngb_update_vbounds_groups(int i, int *nchanged, int *nodelist);
static void ngb_finish_vounds_update_groups(int nchanged, int *nodelist);

static int *Ngb_Node_Tmp_Sibling_groups;

/*! \brief This function is a driver routine for constructing the neighbor
 *         oct-tree, which is done by calling a small number of other
 *         functions.
 *
 *  Does not build a tree if TotNumGroups == 0.
 *
 *  \param[in] ngroups Number of groups in tree.
 *
 *  \return Number of nodes in the tree.
 */
int ngb_treebuild_groups(int ngroups)
{

  TIMER_START(CPU_NGBTREEBUILD);

  mpi_printf("NGBTREE: Ngb-tree construction.  (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));

  double t0 = second();

  int flag;
  do
    {
      int flag_single = ngb_treebuild_construct_groups(ngroups);

      MPI_Allreduce(&flag_single, &flag, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      if(flag == -1)
        {
          Ngb_Node_Tmp_Sibling_groups += Ngb_MaxPart_groups;
          myfree(Ngb_Node_Tmp_Sibling_groups);
          ngb_treefree_groups();

          All.NgbTreeAllocFactor *= 1.15;
          mpi_printf("Increasing NgbTreeAllocFactor, new value=%g\n", All.NgbTreeAllocFactor);

          ngb_treeallocate_groups();
        }
    }
  while(flag == -1);

  int ntopleaves = DomainNLocalTopleave_groups[ThisTask];
  int *list      = DomainListOfLocalTopleaves_groups + DomainFirstLocTopleave_groups[ThisTask];


#pragma omp parallel for schedule(dynamic)
  for(int i = 0; i < ntopleaves; i++)
    {
      int last = -1;
      int no   = Ngb_DomainNodeIndex_groups[list[i]];

      if(no < Ngb_MaxPart_groups || no >= Ngb_MaxPart_groups + Ngb_MaxNodes_groups)
        terminate("i=%d no=%d  task=%d \n", i, no, DomainTask_groups[list[i]]);

      ngb_update_node_recursive_groups(no, Ngb_Node_Tmp_Sibling_groups[no], no, &last, 0);

      /* if there was no particle in the node, we need to initialize nextnode of the node */
      if(no == last)
        Ngb_Nodes_groups[no].u.d.nextnode = -1;

      Ngb_Nodes_groups[no].u.d.sibling = last; /* we temporarily store this here and will later restore this sibling pointer,
                                           which is anyway equal to Ngb_Node_Tmp_Sibling_groups[index] */
    }

  ngb_exchange_topleafdata_groups();

  /* now put in "pseudo" particles as nextnode in non-local topleaves */
  for(int i = 0; i < NTopleaves_groups; i++)
    {
      if(DomainTask_groups[i] != ThisTask)
        {
          int index                     = Ngb_DomainNodeIndex_groups[i];
          Ngb_Nodes_groups[index].u.d.nextnode = Ngb_MaxPart_groups + Ngb_MaxNodes_groups + i;
        }
    }

  /* now update the top-level tree nodes */
  int last = -1;
  ngb_update_node_recursive_groups(Ngb_MaxPart_groups, -1, -1, &last, 1);

  if(last >= Ngb_MaxPart_groups)
    {
      if(last >= Ngb_MaxPart_groups + Ngb_MaxNodes_groups) /* a pseudo-particle */
        Ngb_Nextnode_groups[last - Ngb_MaxNodes_groups] = -1;
      else
        Ngb_Nodes_groups[last].u.d.nextnode = -1;
    }
  else
    Ngb_Nextnode_groups[last] = -1;

  TIMER_STOPSTART(CPU_NGBTREEBUILD, CPU_LOGS);

  double numnodes = Ngb_NumNodes_groups, tot_numnodes;
  MPI_Reduce(&numnodes, &tot_numnodes, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  double t1 = second();
  mpi_printf("NGBTREE: Ngb-tree construction done for groups. took %g sec  <numnodes>=%g  NTopnodes_groups=%d NTopleaves_groups=%d\n", timediff(t0, t1),
             tot_numnodes / NTask, NTopnodes_groups, NTopleaves_groups);

  Ngb_Node_Tmp_Sibling_groups += Ngb_MaxPart_groups;
  myfree(Ngb_Node_Tmp_Sibling_groups);

  Ngb_MarkerValue_groups = 0;
  memset(Ngb_Marker_groups, 0, (Ngb_MaxPart_groups + Ngb_NumNodes_groups) * sizeof(int));

  TIMER_STOP(CPU_LOGS);

  return Ngb_NumNodes_groups;
}

/*! \brief Converts double precision coordinate to unsigned long long int.
 *
 *  \param[in] d Double precision coordinate that is to be converted.
 *
 *  \return Unsigned long long int represenation of d.
 */
static inline unsigned long long ngb_double_to_int(double d)
{
  union
  {
    double d;
    unsigned long long ull;
  } u;
  u.d = d;
  return (u.ull & 0xFFFFFFFFFFFFFllu);
}

/*! Constructs the neighbor oct-tree.
 *
 *  The index convention for accessing tree nodes is the following:
 *
 *  0...NumPart-1                             reference single particles
 *  Ngb_MaxPart.... Ngb_MaxPart+Numnodes-1    references tree nodes.
 *  Ngb_MaxPart + All.MaxNgb_Nodes....                reference "pseudo particles", i.e. the marker that indicates a top-node lying on
 * another CPU.
 *
 *  `Ngb_Nodes_base' points to the first tree node,
 *  `Ngb_Nodes' is shifted such that Ngb_Nodes_groups[Ngb_MaxPart] gives the first tree node.
 *
 */
int ngb_treebuild_construct_groups(int ngroups)
{
  /* create an empty root node  */
  Ngb_NextFreeNode_groups = Ngb_MaxPart_groups; /* index of first free node */

  for(int i = 0; i < 8; i++)
    Ngb_Nodes_groups[Ngb_NextFreeNode_groups].u.suns[i] = -1;

  Ngb_DomainNodeIndex_groups[0] = Ngb_NextFreeNode_groups;


  Ngb_NumNodes_groups = 1;
  Ngb_NextFreeNode_groups++;

  /* create a set of empty nodes corresponding to the top-level domain
   * grid. We need to generate these nodes first to make sure that we have a
   * complete top-level tree which allows the easy insertion of the
   * pseudo-particles at the right place
   */
  if(ngb_create_empty_nodes_groups(Ngb_MaxPart_groups, 0, 1, 0, 0, 0) < 0)
    return -1;

  Ngb_FirstNonTopLevelNode_groups = Ngb_NextFreeNode_groups;

  Ngb_Node_Tmp_Sibling_groups = (int *)mymalloc("Ngb_Node_Tmp_Sibling_groups", (Ngb_MaxNodes_groups + 1) * sizeof(int));
  Ngb_Node_Tmp_Sibling_groups -= Ngb_MaxPart_groups;

  ngb_record_topnode_siblings_groups(Ngb_MaxPart_groups, -1);

#if NUM_THREADS > 1
  omp_lock_t *node_lock = (omp_lock_t *)mymalloc("node_lock", NTopleaves_groups * sizeof(omp_lock_t));

  for(int i = 0; i < NTopleaves_groups; i++)
    omp_init_lock(&node_lock[i]);
#endif

  unsigned long long *ngbTree_IntPos_list =
      (unsigned long long *)mymalloc("ngbTree_IntPos_list", 3 * ngroups * sizeof(unsigned long long));

  /* now we insert all particles */
#pragma omp parallel
  {
    int out_of_space = 0;

#if NUM_THREADS > 1
    int last_used_no = -1;
#endif

    int threadid = get_thread_num();
    int start, end, size;

    int first_empty_slot = Ngb_NextFreeNode_groups + threadid * TAKE_NSLOTS_IN_ONE_GO;
    int count_empty_slot = TAKE_NSLOTS_IN_ONE_GO;

#pragma omp barrier
    if(threadid == 0)
      Ngb_NextFreeNode_groups += NUM_THREADS * TAKE_NSLOTS_IN_ONE_GO;
#pragma omp barrier
#pragma omp flush(Ngb_NextFreeNode_groups)

    size  = (ngroups - 1) / NUM_THREADS + 1;
    start = threadid * size;
    end   = (threadid + 1) * size - 1;
    if(end >= ngroups)
      end = ngroups - 1;

    for(int i = start; i <= end && out_of_space == 0; i++)
      {

        unsigned long long xxb  = ngb_double_to_int(((Group[i].Pos[0] - DomainCorner[0]) * DomainInverseLen) + 1.0);
        unsigned long long yyb  = ngb_double_to_int(((Group[i].Pos[1] - DomainCorner[1]) * DomainInverseLen) + 1.0);
        unsigned long long zzb  = ngb_double_to_int(((Group[i].Pos[2] - DomainCorner[2]) * DomainInverseLen) + 1.0);
        unsigned long long mask = ((unsigned long long)1) << (52 - 1);
        unsigned char shiftx    = (52 - 1);
        unsigned char shifty    = (52 - 2);
        unsigned char shiftz    = (52 - 3);
        unsigned char levels    = 0;

        ngbTree_IntPos_list[3 * i + 0] = xxb;
        ngbTree_IntPos_list[3 * i + 1] = yyb;
        ngbTree_IntPos_list[3 * i + 2] = zzb;

        int no = 0;
        while(TopNodes_groups[no].Daughter >= 0) /* walk down top tree to find correct leaf */
          {
            unsigned char subnode = (((unsigned char)((xxb & mask) >> (shiftx--))) | ((unsigned char)((yyb & mask) >> (shifty--))) |
                                     ((unsigned char)((zzb & mask) >> (shiftz--))));

            mask >>= 1;
            levels++;

            no = TopNodes_groups[no].Daughter + TopNodes_groups[no].MortonToPeanoSubnode[subnode];
          }

        no = TopNodes_groups[no].Leaf;

        if(DomainTask_groups[no] != ThisTask)
          terminate("STOP!  Something is inserted into task=%d, but should be on task=%d no=%d\n", ThisTask, DomainTask_groups[no], no);

        int th = Ngb_DomainNodeIndex_groups[no];

        signed long long centermask = (0xFFF0000000000000llu) >> levels;

        int parent            = -1; /* note: will not be used below before it is changed */
        unsigned char subnode = 0;

#if NUM_THREADS > 1
        if(last_used_no != no)
          {
            if(last_used_no >= 0)
              {
                omp_unset_lock(&node_lock[last_used_no]);
                last_used_no = -1;
              }
            omp_set_lock(&node_lock[no]);
            last_used_no = no;
          }
#endif

        while(1)
          {
            if(th >= Ngb_MaxPart_groups) /* we are dealing with an internal node */
              {
                subnode = (((unsigned char)((xxb & mask) >> (shiftx--))) | ((unsigned char)((yyb & mask) >> (shifty--))) |
                           ((unsigned char)((zzb & mask) >> (shiftz--))));

                centermask >>= 1;
                mask >>= 1;
                levels++;

                if(levels > MAX_TREE_LEVEL)
                  {
                    /* seems like we're dealing with particles at identical (or extremely close)
                     * locations. Shift subnode index to allow tree construction. Note: Multipole moments
                     * of tree are still correct, but one should MAX_TREE_LEVEL large enough to have
                     *      DomainLen/2^MAX_TREE_LEEL  < gravitational softening length
                     */
                    for(int j = 0; j < 8; j++)
                      {
                        if(Ngb_Nodes_groups[th].u.suns[subnode] < 0)
                          break;

                        subnode++;
                        if(subnode >= 8)
                          subnode = 7;
                      }
                  }

                int nn = Ngb_Nodes_groups[th].u.suns[subnode];

                if(nn >= 0) /* ok, something is in the daughter slot already, need to continue */
                  {
                    parent = th;
                    th     = nn;
                  }
                else
                  {
                    /* here we have found an empty slot where we can attach
                     * the new particle as a leaf.
                     */
                    Ngb_Nodes_groups[th].u.suns[subnode] = i;
                    break; /* done for this particle */
                  }
              }
            else
              {
                /* We try to insert into a leaf with a single particle.  Need
                 * to generate a new internal node at this point.
                 * Then resume trying to insert the new particle at
                 * the newly created internal node
                 */
                int thold = th;

                if(count_empty_slot)
                  {
                    th = first_empty_slot + (TAKE_NSLOTS_IN_ONE_GO - count_empty_slot);
                    count_empty_slot--;
                  }
                else
                  {
#pragma omp atomic capture
                    {
                      th = Ngb_NextFreeNode_groups;
                      Ngb_NextFreeNode_groups += TAKE_NSLOTS_IN_ONE_GO;
                    }

                    first_empty_slot = th;
                    count_empty_slot = (TAKE_NSLOTS_IN_ONE_GO - 1);

                    if(first_empty_slot + TAKE_NSLOTS_IN_ONE_GO - Ngb_MaxPart_groups >= Ngb_MaxNodes_groups)
                      {
                        out_of_space = 1;
                        break;
                      }
                  }

                Ngb_Nodes_groups[parent].u.suns[subnode] = th;
                struct NgbNODE *nfreep            = &Ngb_Nodes_groups[th];

                for(int j = 0; j < 8; j++)
                  nfreep->u.suns[j] = -1;

                unsigned long long *intppos = &ngbTree_IntPos_list[3 * thold];

                subnode = (((unsigned char)((intppos[0] & mask) >> shiftx)) | ((unsigned char)((intppos[1] & mask) >> shifty)) |
                           ((unsigned char)((intppos[2] & mask) >> shiftz)));

                nfreep->u.suns[subnode] = thold;
              }
          }
      }

#if NUM_THREADS > 1
    if(last_used_no >= 0)
      {
        omp_unset_lock(&node_lock[last_used_no]);
        last_used_no = -1;
      }
#endif
  }

  myfree(ngbTree_IntPos_list);

#if NUM_THREADS > 1
  myfree(node_lock);
#endif

  if((Ngb_NumNodes_groups = Ngb_NextFreeNode_groups - Ngb_MaxPart_groups) >= Ngb_MaxNodes_groups)
    {
      if(All.NgbTreeAllocFactor > MAX_TREE_ALLOC_FACTOR)
        {
          dump_particles();
          terminate("task %d: out of space for neighbor tree, stopping with particle dump.\n", ThisTask);
        }
      else
        return -1;
    }

  return 0;
}

/*! This function recursively creates a set of empty tree nodes which
 *  corresponds to the top-level tree for the domain grid. This is done to
 *  ensure that this top-level tree is always "complete" so that we can easily
 *  associate the pseudo-particles of other CPUs with tree-nodes at a given
 *  level in the tree, even when the particle population is so sparse that
 *  some of these nodes are actually empty.
 */
int ngb_create_empty_nodes_groups(int no, int topnode, int bits, int x, int y, int z)
{
  if(TopNodes_groups[topnode].Daughter >= 0)
    {
      for(int i = 0; i < 2; i++)
        for(int j = 0; j < 2; j++)
          for(int k = 0; k < 2; k++)
            {
              if(Ngb_NumNodes_groups >= Ngb_MaxNodes_groups)
                {
                  if(All.NgbTreeAllocFactor > MAX_TREE_ALLOC_FACTOR)
                    {
                      dump_particles();
                      terminate("task %d: looks like a serious problem (NTopnodes_groups=%d), stopping with particle dump.\n", ThisTask,
                                NTopnodes_groups);
                    }
                  return -1;
                }

              int sub = 7 & peano_hilbert_key((x << 1) + i, (y << 1) + j, (z << 1) + k, bits);

              int count = i + 2 * j + 4 * k;

              Ngb_Nodes_groups[no].u.suns[count] = Ngb_NextFreeNode_groups;

              for(int n = 0; n < 8; n++)
                Ngb_Nodes_groups[Ngb_NextFreeNode_groups].u.suns[n] = -1;

              if(TopNodes_groups[TopNodes_groups[topnode].Daughter + sub].Daughter == -1)
                Ngb_DomainNodeIndex_groups[TopNodes_groups[TopNodes_groups[topnode].Daughter + sub].Leaf] = Ngb_NextFreeNode_groups;

              Ngb_NextFreeNode_groups++;
              Ngb_NumNodes_groups++;

              if(ngb_create_empty_nodes_groups(Ngb_NextFreeNode_groups - 1, TopNodes_groups[topnode].Daughter + sub, bits + 1, 2 * x + i, 2 * y + j,
                                        2 * z + k) < 0)
                return -1;
            }
    }

  return 0;
}

/*! this routine determines the node ranges a given internal node
 *  and all its subnodes using a recursive computation.  The result is
 *  stored in the Ngb_Nodes[] structure in the sequence of this tree-walk.
 *  mode = 0: process a leave branch, mode = 1: process top-level nodes
 */
void ngb_update_node_recursive_groups(int no, int sib, int father, int *last, int mode)
{
  int j, jj, k, p, pp, nextsib, suns[8];
  MyNgbTreeFloat range_min[3];
  MyNgbTreeFloat range_max[3];
  MyNgbTreeFloat vertex_vmin[3];
  MyNgbTreeFloat vertex_vmax[3];
  if(no >= Ngb_MaxPart_groups && no < Ngb_MaxPart_groups + Ngb_MaxNodes_groups) /* internal node */
    {
      if(*last >= 0)
        {
          if(*last >= Ngb_MaxPart_groups)
            {
              if(*last == no)
                terminate("as");

              if(*last >= Ngb_MaxPart_groups + Ngb_MaxNodes_groups) /* a pseudo-particle */
                Ngb_Nextnode_groups[*last - Ngb_MaxNodes_groups] = no;
              else
                Ngb_Nodes_groups[*last].u.d.nextnode = no;
            }
          else
            Ngb_Nextnode_groups[*last] = no;
        }

      *last = no;

      int not_interal_top_level = 0;

      if(mode == 1)
        {
          if(!(no >= Ngb_MaxPart_groups && no < Ngb_FirstNonTopLevelNode_groups))
            terminate("can't be");

          if(Ngb_Node_Tmp_Sibling_groups[no] != -2)
            not_interal_top_level = 1;
        }

      if(not_interal_top_level)
        {
          p = Ngb_Nodes_groups[no].u.d.nextnode;

          if(p >= Ngb_MaxPart_groups + Ngb_MaxNodes_groups &&
             p < Ngb_MaxPart_groups + Ngb_MaxNodes_groups + NTopleaves_groups) /* a pseudo-particle, i.e. we are dealing with a non-local top-leave */
            ngb_update_node_recursive_groups(p, sib, no, last, mode);
          else
            {
              /* this is local toplevel node */
              *last = Ngb_Nodes_groups[no].u.d.sibling;
            }

          if(Ngb_Node_Tmp_Sibling_groups[no] != sib)
            terminate("Ngb_Node_Tmp_Sibling_groups[no] != sib");

          /* restore the sibling pointer for local toplevel nodes (we had temporarily stored the last element in this branch */
          Ngb_Nodes_groups[no].u.d.sibling = sib;
          Ngb_Nodes_groups[no].father      = father;
        }
      else
        {
          for(j = 0; j < 8; j++)
            suns[j] = Ngb_Nodes_groups[no].u.suns[j]; /* this "backup" is necessary because the nextnode entry will
                                                  overwrite one element (union!) */

         for(k = 0; k < 3; k++)
            {
              range_min[k] = MAX_NGBRANGE_NUMBER;
              range_max[k] = -MAX_NGBRANGE_NUMBER;

              vertex_vmin[k] = MAX_NGBRANGE_NUMBER;
              vertex_vmax[k] = -MAX_NGBRANGE_NUMBER;
           }

          for(j = 0; j < 8; j++)
            {
              if((p = suns[j]) >= 0)
                {
                  /* check if we have a sibling on the same level */
                  for(jj = j + 1; jj < 8; jj++)
                    if((pp = suns[jj]) >= 0)
                      break;

                  if(jj < 8) /* yes, we do */
                    nextsib = pp;
                  else
                    nextsib = sib;


                  ngb_update_node_recursive_groups(p, nextsib, no, last, mode);

                 if(p >= Ngb_MaxPart_groups) /* an internal node or pseudo particle */
                    {
                      if(p >= Ngb_MaxPart_groups + Ngb_MaxNodes_groups) /* a pseudo particle */
                        {
                          /* nothing to be done here because the mass of the
                           * pseudo-particle is still zero. This will be changed
                           * later.
                           */
                        }
                      else
                        {
                         for(k = 0; k < 3; k++)
                            {
                              if(range_min[k] > Ngb_Nodes_groups[p].u.d.range_min[k])
                                range_min[k] = Ngb_Nodes_groups[p].u.d.range_min[k];

                              if(range_max[k] < Ngb_Nodes_groups[p].u.d.range_max[k])
                                range_max[k] = Ngb_Nodes_groups[p].u.d.range_max[k];

                              if(vertex_vmin[k] > Ngb_Nodes_groups[p].vertex_vmin[k])
                                vertex_vmin[k] = Ngb_Nodes_groups[p].vertex_vmin[k];

                              if(vertex_vmax[k] < Ngb_Nodes_groups[p].vertex_vmax[k])
                                vertex_vmax[k] = Ngb_Nodes_groups[p].vertex_vmax[k];

                            }
                        }
                    }
                  else /* a particle */
                    {
                      for(k = 0; k < 3; k++)
                        {
                          if(range_min[k] > Group[p].Pos[k])
                            range_min[k] = Group[p].Pos[k];

                          if(range_max[k] < Group[p].Pos[k])
                            range_max[k] = Group[p].Pos[k];

                          if(vertex_vmin[k] > Group[p].Vel[k])
                             vertex_vmin[k] = Group[p].Vel[k];

                          if(vertex_vmax[k] < Group[p].Vel[k])
                             vertex_vmax[k] = Group[p].Vel[k];
                        }
                    }
                }
           }

          for(k = 0; k < 3; k++)
            {
              Ngb_Nodes_groups[no].u.d.range_min[k] = range_min[k];
              Ngb_Nodes_groups[no].u.d.range_max[k] = range_max[k];
              Ngb_Nodes_groups[no].vertex_vmin[k]   = vertex_vmin[k];
              Ngb_Nodes_groups[no].vertex_vmax[k]   = vertex_vmax[k];
            }

          Ngb_Nodes_groups[no].u.d.sibling = sib;
          Ngb_Nodes_groups[no].father      = father;

          Ngb_Nodes_groups[no].Ti_Current = All.Ti_Current;
        }
    }
  else /* single particle or pseudo particle */
    {
      if(*last >= 0)
        {
          if(*last >= Ngb_MaxPart_groups)
            {
              if(*last >= Ngb_MaxPart_groups + Ngb_MaxNodes_groups) /* a pseudo-particle */
                Ngb_Nextnode_groups[*last - Ngb_MaxNodes_groups] = no;
              else
                Ngb_Nodes_groups[*last].u.d.nextnode = no;
            }
          else
            {
              Ngb_Nextnode_groups[*last] = no;
            }
        }
      if(no < Ngb_MaxPart_groups) /* only set it for single particles... */
        {
          if(father < Ngb_MaxPart_groups)
            terminate("no=%d father=%d\n", no, father);

          Ngb_Father_groups[no] = father;
        }

      *last = no;
    }
}

void ngb_record_topnode_siblings_groups(int no, int sib)
{
  /* note: when this routine is called, only toplevel tree nodes are present */

  if(Ngb_Nodes_groups[no].u.suns[0] >= 0)
    {
      /* marker value to designate internal nodes in the top-level tree */
      Ngb_Node_Tmp_Sibling_groups[no] = -2;

      if(Ngb_Nodes_groups[no].u.suns[0] >= 0)
        for(int j = 0; j < 8; j++)
          {
            int p = Ngb_Nodes_groups[no].u.suns[j];
            int nextsib;

            if(j < 7)
              nextsib = Ngb_Nodes_groups[no].u.suns[j + 1];
            else
              nextsib = sib;

            ngb_record_topnode_siblings_groups(p, nextsib);
          }
    }
  else
    Ngb_Node_Tmp_Sibling_groups[no] = sib; /* a top-level leave node */
}

void ngb_exchange_topleafdata_groups(void)
{
  struct DomainNODE
  {
    MyNgbTreeFloat range_min[3];
    MyNgbTreeFloat range_max[3];
    MyNgbTreeFloat vertex_vmin[3];
    MyNgbTreeFloat vertex_vmax[3];
 };

  struct DomainNODE *DomainMoment = (struct DomainNODE *)mymalloc("DomainMoment", NTopleaves_groups * sizeof(struct DomainNODE));

  /* share the pseudo-particle data accross CPUs */
  int *recvcounts = (int *)mymalloc("recvcounts", sizeof(int) * NTask);
  int *recvoffset = (int *)mymalloc("recvoffset", sizeof(int) * NTask);
  int *bytecounts = (int *)mymalloc("bytecounts", sizeof(int) * NTask);
  int *byteoffset = (int *)mymalloc("byteoffset", sizeof(int) * NTask);

  for(int task = 0; task < NTask; task++)
    recvcounts[task] = 0;

  for(int n = 0; n < NTopleaves_groups; n++)
    recvcounts[DomainTask_groups[n]]++;

  for(int task = 0; task < NTask; task++)
    bytecounts[task] = recvcounts[task] * sizeof(struct DomainNODE);

  recvoffset[0] = 0, byteoffset[0] = 0;
  for(int task = 1; task < NTask; task++)
    {
      recvoffset[task] = recvoffset[task - 1] + recvcounts[task - 1];
      byteoffset[task] = byteoffset[task - 1] + bytecounts[task - 1];
    }

  struct DomainNODE *loc_DomainMoment =
      (struct DomainNODE *)mymalloc("loc_DomainMoment", recvcounts[ThisTask] * sizeof(struct DomainNODE));

  for(int n = 0, idx = 0; n < NTopleaves_groups; n++)
    {
      if(DomainTask_groups[n] == ThisTask)
        {
          int no = Ngb_DomainNodeIndex_groups[n];

          /* read out the multipole moments from the local base cells */
          for(int k = 0; k < 3; k++)
            {
              loc_DomainMoment[idx].range_min[k]   = Ngb_Nodes_groups[no].u.d.range_min[k];
              loc_DomainMoment[idx].range_max[k]   = Ngb_Nodes_groups[no].u.d.range_max[k];
              loc_DomainMoment[idx].vertex_vmin[k] = Ngb_Nodes_groups[no].vertex_vmin[k];
              loc_DomainMoment[idx].vertex_vmax[k] = Ngb_Nodes_groups[no].vertex_vmax[k];
            }

          idx++;
        }
    }

  MPI_Allgatherv(loc_DomainMoment, bytecounts[ThisTask], MPI_BYTE, DomainMoment, bytecounts, byteoffset, MPI_BYTE, MPI_COMM_WORLD);

  for(int task = 0; task < NTask; task++)
    recvcounts[task] = 0;

  for(int n = 0; n < NTopleaves_groups; n++)
    {
      int task = DomainTask_groups[n];
      if(task != ThisTask)
        {
          int no  = Ngb_DomainNodeIndex_groups[n];
          int idx = recvoffset[task] + recvcounts[task]++;

          for(int k = 0; k < 3; k++)
            {
              Ngb_Nodes_groups[no].u.d.range_min[k] = DomainMoment[idx].range_min[k];
              Ngb_Nodes_groups[no].u.d.range_max[k] = DomainMoment[idx].range_max[k];
              Ngb_Nodes_groups[no].vertex_vmin[k]   = DomainMoment[idx].vertex_vmin[k];
              Ngb_Nodes_groups[no].vertex_vmax[k]   = DomainMoment[idx].vertex_vmax[k];
            }
         Ngb_Nodes_groups[no].Ti_Current = All.Ti_Current;
        }
    }

  myfree(loc_DomainMoment);
  myfree(byteoffset);
  myfree(bytecounts);
  myfree(recvoffset);
  myfree(recvcounts);
  myfree(DomainMoment);
}



void ngb_finish_vounds_update_groups(int nchanged, int *nodelist)
{
  struct DomainNODE
  {
    int node;
    MyNgbTreeFloat vertex_vmin[3];
    MyNgbTreeFloat vertex_vmax[3];
  };

  /* share the pseudo-particle data accross CPUs */
  int *recvcounts = (int *)mymalloc("recvcounts", sizeof(int) * NTask);
  int *bytecounts = (int *)mymalloc("bytecounts", sizeof(int) * NTask);
  int *byteoffset = (int *)mymalloc("byteoffset", sizeof(int) * NTask);

  MPI_Allgather(&nchanged, 1, MPI_INT, recvcounts, 1, MPI_INT, MPI_COMM_WORLD);

  for(int task = 0; task < NTask; task++)
    bytecounts[task] = recvcounts[task] * sizeof(struct DomainNODE);

  byteoffset[0] = 0;
  for(int task = 1; task < NTask; task++)
    byteoffset[task] = byteoffset[task - 1] + bytecounts[task - 1];

  struct DomainNODE *loc_DomainMoment =
      (struct DomainNODE *)mymalloc("loc_DomainMoment", recvcounts[ThisTask] * sizeof(struct DomainNODE));

  for(int i = 0; i < nchanged; i++)
    {
      int no                   = nodelist[i];
      loc_DomainMoment[i].node = no;

      for(int j = 0; j < 3; j++)
        {
          loc_DomainMoment[i].vertex_vmin[j] = Ngb_Nodes_groups[no].vertex_vmin[j];
          loc_DomainMoment[i].vertex_vmax[j] = Ngb_Nodes_groups[no].vertex_vmax[j];
        }
    }

  int tot_nchanged = 0;
  for(int task = 0; task < NTask; task++)
    tot_nchanged += recvcounts[task];

  struct DomainNODE *tot_DomainMoment = (struct DomainNODE *)mymalloc("tot_DomainMoment", tot_nchanged * sizeof(struct DomainNODE));

  MPI_Allgatherv(loc_DomainMoment, bytecounts[ThisTask], MPI_BYTE, tot_DomainMoment, bytecounts, byteoffset, MPI_BYTE, MPI_COMM_WORLD);

  for(int i = 0; i < tot_nchanged; i++)
    {
      int no = tot_DomainMoment[i].node;

      for(int j = 0; j < 3; j++)
        {
          Ngb_Nodes_groups[no].vertex_vmin[j] = tot_DomainMoment[i].vertex_vmin[j];
          Ngb_Nodes_groups[no].vertex_vmax[j] = tot_DomainMoment[i].vertex_vmax[j];
        }

      no = Ngb_Nodes_groups[no].father;

      while(no >= 0)
        {
          int flag_changed = 0;

          for(int j = 0; j < 3; j++)
            {
              if(Ngb_Nodes_groups[no].vertex_vmin[j] > tot_DomainMoment[i].vertex_vmin[j])
                {
                  Ngb_Nodes_groups[no].vertex_vmin[j] = tot_DomainMoment[i].vertex_vmin[j];
                  flag_changed                 = 1;
                }

              if(Ngb_Nodes_groups[no].vertex_vmax[j] < tot_DomainMoment[i].vertex_vmax[j])
                {
                  Ngb_Nodes_groups[no].vertex_vmax[j] = tot_DomainMoment[i].vertex_vmax[j];
                  flag_changed                 = 1;
                }
            }

          if(flag_changed == 0)
            break;

          no = Ngb_Nodes_groups[no].father;
        }
    }

  myfree(tot_DomainMoment);
  myfree(loc_DomainMoment);
  myfree(byteoffset);
  myfree(bytecounts);
  myfree(recvcounts);
}


/*! \brief Allocates arrays for neighbor tree.
 *
 *  \return void
 */
void ngb_treeallocate_groups(void)
{
  if(Ngb_MaxPart_groups == 0)
    {
      Ngb_MaxPart_groups  = MaxNgroups;
      Ngb_MaxNodes_groups = (int)(All.NgbTreeAllocFactor * (MaxNgroups + BASENUMBER)) + NTopnodes_groups;
    }

  if(TotNgroups == 0)
    return;

  if(Ngb_Nodes_groups)
    terminate("already allocated");

  Ngb_DomainNodeIndex_groups = (int *)mymalloc_movable(&Ngb_DomainNodeIndex_groups, "Ngb_DomainNodeIndex_groups", NTopleaves_groups * sizeof(int));
  Ngb_Nodes_groups = (struct NgbNODE *)mymalloc_movable(&Ngb_Nodes_groups, "Ngb_Nodes_groups", (Ngb_MaxNodes_groups + 1) * sizeof(struct NgbNODE));
  Ngb_Nodes_groups -= Ngb_MaxPart_groups;

  Ngb_Nextnode_groups = (int *)mymalloc_movable(&Ngb_Nextnode_groups, "Ngb_Nextnode_groups", (Ngb_MaxPart_groups + NTopleaves_groups) * sizeof(int));
  Ngb_Father_groups   = (int *)mymalloc_movable(&Ngb_Father_groups, "Ngb_Father_groups", Ngb_MaxPart_groups * sizeof(int));

  Ngb_Marker_groups = (int *)mymalloc_movable(&Ngb_Marker_groups, "Ngb_Marker_groups", (Ngb_MaxNodes_groups + Ngb_MaxPart_groups) * sizeof(int));

#if NUM_THREADS > 1
  Ngb_NodeLocks_groups = (omp_lock_t *)mymalloc("Ngb_NodeLocks_groups", Ngb_MaxNodes_groups * sizeof(omp_lock_t));
  int i;
  for(i = 0; i < Ngb_MaxNodes_groups; i++)
    omp_init_lock(&Ngb_NodeLocks_groups[i]);
  Ngb_NodeLocks_groups -= Ngb_MaxPart_groups;
#endif
}

/*! \brief This function frees the memory allocated for the neighbor tree.
 *
 *  \return void
 */
void ngb_treefree_groups(void)
{
  if(TotNgroups == 0)
    return;

  if(Ngb_Nodes_groups)
    {
#if NUM_THREADS > 1
      Ngb_NodeLocks_groups += Ngb_MaxPart_groups;
      myfree(Ngb_NodeLocks_groups);
#endif
      myfree_movable(Ngb_Marker_groups);
      myfree_movable(Ngb_Father_groups);
      myfree_movable(Ngb_Nextnode_groups);
      Ngb_Nodes_groups += Ngb_MaxPart_groups;
      myfree_movable(Ngb_Nodes_groups);
      myfree_movable(Ngb_DomainNodeIndex_groups);

      Ngb_MaxPart_groups  = 0;
      Ngb_MaxNodes_groups = 0;
    }
  else
    terminate("trying to free the tree even though it's not allocated");
}

/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/forcetree_walk.c
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
#include <time.h>

#include "allvars.h"
#include "proto.h"
#if defined(RADCOOL) && !defined(GFM)
#error "Need to compile with GFM if RADCOOL option is chosen"
#endif

#ifdef RADCOOL_HOTHALO
static double effectivemass;
#endif

/*! \file forcetree_walk.c
 *  \brief Gravitational tree walk code
 *
 *  This file contains the various gravitational tree walks.
 */

/*! \brief variable for short-range lookup table
 *
 *  contains the factor needed for the short range
 *  contribution of the tree to the gravity force
 */

static float shortrange_table[NTAB + 1];

/*! \brief variable for short-range lookup table
 *
 *  contains the factor needed for the short range
 *  contribution of the tree to the potential energy
 */
static float shortrange_table_potential[NTAB + 1];

/*! \brief Initializes the short range table.
 *
 *  The short range table contains the complementary error function
 *  needed for the computation of the short range part of the gravity
 *  force/potential in case of the TreePM algorithm.
 *
 *  \return void
 */


/*! \brief This routine calculates the (short range) force contribution
 *   for a given particle in case the Tree(PM) algorithm is used.
 *
 *  In the TreePM algorithm, the tree is walked only locally around the
 *  target coordinate.  Tree nodes that fall outside a box of half
 *  side-length Rcut= RCUT*ASMTH*MeshSize can be discarded. The short-range
 *  potential is modified by a complementary error function, multiplied
 *  with the Newtonian form. The resulting short-range suppression compared
 *  to the Newtonian force is tabulated, because looking up from this table
 *  is faster than recomputing the corresponding factor, despite the
 *  memory-access penalty (which reduces cache performance) incurred by the
 *  table.
 *
 *  Depending on the value of TypeOfOpeningCriterion, either the geometrical BH
 *  cell-opening criterion, or the `relative' opening criterion is used.
 *
 *  \param[in] in Gravdata communicated into function.
 *  \param[in, out] out Gravdata communicated from function.
 *  \param[in] target Index of the particle to be processed.
 *  \param[in] mode 0: process local particle (phase 1), 1: process imported
 *             particle (phase 2).
 *  \param[in] thread_id Id of this thread.
 *  \param[in, out] firstnode First node involved in this algorithm.
 *  \param[in] measure_cost_flag Whether the cost of the tree walk should be
 *             measured.
 *
 *  \return Number of interactions processed for particle i.
 */
int tracerBH_neighbors_force_treeevaluate(gravdata_in *in, gravdata_out *out, int target, int mode, int thread_id, int numnodes, int *firstnode,
                       int measure_cost_flag)
{
  struct NODE *nop = NULL;

  double acc_x = 0;
  double acc_y = 0;
  double acc_z = 0;

  int no_of_BHs_ngb = 0;

  int ninteractions = 0;

  double pos_x = in->Pos[0];
  double pos_y = in->Pos[1];
  double pos_z = in->Pos[2];
  double aold  = All.ErrTolForceAcc * in->OldAcc;
  MyDouble HostHaloMass  = in->HostHaloMass;

  double vvir  = pow(10 * All.G * All.cf_H * HostHaloMass, 1.0 / 3);
  double rvir_phys = vvir / (10 * All.cf_H);     /* physical */
  double rvir      = (rvir_phys / All.cf_atime); /* comoving */
  double rcut = All.MaximumNeighboringTracerDistance * rvir;
  double rcut2 = rcut * rcut;

  for(int k = 0; k < numnodes; k++)
    {
      int no;

      if(mode == 0)
        no = Tree_MaxPart; /* root node */
      else
        {
          no = firstnode[k];
          no = Nodes[no].u.d.nextnode; /* open it */
        }

      while(no >= 0)
        {
          double dx, dy, dz, r2, mass, hmax;
          int no_of_BHs = 0;
          if(no < Tree_MaxPart) /* single particle */
            {
              dx = GRAVITY_NEAREST_X(Tree_Pos_list[3 * no + 0] - pos_x);
              dy = GRAVITY_NEAREST_Y(Tree_Pos_list[3 * no + 1] - pos_y);
              dz = GRAVITY_NEAREST_Z(Tree_Pos_list[3 * no + 2] - pos_z);
              r2 = dx * dx + dy * dy + dz * dz;
              mass = P[no].Mass;    

 
#if(PTYPE_USED_FOR_ENVIRONMENT_BASED_SEEDING == 0)
              if(!(P[no].Type == 0))
                  terminate("\n Particle type %d was inserted in force tree, but it should only contain gas \n", P[no].Type);
              if(!(SphP[no].IsThisTheDensestCell == 1))
                  terminate("\n One of the gas particles in the force tree is not the densest gas in an FOF \n");
              if(r2 < rcut2)
                 if (SphP[no].HostHaloMass > HostHaloMass)
                   no_of_BHs = 1;
#else
              if(!(P[no].Type == 1))
                  terminate("\n Particle type %d was inserted in force tree, but it should only contain DM \n", P[no].Type);
       	      if(!(P[no].IsThisTheMinPotential == 1))
                  terminate("\n One of the DM particles in the force tree is not the min potential in an FOF \n");
              if(r2 < rcut2)
                 if (P[no].HostHaloMass > HostHaloMass)
                   no_of_BHs = 1;
#endif
              if(measure_cost_flag)
                Thread[thread_id].P_CostCount[no]++;

              double h_j = All.ForceSoftening[P[no].SofteningType];
              no = Nextnode[no];
            }
          else if(no < Tree_MaxPart + Tree_MaxNodes) /* we have an  internal node */
            {
              if(mode == 1)
                {
                  if(no <
                     Tree_FirstNonTopLevelNode) /* we reached a top-level node again, which means that we are done with the branch */
                    {
                      no = -1;
                      continue;
                    }
                }

              nop = &Nodes[no];

              mass = nop->u.d.mass;
              dx   = GRAVITY_NEAREST_X(nop->u.d.s[0] - pos_x);
              dy   = GRAVITY_NEAREST_Y(nop->u.d.s[1] - pos_y);
              dz   = GRAVITY_NEAREST_Z(nop->u.d.s[2] - pos_z);

              int no_of_BHs_internal = nop->u.d.no_of_BHs;
              r2 = dx * dx + dy * dy + dz * dz;

//              if(no_of_BHs_internal > 1) /* check Barnes-Hut opening criterion */
//                {
                     /* open cell */
                      no = nop->u.d.nextnode;
                      continue;
//                }

                /* ok, node can be used */
             if(measure_cost_flag && mass)
                Thread[thread_id].Node_CostCount[no]++;
              no = nop->u.d.sibling;
            }
          else if(no >= Tree_ImportedNodeOffset) /* point from imported nodelist */
            {
              int n = no - Tree_ImportedNodeOffset;

              dx = GRAVITY_NEAREST_X(Tree_Points[n].Pos[0] - pos_x);
              dy = GRAVITY_NEAREST_Y(Tree_Points[n].Pos[1] - pos_y);
              dz = GRAVITY_NEAREST_Z(Tree_Points[n].Pos[2] - pos_z);

              r2 = dx * dx + dy * dy + dz * dz;

              mass = Tree_Points[n].Mass;

#if(PTYPE_USED_FOR_ENVIRONMENT_BASED_SEEDING == 0)
              if(!(Tree_Points[n].Type == 0))
                  terminate("\n Tree_Points type %d was inserted in force tree, but it should only contain gas \n", Tree_Points[n].Type);
              if(!(Tree_Points[n].IsThisTheDensestCell == 1))
                  terminate("\n One of the gas particles in the force tree is not the densest cell in an FOF \n");
#else
              if(!(Tree_Points[n].Type == 1))
                  terminate("\n Tree_Points type %d was inserted in force tree, but it should only contain DM \n", Tree_Points[n].Type);
              if(!(Tree_Points[n].IsThisTheMinPotential == 1))
                  terminate("\n One of the DM particles in the force tree is not the minimum potential in an FOF \n");
#endif
              if(r2 < rcut2)
                 if (Tree_Points[n].HostHaloMass > HostHaloMass)
                    no_of_BHs = Tree_Points[n].no_of_BHs;
              if(measure_cost_flag)
                 Thread[thread_id].TreePoints_CostCount[n]++;

              no = Nextnode[no - Tree_MaxNodes];
            }
          else /* pseudo particle */
            {
              if(mode == 0)
                {
                  tree_treefind_export_node_threads(no, target, thread_id);
                }

              no = Nextnode[no - Tree_MaxNodes];
              continue;
            }
          /* now evaluate the multipole moment */
          if(no_of_BHs)
            {
                no_of_BHs_ngb += no_of_BHs;
            }
        }
    }
 out->no_of_BHs_ngb = no_of_BHs_ngb;
 return ninteractions;
}



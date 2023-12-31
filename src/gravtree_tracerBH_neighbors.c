/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/gravtree.c
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
#include <stdlib.h>
#include <string.h>

#include "allvars.h"
#include "domain.h"
#include "proto.h"

/*! \file gravtree.c
 *  \brief main driver routines for gravitational (short-range) force computation
 *
 *  This file contains the code for the gravitational force computation by
 *  means of the tree algorithm. To this end, a tree force is computed for all
 *  active local particles, and particles are exported to other processors if
 *  needed, where they can receive additional force contributions. If the
 *  TreePM algorithm is enabled, the force computed will only be the
 *  short-range part.
 */

static double ThreadsCosttotal[NUM_THREADS]; /*!< The total cost of the particles/nodes processed by each thread */
static int ThreadFirstExec[NUM_THREADS]; /*!< Keeps track, if a given thread executes the gravity_primary_loop() for the first time */
static int MeasureCostFlag;              /*!< Whether the tree costs are measured for the current time step */

static int tracerBH_neighbors_evaluate(int target, int mode, int threadid);

typedef gravdata_in data_in;

typedef gravdata_out data_out;

#ifdef DETAILEDTIMINGS
static double tstart;
static int current_timebin;
#endif

/* local data structure for collecting particle/cell data that is sent to other processors if needed */
static data_in *DataGet;

static data_out *DataResult;

/*! \brief Routine that fills the relevant particle/cell data into the input
 *         structure defined above. Needed by generic_comm_helpers2.
 *
 *  \param[out] in Data structure to fill.
 *  \param[in] i Index of particle in P and SphP arrays.
 *  \param[in] firstnode First note of communication.
 *
 *  \return void
 */
static void particle2in(data_in *in, int i, int firstnode)
{
  if(i < NumPart)
    {
#ifdef CELL_CENTER_GRAVITY
      if(P[i].Type == 0)
        {
          for(int k = 0; k < 3; k++)
            in->Pos[k] = SphP[i].Center[k];
        }
      else
#endif
        {
          for(int k = 0; k < 3; k++)
            in->Pos[k] = P[i].Pos[k];
        }

#ifdef TREECOLV2_VEL
      for(int k = 0; k < 3; k++)
        in->Vel[k] = P[i].Vel[k];

      if(P[i].Type == 0)
        in->Vth2 = SphP[i].Vth2;
#endif

      in->Type          = P[i].Type;
      in->SofteningType = P[i].SofteningType;
      in->OldAcc        = P[i].OldAcc;
      if (P[i].Type == 5)
        in->HostHaloMass = BPP(i).HostHaloMass; 
      else
#if(PTYPE_USED_FOR_ENVIRONMENT_BASED_SEEDING == 0)
        in->HostHaloMass = SphP[i].HostHaloMass;
#else
        in->HostHaloMass = P[i].HostHaloMass; 
#endif
    }
  else
    {
      i -= Tree_ImportedNodeOffset;

      for(int k = 0; k < 3; k++)
        in->Pos[k] = Tree_Points[i].Pos[k];

#ifdef TREECOLV2_VEL
      for(int k = 0; k < 3; k++)
        in->Vel[k] = Tree_Points[i].Vel[k];
#endif

      in->Type          = Tree_Points[i].Type;
      in->SofteningType = Tree_Points[i].SofteningType;
      in->OldAcc        = Tree_Points[i].OldAcc;

      if (Tree_Points[i].Type == 5)
        in->HostHaloMass = Tree_Points[i].HostHaloMass;
      else
        in->HostHaloMass = Tree_Points[i].HostHaloMass;
    }
  in->Firstnode = firstnode;
}

/*! \brief Routine to store or combine result data. Needed by
 *         generic_comm_helpers2.
 *
 *  \param[in] out Data to be moved to appropriate variables in global
 *  particle and cell data arrays (P, SphP,...)
 *  \param[in] i Index of particle in P and SphP arrays
 *  \param[in] mode Mode of function: local particles or information that was
 *  communicated from other tasks and has to be added locally?
 *
 *  \return void
 */
static void out2particle(data_out *out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES) /* initial store */
    {
      if(i < NumPart)
        {
#ifdef SEED_HALO_ENVIRONMENT_CRITERION
         P[i].no_of_BHs_ngb = out->no_of_BHs_ngb;        
#endif

        }
      else
        {
          int idx                                      = Tree_ResultIndexList[i - Tree_ImportedNodeOffset];
#ifdef SEED_HALO_ENVIRONMENT_CRITERION
          Tree_ResultsActiveImported[idx].no_of_BHs_ngb = out->no_of_BHs_ngb;
#endif

       }
    }
  else /* combine */
    {
      if(i < NumPart)
        {
#ifdef SEED_HALO_ENVIRONMENT_CRITERION
          P[i].no_of_BHs_ngb += out->no_of_BHs_ngb;
#endif
        }
      else
        {
          int idx = Tree_ResultIndexList[i - Tree_ImportedNodeOffset];
#ifdef SEED_HALO_ENVIRONMENT_CRITERION
          Tree_ResultsActiveImported[idx].no_of_BHs_ngb += out->no_of_BHs_ngb;
#endif
        }
    }
}

#define OMIT_GENERIC_COMM_PATTERN_FOR_GIVEN_PARTICLES
#include "generic_comm_helpers2.h"

/*! \brief Primary loop of gravity calculation.
 *
 *  Gravitational interactions between local particles; see gravity_tree(..).
 *  Equivalent to 'kernel_local'.
 *
 *  \return void
 */
static void gravity_primary_loop_tracerBH_neighbors(void)
{
  TIMER_STOPSTART(CPU_TREEBALSNDRCV, CPU_TREEWALK1);

#ifdef DETAILEDTIMINGS
  double t0 = second();
#endif

  int idx;
#ifdef GENERIC_ASYNC
  int flag = 0;
#endif
  /* do local particles */
#pragma omp parallel private(idx)
  {
    int threadid = get_thread_num();
    double cost  = 0;
#ifdef GENERIC_ASYNC
    int count = 0;
#endif

    if(ThreadFirstExec[threadid])
      {
        ThreadFirstExec[threadid] = 0;

        if(MeasureCostFlag)
          {
            memset(Thread[threadid].P_CostCount, 0, NumPart * sizeof(int));
            memset(Thread[threadid].TreePoints_CostCount, 0, Tree_NumPartImported * sizeof(int));
            memset(Thread[threadid].Node_CostCount + Tree_MaxPart, 0, Tree_NumNodes * sizeof(int));
          }
      }

    for(int j = 0; j < NTask; j++)
      Thread[threadid].Exportflag[j] = -1;

    while(1)
      {
        if(Thread[threadid].ExportSpace < MinSpace)
          break;

#ifdef GENERIC_ASYNC
        if(threadid == 0)
          {
            if((count & POLLINGINTERVAL) == 0)
              if(generic_polling_primary(count, Nforces))
                flag = 1;

            count++;
          }

        if(flag)
          break;
#endif

#pragma omp atomic capture
        idx = NextParticle++;

        if(idx >= Nforces)
          break;

        int i = TargetList[idx];

        cost += tracerBH_neighbors_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
      }

    ThreadsCosttotal[threadid] += cost;
  }
  TIMER_STOPSTART(CPU_TREEWALK1, CPU_TREEBALSNDRCV);
}

/*! \brief Secondary loop of gravity calculation.
 *
 *  Gravitational interactions between imported particles; see gravity_tree(.).
 *  Equivalent to 'kernel_imported'.
 *
 *  \return void
 */
void gravity_secondary_loop_tracerBH_neighbors(void)
{
  TIMER_STOPSTART(CPU_TREEBALSNDRCV, CPU_TREEWALK2);

#ifdef DETAILEDTIMINGS
  double t0 = second();
#endif

  /* now do the particles that were sent to us */
  int i, cnt = 0;
#pragma omp parallel private(i)
  {
    int threadid = get_thread_num();
    double cost  = 0;
#ifdef GENERIC_ASYNC
    int count = 0;
#endif

    if(ThreadFirstExec[threadid])
      {
        ThreadFirstExec[threadid] = 0;

        if(MeasureCostFlag)
          {
            memset(Thread[threadid].P_CostCount, 0, NumPart * sizeof(int));
            memset(Thread[threadid].TreePoints_CostCount, 0, Tree_NumPartImported * sizeof(int));
            memset(Thread[threadid].Node_CostCount + Tree_MaxPart, 0, Tree_NumNodes * sizeof(int));
          }
      }

    while(1)
      {
#pragma omp atomic capture
        i = cnt++;

        if(i >= Nimport)
          break;

#ifdef GENERIC_ASYNC
        if(threadid == 0)
          {
            if((count & POLLINGINTERVAL) == 0)
              generic_polling_secondary();
          }

        count++;
#endif

        cost += tracerBH_neighbors_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }

    ThreadsCosttotal[threadid] += cost;
  }

#ifdef DETAILEDTIMINGS
  double t1 = second();

  fprintf(FdDetailed, "%d %d %d %d %g %g\n", All.NumCurrentTiStep, current_timebin, DETAILED_TIMING_GRAVWALK, MODE_IMPORTED_PARTICLES,
          timediff(tstart, t0), timediff(tstart, t1));
#endif

  TIMER_STOPSTART(CPU_TREEWALK2, CPU_TREEBALSNDRCV);
}

/*! \brief This function computes the gravitational forces for all active
 *         particles.
 *
 *  The tree walk is done in two phases: First the local part of the force tree
 *  is processed (gravity_primary_loop() ). Whenever an external node is
 *  encountered during the walk, this node is saved on a list. This node list
 *  along with data about the particles is then exchanged among tasks.
 *  In the second phase (gravity_secondary_loop() ) each task now continues
 *  the tree walk for the imported particles. Finally the resulting partial
 *  forces are send back to the original task and are summed up there to
 *  complete the tree force calculation.
 *
 *  If only the tree algorithm is used in a periodic box, the whole tree walk
 *  is done twice. First a normal tree walk is done as described above, and
 *  afterwards a second tree walk, which adds the needed Ewald corrections is
 *  performed.
 *
 *  Particles are only exported to other processors when really needed,
 *  thereby allowing a good use of the communication buffer. Every particle is
 *  sent at most once to a given processor together with the complete list of
 *  relevant tree nodes to be checked on the other task.
 *
 *  Particles which drifted into the domain of another task are sent to this
 *  task for the force computation. Afterwards the resulting force is sent
 *  back to the originating task.
 *
 *  In order to improve the work load balancing during a domain decomposition,
 *  the work done by each node/particle is measured. The work is measured for
 *  the interaction partners (i.e. the nodes or particles) and not for the
 *  particles itself that require a force computation. This way, work done for
 *  imported particles is accounted for at the task where the work actually
 *  incurred. The cost measurement is only done for the "GRAVCOSTLEVELS"
 *  highest occupied time bins. The variable 'MeasureCostFlag' will state
 *  whether a measurement is done at the present time step.
 *
 *  The tree imbalance can be further reduced using the CHUNKING option. The
 *  particles requiring a force computation are split into chunks of size
 *  'Nchunksize'. A set of every 'Nchunk' -th chunk is processed first.
 *  Then the process is repeated, processing the next set of chunks. This way
 *  the amount of exported particles is more balanced, as communication heavy
 *  regions are mixed with less communication intensive regions.
 *
 * \param[in] timebin Time bin for which gravity should be calculated.
 *
 * \return void
 */
void gravity_tree_tracerBH_neighbors(int timebin)
{
  int idx, i, k, ncount, iter = 0, maxiter;
  int IncludeParticleInTreeWalk;
  struct detailed_timings
  {
    double all, tree1, tree2, tree, commwait;
    double sumnexport, costtotal, numnodes;
  } timer, tisum, timax;
  memset(&timer, 0, sizeof(struct detailed_timings));
  double Costtotal;
  int ngrp;
  int recvTask;

  TIMER_STORE;
  TIMER_START(CPU_TREE);

  /* allocate buffers to arrange communication */
  mpi_printf("GRAVTREE: Begin tree force.  (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));

  TIMER_STOPSTART(CPU_TREE, CPU_TREECOSTMEASURE);

  for(i = 0; i < NUM_THREADS; i++)
    {
      ThreadsCosttotal[i] = 0;
      ThreadFirstExec[i]  = 0;
    }

  /* find the level (if any) for which we measure gravity cost */
  for(i = 0, TakeLevel = -1; i < GRAVCOSTLEVELS; i++)
    {
      if(All.LevelToTimeBin[i] == timebin)
        {
          TakeLevel = i;
          break;
        }
    }

  if(TakeLevel < 0) /* we have not found a matching slot */
    {
      for(i = 0; i < GRAVCOSTLEVELS; i++)
        {
          if(All.LevelToTimeBin[i] < 0)
            {
              All.LevelToTimeBin[i]       = timebin;
              TakeLevel                   = i;
              All.LevelHasBeenMeasured[i] = 0;
              break;
            }
        }

      if(TakeLevel < 0)
        {
          if(All.HighestOccupiedGravTimeBin - timebin < GRAVCOSTLEVELS) /* we should have space */
            {
              /* clear levels that are out of range */
              for(i = 0; i < GRAVCOSTLEVELS; i++)
                {
                  if(All.LevelToTimeBin[i] > All.HighestOccupiedGravTimeBin)
                    {
                      All.LevelToTimeBin[i]       = timebin;
                      TakeLevel                   = i;
                      All.LevelHasBeenMeasured[i] = 0;
                      break;
                    }
                  if(All.LevelToTimeBin[i] < All.HighestOccupiedGravTimeBin - (GRAVCOSTLEVELS - 1))
                    {
                      All.LevelToTimeBin[i]       = timebin;
                      TakeLevel                   = i;
                      All.LevelHasBeenMeasured[i] = 0;
                      break;
                    }
                }

              if(TakeLevel < 0)
                {
                  if(timebin > All.HighestOccupiedGravTimeBin)
                    {
                      for(i = 0; i < GRAVCOSTLEVELS; i++)
                        {
                          if(All.LevelToTimeBin[i] == All.HighestOccupiedGravTimeBin)
                            {
                              All.LevelToTimeBin[i]       = timebin;
                              TakeLevel                   = i;
                              All.LevelHasBeenMeasured[i] = 0;
                              break;
                            }
                        }
                    }
                }

              if(TakeLevel < 0)
                {
                  mpi_printf("All.HighestOccupiedGravTimeBin=%d   timebin=%d\n", All.HighestOccupiedGravTimeBin, timebin);
                  for(i = 0; i < GRAVCOSTLEVELS; i++)
                    {
                      mpi_printf("All.LevelToTimeBin[i=%d]=%d\n", i, All.LevelToTimeBin[i]);
                    }

                  fflush(stdout);
                  MPI_Barrier(MPI_COMM_WORLD);
                  fflush(stdout);
                  MPI_Barrier(MPI_COMM_WORLD);

                  terminate("TakeLevel=%d < 0", TakeLevel);
                }
            }
        }
    }

  MeasureCostFlag = 0;

  if(TakeLevel >= 0)
    if(All.LevelHasBeenMeasured[TakeLevel] == 0)
      {
        MeasureCostFlag = 1;

        Thread[0].P_CostCount          = (int *)mymalloc("Thread_P_CostCount", NumPart * sizeof(int));
        Thread[0].TreePoints_CostCount = (int *)mymalloc("Threads_TreePoints_CostCount", Tree_NumPartImported * sizeof(int));
        Thread[0].Node_CostCount       = (int *)mymalloc("Threads_Node_CostCount", Tree_NumNodes * sizeof(int));

        for(i = 1; i < NUM_THREADS; i++)
          {
            Thread[i].P_CostCount          = (int *)mymalloc("Threads_P_CostCount", NumPart * sizeof(int));
            Thread[i].TreePoints_CostCount = (int *)mymalloc("Threads_TreePoints_CostCount", Tree_NumPartImported * sizeof(int));
            Thread[i].Node_CostCount       = (int *)mymalloc("Threads_Node_CostCount", Tree_NumNodes * sizeof(int));
          }

        for(i = 0; i < NUM_THREADS; i++)
          Thread[i].Node_CostCount -= Tree_MaxPart;

        for(i = 0; i < NUM_THREADS; i++)
          ThreadFirstExec[i] = 1;
      }

  TIMER_STOPSTART(CPU_TREECOSTMEASURE, CPU_TREE);

  /* Create list of targets. We do this here to simplify the treatment of the two possible sources of points */

  TargetList           = (int *)mymalloc("TargetList", (NumPart + Tree_NumPartImported) * sizeof(int));
  Tree_ResultIndexList = (int *)mymalloc("Tree_ResultIndexList", Tree_NumPartImported * sizeof(int));

  Nforces = 0;
  
  for(i = 0; i < NumPart; i++)
    P[i].no_of_BHs_ngb = 0;

  for(idx = 0; idx < TimeBinsGravity.NActiveParticles; idx++)
    {
      i = TimeBinsGravity.ActiveParticleList[idx];
      if(i < 0)
        continue;

#ifdef CREATE_SUBFOFS
      if(All.SubFOF_mode == 1)          
#if(PTYPE_USED_FOR_ENVIRONMENT_BASED_SEEDING==0)
          IncludeParticleInTreeWalk = SphP[i].IsThisTheDensestCell_bFOF;
#else
       	  IncludeParticleInTreeWalk = P[i].IsThisTheMinPotential_bFOF;
#endif
      else if(All.SubFOF_mode == 0)
#if(PTYPE_USED_FOR_ENVIRONMENT_BASED_SEEDING==0)
          IncludeParticleInTreeWalk = SphP[i].IsThisTheDensestCell;
#else
          IncludeParticleInTreeWalk = P[i].IsThisTheMinPotential;
#endif
#endif

#ifndef CREATE_SUBFOFS
#if(PTYPE_USED_FOR_ENVIRONMENT_BASED_SEEDING==0)
      IncludeParticleInTreeWalk = SphP[i].IsThisTheDensestCell;
#else
      IncludeParticleInTreeWalk = P[i].IsThisTheMinPotential;
#endif
#endif


      if(IncludeParticleInTreeWalk == 0)
        continue;

#ifdef TRACER_PARTICLE
      if(P[i].Type == TRACER_PARTICLE)
        continue;
#endif
      if(Tree_Task_list[i] == ThisTask)
        TargetList[Nforces++] = i;
    }

  for(i = 0, ncount = 0; i < Tree_NumPartImported; i++)
#ifndef HIERARCHICAL_GRAVITY
    if(Tree_Points[i].ActiveFlag)
#endif
      {
        Tree_ResultIndexList[i] = ncount++;
        TargetList[Nforces++]   = i + Tree_ImportedNodeOffset;
      }

  Tree_ResultsActiveImported =
      (struct resultsactiveimported_data *)mymalloc("Tree_ResultsActiveImported", ncount * sizeof(struct resultsactiveimported_data));


#ifdef CHUNKING
  permutate_chunks_in_list(Nforces, TargetList);
#endif

  generic_set_MaxNexport();

  /******************************************/
  /* now execute the tree walk calculations */
  /******************************************/

  TIMER_STOPSTART(CPU_TREE, CPU_TREEBALSNDRCV);

#ifdef DETAILEDTIMINGS
  tstart          = second();
  current_timebin = timebin;
#endif

  iter = generic_comm_pattern(Nforces, gravity_primary_loop_tracerBH_neighbors, gravity_secondary_loop_tracerBH_neighbors);

  /* now communicate the forces in Tree_ResultsActiveImported */

  TIMER_STOPSTART(CPU_TREEBALSNDRCV, CPU_TREESENDBACK);

#ifdef DETAILEDTIMINGS
  double tend = second();

  fprintf(FdDetailed, "%d %d %d %d %g %g\n", All.NumCurrentTiStep, current_timebin, DETAILED_TIMING_GRAVWALK, MODE_FINISHED,
          timediff(tstart, tend), timediff(tstart, tend));
  fflush(FdDetailed);
#endif

  for(int j = 0; j < NTask; j++)
    Recv_count[j] = 0;

  int n;
  for(i = 0, n = 0, k = 0; i < NTask; i++)
    for(int j = 0; j < Force_Recv_count[i]; j++, n++)
      {
#ifndef HIERARCHICAL_GRAVITY
        if(Tree_Points[n].ActiveFlag)
#endif
          {
            Tree_ResultsActiveImported[k].index = Tree_Points[n].index;
            Recv_count[i]++;
            k++;
          }
      }

  MPI_Alltoall(Recv_count, 1, MPI_INT, Send_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(i = 0, Nexport = 0, Nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; i < NTask; i++)
    {
      Nexport += Send_count[i];
      Nimport += Recv_count[i];

      if(i > 0)
        {
          Send_offset[i] = Send_offset[i - 1] + Send_count[i - 1];
          Recv_offset[i] = Recv_offset[i - 1] + Recv_count[i - 1];
        }
    }

  struct resultsactiveimported_data *tmp_results =
      (struct resultsactiveimported_data *)mymalloc("tmp_results", Nexport * sizeof(struct resultsactiveimported_data));
  memset(tmp_results, -1, Nexport * sizeof(struct resultsactiveimported_data));

  /* exchange  data */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              MPI_Sendrecv(&Tree_ResultsActiveImported[Recv_offset[recvTask]],
                           Recv_count[recvTask] * sizeof(struct resultsactiveimported_data), MPI_BYTE, recvTask, TAG_FOF_A,
                           &tmp_results[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct resultsactiveimported_data),
                           MPI_BYTE, recvTask, TAG_FOF_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  for(i = 0; i < Nexport; i++)
    {
      int target = tmp_results[i].index;
#ifdef SEED_HALO_ENVIRONMENT_CRITERION
      P[target].no_of_BHs_ngb = tmp_results[i].no_of_BHs_ngb;
#endif
    }

  myfree(tmp_results);

  myfree(Tree_ResultsActiveImported);
  myfree(Tree_ResultIndexList);
  myfree(TargetList);

  TIMER_STOPSTART(CPU_TREESENDBACK, CPU_TREECOSTMEASURE);

  /* assign node cost to particles */
  if(MeasureCostFlag)
    {
      for(int threadid = 0; threadid < NUM_THREADS; threadid++)
        if(ThreadFirstExec[threadid])
          {
            /* this could happen if neither the primary nor the secondary loop had anything to do */
            ThreadFirstExec[threadid] = 0;
            memset(Thread[threadid].P_CostCount, 0, NumPart * sizeof(int));
            memset(Thread[threadid].TreePoints_CostCount, 0, Tree_NumPartImported * sizeof(int));
            memset(Thread[threadid].Node_CostCount + Tree_MaxPart, 0, Tree_NumNodes * sizeof(int));
          }

      force_assign_cost_values();
      domain_init_sum_cost();

      All.LevelHasBeenMeasured[TakeLevel] = 1;

      if(All.TypeOfOpeningCriterion == 1 && All.Ti_Current == 0)
        All.LevelHasBeenMeasured[TakeLevel] = 0;

      for(i = 0; i < NUM_THREADS; i++)
        Thread[i].Node_CostCount += Tree_MaxPart;

      for(i = NUM_THREADS - 1; i >= 1; i--)
        {
          myfree(Thread[i].Node_CostCount);
          myfree(Thread[i].TreePoints_CostCount);
          myfree(Thread[i].P_CostCount);
        }

      myfree(Thread[0].Node_CostCount);
      myfree(Thread[0].TreePoints_CostCount);
      myfree(Thread[0].P_CostCount);
    }

  TIMER_STOPSTART(CPU_TREECOSTMEASURE, CPU_TREE);

  if(All.TypeOfOpeningCriterion == 1)
    /* This will switch to the relative opening criterion for the following
     * force computations */
    All.ErrTolTheta = 0;

  mpi_printf("GRAVTREE: tree-force is done.\n");

  /*  gather some diagnostic information */

  TIMER_STOPSTART(CPU_TREE, CPU_LOGS);

  Costtotal = 0;
  for(i = 0; i < NUM_THREADS; i++)
    Costtotal += ThreadsCosttotal[i];

  timer.tree1      = TIMER_DIFF(CPU_TREEWALK1);
  timer.tree2      = TIMER_DIFF(CPU_TREEWALK2);
  timer.tree       = timer.tree1 + timer.tree2;
  timer.commwait   = TIMER_DIFF(CPU_TREEBALSNDRCV) + TIMER_DIFF(CPU_TREESENDBACK);
  timer.all        = timer.tree + timer.commwait + TIMER_DIFF(CPU_TREE) + TIMER_DIFF(CPU_TREECOSTMEASURE);
  timer.sumnexport = SumNexport;
  timer.costtotal  = Costtotal;
  timer.numnodes   = Tree_NumNodes;

  MPI_Reduce(&iter, &maxiter, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce((double *)&timer, (double *)&tisum, (int)(sizeof(struct detailed_timings) / sizeof(double)), MPI_DOUBLE, MPI_SUM, 0,
             MPI_COMM_WORLD);
  MPI_Reduce((double *)&timer, (double *)&timax, (int)(sizeof(struct detailed_timings) / sizeof(double)), MPI_DOUBLE, MPI_MAX, 0,
             MPI_COMM_WORLD);

  All.TotNumOfForces += TimeBinsGravity.GlobalNActiveParticles;

#ifdef TGSET
  if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin)
#endif
    if(ThisTask == 0)
      {
        fprintf(FdTimings, "Nf=%9lld  timebin=%d  total-Nf=%lld\n", TimeBinsGravity.GlobalNActiveParticles, timebin,
                All.TotNumOfForces);

        fprintf(FdTimings, "   work-load balance: %g (%g %g), rel1to2: %g\n", timax.tree / ((tisum.tree + 1e-20) / NTask),
                timax.tree1 / ((tisum.tree1 + 1e-20) / NTask), timax.tree2 / ((tisum.tree2 + 1e-20) / NTask),
                tisum.tree1 / (tisum.tree1 + tisum.tree2 + 1e-20));
        fprintf(FdTimings, "   number of iterations:  max=%d, exported fraction: %g\n", maxiter,
                tisum.sumnexport / (TimeBinsGravity.GlobalNActiveParticles + 1e-20));
        fprintf(FdTimings, "   part/sec: raw=%g, effective=%g     ia/part: avg=%g\n",
                TimeBinsGravity.GlobalNActiveParticles / (tisum.tree + 1.0e-20),
                TimeBinsGravity.GlobalNActiveParticles / ((timax.tree + 1.0e-20) * NTask),
                tisum.costtotal / (TimeBinsGravity.GlobalNActiveParticles + 1.0e-20));

        fprintf(FdTimings, "   maximum number of nodes: %g, filled: %g\n", timax.numnodes, timax.numnodes / Tree_MaxNodes);

        fprintf(FdTimings, "   avg times: all=%g  tree1=%g  tree2=%g  commwait=%g sec\n", tisum.all / NTask, tisum.tree1 / NTask,
                tisum.tree2 / NTask, tisum.commwait / NTask);

        myflush(FdTimings);
      }

  TIMER_STOP(CPU_LOGS);
}

static int tracerBH_neighbors_evaluate(int target, int mode, int threadid)
{
  int cost, numnodes, *firstnode;
  data_in local, *target_data;
  data_out out, *target_result;

  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      target_data   = &local;
      target_result = &out;

      numnodes  = 1;
      firstnode = NULL;
    }
  else
    {
      target_data   = &DataGet[target];
      target_result = &DataResult[target];
      generic_get_numnodes(target, &numnodes, &firstnode);
    }

  cost = tracerBH_neighbors_force_treeevaluate(target_data, target_result, target, mode, threadid, numnodes, firstnode, MeasureCostFlag);

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  /*
     else
     DataResult[target] = out;
   */

  /* note: for imported particles, we already have the result places into DataResult[target] */

  return cost;
}

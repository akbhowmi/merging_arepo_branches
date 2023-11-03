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

static MyIDType *MinID;
static int *Head, *Len, *Next, *Tail, *MinIDTask;

/*! \brief Main routine to execute the friend of friends group finder.
 *
 *  If called with num == -1 as argument, only FOF is carried out and no group
 *  catalogs are saved to disk. If num >= 0, the code will store the
 *  group/subgroup catalogs, and bring the particles into output order.
 *  In this case, the calling routine (which is normally savepositions()) will
 *  need to free PS[] and bring the particles back into the original order,
 *  as well as reestablished the mesh.
 *
 *  \param[in] num Index of output; if negative, no output written.
 *
 *  \return void
 */
void fof_fof(int num)
{
  int i, start, lenloc, largestgroup;
  double t0, t1, cputime;

  TIMER_START(CPU_FOF);

  mpi_printf("FOF: Begin to compute FoF group catalogue...  (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));

  if(num >= 0 && RestartFlag != RESTART_FOF_SUBFIND && RestartFlag != RESTART_SNAP_CONVERSION)
    {
      /* let's discard an existing mesh - we do this here to reduce the peak memory usage, even at the price of
       * having to recreate it later */
#ifdef CREATE_SUBFOFS
      if(All.SubFOF_mode == 0)
#endif
         free_mesh();
    }

  if(RestartFlag != RESTART_SNAP_CONVERSION)
    {
      ngb_treefree();

      domain_free();
    }

#ifdef ADD_GROUP_PROPERTIES
  int ngroups_cat    = get_number_of_groups_in_catalogue(num);
  int nsubgroups_cat = get_number_of_subgroups_in_catalogue(num);
  subfind_add_grp_props_read_catalogue(num, ngroups_cat, nsubgroups_cat);
#endif

  domain_Decomposition();

#ifdef BUILD_NGBTREE_FOR_ALL_PARTICLES
  printf("\n Building NGBTREE with all particles");
  ngb_treeallocate();
  ngb_treebuild(NumPart);
  ngb_treefree();
  printf("\n Successfully Built NGBTREE with all particles");
#endif

  ngb_treeallocate();
  ngb_treebuild(NumGas);

  /* check */
  for(i = 0; i < NumPart; i++)
#ifdef GFM_RPROCESS_CHANNELS_NS_KICKS
    if((P[i].Mass == 0 && P[i].ID == 0) || (P[i].Type == 4 && P[i].Mass == 0 && STP(i).NSNS_channel == -1))
#else
    if((P[i].Mass == 0 && P[i].ID == 0) || (P[i].Type == 4 && P[i].Mass == 0))
#endif
      terminate("this should not happen");

  /* this structure will hold auxiliary information for each particle, needed only during group finding */
  PS = (struct subfind_data *)mymalloc_movable(&PS, "PS", All.MaxPart * sizeof(struct subfind_data));

  memset(PS, 0, NumPart * sizeof(struct subfind_data));

  /* First, we save the original location of the particles, in order to be able to revert to this layout later on */
  for(i = 0; i < NumPart; i++)
    {
      PS[i].OriginTask  = ThisTask;
      PS[i].OriginIndex = i;
    }

  fof_OldMaxPart    = All.MaxPart;
  fof_OldMaxPartSph = All.MaxPartSph;
#ifdef GFM
  fof_OldMaxPartStar = All.MaxPartStar;
#endif
#ifdef BLACK_HOLES
  fof_OldMaxPartBHs = All.MaxPartBHs;
#endif
#ifdef DUST_LIVE
  fof_OldMaxPartDust = All.MaxPartDust;
#endif

  LinkL = fof_get_comoving_linking_length();
#ifdef CREATE_SUBFOFS
  if(All.SubFOF_mode == 1)
     LinkL /= All.LinkingLengthReductionFactor;
#endif
  mpi_printf("FOF: Comoving linking length: %g    (presently allocated=%g MB)\n", LinkL, AllocatedBytes / (1024.0 * 1024.0));

  MinID     = (MyIDType *)mymalloc("MinID", NumPart * sizeof(MyIDType));
  MinIDTask = (int *)mymalloc("MinIDTask", NumPart * sizeof(int));

  Head = (int *)mymalloc("Head", NumPart * sizeof(int));
  Len  = (int *)mymalloc("Len", NumPart * sizeof(int));
  Next = (int *)mymalloc("Next", NumPart * sizeof(int));
  Tail = (int *)mymalloc("Tail", NumPart * sizeof(int));

#ifdef HIERARCHICAL_GRAVITY
  timebin_make_list_of_active_particles_up_to_timebin(&TimeBinsGravity, All.HighestOccupiedTimeBin);
#endif

#ifdef PREVENT_SPURIOUS_RESEEDING2
  GasNeighbor();  
#endif

#ifdef PREVENT_SEEDING_AROUND_BLACKHOLE_NEIGHBORS2
  construct_forcetree(0, 0, 0, All.HighestOccupiedTimeBin);
  BH_Neighbor();
  myfree(Father);
  myfree(Nextnode);
#ifdef BLACK_HOLES
  myfree(Tree_AuxBH_Points);
#endif
#ifdef SINKS
  myfree(Tree_AuxSinks_Points);
#endif
  myfree(Tree_Points);
  force_treefree();
#endif

  construct_forcetree(0, 0, 1, All.HighestOccupiedTimeBin); /* build tree for all particles */ 


#if defined(CREATE_SUBFOFS) && defined(BFOFS_AS_GASCLUMPS)
 if(All.SubFOF_mode == 0) 
#endif
#if defined(SUBFIND) || (defined(GFM_WINDS_VARIABLE) && (GFM_WINDS_VARIABLE == 1)) || defined(GFM_WINDS_LOCAL)
   subfind_density_hsml_guess();
#endif



  /* initialize link-lists */
  for(i = 0; i < NumPart; i++)
    {
      Head[i] = Tail[i] = i;
      Len[i]            = 1;
      Next[i]           = -1;
      MinID[i]          = P[i].ID;
      MinIDTask[i]      = ThisTask;
    }

#ifndef ADD_GROUP_PROPERTIES

  /* call routine to find primary groups */
  cputime = fof_find_groups(MinID, Head, Len, Next, Tail, MinIDTask);
  mpi_printf("FOF: group finding took = %g sec\n", cputime);

#endif

#ifdef FOF_SECONDARY_LINK_TARGET_TYPES
  myfree(Father);
  myfree(Nextnode);
#ifdef BLACK_HOLES
  myfree(Tree_AuxBH_Points);
#endif
#ifdef SINKS
  myfree(Tree_AuxSinks_Points);
#endif
  myfree(Tree_Points);

  /* now rebuild the tree with all the types selected as secondary link targets */
  construct_forcetree(0, 0, 2, All.HighestOccupiedTimeBin);
#endif

#ifdef HIERARCHICAL_GRAVITY
  timebin_make_list_of_active_particles_up_to_timebin(&TimeBinsGravity, All.HighestActiveTimeBin);
#endif

#ifndef ADD_GROUP_PROPERTIES

  /* call routine to attach secondary particles/cells to primary groups */
  cputime = fof_find_nearest_dmparticle(MinID, Head, Len, Next, Tail, MinIDTask);

  mpi_printf("FOF: attaching gas and star particles to nearest dm particles took = %g sec\n", cputime);

#endif

  /* calculate velocity dispersion etc. */

#if defined(CREATE_SUBFOFS) && defined(BFOFS_AS_GASCLUMPS)
 if(All.SubFOF_mode == 0)
  {
#endif
#if(defined(GFM_WINDS_VARIABLE) && (GFM_WINDS_VARIABLE == 1)) || defined(GFM_WINDS_LOCAL) || \
    (defined(GFM_VARIABLE_IMF) && (GFM_VARIABLE_IMF == 0))
  if(num == -2)
    {
      /* at this point, TotNgroups/PS[].GrNr are not defined yet. The following selects just gas particles
         for processing in subfind_density() */

      TotNgroups = 1;
      for(i = 0; i < NumPart; i++)
        {
          if(P[i].Type == 0)
            PS[i].GrNr = 0;
          else
            PS[i].GrNr = 1;
        }

      cputime = subfind_density(0);
      mpi_printf("GFM_WINDS_VARIABLE/GFM_WINDS_LOCAL/GFM_VARIABLE_IMF: velocity dispersion calculation took %g sec\n", cputime);
      for(i = 0; i < NumPart; i++)
        {
          if(P[i].Type == 0)
            SphP[i].w.DMVelDisp = PS[i].SubfindVelDisp;
        }
    }
#endif  // #if(defined(GFM_WINDS_VARIABLE) && (GFM_WINDS_VARIABLE == 1)) || defined(GFM_WINDS_LOCAL) || (defined(GFM_VARIABLE_IMF) &&
        // (GFM_VARIABLE_IMF == 0))

#if defined(CREATE_SUBFOFS) && defined(BFOFS_AS_GASCLUMPS)
   }
#endif


  myfree(Father);
  myfree(Nextnode);
#ifdef BLACK_HOLES
  myfree(Tree_AuxBH_Points);
#endif
#ifdef SINKS
  myfree(Tree_AuxSinks_Points);
#endif
  myfree(Tree_Points);
  force_treefree();

  myfree(Tail);
  myfree(Next);
  myfree(Len);

  t0 = second();

  FOF_PList = (struct fof_particle_list *)mymalloc_movable(&FOF_PList, "FOF_PList", NumPart * sizeof(struct fof_particle_list));

  for(i = 0; i < NumPart; i++)
    {
#ifdef ADD_GROUP_PROPERTIES
      FOF_PList[i].MinID        = P[i].MinID;
      FOF_PList[i].MinIDTask    = P[i].MinIDTask;
      FOF_PList[i].OriginalGrNr = P[i].OriginalGrNr;
#else
      FOF_PList[i].MinID     = MinID[Head[i]];
      FOF_PList[i].MinIDTask = MinIDTask[Head[i]];
#endif
      FOF_PList[i].Pindex = i;
    }

  myfree_movable(Head);
  myfree_movable(MinIDTask);
  myfree_movable(MinID);

  FOF_GList = (struct fof_group_list *)mymalloc_movable(&FOF_GList, "FOF_GList", sizeof(struct fof_group_list) * NumPart);

  fof_compile_catalogue();

  t1 = second();
  mpi_printf("FOF: compiling local group data and catalogue took = %g sec\n", timediff(t0, t1));

  MPI_Allreduce(&Ngroups, &TotNgroups, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  sumup_large_ints(1, &Nids, &TotNids);

  if(TotNgroups > 0)
    {
      int largestloc = 0;

      for(i = 0; i < NgroupsExt; i++)
        if(FOF_GList[i].LocCount + FOF_GList[i].ExtCount > largestloc)
          largestloc = FOF_GList[i].LocCount + FOF_GList[i].ExtCount;
      MPI_Allreduce(&largestloc, &largestgroup, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    }
  else
    largestgroup = 0;

  mpi_printf("FOF: Total number of FOF groups with at least %d particles: %d\n", FOF_GROUP_MIN_LEN, TotNgroups);
  mpi_printf("FOF: Largest FOF group has %d particles.\n", largestgroup);
  mpi_printf("FOF: Total number of particles in FOF groups: %lld\n", TotNids);

  t0 = second();

  MaxNgroups = 2 * imax(NgroupsExt, TotNgroups / NTask + 1);

  Group = (struct group_properties *)mymalloc_movable(&Group, "Group", sizeof(struct group_properties) * MaxNgroups);

  mpi_printf("FOF: group properties are now allocated.. (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));

#ifdef PROBABILISTIC_SEEDING
  srand(time(0));
#endif

  printf("\n Value of P[0].ID and FOF_GList[0].MinID before  %d %d at task %d \n",P[0].ID,FOF_GList[0].MinIDTask,ThisTask);
  for(i = 0, start = 0; i < NgroupsExt; i++)
    {
      while(FOF_PList[start].MinID < FOF_GList[i].MinID)
        {
          start++;
          if(start > NumPart)
            terminate("start > NumPart");
        }

      if(FOF_PList[start].MinID != FOF_GList[i].MinID)
        terminate("ID mismatch");

      for(lenloc = 0; start + lenloc < NumPart;)
        if(FOF_PList[start + lenloc].MinID == FOF_GList[i].MinID)
          lenloc++;
        else
          break;

      Group[i].MinID     = FOF_GList[i].MinID;
      Group[i].MinIDTask = FOF_GList[i].MinIDTask;

      fof_compute_group_properties(i, start, lenloc);
      start += lenloc;
    }
  fof_exchange_group_data();

  fof_finish_group_properties();
  t1 = second();
  mpi_printf("FOF: computation of group properties took = %g sec\n", timediff(t0, t1));

#if defined(CALCULATE_SPIN_STARFORMINGGAS) || (defined(GFM_BIPOLAR_WINDS) && (GFM_BIPOLAR_WINDS == 3))
  fof_spin_measurement();
#endif


#ifdef BLACK_HOLES
#ifdef SEED_HALO_ENVIRONMENT_CRITERION
if(TotNgroups > 0)
{
  fof_tag_densest_gas_cell();
//#ifdef CREATE_SUBFOFS
//  if(All.SubFOF_mode == All.SeedingModeFOForSubFOF)
//   {
//#endif
    int timebin = All.HighestSynchronizedTimeBin;
    All.Build_tree_with_only_BHs = 1;
    construct_forcetree(0, 0, 0, timebin); /* build force tree with all particles */
    gravity_tree_tracerBH_neighbors(timebin);

    myfree(Father);
    myfree(Nextnode);
#ifdef BLACK_HOLES
    myfree(Tree_AuxBH_Points);
#endif /* #ifdef BLACK_HOLES */
#ifdef SINKS
    myfree(Tree_AuxSinks_Points);
#endif /* #ifdef SINKS */
    myfree(Tree_Points);
    force_treefree();
    All.Build_tree_with_only_BHs = 0;
//#ifdef CREATE_SUBFOFS
//   }
//#endif

 Group2 = (struct group_properties *)mymalloc_movable(&Group2, "Group2", sizeof(struct group_properties) * MaxNgroups);
// printf("\n Value of P[0].ID and FOF_GList[0].MinID after  %d %d at task %d \n",P[0].ID,FOF_GList[0].MinIDTask,ThisTask);
 for(i = 0, start = 0; i < NgroupsExt; i++)
    {
      while(FOF_PList[start].MinID < FOF_GList[i].MinID)
        {
          start++;
          if(start > NumPart)
            terminate("start > NumPart");
        }

      if(FOF_PList[start].MinID != FOF_GList[i].MinID)
        terminate("ID mismatch");

      for(lenloc = 0; start + lenloc < NumPart;)
        if(FOF_PList[start + lenloc].MinID == FOF_GList[i].MinID)
          lenloc++;
        else
          break;

      Group2[i].MinID     = FOF_GList[i].MinID;
      Group2[i].MinIDTask = FOF_GList[i].MinIDTask;

      fof_compute_halo_environment(i, start, lenloc);

      start += lenloc;
    }
//   printf("\n Value of NgroupsExt and start at the start of exchange  %d, %d at task %d \n",NgroupsExt,start,ThisTask);
   fof_exchange_halo_environment_data();
   fof_finish_halo_environment();
 }
#endif




#ifndef SEED_BLACKHOLES_IN_SUBHALOS
  if(num < 0)
   {
    MPI_Barrier(MPI_COMM_WORLD);
#ifdef CREATE_SUBFOFS
    if (All.SubFOF_mode == All.SeedingModeFOForSubFOF)
     {
#endif

       fof_prepare_to_seed_black_holes();
       fof_make_black_holes();
#ifdef CREATE_SUBFOFS
     }
#endif

    MPI_Barrier(MPI_COMM_WORLD);
   }
#endif
#endif

#if(defined(GFM_WINDS_VARIABLE) && (GFM_WINDS_VARIABLE == 0)) || defined(GFM_BIPOLAR_WINDS) || defined(GFM_AGN_RADIATION) || \
    defined(MASSIVE_SEEDS_MERGER) || defined(BH_NF_RADIO)
#ifndef SEED_HALO_ENVIRONMENT_CRITERION
  if(num < 0)
#endif
    fof_assign_HostHaloMass();
#endif

  fof_assign_group_numbers();

#ifdef HCOUTPUT
  {
    if(!All.HCOutput_Halo_Initialized)
      {
        int count = 0;

        int i;
        for(i = 0; i < NumPart; i++)
          {
            if(P[i].Type == 1 && PS[i].GrNr == 0)
              {
                P[i].HCOutput_Center_Marker = 1;
                count++;
              }
            else
              P[i].HCOutput_Center_Marker = 0;
          }

        int countall;
        MPI_Allreduce(&count, &countall, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        mpi_printf("FOF: Marked %d particles for HCOUTPUT.\n", countall);
      }
  }
#endif

#ifdef AURIGA_MOVIE
  {
    if(!All.Auriga_Movie_Halo_Initialized)
      {
        int count = 0;

        int i;
        for(i = 0; i < NumPart; i++)
          {
            if(P[i].Type == 1 && PS[i].GrNr == 0)
              {
                P[i].Auriga_Movie_Center_Marker = 1;
                count++;
              }
            else
              P[i].Auriga_Movie_Center_Marker = 0;
          }

        int countall;
        MPI_Allreduce(&count, &countall, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        mpi_printf("FOF: Marked %d particles for AURIGA MOVIE.\n", countall);
      }

    /*
    int k;
    for(k = 0; k < 3; k++)
      All.Auriga_Movie_Halo_Center[k] = 0;

    int i;
    for(i = 0; i < Ngroups; i++)
      {
        if(Group[i].GrNr == 0)
          {
            for(k = 0; k < 3; k++)
              All.Auriga_Movie_Halo_Center[k] = Group[i].CM[k];
            break;
          }
      }

    MPI_Allreduce(MPI_IN_PLACE, All.Auriga_Movie_Halo_Center, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    */
  }
#endif

  mpi_printf("FOF: Finished computing FoF groups.  (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));

  myfree_movable(FOF_GList);
  myfree_movable(FOF_PList);

#ifdef SUBFIND
  if(num >= 0)
    {
      TIMER_STOP(CPU_FOF);

      subfind(num);
//#ifdef SEED_BLACKHOLES_IN_SUBHALOS
//      subfind_make_black_holes();
//#endif
#ifdef STORE_MERGERS_IN_SNAPSHOT
      mpi_printf("ATTEMPTING TO STORE MERGERS.\n");
      save_mergers(num);
#endif

      TIMER_START(CPU_FOF);
    }

#ifdef BLACK_HOLES
#ifdef SEED_BLACKHOLES_IN_SUBHALOS
  if(num < 0)
   {
    MPI_Barrier(MPI_COMM_WORLD);
    mpi_printf_task(ThisTask, "\n Trying to make subhalos \n");
    subfind(num);
    mpi_printf_task(ThisTask, "\n Trying to make black holes \n");
    subfind_make_black_holes();
    MPI_Barrier(MPI_COMM_WORLD);
   }
#endif
#endif

#else

  Nsubgroups = 0;
  TotNsubgroups = 0;
  if(num >= 0)
    {
      TIMER_STOP(CPU_FOF);
      TIMER_START(CPU_SNAPSHOT);
      fof_save_groups(num);
#ifdef STORE_MERGERS_IN_SNAPSHOT
      mpi_printf("ATTEMPTING TO STORE MERGERS.\n");
      save_mergers(num); 
#endif

      TIMER_STOP(CPU_SNAPSHOT);
      TIMER_START(CPU_FOF);
    }
#endif

#ifdef SEED_HALO_ENVIRONMENT_CRITERION
if(TotNgroups > 0)
  myfree_movable(Group2);
#endif

  myfree_movable(Group);

  mpi_printf("FOF: All FOF related work finished.  (presently allocated=%g MB)\n", AllocatedBytes / (1024.0 * 1024.0));

#ifndef FOF_STOREIDS
  if(num >= 0)
    {
      TIMER_STOP(CPU_FOF);
      TIMER_START(CPU_SNAPSHOT);

      /* now distribute the particles into output order */
#ifdef CREATE_SUBFOFS
      if(All.SubFOF_mode == All.ParticleFOFdatacorrespondence)
        {
#endif
          t0 = second();
          fof_prepare_output_order();
          fof_subfind_exchange(MPI_COMM_WORLD); /* distribute particles such that FOF groups will appear in coherent way in snapshot files */
          t1 = second();
          mpi_printf("FOF: preparing output order of particles took %g sec\n", timediff(t0, t1));
#ifdef CREATE_SUBFOFS
        }
      if(All.SubFOF_mode == 0)
         myfree(PS);
#endif
      TIMER_STOP(CPU_SNAPSHOT);
      TIMER_START(CPU_FOF);
    } 
  else
    myfree(PS);
#else
  myfree(PS);
#endif

  TIMER_STOP(CPU_FOF);
}

/*! \brief Sorts groups by the desired output order.
 *
 *  \return void
 */
void fof_prepare_output_order(void)
{
  int i, off, ntype[NTYPES];

  struct data_aux_sort *aux_sort = (struct data_aux_sort *)mymalloc("aux_sort", sizeof(struct data_aux_sort) * NumPart);

  for(i = 0; i < NTYPES; i++)
    ntype[i] = 0;

  for(i = 0; i < NumPart; i++)
    {
      aux_sort[i].OriginTask  = ThisTask;
      aux_sort[i].OriginIndex = i;
      aux_sort[i].GrNr        = PS[i].GrNr;
#ifdef SUBFIND
      aux_sort[i].SubNr            = PS[i].SubNr;
      aux_sort[i].DM_BindingEnergy = PS[i].BindingEnergy;
#endif
      aux_sort[i].Type = P[i].Type;
#ifdef FOF_FUZZ_SORT_BY_NEAREST_GROUP
      aux_sort[i].key = PS[i].GroupNr;
#else
      aux_sort[i].ID = P[i].ID;
#endif
#if defined(RECOMPUTE_POTENTIAL_IN_SNAPSHOT) || defined(CALCULATE_QUANTITIES_IN_POSTPROCESS)
      aux_sort[i].FileOrder = P[i].FileOrder;
#endif

      ntype[P[i].Type]++;
    }

  qsort(aux_sort, NumPart, sizeof(struct data_aux_sort), fof_compare_aux_sort_Type);

  if(RestartFlag == RESTART_RECALC_POTENTIAL || RestartFlag == RESTART_CALC_ADDITIONAL)
    {
#if defined(RECOMPUTE_POTENTIAL_IN_SNAPSHOT) || defined(CALCULATE_QUANTITIES_IN_POSTPROCESS)
      for(i = 0, off = 0; i < NTYPES; off += ntype[i], i++)
        parallel_sort(aux_sort + off, ntype[i], sizeof(struct data_aux_sort), fof_compare_aux_sort_FileOrder);
#endif
    }
  else
    {
      for(i = 0, off = 0; i < NTYPES; off += ntype[i], i++)
        parallel_sort(aux_sort + off, ntype[i], sizeof(struct data_aux_sort), fof_compare_aux_sort_GrNr);
    }

  for(i = 0; i < NumPart; i++)
    {
      aux_sort[i].TargetTask  = ThisTask;
      aux_sort[i].TargetIndex = i;
    }

  /* now bring back into starting order */
  parallel_sort(aux_sort, NumPart, sizeof(struct data_aux_sort), fof_compare_aux_sort_OriginTask_OriginIndex);

  for(i = 0; i < NumPart; i++)
    {
      PS[i].TargetTask  = aux_sort[i].TargetTask;
      PS[i].TargetIndex = aux_sort[i].TargetIndex;
    }

  myfree(aux_sort);
}

/*! \brief Calculate linking length based on mean particle separation.
 *
 *  \return Linking length.
 */
double fof_get_comoving_linking_length(void)
{
  int i, ndm;
  long long ndmtot;
  double mass, masstot, rhodm;
#if defined(CREATE_SUBFOFS) && defined(BFOFS_AS_GASCLUMPS)                                         
  int Link_Types;
#endif

  for(i = 0, ndm = 0, mass = 0; i < NumPart; i++)
   {
#if defined(CREATE_SUBFOFS) && defined(BFOFS_AS_GASCLUMPS)
    if(All.SubFOF_mode == 1)
       Link_Types = GAS_CLUMP_PRIMARY_LINK_TYPES;
    else if(All.SubFOF_mode == 0)
       Link_Types = FOF_PRIMARY_LINK_TYPES;
    if(((1 << P[i].Type) & (Link_Types)))
#else
    if(((1 << P[i].Type) & (FOF_PRIMARY_LINK_TYPES)))
#endif
      {
        ndm++;
        mass += P[i].Mass;
      }
   }
  sumup_large_ints(1, &ndm, &ndmtot);
  MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  rhodm = (All.Omega0 - All.OmegaBaryon) * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G);

  return FOF_LINKLENGTH * pow(masstot / ndmtot / rhodm, 1.0 / 3);
}

/*! \brief Compiles the group catalogue.
 *
 *  Combines results from all tasks.
 *
 *  \return void
 */
void fof_compile_catalogue(void)
{
  int i, j, start, nimport, ngrp, recvTask;
  struct fof_group_list *get_FOF_GList;
  int MinParticlesPerGroup;

  /* sort according to MinID */
  mysort(FOF_PList, NumPart, sizeof(struct fof_particle_list), fof_compare_FOF_PList_MinID);

  for(i = 0; i < NumPart; i++)
    {
      FOF_GList[i].MinID     = FOF_PList[i].MinID;
      FOF_GList[i].MinIDTask = FOF_PList[i].MinIDTask;
      if(FOF_GList[i].MinIDTask == ThisTask)
        {
          FOF_GList[i].LocCount = 1;
          FOF_GList[i].ExtCount = 0;
        }
      else
        {
          FOF_GList[i].LocCount = 0;
          FOF_GList[i].ExtCount = 1;
        }

#ifdef ADD_GROUP_PROPERTIES
      FOF_GList[i].OriginalGrNr = FOF_PList[i].OriginalGrNr;
#endif

#ifdef TRACER_PARTICLE
      /* count tracer particles independently so we can enforce FOF_GROUP_MIN_LEN against only non-tracers */
      if(P[FOF_PList[i].Pindex].Type == TRACER_PARTICLE)
        {
          if(FOF_GList[i].MinIDTask == ThisTask)
            {
              FOF_GList[i].LocTrCount = 1;
              FOF_GList[i].ExtTrCount = 0;
            }
          else
            {
              FOF_GList[i].LocTrCount = 0;
              FOF_GList[i].ExtTrCount = 1;
            }
        }
      else
        {
          FOF_GList[i].LocTrCount = 0;
          FOF_GList[i].ExtTrCount = 0;
        }
#endif
    }

  /* eliminate duplicates in FOF_GList with respect to MinID */

  if(NumPart)
    NgroupsExt = 1;
  else
    NgroupsExt = 0;

  for(i = 1, start = 0; i < NumPart; i++)
    {
      if(FOF_GList[i].MinID == FOF_GList[start].MinID)
        {
          FOF_GList[start].LocCount += FOF_GList[i].LocCount;
          FOF_GList[start].ExtCount += FOF_GList[i].ExtCount;
#ifdef TRACER_PARTICLE
          FOF_GList[start].LocTrCount += FOF_GList[i].LocTrCount;
          FOF_GList[start].ExtTrCount += FOF_GList[i].ExtTrCount;
#endif
        }
      else
        {
          start            = NgroupsExt;
          FOF_GList[start] = FOF_GList[i];
          NgroupsExt++;
        }
    }

  /* sort the remaining ones according to task */
  mysort(FOF_GList, NgroupsExt, sizeof(struct fof_group_list), fof_compare_FOF_GList_MinIDTask);

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

  get_FOF_GList = (struct fof_group_list *)mymalloc("get_FOF_GList", nimport * sizeof(struct fof_group_list));

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the group info */
              MPI_Sendrecv(&FOF_GList[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct fof_group_list), MPI_BYTE, recvTask,
                           TAG_DENS_A, &get_FOF_GList[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct fof_group_list),
                           MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  for(i = 0; i < nimport; i++)
    get_FOF_GList[i].MinIDTask = i;

  /* sort the groups according to MinID */
  mysort(FOF_GList, NgroupsExt, sizeof(struct fof_group_list), fof_compare_FOF_GList_MinID);
  mysort(get_FOF_GList, nimport, sizeof(struct fof_group_list), fof_compare_FOF_GList_MinID);

  /* merge the imported ones with the local ones */
  for(i = 0, start = 0; i < nimport; i++)
    {
      while(FOF_GList[start].MinID < get_FOF_GList[i].MinID)
        {
          start++;
          if(start >= NgroupsExt)
            terminate("start >= NgroupsExt");
        }

      if(get_FOF_GList[i].LocCount != 0)
        terminate("get_FOF_GList[i].LocCount != 0");

      if(FOF_GList[start].MinIDTask != ThisTask)
        terminate("FOF_GList[start].MinIDTask != ThisTask");

      if(FOF_GList[start].MinID != get_FOF_GList[i].MinID)
        terminate("FOF_GList[start].MinID != get_FOF_GList[i].MinID start=%d i=%d FOF_GList[start].MinID=%" MYIDTYPE_PRI
                  " get_FOF_GList[i].MinID=%" MYIDTYPE_PRI "\n",
                  start, i, FOF_GList[start].MinID, get_FOF_GList[i].MinID);

      FOF_GList[start].ExtCount += get_FOF_GList[i].ExtCount;
#ifdef TRACER_PARTICLE
      FOF_GList[start].ExtTrCount += get_FOF_GList[i].ExtTrCount;

      if(FOF_GList[start].LocTrCount > FOF_GList[start].LocCount || FOF_GList[start].ExtTrCount > FOF_GList[start].ExtCount)
        terminate("Error: FOF_GList tracer counts exceed total counts!");
#endif
    }

  /* copy the size information back into the list, to inform the others */
  for(i = 0, start = 0; i < nimport; i++)
    {
      while(FOF_GList[start].MinID < get_FOF_GList[i].MinID)
        {
          start++;
          if(start >= NgroupsExt)
            terminate("start >= NgroupsExt");
        }

      get_FOF_GList[i].ExtCount = FOF_GList[start].ExtCount;
      get_FOF_GList[i].LocCount = FOF_GList[start].LocCount;
#ifdef TRACER_PARTICLE
      get_FOF_GList[i].ExtTrCount = FOF_GList[start].ExtTrCount;
      get_FOF_GList[i].LocTrCount = FOF_GList[start].LocTrCount;
#endif
    }

  /* sort the imported/exported list according to MinIDTask */
  mysort(get_FOF_GList, nimport, sizeof(struct fof_group_list), fof_compare_FOF_GList_MinIDTask);
  mysort(FOF_GList, NgroupsExt, sizeof(struct fof_group_list), fof_compare_FOF_GList_MinIDTask);

  for(i = 0; i < nimport; i++)
    get_FOF_GList[i].MinIDTask = ThisTask;

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the group info */
              MPI_Sendrecv(&get_FOF_GList[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct fof_group_list), MPI_BYTE,
                           recvTask, TAG_DENS_A, &FOF_GList[Send_offset[recvTask]],
                           Send_count[recvTask] * sizeof(struct fof_group_list), MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD,
                           MPI_STATUS_IGNORE);
            }
        }
    }

  myfree(get_FOF_GList);

  /* eliminate all groups that are too small, and count local groups */

  
  MinParticlesPerGroup = FOF_GROUP_MIN_LEN;
#if defined(CREATE_SUBFOFS) && defined(BFOFS_AS_GASCLUMPS)
  if (All.SubFOF_mode == 1)
       MinParticlesPerGroup = 3;  
#endif

  for(i = 0, Ngroups = 0, Nids = 0; i < NgroupsExt; i++)
    {
        if(FOF_GList[i].LocCount + FOF_GList[i].ExtCount
#ifdef TRACER_PARTICLE
             - FOF_GList[i].LocTrCount - FOF_GList[i].ExtTrCount
#endif
         < MinParticlesPerGroup)
        {
          FOF_GList[i] = FOF_GList[NgroupsExt - 1];
          NgroupsExt--;
          i--;
        }

      else
        {
          if(FOF_GList[i].MinIDTask == ThisTask)
            {
              Ngroups++;
              Nids += FOF_GList[i].LocCount + FOF_GList[i].ExtCount;
            }
        }
    }

  /* sort the group list according to MinID */
  mysort(FOF_GList, NgroupsExt, sizeof(struct fof_group_list), fof_compare_FOF_GList_MinID);
}

/*! \brief Assigns each group a global group number.
 *
 *  \return void
 */
void fof_assign_group_numbers(void)
{
  int i, j, ngr, start, lenloc;
  long long totNids;
  double t0, t1;

  mpi_printf("FOF: start assigning group numbers\n");

  t0 = second();

  /* assign group numbers (at this point, both Group and FOF_GList are sorted by MinID) */
  for(i = 0; i < NgroupsExt; i++)
    {
      FOF_GList[i].LocCount += FOF_GList[i].ExtCount; /* total length */
      FOF_GList[i].ExtCount = ThisTask;               /* original task */
    }

  parallel_sort(FOF_GList, NgroupsExt, sizeof(struct fof_group_list), fof_compare_FOF_GList_LocCountTaskDiffMinID);

  for(i = 0, ngr = 0; i < NgroupsExt; i++)
    {
      if(FOF_GList[i].ExtCount == FOF_GList[i].MinIDTask)
        ngr++;

      FOF_GList[i].GrNr = ngr - 1;
    }

  MPI_Allgather(&ngr, 1, MPI_INT, Send_count, 1, MPI_INT, MPI_COMM_WORLD);

  /* count how many groups there are on earlier CPUs */
  long long ngr_sum;
  for(j = 0, ngr_sum = 0; j < ThisTask; j++)
    ngr_sum += Send_count[j];

  for(i = 0; i < NgroupsExt; i++)
    FOF_GList[i].GrNr += ngr_sum;

  sumup_large_ints(1, &ngr, &ngr_sum);
  if(ngr_sum != TotNgroups)
    terminate("inconsistency (ngr_sum=%lld)", ngr_sum);

  /* bring the group list back into the original order */
  parallel_sort(FOF_GList, NgroupsExt, sizeof(struct fof_group_list), fof_compare_FOF_GList_ExtCountMinID);

  /* Assign the group numbers to the group properties array */
  for(i = 0, start = 0; i < Ngroups; i++)
    {
      while(FOF_GList[start].MinID < Group[i].MinID)
        {
          start++;
          if(start >= NgroupsExt)
            terminate("start >= NgroupsExt");
        }
      Group[i].GrNr = FOF_GList[start].GrNr;
#ifdef SEED_HALO_ENVIRONMENT_CRITERION
      Group2[i].GrNr = FOF_GList[start].GrNr;      
#endif

    }

  /* sort the groups according to group-number */
  parallel_sort(Group, Ngroups, sizeof(struct group_properties), fof_compare_Group_GrNr);

#ifdef SEED_HALO_ENVIRONMENT_CRITERION
  parallel_sort(Group2, Ngroups, sizeof(struct group_properties), fof_compare_Group_GrNr);
#endif


  for(i = 0; i < NumPart; i++)
    PS[i].GrNr = TotNgroups + 1; /* this marks all particles that are not in any group */

  for(i = 0, start = 0, Nids = 0; i < NgroupsExt; i++)
    {
      while(FOF_PList[start].MinID < FOF_GList[i].MinID)
        {
          start++;
          if(start > NumPart)
            terminate("start > NumPart");
        }

      if(FOF_PList[start].MinID != FOF_GList[i].MinID)
        terminate("FOF_PList[start=%d].MinID=%" MYIDTYPE_PRI " != FOF_GList[i=%d].MinID=%" MYIDTYPE_PRI, start, FOF_PList[start].MinID,
                  i, FOF_GList[i].MinID);

      for(lenloc = 0; start + lenloc < NumPart;)
        if(FOF_PList[start + lenloc].MinID == FOF_GList[i].MinID)
          {
            PS[FOF_PList[start + lenloc].Pindex].GrNr = FOF_GList[i].GrNr;
            Nids++;
            lenloc++;
          }
        else
          break;

      start += lenloc;
    }

  sumup_large_ints(1, &Nids, &totNids);

  if(totNids != TotNids)
    terminate("Task=%d Nids=%d totNids=%lld TotNids=%lld\n", ThisTask, Nids, totNids, TotNids);

  t1 = second();

  mpi_printf("FOF: Assigning of group numbers took = %g sec\n", timediff(t0, t1));
}

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
void fof_compute_group_properties(int gr, int start, int len)
{
  int j, k, index, type, start_index = FOF_PList[start].Pindex;
  double xyz[3];

  Group[gr].Len  = 0;
  double gr_Mass = 0;
#ifdef USE_SFR
  double gr_Sfr = 0;
#endif
#ifdef PROBABILISTIC_SEEDING
//  Group[gr].PlaceSeedIfCriterionSatisfied = 0; 
  Group[gr].RandomFractionForSeed = get_random_fraction(3); 
  mpi_printf("For group %d, generated a random seed of %.3f", gr, Group[gr].RandomFractionForSeed);
#endif

#ifdef GFM_STELLAR_EVOLUTION
  double *GroupMassMetallicity = NULL, *GroupMassMetals = NULL; /* FIXME */
  double gr_GasMassMetallicity = 0, gr_StellarMassMetallicity = 0;
  double gr_GasMassMetals[GFM_N_CHEM_ELEMENTS], gr_StellarMassMetals[GFM_N_CHEM_ELEMENTS];
  for(k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
    {
      gr_GasMassMetals[k]     = 0;
      gr_StellarMassMetals[k] = 0;
    }
#ifdef GFM_DUST
  double *GroupMassDustMetallicity = NULL;
  double gr_GasMassDustMetallicity = 0;
  int chan                         = 0;
#endif
#endif
#ifdef BLACK_HOLES
  double gr_BH_Mass = 0, gr_BH_Mdot = 0, gr_MaxDens = 0;

#ifdef CHECK_FOR_ENOUGH_GAS_MASS_IN_DCBH_FORMING_POCKETS
  MyFloat gr_MaxNeighboringDCBHFormingGasMass = 0;
#endif

#ifdef OUTPUT_LOG_FILES_FOR_SEEDING
  MyFloat gr_Metallicity_maxdens = 0, gr_Sfr_maxdens = 0, gr_Mass_maxdens = 0;
  MyIDType gr_ID_maxdens = 0;
  MyFloat gr_StarFormingGasMass = 0;
  MyFloat gr_StarFormingGasMassMetallicity = 0;
  MyFloat gr_StarFormingMetalFreeGasMass = 0;
#ifdef CALCULATE_LYMAN_WERNER_INTENSITY_LOCAL_STARFORMINGGAS
  MyFloat gr_StarFormingMetalFreeLymanWernerGasMass = 0;
  MyFloat gr_LymanWernerGasMass = 0;
  MyFloat gr_MaxLymanWernerIntensityInDenseMetalPoorGas = 0;
#endif
#endif

#ifdef SEED_BASED_ON_PROBABLISTIC_HALO_PROPERTIES
#ifdef PROBABILISTIC_SEED_MASS_HALO_MASS_RATIO_CRITERION
  MyFloat gr_RandomMinHaloMassForSeeding_maxdens = 0;
#endif
#endif

#ifdef EVOLVING_SEEDING_PROBABILITY
 MyFloat gr_SecondRandomNumberForSeeding_maxdens = 0;
 MyFloat gr_SecondRandomNumberForSeeding_average = 0;
#endif

#ifdef SEED_HALO_ENVIRONMENT_CRITERION2
 MyFloat gr_ThirdRandomNumberForSeeding_maxdens = 0;
#endif

#ifdef UNIFORM_SEEDMASS_DISTRIBUTION
 MyFloat gr_DrawnSeedBlackHoleMass_maxdens = 0;
#endif

#ifdef PREVENT_SPURIOUS_RESEEDING
  MyFloat gr_SeedMass_maxdens = 0;
#endif



#ifdef PREVENT_SPURIOUS_RESEEDING
  MyFloat gr_TotalGasSeedMass=0;
#endif

  int gr_CouldHaveBeenABlackHole_sum = 0;


#ifdef PREVENT_SEEDING_AROUND_BLACKHOLE_NEIGHBORS2
  int gr_BHNeighborExists_maxdens=0;
#endif

#ifdef CALCULATE_LYMAN_WERNER_INTENSITY_LOCAL_SOURCES
  MyFloat gr_LymanWernerIntensityLocalSources_maxdens_type2 = 0;
  MyFloat gr_LymanWernerIntensityLocalSources_maxdens_type3 = 0;
#endif


#ifdef CALCULATE_LYMAN_WERNER_INTENSITY_LOCAL_STARFORMINGGAS
  MyFloat gr_LymanWernerIntensityLocalStarFormingGas_maxdens_type2 = 0;
  MyFloat gr_LymanWernerIntensityLocalStarFormingGas_maxdens_type3 = 0;
#endif


#ifdef CALCULATE_LYMAN_WERNER_INTENSITY_ALL_SOURCES
  MyFloat gr_LymanWernerIntensityAllSources_maxdens_type2 = 0;
  MyFloat gr_LymanWernerIntensityAllSources_maxdens_type3 = 0;
#endif

#ifdef PREVENT_SPURIOUS_RESEEDING2
  int gr_NeighborOfBlackhole_maxdens=0;
#endif

  Group[gr].index_maxdens = Group[gr].task_maxdens = -1;
#ifdef SEED_HALO_ENVIRONMENT_CRITERION
  Group[gr].index_MinPot = Group[gr].task_MinPot = -1;
#endif

#ifdef BH_NF_RADIO
  double gr_Min_BH_Potential    = MAX_FLOAT_NUMBER;
  Group[gr].ID_Min_BH_Potential = 0;
  double gr_XrayLum = 0, gr_RadioLum = 0;
#endif
#endif
#ifdef GFM_WINDS
  double gr_WindMass = 0;
#endif
  double gr_CM[3], gr_Vel[3];
#ifdef CONSTRUCT_FOF_NGBTREE
  double gr_CM_unwrapped[3];
#endif

#if defined(GFM_BIPOLAR_WINDS) && GFM_BIPOLAR_WINDS != 3
  double gr_GravAcc[3];
#endif
  for(k = 0; k < 3; k++)
    {
      gr_CM[k]  = 0;
      gr_Vel[k] = 0;
#ifdef CONSTRUCT_FOF_NGBTREE
      gr_CM_unwrapped[k]  = 0;
#endif

#if defined(GFM_BIPOLAR_WINDS) && GFM_BIPOLAR_WINDS != 3
      gr_GravAcc[k] = 0;
#endif
      Group[gr].FirstPos[k] = P[start_index].Pos[k];
    }

#if defined(CALCULATE_SPIN_STARFORMINGGAS) || (defined(GFM_BIPOLAR_WINDS) && (GFM_BIPOLAR_WINDS == 3))
  double gr_MinPotential = MAX_FLOAT_NUMBER;
#endif

#ifdef SEED_HALO_ENVIRONMENT_CRITERION
  double gr_DMMinPotential = MAX_FLOAT_NUMBER;
#endif


  double gr_MassType[NTYPES];
  for(k = 0; k < NTYPES; k++)
    {
      Group[gr].LenType[k] = 0;
      gr_MassType[k]       = 0;
    }
  /* calculate */
  for(k = 0; k < len; k++)
    {
      index = FOF_PList[start + k].Pindex;

      Group[gr].Len++;
      gr_Mass += P[index].Mass;
      type = P[index].Type;

#ifdef TRACER_MC
      Group[gr].LenType[TRACER_MC] += get_number_of_tracers(index);
#endif

#ifdef GFM_WINDS_SAVE_PARTTYPE
      /* count wind as new particle type for LenType since we save these separately */
      if(P[index].Type == PTYPE_STARS && STP(index).BirthTime <= 0)
        Group[gr].LenType[GFM_WINDS_SAVE_PARTTYPE]++;
      else
#endif
        Group[gr].LenType[type]++;

#ifdef GFM_WINDS
      /* count wind as gas for mass, but not for LenType since we use this to construct offset tables */
      if(P[index].Type == PTYPE_STARS && STP(index).BirthTime <= 0)
        type = 0;
#endif
      gr_MassType[type] += P[index].Mass;

#ifdef BH_NF_RADIO
      if(P[index].Type == PTYPE_GAS)
        gr_XrayLum += get_cooling_luminosity(index);
#endif

#ifdef USE_SFR
      if(P[index].Type == PTYPE_GAS)
        gr_Sfr += SphP[index].Sfr;
#endif
#ifdef GFM_STELLAR_EVOLUTION
      if(P[index].Type == PTYPE_GAS)
        {
          GroupMassMetallicity = &(gr_GasMassMetallicity);
          GroupMassMetals      = gr_GasMassMetals;
        }
      if(P[index].Type == PTYPE_STARS)
        {
          GroupMassMetallicity = &(gr_StellarMassMetallicity);
          GroupMassMetals      = gr_StellarMassMetals;
        }
#ifdef GFM_WINDS
      if(P[index].Type == PTYPE_STARS && STP(index).BirthTime <= 0)
        {
          GroupMassMetallicity = &(gr_GasMassMetallicity);
          GroupMassMetals      = gr_GasMassMetals;
        }
#endif
      if(P[index].Type == PTYPE_GAS)
        {
#if GFM_STELLAR_EVOLUTION == 1
          *GroupMassMetallicity += SphP[index].Metallicity * P[index].Mass;
#else
          *GroupMassMetallicity += SphP[index].MassMetallicity;
#endif
          for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
            {
#if GFM_STELLAR_EVOLUTION == 1
              GroupMassMetals[j] += SphP[index].MetalsFraction[j] * P[index].Mass;
#else
              GroupMassMetals[j] += SphP[index].MassMetals[j];
#endif
            }
        }
      if(P[index].Type == PTYPE_STARS)
        {
          *GroupMassMetallicity += STP(index).Metallicity * P[index].Mass;
          for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
            {
              GroupMassMetals[j] += STP(index).MassMetals[j];
            }
        }
#endif
#if defined(GFM_STELLAR_EVOLUTION) && defined(GFM_DUST)
      if(P[index].Type == PTYPE_GAS)
        {
          GroupMassDustMetallicity = &(gr_GasMassDustMetallicity);
        }
#ifdef GFM_WINDS
      if(P[index].Type == PTYPE_STARS && STP(index).BirthTime <= 0)
        {
          GroupMassDustMetallicity = &(gr_GasMassDustMetallicity);
        }
#endif
      if(P[index].Type == PTYPE_GAS)
        {
          for(chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
            for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
              *GroupMassDustMetallicity += SphP[index].MetalsDustFraction[chan][j] * P[index].Mass;
        }
      if(P[index].Type == PTYPE_STARS && STP(index).BirthTime <= 0)
        {
          for(chan = 0; chan < GFM_DUST_N_CHANNELS; chan++)
            for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
              *GroupMassDustMetallicity += STP(index).MassMetals[j] * STP(index).InitialDustFractions[chan][j];
        }
#endif

#ifdef BLACK_HOLES
#ifdef PREVENT_SPURIOUS_RESEEDING
      if(P[index].Type == 0)
        {
          gr_TotalGasSeedMass += SphP[index].SeedMass; 
        }
      else if(P[index].Type == 4)
        {
          gr_TotalGasSeedMass += STP(index).SeedMass;
        }
#endif     
      if(P[index].Type == 0)
         gr_CouldHaveBeenABlackHole_sum += SphP[index].CouldHaveBeenABlackHole;


#ifdef CALCULATE_LYMAN_WERNER_INTENSITY_LOCAL_STARFORMINGGAS
      if(P[index].Type == 0)
        {
          if ((SphP[index].StarFormingGasLymanWernerIntensity_type2 + SphP[index].StarFormingGasLymanWernerIntensity_type3) > All.MinLymanWernerFluxForNewSeed)        
            {
                 gr_LymanWernerGasMass += P[index].Mass;
            }
        }
#endif

#ifdef EVOLVING_SEEDING_PROBABILITY
      if(P[index].Type == 0)
          gr_SecondRandomNumberForSeeding_average += P[index].Mass * SphP[index].SecondRandomNumberForSeeding;    
#endif

#ifdef OUTPUT_LOG_FILES_FOR_SEEDING
#ifdef SUPPRESS_STARFORMATION_ABOVE_CRITICAL_LYMANWERNERFLUX
      if((P[index].Type == 0) && (SphP[index].GasIsDense == 1))
#else
      if((P[index].Type == 0) && (SphP[index].Sfr > 0))
#endif
      {
         gr_StarFormingGasMass += P[index].Mass;
         gr_StarFormingGasMassMetallicity += P[index].Mass * SphP[index].Metallicity;
         if (SphP[index].MassMetallicity < All.MaxMetallicityForAssumingMetalFree * P[index].Mass * GFM_SOLAR_METALLICITY)
           {
#ifdef REFINEMENT_HIGH_RES_GAS
             gr_StarFormingMetalFreeGasMass += SphP[index].HighResMass;
#else
             gr_StarFormingMetalFreeGasMass += P[index].Mass;
#endif
#ifdef CALCULATE_LYMAN_WERNER_INTENSITY_LOCAL_STARFORMINGGAS
             if ((SphP[index].StarFormingGasLymanWernerIntensity_type2 + SphP[index].StarFormingGasLymanWernerIntensity_type3) > All.MinLymanWernerFluxForNewSeed)
#ifdef REFINEMENT_HIGH_RES_GAS 
                 gr_StarFormingMetalFreeLymanWernerGasMass += SphP[index].HighResMass;
#else
                 gr_StarFormingMetalFreeLymanWernerGasMass += P[index].Mass;
#endif
#endif
           }
      }
#endif

      if(P[index].Type == 5)
        {
          gr_BH_Mdot += BPP(index).BH_Mdot;
          gr_BH_Mass += BPP(index).BH_Mass;
#ifdef BH_NF_RADIO
          if(P[index].Potential < gr_Min_BH_Potential)
            {
              gr_Min_BH_Potential           = P[index].Potential;
              Group[gr].ID_Min_BH_Potential = P[index].ID;
            }
#endif
        }
      if(P[index].Type == PTYPE_GAS)
        {
#ifdef CALCULATE_LYMAN_WERNER_INTENSITY_LOCAL_STARFORMINGGAS
#ifdef SUPPRESS_STARFORMATION_ABOVE_CRITICAL_LYMANWERNERFLUX
          if(SphP[index].GasIsDense == 1)
#else
          if(SphP[index].Sfr > 0)
#endif
            if (SphP[index].MassMetallicity < All.MaxMetallicityForAssumingMetalFree * P[index].Mass * GFM_SOLAR_METALLICITY)
               if(SphP[index].StarFormingGasLymanWernerIntensity_type2 + SphP[index].StarFormingGasLymanWernerIntensity_type3 > gr_MaxLymanWernerIntensityInDenseMetalPoorGas)
                  gr_MaxLymanWernerIntensityInDenseMetalPoorGas = SphP[index].StarFormingGasLymanWernerIntensity_type2 + SphP[index].StarFormingGasLymanWernerIntensity_type3;   
#endif

#ifdef CHECK_FOR_ENOUGH_GAS_MASS_IN_DCBH_FORMING_POCKETS
          if(SphP[index].NeighboringDCBHFormingGasMass > gr_MaxNeighboringDCBHFormingGasMass)
              gr_MaxNeighboringDCBHFormingGasMass = SphP[index].NeighboringDCBHFormingGasMass;
#endif
#ifdef SEED_BHS_FROM_ELIGIBLE_GAS_CELLS
#if defined(SEED_STARFORMINGMETALFREEGASMASS_CRITERION) || defined(SEED_STARFORMINGMETALFREELYMANWERNERGASMASS_CRITERION)
#ifdef SUPPRESS_STARFORMATION_ABOVE_CRITICAL_LYMANWERNERFLUX
             if(SphP[index].GasIsDense == 1)
#else
             if(SphP[index].Sfr > 0)
#endif
             if (SphP[index].MassMetallicity < All.MaxMetallicityForAssumingMetalFree * P[index].Mass * GFM_SOLAR_METALLICITY)
#ifdef SEED_STARFORMINGMETALFREELYMANWERNERGASMASS_CRITERION
             if ((SphP[index].StarFormingGasLymanWernerIntensity_type2 + SphP[index].StarFormingGasLymanWernerIntensity_type3) > All.MinLymanWernerFluxForNewSeed)          
#endif
#endif
#endif
             if(SphP[index].Density > gr_MaxDens)
            {
              gr_MaxDens              = SphP[index].Density;
#ifdef OUTPUT_LOG_FILES_FOR_SEEDING
              gr_Metallicity_maxdens  = SphP[index].Metallicity;
              gr_Sfr_maxdens =  SphP[index].Sfr;
              gr_Mass_maxdens = P[index].Mass;
              gr_ID_maxdens = P[index].ID;
#endif

#ifdef SEED_BASED_ON_PROBABLISTIC_HALO_PROPERTIES
#ifdef PROBABILISTIC_SEED_MASS_HALO_MASS_RATIO_CRITERION
              if (SphP[index].RandomMinHaloMassForSeeding < 0)
                  SphP[index].RandomMinHaloMassForSeeding = - SphP[index].RandomMinHaloMassForSeeding;
              gr_RandomMinHaloMassForSeeding_maxdens  = SphP[index].RandomMinHaloMassForSeeding;
#endif
#endif

#ifdef EVOLVING_SEEDING_PROBABILITY
             gr_SecondRandomNumberForSeeding_maxdens  = SphP[index].SecondRandomNumberForSeeding;
#endif

#ifdef SEED_HALO_ENVIRONMENT_CRITERION2
             gr_ThirdRandomNumberForSeeding_maxdens  = SphP[index].ThirdRandomNumberForSeeding;
#endif

#ifdef UNIFORM_SEEDMASS_DISTRIBUTION
             gr_DrawnSeedBlackHoleMass_maxdens  = SphP[index].DrawnSeedBlackHoleMass;
#endif


#ifdef PREVENT_SPURIOUS_RESEEDING
              gr_SeedMass_maxdens  = SphP[index].SeedMass;
#endif

#ifdef CALCULATE_LYMAN_WERNER_INTENSITY_LOCAL_SOURCES
              gr_LymanWernerIntensityLocalSources_maxdens_type2  = SphP[index].StellarLymanWernerIntensity_type2; 
              gr_LymanWernerIntensityLocalSources_maxdens_type3  = SphP[index].StellarLymanWernerIntensity_type3;
#endif

#ifdef CALCULATE_LYMAN_WERNER_INTENSITY_LOCAL_STARFORMINGGAS
              gr_LymanWernerIntensityLocalStarFormingGas_maxdens_type2  = SphP[index].StarFormingGasLymanWernerIntensity_type2;
              gr_LymanWernerIntensityLocalStarFormingGas_maxdens_type3  = SphP[index].StarFormingGasLymanWernerIntensity_type3;
#endif

#ifdef CALCULATE_LYMAN_WERNER_INTENSITY_ALL_SOURCES
      	      gr_LymanWernerIntensityAllSources_maxdens_type2  = P[index].StellarAllLymanWernerIntensity_type2;
              gr_LymanWernerIntensityAllSources_maxdens_type3  = P[index].StellarAllLymanWernerIntensity_type3;
#endif

#ifdef PREVENT_SEEDING_AROUND_BLACKHOLE_NEIGHBORS2
              gr_BHNeighborExists_maxdens  = SphP[index].BHNeighborExists;
#endif
#ifdef PREVENT_SPURIOUS_RESEEDING2
              gr_NeighborOfBlackhole_maxdens = SphP[index].NeighborOfBlackhole;
#endif
              Group[gr].index_maxdens = index;
              Group[gr].task_maxdens  = ThisTask;
            }
        }
#endif

#ifdef GFM_WINDS
      if(P[index].Type == PTYPE_STARS && STP(index).BirthTime <= 0)
        gr_WindMass += P[index].Mass;
#endif

#if defined(CALCULATE_SPIN_STARFORMINGGAS) || (defined(GFM_BIPOLAR_WINDS) && (GFM_BIPOLAR_WINDS == 3))
      if(P[index].Potential < gr_MinPotential)
        {
          gr_MinPotential               = P[index].Potential;
          Group[gr].Pos_MinPotential[0] = P[index].Pos[0];
          Group[gr].Pos_MinPotential[1] = P[index].Pos[1];
          Group[gr].Pos_MinPotential[2] = P[index].Pos[2];
        }
#endif

#ifdef SEED_HALO_ENVIRONMENT_CRITERION
      if(P[index].Type == 1)
          if(P[index].Potential < gr_DMMinPotential) 
            {
               gr_DMMinPotential               = P[index].Potential;
               Group[gr].index_MinPot = index;
               Group[gr].task_MinPot  = ThisTask;
            }
#endif


      for(j = 0; j < 3; j++)
        {
          xyz[j] = fof_periodic_nearest(P[index].Pos[j] - P[start_index].Pos[j]);
          gr_CM[j] += P[index].Mass * xyz[j];

#ifdef CONSTRUCT_FOF_NGBTREE
          gr_CM_unwrapped[j] += P[index].Mass * P[index].Pos[j];
#endif

          gr_Vel[j] += P[index].Mass * P[index].Vel[j];
#if defined(GFM_BIPOLAR_WINDS) && GFM_BIPOLAR_WINDS != 3
          gr_GravAcc[j] += P[index].Mass * P[index].GravAccel[j];
#endif
        }
    }
#ifdef CREATE_SUBFOFS
  if(All.SubFOF_mode == 0)
   {
#endif
#if defined(EVOLVING_SEEDING_PROBABILITY) || defined(SEED_HALO_ENVIRONMENT_CRITERION2) || defined(UNIFORM_SEEDMASS_DISTRIBUTION)
    for(k = 0; k < len; k++)
      {
        index = FOF_PList[start + k].Pindex;
        if(P[index].Type == PTYPE_GAS)
        {
#ifdef EVOLVING_SEEDING_PROBABILITY
           SphP[index].SecondRandomNumberForSeeding = gr_SecondRandomNumberForSeeding_maxdens;
#endif
#ifdef SEED_HALO_ENVIRONMENT_CRITERION2
           SphP[index].ThirdRandomNumberForSeeding = gr_ThirdRandomNumberForSeeding_maxdens;
#endif
#ifdef UNIFORM_SEEDMASS_DISTRIBUTION
           SphP[index].DrawnSeedBlackHoleMass = gr_DrawnSeedBlackHoleMass_maxdens;
#endif
        }
      }
#endif

#ifdef SEED_BASED_ON_PROBABLISTIC_HALO_PROPERTIES
#ifdef PROBABILISTIC_SEED_MASS_HALO_MASS_RATIO_CRITERION
    for(k = 0; k < len; k++)
      {
       	index = FOF_PList[start + k].Pindex;
        if(P[index].Type == PTYPE_GAS)
        {
           SphP[index].RandomMinHaloMassForSeeding = gr_RandomMinHaloMassForSeeding_maxdens;
        }
      }
#endif
#endif
#ifdef CREATE_SUBFOFS
  }
#endif

  /* put values into group struct */
  Group[gr].Mass = gr_Mass;
#ifdef USE_SFR
  Group[gr].Sfr = gr_Sfr;
#endif

#ifdef EVOLVING_SEEDING_PROBABILITY
  Group[gr].SecondRandomNumberForSeeding_average = gr_SecondRandomNumberForSeeding_average; 
#endif

#ifdef GFM_STELLAR_EVOLUTION
  Group[gr].GasMassMetallicity     = gr_GasMassMetallicity;
  Group[gr].StellarMassMetallicity = gr_StellarMassMetallicity;
  for(k = 0; k < GFM_N_CHEM_ELEMENTS; k++)
    {
      Group[gr].GasMassMetals[k]     = gr_GasMassMetals[k];
      Group[gr].StellarMassMetals[k] = gr_StellarMassMetals[k];
    }
#ifdef GFM_DUST
  Group[gr].GasMassDustMetallicity = gr_GasMassDustMetallicity;
#endif
#endif

#ifdef SEED_HALO_ENVIRONMENT_CRITERION
  Group[gr].DMMinPotential = gr_DMMinPotential;
#endif


#ifdef BLACK_HOLES
  Group[gr].BH_Mass = gr_BH_Mass;
  Group[gr].BH_Mdot = gr_BH_Mdot;
  Group[gr].MaxDens = gr_MaxDens;

#ifdef OUTPUT_LOG_FILES_FOR_SEEDING
  Group[gr].Metallicity_maxdens = gr_Metallicity_maxdens;
  Group[gr].Sfr_maxdens = gr_Sfr_maxdens;
  Group[gr].Mass_maxdens = gr_Mass_maxdens;
  Group[gr].ID_maxdens = gr_ID_maxdens;
#endif

#ifdef SEED_BASED_ON_PROBABLISTIC_HALO_PROPERTIES
#ifdef PROBABILISTIC_SEED_MASS_HALO_MASS_RATIO_CRITERION
  Group[gr].RandomMinHaloMassForSeeding_maxdens = gr_RandomMinHaloMassForSeeding_maxdens;
#endif
#endif

#ifdef EVOLVING_SEEDING_PROBABILITY
  Group[gr].SecondRandomNumberForSeeding_maxdens = gr_SecondRandomNumberForSeeding_maxdens;
#endif

#ifdef SEED_HALO_ENVIRONMENT_CRITERION2
  Group[gr].ThirdRandomNumberForSeeding_maxdens = gr_ThirdRandomNumberForSeeding_maxdens;
#endif

#ifdef UNIFORM_SEEDMASS_DISTRIBUTION
  Group[gr].DrawnSeedBlackHoleMass_maxdens = gr_DrawnSeedBlackHoleMass_maxdens;
#endif

#ifdef CHECK_FOR_ENOUGH_GAS_MASS_IN_DCBH_FORMING_POCKETS
  Group[gr].MaxNeighboringDCBHFormingGasMass = gr_MaxNeighboringDCBHFormingGasMass;
#endif

#ifdef CALCULATE_LYMAN_WERNER_INTENSITY_LOCAL_SOURCES
  Group[gr].LymanWernerIntensityLocalSources_maxdens_type2 =  gr_LymanWernerIntensityLocalSources_maxdens_type2;
  Group[gr].LymanWernerIntensityLocalSources_maxdens_type3 =  gr_LymanWernerIntensityLocalSources_maxdens_type3;
  Group[gr].MaxLymanWernerIntensityInDenseMetalPoorGas = gr_MaxLymanWernerIntensityInDenseMetalPoorGas;
#endif

#ifdef CALCULATE_LYMAN_WERNER_INTENSITY_LOCAL_STARFORMINGGAS
  Group[gr].LymanWernerIntensityLocalStarFormingGas_maxdens_type2 =  gr_LymanWernerIntensityLocalStarFormingGas_maxdens_type2;
  Group[gr].LymanWernerIntensityLocalStarFormingGas_maxdens_type3 =  gr_LymanWernerIntensityLocalStarFormingGas_maxdens_type3;
#endif

#ifdef CALCULATE_LYMAN_WERNER_INTENSITY_ALL_SOURCES
  Group[gr].LymanWernerIntensityAllSources_maxdens_type2   =  gr_LymanWernerIntensityAllSources_maxdens_type2;
  Group[gr].LymanWernerIntensityAllSources_maxdens_type3 =  gr_LymanWernerIntensityAllSources_maxdens_type3;
#endif

#ifdef PREVENT_SEEDING_AROUND_BLACKHOLE_NEIGHBORS2
  Group[gr].BHNeighborExists_maxdens =  gr_BHNeighborExists_maxdens;
#endif

#ifdef PREVENT_SPURIOUS_RESEEDING2
  Group[gr].NeighborOfBlackhole_maxdens =  gr_NeighborOfBlackhole_maxdens;
#endif

#ifdef PREVENT_SPURIOUS_RESEEDING
  Group[gr].TotalGasSeedMass = gr_TotalGasSeedMass;
  Group[gr].SeedMass_maxdens = gr_SeedMass_maxdens;
#endif

  Group[gr].CouldHaveBeenABlackHole_sum = gr_CouldHaveBeenABlackHole_sum;

#ifdef OUTPUT_LOG_FILES_FOR_SEEDING
  Group[gr].StarFormingGasMass = gr_StarFormingGasMass;
  Group[gr].StarFormingGasMassMetallicity = gr_StarFormingGasMassMetallicity;
  Group[gr].StarFormingMetalFreeGasMass = gr_StarFormingMetalFreeGasMass;
#ifdef CALCULATE_LYMAN_WERNER_INTENSITY_LOCAL_STARFORMINGGAS
  Group[gr].StarFormingMetalFreeLymanWernerGasMass = gr_StarFormingMetalFreeLymanWernerGasMass;
  Group[gr].LymanWernerGasMass = gr_LymanWernerGasMass;
#endif
#endif

#ifdef BH_NF_RADIO
  Group[gr].Min_BH_Potential = gr_Min_BH_Potential;
  Group[gr].XrayLum          = gr_XrayLum;
  Group[gr].RadioLum         = gr_RadioLum;
#endif
#endif

#ifdef GFM_WINDS
  Group[gr].WindMass = gr_WindMass;
#endif
  for(k = 0; k < 3; k++)
    {
      Group[gr].CM[k]  = gr_CM[k];
      Group[gr].Vel[k] = gr_Vel[k];
#ifdef CONSTRUCT_FOF_NGBTREE
      Group[gr].CM_unwrapped[k]  = gr_CM_unwrapped[k];
#endif


#if defined(GFM_BIPOLAR_WINDS) && GFM_BIPOLAR_WINDS != 3
      Group[gr].GravAcc[k] = gr_GravAcc[k];
#endif
    }

#if defined(GFM_BIPOLAR_WINDS) && GFM_BIPOLAR_WINDS == 3
  Group[gr].MinPotential = gr_MinPotential;
#endif

  for(k = 0; k < NTYPES; k++)
    Group[gr].MassType[k] = gr_MassType[k];
}

/*---------------------------------*/
#ifdef BH_NF_RADIO
double get_cooling_luminosity(int i)
{
#ifndef SMUGGLE_SFR
  double lum                = 0;
  double eos_dens_threshold = All.PhysDensThresh;

#ifdef MODIFIED_EOS
  eos_dens_threshold *= All.FactorDensThresh;
#endif /* MODIFIED_EOS */

#else
  double lum = 0;
  double eos_dens_threshold = All.NFDensThreshold;
#endif

  double rho = SphP[i].Density * All.cf_a3inv;

  if(rho < eos_dens_threshold)
    {
      double ne = SphP[i].Ne;

#if defined(GFM_AGN_RADIATION) || defined(GFM_UVB_CORRECTIONS)
#ifdef GFM_AGN_RADIATION
      update_radiation_state(rho, SphP[i].MetalsFraction[element_index_Hydrogen], SphP[i].AGNBolIntensity);
#else
      update_radiation_state(rho, SphP[i].MetalsFraction[element_index_Hydrogen], 0);
#endif
#endif

#ifdef GFM_COOLING_METAL
      update_gas_state(rho, SphP[i].MetalsFraction[element_index_Hydrogen], SphP[i].Metallicity);
#endif

#ifdef SMUGGLE_SFR
      double dens = rho * All.UnitDensity_in_cgs * All.HubbleParam * All.HubbleParam; /* convert to physical cgs units */
      double u    = SphP[i].Utherm * All.UnitPressure_in_cgs / All.UnitDensity_in_cgs;
      double temp = convert_u_to_temp(u, dens, &ne);

      ne = SphP[i].Ne;

      /*only hot enough gas contributes to Xray luminosity */
      if(temp < 1.0e4)
        return 0.0;
#endif

      double tcool = GetCoolingTime(SphP[i].Utherm, rho, &ne);

      if(tcool > 0)
        lum = P[i].Mass * SphP[i].Utherm / tcool;
    }

  return lum;
}
#endif

/*! \brief Global exchange of identified groups to their appropriate task.
 *
 *  \return void
 */
void fof_exchange_group_data(void)
{
  struct group_properties *get_Group;
  int i, j, ngrp, recvTask, nimport, start;
  double xyz[3];

  /* sort the groups according to task */
  mysort(Group, NgroupsExt, sizeof(struct group_properties), fof_compare_Group_MinIDTask);

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
              MPI_Sendrecv(&Group[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct group_properties), MPI_BYTE, recvTask,
                           TAG_DENS_A, &get_Group[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(struct group_properties),
                           MPI_BYTE, recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }

  /* sort the groups again according to MinID */
  mysort(Group, NgroupsExt, sizeof(struct group_properties), fof_compare_Group_MinID);
  mysort(get_Group, nimport, sizeof(struct group_properties), fof_compare_Group_MinID);

  /* now add in the partial imported group data to the main ones */
  printf("\n Value of Group[0].MinID %d at task %d",Group[0].MinID,ThisTask);
 
  for(i = 0, start = 0; i < nimport; i++)
    {
      while(Group[start].MinID < get_Group[i].MinID)
        {
          start++;
          if(start >= NgroupsExt)
            terminate("start >= NgroupsExt");
        }

      Group[start].Len += get_Group[i].Len;
      Group[start].Mass += get_Group[i].Mass;

      for(j = 0; j < NTYPES; j++)
        {
          Group[start].LenType[j] += get_Group[i].LenType[j];
          Group[start].MassType[j] += get_Group[i].MassType[j];
        }

#ifdef USE_SFR
      Group[start].Sfr += get_Group[i].Sfr;
#endif
#ifdef BH_NF_RADIO
      Group[start].XrayLum += get_Group[i].XrayLum;
#endif

#ifdef GFM_STELLAR_EVOLUTION
      Group[start].GasMassMetallicity += get_Group[i].GasMassMetallicity;
      Group[start].StellarMassMetallicity += get_Group[i].StellarMassMetallicity;
      for(j = 0; j < GFM_N_CHEM_ELEMENTS; j++)
        {
          Group[start].GasMassMetals[j] += get_Group[i].GasMassMetals[j];
          Group[start].StellarMassMetals[j] += get_Group[i].StellarMassMetals[j];
        }
#ifdef GFM_DUST
      Group[start].GasMassDustMetallicity += get_Group[i].GasMassDustMetallicity;
#endif
#endif
#ifdef BLACK_HOLES
      Group[start].BH_Mdot += get_Group[i].BH_Mdot;
      Group[start].BH_Mass += get_Group[i].BH_Mass;

#ifdef EVOLVING_SEEDING_PROBABILITY
      Group[start].SecondRandomNumberForSeeding_average += get_Group[i].SecondRandomNumberForSeeding_average;
#endif

#ifdef PREVENT_SPURIOUS_RESEEDING
      Group[start].TotalGasSeedMass += get_Group[i].TotalGasSeedMass;
#endif

      Group[start].CouldHaveBeenABlackHole_sum += get_Group[i].CouldHaveBeenABlackHole_sum;

#ifdef OUTPUT_LOG_FILES_FOR_SEEDING
      Group[start].StarFormingGasMass += get_Group[i].StarFormingGasMass;
      Group[start].StarFormingGasMassMetallicity += get_Group[i].StarFormingGasMassMetallicity;
      Group[start].StarFormingMetalFreeGasMass += get_Group[i].StarFormingMetalFreeGasMass;
#ifdef CALCULATE_LYMAN_WERNER_INTENSITY_LOCAL_STARFORMINGGAS
      Group[start].StarFormingMetalFreeLymanWernerGasMass += get_Group[i].StarFormingMetalFreeLymanWernerGasMass;
      Group[start].LymanWernerGasMass += get_Group[i].LymanWernerGasMass;
#endif
#endif

#ifdef CHECK_FOR_ENOUGH_GAS_MASS_IN_DCBH_FORMING_POCKETS
      if(get_Group[i].MaxNeighboringDCBHFormingGasMass > Group[start].MaxNeighboringDCBHFormingGasMass)
         Group[start].MaxNeighboringDCBHFormingGasMass = get_Group[i].MaxNeighboringDCBHFormingGasMass;
#endif
#ifdef CALCULATE_LYMAN_WERNER_INTENSITY_LOCAL_STARFORMINGGAS
      if (get_Group[i].MaxLymanWernerIntensityInDenseMetalPoorGas > Group[start].MaxLymanWernerIntensityInDenseMetalPoorGas)
         Group[start].MaxLymanWernerIntensityInDenseMetalPoorGas = get_Group[i].MaxLymanWernerIntensityInDenseMetalPoorGas;
#endif

#ifdef SEED_HALO_ENVIRONMENT_CRITERION
      if(get_Group[i].DMMinPotential < Group[start].DMMinPotential)
        {
          Group[start].DMMinPotential = get_Group[i].DMMinPotential;
          Group[start].index_MinPot = get_Group[i].index_MinPot;
          Group[start].task_MinPot  = get_Group[i].task_MinPot;
        }

#endif


      if(get_Group[i].MaxDens > Group[start].MaxDens)
        {
          Group[start].MaxDens       = get_Group[i].MaxDens;
          Group[start].index_maxdens = get_Group[i].index_maxdens;
          Group[start].task_maxdens  = get_Group[i].task_maxdens;
#ifdef OUTPUT_LOG_FILES_FOR_SEEDING
          Group[start].Metallicity_maxdens = get_Group[i].Metallicity_maxdens;
          Group[start].Sfr_maxdens = get_Group[i].Sfr_maxdens;
          Group[start].Mass_maxdens = get_Group[i].Mass_maxdens;
          Group[start].ID_maxdens = get_Group[i].ID_maxdens;
#endif

#ifdef SEED_BASED_ON_PROBABLISTIC_HALO_PROPERTIES
#ifdef PROBABILISTIC_SEED_MASS_HALO_MASS_RATIO_CRITERION
         Group[start].RandomMinHaloMassForSeeding_maxdens = get_Group[i].RandomMinHaloMassForSeeding_maxdens;
#endif
#endif

#ifdef EVOLVING_SEEDING_PROBABILITY
         Group[start].SecondRandomNumberForSeeding_maxdens = get_Group[i].SecondRandomNumberForSeeding_maxdens;
#endif

#ifdef SEED_HALO_ENVIRONMENT_CRITERION2
         Group[start].ThirdRandomNumberForSeeding_maxdens = get_Group[i].ThirdRandomNumberForSeeding_maxdens;
#endif

#ifdef UNIFORM_SEEDMASS_DISTRIBUTION
         Group[start].DrawnSeedBlackHoleMass_maxdens = get_Group[i].DrawnSeedBlackHoleMass_maxdens;
#endif

#ifdef CALCULATE_LYMAN_WERNER_INTENSITY_LOCAL_SOURCES
          Group[start].LymanWernerIntensityLocalSources_maxdens_type2 =  get_Group[i].LymanWernerIntensityLocalSources_maxdens_type2;
          Group[start].LymanWernerIntensityLocalSources_maxdens_type3 =  get_Group[i].LymanWernerIntensityLocalSources_maxdens_type3;
#endif

#ifdef CALCULATE_LYMAN_WERNER_INTENSITY_LOCAL_STARFORMINGGAS
          Group[start].LymanWernerIntensityLocalStarFormingGas_maxdens_type2 =  get_Group[i].LymanWernerIntensityLocalStarFormingGas_maxdens_type2;
          Group[start].LymanWernerIntensityLocalStarFormingGas_maxdens_type3 =  get_Group[i].LymanWernerIntensityLocalStarFormingGas_maxdens_type3;
#endif

#ifdef CALCULATE_LYMAN_WERNER_INTENSITY_ALL_SOURCES
          Group[start].LymanWernerIntensityAllSources_maxdens_type2   =  get_Group[i].LymanWernerIntensityAllSources_maxdens_type2;
          Group[start].LymanWernerIntensityAllSources_maxdens_type3 =  get_Group[i].LymanWernerIntensityAllSources_maxdens_type3;
#endif

#ifdef PREVENT_SPURIOUS_RESEEDING
          Group[start].SeedMass_maxdens = get_Group[i].SeedMass_maxdens;
#endif

#ifdef PREVENT_SPURIOUS_RESEEDING2
          Group[start].NeighborOfBlackhole_maxdens = get_Group[i].NeighborOfBlackhole_maxdens;
#endif

#ifdef PREVENT_SEEDING_AROUND_BLACKHOLE_NEIGHBORS2
          Group[start].BHNeighborExists_maxdens = get_Group[i].BHNeighborExists_maxdens;
#endif
        }

#if defined(BH_NF_RADIO)
      if(get_Group[i].Min_BH_Potential < Group[start].Min_BH_Potential)
        {
          Group[start].Min_BH_Potential    = get_Group[i].Min_BH_Potential;
          Group[start].ID_Min_BH_Potential = get_Group[i].ID_Min_BH_Potential;
        }
#endif

#endif

#ifdef GFM_WINDS
      Group[start].WindMass += get_Group[i].WindMass;
#endif

#if defined(CALCULATE_SPIN_STARFORMINGGAS) || (defined(GFM_BIPOLAR_WINDS) && (GFM_BIPOLAR_WINDS == 3))
      if(get_Group[i].MinPotential < Group[start].MinPotential)
        {
          Group[start].MinPotential        = get_Group[i].MinPotential;
          Group[start].Pos_MinPotential[0] = get_Group[i].Pos_MinPotential[0];
          Group[start].Pos_MinPotential[1] = get_Group[i].Pos_MinPotential[1];
          Group[start].Pos_MinPotential[2] = get_Group[i].Pos_MinPotential[2];
        }
#endif

      for(j = 0; j < 3; j++)
        {
          xyz[j] = fof_periodic_nearest(get_Group[i].CM[j] / get_Group[i].Mass + get_Group[i].FirstPos[j] - Group[start].FirstPos[j]);
          Group[start].CM[j] += get_Group[i].Mass * xyz[j];
          Group[start].Vel[j] += get_Group[i].Vel[j];
#ifdef CONSTRUCT_FOF_NGBTREE
          Group[start].CM_unwrapped[j] += get_Group[i].CM_unwrapped[j];
#endif
#if defined(GFM_BIPOLAR_WINDS) && GFM_BIPOLAR_WINDS != 3
          Group[start].GravAcc[j] += get_Group[i].GravAcc[j];
#endif
        }
    }

  myfree(get_Group);
}

/*! \brief Finalizes group property calculation.
 *
 *  Called after a loop over all particles of a group is already completed.
 *
 *  \return void
 */
void fof_finish_group_properties(void)
{
  for(int i = 0; i < NgroupsExt; i++)
    {
#ifdef BLACK_HOLES
#ifdef EVOLVING_SEEDING_PROBABILITY
      Group[i].SecondRandomNumberForSeeding_average /= Group[i].MassType[0];
#endif   
#endif
      if(Group[i].MinIDTask == (MyIDType)ThisTask)
        {
          for(int j = 0; j < 3; j++)
            {
              Group[i].Vel[j] /= Group[i].Mass;
#if defined(GFM_BIPOLAR_WINDS) && GFM_BIPOLAR_WINDS != 3
              Group[i].GravAcc[j] /= Group[i].Mass;
#endif
              Group[i].CM[j] = fof_periodic_wrap(Group[i].CM[j] / Group[i].Mass + Group[i].FirstPos[j]);
#ifdef CONSTRUCT_FOF_NGBTREE
              Group[i].CM_unwrapped[j] /= Group[i].Mass;
#endif
            }

#ifdef BH_NF_RADIO
          if(Group[i].LenType[PTYPE_BNDRY] > 0) /* only if a BH is actually present */
            {
              double vvir       = pow(10.0 * All.G * All.cf_H * Group[i].Mass, 1.0 / 3);
              Group[i].RadioLum = (All.Hubble / All.cf_H) * blackhole_get_radio_efficiency(vvir) * Group[i].XrayLum;
            }
          else
            Group[i].RadioLum = 0;
#endif
        }
    }

  /* eliminate the non-local groups */
  int ngr = NgroupsExt;
  for(int i = 0; i < ngr; i++)
    {
      if(Group[i].MinIDTask != (MyIDType)ThisTask)
        {
          Group[i] = Group[ngr - 1];
          i--;
          ngr--;
        }
    }

  if(ngr != Ngroups)
    terminate("ngr != Ngroups");

  mysort(Group, Ngroups, sizeof(struct group_properties), fof_compare_Group_MinID);
}

/*! \brief Do periodic wrap for distances.
 *
 *         Note that this only works for a cubic box.
 *
 *  \param[in] x Coordinate.
 *
 *  \return distance within [-0.5*BoxSize, 0.5*BoxSize].
 */
double fof_periodic_nearest(const double x) { return coord_nearest(All.BoxSize, 0.5 * All.BoxSize, x); }

/*! \brief Do periodic wrap for coordinates.
 *
 *         Note that this only works for a cubic box.
 *
 *  \param[in] x Coordinate.
 *
 *  \return coordinate within [0, BoxSize).
 */
double fof_periodic_wrap(double x)
{
  while(x >= All.BoxSize)
    x -= All.BoxSize;
  while(x < 0)
    x += All.BoxSize;
  return x;
}
#ifdef PROBABILISTIC_SEEDING
float get_random_fraction(int number_of_decimal_places)
{
    int upper = pow(10 , number_of_decimal_places), lower = 0, num;
    float random_fraction;
    num =  (rand() % (upper - lower + 1)) + lower;
    random_fraction = (float)num / upper;
    return random_fraction;
}
#endif	
#endif /* of FOF */

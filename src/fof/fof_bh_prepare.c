/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/fof/fof_bh.c
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
#include <stdbool.h>

#include "../allvars.h"
#include "../domain.h"
#include "../proto.h"
#include "../subfind/subfind.h"
#include "fof.h"

#ifdef FOF
#ifdef BLACK_HOLES

void fof_prepare_to_seed_black_holes(void)
{
  int idx, i, j, n, ntot;
  int nexport, nimport, recvTask, level;
  int *import_indices, *export_indices;
  double csnd = 0, rad = 0, dt_courant = 0;
  char msg[200];
  
#ifdef BH_BASED_CGM_ZOOM
  int local_bh_on_task = -1, new_fof_bh_flag = 0;
#endif

#ifdef SEED_HALO_ENVIRONMENT_CRITERION2
  MyFloat EnvironmentBasedSeedProbability;  
#endif

  MyFloat SeedBlackHoleMass;


#ifdef SEED_BASED_ON_PROBABLISTIC_HALO_PROPERTIES
#ifdef PROBABILISTIC_SEED_MASS_HALO_MASS_RATIO_CRITERION
    int SeedFromLowerEndLocal,SeedFromLowerEndLocalTemp; 
    int StartFromLeftLocal,StartFromLeftLocalTemp;
    MyFloat MinHaloMassForSeedingLowerEndLocal,MinHaloMassForSeedingLowerEndLocalTemp;
    int SeedTheBH;
    MPI_Status status;
    MyFloat AppliedLogMinHaloMassSeedingThresholdAverage, AppliedLogMinHaloMassSeedingThresholdStdev, CurrentRedshift;
    CurrentRedshift = 1./All.Time - 1.;
    AppliedLogMinHaloMassSeedingThresholdAverage = All.LogMinHaloMassSeedingThresholdAverage;
    AppliedLogMinHaloMassSeedingThresholdStdev = All.LogMinHaloMassSeedingThresholdStdev;
#ifdef EVOLVING_SEEDHALOMASS_DISTRIBUTION
    MyFloat slope = 0.06; 
    AppliedLogMinHaloMassSeedingThresholdAverage += fmax(slope * (CurrentRedshift - 12.) , 0.); 
#endif 
#ifdef EVOLVING_SEEDHALOMASS_DISTRIBUTION_METALENRICHMENT
    AppliedLogMinHaloMassSeedingThresholdAverage = (AppliedLogMinHaloMassSeedingThresholdAverage + 10.) * fmax((CurrentRedshift - All.RedshiftOnsetMetalEnrichmentImpact) * All.SlopeOnsetMetalEnrichmentImpact + 1. , 1.)  - 10.;
#endif

#ifdef EVOLVING_SEEDHALOMASS_DISTRIBUTION_QUADRATIC_MODEL
    AppliedLogMinHaloMassSeedingThresholdAverage = All.Coffmean_2 * CurrentRedshift * CurrentRedshift + All.Coffmean_1 * CurrentRedshift + All.Coffmean_0 - 10.;
#endif

#ifdef EVOLVING_SEEDHALOMASS_DISTRIBUTION_DOUBLE_POWERLAW_MODEL
    AppliedLogMinHaloMassSeedingThresholdAverage = fmax(All.Slope_lowz * CurrentRedshift + All.Intercept_lowz , All.Slope_highz * CurrentRedshift + All.Intercept_highz) - 10.;
#endif

#ifdef EVOLVING_SEEDHALOMASS_DISTRIBUTION_DOUBLE_POWERLAW_MODEL2
    AppliedLogMinHaloMassSeedingThresholdAverage = fmax(All.Slope_lowz * (CurrentRedshift - All.TransitionRedshift) , All.Slope_highz * (CurrentRedshift - All.TransitionRedshift)) + All.LogMinHaloMassSeedingThresholdAverageTransition - 10.;
#endif

    MyFloat SelectedRandomNumberForSeeding;
    double SenderTime;
    MPI_Barrier(MPI_COMM_WORLD);
    if(ThisTask == 0)
    {
        SenderTime = 0;
        SeedFromLowerEndLocal = All.SeedFromLowerEnd;
        StartFromLeftLocal = All.StartFromLeft;
        MinHaloMassForSeedingLowerEndLocal = All.MinHaloMassForSeedingLowerEnd;
        SeedFromLowerEndLocalTemp = SeedFromLowerEndLocal;
        StartFromLeftLocalTemp = StartFromLeftLocal;
        MinHaloMassForSeedingLowerEndLocalTemp = MinHaloMassForSeedingLowerEndLocal;
    }
    else
    {
        //printf("Waiting to recv at task %d", ThisTask);
        MPI_Recv(&SenderTime, 1, MPI_DOUBLE, ThisTask-1, TAG_BHSEED4, MPI_COMM_WORLD, &status);
        MPI_Recv(&SeedFromLowerEndLocal, 1, MPI_INT, ThisTask-1, TAG_BHSEED1, MPI_COMM_WORLD, &status);
        MPI_Recv(&StartFromLeftLocal, 1, MPI_INT, ThisTask-1, TAG_BHSEED2, MPI_COMM_WORLD, &status);
        MPI_Recv(&MinHaloMassForSeedingLowerEndLocal, 1, MPI_DOUBLE, ThisTask-1, TAG_BHSEED3, MPI_COMM_WORLD, &status);
        SeedFromLowerEndLocalTemp = SeedFromLowerEndLocal;
        StartFromLeftLocalTemp = StartFromLeftLocal;
        MinHaloMassForSeedingLowerEndLocalTemp = MinHaloMassForSeedingLowerEndLocal;
//        mpi_printf_task(ThisTask,"Recieved %.8f at task %d at time %.7f", MinHaloMassForSeedingLowerEndLocal, ThisTask, All.Time);
    }

//    mpi_printf_task(ThisTask,"\nVery Ready to make black holes %d at time %.7f, SenderTime %.7f, !!\n",ThisTask,All.Time, SenderTime);

#endif
#endif

  for(n = 0; n < NTask; n++)
    Send_count[n] = 0;

  for(i = 0; i < Ngroups; i++)
    {      
#ifdef UNIFORM_SEEDMASS_DISTRIBUTION
     SeedBlackHoleMass = Group[i].DrawnSeedBlackHoleMass_maxdens;
#else
     SeedBlackHoleMass = All.SeedBlackHoleMass;
#endif

#ifdef SEED_HALO_ENVIRONMENT_CRITERION2
     if((Group2[i].NumberOfMajorNeighbors == 0) && (All.TotNumBHs > All.MinNumberOfBlackHolesForEnvironmentBasedSeeding))
         EnvironmentBasedSeedProbability = All.SeedProbabilityNoNeighbors;
     else if((Group2[i].NumberOfMajorNeighbors == 1) && (All.TotNumBHs > All.MinNumberOfBlackHolesForEnvironmentBasedSeeding))
         EnvironmentBasedSeedProbability = All.SeedProbabilityTwentyNeighbors;
     else
         EnvironmentBasedSeedProbability = 1;
#endif

#if defined(CORRECT_FOR_HALO_MASS_BIAS_IN_ENVIRONMENT_BASED_SEEDING) && defined(SEED_HALO_ENVIRONMENT_CRITERION2)
     if(EnvironmentBasedSeedProbability < 0.999)
       	 EnvironmentBasedSeedProbability *= All.SlopeForEnvironmentBasedSeedProbabilityvsbFOFmass * (log10(Group[i].Mass * 1e10) - (10 + AppliedLogMinHaloMassSeedingThresholdAverage)) + 1;
#endif

#ifdef PREVENT_SPURIOUS_RESEEDING
    if(Group[i].TotalGasSeedMass < All.MaxPaintedGasFractionForReseeding * SeedBlackHoleMass)
    if(Group[i].SeedMass_maxdens == 0)
     {
#endif
#ifdef PREVENT_SPURIOUS_RESEEDING2
    if(Group[i].NeighborOfBlackhole_maxdens == 1)
     {
#endif

#ifdef SEED_BASED_ON_PROBABLISTIC_HALO_PROPERTIES
#ifdef PROBABILISTIC_SEED_MASS_HALO_MASS_RATIO_CRITERION
      SeedTheBH = 0;
      if((Group[i].LenType[5] == 0) && (Group[i].CouldHaveBeenABlackHole_sum == 0))
       {
        if ((StartFromLeftLocal == 0))
           SelectedRandomNumberForSeeding = Group[i].RandomMinHaloMassForSeeding_maxdens;
        else if ((StartFromLeftLocal == 1))
           SelectedRandomNumberForSeeding = -Group[i].RandomMinHaloMassForSeeding_maxdens;
        if ((SeedFromLowerEndLocal == 0) & (Group[i].LenType[1] * All.massDMpart >= (All.Omega0 - All.OmegaBaryon) / All.Omega0 * pow(10,SelectedRandomNumberForSeeding * AppliedLogMinHaloMassSeedingThresholdStdev + AppliedLogMinHaloMassSeedingThresholdAverage + All.LogMinHaloMassSeedingThresholdOffset)))
            {
              SeedFromLowerEndLocal = 1;
              if (StartFromLeftLocal == 0)
                  MinHaloMassForSeedingLowerEndLocal = -SelectedRandomNumberForSeeding;
       	      if (StartFromLeftLocal == 1)
                  MinHaloMassForSeedingLowerEndLocal = fmin(-SelectedRandomNumberForSeeding,1);
              SeedTheBH =	1;
       	    }
         else if ((SeedFromLowerEndLocal == 1) & (Group[i].LenType[1] * All.massDMpart >= (All.Omega0 - All.OmegaBaryon) / All.Omega0 * pow(10,MinHaloMassForSeedingLowerEndLocal * AppliedLogMinHaloMassSeedingThresholdStdev + AppliedLogMinHaloMassSeedingThresholdAverage - All.LogMinHaloMassSeedingThresholdOffset)))
       	    {
              if (StartFromLeftLocal == 1) StartFromLeftLocal = 0;
              else if (StartFromLeftLocal == 0) StartFromLeftLocal = 1;
              SeedFromLowerEndLocal = 0;
       	      SeedTheBH =	1;
       	    }
       }
#endif
#endif

#ifdef PREVENT_SPURIOUS_RESEEDING 
     }
#endif
#ifdef PREVENT_SPURIOUS_RESEEDING2
    }    
#endif

#ifdef BH_BASED_CGM_ZOOM
      if((new_fof_bh_flag == 0) & (All.TotNumBHs == 0))
#endif
#ifdef GAS_BASED_SEED_MODEL
#ifdef ONLY_SEED_IN_VALID_FOFS
        if(Group[i].LenType[1] >= FOF_GROUP_MIN_LEN) 
#endif
#ifdef SEED_STARFORMINGGASMASS_CRITERION
        if(Group[i].StarFormingGasMass > SeedBlackHoleMass * All.MinStarFormingGasParamForNewSeed)
#endif
#ifdef SEED_STARFORMINGMETALFREEGASMASS_CRITERION
        if(Group[i].StarFormingMetalFreeGasMass > SeedBlackHoleMass * All.MinStarFormingMetalFreeGasParamForNewSeed)
#endif
#ifdef SEED_STARFORMINGMETALFREELYMANWERNERGASMASS_CRITERION
        if(Group[i].StarFormingMetalFreeLymanWernerGasMass > SeedBlackHoleMass * All.MinStarFormingMetalFreeLymanWernerGasParamForNewSeed)
#ifdef CHECK_FOR_ENOUGH_GAS_MASS_IN_DCBH_FORMING_POCKETS
        if(Group[i].MaxNeighboringDCBHFormingGasMass > SeedBlackHoleMass * All.MinStarFormingMetalFreeLymanWernerGasParamForNewSeed)
#endif
#endif
#ifdef SEED_LYMANWERNERGASMASS_CRITERION
        if(Group[i].LymanWernerGasMass > SeedBlackHoleMass * All.MinLymanWernerGasParamForNewSeed)
#endif
#ifdef SEED_GASSPIN_CRITERION
        if(Group[i].DensGasDimensionlessSpin < Group[i].DensGasDimensionlessSpin_Max)
        if(Group[i].Mass > SeedBlackHoleMass * 20./sqrt(1-Group[i].DensGasDimensionlessSpin / Group[i].DensGasDimensionlessSpin_Max))         
#endif
#ifdef SEED_STARFORMINGGASMETALLICITY_CRITERION
        if(Group[i].StarFormingGasMassMetallicity < All.MaxStarFormingGasMetallicityForNewSeed * Group[i].StarFormingGasMass * GFM_SOLAR_METALLICITY)
#endif
#ifdef SEED_GAS_DENSITY_CRITERION
#if(SEED_GAS_DENSITY_CRITERION == 0)
        if ((Group[i].MaxDens * All.cf_a3inv >= All.PhysDensThresh && Group[i].MaxDens >= All.OverDensThresh) && (Group[i].Mass_maxdens > 0.0))            
#else
        if(Group[i].Sfr_maxdens > 0)
#endif
#endif
#ifdef SEED_GAS_METALLICITY_CRITERION
        if (Group[i].Metallicity_maxdens  >= (All.MinMetallicityForNewSeed*GFM_SOLAR_METALLICITY) && Group[i].Metallicity_maxdens  <= (All.MaxMetallicityForNewSeed*GFM_SOLAR_METALLICITY))
#endif
#ifdef SEED_MASS_HALO_MASS_RATIO_CRITERION
#if(SEED_MASS_HALO_MASS_RATIO_CRITERION == 0)
        if(Group[i].LenType[1] * All.massDMpart >= (All.Omega0 - All.OmegaBaryon) / All.Omega0 * All.MinSeedMassHaloMassRatioForNewSeed * SeedBlackHoleMass)
#else
        if(Group[i].LenType[1] * All.massDMpart >= (All.Omega0 - All.OmegaBaryon) / All.Omega0 * All.MinSeedMassHaloMassRatioForNewSeed * SeedBlackHoleMass * (1 + All.MinSeedMassHaloMassLymanWernerParam * (Group[i].LymanWernerIntensityLocalStarFormingGas_maxdens_type2 + Group[i].LymanWernerIntensityLocalStarFormingGas_maxdens_type3)))
#endif
#endif
#ifdef NO_SEEDING_IN_LOW_RESOLUTION_CONTAMINATED_HALOS
       if(Group[i].MassType[2] / (Group[i].MassType[1] + Group[i].MassType[2]) < All.MaximumLowResolutionContaminationForSeeding)
#endif
#ifdef SEED_LYMAN_WERNER_INTENSITY_CRITERION
        if(Group[i].LymanWernerIntensityLocalSources_maxdens_type2 + Group[i].LymanWernerIntensityLocalSources_maxdens_type3 > All.MinLymanWernerFluxForNewSeed)
#endif
#ifdef SEED_HALO_ENVIRONMENT_CRITERION1
        if(Group2[i].NumberOfMajorNeighbors >= All.MinNumberOfNeighborsForSeeding)
#endif
#ifdef SEED_HALO_ENVIRONMENT_CRITERION2
        if(Group[i].ThirdRandomNumberForSeeding_maxdens < EnvironmentBasedSeedProbability)
#endif
#ifdef SEED_BASED_ON_PROBABLISTIC_HALO_PROPERTIES
#ifdef PROBABILISTIC_SEED_MASS_HALO_MASS_RATIO_CRITERION
        if(SeedTheBH == 1)                                                                                                                            
#endif
#endif
#else 
        if(Group[i].LenType[1] * All.massDMpart >= (All.Omega0 - All.OmegaBaryon) / All.Omega0 * All.MinFoFMassForNewSeed)
#endif
#ifdef PROBABILISTIC_SEEDING
        if(Group[i].RandomFractionForSeed < All.ConstantSeedProbability) 
#endif
#ifdef PREVENT_SPURIOUS_RESEEDING
        if(Group[i].TotalGasSeedMass < All.MaxPaintedGasFractionForReseeding * SeedBlackHoleMass)
        if(Group[i].SeedMass_maxdens == 0)
#endif
#ifdef PREVENT_SEEDING_AROUND_BLACKHOLE_NEIGHBORS2
        if(Group[i].BHNeighborExists_maxdens == 1) 
#endif
#ifdef PREVENT_SPURIOUS_RESEEDING2
        if(Group[i].NeighborOfBlackhole_maxdens == 1)
#endif
            if((Group[i].LenType[5] == 0) && (Group[i].CouldHaveBeenABlackHole_sum == 0))
	      {
		if(Group[i].index_maxdens >= 0)
#ifdef GAS_BASED_SEED_MODEL
                  if(Group[i].Mass_maxdens != 0 || Group[i].ID_maxdens != 0)
#else
                  if(P[Group[i].index_maxdens].Mass != 0 || P[Group[i].index_maxdens].ID != 0)
#endif


                    Group[i].AllSeedingCriteriaSatisfied = 1;
		    Send_count[Group[i].task_maxdens]++;
#ifdef BH_BASED_CGM_ZOOM
		new_fof_bh_flag = 1;
#endif
	      }
    }


#ifdef SEED_BASED_ON_PROBABLISTIC_HALO_PROPERTIES
#ifdef PROBABILISTIC_SEED_MASS_HALO_MASS_RATIO_CRITERION
   if(ThisTask == NTask-1)
   {
        All.SeedFromLowerEnd  = SeedFromLowerEndLocal;
        All.StartFromLeft = StartFromLeftLocal;
        All.MinHaloMassForSeedingLowerEnd = MinHaloMassForSeedingLowerEndLocal;
        MPI_Send(&All.SeedFromLowerEnd, 1, MPI_INT, 0, TAG_BHSEED5, MPI_COMM_WORLD);
        MPI_Send(&All.StartFromLeft, 1, MPI_INT, 0, TAG_BHSEED6, MPI_COMM_WORLD);
        MPI_Send(&All.MinHaloMassForSeedingLowerEnd, 1, MPI_DOUBLE, 0, TAG_BHSEED7, MPI_COMM_WORLD);
   }
   else
   {
        MPI_Send(&All.Time, 1, MPI_DOUBLE, ThisTask + 1, TAG_BHSEED4, MPI_COMM_WORLD);
        MPI_Send(&SeedFromLowerEndLocal, 1, MPI_INT, ThisTask + 1, TAG_BHSEED1, MPI_COMM_WORLD);
        MPI_Send(&StartFromLeftLocal, 1, MPI_INT, ThisTask + 1, TAG_BHSEED2, MPI_COMM_WORLD);
        MPI_Send(&MinHaloMassForSeedingLowerEndLocal, 1, MPI_DOUBLE, ThisTask + 1, TAG_BHSEED3, MPI_COMM_WORLD);
//        mpi_printf_task(ThisTask,"Sent %.8f by task %d at time %.7f", MinHaloMassForSeedingLowerEndLocal, ThisTask, All.Time);
   }

   if(ThisTask == 0)
   {
        MPI_Recv(&All.SeedFromLowerEnd, 1, MPI_INT, NTask-1, TAG_BHSEED5, MPI_COMM_WORLD, &status);
        MPI_Recv(&All.StartFromLeft, 1, MPI_INT, NTask-1, TAG_BHSEED6, MPI_COMM_WORLD, &status);
        MPI_Recv(&All.MinHaloMassForSeedingLowerEnd, 1, MPI_DOUBLE, NTask-1, TAG_BHSEED7, MPI_COMM_WORLD, &status);
   }

#endif
#endif

 MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);
 for(j = 0, nimport = nexport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nexport += Send_count[j];
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  import_indices = (int *)mymalloc("import_indices", nimport * sizeof(int));
  export_indices = (int *)mymalloc("export_indices", nexport * sizeof(int));

  for(n = 0; n < NTask; n++)
    Send_count[n] = 0;

#ifdef BH_BASED_CGM_ZOOM
  new_fof_bh_flag = 0;
#endif


  for(i = 0; i < Ngroups; i++)
    {
         if (Group[i].AllSeedingCriteriaSatisfied == 1)
           {
               export_indices[Send_offset[Group[i].task_maxdens] + Send_count[Group[i].task_maxdens]] = Group[i].index_maxdens;
               Send_count[Group[i].task_maxdens]++;
           }
    }
    

  memcpy(&import_indices[Recv_offset[ThisTask]], &export_indices[Send_offset[ThisTask]], Send_count[ThisTask] * sizeof(int));


  for(level = 1; level < (1 << PTask); level++)
    {
      recvTask = ThisTask ^ level;

      if(recvTask < NTask)
        MPI_Sendrecv(&export_indices[Send_offset[recvTask]], Send_count[recvTask] * sizeof(int), MPI_BYTE, recvTask, TAG_FOF_E,
                     &import_indices[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(int), MPI_BYTE, recvTask, TAG_FOF_E,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

  MPI_Allreduce(&nimport, &ntot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  mpi_printf("BLACK_HOLES: Making %d new black hole particles\n", ntot);

  All.TotNumBHs += ntot;

  for(n = 0; n < nimport; n++)
    {
     if(P[import_indices[n]].Type != 0)
        {
          sprintf(msg, "Particle n = %d import_indices[n] = %d is not a gas particle\n", n, import_indices[n]);
          terminate(msg);
        }
      SphP[import_indices[n]].CouldHaveBeenABlackHole = 1;
    }

  myfree(export_indices);
  myfree(import_indices);

}
#endif /* BLACK_HOLES */
#endif

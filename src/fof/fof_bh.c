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

void fof_make_black_holes(void)
{
  int idx, i, j, n, ntot;
  int nexport, nimport, recvTask, level;
  int *import_indices, *export_indices;
#ifdef OUTPUT_LOG_FILES_FOR_SEEDING
  MyFloat *import_parentFOFDMMass,*export_parentFOFDMMass;
  int *import_parentFOFTask,*export_parentFOFTask;
#endif
  double csnd = 0, rad = 0, dt_courant = 0;
  char msg[200];

#ifdef EVOLVING_SEEDING_PROBABILITY
  MyFloat Current_Redshift = 1/All.Time - 1;
  All.RedshiftDependentSeedingProbability = pow(10,Current_Redshift * All.SlopeLinearModelForSeedingProbability + All.InterceptLinearModelForSeedingProbability);
#endif
  
#ifdef BH_BASED_CGM_ZOOM
  int local_bh_on_task = -1, new_fof_bh_flag = 0;
#endif

  for(n = 0; n < NTask; n++)
    Send_count[n] = 0;

  for(i = 0; i < Ngroups; i++)
    {
#ifdef BH_BASED_CGM_ZOOM
      if((new_fof_bh_flag == 0) & (All.TotNumBHs == 0))
#endif      
         if(Group[i].AllSeedingCriteriaSatisfied == 1)
#ifdef EVOLVING_SEEDING_PROBABILITY
         if (Group[i].SecondRandomNumberForSeeding_maxdens < All.RedshiftDependentSeedingProbability)
#endif
	      {
		    Send_count[Group[i].task_maxdens]++;
#ifdef BH_BASED_CGM_ZOOM
		    new_fof_bh_flag = 1;
#endif
	      }
    }

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

#ifdef OUTPUT_LOG_FILES_FOR_SEEDING
  import_parentFOFDMMass = (MyFloat *)mymalloc("import_parentFOFDMMass", nimport * sizeof(MyFloat));
  export_parentFOFDMMass = (MyFloat *)mymalloc("export_parentFOFDMMass", nexport * sizeof(MyFloat));
  import_parentFOFTask = (int *)mymalloc("import_parentFOFTask", nimport * sizeof(int));
  export_parentFOFTask = (int *)mymalloc("export_parentFOFTask", nexport * sizeof(int));
#endif

  for(n = 0; n < NTask; n++)
    Send_count[n] = 0;

#ifdef BH_BASED_CGM_ZOOM
  new_fof_bh_flag = 0;
#endif
 
  int extra_output1 = 0;
  double extra_output2 = 0;
  
  for(i = 0; i < Ngroups; i++)
    {
#ifdef SEED_HALO_ENVIRONMENT_CRITERION
     extra_output1 = Group2[i].NumberOfMajorNeighbors;
#endif

#ifdef UNIFORM_SEEDMASS_DISTRIBUTION
     extra_output2 = Group[i].DrawnSeedBlackHoleMass_maxdens;
#endif

#ifdef BH_BASED_CGM_ZOOM
      if((new_fof_bh_flag == 0) & (All.TotNumBHs == 0))
#endif
         if(Group[i].AllSeedingCriteriaSatisfied == 1)
#ifdef EVOLVING_SEEDING_PROBABILITY
         if (Group[i].SecondRandomNumberForSeeding_maxdens < All.RedshiftDependentSeedingProbability)
#endif
            {                  
#ifdef OUTPUT_LOG_FILES_FOR_SEEDING
                     fprintf(FdBlackHolesSeeding2,
"%d %g %llu %g %g %g %g %g %g %g %g %d %g\n"
, ThisTask, All.Time, (long long)Group[i].ID_maxdens,Group[i].MassType[4],
Group[i].StellarMassMetallicity/Group[i].MassType[4], Group[i].Sfr, Group[i].MassType[1], Group[i].GasMassMetallicity/Group[i].MassType[0], Group[i].StarFormingGasMass, 
Group[i].StarFormingGasMassMetallicity/Group[i].StarFormingGasMass,Group[i].StarFormingMetalFreeGasMass,extra_output1,extra_output2);
                     myflush(FdBlackHolesSeeding2);
#endif
                     export_indices[Send_offset[Group[i].task_maxdens] + Send_count[Group[i].task_maxdens]] = Group[i].index_maxdens;
#ifdef OUTPUT_LOG_FILES_FOR_SEEDING
                     export_parentFOFDMMass[Send_offset[Group[i].task_maxdens] + Send_count[Group[i].task_maxdens]] = Group[i].MassType[1];
                     export_parentFOFTask[Send_offset[Group[i].task_maxdens] + Send_count[Group[i].task_maxdens]] = ThisTask;
#endif
                     Send_count[Group[i].task_maxdens]++;
                  
#ifdef BH_BASED_CGM_ZOOM
              new_fof_bh_flag = 1;
#endif
            }
    }

  memcpy(&import_indices[Recv_offset[ThisTask]], &export_indices[Send_offset[ThisTask]], Send_count[ThisTask] * sizeof(int));

#ifdef OUTPUT_LOG_FILES_FOR_SEEDING
  memcpy(&import_parentFOFDMMass[Recv_offset[ThisTask]], &export_parentFOFDMMass[Send_offset[ThisTask]], Send_count[ThisTask] * sizeof(MyFloat));
  memcpy(&import_parentFOFTask[Recv_offset[ThisTask]], &export_parentFOFTask[Send_offset[ThisTask]], Send_count[ThisTask] * sizeof(int));
#endif

  for(level = 1; level < (1 << PTask); level++)
    {
      recvTask = ThisTask ^ level;

      if(recvTask < NTask)
        MPI_Sendrecv(&export_indices[Send_offset[recvTask]], Send_count[recvTask] * sizeof(int), MPI_BYTE, recvTask, TAG_FOF_E,
                     &import_indices[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(int), MPI_BYTE, recvTask, TAG_FOF_E,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#ifdef OUTPUT_LOG_FILES_FOR_SEEDING
        MPI_Sendrecv(&export_parentFOFDMMass[Send_offset[recvTask]], Send_count[recvTask] * sizeof(MyFloat), MPI_BYTE, recvTask, TAG_FOF_DMMASS,
                     &import_parentFOFDMMass[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(MyFloat), MPI_BYTE, recvTask, TAG_FOF_DMMASS,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Sendrecv(&export_parentFOFTask[Send_offset[recvTask]], Send_count[recvTask] * sizeof(int), MPI_BYTE, recvTask, TAG_FOF_TASK, 
                     &import_parentFOFTask[Recv_offset[recvTask]], Recv_count[recvTask] * sizeof(int), MPI_BYTE, recvTask, TAG_FOF_TASK,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#endif
    }

  MPI_Allreduce(&nimport, &ntot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  mpi_printf("BLACK_HOLES: Making %d new black hole particles\n", ntot);

  All.TotNumBHs += ntot;

  for(n = 0; n < nimport; n++)
    {
#ifdef OUTPUT_LOG_FILES_FOR_SEEDING
      fprintf(FdBlackHolesSeeding, "%d %g %llu %g %g %g %g %d %d %g %g\n", ThisTask, All.Time, (long long)P[import_indices[n]].ID, SphP[import_indices[n]].Density, SphP[import_indices[n]].Metallicity, SphP[import_indices[n]].Sfr, import_parentFOFDMMass[n], import_indices[n], import_parentFOFTask[n], SphP[import_indices[n]].Hsml, get_default_softening_of_particletype(5));
      myflush(FdBlackHolesSeeding);
#endif

      if(P[import_indices[n]].Type != 0)
        {
          sprintf(msg, "Particle n = %d import_indices[n] = %d is not a gas particle\n", n, import_indices[n]);
          terminate(msg);
        }

      P[import_indices[n]].Type          = 5; /* make it a black hole particle */
      P[import_indices[n]].SofteningType = All.SofteningTypeOfPartType[5];

#ifdef BH_BASED_CGM_ZOOM
      All.BlackHolePosition[0] = P[import_indices[n]].Pos[0];
      All.BlackHolePosition[1] = P[import_indices[n]].Pos[1];
      All.BlackHolePosition[2] = P[import_indices[n]].Pos[2];
      All.BlackHoleTask        = ThisTask;
      local_bh_on_task         = ThisTask;
#endif

#ifdef INDIVIDUAL_GRAVITY_SOFTENING
      if(((1 << P[import_indices[n]].Type) & (INDIVIDUAL_GRAVITY_SOFTENING)))
        P[import_indices[n]].SofteningType = get_softening_type_from_mass(P[import_indices[n]].Mass);
#endif

      if(NumBHs + 1 > All.MaxPartBHs)
        {
          sprintf(msg, "On Task=%d with NumPart=%d I tried to make too many BH. Sorry, no space left...(All.MaxPartBHs=%d)\n",
                  ThisTask, NumPart, All.MaxPartBHs);
          terminate(msg);
        }

      /* clear all fields to zero in the new BH structure */
      memset(&BHP[NumBHs], 0, sizeof(struct bh_particle_data));

      P[import_indices[n]].AuxDataID = NumBHs;
#ifdef MEASURE_POTMIN_AROUND_BH
      BHP[NumBHs].BH_MinPot = 0; /* set potential to zero to avoid first repositioning step */
#endif
      BHP[NumBHs].PID = import_indices[n];
      NumBHs++;

#ifdef STELLARAGE
      P[import_indices[n]].StellarAge = All.Time;
#endif

#ifdef UNIFORM_SEEDMASS_DISTRIBUTION
      BPP(import_indices[n]).BH_Mass =  SphP[import_indices[n]].DrawnSeedBlackHoleMass;
#else
      BPP(import_indices[n]).BH_Mass = All.SeedBlackHoleMass;
#endif


      BPP(import_indices[n]).BH_Mdot       = 0;
      BPP(import_indices[n]).BH_CumMass_QM = 0;
      BPP(import_indices[n]).BH_CumEgy_QM  = 0;
      BPP(import_indices[n]).BH_CountProgs = 1;

#ifdef PREVENT_SPURIOUS_RESEEDING
      BPP(import_indices[n]).Time_Of_Seeding = All.Time;
      BPP(import_indices[n]).NeighborsHaveBeenPainted = 1;
#ifdef ACCOUNT_FOR_SWALLOWED_PAINTED_GAS
      BPP(import_indices[n]).SeedMass = 0;
#endif
#endif      

      /* the new blackhole is not active for draining gas now. Mdot is zero anyways */
      BPP(import_indices[n]).SwallowID = P[import_indices[n]].ID;
	
      /* assign Courant timestep of the cell to BH_DtGasNeighbor */

      rad = csnd = 0;

      rad  = get_cell_radius(import_indices[n]);
      csnd = get_sound_speed(import_indices[n]);
#ifdef VORONOI_STATIC_MESH
      csnd +=
          sqrt(P[import_indices[n]].Vel[0] * P[import_indices[n]].Vel[0] + P[import_indices[n]].Vel[1] * P[import_indices[n]].Vel[1] +
               P[import_indices[n]].Vel[2] * P[import_indices[n]].Vel[2]);
#endif

      if(csnd <= 0)
        csnd = 1.0e-30;

      dt_courant = rad / csnd;

#ifdef TREE_BASED_TIMESTEPS
      if(dt_courant > SphP[import_indices[n]].CurrentMaxTiStep)
        dt_courant = SphP[import_indices[n]].CurrentMaxTiStep;
#endif

      BPP(import_indices[n]).BH_DtGasNeighbor = dt_courant;

      BPP(import_indices[n]).BH_Hsml = get_default_softening_of_particletype(5);
#ifdef BH_BUBBLES
      BPP(import_indices[n]).BH_Mass_bubbles = All.SeedBlackHoleMass;
      BPP(import_indices[n]).BH_Mass_ini     = All.SeedBlackHoleMass;
#endif

#ifdef DRAINGAS
      BPP(import_indices[n]).NearestDist = MAX_REAL_NUMBER;
#endif

#ifdef BH_FRICTION
      BPP(import_indices[n]).BH_MinPotTime          = -1;
      BPP(import_indices[n]).BH_MinPotTime_Previous = -1;
      BPP(import_indices[n]).BH_MinPotCumAvgTime    = 0;
#endif

#ifdef BH_NF_RADIO
      BPP(import_indices[n]).BH_RadioEgyFeedback = 0;
      BPP(import_indices[n]).BH_Mdot_quasar      = 0;
      BPP(import_indices[n]).BH_Mdot_radio       = 0;
      BPP(import_indices[n]).BH_HaloVvir         = 0;
      BPP(import_indices[n]).BH_XrayLum          = 0;
      BPP(import_indices[n]).BH_RadioLum         = 0;
      BPP(import_indices[n]).BH_CumMass_QM       = 0;
      BPP(import_indices[n]).BH_CumEgy_QM        = 0;
#endif

#ifdef BH_SPIN_EVOLUTION
      double Phi, CosTheta, SinTheta;

      BPP(import_indices[n]).BH_SpinParameter             = All.BHInitialSpin;
      BPP(import_indices[n]).BlackHoleRadiativeEfficiency = All.BlackHoleRadiativeEfficiency;
      BPP(import_indices[n]).BH_FlagOngAccEpis            = 0;

      Phi                                          = 2 * M_PI * get_random_number();
      CosTheta                                     = 1 - 2 * get_random_number();
      SinTheta                                     = sqrt(1 - CosTheta * CosTheta);
      BPP(import_indices[n]).BH_SpinOrientation[0] = SinTheta * cos(Phi);
      BPP(import_indices[n]).BH_SpinOrientation[1] = SinTheta * sin(Phi);
      BPP(import_indices[n]).BH_SpinOrientation[2] = CosTheta;
#endif

#ifdef USE_SFR
      Stars_converted++;
#endif

#ifdef VORONOI_DYNAMIC_UPDATE
      voronoi_remove_connection(import_indices[n]);
#endif

      timebin_add_particle(&TimeBinsBHAccretion, import_indices[n], -1, P[import_indices[n]].TimeBinHydro,
                           TimeBinSynchronized[P[import_indices[n]].TimeBinHydro]);
    }

#ifdef BH_BASED_CGM_ZOOM
  MPI_Allreduce(&local_bh_on_task, &All.BlackHoleTask, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if(All.BlackHoleTask > -1)
    MPI_Bcast(All.BlackHolePosition, 3, MPI_DOUBLE, All.BlackHoleTask, MPI_COMM_WORLD);
#endif


  /* remove the cells that we converted to black holes from the list of active cells */
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      n = TimeBinsHydro.ActiveParticleList[idx];
      if(n < 0)
        continue;

      if(P[n].Type == 5)
        timebin_remove_particle(&TimeBinsHydro, idx, P[n].TimeBinHydro);
    }

  All.TotNumGas -= ntot;


#ifdef MASSIVE_SEEDS
  blackhole_massiveseeds(nimport, import_indices);
#endif

#ifdef OUTPUT_LOG_FILES_FOR_SEEDING
  myfree(export_parentFOFTask);
  myfree(import_parentFOFTask);

  myfree(export_parentFOFDMMass);
  myfree(import_parentFOFDMMass);
#endif
  myfree(export_indices);
  myfree(import_indices);

}
#endif /* BLACK_HOLES */
#endif

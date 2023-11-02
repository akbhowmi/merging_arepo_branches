/*!
 * \copyright   This file is part of the AREPO code developed by Aklant Kumar Bhowmick.
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


#ifdef CALCULATE_LYMAN_WERNER_INTENSITY_ALL_STARFORMINGGAS
/*! \brief Routine that computes gravitational force by direct summation.
 *
 *  Called by gravity() (in accel.c).
 *
 *  \return void
 */
void calc_lyman_werner_intensity_for_starforminggas(void)
{
  int i, idx;
  MyFloat normalization;

for(idx = 0; idx < NumPart; idx++)
    {
      i = idx;
      if(i < 0)
        continue;

      double fac, wp;
      double dx, dy, dz, r, r2;
      double h, h_inv, h3_inv, u;
      int k;

      /* set softening to corresponding particle's softening length */
      h = All.ForceSoftening[All.SofteningTypeOfPartType[4]];
      P[i].StarFormingGasAllLymanWernerIntensity_type2 = 0; 
      P[i].StarFormingGasAllLymanWernerIntensity_type3 = 0;
      for(k = 0; k < All.TotPartStarFormingGas; k++)
        {
          if(PartLymanWernerStarFormingGasListGlobal[k].ID == P[i].ID)
            continue;
          if(PartLymanWernerStarFormingGasListGlobal[k].ID == 0)
            continue;
          dx = P[i].Pos[0] - PartLymanWernerStarFormingGasListGlobal[k].pos[0];
          dy = P[i].Pos[1] - PartLymanWernerStarFormingGasListGlobal[k].pos[1];
          dz = P[i].Pos[2] - PartLymanWernerStarFormingGasListGlobal[k].pos[2];

          r2 = dx * dx + dy * dy + dz * dz;
          r  = sqrt(r2);

          // using spline softening
          if(r >= h)
            {
              fac = 1 / (r2 * r);
              wp  = -1 / r;
            }
          else
            {
              wp = -1 / r; 
            }

          if (PartLymanWernerStarFormingGasListGlobal[k].metallicity < 0.001 * GFM_SOLAR_METALLICITY) /* This is a Pop III star */
            {
              normalization = 15. * PartLymanWernerStarFormingGasListGlobal[k].mass * 1.0e7 * All.HubbleParam;
              P[i].StarFormingGasAllLymanWernerIntensity_type3 += normalization * wp * wp;
            }
          else if((PartLymanWernerStarFormingGasListGlobal[k].metallicity >= 0.001 * GFM_SOLAR_METALLICITY) 
                  & (PartLymanWernerStarFormingGasListGlobal[k].metallicity < 0.1 * GFM_SOLAR_METALLICITY)) /* This is a Pop II star */
            {
              normalization = 3. * PartLymanWernerStarFormingGasListGlobal[k].mass * 1.0e7 * All.HubbleParam;
              P[i].StarFormingGasAllLymanWernerIntensity_type2 += normalization * wp * wp;
            }
        }
    }
}

//void special_particle_create_list()
void lyman_werner_starforminggas_create_list()
{
  struct lyman_werner_starforminggas_data *LymanWernerPartList;
  int i, j, nsrc, nimport, ngrp;
  All.MaxPartStarFormingGas = 0;
//  for(i = 0, nsrc = 0; i < NumPart; i++)
//    {
//       if((P[i].Type == 0) && (SphP[i].Sfr > 0))
//             All.MaxPartStarFormingGas++;
//    }
  MPI_Allreduce(&All.MaxPartSph, &All.TotPartStarFormingGas, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  LymanWernerPartList =
      (struct lyman_werner_starforminggas_data *)mymalloc("LymanWernerPartList", All.TotPartStarFormingGas * sizeof(struct lyman_werner_starforminggas_data));

  for(i = 0, nsrc = 0; i < NumPart; i++)
    {
      if((P[i].Type == 0) && (SphP[i].Sfr > 0))
        {
          LymanWernerPartList[nsrc].ID = P[i].ID;
          LymanWernerPartList[nsrc].pos[0] = P[i].Pos[0];
          LymanWernerPartList[nsrc].pos[1] = P[i].Pos[1];
          LymanWernerPartList[nsrc].pos[2] = P[i].Pos[2];
          LymanWernerPartList[nsrc].mass = SphP[i].Sfr * 5e6 * 1e-10 * All.HubbleParam;
          LymanWernerPartList[nsrc++].metallicity = SphP[i].Metallicity;
        }
    }

  for(j = 0; j < NTask; j++)
    Send_count[j] = nsrc;

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = 0;
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  /* exchange particle data */
  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&LymanWernerPartList[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct lyman_werner_starforminggas_data),
                           MPI_BYTE, recvTask, TAG_LYMAN_B, &PartLymanWernerStarFormingGasListGlobal[Recv_offset[recvTask]],
                           Recv_count[recvTask] * sizeof(struct lyman_werner_starforminggas_data), MPI_BYTE, recvTask, TAG_LYMAN_B, MPI_COMM_WORLD,
                           MPI_STATUS_IGNORE);
            }
        }
    }

  myfree(LymanWernerPartList);
}

/*! \brief Updates list of special particles, i.e. particles for which gravity
 *  is calculated by direct summation.
 *
 *  Called in run() (run.c).
 *
 *  \return void
 */
void lyman_werner_starforminggas_update_list()
{
  struct lyman_werner_starforminggas_data *LymanWernerPartList;

  int i, j, nsrc, nimport, ngrp;
//  All.MaxPartStarFormingGas = 0;
//  for(i = 0, nsrc = 0; i < NumPart; i++)
 //   {
 //      if((P[i].Type == 0) && (SphP[i].Sfr > 0))
 //            All.MaxPartStarFormingGas++;
//    }
  MPI_Allreduce(&All.MaxPartSph, &All.TotPartStarFormingGas, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  LymanWernerPartList =
      (struct lyman_werner_starforminggas_data *)mymalloc("LymanWernerPartList", All.TotPartStarFormingGas * sizeof(struct lyman_werner_starforminggas_data));

  for(i = 0, nsrc = 0; i < NumPart; i++)
    {
      if((P[i].Type == 0) && (SphP[i].Sfr > 0))
        {
          LymanWernerPartList[nsrc].ID = P[i].ID;
          LymanWernerPartList[nsrc].pos[0] = P[i].Pos[0];
          LymanWernerPartList[nsrc].pos[1] = P[i].Pos[1];
          LymanWernerPartList[nsrc].pos[2] = P[i].Pos[2];
          LymanWernerPartList[nsrc].mass = SphP[i].Sfr * 5e6 * 1e-10 * All.HubbleParam;
          LymanWernerPartList[nsrc++].metallicity = STP(i).Metallicity;

        }
    }

  for(j = 0; j < NTask; j++)
    Send_count[j] = nsrc;

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
    {
      nimport += Recv_count[j];

      if(j > 0)
        {
          Send_offset[j] = 0;
          Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
        }
    }

  /* exchange particle data */
  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
        {
          if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
            {
              /* get the particles */
              MPI_Sendrecv(&LymanWernerPartList[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct lyman_werner_starforminggas_data),
                           MPI_BYTE, recvTask, TAG_LYMAN_B, &PartLymanWernerStarFormingGasListGlobal[Recv_offset[recvTask]],
                           Recv_count[recvTask] * sizeof(struct lyman_werner_starforminggas_data), MPI_BYTE, recvTask, TAG_LYMAN_B, MPI_COMM_WORLD,
                           MPI_STATUS_IGNORE);
            }
        }
    }
  myfree(LymanWernerPartList);
}
#endif




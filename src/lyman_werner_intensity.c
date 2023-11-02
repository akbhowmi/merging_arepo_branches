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


#ifdef CALCULATE_LYMAN_WERNER_INTENSITY_ALL_SOURCES
/*! \brief Routine that computes gravitational force by direct summation.
 *
 *  Called by gravity() (in accel.c).
 *
 *  \return void
 */
void calc_lyman_werner_intensity_for_stars(void)
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
      P[i].StellarAllLymanWernerIntensity_type2 = 0; 
      P[i].StellarAllLymanWernerIntensity_type3 = 0;
      for(k = 0; k < All.TotPartStar; k++)
        {
          if(PartLymanWernerListGlobal[k].ID == P[i].ID)
            continue;
          if(PartLymanWernerListGlobal[k].ID == 0)
            continue;
          if(PartLymanWernerListGlobal[k].age > 0.005)  /*Only stars with age less than 5 Myr should contribute*/ 
            continue;

          dx = P[i].Pos[0] - PartLymanWernerListGlobal[k].pos[0];
          dy = P[i].Pos[1] - PartLymanWernerListGlobal[k].pos[1];
          dz = P[i].Pos[2] - PartLymanWernerListGlobal[k].pos[2];

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

          if (PartLymanWernerListGlobal[k].metallicity < 0.001 * GFM_SOLAR_METALLICITY) /* This is a Pop III star */
            {
              normalization = 15. * PartLymanWernerListGlobal[k].mass * 1.0e7 * All.HubbleParam;
              P[i].StellarAllLymanWernerIntensity_type3 += normalization * wp * wp;
            }
          else if((PartLymanWernerListGlobal[k].metallicity >= 0.001 * GFM_SOLAR_METALLICITY) 
                  & (PartLymanWernerListGlobal[k].metallicity < 0.1 * GFM_SOLAR_METALLICITY)) /* This is a Pop II star */
            {
              normalization = 3. * PartLymanWernerListGlobal[k].mass * 1.0e7 * All.HubbleParam;
              P[i].StellarAllLymanWernerIntensity_type2 += normalization * wp * wp;
            }
        }
    }
}

//void special_particle_create_list()
void lyman_werner_source_create_list()
{
  int tot_stars;
  struct lyman_werner_particle_data *LymanWernerPartList;
  MPI_Allreduce(&All.MaxPartStar, &tot_stars, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  LymanWernerPartList =
      (struct lyman_werner_particle_data *)mymalloc("LymanWernerPartList", tot_stars * sizeof(struct lyman_werner_particle_data));

  int i, j, nsrc, nimport, ngrp;
  for(i = 0, nsrc = 0; i < NumPart; i++)
    {
      if((P[i].Type == 4) && (STP(i).BirthTime > 0))
        {
          LymanWernerPartList[nsrc].ID = P[i].ID;
          LymanWernerPartList[nsrc].pos[0] = P[i].Pos[0];
          LymanWernerPartList[nsrc].pos[1] = P[i].Pos[1];
          LymanWernerPartList[nsrc].pos[2] = P[i].Pos[2];
          LymanWernerPartList[nsrc].mass = P[i].Mass;
          LymanWernerPartList[nsrc].age = STP(i).StellarAgeGyr;
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
              MPI_Sendrecv(&LymanWernerPartList[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct lyman_werner_particle_data),
                           MPI_BYTE, recvTask, TAG_LYMAN_A, &PartLymanWernerListGlobal[Recv_offset[recvTask]],
                           Recv_count[recvTask] * sizeof(struct lyman_werner_particle_data), MPI_BYTE, recvTask, TAG_LYMAN_A, MPI_COMM_WORLD,
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
void lyman_werner_source_update_list()
{
  int tot_stars;
  struct lyman_werner_particle_data *LymanWernerPartList;
  MPI_Allreduce(&All.MaxPartStar, &tot_stars, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  LymanWernerPartList =
      (struct lyman_werner_particle_data *)mymalloc("LymanWernerPartList", tot_stars * sizeof(struct lyman_werner_particle_data));

  int i, j, nsrc, nimport, ngrp;
  for(i = 0, nsrc = 0; i < NumPart; i++)
    {
      if((P[i].Type == 4) && (STP(i).BirthTime > 0))
        {
          LymanWernerPartList[nsrc].ID = P[i].ID;
          LymanWernerPartList[nsrc].pos[0] = P[i].Pos[0];
          LymanWernerPartList[nsrc].pos[1] = P[i].Pos[1];
          LymanWernerPartList[nsrc].pos[2] = P[i].Pos[2];
          LymanWernerPartList[nsrc].mass = P[i].Mass;
          LymanWernerPartList[nsrc].age = STP(i).StellarAgeGyr;
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
              MPI_Sendrecv(&LymanWernerPartList[Send_offset[recvTask]], Send_count[recvTask] * sizeof(struct lyman_werner_particle_data),
                           MPI_BYTE, recvTask, TAG_LYMAN_A, &PartLymanWernerListGlobal[Recv_offset[recvTask]],
                           Recv_count[recvTask] * sizeof(struct lyman_werner_particle_data), MPI_BYTE, recvTask, TAG_LYMAN_A, MPI_COMM_WORLD,
                           MPI_STATUS_IGNORE);
            }
        }
    }
  myfree(LymanWernerPartList);
}
#endif




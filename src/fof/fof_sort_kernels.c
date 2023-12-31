/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/fof/fof_sort_kernels.c
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

#include "../allvars.h"
#include "../domain.h"
#include "../proto.h"
#include "../subfind/subfind.h"
#include "fof.h"

/*! \file fof_sort_kernels.c
 *  \brief here are various sort kernels used by the parallel FoF group finder
 */

#ifdef FOF

/*! \brief Comparison function for fof_local_sort_data objects.
 *
 *  Sorting kernel comparing element targetindex.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b.
 */
int fof_compare_local_sort_data_targetindex(const void *a, const void *b)
{
  if(((struct fof_local_sort_data *)a)->targetindex < ((struct fof_local_sort_data *)b)->targetindex)
    return -1;

  if(((struct fof_local_sort_data *)a)->targetindex > ((struct fof_local_sort_data *)b)->targetindex)
    return +1;

  return 0;
}

/*! \brief Comparison function for data_aux_sort objects.
 *
 *  Sorting kernel comparing element Type.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b.
 */
int fof_compare_aux_sort_Type(const void *a, const void *b)
{
  if(((struct data_aux_sort *)a)->Type < ((struct data_aux_sort *)b)->Type)
    return -1;

  if(((struct data_aux_sort *)a)->Type > ((struct data_aux_sort *)b)->Type)
    return +1;

  return 0;
}

#if defined(RECOMPUTE_POTENTIAL_IN_SNAPSHOT) || defined(CALCULATE_QUANTITIES_IN_POSTPROCESS)
/*! \brief Comparison function for data_aux_sort objects.
 *
 *  Sorting kernel comparing element FileOrder.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b.
 */
int fof_compare_aux_sort_FileOrder(const void *a, const void *b)
{
  if(((struct data_aux_sort *)a)->FileOrder < ((struct data_aux_sort *)b)->FileOrder)
    return -1;

  if(((struct data_aux_sort *)a)->FileOrder > ((struct data_aux_sort *)b)->FileOrder)
    return +1;

  return 0;
}
#endif

/*! \brief Comparison function for data_aux_sort objects.
 *
 *  Sorting kernel comparing elements (most important fist):
 *   GrNr, SubNr, DM_BindingEnergy
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b.
 */
int fof_compare_aux_sort_GrNr(const void *a, const void *b)
{
  if(((struct data_aux_sort *)a)->GrNr < ((struct data_aux_sort *)b)->GrNr)
    return -1;

  if(((struct data_aux_sort *)a)->GrNr > ((struct data_aux_sort *)b)->GrNr)
    return +1;

#ifdef SUBFIND
  if(((struct data_aux_sort *)a)->SubNr < ((struct data_aux_sort *)b)->SubNr)
    return -1;

  if(((struct data_aux_sort *)a)->SubNr > ((struct data_aux_sort *)b)->SubNr)
    return +1;

  if(((struct data_aux_sort *)a)->DM_BindingEnergy < ((struct data_aux_sort *)b)->DM_BindingEnergy)
    return -1;

  if(((struct data_aux_sort *)a)->DM_BindingEnergy > ((struct data_aux_sort *)b)->DM_BindingEnergy)
    return +1;
#endif

#if defined(FOF_FUZZ_SORT_BY_NEAREST_GROUP)
  if(((struct data_aux_sort *)a)->key < ((struct data_aux_sort *)b)->key)
    return -1;

  if(((struct data_aux_sort *)a)->key > ((struct data_aux_sort *)b)->key)
    return +1;
#else
  if(((struct data_aux_sort *)a)->ID < ((struct data_aux_sort *)b)->ID)
    return -1;

  if(((struct data_aux_sort *)a)->ID > ((struct data_aux_sort *)b)->ID)
    return +1;
#endif

  return 0;
}

/*! \brief Comparison function for data_aux_sort objects.
 *
 *  Sorting kernel comparing elements (most important first):
 *   OriginTask, OriginIndex
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b.
 */
int fof_compare_aux_sort_OriginTask_OriginIndex(const void *a, const void *b)
{
  if(((struct data_aux_sort *)a)->OriginTask < ((struct data_aux_sort *)b)->OriginTask)
    return -1;

  if(((struct data_aux_sort *)a)->OriginTask > ((struct data_aux_sort *)b)->OriginTask)
    return +1;

  if(((struct data_aux_sort *)a)->OriginIndex < ((struct data_aux_sort *)b)->OriginIndex)
    return -1;

  if(((struct data_aux_sort *)a)->OriginIndex > ((struct data_aux_sort *)b)->OriginIndex)
    return +1;

  return 0;
}

/*! \brief Comparison function for fof_particle_list objects.
 *
 *  Sorting kernel comparing element MinID.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b.
 */
int fof_compare_FOF_PList_MinID(const void *a, const void *b)
{
  if(((struct fof_particle_list *)a)->MinID < ((struct fof_particle_list *)b)->MinID)
    return -1;

  if(((struct fof_particle_list *)a)->MinID > ((struct fof_particle_list *)b)->MinID)
    return +1;

  return 0;
}

/*! \brief Comparison function for fof_group_list objects.
 *
 *  Sorting kernel comparing element MinID.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b.
 */
int fof_compare_FOF_GList_MinID(const void *a, const void *b)
{
  if(((struct fof_group_list *)a)->MinID < ((struct fof_group_list *)b)->MinID)
    return -1;

  if(((struct fof_group_list *)a)->MinID > ((struct fof_group_list *)b)->MinID)
    return +1;

  return 0;
}

/*! \brief Comparison function for fof_group_list objects.
 *
 *  Sorting kernel comparing element MinIDTask.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b.
 */
int fof_compare_FOF_GList_MinIDTask(const void *a, const void *b)
{
  if(((struct fof_group_list *)a)->MinIDTask < ((struct fof_group_list *)b)->MinIDTask)
    return -1;

  if(((struct fof_group_list *)a)->MinIDTask > ((struct fof_group_list *)b)->MinIDTask)
    return +1;

  return 0;
}

/*! \brief Comparison function for fof_group_list objects.
 *
 *  Sorting kernel comparing elements (most important first):
 *   MinIDTask, MinID.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b.
 */
int fof_compare_FOF_GList_MinIDTask_MinID(const void *a, const void *b)
{
  if(((struct fof_group_list *)a)->MinIDTask < ((struct fof_group_list *)b)->MinIDTask)
    return -1;

  if(((struct fof_group_list *)a)->MinIDTask > ((struct fof_group_list *)b)->MinIDTask)
    return +1;

  if(((struct fof_group_list *)a)->MinID < ((struct fof_group_list *)b)->MinID)
    return -1;

  if(((struct fof_group_list *)a)->MinID > ((struct fof_group_list *)b)->MinID)
    return +1;

  return 0;
}

/*! \brief Comparison function for fof_group_list objects.
 *
 *  Sorting kernel comparing elements (most important first):.
 *   LocCount, MinID, ExtCount.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b, except for LocCount where -1 if a > b
 */
int fof_compare_FOF_GList_LocCountTaskDiffMinID(const void *a, const void *b)
{
  if(((struct fof_group_list *)a)->LocCount > ((struct fof_group_list *)b)->LocCount)
    return -1;

  if(((struct fof_group_list *)a)->LocCount < ((struct fof_group_list *)b)->LocCount)
    return +1;

#ifdef ADD_GROUP_PROPERTIES
  if(((struct fof_group_list *)a)->OriginalGrNr < ((struct fof_group_list *)b)->OriginalGrNr)
    return -1;

  if(((struct fof_group_list *)a)->OriginalGrNr > ((struct fof_group_list *)b)->OriginalGrNr)
    return +1;
#endif

  if(((struct fof_group_list *)a)->MinID < ((struct fof_group_list *)b)->MinID)
    return -1;

  if(((struct fof_group_list *)a)->MinID > ((struct fof_group_list *)b)->MinID)
    return +1;

  if(labs(((struct fof_group_list *)a)->ExtCount - ((struct fof_group_list *)a)->MinIDTask) <
     labs(((struct fof_group_list *)b)->ExtCount - ((struct fof_group_list *)b)->MinIDTask))
    return -1;

  if(labs(((struct fof_group_list *)a)->ExtCount - ((struct fof_group_list *)a)->MinIDTask) >
     labs(((struct fof_group_list *)b)->ExtCount - ((struct fof_group_list *)b)->MinIDTask))
    return +1;

  return 0;
}

/*! \brief Comparison function for fof_group_list objects.
 *
 *  Sorting kernel comparing elements (most important first):
 *   ExtCount, MinID.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b.
 */
int fof_compare_FOF_GList_ExtCountMinID(const void *a, const void *b)
{
  if(((struct fof_group_list *)a)->ExtCount < ((struct fof_group_list *)b)->ExtCount)
    return -1;

  if(((struct fof_group_list *)a)->ExtCount > ((struct fof_group_list *)b)->ExtCount)
    return +1;

  if(((struct fof_group_list *)a)->MinID < ((struct fof_group_list *)b)->MinID)
    return -1;

  if(((struct fof_group_list *)a)->MinID > ((struct fof_group_list *)b)->MinID)
    return +1;

  return 0;
}

/*! \brief Comparison function for group_properties objects.
 *
 *  Sorting kernel comparing element MinID.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b.
 */
int fof_compare_Group_MinID(const void *a, const void *b)
{
  if(((struct group_properties *)a)->MinID < ((struct group_properties *)b)->MinID)
    return -1;

  if(((struct group_properties *)a)->MinID > ((struct group_properties *)b)->MinID)
    return +1;

  return 0;
}

/*! \brief Comparison function for group_properties objects.
 *
 *  Sorting kernel comparing element GrNr.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b.
 */
int fof_compare_Group_GrNr(const void *a, const void *b)
{
  if(((struct group_properties *)a)->GrNr < ((struct group_properties *)b)->GrNr)
    return -1;

  if(((struct group_properties *)a)->GrNr > ((struct group_properties *)b)->GrNr)
    return +1;

  return 0;
}

/*! \brief Comparison function for group_properties objects.
 *
 *  Sorting kernel comparing element MinIDTask.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b.
 */
int fof_compare_Group_MinIDTask(const void *a, const void *b)
{
  if(((struct group_properties *)a)->MinIDTask < ((struct group_properties *)b)->MinIDTask)
    return -1;

  if(((struct group_properties *)a)->MinIDTask > ((struct group_properties *)b)->MinIDTask)
    return +1;

  return 0;
}

/*! \brief Comparison function for group_properties objects.
 *
 *  Sorting kernel comparing elements (most imporant first):
 *   MinIDTask, MinID.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b.
 */
int fof_compare_Group_MinIDTask_MinID(const void *a, const void *b)
{
  if(((struct group_properties *)a)->MinIDTask < ((struct group_properties *)b)->MinIDTask)
    return -1;

  if(((struct group_properties *)a)->MinIDTask > ((struct group_properties *)b)->MinIDTask)
    return +1;

  if(((struct group_properties *)a)->MinID < ((struct group_properties *)b)->MinID)
    return -1;

  if(((struct group_properties *)a)->MinID > ((struct group_properties *)b)->MinID)
    return +1;

  return 0;
}

/*! \brief Comparison function for id_list objects.
 *
 *  Sorting kernel comparing elements (most important first):
 *   GrNr, Type, ID.
 *
 *  \param[in] a First object to compare.
 *  \param[in] b Second object to compare.
 *
 *  \return (-1,0,1), -1 if a < b.
 */
int fof_compare_ID_list_GrNrID(const void *a, const void *b)
{
  if(((struct id_list *)a)->GrNr < ((struct id_list *)b)->GrNr)
    return -1;

  if(((struct id_list *)a)->GrNr > ((struct id_list *)b)->GrNr)
    return +1;

  if(((struct id_list *)a)->Type < ((struct id_list *)b)->Type)
    return -1;

  if(((struct id_list *)a)->Type > ((struct id_list *)b)->Type)
    return +1;

  if(((struct id_list *)a)->ID < ((struct id_list *)b)->ID)
    return -1;

  if(((struct id_list *)a)->ID > ((struct id_list *)b)->ID)
    return +1;

  return 0;
}

#endif /* of FOF */

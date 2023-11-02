/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/fof/fof_vars.c
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

#include "../../allvars.h"
#include "mergers_io.h"

/*! \file merger_vars.c
 *  \brief In this file we create instances for the global variables used by FOF, which are declared in fof.h
 */

#ifdef STORE_MERGERS_IN_SNAPSHOT

int Nmergers=0, TotNmergers;
struct merger_properties *MergerEvents;

struct merger_properties *MergerEventsAll;
struct merger_catalogue *MergerCat, *MergerCatLocal, *MergerCatSend;

struct merger_header merger_catalogue_header;

#endif 

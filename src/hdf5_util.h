/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/hdf5_util.h
 * \date        11/2019
 * \author      Simon May
 * \brief
 * \details
 *
 */

#ifndef HDF5_UTIL_H
#define HDF5_UTIL_H

#include <stddef.h>

#ifdef HAVE_HDF5
#include <hdf5.h>

hid_t my_H5Fcreate(const char *fname, unsigned flags, hid_t fcpl_id, hid_t fapl_id);
hid_t my_H5Gcreate(hid_t loc_id, const char *groupname, size_t size_hint);
hid_t my_H5Dcreate(hid_t loc_id, const char *datasetname, hid_t type_id, hid_t space_id, hid_t dcpl_id);
hid_t my_H5Acreate(hid_t loc_id, const char *attr_name, hid_t type_id, hid_t space_id, hid_t acpl_id);
hid_t my_H5Screate(H5S_class_t type);
hid_t my_H5Screate_simple(int rank, const hsize_t *current_dims, const hsize_t *maximum_dims);
herr_t my_H5Dwrite(hid_t dataset_id, hid_t mem_type_id, hid_t mem_space_id, hid_t file_space_id, hid_t xfer_plist_id, const void *buf,
                   const char *datasetname);
herr_t my_H5Awrite(hid_t attr_id, hid_t mem_type_id, const void *buf, const char *attr_name);
hid_t my_H5Fopen(const char *fname, unsigned int flags, hid_t fapl_id);
hid_t my_H5Dopen(hid_t file_id, const char *datasetname);
hid_t my_H5Dopen2(hid_t file_id, const char *datasetname, hid_t dapl_id);
hid_t my_H5Dopen_if_existing(hid_t file_id, const char *datasetname);
herr_t my_H5Dread(hid_t dataset_id, hid_t mem_type_id, hid_t mem_space_id, hid_t file_space_id, hid_t xfer_plist_id, void *buf,
                  const char *datasetname);
hid_t my_H5Gopen(hid_t loc_id, const char *groupname);
hid_t my_H5Aopen_name(hid_t loc_id, const char *attr_name);
herr_t my_H5Aread(hid_t attr_id, hid_t mem_type_id, void *buf, const char *attr_name, hssize_t size);

herr_t my_H5Aclose(hid_t attr_id, const char *attr_name);
herr_t my_H5Dclose(hid_t dataset_id, const char *datasetname);
herr_t my_H5Gclose(hid_t group_id, const char *groupname);
herr_t my_H5Fclose(hid_t file_id, const char *fname);
herr_t my_H5Sclose(hid_t dataspace_id, H5S_class_t type);

hid_t my_H5Tcopy(hid_t type_id);
herr_t my_H5Tclose(hid_t type_id);

herr_t my_H5Sselect_hyperslab(hid_t space_id, H5S_seloper_t op, const hsize_t *start, const hsize_t *stride, const hsize_t *count,
                              const hsize_t *block);
size_t my_H5Tget_size(hid_t datatype_id);
herr_t my_H5Tset_size(hid_t datatype_id, size_t size);

herr_t my_H5Sset_extent_simple(hid_t space_id, int rank, const hsize_t *current_size, const hsize_t *maximum_size,
                               const char *attr_name);
hid_t my_H5Dget_space(hid_t dataset_id, const char *datasetname);

#ifdef HDF5_FILTERS
htri_t my_H5Pall_filters_avail(hid_t plist_id);
hid_t my_H5Pcreate(hid_t class_id);
herr_t my_H5Pclose(hid_t plist);
herr_t my_H5Pset_chunk(hid_t plist, int ndims, const hsize_t *dim);
herr_t my_H5Pset_shuffle(hid_t plist_id);
herr_t my_H5Pset_deflate(hid_t plist_id, uint level);
herr_t my_H5Pset_fletcher32(hid_t plist_id);
#endif

#endif

#endif

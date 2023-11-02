#include <stdio.h>
#include "arepoconfig.h"
#ifdef HAVE_HDF5
#include <hdf5.h>
#include "hdf5_util.h"

void write_compile_time_options_in_hdf5(hid_t handle)
{
hid_t hdf5_dataspace, hdf5_attribute;
double val;
hid_t atype = H5Tcopy(H5T_C_S1);
H5Tset_size(atype, 1);
hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "NTYPES", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 6;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "NTYPES");
my_H5Aclose(hdf5_attribute, "NTYPES");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "COOLING", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "COOLING");
my_H5Aclose(hdf5_attribute, "COOLING");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "UVB_SELF_SHIELDING", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "UVB_SELF_SHIELDING");
my_H5Aclose(hdf5_attribute, "UVB_SELF_SHIELDING");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "USE_SFR", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "USE_SFR");
my_H5Aclose(hdf5_attribute, "USE_SFR");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "VORONOI", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "VORONOI");
my_H5Aclose(hdf5_attribute, "VORONOI");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "REGULARIZE_MESH_CM_DRIFT", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "REGULARIZE_MESH_CM_DRIFT");
my_H5Aclose(hdf5_attribute, "REGULARIZE_MESH_CM_DRIFT");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED");
my_H5Aclose(hdf5_attribute, "REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "REGULARIZE_MESH_FACE_ANGLE", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "REGULARIZE_MESH_FACE_ANGLE");
my_H5Aclose(hdf5_attribute, "REGULARIZE_MESH_FACE_ANGLE");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "TREE_BASED_TIMESTEPS", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "TREE_BASED_TIMESTEPS");
my_H5Aclose(hdf5_attribute, "TREE_BASED_TIMESTEPS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "REFINEMENT_SPLIT_CELLS", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "REFINEMENT_SPLIT_CELLS");
my_H5Aclose(hdf5_attribute, "REFINEMENT_SPLIT_CELLS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "REFINEMENT_MERGE_CELLS", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "REFINEMENT_MERGE_CELLS");
my_H5Aclose(hdf5_attribute, "REFINEMENT_MERGE_CELLS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "NODEREFINE_BACKGROUND_GRID", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "NODEREFINE_BACKGROUND_GRID");
my_H5Aclose(hdf5_attribute, "NODEREFINE_BACKGROUND_GRID");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "BH_BASED_CGM_ZOOM", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "BH_BASED_CGM_ZOOM");
my_H5Aclose(hdf5_attribute, "BH_BASED_CGM_ZOOM");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SELFGRAVITY", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SELFGRAVITY");
my_H5Aclose(hdf5_attribute, "SELFGRAVITY");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GRAVITY_NOT_PERIODIC", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "GRAVITY_NOT_PERIODIC");
my_H5Aclose(hdf5_attribute, "GRAVITY_NOT_PERIODIC");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "EVALPOTENTIAL", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "EVALPOTENTIAL");
my_H5Aclose(hdf5_attribute, "EVALPOTENTIAL");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "ADAPTIVE_HYDRO_SOFTENING", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "ADAPTIVE_HYDRO_SOFTENING");
my_H5Aclose(hdf5_attribute, "ADAPTIVE_HYDRO_SOFTENING");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "DOUBLEPRECISION", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 1;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "DOUBLEPRECISION");
my_H5Aclose(hdf5_attribute, "DOUBLEPRECISION");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "DOUBLEPRECISION_FFTW", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "DOUBLEPRECISION_FFTW");
my_H5Aclose(hdf5_attribute, "DOUBLEPRECISION_FFTW");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "BLACK_HOLES", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "BLACK_HOLES");
my_H5Aclose(hdf5_attribute, "BLACK_HOLES");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "DRAINGAS", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 3;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "DRAINGAS");
my_H5Aclose(hdf5_attribute, "DRAINGAS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "BH_EXACT_INTEGRATION", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "BH_EXACT_INTEGRATION");
my_H5Aclose(hdf5_attribute, "BH_EXACT_INTEGRATION");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "BH_BONDI_DEFAULT", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "BH_BONDI_DEFAULT");
my_H5Aclose(hdf5_attribute, "BH_BONDI_DEFAULT");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "BH_DO_NOT_PREVENT_MERGERS", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "BH_DO_NOT_PREVENT_MERGERS");
my_H5Aclose(hdf5_attribute, "BH_DO_NOT_PREVENT_MERGERS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "BH_NEW_LOGS", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "BH_NEW_LOGS");
my_H5Aclose(hdf5_attribute, "BH_NEW_LOGS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "BH_DF_DISCRETE", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "BH_DF_DISCRETE");
my_H5Aclose(hdf5_attribute, "BH_DF_DISCRETE");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "BH_FAST_WIND", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "BH_FAST_WIND");
my_H5Aclose(hdf5_attribute, "BH_FAST_WIND");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "BH_FAST_WIND_STOCHASTIC", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "BH_FAST_WIND_STOCHASTIC");
my_H5Aclose(hdf5_attribute, "BH_FAST_WIND_STOCHASTIC");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "VORONOI_DYNAMIC_UPDATE", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "VORONOI_DYNAMIC_UPDATE");
my_H5Aclose(hdf5_attribute, "VORONOI_DYNAMIC_UPDATE");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "NO_MPI_IN_PLACE", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "NO_MPI_IN_PLACE");
my_H5Aclose(hdf5_attribute, "NO_MPI_IN_PLACE");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "NO_ISEND_IRECV_IN_DOMAIN", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "NO_ISEND_IRECV_IN_DOMAIN");
my_H5Aclose(hdf5_attribute, "NO_ISEND_IRECV_IN_DOMAIN");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG");
my_H5Aclose(hdf5_attribute, "FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "ENLARGE_DYNAMIC_RANGE_IN_TIME", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "ENLARGE_DYNAMIC_RANGE_IN_TIME");
my_H5Aclose(hdf5_attribute, "ENLARGE_DYNAMIC_RANGE_IN_TIME");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "OUTPUT_PRESSURE", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "OUTPUT_PRESSURE");
my_H5Aclose(hdf5_attribute, "OUTPUT_PRESSURE");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "OUTPUTACCELERATION", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "OUTPUTACCELERATION");
my_H5Aclose(hdf5_attribute, "OUTPUTACCELERATION");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "HAVE_HDF5", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "HAVE_HDF5");
my_H5Aclose(hdf5_attribute, "HAVE_HDF5");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "OUTPUT_COOLHEAT", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "OUTPUT_COOLHEAT");
my_H5Aclose(hdf5_attribute, "OUTPUT_COOLHEAT");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "DEBUG", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "DEBUG");
my_H5Aclose(hdf5_attribute, "DEBUG");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "HOST_MEMORY_REPORTING", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "HOST_MEMORY_REPORTING");
my_H5Aclose(hdf5_attribute, "HOST_MEMORY_REPORTING");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "GFM");
my_H5Aclose(hdf5_attribute, "GFM");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM_STELLAR_EVOLUTION", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 0;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "GFM_STELLAR_EVOLUTION");
my_H5Aclose(hdf5_attribute, "GFM_STELLAR_EVOLUTION");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM_SET_METALLICITY", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "GFM_SET_METALLICITY");
my_H5Aclose(hdf5_attribute, "GFM_SET_METALLICITY");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM_NO_METAL_ENRICHMENT", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "GFM_NO_METAL_ENRICHMENT");
my_H5Aclose(hdf5_attribute, "GFM_NO_METAL_ENRICHMENT");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM_COOLING_METAL", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "GFM_COOLING_METAL");
my_H5Aclose(hdf5_attribute, "GFM_COOLING_METAL");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM_OUTPUT_MASK", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 1+2+4+8+16+32+64+128;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "GFM_OUTPUT_MASK");
my_H5Aclose(hdf5_attribute, "GFM_OUTPUT_MASK");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM_NORMALIZED_METAL_ADVECTION", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "GFM_NORMALIZED_METAL_ADVECTION");
my_H5Aclose(hdf5_attribute, "GFM_NORMALIZED_METAL_ADVECTION");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM_DISCRETE_ENRICHMENT", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "GFM_DISCRETE_ENRICHMENT");
my_H5Aclose(hdf5_attribute, "GFM_DISCRETE_ENRICHMENT");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SMUGGLE_SFR", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SMUGGLE_SFR");
my_H5Aclose(hdf5_attribute, "SMUGGLE_SFR");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SMUGGLE_STAR_FEEDBACK", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SMUGGLE_STAR_FEEDBACK");
my_H5Aclose(hdf5_attribute, "SMUGGLE_STAR_FEEDBACK");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SMUGGLE_STAR_FEEDBACK_TIME_LIMITER", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SMUGGLE_STAR_FEEDBACK_TIME_LIMITER");
my_H5Aclose(hdf5_attribute, "SMUGGLE_STAR_FEEDBACK_TIME_LIMITER");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SMUGGLE_COMPUTE_SFR_FROM_H2", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SMUGGLE_COMPUTE_SFR_FROM_H2");
my_H5Aclose(hdf5_attribute, "SMUGGLE_COMPUTE_SFR_FROM_H2");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SMUGGLE_OUTPUT_STELLAR_FEEDBACK", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SMUGGLE_OUTPUT_STELLAR_FEEDBACK");
my_H5Aclose(hdf5_attribute, "SMUGGLE_OUTPUT_STELLAR_FEEDBACK");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SMUGGLE_OUTPUT_MOLECULAR_FRACTION", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SMUGGLE_OUTPUT_MOLECULAR_FRACTION");
my_H5Aclose(hdf5_attribute, "SMUGGLE_OUTPUT_MOLECULAR_FRACTION");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SMUGGLE_OUTPUT_VIRIAL_PARAM", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SMUGGLE_OUTPUT_VIRIAL_PARAM");
my_H5Aclose(hdf5_attribute, "SMUGGLE_OUTPUT_VIRIAL_PARAM");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SMUGGLE_RADIATION_FEEDBACK", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SMUGGLE_RADIATION_FEEDBACK");
my_H5Aclose(hdf5_attribute, "SMUGGLE_RADIATION_FEEDBACK");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SMUGGLE_RADIATION_FEEDBACK_DEBUG", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SMUGGLE_RADIATION_FEEDBACK_DEBUG");
my_H5Aclose(hdf5_attribute, "SMUGGLE_RADIATION_FEEDBACK_DEBUG");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SMUGGLE_OMEGA_WEIGHT_SN", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SMUGGLE_OMEGA_WEIGHT_SN");
my_H5Aclose(hdf5_attribute, "SMUGGLE_OMEGA_WEIGHT_SN");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SMUGGLE_MOLEC_COOLING", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SMUGGLE_MOLEC_COOLING");
my_H5Aclose(hdf5_attribute, "SMUGGLE_MOLEC_COOLING");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SMUGGLE_COSMIC_RAY_HEATING", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SMUGGLE_COSMIC_RAY_HEATING");
my_H5Aclose(hdf5_attribute, "SMUGGLE_COSMIC_RAY_HEATING");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SMUGGLE_PHOTOELECTRIC_HEATING", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SMUGGLE_PHOTOELECTRIC_HEATING");
my_H5Aclose(hdf5_attribute, "SMUGGLE_PHOTOELECTRIC_HEATING");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SMUGGLE_SN_COOLING_RADIUS_BOOST", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SMUGGLE_SN_COOLING_RADIUS_BOOST");
my_H5Aclose(hdf5_attribute, "SMUGGLE_SN_COOLING_RADIUS_BOOST");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SMUGGLE_DISCRETE_SN", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SMUGGLE_DISCRETE_SN");
my_H5Aclose(hdf5_attribute, "SMUGGLE_DISCRETE_SN");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SMUGGLE_AGB_WINDS", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SMUGGLE_AGB_WINDS");
my_H5Aclose(hdf5_attribute, "SMUGGLE_AGB_WINDS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SMUGGLE_STOCHASTIC_HII_PHOTOIONIZATION", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 0;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "SMUGGLE_STOCHASTIC_HII_PHOTOIONIZATION");
my_H5Aclose(hdf5_attribute, "SMUGGLE_STOCHASTIC_HII_PHOTOIONIZATION");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

my_H5Tclose(atype);
}
#endif

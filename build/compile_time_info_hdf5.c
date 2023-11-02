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
hdf5_attribute = my_H5Acreate(handle, "PERIODIC", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "PERIODIC");
my_H5Aclose(hdf5_attribute, "PERIODIC");
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
hdf5_attribute = my_H5Acreate(handle, "MHD", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "MHD");
my_H5Aclose(hdf5_attribute, "MHD");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "MHD_POWELL", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "MHD_POWELL");
my_H5Aclose(hdf5_attribute, "MHD_POWELL");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "MHD_POWELL_LIMIT_TIMESTEP", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "MHD_POWELL_LIMIT_TIMESTEP");
my_H5Aclose(hdf5_attribute, "MHD_POWELL_LIMIT_TIMESTEP");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "MHD_SEEDFIELD", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "MHD_SEEDFIELD");
my_H5Aclose(hdf5_attribute, "MHD_SEEDFIELD");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "RIEMANN_HLLD", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "RIEMANN_HLLD");
my_H5Aclose(hdf5_attribute, "RIEMANN_HLLD");
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
hdf5_attribute = my_H5Acreate(handle, "REFINEMENT_HIGH_RES_GAS", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "REFINEMENT_HIGH_RES_GAS");
my_H5Aclose(hdf5_attribute, "REFINEMENT_HIGH_RES_GAS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SELFGRAVITY", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SELFGRAVITY");
my_H5Aclose(hdf5_attribute, "SELFGRAVITY");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "HIERARCHICAL_GRAVITY", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "HIERARCHICAL_GRAVITY");
my_H5Aclose(hdf5_attribute, "HIERARCHICAL_GRAVITY");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "CELL_CENTER_GRAVITY", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "CELL_CENTER_GRAVITY");
my_H5Aclose(hdf5_attribute, "CELL_CENTER_GRAVITY");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "ALLOW_DIRECT_SUMMATION", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "ALLOW_DIRECT_SUMMATION");
my_H5Aclose(hdf5_attribute, "ALLOW_DIRECT_SUMMATION");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "DIRECT_SUMMATION_THRESHOLD", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 16000;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "DIRECT_SUMMATION_THRESHOLD");
my_H5Aclose(hdf5_attribute, "DIRECT_SUMMATION_THRESHOLD");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "ENFORCE_JEANS_STABILITY_OF_CELLS_EEOS", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "ENFORCE_JEANS_STABILITY_OF_CELLS_EEOS");
my_H5Aclose(hdf5_attribute, "ENFORCE_JEANS_STABILITY_OF_CELLS_EEOS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "ENFORCE_JEANS_STABILITY_OF_CELLS", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "ENFORCE_JEANS_STABILITY_OF_CELLS");
my_H5Aclose(hdf5_attribute, "ENFORCE_JEANS_STABILITY_OF_CELLS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "EVALPOTENTIAL", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "EVALPOTENTIAL");
my_H5Aclose(hdf5_attribute, "EVALPOTENTIAL");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "NSOFTTYPES", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 4;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "NSOFTTYPES");
my_H5Aclose(hdf5_attribute, "NSOFTTYPES");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "MULTIPLE_NODE_SOFTENING", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "MULTIPLE_NODE_SOFTENING");
my_H5Aclose(hdf5_attribute, "MULTIPLE_NODE_SOFTENING");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "INDIVIDUAL_GRAVITY_SOFTENING", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 32;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "INDIVIDUAL_GRAVITY_SOFTENING");
my_H5Aclose(hdf5_attribute, "INDIVIDUAL_GRAVITY_SOFTENING");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "ADAPTIVE_HYDRO_SOFTENING", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "ADAPTIVE_HYDRO_SOFTENING");
my_H5Aclose(hdf5_attribute, "ADAPTIVE_HYDRO_SOFTENING");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "PMGRID", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 256;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "PMGRID");
my_H5Aclose(hdf5_attribute, "PMGRID");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "RCUT", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 5.0;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "RCUT");
my_H5Aclose(hdf5_attribute, "RCUT");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "PLACEHIGHRESREGION", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 2;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "PLACEHIGHRESREGION");
my_H5Aclose(hdf5_attribute, "PLACEHIGHRESREGION");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "ENLARGEREGION", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 1.1;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "ENLARGEREGION");
my_H5Aclose(hdf5_attribute, "ENLARGEREGION");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GRIDBOOST", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 1;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "GRIDBOOST");
my_H5Aclose(hdf5_attribute, "GRIDBOOST");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "PM_ZOOM_OPTIMIZED", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "PM_ZOOM_OPTIMIZED");
my_H5Aclose(hdf5_attribute, "PM_ZOOM_OPTIMIZED");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "CHUNKING", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "CHUNKING");
my_H5Aclose(hdf5_attribute, "CHUNKING");
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
hdf5_attribute = my_H5Acreate(handle, "OUTPUT_COORDINATES_IN_DOUBLEPRECISION", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "OUTPUT_COORDINATES_IN_DOUBLEPRECISION");
my_H5Aclose(hdf5_attribute, "OUTPUT_COORDINATES_IN_DOUBLEPRECISION");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "NGB_TREE_DOUBLEPRECISION", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "NGB_TREE_DOUBLEPRECISION");
my_H5Aclose(hdf5_attribute, "NGB_TREE_DOUBLEPRECISION");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "FOF", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "FOF");
my_H5Aclose(hdf5_attribute, "FOF");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "FOF_PRIMARY_LINK_TYPES", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 2;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "FOF_PRIMARY_LINK_TYPES");
my_H5Aclose(hdf5_attribute, "FOF_PRIMARY_LINK_TYPES");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "FOF_SECONDARY_LINK_TYPES", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 1+4+16+32;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "FOF_SECONDARY_LINK_TYPES");
my_H5Aclose(hdf5_attribute, "FOF_SECONDARY_LINK_TYPES");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SUBFIND", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SUBFIND");
my_H5Aclose(hdf5_attribute, "SUBFIND");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SAVE_HSML_IN_SNAPSHOT", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SAVE_HSML_IN_SNAPSHOT");
my_H5Aclose(hdf5_attribute, "SAVE_HSML_IN_SNAPSHOT");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SUBFIND_CALC_MORE", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SUBFIND_CALC_MORE");
my_H5Aclose(hdf5_attribute, "SUBFIND_CALC_MORE");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SOFTEREQS", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SOFTEREQS");
my_H5Aclose(hdf5_attribute, "SOFTEREQS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "CREATE_SUBFOFS", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "CREATE_SUBFOFS");
my_H5Aclose(hdf5_attribute, "CREATE_SUBFOFS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GAS_BASED_SEED_MODEL", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "GAS_BASED_SEED_MODEL");
my_H5Aclose(hdf5_attribute, "GAS_BASED_SEED_MODEL");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "ONLY_SEED_IN_VALID_FOFS", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "ONLY_SEED_IN_VALID_FOFS");
my_H5Aclose(hdf5_attribute, "ONLY_SEED_IN_VALID_FOFS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SEED_HALO_ENVIRONMENT_CRITERION", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SEED_HALO_ENVIRONMENT_CRITERION");
my_H5Aclose(hdf5_attribute, "SEED_HALO_ENVIRONMENT_CRITERION");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SEED_BASED_ON_PROBABLISTIC_HALO_PROPERTIES", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SEED_BASED_ON_PROBABLISTIC_HALO_PROPERTIES");
my_H5Aclose(hdf5_attribute, "SEED_BASED_ON_PROBABLISTIC_HALO_PROPERTIES");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "PROBABILISTIC_SEED_MASS_HALO_MASS_RATIO_CRITERION", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "PROBABILISTIC_SEED_MASS_HALO_MASS_RATIO_CRITERION");
my_H5Aclose(hdf5_attribute, "PROBABILISTIC_SEED_MASS_HALO_MASS_RATIO_CRITERION");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "EVOLVING_SEEDHALOMASS_DISTRIBUTION_DOUBLE_POWERLAW_MODEL2", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "EVOLVING_SEEDHALOMASS_DISTRIBUTION_DOUBLE_POWERLAW_MODEL2");
my_H5Aclose(hdf5_attribute, "EVOLVING_SEEDHALOMASS_DISTRIBUTION_DOUBLE_POWERLAW_MODEL2");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "INCLUDE_MERGERS_OF_UNRESOLVED_SEED_BHS", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "INCLUDE_MERGERS_OF_UNRESOLVED_SEED_BHS");
my_H5Aclose(hdf5_attribute, "INCLUDE_MERGERS_OF_UNRESOLVED_SEED_BHS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "PTYPE_USED_FOR_ENVIRONMENT_BASED_SEEDING", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 1;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "PTYPE_USED_FOR_ENVIRONMENT_BASED_SEEDING");
my_H5Aclose(hdf5_attribute, "PTYPE_USED_FOR_ENVIRONMENT_BASED_SEEDING");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "STORE_MERGERS_IN_SNAPSHOT", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "STORE_MERGERS_IN_SNAPSHOT");
my_H5Aclose(hdf5_attribute, "STORE_MERGERS_IN_SNAPSHOT");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "OUTPUT_LOG_FILES_FOR_SEEDING", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "OUTPUT_LOG_FILES_FOR_SEEDING");
my_H5Aclose(hdf5_attribute, "OUTPUT_LOG_FILES_FOR_SEEDING");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "OUTPUT_HOST_PROPERTIES_FOR_BH_MERGERS", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "OUTPUT_HOST_PROPERTIES_FOR_BH_MERGERS");
my_H5Aclose(hdf5_attribute, "OUTPUT_HOST_PROPERTIES_FOR_BH_MERGERS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "PREVENT_SPURIOUS_RESEEDING2", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "PREVENT_SPURIOUS_RESEEDING2");
my_H5Aclose(hdf5_attribute, "PREVENT_SPURIOUS_RESEEDING2");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "PREVENT_SPURIOUS_RESEEDING", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "PREVENT_SPURIOUS_RESEEDING");
my_H5Aclose(hdf5_attribute, "PREVENT_SPURIOUS_RESEEDING");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "OUTPUT_STELLAR_AGE", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "OUTPUT_STELLAR_AGE");
my_H5Aclose(hdf5_attribute, "OUTPUT_STELLAR_AGE");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "BLACK_HOLES", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "BLACK_HOLES");
my_H5Aclose(hdf5_attribute, "BLACK_HOLES");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "BH_THERMALFEEDBACK", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "BH_THERMALFEEDBACK");
my_H5Aclose(hdf5_attribute, "BH_THERMALFEEDBACK");
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
hdf5_attribute = my_H5Acreate(handle, "BH_USE_ALFVEN_SPEED_IN_BONDI", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "BH_USE_ALFVEN_SPEED_IN_BONDI");
my_H5Aclose(hdf5_attribute, "BH_USE_ALFVEN_SPEED_IN_BONDI");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "BH_NEW_CENTERING", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "BH_NEW_CENTERING");
my_H5Aclose(hdf5_attribute, "BH_NEW_CENTERING");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "BH_PRESSURE_CRITERION", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "BH_PRESSURE_CRITERION");
my_H5Aclose(hdf5_attribute, "BH_PRESSURE_CRITERION");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "BH_ADIOS_WIND", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "BH_ADIOS_WIND");
my_H5Aclose(hdf5_attribute, "BH_ADIOS_WIND");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "BH_ADIOS_WIND_WITH_QUASARTHRESHOLD", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "BH_ADIOS_WIND_WITH_QUASARTHRESHOLD");
my_H5Aclose(hdf5_attribute, "BH_ADIOS_WIND_WITH_QUASARTHRESHOLD");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "BH_ADIOS_WIND_WITH_VARIABLE_QUASARTHRESHOLD", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "BH_ADIOS_WIND_WITH_VARIABLE_QUASARTHRESHOLD");
my_H5Aclose(hdf5_attribute, "BH_ADIOS_WIND_WITH_VARIABLE_QUASARTHRESHOLD");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "BH_ADIOS_RANDOMIZED", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "BH_ADIOS_RANDOMIZED");
my_H5Aclose(hdf5_attribute, "BH_ADIOS_RANDOMIZED");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "BH_ADIOS_ONLY_ABOVE_MINIMUM_DENSITY", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "BH_ADIOS_ONLY_ABOVE_MINIMUM_DENSITY");
my_H5Aclose(hdf5_attribute, "BH_ADIOS_ONLY_ABOVE_MINIMUM_DENSITY");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "PROCESS_TIMES_OF_OUTPUTLIST", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "PROCESS_TIMES_OF_OUTPUTLIST");
my_H5Aclose(hdf5_attribute, "PROCESS_TIMES_OF_OUTPUTLIST");
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
hdf5_attribute = my_H5Acreate(handle, "OVERRIDE_PEANOGRID_WARNING", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "OVERRIDE_PEANOGRID_WARNING");
my_H5Aclose(hdf5_attribute, "OVERRIDE_PEANOGRID_WARNING");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "OUTPUT_CPU_CSV", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "OUTPUT_CPU_CSV");
my_H5Aclose(hdf5_attribute, "OUTPUT_CPU_CSV");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "OUTPUT_CENTER_OF_MASS", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "OUTPUT_CENTER_OF_MASS");
my_H5Aclose(hdf5_attribute, "OUTPUT_CENTER_OF_MASS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "OUTPUT_PRESSURE", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "OUTPUT_PRESSURE");
my_H5Aclose(hdf5_attribute, "OUTPUT_PRESSURE");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "OUTPUTPOTENTIAL", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "OUTPUTPOTENTIAL");
my_H5Aclose(hdf5_attribute, "OUTPUTPOTENTIAL");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "HAVE_HDF5", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "HAVE_HDF5");
my_H5Aclose(hdf5_attribute, "HAVE_HDF5");
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
hdf5_attribute = my_H5Acreate(handle, "GENERATE_GAS_IN_ICS", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "GENERATE_GAS_IN_ICS");
my_H5Aclose(hdf5_attribute, "GENERATE_GAS_IN_ICS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SPLIT_PARTICLE_TYPE", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 2+4+8;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "SPLIT_PARTICLE_TYPE");
my_H5Aclose(hdf5_attribute, "SPLIT_PARTICLE_TYPE");
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
hdf5_attribute = my_H5Acreate(handle, "GFM_CONST_IMF", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 0;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "GFM_CONST_IMF");
my_H5Aclose(hdf5_attribute, "GFM_CONST_IMF");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM_PREENRICH", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "GFM_PREENRICH");
my_H5Aclose(hdf5_attribute, "GFM_PREENRICH");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM_WINDS", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "GFM_WINDS");
my_H5Aclose(hdf5_attribute, "GFM_WINDS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM_WINDS_VARIABLE", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 1;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "GFM_WINDS_VARIABLE");
my_H5Aclose(hdf5_attribute, "GFM_WINDS_VARIABLE");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM_WINDS_VARIABLE_HUBBLE", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "GFM_WINDS_VARIABLE_HUBBLE");
my_H5Aclose(hdf5_attribute, "GFM_WINDS_VARIABLE_HUBBLE");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM_WIND_ENERGY_METAL_DEPENDENCE", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "GFM_WIND_ENERGY_METAL_DEPENDENCE");
my_H5Aclose(hdf5_attribute, "GFM_WIND_ENERGY_METAL_DEPENDENCE");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM_WINDS_STRIPPING", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "GFM_WINDS_STRIPPING");
my_H5Aclose(hdf5_attribute, "GFM_WINDS_STRIPPING");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM_WINDS_THERMAL_NEWDEF", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "GFM_WINDS_THERMAL_NEWDEF");
my_H5Aclose(hdf5_attribute, "GFM_WINDS_THERMAL_NEWDEF");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM_COOLING_METAL", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "GFM_COOLING_METAL");
my_H5Aclose(hdf5_attribute, "GFM_COOLING_METAL");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM_AGN_RADIATION", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "GFM_AGN_RADIATION");
my_H5Aclose(hdf5_attribute, "GFM_AGN_RADIATION");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM_STELLAR_PHOTOMETRICS", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "GFM_STELLAR_PHOTOMETRICS");
my_H5Aclose(hdf5_attribute, "GFM_STELLAR_PHOTOMETRICS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM_OUTPUT_MASK", H5T_NATIVE_DOUBLE, hdf5_dataspace, H5P_DEFAULT);
val = 1+2+4+8+16+32+64+256;
my_H5Awrite(hdf5_attribute, H5T_NATIVE_DOUBLE, &val, "GFM_OUTPUT_MASK");
my_H5Aclose(hdf5_attribute, "GFM_OUTPUT_MASK");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM_NORMALIZED_METAL_ADVECTION", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "GFM_NORMALIZED_METAL_ADVECTION");
my_H5Aclose(hdf5_attribute, "GFM_NORMALIZED_METAL_ADVECTION");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM_OUTPUT_BIRTH_POS", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "GFM_OUTPUT_BIRTH_POS");
my_H5Aclose(hdf5_attribute, "GFM_OUTPUT_BIRTH_POS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM_CHEMTAGS", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "GFM_CHEMTAGS");
my_H5Aclose(hdf5_attribute, "GFM_CHEMTAGS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM_DISCRETE_ENRICHMENT", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "GFM_DISCRETE_ENRICHMENT");
my_H5Aclose(hdf5_attribute, "GFM_DISCRETE_ENRICHMENT");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM_SPLITFE", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "GFM_SPLITFE");
my_H5Aclose(hdf5_attribute, "GFM_SPLITFE");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "GFM_RPROCESS", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "GFM_RPROCESS");
my_H5Aclose(hdf5_attribute, "GFM_RPROCESS");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

hdf5_dataspace = my_H5Screate(H5S_SCALAR);
hdf5_attribute = my_H5Acreate(handle, "SHOCK_FINDER_BEFORE_OUTPUT", atype, hdf5_dataspace, H5P_DEFAULT);
my_H5Awrite(hdf5_attribute, atype, "", "SHOCK_FINDER_BEFORE_OUTPUT");
my_H5Aclose(hdf5_attribute, "SHOCK_FINDER_BEFORE_OUTPUT");
my_H5Sclose(hdf5_dataspace, H5S_SCALAR);

my_H5Tclose(atype);
}
#endif

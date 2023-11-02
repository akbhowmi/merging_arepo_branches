#include <stdio.h>
void output_compile_time_options(void)
{
printf(
"        NTYPES=6\n"
"        PERIODIC\n"
"        COOLING\n"
"        UVB_SELF_SHIELDING\n"
"        USE_SFR\n"
"        VORONOI\n"
"        MHD\n"
"        MHD_POWELL\n"
"        MHD_POWELL_LIMIT_TIMESTEP\n"
"        MHD_SEEDFIELD\n"
"        RIEMANN_HLLD\n"
"        REGULARIZE_MESH_CM_DRIFT\n"
"        REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED\n"
"        REGULARIZE_MESH_FACE_ANGLE\n"
"        TREE_BASED_TIMESTEPS\n"
"        REFINEMENT_SPLIT_CELLS\n"
"        REFINEMENT_MERGE_CELLS\n"
"        REFINEMENT_HIGH_RES_GAS\n"
"        SELFGRAVITY\n"
"        HIERARCHICAL_GRAVITY\n"
"        CELL_CENTER_GRAVITY\n"
"        ALLOW_DIRECT_SUMMATION\n"
"        DIRECT_SUMMATION_THRESHOLD=16000\n"
"        ENFORCE_JEANS_STABILITY_OF_CELLS_EEOS\n"
"        ENFORCE_JEANS_STABILITY_OF_CELLS\n"
"        EVALPOTENTIAL\n"
"        NSOFTTYPES=4\n"
"        MULTIPLE_NODE_SOFTENING\n"
"        INDIVIDUAL_GRAVITY_SOFTENING=32\n"
"        ADAPTIVE_HYDRO_SOFTENING\n"
"        PMGRID=256\n"
"        RCUT=5.0\n"
"        PLACEHIGHRESREGION=2\n"
"        ENLARGEREGION=1.1\n"
"        GRIDBOOST=1\n"
"        PM_ZOOM_OPTIMIZED\n"
"        CHUNKING\n"
"        DOUBLEPRECISION=1\n"
"        DOUBLEPRECISION_FFTW\n"
"        OUTPUT_COORDINATES_IN_DOUBLEPRECISION\n"
"        NGB_TREE_DOUBLEPRECISION\n"
"        FOF\n"
"        FOF_PRIMARY_LINK_TYPES=2\n"
"        FOF_SECONDARY_LINK_TYPES=1+4+16+32\n"
"        SUBFIND\n"
"        SAVE_HSML_IN_SNAPSHOT\n"
"        SUBFIND_CALC_MORE\n"
"        SOFTEREQS\n"
"        CREATE_SUBFOFS\n"
"        GAS_BASED_SEED_MODEL\n"
"        ONLY_SEED_IN_VALID_FOFS\n"
"        SEED_HALO_ENVIRONMENT_CRITERION\n"
"        SEED_BASED_ON_PROBABLISTIC_HALO_PROPERTIES\n"
"        PROBABILISTIC_SEED_MASS_HALO_MASS_RATIO_CRITERION\n"
"        EVOLVING_SEEDHALOMASS_DISTRIBUTION_DOUBLE_POWERLAW_MODEL2\n"
"        INCLUDE_MERGERS_OF_UNRESOLVED_SEED_BHS\n"
"        PTYPE_USED_FOR_ENVIRONMENT_BASED_SEEDING=1\n"
"        STORE_MERGERS_IN_SNAPSHOT\n"
"        OUTPUT_LOG_FILES_FOR_SEEDING\n"
"        OUTPUT_HOST_PROPERTIES_FOR_BH_MERGERS\n"
"        PREVENT_SPURIOUS_RESEEDING2\n"
"        PREVENT_SPURIOUS_RESEEDING\n"
"        OUTPUT_STELLAR_AGE\n"
"        BLACK_HOLES\n"
"        BH_THERMALFEEDBACK\n"
"        DRAINGAS=3\n"
"        BH_EXACT_INTEGRATION\n"
"        BH_BONDI_DEFAULT\n"
"        BH_DO_NOT_PREVENT_MERGERS\n"
"        BH_USE_ALFVEN_SPEED_IN_BONDI\n"
"        BH_NEW_CENTERING\n"
"        BH_PRESSURE_CRITERION\n"
"        BH_ADIOS_WIND\n"
"        BH_ADIOS_WIND_WITH_QUASARTHRESHOLD\n"
"        BH_ADIOS_WIND_WITH_VARIABLE_QUASARTHRESHOLD\n"
"        BH_ADIOS_RANDOMIZED\n"
"        BH_ADIOS_ONLY_ABOVE_MINIMUM_DENSITY\n"
"        PROCESS_TIMES_OF_OUTPUTLIST\n"
"        VORONOI_DYNAMIC_UPDATE\n"
"        NO_MPI_IN_PLACE\n"
"        NO_ISEND_IRECV_IN_DOMAIN\n"
"        FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG\n"
"        ENLARGE_DYNAMIC_RANGE_IN_TIME\n"
"        OVERRIDE_PEANOGRID_WARNING\n"
"        OUTPUT_CPU_CSV\n"
"        OUTPUT_CENTER_OF_MASS\n"
"        OUTPUT_PRESSURE\n"
"        OUTPUTPOTENTIAL\n"
"        HAVE_HDF5\n"
"        DEBUG\n"
"        HOST_MEMORY_REPORTING\n"
"        GENERATE_GAS_IN_ICS\n"
"        SPLIT_PARTICLE_TYPE=2+4+8\n"
"        GFM\n"
"        GFM_STELLAR_EVOLUTION=0\n"
"        GFM_CONST_IMF=0\n"
"        GFM_PREENRICH\n"
"        GFM_WINDS\n"
"        GFM_WINDS_VARIABLE=1\n"
"        GFM_WINDS_VARIABLE_HUBBLE\n"
"        GFM_WIND_ENERGY_METAL_DEPENDENCE\n"
"        GFM_WINDS_STRIPPING\n"
"        GFM_WINDS_THERMAL_NEWDEF\n"
"        GFM_COOLING_METAL\n"
"        GFM_AGN_RADIATION\n"
"        GFM_STELLAR_PHOTOMETRICS\n"
"        GFM_OUTPUT_MASK=1+2+4+8+16+32+64+256\n"
"        GFM_NORMALIZED_METAL_ADVECTION\n"
"        GFM_OUTPUT_BIRTH_POS\n"
"        GFM_CHEMTAGS\n"
"        GFM_DISCRETE_ENRICHMENT\n"
"        GFM_SPLITFE\n"
"        GFM_RPROCESS\n"
"        SHOCK_FINDER_BEFORE_OUTPUT\n"
"\n");
}
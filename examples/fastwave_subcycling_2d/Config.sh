# examples/fastwave_subcycling_2d
NTYPES=6
NSOFTTYPES=2
TWODIMS
VORONOI
DOUBLEPRECISION=1
GAMMA=1.6666666666666667
REGULARIZE_MESH_CM_DRIFT
REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED
REGULARIZE_MESH_FACE_ANGLE

LONG_Y=0.2

# generic
VORONOI_DYNAMIC_UPDATE
HAVE_HDF5
DEBUG

MHD
RIEMANN_HLLD
#MHD_DIVBCLEANING

MESHRELAX_DENSITY_IN_INPUT

OUTPUT_IN_DOUBLEPRECISION                # snapshot files will be written in double precision
INPUT_IN_DOUBLEPRECISION                 # initial conditions are in double precision

TREE_BASED_TIMESTEPS     # non-local timestep criterion (take 'signal speed' into account)

# Below this can be removed if we do not use braginskii
BRAGINSKII_VISCOSITY
BRAGINSKII_VISCOSITY_SUBCYCLE

TETRA_INDEX_IN_FACE
VORONOI_MESH_KEEP_DT_AND_DTC
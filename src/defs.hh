#include "Ippl.h"

#define DIM 2
#define GUARDCELLSIZE 1

#define indent_l0 ""
#define indent_l1 "  "
#define indent_l2 "    "
#define indent_l3 "      "
#define indent_l4 "        "
#define indent_l5 "          "

typedef unsigned int                                Idx_t;
typedef double                                      Scalar_t;
typedef double                                      SolverPrecision_t;
typedef std::pair<int,int>                          coord_t;

typedef ParticleSpatialLayout<double, DIM>          playout_t;
typedef playout_t::SingleParticlePos_t              Vector_t;
typedef Vektor<double, 3>                           Vector3D_t;

typedef IntCIC                                      IntOp_t;
typedef UniformCartesian<DIM,double>                Mesh_t;
typedef Edge                                        Center_t;
typedef CenteredFieldLayout<DIM, Mesh_t, Center_t>  FieldLayout_t;
typedef Field<Vector_t, DIM, Mesh_t, Center_t>      VField_t;
typedef Field<double, DIM, Mesh_t, Center_t>        SField_t;
typedef GuardCellSizes<DIM>                         GCS_t;

#pragma once

#include <memory>
#include <ceed.h>

// -----------------------------------------------------------------------------
// libCEED Data Structs
// -----------------------------------------------------------------------------

// libCEED data struct for level
typedef struct CeedData_ *CeedData;
struct CeedData_
{
  Ceed ceed;
  CeedBasis basis_x, basis_u;
  CeedElemRestriction elem_restr_x, elem_restr_u, elem_restr_qd, elem_restr_u_i, elem_restr_qd_i;
  CeedQFunction qf_build;
  CeedQFunctionContext build_ctx;
  CeedQFunction qf_apply;
  CeedOperator op_build;
  CeedOperator op_apply, op_restrict;
  CeedVector q_data, x_ceed, y_ceed;
  CeedInt q_data_size;
  CeedVector u, v;
  CeedVector mesh_coords;

  // Constructor for initialization
  CeedData_();
};

typedef struct FEproblemData_ *FEproblemData;
struct FEproblemData_
{
  CeedInt dim;
  CeedInt num_comp;
  CeedInt num_poly;
  CeedInt num_qpts;
  CeedQuadMode quad_mode;
  CeedInt num_xyz[3];
  CeedInt num_dofs;
  CeedInt mesh_size;
  CeedInt sol_size;

  // Constructor for initialization
  FEproblemData_();
};

/// A structure used to pass additional data to f_build_diff
struct BuildContext
{
  CeedInt dim, space_dim;

  BuildContext(CeedInt dim = 0, CeedInt space_dim = 0);
};

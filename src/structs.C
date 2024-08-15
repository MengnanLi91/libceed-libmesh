#include "structs.h"

CeedData_::CeedData_()
    : ceed(nullptr),
      basis_x(nullptr), basis_u(nullptr),
      elem_restr_x(nullptr), elem_restr_u(nullptr),
      elem_restr_qd(nullptr), elem_restr_u_i(nullptr),
      elem_restr_qd_i(nullptr),
      qf_build(nullptr),
      build_ctx(nullptr),
      qf_apply(nullptr),
      op_build(nullptr),
      op_apply(nullptr),
      op_restrict(nullptr),
      q_data(nullptr),
      x_ceed(nullptr),
      y_ceed(nullptr),
      q_data_size(0),
      u(nullptr),
      v(nullptr),
      mesh_coords(nullptr)
{
}

FEproblemData_::FEproblemData_()
    : dim(0),
      num_comp(0),
      num_poly(0),
      num_qpts(0),
      quad_mode(CEED_GAUSS), // Assuming default to Gaussian quadrature
      num_xyz{0, 0, 0},
      num_dofs(0),
      mesh_size(0),
      sol_size(0)
{
}

BuildContext::BuildContext(CeedInt dim, CeedInt space_dim)
    : dim(dim), space_dim(space_dim)
{
}

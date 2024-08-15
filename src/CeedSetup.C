#include "CeedSetup.h"
#include <iostream>
#include <algorithm>
#include <math.h>

#include "../qfunctions/Diffusion.h"

CeedSetup::CeedSetup(const char *resource)
{
  // Initialize libCEED
  CeedInit(resource, &_ceed);
}

CeedSetup::~CeedSetup()
{
  // Destroy libCEED objects
  // Free dynamically allocated memory.
  // CeedVectorDestroy(&ceed_data->u);
  // CeedVectorDestroy(&ceed_data->v);
  // CeedVectorDestroy(&ceed_data->q_data);
  // CeedVectorDestroy(&ceed_data->mesh_coords);
  // CeedOperatorDestroy(&ceed_data->op_apply);
  // CeedQFunctionDestroy(&ceed_data->qf_apply);
  // CeedQFunctionContextDestroy(&ceed_data->build_ctx);
  // CeedOperatorDestroy(&ceed_data->op_build);
  // CeedQFunctionDestroy(&ceed_data->qf_build);
  // CeedElemRestrictionDestroy(&ceed_data->elem_restr_u);
  // CeedElemRestrictionDestroy(&ceed_data->elem_restr_x);
  // CeedElemRestrictionDestroy(&ceed_data->elem_restr_qd);
  // CeedBasisDestroy(&ceed_data->basis_u);
  // CeedBasisDestroy(&ceed_data->basis_x);
  // CeedDestroy(&_ceed);
}

void CeedSetup::createCeedBasis(FEproblemData feproblem_data)
{
  CeedInt dim = feproblem_data->dim;
  CeedInt num_comp = feproblem_data->num_comp;
  CeedInt poly_order = feproblem_data->num_poly;
  CeedInt num_qpts = feproblem_data->num_qpts;
  CeedQuadMode quad_mode = feproblem_data->quad_mode;

  CeedBasisCreateTensorH1Lagrange(
      _ceed, dim, num_comp, poly_order, num_qpts, quad_mode, &_mesh_basis);
  CeedBasisCreateTensorH1Lagrange(
      _ceed, dim, 1, poly_order, num_qpts, quad_mode, &_sol_basis);
}

void CeedSetup::getCartesianMeshSize(FEproblemData feproblem_data)
{
  // Use the approximate formula:
  //    prob_size ~ num_elem * degree^dim
  CeedInt prob_size = feproblem_data->num_dofs;
  CeedInt degree = feproblem_data->num_poly - 1;
  CeedInt dim = feproblem_data->dim;
  CeedInt num_elem = prob_size / CeedIntPow(degree, dim);
  CeedInt s = 0; // find s: num_elem/2 < 2^s <= num_elem
  while (num_elem > 1)
  {
    num_elem /= 2;
    s++;
  }
  CeedInt r = s % dim;
  for (CeedInt d = 0; d < dim; d++)
  {
    CeedInt sd = s / dim;
    if (r > 0)
    {
      sd++;
      r--;
    }
    feproblem_data->num_xyz[d] = 1 << sd;
  }
}

void CeedSetup::buildCartesianRestriction(FEproblemData feproblem_data, CeedInt num_comp, CeedInt *size,
                                          CeedElemRestriction *restriction, CeedElemRestriction *q_data_restriction)
{
  CeedInt dim = feproblem_data->dim;
  CeedInt p = feproblem_data->num_poly;
  CeedInt num_qpts = feproblem_data->num_qpts;
  CeedInt num_nodes = CeedIntPow(p, dim);        // number of scalar nodes per element
  CeedInt elem_qpts = CeedIntPow(num_qpts, dim); // number of qpts per element
  CeedInt nd[3], num_elem = 1, scalar_size = 1;
  for (CeedInt d = 0; d < dim; d++)
  {
    num_elem *= feproblem_data->num_xyz[d];
    nd[d] = feproblem_data->num_xyz[d] * (p - 1) + 1;
    scalar_size *= nd[d];
  }
  *size = scalar_size * num_comp;
  // elem:         0             1                 n-1
  //           |---*-...-*---|---*-...-*---|- ... -|--...--|
  // num_nodes:   0   1    p-1  p  p+1       2*p             n*p
  CeedInt *el_nodes = static_cast<CeedInt *>(malloc(sizeof(CeedInt) * num_elem * num_nodes));
  for (CeedInt e = 0; e < num_elem; e++)
  {
    CeedInt e_xyz[3] = {1, 1, 1}, re = e;
    for (CeedInt d = 0; d < dim; d++)
    {
      e_xyz[d] = re % feproblem_data->num_xyz[d];
      re /= feproblem_data->num_xyz[d];
    }
    CeedInt *local_elem_nodes = el_nodes + e * num_nodes;
    for (CeedInt l_nodes = 0; l_nodes < num_nodes; l_nodes++)
    {
      CeedInt g_nodes = 0, g_nodes_stride = 1, r_nodes = l_nodes;
      for (CeedInt d = 0; d < dim; d++)
      {
        g_nodes += (e_xyz[d] * (p - 1) + r_nodes % p) * g_nodes_stride;
        g_nodes_stride *= nd[d];
        r_nodes /= p;
      }
      local_elem_nodes[l_nodes] = g_nodes;
    }
  }

  if (restriction)
    CeedElemRestrictionCreate(_ceed, num_elem, num_nodes, num_comp, scalar_size, num_comp * scalar_size, CEED_MEM_HOST, CEED_COPY_VALUES, el_nodes,
                              restriction);
  free(el_nodes);

  if (q_data_restriction)
  {
    CeedElemRestrictionCreateStrided(_ceed, num_elem, elem_qpts, num_comp, num_comp * elem_qpts * num_elem, CEED_STRIDES_BACKEND, q_data_restriction);
  }
}

void CeedSetup::SetCartesianMeshCoords(FEproblemData feproblem_data, CeedData ceed_data)
{
  CeedInt dim = feproblem_data->dim;
  CeedInt mesh_size = feproblem_data->mesh_size;

  CeedVectorCreate(_ceed, mesh_size, &ceed_data->mesh_coords);

  CeedInt p = feproblem_data->num_poly;
  CeedInt nd[3], scalar_size = 1;
  for (CeedInt d = 0; d < dim; d++)
  {
    nd[d] = feproblem_data->num_xyz[d] * (p - 1) + 1;
    scalar_size *= nd[d];
  }
  CeedScalar *coords;
  CeedVectorGetArrayWrite(ceed_data->mesh_coords, CEED_MEM_HOST, &coords);
  CeedScalar *nodes = static_cast<CeedScalar *>(malloc(sizeof(CeedScalar) * p));
  // The H1 basis uses Lobatto quadrature points as nodes.
  CeedLobattoQuadrature(p, nodes, NULL); // nodes are in [-1,1]
  for (CeedInt i = 0; i < p; i++)
    nodes[i] = 0.5 + 0.5 * nodes[i];
  for (CeedInt gs_nodes = 0; gs_nodes < scalar_size; gs_nodes++)
  {
    CeedInt r_nodes = gs_nodes;
    for (CeedInt d = 0; d < dim; d++)
    {
      CeedInt d1d = r_nodes % nd[d];
      coords[gs_nodes + scalar_size * d] = ((d1d / (p - 1)) + nodes[d1d % (p - 1)]) / feproblem_data->num_xyz[d];
      r_nodes /= nd[d];
    }
  }
  free(nodes);
  CeedVectorRestoreArray(ceed_data->mesh_coords, &coords);
}

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

CeedScalar CeedSetup::TransformMeshCoords(FEproblemData feproblem_data, CeedData ceed_data)
{
  CeedInt dim = feproblem_data->dim;
  CeedInt mesh_size = feproblem_data->mesh_size;
  CeedVector mesh_coords = ceed_data->mesh_coords;
  CeedScalar exact_surface_area = (dim == 1 ? 2 : dim == 2 ? 4
                                                           : 6);
  CeedScalar *coords;

  CeedVectorGetArray(mesh_coords, CEED_MEM_HOST, &coords);
  for (CeedInt i = 0; i < mesh_size; i++)
  {
    // map [0,1] to [0,1] varying the mesh density
    coords[i] = 0.5 + 1. / sqrt(3.) * sin((2. / 3.) * M_PI * (coords[i] - 0.5));
  }
  CeedVectorRestoreArray(mesh_coords, &coords);

  return exact_surface_area;
}

void CeedSetup::setupQfunction(FEproblemData feproblem_data, CeedData ceed_data)
{
  CeedInt num_comp_x = feproblem_data->num_comp;
  CeedInt dim = feproblem_data->dim;
  // Context data to be passed to the 'build_diff' QFunction.
  struct BuildContext build_ctx_data;
  build_ctx_data.dim = build_ctx_data.space_dim = feproblem_data->dim;
  CeedQFunctionContextCreate(_ceed, &ceed_data->build_ctx);
  CeedQFunctionContextSetData(ceed_data->build_ctx, CEED_MEM_HOST, CEED_USE_POINTER, sizeof(build_ctx_data), &build_ctx_data);

  // Create the QFunction that builds the diffusion operator (i.e. computes its quadrature data) and set its context data.

  // This creates the QFunction directly.
  CeedQFunctionCreateInterior(_ceed, 1, build_diff, build_diff_loc, &ceed_data->qf_build);
  CeedQFunctionAddInput(ceed_data->qf_build, "dx", num_comp_x * dim, CEED_EVAL_GRAD);
  CeedQFunctionAddInput(ceed_data->qf_build, "weights", 1, CEED_EVAL_WEIGHT);
  CeedQFunctionAddOutput(ceed_data->qf_build, "qdata", dim * (dim + 1) / 2, CEED_EVAL_NONE);
  CeedQFunctionSetContext(ceed_data->qf_build, ceed_data->build_ctx);
}

void CeedSetup::setupOperator(FEproblemData feproblem_data, CeedData ceed_data)
{
  // Create the operator that builds the quadrature data for the diffusion operator.

  CeedOperatorCreate(_ceed, ceed_data->qf_build, CEED_QFUNCTION_NONE, CEED_QFUNCTION_NONE, &ceed_data->op_build);
  CeedOperatorSetField(ceed_data->op_build, "dx", ceed_data->elem_restr_x, ceed_data->basis_x, CEED_VECTOR_ACTIVE);
  CeedOperatorSetField(ceed_data->op_build, "weights", CEED_ELEMRESTRICTION_NONE, ceed_data->basis_x, CEED_VECTOR_NONE);
  CeedOperatorSetField(ceed_data->op_build, "qdata", ceed_data->elem_restr_qd, CEED_BASIS_NONE, CEED_VECTOR_ACTIVE);

  // Compute the quadrature data for the diffusion operator.
  CeedInt dim = feproblem_data->dim;
  CeedInt elem_qpts = CeedIntPow(feproblem_data->num_qpts, dim);
  CeedInt num_elem = 1;
  for (CeedInt d = 0; d < dim; d++)
    num_elem *= feproblem_data->num_xyz[d];
  CeedVectorCreate(_ceed, num_elem * elem_qpts * dim * (dim + 1) / 2, &ceed_data->q_data);
  CeedOperatorApply(ceed_data->op_build, ceed_data->mesh_coords, ceed_data->q_data, CEED_REQUEST_IMMEDIATE);

  // Create the QFunction that defines the action of the diffusion operator.
  CeedQFunctionCreateInterior(_ceed, 1, apply_diff, apply_diff_loc, &ceed_data->qf_apply);
  CeedQFunctionAddInput(ceed_data->qf_apply, "du", dim, CEED_EVAL_GRAD);
  CeedQFunctionAddInput(ceed_data->qf_apply, "qdata", dim * (dim + 1) / 2, CEED_EVAL_NONE);
  CeedQFunctionAddOutput(ceed_data->qf_apply, "dv", dim, CEED_EVAL_GRAD);
  CeedQFunctionSetContext(ceed_data->qf_apply, ceed_data->build_ctx);

  // Create the diffusion operator.
  CeedOperatorCreate(_ceed, ceed_data->qf_apply, CEED_QFUNCTION_NONE, CEED_QFUNCTION_NONE, &ceed_data->op_apply);
  CeedOperatorSetField(ceed_data->op_apply, "du", ceed_data->elem_restr_u, ceed_data->basis_u, CEED_VECTOR_ACTIVE);
  CeedOperatorSetField(ceed_data->op_apply, "qdata", ceed_data->elem_restr_qd, CEED_BASIS_NONE, ceed_data->q_data);
  CeedOperatorSetField(ceed_data->op_apply, "dv", ceed_data->elem_restr_u, ceed_data->basis_u, CEED_VECTOR_ACTIVE);
}

CeedScalar CeedSetup::solve(FEproblemData feproblem_data, CeedData ceed_data)
{
  // Create auxiliary solution-size vectors.

  CeedInt dim = feproblem_data->dim;
  CeedInt sol_size = feproblem_data->sol_size;

  CeedVectorCreate(_ceed, feproblem_data->sol_size, &ceed_data->u);
  CeedVectorCreate(_ceed, feproblem_data->sol_size, &ceed_data->v);
  // Initialize 'u' with sum of coordinates, x+y+z.
  {
    CeedScalar *u_array;
    const CeedScalar *x_array;
    CeedVectorGetArrayWrite(ceed_data->u, CEED_MEM_HOST, &u_array);
    CeedVectorGetArrayRead(ceed_data->mesh_coords, CEED_MEM_HOST, &x_array);
    for (CeedInt i = 0; i < feproblem_data->sol_size; i++)
    {
      u_array[i] = 0;
      for (CeedInt d = 0; d < dim; d++)
        u_array[i] += x_array[i + d * sol_size];
    }
    CeedVectorRestoreArray(ceed_data->u, &u_array);
    CeedVectorRestoreArrayRead(ceed_data->mesh_coords, &x_array);
  }

  // Compute the mesh surface area using the diff operator: surface_area = 1^T \cdot abs( K \cdot x).
  CeedOperatorApply(ceed_data->op_apply, ceed_data->u, ceed_data->v, CEED_REQUEST_IMMEDIATE);

  // Compute and print the sum of the entries of 'v' giving the mesh surface area.
  CeedScalar surface_area = 0.;
  {
    const CeedScalar *v_array;
    CeedVectorGetArrayRead(ceed_data->v, CEED_MEM_HOST, &v_array);
    for (CeedInt i = 0; i < sol_size; i++)
      surface_area += fabs(v_array[i]);
    CeedVectorRestoreArrayRead(ceed_data->v, &v_array);
  }
  // LCOV_EXCL_START
  printf(" done.\n");
  printf("Computed mesh surface area : % .14g\n", surface_area);
  // LCOV_EXCL_STOP

  return surface_area;
}

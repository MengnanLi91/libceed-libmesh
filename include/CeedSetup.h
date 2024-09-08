#pragma once

#include <ceed.h>

#include "structs.h"

class CeedSetup
{
public:
  CeedSetup(const char * resource);
  ~CeedSetup();

  // Setup basis function for libceed
  void createCeedBasis(FEproblemData & feproblem_data);

  // Determine the mesh size based on the given approximate problem size.
  void getCartesianMeshSize(FEproblemData & feproblem_data);

  // Build CeedElemRestriction objects describing the mesh and solution discrete representations.
  void buildCartesianRestriction(FEproblemData & feproblem_data,
                                 CeedInt num_comp,
                                 CeedInt * size,
                                 CeedElemRestriction * restriction,
                                 CeedElemRestriction * q_data_restriction);
  // Create a CeedVector with the mesh coordinates.
  void SetCartesianMeshCoords(FEproblemData & feproblem_data);

  // Apply a transformation to the mesh.
  CeedScalar TransformMeshCoords(FEproblemData & feproblem_data);

  CeedScalar ComputeExactSurface(FEproblemData & feproblem_data);

  // Method to setup QFunction
  void setupQfunction(FEproblemData & feproblem_data);
  void setupQfunctionSurface(FEproblemData & feproblem_data);

  // Method to setup Operator
  void setupOperator(FEproblemData & feproblem_data);
  void setupOperatorSurface(FEproblemData & feproblem_data);

  CeedScalar solve(FEproblemData & feproblem_data);
  CeedScalar solveSurface(FEproblemData & feproblem_data);
  // Getters for libCEED objects
  Ceed & getCeed() { return _ceed; }

  CeedElemRestriction elem_restr_x = nullptr;
  CeedElemRestriction elem_restr_u = nullptr;
  CeedElemRestriction elem_restr_qd = nullptr;

private:
  Ceed _ceed;
  CeedBasis _mesh_basis = nullptr;
  CeedBasis _sol_basis = nullptr;

  CeedQFunction qf_build = nullptr;
  CeedQFunctionContext build_ctx = nullptr;
  CeedQFunction qf_apply = nullptr;
  CeedOperator op_build = nullptr;
  CeedOperator op_apply = nullptr;
  CeedVector q_data = nullptr;
  CeedVector u = nullptr;
  CeedVector v = nullptr;
  CeedVector mesh_coords = nullptr;
  struct BuildContext build_ctx_data;
};

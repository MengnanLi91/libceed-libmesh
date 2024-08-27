#pragma once

#include <ceed.h>

#include "structs.h"

class CeedSetup
{
public:
  CeedSetup(const char *resource);
  ~CeedSetup();

  // Setup basis function for libceed
  void createCeedBasis(FEproblemData &feproblem_data);

  // Determine the mesh size based on the given approximate problem size.
  void getCartesianMeshSize(FEproblemData &feproblem_data);

  // Build CeedElemRestriction objects describing the mesh and solution discrete representations.
  void buildCartesianRestriction(FEproblemData &feproblem_data, CeedInt num_comp, CeedInt *size,
                                 CeedElemRestriction *restriction, CeedElemRestriction *q_data_restriction);
  // Create a CeedVector with the mesh coordinates.
  void SetCartesianMeshCoords(FEproblemData &feproblem_data);

  // Apply a transformation to the mesh.
  CeedScalar TransformMeshCoords(FEproblemData &feproblem_data);

  // Method to setup QFunction
  void setupQfunction(FEproblemData &feproblem_data);

  // Method to setup Operator
  void setupOperator(FEproblemData &feproblem_data);

  CeedScalar solve(FEproblemData &feproblem_data);
  // Getters for libCEED objects
  Ceed &getCeed()
  {
    return _ceed;
  }

  CeedElemRestriction elem_restr_x;
  CeedElemRestriction elem_restr_u;
  CeedElemRestriction elem_restr_qd;
  // void setRestriction(CeedElemRestriction &mesh_restriction, CeedElemRestriction &sol_restriction, CeedElemRestriction &q_data_restriction)
  // {
  //   elem_restr_x = mesh_restriction;
  //   elem_restr_u = sol_restriction;
  //   elem_restr_qd = q_data_restriction;
  // };

private:
  Ceed _ceed;
  CeedBasis _mesh_basis;
  CeedBasis _sol_basis;

  CeedQFunction qf_build;
  CeedQFunctionContext build_ctx;
  CeedQFunction qf_apply;
  CeedOperator op_build;
  CeedOperator op_apply;
  CeedVector q_data;
  CeedVector u;
  CeedVector v;
  CeedVector mesh_coords;
  struct BuildContext build_ctx_data;
};


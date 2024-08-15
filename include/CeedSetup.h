#pragma once

#include <ceed.h>

#include "structs.h"

class CeedSetup
{
public:
  CeedSetup(const char *resource);
  ~CeedSetup();

  // Setup basis function for libceed
  void createCeedBasis(FEproblemData feproblem_data);

  // Determine the mesh size based on the given approximate problem size.
  void getCartesianMeshSize(FEproblemData feproblem_data);

  // Build CeedElemRestriction objects describing the mesh and solution discrete representations.
  void buildCartesianRestriction(FEproblemData feproblem_data, CeedInt num_comp, CeedInt *size,
                                 CeedElemRestriction *restriction, CeedElemRestriction *q_data_restriction);
  // Create a CeedVector with the mesh coordinates.
  void SetCartesianMeshCoords(FEproblemData feproblem_data, CeedData ceed_data);

  // Apply a transformation to the mesh.
  CeedScalar TransformMeshCoords(FEproblemData feproblem_data, CeedData ceed_data);

  // Method to setup QFunction
  void setupQfunction(FEproblemData feproblem_data, CeedData ceed_data);

  // Method to setup Operator
  void setupOperator(FEproblemData feproblem_data, CeedData ceed_data);

  CeedScalar solve(FEproblemData feproblem_data, CeedData ceed_data);
  // Getters for libCEED objects
  Ceed &getCeed()
  {
    return _ceed;
  }
  void cleanup(CeedData ceed_data);

private:
  CeedData ceed_data;
  Ceed _ceed;
  CeedQFunction qf;
  CeedOperator op;
  CeedVector x, b;
  CeedBasis _mesh_basis, _sol_basis;
  CeedInt _mesh_size, _sol_size;
};

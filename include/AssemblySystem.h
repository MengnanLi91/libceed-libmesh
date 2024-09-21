#pragma once

#include <math.h>

// Basic include files needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/vtk_io.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/equation_systems.h"

// Define the Finite Element object.
#include "libmesh/fe.h"

// Define Gauss quadrature rules.
#include "libmesh/quadrature_gauss.h"

// Define useful datatypes for finite element
// matrix and vector components.
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/elem.h"
#include "libmesh/enum_solver_package.h"

// Define the DofMap, which handles degree of freedom
// indexing.
#include "libmesh/dof_map.h"

// Include libCEED header
#include <ceed.h>

#include "structs.h"
#include "FEProblembase.h"
#include "CeedSetup.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

class CeedSetup;

void assemble_poisson_system(EquationSystems & es, const std::string & system_name);

class AssemblySystem
{
public:
  AssemblySystem();

  void assembleFunc(EquationSystems & equation_systems, const std::string & system_name);
  void assembleCEED(FEproblemData & feproblem_data, CeedSetup & ceedsetup);

private:
  Real exact_solution(const Real x, const Real y, const Real z = 0.);

  FEProblembase _fe;
};

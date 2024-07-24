#include <iostream>
#include <algorithm>
#include <math.h>

// Basic include files needed for the mesh functionality.
#include "ceed.h"
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"

#include "mesh_setup.h"
#include "solver.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Function prototype.  This is the function that will assemble
// the linear system for our Poisson problem.  Note that the
// function will take the  EquationSystems object and the
// name of the system we are assembling as input.  From the
//  EquationSystems object we have access to the  Mesh and
// other objects we might need.
void assemble_poisson(EquationSystems &es,
                      const std::string &system_name);

// Function prototype for the exact solution.
// Real exact_solution(const Real x,
//                     const Real y,
//                     const Real z = 0.);

int main(int argc, char **argv)
{
    // Initialize libMesh
    LibMeshInit init(argc, argv);

    // Create a mesh
    Mesh mesh(init.comm());
    setup_mesh(mesh);

    solve_system(mesh);

    // Initialize libCEED
    Ceed ceed;
    CeedInit(argv[1], &ceed);

    // Create a basis
    CeedBasis basis;
    CeedInt P = 3; // Quadratic elements
    CeedInt Q = 4; // Number of quadrature points
    CeedBasisCreateTensorH1Lagrange(ceed, dim, 1, P, Q, CEED_GAUSS, &basis);

    // Clean up libCEED
    CeedBasisDestroy(&basis);
    CeedDestroy(&ceed);

    return 0;
}

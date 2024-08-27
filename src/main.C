#include <iostream>
#include <algorithm>
#include <math.h>

// Basic include files needed for the mesh functionality.
#include "ceed.h"
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"

#include "SetupMesh.h"
#include "LinearSystem.h"
#include "ProblemBase.h"
#include "CeedSetup.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

int main(int argc, char **argv)
{
    // Initialize libMesh
    LibMeshInit init(argc, argv);

    // Create a CeedSetup object for libCEED initialization and setup
    CeedSetup ceed_setup(argv[1]);

    // Create and setup mesh using SetupMesh class
    SetupMesh setup_mesh(init, 2);
    setup_mesh.create_square_mesh(100); // 10x10 elements
    setup_mesh.prepare_for_use();

    Mesh &mesh = setup_mesh.get_mesh();

    // Setup and solve equation systems
    const std::vector<std::string> &system_names{"Poisson"};

    ProblemBase possion_problem(mesh, system_names, ceed_setup);

    possion_problem.initialSetup();

    possion_problem.solve();

    return 0;
}

// solver.c
#include <libmesh/libmesh.h>
#include <libmesh/equation_systems.h>
#include <libmesh/linear_implicit_system.h>
#include <libmesh/dof_map.h>
#include <libmesh/fe.h>
#include <libmesh/fe_interface.h>
#include <libmesh/quadrature_gauss.h>
#include <libmesh/exodusII_io.h>
#include "solver.h"

void solve_system(libMesh::Mesh &mesh) {
    // Define Equation Systems
    libMesh::EquationSystems equation_systems(mesh);
    libMesh::LinearImplicitSystem &system = equation_systems.add_system<libMesh::LinearImplicitSystem>("Poisson");

    // Add variables
    libMesh::FEType fe_type(libMesh::LAGRANGE, libMesh::FIRST);
    system.add_variable("u", fe_type);

    // Initialize the system
    equation_systems.init();

    // Get the degree of freedom map
    libMesh::DofMap &dof_map = system.get_dof_map();
    const unsigned int u_var = system.variable_number("u");

    // Apply Dirichlet boundary conditions (u = 0 on boundary)
    libMesh::MeshBase::const_node_iterator it = mesh.nodes_begin();
    const libMesh::MeshBase::const_node_iterator end = mesh.nodes_end();
    for (; it != end; ++it) {
        const libMesh::Node &node = **it;
        if (node.on_boundary()) {
            dof_map.add_dirichlet_boundary(libMesh::DirichletBoundary({ node.id() }, { &system }, { u_var }, 0.0));
        }
    }

    // Assemble the system matrix and right-hand side
    system.assemble_before_solve = false;
    system.assemble();

    // Solve the system
    system.solve();

    // Output the results
    libMesh::ExodusII_IO exodus_io(mesh);
    exodus_io.write_equation_systems("output.e", equation_systems);
}


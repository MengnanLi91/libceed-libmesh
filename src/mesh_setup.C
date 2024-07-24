// mesh_setup.c
#include <libmesh/libmesh.h>
#include <libmesh/mesh_tools.h>
#include "mesh_setup.h"

void setup_mesh(libMesh::Mesh &mesh) {
    const unsigned int dim = 2;
    mesh.set_mesh_dimension(dim);
    libMesh::MeshTools::Generation::build_square(mesh, 10, 10, 0.0, 1.0, 0.0, 1.0, libMesh::QUAD4);
    mesh.prepare_for_use();
}


#include "SetupMesh.h"
#include <libmesh/mesh_tools.h>

SetupMesh::SetupMesh(LibMeshInit &init, unsigned int dim)
    : mesh(init.comm()), dimension(dim)
{
    mesh.set_mesh_dimension(dimension);
}

void SetupMesh::create_square_mesh(unsigned int elements_per_side)
{
    libMesh::MeshTools::Generation::build_square(mesh,
                                                 elements_per_side, elements_per_side,
                                                 -1, 1.0, -1, 1.0,
                                                 libMesh::QUAD9);
}

void SetupMesh::prepare_for_use()
{
    mesh.prepare_for_use();
}

libMesh::Mesh &SetupMesh::get_mesh()
{
    return mesh;
}

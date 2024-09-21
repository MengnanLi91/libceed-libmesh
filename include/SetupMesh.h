#pragma once

#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
using namespace libMesh;

class SetupMesh
{
public:
  SetupMesh(LibMeshInit & init, unsigned int dim = 2);
  void create_square_mesh(unsigned int elements_per_side);
  void prepare_for_use();
  Mesh & get_mesh();

private:
  Mesh mesh;
  unsigned int dimension;
};

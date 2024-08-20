
#include <stdexcept> // For exception handling
#include <iostream>

#include "CeedUtils.h"

void print_FEproblemData(const FEproblemData_ &data)
{
  std::cout << "dim: " << data.dim << std::endl;
  std::cout << "num_comp: " << data.num_comp << std::endl;
  std::cout << "num_poly: " << data.num_poly << std::endl;
  std::cout << "num_qpts: " << data.num_qpts << std::endl;

  // Print out num_xyz array
  std::cout << "num_xyz: [" << data.num_xyz[0] << ", " << data.num_xyz[1] << ", " << data.num_xyz[2] << "]" << std::endl;

  std::cout << "num_dofs: " << data.num_dofs << std::endl;
  std::cout << "mesh_size: " << data.mesh_size << std::endl;
  std::cout << "sol_size: " << data.sol_size << std::endl;
}

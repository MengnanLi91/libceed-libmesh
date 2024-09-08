
#include <stdexcept> // For exception handling
#include <iostream>

#include "CeedUtils.h"

void
print_FEproblemData(const FEproblemData_ & data)
{
  std::cout << "dim: " << data.dim << std::endl;
  std::cout << "num_comp: " << data.num_comp << std::endl;
  std::cout << "num_poly: " << data.num_poly << std::endl;
  std::cout << "num_qpts: " << data.num_qpts << std::endl;

  // Print out num_xyz array
  std::cout << "num_xyz: [" << data.num_xyz[0] << ", " << data.num_xyz[1] << ", " << data.num_xyz[2]
            << "]" << std::endl;

  std::cout << "num_dofs: " << data.num_dofs << std::endl;
  std::cout << "mesh_size: " << data.mesh_size << std::endl;
  std::cout << "sol_size: " << data.sol_size << std::endl;
}

void
printCeedVector(CeedVector vec)
{
  const CeedScalar * array;
  CeedSize length;

  // Get the length of the CeedVector
  CeedVectorGetLength(vec, &length);

  // Get read access to the CeedVector array
  CeedVectorGetArrayRead(vec, CEED_MEM_HOST, &array);

  // Print each element of the vector
  for (CeedSize i = 0; i < length; i++)
  {
    printf("Element %td: %f\n", i, array[i]);
  }

  // Restore the CeedVector array after reading
  CeedVectorRestoreArrayRead(vec, &array);
}

void
verifyQFunctionContext(CeedQFunctionContext build_ctx)
{
  // verify the data was set correctly
  BuildContext * retrieved_data;
  CeedQFunctionContextGetData(build_ctx, CEED_MEM_HOST, (void **)&retrieved_data);
  std::cout << "Retrieved dim: " << retrieved_data->dim << std::endl;
  std::cout << "Retrieved space_dim: " << retrieved_data->space_dim << std::endl;
  CeedQFunctionContextRestoreData(build_ctx, (void **)&retrieved_data);
}

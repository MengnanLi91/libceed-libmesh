#include <libmesh/exodusII_io.h>

#include "ProblemBase.h"
#include "AssemblySystem.h"
#include "CeedSetup.h"

using namespace libMesh;

class LinearSystem;

// Define the static member outside the class
AssemblySystem *AssemblerContext::assembler_instance = nullptr;

ProblemBase::ProblemBase(Mesh &mesh, const std::vector<std::string> &system_names, CeedSetup &ceed_setup) : _mesh(mesh), _ln_sys_names(system_names), _num_ln_sys(_ln_sys_names.size()), _ln_sys(_num_ln_sys, nullptr), _ceed_setup(ceed_setup)
{
  // Store a pointer to the assembler in the context
  // AssemblerContext::assembler_instance = &_assembler;

  _es = std::make_shared<EquationSystems>(_mesh);
  libmesh_assert_msg(_num_ln_sys, "should have at least one system to solve");

  if (_num_ln_sys)
    for (const auto i : index_range(_ln_sys_names))
    {
      const auto &sys_name = _ln_sys_names[i];
      _ln_sys[i] = std::make_shared<LinearImplicitSystem>(*(_es.get()), sys_name, i);
      _ln_sys_name_to_num[sys_name] = i;
    }
}
void ProblemBase::initialSetup()
{
  std::cout << "Initializing systems..." << std::endl;
  for (const auto &name : _ln_sys_names)
  {
    _es->add_system<LinearImplicitSystem>(name);
    FEType fe_type(SECOND, LAGRANGE);
    _es->get_system<LinearImplicitSystem>(name).add_variable(name, fe_type);
    _es->get_system<LinearImplicitSystem>(name).attach_assemble_function(assemble_poisson_system);
  }
  std::cout << "Number of systems after adding: " << _es->n_systems() << std::endl;

  _es->init();
  _es->print_info();
}
void ProblemBase::solve()
{
  std::cout << "Solving the systems..." << std::endl;
  for (const auto &name : _ln_sys_names)
  {
    _es->get_system<LinearImplicitSystem>(name).solve();
  }

  std::cout << "Number of systems: " << _es->n_systems() << std::endl;
}

void ProblemBase::printInfo()
{
  for (const auto &name : _ln_sys_names)
  {
    _es->get_system<LinearImplicitSystem>(name).print_info();
  }
}

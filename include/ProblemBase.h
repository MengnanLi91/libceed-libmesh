#pragma once

#include <vector>

#include <libmesh/equation_systems.h>
#include <libmesh/mesh.h>
#include <libmesh/linear_implicit_system.h>

#include "LinearSystem.h"
#include "AssemblySystem.h"

using namespace libMesh;

class LinearSystem;

class ProblemBase
{
public:
  ProblemBase(Mesh &mesh, const std::vector<std::string> &system_names, CeedSetup &ceed_setup);

  void initialSetup();
  void solve();
  void printInfo();

  EquationSystems &es() { return *_es; }

  std::shared_ptr<LinearImplicitSystem> getSystembyName(const std::string &name)
  {

    libmesh_assert_msg(_ln_sys_name_to_num.find(name) != _ln_sys_name_to_num.end(), "The system name should have already be created and stored");

    return _ln_sys[_ln_sys_name_to_num[name]];
  }

protected:
  Mesh &_mesh;
  /// The nonlinear system names
  const std::vector<std::string> _ln_sys_names;

  /// The number of nonlinear systems
  const std::size_t _num_ln_sys;

  /// The nonlinear systems
  std::vector<std::shared_ptr<LinearImplicitSystem>> _ln_sys;

  /// Map from nonlinear system name to number
  std::map<const std::string, unsigned int> _ln_sys_name_to_num;

  std::shared_ptr<EquationSystems> _es;

  // std::shared_ptr<LinearSystem> _current_ln;

  AssemblySystem _assembler;
  CeedSetup _ceed_setup;
};

// Define the context structure with the static assembler instance pointer
struct AssemblerContext
{
  static void assemble(libMesh::EquationSystems &es, const std::string &system_name)
  {
    assembler_instance->assembleFunc(es, system_name);
  }
  static AssemblySystem *assembler_instance;
};

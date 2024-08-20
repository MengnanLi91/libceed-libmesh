#pragma once

#include <libmesh/equation_systems.h>
#include <libmesh/mesh.h>
#include <libmesh/linear_implicit_system.h>
#include "libmesh/nonlinear_implicit_system.h"

#include "ProblemBase.h"
#include "AssemblySystem.h"
#include "structs.h"

using namespace libMesh;

class ProblemBase;

class LinearSystem : public LinearImplicitSystem
{
public:
  LinearSystem(EquationSystems &es, const std::string &name, const unsigned int number);

  LinearSystem getSystem(const std::string &name);
  void initialSetup();
  void solveSystem(const std::string &name, CeedSetup &ceedsetup);
  void outputResults(const std::string &filename);
  void addVariable(const std::string &name, Order order, FEFamily family);
  void printInfo();

private:
  EquationSystems &_es;
  const std::string &_name;

  // Create libCEED vectors for the system solution and rhs
  // CeedData _ceed_data;
  FEproblemData _feproblem_data;

  // std::shared_ptr<AssemblySystem> _assembly;
};

// C++ Includes
#include <math.h>

#include "LinearSystem.h"
#include <libmesh/dof_map.h>
#include <libmesh/fe.h>
#include <libmesh/fe_type.h>
#include <libmesh/quadrature_gauss.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/node.h>
#include <libmesh/linear_implicit_system.h>
#include <libmesh/libmesh_common.h>

#include "ProblemBase.h"
#include "AssemblySystem.h"
#include "CeedUtils.h"

LinearSystem::LinearSystem(EquationSystems & es,
                           const std::string & name,
                           const unsigned int number)
  : LinearImplicitSystem(es, name, number), _es(es), _name(name), _feproblem_data()
{
}

void
LinearSystem::initialSetup()
{
  addVariable("u", SECOND, LAGRANGE);
  addVariable("v", SECOND, LAGRANGE);
  this->attach_assemble_function(assemble_poisson_system);
}

void
LinearSystem::addVariable(const std::string & name, Order order, FEFamily family)
{
  FEType fe_type(order, family);
  this->add_variable(name, fe_type);
}

void
LinearSystem::solveSystem(const std::string & /*name*/, CeedSetup & ceedsetup)
{
  AssemblySystem assembly_system;
  // Get the polynomial order
  _feproblem_data.num_dofs = this->n_dofs();
  _feproblem_data.num_poly = this->variable_type(0).order;
  _feproblem_data.num_comp = this->n_components();
  _feproblem_data.quad_mode = CEED_GAUSS;
  auto dim = _es.get_mesh().mesh_dimension();
  _feproblem_data.dim = dim;
  _feproblem_data.num_qpts = this->variable_type(0).order + 2;
  // qrule.n_points();

  print_FEproblemData(_feproblem_data);
  assembly_system.assembleCEED(_feproblem_data, ceedsetup);
  // this->solve();
}

void
LinearSystem::printInfo()
{
  this->print_info();
}

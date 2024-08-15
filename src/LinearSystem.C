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

LinearSystem::LinearSystem(EquationSystems &es, const std::string &name, const unsigned int number) : LinearImplicitSystem(es, name, number), _es(es), _name(name)
{
}

void LinearSystem::initialSetup(CeedSetup &ceedsetup)
{
    AssemblySystem assembly_system;
    addVariable("u", SECOND, LAGRANGE);
    this->attach_assemble_function(assemble_poisson_system);
    this->init();
    // assembly_system.assembleCEED(_es, _name, ceedsetup);
}

void LinearSystem::addVariable(const std::string &name, Order order, FEFamily family)
{
    FEType fe_type(order, family);
    this->add_variable(name, fe_type);
}

void LinearSystem::solveSystem()
{
    this->solve();
}

void LinearSystem::printInfo()
{
    this->print_info();
}

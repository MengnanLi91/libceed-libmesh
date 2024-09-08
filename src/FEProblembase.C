#include "FEProblembase.h"

FEProblembase::FEProblembase(unsigned int dim, libMesh::FEType fe_type, libMesh::Order order)
  : fe(libMesh::FEBase::build(dim, fe_type)), qrule(dim, order)
{
  fe->attach_quadrature_rule(&qrule);
}

void
FEProblembase::reinit(const libMesh::Elem * elem)
{
  fe->reinit(elem);
}

const std::vector<libMesh::Real> &
FEProblembase::get_JxW() const
{
  return fe->get_JxW();
}

const std::vector<std::vector<libMesh::Real>> &
FEProblembase::get_phi() const
{
  return fe->get_phi();
}

const std::vector<std::vector<libMesh::RealGradient>> &
FEProblembase::get_dphi() const
{
  return fe->get_dphi();
}

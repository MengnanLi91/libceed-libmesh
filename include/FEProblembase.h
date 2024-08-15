#pragma once

#include <libmesh/fe_base.h>
#include <libmesh/quadrature_gauss.h>
#include <libmesh/elem.h>
#include <libmesh/dof_map.h>
#include <libmesh/mesh.h>

class FEProblembase
{
public:
  FEProblembase(unsigned int dim, libMesh::FEType fe_type, libMesh::Order order);

  void reinit(const libMesh::Elem *elem);
  const std::vector<libMesh::Real> &get_JxW() const;
  const std::vector<std::vector<libMesh::Real>> &get_phi() const;
  const std::vector<std::vector<libMesh::RealGradient>> &get_dphi() const;

private:
  std::unique_ptr<libMesh::FEBase> fe;
  libMesh::QGauss qrule;
};

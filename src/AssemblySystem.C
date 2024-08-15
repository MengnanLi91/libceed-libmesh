// C++ Includes
// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <libmesh/quadrature.h>
#include <libmesh/quadrature_gauss.h>

#include "AssemblySystem.h"
#include "LinearSystem.h"

void assemble_poisson_system(EquationSystems &es, const std::string &system_name)
{
    // Create an instance of the Assembly class
    AssemblySystem assembly;
    // Call the assembly function from the Assembly class
    assembly.assembleFunc(es, system_name);
}

AssemblySystem::AssemblySystem() : _fe(2, libMesh::FEType(libMesh::SECOND, libMesh::LAGRANGE), libMesh::SECOND)
{
}

void AssemblySystem::assembleCEED(libMesh::EquationSystems &equation_systems, const std::string &system_name, CeedSetup &ceedsetup)
{
    // It is a good idea to make sure we are assembling
    // the proper system.
    libmesh_assert_equal_to(system_name, "Poisson");

    // Get a constant reference to the mesh object.
    const MeshBase &mesh = equation_systems.get_mesh();

    // The dimension that we are running
    const unsigned int dim = mesh.mesh_dimension();

    // Get a reference to the LinearImplicitSystem we are solving
    LinearImplicitSystem &system = equation_systems.get_system<LinearImplicitSystem>("Poisson");

    // A reference to the  DofMap object for this system.  The  DofMap
    // object handles the index translation from node and element numbers
    // to degree of freedom numbers.  We will talk more about the  DofMap
    // in future examples.
    const DofMap &dof_map = system.get_dof_map();
    std::vector<dof_id_type> dof_indices;

    // Get a constant reference to the Finite Element type
    // for the first (and only) variable in the system.
    FEType fe_type = dof_map.variable_type(0);

    // A 5th order Gauss quadrature rule for numerical integration.
    QGauss qrule(dim, FIFTH);

    // Create libCEED vectors for the system solution and rhs
    CeedData ceed_data = nullptr;
    FEproblemData feproblem_data = nullptr;

    feproblem_data->dim = dim;
    feproblem_data->num_dofs = system.n_dofs();
    feproblem_data->num_poly = static_cast<int>(fe_type.order);
    feproblem_data->num_comp = system.n_components();
    feproblem_data->num_qpts = qrule.n_points();
    feproblem_data->quad_mode = CEED_GAUSS;

    ceedsetup.createCeedBasis(feproblem_data);
    ceedsetup.getCartesianMeshSize(feproblem_data);

    printf("Mesh size: nx = %" CeedInt_FMT, feproblem_data->num_xyz[0]);
    if (dim > 1)
        printf(", ny = %" CeedInt_FMT, feproblem_data->num_xyz[1]);
    if (dim > 2)
        printf(", nz = %" CeedInt_FMT, feproblem_data->num_xyz[2]);
    printf("\n");

    ceedsetup.buildCartesianRestriction(feproblem_data, system.n_components(), &feproblem_data->mesh_size, &ceed_data->elem_restr_x, NULL);
    ceedsetup.buildCartesianRestriction(feproblem_data, dim * (dim + 1) / 2, &feproblem_data->sol_size, NULL, &ceed_data->elem_restr_qd);
    ceedsetup.buildCartesianRestriction(feproblem_data, 1, &feproblem_data->sol_size, &ceed_data->elem_restr_u, NULL);

    // LCOV_EXCL_START
    printf("Number of mesh nodes     : %" CeedInt_FMT "\n", feproblem_data->mesh_size / dim);
    printf("Number of solution nodes : %" CeedInt_FMT "\n", feproblem_data->sol_size);
    // LCOV_EXCL_STOP

    ceedsetup.SetCartesianMeshCoords(feproblem_data, ceed_data);

    CeedScalar exact_surface_area = ceedsetup.TransformMeshCoords(feproblem_data, ceed_data);

    ceedsetup.setupQfunction(feproblem_data, ceed_data);
    ceedsetup.setupOperator(feproblem_data, ceed_data);

    CeedScalar surface_area = ceedsetup.solve(feproblem_data, ceed_data);

    printf("Exact mesh surface area    : % .14g\n", exact_surface_area);
    printf("Surface area error         : % .14g\n", surface_area - exact_surface_area);
}

Real AssemblySystem::exact_solution(const Real x,
                                    const Real y,
                                    const Real z)
{
    static const Real pi = acos(-1.);

    return cos(.5 * pi * x) * sin(.5 * pi * y) * cos(.5 * pi * z);
}

void AssemblySystem::assembleFunc(libMesh::EquationSystems &equation_systems, const std::string &system_name)
{
    // It is a good idea to make sure we are assembling
    // the proper system.
    libmesh_assert_equal_to(system_name, "Poisson");

    // Get a constant reference to the mesh object.
    const MeshBase &mesh = equation_systems.get_mesh();

    // The dimension that we are running
    const unsigned int dim = mesh.mesh_dimension();

    // Get a reference to the LinearImplicitSystem we are solving
    LinearImplicitSystem &system = equation_systems.get_system<LinearImplicitSystem>("Poisson");

    // A reference to the  DofMap object for this system.  The  DofMap
    // object handles the index translation from node and element numbers
    // to degree of freedom numbers.  We will talk more about the  DofMap
    // in future examples.
    const DofMap &dof_map = system.get_dof_map();

    // Get a constant reference to the Finite Element type
    // for the first (and only) variable in the system.
    FEType fe_type = dof_map.variable_type(0);

    // Build a Finite Element object of the specified type.  Since the
    // FEBase::build() member dynamically creates memory we will
    // store the object as a std::unique_ptr<FEBase>.  This can be thought
    // of as a pointer that will clean up after itself.  Introduction Example 4
    // describes some advantages of  std::unique_ptr's in the context of
    // quadrature rules.
    std::unique_ptr<FEBase> fe(FEBase::build(dim, fe_type));

    // A 5th order Gauss quadrature rule for numerical integration.
    QGauss qrule(dim, FIFTH);

    // Tell the finite element object to use our quadrature rule.
    fe->attach_quadrature_rule(&qrule);

    // Declare a special finite element object for
    // boundary integration.
    std::unique_ptr<FEBase> fe_face(FEBase::build(dim, fe_type));

    // Boundary integration requires one quadrature rule,
    // with dimensionality one less than the dimensionality
    // of the element.
    QGauss qface(dim - 1, FIFTH);

    // Tell the finite element object to use our
    // quadrature rule.
    fe_face->attach_quadrature_rule(&qface);

    // Here we define some references to cell-specific data that
    // will be used to assemble the linear system.
    //
    // The element Jacobian * quadrature weight at each integration point.
    const std::vector<Real> &JxW = fe->get_JxW();

    // The physical XY locations of the quadrature points on the element.
    // These might be useful for evaluating spatially varying material
    // properties at the quadrature points.
    const std::vector<Point> &q_point = fe->get_xyz();

    // The element shape functions evaluated at the quadrature points.
    const std::vector<std::vector<Real>> &phi = fe->get_phi();

    // The element shape function gradients evaluated at the quadrature
    // points.
    const std::vector<std::vector<RealGradient>> &dphi = fe->get_dphi();

    // Define data structures to contain the element matrix
    // and right-hand-side vector contribution.  Following
    // basic finite element terminology we will denote these
    // "Ke" and "Fe".  These datatypes are templated on
    //  Number, which allows the same code to work for real
    // or complex numbers.
    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;

    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector<dof_id_type> dof_indices;

    // The global system matrix
    SparseMatrix<Number> &matrix = system.get_system_matrix();

    // Now we will loop over all the elements in the mesh.
    // We will compute the element matrix and right-hand-side
    // contribution.
    //
    // Element ranges are a nice way to iterate through all the
    // elements, or all the elements that have some property.  The
    // range will iterate from the first to the last element on
    // the local processor.
    // It is smart to make this one const so that we don't accidentally
    // mess it up!  In case users later modify this program to include
    // refinement, we will be safe and will only consider the active
    // elements; hence we use a variant of the
    // active_local_element_ptr_range.
    for (const auto &elem : mesh.active_local_element_ptr_range())
    {
        // Get the degree of freedom indices for the
        // current element.  These define where in the global
        // matrix and right-hand-side this element will
        // contribute to.
        dof_map.dof_indices(elem, dof_indices);

        // Cache the number of degrees of freedom on this element, for
        // use as a loop bound later.  We use cast_int to explicitly
        // convert from size() (which may be 64-bit) to unsigned int
        // (which may be 32-bit but which is definitely enough to count
        // *local* degrees of freedom.
        const unsigned int n_dofs =
            cast_int<unsigned int>(dof_indices.size());

        // Compute the element-specific data for the current
        // element.  This involves computing the location of the
        // quadrature points (q_point) and the shape functions
        // (phi, dphi) for the current element.
        fe->reinit(elem);

        // With one variable, we should have the same number of degrees
        // of freedom as shape functions.
        libmesh_assert_equal_to(n_dofs, phi.size());

        // Zero the element matrix and right-hand side before
        // summing them.  We use the resize member here because
        // the number of degrees of freedom might have changed from
        // the last element.  Note that this will be the case if the
        // element type is different (i.e. the last element was a
        // triangle, now we are on a quadrilateral).

        // The  DenseMatrix::resize() and the  DenseVector::resize()
        // members will automatically zero out the matrix  and vector.
        Ke.resize(n_dofs, n_dofs);

        Fe.resize(n_dofs);

        // Now loop over the quadrature points.  This handles
        // the numeric integration.
        for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
        {

            // Now we will build the element matrix.  This involves
            // a double loop to integrate the test functions (i) against
            // the trial functions (j).
            for (unsigned int i = 0; i != n_dofs; i++)
                for (unsigned int j = 0; j != n_dofs; j++)
                {
                    Ke(i, j) += JxW[qp] * (dphi[i][qp] * dphi[j][qp]);
                }

            // This is the end of the matrix summation loop
            // Now we build the element right-hand-side contribution.
            // This involves a single loop in which we integrate the
            // "forcing function" in the PDE against the test functions.
            {
                const Real x = q_point[qp](0);
                const Real y = q_point[qp](1);
                const Real eps = 1.e-3;

                // "fxy" is the forcing function for the Poisson equation.
                // In this case we set fxy to be a finite difference
                // Laplacian approximation to the (known) exact solution.
                //
                // We will use the second-order accurate FD Laplacian
                // approximation, which in 2D is
                //
                // u_xx + u_yy = (u(i,j-1) + u(i,j+1) +
                //                u(i-1,j) + u(i+1,j) +
                //                -4*u(i,j))/h^2
                //
                // Since the value of the forcing function depends only
                // on the location of the quadrature point (q_point[qp])
                // we will compute it here, outside of the i-loop
                const Real fxy = -(exact_solution(x, y - eps) +
                                   exact_solution(x, y + eps) +
                                   exact_solution(x - eps, y) +
                                   exact_solution(x + eps, y) -
                                   4. * exact_solution(x, y)) /
                                 eps / eps;

                for (unsigned int i = 0; i != n_dofs; i++)
                    Fe(i) += JxW[qp] * fxy * phi[i][qp];
            }
        }

        // We have now reached the end of the RHS summation,
        // and the end of quadrature point loop, so
        // the interior element integration has
        // been completed.  However, we have not yet addressed
        // boundary conditions.  For this example we will only
        // consider simple Dirichlet boundary conditions.
        //
        // There are several ways Dirichlet boundary conditions
        // can be imposed.  A simple approach, which works for
        // interpolary bases like the standard Lagrange polynomials,
        // is to assign function values to the
        // degrees of freedom living on the domain boundary. This
        // works well for interpolary bases, but is more difficult
        // when non-interpolary (e.g Legendre or Hierarchic) bases
        // are used.
        //
        // Dirichlet boundary conditions can also be imposed with a
        // "penalty" method.  In this case essentially the L2 projection
        // of the boundary values are added to the matrix. The
        // projection is multiplied by some large factor so that, in
        // floating point arithmetic, the existing (smaller) entries
        // in the matrix and right-hand-side are effectively ignored.
        //
        // This amounts to adding a term of the form (in latex notation)
        //
        // \frac{1}{\epsilon} \int_{\delta \Omega} \phi_i \phi_j = \frac{1}{\epsilon} \int_{\delta \Omega} u \phi_i
        //
        // where
        //
        // \frac{1}{\epsilon} is the penalty parameter, defined such that \epsilon << 1
        {

            // The following loop is over the sides of the element.
            // If the element has no neighbor on a side then that
            // side MUST live on a boundary of the domain.
            for (auto side : elem->side_index_range())
                if (elem->neighbor_ptr(side) == nullptr)
                {
                    // The value of the shape functions at the quadrature
                    // points.
                    const std::vector<std::vector<Real>> &phi_face = fe_face->get_phi();

                    // The Jacobian * Quadrature Weight at the quadrature
                    // points on the face.
                    const std::vector<Real> &JxW_face = fe_face->get_JxW();

                    // The XYZ locations (in physical space) of the
                    // quadrature points on the face.  This is where
                    // we will interpolate the boundary value function.
                    const std::vector<Point> &qface_point = fe_face->get_xyz();

                    // Compute the shape function values on the element
                    // face.
                    fe_face->reinit(elem, side);

                    // Some shape functions will be 0 on the face, but for
                    // ease of indexing and generality of code we loop over
                    // them anyway
                    libmesh_assert_equal_to(n_dofs, phi_face.size());

                    // Loop over the face quadrature points for integration.
                    for (unsigned int qp = 0; qp < qface.n_points(); qp++)
                    {
                        // The location on the boundary of the current
                        // face quadrature point.
                        const Real xf = qface_point[qp](0);
                        const Real yf = qface_point[qp](1);

                        // The penalty value.  \frac{1}{\epsilon}
                        // in the discussion above.
                        const Real penalty = 1.e10;

                        // The boundary value.
                        const Real value = exact_solution(xf, yf);

                        // Matrix contribution of the L2 projection.
                        for (unsigned int i = 0; i != n_dofs; i++)
                            for (unsigned int j = 0; j != n_dofs; j++)
                                Ke(i, j) += JxW_face[qp] * penalty * phi_face[i][qp] * phi_face[j][qp];

                        // Right-hand-side contribution of the L2
                        // projection.
                        for (unsigned int i = 0; i != n_dofs; i++)
                            Fe(i) += JxW_face[qp] * penalty * value * phi_face[i][qp];
                    }
                }
        }

        // We have now finished the quadrature point loop,
        // and have therefore applied all the boundary conditions.

        // If this assembly program were to be used on an adaptive mesh,
        // we would have to apply any hanging node constraint equations
        dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);

        // The element matrix and right-hand-side are now built
        // for this element.  Add them to the global matrix and
        // right-hand-side vector.  The  SparseMatrix::add_matrix()
        // and  NumericVector::add_vector() members do this for us.
        matrix.add_matrix(Ke, dof_indices);
        system.rhs->add_vector(Fe, dof_indices);
    }

    // All done!
}

// LIC// ====================================================================
// LIC// Algebraic direct construction of the Bernadou curved Bell basis.
// LIC// ====================================================================

#include "c1_curved_elements_direct.h"

namespace oomph
{
  /// Constructor. The algebraic transformation is built after
  /// upgrade_element().
  template<unsigned BOUNDARY_ORDER>
  DirectBernadouElementBasis<BOUNDARY_ORDER>::DirectBernadouElementBasis()
    : Base(), Algebraic_transformation_is_built(false)
  {
  }


  /// Upgrade the geometry and construct T and C=T^T algebraically.
  template<unsigned BOUNDARY_ORDER>
  void DirectBernadouElementBasis<BOUNDARY_ORDER>::upgrade_element(
    const VertexList& verts,
    const double& su,
    const double& so,
    const C1PlateHelper::CurvedEdgeEnumeration& curved_edge,
    const C1CurviLine& parametric_curve)
  {
    Base::upgrade_element(verts, su, so, curved_edge, parametric_curve);
    Algebraic_transformation_is_built = false;
    build_algebraic_transformation();
  }


  /// Quintic Hermite basis ordered as
  /// [f(0),f(1),f'(0),f'(1),f''(0),f''(1)].
  template<unsigned BOUNDARY_ORDER>
  void DirectBernadouElementBasis<BOUNDARY_ORDER>::quintic_hermite(
    const double& u, Vector<double>& h)
  {
    h.resize(6);
    h[0] = std::pow(1.0 - u, 3) * (1.0 + 3.0 * u + 6.0 * u * u);
    h[1] = std::pow(u, 3) * (10.0 - 15.0 * u + 6.0 * u * u);
    h[2] = std::pow(1.0 - u, 3) * u * (1.0 + 3.0 * u);
    h[3] = (1.0 - u) * std::pow(u, 3) * (3.0 * u - 4.0);
    h[4] = 0.5 * std::pow(1.0 - u, 3) * u * u;
    h[5] = 0.5 * std::pow(1.0 - u, 2) * std::pow(u, 3);
  }


  /// Cubic Hermite basis ordered as
  /// [g(0),g(1),g'(0),g'(1)].
  template<unsigned BOUNDARY_ORDER>
  void DirectBernadouElementBasis<BOUNDARY_ORDER>::cubic_hermite(
    const double& u, Vector<double>& h)
  {
    h.resize(4);
    h[0] = 1.0 - 3.0 * u * u + 2.0 * u * u * u;
    h[1] = 3.0 * u * u - 2.0 * u * u * u;
    h[2] = u - 2.0 * u * u + u * u * u;
    h[3] = -u * u + u * u * u;
  }


  /// Set a coefficient vector to zero with the correct global-dof length.
  template<unsigned BOUNDARY_ORDER>
  void DirectBernadouElementBasis<BOUNDARY_ORDER>::zero_global_dof_coefficients(
    Vector<double>& coefficients) const
  {
    coefficients = Vector<double>(this->n_basis_functions(), 0.0);
  }


  /// Copy a coefficient vector into one row of a matrix.
  template<unsigned BOUNDARY_ORDER>
  void DirectBernadouElementBasis<BOUNDARY_ORDER>::put_row(
    DenseMatrix<double>& matrix,
    const unsigned& row,
    const Vector<double>& coefficients)
  {
    for (unsigned j = 0; j < coefficients.size(); ++j)
    {
      matrix(row, j) = coefficients[j];
    }
  }


  /// Coefficients of w at a vertex in the global dof vector.
  template<unsigned BOUNDARY_ORDER>
  void DirectBernadouElementBasis<BOUNDARY_ORDER>::vertex_value_coefficients(
    const unsigned& vertex, Vector<double>& coefficients) const
  {
    zero_global_dof_coefficients(coefficients);
    coefficients[vertex] = 1.0;
  }


  /// Coefficients of a physical directional first derivative
  /// direction.grad_x(w) at a vertex.
  template<unsigned BOUNDARY_ORDER>
  void DirectBernadouElementBasis<BOUNDARY_ORDER>::
    vertex_physical_first_derivative_coefficients(
      const unsigned& vertex,
      const Vector<double>& direction,
      Vector<double>& coefficients) const
  {
    zero_global_dof_coefficients(coefficients);
    coefficients[3 + 2 * vertex] = direction[0];
    coefficients[4 + 2 * vertex] = direction[1];
  }


  /// Coefficients of the physical mixed second derivative
  /// a^T H_x(w) b at a vertex.
  template<unsigned BOUNDARY_ORDER>
  void DirectBernadouElementBasis<BOUNDARY_ORDER>::
    vertex_physical_second_derivative_coefficients(
      const unsigned& vertex,
      const Vector<double>& a,
      const Vector<double>& b,
      Vector<double>& coefficients) const
  {
    zero_global_dof_coefficients(coefficients);
    coefficients[9 + 3 * vertex] = a[0] * b[0];
    coefficients[10 + 3 * vertex] = a[0] * b[1] + a[1] * b[0];
    coefficients[11 + 3 * vertex] = a[1] * b[1];
  }


  /// Coefficients of q.grad_s(w) at a vertex.
  template<unsigned BOUNDARY_ORDER>
  void DirectBernadouElementBasis<BOUNDARY_ORDER>::
    vertex_reference_first_derivative_coefficients(
      const unsigned& vertex,
      const Vector<double>& q,
      Vector<double>& coefficients) const
  {
    const double vertex_coordinate[3][2] = {{1.0, 0.0}, {0.0, 1.0}, {0.0, 0.0}};

    Vector<double> s(2);
    s[0] = vertex_coordinate[vertex][0];
    s[1] = vertex_coordinate[vertex][1];

    DenseMatrix<double> jacobian(2, 2, 0.0);
    this->get_basic_jacobian(s, jacobian);

    Vector<double> physical_direction(2, 0.0);
    for (unsigned i = 0; i < 2; ++i)
    {
      physical_direction[i] = jacobian(i, 0) * q[0] + jacobian(i, 1) * q[1];
    }

    vertex_physical_first_derivative_coefficients(
      vertex, physical_direction, coefficients);
  }


  /// Coefficients of a^T H_s(w) b at a vertex, including the Hessian of the
  /// coordinate map.
  template<unsigned BOUNDARY_ORDER>
  void DirectBernadouElementBasis<BOUNDARY_ORDER>::
    vertex_reference_second_derivative_coefficients(
      const unsigned& vertex,
      const Vector<double>& a,
      const Vector<double>& b,
      Vector<double>& coefficients) const
  {
    const double vertex_coordinate[3][2] = {{1.0, 0.0}, {0.0, 1.0}, {0.0, 0.0}};

    Vector<double> s(2);
    s[0] = vertex_coordinate[vertex][0];
    s[1] = vertex_coordinate[vertex][1];

    DenseMatrix<double> jacobian(2, 2, 0.0);
    RankThreeTensor<double> map_hessian(2, 2, 2, 0.0);
    this->get_basic_jacobian(s, jacobian);
    this->get_basic_hessian(s, map_hessian);

    Vector<double> physical_a(2, 0.0);
    Vector<double> physical_b(2, 0.0);
    for (unsigned i = 0; i < 2; ++i)
    {
      physical_a[i] = jacobian(i, 0) * a[0] + jacobian(i, 1) * a[1];
      physical_b[i] = jacobian(i, 0) * b[0] + jacobian(i, 1) * b[1];
    }

    vertex_physical_second_derivative_coefficients(
      vertex, physical_a, physical_b, coefficients);

    // Chain-rule contribution from the Hessian of the coordinate map:
    //
    // a^T H_s(w) b = (J a)^T H_x(w) (J b)
    //                  + grad_x(w) . H_F[a,b]. <- This term
    for (unsigned i = 0; i < 2; ++i)
    {
      double map_second_derivative = 0.0;
      for (unsigned alpha = 0; alpha < 2; ++alpha)
      {
        for (unsigned beta = 0; beta < 2; ++beta)
        {
          map_second_derivative +=
            a[alpha] * b[beta] * map_hessian(i, alpha, beta);
        }
      }
      coefficients[3 + 2 * vertex + i] += map_second_derivative;
    }
  }



  /// Assemble T directly from the global and basic dof definitions, then
  /// transpose it to obtain C.
  template<unsigned BOUNDARY_ORDER>
  void DirectBernadouElementBasis<
    BOUNDARY_ORDER>::build_algebraic_transformation() const
  {
    const unsigned nbasis = this->n_basis_functions();
    const unsigned nbasic = this->n_basic_basis_functions();
    const unsigned ninternal = this->n_internal_dofs();
    const unsigned npoint_per_edge = BOUNDARY_ORDER - 1;
    const unsigned nmidnodes = 3 * npoint_per_edge;

    Basic_dofs_from_global_dofs.resize(nbasic, nbasis, 0.0);

    Vector<double> coefficients;

    // ------------------------------------------------------------------
    // 1. Basic vertex values.
    // ------------------------------------------------------------------
    for (unsigned vertex = 0; vertex < 3; ++vertex)
    {
      vertex_value_coefficients(vertex, coefficients);
      put_row(Basic_dofs_from_global_dofs, vertex, coefficients);
    }

    // ------------------------------------------------------------------
    // 2. Basic first and second reference derivatives at the vertices.
    // ------------------------------------------------------------------
    Vector<double> e0(2, 0.0), e1(2, 0.0);
    e0[0] = 1.0;
    e1[1] = 1.0;

    for (unsigned vertex = 0; vertex < 3; ++vertex)
    {
      vertex_reference_first_derivative_coefficients(vertex, e0, coefficients);
      put_row(Basic_dofs_from_global_dofs, 3 + 2 * vertex, coefficients);

      vertex_reference_first_derivative_coefficients(vertex, e1, coefficients);
      put_row(Basic_dofs_from_global_dofs, 4 + 2 * vertex, coefficients);

      vertex_reference_second_derivative_coefficients(
        vertex, e0, e0, coefficients);
      put_row(Basic_dofs_from_global_dofs, 9 + 3 * vertex, coefficients);

      vertex_reference_second_derivative_coefficients(
        vertex, e0, e1, coefficients);
      put_row(Basic_dofs_from_global_dofs, 10 + 3 * vertex, coefficients);

      vertex_reference_second_derivative_coefficients(
        vertex, e1, e1, coefficients);
      put_row(Basic_dofs_from_global_dofs, 11 + 3 * vertex, coefficients);
    }

    // ------------------------------------------------------------------
    // 3. Basic outward-normal derivatives at the three edge midpoints.
    // ------------------------------------------------------------------
    for (unsigned edge = 0; edge < 3; ++edge)
    {
      // [zdec] Method that fills out the normal derivative blocks
      // for rows
      //   18 + edge
    }

    // ------------------------------------------------------------------
    // 4. Additional basic edge values.
    // ------------------------------------------------------------------
    for (unsigned edge = 0; edge < 3; ++edge)
    {
      for (unsigned point = 0; point < npoint_per_edge; ++point)
      {
        // const double u = // [zdec] method to get the ordinate of point
        // [zdec] Method that fills out the value block for rows
        //   21 + edge * npoint_per_edge + point
      }
    }

    // ------------------------------------------------------------------
    // 5. Additional basic edge outward-normal derivatives.
    // ------------------------------------------------------------------
    for (unsigned edge = 0; edge < 3; ++edge)
    {
      for (unsigned point = 0; point < npoint_per_edge; ++point)
      {
        // const double u = // [zdec] method to get the ordinate of point
        // [zdec] Method that fills out the normal derivative block for
        // the rows
        //   21 + nmidnodes + edge * npoint_per_edge + point
      }
    }

    // ------------------------------------------------------------------
    // 6. Internal values are already global dofs.
    // ------------------------------------------------------------------
    for (unsigned i = 0; i < ninternal; ++i)
    {
      zero_global_dof_coefficients(coefficients);
      coefficients[18 + i] = 1.0;
      put_row(
        Basic_dofs_from_global_dofs, 21 + 2 * nmidnodes + i, coefficients);
    }

    if (21 + 2 * nmidnodes + ninternal != nbasic)
    {
      throw OomphLibError("Internal error in algebraic direct dof count.",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

    // C=T^T.
    Direct_basic_association_matrix.resize(nbasis, nbasic, 0.0);
    for (unsigned i = 0; i < nbasis; ++i)
    {
      for (unsigned j = 0; j < nbasic; ++j)
      {
        Direct_basic_association_matrix(i, j) =
          Basic_dofs_from_global_dofs(j, i);
      }
    }

    Algebraic_transformation_is_built = true;
  }


  // Explicit instantiation for the two supported boundary orders.
  template class DirectBernadouElementBasis<3>;
  template class DirectBernadouElementBasis<5>;

} // namespace oomph

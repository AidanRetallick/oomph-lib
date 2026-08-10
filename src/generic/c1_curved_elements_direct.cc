// ====================================================================
// Algebraic direct construction of the Bernadou curved Bell basis.
// ====================================================================

#include "c1_curved_elements_direct.h"
#include <cmath>

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


  /// Build T and C if required.
  template<unsigned BOUNDARY_ORDER>
  void DirectBernadouElementBasis<
    BOUNDARY_ORDER>::ensure_algebraic_transformation() const
  {
    if (!Algebraic_transformation_is_built)
    {
      build_algebraic_transformation();
    }
  }


  /// Return reference location, reference tangent, outward unit reference
  /// normal, and endpoint vertices for one of the three basic edges.
  ///
  /// The edge parameters are:
  ///   edge 0: s=(0,u),     node 2 -> node 1;
  ///   edge 1: s=(u,0),     node 2 -> node 0;
  ///   edge 2: s=(u,1-u),   node 1 -> node 0.
  template<unsigned BOUNDARY_ORDER>
  void DirectBernadouElementBasis<BOUNDARY_ORDER>::get_reference_edge_geometry(
    const unsigned& edge,
    const double& u,
    Vector<double>& s,
    Vector<double>& tangent,
    Vector<double>& outward_normal,
    unsigned& start_vertex,
    unsigned& end_vertex)
  {
    s.resize(2);
    tangent.resize(2);
    outward_normal.resize(2);

    switch (edge)
    {
      // s0=0, parameter u=s1, node 2 -> node 1.
      case 0:
        s[0] = 0.0;
        s[1] = u;
        tangent[0] = 0.0;
        tangent[1] = 1.0;
        outward_normal[0] = -1.0;
        outward_normal[1] = 0.0;
        start_vertex = 2;
        end_vertex = 1;
        break;

      // s1=0, parameter u=s0, node 2 -> node 0.
      case 1:
        s[0] = u;
        s[1] = 0.0;
        tangent[0] = 1.0;
        tangent[1] = 0.0;
        outward_normal[0] = 0.0;
        outward_normal[1] = -1.0;
        start_vertex = 2;
        end_vertex = 0;
        break;

      // s0+s1=1, parameter u=s0, node 1 -> node 0.
      case 2:
        s[0] = u;
        s[1] = 1.0 - u;
        tangent[0] = 1.0;
        tangent[1] = -1.0;
        outward_normal[0] = 1.0 / std::sqrt(2.0);
        outward_normal[1] = 1.0 / std::sqrt(2.0);
        start_vertex = 1;
        end_vertex = 0;
        break;

      default:
        throw OomphLibError("Reference edge index must be 0, 1 or 2.",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
    }
  }


  /// Basic-edge interpolation abscissa in the ordering used by the full
  /// basic basis (not including midpoint).
  template<unsigned BOUNDARY_ORDER>
  double DirectBernadouElementBasis<BOUNDARY_ORDER>::edge_dof_abscissa(
    const unsigned& edge, const unsigned& point_index)
  {
    const unsigned npoint = BOUNDARY_ORDER - 1;
    if (point_index >= npoint)
    {
      throw OomphLibError("Basic-edge dof point index is out of range.",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

    double ascending_point = 0.0;

    if (BOUNDARY_ORDER == 3)
    {
      const double point[2] = {1.0 / 4.0, 3.0 / 4.0};
      ascending_point = point[point_index];
    }
    else
    {
      const double point[4] = {1.0 / 6.0, 2.0 / 6.0, 4.0 / 6.0, 5.0 / 6.0};
      ascending_point = point[point_index];
    }

    // The basic-dof ordering is descending on edges 0 and 2, and ascending
    // on edge 1.
    if (edge == 1)
    {
      return ascending_point;
    }

    if (edge == 0 || edge == 2)
    {
      return 1.0 - ascending_point;
    }

    throw OomphLibError("Reference edge index must be 0, 1 or 2.",
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
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


  /// Derivative with respect to u of the quintic Hermite basis.
  template<unsigned BOUNDARY_ORDER>
  void DirectBernadouElementBasis<BOUNDARY_ORDER>::d_quintic_hermite(
    const double& u, Vector<double>& dh)
  {
    dh.resize(6);
    dh[0] = -30.0 * u * u * std::pow(1.0 - u, 2);
    dh[1] = 30.0 * u * u * std::pow(1.0 - u, 2);
    dh[2] = -std::pow(1.0 - u, 2) * (3.0 * u - 1.0) * (5.0 * u + 1.0);
    dh[3] = -u * u * (3.0 * u - 2.0) * (5.0 * u - 6.0);
    dh[4] = -0.5 * u * std::pow(1.0 - u, 2) * (5.0 * u - 2.0);
    dh[5] = 0.5 * u * u * (u - 1.0) * (5.0 * u - 3.0);
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


  /// destination += scale * source.
  template<unsigned BOUNDARY_ORDER>
  void DirectBernadouElementBasis<BOUNDARY_ORDER>::add_scaled(
    Vector<double>& destination,
    const Vector<double>& source,
    const double& scale)
  {
    const unsigned n = destination.size();
    for (unsigned i = 0; i < n; ++i)
    {
      destination[i] += scale * source[i];
    }
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
  /// a.grad_x(w) at a vertex.
  template<unsigned BOUNDARY_ORDER>
  void DirectBernadouElementBasis<BOUNDARY_ORDER>::
    vertex_physical_first_derivative_coefficients(
      const unsigned& vertex,
      const Vector<double>& a,
      Vector<double>& coefficients) const
  {
    zero_global_dof_coefficients(coefficients);
    coefficients[3 + 2 * vertex] = a[0];
    coefficients[4 + 2 * vertex] = a[1];
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


  /// Coefficients of a.grad_s(w) at a vertex.
  ///   a.grad_s(w) = a.J.grad_x(w) = (a.J).grad_x(w)
  template<unsigned BOUNDARY_ORDER>
  void DirectBernadouElementBasis<BOUNDARY_ORDER>::
    vertex_reference_first_derivative_coefficients(
      const unsigned& vertex,
      const Vector<double>& a,
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
      physical_direction[i] = jacobian(i, 0) * a[0] + jacobian(i, 1) * a[1];
    }

    vertex_physical_first_derivative_coefficients(
      vertex, physical_direction, coefficients);
  }


  /// Coefficients of a^T H_s(w) b at a vertex, including the Hessian of the
  /// coordinate map.
  ///   a^T H_s(w) b = (J a)^T H_x(w) (J b)
  ///                    + grad_x(w) . H_x[a,b].
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
    //                    + grad_x(w) . H_F[a,b]. <- This term
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


  /// Algebraic coefficients of the quintic value trace at an edge point. This
  /// explicitly imposes the constraint that edge values are interpolated
  /// quintically.
  template<unsigned BOUNDARY_ORDER>
  void DirectBernadouElementBasis<BOUNDARY_ORDER>::edge_value_coefficients(
    const unsigned& edge, const double& u, Vector<double>& coefficients) const
  {
    Vector<double> s(2), tangent(2), normal(2);
    unsigned start_vertex = 0;
    unsigned end_vertex = 0;
    get_reference_edge_geometry(edge, u, s, tangent, normal, start_vertex, end_vertex);

    Vector<double> h(6);
    quintic_hermite(u, h);

    Vector<double> endpoint_functional[6];
    vertex_value_coefficients(start_vertex, endpoint_functional[0]);
    vertex_value_coefficients(end_vertex, endpoint_functional[1]);
    vertex_reference_first_derivative_coefficients(
      start_vertex, tangent, endpoint_functional[2]);
    vertex_reference_first_derivative_coefficients(
      end_vertex, tangent, endpoint_functional[3]);
    vertex_reference_second_derivative_coefficients(
      start_vertex, tangent, tangent, endpoint_functional[4]);
    vertex_reference_second_derivative_coefficients(
      end_vertex, tangent, tangent, endpoint_functional[5]);

    zero_global_dof_coefficients(coefficients);
    for (unsigned i = 0; i < 6; ++i)
    {
      add_scaled(coefficients, endpoint_functional[i], h[i]);
    }
  }


  /// Algebraic coefficients of the derivative with respect to the chosen
  /// edge parameter of the quintic value trace.
  /// (Needed for the normal derivative interpolation on straight edges)
  template<unsigned BOUNDARY_ORDER>
  void DirectBernadouElementBasis<BOUNDARY_ORDER>::
    edge_value_parameter_derivative_coefficients(
      const unsigned& edge, const double& u, Vector<double>& coefficients) const
  {
    Vector<double> s(2), tangent(2), normal(2);
    unsigned start_vertex = 0;
    unsigned end_vertex = 0;
    get_reference_edge_geometry(
      edge, u, s, tangent, normal, start_vertex, end_vertex);

    Vector<double> dh(6);
    d_quintic_hermite(u, dh);

    Vector<double> endpoint_functional[6];
    vertex_value_coefficients(start_vertex, endpoint_functional[0]);
    vertex_value_coefficients(end_vertex, endpoint_functional[1]);
    vertex_reference_first_derivative_coefficients(
      start_vertex, tangent, endpoint_functional[2]);
    vertex_reference_first_derivative_coefficients(
      end_vertex, tangent, endpoint_functional[3]);
    vertex_reference_second_derivative_coefficients(
      start_vertex, tangent, tangent, endpoint_functional[4]);
    vertex_reference_second_derivative_coefficients(
      end_vertex, tangent, tangent, endpoint_functional[5]);

    zero_global_dof_coefficients(coefficients);
    for (unsigned i = 0; i < 6; ++i)
    {
      add_scaled(coefficients, endpoint_functional[i], dh[i]);
    }
  }


  /// Algebraic coefficients of the basic outward-normal derivative at an edge
  /// point. This explicitly imposes that the physical normal on straight edges
  /// and the reference normal on the curved edge are cubically interpolated by
  /// the nodal data.
  template<unsigned BOUNDARY_ORDER>
  void DirectBernadouElementBasis<BOUNDARY_ORDER>::
    edge_basic_normal_derivative_coefficients(
      const unsigned& edge, const double& u, Vector<double>& coefficients) const
  {
    Vector<double> s(2);
    Vector<double> reference_tangent(2);
    Vector<double> reference_normal(2);

    unsigned start_vertex = 0;
    unsigned end_vertex = 0;

    get_reference_edge_geometry(edge,
				u,
				s,
				reference_tangent,
				reference_normal,
				start_vertex,
				end_vertex);

    Vector<double> h3(4);
    cubic_hermite(u, h3);

    // The distinguished curved edge interpolates the reference normal:
    //     g = n . grad_s(w),       n = (1,1)/sqrt(2),
    // cubically.
    if (edge == 2)
    {
      Vector<double> g_functional[4];

      vertex_reference_first_derivative_coefficients(
        start_vertex, reference_normal, g_functional[0]);

      vertex_reference_first_derivative_coefficients(
        end_vertex, reference_normal, g_functional[1]);

      vertex_reference_second_derivative_coefficients(
        start_vertex, reference_normal, reference_tangent, g_functional[2]);

      vertex_reference_second_derivative_coefficients(
        end_vertex, reference_normal, reference_tangent, g_functional[3]);

      zero_global_dof_coefficients(coefficients);

      for (unsigned i = 0; i < 4; ++i)
      {
        add_scaled(coefficients, g_functional[i], h3[i]);
      }

      return;
    }

    // Construct the basic outward-normal derivative DOF on a straight edge.
    //
    // The physical edge normal is constant, so the normal derivative
    //
    //     g(u) = n . grad_x(w)
    //
    // is represented by a cubic Hermite interpolant along the edge. Its four
    // Hermite data are obtained from the physical vertex DOFs:
    //
    //     g(0), g(1)
    //       from the vertex gradients,
    //     g'(0), g'(1)
    //       from n^T H_x(w) t,
    //
    // where t = dx/du is the physical edge tangent. Since the edge mapping is
    // affine, both n and t are constant along the edge.
    //
    // The resulting cubic gives the desired physical normal derivative at the
    // requested edge coordinate u. This is then converted to the derivative in
    // the basic outward-normal direction by decomposing the mapped basic normal
    // J n_hat into its physical tangent and normal components:
    //
    //     J n_hat = alpha t + beta n.
    //
    // Hence
    //
    //     n_hat . grad_s(w)
    //       = alpha (t . grad_x(w)) + beta (n . grad_x(w)),
    //
    // with the tangential contribution obtained from the derivative of the
    // quintic edge-value trace and the normal contribution from the cubic
    // Hermite interpolant above. The returned coefficient vector therefore
    // expresses the basic normal-derivative DOF directly in terms of the global
    // physical Hermite DOFs.

    // Get the Jacobian at the requested edge point.
    DenseMatrix<double> jacobian(2, 2, 0.0);
    this->get_basic_jacobian(s, jacobian);

    // Map the basic edge tangent into physical space.
    Vector<double> physical_tangent(2, 0.0);
    for (unsigned i = 0; i < 2; ++i)
    {
      physical_tangent[i] = jacobian(i, 0) * reference_tangent[0] +
			    jacobian(i, 1) * reference_tangent[1];
    }

    // Use its 90-degree rotation as the fixed physical edge normal.
    Vector<double> physical_normal(2, 0.0);
    physical_normal[0] = -physical_tangent[1];
    physical_normal[1] = physical_tangent[0];

#ifdef PARANOID
    // Check that the edge isn't collapsed compared to the element size
    double jacobian_scale_squared = 0.0;
    for (unsigned i = 0; i < 2; ++i)
    {
      for (unsigned j = 0; j < 2; ++j)
      {
        jacobian_scale_squared += jacobian(i, j) * jacobian(i, j);
      }
    }

    const double tangent_norm_squared =
      physical_tangent[0] * physical_tangent[0] +
      physical_tangent[1] * physical_tangent[1];

    const double tolerance = 100.0 * std::numeric_limits<double>::epsilon();

    if (tangent_norm_squared <= tolerance * jacobian_scale_squared)
    {
      throw OomphLibError("Collapsed or nearly collapsed straight edge in "
                          "curved C1 basis construction.",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Construct the four Hermite data for
    //     g(u) = N . grad_x(w),
    // where N is the fixed physical edge normal.
    //
    // Since N and the physical tangent t are constant,
    //     g'(u) = N^T H_x(w) t.
    Vector<double> g_functional[4];

    vertex_physical_first_derivative_coefficients(
      start_vertex, physical_normal, g_functional[0]);

    vertex_physical_first_derivative_coefficients(
      end_vertex, physical_normal, g_functional[1]);

    vertex_physical_second_derivative_coefficients(
      start_vertex, physical_normal, physical_tangent, g_functional[2]);

    vertex_physical_second_derivative_coefficients(
      end_vertex, physical_normal, physical_tangent, g_functional[3]);

    // Evaluate the cubic physical-normal derivative trace at u.
    Vector<double> physical_normal_trace;
    zero_global_dof_coefficients(physical_normal_trace);
    for (unsigned i = 0; i < 4; ++i)
    {
      add_scaled(physical_normal_trace, g_functional[i], h3[i]);
    }

    // Obtain the physical image of the basic outward-normal direction at the
    // requested point.
    //
    // Unlike the physical edge tangent, J*n_hat need not be constant along the
    // straight edge, so this must use the Jacobian at the actual point u.
    Vector<double> mapped_basic_normal(2, 0.0);
    for (unsigned i = 0; i < 2; ++i)
    {
      mapped_basic_normal[i] = jacobian(i, 0) * reference_normal[0] +
			       jacobian(i, 1) * reference_normal[1];
    }

    // The tangential derivative of the quintic edge-value trace.
    Vector<double> value_trace_derivative;
    edge_value_parameter_derivative_coefficients(
      edge, u, value_trace_derivative);

    // Decompose the mapped basic normal into the physical tangent/normal basis:
    //
    //     J*n_hat = alpha*t + beta*N.
    const double determinant = physical_tangent[0] * physical_normal[1] -
                               physical_tangent[1] * physical_normal[0];

    const double alpha = (mapped_basic_normal[0] * physical_normal[1] -
                          mapped_basic_normal[1] * physical_normal[0]) /
                         determinant;

    const double beta = (physical_tangent[0] * mapped_basic_normal[1] -
                         physical_tangent[1] * mapped_basic_normal[0]) /
                        determinant;

    // Hence
    //
    //     n_hat . grad_s(w)
    //       = alpha * (t . grad_x(w))
    //       + beta  * (N . grad_x(w)).
    //
    zero_global_dof_coefficients(coefficients);
    add_scaled(coefficients, value_trace_derivative, alpha);
    add_scaled(coefficients, physical_normal_trace, beta);
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
      edge_basic_normal_derivative_coefficients(edge, 0.5, coefficients);
      put_row(Basic_dofs_from_global_dofs, 18 + edge, coefficients);
    }

    // ------------------------------------------------------------------
    // 4. Additional basic edge values.
    // ------------------------------------------------------------------
    for (unsigned edge = 0; edge < 3; ++edge)
    {
      for (unsigned point = 0; point < npoint_per_edge; ++point)
      {
        const double u = edge_dof_abscissa(edge, point);
        edge_value_coefficients(edge, u, coefficients);
        put_row(Basic_dofs_from_global_dofs,
                21 + edge * npoint_per_edge + point,
                coefficients);
      }
    }

    // ------------------------------------------------------------------
    // 5. Additional basic edge outward-normal derivatives.
    // ------------------------------------------------------------------
    for (unsigned edge = 0; edge < 3; ++edge)
    {
      for (unsigned point = 0; point < npoint_per_edge; ++point)
      {
        const double u = edge_dof_abscissa(edge, point);
        edge_basic_normal_derivative_coefficients(edge, u, coefficients);
        put_row(Basic_dofs_from_global_dofs,
                21 + nmidnodes + edge * npoint_per_edge + point,
                coefficients);
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


  /// Fill T, where
  ///
  ///          basic_dofs = T * global_dofs.
  ///
  /// T has dimensions n_basic_basis_functions by n_basis_functions.
  template<unsigned BOUNDARY_ORDER>
  void DirectBernadouElementBasis<BOUNDARY_ORDER>::
    fill_in_basic_dof_transformation_matrix(
      DenseMatrix<double>& transformation_matrix) const
  {
    ensure_algebraic_transformation();
    transformation_matrix = Basic_dofs_from_global_dofs;
  }


  /// Fill C=T^T, where
  ///
  ///          global_basis = C * full_basic_basis.
  ///
  /// C has dimensions n_basis_functions by n_basic_basis_functions.
  template<unsigned BOUNDARY_ORDER>
  void DirectBernadouElementBasis<BOUNDARY_ORDER>::
    fill_in_direct_basic_association_matrix(
      DenseMatrix<double>& direct_matrix) const
  {
    ensure_algebraic_transformation();
    direct_matrix = Direct_basic_association_matrix;
  }


  /// Fill the global-basis-to-monomial association matrix.
  /// [zdec] We only need this for comparing to original elements
  template<unsigned BOUNDARY_ORDER>
  void DirectBernadouElementBasis<
    BOUNDARY_ORDER>::fill_in_full_association_matrix(DenseMatrix<double>&
                                                       conversion_matrix) const
  {
    ensure_algebraic_transformation();

    const unsigned nbasis = this->n_basis_functions();
    const unsigned nbasic = this->n_basic_basis_functions();

    DenseMatrix<double> monomial_to_basic(nbasic, nbasic, 0.0);
    this->monomial_to_basic_matrix(monomial_to_basic);

    conversion_matrix.resize(nbasis, nbasic, 0.0);
    for (unsigned i = 0; i < nbasis; ++i)
    {
      for (unsigned k = 0; k < nbasic; ++k)
      {
        const double cik = Direct_basic_association_matrix(i, k);
	// [zdec] HMMMMM! Wrote this at first because some entries should really
	// be zero and untouched... Does this need fp tolerance?
        if (cik == 0.0)
        {
          continue;
        }
        for (unsigned j = 0; j < nbasic; ++j)
        {
          conversion_matrix(i, j) += cik * monomial_to_basic(k, j);
        }
      }
    }
  }


  /// Development/validation helper. The returned matrices are in the
  /// monomial basis. direct_matrix and factorised_matrix should agree up to
  /// roundoff when the algebraic direct construction reproduces the original
  /// D*B construction.
  template<unsigned BOUNDARY_ORDER>
  void DirectBernadouElementBasis<BOUNDARY_ORDER>::
    compare_with_factorised_association_matrix(
      DenseMatrix<double>& direct_matrix,
      DenseMatrix<double>& factorised_matrix,
      DenseMatrix<double>& difference_matrix) const
  {
    const unsigned nbasis = this->n_basis_functions();
    const unsigned nbasic = this->n_basic_basis_functions();

    fill_in_full_association_matrix(direct_matrix);

    // The original routine expects a pre-sized, zero-initialised output.
    factorised_matrix.resize(nbasis, nbasic, 0.0);
    Base::fill_in_full_association_matrix(factorised_matrix);

    difference_matrix.resize(nbasis, nbasic, 0.0);
    for (unsigned i = 0; i < nbasis; ++i)
    {
      for (unsigned j = 0; j < nbasic; ++j)
      {
        difference_matrix(i, j) = direct_matrix(i, j) - factorised_matrix(i, j);
      }
    }
  }


  /// Transform reference derivatives of an arbitrary list of functions to
  /// Eulerian derivatives.
  template<unsigned BOUNDARY_ORDER>
  double DirectBernadouElementBasis<BOUNDARY_ORDER>::
    transform_basis_derivatives_reference_to_eulerian(
      const Vector<double>& s,
      const unsigned& nfunction,
      const DShape& d_reference,
      const DShape& d2_reference,
      DShape& d_eulerian,
      DShape& d2_eulerian,
      const bool& compute_second_derivatives) const
  {
    DenseMatrix<double> jacobian(2, 2, 0.0);
    DenseMatrix<double> inverse_jacobian(2, 2, 0.0);
    this->get_basic_jacobian(s, jacobian);

    const double det =
      jacobian(0, 0) * jacobian(1, 1) - jacobian(0, 1) * jacobian(1, 0);

#ifdef PARANOID
    // Check that the jacobian isn't degenerate compared to the Frobenius norm
    double jacobian_norm_squared = 0.0;

    for (unsigned i = 0; i < 2; ++i)
    {
      for (unsigned j = 0; j < 2; ++j)
      {
        jacobian_norm_squared += jacobian(i, j) * jacobian(i, j);
      }
    }

    if ((jacobian_norm_squared == 0.0) ||
        (std::fabs(det) <= 100.0 * std::numeric_limits<double>::epsilon() *
                             jacobian_norm_squared))
    {
      throw OomphLibError(
        "Singular or nearly singular curved-element Jacobian. det(J) = " +
          std::to_string(det) + ".",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // inverse_jacobian(alpha,i) = d s_alpha / d x_i.
    inverse_jacobian(0, 0) = jacobian(1, 1) / det;
    inverse_jacobian(0, 1) = -jacobian(0, 1) / det;
    inverse_jacobian(1, 0) = -jacobian(1, 0) / det;
    inverse_jacobian(1, 1) = jacobian(0, 0) / det;

    RankThreeTensor<double> map_hessian(2, 2, 2, 0.0);
    if (compute_second_derivatives)
    {
      this->get_basic_hessian(s, map_hessian);
    }

    for (unsigned j = 0; j < nfunction; ++j)
    {
      const double gx = d_reference(j, 0) * inverse_jacobian(0, 0) +
                        d_reference(j, 1) * inverse_jacobian(1, 0);
      const double gy = d_reference(j, 0) * inverse_jacobian(0, 1) +
                        d_reference(j, 1) * inverse_jacobian(1, 1);

      d_eulerian(j, 0) = gx;
      d_eulerian(j, 1) = gy;

      if (compute_second_derivatives)
      {
        // H_ref = J^T H_x J + sum_i w_,i H(F_i).
        // Therefore K = H_ref - sum_i w_,i H(F_i), followed by
        // H_x = J^{-T} K J^{-1}.
        const double k00 = d2_reference(j, 0) - gx * map_hessian(0, 0, 0) -
                           gy * map_hessian(1, 0, 0);
        const double k01 = d2_reference(j, 1) - gx * map_hessian(0, 0, 1) -
                           gy * map_hessian(1, 0, 1);
        const double k11 = d2_reference(j, 2) - gx * map_hessian(0, 1, 1) -
                           gy * map_hessian(1, 1, 1);

        const double i00 = inverse_jacobian(0, 0);
        const double i01 = inverse_jacobian(0, 1);
        const double i10 = inverse_jacobian(1, 0);
        const double i11 = inverse_jacobian(1, 1);

        d2_eulerian(j, 0) =
          i00 * i00 * k00 + 2.0 * i00 * i10 * k01 + i10 * i10 * k11;

        d2_eulerian(j, 1) =
          i00 * i01 * k00 + (i00 * i11 + i10 * i01) * k01 + i10 * i11 * k11;

        d2_eulerian(j, 2) =
          i01 * i01 * k00 + 2.0 * i01 * i11 * k01 + i11 * i11 * k11;
      }
    }

    return det;
  }


  // [zdec] better name for this
  /// Evaluate the final global basis and its reference derivatives from C.
  template<unsigned BOUNDARY_ORDER>
  void DirectBernadouElementBasis<BOUNDARY_ORDER>::
    evaluate_global_basis_wrt_reference(const Vector<double>& s,
                                   Shape& value,
                                   DShape& dvalue,
                                   DShape& d2value,
                                   const bool& compute_second_derivatives) const
  {
    ensure_algebraic_transformation();

    const unsigned nbasis = this->n_basis_functions();
    const unsigned nbasic = this->n_basic_basis_functions();

    Shape basic_value(nbasic);
    DShape basic_dvalue(nbasic, 2);
    DShape basic_d2value(nbasic, 3);

    this->full_basic_polynomials(s, basic_value);
    this->dfull_basic_polynomials(s, basic_dvalue);
    if (compute_second_derivatives)
    {
      this->d2full_basic_polynomials(s, basic_d2value);
    }

    for (unsigned i = 0; i < nbasis; ++i)
    {
      value[i] = 0.0;
      dvalue(i, 0) = 0.0;
      dvalue(i, 1) = 0.0;
      if (compute_second_derivatives)
      {
        d2value(i, 0) = 0.0;
        d2value(i, 1) = 0.0;
        d2value(i, 2) = 0.0;
      }
      for (unsigned j = 0; j < nbasic; ++j)
      {
        const double coefficient = Direct_basic_association_matrix(i, j);
        value[i] += coefficient * basic_value[j];
        dvalue(i, 0) += coefficient * basic_dvalue(j, 0);
        dvalue(i, 1) += coefficient * basic_dvalue(j, 1);
        if (compute_second_derivatives)
        {
          d2value(i, 0) += coefficient * basic_d2value(j, 0);
          d2value(i, 1) += coefficient * basic_d2value(j, 1);
          d2value(i, 2) += coefficient * basic_d2value(j, 2);
        }
      }
    }
  }


  /// Scatter a flat list of global basis values into oomph-lib's nodal and
  /// bubble Shape containers, applying the curved-edge node permutation.
  template<unsigned BOUNDARY_ORDER>
  void DirectBernadouElementBasis<BOUNDARY_ORDER>::scatter_values(
    const Shape& value, Shape& psi, Shape& bpsi) const
  {
    unsigned index_shift = 0;
    this->nodal_index_shift(index_shift);

    for (unsigned inode = 0; inode < 3; ++inode)
    {
      const unsigned global_node = (inode + index_shift) % 3;
      psi(global_node, 0) = value[inode];
      psi(global_node, 1) = value[3 + 2 * inode];
      psi(global_node, 2) = value[4 + 2 * inode];
      psi(global_node, 3) = value[9 + 3 * inode];
      psi(global_node, 4) = value[10 + 3 * inode];
      psi(global_node, 5) = value[11 + 3 * inode];
    }

    for (unsigned i = 0; i < this->n_internal_dofs(); ++i)
    {
      bpsi(i, 0) = value[18 + i];
    }
  }


  /// Scatter first derivatives.
  template<unsigned BOUNDARY_ORDER>
  void DirectBernadouElementBasis<BOUNDARY_ORDER>::scatter_first_derivatives(
    const DShape& dvalue, DShape& dpsi, DShape& dbpsi) const
  {
    unsigned index_shift = 0;
    this->nodal_index_shift(index_shift);

    for (unsigned inode = 0; inode < 3; ++inode)
    {
      const unsigned global_node = (inode + index_shift) % 3;
      const unsigned flat_index[6] = {inode,
                                      3 + 2 * inode,
                                      4 + 2 * inode,
                                      9 + 3 * inode,
                                      10 + 3 * inode,
                                      11 + 3 * inode};

      for (unsigned itype = 0; itype < 6; ++itype)
      {
        dpsi(global_node, itype, 0) = dvalue(flat_index[itype], 0);
        dpsi(global_node, itype, 1) = dvalue(flat_index[itype], 1);
      }
    }

    for (unsigned i = 0; i < this->n_internal_dofs(); ++i)
    {
      dbpsi(i, 0, 0) = dvalue(18 + i, 0);
      dbpsi(i, 0, 1) = dvalue(18 + i, 1);
    }
  }


  /// Scatter second derivatives.
  template<unsigned BOUNDARY_ORDER>
  void DirectBernadouElementBasis<BOUNDARY_ORDER>::scatter_second_derivatives(
    const DShape& d2value, DShape& d2psi, DShape& d2bpsi) const
  {
    unsigned index_shift = 0;
    this->nodal_index_shift(index_shift);

    for (unsigned inode = 0; inode < 3; ++inode)
    {
      const unsigned global_node = (inode + index_shift) % 3;
      const unsigned flat_index[6] = {inode,
                                      3 + 2 * inode,
                                      4 + 2 * inode,
                                      9 + 3 * inode,
                                      10 + 3 * inode,
                                      11 + 3 * inode};

      for (unsigned itype = 0; itype < 6; ++itype)
      {
        d2psi(global_node, itype, 0) = d2value(flat_index[itype], 0);
        d2psi(global_node, itype, 1) = d2value(flat_index[itype], 1);
        d2psi(global_node, itype, 2) = d2value(flat_index[itype], 2);
      }
    }

    for (unsigned i = 0; i < this->n_internal_dofs(); ++i)
    {
      d2bpsi(i, 0, 0) = d2value(18 + i, 0);
      d2bpsi(i, 0, 1) = d2value(18 + i, 1);
      d2bpsi(i, 0, 2) = d2value(18 + i, 2);
    }
  }


  /// Return basis values.
  template<unsigned BOUNDARY_ORDER>
  void DirectBernadouElementBasis<BOUNDARY_ORDER>::shape(
    const Vector<double>& s, Shape& psi, Shape& bpsi) const
  {
    Vector<double> s_basic(s);
    this->permute_shape(s_basic);

    Shape value(this->n_basis_functions());
    DShape dvalue(this->n_basis_functions(), 2);
    DShape d2value(this->n_basis_functions(), 3);
    evaluate_global_basis_wrt_reference(s_basic, value, dvalue, d2value, false);
    scatter_values(value, psi, bpsi);
  }


  /// Return basis values and Eulerian first derivatives.
  template<unsigned BOUNDARY_ORDER>
  double DirectBernadouElementBasis<BOUNDARY_ORDER>::d_shape_dx(
    const Vector<double>& s,
    Shape& psi,
    Shape& bpsi,
    DShape& dpsi,
    DShape& dbpsi) const
  {
    Vector<double> s_basic(s);
    this->permute_shape(s_basic);

    const unsigned nbasis = this->n_basis_functions();
    Shape value(nbasis);
    DShape d_reference(nbasis, 2);
    DShape dummy_d2_reference(nbasis, 3);
    evaluate_global_basis_wrt_reference(
      s_basic, value, d_reference, dummy_d2_reference, false);

    DShape d_eulerian(nbasis, 2);
    DShape dummy_d2_eulerian(nbasis, 3);
    const double det = transform_basis_derivatives_reference_to_eulerian(s_basic,
                                                       nbasis,
                                                       d_reference,
                                                       dummy_d2_reference,
                                                       d_eulerian,
                                                       dummy_d2_eulerian,
                                                       false);

    scatter_values(value, psi, bpsi);
    scatter_first_derivatives(d_eulerian, dpsi, dbpsi);
    return det;
  }


  /// Return basis values and Eulerian first and second derivatives.
  template<unsigned BOUNDARY_ORDER>
  double DirectBernadouElementBasis<BOUNDARY_ORDER>::d2_shape_dx2(
    const Vector<double>& s,
    Shape& psi,
    Shape& bpsi,
    DShape& dpsi,
    DShape& dbpsi,
    DShape& d2psi,
    DShape& d2bpsi) const
  {
    Vector<double> s_basic(s);
    this->permute_shape(s_basic);

    const unsigned nbasis = this->n_basis_functions();
    Shape value(nbasis);
    DShape d_reference(nbasis, 2);
    DShape d2_reference(nbasis, 3);
    evaluate_global_basis_wrt_reference(
      s_basic, value, d_reference, d2_reference, true);

    DShape d_eulerian(nbasis, 2);
    DShape d2_eulerian(nbasis, 3);
    const double det = transform_basis_derivatives_reference_to_eulerian(s_basic,
                                                       nbasis,
                                                       d_reference,
                                                       d2_reference,
                                                       d_eulerian,
                                                       d2_eulerian,
                                                       true);

    scatter_values(value, psi, bpsi);
    scatter_first_derivatives(d_eulerian, dpsi, dbpsi);
    scatter_second_derivatives(d2_eulerian, d2psi, d2bpsi);
    return det;
  }


  // Explicit instantiation for the two supported boundary orders.
  template class DirectBernadouElementBasis<3>;
  template class DirectBernadouElementBasis<5>;

} // namespace oomph

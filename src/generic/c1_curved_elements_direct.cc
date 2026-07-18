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
      // [zdec] Method that fills out identity block for vertices
    }

    // ------------------------------------------------------------------
    // 2. Basic first and second reference derivatives at the vertices.
    // ------------------------------------------------------------------
    Vector<double> e0(2, 0.0), e1(2, 0.0);
    e0[0] = 1.0;
    e1[1] = 1.0;

    for (unsigned vertex = 0; vertex < 3; ++vertex)
    {
      // [zdec] Method that fills out the first derivative blocks
      // for rows 
      //   3 + 2 * vertex
      //   4 + 2 * vertex

      // [zdec] Method that fills out the second derivative blocks
      // for rows 
      //   9 + 3 * vertex
      //  10 + 3 * vertex
      //  11 + 3 * vertex
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
      // [zdec] Method that fills out identity block for internal dofs
    }

#ifdef PARANOID
    if (21 + 2 * nmidnodes + ninternal != nbasic)
    {
      throw OomphLibError("Internal error in algebraic direct dof count.",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
#endif

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

// ====================================================================
// Algebraic direct construction of the Bernadou curved Bell basis.
// ====================================================================

#ifndef OOMPH_DIRECT_BERNADOU_ELEMENT_BASIS_H
#define OOMPH_DIRECT_BERNADOU_ELEMENT_BASIS_H

#include "c1_curved_elements.h"

namespace oomph
{

  /// Algebraic direct construction of the P3/P5 Bernadou curved Bell basis.
  ///
  /// Let d_g denote the vector of physical/global degrees of freedom and d_b
  /// the full vector of basic-triangle degrees of freedom. Similarly, let psi_g
  /// denote the vector of physical/global basis functions and psi_b the vector
  /// of basic-triangle degrees of freedom. The trace restrictions and the chain
  /// rule give the rectangular transformation
  ///
  ///                 d_b = T d_g
  ///
  /// directly and algebraically. In order to preserve the interpolation
  ///
  ///                 w = d_b^T psi_b = d_g^T psi_g,
  ///
  /// the global basis can be written in terms of the basic basis
  ///
  ///                 psi_g = C psi_b
  ///                 C = T^T.
  ///
  /// The entries of T are assembled from:
  ///  * the Jacobian and Hessian of the curved coordinate map at the vertices;
  ///  * quintic Hermite interpolation of the value trace;
  ///  * cubic Hermite interpolation of the physical normal derivative on the
  ///    two straight edges;
  ///  * cubic Hermite interpolation of the existing fixed reference-
  ///    transverse derivative on the distinguished curved edge;
  ///  * identity maps for the internal value degrees of freedom.
  template<unsigned BOUNDARY_ORDER>
  class DirectBernadouElementBasis : public BernadouElementBasis<BOUNDARY_ORDER>
  {
  public:
    typedef BernadouElementBasis<BOUNDARY_ORDER> Base;
    typedef typename Base::VertexList VertexList;

    static_assert(BOUNDARY_ORDER == 3 || BOUNDARY_ORDER == 5,
                  "DirectBernadouElementBasis is implemented only for "
                  "boundary orders 3 and 5.");

    /// Constructor. The algebraic transformation is built after
    /// upgrade_element().
    DirectBernadouElementBasis();

    /// Empty destructor.
    virtual ~DirectBernadouElementBasis() {}

    /// Upgrade the geometry and construct T and C=T^T algebraically.
    void upgrade_element(
      const VertexList& verts,
      const double& su,
      const double& so,
      const C1PlateHelper::CurvedEdgeEnumeration& curved_edge,
      const C1CurviLine& parametric_curve) override;

    /// Return basis values.
    void shape(const Vector<double>& s,
               Shape& nodal_basis,
               Shape& bubble_basis) const override;

    /// Return basis values and Eulerian first derivatives.
    double d_shape_dx(const Vector<double>& s,
                      Shape& psi,
                      Shape& bpsi,
                      DShape& dpsi,
                      DShape& dbpsi) const override;

    /// Return basis values and Eulerian first and second derivatives.
    double d2_shape_dx2(const Vector<double>& s,
                        Shape& psi,
                        Shape& bpsi,
                        DShape& dpsi,
                        DShape& dbpsi,
                        DShape& d2psi,
                        DShape& d2bpsi) const override;

    /// Fill the global-basis-to-monomial association matrix.
    /// [zdec] We only need this for comparing to original elements
    void fill_in_full_association_matrix(
      DenseMatrix<double>& conversion_matrix) const override;

    /// Fill T, where
    ///
    ///          basic_dofs = T * global_dofs.
    ///
    /// T has dimensions n_basic_basis_functions by n_basis_functions.
    /// Exists as a public read interface.
    void fill_in_basic_dof_transformation_matrix(
      DenseMatrix<double>& transformation_matrix) const;

    /// Fill C=T^T, where
    ///
    ///          global_basis = C * full_basic_basis.
    ///
    /// C has dimensions n_basis_functions by n_basic_basis_functions.
    /// Exists as a public read interface.
    void fill_in_direct_basic_association_matrix(
      DenseMatrix<double>& direct_matrix) const;

    /// Development/validation helper. The returned matrices are in the
    /// monomial basis. direct_matrix and factorised_matrix should agree up to
    /// roundoff when the algebraic direct construction reproduces the original
    /// D*B construction.
    void compare_with_factorised_association_matrix(
      DenseMatrix<double>& direct_matrix,
      DenseMatrix<double>& factorised_matrix,
      DenseMatrix<double>& difference_matrix) const;

  private:
    /// No copying.
    DirectBernadouElementBasis(const DirectBernadouElementBasis&);
    void operator=(const DirectBernadouElementBasis&);

    /// T: full basic dofs as algebraic functions of the global dofs.
    mutable DenseMatrix<double> Basic_dofs_from_global_dofs;

    /// C=T^T: global basis functions in the full basic basis.
    mutable DenseMatrix<double> Direct_basic_association_matrix;

    /// Have T and C been built for the current geometry?
    mutable bool Algebraic_transformation_is_built;

    /// Build T and C if required.
    void ensure_algebraic_transformation() const;

    /// Assemble T directly from the global and basic dof definitions, then
    /// transpose it to obtain C.
    void build_algebraic_transformation() const;

    /// Return reference location, reference tangent, outward unit reference
    /// normal, and endpoint vertices for one of the three basic edges.
    ///
    /// The edge parameters are:
    ///   edge 0: s=(0,u),     node 2 -> node 1;
    ///   edge 1: s=(u,0),     node 2 -> node 0;
    ///   edge 2: s=(u,1-u),   node 1 -> node 0.
    static void get_reference_edge_geometry(const unsigned& edge,
					    const double& u,
                                            Vector<double>& s,
                                            Vector<double>& tangent,
                                            Vector<double>& outward_normal,
                                            unsigned& start_vertex,
                                            unsigned& end_vertex);

    /// Basic-edge interpolation abscissa in the ordering used by the full
    /// basic basis.
    static double edge_dof_abscissa(const unsigned& edge,
                                    const unsigned& point_index);

    /// Quintic Hermite basis ordered as
    /// [f(0),f(1),f'(0),f'(1),f''(0),f''(1)].
    static void quintic_hermite(const double& u, Vector<double>& h);

    /// Derivative with respect to u of the quintic Hermite basis.
    static void d_quintic_hermite(const double& u, Vector<double>& dh);

    /// Cubic Hermite basis ordered as
    /// [g(0),g(1),g'(0),g'(1)].
    static void cubic_hermite(const double& u, Vector<double>& h);

    /// Set a coefficient vector to zero with the correct global-dof length.
    void zero_global_dof_coefficients(Vector<double>& coefficients) const;

    /// destination += scale * source.
    static void add_scaled(Vector<double>& destination,
                           const Vector<double>& source,
                           const double& scale);

    /// Copy a coefficient vector into one row of a matrix.
    static void put_row(DenseMatrix<double>& matrix,
                        const unsigned& row,
                        const Vector<double>& coefficients);

    /// Coefficients of w at a vertex in the global dof vector.
    void vertex_value_coefficients(const unsigned& vertex,
                                   Vector<double>& coefficients) const;

    /// Coefficients of a physical directional first derivative
    /// direction.grad_x(w) at a vertex.
    void vertex_physical_first_derivative_coefficients(
      const unsigned& vertex,
      const Vector<double>& direction,
      Vector<double>& coefficients) const;

    /// Coefficients of the physical mixed second derivative
    /// a^T H_x(w) b at a vertex.
    void vertex_physical_second_derivative_coefficients(
      const unsigned& vertex,
      const Vector<double>& a,
      const Vector<double>& b,
      Vector<double>& coefficients) const;

    /// Coefficients of q.grad_s(w) at a vertex.
    void vertex_reference_first_derivative_coefficients(
      const unsigned& vertex,
      const Vector<double>& q,
      Vector<double>& coefficients) const;

    /// Coefficients of a^T H_s(w) b at a vertex, including the Hessian of the
    /// coordinate map.
    void vertex_reference_second_derivative_coefficients(
      const unsigned& vertex,
      const Vector<double>& a,
      const Vector<double>& b,
      Vector<double>& coefficients) const;

    /// Algebraic coefficients of the quintic value trace at an edge point.
    void edge_value_coefficients(const unsigned& edge,
                                 const double& u,
                                 Vector<double>& coefficients) const;

    /// Algebraic coefficients of the derivative with respect to the chosen
    /// edge parameter of the quintic value trace.
    void edge_value_parameter_derivative_coefficients(
      const unsigned& edge,
      const double& u,
      Vector<double>& coefficients) const;

    /// Algebraic coefficients of the basic outward-normal derivative at an
    /// edge point.
    void edge_basic_normal_derivative_coefficients(
      const unsigned& edge,
      const double& u,
      Vector<double>& coefficients) const;

    /// Transform reference derivatives of an arbitrary list of functions to
    /// Eulerian derivatives.
    double transform_basis_derivatives_reference_to_eulerian(
      const Vector<double>& s,
      const unsigned& nfunction,
      const DShape& d_reference,
      const DShape& d2_reference,
      DShape& d_eulerian,
      DShape& d2_eulerian,
      const bool& compute_second_derivatives) const;

    // [zdec] better name for this
    /// Evaluate the final global basis and its reference derivatives from C.
    void evaluate_global_basis_wrt_reference(
      const Vector<double>& s,
      Shape& value,
      DShape& dvalue,
      DShape& d2value,
      const bool& compute_second_derivatives) const;

    /// Scatter a flat list of global basis values into oomph-lib's nodal and
    /// bubble Shape containers, applying the curved-edge node permutation.
    void scatter_values(const Shape& value, Shape& psi, Shape& bpsi) const;

    /// Scatter first derivatives.
    void scatter_first_derivatives(const DShape& dvalue,
                                   DShape& dpsi,
                                   DShape& dbpsi) const;

    /// Scatter second derivatives.
    void scatter_second_derivatives(const DShape& d2value,
                                    DShape& d2psi,
                                    DShape& d2bpsi) const;
  };

} // namespace oomph

#endif

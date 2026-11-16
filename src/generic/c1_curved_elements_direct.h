// ====================================================================
// Algebraic direct construction of the Bernadou curved Bell basis.
// ====================================================================

#ifndef OOMPH_DIRECT_BERNADOU_ELEMENT_BASIS_H
#define OOMPH_DIRECT_BERNADOU_ELEMENT_BASIS_H

#include "Vector.h"
#include "shape.h"
#include "bell_element_basis.h"
#include "c1_plate_helper.h"

namespace oomph
{

  /// Non-templated runtime interface for Bernadou curved Bell bases.
  ///
  /// The owning finite element can therefore store a single pointer while
  /// choosing between the P3 and P5 boundary representations at runtime.
  class TemplateFreeBernadouElementBasis
  {
  public:
    /// Vertices of a two-dimensional triangle.
    typedef Vector<Vector<double>> VertexList;

    /// Virtual destructor.
    virtual ~TemplateFreeBernadouElementBasis() {}

    /// Set up geometry and basis coefficients. Vertices must already be
    /// ordered with the curved edge from vertex 0 to vertex 1; s_ubar and
    /// s_obar are their respective curve parameters.
    /// [zdec] Check this ordering comment
    virtual void upgrade_element(
      const VertexList& vertices,
      const double& s_ubar,
      const double& s_obar,
      const C1PlateHelper::CurvedEdgeEnumeration& curved_edge,
      const C1CurviLine& parametric_curve) = 0;

    /// Evaluate nodal and bubble basis values at element-local coordinates s.
    /// Nodal DOF order is [w, w_x, w_y, w_xx, w_xy, w_yy].
    virtual void shape(const Vector<double>& s,
                       Shape& nodal_basis,
                       Shape& bubble_basis) const = 0;

    /// Evaluate values and physical first derivatives at element-local s.
    /// Return the signed determinant of the coordinate-map Jacobian.
    virtual double d_shape_dx(const Vector<double>& s,
                              Shape& nodal_basis,
                              Shape& bubble_basis,
                              DShape& dnodal_basis,
                              DShape& dbubble_basis) const = 0;

    /// Evaluate values and physical first and second derivatives.
    /// Second derivatives are ordered [xx, xy, yy]; return the signed
    /// determinant of the coordinate-map Jacobian.
    virtual double d2_shape_dx2(const Vector<double>& s,
                                Shape& nodal_basis,
                                Shape& bubble_basis,
                                DShape& dnodal_basis,
                                DShape& dbubble_basis,
                                DShape& d2nodal_basis,
                                DShape& d2bubble_basis) const = 0;

    /// Curve parameter at canonical vertex 0; requires an upgraded basis.
    virtual const double& get_s_ubar() const = 0;

    /// Curve parameter at canonical vertex 1; requires an upgraded basis.
    virtual const double& get_s_obar() const = 0;

    /// Curved edge in the owning element's local numbering, or none.
    virtual C1PlateHelper::CurvedEdgeEnumeration curved_edge() const = 0;

    /// Map element-local coordinates s to physical coordinates x.
    virtual void coordinate_x(const Vector<double>& s,
                              Vector<double>& x) const = 0;

    /// Return J(i,alpha)=d x_i/d s_alpha in element-local coordinates.
    virtual void get_jacobian(const Vector<double>& s,
                              DenseMatrix<double>& jacobian) const = 0;

    /// Number of physical basis functions, including internal values.
    virtual unsigned n_basis_functions() const = 0;

    /// Number of internal value DOFs.
    virtual unsigned n_internal_dofs() const = 0;

    /// Number of reference monomials, equal to the full basic-basis size.
    virtual unsigned n_basic_basis_functions() const = 0;

    /// Fill the coefficients M in
    ///
    ///                 psi_global = M m,
    ///
    /// where m is the vector of reference monomials.
    virtual void get_full_association_matrix(
      DenseMatrix<double>& association_matrix) const = 0;

    /// Return the element-local coordinates of internal value DOF idof.
    virtual void get_internal_dofs_location(
      const unsigned& idof,
      Vector<double>& s) const = 0;
  };


  /// Algebraic direct construction of the P3/P5 Bernadou curved Bell basis.
  ///
  /// Let d_g denote the vector of physical/global degrees of freedom and d_b
  /// the full vector of basic-triangle degrees of freedom. Similarly, let psi_g
  /// denote the vector of physical/global basis functions and psi_b the vector
  /// of basic-triangle basis functions. The trace restrictions and the chain
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
  ///                 psi_g = T^T psi_b = T^T A m = M m,
  ///
  /// where m contains reference monomials and psi_b = A m. Only M is
  /// retained per element. T is temporary during setup, and A is shared
  /// statically between elements of the same boundary order.
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
  class DirectBernadouElementBasis : public TemplateFreeBernadouElementBasis
  {
  public:
    /// Non-templated runtime interface.
    typedef TemplateFreeBernadouElementBasis Base;

    /// Triangle vertices in canonical order.
    typedef Base::VertexList VertexList;

    static_assert(BOUNDARY_ORDER == 3 || BOUNDARY_ORDER == 5,
                  "DirectBernadouElementBasis is implemented only for "
                  "boundary orders 3 and 5.");

    /// Construct an unupgraded basis; upgrade_element() builds geometry and M.
    DirectBernadouElementBasis();

    /// Empty destructor.
    virtual ~DirectBernadouElementBasis() {}

    /// Upgrade the geometry and eagerly construct M=T^T A.
    /// verts must already have the curved edge from vertex 0 to vertex 1;
    /// su and so are their respective curve parameters. curved_edge records
    /// the edge in the owning element's local numbering.
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

    /// Copy the stored M, with dimensions n_basis_functions() by
    /// n_basic_basis_functions(), in canonical basis order.
    void get_full_association_matrix(
      DenseMatrix<double>& conversion_matrix) const override;

    /// Curve parameter at canonical vertex 0; requires an upgraded basis.
    const double& get_s_ubar() const override;
    /// Curve parameter at canonical vertex 1; requires an upgraded basis.
    const double& get_s_obar() const override;

    /// Curved edge in the owning element's numbering, or none before upgrade.
    C1PlateHelper::CurvedEdgeEnumeration curved_edge() const override
    {
      return Curved_edge;
    }

    /// Map element-local coordinates s to physical coordinates x.
    void coordinate_x(const Vector<double>& s,
                      Vector<double>& x) const override;

    /// Return J(i,alpha)=d x_i/d s_alpha in element-local coordinates.
    void get_jacobian(const Vector<double>& s,
                      DenseMatrix<double>& jacobian) const override;

    /// Polynomial degree of the curved-boundary geometry: 3 or 5.
    unsigned boundary_order() const
    {
      return BOUNDARY_ORDER;
    }

    /// Degree of the reference monomials: 7 or 9.
    unsigned basic_basis_order() const
    {
      return BOUNDARY_ORDER + 4;
    }

    /// Number of physical basis functions: 18 nodal plus 3 or 10 internal.
    unsigned n_basis_functions() const override
    {
      return (BOUNDARY_ORDER == 3 ? 21 : 28);
    }

    /// Number of internal value DOFs: 3 for P3, 10 for P5.
    unsigned n_internal_dofs() const override
    {
      return n_basis_functions() - 18;
    }

    /// Number of reference monomials/full basic DOFs: 36 or 55.
    unsigned n_basic_basis_functions() const override
    {
      return (BOUNDARY_ORDER + 5) * (BOUNDARY_ORDER + 6) / 2;
    }

    /// Return the element-local coordinates of internal value DOF idof.
    void get_internal_dofs_location(
      const unsigned& idof,
      Vector<double>& s) const override;

  private:

    /// Copy construction is disabled.
    DirectBernadouElementBasis(const DirectBernadouElementBasis&);

    /// Assignment is disabled.
    void operator=(const DirectBernadouElementBasis&);

    /// Physical vertices in canonical order: curved edge from vertex 0 to 1.
    VertexList Vertices;

    /// Curve parameter at canonical vertex 0.
    double S_ubar;

    /// Curve parameter at canonical vertex 1.
    double S_obar;

    /// Curved edge in element-local numbering; none until upgrade.
    C1PlateHelper::CurvedEdgeEnumeration Curved_edge;

    /// Affine mapping vector vertex_2 - vertex_0.
    double A1_vector[2];

    /// Affine mapping vector vertex_2 - vertex_1.
    double B2_vector[2];

    /// Coefficients of the vector-valued correction polynomial Q_s(s0) in
    ///
    ///  (x0,x1)
    ///    = F(s0,s1) = v2 - A1*s0 - B2*s1
    ///            + 0.5*s0*s1*(Q_t(s1) + Q_s(s0)).
    ///
    /// The first index is physical component, the second power 0..3.
    /// For P3 the quadratic and cubic coefficients are zero.
    double Correction_s[2][4];

    /// Coefficients of Q_t(s1), with the same indexing as Correction_s.
    double Correction_t[2][4];

    /// Coefficients M of the final basis in the monomial basis:
    ///
    ///                 psi_global = M m,
    ///                 M = T^T A.
    DenseMatrix<double> Global_basis_monomial_coefficients;

    /// Store canonical vertices and form the polynomial mapping coefficients.
    void set_up_geometry(const VertexList& vertices,
                         const double& s_ubar,
                         const double& s_obar,
                         const C1PlateHelper::CurvedEdgeEnumeration& curved_edge,
                         const C1CurviLine& parametric_curve);

#ifdef PARANOID
    /// Check curve endpoints and nondegenerate endpoint tangent directions.
    void self_check_geometry(const Vector<double>& curve_start,
                             const Vector<double>& curve_end,
                             const Vector<double>& a2,
                             const Vector<double>& b1) const;
#endif

    /// Evaluate a cubic correction polynomial with coefficients in power order.
    static double cubic_polynomial_value(const double coefficient[4],
                                         const double& u);

    /// Evaluate the first derivative of a cubic correction polynomial.
    static double cubic_polynomial_first_derivative(
      const double coefficient[4],
      const double& u);

    /// Evaluate the second derivative of a cubic correction polynomial.
    static double cubic_polynomial_second_derivative(
      const double coefficient[4],
      const double& u);

    /// Canonical basic-to-physical mapping x=F(s).
    void position_x(const Vector<double>& s,
                    Vector<double>& x) const;

    /// Jacobian J(i,alpha)=d x_i/d s_alpha in canonical basic coordinates.
    void jacobian_dx_ds(const Vector<double>& s,
                        DenseMatrix<double>& jacobian) const;

    /// Hessian H(i,alpha,beta)=d2 x_i/(d s_alpha d s_beta).
    void hessian_d2x_ds2(const Vector<double>& s,
                         RankThreeTensor<double>& hessian) const;

    /// Map element-local coordinates to the canonical curved-edge orientation.
    void permute_to_basic_coordinates(Vector<double>& s) const;

    /// Return the constant Jacobian of the local-to-canonical permutation.
    void jacobian_ds_basic_ds_local(DenseMatrix<double>& jacobian) const;

    /// Return the cyclic node shift from canonical to element-local numbering.
    void nodal_index_shift(unsigned& index_shift) const;

    /// Evaluate monomials in descending total degree, increasing s1 power.
    void monomials(const Vector<double>& s,
                   Shape& monomial) const;

    /// Evaluate monomials and first derivatives in canonical coordinates.
    /// If requested, also evaluate second derivatives in [00, 01, 11] order.
    void monomials_and_derivatives(const Vector<double>& s,
                                   Shape& monomial,
                                   DShape& dmonomial,
                                   DShape& d2monomial,
                                   const bool& compute_second_derivatives) const;

    /// Return the fixed matrix A mapping degree-(BOUNDARY_ORDER+4)
    /// monomials to the full basic-triangle cardinal basis.
    /// The function-local static is initialised once, safely under concurrent
    /// calls, and shared by every element of the same boundary order.
    /// Specialisations are defined in c1_curved_association_matrix_direct.cc.
    static const DenseMatrix<double>& monomial_to_basic_matrix();

    /// Build the temporary rectangular matrix T such that d_b=T d_g.
    void build_basic_dof_transformation_matrix(
      DenseMatrix<double>& transformation_matrix) const;

    /// Build and store M=T^T A. Called eagerly by upgrade_element().
    void build_monomial_association_matrix();

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

    /// Coefficients of w at a vertex in the global dof vector (== 1.0).
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

    /// Algebraic coefficients of the quintic value trace at an edge point. This
    /// explicitly imposes the constraint that edge values are interpolated
    /// quintically.
    void edge_value_coefficients(const unsigned& edge,
                                 const double& u,
                                 Vector<double>& coefficients) const;

    /// Algebraic coefficients of the derivative with respect to the edge
    /// parameter of the quintic value trace.
    /// (Needed for the normal derivative interpolation on straight edges)
    void edge_value_parameter_derivative_coefficients(
      const unsigned& edge,
      const double& u,
      Vector<double>& coefficients) const;

    /// Algebraic coefficients of the basic outward-normal derivative at an edge
    /// point. This explicitly imposes that the physical normal on straight
    /// edges and the reference normal on the curved edge are cubically
    /// interpolated by the nodal data.
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

    /// Evaluate global basis values from M without computing derivatives.
    void evaluate_global_basis_values_wrt_reference(
      const Vector<double>& s,
      Shape& value) const;

    /// Evaluate the global basis and reference derivatives from M.
    /// s is in canonical coordinates; second derivatives are optional.
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

  /// Shared 36-by-36 association matrix for a cubic boundary.
  template<>
  const DenseMatrix<double>&
  DirectBernadouElementBasis<3>::monomial_to_basic_matrix();

  /// Shared 55-by-55 association matrix for a quintic boundary.
  template<>
  const DenseMatrix<double>&
  DirectBernadouElementBasis<5>::monomial_to_basic_matrix();

} // namespace oomph

#endif

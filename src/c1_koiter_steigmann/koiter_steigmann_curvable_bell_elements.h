// Header file for the Koiter Steigmann - Bell/Bernadou elements
#ifndef OOMPH_KS_CURVABLE_BELL_ELEMENTS_HEADER
#define OOMPH_KS_CURVABLE_BELL_ELEMENTS_HEADER

// oomph-lib headers
#include "koiter_steigmann_equations.h"
#include "src/generic/subparametric_Telement.h"
#include "src/generic/oomph_definitions.h"

namespace oomph
{
  //===start of rotation helper class=========================================
  /// Helper class to contain all the rotation information in the element.
  class RotatedBoundaryHelper
  {
  public:
    /// Constructor: just initialise the member data to their defaults (zeros)
    RotatedBoundaryHelper(FiniteElement* const& parent_element_pt)
      : Parent_element_pt(parent_element_pt),
        Nnode(Parent_element_pt->nvertex_node()),
        Boundary_coordinate_of_node(3, 0.0),
        Nodal_boundary_parametrisation_pt(3, 0),
        Rotation_matrix_at_node(3, DenseMatrix<double>(6, 6, 0.0))
    {
    }

    /// Destructor
    ~RotatedBoundaryHelper() {}

    CurvilineGeomObject* nodal_boundary_parametrisation_pt(
      const unsigned& j_node)
    {
      return Nodal_boundary_parametrisation_pt[j_node];
    }


    /// Add a new boundary parametrisation to nodes all the nodes in the
    /// vector node_on_boundary
    void set_nodal_boundary_parametrisation(
      const Vector<unsigned>& node_on_boundary,
      const Vector<double>& boundary_coord_of_node,
      CurvilineGeomObject* const& boundary_parametrisation_pt)
    {
      // Loop over all the nodes in node_on_boundary and add the boundary
      // pointer to their vector of boundaries
      unsigned n_node = node_on_boundary.size();
      for (unsigned j = 0; j < n_node; j++)
      {
        // The j-th node on the boundary
        unsigned j_node = node_on_boundary[j];

        // Set the boundary parametrisation data pointer for this node
        Nodal_boundary_parametrisation_pt[j_node] = boundary_parametrisation_pt;

        // Set the coordinate of node j on this boundary
        Boundary_coordinate_of_node[j_node] = boundary_coord_of_node[j];

        update_rotation_matrices();
      } // end of loop over nodes in node_on_boundary [j]
    } // end of set_nodal_boundary_parametrisation()


    /// Update all rotation matrices (checks if they are needed unless flag is
    /// true)
    void update_rotation_matrices()
    {
      // [zdec] hard coded the three vertex nodes
      unsigned n_vertex = 3;
      // Loop over each vertex
      for (unsigned j_node = 0; j_node < n_vertex; j_node++)
      {
        // If this node does not have a parametrisation (the pointer is still
        // null) skip over it, otherwise we go on to fill out the rotation
        // matrix
        if (!nodal_boundary_parametrisation_pt(j_node))
        {
          continue;
        }

        // Initialise the two basis vectors and their jacobians
        Vector<Vector<double>> bi(2, Vector<double>(2, 0.0));
        Vector<DenseMatrix<double>> dbidx(2, DenseMatrix<double>(2, 2, 0.0));

        // Our new coordinate system:
        //     (l, s)=(normal component, tangent component)
        // which we define in terms of basis vectors (rescaled)
        //     ni=dxi/dl / |n|           <-- Jacobian col 1
        //     ti=dxi/ds / |t|           <-- Jacobian col 2
        // and their derivatives
        //     dnidxj=d/dxj(dxi/dl / |n|) <-- Hessian `col' 1
        //     dtidxj=d/dxj(dxi/ds / |t|) <-- Hessian `col' 2

        // [zdec] we use i and j for brevity
        // but it should be alpha & beta
        // Need to write up how the transformation is done

        // Storage for our basis and derivatives
        Vector<double> ni(2, 0.0);
        Vector<double> ti(2, 0.0);
        Vector<double> dnids(2, 0.0);
        Vector<double> dtids(2, 0.0);

        // All tensors assumed evaluated on the boundary
        // Jacobian of inverse mapping
        DenseMatrix<double> jac_inv(2, 2, 0.0);
        // Hessian of mapping [zdec] (not needed because...)
        Vector<DenseMatrix<double>> hess(2, DenseMatrix<double>(2, 2, 0.0));
        // Hessian of inverse mapping [zdec] (...this can be found by hand)
        Vector<DenseMatrix<double>> hess_inv(2, DenseMatrix<double>(2, 2, 0.0));

        // The basis is defined in terms of the boundary parametrisation
        Vector<double> boundary_coord = {Boundary_coordinate_of_node[j_node]};
        CurvilineGeomObject* boundary_pt =
          Nodal_boundary_parametrisation_pt[j_node];
        Vector<double> x(2, 0.0);
        Vector<double> dxids(2, 0.0);
        Vector<double> d2xids2(2, 0.0);

        // Get position (debug)
        boundary_pt->position(boundary_coord, x);
        // Get tangent vector
        boundary_pt->dposition(boundary_coord, dxids);
        // Get second derivative
        boundary_pt->d2position(boundary_coord, d2xids2);

        double mag_t = sqrt(dxids[0] * dxids[0] + dxids[1] * dxids[1]);
        // ti is the normalised tangent vector
        ti[0] = dxids[0] / mag_t;
        ti[1] = dxids[1] / mag_t;
        // Derivative of (normalised) tangent
        dtids[0] = d2xids2[0] / std::pow(mag_t, 2) -
                   (dxids[0] * d2xids2[0] + dxids[1] * d2xids2[1]) * dxids[0] /
                     std::pow(mag_t, 4);
        dtids[1] = d2xids2[1] / std::pow(mag_t, 2) -
                   (dxids[0] * d2xids2[0] + dxids[1] * d2xids2[1]) * dxids[1] /
                     std::pow(mag_t, 4);
        // n = (t x e_z) implies
        ni[0] = ti[1];
        ni[1] = -ti[0];
        // Same for dnids
        dnids[0] = dtids[1];
        dnids[1] = -dtids[0];

        // Need inverse of mapping to calculate ds/dxi ----------------
        //   /  dx/dl  dx/ds  \ -1  ___  __1__ /  dy/ds -dx/ds \ .
        //   \  dy/dl  dy/ds  /     ---   det  \ -dy/dl  dx/dl /
        //
        //                          ___  /  dl/dx  dl/dy  \ .
        //                          ---  \  ds/dx  ds/dy  /
        //
        // Fill out inverse of Jacobian
        double det = (ni[0] * ti[1] - ni[1] * ti[0]);
        jac_inv(0, 0) = ti[1] / det;
        jac_inv(0, 1) = -ti[0] / det;
        jac_inv(1, 0) = -ni[1] / det;
        jac_inv(1, 1) = ni[0] / det;

        // Fill out the Hessian
        // (unneeded -- can calculate the inverse components by hand)
        for (unsigned alpha = 0; alpha < 2; alpha++)
        {
          // hess[alpha](0,0) = 0.0;
          hess[alpha](0, 1) = dnids[alpha];
          hess[alpha](1, 0) = dnids[alpha];
          hess[alpha](1, 1) = dtids[alpha];
        }

        // Fill out inverse of Hessian
        // H^{-1}abg = J^{-1}ad Hdez J^{-1}eb J^{-1}zg
        for (unsigned alpha = 0; alpha < 2; alpha++)
        {
          for (unsigned beta = 0; beta < 2; beta++)
          {
            for (unsigned gamma = 0; gamma < 2; gamma++)
            {
              for (unsigned alpha2 = 0; alpha2 < 2; alpha2++)
              {
                for (unsigned beta2 = 0; beta2 < 2; beta2++)
                {
                  for (unsigned gamma2 = 0; gamma2 < 2; gamma2++)
                  {
                    hess_inv[alpha](beta, gamma) -=
                      jac_inv(alpha, alpha2) * hess[alpha2](beta2, gamma2) *
                      jac_inv(beta2, beta) * jac_inv(gamma2, gamma);
                  }
                }
              }
            }
          }
        }

        // Fill in the rotation matrix using the new basis
        fill_in_rotation_matrix_at_node_with_basis(j_node, jac_inv, hess_inv);

      } // end loop over vertices
    } // end of update_rotation_matrices()


    /// Access function to fill out rot_mat using rotation matrix
    void get_rotation_matrix_at_node(const unsigned& j_node,
                                     DenseMatrix<double>& rot_mat)
    {
      rot_mat = Rotation_matrix_at_node[j_node];
    }

  private:
    /// Helper function to fill in the rotation matrix for a given basis
    void fill_in_rotation_matrix_at_node_with_basis(
      const unsigned& j_node,
      const DenseMatrix<double>& jac_inv,
      const Vector<DenseMatrix<double>>& hess_inv)
    {
      // Rotation matrix, b constructed using submatrices b1, b12, b22
      DenseMatrix<double> b1(2, 2, 0.0), b22(3, 3, 0.0), b12(2, 3, 0.0);

      // Fill in the submatrices
      // Loop over the rotated first derivatives
      for (unsigned mu = 0; mu < 2; mu++)
      {
        // Loop over the unrotated first derivatives
        for (unsigned alpha = 0; alpha < 2; alpha++)
        {
          // Fill in b1 - the Jacobian
          // Fill in the affine rotation of the first derivatives
          b1(mu, alpha) = jac_inv(mu, alpha);

          // Loop over unrotated second derivatives
          for (unsigned beta = 0; beta < 2; ++beta)
          {
            // Avoid double counting the cross derivative
            if (alpha <= beta)
            {
              // Define column index
              const unsigned col = alpha + beta;

              // Fill in the non-affine part of the rotation of the first
              // derivatives
              b12(mu, col) += hess_inv[mu](alpha, beta);
              // [zdec] debug mixed derivative -- add extra
              if (alpha < beta)
              {
                // b12(mu, col) -= hess_inv[mu](alpha, beta);
              }
              // Loop over the rotated second derivatives
              for (unsigned nu = 0; nu < 2; nu++)
              {
                // // Avoid double counting the cross derivative
                // if (mu <= nu)
                {
                  // Fill in b22 - the Affine part of the Jacobian derivative
                  // Redefine row index for the next submatrix
                  unsigned row_b22 = mu + nu;
                  // Fill in the affine part of the rotation of the second
                  // derivatives [zdec] if( beta>= alpha) ?
                  b22(row_b22, col) += jac_inv(mu, alpha) * jac_inv(nu, beta);
                }
              }
            }
          }
        }
      }

      // Fill in the submatrices to the full (6x6) matrix
      Rotation_matrix_at_node[j_node](0, 0) = 1.0;
      // Fill in b1 --- the affine contribution to rotation of the
      // first derivatives
      for (unsigned i = 0; i < 2; ++i)
      {
        for (unsigned j = 0; j < 2; ++j)
        {
          Rotation_matrix_at_node[j_node](1 + i, 1 + j) = b1(i, j);
        }
      }
      // Fill in b21 --- the non-affine (second derivative dependent)
      // rotation of the first derivatives
      for (unsigned i = 0; i < 2; ++i)
      {
        for (unsigned j = 0; j < 3; ++j)
        {
          Rotation_matrix_at_node[j_node](1 + i, 3 + j) = b12(i, j);
        }
      }
      // Fill in b22 --- the rotation of the second derivatives
      for (unsigned i = 0; i < 3; ++i)
      {
        for (unsigned j = 0; j < 3; ++j)
        {
          Rotation_matrix_at_node[j_node](3 + i, 3 + j) = b22(i, j);
        }
      }
    } // end fill_in_rotation_matrix_at_node_with_basis

    /// Pointer to the `parent' finite element which this is a helper force
    FiniteElement* Parent_element_pt;

    /// The number of nodes (that we store rotation data for) in the fvk element
    /// that uses this helper
    unsigned Nnode;

    /// Vector containing boundary parametrised location for each node
    Vector<double> Boundary_coordinate_of_node;

    /// Vector containing boundary parametrisation at each node
    Vector<CurvilineGeomObject*> Nodal_boundary_parametrisation_pt;

    /// Vector containing <rotation matrix at each node>
    Vector<DenseMatrix<double>> Rotation_matrix_at_node;
  };

  //---end of rotation helper class-------------------------------------------




  //===============================================================================
  /// KoiterSteigmannC1CurvableBellElement elements are a subparametric
  /// scheme with  linear Lagrange interpolation for approximating the geometry
  /// and the C1-functions for approximating variables.
  //==============================================================================

  // [zdec] NNODE_1D(=2) should not be required here due to the fact that all
  // node should be vertex nodes? Is any templating needed???
  // template<unsigned DIM,
  //          unsigned NNODE_1D,
  //          unsigned BOUNDARY_ORDER,
  //          template<unsigned DIM_, unsigned NNODE_1D_>
  //          class PLATE_EQUATIONS>
  class KoiterSteigmannC1CurvableBellElement
    : public virtual CurvableBellElement<2>,
      public virtual KoiterSteigmannEquations
  {
  public:

    //----------------------------------------------------------------------
    // Class construction

    /// Constructor: Call constructors for C1CurvedBellElement and
    /// KoiterSteigmannEquations
    KoiterSteigmannC1CurvableBellElement()
      : CurvableBellElement<Nnode_1D>(Nfield, Field_is_bell_interpolated),
        KoiterSteigmannEquations()
    {
      // Use the higher order integration scheme
      delete this->integral_pt();
      // Do we want something that is order 8 instead?
      TGauss<2, 9>* new_integral_pt = new TGauss<2, 9>;
      this->set_integration_scheme(new_integral_pt);

      // Rotated dof helper
      Rotated_boundary_helper_pt = new RotatedBoundaryHelper(this);
    }

    /// Destructor
    ~KoiterSteigmannC1CurvableBellElement()
    {
      // Delete the rotated bonudary helper we made for this element
      delete Rotated_boundary_helper_pt;
    }

    /// Broken copy constructor
    KoiterSteigmannC1CurvableBellElement(
      const KoiterSteigmannC1CurvableBellElement& dummy)
    {
      BrokenCopy::broken_copy("KoiterSteigmannC1CurvableBellElement");
    }

    /// Broken assignment operator
    void operator=(
      const KoiterSteigmannC1CurvableBellElement&)
    {
      BrokenCopy::broken_assign("KoiterSteigmannC1CurvableBellElement");
    }


    //----------------------------------------------------------------------
    // Output and documentation

    /// \short Output function:
    ///  x,y,u   or    x,y,z,u
    void output(std::ostream& outfile)
    {
      KoiterSteigmannEquations::output(outfile);
    }

    ///  \short Output function:
    ///   x,y,u   or    x,y,z,u at n_plot^DIM plot points
    void output(std::ostream& outfile, const unsigned& n_plot)
    {
      KoiterSteigmannEquations::output(outfile, n_plot);
    }

    ///  \short Output function:
    ///   x,y,u   or    x,y,z,u at n_plot^DIM plot points
    void output_interpolated_exact_soln(
      std::ostream& outfile,
      FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
      const unsigned& n_plot);

    /// \short C-style output function:
    ///  x,y,u   or    x,y,z,u
    void output(FILE* file_pt)
    {
      KoiterSteigmannEquations::output(file_pt);
    }

    ///  \short C-style output function:
    ///   x,y,u   or    x,y,z,u at n_plot^DIM plot points
    void output(FILE* file_pt, const unsigned& n_plot)
    {
      KoiterSteigmannEquations::output(file_pt, n_plot);
    }


    /// \short Output function for an exact solution:
    ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^DIM plot points
    void output_fct(std::ostream& outfile,
                    const unsigned& n_plot,
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
    {
      KoiterSteigmannEquations::output_fct(
        outfile, n_plot, exact_soln_pt);
    }

    /// \short Output function for a time-dependent exact solution.
    ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^DIM plot points
    /// (Calls the steady version)
    void output_fct(std::ostream& outfile,
                    const unsigned& n_plot,
                    const double& time,
                    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
    {
      KoiterSteigmannEquations::output_fct(
        outfile, n_plot, time, exact_soln_pt);
    }


    // /// \short enum to enumerate the possible edges that could be curved
    // typedef typename MyC1CurvedElements::Edge Edge;

    // [zdec]
    // /// \short Get the pointer to the Curved shape class data member
    // const MyC1CurvedElements::BernadouElementBasis<BOUNDARY_ORDER>* curved_shape_pt()
    // {
    //   return &Curved_shape;
    // };

    // /// \short get the coordinate
    // inline void get_coordinate_x(const Vector<double>& s,
    //                              Vector<double>& x) const;

    // /// \short get the coordinate i
    // double interpolated_x(const Vector<double>& s, const unsigned& i) const
    // {
    //   Vector<double> r(2);
    //   get_coordinate_x(s, r);
    //   return r[i];
    // }


    //----------------------------------------------------------------------
    // Jacobian and residual contributions

    /// Add the element's contribution to its residual vector (wrapper) with
    /// cached association matrix
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Store the expensive-to-construct matrix
      this->store_association_matrix();
      // Call the generic routine with the flag set to 1
      KoiterSteigmannEquations::fill_in_contribution_to_residuals(residuals);
      // Remove the expensive-to-construct matrix
      this->delete_association_matrix();
    }

    /// Add the element's contribution to its residual vector and
    /// element Jacobian matrix (wrapper) with caching of association matrix
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      // Store the expensive-to-construct matrix
      this->store_association_matrix();
      // Call the generic routine with the flag set to 1
      KoiterSteigmannEquations::fill_in_contribution_to_jacobian(residuals,
                                                                 jacobian);
      // Remove the expensive-to-construct matrix
      this->delete_association_matrix();
    }


    //----------------------------------------------------------------------
    // Geometry and boundaries

    /// Function useful for setting boundary conditions that streamlines the
    /// boundary condition setting process by pinning the k-th value of the i-th
    /// displacement field at all nodes on boundary b to the value prescribed by
    /// specified_u_ik_pt. This allows the user to set BCs without understanding
    /// the underlying geometry elements of how they are integrated with
    /// KoiterSteigmannEquations.
    ///
    /// As the dofs are Hermite, the user needs to take care to ensure they are
    /// setting values that are implicitly prescribed at nodes by the higher
    /// order continuous information AS WELL AS to be careful to track the
    /// coordinate frame the Hermite data is describing (e.g. (x,y) vs
    /// (normal,tangent)).
    void set_boundary_condition(const unsigned& i_field,
                                const unsigned& k_type,
                                const unsigned& b_boundary,
                                const ScalarFctPt& specified_u_ik_pt);

    /// Get the zeta coordinate
    inline void interpolated_zeta(const Vector<double>& s,
                                  Vector<double>& zeta) const
    {
      // If there is a macro element use it
      if (this->Macro_elem_pt != 0)
      {
        this->get_x_from_macro_element(s, zeta);
      }
      // Otherwise interpolate zeta_nodal using the shape functions
      else
      {
        interpolated_x(s, zeta);
      }
    }

    /// Return true if the element has been upgraded to interpolate a curved
    /// boundary
    bool element_is_curved() const
    {
      return CurvableBellElement<Nnode_1D>::element_is_curved();
    }

    /// Upgrade the Bell element to a curved Bernadou element
    virtual void upgrade_element_to_curved(
      const MyC1CurvedElements::Edge& curved_edge,
      const double& s_ubar,
      const double& s_obar,
      CurvilineGeomObject* parametric_edge,
      const unsigned& boundary_order)
    {
      CurvableBellElement<Nnode_1D>::upgrade_element_to_curved(
        curved_edge, s_ubar, s_obar, parametric_edge, boundary_order);
    }


    //----------------------------------------------------------------------
    // Member data acess functions

    /// Access function to rotated boundary helper object
    RotatedBoundaryHelper* rotated_boundary_helper_pt()
    {
      return Rotated_boundary_helper_pt;
    }

    /// Access the number of fields
    unsigned nfield()
    {
      return Nfield;
    }

    /// Required  # of `values' (pinned or dofs)
    /// at node n
    inline unsigned required_nvalue() const
    {
      return Initial_Nvalue;
    }




  protected:

    /// Transform the shape functions so that they correspond to the new rotated
    /// dofs
    inline void rotate_shape(Shape& shape) const;

    /// Transform the shape functions and first derivatives so that they
    /// correspond to the new rotated dofs
    inline void rotate_shape(Shape& shape, DShape& dshape) const;

    /// Transform the shape functions, first and second derivatives so that they
    /// correspond to the new rotated dofs
    inline void rotate_shape(Shape& shape,
                             DShape& dshape,
                             DShape& d2shape) const;


    /// Required # of `values' (pinned or dofs) at node n
    inline unsigned required_nvalue(const unsigned& n) const
    {
      return Initial_Nvalue;
    }


    //----------------------------------------------------------------------------
    // Interface to KoiterSteigmannEquations (can this all be (static) data?)

    /// Interface to return a vector of the index of each displacement unkonwn
    /// in the grander scheme of unknowns
    virtual Vector<unsigned> u_field_indices() const
    {
      return {0, 1, 2};
    }

    /// Interface to return a vector of the index of the u_i displacement
    /// unkonwn in the grander scheme of unknowns
    virtual unsigned u_i_field_index(const unsigned& i_field) const
    {
      return i_field;
    }

    // [zdec] should always be three
    /// Interface to return the number of nodes used by u
    virtual unsigned nu_node() const
    {
      return CurvableBellElement<Nnode_1D>::nnode_for_field(
        u_i_field_index(0));
    }

    // [zdec] should always be 1,2,3? also not used in KoiterSteigmannEquations
    // yet (not very future proof)
    /// Interface to get the local indices of the nodes used by u
    virtual Vector<unsigned> get_u_node_indices() const
    {
      return CurvableBellElement<Nnode_1D>::nodal_indices_for_field(
        u_i_field_index(0));
    }

    /// Interface to get the number of basis types for u at node j
    virtual unsigned nu_type_at_each_node() const
    {
      return CurvableBellElement<Nnode_1D>::nnodal_basis_type_for_field(
        u_i_field_index(0));
    }

    /// Interface to retrieve the value of u_i at node j of type k
    virtual double get_u_i_value_at_node_of_type(
      const unsigned& i_field,
      const unsigned& j_node,
      const unsigned& k_type) const
    {
      unsigned n_u_types = nu_type_at_each_node();
      unsigned nodal_type_index = i_field * n_u_types + k_type;
      return raw_nodal_value(j_node, nodal_type_index);
    }

    /// Interface to retrieve the t-th history value of u_i at node j of type k
    virtual double get_u_i_value_at_node_of_type(
      const unsigned& t_time,
      const unsigned& i_field,
      const unsigned& j_node,
      const unsigned& k_type) const
    {
      unsigned n_u_types = nu_type_at_each_node();
      unsigned nodal_type_index = i_field * n_u_types + k_type;
      return raw_nodal_value(t_time, j_node, nodal_type_index);
    }

    /// Interface to get the pointer to the internal data used to interpolate
    /// the i-th displacement field
    virtual Data* u_i_internal_data_pt(const unsigned& i_field) const
    {
      unsigned index =
        CurvableBellElement<Nnode_1D>::index_of_internal_data_for_field(
          i_field);
      return internal_data_pt(index);
    }

    /// Interface to get the number of internal types for the u fields
    virtual unsigned nu_type_internal() const
    {
      return CurvableBellElement<Nnode_1D>::ninternal_basis_type_for_field(
        u_field_indices()[0]);
    }

    /// (pure virtual) interface to retrieve the value of u_alpha of internal
    /// type k
    virtual double get_u_i_internal_value_of_type(
      const unsigned& i_field,
      const unsigned& k_type) const
    {
      unsigned index = index_of_internal_data_for_field(i_field);
      return CurvableBellElement<Nnode_1D>::internal_value_for_field_of_type(
        index, k_type);
    }

    /// (pure virtual) interface to retrieve the t-th history value of u_alpha
    /// of internal type k
    virtual double get_u_i_internal_value_of_type(
      const unsigned& t_time,
      const unsigned& i_field,
      const unsigned& k_type) const
    {
      unsigned index = index_of_internal_data_for_field(i_field);
      return CurvableBellElement<Nnode_1D>::internal_value_for_field_of_type(
        t_time, index, k_type);
    }


    /// (pure virtual) Out-of-plane basis functions at local coordinate s
    virtual void basis_u_koiter_steigmann(
      const Vector<double>& s,
      Shape& psi_n,
      Shape& psi_i) const;

    /// (pure virtual) Out-of-plane basis functions and derivs w.r.t. global
    /// coords at local coordinate s; return det(Jacobian of mapping)
    virtual double d2basis_u_eulerian_koiter_steigmann(
      const Vector<double>& s,
      Shape& psi_n,
      Shape& psi_i,
      DShape& dpsi_n_dx,
      DShape& dpsi_i_dx,
      DShape& d2psi_n_dx2,
      DShape& d2psi_i_dx2) const;

    /// (pure virtual) Out-of-plane basis/test functions at local coordinate s
    virtual void basis_and_test_u_koiter_steigmann(
      const Vector<double>& s,
      Shape& psi_n,
      Shape& psi_i,
      Shape& test_n,
      Shape& test_i) const;

    /// (pure virtual) Out-of-plane basis/test functions and first derivs w.r.t.
    /// to global coords at local coordinate s; return det(Jacobian of mapping)
    virtual double dbasis_and_dtest_u_eulerian_koiter_steigmann(
      const Vector<double>& s,
      Shape& psi_n,
      Shape& psi_i,
      DShape& dpsi_n_dx,
      DShape& dpsi_i_dx,
      Shape& test_n,
      Shape& test_i,
      DShape& dtest_n_dx,
      DShape& dtest_i_dx) const;

    /// (pure virtual) Out-of-plane basis/test functions and first/second derivs
    /// w.r.t. to global coords at local coordinate s;
    /// return det(Jacobian of mapping)
    virtual double d2basis_and_d2test_u_eulerian_koiter_steigmann(
      const Vector<double>& s,
      Shape& psi_n,
      Shape& psi_i,
      DShape& dpsi_n_dx,
      DShape& dpsi_i_dx,
      DShape& d2psi_n_dx2,
      DShape& d2psi_i_dx2,
      Shape& test_n,
      Shape& test_i,
      DShape& dtest_n_dx,
      DShape& dtest_i_dx,
      DShape& d2test_n_dx2,
      DShape& d2test_i_dx2) const;



    // Private Data Members
  private:

    /// Pointer to an instance of rotated boundary helper
    RotatedBoundaryHelper* Rotated_boundary_helper_pt;

    /// Number of nodes along each edge of the element (is always 2)
    static const unsigned Nnode_1D = 2;

    /// Static number of fields (is always 3)
    static const unsigned Nfield;

    /// Static bool vector with the Bell interpolation of the fields
    /// (always only w)
    static const std::vector<bool> Field_is_bell_interpolated;

    /// \short Static int that holds the number of variables at
    /// nodes: always the same
    static const unsigned Initial_Nvalue;

  };


  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////


  //==============================================================================
  /// Face geometry for the KoiterSteigmannC1CurvableBellElement elements:
  /// The spatial dimension of the face elements is one lower than that of the
  /// bulk element but they have the same number of points along their 1D edges.
  //==============================================================================
  // template<unsigned DIM,
  //          unsigned NNODE_1D,
  //          unsigned BOUNDARY_ORDER,
  //          template<unsigned DIM_, unsigned NNODE_1D_>
  //          class PLATE_EQUATIONS>
  template<>
  class FaceGeometry<KoiterSteigmannC1CurvableBellElement>
    : public virtual TElement<1, 2>
  {
  public:
    /// \short Constructor: Call the constructor for the
    /// appropriate lower-dimensional TElement
    FaceGeometry() : TElement<1, 2>() {}
  };


  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////



  //============================================================================
  /// Function useful for setting boundary conditions that streamlines the
  /// boundary condition setting process by pinning the k-th value of the i-th
  /// displacement field at all nodes on boundary b to the value prescribed by
  /// specified_u_ik_pt. This allows the user to set BCs without understanding
  /// the underlying geometry elements of how they are integrated with
  /// KoiterSteigmannEquations.
  ///
  /// As the dofs are Hermite, the user needs to take care to ensure they are
  /// setting values that are implicitly prescribed at nodes by the higher
  /// order continuous information AS WELL AS to be careful to track the
  /// coordinate frame the Hermite data is describing (e.g. (x,y) vs
  /// (normal,tangent)).
  //============================================================================
  void KoiterSteigmannC1CurvableBellElement::
  set_boundary_condition(const unsigned& i_field,
                              const unsigned& k_type,
                              const unsigned& b_boundary,
                              const ScalarFctPt& specified_u_ik_pt)
  {
    const unsigned first_nodal_type_index =
      this->first_nodal_type_index_for_field(i_field);
    const unsigned n_vertices = nu_node();

#ifdef PARANOID
    // Check that the dof number is a sensible value
    unsigned n_type = nu_type_at_each_node();
    if (k_type >= n_type)
    {
      throw OomphLibError(
        "Foppl von Karman elements only have 6 Hermite deflection degrees\
of freedom at internal points. They are {w ; w,x ; w,y ; w,xx ; w,xy ; w,yy}",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // These KS elements only have displacement dofs at vertices, we assume
    // these are the first n_vertices nodes
    for (unsigned n = 0; n < n_vertices; ++n)
    {
      // Get node
      Node* nod_pt = this->node_pt(n);
      // Check if it is on the boundary
      bool is_boundary_node = nod_pt->is_on_boundary(b_boundary);
      if (is_boundary_node)
      {
        // Extract nodal coordinates from node:
        Vector<double> x(2);
        x[0] = nod_pt->x(0);
        x[1] = nod_pt->x(1);
        // Get value
        double value;
        specified_u_ik_pt(x, value);
        // Pin and set the value
        nod_pt->pin(first_nodal_type_index + k_type);
        nod_pt->set_value(first_nodal_type_index + k_type, value);
      }
    }
  }



  //======================================================================
  /// (pure virtual) Basis functions at local coordinate s
  //======================================================================
  void KoiterSteigmannC1CurvableBellElement::
    basis_u_koiter_steigmann(const Vector<double>& s,
                               Shape& psi_n,
                               Shape& psi_i) const
  {
    throw OomphLibError("This still needs testing for curved elements.",
                        "void KoiterSteigmannEquations::\
shape_and_test_foeppl_von_karman(...)",
                        OOMPH_EXCEPTION_LOCATION);

    this->c1_basis(s, psi_n, psi_i);

    // Rotate the degrees of freedom
    rotate_shape(psi_n);
  }


  //======================================================================
  /// Define the shape functions and test functions and derivatives
  /// w.r.t. global coordinates and return Jacobian of mapping.
  ///
  /// Galerkin: Test functions = shape functions
  //======================================================================
  void KoiterSteigmannC1CurvableBellElement::
    basis_and_test_u_koiter_steigmann(const Vector<double>& s,
                                      Shape& psi_n,
                                      Shape& psi_i,
                                      Shape& test_n,
                                      Shape& test_i) const
  {
    // Use the c1 basis
    this->c1_basis(s, psi_n, psi_i);

    // Rotate the degrees of freedom
    rotate_shape(psi_n);

    // Galerkin
    // (Shallow) copy the basis functions
    test_n = psi_n;
    test_i = psi_i;
  }


  //======================================================================
  /// Fetch the basis functions and test functions and first derivatives
  /// w.r.t. global coordinates and return Jacobian of mapping.
  ///
  /// Galerkin: Test functions = shape functions
  //======================================================================
  double KoiterSteigmannC1CurvableBellElement::
    dbasis_and_dtest_u_eulerian_koiter_steigmann(const Vector<double>& s,
                                                 Shape& psi_n,
                                                 Shape& psi_i,
                                                 DShape& dpsi_n_dx,
                                                 DShape& dpsi_i_dx,
                                                 Shape& test_n,
                                                 Shape& test_i,
                                                 DShape& dtest_n_dx,
                                                 DShape& dtest_i_dx) const
  {
    // Get the basis
    double J = this->d_c1_basis_eulerian(s, psi_n, psi_i, dpsi_n_dx, dpsi_i_dx);

    // Rotate the degrees of freedom
    rotate_shape(psi_n, dpsi_n_dx);
    // Galerkin
    // (Shallow) copy the basis functions
    test_n = psi_n;
    dtest_n_dx = dpsi_n_dx;
    test_i = psi_i;
    dtest_i_dx = dpsi_i_dx;

    return J;
  }


  //======================================================================
  /// Fetch the basis functions and first/second derivatives
  /// w.r.t. global coordinates and return Jacobian of mapping.
  ///
  /// Galerkin: Test functions = shape functions
  //======================================================================
  double KoiterSteigmannC1CurvableBellElement::
    d2basis_u_eulerian_koiter_steigmann(const Vector<double>& s,
                                        Shape& psi_n,
                                        Shape& psi_i,
                                        DShape& dpsi_n_dx,
                                        DShape& dpsi_i_dx,
                                        DShape& d2psi_n_dx,
                                        DShape& d2psi_i_dx) const
  {
    // Call the geometrical shape functions and derivatives
    double J = CurvableBellElement<Nnode_1D>::d2_c1_basis_eulerian(
      s, psi_n, psi_i, dpsi_n_dx, dpsi_i_dx, d2psi_n_dx, d2psi_i_dx);
    // Rotate the dofs
    rotate_shape(psi_n, dpsi_n_dx, d2psi_n_dx);

    // Return the jacobian
    return J;
  }


  //======================================================================
  /// Fetch the basis functions and test functions and first/second derivatives
  /// w.r.t. global coordinates and return Jacobian of mapping.
  ///
  /// Galerkin: Test functions = shape functions
  //======================================================================
  double KoiterSteigmannC1CurvableBellElement::
    d2basis_and_d2test_u_eulerian_koiter_steigmann(const Vector<double>& s,
                                                   Shape& psi_n,
                                                   Shape& psi_i,
                                                   DShape& dpsi_n_dx,
                                                   DShape& dpsi_i_dx,
                                                   DShape& d2psi_n_dx,
                                                   DShape& d2psi_i_dx,
                                                   Shape& test_n,
                                                   Shape& test_i,
                                                   DShape& dtest_n_dx,
                                                   DShape& dtest_i_dx,
                                                   DShape& d2test_n_dx,
                                                   DShape& d2test_i_dx) const
  {
    // Call the geometrical shape functions and derivatives
    double J = CurvableBellElement<Nnode_1D>::d2_c1_basis_eulerian(
      s, psi_n, psi_i, dpsi_n_dx, dpsi_i_dx, d2psi_n_dx, d2psi_i_dx);
    // Rotate the dofs
    rotate_shape(psi_n, dpsi_n_dx, d2psi_n_dx);

    // Galerkin
    // Set the test functions equal to the shape functions (this is a shallow
    // copy)
    test_n = psi_n;
    dtest_n_dx = dpsi_n_dx;
    d2test_n_dx = d2psi_n_dx;
    test_i = psi_i;
    dtest_i_dx = dpsi_i_dx;
    d2test_i_dx = d2psi_i_dx;

    // Return the jacobian
    return J;
  }


  //======================================================================
  /// Rotate the shape functions according to the specified basis on the
  /// boundary. This function does a DenseDoubleMatrix solve to determine
  /// new basis - which could be speeded up by caching the matrices higher
  /// up and performing the LU decomposition only once
  //======================================================================
  inline void KoiterSteigmannC1CurvableBellElement::rotate_shape(
    Shape& psi) const
  {
    const unsigned n_dof_types = nu_type_at_each_node();

    // Get the nodes that need rotating
    Vector<unsigned> nodes_to_rotate;
    for (unsigned j_node = 0; j_node < 3; j_node++)
    {
      // If the node has had its boundary parametrisation set, its shape
      // functions need rotating
      if (Rotated_boundary_helper_pt->nodal_boundary_parametrisation_pt(j_node))
      {
        nodes_to_rotate.push_back(j_node);
      }
    }

    // Loop over the nodes with rotated dofs
    unsigned n_nodes_to_rotate = nodes_to_rotate.size();
    for (unsigned j = 0; j < n_nodes_to_rotate; j++)
    {
      // Get the nodes
      unsigned j_node = nodes_to_rotate[j];

      // Construct the vectors to hold the shape functions
      Vector<double> psi_vector(n_dof_types);

      // Get the rotation matrix
      DenseDoubleMatrix rotation_matrix(n_dof_types, n_dof_types, 0.0);
      this->Rotated_boundary_helper_pt->get_rotation_matrix_at_node(
        j_node, rotation_matrix);

      // Copy to the vectors
      for (unsigned l = 0; l < n_dof_types; ++l)
      {
        // Copy to the vectors
        for (unsigned k = 0; k < n_dof_types; ++k)
        {
          // Copy over shape functions
          // psi_vector[l]=psi(inode,l);
          psi_vector[l] += psi(j_node, k) * rotation_matrix(l, k);
        }
      }

      // Copy back to shape the rotated vetcors
      for (unsigned l = 0; l < n_dof_types; ++l)
      {
        // Copy over shape functions
        psi(j_node, l) = psi_vector[l];
      }
    }
  }


  //======================================================================
  /// Rotate the shape functions according to the specified basis on the
  /// boundary. This function does a DenseDoubleMatrix solve to determine
  /// new basis - which could be speeded up by caching the matrices higher
  /// up and performing the LU decomposition only once
  //======================================================================
  inline void KoiterSteigmannC1CurvableBellElement::rotate_shape(
    Shape& psi, DShape& dpsidx) const
  {
    const unsigned n_dof_types = nu_type_at_each_node();
    const unsigned n_dim = this->dim();

    // Get the nodes that need rotating
    Vector<unsigned> nodes_to_rotate;
    for (unsigned j_node = 0; j_node < 3; j_node++)
    {
      // If the node has had its boundary parametrisation set, its shape
      // functions need rotating
      if (Rotated_boundary_helper_pt->nodal_boundary_parametrisation_pt(j_node))
      {
        nodes_to_rotate.push_back(j_node);
      }
    }

    // Loop over the nodes with rotated dofs
    unsigned n_nodes_to_rotate = nodes_to_rotate.size();
    for (unsigned j = 0; j < n_nodes_to_rotate; j++)
    {
      // Get the nodes
      unsigned j_node = nodes_to_rotate[j];

      // Construct the vectors to hold the shape functions
      Vector<double> psi_vector(n_dof_types);
      Vector<Vector<double>> dpsi_vector_dxi(n_dim,
                                             Vector<double>(n_dof_types));

      // Get the rotation matrix
      DenseDoubleMatrix rotation_matrix(n_dof_types, n_dof_types, 0.0);
      this->Rotated_boundary_helper_pt->get_rotation_matrix_at_node(
        j_node, rotation_matrix);

      // Copy to the vectors
      for (unsigned l = 0; l < n_dof_types; ++l)
      {
        // Copy to the vectors
        for (unsigned k = 0; k < n_dof_types; ++k)
        {
          // Copy over shape functions
          // psi_vector[l]=psi(inode,l);
          psi_vector[l] += psi(j_node, k) * rotation_matrix(l, k);
          // Copy over first derivatives
          for (unsigned i = 0; i < n_dim; ++i)
          {
            dpsi_vector_dxi[i][l] +=
              dpsidx(j_node, k, i) * rotation_matrix(l, k);
          }
        }
      }

      // Copy back to shape the rotated vetcors
      for (unsigned l = 0; l < n_dof_types; ++l)
      {
        // Copy over shape functions
        psi(j_node, l) = psi_vector[l];
        // Copy over first derivatives
        for (unsigned i = 0; i < n_dim; ++i)
        {
          dpsidx(j_node, l, i) = dpsi_vector_dxi[i][l];
        }
      }
    }
  }


  //======================================================================
  /// Rotate the shape functions according to the specified basis on the
  /// boundary. This function does a DenseDoubleMatrix solve to determine
  /// new basis - which could be speeded up by caching the matrices higher
  /// up and performing the LU decomposition only once
  //======================================================================
  inline void KoiterSteigmannC1CurvableBellElement::rotate_shape(
    Shape& psi, DShape& dpsidx, DShape& d2psidx) const
  {
    const unsigned n_dof_types = nu_type_at_each_node();
    const unsigned n_dim = this->dim();
    // n_dimth triangle number
    const unsigned n_2ndderiv = ((n_dim + 1) * (n_dim)) / 2;

    // Get the nodes that need rotating
    Vector<unsigned> nodes_to_rotate;
    for (unsigned j_node = 0; j_node < 3; j_node++)
    {
      // If the node has had its boundary parametrisation set, its shape
      // functions need rotating
      if (Rotated_boundary_helper_pt->nodal_boundary_parametrisation_pt(j_node))
      {
        nodes_to_rotate.push_back(j_node);
      }
    }

    // Loop over the nodes with rotated dofs
    unsigned n_nodes_to_rotate = nodes_to_rotate.size();
    for (unsigned j = 0; j < n_nodes_to_rotate; j++)
    {
      // Get the nodes
      unsigned j_node = nodes_to_rotate[j];

      // Construct the vectors to hold the shape functions
      Vector<double> psi_vector(n_dof_types);
      Vector<Vector<double>> dpsi_vector_dxi(n_dim,
                                             Vector<double>(n_dof_types));
      Vector<Vector<double>> d2psi_vector_dxidxj(n_2ndderiv,
                                                 Vector<double>(n_dof_types));

      // Get the rotation matrix
      DenseDoubleMatrix rotation_matrix(n_dof_types, n_dof_types, 0.0);
      this->Rotated_boundary_helper_pt->get_rotation_matrix_at_node(
        j_node, rotation_matrix);

      // Copy to the vectors
      for (unsigned l = 0; l < n_dof_types; ++l)
      {
        // Copy to the vectors
        for (unsigned k = 0; k < n_dof_types; ++k)
        {
          // Copy over shape functions
          // psi_vector[l]=psi(inode,l);
          psi_vector[l] += psi(j_node, k) * rotation_matrix(l, k);
          // Copy over first derivatives
          for (unsigned i = 0; i < n_dim; ++i)
          {
            dpsi_vector_dxi[i][l] +=
              dpsidx(j_node, k, i) * rotation_matrix(l, k);
          }
          for (unsigned i = 0; i < n_2ndderiv; ++i)
          {
            d2psi_vector_dxidxj[i][l] +=
              d2psidx(j_node, k, i) * rotation_matrix(l, k);
          }
        }
      }

      // Copy back to shape the rotated vetcors
      for (unsigned l = 0; l < n_dof_types; ++l)
      {
        // Copy over shape functions
        psi(j_node, l) = psi_vector[l];
        // Copy over first derivatives
        for (unsigned i = 0; i < n_dim; ++i)
        {
          dpsidx(j_node, l, i) = dpsi_vector_dxi[i][l];
        }
        // Copy over second derivatives
        for (unsigned i = 0; i < n_2ndderiv; ++i)
        {
          d2psidx(j_node, l, i) = d2psi_vector_dxidxj[i][l];
        }
      }
    }
  }



//   //======================================================================
//   /// Rotate the shape functions according to the specified basis on the
//   /// boundary. This function does a DenseDoubleMatrix solve to determine
//   /// new basis - which could be speeded up by caching the matrices higher
//   /// up and performing the LU decomposition only once
//   //======================================================================
//   // template<unsigned DIM,
//   //          unsigned NNODE_1D,
//   //          unsigned BOUNDARY_ORDER,
//   //          template<unsigned DIM_, unsigned NNODE_1D_>
//   //          class PLATE_EQUATIONS>
//   inline void KoiterSteigmannC1CurvableBellElement<
//     DIM,
//     NNODE_1D,
//     BOUNDARY_ORDER,
//     PLATE_EQUATIONS>::rotate_shape(Shape& psi) const
//   {
//     // Loop over the nodes with rotated dofs
//     for (unsigned i = 0; i < Nnodes_to_rotate; ++i)
//     {
//       // Get the nodes
//       unsigned inode = Nodes_to_rotate[i];

//       // Construct the vectors to hold the shape functions
//       Vector<double> psi_vector(6, 0);

//       // Get the rotation matrix
//       DenseDoubleMatrix rotation_matrix(6, 6, 0.0);
//       this->rotation_matrix_at_node(inode, rotation_matrix);

//       // Copy to the vectors
//       for (unsigned l = 0; l < 6; ++l)
//         for (unsigned k = 0; k < 6; ++k)
//         {
//           psi_vector[l] += psi(inode, k) * rotation_matrix(l, k);
//         }

//       // Copy back to shape the rotated vetcors
//       for (unsigned l = 0; l < 6; ++l)
//       {
//         psi(inode, l) = psi_vector[l];
//       }
//     }
//   }

//   //======================================================================
//   /// Rotate the shape functions according to the specified basis on the
//   /// boundary. This function does a DenseDoubleMatrix solve to determine
//   /// new basis - which could be speeded up by caching the matrices higher
//   /// up and performing the LU decomposition only once
//   //======================================================================
//   // template<unsigned DIM,
//   //          unsigned NNODE_1D,
//   //          unsigned BOUNDARY_ORDER,
//   //          template<unsigned DIM_, unsigned NNODE_1D_>
//   //          class PLATE_EQUATIONS>
//   inline void KoiterSteigmannC1CurvableBellElement<
//     DIM,
//     NNODE_1D,
//     BOUNDARY_ORDER,
//     PLATE_EQUATIONS>::rotate_shape(Shape& psi, DShape& dpsidx) const
//   {
//     // Loop over the nodes with rotated dofs
//     for (unsigned i = 0; i < Nnodes_to_rotate; ++i)
//     {
//       // Get the nodes
//       unsigned inode = Nodes_to_rotate[i];

//       // Construct the vectors to hold the shape functions
//       Vector<double> psi_vector(6, 0);
//       Vector<Vector<double>> dpsi_vector_dxi(2, Vector<double>(6, 0));

//       // Get the rotation matrix
//       DenseDoubleMatrix rotation_matrix(6, 6, 0.0);
//       this->rotation_matrix_at_node(inode, rotation_matrix);

//       // Copy to the vectors
//       for (unsigned l = 0; l < 6; ++l)
//       {
//         // Copy to the vectors
//         for (unsigned k = 0; k < 6; ++k)
//         {
//           // Copy over shape functions
//           // psi_vector[l]=psi(inode,l);
//           psi_vector[l] += psi(inode, k) * rotation_matrix(l, k);
//           // Copy over first derivatives
//           for (unsigned i = 0; i < 2; ++i)
//           {
//             dpsi_vector_dxi[i][l] +=
//               dpsidx(inode, k, i) * rotation_matrix(l, k);
//           }
//         }
//       }

//       // Copy back to shape the rotated vetcors
//       for (unsigned l = 0; l < 6; ++l)
//       {
//         // Copy over shape functions
//         psi(inode, l) = psi_vector[l];
//         // Copy over first derivatives
//         for (unsigned i = 0; i < 2; ++i)
//         {
//           dpsidx(inode, l, i) = dpsi_vector_dxi[i][l];
//         }
//       }
//     }
//   }

//   //======================================================================
//   /// Rotate the shape functions according to the specified basis on the
//   /// boundary. This function does a DenseDoubleMatrix solve to determine
//   /// new basis - which could be speeded up by caching the matrices higher
//   /// up and performing the LU decomposition only once
//   //======================================================================
//   // template<unsigned DIM,
//   //          unsigned NNODE_1D,
//   //          unsigned BOUNDARY_ORDER,
//   //          template<unsigned DIM_, unsigned NNODE_1D_>
//   //          class PLATE_EQUATIONS>
//   inline void KoiterSteigmannC1CurvableBellElement<
//     DIM,
//     NNODE_1D,
//     BOUNDARY_ORDER,
//     PLATE_EQUATIONS>::rotate_shape(Shape& psi,
//                                    DShape& dpsidx,
//                                    DShape& d2psidx) const
//   {
//     // Loop over the nodes with rotated dofs
//     for (unsigned i = 0; i < Nnodes_to_rotate; ++i)
//     {
//       // Get the nodes
//       unsigned inode = Nodes_to_rotate[i];

//       // Construct the vectors to hold the shape functions
//       Vector<double> psi_vector(6, 0);
//       Vector<Vector<double>> dpsi_vector_dxi(2, Vector<double>(6, 0));
//       Vector<Vector<double>> d2psi_vector_dxidxj(3, Vector<double>(6, 0));

//       // Get the rotation matrix
//       DenseDoubleMatrix rotation_matrix(6, 6, 0.0);
//       this->rotation_matrix_at_node(inode, rotation_matrix);

//       // Copy to the vectors
//       for (unsigned l = 0; l < 6; ++l)
//       {
//         // Copy to the vectors
//         for (unsigned k = 0; k < 6; ++k)
//         {
//           // Copy over shape functions
//           // psi_vector[l]=psi(inode,l);
//           psi_vector[l] += psi(inode, k) * rotation_matrix(l, k);
//           // Copy over first derivatives
//           for (unsigned i = 0; i < 2; ++i)
//           {
//             dpsi_vector_dxi[i][l] +=
//               dpsidx(inode, k, i) * rotation_matrix(l, k);
//           }
//           for (unsigned i = 0; i < 3; ++i)
//           {
//             d2psi_vector_dxidxj[i][l] +=
//               d2psidx(inode, k, i) * rotation_matrix(l, k);
//           }
//         }
//       }

//       // Copy back to shape the rotated vetcors
//       for (unsigned l = 0; l < 6; ++l)
//       {
//         // Copy over shape functions
//         psi(inode, l) = psi_vector[l];
//         // Copy over first derivatives
//         for (unsigned i = 0; i < 2; ++i)
//         {
//           dpsidx(inode, l, i) = dpsi_vector_dxi[i][l];
//         }
//         // Copy over second derivatives
//         for (unsigned i = 0; i < 3; ++i)
//         {
//           d2psidx(inode, l, i) = d2psi_vector_dxidxj[i][l];
//         }
//       }
//     }
//   }

//   //==============================================================================
//   /// Get the jth bubble dof at the lth internal point. Deliberately broken for
//   /// the case where there is no curved edge
//   //==============================================================================
//   // template<unsigned DIM,
//   //          unsigned NNODE_1D,
//   //          unsigned BOUNDARY_ORDER,
//   //          template<unsigned DIM_, unsigned NNODE_1D_>
//   //          class PLATE_EQUATIONS>
//   double KoiterSteigmannC1CurvableBellElement<
//     DIM,
//     NNODE_1D,
//     BOUNDARY_ORDER,
//     PLATE_EQUATIONS>::get_u_bubble_dof(const unsigned& l,
//                                        const unsigned& j) const
//   {
//     // Now give the lth internal degree of freedom
//     if (Curved_edge != MyC1CurvedElements::none && j <= 2)
//     {
//       // Internal dofs per displacement
//       unsigned n_internal_dofs = Curved_shape.n_internal_dofs();
//       return this->internal_data_pt(Bubble_u_internal_index)
//         ->value(l + j * n_internal_dofs);
//     }
//     // Deliberately break this function for the below cases
//     // If there is no curved edge then we cannot return anything meaningful
//     else if (Curved_edge == MyC1CurvedElements::none)
//     {
//       throw OomphLibError(
//         "There are no time-dependent internal 'bubble' dofs for \
// this element.",
//         OOMPH_CURRENT_FUNCTION,
//         OOMPH_EXCEPTION_LOCATION);
//       // Return dummy value 0.0
//       return 0;
//     }
//     // For these elements we only have a single dof at each internal point
//     else // if (j>2)
//     {
//       throw OomphLibError(
//         "There are only three degrees of freedom at the internal points in this \
// element.",
//         OOMPH_CURRENT_FUNCTION,
//         OOMPH_EXCEPTION_LOCATION);
//       // Return dummy value 0.0
//       return 0;
//     }
//   }

//   //==============================================================================
//   /// Get the jth bubble dof at the lth internal point. Deliberately broken for
//   /// case when there is no curved edge.
  //   //==============================================================================
  //   // template<unsigned DIM,
  //   //          unsigned NNODE_1D,
  //   //          unsigned BOUNDARY_ORDER,
  //   //          template<unsigned DIM_, unsigned NNODE_1D_>
  //   //          class PLATE_EQUATIONS>
  //   int KoiterSteigmannC1CurvableBellElement<
  //     DIM,
  //     NNODE_1D,
  //     BOUNDARY_ORDER,
  //     PLATE_EQUATIONS>::local_u_bubble_equation(const unsigned& l,
  //                                               const unsigned& j) const
  //   {
  //     // Deliberately break this function for the below cases
  //     // If there is no curved edge then we cannot return anything meaningful
  //     if (Curved_edge == MyC1CurvedElements::none)
  //     {
  //       throw OomphLibError(
  //         "There are no time-dependent internal 'bubble' dofs for this
  //         element.", OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  //       // Return dummy value -2
  //       return -2;
  //     }
  //     // For these elements we only have a single dof at each internal point
  //     else if (j > 2)
  //     {
  //       throw OomphLibError(
  //         "There are only three equations at the internal points in this \
// element.",
  //         OOMPH_CURRENT_FUNCTION,
  //         OOMPH_EXCEPTION_LOCATION);
  //       // Return dummy value -2
  //       return -2;
  //     }
  //     // Now give the lth internal equation number
  //     else
  //     {
  //       // Internal dofs per displacement
  //       unsigned n_internal_dofs = Curved_shape.n_internal_dofs();
  //       return this->internal_local_eqn(Bubble_u_internal_index,
  //                                       l + n_internal_dofs * j);
  //     }
  //   }

  //   //==============================================================================
  //   /// Set up the rotated degrees of freedom: includes a check for the
  //   number of
  //   /// rotation nodes being greater than three.
  //   //==============================================================================
  //   // template<unsigned DIM,
  //   //          unsigned NNODE_1D,
  //   //          unsigned BOUNDARY_ORDER,
  //   //          template<unsigned DIM_, unsigned NNODE_1D_>
  //   //          class PLATE_EQUATIONS>
  //   void KoiterSteigmannC1CurvableBellElement<DIM,
  //                                                  NNODE_1D,
  //                                                  BOUNDARY_ORDER,
  //                                                  PLATE_EQUATIONS>::
  //     set_up_rotated_dofs(const unsigned& nnodes_to_rotate,
  //                         const Vector<unsigned>& nodes_to_rotate,
  //                         const BasisVectorsFctPt& basis_vectors_fct_pt)
  //   {
  //     // Change the member Nnode_to_rotate
  //     Nnodes_to_rotate = nnodes_to_rotate;
  // #ifdef PARANOID
  //     // Check that the number of nodes is smaller than 3
  //     if (nnodes_to_rotate > 3)
  //     {
  //       throw OomphLibError(
  //         "There are only three nodes per element, so we cannot rotate more
  //         than\
//  three ",
  //         OOMPH_CURRENT_FUNCTION,
  //         OOMPH_EXCEPTION_LOCATION);
  //     }
  // #endif

  //     Nodes_to_rotate = nodes_to_rotate;

  //     // Point to the basis vectors function
  //     Rotated_basis_fct_pt = basis_vectors_fct_pt;
  //   }

  //   //==============================================================================
  //   /// Rotate the shape functions according to
  //   /// w.r.t. global coordinates and return Jacobian of mapping.
  //   ///
  //   /// Galerkin: Test functions = shape functions
  //   //==============================================================================
  //   // template<unsigned DIM,
  //   //          unsigned NNODE_1D,
  //   //          unsigned BOUNDARY_ORDER,
  //   //          template<unsigned DIM_, unsigned NNODE_1D_>
  //   //          class PLATE_EQUATIONS>
  //   void KoiterSteigmannC1CurvableBellElement<DIM,
  //                                                  NNODE_1D,
  //                                                  BOUNDARY_ORDER,
  //                                                  PLATE_EQUATIONS>::
  //     rotation_matrix_at_node(const unsigned& inode,
  //                             DenseDoubleMatrix& rotation_matrix) const
  //   {
  //     // Initialise x normal and tangent
  //     Vector<double> x(2, 0.0);

  //     // Get the node pointer
  //     Node* nod_pt = this->node_pt(inode);

  //     // Get the position of the vertex
  //     x[0] = nod_pt->x(0);
  //     x[1] = nod_pt->x(1);

  //     // Initialise the two basis vectors
  //     Vector<Vector<double>> bi(2, Vector<double>(2, 0.0));
  //     Vector<DenseMatrix<double>> Dbi(2, DenseMatrix<double>(2, 2, 0.0));

  //     // Now find the two new basis vectors
  //     // Get the normal - overload if not continuous
  //     (*Rotated_basis_fct_pt)(x, bi[0], bi[1], Dbi[0], Dbi[1]);

  //     // Rotation matrix, B
  //     DenseMatrix<double> b1(2, 2, 0.0), b22(3, 3, 0.0), b21(3, 2, 0.0);

  //     // Fill in the submatrices
  //     for (unsigned alpha = 0; alpha < 2; ++alpha)
  //     {
  //       for (unsigned beta = 0; beta < 2; ++beta)
  //       {
  //         // Fill in b1 - the Jacobian
  //         // Fill in the rotation of the first derivatives
  //         b1(alpha, beta) = bi[beta][alpha];

  //         // Avoid double counting the cross derivative
  //         if (alpha <= beta)
  //         {
  //           // Define row index
  //           const unsigned row = alpha + beta;
  //           for (unsigned gamma = 0; gamma < 2; ++gamma)
  //           {
  //             // Fill in b21 - the non affine part of the Jacobian derivative
  //             // Define column index
  //             unsigned col_b21 = gamma;
  //             // Fill in the non-affine part of the rotation of the second
  //             // derivatives if( beta>= alpha) ?
  //             b21(row, col_b21) += Dbi[gamma](alpha, beta);
  //             for (unsigned delta = 0; delta < 2; ++delta)
  //             {
  //               // Fill in b22 - the Affine part of the Jacobian derivative
  //               // Redefine column index for the next submatrix
  //               unsigned col_b22 = gamma + delta;
  //               // Fill in the affine part of the rotation of the second
  //               // derivatives if( beta>= alpha) ?
  //               b22(row, col_b22) += bi[gamma][alpha] * bi[delta][beta];
  //             }
  //           }
  //         }
  //       }
  //     }

  //     // Fill in the submatrices to the full (6x6) matrix - we need to right
  //     // multiply this matrix so we need the transpose of the Jacobian W dof
  //     // remains the same
  //     rotation_matrix(0, 0) = 1.0;
  //     // Fill in b1
  //     for (unsigned i = 0; i < 2; ++i)
  //     {
  //       for (unsigned j = 0; j < 2; ++j)
  //       {
  //         rotation_matrix(1 + j, 1 + i) = b1(i, j);
  //       }
  //     }
  //     // Fill in b21
  //     for (unsigned i = 0; i < 3; ++i)
  //     {
  //       for (unsigned j = 0; j < 2; ++j)
  //       {
  //         rotation_matrix(1 + j, 3 + i) = b21(i, j);
  //       }
  //     }
  //     // Fill in b22
  //     for (unsigned i = 0; i < 3; ++i)
  //     {
  //       for (unsigned j = 0; j < 3; ++j)
  //       {
  //         rotation_matrix(3 + j, 3 + i) = b22(i, j);
  //       }
  //     }
  //   }
  // /// \short get the coordinate i
  // inline void interpolated_zeta(const Vector<double>& s,
  //                               Vector<double>& zeta) const
  // {
  //   get_coordinate_x(s, zeta);
  // }

  // /// \short get the coordinate i
  //  void interpolated_x (const Vector<double>& s, Vector<double>& r) const
  //  { get_coordinate_x(s,r);}

  // inline void get_internal_dofs_location(const unsigned& s,
  //                                        Vector<double>& x) const;

  // /// Access to number of internal dofs
  // unsigned number_of_internal_dofs() const
  // {
  //   return this->Number_of_internal_dofs;
  // }

  //  //=======================================================================
  //  /// Return FE interpolated position x[] at local coordinate s as Vector
  //  //=======================================================================
  //   void interpolated_x(const Vector<double> &s, Vector<double> &x)
  //    const
  //   {
  //     //Find the number of nodes
  //    const unsigned n_node = this->nnode();
  //    //Find the number of positional types
  //    const unsigned n_position_type = 1;
  //    //Find the dimension stored in the node
  //    const unsigned nodal_dim = this->nodal_dimension();
  //
  //    //Assign storage for the local shape function
  //    Shape psi(n_node,n_position_type);
  //    //Find the values of shape function
  //    this->shape(s,psi);
  //
  //    //Loop over the dimensions
  //    for(unsigned i=0;i<nodal_dim;i++)
  //     {
  //      //Initilialise value of x[i] to zero
  //      x[i] = 0.0;
  //      //Loop over the local nodes
  //      for(unsigned l=0;l<n_node;l++)
  //       {
  //        //Loop over the number of dofs
  //        for(unsigned k=0;k<n_position_type;k++)
  //         {
  //          x[i] += this->nodal_position_gen(l,k,i)*psi(l,k);
  //         }
  //       }
  //     }
  //   }

  // /// \short Set up the rotated degrees of freedom
  // inline void set_up_rotated_dofs(
  //   const unsigned& nnodes_to_rotate,
  //   const Vector<unsigned>& nodes_to_rotate,
  //   const BasisVectorsFctPt& basis_vectors_fct_pt);

  // // HERE wrapper around locate zeta - hacky way to get the interface working
  // // needs FIXING
  // void locate_zeta(const Vector<double>& zeta,
  //                  GeomObject*& geom_object_pt,
  //                  Vector<double>& s,
  //                  const bool& use_coordinate_as_initial_guess)
  // {
  //   // Temporarily set nnodal_position_type to be one
  //   this->set_nnodal_position_type(1);
  //   FiniteElement::locate_zeta(
  //     zeta, geom_object_pt, s, use_coordinate_as_initial_guess);
  //   // Set it back to six
  //   this->set_nnodal_position_type(6);
  // }

  // // Upgrade an element to its curved counterpart
  // inline void upgrade_to_curved_element(const Edge& curved_edge,
  //                                       const double& s_ubar,
  //                                       const double& s_obar,
  //                                       CurvilineGeomObject*
  //                                       parametric_edge);

  // // Precompute the association matrix
  // void precompute_association_matrix(DenseMatrix<double>& m)
  // {
  //   // If the element has been upgraded
  //   if (Curved_edge == MyC1CurvedElements::none)
  //   {
  //   } // Do nothing
  //   else
  //   {
  //     Curved_shape.fill_in_full_association_matrix(m);
  //   }
  // };

  // // Get the number of basis functions, wrapper
  // double n_basis_functions()
  // {
  //   return Curved_shape.n_basis_functions();
  // };

  // // Get the number of basic basis functions, wrapper
  // double n_basic_basis_functions()
  // {
  //   return Curved_shape.n_basic_basis_functions();
  // };

  // /// \short Shape, test functions & derivs. w.r.t. to global coords. at
  // /// integration point ipt. Return Jacobian.
  // inline double d2shape_and_d2test_eulerian_at_knot_biharmonic(const
  // unsigned& ipt,
  //                                                         Shape &psi,
  //                                                         DShape &dpsidx,
  //                                                         DShape &d2psidx,
  //                                                         Shape &test,
  //                                                         DShape &dtestdx,
  //                                                         DShape &d2testdx)
  //  const;

  // inline double dshape_and_dtest_eulerian_at_knot_biharmonic(const unsigned
  // &ipt,
  //                                                         Shape &psi,
  //                                                         DShape &dpsidx,
  //                                                         Shape &test,
  //                                                         DShape &dtestdx)
  //  const;

  // // Inline functions:
  // //==============================================================================
  // /// Get the mapped position in the element. For straight sided elements
  // this
  // /// is and affine mapping.
  // //==============================================================================
  // // template<unsigned DIM,
  // //          unsigned NNODE_1D,
  // //          unsigned BOUNDARY_ORDER,
  // //          template<unsigned DIM_, unsigned NNODE_1D_>
  // //          class PLATE_EQUATIONS>
  // void KoiterSteigmannC1CurvableBellElement::get_coordinate_x(const
  // Vector<double>& s,
  //                                      Vector<double>& x) const
  // {
  //   // If the element has been upgraded
  //   if (Curved_edge == MyC1CurvedElements::none)
  //   {
  //     this->my_interpolated_x(s, x);
  //   }
  //   else
  //   {
  //     Curved_shape.coordinate_x(s, x);
  //   }
  // };

  //   // Inline functions:
  //   //==============================================================================
  //   /// Get the mapped position in the element. For straight sided elements
  //   this
  //   /// is and affine mapping.
  //   //==============================================================================
  //   // template<unsigned DIM,
  //   //          unsigned NNODE_1D,
  //   //          unsigned BOUNDARY_ORDER,
  //   //          template<unsigned DIM_, unsigned NNODE_1D_>
  //   //          class PLATE_EQUATIONS>
  //   void KoiterSteigmannC1CurvableBellElement<
  //     DIM,
  //     NNODE_1D,
  //     BOUNDARY_ORDER,
  //     PLATE_EQUATIONS>::get_internal_dofs_location(const unsigned& dof,
  //                                                  Vector<double>& s) const
  //   {
  //     // If the element has been upgraded
  //     if (Curved_edge == MyC1CurvedElements::none)
  //     {
  //       throw OomphLibError(
  //         "There are no internal dofs for these elements as they have not
  //         been\
// upgraded to curved elements.",
  //         OOMPH_CURRENT_FUNCTION,
  //         OOMPH_EXCEPTION_LOCATION);
  //     }
  //     else
  //     {
  //       Curved_shape.get_internal_dofs_location(dof, s);
  //     }
  //   };

  //   //==============================================================================
  //   /// Upgrade an element to its curved counterpart: this adds internal data
  //   to
  //   /// elements and upgrades the shape class data member.
  //   //==============================================================================
  //   // template<unsigned DIM,
  //   //          unsigned NNODE_1D,
  //   //          unsigned BOUNDARY_ORDER,
  //   //          template<unsigned DIM_, unsigned NNODE_1D_>
  //   //          class PLATE_EQUATIONS>
  //   inline void KoiterSteigmannC1CurvableBellElement<DIM,
  //                                                         NNODE_1D,
  //                                                         BOUNDARY_ORDER,
  //                                                         PLATE_EQUATIONS>::
  //     upgrade_to_curved_element(const Edge& curved_edge,
  //                               const double& s_ubar,
  //                               const double& s_obar,
  //                               CurvilineGeomObject* parametric_edge)
  //   {
  // #ifdef PARANOID
  //     // When upgrading add to count
  //     Curved_edge_counter += 1;
  //     // Check that we haven't upgraded this element already
  //     if (Curved_edge_counter > 1)
  //     {
  //       // SCREAM
  //       throw OomphLibError(
  //         "Cannot upgrade more than a single edge to be curved in C1 Curved
  //         Bell \
// Elements.",
  //         OOMPH_CURRENT_FUNCTION,
  //         OOMPH_EXCEPTION_LOCATION);
  //     }
  // #endif
  //     using namespace MyC1CurvedElements;
  //     // Add the curved edge
  //     Curved_edge = curved_edge;

  //     // Set the integral pointer
  //     // HERE do we need order 16 accuracy rather than 15?
  //     // Do we need specific integration scheme for the p16 term u_i,xj
  //     psi_xj
  //     // terms?
  //     delete this->integral_pt();
  //     if (BOUNDARY_ORDER == 5)
  //     {
  //       TGauss<2, 16>* new_integral_pt = new TGauss<2, 16>;
  //       // Set the Integration scheme
  //       this->set_integration_scheme(new_integral_pt);
  //       // TGauss<2,5>* new_integral_pt = new TGauss<2,5>;
  //     }
  //     else
  //     {
  //       TGauss<2, 13>* new_integral_pt = new TGauss<2, 13>;
  //       // Set the Integration scheme
  //       this->set_integration_scheme(new_integral_pt);
  //     }
  //     // Set the number of internal dofs to 1
  //     this->Number_of_internal_dof_types = 1;
  //     this->Number_of_internal_dofs = Curved_shape.n_internal_dofs();
  //     Bubble_u_internal_index = this->add_internal_data(new Data(
  //       (this->Number_of_internal_dofs) * (this->Number_of_displacements)));

  //     // Set up the data of the element
  //     typename BernadouElementBasis<BOUNDARY_ORDER>::VertexList vertices(
  //       3, Vector<double>(2, 0.0));

  //     // Now switch to upgrade
  //     // The shape functions are designed such that the curved edge is always
  //     edge
  //     // two. So this is where we set that up. This is temporary and not the
  //     final
  //     // solution we want
  //     switch (curved_edge)
  //     {
  //       // Throw an error if an edge is upgraded to none
  //       case none:
  //         throw OomphLibError(
  //           "Cannot upgrade edge 'none'. Curved elements must have\
// one side defined by a parametric function.",
  //           OOMPH_CURRENT_FUNCTION,
  //           OOMPH_EXCEPTION_LOCATION);
  //         break;
  //       case zero:
  //         // Everything cyclicly permutes
  //         for (unsigned i = 0; i < 2; ++i)
  //         {
  //           vertices[2][i] = this->node_pt(0)->x(i);
  //           vertices[0][i] = this->node_pt(1)->x(i);
  //           vertices[1][i] = this->node_pt(2)->x(i);
  //         }
  //         break;
  //       case one:
  //         // Everything cyclicly permutes
  //         for (unsigned i = 0; i < 2; ++i)
  //         {
  //           vertices[2][i] = this->node_pt(1)->x(i);
  //           vertices[0][i] = this->node_pt(2)->x(i);
  //           vertices[1][i] = this->node_pt(0)->x(i);
  //         }
  //         break;
  //       case two:
  //         // Everything is just copied over
  //         for (unsigned i = 0; i < 2; ++i)
  //         {
  //           vertices[2][i] = this->node_pt(2)->x(i);
  //           vertices[0][i] = this->node_pt(0)->x(i);
  //           vertices[1][i] = this->node_pt(1)->x(i);
  //         }
  //         break;
  //     }

  //     // Add the vertices to make the shape functions fully functional
  //     Curved_shape.upgrade_element(
  //       vertices, s_ubar, s_obar, curved_edge, *parametric_edge);
  //   }

  //   //======================================================================
  //   /// Define the shape functions and test functions and derivatives
  //   /// w.r.t. global coordinates and return Jacobian of mapping.
  //   ///
  //   /// Galerkin: Test functions = shape functions
  //   //======================================================================
  //   // template<unsigned DIM,
  //   //          unsigned NNODE_1D,
  //   //          unsigned BOUNDARY_ORDER,
  //   //          template<unsigned DIM_, unsigned NNODE_1D_>
  //   //          class PLATE_EQUATIONS>
  //   void KoiterSteigmannC1CurvableBellElement<
  //     DIM,
  //     NNODE_1D,
  //     BOUNDARY_ORDER,
  //     PLATE_EQUATIONS>::shape_and_test_biharmonic(const Vector<double>& s,
  //                                                 Shape& psi,
  //                                                 Shape& psi_b,
  //                                                 Shape& test,
  //                                                 Shape& test_b) const
  //   {
  //     throw OomphLibError(
  //       "This still needs testing for curved elements.",
  //       "void
  //       KoiterSteigmannC1CurvableBellElement<DIM,NNODE_1D,BOUNDARY_ORDER,PLATE_EQUATIONS>::\
// shape_and_test_biharmonic(...)",
  //       OOMPH_EXCEPTION_LOCATION); // HERE

  //     // Get dummy shape functions for the Bell call
  //     DShape dpsidx(3, 6, 2);
  //     DShape d2psidx(3, 6, 3);

  //     // Vertices
  //     Vector<Vector<double>> v(3, Vector<double>(2));
  //     for (unsigned inode = 0; inode < 3; ++inode)
  //     {
  //       // Get the position vector
  //       Node* nod_pt = this->node_pt(inode);
  //       v[inode][0] = nod_pt->x(0);
  //       v[inode][1] = nod_pt->x(1);
  //     }

  //     // If the element has not been upgraded
  //     if (Curved_edge == MyC1CurvedElements::none)
  //     {
  //       // Get J
  //       this->J_eulerian1(s);
  //       Bell_basis.d2_basis_eulerian(s, v, psi, dpsidx, d2psidx);
  //     }
  //     else // i.e if has curved edge
  //     {
  //       Curved_shape.shape(s, psi, psi_b);
  //     }

  //     // Rotate the degrees of freedom
  //     rotate_shape(psi);

  //     // Galerkin
  //     // (Shallow) copy the basis functions
  //     test = psi;
  //     test_b = psi_b;
  //   }


  //   //======================================================================
  //   /// Define the shape functions and test functions and derivatives
  //   /// w.r.t. global coordinates and return Jacobian of mapping.
  //   ///
  //   /// Galerkin: Test functions = shape functions
  //   //======================================================================
  //   // template<unsigned DIM,
  //   //          unsigned NNODE_1D,
  //   //          unsigned BOUNDARY_ORDER,
  //   //          template<unsigned DIM_, unsigned NNODE_1D_>
  //   //          class PLATE_EQUATIONS>
  //   double KoiterSteigmannC1CurvableBellElement<DIM,
  //                                                    NNODE_1D,
  //                                                    BOUNDARY_ORDER,
  //                                                    PLATE_EQUATIONS>::
  //     dshape_and_dtest_eulerian_biharmonic(const Vector<double>& s,
  //                                          Shape& psi,
  //                                          Shape& psi_b,
  //                                          DShape& dpsidx,
  //                                          DShape& dpsi_b_dx,
  //                                          Shape& test,
  //                                          Shape& test_b,
  //                                          DShape& dtestdx,
  //                                          DShape& dtest_b_dx) const
  //   {
  //     // Throw if called
  //     throw OomphLibError(
  //       "This still needs testing for curved elements.",
  //       "void
  //       KoiterSteigmannC1CurvableBellElement<DIM,NNODE_1D,BOUNDARY_ORDER,PLATE_EQUATIONS>::\
// dshape_and_dtest_biharmonic(...)",
  //       OOMPH_EXCEPTION_LOCATION); // HERE

  //     // Now set up dummy DShape so we can call Bell
  //     DShape d2psidx(3, 6, 3);
  //     double J = this->J_eulerian1(s);

  //     // Vertices
  //     Vector<Vector<double>> v(3, Vector<double>(2));
  //     for (unsigned inode = 0; inode < 3; ++inode)
  //     {
  //       // Get the position vector
  //       Node* nod_pt = this->node_pt(inode);
  //       v[inode][0] = nod_pt->x(0);
  //       v[inode][1] = nod_pt->x(1);
  //     }

  //     // If the element has been upgraded
  //     if (Curved_edge == MyC1CurvedElements::none)
  //     {
  //       // Get J
  //       J = this->J_eulerian1(s);
  //       Bell_basis.d2_basis_eulerian(s, v, psi, dpsidx, d2psidx);
  //     }
  //     else // i.e if has curved edge
  //     {
  //       J = Curved_shape.d_shape_dx(s, psi, psi_b, dpsidx, dpsi_b_dx);
  //     }

  //     // Rotate the degrees of freedom
  //     rotate_shape(psi, dpsidx);
  //     // Galerkin
  //     // (Shallow) copy the basis functions
  //     test = psi;
  //     dtestdx = dpsidx;
  //     test_b = psi_b;
  //     dtest_b_dx = dpsi_b_dx;

  //     return J;
  //   }

  //   // template<unsigned DIM,
  //   //          unsigned NNODE_1D,
  //   //          unsigned BOUNDARY_ORDER,
  //   //          template<unsigned DIM_, unsigned NNODE_1D_>
  //   //          class PLATE_EQUATIONS>
  //   double KoiterSteigmannC1CurvableBellElement<DIM,
  //                                                    NNODE_1D,
  //                                                    BOUNDARY_ORDER,
  //                                                    PLATE_EQUATIONS>::
  //     d2shape_and_d2test_eulerian_biharmonic(const Vector<double>& s,
  //                                            Shape& psi,
  //                                            Shape& psi_b,
  //                                            DShape& dpsidx,
  //                                            DShape& dpsi_bdx,
  //                                            DShape& d2psidx,
  //                                            DShape& d2psi_bdx,
  //                                            Shape& test,
  //                                            Shape& test_b,
  //                                            DShape& dtestdx,
  //                                            DShape& dtest_bdx,
  //                                            DShape& d2testdx,
  //                                            DShape& d2test_bdx) const
  //   {
  //     // Call the geometrical shape functions and derivatives
  //     double J = this->J_eulerian1(s);
  //     // Vertices
  //     Vector<Vector<double>> v(3, Vector<double>(2));
  //     for (unsigned inode = 0; inode < 3; ++inode)
  //     {
  //       // Get the position vector
  //       Node* nod_pt = this->node_pt(inode);
  //       v[inode][0] = nod_pt->x(0);
  //       v[inode][1] = nod_pt->x(1);
  //     }

  //     // If the element has been upgraded
  //     if (Curved_edge == MyC1CurvedElements::none)
  //     {
  //       // Get J
  //       J = this->J_eulerian1(s);
  //       Bell_basis.d2_basis_eulerian(s, v, psi, dpsidx, d2psidx);
  //     }
  //     // i.e if has curved edge and precomputed matrix
  //     else if (this->get_association_matrix_pt() != 0)
  //     {
  //       J = Curved_shape.d2_shape_dx2(s,
  //                                     psi,
  //                                     psi_b,
  //                                     dpsidx,
  //                                     dpsi_bdx,
  //                                     d2psidx,
  //                                     d2psi_bdx,
  //                                     *(this->get_association_matrix_pt()));
  //     }
  //     // i.e if has curved edge but no precomputed matrix
  //     else
  //     {
  //       J = Curved_shape.d2_shape_dx2(
  //         s, psi, psi_b, dpsidx, dpsi_bdx, d2psidx, d2psi_bdx);
  //     }

  //     // Rotate the dofs
  //     rotate_shape(psi, dpsidx, d2psidx);

  //     // Galerkin
  //     // Set the test functions equal to the shape functions (this is a
  //     shallow
  //     // copy)
  //     test = psi;
  //     dtestdx = dpsidx;
  //     d2testdx = d2psidx;
  //     test_b = psi_b;
  //     dtest_bdx = dpsi_bdx;
  //     d2test_bdx = d2psi_bdx;

  //     // Return the jacobian
  //     return J;
  //   }
  ////==============================================================================
  ///// Define the shape functions and test functions and derivatives
  ///// w.r.t. global coordinates and return Jacobian of mapping.
  /////
  ///// Galerkin: Test functions = shape functions
  ////==============================================================================
  // template<unsigned DIM, unsigned NNODE_1D, unsigned BOUNDARY_ORDER>
  // double
  // KoiterSteigmannC1CurvableBellElement<DIM,NNODE_1D,BOUNDARY_ORDER>::
  //  dshape_and_dtest_eulerian_at_knot_biharmonic(
  //   const unsigned &ipt,
  //   Shape &psi,
  //   DShape &dpsidx,
  //   Shape &test,
  //   DShape &dtestdx) const
  //{
  //  const double J =
  //  this->dshape_and_dtest_eulerian_at_knot(ipt,psi,dpsidx,test
  //   ,dtestdx);
  //
  //  return J;
  // }
  //
  // template<unsigned DIM, unsigned NNODE_1D, unsigned BOUNDARY_ORDER>
  // double
  // KoiterSteigmannC1CurvableBellElement<DIM,NNODE_1D,BOUNDARY_ORDER>::
  //  d2shape_and_d2test_eulerian_at_knot_biharmonic(
  //   const unsigned &ipt,
  //   Shape &psi,
  //   DShape &dpsidx,
  //   DShape &d2psidx,
  //   Shape &test,
  //   DShape &dtestdx,
  //   DShape &d2testdx) const
  //{
  //  //Call the geometrical shape functions and derivatives
  //  const double J = this->d2shape_and_d2test_eulerian_at_knot(ipt,psi,dpsidx,
  //   d2psidx,test,dtestdx,d2testdx);
  //
  //  //Return the jacobian
  //  return J;
  // }

  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////


  //========= start_of_duplicate_node_constraint_element ==================
  /// Non-geometric element used to constrain dofs between duplicated
  /// vertices where the Hemite data at each node is different but must
  /// be compatible.
  ///
  /// If the first (left) node uses coordinates (s_1,s_2) for the fields
  /// (U,V,W) and the second (right) uses coordinates (t_1, t_2) for the fields
  /// (u,v,w) then enforcing (U,V,W)=(u,v,w), using the chain rule we arrive at
  /// three equations for displacement (alpha=1,2):
  ///     0 = (U_\alpha - u_\alpha)
  ///     0 = (W-w)
  /// two equations constraining gradient (alpha=1,2):
  ///     0 = (dW_1/ds_\alpha - dw_2/dt_\beta J_{\beta\alpha})
  /// and three equations constraining curvature (alpha,beta=1,2; beta>=alpha):
  ///     0 = (d^2W_1/ds_\alpha ds_\beta
  ///          - J_{\alpha\gamma} * J_{\beta\delta} * d^2w_2/dt_\gamma dt_\delta
  ///          - H_{\gamma\alpha\beta} * dw_2/dt_gamma)
  /// where L_i, i=0,..,7, are Lagrange multipliers -- dofs which are
  /// stored in the internal data of this element.
  //=======================================================================
  class DuplicateNodeConstraintElement : public virtual GeneralisedElement
  {
  public:
    /// Construcor. Needs the two node pointers so that we can retrieve the
    /// boundary data at solve time
    DuplicateNodeConstraintElement(
      Node* const& left_node_pt,
      Node* const& right_node_pt,
      CurvilineGeomObject* const& left_boundary_pt,
      CurvilineGeomObject* const& right_boundary_pt,
      Vector<double> const& left_coord,
      Vector<double> const& right_coord)
      : Left_node_pt(left_node_pt),
        Right_node_pt(right_node_pt),
        Left_boundary_pt(left_boundary_pt),
        Right_boundary_pt(right_boundary_pt),
        Left_node_coord(left_coord),
        Right_node_coord(right_coord)
    {
      // Add internal data which stores the eight Lagrange multipliers
      Index_of_lagrange_data = add_internal_data(new Data(18));

      // Add each node as external data
      Index_of_left_data = add_external_data(Left_node_pt);
      Index_of_right_data = add_external_data(Right_node_pt);
    }

    /// Destructor
    ~DuplicateNodeConstraintElement()
    {
      // Must remove Lagrange multiplier data?
    }

    /// Add the contribution to the residuals from the Lagrange multiplier
    /// constraining equations
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      fill_in_generic_residual_contribution_constraint(
        residuals, GeneralisedElement::Dummy_matrix, 0);
    }

    /// Add the contribution to the Jacobian from the Lagrange multiplier
    /// constraining equations
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      fill_in_generic_residual_contribution_constraint(
        residuals, jacobian, 1);
    }

    /// Validate constraints which contain no unpinned dofs and pin their
    /// corrosponding lagrange multiplier as it is used in no equations and
    /// it's own equation is trivially satisfied (Jacobian has a zero column
    /// and row if unpinned => singular)
    // [zdec] Do we want a bool in the element to determine whether we enforce
    // constraints that are already fully pinned (tears may be desired in some
    // dofs?)
    void validate_and_pin_redundant_constraints()
    {
      // [zdec] debug
      std::cout << "Constraint Element" << std::endl;

      // Start by unpinning all lagrange multipliers in case the boundary
      // conditions are less restrictive than previously
      internal_data_pt(Index_of_lagrange_data)->unpin_all();


      // [zdec] This full description might be overkill for the code but it will
      // go in my thesis.

      // We need to keep trak of which fvk dofs are already 'used' by Lagrange
      // constraints.  If dofs 3 and 4 in the right node (dw/dl_1, dw/dl_2) are
      // the only unpinned dofs between three lagrange constraints (e.g. 3,4,5),
      // then including all three constraints will result in a three
      // (consistent) linearly dependent equations and hence a singular matrix.
      // Therefore, each time we apply a constraint we must 'use' a dof by
      // marking it as effectively pinned by the lagrange constraint. Generally,
      // checking we are maximally constraining our duplicated nodes without
      // introducing linear dependent equations can be a tedious problem (we may
      // mark dof A as 'used' when choosing between A and B only for the next
      // constraint to contain only dof A) but we can safely mark the first free
      // dof provided we choose a constraint and dof order that prioritise
      // marking dofs which aren't used again (i.e. right dofs).

      // The number of fileds we are coupling
      unsigned n_field = 3;
      // The number of nodal types per field
      unsigned n_type = 6;
      // Total number of values per node
      unsigned n_nodal_val = n_field * n_type;

      // We use a vector of booleans to keep track of dofs that might be reused
      // (no need to track right dofs which are used once)
      std::vector<bool> right_data_used(n_nodal_val, false);
      std::vector<bool> left_data_used(n_nodal_val, false);

      // Store each data
      Data* left_data_pt = external_data_pt(Index_of_left_data);
      Data* right_data_pt = external_data_pt(Index_of_right_data);

      // We also want to store the jacobian and the hessian of the mapping
      DenseMatrix<double> jac_of_transform(2, 2, 0.0);
      Vector<DenseMatrix<double>> hess_of_transform(
        2, DenseMatrix<double>(2, 2, 0.0));
      get_jac_and_hess_of_coordinate_transform(jac_of_transform,
                                               hess_of_transform);


      // Each displacement uses the same six* constraints to constrain its C1
      // continuity between the left and right node - here we loop over the
      // displacement fields and apply these constraints to each.
      //
      //           * up to six, some may be unnecessary due to pinned values.
      for (unsigned i_field = 0; i_field < 3; i_field++)
      {
        // Constraint 0 uses dof 0 in each node (u_i_left = u_i_right)
        {
          // Index of the condition applied to this field, 0 corresponds to the
          // pure displacement dof type
          unsigned k_type = 0;
          // Index of the condition on the element (amongst all fields)
          unsigned condition_index = i_field * n_type + k_type;
          // Index of the val associated with displacement in the right node
          unsigned right_ui_index = i_field * n_type + k_type;
          // Index of the val associated with displacement in the left node
          unsigned left_ui_index = i_field * n_type + k_type;


          // Get whether each value is pinned
          bool right_ui_pinned = right_data_pt->is_pinned(right_ui_index);
          bool left_ui_pinned = left_data_pt->is_pinned(left_ui_index);

          // If anything is free, mark it as used and continue without doing
          // anything else
          if (!right_ui_pinned && !right_data_used[right_ui_index])
          {
            // [zdec] debug
            std::cout << "eqn " << condition_index << " depends on R dof"
                      << k_type << std::endl;
            right_data_used[right_ui_index] = true;
          }
          else if (!left_ui_pinned && !left_data_used[left_ui_index])
          {
            // [zdec] debug
            std::cout << "eqn " << condition_index << " depends on L dof"
                      << k_type << std::endl;
            left_data_used[left_ui_index] = true;
          }
          else
          {
            // --------------------------------------------------------------

            // [zdec] Below has to be true before solving here as there is no
            // possibility that free dofs are being used elsewhere

            // If we made it here, it is because all dofs in the constraint are
            // pinned so we need to check the constraint is satisfied manually
            // and then remove it by pinning the corresponding lagrange
            // multiplier

            // // Calculate the residual of the constraint
            // double constraint_residual =
            //   right_data_pt->value(condition_index) -
            //   left_data_pt->value(condition_index);
            // // Check that the constraint is met and we don't have a tear
            // if(constraint_residual > Constraint_tolerance)
            // {
            //   throw_unsatisfiable_constraint_error(condition_index,
            //   constraint_residual);
            // }

            // If it is met, we pin the lagrange multiplier that corresponds to
            // this constraint as it is redundant and results in a zero
            // row/column
            internal_data_pt(Index_of_lagrange_data)->pin(condition_index);
          }
        }

        // Constraints 1,2 use dofs 1,2 respectively from the right node and
        // both constraints use both dofs 1 and 2 in the left
        for (unsigned k_type = 1; k_type < 3; k_type++)
        {
          // The index of the derivative
          unsigned alpha = k_type - 1;
          // Index of the condition on the element
          unsigned condition_index = i_field * n_type + k_type;
          // Index of the right nodes alpha-th derivative value
          unsigned right_uida_index = i_field * n_type + k_type;
          // Index of the left nodes first derivative value
          unsigned left_uid1_index = i_field * n_type + 1;
          // Index of the left nodes second derivative value
          unsigned left_uid2_index = i_field * n_type + 2;

          // Get whether each value is pinned
          bool right_duida_pinned = right_data_pt->is_pinned(right_uida_index);
          bool left_duid1_pinned = left_data_pt->is_pinned(left_uid1_index);
          bool left_duid2_pinned = left_data_pt->is_pinned(left_uid2_index);

          // If anything is free, mark it as used and continue without doing
          // anything else. We also need to check that each dof hasn't become
          // decoupled from this constraint by ensuring that its coefficient (if
          // it has one) is sufficiently large (> Orthogonality_tolerance)
          if (!right_duida_pinned && !right_data_used[right_uida_index])
          {
            // // [zdec] debug
            // std::cout << "eqn " << condition_index << " depends on dof R"
            //        << right_duida_index << std::endl;
            right_data_used[right_uida_index] = true;
            continue;
          }
          if (!left_duid1_pinned && !left_data_used[left_uid1_index])
          {
            // // [zdec] debug
            // std::cout << "eqn " << condition_index << " depends on dof L"
            //        << left_duid1_index << std::endl;
            double coeff = jac_of_transform(0, alpha);
            if (fabs(coeff) > Orthogonality_tolerance)
            {
              left_data_used[left_uid1_index] = true;
              continue;
            }
          }
          if (!left_duid2_pinned && !left_data_used[left_uid2_index])
          {
            // [zdec] debug
            // std::cout << "eqn " << condition_index << " depends on dof L"
            //        << left_duid2_index << std::endl;
            double coeff = jac_of_transform(1, alpha);
            if (fabs(coeff) > Orthogonality_tolerance)
            {
              left_data_used[left_uid2_index] = true;
              continue;
            }
          }
          // ---------------------------------------------------------------------
          // If we made it here, it is because all dofs in the constraint are
          // pinned so we need to check the constraint is satisfied manually and
          // then remove it by pinning the corresponding lagrange multiplier

          // // Calculate the residual of the constraint
          // double constraint_residual = right_data_pt->value(condition_index);
          // for(unsigned beta = 0; beta < 2; beta++)
          // {
          //   constraint_residual +=
          //     - left_data_pt->value(3+beta) * jac_of_transform(beta,alpha);
          // }
          // // Check that the constraint is met and we don't have a tear
          // if(constraint_residual > Constraint_tolerance)
          // {
          //   throw_unsatisfiable_constraint_error(condition_index,
          //   constraint_residual);
          // }

          // If it is met, we pin the lagrange multiplier that corresponds to
          // this constraint as it is redundant and results in a zero row/column
          internal_data_pt(Index_of_lagrange_data)->pin(condition_index);
        }

        // Constraints 3-5 use dofs 3-5 respectively from the right node and
        // all use dofs 1-5 from the left node
        for (unsigned k_type = 3; k_type < 6; k_type++)
        {
          // Index of the condition on the element
          unsigned condition_index = i_field * n_type + k_type;
          // Index of the right nodes alpha-th derivative value
          unsigned right_uidadb_index = i_field * n_type + k_type;
          // Index of the left nodes d1 derivative value
          unsigned left_uid1_index = i_field * n_type + 1;
          // Index of the left nodes d2 derivative value
          unsigned left_uid2_index = i_field * n_type + 2;
          // Index of the left nodes d1d1 derivative value
          unsigned left_uid1d1_index = i_field * n_type + 3;
          // Index of the left nodes d1d2 derivative value
          unsigned left_uid1d2_index = i_field * n_type + 4;
          // Index of the left nodes d2d2 derivative value
          unsigned left_uid2d2_index = i_field * n_type + 5;

          // The index of the derivatives
          unsigned alpha, beta;
          switch (k_type)
          {
              // The third index is when both derivatives are first
            case 3:
              alpha = 0;
              beta = 0;
              break;
              // The fourth index is when the derivatives are mixed
            case 4:
              alpha = 0;
              beta = 1;
              break;
              // The fifth index is when both derivatives are second
            case 5:
              alpha = 1;
              beta = 1;
              break;
          }

          // Get whether each value is pinned
          bool right_uidadb_pinned =
            right_data_pt->is_pinned(right_uidadb_index);
          bool left_uid1_pinned = left_data_pt->is_pinned(left_uid1_index);
          bool left_uid2_pinned = left_data_pt->is_pinned(left_uid2_index);
          bool left_uid1d1_pinned = left_data_pt->is_pinned(left_uid1d1_index);
          bool left_uid1d2_pinned = left_data_pt->is_pinned(left_uid1d2_index);
          bool left_uid2d2_pinned = left_data_pt->is_pinned(left_uid2d2_index);

          // If anything is free, mark it as used and continue without doing
          // anything else. We also need to check that each dof hasn't become
          // decoupled from this constraint by ensuring that its coefficient
          // (if it has one) is sufficiently large (> Orthogonality_tolerance)
          if (!right_uidadb_pinned && !right_data_used[right_uidadb_index])
          {
            // [zdec] debug
            std::cout << "eqn " << condition_index << " depends on dof R"
                      << condition_index << std::endl;
            right_data_used[right_uidadb_index] = true;
            continue;
          }
          if (!left_uid1_pinned && !left_data_used[left_uid1_index])
          {
            double coeff = hess_of_transform[0](alpha, beta);
            if (fabs(coeff) > Orthogonality_tolerance)
            {
              // [zdec] debug
              std::cout << "eqn " << condition_index << " depends on dof L3"
                        << std::endl;
              left_data_used[left_uid1_index] = true;
              continue;
            }
          }
          if (!left_uid2_pinned && !left_data_used[left_uid2_index])
          {
            double coeff = hess_of_transform[1](alpha, beta);
            if (fabs(coeff) > Orthogonality_tolerance)
            {
              // [zdec] debug
              std::cout << "eqn " << condition_index << " depends on dof L4"
                        << std::endl;
              left_data_used[left_uid2_index] = true;
              continue;
            }
          }
          if (!left_uid1d1_pinned && !left_data_used[left_uid1d1_pinned])
          {
            double coef =
              jac_of_transform(0, alpha) * jac_of_transform(0, beta);
            if (fabs(coef) > Orthogonality_tolerance)
            {
              // [zdec] debug
              std::cout << "eqn " << condition_index << " depends on dof L5"
                        << std::endl;
              left_data_used[left_uid1d1_index] = true;
              continue;
            }
          }
          if (!left_uid1d2_pinned && !left_data_used[left_uid1d2_index])
          {
            double coef =
              jac_of_transform(0, alpha) * jac_of_transform(1, beta) +
              jac_of_transform(1, alpha) * jac_of_transform(0, beta);
            if (fabs(coef) > Orthogonality_tolerance)
            {
              // [zdec] debug
              std::cout << "eqn " << condition_index << " depends on dof L6"
                        << std::endl;
              left_data_used[left_uid1d2_index] = true;
              continue;
            }
          }
          if (!left_uid2d2_pinned && !left_data_used[left_uid2d2_index])
          {
            double coef =
              jac_of_transform(1, alpha) * jac_of_transform(1, beta);
            if (fabs(coef) > Orthogonality_tolerance)
            {
              // [zdec] debug
              std::cout << "eqn " << condition_index << " depends on dof L7"
                        << std::endl;
              left_data_used[left_uid2d2_index] = true;
              continue;
            }
          }
          // -------------------------------------------------------------------
          // if we made it here, it is because all dofs in the constraint are
          // pinned so we need to check the constraint is satisfied manually
          // and then remove it by pinning the corresponding lagrange
          // multiplier

          // // calculate the residual of the constraint
          // double constraint_residual = right_data_pt->value(condition_index);
          // for(unsigned gamma = 0; gamma < 2; gamma++)
          // {
          //   constraint_residual +=
          //     - left_data_pt->value(3+gamma) *
          //     hess_of_transform[gamma](alpha,beta);
          //   for(unsigned delta = 0; delta < 2; delta++)
          //   {
          //     constraint_residual +=
          //       - left_data_pt->value(5+gamma+delta)
          //       * jac_of_transform(gamma,alpha)
          //       * jac_of_transform(delta,beta);
          //   }
          // }
          // // check that the constraint is met and we don't have a tear
          // if(constraint_residual > constraint_tolerance)
          // {
          //   throw_unsatisfiable_constraint_error(condition_index,
          //   constraint_residual);
          // }

          // if it is met, we pin the lagrange multiplier that corresponds to
          // this constraint as it is redundant and results in a zero
          // row/column
          internal_data_pt(Index_of_lagrange_data)->pin(condition_index);
        }
      }
    } // End validate_and_pin_redundant_constraints()


  private:
    /// Throw an error about a constraint that cannot be satisfied as it has no
    /// free variables but still has a residual greater than a requested error
    /// tokerabce. Takes the index and the residual of the offending constraint
    void throw_unsatisfiable_constraint_error(const unsigned& i,
                                              const double& res)
    {
      // Get the position of the nodes so we can be a little helpful about
      // where the boundary conditions are contradictory.
      Vector<double> x(2, 0.0);
      Left_boundary_pt->position(Left_node_coord, x);
      std::string error_string =
        "Constraint " + std::to_string(i) + " on the nodes at x = (" +
        std::to_string(x[0]) + ", " + std::to_string(x[1]) +
        ") has no free variables but is not satisfied to within the " +
        "tolerance (" + std::to_string(Constraint_tolerance) + ")." +
        "The residual of the constraint is: C_" +
        std::to_string(Constraint_tolerance) + " = " + std::to_string(res) +
        "\n";
      throw OomphLibError(
        error_string, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    } // End of throw_unsatisfiable_constraint_error


    /// Function to calculate Jacobian and Hessian of the coordinate mapping
    void get_jac_and_hess_of_coordinate_transform(
      DenseMatrix<double>& jac_of_transform,
      Vector<DenseMatrix<double>>& hess_of_transform)
    {
      //----------------------------------------------------------------------
      // We need the parametrisations either side of the vertex which define
      // the coordinates each node uses for its Hermite dofs.
      Vector<double> left_x(2, 0.0); // [zdec] debug
      Vector<double> right_x(2, 0.0); // [zdec] debug
      Vector<double> left_dxids(2, 0.0);
      Vector<double> left_d2xids2(2, 0.0);
      Vector<double> right_dxids(2, 0.0);
      Vector<double> right_d2xids2(2, 0.0);
      Left_boundary_pt->position(Left_node_coord, left_x); // [zdec] debug
      Right_boundary_pt->position(Right_node_coord, right_x); // [zdec] debug
      Left_boundary_pt->dposition(Left_node_coord, left_dxids);
      Left_boundary_pt->d2position(Left_node_coord, left_d2xids2);
      Right_boundary_pt->dposition(Right_node_coord, right_dxids);
      Right_boundary_pt->d2position(Right_node_coord, right_d2xids2);

      // Get the speed of each parametrisation
      double left_mag =
        sqrt(left_dxids[0] * left_dxids[0] + left_dxids[1] * left_dxids[1]);
      double right_mag =
        sqrt(right_dxids[0] * right_dxids[0] + right_dxids[1] * right_dxids[1]);

      //----------------------------------------------------------------------
      // Normalise dxids to find the tangent vectors and their
      // derivatives either side of the vertex
      Vector<double> left_ti(2, 0.0);
      Vector<double> left_ni(2, 0.0);
      Vector<double> left_dtids(2, 0.0);
      Vector<double> left_dnids(2, 0.0);
      Vector<double> right_ti(2, 0.0);
      Vector<double> right_ni(2, 0.0);
      Vector<double> right_dtids(2, 0.0);
      Vector<double> right_dnids(2, 0.0);
      for (unsigned alpha = 0; alpha < 2; alpha++)
      {
        // Fill in the tangents either side of the vertex
        left_ti[alpha] = left_dxids[alpha] / left_mag;
        right_ti[alpha] = right_dxids[alpha] / right_mag;
        // Fill in the derivatives of the (normalised) tangents either side of
        // the vertex
        left_dtids[alpha] =
          left_d2xids2[alpha] / std::pow(left_mag, 2) -
          (left_dxids[0] * left_d2xids2[0] + left_dxids[1] * left_d2xids2[1]) *
            left_dxids[alpha] / std::pow(left_mag, 4);
        right_dtids[alpha] = right_d2xids2[alpha] / std::pow(right_mag, 2) -
                             (right_dxids[0] * right_d2xids2[0] +
                              right_dxids[1] * right_d2xids2[1]) *
                               right_dxids[alpha] / std::pow(right_mag, 4);
        // Use these to fill out the corresponding vectors for the normal
        // direction (nx,ny) = (ty,-tx)
      }
      // Use orthogonality to fill in normals and their derivatives
      for (unsigned alpha = 0; alpha < 2; alpha++)
      {
        left_ni[alpha] = pow(-1, alpha) * left_ti[(alpha + 1) % 2];
        right_ni[alpha] = pow(-1, alpha) * right_ti[(alpha + 1) % 2];
        left_dnids[alpha] = pow(-1, alpha) * left_dtids[(alpha + 1) % 2];
        right_dnids[alpha] = pow(-1, alpha) * right_dtids[(alpha + 1) % 2];
      }

      //----------------------------------------------------------------------
      // We need to fill out the Jacobians and Hessians of the boundary
      // coordinates either side of the vertex
      DenseMatrix<double> left_jac(2, 2, 0.0);
      DenseMatrix<double> right_jac(2, 2, 0.0);
      Vector<DenseMatrix<double>> left_hess(2, DenseMatrix<double>(2, 2, 0.0));
      Vector<DenseMatrix<double>> right_hess(2, DenseMatrix<double>(2, 2, 0.0));
      for (unsigned alpha = 0; alpha < 2; alpha++)
      {
        // Fill in Jacobians {{nx,tx},{ny,ty}}
        left_jac(alpha, 0) = left_ni[alpha];
        left_jac(alpha, 1) = left_ti[alpha];
        right_jac(alpha, 0) = right_ni[alpha];
        right_jac(alpha, 1) = right_ti[alpha];
        // Fill in Hessians
        // left_hess[alpha](0,0) = 0.0;
        left_hess[alpha](0, 1) = left_dnids[alpha];
        left_hess[alpha](1, 0) = left_dnids[alpha];
        left_hess[alpha](1, 1) = left_dtids[alpha];
        // right_hess[alpha](0,0) = 0.0;
        right_hess[alpha](0, 1) = right_dnids[alpha];
        right_hess[alpha](1, 0) = right_dnids[alpha];
        right_hess[alpha](1, 1) = right_dtids[alpha];
      }

      //----------------------------------------------------------------------
      // We need the inverse Jacobian and Hessian for the left parametrisation
      DenseMatrix<double> left_jac_inv(2, 2, 0.0);
      Vector<DenseMatrix<double>> left_hess_inv(2,
                                                DenseMatrix<double>(2, 2, 0.0));
      left_jac_inv(0, 0) = left_jac(1, 1);
      left_jac_inv(0, 1) = -left_jac(0, 1);
      left_jac_inv(1, 0) = -left_jac(1, 0);
      left_jac_inv(1, 1) = left_jac(0, 0);
      // Fill out inverse of Hessian
      // H^{-1}abg = J^{-1}ad Hdez J^{-1}eb J^{-1}zg
      for (unsigned alpha = 0; alpha < 2; alpha++)
      {
        for (unsigned beta = 0; beta < 2; beta++)
        {
          for (unsigned gamma = 0; gamma < 2; gamma++)
          {
            for (unsigned alpha2 = 0; alpha2 < 2; alpha2++)
            {
              for (unsigned beta2 = 0; beta2 < 2; beta2++)
              {
                for (unsigned gamma2 = 0; gamma2 < 2; gamma2++)
                {
                  left_hess_inv[alpha](beta, gamma) -=
                    left_jac_inv(alpha, alpha2) *
                    left_hess[alpha2](beta2, gamma2) *
                    left_jac_inv(beta2, beta) * left_jac_inv(gamma2, gamma);
                }
              }
            }
          }
        }
      }

      //----------------------------------------------------------------------
      //----------------------------------------------------------------------
      // Use these to calculate the Jacobian of the left->right transform
      //     J = J_{left}^{-1}J_{right}
      // and the Hessian of the left->right transform
      //     H = H_{left}^{-1}J_{right}J_{right} + J_{left}^{-1}H_{right}
      for (unsigned alpha = 0; alpha < 2; alpha++)
      {
        for (unsigned beta = 0; beta < 2; beta++)
        {
          for (unsigned gamma = 0; gamma < 2; gamma++)
          {
            // Add contribution to J
            jac_of_transform(alpha, beta) +=
              left_jac_inv(alpha, gamma) * right_jac(gamma, beta);
            for (unsigned mu = 0; mu < 2; mu++)
            {
              // Add second term contribution to H
              hess_of_transform[alpha](beta, gamma) +=
                left_jac_inv(alpha, mu) * right_hess[mu](beta, gamma);
              for (unsigned nu = 0; nu < 2; nu++)
              {
                // Add first term contribution to H
                hess_of_transform[alpha](beta, gamma) +=
                  left_hess_inv[alpha](mu, nu) * right_jac(mu, beta) *
                  right_jac(nu, gamma);
              }
            }
          }
        }
      }

      // // [zdec] debug
      // std::ofstream jac_and_hess;

      // jac_and_hess.open("corner_jac_and_hess_new.csv", std::ios_base::app);
      // jac_and_hess << "Jacobian :" << std::endl
      //              << jac_of_transform(0, 0) << " " << jac_of_transform(0, 1)
      //              << std::endl
      //              << jac_of_transform(1, 0) << " " << jac_of_transform(1, 1)
      //              << std::endl
      //              << "Hessian [x]:" << std::endl
      //              << hess_of_transform[0](0, 0) << " " <<
      //              hess_of_transform[0](0, 1)
      //              << std::endl
      //              << hess_of_transform[0](1, 0) << " " <<
      //              hess_of_transform[0](1, 1)
      //              << std::endl
      //              << "Hessian [y]:" << std::endl
      //              << hess_of_transform[1](0, 0) << " " <<
      //              hess_of_transform[1](0, 1)
      //              << std::endl
      //              << hess_of_transform[1](1, 0) << " " <<
      //              hess_of_transform[1](1, 1)
      //              << std::endl
      //              << std::endl;
      // jac_and_hess.close();


      // jac_and_hess.open("invleft_jac_and_hess_new.csv", std::ios_base::app);
      // jac_and_hess << "Jacobian :" << std::endl
      //              << left_jac_inv(0, 0) << " " << left_jac_inv(0, 1) <<
      //              std::endl
      //              << left_jac_inv(1, 0) << " " << left_jac_inv(1, 1) <<
      //              std::endl
      //              << "Hessian [x]:" << std::endl
      //              << left_hess_inv[0](0, 0) << " " << left_hess_inv[0](0, 1)
      //              << std::endl
      //              << left_hess_inv[0](1, 0) << " " << left_hess_inv[0](1, 1)
      //              << std::endl
      //              << "Hessian [y]:" << std::endl
      //              << left_hess_inv[1](0, 0) << " " << left_hess_inv[1](0, 1)
      //              << std::endl
      //              << left_hess_inv[1](1, 0) << " " << left_hess_inv[1](1, 1)
      //              << std::endl
      //              << std::endl;
      // jac_and_hess.close();

      // jac_and_hess.open("left_jac_and_hess_new.csv", std::ios_base::app);
      // jac_and_hess << "Jacobian :" << std::endl
      //              << left_jac(0, 0) << " " << left_jac(0, 1) << std::endl
      //              << left_jac(1, 0) << " " << left_jac(1, 1) << std::endl
      //              << "Hessian [x]:" << std::endl
      //              << left_hess[0](0, 0) << " " << left_hess[0](0, 1)
      //              << std::endl
      //              << left_hess[0](1, 0) << " " << left_hess[0](1, 1)
      //              << std::endl
      //              << "Hessian [y]:" << std::endl
      //              << left_hess[1](0, 0) << " " << left_hess[1](0, 1)
      //              << std::endl
      //              << left_hess[1](1, 0) << " " << left_hess[1](1, 1)
      //              << std::endl
      //              << std::endl;
      // jac_and_hess.close();

      // jac_and_hess.open("right_jac_and_hess_new.csv", std::ios_base::app);
      // jac_and_hess << "Jacobian :" << std::endl
      //              << right_jac(0, 0) << " " << right_jac(0, 1) << std::endl
      //              << right_jac(1, 0) << " " << right_jac(1, 1) << std::endl
      //              << "Hessian [x]:" << std::endl
      //              << right_hess[0](0, 0) << " " << right_hess[0](0, 1)
      //              << std::endl
      //              << right_hess[0](1, 0) << " " << right_hess[0](1, 1)
      //              << std::endl
      //              << "Hessian [y]:" << std::endl
      //              << right_hess[1](0, 0) << " " << right_hess[1](0, 1)
      //              << std::endl
      //              << right_hess[1](1, 0) << " " << right_hess[1](1, 1)
      //              << std::endl
      //              << std::endl;
      // jac_and_hess.close();


      // // [zdec] debug
      // std::ofstream debug_stream;
      // debug_stream.open("left_norm_and_tan.dat", std::ios_base::app);
      // debug_stream << left_x[0] << " " << left_x[1] << " " << left_ni[0] << "
      // "
      //              << left_ni[1] << " " << left_ti[0] << " " << left_ti[1] <<
      //              " "
      //              << left_dnids[0] << " " << left_dnids[1] << " " <<
      //              left_dtids[0]
      //              << " " << left_dtids[1] << " " << left_d2xids2[0] << " "
      //              << left_d2xids2[1] << std::endl;
      // debug_stream.close();
      // debug_stream.open("right_norm_and_tan.dat", std::ios_base::app);
      // debug_stream << right_x[0] << " " << right_x[1] << " " << right_ni[0]
      // << " "
      //              << right_ni[1] << " " << right_ti[0] << " " << right_ti[1]
      //              << " "
      //              << right_dnids[0] << " " << right_dnids[1] << " " <<
      //              right_dtids[0]
      //              << " " << right_dtids[1] << " " << right_d2xids2[0] << " "
      //              << right_d2xids2[1] << std::endl;
      // debug_stream.close();

    } // End get_jac_and_hess_of_coordinate_transform


    /// Add the contribution to the residuals (and jacobain if flag is 1) from
    /// the Lagrange multiplier constraining equations
    void fill_in_generic_residual_contribution_constraint(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      const unsigned& flag)
    {
      // Number of fields we are constraining
      unsigned n_field = 3;
      // Number of nodal types per field
      unsigned n_type = 6;
      // Total number of nodal values
      unsigned n_val = n_field * n_type;
      // Dimension of constrained elements
      unsigned dim = 2;

      // [zdec] debug
      std::cout << std::endl
                << std::endl
                << "ADD CONTRIBUTION FROM CONSTRAINTS" << std::endl
                << "=============================================" << std::endl;
      //----------------------------------------------------------------------
      //----------------------------------------------------------------------
      // Calculate Jacobian and Hessian of coordinate transform between
      // each boundary coordinate
      DenseMatrix<double> jac_of_transform(dim, dim, 0.0);
      Vector<DenseMatrix<double>> hess_of_transform(
        dim, DenseMatrix<double>(dim, dim, 0.0));
      get_jac_and_hess_of_coordinate_transform(jac_of_transform,
                                               hess_of_transform);

      //----------------------------------------------------------------------
      //----------------------------------------------------------------------
      // Use the jac and hess of transform to add the residual
      // contributions from the constraint
      // [zdec]::TODO make indexing (alpha,beta,gamma,...) consistent

      // Store the internal data pointer which stores the Lagrange multipliers
      Vector<double> lagrange_value(n_val, 0.0);
      internal_data_pt(Index_of_lagrange_data)->value(lagrange_value);

      // Store the left and right nodal dofs with indices
      // 0: u0                6: u1
      // 1: du0/dx0           7: du1/dx0
      // 2: du0/dx1               :
      // 3: d^2u0/dx0^2           :
      // 4: d^2u0/dx0dx1          :
      // 5: d^2u0/dx1^2      17: d^2u2/dx1^2
      Vector<double> left_value(n_val, 0.0);
      Vector<double> right_value(n_val, 0.0);
      Left_node_pt->value(left_value);
      Right_node_pt->value(right_value);


      //----------------------------------------------------------------------
      // First the contributions to the right node external equations as these
      // are independent of both field and value type
      for (unsigned i_val = 0; i_val < n_val; i_val++)
      {
        int right_eqn_number = external_local_eqn(Index_of_right_data, i_val);

        // If this dof isn't pinned we add to the residual
        if (right_eqn_number >= 0)
        {
          // Right dof term in the constraint always lambda_i*W_i
          residuals[right_eqn_number] += lagrange_value[i_val];

          // If flag, then add the jacobian contribution
          if (flag)
          {
            // The contributions to the right node's equations are just
            // r_i * L_i
            int lagrange_dof_number =
              internal_local_eqn(Index_of_lagrange_data, i_val);
            // If this dof isn't pinned then add the contributions
            // [zdec] should never be pinned if right value is not
            if (lagrange_dof_number >= 0)
            {
              // Add the contribution to the jacobian
              jacobian(right_eqn_number, lagrange_dof_number) += 1.0;
              // And by symmetry, we can add the transpose contribution to the
              // jacobian
              jacobian(lagrange_dof_number, right_eqn_number) += 1.0;
            } // End pinned check
          } // End Jacobian contribution [if (flag)]
        }
      } // End for loop adding contributions to right nodal equations


      // Now, loop over the three displacement fields as the left node residual
      // contributions as well as the lagrange residual contributions are the
      // same for each field
      for (unsigned i_field = 0; i_field < n_field; i_field++)
      {
        // ---------------------------------------------------------------------
        // The contributions to the left node external equations First
        // is the displacement: - lambda_i*(u_i)
        {
          // Index of the left nodes i-th displacement (0-th type) value
          unsigned left_ui_index = i_field * n_type + 0;
          // External equation index
          int left_eqn_number =
            external_local_eqn(Index_of_left_data, left_ui_index);
          // If this dof isn't pinned we add to the residual
          if (left_eqn_number >= 0)
          {
	    // Add residual contribution which comes from lagrange value
	    unsigned lagrange_index = i_field * n_type + 0;
            residuals[left_eqn_number] += - lagrange_value[lagrange_index];

	    // Get the local equation number for this dof and add the jacobian
	    // contribution if unpinned (and we are making the jacobian)
            // [zdec] should never be pinned if right value is not
            int lagrange_dof_number =
              internal_local_eqn(Index_of_lagrange_data, lagrange_index);
            if (flag && (lagrange_dof_number >= 0))
            {
	      // Add the contribution to the jacobian
	      jacobian(left_eqn_number, lagrange_dof_number) += -1.0;
              // And by symmetry, we can add the transpose contribution to the
              // jacobian
              jacobian(lagrange_dof_number, left_eqn_number) += -1.0;
            }
          }
        }

        // Next two are from left gradient of u_i:
        //     - lambda_{1+\beta} * ui_{1+\alpha} * J_{\alpha\beta}
        //     - lambda_{3+\beta+\gamma} * ui_{1+\alpha} * H_{\alpha\beta\gamma}
        // gamma>=beta so we don't double count lambda_6 condition
        for (unsigned alpha = 0; alpha < 2; alpha++)
        {
          // Index of the left nodes d1 derivative value
          unsigned left_uida_index = i_field * n_type + (1 + alpha);
          // Eqn number is the index of the alpha-th derivative of w which is
          // the 1+alpha-th dof (alpha=0,1)
          int left_eqn_number =
            external_local_eqn(Index_of_left_data, left_uida_index);
          // If this dof isn't pinned we add to the residual
          if (left_eqn_number >= 0)
          {
            // Loop over the lagrange multipliers associated with the right
            // first derivatives
            for (unsigned beta = 0; beta < 2; beta++)
            {
              // Add residual contribution from the lagrange value associated
              // with the right beta-th derivative
              unsigned lagrange_uidb_index = i_field * n_type + (1 + beta);
              residuals[left_eqn_number] +=
                - lagrange_value[lagrange_uidb_index]
		* jac_of_transform(alpha, beta);

              // Get the local equation number for this dof and add the jacobian
              // contribution if unpinned (and we are making the jacobian)
              int lagrange_dof_number =
                internal_local_eqn(Index_of_lagrange_data, lagrange_uidb_index);
              if (flag && (lagrange_dof_number >= 0))
              {
                double jac_term = - jac_of_transform(alpha, beta);
                // Orthogonality check (for jacobian cleanliness)
                if (fabs(jac_term) > Orthogonality_tolerance)
                {
                  // Add the contribution to the jacobian
                  jacobian(left_eqn_number, lagrange_dof_number) += jac_term;
                  // And by symmetry, we can add the transpose contribution to
                  // the jacobian
                  jacobian(lagrange_dof_number, left_eqn_number) += jac_term;
                } // End orthogonality check
              } // End jacobian contributions for first derivative terms

              // Loop over the lagrange multipliers associated with the right
              // second derivatives. gamma>=beta so we don't double count the
              // mixed derivative constraint
              for (unsigned gamma = beta; gamma < 2; gamma++)
              {
		// Add residual contribution from the lagrange value associated
		// with the right beta+gamma-th second derivative
                unsigned lagrange_uidbdg_index =
                  i_field * n_type + (3 + beta + gamma);
                residuals[left_eqn_number] +=
                  - lagrange_value[lagrange_uidbdg_index]
		  * hess_of_transform[alpha](beta, gamma);

		// Get the local equation number for this dof and add the
		// jacobian contribution if unpinned (and we are making the
		// jacobian)
                int lagrange_dof_number = internal_local_eqn(
                  Index_of_lagrange_data, lagrange_uidbdg_index);
                if (flag && (lagrange_dof_number >= 0))
                {
                  double jac_term = - hess_of_transform[alpha](beta, gamma);
                  // Orthogonality check
                  if (fabs(jac_term) > Orthogonality_tolerance)
                  {
                    // Add the contribution to the jacobian
                    jacobian(left_eqn_number, lagrange_dof_number) += jac_term;
                    // And by symmetry, we can add the transpose contribution
                    // to the jacobian
                    jacobian(lagrange_dof_number, left_eqn_number) += jac_term;
                  } // End of orthogonality check
                } // End of jacobian contribution for second derivative terms
              } // End loop over second derivative Lagrange multipliers [gamma]
            } // End loop over first derivative Lagrange multipliers [beta]
          } // End of if unpinned
        } // End loop adding contributions to the left nodal gradient equations

        // Last three are the left second derivatives of u_i (delta>gamma):
        //     - lambda_{3+\gamma+\delta} * w_{3+\alpha+\beta}
        //       * J_{\alpha\gamma} * J_{\beta\delta}
        // Index second derivative (equation) using alpha & beta
        for (unsigned alpha = 0; alpha < 2; alpha++)
        {
          // Note that d^2w/ds_1ds_2 is counted twice in the summation so we
          // allow alpha and beta to loop over both indices (unlike gamma+delta)
          for (unsigned beta = 0; beta < 2; beta++)
          {
            // Index of lagrange value associated with the right alpha+beta-th
            // second derivatives
            unsigned left_uidadb_index = i_field * n_type + (3 + alpha + beta);
            // Eqn number is the index of the second derivative of w which is
            // the 3+alpha+beta-th dof
            int left_eqn_number =
              external_local_eqn(Index_of_left_data, left_uidadb_index);
            // If this dof isn't pinned we add to the residual
            if (left_eqn_number >= 0)
            {
              // Index lagrange multipliers using gamma and delta
              for (unsigned gamma = 0; gamma < 2; gamma++)
              {
                // delta>=gamma so we don't double count the lagrange_value
                // associated with the mixed derivative
                for (unsigned delta = gamma; delta < 2; delta++)
                {
                  // Add residual contribution
                  unsigned lagrange_uidgdd_index =
                    i_field * n_type + (3 + gamma + delta);
                  residuals[left_eqn_number] +=
                    -lagrange_value[lagrange_uidgdd_index] *
                    jac_of_transform(alpha, gamma) *
                    jac_of_transform(beta, delta);

		  // Get the local equation number for this dof and add the
		  // jacobian contribution if unpinned (and we are making the
		  // jacobian)
                  int lagrange_dof_number = internal_local_eqn(
                    Index_of_lagrange_data, lagrange_uidgdd_index);
                  if (flag && (lagrange_dof_number >= 0))
                  {
                    // Find the jacobian matrix contribution
                    double jac_term = -jac_of_transform(alpha, gamma) *
                                      jac_of_transform(beta, delta);
                    // Orthogonality check
                    if (fabs(jac_term) > Orthogonality_tolerance)
                    {
                      // Add the contribution to the jacobian
                      jacobian(left_eqn_number, lagrange_dof_number) +=
                        jac_term;
                      // And by symmetry, we can add the transpose contribution
                      // to the jacobian
                      jacobian(lagrange_dof_number, left_eqn_number) +=
                        jac_term;
                    }
                  } // End jacobian

                }
              } // End loops over the conditions (gamma,delta)
            } // End if dof isn't pinned
          }
        } // End loops adding contributions to the left nodal curvature
          // equations (alpha,beta)


        //----------------------------------------------------------------------
        // Now add contributions to the internal (lagrange multiplier) equations
        // (note jacobian contributions will have already been thanks to use of
        // symmetry)

        // First pair of dofs are equal
        {
          // Index of the condition on the nodes (displacement constraint
          // corresponds to zeroth type for each field)
          unsigned lagrange_index = i_field * n_type + 0;
          // Get the internal data eqn number for this constraint
          int lagrange_eqn_number =
            internal_local_eqn(Index_of_lagrange_data, lagrange_index);
          // If this dof isn't pinned we add to the residual
          if (lagrange_eqn_number >= 0)
          {
	    // Add contributions from right and left nodes displacement
	    // (zero-th) type
	    unsigned right_ui_index = i_field * n_type + 0;
	    unsigned left_ui_index = i_field * n_type + 0;
            residuals[lagrange_eqn_number] +=
              (right_value[right_ui_index] - left_value[left_ui_index]);
          }
	  // Jacobian is added during right and left nodal equations using
	  // symmetry
        }

        // Next two (first derivatives of w) are related by
        //     grad_r(w) = grad_l(w)*J
        // where  J is the Jacobian grad_r(left coords)
        for (unsigned alpha = 0; alpha < 2; alpha++)
        {
	  // Index of the condition on the nodes
	  unsigned lagrange_index = i_field * n_type + (1 + alpha);
          // Get the internal data eqn number for this constraint
          int lagrange_eqn_number =
            internal_local_eqn(Index_of_lagrange_data, lagrange_index);
          // If this dof isn't pinned we add to the residual
          if (lagrange_eqn_number >= 0)
          {
	    // Add contribution from right node
	    unsigned right_uida_index = i_field * n_type + (1 + alpha);
            residuals[lagrange_eqn_number] +=
              (right_value[right_uida_index]);
	    // Add contribuions from left node
            for (unsigned beta = 0; beta < 2; beta++)
            {
	      unsigned left_uidb_index = i_field * n_type + (1 + beta);
              residuals[lagrange_eqn_number] +=
                - left_value[left_uidb_index] * jac_of_transform(beta, alpha);
            }
	    // Jacobian is added during right and left nodal equations using
	    // symmetry
          }
        }

        // Final three (second derivatives of w) are related by:
        //     grad_r(grad_r(w)) = grad_l(grad_l(w))*J*J + grad_l(w)*H
        // where H is the Hessian: grad_r(grad_r(left coords))
        // Loop over index of first derivative (0 or 1)
        for (unsigned alpha = 0; alpha < 2; alpha++)
        {
          // Loop over index of second derivative
          // (>=alpha to prevent double counting mixed deriv)
          for (unsigned beta = alpha; beta < 2; beta++)
          {
	    // Index of the condition on the nodes
	    unsigned lagrange_index = i_field * n_type + (3 + alpha + beta);
	    // Get the internal data eqn number for this constraint
	    int lagrange_eqn_number =
	      internal_local_eqn(Index_of_lagrange_data, lagrange_index);
	    // If this dof isn't pinned we add to the residual
	    if (lagrange_eqn_number >= 0)
	    {
	      // Add contributions from right node
	      unsigned right_uidadb_index = i_field * n_type + (3 + alpha + beta);
	      residuals[lagrange_eqn_number] +=
		right_value[right_uidadb_index];
	      // Loop over the left node derivatives
	      for (unsigned gamma = 0; gamma < 2; gamma++)
	      {
		// Add contributions from left node first derivatives
		unsigned left_uidg_index = i_field * n_type + (1 + gamma);
		residuals[lagrange_eqn_number] +=
		  - left_value[left_uidg_index]
		  * hess_of_transform[gamma](alpha, beta);
		// Loop over the left derivatives again to get second
		// derivatives
		for (unsigned delta = 0; delta < 2; delta++)
		{
		  // Add contributions from left node second derivatives
		  unsigned left_uidgdd_index =
		    i_field * n_type + (3 + gamma + delta);
		  residuals[lagrange_eqn_number] +=
		    - left_value[left_uidgdd_index]
		    * jac_of_transform(gamma, alpha)
		    * jac_of_transform(delta, beta);
		}
	      }
	      // Jacobian is added during right and left nodal equations using
	      // symmetry
	    } // End if eqn not pinned
          }
        }
      } // End loop over displacements [i_field]
    } // End fill_in_generic_residual_contribution_constraint

    /// Store the index of the internal data keeping the Lagrange multipliers
    unsigned Index_of_lagrange_data;

    /// Store the index of the external data for the left node
    unsigned Index_of_left_data;

    /// Store the index of the external data for the right node
    unsigned Index_of_right_data;

    /// Pointer to the left node (before the vertex when anticlockwise)
    Node* Left_node_pt;

    /// Pointer to the right node (after the vertex when anticlockwise)
    Node* Right_node_pt;

    /// Pointer to the left node's boundary parametrisation
    CurvilineGeomObject* Left_boundary_pt;

    /// Pointer to the right node's boundary parametrisation
    CurvilineGeomObject* Right_boundary_pt;

    /// Coordinate of the left node on the left boundary
    Vector<double> Left_node_coord;

    /// Coordinate of the left node on the left boundary
    Vector<double> Right_node_coord;

    /// Tolerance for validating fully pinned constraints
    // [zdec] does this wnat to be the problem residual tolerance?
    double Constraint_tolerance = 1.0e-10;

    /// Tolerance for checking whether a dof has become decoupled from an
    /// equation.
    /// i.e. in the equation y=Ax, how small does A have to be before y no
    /// longer /numerically/ depends on x? This becomes relevant when derivative
    /// directions become orthogonal, we need to ensure they aren't considered
    /// linearly dependent. (We choose this to be slightly larger than machine
    /// precision and it shouldn't generally need to be touched)
    double Orthogonality_tolerance = 1.0e-15;
  }; // End of DuplicateNodeConstraintElement class definition


} // namespace oomph


#endif

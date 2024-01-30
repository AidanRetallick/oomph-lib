// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC//    Version 1.0; svn revision $LastChangedRevision: 1097 $
// LIC//
// LIC// $LastChangedDate: 2015-12-17 11:53:17 +0000 (Thu, 17 Dec 2015) $
// LIC//
// LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
// LIC//
// LIC// This library is free software; you can redistribute it and/or
// LIC// modify it under the terms of the GNU Lesser General Public
// LIC// License as published by the Free Software Foundation; either
// LIC// version 2.1 of the License, or (at your option) any later version.
// LIC//
// LIC// This library is distributed in the hope that it will be useful,
// LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
// LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// LIC// Lesser General Public License for more details.
// LIC//
// LIC// You should have received a copy of the GNU Lesser General Public
// LIC// License along with this library; if not, write to the Free Software
// LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// LIC// 02110-1301  USA.
// LIC//
// LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
// LIC//
// LIC//====================================================================

#ifndef OOMPH_FVK_CURVABLE_BELL_ELEMENTS_HEADER
#define OOMPH_FVK_CURVABLE_BELL_ELEMENTS_HEADER

#include "foeppl_von_karman_equations.h"
#include "src/generic/bell_element_basis.h"
#include "src/generic/c1_curved_elements.h"
#include "src/generic/my_geom_object.h"
#include "src/generic/subparametric_Telement.h"
#include "src/generic/oomph_definitions.h"

// [zdec] TODO:
// -- Move function definitions to .cc
// -- Optimise basis and test functions by not retrieving d and d2 basis every
// time
// --

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
        // Hessian of inverse mapping [zdec] (...this can be found by
        // hand)
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


        // // [zdec] debug
        // std::ofstream jac_and_hess;
        // jac_and_hess.open("jac_and_hess_new.csv", std::ios_base::app);
        // jac_and_hess << "Jacobian inverse:" << std::endl
        // 		   << bi[0][0] << " " << bi[0][1] << std::endl
        // 		   << bi[1][0] << " " << bi[1][1] << std::endl
        // 		   << "Hessian inverse [x]:" << std::endl
        // 		   << Dbi[0](0,0) << " " << Dbi[0](0,1) << std::endl
        // 		   << Dbi[0](1,0) << " " << Dbi[0](1,1) << std::endl
        // 		   << "Hessian inverse [y]:" << std::endl
        // 		   << Dbi[1](0,0) << " " << Dbi[1](0,1) << std::endl
        // 		   << Dbi[1](1,0) << " " << Dbi[1](1,1) << std::endl <<
        // std::endl;


        // // [zdec] debug
        // std::ofstream jac_and_hess;
        // jac_and_hess.open("jac_and_hess_new.csv", std::ios_base::app);
        // jac_and_hess << "Jacobian inverse:" << std::endl
        //              << jac_inv(0, 0) << " " << jac_inv(0, 1) << std::endl
        //              << jac_inv(1, 0) << " " << jac_inv(1, 1) << std::endl
        //              << "Hessian inverse [x]:" << std::endl
        //              << hess_inv[0](0, 0) << " " << hess_inv[0](0, 1)
        //              << std::endl
        //              << hess_inv[0](1, 0) << " " << hess_inv[0](1, 1)
        //              << std::endl
        //              << "Hessian inverse [y]:" << std::endl
        //              << hess_inv[1](0, 0) << " " << hess_inv[1](0, 1)
        //              << std::endl
        //              << hess_inv[1](1, 0) << " " << hess_inv[1](1, 1)
        //              << std::endl
        //              << std::endl;
        // jac_and_hess.close();

        // // [zdec] debug
        // std::ofstream debug_stream;
        // debug_stream.open("norm_and_tan.dat", std::ios_base::app);
        // debug_stream << x[0] << " " << x[1] << " " << ni[0] << " " << ni[1]
        //              << " " << ti[0] << " " << ti[1] << " " << dnids[0] << "
        //              "
        //              << dnids[1] << " " << dtids[0] << " " << dtids[1] << " "
        //              << d2xids2[0] << " " << d2xids2[1] << std::endl;
        // debug_stream.close();

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


  //============================================================================
  /// FoepplVonKarmanC1CurvableBellElement elements are a subparametric scheme
  /// with  linear Lagrange interpolation for approximating the geometry and
  /// the C1-functions for approximating variables.
  /// These elements use TElement<2,NNODE_1D> with Lagrangian bases to
  /// interpolate the in-plane unknowns u_x, u_y although with sub-parametric
  /// simplex shape interpolation.

  /// The Bell Hermite basis is used to interpolate the out-of-plane unknown w
  /// for straight sided elements and is upgraded to Bernadou basis if the
  /// element is on a curved boundary. This upgrade similarly corrects the shape
  /// interpolation.
  ///
  /// As a result, all nodes contain the two in-plane Lagrangian dofs, vertex
  /// nodes (i.e. nodes 0,1,2) each contain an additional six out-of-plane
  /// Hermite dofs.
  ///
  ///  e.g FeopplVonKarmanC1CurvableBellElement<4> (unupgraded)
  ///                        | 0 | 1 | 2 | 3  | 4  |  5   |  6    |  7   |
  ///          o <---------- |u_x|u_y| w |dwdx|dwdy|d2wdx2|d2wdxdy|d2wdy2|
  ///         / \            .
  ///        x   x <-------- |u_x|u_y|
  ///       /     \          .
  ///      x   x   x         .
  ///     /         \        .
  ///    o---x---x---o       .
  ///
  ///  e.g FeopplVonKarmanC1CurvableBellElement<3> (upgraded)
  ///                         | 0 | 1 | 2 | 3  | 4  |  5   |  6    |  7   |
  ///          o_ <---------- |u_x|u_y| w |dwdn|dwdt|d2wdn2|d2wdndt|d2wdt2|
  ///         /  \            .
  ///        / .  |           .
  ///       x     x <-------- |u_x|u_y|
  ///      / .  . |           .
  ///     /        \          .
  ///    o-----x----o         .
  ///
  //============================================================================

  template<unsigned NNODE_1D>
  class FoepplVonKarmanC1CurvableBellElement
    : public virtual CurvableBellElement<NNODE_1D>,
      public virtual FoepplVonKarmanEquations
  {
  public:

    //----------------------------------------------------------------------
    // Class construction

    /// Constructor: Call constructors for C1CurvableBellElement and
    /// FoepplVonKarmanEquations
    FoepplVonKarmanC1CurvableBellElement()
      : CurvableBellElement<NNODE_1D>(Nfield, Field_is_bell_interpolated),
        FoepplVonKarmanEquations()
    {
      // Use the higher order integration scheme
      delete this->integral_pt();
      // Use the higher order integration scheme
      TGauss<2, 4>* new_integral_pt = new TGauss<2, 4>;
      this->set_integration_scheme(new_integral_pt);

      // Rotated dof helper
      Rotated_boundary_helper_pt = new RotatedBoundaryHelper(this);
    }

    /// Destructor: clean up alloacations
    ~FoepplVonKarmanC1CurvableBellElement()
    {
      // [zdec] dont we need to delete the integral pt?
      delete Rotated_boundary_helper_pt;
    }

    /// Broken copy constructor
    FoepplVonKarmanC1CurvableBellElement(
      const FoepplVonKarmanC1CurvableBellElement<NNODE_1D>& dummy)
    {
      BrokenCopy::broken_copy("FoepplVonKarmanC1CurvableBellElement");
    }

    /// Broken assignment operator
    void operator=(const FoepplVonKarmanC1CurvableBellElement<NNODE_1D>&)
    {
      BrokenCopy::broken_assign("FoepplVonKarmanC1CurvableBellElement");
    }


    //----------------------------------------------------------------------
    // Output and documentation

    /// Output function:
    ///  x, y, ux, uy, w
    void output(std::ostream& outfile)
    {
      FoepplVonKarmanEquations::output(outfile);
    }

    /// Full output function with a rich set of unknowns:
    ///  x, y, ux, uy, w, dw, ddw, du, strain, stress, principal stress
    void full_output(std::ostream& outfile)
    {
      FoepplVonKarmanEquations::full_output(outfile);
    }

    /// Output function:
    ///   x,y,u   or    x,y,z,u at n_plot*(n_plot+1)/2 plot points
    void output(std::ostream& outfile, const unsigned& n_plot)
    {
      FoepplVonKarmanEquations::output(outfile, n_plot);
    }

    /// Output function:
    ///   x,y,u   or    x,y,z,u at n_plot*(n_plot+1)/2plot points
    void output_interpolated_exact_soln(
      std::ostream& outfile,
      FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
      const unsigned& n_plot);

    /// C-style output function:
    ///  x,y,u   or    x,y,z,u
    void output(FILE* file_pt)
    {
      FoepplVonKarmanEquations::output(file_pt);
    }

    /// C-style output function:
    ///   x,y,u   or    x,y,z,u at n_plot*(n_plot+1)/2 plot points
    void output(FILE* file_pt, const unsigned& n_plot)
    {
      FoepplVonKarmanEquations::output(file_pt, n_plot);
    }


    /// Output function for an exact solution:
    ///  x,y,u_exact   or    x,y,z,u_exact at n_plot*(n_plot+1)/2 plot points
    void output_fct(std::ostream& outfile,
                    const unsigned& n_plot,
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
    {
      FoepplVonKarmanEquations::output_fct(outfile, n_plot, exact_soln_pt);
    }

    /// Output function for a time-dependent exact solution.
    ///  x,y,u_exact   or    x,y,z,u_exact at n_plot*(n_plot+1)/2 points
    /// (Calls the steady version)
    void output_fct(std::ostream& outfile,
                    const unsigned& n_plot,
                    const double& time,
                    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
    {
      FoepplVonKarmanEquations::output_fct(
        outfile, n_plot, time, exact_soln_pt);
    }


    //----------------------------------------------------------------------
    // Jacobian and residual contributions

    /// Add the element's contribution to its residual vector (wrapper) with
    /// cached association matrix
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Store the expensive-to-construct matrix
      this->store_association_matrix();
      // Call the generic routine with the flag set to 1
      FoepplVonKarmanEquations::fill_in_contribution_to_residuals(residuals);
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
      FoepplVonKarmanEquations::fill_in_contribution_to_jacobian(residuals,
                                                                 jacobian);
      // Remove the expensive-to-construct matrix
      this->delete_association_matrix();
    }


    //----------------------------------------------------------------------
    // Geometry and boundaries

    /// Function to pin all deflection dofs
    void pin_all_deflection_dofs() const;

    /// Function to pin the j-th in-plane displacement dof at all nodes along
    /// boundary b to the value prescribed by specified_u_j_pt
    void fix_in_plane_displacement_dof(const unsigned& j_type,
                                       const unsigned& b,
                                       const ScalarFctPt& specified_u_j_pt);

    /// Function to pin the j-th out-of-plane displacement dof at all nodes
    /// along boundary b to the value prescribed by specified_w_j_pt
    void fix_out_of_plane_displacement_dof(const unsigned& dof_number,
                                           const unsigned& b,
                                           const ScalarFctPt& specified_w_j_pt);


    // [zdec] I think i misnamed this, should it be interpolated_x?
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
      return CurvableBellElement<NNODE_1D>::element_is_curved();
    }

    /// Upgrade the Bell element to a curved Bernadou element. Expects, in
    /// order, the unsigned enumeration of the edge that is on the boundary
    /// (curved_edge) as well as the coordinates of the start and end of that
    /// edge on the boundary (s_ubar,s_obar), the parametric description of the
    /// curved edge (parametric_edge) and lastly the polynomial order of the
    /// boundary interpolation (boundary_order) which can be either 3 or 5.
    virtual void upgrade_element_to_curved(
      const MyC1CurvedElements::Edge& curved_edge,
      const double& s_ubar,
      const double& s_obar,
      CurvilineGeomObject* parametric_edge,
      const unsigned& boundary_order)
    {
      CurvableBellElement<NNODE_1D>::upgrade_element_to_curved(
        curved_edge, s_ubar, s_obar, parametric_edge, boundary_order);
    }


    //----------------------------------------------------------------------
    // Member data access functions

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

    /// Access the status of the interpolation of each field
    std::vector<bool> field_is_bell_interpolated()
    {
      return Field_is_bell_interpolated;
    }

    /// Required  # of `values' (pinned or dofs)
    /// at node n
    inline unsigned required_nvalue(const unsigned& n) const
    {
      return Initial_Nvalue[n];
    }



  protected:
    //----------------------------------------------------------------------------
    // Interface to FoepplVonKarmanEquations (can this all be (static) data?)

    /// Function to return a vector of the indices of the in-plane fvk
    /// displacement unkonwns in the grander scheme of unknowns
    virtual Vector<unsigned> u_field_indices() const
    {
      return {0, 1};
    }

    /// Function to return the index of the alpha-th in-plane displacement
    /// unkonwns in the grander scheme of unknowns
    virtual unsigned u_alpha_field_index(const unsigned& alpha) const
    {
      return alpha;
    }

    /// Field index for w
    virtual unsigned w_field_index() const
    {
      return 2;
    }

    /// Interface to return the number of nodes used by u
    virtual unsigned nu_node() const
    {
      return CurvableBellElement<NNODE_1D>::nnode_for_field(
        u_alpha_field_index(0));
    }

    /// Interface to return the number of nodes used by w
    virtual unsigned nw_node() const
    {
      return CurvableBellElement<NNODE_1D>::nnode_for_field(w_field_index());
    }


    /// Interface to get the local indices of the nodes used by u
    virtual Vector<unsigned> get_u_node_indices() const
    {
      return CurvableBellElement<NNODE_1D>::nodal_indices_for_field(
        u_alpha_field_index(0));
    }

    /// Interface to get the local indices of the nodes used by w
    virtual Vector<unsigned> get_w_node_indices() const
    {
      return CurvableBellElement<NNODE_1D>::nodal_indices_for_field(
        w_field_index());
    }


    /// Interface to get the number of basis types for u at node j
    virtual unsigned nu_type_at_each_node() const
    {
      return CurvableBellElement<NNODE_1D>::nnodal_basis_type_for_field(
        u_alpha_field_index(0));
    }

    /// Interface to get the number of basis types for w at node j
    virtual unsigned nw_type_at_each_node() const
    {
      return CurvableBellElement<NNODE_1D>::nnodal_basis_type_for_field(
        w_field_index());
    }

    /// Interface to retrieve the value of u_alpha at node j of
    /// type k
    virtual double get_u_alpha_value_at_node_of_type(
      const unsigned& alpha,
      const unsigned& j_node,
      const unsigned& k_type) const
    {
      unsigned n_u_types = nu_type_at_each_node();
      unsigned nodal_type_index = alpha * n_u_types + k_type;
      return raw_nodal_value(j_node, nodal_type_index);
    }

    /// Interface to retrieve the t-th history value of
    /// u_alpha at node j of type k
    virtual double get_u_alpha_value_at_node_of_type(
      const unsigned& t_time,
      const unsigned& alpha,
      const unsigned& j_node,
      const unsigned& k_type) const
    {
      unsigned n_u_types = nu_type_at_each_node();
      unsigned nodal_type_index = alpha * n_u_types + k_type;
      return raw_nodal_value(t_time, j_node, nodal_type_index);
    }

    /// Interface to retrieve the value of w at node j of type k
    virtual double get_w_value_at_node_of_type(const unsigned& j_node,
                                               const unsigned& k_type) const
    {
      unsigned n_u_fields = 2;
      unsigned n_u_types = nu_type_at_each_node();
      unsigned nodal_type_index = n_u_fields * n_u_types + k_type;
      return raw_nodal_value(j_node, nodal_type_index);
    }

    /// Interface to retrieve the t-th history value of w at node j
    /// of type k
    virtual double get_w_value_at_node_of_type(const unsigned& t_time,
                                               const unsigned& j_node,
                                               const unsigned& k_type) const
    {
      unsigned n_u_fields = 2;
      unsigned n_u_types = nu_type_at_each_node();
      unsigned nodal_type_index = n_u_fields * n_u_types + k_type;
      return raw_nodal_value(t_time, j_node, nodal_type_index);
    }

    // [IN-PLANE-INTERNAL]
    // /// (pure virtual) interface to get the pointer to the internal data used
    // to
    // /// interpolate u (NOTE: assumes each u field has exactly one internal
    // data) virtual Vector<Data*> u_internal_data_pts()
    //  {
    //   // Write me
    //  }

    /// Interface to get the pointer to the internal data used to
    /// interpolate w (NOTE: assumes w field has exactly one internal data)
    virtual Data* w_internal_data_pt() const
    {
      unsigned i_field = w_field_index();
      unsigned index =
        CurvableBellElement<NNODE_1D>::index_of_internal_data_for_field(
          i_field);
      return internal_data_pt(index);
    }


    /// Interface to get the number of internal types for the u
    /// fields (left in tact to ensure that this always returns zero)
    virtual unsigned nu_type_internal() const
    {
      return CurvableBellElement<NNODE_1D>::ninternal_basis_type_for_field(
        u_field_indices()[0]);
    }

    /// Interface to get the number of internal types for the w
    /// fields
    virtual unsigned nw_type_internal() const
    {
      // // [zdec] debug
      // oomph_info
      //   << "We have "
      //   << CurvableBellElement<NNODE_1D>::ninternal_basis_type_for_field(
      // 	  w_field_index())
      //   << " in here." << std::endl;

      return CurvableBellElement<NNODE_1D>::ninternal_basis_type_for_field(
        w_field_index());
    }


    // [IN-PLANE-INTERNAL]
    // /// (pure virtual) interface to retrieve the value of u_alpha of internal
    // /// type k
    // virtual double get_u_alpha_internal_value_of_type(const unsigned& alpha,
    // 						    const unsigned& k_type) const = 0;

    // /// (pure virtual) interface to retrieve the t-th history value of
    // u_alpha of
    // /// internal type k
    // virtual double get_u_alpha_internal_value_of_type(const unsigned& time,
    // 						    const unsigned& alpha,
    // 						    const unsigned& k_type) const = 0;


    /// Interface to retrieve the value of w of internal type k
    virtual double get_w_internal_value_of_type(const unsigned& k_type) const
    {
      unsigned index = w_field_index();
      return CurvableBellElement<NNODE_1D>::internal_value_for_field_of_type(
        index, k_type);
    }

    /// Interface to retrieve the t-th history value of w of
    /// internal type k
    virtual double get_w_internal_value_of_type(const unsigned& t_time,
                                                const unsigned& k_type) const
    {
      unsigned index = w_field_index();
      return CurvableBellElement<NNODE_1D>::internal_value_for_field_of_type(
        t_time, index, k_type);
    }


    /// In-plane basis functions and derivatives w.r.t. global
    /// coords at local coordinate s; return det(Jacobian of mapping)
    virtual double basis_u_foeppl_von_karman(const Vector<double>& s,
                                             Shape& psi_n) const;

    /// In-plane basis functions and derivatives w.r.t. global
    /// coords at local coordinate s; return det(Jacobian of mapping)
    virtual double dbasis_u_eulerian_foeppl_von_karman(const Vector<double>& s,
                                                       Shape& psi_n,
                                                       DShape& dpsi_n_dx) const;

    /// In-plane basis/test functions at and derivatives w.r.t
    /// global coords at local coordinate s; return det(Jacobian of mapping)
    virtual double dbasis_and_dtest_u_eulerian_foeppl_von_karman(
      const Vector<double>& s,
      Shape& psi_n,
      DShape& dpsi_n_dx,
      Shape& test_n,
      DShape& dtest_n_dx) const;

    /// Out-of-plane basis functions at local coordinate s
    virtual void basis_w_foeppl_von_karman(const Vector<double>& s,
                                           Shape& psi_n,
                                           Shape& psi_i) const;

    /// Out-of-plane basis functions and derivs w.r.t. global
    /// coords at local coordinate s; return det(Jacobian of mapping)
    virtual double d2basis_w_eulerian_foeppl_von_karman(
      const Vector<double>& s,
      Shape& psi_n,
      Shape& psi_i,
      DShape& dpsi_n_dx,
      DShape& dpsi_i_dx,
      DShape& d2psi_n_dx2,
      DShape& d2psi_i_dx2) const;

    /// Out-of-plane basis/test functions at local coordinate s
    virtual void basis_and_test_w_foeppl_von_karman(const Vector<double>& s,
                                                    Shape& psi_n,
                                                    Shape& psi_i,
                                                    Shape& test_n,
                                                    Shape& test_i) const;

    /// Out-of-plane basis/test functions and first derivs w.r.t.
    /// to global coords at local coordinate s; return det(Jacobian of mapping)
    virtual double dbasis_and_dtest_w_eulerian_foeppl_von_karman(
      const Vector<double>& s,
      Shape& psi_n,
      Shape& psi_i,
      DShape& dpsi_n_dx,
      DShape& dpsi_i_dx,
      Shape& test_n,
      Shape& test_i,
      DShape& dtest_n_dx,
      DShape& dtest_i_dx) const;

    /// Out-of-plane basis/test functions and first/second derivs
    /// w.r.t. to global coords at local coordinate s;
    /// return det(Jacobian of mapping)
    virtual double d2basis_and_d2test_w_eulerian_foeppl_von_karman(
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

    // End of FoepplVonKarmanEquations interface functions
    //----------------------------------------------------------------------------


    //----------------------------------------------------------------------
    // Rotation handling helpers

    /// Transform the shape functions so that they correspond to
    /// the new rotated dofs
    inline void rotate_shape(Shape& shape) const;

    /// Transform the shape functions and first derivatives so that they
    /// correspond to the new rotated dofs
    inline void rotate_shape(Shape& shape, DShape& dshape) const;

    /// Transform the shape functions, first and second   derivatives so
    /// that they correspond to the new rotated dofs
    inline void rotate_shape(Shape& shape,
                             DShape& dshape,
                             DShape& d2shape) const;



    // All member data is private
  private:

    /// Pointer to an instance of rotated boundary helper
    RotatedBoundaryHelper* Rotated_boundary_helper_pt;

    /// Static number of fields (is always 3)
    static const unsigned Nfield;

    /// Static bool vector with the Bell interpolation of the fields
    /// (always only w)
    static const std::vector<bool> Field_is_bell_interpolated;

    /// Array of static ints that holds the number of variables at
    /// nodes: constant across elements but varies across nodes
    static const unsigned Initial_Nvalue[];

  };


  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////



  //==============================================================================
  /// Face geometry for the FoepplVonKarmanC1CurvableBellElement elements: The
  /// spatial dimension of the face elements is one lower than that of the bulk
  /// element but they have the same number of points along their 1D edges.
  //==============================================================================
  template<unsigned NNODE_1D>
  class FaceGeometry<FoepplVonKarmanC1CurvableBellElement<NNODE_1D>>
    : public virtual TElement<1, NNODE_1D>
  {
  public:
    /// Constructor: Call the constructor for the
    /// appropriate lower-dimensional TElement
    FaceGeometry() : TElement<1, NNODE_1D>() {}
  };



  // //[zdec] old rotation
  // //==============================================================================
  // /// Set up the rotated degrees of freedom: includes a check for the number
  // of
  // /// rotation nodes being greater than three.
  // //==============================================================================
  // template<unsigned NNODE_1D>
  // void FoepplVonKarmanC1CurvableBellElement<NNODE_1D>::set_up_rotated_dofs(
  //   const unsigned& nnodes_to_rotate,
  //   const Vector<unsigned>& nodes_to_rotate,
  //   const BasisVectorsFctPt& basis_vectors_fct_pt)
  // {
  //   // Change the member Nnode_to_rotate
  //   Nnodes_to_rotate = nnodes_to_rotate;
  //   #ifdef PARANOID
  //   // Check that the number of nodes is smaller than 3
  //   if (nnodes_to_rotate > 3)
  //   {
  //     throw OomphLibError(
  // 	"There are only three nodes per element, so we cannot rotate more than
  // three ", 	OOMPH_CURRENT_FUNCTION, 	OOMPH_EXCEPTION_LOCATION);
  //   }
  //   #endif

  //   Nodes_to_rotate = nodes_to_rotate;

  //   // Point to the basis vectors function
  //   Rotated_basis_fct_pt = basis_vectors_fct_pt;
  // }


  // // [zdec] old rotation -- delete
  // //==============================================================================
  // /// Rotate the shape functions according to
  // /// w.r.t. global coordinates and return Jacobian of mapping.
  // ///
  // /// Galerkin: Test functions = shape functions
  // //==============================================================================
  // template<unsigned NNODE_1D>
  // void
  // FoepplVonKarmanC1CurvableBellElement<NNODE_1D>::rotation_matrix_at_node(
  //   const unsigned& inode, DenseDoubleMatrix& rotation_matrix) const
  // {
  //   // Initialise x normal and tangent
  //   Vector<double> x(2, 0.0);

  //   // Get the node pointer
  //   Node* nod_pt = this->node_pt(inode);

  //   // Get the position of the vertex
  //   x[0] = nod_pt->x(0);
  //   x[1] = nod_pt->x(1);

  //   // Initialise the two basis vectors
  //   Vector<Vector<double>> bi(2, Vector<double>(2, 0.0));
  //   Vector<DenseMatrix<double>> dbi(2, DenseMatrix<double>(2, 2, 0.0));

  //   (*Rotated_basis_fct_pt)(x,bi[0],bi[1],dbi[0],dbi[1]);

  //   // Rotation matrix, B
  //   DenseMatrix<double> b1(2, 2, 0.0), b22(3, 3, 0.0), b21(3, 2, 0.0);

  //   // Fill in the submatrices
  //   for (unsigned alpha = 0; alpha < 2; ++alpha)
  //   {
  //     for (unsigned beta = 0; beta < 2; ++beta)
  //     {
  //       // Fill in b1 - the Jacobian
  //       // Fill in the rotation of the first derivatives
  //       b1(alpha, beta) = bi[beta][alpha];

  //       // Avoid double counting the cross derivative
  //       if (alpha <= beta)
  //       {
  //         // Define row index
  //         const unsigned row = alpha + beta;
  //         for (unsigned gamma = 0; gamma < 2; ++gamma)
  //         {
  //           // Fill in b21 - the non affine part of the Jacobian derivative
  //           // Define column index
  //           unsigned col_b21 = gamma;
  //           // Fill in the non-affine part of the rotation of the second
  //           // derivatives if( beta>= alpha) ?
  //           b21(row, col_b21) += dbi[gamma](alpha, beta);
  //           for (unsigned delta = 0; delta < 2; ++delta)
  //           {
  //             // Fill in b22 - the Affine part of the Jacobian derivative
  //             // Redefine column index for the next submatrix
  //             unsigned col_b22 = gamma + delta;
  //             // Fill in the affine part of the rotation of the second
  //             // derivatives if( beta>= alpha) ?
  //             b22(row, col_b22) += bi[gamma][alpha] * bi[delta][beta];
  //           }
  //         }
  //       }
  //     }
  //   }


  //   // Fill in the submatrices to the full (6x6) matrix - we need to right
  //   // multiply this matrix so we need the transpose of the Jacobian W dof
  //   // remains the same
  //   rotation_matrix(0, 0) = 1.0;
  //   // Fill in b1
  //   for (unsigned i = 0; i < 2; ++i)
  //   {
  //     for (unsigned j = 0; j < 2; ++j)
  //     {
  //       rotation_matrix(1 + j, 1 + i) = b1(i, j);
  //     }
  //   }
  //   // Fill in b21
  //   for (unsigned i = 0; i < 3; ++i)
  //   {
  //     for (unsigned j = 0; j < 2; ++j)
  //     {
  //       rotation_matrix(1 + j, 3 + i) = b21(i, j);
  //     }
  //   }
  //   // Fill in b22
  //   for (unsigned i = 0; i < 3; ++i)
  //   {
  //     for (unsigned j = 0; j < 3; ++j)
  //     {
  //       rotation_matrix(3 + j, 3 + i) = b22(i, j);
  //     }
  //   }
  // }


  //============================================================================
  /// In-plane basis functions and derivatives w.r.t. global
  /// coords at local coordinate s; return det(Jacobian of mapping)
  //============================================================================
  template<unsigned NNODE_1D>
  double FoepplVonKarmanC1CurvableBellElement<
    NNODE_1D>::basis_u_foeppl_von_karman(const Vector<double>& s,
                                         Shape& psi_n) const
  {
    // Dimension
    unsigned dim = this->dim();

    // Initialise and get dpsi w.r.t local coord
    const unsigned n_node = this->nnode();
    DShape dummy_dpsids(n_node, dim);
    SubparametricTriangleElement<NNODE_1D>::
      dshape_local(s, psi_n, dummy_dpsids);
    double J;

    // Get the Jacobian of the mapping
    DenseMatrix<double> jacobian(dim, dim);
    DenseMatrix<double> inverse_jacobian(dim, dim);
    J = CurvableBellElement<NNODE_1D>::local_to_eulerian_mapping(
      s, jacobian, inverse_jacobian);

    return J;
  }


  //============================================================================
  /// Fetch the in-plane basis functions and their derivatives w.r.t. global
  /// coordinates at s and return Jacobian of mapping.
  //============================================================================
  template<unsigned NNODE_1D>
  double FoepplVonKarmanC1CurvableBellElement<
    NNODE_1D>::dbasis_u_eulerian_foeppl_von_karman(const Vector<double>& s,
                                                   Shape& psi_n,
                                                   DShape& dpsi_n_dx) const
  {
    // Dimension
    unsigned dim = this->dim();

    // Initialise and get dpsi w.r.t local coord
    const unsigned n_node = this->nnode();
    DShape dpsids(n_node, dim);
    SubparametricTriangleElement<NNODE_1D>::dshape_local(s, psi_n, dpsids);
    double J;

    // Get the Jacobian of the mapping
    DenseMatrix<double> jacobian(dim, dim);
    DenseMatrix<double> inverse_jacobian(dim, dim);
    J = CurvableBellElement<NNODE_1D>::local_to_eulerian_mapping(
      s, jacobian, inverse_jacobian);

    // Now find the global derivatives
    for (unsigned l = 0; l < n_node; ++l)
    {
      for (unsigned i = 0; i < this->dim(); ++i)
      {
        // Initialise to zero
        dpsi_n_dx(l, i) = 0.0;
        for (unsigned j = 0; j < this->dim(); ++j)
        {
          // Convert to local coordinates
          dpsi_n_dx(l, i) += inverse_jacobian(i, j) * dpsids(l, j);
        }
      }
    }
    return J;
  }


  //============================================================================
  /// Fetch the in-plane basis functions and test functions and derivatives
  /// w.r.t. global coordinates at s and return Jacobian of mapping.
  ///
  /// Galerkin: Test functions = shape functions
  //=============================================================================
  template<unsigned NNODE_1D>
  double FoepplVonKarmanC1CurvableBellElement<NNODE_1D>::
    dbasis_and_dtest_u_eulerian_foeppl_von_karman(const Vector<double>& s,
                                                  Shape& psi_n,
                                                  DShape& dpsi_n_dx,
                                                  Shape& test_n,
                                                  DShape& dtest_n_dx) const
  {
    // Initialise and get dpsi w.r.t local coord
    const unsigned n_node = this->nnode();
    DShape dpsids(n_node, this->dim());
    SubparametricTriangleElement<NNODE_1D>::dshape_local(s, psi_n, dpsids);
    double J;

    // Get the Jacobian of the mapping
    DenseMatrix<double> jacobian(this->dim(), this->dim()),
      inverse_jacobian(this->dim(), this->dim());
    J = CurvableBellElement<NNODE_1D>::local_to_eulerian_mapping(
      s, jacobian, inverse_jacobian);

    // Now find the global derivatives
    for (unsigned l = 0; l < n_node; ++l)
    {
      for (unsigned i = 0; i < this->dim(); ++i)
      {
        // Initialise to zero
        dpsi_n_dx(l, i) = 0.0;
        for (unsigned j = 0; j < this->dim(); ++j)
        {
          // Convert to local coordinates
          dpsi_n_dx(l, i) += inverse_jacobian(i, j) * dpsids(l, j);
        }
      }
    }
    test_n = psi_n;
    dtest_n_dx = dpsi_n_dx;
    return J;
  }


  //======================================================================
  /// (pure virtual) Out-of-plane basis functions at local coordinate s
  //======================================================================
  template<unsigned NNODE_1D>
  void FoepplVonKarmanC1CurvableBellElement<
    NNODE_1D>::basis_w_foeppl_von_karman(const Vector<double>& s,
                                         Shape& psi_n,
                                         Shape& psi_i) const
  {
    throw OomphLibError("This still needs testing for curved elements.",
                        "void FoepplVonKarmanC1CurvableBellElement<NNODE_1D>::\
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
  template<unsigned NNODE_1D>
  void FoepplVonKarmanC1CurvableBellElement<
    NNODE_1D>::basis_and_test_w_foeppl_von_karman(const Vector<double>& s,
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
  template<unsigned NNODE_1D>
  double FoepplVonKarmanC1CurvableBellElement<NNODE_1D>::
    dbasis_and_dtest_w_eulerian_foeppl_von_karman(const Vector<double>& s,
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
  template<unsigned NNODE_1D>
  double FoepplVonKarmanC1CurvableBellElement<
    NNODE_1D>::d2basis_w_eulerian_foeppl_von_karman(const Vector<double>& s,
                                                    Shape& psi_n,
                                                    Shape& psi_i,
                                                    DShape& dpsi_n_dx,
                                                    DShape& dpsi_i_dx,
                                                    DShape& d2psi_n_dx,
                                                    DShape& d2psi_i_dx) const
  {
    // Call the geometrical shape functions and derivatives
    double J = CurvableBellElement<NNODE_1D>::d2_c1_basis_eulerian(
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
  template<unsigned NNODE_1D>
  double FoepplVonKarmanC1CurvableBellElement<NNODE_1D>::
    d2basis_and_d2test_w_eulerian_foeppl_von_karman(const Vector<double>& s,
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
    double J = CurvableBellElement<NNODE_1D>::d2_c1_basis_eulerian(
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
  template<unsigned NNODE_1D>
  inline void FoepplVonKarmanC1CurvableBellElement<NNODE_1D>::rotate_shape(
    Shape& psi) const
  {
    const unsigned n_dof_types = nw_type_at_each_node();

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
  template<unsigned NNODE_1D>
  inline void FoepplVonKarmanC1CurvableBellElement<NNODE_1D>::rotate_shape(
    Shape& psi, DShape& dpsidx) const
  {
    const unsigned n_dof_types = nw_type_at_each_node();
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
  template<unsigned NNODE_1D>
  inline void FoepplVonKarmanC1CurvableBellElement<NNODE_1D>::rotate_shape(
    Shape& psi, DShape& dpsidx, DShape& d2psidx) const
  {
    const unsigned n_dof_types = nw_type_at_each_node();
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


  //=============================================================================
  /// Function to pin all deflection dofs
  //=============================================================================
  template<unsigned NNODE_1D>
  void FoepplVonKarmanC1CurvableBellElement<NNODE_1D>::pin_all_deflection_dofs()
    const
  {
    const unsigned w_index = w_field_index();

    // Get nodes
    const unsigned n_node = nw_node();
    const Vector<unsigned> nodes = get_w_node_indices();

    // Curved Bell elements only have deflection dofs at vertices
    for (unsigned j_nodei = 0; j_nodei < n_node; j_nodei++)
    {
      // Get the j_nodei-th node used by i_field
      unsigned j_node = nodes[j_nodei];
      Node* nod_pt = this->node_pt(j_node);

      // Get the number of types at the current node
      unsigned n_type = nw_type_at_each_node();

      // Check if it is on the boundary
      for (unsigned k_type = 0; k_type < n_type; k_type++)
      {
        // Pin and set the value
        nod_pt->pin(2 + k_type);
        nod_pt->set_value(2 + k_type, 0.0);
      }
    }

    // Now fix internal dofs
    unsigned n_internal = nw_type_internal();
    for (unsigned k_type = 0; k_type < n_internal; k_type++)
    {
      // Get node
      // Pin and set the value
      this->internal_data_for_field_pt(w_index)->pin(k_type);
    }
  }


  //=============================================================================
  /// Function to pin the j-th in-plane displacement dof at all nodes along
  /// boundary b to the value prescribed by specified_u_j_pt
  //=============================================================================
  template<unsigned NNODE_1D>
  void FoepplVonKarmanC1CurvableBellElement<NNODE_1D>::
  fix_in_plane_displacement_dof(const unsigned& alpha,
				const unsigned& b,
				const ScalarFctPt& specified_u_j_pt)
  {
    // Initialise constants that we use in this function
    const unsigned field_index = u_alpha_field_index(alpha);
    const unsigned nodal_type_index =
      this->first_nodal_type_index_for_field(field_index);
    const unsigned n_node = nu_node();
    const unsigned dim = this->dim();

#ifdef PARANOID
    // Check that the dof number is a sensible value
    const unsigned n_type = 2 * nu_type_at_each_node();
    if (alpha >= n_type)
    {
      throw OomphLibError(
        "Foppl von Karman elements only have 2 in-plane displacement degrees\
of freedom at internal points. They are {ux, uy}",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Bell elements only have deflection dofs at vertices
    for (unsigned n = 0; n < n_node; ++n)
    {
      // Get node
      Node* nod_pt = this->node_pt(n);
      // Check if it is on the boundary
      bool is_boundary_node = nod_pt->is_on_boundary(b);
      if (is_boundary_node)
      {
        // Extract nodal coordinates from node:
        // Since the element isn't necessarily isoparametric the nodes
        // 'position' is not necessarily correct?
        Vector<double> x(dim), s(dim);
        this->local_coordinate_of_node(n, s);
        interpolated_x(s, x);
        // Fill in value
        double value;
        specified_u_j_pt(x, value);
        // Pin and set the value
        nod_pt->pin(alpha);
        nod_pt->set_value(nodal_type_index, value);
      }
    }
  }


  //=============================================================================
  /// Function to pin particular out-of-plane displacement dofs along boundary
  /// b to the value prescribed by specified_w_j_pt
  //=============================================================================
  template<unsigned NNODE_1D>
  void FoepplVonKarmanC1CurvableBellElement<NNODE_1D>::
  fix_out_of_plane_displacement_dof(const unsigned& k_type,
				    const unsigned& b_boundary,
				    const ScalarFctPt& specified_w_j_pt)
  {
    const unsigned w_index = w_field_index();
    const unsigned first_nodal_type_index =
      this->first_nodal_type_index_for_field(w_index);
    const unsigned n_vertices = nw_node();

#ifdef PARANOID
    // Check that the dof number is a sensible value
    unsigned n_type = nw_type_at_each_node();
    if (k_type >= n_type)
    {
      throw OomphLibError(
        "Foppl von Karman elements only have 6 Hermite deflection degrees\
of freedom at internal points. They are {w ; w,x ; w,y ; w,xx ; w,xy ; w,yy}",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Bell elements only have deflection dofs at vertices
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
        specified_w_j_pt(x, value);
        // Pin and set the value
        nod_pt->pin(first_nodal_type_index + k_type);
        nod_pt->set_value(first_nodal_type_index + k_type, value);
      }
    }
  }


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
      Index_of_lagrange_data = add_internal_data(new Data(8));

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
      fill_in_generic_residual_contribution_constraint(residuals, jacobian, 1);
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

      // We use a vector of booleans to keep track of dofs that might be reused
      // (no need to track right dofs which are used once)
      std::vector<bool> right_data_used(8, false);
      std::vector<bool> left_data_used(8, false);

      // Store each data
      Data* left_data_pt = external_data_pt(Index_of_left_data);
      Data* right_data_pt = external_data_pt(Index_of_right_data);

      // We also want to store the jacobian and the hessian of the mapping
      DenseMatrix<double> jac_of_transform(2, 2, 0.0);
      Vector<DenseMatrix<double>> hess_of_transform(
        2, DenseMatrix<double>(2, 2, 0.0));
      get_jac_and_hess_of_coordinate_transform(jac_of_transform,
                                               hess_of_transform);

      // Constraints 0-2 use dofs 0-2 respectively in each node
      for (unsigned i_con = 0; i_con < 3; i_con++)
      {
        // Get whether each value is pinned
        bool right_ui_pinned = right_data_pt->is_pinned(i_con);
        bool left_ui_pinned = left_data_pt->is_pinned(i_con);
        // If anything is free, mark it as used and continue without doing
        // anything else
        if (!right_ui_pinned && !right_data_used[i_con])
        {
          // [zdec] debug
          std::cout << "eqn " << i_con << " depends on dof R" << i_con
                    << std::endl;
          right_data_used[i_con] = true;
          continue;
        }
        if (!left_ui_pinned && !left_data_used[i_con])
        {
          // [zdec] debug
          std::cout << "eqn " << i_con << " depends on dof L" << i_con
                    << std::endl;
          left_data_used[i_con] = true;
          continue;
        }
        // ---------------------------------------------------------------------
        // If we made it here, it is because all dofs in the constraint are
        // pinned so we need to check the constraint is satisfied manually and
        // then remove it by pinning the corresponding lagrange multiplier

        // // Calculate the residual of the constraint
        // double constraint_residual =
        //   right_data_pt->value(i_con) - left_data_pt->value(i_con);
        // // Check that the constraint is met and we don't have a tear
        // if(constraint_residual > Constraint_tolerance)
        // {
        //   throw_unsatisfiable_constraint_error(i_con, constraint_residual);
        // }

        // If it is met, we pin the lagrange multiplier that corresponds to
        // this constraint as it is redundant and results in a zero row/column
        internal_data_pt(Index_of_lagrange_data)->pin(i_con);
      }

      // Constraints 3-4 use dofs 3-4 respectively from the right node and
      // both in the left
      for (unsigned alpha = 0; alpha < 2; alpha++)
      {
        // Get the condition associated with this right derivative
        unsigned i_con = 3 + alpha;

        // Get whether each value is pinned
        bool right_ui_pinned = right_data_pt->is_pinned(i_con);
        bool left_u3_pinned = left_data_pt->is_pinned(3);
        bool left_u4_pinned = left_data_pt->is_pinned(4);
        // If anything is free, mark it as used and continue without doing
        // anything else. We also need to check that each dof hasn't become
        // decoupled from this constraint by ensuring that its coefficient (if
        // it has one) is sufficiently large (> Orthogonality_tolerance)
        if (!right_ui_pinned && !right_data_used[i_con])
        {
          // [zdec] debug
          std::cout << "eqn " << i_con << " depends on dof R" << i_con
                    << std::endl;
          right_data_used[i_con] = true;
          continue;
        }
        if (!left_u3_pinned && !left_data_used[3])
        {
          // [zdec] debug
          std::cout << "eqn " << i_con << " depends on dof L3" << std::endl;
          double coeff = jac_of_transform(0, alpha);
          if (fabs(coeff) > Orthogonality_tolerance)
          {
            left_data_used[3] = true;
            continue;
          }
        }
        if (!left_u4_pinned && !left_data_used[4])
        {
          // [zdec] debug
          std::cout << "eqn " << i_con << " depends on dof L4" << std::endl;
          double coeff = jac_of_transform(1, alpha);
          if (fabs(coeff) > Orthogonality_tolerance)
          {
            left_data_used[4] = true;
            continue;
          }
        }
        // ---------------------------------------------------------------------
        // If we made it here, it is because all dofs in the constraint are
        // pinned so we need to check the constraint is satisfied manually and
        // then remove it by pinning the corresponding lagrange multiplier

        // // Calculate the residual of the constraint
        // double constraint_residual = right_data_pt->value(i_con);
        // for(unsigned beta = 0; beta < 2; beta++)
        // {
        //   constraint_residual +=
        //     - left_data_pt->value(3+beta) * jac_of_transform(beta,alpha);
        // }
        // // Check that the constraint is met and we don't have a tear
        // if(constraint_residual > Constraint_tolerance)
        // {
        //   throw_unsatisfiable_constraint_error(i_con, constraint_residual);
        // }

        // If it is met, we pin the lagrange multiplier that corresponds to
        // this constraint as it is redundant and results in a zero row/column
        internal_data_pt(Index_of_lagrange_data)->pin(i_con);
      }

      // Constraints 5-7 use dofs 5-7 respectively from the right node and
      // all use dofs 3-7 from the left node
      for (unsigned alpha = 0; alpha < 2; alpha++)
      {
        // beta>=alpha so we dont double count constraint 6
        for (unsigned beta = alpha; beta < 2; beta++)
        {
          // The index of the constraint
          unsigned i_con = 5 + alpha + beta;

          // Get whether each value is pinned
          bool right_ui_pinned = right_data_pt->is_pinned(i_con);
          bool left_u3_pinned = left_data_pt->is_pinned(3);
          bool left_u4_pinned = left_data_pt->is_pinned(4);
          bool left_u5_pinned = left_data_pt->is_pinned(5);
          bool left_u6_pinned = left_data_pt->is_pinned(6);
          bool left_u7_pinned = left_data_pt->is_pinned(7);
          // If anything is free, mark it as used and continue without doing
          // anything else. We also need to check that each dof hasn't become
          // decoupled from this constraint by ensuring that its coefficient (if
          // it has one) is sufficiently large (> Orthogonality_tolerance)
          if (!right_ui_pinned && !right_data_used[i_con])
          {
            // [zdec] debug
            std::cout << "eqn " << i_con << " depends on dof R" << i_con
                      << std::endl;
            right_data_used[i_con] = true;
            continue;
          }
          if (!left_u3_pinned && !left_data_used[3])
          {
            double coeff = hess_of_transform[0](alpha, beta);
            if (fabs(coeff) > Orthogonality_tolerance)
            {
              // [zdec] debug
              std::cout << "eqn " << i_con << " depends on dof L3" << std::endl;
              left_data_used[3] = true;
              continue;
            }
          }
          if (!left_u4_pinned && !left_data_used[4])
          {
            double coeff = hess_of_transform[1](alpha, beta);
            if (fabs(coeff) > Orthogonality_tolerance)
            {
              // [zdec] debug
              std::cout << "eqn " << i_con << " depends on dof L4" << std::endl;
              left_data_used[4] = true;
              continue;
            }
          }
          if (!left_u5_pinned && !left_data_used[5])
          {
            double coef =
              jac_of_transform(0, alpha) * jac_of_transform(0, beta);
            if (fabs(coef) > Orthogonality_tolerance)
            {
              // [zdec] debug
              std::cout << "eqn " << i_con << " depends on dof L5" << std::endl;
              left_data_used[5] = true;
              continue;
            }
          }
          if (!left_u6_pinned && !left_data_used[6])
          {
            double coef =
              jac_of_transform(0, alpha) * jac_of_transform(1, beta) +
              jac_of_transform(1, alpha) * jac_of_transform(0, beta);
            if (fabs(coef) > Orthogonality_tolerance)
            {
              // [zdec] debug
              std::cout << "eqn " << i_con << " depends on dof L6" << std::endl;
              left_data_used[6] = true;
              continue;
            }
          }
          if (!left_u7_pinned && !left_data_used[7])
          {
            double coef =
              jac_of_transform(1, alpha) * jac_of_transform(1, beta);
            if (fabs(coef) > Orthogonality_tolerance)
            {
              // [zdec] debug
              std::cout << "eqn " << i_con << " depends on dof L7" << std::endl;
              left_data_used[7] = true;
              continue;
            }
          }
          // -------------------------------------------------------------------
          // If we made it here, it is because all dofs in the constraint are
          // pinned so we need to check the constraint is satisfied manually and
          // then remove it by pinning the corresponding lagrange multiplier

          // // Calculate the residual of the constraint
          // double constraint_residual = right_data_pt->value(i_con);
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
          // // Check that the constraint is met and we don't have a tear
          // if(constraint_residual > Constraint_tolerance)
          // {
          //   throw_unsatisfiable_constraint_error(i_con, constraint_residual);
          // }

          // If it is met, we pin the lagrange multiplier that corresponds to
          // this constraint as it is redundant and results in a zero row/column
          internal_data_pt(Index_of_lagrange_data)->pin(i_con);
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
      // [zdec] debug
      std::cout << std::endl
                << std::endl
                << "ADD CONTRIBUTION FROM CONSTRAINTS" << std::endl
                << "=============================================" << std::endl;
      //----------------------------------------------------------------------
      //----------------------------------------------------------------------
      // Calculate Jacobian and Hessian of coordinate transform between
      // each boundary coordinate
      DenseMatrix<double> jac_of_transform(2, 2, 0.0);
      Vector<DenseMatrix<double>> hess_of_transform(
        2, DenseMatrix<double>(2, 2, 0.0));
      get_jac_and_hess_of_coordinate_transform(jac_of_transform,
                                               hess_of_transform);

      //----------------------------------------------------------------------
      //----------------------------------------------------------------------
      // Use the jac and hess of transform to add the residual
      // contributions from the constraint
      // [zdec]::TODO make indexing (alpha,beta,gamma,...) consistent

      // Store the internal data pointer which stores the Lagrange multipliers
      Vector<double> lagrange_value(8, 0.0);
      internal_data_pt(Index_of_lagrange_data)->value(lagrange_value);

      // Store the left and right nodal dofs
      // 0: u_1        4: dw/ds_2
      // 1: u_2        5: d^2w/ds_1^2
      // 2: w          6: d^2w/ds_1ds_2
      // 3: dw/ds_1    7: d^2w/ds_2^2
      Vector<double> left_value(8, 0.0);
      Vector<double> right_value(8, 0.0);
      Left_node_pt->value(left_value);
      Right_node_pt->value(right_value);


      //----------------------------------------------------------------------
      // First the contributions to the right node external equations
      unsigned n_external_type = 8;
      for (unsigned k_type = 0; k_type < n_external_type; k_type++)
      {
        int right_eqn_number = external_local_eqn(Index_of_right_data, k_type);

        // If this dof isn't pinned we add to the residual
        if (right_eqn_number >= 0)
        {
          // Right dof term in the constraint always lambda_i*W_i
          residuals[right_eqn_number] += lagrange_value[k_type];
        }
      } // End for loop adding contributions to right nodal equations


      //----------------------------------------------------------------------
      // Next, the contributions to the left node external equations
      // First three are displacements:
      //     - lambda_i*(U_alpha or W_0)
      for (unsigned i = 0; i < 3; i++)
      {
        // Eqn number is just the index of (u_x,u_y,w) which is the ith dof
        // (i=0,1,2)
        int left_eqn_number = external_local_eqn(Index_of_left_data, i);
        // If this dof isn't pinned we add to the residual
        if (left_eqn_number >= 0)
        {
          residuals[left_eqn_number] += -lagrange_value[i];
        }
      } // End loop adding contribution to the left nodal displacement equations

      // Next two are from gradient of w:
      //     - lambda_{3+\beta} * w_{1+\alpha} * J_{\alpha\beta}
      //     - lambda_{5+\beta+\gamma} * w_{1+\alpha} * H_{\alpha\beta\gamma}
      // gamma>=beta so we don't double count lambda_6 condition
      for (unsigned alpha = 0; alpha < 2; alpha++)
      {
        // Eqn number is the index of the alpha-th derivative of w which is the
        // 3+alpha-th dof (alpha=0,1)
        int left_eqn_number = external_local_eqn(Index_of_left_data, 3 + alpha);
        // If this dof isn't pinned we add to the residual
        if (left_eqn_number >= 0)
        {
          for (unsigned beta = 0; beta < 2; beta++)
          {
            residuals[left_eqn_number] +=
              -lagrange_value[3 + beta] * jac_of_transform(alpha, beta);
            // gamma>=beta so we don't double count the
            // lagrange_value[6] constraint
            for (unsigned gamma = beta; gamma < 2; gamma++)
            {
              residuals[left_eqn_number] +=
                -lagrange_value[5 + beta + gamma] *
                hess_of_transform[alpha](beta, gamma);
            }
          }
        } // End of if unpinned
      } // End loop adding contributions to the left nodal gradient equations

      // Last three are the second derivatives of w (delta>gamma):
      //     - lambda_{5+\gamma+\delta} * w_{3+\alpha+\beta}
      //       * J_{\alpha\gamma} * J_{\beta\delta}
      // Index second derivative (equation) using alpha & beta
      for (unsigned alpha = 0; alpha < 2; alpha++)
      {
        // Note that d^2w/ds_1ds_2 is counted twice in the summation
        for (unsigned beta = 0; beta < 2; beta++)
        {
          int left_eqn_number =
            external_local_eqn(Index_of_left_data, 5 + alpha + beta);
          // If this dof isn't pinned we add to the residual
          if (left_eqn_number >= 0)
          {
            // Index constraint using gamma and delta
            for (unsigned gamma = 0; gamma < 2; gamma++)
            {
              // delta>=gamma so we don't double count the
              // lagrange_value[6] constraint
              for (unsigned delta = gamma; delta < 2; delta++)
              {
                residuals[left_eqn_number] +=
                  -lagrange_value[5 + gamma + delta] *
                  jac_of_transform(alpha, gamma) *
                  jac_of_transform(beta, delta);
              }
            } // End loops over the conditions (gamma,delta)
          } // End if dof isn't pinned
        }
      } // End loops adding contributions to the left nodal curvature equations
        // (alpha,beta)


      //----------------------------------------------------------------------
      // Now add contributions to the internal (lagrange multiplier) equations

      // Storage for tensor products so they can be reused in the Jacobians
      // (upper halves of symmetric 2x2 matrices stored in vectors of length
      //  three: {A11,A12,A22} )
      // Left gradient of w multiplied by the jacobian of the coordiante
      // transform
      Vector<double> DwJ(2, 0.0);
      // Left Hessian of w multiplied in each index by the jacobian of the
      // coordinate transform
      Vector<double> D2wJJ(3, 0.0);
      // Left gradient of w multipied in the first index of the Hessian of the
      // coordinate transform
      Vector<double> DwH(3, 0.0);

      // First three (u,v,w) dofs are equal
      for (unsigned i_dof = 0; i_dof < 3; i_dof++)
      {
        // Get the internal data eqn number for this constraint
        int internal_eqn_number =
          internal_local_eqn(Index_of_lagrange_data, i_dof);
        // If this dof isn't pinned we add to the residual
        if (internal_eqn_number >= 0)
        {
          // || Add constraining residual ||
          residuals[internal_eqn_number] +=
            (right_value[i_dof] - left_value[i_dof]);
          // [zdec] WHAT IF BOTH ARE PINNED?
          // More generally, what if we have no dofs to satisfy our constraints?
          // Do we check (if all dofs are pinned) that a constraint is met,
          // before pinning the relevant dofs
        }
      }

      // Next two (first derivatives of w) are related by
      //     grad_r(w) = grad_l(w)*J
      // where  J is the Jacobian grad_r(left coords)
      for (unsigned alpha = 0; alpha < 2; alpha++)
      {
        // Get the internal data eqn number for this constraint
        int internal_eqn_number =
          internal_local_eqn(Index_of_lagrange_data, 3 + alpha);
        // If this dof isn't pinned we add to the residual
        if (internal_eqn_number >= 0)
        {
          // Calculate grad_l(w)*J
          for (unsigned beta = 0; beta < 2; beta++)
          {
            DwJ[alpha] += left_value[3 + beta] * jac_of_transform(beta, alpha);
          }
          // || Add constraining residual ||
          residuals[internal_eqn_number] +=
            (right_value[3 + alpha] - DwJ[alpha]);
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
          // Get the internal data eqn number for this constraint
          int internal_eqn_number =
            internal_local_eqn(Index_of_lagrange_data, 5 + alpha + beta);
          // If this dof isn't pinned we add to the residual
          if (internal_eqn_number >= 0)
          {
            // Loops to calculate D(D(w))*J*J and D(w)*H
            for (unsigned gamma = 0; gamma < 2; gamma++)
            {
              for (unsigned delta = 0; delta < 2; delta++)
              {
                // Add contribution to D(D(w))*J*J
                D2wJJ[alpha + beta] += left_value[5 + gamma + delta] *
                                       jac_of_transform(gamma, alpha) *
                                       jac_of_transform(delta, beta);
                // // [zdec] debug
                // std::cout << alpha << " "
                // 	  << beta << " "
                // 	  << gamma << " "
                // 	  << delta << ": D2wJJ += "
                // 	  << left_value[5+gamma+delta] << " * "
                // 	  << jac_of_transform(gamma,alpha) << " * "
                // 	  << jac_of_transform(delta,beta) << " = "
                // 	  << D2wJJ[alpha+beta] << std::endl;
              }
              // Add contributions to D(w)*H
              DwH[alpha + beta] +=
                left_value[3 + gamma] * hess_of_transform[gamma](alpha, beta);
              // // [zdec] debug
              // std::cout << alpha << " "
              // 		<< beta << " "
              // 		<< gamma << " "
              // 		<< ": DwH += "
              // 		<< left_value[3+gamma] << " * "
              // 		<< hess_of_transform[gamma](alpha,beta) << " = "
              // 		<< DwH[alpha+beta] << std::endl;
            }
            // // [zdec] debug
            // std::cout << alpha << " "
            // 	      << beta << " "
            // 	      << ": res += " << right_value[5+alpha+beta]
            // 	      << " - " << D2wJJ[alpha+beta]
            // 	      << " - " << DwH[alpha+beta] << std::endl;

            // || Add constraining residual ||
            residuals[internal_eqn_number] +=
              (right_value[5 + alpha + beta] - D2wJJ[alpha + beta] -
               DwH[alpha + beta]);
          } // End if eqn not pinned
        }
      }

      //----------------------------------------------------------------------
      //----------------------------------------------------------------------
      // If flag, then add the jacobian contribution
      //
      if (flag)
      {
        //----------------------------------------------------------------------
        // First the contributions to the right node's equations which is just
        // r_i * L_i
        for (unsigned k_type = 0; k_type < 8; k_type++)
        {
          // Get the rhs equation number
          int right_eqn_number =
            external_local_eqn(Index_of_right_data, k_type);
          // Get the lagrange dof number
          int lagrange_dof_number =
            internal_local_eqn(Index_of_lagrange_data, k_type);
          // If this equation and dof aren't pinned then add their contribution
          if (right_eqn_number >= 0 && lagrange_dof_number >= 0)
          {
            // Add the contribution to the jacobian
            jacobian(right_eqn_number, lagrange_dof_number) += 1.0;
            // And by symmetry, we can add the transpose contribution to the
            // jacobian
            jacobian(lagrange_dof_number, right_eqn_number) += 1.0;
          } // End pinned check
        } // End loop over right node equations

        //----------------------------------------------------------------------
        // Next the contributions to the left node's equations which differ
        // between dof types
        // The first three are just -l_i * L_i
        for (unsigned k_type = 0; k_type < 3; k_type++)
        {
          // Get the lhs equation number
          int left_eqn_number = external_local_eqn(Index_of_left_data, k_type);
          // Get the lagrange dof number
          int lagrange_dof_number =
            internal_local_eqn(Index_of_lagrange_data, k_type);
          // If this equation and dof aren't pinned then add their contribution
          if (left_eqn_number >= 0 && lagrange_dof_number >= 0)
          {
            // Add the contribution to the jacobian
            jacobian(left_eqn_number, lagrange_dof_number) += -1.0;
            // And by symmetry, we can add the transpose contribution to the
            // jacobian
            jacobian(lagrange_dof_number, left_eqn_number) += -1.0;
          }
        } // End of three displacements

        // The next two equations come from the l_\alpha derivative of:
        // -l_\alpha * ( J_{\alpha\beta} * L_{3+\beta}
        //              + H_{\alpha\beta\gamma} * L_{5+\beta+gamma} )
        for (unsigned alpha = 0; alpha < 2; alpha++)
        {
          // Get the index of this equation in the left node
          unsigned left_type = 3 + alpha;
          // Get the local equation number for this equation
          int left_eqn_number =
            external_local_eqn(Index_of_left_data, left_type);
          // Check this equation isn't pinned
          if (left_eqn_number >= 0)
          {
            // Loop over the Lagrange dofs
            for (unsigned beta = 0; beta < 2; beta++)
            {
              // Get the index of this Lagrange dof
              unsigned lagrange_type = 3 + beta;
              // Get the local equation number for this dof
              int lagrange_dof_number =
                internal_local_eqn(Index_of_lagrange_data, lagrange_type);
              // Check this dof isn't pinned
              if (lagrange_dof_number >= 0)
              {
                double jac_term = -jac_of_transform(alpha, beta);
                // Orthogonality check (for jacobian conditioning)
                if (fabs(jac_term) > Orthogonality_tolerance)
                {
                  // Add the contribution to the jacobian
                  jacobian(left_eqn_number, lagrange_dof_number) += jac_term;
                  // And by symmetry, we can add the transpose contribution to
                  // the jacobian
                  jacobian(lagrange_dof_number, left_eqn_number) += jac_term;
                } // End orthogonality check
              } // End dof pinned check

              // Add the Hessian terms (gamma>=beta so we don't recount fifth
              // constraint)
              for (unsigned gamma = beta; gamma < 0; gamma++)
              {
                // Get the index of this Lagrange dof
                lagrange_type = 5 + beta + gamma;
                // Get the local equation number for this dof
                lagrange_dof_number =
                  internal_local_eqn(Index_of_lagrange_data, lagrange_type);
                // Check that this dof isn't pinned
                if (lagrange_dof_number)
                {
                  double jac_term = -hess_of_transform[alpha](beta, gamma);
                  // Orthogonality check
                  if (fabs(jac_term) > Orthogonality_tolerance)
                  {
                    // Add the contribution to the jacobian
                    jacobian(left_eqn_number, lagrange_dof_number) += jac_term;
                    // And by symmetry, we can add the transpose contribution to
                    // the jacobian
                    jacobian(lagrange_dof_number, left_eqn_number) += jac_term;
                  } // End of orthogonality check
                } // End dof pinned check
              } // End loop over second derivative index gamma
            } // End for loop over Lagrange dofs
          } // End equation pinned check
        } // End for loop over left node equations

        //----------------------------------------------------------------------
        // Lastly, fill in the jacobian entries for the three second derivative
        // equations in the left node. The two derivatives
        // -l_{5+\alpha+\beta}
        //     * J_{\alpha\gamma} * J_{\beta\delta} * L_{5+\gamma+\delta}
        for (unsigned alpha = 0; alpha < 2; alpha++)
        {
          for (unsigned beta = 0; beta < 2; beta++)
          {
            // Get the index of this equation in the left node
            unsigned left_type = 5 + alpha + beta;
            // Get the local equation number for this equation
            int left_eqn_number =
              external_local_eqn(Index_of_left_data, left_type);
            // Check this equation isn't pinned
            if (left_eqn_number >= 0)
            {
              // Loop over the lagrange dofs
              for (unsigned gamma = 0; gamma < 2; gamma++)
              {
                // The second derivative must have index >= gamma to prevent
                // double counting the mixed derivative condition.
                for (unsigned delta = gamma; delta < 2; delta++)
                {
                  // Get the index of this Lagrange dof
                  unsigned lagrange_type = 5 + gamma + delta;
                  // Get the local equation number for this dof
                  int lagrange_dof_number =
                    internal_local_eqn(Index_of_lagrange_data, lagrange_type);
                  // Check this dof isn't pinned
                  if (lagrange_dof_number >= 0)
                  {
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
                  } // End dof pinned check
                } // [delta]
              } // End lagrange loops [gamma]
            } // End equation pinned check
          } // [beta]
        } // End equation loops [alpha]
      } // End if(flag) check for jacobian
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
  };


} // namespace oomph


#endif

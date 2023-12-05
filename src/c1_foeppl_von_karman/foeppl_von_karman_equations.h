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
// Header file for the Foeppl-von Karman equation elements
#ifndef OOMPH_C1FVK_EQUATIONS_HEADER
#define OOMPH_C1FVK_EQUATIONS_HEADER


#include <sstream>

// OOMPH-LIB headers
#include "../generic/nodes.h"
#include "../generic/oomph_utilities.h"
#include "../generic/Telements.h"

// [zdec] TODO:
// -- Move function definitions to .cc

namespace oomph
{
  //=============================================================================
  /// A class for all subparametric elements that solve the
  /// Foeppl-von Karman equations (with artificial damping).
  /// \f[
  /// \nabla^4 w - \eta \frac{\partial}{\partial x_\alpha}
  /// \left(\sigma_{\alpha\beta}\frac{\partial w}{\partial x_\beta}\right)
  /// = p(\vec{x}) - \mu\frac{\partial w}{\partial t}
  /// \f]
  /// \f[
  /// \frac{\partial \sigma_{\alpha\beta}}{\partial x_\beta} = t_\alpha
  /// \f]
  ///
  /// This class contains the generic maths. Shape functions, geometric
  /// mapping etc. must get implemented in derived class.
  ///
  /// The equations allow for different basis and test functions
  /// for in-plane and out-of-plane unknowns, however it does
  /// assume that in-plane u_x and u_y have the same interpolation.
  ///
  /// In general, field interpolation may rely on nodal and internal data and
  /// these have been considered separately. To date, only the out-of-plane
  /// field has been written with the possibility of internal data in mind,
  /// however this can easily be added for in-plane dofs if required.
  //=============================================================================
  class FoepplVonKarmanEquations : public virtual FiniteElement
  {
  public:
    /// A pointer to a scalar function of the position. Can be used for
    /// out-of-plane forcing, swelling, isotropic-prestrain, etc.
    typedef void (*ScalarFctPt)(const Vector<double>& x,
				double& f);

    /// A pointer to a vector function of the position. Can be used for
    /// in-of-plane forcing, anisotropic-prestrain, etc.
    typedef void (*VectorFctPt)(const Vector<double>& x,
                                Vector<double>& f);

    /// Function pointer to the Error Metric we are using
    /// e.g could be that we are just interested in error on w etc.
    typedef void (*ErrorMetricFctPt)(const Vector<double>& x,
                                     const Vector<double>& u,
                                     const Vector<double>& u_exact,
                                     double& error,
                                     double& norm);

    /// Function pointer to the Error Metric we are using if we want multiple
    /// errors.  e.g could be we want errors seperately on each displacment
    typedef void (*MultipleErrorMetricFctPt)(const Vector<double>& x,
                                             const Vector<double>& u,
                                             const Vector<double>& u_exact,
                                             Vector<double>& error,
                                             Vector<double>& norm);


    //--------------------------------------------------------------------------
    // Pure virtual interfaces which must be implemented when geometry and bases
    // are added in the derived class

    /// (pure virtual) interface to return the number of nodes used by u
    virtual unsigned nu_node() const = 0;

    /// (pure virtual) interface to return the number of nodes used by w
    virtual unsigned nw_node() const = 0;


    /// (pure virtual) interface to get the local indices of the nodes used by u
    virtual Vector<unsigned> get_u_node_indices() const = 0;

    /// (pure virtual) interface to get the local indices of the nodes used by w
    virtual Vector<unsigned> get_w_node_indices() const = 0;


    /// (pure virtual) interface to get the number of basis types for u at node
    /// j
    virtual unsigned nu_type_at_each_node() const = 0;

    /// (pure virtual) interface to get the number of basis types for w at node
    /// j
    virtual unsigned nw_type_at_each_node() const = 0;


    /// (pure virtual) interface to retrieve the value of u_alpha at node j of
    /// type k
    virtual double get_u_alpha_value_at_node_of_type(
      const unsigned& alpha,
      const unsigned& j_node,
      const unsigned& k_type) const = 0;

    /// (pure virtual) interface to retrieve the t-th history value value of
    /// u_alpha at node j of type k
    virtual double get_u_alpha_value_at_node_of_type(
      const unsigned& t_time,
      const unsigned& alpha,
      const unsigned& j_node,
      const unsigned& k_type) const = 0;

    /// (pure virtual) interface to retrieve the value of w at node j of type k
    virtual double get_w_value_at_node_of_type(
      const unsigned& j_node, const unsigned& k_type) const = 0;

    /// (pure virtual) interface to retrieve the t-th history value of w at node
    /// j of type k
    virtual double get_w_value_at_node_of_type(
      const unsigned& t_time,
      const unsigned& j_node,
      const unsigned& k_type) const = 0;

    // [IN-PLANE-INTERNAL]
    // /// (pure virtual) interface to get the pointer to the internal data used
    // to
    // /// interpolate u (NOTE: assumes each u field has exactly one internal
    // data) virtual Vector<Data*> u_internal_data_pts() const = 0;

    /// (pure virtual) interface to get the pointer to the internal data used to
    /// interpolate w (NOTE: assumes w field has exactly one internal data)
    virtual Data* w_internal_data_pt() const = 0;


    /// (pure virtual) interface to get the number of internal types for the u
    /// fields
    virtual unsigned nu_type_internal() const = 0;

    /// (pure virtual) interface to get the number of internal types for the w
    /// fields
    virtual unsigned nw_type_internal() const = 0;


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

    /// (pure virtual) interface to retrieve the value of w of internal type k
    virtual double get_w_internal_value_of_type(
      const unsigned& k_type) const = 0;

    /// (pure virtual) interface to retrieve the t-th history value of w of
    /// internal type k
    virtual double get_w_internal_value_of_type(
      const unsigned& time, const unsigned& k_type) const = 0;


  protected:
    /// (pure virtual) In-plane basis functions and derivatives w.r.t. global
    /// coords at local coordinate s; return det(Jacobian of mapping)
    virtual double basis_u_foeppl_von_karman(
      const Vector<double>& s,
      Shape& psi_n) const = 0;

    /// (pure virtual) In-plane basis functions and derivatives w.r.t. global
    /// coords at local coordinate s; return det(Jacobian of mapping)
    virtual double dbasis_u_eulerian_foeppl_von_karman(
      const Vector<double>& s,
      Shape& psi_n,
      DShape& dpsi_n_dx) const = 0;

    /// (pure virtual) In-plane basis/test functions at and derivatives w.r.t
    /// global coords at local coordinate s; return det(Jacobian of mapping)
    virtual double dbasis_and_dtest_u_eulerian_foeppl_von_karman(
      const Vector<double>& s,
      Shape& psi_n,
      DShape& dpsi_n_dx,
      Shape& test_n,
      DShape& dtest_n_dx) const = 0;

    /// (pure virtual) Out-of-plane basis functions at local coordinate s
    virtual void basis_w_foeppl_von_karman(
      const Vector<double>& s,
      Shape& psi_n,
      Shape& psi_i) const = 0;

    /// (pure virtual) Out-of-plane basis functions and derivs w.r.t. global
    /// coords at local coordinate s; return det(Jacobian of mapping)
    virtual double d2basis_w_eulerian_foeppl_von_karman(
      const Vector<double>& s,
      Shape& psi_n,
      Shape& psi_i,
      DShape& dpsi_n_dx,
      DShape& dpsi_i_dx,
      DShape& d2psi_n_dx2,
      DShape& d2psi_i_dx2) const = 0;

    /// (pure virtual) Out-of-plane basis/test functions at local coordinate s
    virtual void basis_and_test_w_foeppl_von_karman(
      const Vector<double>& s,
      Shape& psi_n,
      Shape& psi_i,
      Shape& test_n,
      Shape& test_i) const = 0;

    /// (pure virtual) Out-of-plane basis/test functions and first derivs w.r.t.
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
      DShape& dtest_i_dx) const = 0;

    /// (pure virtual) Out-of-plane basis/test functions and first/second derivs
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
      DShape& d2test_i_dx2) const = 0;

    // End of pure virtual interface functions
    //----------------------------------------------------------------------------


    //----------------------------------------------------------------------------
    // Start of established FvK equation functions

  public:
    /// Constructor
    FoepplVonKarmanEquations()
      : Pressure_fct_pt(0),
        In_plane_forcing_fct_pt(0),
        Swelling_fct_pt(0),
	Solve_u_exact(false),
        Error_metric_fct_pt(0),
        Multiple_error_metric_fct_pt(0),
        U_is_damped(false),
        W_is_damped(true),
        Association_matrix_pt(0)
    {
      Eta_pt = &Default_Eta_Value;
      Nu_pt = &Default_Nu_Value;
      Mu_pt = &Default_Mu_Value;
    }

    /// Broken copy constructor
    FoepplVonKarmanEquations(const FoepplVonKarmanEquations& dummy)
    {
      BrokenCopy::broken_copy("FoepplVonKarmanEquations");
    }

    /// Broken assignment operator
    void operator=(const FoepplVonKarmanEquations&)
    {
      BrokenCopy::broken_assign("FoepplVonKarmanEquations");
    }


    /// Get pointer to association matrix
    DenseMatrix<double>* get_association_matrix_pt() const
    {
      return Association_matrix_pt;
    };

    /// Eta [zdec] why the ampersand?
    const double& get_eta() const
    {
      return *Eta_pt;
    }

    /// Pointer to eta
    const double*& eta_pt()
    {
      return Eta_pt;
    }

    /// Output with default number of plot points
    void output(std::ostream& outfile)
    {
      const unsigned n_plot = 5;
      FoepplVonKarmanEquations::output(outfile, n_plot);
    }

    /// Output with default number of plot points
    void full_output(std::ostream& outfile)
    {
      const unsigned n_plot = 5;
      FoepplVonKarmanEquations::full_output(outfile, n_plot);
    }

    /// Output FE representation of soln: x,y,u or x,y,z,u at
    /// n_plot*(n_plot+1)/2 plot points
    void output(std::ostream& outfile, const unsigned& n_plot);

    /// Full output function with a rich set of unknowns:
    ///  x, y, ux, uy, w, du, dw, d2w, strain, stress, principal stress
    /// at n_plot*(n_plot+1)/2 plot points
    void full_output(std::ostream& outfile, const unsigned& n_plot);

    /// C_style output with default number of plot points
    void output(FILE* file_pt)
    {
      const unsigned n_plot = 5;
      FoepplVonKarmanEquations::output(file_pt, n_plot);
    }

    /// C-style output FE representation of soln: x,y,u or x,y,z,u at
    /// n_plot*(n_plot+1)/2 plot points
    void output(FILE* file_pt, const unsigned& n_plot);

    /// Output exact soln: x,y,u_exact at n_plot*(n_plot_1)/2 plot points
    void output_fct(std::ostream& outfile,
                    const unsigned& n_plot,
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt);

    /// Output exact soln: x,y,u_exact or x,y,z,u_exact at n_plot^DIM plot
    /// points (dummy time-dependent version to keep intel compiler happy)
    virtual void output_fct(
      std::ostream& outfile,
      const unsigned& n_plot,
      const double& time,
      FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
    {
      throw OomphLibError(
        "There is no time-dependent output_fct() for these elements ",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }


    /// (pure virtual) interface to return a vector of the indices
    virtual Vector<unsigned> u_field_indices() const = 0;

    /// (pure virtual) interface to return the number of nodes used by w
    virtual unsigned w_field_index() const = 0;


    /// Get error against and norm of exact solution
    void compute_error(std::ostream& outfile,
                       FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
                       double& error,
                       double& norm);

    /// Get error against and norm of exact solution
    void compute_error_in_deflection(
      std::ostream& outfile,
      FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
      double& error,
      double& norm);

    /// Dummy, time dependent error checker
    void compute_error(std::ostream& outfile,
                       FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,
                       const double& time,
                       double& error,
                       double& norm)
    {
      throw OomphLibError(
        "There is no time-dependent compute_error() for these elements",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }

    /// Fill in the strain tensor from displacement gradients
    void get_epsilon(DenseMatrix<double>& epsilon,
                     const DenseMatrix<double>& grad_u,
                     const DenseMatrix<double>& grad_w,
                     const double& c_swell) const
    {
      // Truncated Green Lagrange strain tensor
      DenseMatrix<double> dummy_epsilon(this->dim(), this->dim(), 0.0);
      for (unsigned alpha = 0; alpha < this->dim(); ++alpha)
      {
        for (unsigned beta = 0; beta < this->dim(); ++beta)
        {
          // Truncated Green Lagrange strain tensor
          dummy_epsilon(alpha, beta) +=
            0.5 * grad_u(alpha, beta) + 0.5 * grad_u(beta, alpha) +
            0.5 * grad_w(0, alpha) * grad_w(0, beta);
        }
        // Swelling slack
        dummy_epsilon(alpha, alpha) -= c_swell;
      }
      epsilon = dummy_epsilon;
    }

    /// Fill in the stress tensor from displacement gradients
    void get_sigma(DenseMatrix<double>& sigma,
                   const DenseMatrix<double>& grad_u,
                   const DenseMatrix<double>& grad_w,
                   const double c_swell) const
    {
      // Get the Poisson ratio
      double nu(get_nu());

      // Truncated Green Lagrange strain tensor
      DenseMatrix<double> epsilon(this->dim(), this->dim(), 0.0);
      get_epsilon(epsilon, grad_u, grad_w, c_swell);

      // Empty sigma
      sigma(0, 0) = 0.0;
      sigma(0, 1) = 0.0;
      sigma(1, 0) = 0.0;
      sigma(1, 1) = 0.0;

      // Now construct the Stress
      for (unsigned alpha = 0; alpha < this->dim(); ++alpha)
      {
        for (unsigned beta = 0; beta < this->dim(); ++beta)
        {
          // The Laplacian term: Trace[ \epsilon ] I
          // \nu * \epsilon_{\alpha \beta} delta_{\gamma \gamma}
          sigma(alpha, alpha) += nu * epsilon(beta, beta) / (1 - nu * nu);

          // The scalar transform term: \epsilon
          // (1-\nu) * \epsilon_{\alpha \beta}
          sigma(alpha, beta) += (1 - nu) * epsilon(alpha, beta) / (1 - nu * nu);
        }
      }
    }

    /// Fill in the stress tensor using a precalculated strain tensor
    void get_sigma_from_epsilon(DenseMatrix<double>& sigma,
                                const DenseMatrix<double>& epsilon) const
    {
      // Get the Poisson ratio
      double nu(get_nu());

      // Now construct the Stress
      sigma(0, 0) = (epsilon(0, 0) + nu * epsilon(1, 1)) / (1.0 - nu * nu);
      sigma(1, 1) = (epsilon(1, 1) + nu * epsilon(0, 0)) / (1.0 - nu * nu);
      sigma(0, 1) = epsilon(0, 1) / (1.0 + nu);
      sigma(1, 0) = sigma(0, 1);
    }

    /// Get the principal stresses from the stress tensor
    void get_principal_stresses(const DenseMatrix<double>& sigma,
                                Vector<double>& eigenvals,
                                DenseMatrix<double>& eigenvecs) const
    {
      // Ensure that our eigenvectors are the right size
      eigenvals.resize(2);
      eigenvecs.resize(2);

      // Store the axial and shear stresses
      double s00 = sigma(0, 0);
      double s01 = sigma(0, 1);
      double s11 = sigma(1, 1);

      // Calculate the principal stress magnitudes
      eigenvals[0] = 0.5 * ((s00 + s11) + sqrt((s00 + s11) * (s00 + s11) -
                                               4.0 * (s00 * s11 - s01 * s01)));
      eigenvals[1] = 0.5 * ((s00 + s11) - sqrt((s00 + s11) * (s00 + s11) -
                                               4.0 * (s00 * s11 - s01 * s01)));

      // Handle the shear free case
      if (s01 == 0.0)
      {
        eigenvecs(0, 0) = 1.0;
        eigenvecs(1, 0) = 0.0;
        eigenvecs(0, 1) = 0.0;
        eigenvecs(1, 1) = 1.0;
      }

      else
      {
        // TODO: better (more general) sign choice for streamlines

        // For max eval we choose y-positive evecs (suited to swelling sheet
        // problem)
        double sign = (eigenvals[0] - s00 < 0.0) ? -1.0 : 1.0;
        // Calculate the normalised principal stress direction for eigenvals[0]
        eigenvecs(0, 0) =
          sign *
          (s01 / sqrt(s01 * s01 + (eigenvals[0] - s00) * (eigenvals[0] - s00)));
        eigenvecs(1, 0) =
          sign *
          ((eigenvals[0] - s00) /
           sqrt(s01 * s01 + (eigenvals[0] - s00) * (eigenvals[0] - s00)));

        // For min eval we choose x-positive evecs (suited to swelling sheet
        // problem)
        sign = (s01 < 0.0) ? -1.0 : 1.0;
        // Calculate the normalised principal stress direction for eigenvals[1]
        eigenvecs(0, 1) =
          sign *
          (s01 / sqrt(s01 * s01 + (eigenvals[1] - s00) * (eigenvals[1] - s00)));
        eigenvecs(1, 1) =
          sign *
          ((eigenvals[1] - s00) /
           sqrt(s01 * s01 + (eigenvals[1] - s00) * (eigenvals[1] - s00)));
      }
    }

    /// Access function: Pointer to pressure function
    ScalarFctPt& pressure_fct_pt()
    {
      return Pressure_fct_pt;
    }

    /// Access function: Pointer to pressure function. Const version
    ScalarFctPt pressure_fct_pt() const
    {
      return Pressure_fct_pt;
    }

    /// Access function: Pointer to in plane forcing function
    VectorFctPt& in_plane_forcing_fct_pt()
    {
      return In_plane_forcing_fct_pt;
    }

    /// Access function: Pointer to in plane forcing function. Const version
    VectorFctPt in_plane_forcing_fct_pt() const
    {
      return In_plane_forcing_fct_pt;
    }

    /// Access function: Pointer to swelling function
    ScalarFctPt& swelling_fct_pt()
    {
      return Swelling_fct_pt;
    }

    /// Access function: Pointer to swelling function. Const version
    ScalarFctPt swelling_fct_pt() const
    {
      return Swelling_fct_pt;
    }

    /// Access function: Pointer to error metric function
    ErrorMetricFctPt& error_metric_fct_pt()
    {
      return Error_metric_fct_pt;
    }

    /// Access function: Pointer to multiple error metric function
    MultipleErrorMetricFctPt& multiple_error_metric_fct_pt()
    {
      return Multiple_error_metric_fct_pt;
    }

    /// Access function: Pointer to error metric function function
    ErrorMetricFctPt error_metric_fct_pt() const
    {
      return Error_metric_fct_pt;
    }

    /// Access function: Pointer to multiple error metric function
    MultipleErrorMetricFctPt multiple_error_metric_fct_pt() const
    {
      return Multiple_error_metric_fct_pt;
    }

    /// Access function: Bool that swaps element residuals to solve
    ///  (ux,uy,w) = (ux_ex,uy_ex,w_ex)
    /// which is useful for applying initial guesses or for debugging
    bool get_solve_u_exact() const
    {
      return Solve_u_exact;
    }

    /// Access function: Get pointer to u_exact function
    void get_u_exact_pt(VectorFctPt& u_exact_pt)
    {
      u_exact_pt = U_exact_pt;
    }

    /// Function to get the u_exact displacement field at global coordinates x
    void get_u_exact(const Vector<double>& x, Vector<double>& u_exact)
    {
      // Get the exact function pointer
      VectorFctPt u_exact_function_pt = 0;
      get_u_exact_pt(u_exact_function_pt);

      // Call the exact function at x and assign it to u
      (*u_exact_function_pt)(x,u_exact);
    }

    /// Activate the alternative residuals which are used to solve for an exact
    /// displacement field assignment (doesn't set exact function pointer)
    void activate_u_exact_solve()
    {
      // Throw an error if we haven't set a function pointer for u_exact yet
      if(!U_exact_pt)
      {
	throw OomphLibError("You need to set U_exact_pt before activating\
 Solve_u_exact",
			    OOMPH_EXCEPTION_LOCATION,
			    OOMPH_CURRENT_FUNCTION);
      }
      // Set the alternative equations flag to be true
      Solve_u_exact = true;
    }

    /// Activate the alternative residuals and set the exact function pointer
    void activate_u_exact_solve(VectorFctPt u_exact_pt)
    {
      // Update the exact displacement function pointer
      U_exact_pt = u_exact_pt;
      // Call the basic function which does the rest of the work
      activate_u_exact_solve();
    }

    /// Deactivate the alternative residuals
    void deactivate_u_exact_solve()
    {
      // Set the alternative equations flag to be false
      Solve_u_exact = false;
    }

    /// Access value of the damping flag for u
    virtual bool u_is_damped() const
    {
      return U_is_damped;
    }

    /// Access value of the damping flag for w
    virtual bool w_is_damped() const
    {
      return W_is_damped;
    }

    /// Access function to the Poisson ratio.
    const double*& nu_pt()
    {
      return Nu_pt;
    }

    /// Access function to the Poisson ratio (const version)
    const double& get_nu() const
    {
      return *Nu_pt;
    }

    /// Access function to the dampening coefficient.
    const double*& mu_pt()
    {
      return Mu_pt;
    }

    /// Access function to the dampening coefficient (const version)
    const double& get_mu() const
    {
      return *Mu_pt;
    }

    /// Get pressure term at (Eulerian) position x. This function is
    /// virtual to allow overloading in multi-physics problems where
    /// the strength of the pressure function might be determined by
    /// another system of equations.
    inline virtual void get_pressure_foeppl_von_karman(const unsigned& ipt,
                                                       const Vector<double>& x,
                                                       double& pressure) const
    {
      // If no pressure function has been set, return zero
      if (Pressure_fct_pt == 0)
      {
        pressure = 0.0;
      }
      else
      {
        // Get pressure strength
        (*Pressure_fct_pt)(x, pressure);
      }
    }

    /// Get pressure term at (Eulerian) position x. This function is
    /// virtual to allow overloading in multi-physics problems where
    /// the strength of the pressure function might be determined by
    /// another system of equations.
    inline virtual void get_in_plane_forcing_foeppl_von_karman(
      const unsigned& ipt,
      const Vector<double>& x,
      Vector<double>& traction) const
    {
      // In plane is same as DIM of problem (2)
      traction.resize(this->dim());
      // If no pressure function has been set, return zero
      if (In_plane_forcing_fct_pt == 0)
      {
        traction[0] = 0.0;
        traction[1] = 0.0;
      }
      else
      {
        // Get pressure strength
        (*In_plane_forcing_fct_pt)(x, traction);
      }
    }


    /// Get swelling at (Eulerian) position x. This function is
    /// virtual to allow overloading. [zdec] add ipt
    inline virtual void get_swelling_foeppl_von_karman(const Vector<double>& x,
                                                       double& swelling) const
    {
      // If no swelling function has been set, return zero
      if (Swelling_fct_pt == 0)
      {
        swelling = 0.0;
      }
      else
      {
        // Get swelling magnitude
        (*Swelling_fct_pt)(x, swelling);
      }
    }


    /// Add the element's contribution to its residual vector (wrapper)
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Call the generic residuals function with flag set to 0
      // using a dummy matrix argument
      fill_in_generic_residual_contribution_foeppl_von_karman(
        residuals, GeneralisedElement::Dummy_matrix, 0);
    }


    /// Add the element's contribution to its residual vector and
    /// element Jacobian matrix (wrapper)
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      // Call the generic routine with the flag set to 1
      fill_in_generic_residual_contribution_foeppl_von_karman(
        residuals, jacobian, 1);
    }


    /// Add the element's contribution to its residual vector and
    /// element Jacobian matrix (wrapper)
    void fill_in_contribution_to_jacobian_and_mass_matrix(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& mass_matrix)
    {
      // Call fill in Jacobian
      fill_in_contribution_to_jacobian(residuals, jacobian);
      // There is no mass matrix: we will just want J w = 0

      // -- COPIED FROM DISPLACMENT FVK EQUATIONS --
      // Dummy diagonal (won't result in global unit matrix but
      // doesn't matter for zero eigenvalue/eigenvector
      unsigned ndof = mass_matrix.nrow();
      for (unsigned i = 0; i < ndof; i++)
      {
        mass_matrix(i, i) += 1.0;
      }
    }


    /// Return FE representation of the three displacements at local coordinate
    /// s
    ///   0: u_x
    ///   1: u_y
    ///   2: w
    virtual inline Vector<double> interpolated_fvk_disp(
      const Vector<double>& s) const
    {
      // The indices of in-plane and out-of-plane unknowns
      const unsigned n_u_fields = 2;
      const unsigned n_w_fields = 1;

      // Find the dimension of the element
      const unsigned dim = this->dim();

      // Find out how many nodes there are for each field
      const unsigned n_u_node = nu_node();
      const unsigned n_w_node = nw_node();

      // Get the vector of nodes used for each field
      const Vector<unsigned> u_nodes = get_u_node_indices();
      const Vector<unsigned> w_nodes = get_w_node_indices();

      // Find out how many basis types there are at each node
      const unsigned n_u_nodal_type = nu_type_at_each_node();
      const unsigned n_w_nodal_type = nw_type_at_each_node();

      // Find out how many basis types there are internally
      // NOTE: In general, the in-plane basis may require internal basis/test
      // functions however, no such basis has been used so far and so this is
      // not implemented. Comments marked with "[IN-PLANE-INTERNAL]" indicate
      // locations where changes must be made in the case that internal data is
      // used for the in-plane unknowns.
      const unsigned n_u_internal_type = nu_type_internal();
      const unsigned n_w_internal_type = nw_type_internal();

#ifdef PARANOID
      // [IN-PLANE-INTERNAL]
      // This PARANOID block should be deleted if/when internal in-plane
      // contributions have been implemented

      // Throw an error if the number of internal in-plane basis types is
      // non-zero
      if (n_u_internal_type != 0)
      {
        throw OomphLibError(
          "The number of internal basis types for u is non-zero\
 but this functionality is not yet implemented. If you want to implement this\
 look for comments containing the tag [IN-PLANE-INTERNAL].",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // In-plane local basis & test functions
      // ------------------------------------------
      // Local in-plane nodal basis and test functions
      Shape psi_n_u(n_u_node, n_u_nodal_type);
      // [IN-PLANE-INTERNAL]
      // // Local in-plane internal basis and test functions
      // Shape psi_i_u(n_u_internal_type);

      // Out-of-plane local basis & test functions
      // ------------------------------------------
      // Nodal basis & test functions
      Shape psi_n_w(n_w_node, n_w_nodal_type);
      Shape test_n_w(n_w_node, n_w_nodal_type);
      // Internal basis & test functions
      Shape psi_i_w(n_w_internal_type);
      Shape test_i_w(n_w_internal_type);

      // Allocate and find global x
      Vector<double> interp_x(dim, 0.0);
      interpolated_x(s, interp_x);

      // Call the derivatives of the shape and test functions for the out of
      // plane unknown
      basis_w_foeppl_von_karman(s, psi_n_w, psi_i_w);

      // [IN-PLANE-INTERNAL]
      // This does not include internal basis functions for u. If they are
      // needed they must be added (as above for w)
      basis_u_foeppl_von_karman(s, psi_n_u);


      //======================= INTERPOLATION =================================

      // Create space for in-plane interpolated unknowns
      Vector<double> interpolated_u(n_u_fields, 0.0);
      // Create space for out-of-plane interpolated unknowns
      Vector<double> interpolated_w(n_w_fields, 0.0);

      //---Nodal contribution to the in-plane unknowns--------------------------
      // Loop over nodes used by in-plane fields
      for (unsigned j_node = 0; j_node < n_u_node; j_node++)
      {
        // Get the j-th node used by in-plane fields
        unsigned j_node_local = u_nodes[j_node];

        // Loop over types
        for (unsigned k_type = 0; k_type < n_u_nodal_type; k_type++)
        {
          // Loop over in-plane unknowns
          for (unsigned alpha = 0; alpha < n_u_fields; alpha++)
          {
            // Get the nodal value of type k
            double u_value =
              get_u_alpha_value_at_node_of_type(alpha, j_node_local, k_type);

            // --- Displacement ---
            // Add nodal contribution of type k to the interpolated displacement
            interpolated_u[alpha] += u_value * psi_n_u(j_node, k_type);

          } // End of loop over the index of u -- alpha
        } // End of loop over the types -- k_type
      } // End of loop over the nodes -- j_node


      //---Internal contribution to the in-plane
      // unknowns-------------------------
      // [IN-PLANE-INTERNAL]
      // Internal contributions to in-plane interpolation not written


      //---Nodal contribution to the out-of-plane
      // unknowns------------------------
      for (unsigned j_node = 0; j_node < n_w_node; j_node++)
      {
        // Get the j-th node used by in-plane fields
        unsigned j_node_local = w_nodes[j_node];

        // Loop over types
        for (unsigned k_type = 0; k_type < n_w_nodal_type; k_type++)
        {
          // Get the nodal value of type k
          double w_value = get_w_value_at_node_of_type(j_node_local, k_type);

          // --- Displacement ---
          // Add nodal contribution of type k to the interpolated displacement
          interpolated_w[0] += w_value * psi_n_w(j_node, k_type);

        } // End of loop over the types -- k_type
      } // End of loop over the nodes -- j_node

      //---Internal contribution to the out-of-plane field----------------------
      // Loop over the internal data types
      for (unsigned k_type = 0; k_type < n_w_internal_type; k_type++)
      {
        // Get the nodal value of type k
        double w_value = get_w_internal_value_of_type(k_type);

        // --- Displacement ---
        // Add nodal contribution of type k to the interpolated displacement
        interpolated_w[0] += w_value * psi_i_w(k_type);

      } // End of loop over internal types -- k_type

      //====================== END OF INTERPOLATION ============================

      // Copy our interpolated fields into the order we want and return them.
      Vector<double> interpolated_vals(12, 0.0);
      interpolated_vals[0] = interpolated_u[0]; // ux
      interpolated_vals[1] = interpolated_u[1]; // uy
      interpolated_vals[2] = interpolated_w[0]; // w
      return (interpolated_vals);
    } // end of interpolated_fvk_disp


    /// Return FE representation of the three displacements and their
    /// derivatives at local coordinate s:
    ///   0: u_x          6: du_y/dy
    ///   1: u_y          7: dw/dx
    ///   2: w            8: dw/dy
    ///   3: du_x/dx      9: d2w/dx2
    ///   4: du_x/dy     10: d2w/dxdy
    ///   5: du_y/dx     11: d2w/dy2
    virtual inline Vector<double> interpolated_fvk_disp_and_deriv(
      const Vector<double>& s) const
    {
      // The indices of in-plane and out-of-plane unknowns
      const unsigned n_u_fields = 2;
      const unsigned n_w_fields = 1;

      // Find the dimension of the element
      const unsigned dim = this->dim();
      // The number of first derivatives is the dimension of the element
      const unsigned n_deriv = dim;
      // The number of second derivatives is the triangle number of the
      // dimension
      const unsigned n_2deriv = dim * (dim + 1) / 2;

      // Find out how many nodes there are for each field
      const unsigned n_u_node = nu_node();
      const unsigned n_w_node = nw_node();

      // Get the vector of nodes used for each field
      const Vector<unsigned> u_nodes = get_u_node_indices();
      const Vector<unsigned> w_nodes = get_w_node_indices();

      // Find out how many basis types there are at each node
      const unsigned n_u_nodal_type = nu_type_at_each_node();
      const unsigned n_w_nodal_type = nw_type_at_each_node();

      // Find out how many basis types there are internally
      // NOTE: In general, the in-plane basis may require internal basis/test
      // functions however, no such basis has been used so far and so this is
      // not implemented. Comments marked with "[IN-PLANE-INTERNAL]" indicate
      // locations where changes must be made in the case that internal data is
      // used for the in-plane unknowns.
      const unsigned n_u_internal_type = nu_type_internal();
      const unsigned n_w_internal_type = nw_type_internal();

#ifdef PARANOID
      // [IN-PLANE-INTERNAL]
      // This PARANOID block should be deleted if/when internal in-plane
      // contributions have been implemented

      // Throw an error if the number of internal in-plane basis types is
      // non-zero
      if (n_u_internal_type != 0)
      {
        throw OomphLibError(
          "The number of internal basis types for u is non-zero\
 but this functionality is not yet implemented. If you want to implement this\
 look for comments containing the tag [IN-PLANE-INTERNAL].",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // In-plane local basis & test functions
      // ------------------------------------------
      // Local in-plane nodal basis and test functions
      Shape psi_n_u(n_u_node, n_u_nodal_type);
      Shape test_n_u(n_u_node, n_u_nodal_type);
      DShape dpsi_n_udxi(n_u_node, n_u_nodal_type, n_deriv);
      DShape dtest_n_udxi(n_u_node, n_u_nodal_type, n_deriv);
      // [IN-PLANE-INTERNAL]
      // Uncomment this block of Shape/DShape functions to use for the internal
      // in-plane basis
      // // Local in-plane internal basis and test functions
      // Shape psi_i_u(n_u_internal_type);
      // Shape test_i_u(n_u_internal_type);
      // DShape dpsi_i_udxi(n_u_internal_type, n_deriv);
      // DShape dtest_i_udxi(n_u_internal_type, n_deriv);

      // Out-of-plane local basis & test functions
      // ------------------------------------------
      // Nodal basis & test functions
      Shape psi_n_w(n_w_node, n_w_nodal_type);
      Shape test_n_w(n_w_node, n_w_nodal_type);
      DShape dpsi_n_wdxi(n_w_node, n_w_nodal_type, n_deriv);
      DShape dtest_n_wdxi(n_w_node, n_w_nodal_type, n_deriv);
      DShape d2psi_n_wdxi2(n_w_node, n_w_nodal_type, n_2deriv);
      DShape d2test_n_wdxi2(n_w_node, n_w_nodal_type, n_2deriv);
      // Internal basis & test functions
      Shape psi_i_w(n_w_internal_type);
      Shape test_i_w(n_w_internal_type);
      DShape dpsi_i_wdxi(n_w_internal_type, n_deriv);
      DShape dtest_i_wdxi(n_w_internal_type, n_deriv);
      DShape d2psi_i_wdxi2(n_w_internal_type, n_2deriv);
      DShape d2test_i_wdxi2(n_w_internal_type, n_2deriv);

      // Allocate and find global x
      Vector<double> interp_x(dim, 0.0);
      interpolated_x(s, interp_x);

      // Call the derivatives of the shape and test functions for the out of
      // plane unknown
      d2basis_w_eulerian_foeppl_von_karman(s,
                                           psi_n_w,
                                           psi_i_w,
                                           dpsi_n_wdxi,
                                           dpsi_i_wdxi,
                                           d2psi_n_wdxi2,
                                           d2psi_i_wdxi2);

      // [IN-PLANE-INTERNAL]
      // This does not include internal basis functions for u. If they are
      // needed they must be added (as above for w)
      dbasis_u_eulerian_foeppl_von_karman(s, psi_n_u, dpsi_n_udxi);


      //======================= INTERPOLATION =================================

      // Create space for in-plane interpolated unknowns
      Vector<double> interpolated_u(n_u_fields, 0.0);
      Vector<double> interpolated_dudt(n_u_fields, 0.0);
      DenseMatrix<double> interpolated_dudxi(n_u_fields, n_deriv, 0.0);
      // Create space for out-of-plane interpolated unknowns
      Vector<double> interpolated_w(n_w_fields, 0.0);
      Vector<double> interpolated_dwdt(n_w_fields, 0.0);
      DenseMatrix<double> interpolated_dwdxi(n_w_fields, n_deriv, 0.0);
      DenseMatrix<double> interpolated_d2wdxi2(n_w_fields, n_2deriv, 0.0);

      //---Nodal contribution to the in-plane unknowns--------------------------
      // Loop over nodes used by in-plane fields
      for (unsigned j_node = 0; j_node < n_u_node; j_node++)
      {
        // Get the j-th node used by in-plane fields
        unsigned j_node_local = u_nodes[j_node];

        // By default, we have no history values (no damping)
        unsigned n_time = 0;
        // Turn on damping if this field requires it AND we are not doing a
        // steady solve
        TimeStepper* timestepper_pt =
          this->node_pt(j_node_local)->time_stepper_pt();
        bool damping = u_is_damped() && !(timestepper_pt->is_steady());
        if (damping)
        {
          n_time = timestepper_pt->ntstorage();
        }

        // Loop over types
        for (unsigned k_type = 0; k_type < n_u_nodal_type; k_type++)
        {
          // Loop over in-plane unknowns
          for (unsigned alpha = 0; alpha < n_u_fields; alpha++)
          {
            // --- Time derivative ---
            double nodal_dudt_value = 0.0;
            // Loop over the history values (if damping, then n_time>0) and
            // add history contribution to nodal contribution to time derivative
            for (unsigned t_time = 0; t_time < n_time; t_time++)
            {
              nodal_dudt_value += get_u_alpha_value_at_node_of_type(
                                    t_time, alpha, j_node_local, k_type) *
                                  timestepper_pt->weight(1, t_time);
            }
            interpolated_dudt[alpha] +=
              nodal_dudt_value * psi_n_u(j_node, k_type);

            // Get the nodal value of type k
            double u_value =
              get_u_alpha_value_at_node_of_type(alpha, j_node_local, k_type);

            // --- Displacement ---
            // Add nodal contribution of type k to the interpolated displacement
            interpolated_u[alpha] += u_value * psi_n_u(j_node, k_type);

            // --- First derivatives ---
            for (unsigned l_deriv = 0; l_deriv < n_deriv; l_deriv++)
            {
              // Add the nodal contribution of type k to the derivative of the
              // displacement
              interpolated_dudxi(alpha, l_deriv) +=
                u_value * dpsi_n_udxi(j_node, k_type, l_deriv);
            }

            // --- No second derivatives for in-plane ---

          } // End of loop over the index of u -- alpha
        } // End of loop over the types -- k_type
      } // End of loop over the nodes -- j_node


      //---Internal contribution to the in-plane
      // unknowns-------------------------
      // [IN-PLANE-INTERNAL]
      // Internal contributions to in-plane interpolation not written


      //---Nodal contribution to the out-of-plane
      // unknowns------------------------
      for (unsigned j_node = 0; j_node < n_w_node; j_node++)
      {
        // Get the j-th node used by in-plane fields
        unsigned j_node_local = w_nodes[j_node];

        // By default, we have no history values (no damping)
        unsigned n_time = 0;
        // Turn on damping if this field requires it AND we are not doing a
        // steady solve
        TimeStepper* timestepper_pt =
          this->node_pt(j_node_local)->time_stepper_pt();
        bool damping = w_is_damped() && !(timestepper_pt->is_steady());
        if (damping)
        {
          n_time = timestepper_pt->ntstorage();
        }

        // Loop over types
        for (unsigned k_type = 0; k_type < n_w_nodal_type; k_type++)
        {
          // --- Time derivative ---
          double nodal_dwdt_value = 0.0;
          // Loop over the history values (if damping, then n_time>0) and
          // add history contribution to nodal contribution to time derivative
          for (unsigned t_time = 0; t_time < n_time; t_time++)
          {
            nodal_dwdt_value +=
              get_w_value_at_node_of_type(t_time, j_node_local, k_type) *
              timestepper_pt->weight(1, t_time);
          }
          interpolated_dwdt[0] += nodal_dwdt_value * psi_n_w(j_node, k_type);

          // Get the nodal value of type k
          double w_value = get_w_value_at_node_of_type(j_node_local, k_type);

          // --- Displacement ---
          // Add nodal contribution of type k to the interpolated displacement
          interpolated_w[0] += w_value * psi_n_w(j_node, k_type);

          // --- First derivatives ---
          for (unsigned l_deriv = 0; l_deriv < n_deriv; l_deriv++)
          {
            // Add the nodal contribution of type k to the derivative of the
            // displacement
            interpolated_dwdxi(0, l_deriv) +=
              w_value * dpsi_n_wdxi(j_node, k_type, l_deriv);
          }

          // --- Second derivatives ---
          for (unsigned l_2deriv = 0; l_2deriv < n_2deriv; l_2deriv++)
          {
            // Add the nodal contribution of type k to the derivative of the
            // displacement
            interpolated_d2wdxi2(0, l_2deriv) +=
              w_value * d2psi_n_wdxi2(j_node, k_type, l_2deriv);
          }

        } // End of loop over the types -- k_type
      } // End of loop over the nodes -- j_node

      //---Internal contribution to the out-of-plane field----------------------
      // --- Set up damping ---
      // By default, we have no history values (no damping)
      unsigned n_time = 0;
      // Turn on damping if this field requires it AND we are not doing a steady
      // solve
      TimeStepper* timestepper_pt =
        this->w_internal_data_pt()->time_stepper_pt();
      bool damping = w_is_damped() && !(timestepper_pt->is_steady());
      if (damping)
      {
        n_time = timestepper_pt->ntstorage();
      }
      // Loop over the internal data types
      for (unsigned k_type = 0; k_type < n_w_internal_type; k_type++)
      {
        // --- Time derivative ---
        double nodal_dwdt_value = 0.0;
        // Loop over the history values (if damping, then n_time>0) and
        // add history contribution to nodal contribution to time derivative
        for (unsigned t_time = 0; t_time < n_time; t_time++)
        {
          nodal_dwdt_value += get_w_internal_value_of_type(t_time, k_type) *
                              timestepper_pt->weight(1, t_time);
        }
        interpolated_dwdt[0] += nodal_dwdt_value * psi_i_w(k_type);

        // Get the nodal value of type k
        double w_value = get_w_internal_value_of_type(k_type);

        // --- Displacement ---
        // Add nodal contribution of type k to the interpolated displacement
        interpolated_w[0] += w_value * psi_i_w(k_type);

        // --- First derivatives ---
        for (unsigned l_deriv = 0; l_deriv < n_deriv; l_deriv++)
        {
          // Add the nodal contribution of type k to the derivative of the
          // displacement
          interpolated_dwdxi(0, l_deriv) +=
            w_value * dpsi_i_wdxi(k_type, l_deriv);
        }

        // --- Second derivatives ---
        for (unsigned l_2deriv = 0; l_2deriv < n_2deriv; l_2deriv++)
        {
          // Add the nodal contribution of type k to the derivative of the
          // displacement
          interpolated_d2wdxi2(0, l_2deriv) +=
            w_value * d2psi_i_wdxi2(k_type, l_2deriv);
        }

      } // End of loop over internal types -- k_type

      //====================== END OF INTERPOLATION ============================

      // Copy our interpolated fields into f in the order we want and return it.
      Vector<double> interpolated_vals(12, 0.0);
      interpolated_vals[0] = interpolated_u[0]; // ux
      interpolated_vals[1] = interpolated_u[1]; // uy
      interpolated_vals[2] = interpolated_w[0]; // w
      interpolated_vals[3] = interpolated_dudxi(0, 0); // duxdx
      interpolated_vals[4] = interpolated_dudxi(0, 1); // duxdy
      interpolated_vals[5] = interpolated_dudxi(1, 0); // duydx
      interpolated_vals[6] = interpolated_dudxi(1, 1); // duydy
      interpolated_vals[7] = interpolated_dwdxi(0, 0); // dwdx
      interpolated_vals[8] = interpolated_dwdxi(0, 1); // dwdy
      interpolated_vals[9] = interpolated_d2wdxi2(0, 0); // d2wdx2
      interpolated_vals[10] = interpolated_d2wdxi2(0, 1); // d2wdxdy
      interpolated_vals[11] = interpolated_d2wdxi2(0, 2); // d2wdy2
      return (interpolated_vals);
    }
    // END OF INTERPOLATED U FVK [zdec]

    /// Self-test: Return 0 for OK
    unsigned self_test();


  protected:
    /// Compute element residual Vector only (if flag=and/or element
    /// Jacobian matrix
    virtual void fill_in_generic_residual_contribution_foeppl_von_karman(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      const unsigned& flag);

    /// Add the element's contribution to its residual vector for the
    /// alternative problem of (ux,uy,w) = (ux_exact,uy_exact,w_exact)
    void fill_in_contribution_to_residual_u_exact(Vector<double> &residuals);

    /// Pointer to pressure function:
    ScalarFctPt Pressure_fct_pt;

    /// Pointer to in plane forcing function (i.e. the shear force applied to
    /// the face)
    VectorFctPt In_plane_forcing_fct_pt;

    /// Pointer to swelling function:
    ScalarFctPt Swelling_fct_pt;

    /// Pointer to Poisson ratio, which this element cannot modify
    const double* Nu_pt;

    /// Pointer to the dampening coefficient, which this element cannot modify
    const double* Mu_pt;

    /// Pointer to global eta
    const double* Eta_pt;

    /// Default value for physical constant: Poisson ratio.
    static const double Default_Nu_Value;

    /// Default value for constant: dampening coefficient.
    static const double Default_Mu_Value;

    /// Default eta value so that we use 'natural' nondim and have no h
    /// dependence.
    static const double Default_Eta_Value;

    /// Flag to swap to alternative solve
    bool Solve_u_exact;

    /// Pointer to exact displacement field function pointer
    VectorFctPt U_exact_pt;

    /// Pointer to error metric
    ErrorMetricFctPt Error_metric_fct_pt;

    /// Pointer to error metric when we want multiple errors
    MultipleErrorMetricFctPt Multiple_error_metric_fct_pt;

    /// Flag to control damping of the in-plane variables
    bool U_is_damped;

    /// Flag to control damping of the out-of-plane variable
    bool W_is_damped;

    /// Pointer to precomputed matrix that associates shape functions to
    /// monomials
    DenseMatrix<double>* Association_matrix_pt;
  };


} // end namespace oomph
#endif

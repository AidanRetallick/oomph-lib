// Header file for the Koiter Steigmann Equations class
#ifndef OOMPH_C1_KOITER_STEIGMANN_EQUATIONS_HEADER
#define OOMPH_C1_KOITER_STEIGMANN_EQUATIONS_HEADER

// std headers
#include <sstream>

// oomph-lib headers
#include "src/generic/nodes.h"
#include "src/generic/oomph_utilities.h"
#include "src/generic/Telements.h"

namespace oomph
{
  //=============================================================
  /// A class for all subparametric elements that solve the 2D-
  /// Biharmonic equations.
  /// \[
  /// \frac{\partial^4 u}{\partial x_i^4} = f(x_j)
  /// \f]
  /// This contains the generic maths. Shape functions, geometric
  /// mapping etc. must get implemented in derived class.
  //=============================================================
  class KoiterSteigmannEquations : public virtual FiniteElement
  {
  public:
    //----------------------------------------------------------------------
    // Types defined for class

    /// A pointer to a scalar function of the position. Can be used for
    /// out-of-plane forcing, swelling, isotropic-prestrain, etc.
    typedef void (*ScalarFctPt)(const Vector<double>& x, double& f);

    /// A pointer to a vector function of the position. Can be used for
    /// in-of-plane forcing, anisotropic-prestrain, etc.
    typedef void (*VectorFctPt)(const Vector<double>& x, Vector<double>& f);

    /// \short Function pointer to prestress function  e_0(x) --
    /// x is a Vector!
    typedef void (*PrestrainFctPt)(const Vector<double>& x,
                                   DenseMatrix<double>& prestrain);

    /// \short Function pointer to pressure function fct(x,f(x)) --
    /// x is a Vector!
    // typedef void (*PressureFctPt)(const Vector<double>& x, double& f);
    typedef void (*StressFctPt)(const Vector<double>& x,
                                const Vector<double>& u,
                                const DenseMatrix<double>& strain,
                                const DenseMatrix<double>& g_tensor,
                                DenseMatrix<double>& stress);

    typedef void (*DStressFctPt)(const Vector<double>& x,
                                 const Vector<double>& u,
                                 const DenseMatrix<double>& strain,
                                 const DenseMatrix<double>& g_tensor,
                                 RankThreeTensor<double>& d_stress_dui,
                                 RankFourTensor<double>& d_stress_dstrain);

    /// \short Function pointer to gradient of pressure function  fct(x,g(x)) --
    /// x is a Vector!
    typedef void (*PressureVectorFctPt)(const Vector<double>& x,
                                        const Vector<double>& u,
                                        const DenseMatrix<double>& grad_u,
                                        const Vector<double>& n,
                                        Vector<double>& pressure_vector);

    /// \short Function pointer to gradient of pressure function  fct(x,g(x)) --
    /// x is a Vector!
    typedef void (*DPressureVectorFctPt)(
      const Vector<double>& x,
      const Vector<double>& u,
      const DenseMatrix<double>& grad_u,
      const Vector<double>& n,
      DenseMatrix<double>& d_pressure_vector);

    /// \short Function pointer to gradient of pressure function  fct(x,g(x)) --
    /// x is a Vector!
    typedef void (*DPressureVectorDMatrixFctPt)(
      const Vector<double>& x,
      const Vector<double>& u,
      const DenseMatrix<double>& grad_u,
      const Vector<double>& n,
      RankThreeTensor<double>& d_pressure_vector);

    /// \short Function pointer to the Error Metric we are using
    ///  e.g could be that we are just interested in error on w etc.
    typedef void (*ErrorMetricFctPt)(const Vector<double>& x,
                                     const Vector<double>& u,
                                     const Vector<double>& u_exact,
                                     double& error,
                                     double& norm);

    /// \shot Function pointer to the Error Metric we are using if we want
    /// multiple
    ///  errors.  e.g could be we want errors seperately on each displacment
    typedef void (*MultipleErrorMetricFctPt)(const Vector<double>& x,
                                             const Vector<double>& u,
                                             const Vector<double>& u_exact,
                                             Vector<double>& error,
                                             Vector<double>& norm);

    // End of types defined for class
    //----------------------------------------------------------------------


    //--------------------------------------------------------------------------
    // Class construction

    /// Constructor
    KoiterSteigmannEquations()
      : Use_finite_difference_jacobian(false),
	Ui_is_damped{false,false,false},
	Pressure_fct_pt(0),
        D_pressure_dr_fct_pt(0),
        D_pressure_dn_fct_pt(0),
        D_pressure_d_grad_u_fct_pt(0),
	Prestrain_fct_pt(0),
        Stress_fct_pt(0),
        D_stress_fct_pt(0),
        Error_metric_fct_pt(0),
        Multiple_error_metric_fct_pt(0)
    {
      // Poisson ratio is 0.5 (incompressible) by default.
      Nu_pt = &Default_Nu_Value;
      Mu_pt = &Default_Mu_Value;
      Eta_u_pt = &Default_Eta_Value;
      Eta_sigma_pt = &Default_Eta_Value;

      // Thickness ratio small (0.001). Might expect thickness independence in
      // this limit depending on the magnitude of the external forcing
      Thickness_pt = &Default_Thickness_Value;

      // [zdec] check these comments
      // Default parameters should lead to membrane like behaviour by default:
      // order 1 displacements for small thickness sheet
      // By default the displacements are asumed to be of order 1
    }

    /// Broken copy constructor
    KoiterSteigmannEquations(const KoiterSteigmannEquations& dummy)
    {
      BrokenCopy::broken_copy("KoiterSteigmannEquations");
    }

    /// Broken assignment operator
    void operator=(const KoiterSteigmannEquations&)
    {
      BrokenCopy::broken_assign("KoiterSteigmannEquations");
    }


    //----------------------------------------------------------------------
    // Output and documentation

    /// Output with default number of plot points
    void output(std::ostream& outfile)
    {
      const unsigned n_plot = 5;
      KoiterSteigmannEquations::output(outfile, n_plot);
    }

    /// \short Output FE representation of soln: x,y,u or x,y,z,u at
    /// n_plot^DIM plot points
    void output(std::ostream& outfile, const unsigned& n_plot);

    /// C_style output with default number of plot points
    void output(FILE* file_pt)
    {
      const unsigned n_plot = 5;
      KoiterSteigmannEquations::output(file_pt, n_plot);
    }

    /// \short C-style output FE representation of soln: x,y,u or x,y,z,u at
    /// n_plot^DIM plot points
    void output(FILE* file_pt, const unsigned& n_plot);

    /// Output exact soln: x,y,u_exact or x,y,z,u_exact at n_plot^DIM plot
    /// points
    void output_fct(std::ostream& outfile,
                    const unsigned& n_plot,
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt);

    /// \short Output exact soln: x,y,u_exact or x,y,z,u_exact at
    /// n_plot^DIM plot points (dummy time-dependent version to
    /// keep intel compiler happy)
    void output_fct(std::ostream& outfile,
                    const unsigned& n_plot,
                    const double& time,
                    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
    {
      throw OomphLibError(
        "There is no time-dependent output_fct() for these elements ",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
    }


    //----------------------------------------------------------------------
    // Error and norms

    /// Get error against and norm of exact solution
    void compute_error(
      std::ostream& outfile,
      FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
      double& error,
      double& norm);

    /// Get error against and norm of exact solution
    void compute_error(
      std::ostream& outfile,
      FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
      Vector<double>& error,
      Vector<double>& norm);

    /// Dummy, time dependent error checker
    inline void compute_error(
      std::ostream& outfile,
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


    //----------------------------------------------------------------------
    // Dependent variables

    // Return the determinant of the metric tensor
    inline double two_by_two_determinant(const DenseMatrix<double>& g_tensor)
    {
      // Now return determinant
      return g_tensor(0, 0) * g_tensor(1, 1) - g_tensor(0, 1) * g_tensor(1, 0);
    }

    /// Return FE representation of unknown values u(s) at local coordinate s
    inline void interpolated_koiter_steigmann_disp(
      const Vector<double>& s, Vector<Vector<double>>& interpolated_u) const
    {
      // Number of displacement fields we are interpolating
      const unsigned n_displacements = this->Number_of_displacements;

      // Dimension of the plate domain
      const unsigned dim = this->dim();
      // The number of first derivatives is the dimension of the element
      const unsigned n_deriv = dim;
      // The number of second derivatives is the triangle number of the
      // dimension
      const unsigned n_2deriv = dim * (dim + 1) / 2;

      // Find out how many nodes are used to interpolate each displacement field
      const unsigned n_u_node = nu_node();

      // Vector of nodal indices used in interpolating displacement fields
      const Vector<unsigned> u_nodes = get_u_node_indices();

      // Find number of position dofs at each node per displacement field
      const unsigned n_u_nodal_type = nu_type_at_each_node();

      // Find the number of internal dofs per displacement field
      const unsigned n_u_internal_type = nu_type_internal();

      // [zdec] REMOVE TEST FUNCTIONS
      // Local c1-shape funtion
      Shape psi_n(n_u_node, n_u_nodal_type);
      Shape psi_i(n_u_internal_type);
      DShape dpsi_n_dxi(n_u_node, n_u_nodal_type, n_deriv);
      DShape dpsi_i_dxi(n_u_internal_type, n_deriv);
      DShape d2psi_n_dxi2(n_u_node, n_u_nodal_type, n_2deriv);
      DShape d2psi_i_dxi2(n_u_internal_type, n_2deriv);

      // Find values of c1-shape function
      d2basis_u_eulerian_koiter_steigmann(s,
					  psi_n,
					  psi_i,
					  dpsi_n_dxi,
					  dpsi_i_dxi,
					  d2psi_n_dxi2,
					  d2psi_i_dxi2);

      // Interpolated unknown
      // Loop over displacements
      for (unsigned i_field = 0; i_field < n_displacements; i_field++)
      {
        // Loop over nodes
        for (unsigned j_node = 0; j_node < n_u_node; j_node++)
        {
          // Loop over hermite dofs
          for (unsigned k_type = 0; k_type < n_u_nodal_type; k_type++)
          {
            // Get the kth nodal value at node j for displacement i
            double u_value =
	      get_u_i_value_at_node_of_type(i_field, j_node, k_type);
            interpolated_u[i_field][0] += u_value * psi_n(j_node, k_type);
            interpolated_u[i_field][1] += u_value * dpsi_n_dxi(j_node, k_type, 0);
            interpolated_u[i_field][2] += u_value * dpsi_n_dxi(j_node, k_type, 1);
            interpolated_u[i_field][3] += u_value * d2psi_n_dxi2(j_node, k_type, 0);
            interpolated_u[i_field][4] += u_value * d2psi_n_dxi2(j_node, k_type, 1);
            interpolated_u[i_field][5] += u_value * d2psi_n_dxi2(j_node, k_type, 2);
          }
        }
        // Loop over internal dofs
	for (unsigned k_type = 0; k_type < n_u_internal_type; k_type++)
	{
	  double u_value = get_u_i_internal_value_of_type(i_field, k_type);
	  interpolated_u[i_field][0] += u_value * psi_i(k_type);
	  interpolated_u[i_field][1] += u_value * dpsi_i_dxi(k_type, 0);
	  interpolated_u[i_field][2] += u_value * dpsi_i_dxi(k_type, 1);
	  interpolated_u[i_field][3] += u_value * d2psi_i_dxi2(k_type, 0);
	  interpolated_u[i_field][4] += u_value * d2psi_i_dxi2(k_type, 1);
	  interpolated_u[i_field][5] += u_value * d2psi_i_dxi2(k_type, 2);
	}
      }
    }

    // Fill in the Kirchhoff St Venant stress tensor
    inline void fill_in_kirchhoff_st_venant_stress(
      const DenseMatrix<double>& strain, DenseMatrix<double>& stress)
    {
      // Zero the stress
      for (unsigned alpha = 0; alpha < 2; ++alpha)
      {
        for (unsigned beta = 0; beta < 2; ++beta)
        {
          stress(alpha, beta) = 0.0;
        }
      }

      // Get Poisson ratio
      double nu = get_nu();
      // Loop over alpha and beta
      for (unsigned alpha = 0; alpha < 2; ++alpha)
      {
        for (unsigned beta = 0; beta < 2; ++beta)
        {
          stress(alpha, beta) += (1 - nu) * strain(alpha, beta) / (1 - nu * nu);
          stress(alpha, alpha) += nu * strain(beta, beta) / (1 - nu * nu);
        }
      }
    }

    // Fill in the Kirchhoff St Venant stress tensor
    inline void fill_in_d_kirchhoff_st_venant_stress_du_unknown(
      const RankFourTensor<double>& d_strain_du_unknown,
      RankFourTensor<double>& d_stress_du_unknown)
    {
      // Zero the vector
      for (unsigned alpha = 0; alpha < 2; ++alpha)
      {
        for (unsigned beta = 0; beta < 2; ++beta)
        {
          for (unsigned i = 0; i < Number_of_displacements; ++i)
          {
            d_stress_du_unknown(alpha, beta, i, 0) = 0.0;
            d_stress_du_unknown(alpha, beta, i, 1) = 0.0;
            d_stress_du_unknown(alpha, beta, i, 2) = 0.0;
          }
        }
      }

      // Get Poisson ratio
      double nu = get_nu();
      // Loop over alpha and beta
      for (unsigned alpha = 0; alpha < 2; ++alpha)
      {
        for (unsigned beta = 0; beta < 2; ++beta)
        {
          for (unsigned i = 0; i < Number_of_displacements; ++i)
          {
            for (unsigned gamma = 0; gamma < 2; ++gamma)
            {
              d_stress_du_unknown(alpha, beta, i, 1 + gamma) +=
                (1 - nu) * d_strain_du_unknown(alpha, beta, i, gamma) /
                (1 - nu * nu);
              d_stress_du_unknown(alpha, alpha, i, 1 + gamma) +=
                nu * d_strain_du_unknown(beta, beta, i, gamma) / (1 - nu * nu);
            }
          }
        }
      }
    }

    // Get the (Green Lagrange) strain tensor
    // \Gamma_{\alpha\beta\gamma} = ( E_{\alpha\beta,gamma} +
    // E_{\gamma\alpha,\beta}
    //   - E_{\beta\gamma,\alpha})
    inline void fill_in_stress_tensor(const Vector<double>& x,
                                      const Vector<double>& r,
                                      const DenseMatrix<double>& strain,
                                      const DenseMatrix<double>& g_tensor,
                                      DenseMatrix<double>& stress)
    {
      // Zero the stress
      for (unsigned alpha = 0; alpha < 2; ++alpha)
      {
        for (unsigned beta = 0; beta < 2; ++beta)
        {
          stress(alpha, beta) = 0.0;
        }
      }
      // IF no stress function has been set use Kirchhoff st venant
      if (Stress_fct_pt == 0)
      {
        fill_in_kirchhoff_st_venant_stress(strain, stress);
      }
      // Use the user defined one
      else
      {
        (*Stress_fct_pt)(x, r, strain, g_tensor, stress);
      }
    }

    // Get the (Green Lagrange) strain tensor
    // \Gamma_{\alpha\beta\gamma} = ( E_{\alpha\beta,gamma} +
    // E_{\gamma\alpha,\beta}
    //   - E_{\beta\gamma,\alpha})
    inline void fill_in_d_stress_tensor_du_unknown(
      const Vector<double>& x,
      const Vector<double>& r,
      const DenseMatrix<double>& strain,
      const DenseMatrix<double>& g_tensor,
      DenseMatrix<double>& stress,
      const RankFourTensor<double>& d_strain_dunknown,
      RankFourTensor<double>& d_stress_dunknown)
    {
      // Zero the vector
      for (unsigned alpha = 0; alpha < 2; ++alpha)
      {
        for (unsigned beta = 0; beta < 2; ++beta)
        {
          for (unsigned i = 0; i < Number_of_displacements; ++i)
          {
            // Flatpacked array ui ui,alpha
            for (unsigned k = 0; k < 3; ++k)
            {
              d_stress_dunknown(alpha, beta, i, k) = 0.0;
            }
          }
        }
      }

      // IF not set use Kirchhoff st venant
      if (Stress_fct_pt == 0 && D_stress_fct_pt == 0)
      {
        fill_in_d_kirchhoff_st_venant_stress_du_unknown(d_strain_dunknown,
                                                        d_stress_dunknown);
      }
      // Use the user defined one
      else if (D_stress_fct_pt != 0)
      {
        // Initialise the tensor (SLOW)
        RankFourTensor<double> d_stress_d_epsilon(2, 2, 2, 2, 0.0);
        RankThreeTensor<double> d_stress_du(2, 2, Number_of_displacements, 0.0);
        // Get the user defined tensor
        (*D_stress_fct_pt)(
          x, r, strain, g_tensor, d_stress_du, d_stress_d_epsilon);
        for (unsigned i = 0; i < Number_of_displacements; ++i)
        {
          for (unsigned alpha = 0; alpha < 2; ++alpha)
          {
            for (unsigned beta = 0; beta < 2; ++beta)
            {
              d_stress_dunknown(alpha, beta, i, 0) +=
                d_stress_du(alpha, beta, i);
              for (unsigned gamma = 0; gamma < 2; ++gamma)
              {
                for (unsigned delta = 0; delta < 2; ++delta)
                {
                  for (unsigned mu = 0; mu < 2; ++mu)
                  {
                    // Fill in using the rank four tensor (SLOW)
                    d_stress_dunknown(alpha, beta, i, 1 + mu) +=
                      d_stress_d_epsilon(alpha, beta, gamma, delta) *
                      d_strain_dunknown(gamma, delta, i, mu);
                  }
                }
              }
            }
          }
        }
      }
      else
      {
        throw OomphLibError(
          "You have defined a function for stress but not d_stress_depsilon. A finite\
difference stress is not yet defined so an error has been thrown.",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
        // USE THE USER DEFINED ONE FILL IN
        //     (*D_stress_fct_pt)(strain,stress);
        /// HERE FINITE DIFFERENCE + WARNING
      }
    }

    // Return the interpolated unit normal
    inline void fill_in_metric_tensor(
      const DenseMatrix<double>& interpolated_drdxi,
      DenseMatrix<double>& g_tensor)
    {
      // Zero the metric tensor
      for (unsigned alpha = 0; alpha < 2; ++alpha)
      {
        for (unsigned beta = 0; beta < 2; ++beta)
        {
          g_tensor(alpha, beta) = 0.0;
        }
      }
      // Fill in metric tensor
      // Loop over plane tangent vectors
      for (unsigned alpha = 0; alpha < 2; ++alpha)
      {
        // Loop over plane tangent vectors
        for (unsigned beta = 0; beta < 2; ++beta)
        {
          // LSum over components of tangent vectors
          for (unsigned i = 0; i < this->Number_of_displacements; ++i)
          {
            g_tensor(alpha, beta) +=
              interpolated_drdxi(i, alpha) * interpolated_drdxi(i, beta);
          }
        }
      }
    }

    // Get the tangent vectors
    inline void fill_in_tangent_vectors(
      const DenseMatrix<double>& interpolated_dudxi,
      DenseMatrix<double>& interpolated_drdxi)
    {
      // Zero the tangent vectors
      for (unsigned i = 0; i < this->Number_of_displacements; ++i)
      {
        for (unsigned beta = 0; beta < 2; ++beta)
        {
          interpolated_drdxi(i, beta) = 0.0;
        }
      }
      // Fill in the tangent vectors
      for (unsigned i = 0; i < this->Number_of_displacements; ++i)
      {
        for (unsigned beta = 0; beta < 2; ++beta)
        {
          // Scale the displacements
          const double eta_u = this->eta_u();
          // The displacement part
          interpolated_drdxi(i, beta) += eta_u * interpolated_dudxi(i, beta);
          // The coordinate part
          interpolated_drdxi(i, beta) += (double)(i == beta);
        }
      }
    }

    // Return the determinant of the metric tensor
    inline double metric_tensor_determinant(
      const DenseMatrix<double>& interpolated_drdxi)
    {
      // Intitialise
      DenseMatrix<double> g_tensor(2, 2, 0.0);
      // Fill in metric tensor
      fill_in_metric_tensor(interpolated_drdxi, g_tensor);
      // Now return determinant
      return g_tensor(0, 0) * g_tensor(1, 1) - g_tensor(0, 1) * g_tensor(1, 0);
    }

    // Return the determinant of the metric tensor
    inline void d_metric_tensor_determinant_du_unknown(
      const DenseMatrix<double>& g_tensor,
      const RankFourTensor<double>& d_g_tensor_du_unknown,
      DenseMatrix<double>& d_metric_tensor_determinant_du_unknown)
    {
      // // Adjugate metric tensor

      // adj_g(0,0)= g_tensor(1,1);
      // adj_g(1,1)= g_tensor(0,0);
      // adj_g(0,1)=-g_tensor(0,1);
      // adj_g(1,0)=-g_tensor(1,0);

      // Loop over derivatives
      for (unsigned i = 0; i < this->Number_of_displacements; ++i)
      {
        d_metric_tensor_determinant_du_unknown(i, 0) = 0.0;
        d_metric_tensor_determinant_du_unknown(i, 1) = 0.0;
        for (unsigned gamma = 0; gamma < 2; ++gamma)
        {
          for (unsigned delta = 0; delta < 2; ++delta)
          {
            // Derivative: D(det A)(t) = adj(A) : D(A)(t)
            // d_metric_tensor_determinant_du_unknown[i]+=adj_g(gamma,delta)
            //   * d_g_tensor_du_unknown(gamma,delta,i) ;
            // But we can write the adjugate like this (as it is symmetric)
            for (unsigned mu = 0; mu < 2; ++mu)
            {
              d_metric_tensor_determinant_du_unknown(i, mu) +=
                g_tensor((gamma + 1) % 2, (delta + 1) % 2) *
                d_g_tensor_du_unknown(gamma, delta, i, mu) *
                (gamma != delta ? -1 : 1);
            }
          }
        }
      }
    }

    // Return the interpolated unit normal
    inline void fill_in_unit_normal(
      const DenseMatrix<double>& interpolated_drdxi, Vector<double>& normal)
    {
      // Determinant of metric tensor
      const double g = metric_tensor_determinant(interpolated_drdxi);
      const double sqrt_g = sqrt(g);

      // Zero the normal
      for (unsigned i = 0; i < this->Number_of_displacements; ++i)
      {
        normal[i] = 0.0;
      }

      // Loop over plane tangent vectors
      for (unsigned alpha = 0; alpha < 2; ++alpha)
      {
        // Beta != Alpha
        const unsigned beta = (alpha + 1) % 2;
        Vector<double> g_alpha(3), g_beta(3);

        // Loop over displacement components
        for (unsigned i = 0; i < this->Number_of_displacements; ++i)
        {
          // Now fill in the tangent vectors
          g_alpha[i] = interpolated_drdxi(i, alpha);
          g_beta[i] = interpolated_drdxi(i, beta);
        }
        // Compute the cross product
        Vector<double> tmp(this->Number_of_displacements, 0.0);
        VectorHelpers::cross(g_alpha, g_beta, tmp);
        // And fill in the normal
        for (unsigned i = 0; i < this->Number_of_displacements; ++i)
        {
          normal[i] += tmp[i] / (2 * sqrt_g) * (alpha == 0 ? 1 : -1);
        }
      }
    }

    // Return the interpolated unit normal
    inline void fill_in_d_unit_normal_du_unknown(
      const DenseMatrix<double>& interpolated_drdxi,
      const DenseMatrix<double>& g_tensor,
      const RankFourTensor<double>& d_g_tensor_du_unknown,
      RankThreeTensor<double>& d_normal_du_unknown)
    {
      // Loop over plane tangent vectors
      for (unsigned i = 0; i < this->Number_of_displacements; ++i)
      {
        for (unsigned j = 0; j < this->Number_of_displacements; ++j)
        {
          d_normal_du_unknown(i, j, 0) = 0.0;
          d_normal_du_unknown(i, j, 1) = 0.0;
        }
      }

      // Determinant of metric tensor
      const double g(
        two_by_two_determinant(
          g_tensor));
      const double sqrt_g_inv(1 / sqrt(g));
      const double sqrt_g_cubed_inv(1 / std::pow(g, 1.5));

      // Get the metric tensor determinant
      DenseMatrix<double> d_g_du_unknown(this->Number_of_displacements, 2, 0.0);
      d_metric_tensor_determinant_du_unknown(
        g_tensor, d_g_tensor_du_unknown, d_g_du_unknown);

      // Loop over plane tangent vectors
      double tmp;
      for (unsigned j = 0; j < this->Number_of_displacements; ++j)
      {
        // Eta parameter
        const double eta_u = this->eta_u();
        for (unsigned alpha = 0; alpha < 2; ++alpha)
        {
          // Beta != Alpha
          const unsigned beta = (alpha + 1) % 2, k = (j + 1) % 3,
                         i = (j + 2) % 3;

          // Cross product between g_alpha and d_g_beta_du_unknown
          // NB:  dg_beta_duj_unknown[j]  = d2uidx_du_unknown[beta]
          // Error could be here somewhere
          tmp =
            interpolated_drdxi(i, alpha) * eta_u; // * d2uidx_du_unknown[beta];
          d_normal_du_unknown(k, j, beta) +=
            tmp * (sqrt_g_inv) * (alpha == 0 ? 1 : -1);

          tmp =
            interpolated_drdxi(k, beta) * eta_u; // *  d2uidx_du_unknown[alpha];
          d_normal_du_unknown(i, j, alpha) +=
            tmp * (sqrt_g_inv) * (alpha == 0 ? 1 : -1);

          // Loop over displacement components
          // And fill in the normal
          for (unsigned ii = 0; ii < this->Number_of_displacements; ++ii)
          {
            // The cross product will be epsilon_ijk vj vk
            const unsigned jj = (ii + 1) % 3, kk = (ii + 2) % 3;
            // Error could be in 'tmp'
            // Cross product between g_alpha and g_beta
            tmp = interpolated_drdxi(jj, alpha) * interpolated_drdxi(kk, beta) -
                  interpolated_drdxi(kk, alpha) * interpolated_drdxi(jj, beta);
            for (unsigned gamma = 0; gamma < 2; ++gamma)
            {
              d_normal_du_unknown(ii, j, gamma) -=
                tmp / (4) * sqrt_g_cubed_inv * d_g_du_unknown(j, gamma) *
                (alpha == 0 ? 1 : -1);
            }
          }
        }
      }
    }


    /// Get the (Green Lagrange) strain tensor and add a prestrain
    inline void fill_in_strain_tensor(
      const DenseMatrix<double>& interpolated_dudxi,
      const DenseMatrix<double>& prestrain,
      DenseMatrix<double>& strain)
    {
      // Fill in the strain as usual
      fill_in_strain_tensor(interpolated_dudxi, strain);

      // Add the prestrain contributions
      for (unsigned alpha = 0; alpha < 2; ++alpha)
      {
        for (unsigned beta = 0; beta < 2; ++beta)
        {
          strain(alpha, beta) += prestrain(alpha,beta);
        }
      }
    }


    // Get the (Green Lagrange) strain tensor
    inline void fill_in_strain_tensor(
      const DenseMatrix<double>& interpolated_dudxi,
      DenseMatrix<double>& strain)
    {
      // Zero the strain
      for (unsigned alpha = 0; alpha < 2; ++alpha)
      {
        for (unsigned beta = 0; beta < 2; ++beta)
        {
          strain(alpha, beta) = 0.0;
        }
      }
      // Loop over alpha and beta
      for (unsigned alpha = 0; alpha < 2; ++alpha)
      {
        for (unsigned beta = 0; beta < 2; ++beta)
        {
          // Fill in linear terms
          strain(alpha, beta) += interpolated_dudxi(alpha, beta) / 2.;
          strain(alpha, beta) += interpolated_dudxi(beta, alpha) / 2.;

          // Loop over displacements
          for (unsigned i = 0; i < this->Number_of_displacements; ++i)
          {
            // Scale the displacements
            const double eta_u = this->eta_u();
            // Nonlinear terms
            strain(alpha, beta) += eta_u * interpolated_dudxi(i, alpha) *
                                   interpolated_dudxi(i, beta) / 2.;
          }
        }
      }
    }

    // Get the (Green Lagrange) strain tensor
    inline void fill_in_d_g_tensor_du_unknown(
      const DenseMatrix<double>& interpolated_dudxi,
      RankFourTensor<double>& d_g_tensor_du_unknown)
    {
      // Loop over alpha and beta
      for (unsigned alpha = 0; alpha < 2; ++alpha)
      {
        for (unsigned beta = 0; beta < 2; ++beta)
        {
          // Zero the vector
          for (unsigned i = 0; i < this->Number_of_displacements; ++i)
          {
            d_g_tensor_du_unknown(alpha, beta, i, 0) = 0.0;
            d_g_tensor_du_unknown(alpha, beta, i, 1) = 0.0;
          }
          // Fill in linear terms relating to u_\alpha,\beta
          const double eta_u = this->eta_u();
          d_g_tensor_du_unknown(alpha, beta, beta, alpha) += eta_u;
          d_g_tensor_du_unknown(alpha, beta, alpha, beta) += eta_u;

          // Loop over displacements
          for (unsigned i = 0; i < this->Number_of_displacements; ++i)
          {
            // Scale the displacements
            const double eta_u = this->eta_u();
            // Nonlinear terms
            d_g_tensor_du_unknown(alpha, beta, i, alpha) +=
              eta_u * eta_u * interpolated_dudxi(i, beta);
            d_g_tensor_du_unknown(alpha, beta, i, beta) +=
              eta_u * eta_u * interpolated_dudxi(i, alpha);
          }
        }
      }
    }

    // Get the (Green Lagrange) strain tensor
    inline void fill_in_d_strain_tensor_du_unknown(
      const DenseMatrix<double>& interpolated_dudxi,
      RankFourTensor<double>& d_epsilon_tensor_du_unknown)
    {
      // Loop over alpha and beta
      for (unsigned alpha = 0; alpha < 2; ++alpha)
      {
        for (unsigned beta = 0; beta < 2; ++beta)
        {
          // Zero the vector
          for (unsigned i = 0; i < this->Number_of_displacements; ++i)
          {
            d_epsilon_tensor_du_unknown(alpha, beta, i, 0) = 0.0;
            d_epsilon_tensor_du_unknown(alpha, beta, i, 1) = 0.0;
          }
          // Fill in linear terms
          d_epsilon_tensor_du_unknown(alpha, beta, beta, alpha) += 1 / 2.;
          d_epsilon_tensor_du_unknown(alpha, beta, alpha, beta) += 1 / 2.;

          // Loop over displacements
          for (unsigned i = 0; i < this->Number_of_displacements; ++i)
          {
            // Scale the displacements
            const double eta_u = this->eta_u();
            // Nonlinear terms
            d_epsilon_tensor_du_unknown(alpha, beta, i, alpha) +=
              eta_u * interpolated_dudxi(i, beta) / 2.;
            d_epsilon_tensor_du_unknown(alpha, beta, i, beta) +=
              eta_u * interpolated_dudxi(i, alpha) / 2.;
          }
        }
      }
    }

    // Fill in the curvature tensor
    inline void fill_in_curvature_tensor(
      const Vector<double>& unit_normal,
      const DenseMatrix<double>& interpolated_d2rdxi2,
      DenseMatrix<double>& curvature)
    {
      // Zero the curvature
      for (unsigned alpha = 0; alpha < 2; ++alpha)
      {
        for (unsigned beta = 0; beta < 2; ++beta)
        {
          curvature(alpha, beta) = 0.0;
        }
      }
      // Loop over displacements and then coordinates
      for (unsigned i = 0; i < this->Number_of_displacements; ++i)
      {
        for (unsigned alpha = 0; alpha < 2; ++alpha)
        {
          for (unsigned beta = 0; beta < 2; ++beta)
          {
            // Scale the displacements
            curvature(alpha, beta) +=
              unit_normal[i] * interpolated_d2rdxi2(i, alpha + beta);
          }
        }
      }
    }

    // Fill in the curvature tensor
    inline void fill_in_d_curvature_tensor_du_unknown(
      const Vector<double>& unit_normal,
      const DenseMatrix<double>& interpolated_d2rdxi2,
      const RankThreeTensor<double>& d_unit_normal_dunknown,
      RankFourTensor<double>& d_curvature_du_unknown)
    {
      // Zero the vector
      for (unsigned alpha = 0; alpha < 2; ++alpha)
      {
        for (unsigned beta = 0; beta < 2; ++beta)
        {
          for (unsigned i = 0; i < this->Number_of_displacements; ++i)
          {
            // Flatpacked data
            for (unsigned k = 0; k < 5; ++k)
            {
              d_curvature_du_unknown(alpha, beta, i, k) = 0.0;
            }
          }
        }
      }
      // Loop over displacements and then coordinates
      for (unsigned i = 0; i < this->Number_of_displacements; ++i)
      {
        // Scale the displacements
        for (unsigned alpha = 0; alpha < 2; ++alpha)
        {
          for (unsigned beta = 0; beta < 2; ++beta)
          {
            d_curvature_du_unknown(alpha, beta, i, 2 + alpha + beta) +=
              unit_normal[i];
            for (unsigned j = 0; j < this->Number_of_displacements; ++j)
            {
              // Note here i is the component (that and we have eta_u_i) and j
              // is unknown
              for (unsigned gamma = 0; gamma < 2; ++gamma)
                d_curvature_du_unknown(alpha, beta, j, gamma) +=
                  d_unit_normal_dunknown(i, j, gamma) *
                  interpolated_d2rdxi2(i, alpha + beta);
            }
          }
        }
      }
    }

    // Fill in the moment tensor
    inline void fill_in_moment_tensor(
      const DenseMatrix<double>& interpolated_drxi,
      const Vector<double>& unit_normal,
      const DenseMatrix<double>& curvature,
      RankThreeTensor<double>& moment)
    {
      // Zero the stress
      for (unsigned i = 0; i < this->Number_of_displacements; ++i)
      {
        for (unsigned alpha = 0; alpha < 2; ++alpha)
        {
          for (unsigned beta = 0; beta < 2; ++beta)
          {
            moment(i, alpha, beta) = 0.0;
          }
        }
      }
      // Get thickness
      const double h = this->get_thickness(), nu = this->get_nu();
      // Loop over displacements and then coordinates
      for (unsigned i = 0; i < this->Number_of_displacements; ++i)
      {
        for (unsigned alpha = 0; alpha < 2; ++alpha)
        {
          for (unsigned beta = 0; beta < 2; ++beta)
          {
            // Get moment i
            moment(i, alpha, beta) += h * h * (1 - nu) * unit_normal[i] *
                                      curvature(alpha, beta) /
                                      (12. * (1 - nu * nu));
            moment(i, alpha, alpha) += h * h * nu * unit_normal[i] *
                                       curvature(beta, beta) /
                                       (12. * (1 - nu * nu));
          }
        }
      }
    }

    // Fill in the moment tensor
    inline void fill_in_d_moment_tensor_du_unknown(
      const DenseMatrix<double>& interpolated_d2rxi2,
      const Vector<double>& unit_normal,
      const DenseMatrix<double>& curvature,
      const RankThreeTensor<double>& d_unit_normal_du_unknown,
      const RankFourTensor<double>& d_curvature_du_unknown,
      RankFiveTensor<double>& d_moment_du_unknown)
    {
      // Get thickness
      const double h = this->get_thickness();
      // Get Poisson ratio
      const double nu = this->get_nu();
      // Zero the moment derivative tensor
      for (unsigned i = 0; i < this->Number_of_displacements; ++i)
      {
        for (unsigned alpha = 0; alpha < 2; ++alpha)
        {
          for (unsigned beta = 0; beta < 2; ++beta)
          {
            for (unsigned j = 0; j < this->Number_of_displacements; ++j)
            {
              // Flatpacked first and second derivatives
              for (unsigned k = 0; k < 5; ++k)
              {
                d_moment_du_unknown(i, alpha, beta, j, k) = 0.0;
              }
            }
          }
        }
      }
      // Loop over displacements and then coordinates
      for (unsigned i = 0; i < this->Number_of_displacements; ++i)
      {
        for (unsigned alpha = 0; alpha < 2; ++alpha)
        {
          for (unsigned beta = 0; beta < 2; ++beta)
          {
            // Loop over displacements and then coordinates
            for (unsigned j = 0; j < this->Number_of_displacements; ++j)
            {
              // Get moment i
              // Construct the tensor
              for (unsigned k = 0; k < 5; ++k)
              {
                // First terms
                d_moment_du_unknown(i, alpha, beta, j, k) +=
                  (1 - nu) * h * h * unit_normal[i] *
                  d_curvature_du_unknown(alpha, beta, j, k) /
                  (12. * (1 - nu * nu));
                d_moment_du_unknown(i, alpha, alpha, j, k) +=
                  nu * h * h * unit_normal[i] *
                  d_curvature_du_unknown(beta, beta, j, k) /
                  (12. * (1 - nu * nu));
              }
              // Second terms
              for (unsigned gamma = 0; gamma < 2; ++gamma)
              {
                d_moment_du_unknown(i, alpha, beta, j, gamma) +=
                  (1 - nu) * h * h * d_unit_normal_du_unknown(i, j, gamma) *
                  curvature(alpha, beta) / (12. * (1 - nu * nu));
                d_moment_du_unknown(i, alpha, alpha, j, gamma) +=
                  nu * h * h * d_unit_normal_du_unknown(i, j, gamma) *
                  curvature(beta, beta) / (12. * (1 - nu * nu));
              }
            }
          }
        }
      }
    }

    inline void fill_in_total_tension(
      const DenseMatrix<double>& stress,
      const RankThreeTensor<double>& christoffel_tensor,
      const RankThreeTensor<double>& moment_tensors,
      const DenseMatrix<double>& interpolated_dudxi,
      const DenseMatrix<double>& interpolated_d2rdxi2,
      DenseMatrix<double>& tension_vectors)
    {
      // Zero the tension
      for (unsigned i = 0; i < this->Number_of_displacements; ++i)
      {
        // Loop over inplane components
        for (unsigned alpha = 0; alpha < 2; ++alpha)
        {
          tension_vectors(i, alpha) = 0.0;
        }
      }

      double delta_ibeta;
      // Loop over displacement components
      for (unsigned i = 0; i < this->Number_of_displacements; ++i)
      {
        // Loop over inplane components
        for (unsigned alpha = 0; alpha < 2; ++alpha)
        {
          // Loop over inplane components
          for (unsigned beta = 0; beta < 2; ++beta)
          {
            // Kronecker Delta: delta_{i\beta}
            delta_ibeta = (i == beta ? 1.0 : 0.0);
            // Local copy of eta_u and eta_sigma
            const double eta_u = this->eta_u();
            const double eta_sigma = this->eta_sigma();
            // The diagonal parts of the tangent matrix
            tension_vectors(i, alpha) +=
              (delta_ibeta + eta_u * interpolated_dudxi(i, beta)) * eta_sigma *
              stress(alpha, beta);

            for (unsigned gamma = 0; gamma < 2; ++gamma)
            {
              // Shouldn't this be -M_{i \al \be} \Ga^\ga_{\al \be} ?
              tension_vectors(i, alpha) -=
                eta_u * moment_tensors(i, beta, gamma) *
                christoffel_tensor(alpha, beta, gamma);
            }
          }
        }
      }
    }

    inline void d_fill_in_total_tension_du_unknown(
      const DenseMatrix<double>& stress,
      const RankThreeTensor<double>& christoffel_tensor,
      const RankThreeTensor<double>& moment_tensors,
      const DenseMatrix<double>& interpolated_dudxi,
      const DenseMatrix<double>& interpolated_d2udxi2,
      RankFourTensor<double>& d_stress_du_unknown,
      const RankFiveTensor<double>& d_christoffel_tensor_du_unknown,
      const RankFiveTensor<double>& d_moment_tensors_du_unknown,
      RankFourTensor<double>& d_tension_vectors_du_unknown)
    {
      // Zero the derivatives of tension vectors
      for (unsigned i = 0; i < this->Number_of_displacements; ++i)
      {
        for (unsigned beta = 0; beta < 2; ++beta)
        {
          for (unsigned j = 0; j < this->Number_of_displacements; ++j)
          {
            for (unsigned k = 0; k < 6; ++k)
            {
              d_tension_vectors_du_unknown(i, beta, j, k) = 0.0;
            }
          }
        }
      }

      // Loop over displacement components
      for (unsigned i = 0; i < this->Number_of_displacements; ++i)
      {
        // Because we don't use a d_tangent_dunknown
        //  const double eta_g =(i==2 ? eta_g_z() : eta_g_xy());
        // Loop over inplane components
        double delta_ibeta;
        for (unsigned gamma = 0; gamma < 2; ++gamma)
        {
          // Loop over inplane components
          for (unsigned beta = 0; beta < 2; ++beta)
          {
            // Kronecker Delta: delta_{i\beta}
            delta_ibeta = (i == beta ? 1.0 : 0.0);
            // Local copy
            const double eta_u = this->eta_u();
            const double eta_sigma = this->eta_sigma();
            // Fill in dtension_du_unknown
            d_tension_vectors_du_unknown(i, gamma, i, 1 + beta) +=
              eta_sigma * eta_u * stress(gamma, beta);
            for (unsigned j = 0; j < this->Number_of_displacements; ++j)
            {
              // The diagonal parts of the tangent matrix
              for (unsigned mu = 0; mu < 3; ++mu)
              {
                d_tension_vectors_du_unknown(i, gamma, j, mu) +=
                  eta_sigma * d_stress_du_unknown(gamma, beta, j, mu) *
                  (delta_ibeta + eta_u * interpolated_dudxi(i, beta));
              }
            }
            for (unsigned alpha = 0; alpha < 2; ++alpha)
            {
              for (unsigned j = 0; j < this->Number_of_displacements; ++j)
              {
                for (unsigned k = 0; k < 5; ++k)
                {
                  d_tension_vectors_du_unknown(i, gamma, j, 1 + k) -=
                    eta_u * d_moment_tensors_du_unknown(i, alpha, beta, j, k) *
                    christoffel_tensor(gamma, alpha, beta);
                  d_tension_vectors_du_unknown(i, gamma, j, 1 + k) -=
                    eta_u * moment_tensors(i, alpha, beta) *
                    d_christoffel_tensor_du_unknown(gamma, alpha, beta, j, k);
                }
              }
            }
          }
        }
      }
    }

    // Get the Christtoffel tensor of the second kind
    // \Gamma_{\alpha\beta\gamma} = ( E_{\alpha\beta,gamma} +
    // E_{\gamma\alpha,\beta}
    //   - E_{\beta\gamma,\alpha})
    inline void fill_in_second_christoffel_tensor(
      const DenseMatrix<double>& interpolated_dudxi,
      const DenseMatrix<double>& interpolated_d2udxi2,
      RankThreeTensor<double>& gamma_tensor)
    {
      // Zero the christoffel
      for (unsigned gamma = 0; gamma < 2; ++gamma)
      {
        for (unsigned alpha = 0; alpha < 2; ++alpha)
        {
          for (unsigned beta = 0; beta < 2; ++beta)
          {
            gamma_tensor(alpha, beta, gamma) = 0.0;
          }
        }
      }

      // An intermediate step
      RankThreeTensor<double> strain_gradient(2, 2, 2, 0.0);

      // Loop over alpha and beta
      for (unsigned gamma = 0; gamma < 2; ++gamma)
      {
        for (unsigned alpha = 0; alpha < 2; ++alpha)
        {
          for (unsigned beta = 0; beta < 2; ++beta)
          {
            // Fill in linear terms
            strain_gradient(alpha, beta, gamma) +=
              interpolated_d2udxi2(alpha, beta + gamma) / 2. +
              interpolated_d2udxi2(beta, alpha + gamma) / 2.;
            // Loop over displacements
            for (unsigned i = 0; i < this->Number_of_displacements; ++i)
            {
              // Scale the displacements
              const double eta_u = this->eta_u();
              // Nonlinear terms
              strain_gradient(alpha, beta, gamma) +=
                eta_u *
                (interpolated_dudxi(i, alpha) *
                 interpolated_d2udxi2(i, beta + gamma)) /
                2.;
              strain_gradient(alpha, beta, gamma) +=
                eta_u *
                (interpolated_dudxi(i, beta) *
                 interpolated_d2udxi2(i, alpha + gamma)) /
                2.;
            }
          }
        }
      }

      // Fill in Christoffel
      for (unsigned alpha = 0; alpha < 2; ++alpha)
      {
        for (unsigned beta = 0; beta < 2; ++beta)
        {
          for (unsigned gamma = 0; gamma < 2; ++gamma)
          {
            gamma_tensor(alpha, beta, gamma) +=
              strain_gradient(alpha, beta, gamma) +
              strain_gradient(gamma, alpha, beta) -
              strain_gradient(gamma, beta, alpha);
          }
        }
      }
    }

    // Get the Christtoffel tensor of the second kind
    // \Gamma_{\alpha\beta\gamma} = ( E_{\alpha\beta,gamma} +
    // E_{\gamma\alpha,\beta}
    //   - E_{\beta\gamma,\alpha})
    inline void fill_in_d_second_christoffel_tensor_dui_unknown(
      const DenseMatrix<double>& interpolated_dudxi,
      const DenseMatrix<double>& interpolated_d2udxi2,
      RankFiveTensor<double>& d_gamma_tensor_du_unknown)
    {
      // An intermediate step
      RankFiveTensor<double> d_strain_gradient_du_unknown(2, 2, 2, 3, 5, 0.0);

      // Zero
      for (unsigned gamma = 0; gamma < 2; ++gamma)
      {
        for (unsigned alpha = 0; alpha < 2; ++alpha)
        {
          for (unsigned beta = 0; beta < 2; ++beta)
          {
            for (unsigned i = 0; i < this->Number_of_displacements; ++i)
            {
              for (unsigned k = 0; k < 5; ++k)
              {
                d_gamma_tensor_du_unknown(alpha, beta, gamma, i, k) = 0.0;
              }
            }
          }
        }
      }

      // Loop over alpha and beta
      for (unsigned gamma = 0; gamma < 2; ++gamma)
      {
        for (unsigned alpha = 0; alpha < 2; ++alpha)
        {
          for (unsigned beta = 0; beta < 2; ++beta)
          {
            // Fill in linear terms
            d_strain_gradient_du_unknown(
              alpha, beta, gamma, alpha, 2 + beta + gamma) += 1 / 2.;
            d_strain_gradient_du_unknown(
              alpha, beta, gamma, beta, 2 + alpha + gamma) += 1 / 2.;
            // Loop over displacements
            for (unsigned i = 0; i < this->Number_of_displacements; ++i)
            {
              // Scale the displacements
              const double eta_u = this->eta_u();
              // Nonlinear terms
              d_strain_gradient_du_unknown(alpha, beta, gamma, i, alpha) +=
                eta_u * (interpolated_d2udxi2(i, beta + gamma)) / 2.;
              d_strain_gradient_du_unknown(alpha, beta, gamma, i, beta) +=
                eta_u * (interpolated_d2udxi2(i, alpha + gamma)) / 2.;
              d_strain_gradient_du_unknown(
                alpha, beta, gamma, i, 2 + beta + gamma) +=
                eta_u * (interpolated_dudxi(i, alpha)) / 2.;
              d_strain_gradient_du_unknown(
                alpha, beta, gamma, i, 2 + alpha + gamma) +=
                eta_u * (interpolated_dudxi(i, beta)) / 2.;
            }
          }
        }
      }
      // Fill in Christoffel
      for (unsigned i = 0; i < this->Number_of_displacements; ++i)
      {
        for (unsigned alpha = 0; alpha < 2; ++alpha)
        {
          for (unsigned beta = 0; beta < 2; ++beta)
          {
            for (unsigned gamma = 0; gamma < 2; ++gamma)
            {
              for (unsigned k = 0; k < 5; ++k)
              {
                d_gamma_tensor_du_unknown(alpha, beta, gamma, i, k) +=
                  +d_strain_gradient_du_unknown(alpha, beta, gamma, i, k) +
                  d_strain_gradient_du_unknown(gamma, alpha, beta, i, k) -
                  d_strain_gradient_du_unknown(gamma, beta, alpha, i, k);
              }
            }
          }
        }
      }
    }


    //----------------------------------------------------------------------
    // Control parameters/forcing

    /// Enable damping in all displacements
    inline void enable_damping()
    {
      Ui_is_damped[0] = true;
      Ui_is_damped[1] = true;
      Ui_is_damped[2] = true;
    }

    /// Enable damping in the i-th displacement
    inline void enable_damping(const unsigned i_field)
    {
      // [zdec] Range checking?
      Ui_is_damped[i_field] = true;
    }

    /// Disable damping in all displacements
    inline void disable_damping()
    {
      Ui_is_damped[0] = false;
      Ui_is_damped[1] = false;
      Ui_is_damped[2] = false;
    }

    /// Disable damping in the i-th displacement
    inline void disable_damping(const unsigned i_field)
    {
      // [zdec] Range checking?
      Ui_is_damped[i_field] = false;
    }

    /// Get prestrain at (Eulerian) position x. This function is virtual to
    /// allow overloading in multi-physics problems where the strength of the
    /// prestrain function might be determined by another system of equations.
    inline virtual void get_prestrain(const unsigned& ipt,
				      const Vector<double>& x,
				      DenseMatrix<double>& prestrain) const
    {
      // Dimension of the plate
      unsigned dim = this->dim();

      // If no prestrain function has been set, return zero
      // Zero the pressure (as a precaution)
      for (unsigned alpha = 0; alpha < dim; ++alpha)
      {
	for (unsigned beta = 0; beta < dim; ++beta)
	{
	  prestrain(alpha, beta) = 0.0;
	}
      }
      if (Prestrain_fct_pt != 0)
      {
        // Get prestrain strength
        (*Prestrain_fct_pt)(x, prestrain);
      }
    }

    /// Get pressure term at (Eulerian) position x. This function is
    /// virtual to allow overloading in multi-physics problems where
    /// the strength of the pressure function might be determined by
    /// another system of equations.
    inline virtual void get_pressure(const unsigned& ipt,
                                     const Vector<double>& x,
				     const Vector<double>& u,
				     const DenseMatrix<double>& grad_u,
                                     const Vector<double>& n,
                                     Vector<double>& pressure) const
    {
      // If no pressure function has been set, return zero
      // Zero the pressure (as a precaution)
      pressure[0] = 0.0;
      pressure[1] = 0.0;
      pressure[2] = 0.0;
      if (Pressure_fct_pt != 0)
      {
        // Get pressure strength
        (*Pressure_fct_pt)(x, u, grad_u, n, pressure);
      }
    }

    /// Get pressure term at (Eulerian) position x. This function is
    /// virtual to allow overloading in multi-physics problems where
    /// the strength of the pressure function might be determined by
    /// another system of equations.
    inline virtual void get_d_pressure_dn(
      const unsigned& ipt,
      const Vector<double>& x,
      const Vector<double>& u,
      const DenseMatrix<double>& grad_u,
      const Vector<double>& n,
      DenseMatrix<double>& d_pressure_dn) const
    {
      // If no pressure function has been set, return zero
      // Zero the pressure (as a precaution)
      for (unsigned i = 0; i < Number_of_displacements; ++i)
      {
	for (unsigned j = 0; j < Number_of_displacements; ++j)
	{
	  d_pressure_dn(i, j) = 0.0;
	}
      }
      if (D_pressure_dn_fct_pt != 0)
      {
        // Get d_p_dn
        (*D_pressure_dn_fct_pt)(x, u, grad_u, n, d_pressure_dn);
      }
    }

    /// Get pressure term at (Eulerian) position x. This function is
    /// virtual to allow overloading in multi-physics problems where
    /// the strength of the pressure function might be determined by
    /// another system of equations.
    inline virtual void get_d_pressure_dr(
      const unsigned& ipt,
      const Vector<double>& x,
      const Vector<double>& u,
      const DenseMatrix<double>& grad_u,
      const Vector<double>& n,
      DenseMatrix<double>& d_pressure_dr) const
    {
      // If no pressure function has been set, return zero
      // Zero the pressure (as a precaution)
      for (unsigned i = 0; i < Number_of_displacements; ++i)
      {
	for (unsigned j = 0; j < Number_of_displacements; ++j)
	{
	  d_pressure_dr(i, j) = 0.0;
	}
      }
      if (D_pressure_dr_fct_pt != 0)
      {
        // Get d_p_dn
        (*D_pressure_dr_fct_pt)(x, u, grad_u, n, d_pressure_dr);
      }
    }

    /// Get pressure term at (Eulerian) position x. This function is
    /// virtual to allow overloading in multi-physics problems where
    /// the strength of the pressure function might be determined by
    /// another system of equations.
    inline void get_d_pressure_d_grad_u(
      const unsigned& ipt,
      const Vector<double>& x,
      const Vector<double>& u,
      const DenseMatrix<double>& grad_u,
      const Vector<double>& n,
      RankThreeTensor<double>& d_pressure_d_grad_u) const
    {
      // [zdec] Preserve just to check whether I am stupid
      // // If no pressure function has been set, return zero
      // if (D_pressure_dn_fct_pt == 0)
      // {
      //   for (unsigned i = 0; i < Number_of_displacements; ++i)
      //   {
      //     for (unsigned j = 0; j < Number_of_displacements; ++j)
      //     {
      //       for (unsigned alpha = 0; alpha < 2; ++alpha)
      //       {
      //         d_pressure_d_grad_u(i, j, alpha) = 0.0;
      //       }
      //     }
      //   }
      // }
      // else
      // {
      //   // Zero the pressure (as a precaution)
      //   for (unsigned i = 0; i < Number_of_displacements; ++i)
      //   {
      //     for (unsigned j = 0; j < Number_of_displacements; ++j)
      //     {
      //       for (unsigned alpha = 0; alpha < 2; ++alpha)
      //       {
      //         d_pressure_d_grad_u(i, j, alpha) = 0.0;
      //       }
      //     }
      //   }
      //   // Get d_p_dn
      //   (*D_pressure_d_grad_u_fct_pt)(x, u, grad_u, n, d_pressure_d_grad_u);
      // }

      // If no pressure function has been set, return zero
      // Zero the pressure (as a precaution)
      for (unsigned i = 0; i < Number_of_displacements; ++i)
      {
	for (unsigned j = 0; j < Number_of_displacements; ++j)
	{
	  for (unsigned alpha = 0; alpha < 2; ++alpha)
	  {
	    d_pressure_d_grad_u(i, j, alpha) = 0.0;
	  }
	}
      }
      if (D_pressure_d_grad_u_fct_pt != 0)
      {
        // Get d_p_dn
        (*D_pressure_d_grad_u_fct_pt)(x, u, grad_u, n, d_pressure_d_grad_u);
      }
    }


    //----------------------------------------------------------------------
    // Jacobian and residual contributions

    /// Enable Finite difference Jacobian
    void enable_finite_difference_jacobian()
    {
      Use_finite_difference_jacobian = true;
    }

    /// Disable Finite difference Jacobian
    void disable_finite_difference_jacobian()
    {
      Use_finite_difference_jacobian = false;
    }

    /// Return the index at which the ith displacement kth type (unknown) value
    /// is stored. In derived multi-physics elements, this function should be
    /// overloaded to reflect the chosen storage scheme. Note that these
    /// equations require that the unknown is always stored at the same index at
    /// each node, the inteface may .
    inline virtual unsigned nodal_index_of_u_i(const unsigned& i_field,
					       const unsigned& k_type) const
    {
      return i_field * (this->nu_type_at_each_node()) + k_type;
    }

    /// Add the element's contribution to its residual vector (wrapper)
    void fill_in_contribution_to_residuals(Vector<double>& residuals)
    {
      // Call the generic residuals function with flag set to 0
      // using a dummy matrix argument
      fill_in_generic_residual_contribution_koiter_steigmann(
        residuals, GeneralisedElement::Dummy_matrix, 0);
    }


    /// Add the element's contribution to its residual vector and
    /// element Jacobian matrix (wrapper)
    void fill_in_contribution_to_jacobian(Vector<double>& residuals,
                                          DenseMatrix<double>& jacobian)
    {
      // Call the generic routine with the flag set to 1
      if (!Use_finite_difference_jacobian)
      {
        fill_in_generic_residual_contribution_koiter_steigmann(
          residuals, jacobian, 1);
      }
      else
      {
        // Otherwise call the default
        FiniteElement::fill_in_contribution_to_jacobian(residuals, jacobian);
      }
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


    //----------------------------------------------------------------------
    // Misc

    /// Self-test: Return 0 for OK
    unsigned self_test();


    //--------------------------------------------------------------------------
    // Member data access functions

    /// Access function: Pointer to pressure function
    inline PressureVectorFctPt& pressure_fct_pt()
    {
      return Pressure_fct_pt;
    }

    /// Access function: Pointer to prestrain function
    inline PrestrainFctPt& prestrain_fct_pt()
    {
      return Prestrain_fct_pt;
    }

    /// Access function: Pointer to pressure function
    inline DPressureVectorFctPt& d_pressure_dn_fct_pt()
    {
      return D_pressure_dn_fct_pt;
    }

    /// Access function: Pointer to pressure function
    inline DPressureVectorFctPt& d_pressure_dr_fct_pt()
    {
      return D_pressure_dr_fct_pt;
    }

    /// Access function: Pointer to pressure function
    inline DPressureVectorDMatrixFctPt& d_pressure_d_grad_u_fct_pt()
    {
      return D_pressure_d_grad_u_fct_pt;
    }

    /// Access function: Pointer to Stress function
    inline StressFctPt& stress_fct_pt()
    {
      return Stress_fct_pt;
    }

    /// Access function: Pointer to Stress function
    inline DStressFctPt& d_stress_fct_pt()
    {
      return D_stress_fct_pt;
    }

    /// Access function: Pointer to error metric function
    inline ErrorMetricFctPt& error_metric_fct_pt()
    {
      return Error_metric_fct_pt;
    }

    /// Access function: Pointer to multiple error metric function
    inline MultipleErrorMetricFctPt& multiple_error_metric_fct_pt()
    {
      return Multiple_error_metric_fct_pt;
    }

    /// Access function: Pointer to pressure function. Const version
    inline PressureVectorFctPt pressure_fct_pt() const
    {
      return Pressure_fct_pt;
    }

    /// Access function: Pointer to pressure function
    inline DPressureVectorFctPt d_pressure_dn_fct_pt() const
    {
      return D_pressure_dn_fct_pt;
    }

    /// Access function: Pointer to pressure function
    inline DPressureVectorFctPt d_pressure_dr_fct_pt() const
    {
      return D_pressure_dr_fct_pt;
    }

    /// Access function: Pointer to pressure function
    inline DPressureVectorDMatrixFctPt d_pressure_d_grad_u_fct_pt() const
    {
      return D_pressure_d_grad_u_fct_pt;
    }

    /// Access function: Pointer to stress function. Const version
    inline StressFctPt stress_fct_pt() const
    {
      return Stress_fct_pt;
    }

    /// Access function: Pointer to stress function. Const version
    inline DStressFctPt d_stress_fct_pt() const
    {
      return D_stress_fct_pt;
    }

    /// Access function: Pointer to error metric function function
    inline ErrorMetricFctPt error_metric_fct_pt() const
    {
      return Error_metric_fct_pt;
    }

    /// Access function: Pointer to multiple error metric function
    inline MultipleErrorMetricFctPt multiple_error_metric_fct_pt() const
    {
      return Multiple_error_metric_fct_pt;
    }

    /// Access function to the Poisson ratio pointer.
    inline const double*& nu_pt()
    {
      return Nu_pt;
    }

    /// Access function to the Poisson ratio (const version)
    inline const double& get_nu() const
    {
      return *Nu_pt;
    }

    /// Access function to the coefficient of damping pointer
    inline const double*& mu_pt()
    {
      return Mu_pt;
    }

    /// Access function to the coefficient of damping
    inline const double& get_mu()
    {
      return *Mu_pt;
    }

    /// Access function to the thickness pointer.
    inline const double*& thickness_pt()
    {
      return Thickness_pt;
    }

    /// Access function to the thickness
    inline const double& get_thickness() const
    {
      return *Thickness_pt;
    }

    /// Access function to the displacment scaling.
    virtual const double*& eta_u_pt()
    {
      return Eta_u_pt;
    }

    /// Access function to the in plane displacment scaling in the displacement.
    virtual const double*& eta_sigma_pt()
    {
      return Eta_sigma_pt;
    }

    /// Access function to the displacement scaling (HERE always 1. - until
    /// valid.).
    virtual const double& eta_u() const
    {
      return *Eta_u_pt;
    }

    /// Access function to the in plane displacment scaling.
    virtual const double& eta_sigma() const
    {
      return *Eta_sigma_pt;
    }



  protected:
    //--------------------------------------------------------------------------
    // Pure virtual interfaces which must be implemented when geometry and bases
    // are added in the derived class

    /// (pure virtual) interface to return a vector of the index of each
    /// displacement unkonwn in the grander scheme of unknowns
    virtual Vector<unsigned> u_field_indices() const = 0;

    /// (pure virtual) interface to return a vector of the index of the u_i
    /// displacement unkonwn in the grander scheme of unknowns
    virtual unsigned u_i_field_index(const unsigned& i_field) const
    {
      return i_field;
    }

    /// (pure virtual) interface to return the number of nodes used to
    /// interpolate each displacement field
    virtual unsigned nu_node() const = 0;

    /// (pure virtual) interface to get the local indices of the nodes used by u
    virtual Vector<unsigned> get_u_node_indices() const = 0;

    /// (pure virtual) interface to get the number of basis types for u at node
    /// j
    virtual unsigned nu_type_at_each_node() const = 0;

    /// (pure virtual) interface to retrieve the value of u_alpha at node j of
    /// type k
    virtual double get_u_i_value_at_node_of_type(
      const unsigned& i_field,
      const unsigned& j_node,
      const unsigned& k_type) const = 0;

    /// (pure virtual) interface to retrieve the t-th history value value of
    /// u_alpha at node j of type k
    virtual double get_u_i_value_at_node_of_type(
      const unsigned& t_time,
      const unsigned& i_field,
      const unsigned& j_node,
      const unsigned& k_type) const = 0;

    /// (pure virtual) interface to get the pointer to the internal data used to
    /// interpolate u_i (NOTE: assumes each u field has exactly one internal
    /// data)
    virtual Data* u_i_internal_data_pt(
      const unsigned& i_field) const = 0;

    /// (pure virtual) interface to get the number of internal types for the u
    /// fields
    virtual unsigned nu_type_internal() const = 0;

    /// (pure virtual) interface to retrieve the value of u_alpha of internal
    /// type k
    virtual double get_u_i_internal_value_of_type(
      const unsigned& i_field,
      const unsigned& k_type) const = 0;

    /// (pure virtual) interface to retrieve the t-th history value of u_alpha
    /// of internal type k
    virtual double get_u_i_internal_value_of_type(
      const unsigned& t_time,
      const unsigned& i_field,
      const unsigned& k_type) const = 0;


    /// (pure virtual) Out-of-plane basis functions at local coordinate s
    virtual void basis_u_koiter_steigmann(const Vector<double>& s,
                                           Shape& psi_n,
                                           Shape& psi_i) const = 0;

    /// (pure virtual) Out-of-plane basis functions and derivs w.r.t. global
    /// coords at local coordinate s; return det(Jacobian of mapping)
    virtual double d2basis_u_eulerian_koiter_steigmann(
      const Vector<double>& s,
      Shape& psi_n,
      Shape& psi_i,
      DShape& dpsi_n_dx,
      DShape& dpsi_i_dx,
      DShape& d2psi_n_dx2,
      DShape& d2psi_i_dx2) const = 0;

    /// (pure virtual) Out-of-plane basis/test functions at local coordinate s
    virtual void basis_and_test_u_koiter_steigmann(
      const Vector<double>& s,
      Shape& psi_n,
      Shape& psi_i,
      Shape& test_n,
      Shape& test_i) const = 0;

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
      DShape& dtest_i_dx) const = 0;

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
      DShape& d2test_i_dx2) const = 0;

    // End of pure virtual interface functions
    //----------------------------------------------------------------------------

    /// Compute element residual Vector only (if flag=and/or element
    /// Jacobian matrix
    virtual void fill_in_generic_residual_contribution_koiter_steigmann(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      const unsigned& flag);



    // All member data is private
  private:
    /// Flag to use finite difference jacobian
    bool Use_finite_difference_jacobian;

    /// Flags to determine which displacement fields are damped
    bool Ui_is_damped[3];

    /// Pointer to pressure function:
    PressureVectorFctPt Pressure_fct_pt;

    /// Pointer to pressure function:
    DPressureVectorFctPt D_pressure_dr_fct_pt;

    /// Pointer to pressure function:
    DPressureVectorFctPt D_pressure_dn_fct_pt;

    /// Pointer to pressure function:
    DPressureVectorDMatrixFctPt D_pressure_d_grad_u_fct_pt;

    /// Pointer to prestrain function:
    PrestrainFctPt Prestrain_fct_pt;

    /// Pointer to stress function
    StressFctPt Stress_fct_pt;

    /// Pointer to epsilon derivative of stress function
    DStressFctPt D_stress_fct_pt;

    /// Pointer to Poisson ratio, which this element cannot modify
    const double* Nu_pt;

    /// Pointer to Mu, the damping coefficient
    const double* Mu_pt;

    /// Pointer to Non dimensional thickness parameter
    const double* Thickness_pt;

    /// Pointer to displacement scaling
    const double* Eta_u_pt;

    /// Pointer to stress scaling
    const double* Eta_sigma_pt;

    /// Pointer to error metric
    ErrorMetricFctPt Error_metric_fct_pt;

    /// Pointer to error metric when we want multiple errors
    MultipleErrorMetricFctPt Multiple_error_metric_fct_pt;

    /// unsigned that holds the number if displacments the element has
    static const unsigned Number_of_displacements;

    /// Default value of Nu incase not user set
    static const double Default_Nu_Value;

    /// Default value of Eta incase not user set
    static const double Default_Eta_Value;

    /// Default value of Mu incase not user set
    static const double Default_Mu_Value;

    /// Default value of Thickness incase not user set
    static const double Default_Thickness_Value;
  };

} // end namespace oomph
#endif

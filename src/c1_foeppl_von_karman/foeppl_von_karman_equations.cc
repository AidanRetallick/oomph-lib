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
// Non--inline functions for the foeppl_von_karman equations
#include "foeppl_von_karman_equations.h"

namespace oomph
{
  /// Default to incompressible sheet Poisson ratio
  const double FoepplVonKarmanEquations::Default_Nu_Value = 0.5;

  /// Default to no damping
  const double FoepplVonKarmanEquations::Default_Mu_Value = 0.0;

  /// Default to sheet with t/L=0.01, nu=0.5: FvK eta=12*(1.0-nu^2)*(L/t)^2
  const double FoepplVonKarmanEquations::Default_Eta_Value = 90.0e3;


  //======================================================================
  /// Fill in generic residual and jacobian contribution
  // [zdec] Local equation calls are hard coded (eg k_type+2) fix later
  //======================================================================
  void FoepplVonKarmanEquations::
    fill_in_generic_residual_contribution_foeppl_von_karman(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      const unsigned& flag)
  {
    // The indices of in-plane and out-of-plane unknowns
    const unsigned n_u_displacements = 2;
    const unsigned n_w_displacements = 1;
    const Vector<unsigned> u_indices = u_field_indices();
    const unsigned w_index = w_field_index();

    // Find the dimension of the element [zdec] will this ever not be 2?
    const unsigned dim = this->dim();
    // The number of first derivatives is the dimension of the element
    const unsigned n_deriv = dim;
    // The number of second derivatives is the triangle number of the dimension
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
    // functions however, no such basis has been used so far and so this is not
    // implemented. Comments marked with "[IN-PLANE-INTERNAL]" indicate
    // locations where changes must be made in the case that internal data is
    // used for the in-plane unknowns.
    const unsigned n_u_internal_type = nu_type_internal();
    const unsigned n_w_internal_type = nw_type_internal();

    // [IN-PLANE-INTERNAL]
    // Throw an error if the number of internal in-plane basis types is non-zero
    if (n_u_internal_type != 0)
    {
      throw OomphLibError("The number of internal basis types for u is non-zero\
 but this functionality is not yet implemented. If you want to implement this\
 look for comments containing the tag [IN-PLANE-INTERNAL].",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

    // Get the state of the flag to determine whether we are assigning a
    // displacement field to our elements rather than solving FvK equations
    const bool solve_u_exact = get_solve_u_exact();

    // Get the Poisson ratio of the plate
    const double nu = get_nu();

    // Get the plate parameters
    const double eta = get_eta();

    // Get the virtual damping coefficient
    const double mu = get_mu();

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

    // Set the value of n_intpt
    const unsigned n_intpt = this->integral_pt()->nweight();
    // Integers to store the local equation and unknown numbers
    int local_eqn = 0;
    int local_unknown = 0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the integral weight
      double w = this->integral_pt()->weight(ipt);

      // Allocate and initialise to zero
      Vector<double> interp_x(dim, 0.0);
      Vector<double> s(dim);
      s[0] = this->integral_pt()->knot(ipt, 0);
      s[1] = this->integral_pt()->knot(ipt, 1);
      interpolated_x(s, interp_x);

      // Call the derivatives of the shape and test functions for the out of
      // plane unknown
      double J =
        d2basis_and_d2test_w_eulerian_foeppl_von_karman(s,
                                                        psi_n_w,
                                                        psi_i_w,
                                                        dpsi_n_wdxi,
                                                        dpsi_i_wdxi,
                                                        d2psi_n_wdxi2,
                                                        d2psi_i_wdxi2,
                                                        test_n_w,
                                                        test_i_w,
                                                        dtest_n_wdxi,
                                                        dtest_i_wdxi,
                                                        d2test_n_wdxi2,
                                                        d2test_i_wdxi2);

      // [IN-PLANE-INTERNAL]
      // This does not include internal basis functions for u. If they are
      // needed they must be added (as above for w)
      dbasis_and_dtest_u_eulerian_foeppl_von_karman(
        s, psi_n_u, dpsi_n_udxi, test_n_u, dtest_n_udxi);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      //==========================INTERPOLATION=================================
      // Create space for in-plane interpolated unknowns
      Vector<double> interpolated_u(n_u_displacements, 0.0);
      Vector<double> interpolated_dudt(n_u_displacements, 0.0);
      DenseMatrix<double> interpolated_dudxi(n_u_displacements, n_deriv, 0.0);
      // Create space for out-of-plane interpolated unknowns
      Vector<double> interpolated_w(n_w_displacements, 0.0);
      Vector<double> interpolated_dwdt(n_w_displacements, 0.0);
      DenseMatrix<double> interpolated_dwdxi(n_w_displacements, n_deriv, 0.0);
      DenseMatrix<double> interpolated_d2wdxi2(n_w_displacements, n_2deriv, 0.0);

      //---Nodal contribution to the in-plane unknowns-------------------------
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
          for (unsigned alpha = 0; alpha < n_u_displacements; alpha++)
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


      //---Internal contribution to the in-plane unknowns----------------------
      // [IN-PLANE-INTERNAL]
      // Internal contributions to in-plane interpolation not written


      //---Nodal contribution to the out-of-plane unknowns---------------------
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

      //---Internal contribution to the out-of-plane
      // field------------------------
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

      //==================== END OF INTERPOLATION ==========================


      //----------------- ALTERNATIVE RESIDUALS IF SET ---------------------
      // If the flag Solve_u_exact is set, then we solve for a displacement
      // field assignement instead of the fvk equations.
      //   ux = ux_exact
      //   uy = uy_exact
      //    w = w_exact
      // We run this if block and then break so that we skip over the remaining
      // (normal) FvK equations residuals
      if (solve_u_exact)
      {
        // Get the exact displacement field we are trying to assign:
        //   u_exact = (ux_exact,uy_exact,w_exact)
        Vector<double> u_exact(3, 0.0);
        get_u_exact(interp_x, u_exact);

        //--- In-plane equations --------
        // Loop over the in-plane displacements alpha = 0,1
        for (unsigned alpha = 0; alpha < 2; alpha++)
        {
          // Loop over nodal equations
          for (unsigned j_node = 0; j_node < n_u_node; j_node++)
          {
            for (unsigned k_type = 0; k_type < n_u_nodal_type; k_type++)
            {
              // Get the local equation number associated with this dof, if it
              // is not pinned, add its contribution to the residuals
              int eqn_no = nodal_local_eqn(j_node, k_type);
              if (eqn_no >= 0)
              {
                // u_alpha - u_alpha_exact = 0
                residuals[eqn_no] +=
                  (interpolated_u[alpha] - u_exact[alpha]) *
                  test_n_u(j_node, (alpha * n_u_nodal_type) + k_type) * W;

                // Now the Jacobian
                if (flag)
                {
                  // Loop over the nodal dofs
                  for (unsigned jj_node = 0; jj_node < n_u_node; jj_node++)
                  {
                    for (unsigned kk_type = 0; kk_type < n_u_nodal_type;
                         kk_type++)
                    {
                      // Get the eqn number associated with this dof, if it is
                      // not pinned, add its contribution to the jacobian
                      int dof_no = nodal_local_eqn(
                        jj_node, (alpha * n_u_nodal_type) + kk_type);
                      if (dof_no >= 0)
                      {
                        jacobian(eqn_no, dof_no) +=
                          psi_n_u(jj_node, (alpha * n_u_nodal_type) + kk_type) *
                          test_n_u(j_node, (alpha * n_u_nodal_type) + k_type) *
                          W;
                      }
                    }
                  }

                  // [IN-PLANE-INTERNAL]
                  // Internal dof contributions to u would need to be added

                } // End jacobian contributions
              }
            }
          } // End nodal equations

          // [IN-PLANE-INTERNAL]
          // If we ever add internal dof contributions to in-plane, we need to
          // add their equations here

        } // End in-plane equations

        //--- Out-of-plane equations ----
        // Loop over nodal equations
        for (unsigned j_node = 0; j_node < n_w_node; j_node++)
        {
          for (unsigned k_type = 0; k_type < n_w_nodal_type; k_type++)
          {
            // Get the local equation number associated with this dof, if it
            // is not pinned, add its contribution to the residuals
            int eqn_no = nodal_local_eqn(j_node, k_type);
            if (eqn_no >= 0)
            {
              // w - w_exact = 0
              residuals[eqn_no] +=
                (interpolated_w[0] - u_exact[2]) * test_n_w(j_node, k_type) * W;

              // Now the Jacobian
              if (flag)
              {
                // Jac contribution from the nodal dofs
                for (unsigned jj_node = 0; jj_node < n_w_node; jj_node++)
                {
                  for (unsigned kk_type = 0; kk_type < n_w_nodal_type;
                       kk_type++)
                  {
                    // If this dof isn't pinned add its contribution to the
                    // jacobian
                    int dof_no = nodal_local_eqn(jj_node, kk_type + 2);
                    if (dof_no >= 0)
                    {
                      jacobian(eqn_no, dof_no) += psi_n_w(jj_node, kk_type) *
                                                  test_n_w(j_node, k_type) * W;
                    }
                  }
                } // End nodal contribution to nodal equation jac
                // Jac contribution from the internal dofs
                for (unsigned kk_type = 0; kk_type < n_w_internal_type;
                     kk_type++)
                {
                  // If this dof isn't pinned add its contribution to the
                  // jacobian
                  int dof_no = internal_local_eqn(w_index, kk_type);
                  if (dof_no >= 0)
                  {
                    jacobian(eqn_no, dof_no) +=
                      psi_i_w(kk_type) * test_n_w(j_node, k_type) * W;
                  }
                } // End internal contribution to nodal equation jac
              } // End nodal jac
            }
          }
        } // End out-of-plane nodal equations

        // Loop over out-of-plane internal equations
        for (unsigned k_type = 0; k_type < n_w_internal_type; k_type++)
        {
          // Get the local equation number associated with this dof, if it
          // is not pinned, add its contribution to the residuals
          int eqn_no = internal_local_eqn(w_index, k_type);
          if (eqn_no >= 0)
          {
            // w - w_exact = 0
            residuals[eqn_no] +=
              (interpolated_w[0] - u_exact[2]) * test_i_w(k_type) * W;

            // Now the Jacobian
            if (flag)
            {
              // Jac contribution from the nodal dofs
              for (unsigned jj_node = 0; jj_node < n_w_node; jj_node++)
              {
                for (unsigned kk_type = 0; kk_type < n_w_nodal_type; kk_type++)
                {
                  // If this dof isn't pinned add its contribution to the
                  // jacobian
                  int dof_no = nodal_local_eqn(jj_node, kk_type + 2);
                  if (dof_no >= 0)
                  {
                    jacobian(eqn_no, dof_no) +=
                      psi_n_w(jj_node, kk_type) * test_i_w(k_type) * W;
                  }
                }
              }
              // Jac contribution from the internal dofs
              for (unsigned kk_type = 0; kk_type < n_w_internal_type; kk_type++)
              {
                // If this dof isn't pinned add its contribution to the
                // jacobian
                int dof_no = internal_local_eqn(w_index, kk_type);
                if (dof_no >= 0)
                {
                  jacobian(eqn_no, dof_no) +=
                    psi_i_w(kk_type) * test_i_w(k_type) * W;
                }
              } // End internal contribution to internal equation jac
            } // End internal jac
          }
        } // End out-of-plane internal equations

        // Now break so that we don't do the rest of
        break;
      }


      //=============== FVK RESIDUALS AND JACOBIAN =========================
      // Get the swelling function
      double c_swell(0.0);
      get_swelling_foeppl_von_karman(interp_x, c_swell);

      // Get the stress
      DenseMatrix<double> sigma(2, 2, 0.0);
      get_sigma(sigma, interpolated_dudxi, interpolated_dwdxi, c_swell);

      // Get pressure function
      double pressure(0.0);
      Vector<double> traction(2, 0.0);
      get_pressure_foeppl_von_karman(ipt, interp_x, pressure);
      get_in_plane_forcing_foeppl_von_karman(ipt, interp_x, traction);


      //--------------------------------------------------------------------------
      // Outline of residual and jacobian calculation. Note, it is assumed only
      // w requires internal dofs here. Hence, only the lines with a leading 'o'
      // have been written. If u uses internal dofs the parts marked with a
      // leading x need adding.
      //
      // o w_nodal_eqn
      // o   w_nodal_eqn/u_nodal_dofs
      // x   w_nodal_eqn/u_internal_dofs
      // o   w_nodal_eqn/w_nodal_dofs
      // o   w_nodal_eqn/w_internal_dofs
      // o w_internal_eqn
      // o   w_internal_eqn/u_nodal_dofs
      // x   w_internal_eqn/u_internal_dofs
      // o   w_internal_eqn/w_nodal_dofs
      // o   w_internal_eqn/w_internal_dofs
      // o u_nodal_eqn
      // o   u_nodal_eqn/u_nodal_dofs
      // x   u_nodal_eqn/u_internal_dofs
      // o   u_nodal_eqn/w_nodal_dofs
      // o   u_nodal_eqn/w_internal_dofs
      // x u_internal_eqn
      // x   u_internal_eqn/u_nodal_dofs
      // x   u_internal_eqn/u_internal_dofs
      // x   u_internal_eqn/w_nodal_dofs
      // x   u_internal_eqn/w_internal_dofs


      // w_nodal_eqn
      //--------------------------------------------------------------------------
      // Loop over the nodal test functions
      for (unsigned j_node = 0; j_node < n_w_node; j_node++)
      {
        for (unsigned k_type = 0; k_type < n_w_nodal_type; k_type++)
        {
          // Get the local equation
          local_eqn = nodal_local_eqn(j_node, k_type + 2);
          // IF it's not a boundary condition
          if (local_eqn >= 0)
          {
            // Add virtual time derivative term for damped buckling
            residuals[local_eqn] +=
              mu * interpolated_dwdt[0] * test_n_w(j_node, k_type) * W;
            // Add body force/pressure term here
            residuals[local_eqn] += -pressure * test_n_w(j_node, k_type) * W;

            for (unsigned alpha = 0; alpha < 2; ++alpha)
            {
              for (unsigned beta = 0; beta < 2; ++beta)
              {
                // w_{,\alpha\beta} \delta \kappa_\alpha\beta
                residuals[local_eqn] +=
                  (1.0 - nu) * interpolated_d2wdxi2(0, alpha + beta) *
                  d2test_n_wdxi2(j_node, k_type, alpha + beta) * W;
                // w_{,\alpha\alpha} \delta \kappa_\beta\beta
                residuals[local_eqn] +=
                  (nu)*interpolated_d2wdxi2(0, beta + beta) *
                  d2test_n_wdxi2(j_node, k_type, alpha + alpha) * W;
                // sigma_{\alpha\beta} w_{,\alpha} \delta w_{,\beta}
                residuals[local_eqn] += eta * sigma(alpha, beta) *
                                        interpolated_dwdxi(0, alpha) *
                                        dtest_n_wdxi(j_node, k_type, beta) * W;
              }
            }

            // Calculate the jacobian
            //--------------------------------------------------------------------
            if (flag)
            {
              // d(w_nodal_equation)/d(u_nodal_dofs)
              //------------------------------------------------------------------
              // Loop over the in-plane nodes
              for (unsigned jj_node = 0; jj_node < n_u_node; jj_node++)
              {
                // Loop over the in-plane types
                for (unsigned kk_type = 0; kk_type < n_u_nodal_type; kk_type++)
                {
                  // Loop over in-plane displacements
                  for (unsigned gamma = 0; gamma < 2; ++gamma)
                  {
                    local_unknown = this->nodal_local_eqn(jj_node, gamma);
                    // If at a non-zero degree of freedom add in the entry
                    if (local_unknown >= 0)
                    {
                      // Loop over in-plane displacements
                      for (unsigned alpha = 0; alpha < 2; ++alpha)
                      {
                        // Loop over dimensions
                        for (unsigned beta = 0; beta < 2; ++beta)
                        {
                          // The Nonlinear Terms
                          if (alpha == beta)
                          {
                            jacobian(local_eqn, local_unknown) +=
                              eta * nu * dpsi_n_udxi(jj_node, kk_type, gamma) *
                              interpolated_dwdxi(0, alpha) *
                              dtest_n_wdxi(j_node, k_type, beta) /
                              (1.0 - nu * nu) * W;
                          }
                          if (alpha == gamma)
                          {
                            jacobian(local_eqn, local_unknown) +=
                              eta * (1.0 - nu) / 2.0 *
                              dpsi_n_udxi(jj_node, kk_type, beta) *
                              interpolated_dwdxi(0, alpha) *
                              dtest_n_wdxi(j_node, k_type, beta) /
                              (1.0 - nu * nu) * W;
                          }
                          if (beta == gamma)
                          {
                            jacobian(local_eqn, local_unknown) +=
                              eta * (1.0 - nu) / 2.0 *
                              dpsi_n_udxi(jj_node, kk_type, alpha) *
                              interpolated_dwdxi(0, alpha) *
                              dtest_n_wdxi(j_node, k_type, beta) /
                              (1.0 - nu * nu) * W;
                          }
                        } // End loop over dimensions --beta
                      } // End loop over in plane displacements -- alpha
                    } // End if dof not pinned -- local_unknown >= 0
                  } // End loop over displacements -- gamma
                } // End loop over in-plane types -- kk_type
              } // End loop over in-plane nodes -- jj_node


              // d(w_nodal_eqn)/d(u_internal_dofs)
              //------------------------------------------------------------------
              // [IN-PLANE-INTERNAL] Not implemented, if needed, write it


              // d(w_nodal_eqn)/d(w_nodal_dofs)
              //------------------------------------------------------------------
              // Loop over the w basis nodes
              for (unsigned j_node2 = 0; j_node2 < n_w_node; j_node2++)
              {
                // Loop over w basis function types
                for (unsigned k_type2 = 0; k_type2 < n_w_nodal_type; k_type2++)
                {
                  local_unknown = this->nodal_local_eqn(j_node2, k_type2 + 2);
                  // If at a non-zero degree of freedom add in the entry
                  if (local_unknown >= 0)
                  {
                    // Contribution from damping
                    jacobian(local_eqn, local_unknown) +=
                      mu * psi_n_w(j_node2, k_type2) *
                      node_pt(j_node2)->time_stepper_pt()->weight(1, 0) *
                      test_n_w(j_node, k_type) * W;
                    // Loop over dimensions
                    for (unsigned alpha = 0; alpha < 2; ++alpha)
                    {
                      for (unsigned beta = 0; beta < 2; ++beta)
                      {
                        // Linear Terms
                        // w_{,\alpha\beta} \delta \kappa_\alpha\beta
                        jacobian(local_eqn, local_unknown) +=
                          (1.0 - nu) *
                          d2psi_n_wdxi2(j_node2, k_type2, alpha + beta) *
                          d2test_n_wdxi2(j_node, k_type, alpha + beta) * W;
                        // w_{,\alpha\alpha} \delta \kappa_\beta\beta
                        jacobian(local_eqn, local_unknown) +=
                          nu * d2psi_n_wdxi2(j_node2, k_type2, beta + beta) *
                          d2test_n_wdxi2(j_node, k_type, alpha + alpha) * W;
                        // Nonlinear terms
                        jacobian(local_eqn, local_unknown) +=
                          eta * sigma(alpha, beta) *
                          dpsi_n_wdxi(j_node2, k_type2, alpha) *
                          dtest_n_wdxi(j_node, k_type, beta) * W;
                        // Nonlinear terms
                        jacobian(local_eqn, local_unknown) +=
                          eta * nu / (1.0 - nu * nu) *
                          interpolated_dwdxi(0, alpha) *
                          dpsi_n_wdxi(j_node2, k_type2, alpha) *
                          interpolated_dwdxi(0, beta) *
                          dtest_n_wdxi(j_node, k_type, beta) * W;
                        jacobian(local_eqn, local_unknown) +=
                          eta / 2.0 * (1.0 - nu) / (1.0 - nu * nu) *
                          interpolated_dwdxi(0, alpha) *
                          dpsi_n_wdxi(j_node2, k_type2, beta) *
                          interpolated_dwdxi(0, alpha) *
                          dtest_n_wdxi(j_node, k_type, beta) * W;
                        jacobian(local_eqn, local_unknown) +=
                          eta / 2.0 * (1.0 - nu) / (1.0 - nu * nu) *
                          interpolated_dwdxi(0, beta) *
                          dpsi_n_wdxi(j_node2, k_type2, alpha) *
                          interpolated_dwdxi(0, alpha) *
                          dtest_n_wdxi(j_node, k_type, beta) * W;
                      } // End loop over beta
                    } // End loop over alpha
                  } // End if dof not pinned -- local_unknown>=0
                } // End loop over the nodal basis functions -- k_type2
              } // End loop over the nodes used by w -- j_node2


              // d(w_nodal_eqn)/d(w_internal_dofs)
              //------------------------------------------------------------------
              // Loop over the internal test functions for w
              for (unsigned k_type2 = 0; k_type2 < n_w_internal_type; k_type2++)
              {
                local_unknown = internal_local_eqn(w_index, k_type2);
                // If at a non-zero degree of freedom add in the entry
                if (local_unknown >= 0)
                {
                  // Contribution from buckle stabilising drag
                  jacobian(local_eqn, local_unknown) +=
                    psi_i_w(k_type2) *
                    w_internal_data_pt()->time_stepper_pt()->weight(1, 0) * mu *
                    test_n_w(j_node, k_type) * W;
                  // Loop over dimensions
                  for (unsigned alpha = 0; alpha < 2; ++alpha)
                  {
                    for (unsigned beta = 0; beta < 2; ++beta)
                    {
                      // The Linear Terms
                      // w_{,\alpha\beta} \delta \kappa_\alpha\beta
                      jacobian(local_eqn, local_unknown) +=
                        (1.0 - nu) * d2psi_i_wdxi2(k_type2, alpha + beta) *
                        d2test_n_wdxi2(j_node, k_type, alpha + beta) * W;
                      // w_{,\alpha\alpha} \delta \kappa_\beta\beta
                      jacobian(local_eqn, local_unknown) +=
                        nu * d2psi_i_wdxi2(k_type2, beta + beta) *
                        d2test_n_wdxi2(j_node, k_type, alpha + alpha) * W;
                      // Nonlinear terms
                      jacobian(local_eqn, local_unknown) +=
                        eta * sigma(alpha, beta) * dpsi_i_wdxi(k_type2, alpha) *
                        dtest_n_wdxi(j_node, k_type, beta) * W;
                      // Nonlinear terms
                      jacobian(local_eqn, local_unknown) +=
                        eta * nu / (1.0 - nu * nu) *
                        interpolated_dwdxi(0, alpha) *
                        dpsi_i_wdxi(k_type2, alpha) *
                        interpolated_dwdxi(0, beta) *
                        dtest_n_wdxi(j_node, k_type, beta) * W;
                      jacobian(local_eqn, local_unknown) +=
                        eta / 2.0 * (1.0 - nu) / (1.0 - nu * nu) *
                        interpolated_dwdxi(0, alpha) *
                        dpsi_i_wdxi(k_type2, beta) *
                        interpolated_dwdxi(0, alpha) *
                        dtest_n_wdxi(j_node, k_type, beta) * W;
                      jacobian(local_eqn, local_unknown) +=
                        eta / 2.0 * (1.0 - nu) / (1.0 - nu * nu) *
                        interpolated_dwdxi(0, beta) *
                        dpsi_i_wdxi(k_type2, alpha) *
                        interpolated_dwdxi(0, alpha) *
                        dtest_n_wdxi(j_node, k_type, beta) * W;
                    } // End of loop over beta
                  } // End of loop over alpha
                } // End of if local dof is not pinned -- local_unknown>0
              } // End of loop over the internal types -- k_type2

            } // End of if making the jacobian -- flag
          } // End of if local eqn not pinned -- local_eqn>0
        } // End of loop over w nodal type -- k_type
      } // End of loop over w nodes -- j_node


      // w_internal_eqn
      //--------------------------------------------------------------------------
      // Loop over the internal test functions
      for (unsigned k_type = 0; k_type < n_w_internal_type; k_type++)
      {
        // Get the local equation
        local_eqn = internal_local_eqn(w_index, k_type);
        // IF it's not a boundary condition
        if (local_eqn >= 0)
        {
          // Add virtual time derivative term for damped buckling
          residuals[local_eqn] +=
            mu * interpolated_dwdt[0] * test_i_w(k_type) * W;
          // Add body force/pressure term here
          residuals[local_eqn] += -pressure * test_i_w(k_type) * W;

          for (unsigned alpha = 0; alpha < 2; ++alpha)
          {
            for (unsigned beta = 0; beta < 2; ++beta)
            {
              // w_{,\alpha\beta} \delta \kappa_\alpha\beta
              residuals[local_eqn] += (1.0 - nu) *
                                      interpolated_d2wdxi2(0, alpha + beta) *
                                      d2test_i_wdxi2(k_type, alpha + beta) * W;
              // w_{,\alpha\alpha} \delta \kappa_\beta\beta
              residuals[local_eqn] +=
                (nu)*interpolated_d2wdxi2(0, beta + beta) *
                d2test_i_wdxi2(k_type, alpha + alpha) * W;
              // sigma_{\alpha\beta} w_\alpha \delta w_{,\beta}
              residuals[local_eqn] += eta * sigma(alpha, beta) *
                                      interpolated_dwdxi(0, alpha) *
                                      dtest_i_wdxi(k_type, beta) * W;
            }
          }

          // Calculate the jacobian
          //-----------------------
          if (flag)
          {
            // d(w_internal_eqn)/d(u_nodal_dofs)
            //------------------------------------------------------------------
            // Loop over the in--plane unknowns again
            for (unsigned jj_node = 0; jj_node < n_u_node; jj_node++)
            {
              for (unsigned kk_type = 0; kk_type < n_u_nodal_type; kk_type++)
              {
                // Loop over displacements
                for (unsigned gamma = 0; gamma < 2; ++gamma)
                {
                  local_unknown = this->nodal_local_eqn(jj_node, gamma);
                  // If at a non-zero degree of freedom add in the entry
                  if (local_unknown >= 0)
                  {
                    // Loop over displacements
                    for (unsigned alpha = 0; alpha < 2; ++alpha)
                    {
                      // Loop over dimensions
                      for (unsigned beta = 0; beta < 2; ++beta)
                      {
                        // The Nonlinear Terms
                        if (alpha == beta)
                        {
                          jacobian(local_eqn, local_unknown) +=
                            eta * nu * dpsi_n_udxi(jj_node, kk_type, gamma) *
                            interpolated_dwdxi(0, alpha) *
                            dtest_i_wdxi(k_type, beta) / (1.0 - nu * nu) * W;
                        }
                        if (alpha == gamma)
                        {
                          jacobian(local_eqn, local_unknown) +=
                            eta * (1.0 - nu) / 2.0 *
                            dpsi_n_udxi(jj_node, kk_type, beta) *
                            interpolated_dwdxi(0, alpha) *
                            dtest_i_wdxi(k_type, beta) / (1.0 - nu * nu) * W;
                        }
                        if (beta == gamma)
                        {
                          jacobian(local_eqn, local_unknown) +=
                            eta * (1.0 - nu) / 2.0 *
                            dpsi_n_udxi(jj_node, kk_type, alpha) *
                            interpolated_dwdxi(0, alpha) *
                            dtest_i_wdxi(k_type, beta) / (1.0 - nu * nu) * W;
                        }
                      } // End of loop over beta
                    } // End of loop over alpha
                  } // End of if local dof not pinned -- local_unknown>0
                } // End of loop over in-plane fields -- gamma
              } // End of loop over in-plane nodal types -- kk_type
            } // End of loop over in-plane nodes -- jj_node


            // d(w_internal_eqn)/d(u_internal_dofs)
            //------------------------------------------------------------------
            // [IN-PLANE-INTERNAL] Not implemented, if needed, write it


            // d(w_internal_eqn)/d(w_nodal_dofs)
            //------------------------------------------------------------------
            // Loop over the test functions again
            for (unsigned j_node2 = 0; j_node2 < n_w_node; j_node2++)
            {
              // Loop over position dofs
              for (unsigned k_type2 = 0; k_type2 < n_w_nodal_type; k_type2++)
              {
                local_unknown = this->nodal_local_eqn(
                  j_node2, k_type2 + 2); // [zdec] manually indexed
                // If at a non-zero degree of freedom add in the entry
                if (local_unknown >= 0)
                {
                  // Contribution from buckle stabilising drag
                  jacobian(local_eqn, local_unknown) +=
                    mu * psi_n_w(j_node2, k_type2) *
                    node_pt(j_node2)->time_stepper_pt()->weight(1, 0) *
                    test_i_w(k_type) * W;
                  // Loop over dimensions
                  for (unsigned alpha = 0; alpha < 2; ++alpha)
                  {
                    for (unsigned beta = 0; beta < 2; ++beta)
                    {
                      // The Linear Terms
                      // w_{,\alpha\beta} \delta \kappa_\alpha\beta
                      jacobian(local_eqn, local_unknown) +=
                        (1.0 - nu) *
                        d2psi_n_wdxi2(j_node2, k_type2, alpha + beta) *
                        d2test_i_wdxi2(k_type, alpha + beta) * W;
                      // w_{,\alpha\alpha} \delta \kappa_\beta\beta
                      jacobian(local_eqn, local_unknown) +=
                        nu * d2psi_n_wdxi2(j_node2, k_type2, beta + beta) *
                        d2test_i_wdxi2(k_type, alpha + alpha) * W;
                      // Nonlinear terms
                      jacobian(local_eqn, local_unknown) +=
                        eta * sigma(alpha, beta) *
                        dpsi_n_wdxi(j_node2, k_type2, alpha) *
                        dtest_i_wdxi(k_type, beta) * W;
                      // Nonlinear terms
                      jacobian(local_eqn, local_unknown) +=
                        eta * nu / (1.0 - nu * nu) *
                        interpolated_dwdxi(0, alpha) *
                        dpsi_n_wdxi(j_node2, k_type2, alpha) *
                        interpolated_dwdxi(0, beta) *
                        dtest_i_wdxi(k_type, beta) * W;
                      jacobian(local_eqn, local_unknown) +=
                        eta / 2.0 * (1.0 - nu) / (1.0 - nu * nu) *
                        interpolated_dwdxi(0, alpha) *
                        dpsi_n_wdxi(j_node2, k_type2, beta) *
                        interpolated_dwdxi(0, alpha) *
                        dtest_i_wdxi(k_type, beta) * W;
                      jacobian(local_eqn, local_unknown) +=
                        eta / 2.0 * (1.0 - nu) / (1.0 - nu * nu) *
                        interpolated_dwdxi(0, beta) *
                        dpsi_n_wdxi(j_node2, k_type2, alpha) *
                        interpolated_dwdxi(0, alpha) *
                        dtest_i_wdxi(k_type, beta) * W;
                    } // End loop over beta
                  } // End loop over alpha
                } // End of if dof is pinned -- local_unknown>0
              } // End of loop over w nodal types -- k_type2
            } // End of loop over w nodes -- j_node2

            // d(w_internal_eqn)/d(w_internal_dofs)
            // Loop over the test functions again
            for (unsigned k_type2 = 0; k_type2 < n_w_internal_type; k_type2++)
            {
              local_unknown =
                internal_local_eqn(w_index, k_type2); // [zdec] manually indexed
              // If at a non-zero degree of freedom add in the entry
              if (local_unknown >= 0)
              {
                // Contribution from buckle stabilising drag
                jacobian(local_eqn, local_unknown) +=
                  mu * psi_i_w(k_type2) *
                  w_internal_data_pt()->time_stepper_pt()->weight(1, 0) *
                  test_i_w(k_type) * W;
                // Loop over dimensions
                for (unsigned alpha = 0; alpha < 2; alpha++)
                {
                  for (unsigned beta = 0; beta < 2; beta++)
                  {
                    // The Linear Terms
                    // w_{,\alpha\beta} \delta \kappa_\alpha\beta
                    jacobian(local_eqn, local_unknown) +=
                      (1.0 - nu) * d2psi_i_wdxi2(k_type2, alpha + beta) *
                      d2test_i_wdxi2(k_type, alpha + beta) * W;
                    // w_{,\alpha\alpha} \delta \kappa_\beta\beta
                    jacobian(local_eqn, local_unknown) +=
                      nu * d2psi_i_wdxi2(k_type2, beta + beta) *
                      d2test_i_wdxi2(k_type, alpha + alpha) * W;
                    // Nonlinear terms
                    jacobian(local_eqn, local_unknown) +=
                      eta * sigma(alpha, beta) * dpsi_i_wdxi(k_type2, alpha) *
                      dtest_i_wdxi(k_type, beta) * W;
                    // Nonlinear terms
                    jacobian(local_eqn, local_unknown) +=
                      eta * nu / (1.0 - nu * nu) *
                      interpolated_dwdxi(0, alpha) *
                      dpsi_i_wdxi(k_type2, alpha) *
                      interpolated_dwdxi(0, beta) * dtest_i_wdxi(k_type, beta) *
                      W;
                    jacobian(local_eqn, local_unknown) +=
                      eta / 2.0 * (1.0 - nu) / (1.0 - nu * nu) *
                      interpolated_dwdxi(0, alpha) *
                      dpsi_i_wdxi(k_type2, beta) *
                      interpolated_dwdxi(0, alpha) *
                      dtest_i_wdxi(k_type, beta) * W;
                    jacobian(local_eqn, local_unknown) +=
                      eta / 2.0 * (1.0 - nu) / (1.0 - nu * nu) *
                      interpolated_dwdxi(0, beta) *
                      dpsi_i_wdxi(k_type2, alpha) *
                      interpolated_dwdxi(0, alpha) *
                      dtest_i_wdxi(k_type, beta) * W;
                  } // End of loop over beta
                } // End of loop over alpha
              } // End of if local dof is not pinned -- local_unknown>0
            } // End of loop over internal types -- k_type2

          } // End of if jacobian -- flag
        } // End of if not pinned -- local_eqn>0
      } // End loop over w internal types


      // u_nodal_eqn
      //--------------------------------------------------------------------------
      for (unsigned j_node = 0; j_node < n_u_node; j_node++)
      {
        // Loop over nodal types
        for (unsigned k_type = 0; k_type < n_u_nodal_type; k_type++)
        {
          // Now loop over displacement equations
          for (unsigned alpha = 0; alpha < 2; alpha++)
          {
            // Get the local equation
            local_eqn = this->nodal_local_eqn(j_node, alpha);
            // IF it's not a boundary condition
            if (local_eqn >= 0)
            {
              // Now add on the stress terms
              for (unsigned beta = 0; beta < 2; beta++)
              {
                // Variations in u_alpha
                residuals[local_eqn] +=
                  sigma(alpha, beta) * dtest_n_udxi(j_node, k_type, beta) * W;
              }
              // Add the forcing term
              residuals[local_eqn] +=
                traction[alpha] * test_n_u(j_node, k_type) * W;

              // Now loop over Jacobian
              if (flag)
              {
                // d(u_nodal_eqn)/d(u_nodal_dofs)
                //----------------------------------------------------------------
                // Loop over in-plane nodes again
                for (unsigned jj_node = 0; jj_node < n_u_node; jj_node++)
                {
                  for (unsigned kk_type = 0; kk_type < n_u_nodal_type;
                       kk_type++)
                  {
                    // Now loop over displacement unknowns
                    for (unsigned gamma = 0; gamma < 2; ++gamma)
                    {
                      // Get the local equation
                      local_unknown = this->nodal_local_eqn(jj_node, gamma);
                      // If the unknown is a degree of freedom
                      if (local_unknown >= 0)
                      {
                        // Now add on the stress terms
                        for (unsigned beta = 0; beta < 2; ++beta)
                        {
                          // The Linear Terms
                          if (alpha == beta)
                          {
                            jacobian(local_eqn, local_unknown) +=
                              nu * dpsi_n_udxi(jj_node, kk_type, gamma) *
                              dtest_n_udxi(j_node, k_type, beta) /
                              (1.0 - nu * nu) * W;
                          }
                          if (alpha == gamma)
                          {
                            jacobian(local_eqn, local_unknown) +=
                              (1.0 - nu) / 2.0 *
                              dpsi_n_udxi(jj_node, kk_type, beta) *
                              dtest_n_udxi(j_node, k_type, beta) /
                              (1.0 - nu * nu) * W;
                          }
                          if (beta == gamma)
                          {
                            jacobian(local_eqn, local_unknown) +=
                              (1.0 - nu) / 2.0 *
                              dpsi_n_udxi(jj_node, kk_type, alpha) *
                              dtest_n_udxi(j_node, k_type, beta) /
                              (1.0 - nu * nu) * W;
                          }
                        } // End loop beta
                      } // End if local unknown not pinned -- local_unknown>=0
                    } // End loop in-plane unknowns -- gamma
                  } // End loop u nodal types -- kk_type
                } // End loop u nodes -- jj_node


                // d(u_nodal_eqn)/d(u_internal_dofs)
                //------------------------------------------------------------------
                // [IN-PLANE-INTERNAL] Not implemented, if needed, write it


                // d(u_nodal_eqn)/d(w_nodal_dofs)
                // Loop over the w nodal dofs
                for (unsigned j_node2 = 0; j_node2 < n_w_node; j_node2++)
                {
                  // Loop over position dofs
                  for (unsigned k_type2 = 0; k_type2 < n_w_nodal_type;
                       k_type2++)
                  {
                    // If at a non-zero degree of freedom add in the entry
                    // Get the local equation
                    local_unknown = this->nodal_local_eqn(
                      j_node2, k_type2 + 2); // [zdec] manually indexed
                    // If the unknown is a degree of freedom
                    if (local_unknown >= 0)
                    {
                      // Now add on the stress terms
                      for (unsigned beta = 0; beta < 2; ++beta)
                      {
                        // The NonLinear Terms
                        jacobian(local_eqn, local_unknown) +=
                          nu * dpsi_n_wdxi(j_node2, k_type2, beta) *
                          interpolated_dwdxi(0, beta) *
                          dtest_n_udxi(j_node, k_type, alpha) /
                          (1.0 - nu * nu) * W;
                        jacobian(local_eqn, local_unknown) +=
                          (1.0 - nu) / 2.0 *
                          dpsi_n_wdxi(j_node2, k_type2, alpha) *
                          interpolated_dwdxi(0, beta) *
                          dtest_n_udxi(j_node, k_type, beta) / (1.0 - nu * nu) *
                          W;
                        jacobian(local_eqn, local_unknown) +=
                          (1.0 - nu) / 2.0 *
                          dpsi_n_wdxi(j_node2, k_type2, beta) *
                          interpolated_dwdxi(0, alpha) *
                          dtest_n_udxi(j_node, k_type, beta) / (1 - nu * nu) *
                          W;
                      } // End loop beta
                    } // End if local dof is not pinned -- local_unknown>=0
                  } // End loop w nodal type -- kk_type
                } // End loop w node -- jj_node

                // d(u_nodal_eqn)/d(w_internal_dofs)
                //----------------------------------------------------------------
                // Loop over the internal types for w
                for (unsigned k_type2 = 0; k_type2 < n_w_internal_type;
                     k_type2++)
                {
                  // Get the local equation
                  local_unknown = internal_local_eqn(
                    w_index, k_type2); // [zdec] manually indexed
                  // If at a non-zero degree of freedom add in the entry
                  if (local_unknown >= 0)
                  {
                    // Now add on the stress terms
                    for (unsigned beta = 0; beta < 2; ++beta)
                    {
                      // The NonLinear Terms
                      jacobian(local_eqn, local_unknown) +=
                        nu * dpsi_i_wdxi(k_type2, beta) *
                        interpolated_dwdxi(0, beta) *
                        dtest_n_udxi(j_node, k_type, alpha) / (1.0 - nu * nu) *
                        W;
                      jacobian(local_eqn, local_unknown) +=
                        (1.0 - nu) / 2.0 * dpsi_i_wdxi(k_type2, alpha) *
                        interpolated_dwdxi(0, beta) *
                        dtest_n_udxi(j_node, k_type, beta) / (1.0 - nu * nu) *
                        W;
                      jacobian(local_eqn, local_unknown) +=
                        (1.0 - nu) / 2.0 * dpsi_i_wdxi(k_type2, beta) *
                        interpolated_dwdxi(0, alpha) *
                        dtest_n_udxi(j_node, k_type, beta) / (1.0 - nu * nu) *
                        W;
                    } // End loop beta
                  } // End if local dof is unpinned -- local_unknown>=0
                } // End loop w internal type -- k_type2

              } // End if make jacobian -- flag
            } // End if local eqn not pinned -- local_equation>=0
          } // End loop over in-plane fields -- alpha
        } // End of loop over u nodal types -- k_type
      } // End loop over u nodes -- j_node


      // u_internal_eqn
      //------------------------------------------------------------------
      // [IN-PLANE-INTERNAL] Not implemented, if needed, write it
      {
        // d(u_internal_eqn)/d(u_internal_dofs)
        //------------------------------------------------------------------
        // [IN-PLANE-INTERNAL] Not implemented, if needed, write it

        // d(u_internal_eqn)/d(u_internal_dofs)
        //------------------------------------------------------------------
        // [IN-PLANE-INTERNAL] Not implemented, if needed, write it

        // d(u_internal_eqn)/d(u_internal_dofs)
        //------------------------------------------------------------------
        // [IN-PLANE-INTERNAL] Not implemented, if needed, write it

        // d(u_internal_eqn)/d(u_internal_dofs)
        //------------------------------------------------------------------
        // [IN-PLANE-INTERNAL] Not implemented, if needed, write it
      }


    } // End of loop over integration points
  } // End of fill_in_generic_residual_contribution_foeppl_von_karman


  //======================================================================
  /// Self-test:  Return 0 for OK
  //======================================================================
  unsigned FoepplVonKarmanEquations::self_test()
  {
    bool passed = true;

    // Check lower-level stuff
    if (FiniteElement::self_test() != 0)
    {
      passed = false;
    }

    // Return verdict
    if (passed)
    {
      return 0;
    }
    else
    {
      return 1;
    }
  }


  //======================================================================
  /// Full output function (large - 26 values per plot point):
  ///
  ///   - : x
  ///   - : y
  ///   0 : ux           13: strain_xx
  ///   1 : uy	         14: strain_xy
  ///   2 : w	         15: strain_yy
  ///   3 : du_x/dx      16: stress_xx
  ///   4 : du_x/dy      17: stress_xy
  ///   5 : du_y/dx      18: stress_yy
  ///   6 : du_y/dy      19: principal_stress_val_1
  ///   7 : dw/dx        20: principal_stress_val_2
  ///   8 : dw/dy        21: principal_stress_vec_1x
  ///   9 : d2w/dx2      22: principal_stress_vec_1y
  ///   10: d2w/dxd      23: principal_stress_vec_2x
  ///   11: d2w/dy2      24: principal_stress_vec_2y
  ///
  /// nplot points in each coordinate direction
  //======================================================================
  void FoepplVonKarmanEquations::full_output(std::ostream& outfile,
                                             const unsigned& nplot)
  {
    unsigned dim = this->dim();

    // Vector of local coordinates
    Vector<double> s(dim), x(dim);

    // Tecplot header info
    outfile << this->tecplot_zone_string(nplot);

    // Storage for variables
    double c_swell(0.0);
    Vector<double> interpolated_vals(12, 0.0);
    DenseMatrix<double> interpolated_dwdxj(1, dim, 0.0);
    DenseMatrix<double> interpolated_duidxj(2, dim, 0.0);
    DenseMatrix<double> epsilon(2, 2, 0.0);
    DenseMatrix<double> sigma(2, 2, 0.0);
    DenseMatrix<double> sigma_eigenvecs(2, 2, 0.0);
    Vector<double> sigma_eigenvals(2, 0.0);
    unsigned num_plot_points = this->nplot_points(nplot);
    Vector<double> r(3);

    // Loop over plot points
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local and global coordinates of plot point
      this->get_s_plot(iplot, nplot, s);
      interpolated_x(s, x);

      // Get interpolated unknowns
      interpolated_vals = interpolated_fvk_disp_and_deriv(s);

      // Get degree of swelling for use in the strain tensor
      this->get_swelling_foeppl_von_karman(x, c_swell);

      // TODO: make the output customisable with flags e.g. output_axial_strain,
      //       output_principal_stress
      // TODO: make the indexing below more general and not hard coded.
      // Copy gradients from u into interpolated gradient matrices...
      interpolated_dwdxj(0, 0) = interpolated_vals[7]; // dw/dx
      interpolated_dwdxj(0, 1) = interpolated_vals[8]; // dw/dy
      interpolated_duidxj(0, 0) = interpolated_vals[3]; // du_x/dx
      interpolated_duidxj(0, 1) = interpolated_vals[4]; // du_x/dy
      interpolated_duidxj(1, 0) = interpolated_vals[5]; // du_y/dx
      interpolated_duidxj(1, 1) = interpolated_vals[6]; // du_y/dy
      // ...which are used to retrieve the strain tensor epsilon.
      get_epsilon(epsilon, interpolated_duidxj, interpolated_dwdxj, c_swell);

      // Use epsilon to find the stress sigma.
      get_sigma_from_epsilon(sigma, epsilon);

      // Get the principal values and the corresponding directions of stress.
      get_principal_stresses(sigma, sigma_eigenvals, sigma_eigenvecs);


      for (unsigned i = 0; i < this->dim(); i++)
      {
        outfile << x[i] << " ";
      }

      // Loop for variables
      for (Vector<double>::iterator it = interpolated_vals.begin();
           it != interpolated_vals.end();
           ++it)
      {
        outfile << *it << " ";
      }

      // Output axial strains
      outfile << epsilon(0, 0) << " " << epsilon(0, 1) << " " << epsilon(1, 1)
              << " ";

      // Output axial stress
      outfile << sigma(0, 0) << " " << sigma(0, 1) << " " << sigma(1, 1) << " ";

      // Output principal stresses
      outfile << sigma_eigenvals[0] << " " << sigma_eigenvals[1] << " ";

      // Output principal stress directions
      outfile << sigma_eigenvecs(0, 0) << " " << sigma_eigenvecs(1, 0) << " "
              << sigma_eigenvecs(0, 1) << " " << sigma_eigenvecs(1, 1) << " ";

      // End output line
      outfile << "\n";
    }
    // Write tecplot footer (e.g. FE connectivity lists)
    this->write_tecplot_zone_footer(outfile, nplot);
  }

  //======================================================================
  /// Default output function (small - 5 values per plot point):
  ///
  ///   - : x
  ///   - : y
  ///   0 : ux
  ///   1 : uy
  ///   2 : w
  ///
  /// nplot points in each coordinate direction
  //======================================================================
  void FoepplVonKarmanEquations::output(std::ostream& outfile,
                                        const unsigned& nplot)
  {
    // Tecplot header info
    outfile << this->tecplot_zone_string(nplot);

    // Number of coordinates
    unsigned dim = this->dim();
    // Number of unknowns
    unsigned n_unknown = 3;

    // Loop over plot points
    unsigned num_plot_points = this->nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      Vector<double> s(dim);
      this->get_s_plot(iplot, nplot, s);

      // Output the interpolated global coords
      Vector<double> x(dim, 0.0);
      interpolated_x(s, x);
      for (unsigned i = 0; i < dim; i++)
      {
        outfile << x[i] << " ";
      }

      // Output the interpolated unknowns
      Vector<double> interpolated_vals(n_unknown, 0.0);
      interpolated_vals = interpolated_fvk_disp_and_deriv(s);
      for (Vector<double>::iterator it = interpolated_vals.begin();
           it != interpolated_vals.end();
           ++it)
      {
        outfile << *it << " ";
      }

      // End plot point line
      outfile << "\n";

    }

    // Write tecplot footer (e.g. FE connectivity lists)
    this->write_tecplot_zone_footer(outfile, nplot);
  }


  //======================================================================
  /// C-style output function:
  ///
  ///   -: x
  ///   -: y
  ///   0: u_x
  ///   1: u_y
  ///   2: w
  ///
  /// nplot points in each coordinate direction
  //======================================================================
  void FoepplVonKarmanEquations::output(FILE* file_pt, const unsigned& nplot)
  {
    // Tecplot header info
    fprintf(file_pt, "%s", this->tecplot_zone_string(nplot).c_str());

    // Number of coordinates
    unsigned dim = this->dim();
    // Number of unknowns
    unsigned n_unknown = 3;

    // Loop over plot points
    unsigned num_plot_points = this->nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      Vector<double> s(dim);
      this->get_s_plot(iplot, nplot, s);

      // Output the interpolated global coords
      Vector<double> x(dim, 0.0);
      interpolated_x(s, x);
      for (unsigned i = 0; i < dim; i++)
      {
        fprintf(file_pt, "%g ", x[i]);
      }

      // Output the interpolated unknowns
      Vector<double> interpolated_vals(n_unknown, 0.0);
      interpolated_vals = interpolated_fvk_disp_and_deriv(s);
      for (unsigned i = 0; i < n_unknown; i++)
      {
        fprintf(file_pt, "%g \n", interpolated_vals[i]);
      }
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    this->write_tecplot_zone_footer(file_pt, nplot);
  }


  //======================================================================
  /// Output exact solution
  ///
  /// Solution is provided via function pointer.
  /// Plot at a given number of plot points.
  ///
  ///   x,y,u_exact    or    x,y,z,u_exact
  //======================================================================
  void FoepplVonKarmanEquations::output_fct(
    std::ostream& outfile,
    const unsigned& nplot,
    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
  {
    // Number of coordinates
    unsigned dim = this->dim();

    // Vector of local coordinates
    Vector<double> s(dim);

    // Vector for coordintes
    Vector<double> x(dim);

    // Tecplot header info
    outfile << this->tecplot_zone_string(nplot);

    // Exact solution Vector (here a scalar)
    //
    // [zdec] Why was this->required_nvalue(0) called here? Caused range
    // checking errors due to required_nvalue(0) being 8 but exact_soln_pt
    // resizing to a vector of length one. Replaced with 1 for now as the above
    // comment states we already know this function to be scalar.
    //
    // Vector<double> exact_soln(this->required_nvalue(0),0.0);
    unsigned nvalue = 1;
    Vector<double> exact_soln(nvalue, 0.0);

    // Loop over plot points
    unsigned num_plot_points = this->nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      this->get_s_plot(iplot, nplot, s);

      // Get x position as Vector
      interpolated_x(s, x);

      // Get exact solution at this point
      (*exact_soln_pt)(x, exact_soln);

      // Output x,y,...,u_exact
      for (unsigned i = 0; i < dim; i++)
      {
        outfile << x[i] << " ";
      }
      // Loop over variables
      for (unsigned j = 0; j < nvalue; j++)
      {
        outfile << exact_soln[j] << " ";
      }
      outfile << std::endl;
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    this->write_tecplot_zone_footer(outfile, nplot);
  } // End of output_fct(...)


  //======================================================================
  /// Validate against exact solution
  ///
  /// Solution is provided via function pointer.
  /// Plot error at a given number of plot points.
  ///
  /// [zdec] Is this function necessary?
  //======================================================================
  void FoepplVonKarmanEquations::compute_error_in_deflection(
    std::ostream& outfile,
    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
    double& error,
    double& norm)
  {
    // Initialise
    error = 0.0;
    norm = 0.0;
    // Find out how many nodes there are
    // const unsigned n_u_node = this->nnode();
    const unsigned n_w_node = nw_node();
    // Find out how many nodes positional dofs there are [zdec] assumes each
    // node has the same dofs
    const unsigned n_w_nodal_type = nw_type_at_each_node();
    // Find out how many internal types there are
    const unsigned n_w_internal_type = nw_type_internal();

    // Find the dimension of the element [zdec] will this ever not be 2?
    const unsigned dim = this->dim();
    // The number of first derivatives is the dimension of the element
    const unsigned n_deriv = dim;
    // The number of second derivatives is the triangle number of the dimension
    const unsigned n_2deriv = dim * (dim + 1) / 2;

    // Vector of local coordinates
    Vector<double> s(dim);

    // Vector for coordintes
    Vector<double> x(dim);

    // Set the value of n_intpt
    unsigned n_intpt = this->integral_pt()->nweight();

    // Exact solution Vector (here a scalar)
    Vector<double> exact_soln(this->required_nvalue(0), 0.0);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < this->dim(); i++)
      {
        s[i] = this->integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      double w = this->integral_pt()->weight(ipt);

      // Local out-of-plane nodal basis and test funtions
      Shape psi_n_w(n_w_node, n_w_nodal_type);
      DShape dpsi_n_wdxi(n_w_node, n_w_nodal_type, n_deriv);
      DShape d2psi_n_wdxi2(n_w_node, n_w_nodal_type, n_2deriv);
      // Local out-of-plane internal basis and test functions
      Shape psi_i_w(n_w_internal_type);
      DShape dpsi_i_wdxi(n_w_internal_type, n_deriv);
      DShape d2psi_i_wdxi2(n_w_internal_type, n_2deriv);

      // Call the derivatives of the shape and test functions for the out of
      // plane unknown
      double J = d2basis_w_eulerian_foeppl_von_karman(s,
                                                      psi_n_w,
                                                      psi_i_w,
                                                      dpsi_n_wdxi,
                                                      dpsi_i_wdxi,
                                                      d2psi_n_wdxi2,
                                                      d2psi_i_wdxi2);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Get x position as Vector
      interpolated_x(s, x);

      // Get FE function value
      Vector<double> u_fe(this->required_nvalue(0), 0.0);
      u_fe = interpolated_fvk_disp_and_deriv(s);

      // Get exact solution at this point
      (*exact_soln_pt)(x, exact_soln);

      // Output x,y,...,error
      for (unsigned i = 0; i < this->dim(); i++)
      {
        outfile << x[i] << " ";
      }
      for (unsigned ii = 0; ii < this->required_nvalue(0); ii++)
      {
        outfile << exact_soln[ii] << " " << exact_soln[ii] - u_fe[ii] << " ";
      }
      outfile << std::endl;

      // Loop over variables
      double tmp1 = 0.0, tmp2 = 0.0;
      for (unsigned ii = 0; ii < 1; ii++)
      {
        // Add to error and norm
        tmp1 = (exact_soln[ii] * exact_soln[ii] * W);
        tmp2 = ((exact_soln[ii] - u_fe[ii]) * (exact_soln[ii] - u_fe[ii]) * W);
        norm += tmp1;
        error += tmp2;
      }
    } // End of loop over integration pts
  }
  // HERE OVERLOAD COMPUTE_ERROR SO WE CAN GET A VECTOR OF NORMS


  //======================================================================
  /// Validate against exact solution
  ///
  /// Solution is provided via function pointer.
  /// Plot error at a given number of plot points.
  //======================================================================
  void FoepplVonKarmanEquations::compute_error(
    std::ostream& outfile,
    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
    double& error,
    double& norm)
  {
    // Initialise
    error = 0.0;
    norm = 0.0;

    // Find out how many nodes there are
    // const unsigned n_u_node = this->nnode();
    const unsigned n_w_node = nw_node();
    // Find out how many nodes positional dofs there are
    unsigned n_w_nodal_type = nw_type_at_each_node(); // [zdec]
    // Find out how many bubble nodes there are
    const unsigned n_w_internal_type = nw_type_internal();

    // Find the dimension of the element [zdec] will this ever not be 2?
    const unsigned dim = this->dim();
    // The number of first derivatives is the dimension of the element
    const unsigned n_deriv = dim;
    // The number of second derivatives is the triangle number of the dimension
    const unsigned n_2deriv = dim * (dim + 1) / 2;

    // Vector of local coordinates

    // Vector of local coordinates
    Vector<double> s(this->dim());

    // Vector for coordintes
    Vector<double> x(this->dim());

    // Set the value of n_intpt
    unsigned n_intpt = this->integral_pt()->nweight();

    // Tecplot
    // outfile << "ZONE" << std::endl;

    // Exact solution Vector (here a scalar)
    Vector<double> exact_soln(this->required_nvalue(0), 0.0);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Assign values of s
      for (unsigned i = 0; i < this->dim(); i++)
      {
        s[i] = this->integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      double w = this->integral_pt()->weight(ipt);

      // Local out-of-plane nodal basis and test funtions
      Shape psi_n_w(n_w_node, n_w_nodal_type);
      DShape dpsi_n_wdxi(n_w_node, n_w_nodal_type, n_deriv);
      DShape d2psi_n_wdxi2(n_w_node, n_w_nodal_type, n_2deriv);
      // Local out-of-plane internal basis and test functions
      Shape psi_i_w(n_w_internal_type);
      DShape dpsi_i_wdxi(n_w_internal_type, n_deriv);
      DShape d2psi_i_wdxi2(n_w_internal_type, n_2deriv);

      // Call the derivatives of the shape and test functions for the out of
      // plane unknown
      double J = d2basis_w_eulerian_foeppl_von_karman(s,
                                                      psi_n_w,
                                                      psi_i_w,
                                                      dpsi_n_wdxi,
                                                      dpsi_i_wdxi,
                                                      d2psi_n_wdxi2,
                                                      d2psi_i_wdxi2);

      // Premultiply the weights and the Jacobian
      double W = w * J;
      // double Wlin = w*Jlin;

      // Get x position as Vector
      interpolated_x(s, x);

      // Get FE function value
      Vector<double> u_fe(this->required_nvalue(0), 0.0);
      u_fe = interpolated_fvk_disp_and_deriv(s);

      // Get exact solution at this point
      (*exact_soln_pt)(x, exact_soln);

      // Loop over variables
      double tmp1 = 0.0, tmp2 = 0.0;
      if (Error_metric_fct_pt == 0)
      {
        // Output x,y,...,error
        for (unsigned i = 0; i < this->dim(); i++)
        {
          outfile << x[i] << " ";
        }
        for (unsigned ii = 0; ii < this->required_nvalue(0); ii++)
        {
          outfile << exact_soln[ii] << " " << exact_soln[ii] - u_fe[ii] << " ";
        }
        outfile << std::endl;

        for (unsigned ii = 0; ii < 1; ii++)
        {
          // Add to error and norm
          tmp1 = (exact_soln[ii] * exact_soln[ii]);
          tmp2 = ((exact_soln[ii] - u_fe[ii]) * (exact_soln[ii] - u_fe[ii]));
        }
      }
      else
      {
        // Get the metric
        (*Error_metric_fct_pt)(x, u_fe, exact_soln, tmp2, tmp1);
      }
      norm += tmp1 * W;
      error += tmp2 * W;
    } // End of loop over integration pts
  } // End of compute_error(...)


  //============================================================================
  /// Return FE representation of the three displacements at local coordinate
  /// s
  ///   0: u_x
  ///   1: u_y
  ///   2: w
  //============================================================================
  Vector<double> FoepplVonKarmanEquations::interpolated_fvk_disp(
    const Vector<double>& s) const
  {
    // The number of in-plane and out-of-plane unknowns
    const unsigned n_u_displacements = 2;
    const unsigned n_w_displacements = 1;

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
      throw OomphLibError("The number of internal basis types for u is non-zero\
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
    Vector<double> interpolated_u(n_u_displacements, 0.0);
    // Create space for out-of-plane interpolated unknowns
    Vector<double> interpolated_w(n_w_displacements, 0.0);

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
        for (unsigned alpha = 0; alpha < n_u_displacements; alpha++)
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
  } // end of interpolated_fvk_disp(...)



  // [zdec] should this be in the .cc?
  //============================================================================
  /// Return FE representation of the three displacements and their
  /// derivatives at local coordinate s:
  ///   0: u_x          6: du_y/dy
  ///   1: u_y          7: dw/dx
  ///   2: w            8: dw/dy
  ///   3: du_x/dx      9: d2w/dx2
  ///   4: du_x/dy     10: d2w/dxdy
  ///   5: du_y/dx     11: d2w/dy2
  //============================================================================
  Vector<double> FoepplVonKarmanEquations::interpolated_fvk_disp_and_deriv(
    const Vector<double>& s) const
  {
    // The number of in-plane and out-of-plane unknowns
    const unsigned n_u_displacements = 2;
    const unsigned n_w_displacements = 1;

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
      throw OomphLibError("The number of internal basis types for u is non-zero\
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
    Vector<double> interpolated_u(n_u_displacements, 0.0);
    Vector<double> interpolated_dudt(n_u_displacements, 0.0);
    DenseMatrix<double> interpolated_dudxi(n_u_displacements, n_deriv, 0.0);
    // Create space for out-of-plane interpolated unknowns
    Vector<double> interpolated_w(n_w_displacements, 0.0);
    Vector<double> interpolated_dwdt(n_w_displacements, 0.0);
    DenseMatrix<double> interpolated_dwdxi(n_w_displacements, n_deriv, 0.0);
    DenseMatrix<double> interpolated_d2wdxi2(n_w_displacements, n_2deriv, 0.0);

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
        for (unsigned alpha = 0; alpha < n_u_displacements; alpha++)
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
    TimeStepper* timestepper_pt = this->w_internal_data_pt()->time_stepper_pt();
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
  // End of interpolated_fvk_disp_and_deriv(...)

} // namespace oomph

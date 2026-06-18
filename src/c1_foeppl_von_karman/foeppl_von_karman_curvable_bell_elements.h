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
#include "src/generic/c1_plate_helper.h"
#include "src/generic/subparametric_Telement.h"
#include "src/generic/oomph_definitions.h"

// [zdec] TODO:
// -- Move function definitions to .cc
// -- Optimise basis and test functions by not retrieving d and d2 basis every
// time
// --

namespace oomph
{

 
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
  class FvKDuplicateNodeConstraintElement : public virtual DuplicateNodeConstraintElement
  {
  public:
   
    /// Construcor. Needs the two node pointers so that we can retrieve the
    /// boundary data at solve time
    FvKDuplicateNodeConstraintElement(
      Node* const& left_node_pt,
      Node* const& right_node_pt,
      C1CurviLine* const& left_boundary_pt,
      C1CurviLine* const& right_boundary_pt,
      Vector<double> const& left_coord,
      Vector<double> const& right_coord)
     : DuplicateNodeConstraintElement(
      left_node_pt,
      right_node_pt,
      left_boundary_pt,
      right_boundary_pt,
      left_coord,
      right_coord)      
    {
      // Add internal data which stores the eight Lagrange multipliers
      Index_of_lagrange_data = add_internal_data(new Data(8));

      // Add each node as external data
      Index_of_left_data = add_external_data(Left_node_pt);
      Index_of_right_data = add_external_data(Right_node_pt);
    }

    /// Destructor
    ~FvKDuplicateNodeConstraintElement()
    {
     unsigned n_internal=ninternal_data();
     for (unsigned i=0;i<n_internal;i++)
      {
       delete this->internal_data_pt(i);
      }
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
    void pin_redundant_constraints()
    {
      // Start by unpinning all lagrange multipliers in case the boundary
      // conditions are less restrictive than previously
      internal_data_pt(Index_of_lagrange_data)->unpin_all();


      // [zdec] This full description might be overkill for the code but it will
      // go in my thesis.

      // We need to keep track of which fvk dofs are already 'used' by Lagrange
      // constraints. If dofs 3 and 4 in the right node (dw/dl_1, dw/dl_2) are
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

      // The number of nodal types per field
      unsigned n_type = 6;

      // We use a vector of booleans to keep track of dofs that might be reused
      // (no need to track right dofs which are used once)
      std::vector<bool> right_data_used(n_type, false);
      std::vector<bool> left_data_used(n_type, false);

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
      for (unsigned k_type = 0; k_type < 3; k_type++)
      {
	// Index of the condition on the element
	unsigned condition_index = k_type;
        // Index of the val associated with displacement in the right node
        unsigned right_ui_index = k_type;
        // Index of the val associated with displacement in the left node
        unsigned left_ui_index = k_type;

	// Get whether each value is pinned
	bool right_ui_pinned = right_data_pt->is_pinned(right_ui_index);
        bool left_ui_pinned = left_data_pt->is_pinned(left_ui_index);

        // If anything is free, mark it as used and continue without doing
        // anything else
        if (!right_ui_pinned && !right_data_used[right_ui_index])
        {
          // // [zdec] debug
          // std::cout << "eqn " << condition_index
	  //           << " depends on dof R" << right_ui_index
          //           << std::endl;
          right_data_used[right_ui_index] = true;
        }
        else if (!left_ui_pinned && !left_data_used[left_ui_index])
        {
          // // [zdec] debug
          // std::cout << "eqn " << condition_index
	  //           << " depends on dof L" << left_ui_index
          //           << std::endl;
          left_data_used[left_ui_index] = true;
        }
	else
	{
	  // ---------------------------------------------------------------------
	  // If we made it here, it is because all dofs in the constraint are
	  // pinned so we need to check the constraint is satisfied manually and
	  // then remove it by pinning the corresponding lagrange multiplier
         
         // hierher Aidan: should this remain alive?

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
	  internal_data_pt(Index_of_lagrange_data)->pin(condition_index);
	}
      } // End for loop over first three conditions [k_type]


      // Constraints 3-4 use dofs 3-4 (first derivatives of w) respectively from
      // the right node and both in the left
      for (unsigned alpha = 0; alpha < 2; alpha++)
      {
	// Index of the condition on the element
	unsigned condition_index = 3 + alpha;
	// Index of the right nodes alpha-th derivative value
        unsigned right_dwda_index = 3 + alpha;
	// Index of the left nodes first derivative value
	unsigned left_dwd1_index = 3;
        // Index of the left nodes second derivative value
        unsigned left_dwd2_index = 4;

        // Get whether each nodal value is pinned
        bool right_dwda_pinned = right_data_pt->is_pinned(right_dwda_index);
        bool left_dwd1_pinned = left_data_pt->is_pinned(left_dwd1_index);
        bool left_dwd2_pinned = left_data_pt->is_pinned(left_dwd2_index);

        // If anything is free, mark it as used and continue without doing
        // anything else. We also need to check that each dof hasn't become
        // decoupled from this constraint by ensuring that its coefficient (if
        // it has one) is sufficiently large (> Orthogonality_tolerance)
        if (!right_dwda_pinned && !right_data_used[right_dwda_index])
        {
          // // [zdec] debug
          // std::cout << "eqn " << condition_index
	  //           << " depends on dof R" << right_dwda_index
          //           << std::endl;
          right_data_used[right_dwda_index] = true;
          continue;
        }
        if (!left_dwd1_pinned && !left_data_used[left_dwd1_index])
        {
          // // [zdec] debug
          // std::cout << "eqn " << condition_index
	  //           << " depends on dof L" << left_dwd1_index
	  //           << std::endl;
          double coeff = jac_of_transform(0, alpha);
          if (fabs(coeff) > Orthogonality_tolerance)
          {
            left_data_used[left_dwd1_index] = true;
            continue;
          }
        }
        if (!left_dwd2_pinned && !left_data_used[left_dwd2_index])
        {
          // // [zdec] debug
          // std::cout << "eqn " << condition_index
	  //           << " depends on dof L" << left_dwd2_index
	  //           << std::endl;
          double coeff = jac_of_transform(1, alpha);
          if (fabs(coeff) > Orthogonality_tolerance)
          {
            left_data_used[left_dwd2_index] = true;
            continue;
          }
        }
        // ---------------------------------------------------------------------
        // If we made it here, it is because all dofs in the constraint are
        // pinned so we need to check the constraint is satisfied manually and
        // then remove it by pinning the corresponding lagrange multiplier

        // hierher Aidan: should this remain alive?

        
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
        internal_data_pt(Index_of_lagrange_data)->pin(condition_index);
      }

      // Constraints 5-7 use dofs 5-7 respectively from the right node and
      // all use dofs 3-7 from the left node
      for (unsigned alpha = 0; alpha < 2; alpha++)
      {
        // beta>=alpha so we dont double count constraint 6
        for (unsigned beta = alpha; beta < 2; beta++)
        {
          // The index of the constraint
          unsigned condition_index = 5 + alpha + beta;
          // Index of the right nodes alpha-beta-th derivative value
          unsigned right_dwdadb_index = 5 + alpha + beta;
          // Index of the left nodes d1 derivative value
          unsigned left_dwd1_index = 3;
          // Index of the left nodes d2 derivative value
          unsigned left_dwd2_index = 4;
          // Index of the left nodes d1d1 derivative value
          unsigned left_dwd1d1_index = 5;
          // Index of the left nodes d1d2 derivative value
          unsigned left_dwd1d2_index = 6;
          // Index of the left nodes d2d2 derivative value
          unsigned left_dwd2d2_index = 7;

          // Get whether each nodal value is pinned
          bool right_dwdadb_pinned =
	    right_data_pt->is_pinned(right_dwdadb_index);
          bool left_dwd1_pinned = left_data_pt->is_pinned(left_dwd1_index);
          bool left_dwd2_pinned = left_data_pt->is_pinned(left_dwd2_index);
          bool left_dwd1d1_pinned = left_data_pt->is_pinned(left_dwd1d1_index);
          bool left_dwd1d2_pinned = left_data_pt->is_pinned(left_dwd1d2_index);
          bool left_dwd2d2_pinned = left_data_pt->is_pinned(left_dwd2d2_index);

          // If anything is free, mark it as used and continue without doing
          // anything else. We also need to check that each dof hasn't become
          // decoupled from this constraint by ensuring that its coefficient (if
          // it has one) is sufficiently large (> Orthogonality_tolerance)
          if (!right_dwdadb_pinned && !right_data_used[right_dwdadb_index])
          {
            // // [zdec] debug
            // std::cout << "eqn " << condition_index
	    //           << " depends on dof R" << right_dwdadb_index
            //           << std::endl;
            right_data_used[right_dwdadb_index] = true;
            continue;
          }
          if (!left_dwd1_pinned && !left_data_used[left_dwd1_index])
          {
            double coeff = hess_of_transform[0](alpha, beta);
            if (fabs(coeff) > Orthogonality_tolerance)
            {
              // // [zdec] debug
              // std::cout << "eqn " << condition_index
	      //   	<< " depends on dof L" << left_dwd1_index
	      //   	<< std::endl;
              left_data_used[left_dwd1_index] = true;
              continue;
            }
          }
          if (!left_dwd2_pinned && !left_data_used[left_dwd2_index])
          {
            double coeff = hess_of_transform[1](alpha, beta);
            if (fabs(coeff) > Orthogonality_tolerance)
            {
              // // [zdec] debug
              // std::cout << "eqn " << condition_index
	      //   	<< " depends on dof L" << left_dwd2_index
	      //   	<< std::endl;
              left_data_used[left_dwd2_index] = true;
              continue;
            }
          }
          if (!left_dwd1d1_pinned && !left_data_used[left_dwd1d1_index])
          {
            double coef =
              jac_of_transform(0, alpha) * jac_of_transform(0, beta);
            if (fabs(coef) > Orthogonality_tolerance)
            {
              // // [zdec] debug
              // std::cout << "eqn " << condition_index
	      //   	<< " depends on dof L" << left_dwd1d1_index
	      //   	<< std::endl;
              left_data_used[left_dwd1d1_index] = true;
              continue;
            }
          }
          if (!left_dwd1d2_pinned && !left_data_used[left_dwd1d2_index])
          {
            double coef =
              jac_of_transform(0, alpha) * jac_of_transform(1, beta) +
              jac_of_transform(1, alpha) * jac_of_transform(0, beta);
            if (fabs(coef) > Orthogonality_tolerance)
            {
              // // [zdec] debug
              // std::cout << "eqn " << condition_index
	      //   	<< " depends on dof L" << left_dwd1d2_index
	      //   	<< std::endl;
              left_data_used[left_dwd1d2_index] = true;
              continue;
            }
          }
          if (!left_dwd2d2_pinned && !left_data_used[left_dwd2d2_index])
          {
            double coef =
              jac_of_transform(1, alpha) * jac_of_transform(1, beta);
            if (fabs(coef) > Orthogonality_tolerance)
            {
              // // [zdec] debug
              // std::cout << "eqn " << condition_index
	      //   	<< " depends on dof L" << left_dwd2d2_index
	      //   	<< std::endl;
              left_data_used[left_dwd2d2_index] = true;
              continue;
            }
          }
          // -------------------------------------------------------------------
          // If we made it here, it is because all dofs in the constraint are
          // pinned so we need to check the constraint is satisfied manually and
          // then remove it by pinning the corresponding lagrange multiplier

          // hierher Aidan: should this remain alive?

          
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
          internal_data_pt(Index_of_lagrange_data)->pin(condition_index);
        }
      }
    } // End pin_redundant_constraints()


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




    /// Add the contribution to the residuals (and jacobain if flag is 1) from
    /// the Lagrange multiplier constraining equations
    void fill_in_generic_residual_contribution_constraint(
      Vector<double>& residuals,
      DenseMatrix<double>& jacobian,
      const unsigned& flag)
    {
      // // [zdec] debug
      // std::cout << std::endl
      //           << std::endl
      //           << "ADD CONTRIBUTION FROM CONSTRAINTS" << std::endl
      //           << "=============================================" << std::endl;
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

          // If flag, then add the jacobian contribution
          if (flag)
          {
            // The contributions to the right node's equations are just
            // r_i * L_i
            int lagrange_dof_number =
              internal_local_eqn(Index_of_lagrange_data, k_type);
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


      //----------------------------------------------------------------------
      // Next, the contributions to the left node external equations
      // First three are displacements:
      //     - lambda_i*(U_alpha or W_0)
      for (unsigned k_type = 0; k_type < 3; k_type++)
      {
	// Index of the left nodes i-th displacement (0-th type) value
	unsigned left_ui_index = k_type;
        // External equation index
        int left_eqn_number =
          external_local_eqn(Index_of_left_data, left_ui_index);
        // If this dof isn't pinned we add to the residual
        if (left_eqn_number >= 0)
        {
	  // Add residual contribution which comes from lagrange value
	  unsigned lagrange_index = k_type;
          residuals[left_eqn_number] += -lagrange_value[lagrange_index];

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
      } // End loop adding contribution to the left nodal displacement equations

      // Next two are from gradient of w:
      //     - lambda_{3+\beta} * w_{1+\alpha} * J_{\alpha\beta}
      //     - lambda_{5+\beta+\gamma} * w_{1+\alpha} * H_{\alpha\beta\gamma}
      // gamma>=beta so we don't double count lambda_6 condition
      for (unsigned alpha = 0; alpha < 2; alpha++)
      {
	// Index of the left nodes d1 derivative value
	unsigned left_wda_index = 3 + alpha;
        // Eqn number is the index of the alpha-th derivative of w which is
        // the 1+alpha-th dof (alpha=0,1)
        int left_eqn_number =
          external_local_eqn(Index_of_left_data, left_wda_index);
        // If this dof isn't pinned we add to the residual
        if (left_eqn_number >= 0)
        {
	  // Loop over the lagrange multipliers associated with the right
	  // first derivatives
          for (unsigned beta = 0; beta < 2; beta++)
          {
	    // Add residual contribution from the lagrange value associated
	    // with the right beta-th derivative
            unsigned lagrange_wdb_index = 3 + beta;
            residuals[left_eqn_number] +=
              -lagrange_value[lagrange_wdb_index]
	      * jac_of_transform(alpha, beta);

	    // Get the local equation number for this dof and add the jacobian
	    // contribution if unpinned (and we are making the jacobian)
            int lagrange_dof_number =
              internal_local_eqn(Index_of_lagrange_data, lagrange_wdb_index);
            if (flag && (lagrange_dof_number >= 0))
            {
              double jac_term = -jac_of_transform(alpha, beta);
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

            // gamma>=beta so we don't double count the
            // lagrange_value[6] constraint
            for (unsigned gamma = beta; gamma < 2; gamma++)
            {
	      // Add residual contribution from the lagrange value associated
	      // with the right beta+gamma-th second derivative
	      unsigned lagrange_wdbdg_index = 5 + beta + gamma;
	      residuals[left_eqn_number] +=
		- lagrange_value[lagrange_wdbdg_index]
		* hess_of_transform[alpha](beta, gamma);

	      // Get the local equation number for this dof and add the
	      // jacobian contribution if unpinned (and we are making the
	      // jacobian)
	      int lagrange_dof_number =
		internal_local_eqn(Index_of_lagrange_data, lagrange_wdbdg_index);
	      if (flag && (lagrange_dof_number >= 0))
	      {
		double jac_term = -hess_of_transform[alpha](beta, gamma);
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

      // Last three are the second derivatives of w (delta>gamma):
      //     - lambda_{5+\gamma+\delta} * w_{3+\alpha+\beta}
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
          unsigned left_wdadb_index = 5 + alpha + beta;
          // Eqn number is the index of the second derivative of w
          int left_eqn_number =
            external_local_eqn(Index_of_left_data, left_wdadb_index);
          // If this dof isn't pinned we add to the residual
          if (left_eqn_number >= 0)
          {
            // Index constraint using gamma and delta
            for (unsigned gamma = 0; gamma < 2; gamma++)
            {
              // delta>=gamma so we don't double count the lagrange_value
	      // associated with the mixed derivative
              for (unsigned delta = gamma; delta < 2; delta++)
              {
		// Add residual contribution
		unsigned lagrange_wdgdd_index = 5 + gamma + delta;
                residuals[left_eqn_number] +=
		  -lagrange_value[lagrange_wdgdd_index] *
		  jac_of_transform(alpha, gamma) *
		  jac_of_transform(beta, delta);

                // Get the local equation number for this dof and add the
                // jacobian contribution if unpinned (and we are making the
                // jacobian)
                int lagrange_dof_number = internal_local_eqn(
                  Index_of_lagrange_data, lagrange_wdgdd_index);
                if (flag && (lagrange_dof_number >= 0))
                {
                  // Find the jacobian matrix contribution
                  double jac_term = -jac_of_transform(alpha, gamma) *
                                    jac_of_transform(beta, delta);
                  // Orthogonality check
                  if (fabs(jac_term) > Orthogonality_tolerance)
                  {
                    // Add the contribution to the jacobian
                    jacobian(left_eqn_number, lagrange_dof_number) += jac_term;
                    // And by symmetry, we can add the transpose contribution
                    // to the jacobian
                    jacobian(lagrange_dof_number, left_eqn_number) += jac_term;
                  }
                  } // End jacobian
              }
            } // End loops over the conditions (gamma,delta)
          } // End if dof isn't pinned
        }
      } // End loops adding contributions to the left nodal curvature equations
      // (alpha,beta)


      //----------------------------------------------------------------------
      // Now add contributions to the internal (lagrange multiplier) equations
      // (note jacobian contributions will have already been thanks to use of
      // symmetry)

      // First three (u,v,w) dofs are equal
      for (unsigned i_dof = 0; i_dof < 3; i_dof++)
      {
	// Index of the condition on the nodes (displacement constraint
	// corresponds to zeroth type for each field)
        unsigned lagrange_index = i_dof;
        // Get the internal data eqn number for this constraint
        int lagrange_eqn_number =
          internal_local_eqn(Index_of_lagrange_data, lagrange_index);
	// If this dof isn't pinned we add to the residual
	if (lagrange_eqn_number >= 0)
        {
          // Add contributions from right and left nodes displacement
          // (zero-th) type
          unsigned right_ui_index = i_dof;
          unsigned left_ui_index = i_dof;
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
        unsigned lagrange_index = 3 + alpha;
        // Get the internal data eqn number for this constraint
        int lagrange_eqn_number =
          internal_local_eqn(Index_of_lagrange_data, lagrange_index);
        // If this dof isn't pinned we add to the residual
        if (lagrange_eqn_number >= 0)
        {
          // Add contribution from right node
          unsigned right_wda_index = 3 + alpha;
          residuals[lagrange_eqn_number] += (right_value[right_wda_index]);
          // Add contribuions from left node
          for (unsigned beta = 0; beta < 2; beta++)
          {
            unsigned left_wdb_index = 3 + beta;
            residuals[lagrange_eqn_number] +=
              -left_value[left_wdb_index] * jac_of_transform(beta, alpha);
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
          unsigned lagrange_index = 5 + alpha + beta;
          // Get the internal data eqn number for this constraint
          int lagrange_eqn_number =
            internal_local_eqn(Index_of_lagrange_data, lagrange_index);
          // If this dof isn't pinned we add to the residual
          if (lagrange_eqn_number >= 0)
          {
            // Add contributions from right node
            unsigned right_wdadb_index = 5 + alpha + beta;
            residuals[lagrange_eqn_number] += right_value[right_wdadb_index];
            // Loop over the left node derivatives
            for (unsigned gamma = 0; gamma < 2; gamma++)
            {
              // Add contributions from left node first derivatives
              unsigned left_wdg_index = 3 + gamma;
              residuals[lagrange_eqn_number] +=
                -left_value[left_wdg_index] *
                hess_of_transform[gamma](alpha, beta);
              // Loop over the left derivatives again to get second
              // derivatives
              for (unsigned delta = 0; delta < 2; delta++)
              {
                // Add contributions from left node second derivatives
		unsigned left_wdgdd_index = 5 + gamma + delta;
		residuals[lagrange_eqn_number] +=
                  -left_value[left_wdgdd_index] *
                  jac_of_transform(gamma, alpha) *
                  jac_of_transform(delta, beta);
              }
            }
            // Jacobian is added during right and left nodal equations using
            // symmetry
          } // End if eqn not pinned
	}
      }
    } // End fill_in_generic_residual_contribution_constraint


  }; // End of FvKDuplicateNodeConstraintElement class definition



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

      // hierher Aidan: Can this go; is this handled elsewhere?

      // // Use the higher order integration scheme
      // delete this->integral_pt();
      // // Use the higher order integration scheme
      // TGauss<2, 4>* new_integral_pt = new TGauss<2, 4>;
      // this->set_integration_scheme(new_integral_pt);

      // Rotated dof helper
      Rotated_boundary_helper_pt = new RotatedBoundaryHelper;
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

    /// Full output function with a rich set of unknowns:
    ///  x, y, ux, uy, w, dw, ddw, du, strain, stress, principal stress
   void full_output(std::ostream& outfile, const unsigned& nplot)
    {
     FoepplVonKarmanEquations::full_output(outfile,nplot);
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

   
   /// Clamp: i.e. pin the in-plane displacements and pin the out-of-plane
   /// displacement and its normal derivative. We also apply implied
   /// boundary conditions (e.g. specification of dw/dn also implies
   /// d^2w/dn/dzeta etc. boundary_values_pt[i] describes boundary conditions
   /// for the three displacement components (i=0,1 in plane (x,y);
   /// i=2: out-of-plane (z)).
   /// curviline_pt provides a pointer to the representation of the curvilinear
   /// boundary in the triangle mesh.
   virtual void fully_clamp_specified_boundary(
    const unsigned& b,
    const Vector<BoundaryConditionForC1PlateBending*>& boundary_values_pt,
    TriangleMeshCurviLine* curviline_pt);
   
   
   /// Pin i.e. pin the in-plane and out of plane displacements only.
   /// We also apply implied boundary conditions (e.g. specification of w
   /// also implies dw/dt etc. boundary_values_pt[i] describes boundary conditions
   /// for the three displacement components (i=0,1 in plane (x,y);
   /// i=2: out-of-plane (z)).
   /// curviline_pt provides a pointer to the representation of the curvilinear
   /// boundary in the triangle mesh.
   virtual void pin_specified_boundary(
    const unsigned& b,
    const Vector<BoundaryConditionForC1PlateBending*>& boundary_values_pt,
    TriangleMeshCurviLine* curviline_pt);

   

   
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

   /// Factory to create DuplicateNodeConstraintElement.
   /// which ensures that the deformation is sufficiently smooth
   /// between different parts of a boundary (which may contain isolated
   /// kinks which make it C0 rather than C1). Pass:
   /// -- Pointers to nodes on the "left" and "right" boundary
   /// -- pointers to the C1Curvilines that describe the smooth
   ///    parts of the boundary on either side
   /// -- the boundary coordinates that identifies the corner point
   ///    relative the right and left boundary parametrisation
   ///    (specified via a one-sized vector). 
   DuplicateNodeConstraintElement* duplicate_constraint_element_factory(
    Node* const& left_node_pt,
    Node* const& right_node_pt,
    C1CurviLine* const& left_boundary_pt,
    C1CurviLine* const& right_boundary_pt,
    Vector<double> const& left_coord,
    Vector<double> const& right_coord) override
    {
     return new FvKDuplicateNodeConstraintElement(
      left_node_pt,
      right_node_pt,
      left_boundary_pt,
      right_boundary_pt,
      left_coord,
      right_coord);
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

    /// Required  # of `values' (pinned or dofs)
    /// at node n
    inline unsigned required_nvalue(const unsigned& n) const
    {
      return Initial_Nvalue[n];
    }



   // protected:
  public: // hierher change back to protected once we've dropped the interpolation tests back to lower classes

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

    /// Get the jacobian of the local to global mapping
    virtual double local_to_eulerian_mapping(
      const Vector<double>& s,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& inverse_jacobian) const;

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


  private:

   /// Helper function to impose alpha-th in-plane displacements
   ///(0 or 1 for x or y displacements) according to scalar function
   /// specified in boundary_values_pt. curvline_pt provides the
   /// geometry of the boundary in terms of a not-necessarily-arclength
   /// coordinate zeta.
   void pin_and_impose_specified_in_plane_displacement_along_specified_boundary(
    const unsigned& alpha,
    const unsigned& b,
    BoundaryConditionForC1PlateBending* boundary_values_pt,
    TriangleMeshCurviLine* curviline_pt);

   /// Helper function to impose i-th displacements
   /// (0 or 1 for x or y displacements; 2 for out of plane)
   /// according to scalar function specified in boundary_values_pt.
   /// curvline_pt provides the geometry of the boundary in terms of
   /// a not-necessarily-arclength coordinate zeta.
   void pin_and_impose_specified_displacement_along_specified_boundary(
    const unsigned& i,
    const unsigned& b,
    BoundaryConditionForC1PlateBending* boundary_values_pt,
    TriangleMeshCurviLine* curviline_pt);
     
   /// Helper function to impose out of plane displacements
   /// according to scalar function specified in boundary_values_pt.
   /// curvline_pt provides the geometry of the boundary in terms of
   /// a not-necessarily-arclength coordinate zeta.
   void clamp_and_impose_specified_out_of_plane_displacement_along_specified_boundary(
    const unsigned& b,
    BoundaryConditionForC1PlateBending* boundary_values_pt,
    TriangleMeshCurviLine* curviline_pt);

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
    TElement<2, NNODE_1D>::
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
    TElement<2, NNODE_1D>::dshape_local(s, psi_n, dpsids);
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
    TElement<2, NNODE_1D>::dshape_local(s, psi_n, dpsids);
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
  /// Out-of-plane basis functions at local coordinate s
  //======================================================================
  template<unsigned NNODE_1D>
  void FoepplVonKarmanC1CurvableBellElement<
    NNODE_1D>::basis_w_foeppl_von_karman(const Vector<double>& s,
                                         Shape& psi_n,
                                         Shape& psi_i) const
  {

   // hierher Aidan: Kill this commented out bit?
//     throw OomphLibError("This still needs testing for curved elements.",
//                         "void FoepplVonKarmanC1CurvableBellElement<NNODE_1D>::shape_and_test_foeppl_von_karman(...)",
//                         OOMPH_EXCEPTION_LOCATION);

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
  /// Get the jacobian of the local to global mapping
  //======================================================================
  template<unsigned NNODE_1D>
  double FoepplVonKarmanC1CurvableBellElement<NNODE_1D>::
    local_to_eulerian_mapping(const Vector<double>& s,
			      DenseMatrix<double>& jacobian,
			      DenseMatrix<double>& inverse_jacobian) const
  {
    return CurvableBellElement<NNODE_1D>::
      local_to_eulerian_mapping(s,
				jacobian,
				inverse_jacobian);
  }

  //======================================================================
  /// Rotate the shape functions according to the specified basis on the
  /// boundary. This function does a DenseDoubleMatrix solve to determine
  /// new basis - which could be speeded up by caching the matrices higher
  /// up and performing the LU decomposition only once
 //
 // hierher shouldn't this live further down?
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



 //==========================================================================
 /// Clamp: i.e. pin the in-plane displacements and pin the out-of-plane
 /// displacement and its normal derivative. We also apply implied
 /// boundary conditions (e.g. specification of dw/dn also implies
 /// d^2w/dn/dzeta etc. boundary_values_pt[i] describes boundary conditions
 /// for the three displacement components (i=0,1 in plane (x,y);
 /// i=2: out-of-plane (z)).
 /// curviline_pt provides a pointer to the representation of the curvilinear
 /// boundary in the triangle mesh.
 //==========================================================================
 template<unsigned NNODE_1D>
 void FoepplVonKarmanC1CurvableBellElement<NNODE_1D>::fully_clamp_specified_boundary(
  const unsigned& b,
  const Vector<BoundaryConditionForC1PlateBending*>& boundary_values_pt,
  TriangleMeshCurviLine* curviline_pt)
 {

  // Deal with in plane displacements
  for (unsigned i=0;i<2;i++)
   {
    pin_and_impose_specified_in_plane_displacement_along_specified_boundary(
     i,b,boundary_values_pt[i],curviline_pt);
   }
  
  // Deal with out of plane displacements
  clamp_and_impose_specified_out_of_plane_displacement_along_specified_boundary(
   b,boundary_values_pt[2],curviline_pt);
 }


 //=============================================================================
 /// Pin i.e. pin the in-plane and out of plane displacements only.
 /// We also apply implied boundary conditions (e.g. specification of w
 /// also implies dw/dt etc. boundary_values_pt[i] describes boundary conditions
 /// for the three displacement components (i=0,1 in plane (x,y);
 /// i=2: out-of-plane (z)).
 /// curviline_pt provides a pointer to the representation of the curvilinear
 /// boundary in the triangle mesh.
//=============================================================================
 template<unsigned NNODE_1D>
 void FoepplVonKarmanC1CurvableBellElement<NNODE_1D>::pin_specified_boundary(
  const unsigned& b,
  const Vector<BoundaryConditionForC1PlateBending*>& boundary_values_pt,
  TriangleMeshCurviLine* curviline_pt)
 {

  for (unsigned i=0;i<3;i++)
   {
    pin_and_impose_specified_displacement_along_specified_boundary
     (i,b,boundary_values_pt[i],curviline_pt);
   }

 }
 

 //=============================================================================
 /// Helper function to impose alpha-th in-plane displacements
 ///(0 or 1 for x or y displacements) according to scalar function
 /// specified in boundary_values_pt. curvline_pt provides the
 /// geometry of the boundary in terms of a not-necessarily-arclength
 /// coordinate zeta.
 //=============================================================================
 template<unsigned NNODE_1D>
 void FoepplVonKarmanC1CurvableBellElement<NNODE_1D>::
 pin_and_impose_specified_in_plane_displacement_along_specified_boundary(
  const unsigned& alpha,
  const unsigned& b,
  BoundaryConditionForC1PlateBending* boundary_values_pt,
  TriangleMeshCurviLine* curviline_pt)
 {
    
  // Initialise constants that we use in this function
  const unsigned field_index = u_alpha_field_index(alpha);
  const unsigned nodal_type_index =
   this->first_nodal_type_index_for_field(field_index);
  const unsigned n_node = nu_node();
  
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
    // Get boundary node
    BoundaryNode<Node>* nod_pt =
     dynamic_cast<BoundaryNode<Node>*>(this->node_pt(n));
    
    // Check if it is on the boundary
    if (nod_pt!=0)
     {
      if (nod_pt->is_on_boundary(b))
       {
        // We should only have one coordinate on this boundary
        unsigned nzeta=nod_pt->ncoordinates_on_boundary(b);
#ifdef PARANOID
        if (nzeta!=1)
         {
          std::stringstream error_message;
          error_message << "boundary coordinate must be 1D!"
                        << std::endl;
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
         }
#endif
        
        // Get value itself from boundary condition object
        Vector<double> zeta(nzeta);
        nod_pt->get_coordinates_on_boundary(b,zeta);
        double value=boundary_values_pt->f(zeta[0]);
                
        // Pin and set the value
        nod_pt->pin(alpha);
        nod_pt->set_value(nodal_type_index, value); // hierher Aidan should indices in pin and set_value be the same? 
       }
     }
   }
 }

 //=============================================================================
 /// Helper function to impose i-th displacements
 /// (0 or 1 for x or y displacements; 2 for out of plane)
 /// according to scalar function specified in boundary_values_pt.
 /// curvline_pt provides the geometry of the boundary in terms of
 /// a not-necessarily-arclength coordinate zeta.
 //=============================================================================
 template<unsigned NNODE_1D>
 void FoepplVonKarmanC1CurvableBellElement<NNODE_1D>::
 pin_and_impose_specified_displacement_along_specified_boundary(
  const unsigned& i,
  const unsigned& b,
  BoundaryConditionForC1PlateBending* boundary_values_pt,
  TriangleMeshCurviLine* curviline_pt)
 {

  // In plane
  if (i<2)
   {
    pin_and_impose_specified_in_plane_displacement_along_specified_boundary
     (i,b,boundary_values_pt,curviline_pt);
   }
  // Out of plane
  else
   {
    const unsigned w_index = w_field_index();
    const unsigned first_nodal_type_index =
     this->first_nodal_type_index_for_field(w_index);
    const unsigned n_vertices = nw_node();
    
    // Bell elements only have deflection dofs at vertices
    for (unsigned n = 0; n < n_vertices; ++n)
     {      
      // Get boundary node
      BoundaryNode<Node>* nod_pt =
       dynamic_cast<BoundaryNode<Node>*>(this->node_pt(n));
      
      // Check if it is on the boundary
      if (nod_pt!=0)
       {
        if (nod_pt->is_on_boundary(b))
         {
          // We should only have one coordinate on this boundary
          unsigned nzeta=nod_pt->ncoordinates_on_boundary(b);
#ifdef PARANOID
          if (nzeta!=1)
           {
            std::stringstream error_message;
            error_message << "boundary coordinate must be 1D!"
                          << std::endl;
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
           }
#endif
        
        // Get value itself from boundary condition object
        Vector<double> zeta(nzeta);
        nod_pt->get_coordinates_on_boundary(b,zeta);

        // Get the Jacobian from the GeomObject parametrising the boundary
        GeomObject* geom_obj_pt=curviline_pt->geom_object_pt();
        DenseMatrix<double> drdzeta(nzeta,2);
        geom_obj_pt->dposition(zeta,drdzeta);
        RankThreeTensor<double> d2rdzeta2(nzeta,2,2);
        geom_obj_pt->d2position(zeta,d2rdzeta2);
        double dtdzeta=sqrt(drdzeta(0,0)*drdzeta(0,0)+drdzeta(0,1)*drdzeta(0,1));

        // Foppl von Karman elements only have 6 Hermite deflection degrees
        // of freedom. They are {w ; w,x ; w,y ; w,xx ; w,xy ; w,yy}
        // or their rotated counterparts {w ; w,n ; w,t ; w,nn ; w,nt ; w,tt}
        unsigned nw_type=nw_type_at_each_node();
        for (unsigned k_type=0;k_type<nw_type;k_type++)
         {
          double value=0.0;
          switch (k_type)
           {
           case 0:
            // f itself
            value=boundary_values_pt->f(zeta[0]);
            break;
            
           case 1:
            // df/dn: Normal derivative left free
            // Calls broken virtual function; dies if not implemented
            //value=boundary_values_pt->dfdn(zeta[0]);
            break;
            
           case 2:
            // df/dt, including Jacobian
            value=boundary_values_pt->dfdzeta(zeta[0])/dtdzeta;
            break;
            
           case 3:
            // d^2f/dn^2: Second normal derivative shouldn't be set!
            break;
            
           case 4:
            // d^2f/dndt mixed second normal derivative: left free
            // value=boundary_values_pt->d2fdndzeta(zeta[0])/dtdzeta;
            break;
            
           case 5:
            // d^2f/dt^2  including Jacobian (twice; chain rule yourself to death)
            value=
             (boundary_values_pt->d2fdzeta2(zeta[0])*dtdzeta*dtdzeta-
              boundary_values_pt->dfdzeta(zeta[0])*
              (drdzeta(0,0)*d2rdzeta2(0,0,0)+drdzeta(0,1)*d2rdzeta2(0,0,1)))/
             (pow(dtdzeta,4));
            
            break;

           default:
            std::stringstream error_message;
            error_message << "never get here" << std::endl;
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
           }

          // Skip anything related to normal derivatives
          if ( ( k_type!=1) &&  // dw/dn
               ( k_type!=3) &&  // d^2w/dn^2
               ( k_type!=4)  )  // d^2w/dtdn
           
           {
            // Pin and set the value
            nod_pt->pin(first_nodal_type_index + k_type);
            nod_pt->set_value(first_nodal_type_index + k_type, value);
           }
         }
        
       }
     }

   }
 }


 }

 
 
 //=============================================================================
 /// Helper function to impose out of plane displacements
 /// according to scalar function specified in boundary_values_pt.
 /// curvline_pt provides the geometry of the boundary in terms of
 /// a not-necessarily-arclength coordinate zeta.
 //=============================================================================
 template<unsigned NNODE_1D>
 void FoepplVonKarmanC1CurvableBellElement<NNODE_1D>::
 clamp_and_impose_specified_out_of_plane_displacement_along_specified_boundary(
  const unsigned& b,
  BoundaryConditionForC1PlateBending* boundary_values_pt,
  TriangleMeshCurviLine* curviline_pt)
 {

  const unsigned w_index = w_field_index();
  const unsigned first_nodal_type_index =
   this->first_nodal_type_index_for_field(w_index);
  const unsigned n_vertices = nw_node();
  
  // Bell elements only have deflection dofs at vertices
  for (unsigned n = 0; n < n_vertices; ++n)
   {
    
    // Get boundary node
    BoundaryNode<Node>* nod_pt =
     dynamic_cast<BoundaryNode<Node>*>(this->node_pt(n));
    
    // Check if it is on the boundary
    if (nod_pt!=0)
     {
      if (nod_pt->is_on_boundary(b))
       {
        // We should only have one coordinate on this boundary
        unsigned nzeta=nod_pt->ncoordinates_on_boundary(b);
#ifdef PARANOID
        if (nzeta!=1)
         {
          std::stringstream error_message;
          error_message << "boundary coordinate must be 1D!"
                        << std::endl;
          throw OomphLibError(error_message.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
         }
#endif
        
        // Get value itself from boundary condition object
        Vector<double> zeta(nzeta);
        nod_pt->get_coordinates_on_boundary(b,zeta);

        oomph_info << "hierher FULLY CLAMPING AT "
                   << nod_pt->x(0) << " "
                   << nod_pt->x(1) << " "
                   << b << " "
                   << zeta[0] << " " 
                   << std::endl;

        // Get the Jacobian from the GeomObject parametrising the boundary
        GeomObject* geom_obj_pt=curviline_pt->geom_object_pt();
        DenseMatrix<double> drdzeta(nzeta,2);
        geom_obj_pt->dposition(zeta,drdzeta);
        RankThreeTensor<double> d2rdzeta2(nzeta,2,2);
        geom_obj_pt->d2position(zeta,d2rdzeta2);
        double dtdzeta=sqrt(drdzeta(0,0)*drdzeta(0,0)+drdzeta(0,1)*drdzeta(0,1));

        // Foppl von Karman elements only have 6 Hermite deflection degrees
        // of freedom. They are {w ; w,x ; w,y ; w,xx ; w,xy ; w,yy}
        // or their rotated counterparts {w ; w,n ; w,t ; w,nn ; w,nt ; w,tt}
        unsigned nw_type=nw_type_at_each_node();
        for (unsigned k_type=0;k_type<nw_type;k_type++)
         {
          double value=0.0;
          switch (k_type)
           {
           case 0:
            // f itself
            value=boundary_values_pt->f(zeta[0]);
            break;
            
           case 1:
            // df/dn: Normal derivative
            // Calls broken virtual function; dies if not implemented
            value=boundary_values_pt->dfdn(zeta[0]);
            break;
            
           case 2:
            // df/dt, including Jacobian
            value=boundary_values_pt->dfdzeta(zeta[0])/dtdzeta;
            break;
            
           case 3:
            // d^2f/dn^2: Second normal derivative shouldn't be set!
            break;
            
           case 4:
            // d^2f/dndt mixed second normal derivative: 
            value=boundary_values_pt->d2fdndzeta(zeta[0])/dtdzeta;
            break;
            
           case 5:
            // d^2f/dt^2  including Jacobian (twice; chain rule yourself to death)
            value=
             (boundary_values_pt->d2fdzeta2(zeta[0])*dtdzeta*dtdzeta-
              boundary_values_pt->dfdzeta(zeta[0])*
              (drdzeta(0,0)*d2rdzeta2(0,0,0)+drdzeta(0,1)*d2rdzeta2(0,0,1)))/
             (pow(dtdzeta,4));
            
            break;

           default:
              std::stringstream error_message;
              error_message << "never get here" << std::endl;
              throw OomphLibError(error_message.str(),
                                  OOMPH_CURRENT_FUNCTION,
                                  OOMPH_EXCEPTION_LOCATION);
           }

          // Skip second normal derivative
          if (k_type!=3)
           {            
            // Pin and set the value
            nod_pt->pin(first_nodal_type_index + k_type);
            nod_pt->set_value(first_nodal_type_index + k_type, value);
           }
         }
        
       }
     }

   }
 }


  

 ///////////////////////////////////////////////////////////////////////////
 ///////////////////////////////////////////////////////////////////////////
 ///////////////////////////////////////////////////////////////////////////
 
 


} // namespace oomph


#endif

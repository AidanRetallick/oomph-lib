// Non--inline functions for the kirchhoff_plate_bending equations
#include "koiter_steigmann_equations.h"


namespace oomph
{
  /// Unsigned that holds the number if displacments the element has
  const unsigned KoiterSteigmannEquations::Number_of_displacements = 3;

  /// Default Nu value -- the Poisson ratio
  const double KoiterSteigmannEquations::Default_Nu_Value = 0.5;

  /// Default Eta value -- a placeholder for dimensional coefficients
  const double KoiterSteigmannEquations::Default_Eta_Value = 1.0;

  /// Default Mu value -- the coefficient of damping (may be poorly scaled)
  const double KoiterSteigmannEquations::Default_Mu_Value = 1.0;

  /// Default plate thickness
  const double KoiterSteigmannEquations::Default_Thickness_Value = 0.01;

  /// Output to mathematica
  template<typename T>
  void output_mathematica(std::ostream& outfile,
                          const RankThreeTensor<T>& tensor)
  {
    // Loop over outer index
    const unsigned n1 = tensor.nindex1();
    for (unsigned i = 0; i < n1; ++i)
    {
      const unsigned n2 = tensor.nindex2();
      outfile << (i == 0 ? "{" : "");
      // Loop over mid index
      for (unsigned j = 0; j < n2; ++j)
      {
        const unsigned n3 = tensor.nindex3();
        outfile << (j == 0 ? "{" : "");
        // Loop over inner index
        for (unsigned k = 0; k < n3; ++k)
        {
          outfile << (k == 0 ? "{" : "") << tensor(i, j, k)
                  << (k == n3 - 1 ? "}" : ",");
        }
        outfile << (j == n2 - 1 ? "}" : ",");
      }
      outfile << (i == n1 - 1 ? "}" : ",");
    }
  }

  //======================================================================
  void KoiterSteigmannEquations::
  fill_in_generic_residual_contribution_koiter_steigmann(
    Vector<double>& residuals,
    DenseMatrix<double>& jacobian,
    const unsigned& flag)
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

    // [zdec] INITIALISE ALL TO 0.0
    // Basis funtion (assumed to be the same for all displacements)
    Shape psi_n(n_u_node, n_u_nodal_type);
    Shape psi_i(n_u_internal_type);
    DShape dpsi_n_dxi(n_u_node, n_u_nodal_type, n_deriv);
    DShape dpsi_i_dxi(n_u_internal_type, n_deriv);
    DShape d2psi_n_dxi2(n_u_node, n_u_nodal_type, n_2deriv);
    DShape d2psi_i_dxi2(n_u_internal_type, n_2deriv);

    // Test funcitons (assumed to be the same for all displacements)
    Shape test_n(n_u_node, n_u_nodal_type);
    Shape test_i(n_u_internal_type);
    DShape dtest_n_dxi(n_u_node, n_u_nodal_type, n_deriv);
    DShape dtest_i_dxi(n_u_internal_type, n_deriv);
    DShape d2test_n_dxi2(n_u_node, n_u_nodal_type, n_2deriv);
    DShape d2test_i_dxi2(n_u_internal_type, n_2deriv);

    // Initialise the matrices and Vectors for the basic interpolated variables
    // Calculate values of unknown
    Vector<double> interpolated_u(n_displacements, 0.0);
    DenseMatrix<double> interpolated_dudxi(n_displacements, n_deriv, 0.0);
    DenseMatrix<double> interpolated_d2udxi2(n_displacements, n_2deriv, 0.0);
    Vector<double> interpolated_dudt(n_displacements, 0.0);

    // Allocate and initialise to zero
    Vector<double> interpolated_x(dim);
    Vector<double> s(dim);

    // Set up containers for the normal vector and pressure vector
    Vector<double> n_vector(n_displacements);
    Vector<double> pressure(n_displacements);
    DenseMatrix<double> d_pressure_dn(n_displacements, n_displacements);
    DenseMatrix<double> d_pressure_dr(n_displacements, n_displacements);
    RankThreeTensor<double> d_pressure_d_grad_u(
      n_displacements, n_displacements, dim);

    // Set up matrices for the metric, strain, stress and curvature tensors,
    // along with matrices to hold the two tangent vectors and tensions
    DenseMatrix<double> g_matrix(dim, dim);
    DenseMatrix<double> e_tensor(dim, dim);
    DenseMatrix<double> e0_tensor(dim, dim);
    DenseMatrix<double> s_tensor(dim, dim);
    DenseMatrix<double> b_tensor(dim, dim);
    DenseMatrix<double> g_vectors(n_displacements, dim);
    DenseMatrix<double> t_vectors(n_displacements, dim);

    // Set up Tensors to hold the christoffel symbols (of the deformed
    // configuration) and to hold the 3 moment tensors
    RankThreeTensor<double> gamma_tensor(dim, dim, dim);
    RankThreeTensor<double> m_tensors(n_displacements, dim, dim);

    // [zdec] David didn't think this was good enough clearly
    // IF we are constructing the Jacobian find someway of avoiding the overhead
    // HERE
    // Set up containers for the d2_interpolated_u_dx_dunknown and second deriv.
    Vector<double> du_dui_unknown(6);
    // Set up container for the d_normal vector
    RankThreeTensor<double> d_n_vector_du_unknown(
      n_displacements, n_displacements, dim);
    // Set up tensors for the derivatives of metric, strain, stress and
    // curvature tensors, along  with tensor to hold the derivatives of the two
    // tangent vectors and tensions
    RankFourTensor<double> d_e_tensor_du_unknown(2, 2, 3, 2);
    RankFourTensor<double> d_s_tensor_du_unknown(2, 2, 3, 3);
    RankFourTensor<double> d_g_tensor_du_unknown(2, 2, 3, 2);
    RankFourTensor<double> d_b_tensor_du_unknown(2, 2, 3, 5);
    RankFourTensor<double> d_t_vectors_du_unknown(3, 2, 3, 6);
    // Set up Tensors to hold the derivatives of christoffel symbols (of the
    // deformed configuration) and to hold the derivatives of the moment tensors
    RankFiveTensor<double> d_gamma_tensor_du_unknown(2, 2, 2, 3, 5);
    RankFiveTensor<double> d_m_tensors_du_unknown(3, 2, 2, 3, 5);

    // Storage for the value of the damping coefficient
    const double mu = get_mu();

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

      s[0] = this->integral_pt()->knot(ipt, 0);
      s[1] = this->integral_pt()->knot(ipt, 1);
      this->interpolated_x(s, interpolated_x);
      // Call the derivatives of the shape and test functions for the unknown
      double J = d2basis_and_d2test_u_eulerian_koiter_steigmann(s,
                                                                psi_n,
                                                                psi_i,
                                                                dpsi_n_dxi,
                                                                dpsi_i_dxi,
                                                                d2psi_n_dxi2,
                                                                d2psi_i_dxi2,
                                                                test_n,
                                                                test_i,
                                                                dtest_n_dxi,
                                                                dtest_i_dxi,
                                                                d2test_n_dxi2,
                                                                d2test_i_dxi2);

      // Premultiply the weights and the Jacobian
      double W = w * J;


      // Interpolated displacements and derivatives
      //-----------------------------------------
      for (unsigned i_field = 0; i_field < n_displacements; i_field++)
      {
        // Reset interpolated values
        interpolated_u[i_field] = 0.0;
        interpolated_dudt[i_field] = 0.0;
        for(unsigned alpha = 0; alpha < n_deriv; alpha++)
        {
          interpolated_dudxi(i_field, alpha) = 0.0;
        }
        for(unsigned alphabeta = 0; alphabeta < n_2deriv; alphabeta++)
        {
          interpolated_d2udxi2(i_field, alphabeta) = 0.0;
        }

        // Loop over nodes
        for (unsigned j_node = 0; j_node < n_u_node; j_node++)
        {
          // --------- Damping setup
          // By default, we have no history values (no damping)
          unsigned n_time = 0;
          // Turn on damping if this field requires it AND we are not doing a
          // steady solve
          TimeStepper* timestepper_pt =
            this->node_pt(j_node)->time_stepper_pt();
          bool damping =
            Ui_is_damped[i_field] && !(timestepper_pt->is_steady());
          if (damping)
          {
            n_time = timestepper_pt->ntstorage();
          }
          // ---------

          // Loop over types at each node and add nodal contributions
          for (unsigned k_type = 0; k_type < n_u_nodal_type; k_type++)
          {
            // Get the nodal value of the unknown
            double u_value =
              get_u_i_value_at_node_of_type(i_field, j_node, k_type);
            interpolated_u[i_field] += u_value * psi_n(j_node, k_type);
            // Loop over directions
            for (unsigned alpha = 0; alpha < n_deriv; alpha++)
            {
              interpolated_dudxi(i_field, alpha) +=
                u_value * dpsi_n_dxi(j_node, k_type, alpha);
            }
            for (unsigned alphabeta = 0; alphabeta < n_2deriv; alphabeta++)
            {
              interpolated_d2udxi2(i_field, alphabeta) +=
                u_value * d2psi_n_dxi2(j_node, k_type, alphabeta);
            }
            // Loop over the history values (if damping, then n_time>0) and
            // add history contribution to nodal contribution to time derivative
            double nodal_dudt_value = 0.0;
            for (unsigned t_time = 0; t_time < n_time; t_time++)
            {
              nodal_dudt_value +=
                get_u_i_value_at_node_of_type(t_time, i_field, j_node, k_type) *
                timestepper_pt->weight(1, t_time);
            }
            interpolated_dudt[i_field] +=
              nodal_dudt_value * psi_n(j_node, k_type);
          }
        }

        // Loop over internal dof contributions
        // By default, we have no history values (no damping)
        unsigned n_time = 0;
        // Turn on damping if this field requires it AND we are not doing a
        // steady solve
        TimeStepper* timestepper_pt =
          this->u_i_internal_data_pt(i_field)->time_stepper_pt();
        bool damping = Ui_is_damped[i_field] && !(timestepper_pt->is_steady());
        if (damping)
        {
          n_time = timestepper_pt->ntstorage();
        }
        for (unsigned k_type = 0; k_type < n_u_internal_type; k_type++)
        {
          double u_value = get_u_i_internal_value_of_type(i_field, k_type);
          interpolated_u[i_field] += u_value * psi_i(k_type);
          // Loop over directions
          for (unsigned alpha = 0; alpha < n_deriv; alpha++)
          {
            interpolated_dudxi(i_field, alpha) +=
              u_value * dpsi_i_dxi(k_type, alpha);
          }
          for (unsigned alphabeta = 0; alphabeta < n_2deriv; alphabeta++)
          {
            interpolated_d2udxi2(i_field, alphabeta) +=
              u_value * d2psi_i_dxi2(k_type, alphabeta);
          }
          // Loop over the history values (if damping, then n_time>0) and
          // add history contribution to nodal contribution to time derivative
          double internal_dudt_value = 0.0;
          for (unsigned t_time = 0; t_time < n_time; t_time++)
          {
            internal_dudt_value +=
              get_u_i_internal_value_of_type(t_time, i_field, k_type) *
              timestepper_pt->weight(1, t_time);
          }
          interpolated_dudt[i_field] += internal_dudt_value * psi_i(k_type);
        }
      }

      // Initialise
      // Get tangent vectors
      fill_in_tangent_vectors(interpolated_dudxi, g_vectors);

      // Get metric tensor
      fill_in_metric_tensor(g_vectors, g_matrix);

      // Construct the tensor and vector quantities
      // Get unit normal
      fill_in_unit_normal(g_vectors, n_vector);

      // Get prestrain tensor [zdec] do we need the ipt argument?
      get_prestrain(ipt, interpolated_x, e0_tensor);

      // Get strain tensor
      // fill_in_strain_tensor(interpolated_dudxi, e_tensor);
      fill_in_strain_tensor(interpolated_dudxi, e0_tensor, e_tensor);

      // Get strain tensor
      fill_in_stress_tensor(
        interpolated_x, interpolated_u, e_tensor, g_matrix, s_tensor);

      // Get Curvature tensors
      fill_in_curvature_tensor(n_vector, interpolated_d2udxi2, b_tensor);

      // Get Moment tensors
      fill_in_moment_tensor(
        interpolated_d2udxi2, n_vector, b_tensor, m_tensors);

      // Get second christoffel tensor
      fill_in_second_christoffel_tensor(
        interpolated_dudxi, interpolated_d2udxi2, gamma_tensor);

      // Get total tensions
      fill_in_total_tension(s_tensor,
                            gamma_tensor,
                            m_tensors,
                            interpolated_dudxi,
                            interpolated_d2udxi2,
                            t_vectors);

      // Get pressure function
      get_pressure(ipt,
                   interpolated_x,
                   interpolated_u,
                   interpolated_dudxi,
                   n_vector,
                   pressure);

      if (flag)
      {
        // The derivative of the strain tensor
        fill_in_d_strain_tensor_du_unknown(interpolated_dudxi,
                                           d_e_tensor_du_unknown);

        fill_in_d_g_tensor_du_unknown(interpolated_dudxi,
                                      d_g_tensor_du_unknown);

        fill_in_d_stress_tensor_du_unknown(interpolated_x,
                                           interpolated_u,
                                           e_tensor,
                                           g_matrix,
                                           s_tensor,
                                           d_e_tensor_du_unknown,
                                           d_s_tensor_du_unknown);

        // Fill in normal
        fill_in_d_unit_normal_du_unknown(
          g_vectors, g_matrix, d_g_tensor_du_unknown, d_n_vector_du_unknown);

        // Fill in curvature
        fill_in_d_curvature_tensor_du_unknown(n_vector,
                                              interpolated_d2udxi2,
                                              d_n_vector_du_unknown,
                                              d_b_tensor_du_unknown);

        // Fill in curvature
        fill_in_d_second_christoffel_tensor_dui_unknown(
          interpolated_dudxi, interpolated_d2udxi2, d_gamma_tensor_du_unknown);

        // Fill in moment
        fill_in_d_moment_tensor_du_unknown(interpolated_d2udxi2,
                                           n_vector,
                                           b_tensor,
                                           d_n_vector_du_unknown,
                                           d_b_tensor_du_unknown,
                                           d_m_tensors_du_unknown);

        // Total tension
        d_fill_in_total_tension_du_unknown(s_tensor,
                                           gamma_tensor,
                                           m_tensors,
                                           interpolated_dudxi,
                                           interpolated_d2udxi2,
                                           d_s_tensor_du_unknown,
                                           d_gamma_tensor_du_unknown,
                                           d_m_tensors_du_unknown,
                                           d_t_vectors_du_unknown);

        get_d_pressure_dn(ipt,
                          interpolated_x,
                          interpolated_u,
                          interpolated_dudxi,
                          n_vector,
                          d_pressure_dn);

        get_d_pressure_dr(ipt,
                          interpolated_x,
                          interpolated_u,
                          interpolated_dudxi,
                          n_vector,
                          d_pressure_dr);

        get_d_pressure_d_grad_u(ipt,
                                interpolated_x,
                                interpolated_u,
                                interpolated_dudxi,
                                n_vector,
                                d_pressure_d_grad_u);
      }

      // Loop over the displacement fields which are all treated equally
      for (unsigned i_field = 0; i_field < n_displacements; i_field++)
      {
        // Loop over the nodal test functions by visiting each node and then
        // each type at that node
        for (unsigned j_node = 0; j_node < n_u_node; j_node++)
        {
          for (unsigned k_type = 0; k_type < n_u_nodal_type; k_type++)
          {
            // Get the local equation
            unsigned u_ik_index = nodal_index_of_u_i(i_field, k_type);
            local_eqn = this->nodal_local_eqn(j_node, u_ik_index);
            // IF it's not a boundary condition
            if (local_eqn >= 0)
            {
              // Add body force/pressure term here as well as damping
              residuals[local_eqn] +=
                W * test_n(j_node, k_type) *
                (mu * interpolated_dudt[i_field] - pressure[i_field]);
              for (unsigned alpha = 0; alpha < n_deriv; alpha++)
              {
                // w_{,\alpha\beta} \delta \kappa_\alpha\beta
                residuals[local_eqn] += t_vectors(i_field, alpha) *
		  dtest_n_dxi(j_node, k_type, alpha) * W;

                for (unsigned beta = 0; beta < n_deriv; beta++)
                {
                  // w_{,\alpha\beta} \delta \kappa_\alpha\beta
                  residuals[local_eqn] +=
                    m_tensors(i_field, alpha, beta) *
                    d2test_n_dxi2(j_node, k_type, alpha + beta) * W;
                }
              }
              // Calculate the jacobian
              //-----------------------
              if (flag)
              {
                // Loop over the test functions again
                for (unsigned j_node2 = 0; j_node2 < n_u_node; j_node2++)
                {
                  // Loop over position dofs
                  for (unsigned k_type2 = 0; k_type2 < n_u_nodal_type;
                       k_type2++)
                  {
                    // Fill in the derivatives of basis
                    // The derivatives of u, D u(x,y) and D2 u(x,y)(x,y)  wrt to
                    // the (i_field2,j_node2,i_field2)th unknown
                    du_dui_unknown[0] = psi_n(j_node2, k_type2);
                    // Loop over inplane coordinates
                    for (unsigned alpha = 0; alpha < 2; ++alpha)
                    {
                      // Fill in first two derivatives of basis
                      du_dui_unknown[1 + alpha] =
                        dpsi_n_dxi(j_node2, k_type2, alpha);
                      for (unsigned beta = 0; beta < 2; ++beta)
                      {
                        // Fill in second three derivatives of basis
                        du_dui_unknown[3 + alpha + beta] =
                          d2psi_n_dxi2(j_node2, k_type2, alpha + beta);
                      }
                    }

                    // Loop over displacement dofs
                    for (unsigned i_field2 = 0; i_field2 < n_displacements;
                         i_field2++)
                    {
                      // Get the local equation
                      unsigned u_index2 = nodal_index_of_u_i(i_field2, k_type2);
                      local_unknown = this->nodal_local_eqn(j_node2, u_index2);
                      // If at a non-zero degree of freedom add in the entry
                      if (local_unknown >= 0)
                      {
                        // Damping scope
			if (i_field2 == i_field)
                        {
                          // If damping is on at this node, add contribution
                          TimeStepper* timestepper_pt =
                            this->node_pt(j_node2)->time_stepper_pt();
                          bool damped = Ui_is_damped[i_field2] &&
			    !(timestepper_pt->is_steady());
                          if (damped)
                          {
                            jacobian(local_eqn, local_unknown) +=
                              mu * psi_n(j_node2, k_type2) *
                              timestepper_pt->weight(1, 0) *
                              test_n(j_node, k_type) * W;
                          }
                        }
                        // Add body force/pressure term here
                        jacobian(local_eqn, local_unknown) -=
                          d_pressure_dr(i_field, i_field2) * du_dui_unknown[0] *
                          test_n(j_node, k_type) * W;
                        // Loop over first derivatives of basis
                        for (unsigned mu = 0; mu < 2; ++mu)
                        {
                          jacobian(local_eqn, local_unknown) -=
                            d_pressure_d_grad_u(i_field, i_field2, mu) *
                            du_dui_unknown[1 + mu] * test_n(j_node, k_type) * W;
                        }
                        for (unsigned j = 0; j < n_displacements; ++j)
                        {
                          // Loop over first derivatives of basis
                          for (unsigned mu = 0; mu < 2; ++mu)
                          {
                            jacobian(local_eqn, local_unknown) -=
                              d_pressure_dn(i_field, j) *
                              d_n_vector_du_unknown(j, i_field2, mu) *
                              du_dui_unknown[1 + mu] * test_n(j_node, k_type) *
                              W;
                          }
                        }
                        // Loop over inplane coordinates
                        for (unsigned alpha = 0; alpha < 2; ++alpha)
                        {
                          // Tension may depend on u through the stress
                          jacobian(local_eqn, local_unknown) +=
                            d_t_vectors_du_unknown(
                              i_field, alpha, i_field2, 0) *
                            du_dui_unknown[0] *
                            dtest_n_dxi(j_node, k_type, alpha) * W;
                          // Loop over first and second derivatives of basis
                          for (unsigned m2 = 0; m2 < 5; ++m2)
                          {
                            // w_{,\alpha\beta} \delta \kappa_\alpha\beta
                            jacobian(local_eqn, local_unknown) +=
                              d_t_vectors_du_unknown(
                                i_field, alpha, i_field2, 1 + m2) *
                              du_dui_unknown[1 + m2] *
                              dtest_n_dxi(j_node, k_type, alpha) * W;
                            // Loop over inplane coordinates
                            for (unsigned beta = 0; beta < 2; ++beta)
                            {
                              // w_{,\alpha\beta} \delta \kappa_\alpha\beta
                              jacobian(local_eqn, local_unknown) +=
                                d_m_tensors_du_unknown(
                                  i_field, alpha, beta, i_field2, m2) *
                                du_dui_unknown[1 + m2] *
                                d2test_n_dxi2(j_node, k_type, alpha + beta) * W;
                            }
                          }
                        }
                      }
                    }
                  }
                }
                // Loop over the internal test functions
                for (unsigned k_type2 = 0; k_type2 < n_u_internal_type;
                     k_type2++)
                {
                  // Fill in the tensor derivatives
                  // The derivatives of u, D u(x,y) and D2 u(x,y)(x,y)  wrt to
                  // the (i_field2,j_node2,i_field2)th unknown
                  du_dui_unknown[0] = psi_i(k_type2);
                  // Loop over inplane coordinates
                  for (unsigned alpha = 0; alpha < 2; ++alpha)
                  {
                    // Fill in first two derivatives of basis
                    du_dui_unknown[1 + alpha] = dpsi_i_dxi(k_type2, alpha);
                    // Loop over inplane coordinates
                    for (unsigned beta = 0; beta < 2; ++beta)
                    {
                      // Fill in second three derivatives of basis
                      du_dui_unknown[3 + alpha + beta] =
                        d2psi_i_dxi2(k_type2, alpha + beta);
                    }
                  }
                  // Loop over unknown displacement components i_field2
                  for (unsigned i_field2 = 0; i_field2 < n_displacements;
                       i_field2++)
                  {
                    local_unknown = internal_local_eqn(i_field2, k_type2);
                    // If at a non-zero degree of freedom add in the entry
                    if (local_unknown >= 0)
                    {
                      // Damping scope
		      if (i_field2 == i_field)
                      {
                        // If damping is on at this node, add contribution
                        TimeStepper* timestepper_pt =
                          this->u_i_internal_data_pt(i_field)
			      ->time_stepper_pt();
                        bool damped = Ui_is_damped[i_field2] &&
			  !(timestepper_pt->is_steady());
                        if (damped)
                        {
                          jacobian(local_eqn, local_unknown) +=
                            mu * psi_i(k_type2) * timestepper_pt->weight(1, 0) *
                            test_n(j_node, k_type) * W;
                        }
                      }
                      // Add body force/pressure term here
                      jacobian(local_eqn, local_unknown) -=
                        d_pressure_dr(i_field, i_field2) * du_dui_unknown[0] *
                        test_n(j_node, k_type) * W;
                      // Loop over first derivatives of basis
                      for (unsigned mu = 0; mu < 2; ++mu)
                      {
                        jacobian(local_eqn, local_unknown) -=
                          d_pressure_d_grad_u(i_field, i_field2, mu) *
                          du_dui_unknown[1 + mu] * test_n(j_node, k_type) * W;
                      }
                      // Loop over displacement components
                      for (unsigned j = 0; j < n_displacements; ++j)
                      {
                        // Loop over first derivatives of basis
                        for (unsigned mu = 0; mu < 2; ++mu)
                        {
                          jacobian(local_eqn, local_unknown) -=
                            d_pressure_dn(i_field, j) *
                            d_n_vector_du_unknown(j, i_field2, mu) *
                            du_dui_unknown[1 + mu] * test_n(j_node, k_type) * W;
                        }
                      }
                      // Loop over inplane coordinates
                      for (unsigned alpha = 0; alpha < 2; ++alpha)
                      {
                        // Tension may depend on u through the stress
                        jacobian(local_eqn, local_unknown) +=
                          d_t_vectors_du_unknown(i_field, alpha, i_field2, 0) *
                          du_dui_unknown[0] *
                          dtest_n_dxi(j_node, k_type, alpha) * W;
                        // Loop over
                        for (unsigned m2 = 0; m2 < 5; ++m2)
                        {
                          // w_{,\alpha\beta} \delta \kappa_\alpha\beta
                          jacobian(local_eqn, local_unknown) +=
                            d_t_vectors_du_unknown(
                              i_field, alpha, i_field2, 1 + m2) *
                            du_dui_unknown[1 + m2] *
                            dtest_n_dxi(j_node, k_type, alpha) * W;

                          // Loop over inplane coordinates
                          for (unsigned beta = 0; beta < 2; ++beta)
                          {
                            // w_{,\alpha\beta} \delta \kappa_\alpha\beta
                            jacobian(local_eqn, local_unknown) +=
                              d_m_tensors_du_unknown(
                                i_field, alpha, beta, i_field2, m2) *
                              du_dui_unknown[1 + m2] *
                              d2test_n_dxi2(j_node, k_type, alpha + beta) * W;
                          }
                        }
                      }
                    }
                  }
                }
              } // End of flag
            }
          }
        }

        // Loop over the internal test functions
        for (unsigned k_type = 0; k_type < n_u_internal_type; k_type++)
        {
          // Get the local equation
          local_eqn = internal_local_eqn(i_field, k_type);
          // IF it's not a boundary condition
          if (local_eqn >= 0)
          {
            // Add body force/pressure term here as well as damping
            residuals[local_eqn] +=
              W * test_i(k_type) *
              (mu * interpolated_dudt[i_field] - pressure[i_field]);

            for (unsigned alpha = 0; alpha < 2; ++alpha)
            {
              // w_{,\alpha\beta} \delta \kappa_\alpha\beta
              residuals[local_eqn] +=
                t_vectors(i_field, alpha) * dtest_i_dxi(k_type, alpha) * W;

              for (unsigned beta = 0; beta < 2; ++beta)
              {
                // w_{,\alpha\beta} \delta \kappa_\alpha\beta
                residuals[local_eqn] += m_tensors(i_field, alpha, beta) *
		  d2test_i_dxi2(k_type, alpha + beta) * W;
              }
            }
            // Calculate the jacobian
            //-----------------------
            if (flag)
            {
              // Loop over the test functions again
              for (unsigned j_node2 = 0; j_node2 < n_u_node; j_node2++)
              {
                // Loop over position dofs
                for (unsigned k_type2 = 0; k_type2 < n_u_nodal_type; k_type2++)
                {
                  // Fill in the tensor derivatives
                  // The derivatives of u, D u(x,y) and D2 u(x,y)(x,y)  wrt to
                  // the (i_field2,j_node2,i_field2)th unknown
                  du_dui_unknown[0] = psi_n(j_node2, k_type2);
                  // Loop over inplane coordinates
                  for (unsigned alpha = 0; alpha < 2; ++alpha)
                  {
                    // Fill in first two derivatives of basis
                    du_dui_unknown[1 + alpha] =
                      dpsi_n_dxi(j_node2, k_type2, alpha);
                    // Loop over inplane coordinates
                    for (unsigned beta = 0; beta < 2; ++beta)
                    {
                      // Fill in second three derivatives of basis
                      du_dui_unknown[3 + alpha + beta] =
                        d2psi_n_dxi2(j_node2, k_type2, alpha + beta);
                    }
                  }
                  // Loop over displacement dofs
                  for (unsigned i_field2 = 0; i_field2 < n_displacements;
                       i_field2++)
                  {
                    // Get the local equation
                    unsigned u_index2 = nodal_index_of_u_i(i_field2, k_type2);
                    local_unknown = this->nodal_local_eqn(j_node2, u_index2);
                    // If at a non-zero degree of freedom add in the entry
                    if (local_unknown >= 0)
		    {
		      // Damping scope
		      if (i_field2 == i_field)
		      {
			// If damping is on at this node, add contribution
			TimeStepper* timestepper_pt =
			  this->node_pt(j_node2)->time_stepper_pt();
			bool damped = Ui_is_damped[i_field2] &&
			  !(timestepper_pt->is_steady());
			if (damped)
			{
			  jacobian(local_eqn, local_unknown) +=
			    mu * psi_n(j_node2, k_type2) *
			    timestepper_pt->weight(1, 0) * test_i(k_type) * W;
			}
		      }
		      // Add body force/pressure term here
		      jacobian(local_eqn, local_unknown) -=
			d_pressure_dr(i_field, i_field2) * du_dui_unknown[0] *
			test_i(k_type) * W;
		      for (unsigned j = 0; j < n_displacements; ++j)
		      {
			// Loop over first derivatives of basis
			for (unsigned mu = 0; mu < 2; ++mu)
			{
			  jacobian(local_eqn, local_unknown) -=
			    d_pressure_dn(i_field, j) *
			    d_n_vector_du_unknown(j, i_field2, mu) *
			    du_dui_unknown[1 + mu] * test_i(k_type) * W;
			}
		      }
		      // Loop over inplane coordinates
		      for (unsigned alpha = 0; alpha < 2; ++alpha)
		      {
			// Tension may depend on u through the stress
			jacobian(local_eqn, local_unknown) +=
			  d_t_vectors_du_unknown(i_field, alpha, i_field2, 0) *
			  du_dui_unknown[0] * dtest_i_dxi(k_type, alpha) * W;
			// Loop over first and second derivatives of basis
			for (unsigned m2 = 0; m2 < 5; ++m2)
			{
			  // w_{,\alpha\beta} \delta \kappa_\alpha\beta
			  jacobian(local_eqn, local_unknown) +=
			    d_t_vectors_du_unknown(
			      i_field, alpha, i_field2, 1 + m2) *
			    du_dui_unknown[1 + m2] *
			    dtest_i_dxi(k_type, alpha) * W;

                          // Loop over inplane coordinates
                          for (unsigned beta = 0; beta < 2; ++beta)
                          {
                            // w_{,\alpha\beta} \delta \kappa_\alpha\beta
                            jacobian(local_eqn, local_unknown) +=
                              d_m_tensors_du_unknown(
                                i_field, alpha, beta, i_field2, m2) *
                              du_dui_unknown[1 + m2] *
                              d2test_i_dxi2(k_type, alpha + beta) * W;
                          }
                        }
                      }
                    }
                  }
                }
              }
              // Loop over the internal test functions
              for (unsigned k_type2 = 0; k_type2 < n_u_internal_type; k_type2++)
              {
                // Fill in the derivatives of basis
                // The derivatives of u, D u(x,y) and D2 u(x,y)(x,y)  wrt to
                // the (i_field2,j_node2,i_field2)th unknown
                du_dui_unknown[0] = psi_i(k_type2);
                // Loop over inplane coordinates
                for (unsigned alpha = 0; alpha < 2; ++alpha)
                {
                  // Fill in first two derivatives of basis
                  du_dui_unknown[1 + alpha] = dpsi_i_dxi(k_type2, alpha);
                  // Loop over inplane coordinates
                  for (unsigned beta = 0; beta < 2; ++beta)
                  {
                    // Fill in second three derivatives of basis
                    du_dui_unknown[3 + alpha + beta] =
                      d2psi_i_dxi2(k_type2, alpha + beta);
                  }
                }
                // Loop over displacement dofs
                for (unsigned i_field2 = 0; i_field2 < n_displacements;
                     i_field2++)
                {
                  local_unknown = internal_local_eqn(i_field2, k_type2);
                  // If at a non-zero degree of freedom add in the entry
                  if (local_unknown >= 0)
                  {
                    // Damping scope
		    if (i_field2 == i_field)
                    {
                      // If damping is on at this node, add contribution
                      TimeStepper* timestepper_pt =
                        this->u_i_internal_data_pt(i_field)->time_stepper_pt();
                      bool damped = Ui_is_damped[i_field2] &&
			!(timestepper_pt->is_steady());
                      if (damped)
                      {
                        jacobian(local_eqn, local_unknown) +=
                          mu * psi_i(k_type2) * timestepper_pt->weight(1, 0) *
                          test_i(k_type) * W;
                      }
                    }
                    // Add body force/pressure term here
                    jacobian(local_eqn, local_unknown) -=
                      d_pressure_dr(i_field, i_field2) * du_dui_unknown[0] *
                      test_i(k_type) * W;
                    // Loop over first derivatives of basis
                    for (unsigned mu = 0; mu < 2; ++mu)
                    {
                      jacobian(local_eqn, local_unknown) -=
                        d_pressure_d_grad_u(i_field, i_field2, mu) *
                        du_dui_unknown[1 + mu] * test_i(k_type) * W;
                    }
                    for (unsigned j = 0; j < n_displacements; ++j)
                    {
                      for (unsigned mu = 0; mu < 2; ++mu)
                      {
                        jacobian(local_eqn, local_unknown) -=
                          d_pressure_dn(i_field, j) *
                          d_n_vector_du_unknown(j, i_field2, mu) *
                          du_dui_unknown[1 + mu] * test_i(k_type) * W;
                      }
                    }
                    for (unsigned alpha = 0; alpha < 2; ++alpha)
                    {
                      // Tension may depend on u through the stress
                      jacobian(local_eqn, local_unknown) +=
                        d_t_vectors_du_unknown(i_field, alpha, i_field2, 0) *
                        du_dui_unknown[0] * dtest_i_dxi(k_type, alpha) * W;
                      for (unsigned m2 = 0; m2 < 5; ++m2)
                      {
                        // w_{,\alpha\beta} \delta \kappa_\alpha\beta
                        //   alpha,i_field2)*dtest_i_dxi(l,k,alpha)*W;
                        jacobian(local_eqn, local_unknown) +=
                          d_t_vectors_du_unknown(
                            i_field, alpha, i_field2, 1 + m2) *
                          du_dui_unknown[1 + m2] * dtest_i_dxi(k_type, alpha) *
                          W;

                        for (unsigned beta = 0; beta < 2; ++beta)
                        {
                          // w_{,\alpha\beta} \delta \kappa_\alpha\beta
                          jacobian(local_eqn, local_unknown) +=
                            d_m_tensors_du_unknown(
                              i_field, alpha, beta, i_field2, m2) *
                            du_dui_unknown[1 + m2] *
                            d2test_i_dxi2(k_type, alpha + beta) * W;
                        }
                      }
                    }
                  }
                }
              }
	    } // End of flag
	  }
	}
      } // End loop over the displacements [i_field]
    } // End of loop over integration points
  } // End of fill in generic residual contribution


    //======================================================================
    /// Self-test:  Return 0 for OK
    //======================================================================
  unsigned KoiterSteigmannEquations::self_test()
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
  /// Output function:
  ///
  ///   x,y,u   or    x,y,z,u
  ///
  /// nplot points in each coordinate direction
  //======================================================================
  void KoiterSteigmannEquations::output(std::ostream & outfile,
					const unsigned& nplot)
  {
    // Dimension of the element
    const unsigned dim = this->dim();

    // Vector of local coordinates
    Vector<double> s(dim, 0.0), x(dim, 0.0);

    // Tecplot header info
    outfile << this->tecplot_zone_string(nplot);

    // Loop over plot points
    unsigned num_plot_points = this->nplot_points(nplot);
    // Vector<double> r(3);

    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      this->get_s_plot(iplot, nplot, s);
      Vector<Vector<double>> u(Number_of_displacements,
			       Vector<double>(6, 0.0));
      interpolated_koiter_steigmann_disp(s, u);

      // Get x position as Vector
      this->interpolated_x(s, x);

      for (unsigned i = 0; i < dim; i++)
      {
	outfile << x[i] << " ";
      }

      // Loop for variables
      for (unsigned i = 0; i < Number_of_displacements; i++)
      {
	for (unsigned j = 0; j < 6; j++)
	{
	  outfile << u[i][j] << " ";
	}
      }

      outfile << std::endl;
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    this->write_tecplot_zone_footer(outfile, nplot);
  }


  //======================================================================
  /// Output function:
  ///
  ///   x,y,u,stress
  ///
  /// nplot points in each coordinate direction
  //======================================================================
  void KoiterSteigmannEquations::output_stress(std::ostream & outfile,
					       const unsigned& nplot)
  {
    // Dimension of the element
    const unsigned dim = this->dim();

    // Vector of local and global coordinates
    Vector<double> s(dim, 0.0);
    Vector<double> x(dim, 0.0);

    // Vector of displaced position cordinates
    Vector<double> r(dim, 0.0);
    DenseMatrix<double> drdx(Number_of_displacements, dim, 0.0);

    // Storage for prestrain, strain, metric and stress
    DenseMatrix<double> prestrain(dim, dim, 0.0);
    DenseMatrix<double> strain(dim, dim, 0.0);
    DenseMatrix<double> metric(dim, dim, 0.0);
    DenseMatrix<double> stress(dim, dim, 0.0);

    // Tecplot header info
    outfile << this->tecplot_zone_string(nplot);

    // Loop over plot points
    unsigned num_plot_points = this->nplot_points(nplot);
    // Vector<double> r(3);

    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      this->get_s_plot(iplot, nplot, s);

      // Get the displacement
      Vector<Vector<double>> u(Number_of_displacements,
			       Vector<double>(6, 0.0));
      interpolated_koiter_steigmann_disp(s, u);

      // Extract the gradient and displaced position from the displacements
      DenseMatrix<double> dudx(Number_of_displacements, dim, 0.0);
      for(unsigned i = 0; i < Number_of_displacements; i++)
      {
	r[i] = x[i] + u[i][0];
	for(unsigned j = 0; j < dim; j++)
	{
	  dudx(i,j) = u[i][1+j];
	  drdx(i,j) = dudx(i,j) + (double)(i==j);
	}
      }

      // Get x position as Vector
      this->interpolated_x(s, x);

      // Get the strain and stress
      fill_in_metric_tensor(drdx, metric);
      get_prestrain(iplot, x, prestrain);
      fill_in_strain_tensor(dudx, prestrain, strain);
      fill_in_stress_tensor(x, r, strain, metric, stress);

      for (unsigned i = 0; i < dim; i++)
      {
	outfile << x[i] << " ";
      }

      // Loop for displacements
      for (unsigned i = 0; i < Number_of_displacements; i++)
      {
	for (unsigned j = 0; j < 6; j++)
	{
	  outfile << u[i][j] << " ";
	}
      }

      // Loop for stress
      for (unsigned i = 0; i < 2; i++)
      {
	for (unsigned j = 0; j < 2; j++)
	{
	  outfile << stress(i,j) << " ";
	}
      }

      outfile << std::endl;
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    this->write_tecplot_zone_footer(outfile, nplot);
  }


  //======================================================================
  /// C-style output function:
  ///
  ///   x,y,u   or    x,y,z,u
  ///
  /// nplot points in each coordinate direction
  //======================================================================
  void KoiterSteigmannEquations::output(FILE * file_pt, const unsigned& nplot)
  {
    // Store the dimension of the element
    const unsigned dim = this->dim();

    // Vector of local coordinates
    Vector<double> s(dim), x(dim);
    ;

    // Tecplot header info
    fprintf(file_pt, "%s", this->tecplot_zone_string(nplot).c_str());

    // Loop over plot points

    unsigned num_plot_points = this->nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      this->get_s_plot(iplot, nplot, s);

      // Get x position as Vector
      this->interpolated_x(s, x);

      for (unsigned i = 0; i < dim; i++)
      {
	fprintf(file_pt, "%g ", x[i]);
      }
      Vector<Vector<double>> u(Number_of_displacements,
			       Vector<double>(6, 0.0));
      interpolated_koiter_steigmann_disp(s, u);
      for (unsigned ii = 0; ii < Number_of_displacements; ++ii)
      {
	fprintf(file_pt, "%g \n", u[ii][0]);
      } // interpolated_u_poisson(s));
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
  void KoiterSteigmannEquations::output_fct(
    std::ostream & outfile,
    const unsigned& nplot,
    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
  {
    // Store the dimension of the element
    const unsigned dim = this->dim();

    // Vector of local coordinates
    Vector<double> s(dim);

    // Vector for coordintes
    Vector<double> x(dim);

    // Tecplot header info
    outfile << this->tecplot_zone_string(nplot);

    // Exact solution Vector (here a scalar)
    Vector<double> exact_soln(this->required_nvalue(0), 0.0);

    // Loop over plot points
    unsigned num_plot_points = this->nplot_points(nplot);
    for (unsigned iplot = 0; iplot < num_plot_points; iplot++)
    {
      // Get local coordinates of plot point
      this->get_s_plot(iplot, nplot, s);

      // Get x position as Vector
      this->interpolated_x(s, x);

      // Get exact solution at this point
      (*exact_soln_pt)(x, exact_soln);

      // Output x,y,...,u_exact
      for (unsigned i = 0; i < dim; i++)
      {
	outfile << x[i] << " ";
      }
      // Loop over variables
      for (unsigned j = 0; j < this->required_nvalue(0); j++)
      {
	outfile << exact_soln[j] << " ";
      }
      outfile << std::endl;
    }

    // Write tecplot footer (e.g. FE connectivity lists)
    this->write_tecplot_zone_footer(outfile, nplot);
  }

  //======================================================================
  /// Validate against exact solution
  ///
  /// Solution is provided via function pointer.
  /// Plot error at a given number of plot points.
  ///
  //======================================================================
  void KoiterSteigmannEquations::compute_error(
    std::ostream & outfile,
    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
    Vector<double> & error,
    Vector<double> & norm)
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


    // Vector of local coordinates
    Vector<double> s(dim);

    // Vector for coordintes
    Vector<double> x(dim);

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
      for (unsigned i = 0; i < dim; i++)
      {
	s[i] = this->integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      double w = this->integral_pt()->weight(ipt);

      // Get jacobian of mapping
      double J;
      // [zdec] Get this J without all the stupid functions
      // double Jlin = this->J_eulerian1(s);// Nope
      {
	// Find values of c1-shape function
	J = d2basis_u_eulerian_koiter_steigmann(s,
						psi_n,
						psi_i,
						dpsi_n_dxi,
						dpsi_i_dxi,
						d2psi_n_dxi2,
						d2psi_i_dxi2);
      }
      // Premultiply the weights and the Jacobian
      double W = w * J;
      // double Wlin = w*Jlin;

      // Get x position as Vector
      this->interpolated_x(s, x);

      // Get FE function value
      Vector<Vector<double>> u_fe(n_displacements, Vector<double>(6, 0.0));
      interpolated_koiter_steigmann_disp(s, u_fe);

      // Get exact solution at this point
      (*exact_soln_pt)(x, exact_soln);

      // Check if we have defined a new metric
      Vector<double> tmp1(error.size(), 0.0), tmp2(norm.size(), 0.0);
      if (Multiple_error_metric_fct_pt == 0)
      {
	// HERE tidy this up
	// Do nothing! Maybe add a static warning or something
      }
      else
      {
	// Tmp storage
	Vector<double> u_fe_tmp(18, 0.0);
	// Flatpack
	for (unsigned ii = 0; ii < 18; ++ii)
	{
	  u_fe_tmp[ii] = u_fe[ii / 6][ii % 6];
	}
	// Get the metric
	(*Multiple_error_metric_fct_pt)(x, u_fe_tmp, exact_soln, tmp1, tmp2);
      }

      // Add on to the error and norm
      for (unsigned i = 0; i < norm.size(); ++i)
      {
	error[i] += tmp1[i] * W;
      }
      for (unsigned i = 0; i < norm.size(); ++i)
      {
	norm[i] += tmp2[i] * W;
      }
    } // End of loop over integration pts
  }


  //======================================================================
  /// Validate against exact solution
  ///
  /// Solution is provided via function pointer.
  /// Plot error at a given number of plot points.
  ///
  //======================================================================
  void KoiterSteigmannEquations::compute_error(
    std::ostream & outfile,
    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
    double& error,
    double& norm)
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

    // Initialise
    error = 0.0;
    norm = 0.0;

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
      for (unsigned i = 0; i < dim; i++)
      {
	s[i] = this->integral_pt()->knot(ipt, i);
      }

      // Get the integral weight
      double w = this->integral_pt()->weight(ipt);

      // Get jacobian of mapping
      double J;
      // [zdec] Do this without all the basis and test
      {
	// Find values of c1-shape function
	J = d2basis_u_eulerian_koiter_steigmann(s,
						psi_n,
						psi_i,
						dpsi_n_dxi,
						dpsi_i_dxi,
						d2psi_n_dxi2,
						d2psi_i_dxi2);
      }

      // Premultiply the weights and the Jacobian
      double W = w * J;
      // double Wlin = w*Jlin;

      // Get x position as Vector
      this->interpolated_x(s, x);

      // Get FE function value
      Vector<Vector<double>> u_fe(n_displacements, Vector<double>(6, 0.0));
      interpolated_koiter_steigmann_disp(s, u_fe);

      // Get exact solution at this point
      (*exact_soln_pt)(x, exact_soln);

      // Check if we have defined a new metric
      double tmp1 = 0.0, tmp2 = 0.0;
      if (Error_metric_fct_pt == 0)
      {
	// Output x,y,...,error
	for (unsigned i = 0; i < dim; i++)
	{
	  outfile << x[i] << " ";
	}
	for (unsigned ii = 0; ii < 6 * n_displacements; ii++)
	{
	  outfile << exact_soln[ii % 6] << " "
		  << exact_soln[ii] - u_fe[ii / 6][ii % 6] << " ";
	}
	outfile << std::endl;

	// Loop over variables
	for (unsigned ii = 0; ii < n_displacements; ii++)
	{
	  // Add to error and norm
	  // Size of solution r_exact . r_exact
	  tmp1 += (exact_soln[ii * 6] * exact_soln[ii * 6]);
	  // Default norm is just (r - r_exact) . (r - r_exact)
	  tmp2 += (pow(exact_soln[ii * 6] - u_fe[ii][0], 2));
	}
      }
      // Use a user defined error metric
      else
      {
	// Tmp storage
	Vector<double> u_fe_tmp(18, 0.0);
	// Flatpack
	for (unsigned ii = 0; ii < 18; ++ii)
	{
	  u_fe_tmp[ii] = u_fe[ii / 6][ii % 6];
	}
	// Get the metric
	(*Error_metric_fct_pt)(x, u_fe_tmp, exact_soln, tmp2, tmp1);
      }
      norm += tmp1 * W;
      error += tmp2 * W;
    } // End of loop over integration pts
      // Reset to zero
  }

} // namespace oomph

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
// oomph-lib headers


// test2

#include "c1_curved_elements.h"

namespace oomph
{
  namespace MyC1CurvedElements
  {
    // Declaration of a function protypes  (explicit instantiation)
    /// Fill in matrix to convert the monomials to the basic shape functions
    template<>
    void oomph::MyC1CurvedElements::BernadouElementBasis<
      3>::monomial_to_basic_matrix(oomph::DenseMatrix<double>& m) const;

    /// \short Fill in inverse matrix to convert the monomials to the basic
    /// shape functions (BROKEN)
    template<>
    void oomph::MyC1CurvedElements::BernadouElementBasis<
      3>::inverse_monomial_to_basic_matrix(DenseDoubleMatrix& ib2l) const;

    /// Fill in the (precomputed) full basic polynomials
    template<>
    void BernadouElementBasis<3>::full_basic_polynomials(
      const Vector<double>& s, Shape& phi) const;

    /// Fill in the derivatives of the (precomputed) full basic polynomials
    template<>
    void BernadouElementBasis<3>::dfull_basic_polynomials(
      const Vector<double>& s, DShape& phi) const;

    /// Fill in second derivatives of the (precomputed) full basic polynomials
    template<>
    void BernadouElementBasis<3>::d2full_basic_polynomials(
      const Vector<double>& s, DShape& phi) const;

    /// Fill in matrix to convert the monomials to the basic shape functions
    template<>
    void oomph::MyC1CurvedElements::BernadouElementBasis<
      5>::monomial_to_basic_matrix(oomph::DenseMatrix<double>& m) const;

    /// Fill in inverse matrix to convert the monomials to the basic shape
    /// functions
    template<>
    void oomph::MyC1CurvedElements::BernadouElementBasis<
      5>::inverse_monomial_to_basic_matrix(DenseDoubleMatrix& ib2l) const;

    /// Fill in the (precomputed) full basic polynomials
    template<>
    void BernadouElementBasis<5>::full_basic_polynomials(
      const Vector<double>& s, Shape& phi) const;

    /// Fill in the derivatives of the (precomputed) full basic polynomials
    template<>
    void BernadouElementBasis<5>::dfull_basic_polynomials(
      const Vector<double>& s, DShape& phi) const;

    /// Fill in second derivatives of the (precomputed) full basic polynomials
    template<>
    void BernadouElementBasis<5>::d2full_basic_polynomials(
      const Vector<double>& s, DShape& phi) const;

    /// Parametric function describing curved edge as a function of local coords
    template<>
    void BernadouElementBasis<3>::f_k(const Vector<double>& s,
                                      Vector<double>& fk) const
    {
      // Some shorthands
      Vector<double> a1(2, 0.0), a2(2, 0.0), b1(2, 0.0), b2(2, 0.0);
      // Fill in
      A1(a1);
      A2(a2);
      B1(b1);
      B2(b2);
      // Fill in the vector
      for (unsigned i = 0; i < 2; ++i)
      {
        // The affine part
        fk[i] = vertices[2][i] - a1[i] * s[0] -
                b2[i] * s[1]
                // The additional curved part
                + 0.5 * s[0] * s[1] *
                    ((2 * (a1[i] - b2[i]) - (a2[i] - b1[i])) * (s[1] - s[0]) +
                     a2[i] + b1[i]);
      }
    }

    /// Parametric function describing curved edge as a function of local coords
    template<>
    void BernadouElementBasis<5>::f_k(const Vector<double>& s,
                                      Vector<double>& fk) const
    {
      // Some shorthands
      Vector<double> a1(2, 0.0), a2(2, 0.0), b1(2, 0.0), b2(2, 0.0), d1(2, 0.0),
        d2(2, 0.0);
      // Fill in
      A1(a1);
      A2(a2);
      B1(b1);
      B2(b2);
      D1(d1), D2(d2);
      // Fill in the vector
      for (unsigned i = 0; i < 2; ++i)
      {
        // Some shorthands
        double beta0 = b2[i] - a1[i] + a2[i];
        double beta1 = b2[i] - a1[i] + a2[i] + 0.5 * d1[i];
        double beta2 =
          9 * (a1[i] - b2[i]) - 5 * a2[i] + 4 * b1[i] - d1[i] + 0.5 * d2[i];
        double beta3 =
          -6 * (a1[i] - b2[i]) + 3 * a2[i] - 3 * b1[i] + 0.5 * (d1[i] - d2[i]);

        double beta0t = -b2[i] + a1[i] + b1[i];
        double beta1t = -b2[i] + a1[i] + b1[i] + 0.5 * d2[i];
        double beta2t =
          -9 * (a1[i] - b2[i]) - 5 * b1[i] + 4 * a2[i] - d2[i] + 0.5 * d1[i];
        double beta3t =
          6 * (a1[i] - b2[i]) + 3 * b1[i] - 3 * a2[i] + 0.5 * (d2[i] - d1[i]);

        // The affine part
        fk[i] = vertices[2][i] - a1[i] * s[0] -
                b2[i] * s[1]
                // The additional curved part
                + 0.5 * s[0] * s[1] *
                    (beta3 * pow(s[1], 3) + beta2 * pow(s[1], 2) +
                     beta1 * s[1] + beta0 + beta3t * pow(s[0], 3) +
                     beta2t * pow(s[0], 2) + beta1t * s[0] + beta0t);
      }
    }

    /// Get the physical coordinate
    template<unsigned BOUNDARY_ORDER>
    void BernadouElementBasis<BOUNDARY_ORDER>::coordinate_x(
      const Vector<double>& s, Vector<double>& fk) const
    {
      Vector<double> s_basic(s);
      permute_shape(s_basic);
      f_k(s_basic, fk);
    }

    /// The approximated boundary polynomial at point along arclength s2
    template<unsigned BOUNDARY_ORDER>
    void BernadouElementBasis<BOUNDARY_ORDER>::psi_h(const double& s1,
                                                     Vector<double>& psih) const
    {
// Check the construction of the elements is complete
#ifdef PARANOID
      if (Curved_edge == none)
      {
        throw OomphLibError("The element has not been upgraded yet. Did \
  you forget to set upe the Curved_edge?",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
      // F_k at a point along s0 = 1 -s1
      Vector<double> s(2);
      s[0] = 1 - s1;
      s[1] = s1;
      f_k(s, psih);
    }

    /// The Jacobian of the mapping from basic to global
    template<>
    void BernadouElementBasis<3>::get_basic_jacobian(
      const Vector<double> s, DenseMatrix<double>& jacobian) const
    {
      // Initialise the mapping
      // DenseMatrix<double> jacobian(2,2,0.0);
      Vector<double> a1(2, 0.0), a2(2, 0.0), b1(2, 0.0), b2(2, 0.0);
      // Fill in the vector
      A1(a1);
      A2(a2);
      B1(b1);
      B2(b2);
      // Fill in the vector
      for (unsigned i = 0; i < 2; ++i)
      {
        // Some shorthands
        //  const double a1(A1(i)), a2(A2(i)), b1(B1(i)), b2(B2(i));

        // The affine part
        jacobian(i, 0) =
          -a1[i]
          // The additional curved part
          + 0.5 * s[1] *
              ((2 * (a1[i] - b2[i]) - (a2[i] - b1[i])) * s[1] + a2[i] + b1[i]) -
          s[0] * s[1] * ((2 * (a1[i] - b2[i]) - (a2[i] - b1[i])));

        // The affine part
        jacobian(i, 1) = -b2[i]
                         // The additional curved part
                         + 0.5 * s[0] *
                             (-(2 * (a1[i] - b2[i]) - (a2[i] - b1[i])) * s[0] +
                              a2[i] + b1[i]) +
                         s[0] * s[1] * (2 * (a1[i] - b2[i]) - (a2[i] - b1[i]));
      }
    }

    /// Get the Jacobian of the mapping from basic to global
    template<>
    void BernadouElementBasis<5>::get_basic_jacobian(
      const Vector<double> s, DenseMatrix<double>& jacobian) const
    {
      // Some shorthands
      Vector<double> a1(2, 0.0), a2(2, 0.0), b1(2, 0.0), b2(2, 0.0), d1(2, 0.0),
        d2(2, 0.0);
      // Fill in
      A1(a1);
      A2(a2);
      B1(b1);
      B2(b2);
      D1(d1), D2(d2);
      // Fill in the vector
      for (unsigned i = 0; i < 2; ++i)
      {
        // Some shorthands
        double beta0 = b2[i] - a1[i] + a2[i];
        double beta1 = b2[i] - a1[i] + a2[i] + 0.5 * d1[i];
        double beta2 =
          9 * (a1[i] - b2[i]) - 5 * a2[i] + 4 * b1[i] - d1[i] + 0.5 * d2[i];
        double beta3 =
          -6 * (a1[i] - b2[i]) + 3 * a2[i] - 3 * b1[i] + 0.5 * (d1[i] - d2[i]);

        double beta0t = -b2[i] + a1[i] + b1[i];
        double beta1t = -b2[i] + a1[i] + b1[i] + 0.5 * d2[i];
        double beta2t =
          -9 * (a1[i] - b2[i]) - 5 * b1[i] + 4 * a2[i] - d2[i] + 0.5 * d1[i];
        double beta3t =
          6 * (a1[i] - b2[i]) + 3 * b1[i] - 3 * a2[i] + 0.5 * (d2[i] - d1[i]);

        // The affine part
        jacobian(i, 0) =
          -a1[i]
          // The additional curved part
          + 0.5 * s[1] *
              (beta3 * pow(s[1], 3) + beta2 * pow(s[1], 2) + beta1 * s[1] +
               beta0 + 4 * beta3t * pow(s[0], 3) + 3 * beta2t * pow(s[0], 2) +
               2 * beta1t * s[0] + beta0t);
        // The affine part
        jacobian(i, 1) =
          -b2[i]
          // The additional curved part
          + 0.5 * s[0] *
              (4 * beta3 * pow(s[1], 3) + 3 * beta2 * pow(s[1], 2) +
               2 * beta1 * s[1] + beta0 + beta3t * pow(s[0], 3) +
               beta2t * pow(s[0], 2) + beta1t * s[0] + beta0t);
      }
    }

    /// The Hessian of the global coordinate (of the vector mapping) - like a
    /// second order Jacobian.
    ///        d^2 x_i
    /// or:   ----------    (rank 3) with x the global coordinate and s the
    /// basic.
    //       d s_i ds_j
    template<>
    void BernadouElementBasis<3>::get_basic_hessian(
      const Vector<double>& s, RankThreeTensor<double>& hessian) const
    {
      // DenseMatrix<double> jacobian(2,2,0.0);
      Vector<double> a1(2, 0.0), a2(2, 0.0), b1(2, 0.0), b2(2, 0.0);
      // Fill in the vector
      A1(a1);
      A2(a2);
      B1(b1);
      B2(b2);
      // Fill in the vector
      for (unsigned i = 0; i < 2; ++i)
      {
        // This will all be from curved part
        hessian(i, 0, 0) = -s[1] * (2 * (a1[i] - b2[i]) - (a2[i] - b1[i]));

        hessian(i, 0, 1) = 0.5 * ((2 * (a1[i] - b2[i]) - (a2[i] - b1[i])) *
                                    (2 * s[1] - 2 * s[0]) +
                                  a2[i] + b1[i]);

        hessian(i, 1, 0) = hessian(i, 0, 1);

        hessian(i, 1, 1) = s[0] * (2 * (a1[i] - b2[i]) - (a2[i] - b1[i]));
      }
    }

    /// The Hessian of the global coordinate (of the vector mapping) - like a
    /// second order Jacobian.
    ///        d^2 x_i
    /// or:   ----------    (rank 3) with x the global coordinate and s the
    /// basic.
    //       d s_i ds_j
    template<>
    void BernadouElementBasis<5>::get_basic_hessian(
      const Vector<double>& s, RankThreeTensor<double>& hessian) const
    {
      // Some shorthands
      Vector<double> a1(2, 0.0), a2(2, 0.0), b1(2, 0.0), b2(2, 0.0), d1(2, 0.0),
        d2(2, 0.0);
      // Fill in
      A1(a1);
      A2(a2);
      B1(b1);
      B2(b2);
      D1(d1), D2(d2);
      // Fill in the vector
      for (unsigned i = 0; i < 2; ++i)
      {
        // Some shorthands
        double beta0 = b2[i] - a1[i] + a2[i];
        double beta1 = b2[i] - a1[i] + a2[i] + 0.5 * d1[i];
        double beta2 =
          9 * (a1[i] - b2[i]) - 5 * a2[i] + 4 * b1[i] - d1[i] + 0.5 * d2[i];
        double beta3 =
          -6 * (a1[i] - b2[i]) + 3 * a2[i] - 3 * b1[i] + 0.5 * (d1[i] - d2[i]);

        double beta0t = -b2[i] + a1[i] + b1[i];
        double beta1t = -b2[i] + a1[i] + b1[i] + 0.5 * d2[i];
        double beta2t =
          -9 * (a1[i] - b2[i]) - 5 * b1[i] + 4 * a2[i] - d2[i] + 0.5 * d1[i];
        double beta3t =
          6 * (a1[i] - b2[i]) + 3 * b1[i] - 3 * a2[i] + 0.5 * (d2[i] - d1[i]);

        hessian(i, 0, 0) =
          +0.5 * s[1] *
          (12 * beta3t * pow(s[0], 2) + 6 * beta2t * s[0] + 2 * beta1t);

        hessian(i, 0, 1) =
          +0.5 * (4 * beta3 * pow(s[1], 3) + 3 * beta2 * pow(s[1], 2) +
                  2 * beta1 * s[1] + beta0 + 4 * beta3t * pow(s[0], 3) +
                  3 * beta2t * pow(s[0], 2) + 2 * beta1t * s[0] + beta0t);

        hessian(i, 1, 0) = hessian(i, 0, 1);

        hessian(i, 1, 1) =
          +0.5 * s[0] *
          (12 * beta3 * pow(s[1], 2) + 6 * beta2 * s[1] + 2 * beta1);
      }
    }

    /// The matrix that transforms the global dofs to the local dofs
    template<unsigned BOUNDARY_ORDER>
    void BernadouElementBasis<BOUNDARY_ORDER>::local_to_global_matrix(
      DenseMatrix<double>& d) const
    {
      // Converts global dofs:
      // w(ai) Dw(xj)(ai) D2 w(xj,xk)(ai) w(ei) with ai in {a1,a2,a3} and xj,xk
      // in {x1,x2} to w(ai) Dw(t[i]j)(ai) D2 w(t[i]j,t[i]j)(ai)
      // D2(t[i+1]2,t[i+2]1) w(ei) with with i a cyclic permuter in base 3 and
      // with t[1]j in {A1,A2}, t[2]j in {B1,B2} and t[2]j in {C1,C2} labelling
      // the tangents two at each node

      // Matrix should be
      // DenseMatrix<double> d(21,21,0.0);

      // Initialise sub matrices
      DenseMatrix<double> d1(2, 2, 0.0), d2(2, 2, 0.0), d3(2, 2, 0.0),
        d4(9, 9, 0.0);

      // Construct sub matrices
      for (unsigned i = 0; i < 2; ++i)
      {
        // To transform to new  local basis at node 0
        d1(i, 0) = A1(i);
        d1(i, 1) = A2(i);
        // To transform to new  local basis at node 1
        d2(i, 0) = B1(i);
        d2(i, 1) = B2(i);
        // To transform to new  local basis at node 2
        d3(i, 0) = -B2(i); // C1=-B2
        d3(i, 1) = -A1(i); // C2=-B1

        //  Matrix for second derivatives
        for (unsigned j = 0; j < 2; ++j)
        {
          // at node 0
          d4(i + j, 0) += A1(i) * A1(j);
          d4(i + j, 1) += A2(i) * A2(j);
          d4(i + j, 6) += B2(i) * B2(j); // cross deriv. of tangents on opp.
                                         // face

          // at node 1
          d4(i + j + 3, 2) += B1(i) * B1(j);
          d4(i + j + 3, 3) += B2(i) * B2(j);
          d4(i + j + 3, 7) +=
            A1(i) * A1(j); // cross deriv. of tangents on opp. face

          // at node 2
          d4(i + j + 6, 5) += A1(i) * A1(j);
          d4(i + j + 6, 4) += B2(i) * B2(j);
          d4(i + j + 6, 8) -=
            A2(i) * B1(j); // cross deriv. of tangents on opp. face
        }
      }

      // Now Construct the full matrix
      // Fill in 3x3 identity
      for (unsigned i = 0; i < 3; ++i)
      {
        // Identity matrix at start
        d(i, i) = 1.0;
      }
      // Fill in n_internal_dof x n_internal_dof identity
      for (unsigned i = 0; i < n_internal_dofs(); ++i)
      {
        // Identity matrix at end
        d(i + 18, i + 18) = 1.0;
      }

      // Fill in 2x2 's
      for (unsigned i = 0; i < 2; ++i)
      {
        for (unsigned j = 0; j < 2; ++j)
        {
          // Fill in d1 ; d2 ; d3
          d(i + 3, j + 3) = d1(i, j);
          d(i + 5, j + 5) = d2(i, j);
          d(i + 7, j + 7) = d3(i, j);
        }
      }

      // Fill in d4 (9 x 9)
      for (unsigned i = 0; i < 9; ++i)
      {
        for (unsigned j = 0; j < 9; ++j)
        {
          d(i + 9, j + 9) = d4(i, j);
        }
      }
    }

    /// Get 5th order, 1D hermite shape functions
    /// Dofs: w(0) w(1) w'(0) w'(1) w''(0) w''(1)
    template<unsigned BOUNDARY_ORDER>
    void BernadouElementBasis<BOUNDARY_ORDER>::hermite_shape_1d_5(
      const double& s, Shape& psi) const
    {
      // These can be determined by a simple matrix inversion
      psi[0] = pow(1 - s, 3) * (1 + 3 * s + 6 * s * s);
      psi[1] = pow(s, 3) * (10 - 15 * s + 6 * s * s);
      psi[2] = pow(1 - s, 3) * s * (1 + 3 * s);
      psi[3] = (1 - s) * pow(s, 3) * (3 * s - 4);
      psi[4] = 0.5 * pow(1 - s, 3) * pow(s, 2);
      psi[5] = 0.5 * pow(1 - s, 2) * pow(s, 3);
    }

    /// Get 5th order, 1D hermite shape functions
    /// Dofs: w(0) w(1) w'(0) w'(1) w''(0) w''(1)
    template<unsigned BOUNDARY_ORDER>
    void BernadouElementBasis<BOUNDARY_ORDER>::hermite_shape_1d_5(
      const s_basic_node& s, Shape& psi) const
    {
      // These can be determined by a simple matrix inversion
      switch (s)
      {
        // We don't need this at every point so instead use a switch
        case one_quarter:
          psi[0] = 459. / 512.;
          psi[1] = 53. / 512.;
          psi[2] = 189. / 1024.;
          psi[3] = -39. / 1024.;
          psi[4] = 27. / 2048.;
          psi[5] = 9 / 2048.;
          break;
        case one_half:
          psi[0] = 1. / 2.;
          psi[1] = 1. / 2.;
          psi[2] = 5. / 32.;
          psi[3] = -5. / 32.;
          psi[4] = 1. / 64.;
          psi[5] = 1. / 64.;
          break;
        case three_quarters:
          psi[0] = 53. / 512.;
          psi[1] = 459. / 512.;
          psi[2] = 39. / 1024.;
          psi[3] = -189. / 1024.;
          psi[4] = 9. / 2048.;
          psi[5] = 27. / 2048.;
          break;
      }
    }

    /// Get the derivatives of 5th order 1d shape
    template<unsigned BOUNDARY_ORDER>
    void BernadouElementBasis<BOUNDARY_ORDER>::d_hermite_shape_1d_5(
      const double& s, DShape& dpsi) const
    {
      // These can be determined by a simple matrix inversion
      dpsi(0, 0) = -30 * pow(1 - s, 2) * pow(s, 2);
      dpsi(1, 0) = 30 * pow(1 - s, 2) * pow(s, 2);
      dpsi(2, 0) = pow(1 - s, 2) * (1 - 3 * s) * (1 + 5 * s);
      dpsi(3, 0) = pow(s, 2) * (3 * s - 2) * (6 - 5 * s);
      dpsi(4, 0) = 0.5 * pow(1 - s, 2) * s * (2 - 5 * s);
      dpsi(5, 0) = 0.5 * (1 - s) * s * s * (3 - 5 * s);
    }

    /// Get the derivatives of 5th order 1d shape at the basic nodes
    template<unsigned BOUNDARY_ORDER>
    void BernadouElementBasis<BOUNDARY_ORDER>::d_hermite_shape_1d_5(
      const s_basic_node& s, DShape& dpsi) const
    {
      // These can be determined by a simple matrix inversion
      switch (s)
      {
        // We don't need this at every point so instead use a switch
        case one_quarter:
          dpsi(0, 0) = -135. / 128.;
          dpsi(1, 0) = 135. / 128.;
          dpsi(2, 0) = 81. / 256.;
          dpsi(3, 0) = -95. / 256.;
          dpsi(4, 0) = 27. / 512.;
          dpsi(5, 0) = 21. / 512.;
          break;
        case one_half:
          dpsi(0, 0) = -15. / 8.;
          dpsi(1, 0) = 15. / 8.;
          dpsi(2, 0) = -7. / 16.;
          dpsi(3, 0) = -7. / 16.;
          dpsi(4, 0) = -1. / 32.;
          dpsi(5, 0) = 1. / 32.;
          break;
        case three_quarters:
          dpsi(0, 0) = -135. / 128.;
          dpsi(1, 0) = 135. / 128.;
          dpsi(2, 0) = -95. / 256.;
          dpsi(3, 0) = 81. / 256.;
          dpsi(4, 0) = -21. / 512.;
          dpsi(5, 0) = -27. / 512.;
          break;
      }
    }

    /// Get 3rd order, 1D hermite shape functions
    /// Dofs: w(0) w(1) w'(0) w'(1)
    template<unsigned BOUNDARY_ORDER>
    void BernadouElementBasis<BOUNDARY_ORDER>::hermite_shape_1d_3(
      const double& s, Shape& psi) const
    {
      // These can be determined by a simple matrix inversion
      psi[0] = pow(s - 1, 2) * (1 + 2 * s);
      psi[1] = pow(s, 2) * (3 - 2 * s);
      psi[2] = pow(s - 1, 2) * s;
      psi[3] = (s - 1) * pow(s, 2);
    }

    /// Get 3rd order, 1D hermite shape functions
    /// Dofs: w(0) w(1) w'(0) w'(1)
    template<unsigned BOUNDARY_ORDER>
    void BernadouElementBasis<BOUNDARY_ORDER>::hermite_shape_1d_3(
      const s_basic_node& s, Shape& psi) const
    {
      // These can be determined by a simple matrix inversion
      switch (s)
      {
        // We don't need this at every point so instead use a switch
        case one_quarter:
          psi[0] = 27. / 32.;
          psi[1] = 5. / 32.;
          psi[2] = 9. / 64.;
          psi[3] = -3. / 64.;
          break;
        case one_half:
          psi[0] = 1. / 2.;
          psi[1] = 1. / 2.;
          psi[2] = 1. / 8.;
          psi[3] = -1. / 8.;
          break;
        case three_quarters:
          psi[0] = 5. / 32.;
          psi[1] = 27. / 32.;
          psi[2] = 3. / 64.;
          psi[3] = -9. / 64.;
          break;
      }
    }

    /// The representation on trace 1 will be identical regardless of the
    /// boundary interpolation - as it has to be C1 continuous with straight
    /// sided Bell elements Now define the w trace column vectors fi i in
    /// {1,2,3} -these vectors are effectively the 'basis' functions on the
    /// trace.
    template<unsigned BOUNDARY_ORDER>
    Vector<double> BernadouElementBasis<BOUNDARY_ORDER>::f_1(
      const double& s0) const
    {
      // Get the p5 1d basis polynomials
      Shape psi_5(6);
      hermite_shape_1d_5(s0, psi_5);

      // Now construct f1 (affine mapping on this edge)
      Vector<double> f1(n_basis_functions(), 0.0);
      f1[0] = psi_5[1];
      f1[2] = psi_5[0];
      // First tangential derivatives
      f1[3] = -psi_5[3];
      f1[8] = psi_5[2];
      // Second tangential derivatives
      f1[9] = psi_5[5];
      f1[14] = psi_5[4];

      return f1;
    }

    /// The representation on trace 1 will be identical regardless of the
    /// boundary interpolation - as it has to be C1 continuous with straight
    /// sided Bell elements. First (tangent) derivative of trace
    template<unsigned BOUNDARY_ORDER>
    Vector<double> BernadouElementBasis<BOUNDARY_ORDER>::df_1_ds(
      const double& s0) const
    {
      // Get the p5 1d basis polynomials
      DShape dpsi_5(6, 1);
      d_hermite_shape_1d_5(s0, dpsi_5);

      // Now construct f1 (affine mapping on this edge)
      Vector<double> df1(n_basis_functions(), 0.0);
      df1[0] = dpsi_5(1, 0);
      df1[2] = dpsi_5(0, 0);
      // First tangential derivatives
      df1[3] = -dpsi_5(3, 0);
      df1[8] = dpsi_5(2, 0);
      // Second tangential derivatives
      df1[9] = dpsi_5(5, 0);
      df1[14] = dpsi_5(4, 0);

      return df1;
    }

    /// The representation on trace 2 will be identical regardless of the
    /// boundary interpolation - as it has to be C1 continuous with straight
    /// sided Bell elements. First (tangent) derivative of trace
    template<unsigned BOUNDARY_ORDER>
    Vector<double> BernadouElementBasis<BOUNDARY_ORDER>::f_2(
      const double& s1) const
    {
      // Get the p5 1d basis polynomials
      Shape psi_5(6);
      hermite_shape_1d_5(s1, psi_5);

      // Now construct f2 (affine mapping on this edge)
      Vector<double> f2(n_basis_functions(), 0.0);
      f2[1] = psi_5[1];
      f2[2] = psi_5[0];
      // First tangential derivatives
      f2[6] = -psi_5[3];
      f2[7] = psi_5[2];
      // Second tangential derivatives
      f2[12] = psi_5[5];
      f2[13] = psi_5[4];

      return f2;
    }

    /// The representation on trace 2 will be identical regardless of the
    /// boundary interpolation - as it has to be C1 continuous with straight
    /// sided Bell elements.
    template<unsigned BOUNDARY_ORDER>
    Vector<double> BernadouElementBasis<BOUNDARY_ORDER>::df_2_ds(
      const double& s1) const
    {
      // Get the p5 1d basis polynomials
      DShape dpsi_5(6, 1);
      d_hermite_shape_1d_5(s1, dpsi_5);

      // Now construct f2 (affine mapping on this edge)
      Vector<double> f2(n_basis_functions(), 0.0);
      f2[1] = dpsi_5(1, 0);
      f2[2] = dpsi_5(0, 0);
      // First tangential derivatives
      f2[6] = -dpsi_5(3, 0);
      f2[7] = dpsi_5(2, 0);
      // Second tangential derivatives
      f2[12] = dpsi_5(5, 0);
      f2[13] = dpsi_5(4, 0);

      return f2;
    }

    /// This is the curved edge trace - and will therefore depend on boundary
    /// representation.
    template<>
    Vector<double> BernadouElementBasis<3>::f_3(const double& s0) const
    {
      // Get the p5 1d basis polynomials
      Shape psi_5(6);
      hermite_shape_1d_5(s0, psi_5);

      // Now construct f1
      Vector<double> f3(21, 0.0);
      f3[0] = psi_5[1];
      f3[1] = psi_5[0];

      f3[4] = -psi_5[3];
      f3[5] = psi_5[2];
      // Now do second derivatives
      f3[10] = psi_5[5];
      f3[11] = psi_5[4];

      // Extra depedance on first tangent derivatives (Non affine mapping)
      f3[3] -= (6 * a_tilde_1() + 2 * a_tildetilde_1()) * psi_5[5];
      f3[4] -= (6 * a_tilde_2() + 2 * a_tildetilde_2() + 4) * psi_5[5];
      f3[5] -= (6 * b_tilde_1() - 2 * b_tildetilde_1() + 4) * psi_5[4];
      f3[6] -= (6 * b_tilde_2() - 2 * b_tildetilde_2()) * psi_5[4];

      return f3;
    }

    /// This is the curved edge trace - and will therefore depend on boundary
    /// representation.
    template<>
    Vector<double> BernadouElementBasis<5>::f_3(const double& s1) const
    {
      // Get the p5 1d basis polynomials
      Shape psi_5(6);
      hermite_shape_1d_5(s1, psi_5);
      Vector<double> f3(28, 0.0);
      f3[0] = psi_5[1];
      f3[1] = psi_5[0];

      f3[4] = -psi_5[3];
      f3[5] = psi_5[2];
      // Now do second derivatives
      f3[10] = psi_5[5];
      f3[11] = psi_5[4];

      // Extra depedance on first tangent derivatives (Non affine mapping)
      f3[3] = a_utilde_1() * psi_5[5];
      f3[4] += a_utilde_2() * psi_5[5];
      f3[5] += b_utilde_1() * psi_5[4];
      f3[6] = b_utilde_2() * psi_5[4];

      // Now construct f1
      return f3;
    }

    /// This is the curved edge trace - and will therefore depend on boundary
    /// representation.
    /// Tangential derivative of basis.
    template<>
    Vector<double> BernadouElementBasis<3>::df_3_ds(const double& s0) const
    {
      // Get the p5 1d basis polynomials
      DShape dpsi_5(6, 1);
      d_hermite_shape_1d_5(s0, dpsi_5);

      // Now construct f1
      Vector<double> f3(21, 0.0);
      f3[0] = dpsi_5(1, 0);
      f3[1] = dpsi_5(0, 0);

      f3[4] = -dpsi_5(3, 0);
      f3[5] = dpsi_5(2, 0);
      // Now do second derivatives
      f3[10] = dpsi_5(5, 0);
      f3[11] = dpsi_5(4, 0);

      // Extra depedance on first tangent derivatives (Non affine mapping)
      f3[3] -= (6 * a_tilde_1() + 2 * a_tildetilde_1()) * dpsi_5(5, 0);
      f3[4] -= (6 * a_tilde_2() + 2 * a_tildetilde_2() + 4) * dpsi_5(5, 0);
      f3[5] -= (6 * b_tilde_1() - 2 * b_tildetilde_1() + 4) * dpsi_5(6, 0);
      f3[6] -= (6 * b_tilde_2() - 2 * b_tildetilde_2()) * dpsi_5(6, 0);

      return f3;
    }

    /// This is the curved edge trace - and will therefore depend on boundary
    /// representation.
    /// Tangential derivative of basis.
    template<>
    Vector<double> BernadouElementBasis<5>::df_3_ds(const double& s0) const
    {
      // Get the p5 1d basis polynomials
      DShape dpsi_5(6, 1);
      d_hermite_shape_1d_5(s0, dpsi_5);

      // Now construct f1
      Vector<double> f3(28, 0.0);
      f3[0] = dpsi_5(1, 0);
      f3[1] = dpsi_5(0, 0);

      f3[4] = -dpsi_5(3, 0);
      f3[5] = dpsi_5(2, 0);
      // Now do second derivatives
      f3[10] = dpsi_5(5, 0);
      f3[11] = dpsi_5(4, 0);

      // Extra depedance on first tangent derivatives (Non affine mapping)
      f3[3] += a_utilde_1() * dpsi_5(5, 0);
      f3[4] += a_utilde_2() * dpsi_5(5, 0);
      f3[5] += b_utilde_1() * dpsi_5(4, 0);
      f3[6] += b_utilde_2() * dpsi_5(4, 0);
      return f3;
    }

    /// Now define the  w,n trace column vectors gi i in {1,2,3}
    /// The representation on trace 1 will be identical regardless of the
    /// boundary interpolation - as it has to be C1 continuous with straight
    /// sided Bell elements.
    template<unsigned BOUNDARY_ORDER>
    Vector<double> BernadouElementBasis<BOUNDARY_ORDER>::g_1(
      const double& s0) const
    {
      // Get the p5 1d basis polynomials
      Shape psi_3(4);
      hermite_shape_1d_3(s0, psi_3);
      const double atilde1 = a_tilde_1();
      const double atilde2 = a_tilde_2();

      // Now construct f1
      Vector<double> g1(n_basis_functions(), 0.0);
      g1[3] = 0.5 * (eta_2() - 1 - 2 * atilde1) * psi_3[1];
      g1[4] = -atilde2 * psi_3[1];

      g1[7] = 1.0 * psi_3[0];
      g1[8] = -0.5 * (1 + eta_2()) * psi_3[0];

      // Second derivatives
      g1[9] = 0.5 * (atilde1 - eta_2()) * psi_3[3];
      g1[10] = -pow(atilde2, 2) / (2. + 2. * atilde1) * psi_3[3];
      g1[15] = 1.0 / (2 + 2 * atilde1) * psi_3[3];

      // For convenience
      double d2f_denom(c_tilde_2() * c_tildetilde_1() +
                       c_tilde_1() * c_tildetilde_2());

      g1[13] = -c_tilde_1() * c_tildetilde_1() / d2f_denom * psi_3[2];
      g1[14] =
        -(0.5 + 0.5 * eta_2() + c_tilde_2() * c_tildetilde_2() / d2f_denom) *
        psi_3[2];
      g1[17] = -1.0 / d2f_denom * psi_3[2];
      return g1;
    }

    /// Now define the  w,n trace column vectors gi i in {1,2,3}
    /// The representation on trace 1 will be identical regardless of the
    /// boundary interpolation - as it has to be C1 continuous with straight
    /// sided Bell elements.
    template<unsigned BOUNDARY_ORDER>
    Vector<double> BernadouElementBasis<BOUNDARY_ORDER>::g_2(
      const double& s0) const
    {
      // Get the p5 1d basis polynomials
      Shape psi_3(4);
      hermite_shape_1d_3(s0, psi_3);
      const double btilde1 = b_tilde_1();
      const double btilde2 = b_tilde_2();

      // Now construct f1
      Vector<double> g1(n_basis_functions(), 0.0);
      g1[6] = 0.5 * (-eta_1() - 1 - 2 * btilde2) * psi_3[1];
      g1[5] = -btilde1 * psi_3[1];

      g1[8] = 1.0 * psi_3[0];
      g1[7] = -0.5 * (1 - eta_1()) * psi_3[0]; // This might be wrong

      // Second derivatives
      g1[12] = 0.5 * (btilde2 + eta_1()) * psi_3[3];
      g1[11] = -btilde1 * btilde1 / (2 + 2 * btilde2) * psi_3[3];
      g1[16] = 1.0 / (2 + 2 * btilde2) * psi_3[3];

      // For convenience
      double d2f_denom =
        (c_tilde_2() * c_tildetilde_1() + c_tilde_1() * c_tildetilde_2());

      g1[14] = -c_tilde_2() * c_tildetilde_2() / d2f_denom * psi_3[2];
      g1[13] =
        (0.5 * (-1 + eta_1()) - c_tilde_1() * c_tildetilde_1() / d2f_denom) *
        psi_3[2];
      g1[17] = -1.0 / d2f_denom * psi_3[2];
      return g1;
    }

    /// This is the curved edge trace - and will therefore depend on boundary
    /// representation.
    template<>
    Vector<double> BernadouElementBasis<3>::g_3(const double& s0) const
    {
      // Get the p5 1d basis polynomials
      Shape psi_3(4);
      hermite_shape_1d_3(s0, psi_3);
      const double atilde1 = a_tilde_1();
      const double atilde2 = a_tilde_2();
      const double btilde1 = b_tilde_1();
      const double btilde2 = b_tilde_2();

      // Now construct f1
      Vector<double> g3(21, 0.0);
      g3[3] = 1.0 * psi_3[1];
      g3[4] = -0.5 * psi_3[1];

      g3[5] = -0.5 * psi_3[0];
      g3[6] = 1.0 * psi_3[0];

      // Second derivatives
      g3[3] -= (atilde1 + 0.5 * a_tildetilde_1()) * psi_3[3];
      g3[4] -= (0.5 + atilde2 + 0.5 * a_tildetilde_2()) * psi_3[3];

      g3[9] = (1.0 + atilde1) / (2 * atilde2) * psi_3[3];
      g3[10] = 0.5 * (1 + atilde2 / (1 + atilde1)) * psi_3[3];
      g3[15] = -1.0 / (2 * atilde2 * (1 + atilde1)) * psi_3[3];

      g3[5] += (0.5 + btilde1 - 0.5 * b_tildetilde_1()) * psi_3[2];
      g3[6] += (btilde2 - 0.5 * b_tildetilde_2()) * psi_3[2];

      g3[11] = -0.5 * (1 + btilde1 / (1 + btilde2)) * psi_3[2];
      g3[12] = -(1.0 + btilde2) / (2 * btilde1) * psi_3[2];
      g3[16] = 1.0 / (2 * btilde1 * (1 + btilde2)) * psi_3[2];
      return g3;
    }

    /// This is the curved edge trace - and will therefore depend on boundary
    /// representation.
    template<>
    Vector<double> BernadouElementBasis<5>::g_3(const double& s0) const
    {
      // Get the p5 1d basis polynomials
      Shape psi_3(4);
      hermite_shape_1d_3(s0, psi_3);

      // Now construct f1
      const double atilde1 = a_tilde_1();
      const double atilde2 = a_tilde_2();
      const double btilde1 = b_tilde_1();
      const double btilde2 = b_tilde_2();

      Vector<double> g3(28, 0.0);
      g3[3] = 1.0 * psi_3[1];
      g3[4] = -0.5 * psi_3[1];

      g3[5] = -0.5 * psi_3[0];
      g3[6] = 1.0 * psi_3[0];

      // Second derivatives
      g3[3] += (atilde1 + a_utilde_1() / 2.) / 2. * psi_3[3];
      g3[4] += (1 + atilde2 + a_utilde_2() / 2.) / 2. * psi_3[3];

      g3[9] = (1.0 + atilde1) / (2. * atilde2) * psi_3[3];
      g3[10] = 0.5 * (1 + atilde2 / (1 + atilde1)) * psi_3[3];
      g3[15] = -1.0 / (2 * atilde2 * (1 + atilde1)) * psi_3[3];

      g3[5] -= (1 + btilde1 + b_utilde_1() / 2.) / 2. * psi_3[2];
      g3[6] -= (btilde2 + b_utilde_2() / 2.) / 2. * psi_3[2];

      g3[11] = -0.5 * (1 + btilde1 / (1 + btilde2)) * psi_3[2];
      g3[12] = -(1.0 + btilde2) / (2. * btilde1) * psi_3[2];
      g3[16] = 1.0 / (2 * btilde1 * (1 + btilde2)) * psi_3[2];
      return g3;
    }

    /// Get 36 (54) basis monomials for generic p7(9) polynomial on basic
    /// triangle Rename to something more generic and template it

    template<unsigned BOUNDARY_ORDER>
    void BernadouElementBasis<BOUNDARY_ORDER>::full_basis_monomials(
      const Vector<double>& s, Shape& p7) const
    {
      // Initialise the basis monomials
      // Initialise a counter - this is a bit clumsy but does the job
      unsigned counter = 0;
      // This looks complicated but it is only to remain consistent with
      // Bernadou and Boisserie (1994)
      for (int l = basic_basis_order(); l >= 0; --l)
      {
        for (int m = 0; m <= l; ++m)
        {
          // pow(0,0)=1, so this is fine
          p7[counter] = pow(s[0], l - m) * pow(s[1], m);
          ++counter;
        }
      }
    }

    /// Get first derivatives of the 36(54) basis monomials for generic p7(9)
    /// polynomial
    template<unsigned BOUNDARY_ORDER>
    void BernadouElementBasis<BOUNDARY_ORDER>::dfull_basis_monomials(
      const Vector<double>& s, DShape& dp7) const
    {
      // Now construct
      // Initialise the basis monomials
      // Initialise a counter - this is a bit clumsy but does the job
      unsigned counter = 0;
      // This looks complicated but it is only to remain consistent with
      // Bernadou and Boisserie (1994)
      for (int l = basic_basis_order(); l >= 0; --l)
      {
        for (int m = 0; m <= l; ++m)
        {
          // pow potentially undefined for pow(0,0) in future versions?
          // Do the first derivative wrt s0
          if (l - m != 0)
            dp7(counter, 0) = (l - m) * pow(s[0], l - m - 1) * pow(s[1], m);
          else
            dp7(counter, 0) = 0.0;
          // Do the first derivative wrt s0
          if (m != 0) dp7(counter, 1) = (m)*pow(s[0], l - m) * pow(s[1], m - 1);
          else
            dp7(counter, 1) = 0.0;
          ++counter;
        }
      }
    }

    /// Get second derivatives of the 36(55) basis monomials for generic p7(9)
    /// polynomial
    template<unsigned BOUNDARY_ORDER>
    void BernadouElementBasis<BOUNDARY_ORDER>::d2full_basis_monomials(
      const Vector<double>& s, DShape& d2p7) const
    {
      // Now construct
      // Initialise the basis monomials
      // Initialise a counter - this is a bit clumsy but does the job
      unsigned counter = 0;
      // This looks complicated but it is only to remain consistent with
      // Bernadou and Boisserie (1994)
      for (int l = basic_basis_order(); l >= 0; --l)
      {
        for (int m = 0; m <= l; ++m)
        {
          // pow potentially undefined for pow(0,0) in future versions?
          // Do the second derivative wrt s0
          if (l - m != 0 && l - m != 1)
            d2p7(counter, 0) =
              (l - m) * (l - m - 1) * pow(s[0], l - m - 2) * pow(s[1], m);
          else
            d2p7(counter, 0) = 0.0;

          // Do the cross derivative
          if (l - m != 0 && m != 0)
            d2p7(counter, 1) =
              (l - m) * m * pow(s[0], l - m - 1) * pow(s[1], m - 1);
          else
            d2p7(counter, 1) = 0.0;

          // Do the second derivative wrt s1
          if (m != 0 && m != 1)
            d2p7(counter, 2) =
              (m) * (m - 1) * pow(s[0], l - m) * pow(s[1], m - 2);
          else
            d2p7(counter, 2) = 0.0;
          ++counter;
        }
      }
    }


    /// Submatrix M1 (B_1 in Bernadou and Boisserie 1997)
    /// This transforms w(ai) onto basic element (identity matrix)
    template<unsigned BOUNDARY_ORDER>
    void BernadouElementBasis<BOUNDARY_ORDER>::basic_to_local_submatrix_1(
      DenseMatrix<double>& m1) const
    {
      // This should have been passed in as an empty (initialised) 24x3 matrix
      for (unsigned i = 0; i < 3; ++i) m1(i, i) = 1.0; // Identity matrix
    }

    /// Submatrix M2 (B_2 in Bernadou and Boisserie 1997)
    /// This transforms w,j(ai) onto local dofs (which can be converted to basic
    /// dofs)
    template<unsigned BOUNDARY_ORDER>
    void BernadouElementBasis<BOUNDARY_ORDER>::basic_to_local_submatrix_2(
      DenseMatrix<double>& m2) const
    {
      // This should have been passed in as an empty (initialised) 24x6 matrix
      DenseMatrix<double> b2(6, 6, 0.0);

      // Transforms tangential Derivatives at Node 0 to local dofs
      b2(0, 0) = -1; // D w(A1) is -D w (-ex)
      b2(0, 1) = -1;
      b2(1, 1) = +1; // D w (A2) is D w(-ex + ey)

      // Transforms tangential Derivatives at Node 1 to basic dofs
      b2(2, 2) = +1;
      b2(3, 3) = -1; // -D w (B1) is D w (ex - ey)
      b2(3, 2) = -1; // D w (B2) is D w (-ey)

      // Transforms tangential Derivatives at Node 2 to basic dofs
      b2(4, 5) = +1; // D w (C1) is D w (ey)
      b2(5, 4) = +1; // D w (C2) is D w (ex)

      // Fill the Submatrix in using this
      for (unsigned i = 0; i < 6; ++i)
        for (unsigned j = 0; j < 6; ++j) m2(i + 3, j) = b2(i, j);
    }

    /// Submatrix M3 (B_3 in Bernadou and Boisserie 1997)
    /// This transforms w,jk(ai) onto basic element. There will be dependance on
    /// first derivatives due to the (in general) non affine mapping
    /// This will be different between element orders.
    template<>
    void BernadouElementBasis<3>::basic_to_local_submatrix_3(
      DenseMatrix<double>& m3) const
    {
      // This should have been passed in as an empty (initialised) 24x6 matrix
      DenseMatrix<double> b31(6, 9, 0.0);
      DenseMatrix<double> b32(9, 9, 0.0);
      const double atilde1 = a_tilde_1();
      const double atilde2 = a_tilde_2();
      const double btilde1 = b_tilde_1();
      const double btilde2 = b_tilde_2();


      // Node 0
      // The Second tangent derivative in terms of basic dofs
      // D2 w(a0)(A1,A1) transforms to D2 w(Fk(a0))(s0,s0)
      b32(0, 0) = 1.0;

      //  D2 w(Fk(a0))(s0,s1) is more involved
      // First derivatives
      b31(0, 1) = 2 * atilde1 + a_tildetilde_1() / 2.;
      b31(1, 1) = 3 / 2. + 2 * atilde2 + a_tildetilde_2() / 2.;

      b31(0, 2) = -2 * atilde1 - a_tildetilde_1();
      b31(1, 2) = -1.0 - 2 * atilde2 - a_tildetilde_2();
      // Second derivatives
      b32(0, 1) = 1 + (1 + atilde1) / (2 * atilde2);
      b32(1, 1) = atilde2 / (2 + 2 * atilde1);
      b32(6, 1) = -1.0 / (2 * atilde2 * (1 + atilde1));

      b32(0, 2) = 1 + (1 + atilde1) / (atilde2);
      b32(1, 2) = 1 + atilde2 / (1 + atilde1);
      b32(6, 2) = -1.0 / (atilde2 * (1 + atilde1));

      // Node 1
      //  D2 w(Fk(a1))(s1,s1) is more involved
      // First derivatives
      b31(2, 3) = -1.0 - 2 * btilde1 + b_tildetilde_1();
      b31(3, 3) = -2 * btilde2 + b_tildetilde_2();

      b31(2, 4) = 3 / 2. + 2 * btilde1 - b_tildetilde_1() / 2.;
      b31(3, 4) = 2 * btilde2 - b_tildetilde_2() / 2.;
      // Second derivatives
      b32(2, 3) = 1 + btilde1 / (1 + btilde2);
      b32(3, 3) = 1 + (1 + btilde2) / (btilde1);
      b32(7, 3) = -1.0 / (btilde1 * (1 + btilde2));

      b32(2, 4) = btilde1 / (2 + 2 * btilde2);
      b32(3, 4) = 1 + (1 + btilde2) / (2 * btilde1);
      b32(7, 4) = -1.0 / (2 * btilde1 * (1 + btilde2));

      // D2 w(a2)(C1,C1) transforms to D2 w(Fk(a2))(s1,s1)
      b32(3, 5) = 1.0;

      // Node 3
      // D2 w(a1)(A1,A1) transforms to D2 w(Fk(a1))(s0,s0)
      b32(5, 6) = 1.0;

      //  D2 w(Fk(a0))(s0,s1) is more involved
      // First derivatives
      b31(4, 7) = (c_tilde_1() + c_tildetilde_1()) / 2.;
      b31(5, 7) = (c_tilde_2() + c_tildetilde_2()) / 2.;

      b32(4, 7) =
        -(c_tilde_1() * c_tildetilde_1()) /
        (c_tilde_2() * c_tildetilde_1() + c_tilde_1() * c_tildetilde_2());
      b32(5, 7) =
        -(c_tilde_2() * c_tildetilde_2()) /
        (c_tilde_2() * c_tildetilde_1() + c_tilde_1() * c_tildetilde_2());
      b32(8, 7) = -1.0 / (c_tilde_2() * c_tildetilde_1() +
                          c_tilde_1() * c_tildetilde_2());

      // D2 w(a2)(C2,C2) transforms to D2 w(Fk(a2))(s0,s0)
      b32(4, 8) = 1.0;

      // Now construct the submatrix
      for (unsigned j = 0; j < 9; ++j)
      {
        for (unsigned i = 0; i < 6; ++i)
        {
          // Fill in b31
          m3(i + 3, j) = b31(i, j);
          // Fill in b32 (up to 6 so we don't need a seperate loop)
          m3(i + 9, j) = b32(i, j);
        }
        // Carry on the for loop
        for (unsigned i = 6; i < 9; ++i)
        {
          // Fill in b32 (6 to 9)
          m3(i + 9, j) = b32(i, j);
        }
      }
    }

    /// Submatrix M3 (B_3 in Bernadou and Boisserie 1997)
    /// This transforms w,jk(ai) onto basic element. There will be dependance on
    /// first derivatives due to the (in general) non affine mapping
    /// This will be different between element orders.
    template<>
    void BernadouElementBasis<5>::basic_to_local_submatrix_3(
      DenseMatrix<double>& m3) const
    {
      // This should have been passed in as an empty (initialised) 28x6 matrix
      // b31 differs from the b31<3> but b32 stays the same (because the
      // Jacobian differs)
      DenseMatrix<double> b31(6, 9, 0.0);
      DenseMatrix<double> b32(9, 9, 0.0);
      const double atilde1 = a_tilde_1();
      const double atilde2 = a_tilde_2();
      const double btilde1 = b_tilde_1();
      const double btilde2 = b_tilde_2();
      const double autilde1 = a_utilde_1();
      const double autilde2 = a_utilde_2();
      const double butilde1 = b_utilde_1();
      const double butilde2 = b_utilde_2();


      // Node 0
      // The Second tangent derivative in terms of basic dofs
      // D2 w(a0)(A1,A1) transforms to D2 w(Fk(a0))(s0,s0)
      b32(0, 0) = 1.0;

      //  D2 w(Fk(a0))(s0,s1) is more involved
      // First derivatives
      b31(0, 1) = (2 * atilde1 - autilde1) / 4.;
      b31(1, 1) = (2. + 2 * atilde2 - autilde2) / 4.;

      b31(0, 2) = atilde1 + autilde1 / 2.;
      b31(1, 2) = 1.0 + atilde2 + autilde2 / 2.;
      // Second derivatives
      b32(0, 1) = 1 + (1 + atilde1) / (2 * atilde2);
      b32(1, 1) = atilde2 / (2 + 2 * atilde1);
      b32(6, 1) = -1.0 / (2 * atilde2 * (1 + atilde1));

      b32(0, 2) = 1 + (1 + atilde1) / (atilde2);
      b32(1, 2) = 1 + atilde2 / (1 + atilde1);
      b32(6, 2) = -1.0 / (atilde2 * (1 + atilde1));

      // Node 1
      //  D2 w(Fk(a1))(s1,s1) is more involved
      // First derivatives
      b31(2, 3) = 1.0 + btilde1 + butilde1 / 2.;
      b31(3, 3) = btilde2 + butilde2 / 2.;

      b31(2, 4) = (2 + 2 * btilde1 - butilde1) / 4.;
      b31(3, 4) = (2 * btilde2 - butilde2) / 4.;
      // Second derivatives
      b32(2, 3) = 1 + btilde1 / (1 + btilde2);
      b32(3, 3) = 1 + (1 + btilde2) / (btilde1);
      b32(7, 3) = -1.0 / (btilde1 * (1 + btilde2));

      b32(2, 4) = btilde1 / (2 + 2 * btilde2);
      b32(3, 4) = 1 + (1 + btilde2) / (2 * btilde1);
      b32(7, 4) = -1.0 / (2 * btilde1 * (1 + btilde2));

      // D2 w(a2)(C1,C1) transforms to D2 w(Fk(a2))(s1,s1)
      b32(3, 5) = 1.0;

      // Node 3
      // D2 w(a1)(A1,A1) transforms to D2 w(Fk(a1))(s0,s0)
      b32(5, 6) = 1.0;

      //  D2 w(Fk(a0))(s0,s1) is more involved
      // First derivatives
      b31(4, 7) = (c_tilde_1() + c_tildetilde_1()) / 2.;
      b31(5, 7) = (c_tilde_2() + c_tildetilde_2()) / 2.;

      b32(4, 7) =
        -(c_tilde_1() * c_tildetilde_1()) /
        (c_tilde_2() * c_tildetilde_1() + c_tilde_1() * c_tildetilde_2());
      b32(5, 7) =
        -(c_tilde_2() * c_tildetilde_2()) /
        (c_tilde_2() * c_tildetilde_1() + c_tilde_1() * c_tildetilde_2());
      b32(8, 7) = -1.0 / (c_tilde_2() * c_tildetilde_1() +
                          c_tilde_1() * c_tildetilde_2());

      // D2 w(a2)(C2,C2) transforms to D2 w(Fk(a2))(s0,s0)
      b32(4, 8) = 1.0;

      // Now construct the submatrix
      for (unsigned j = 0; j < 9; ++j)
      {
        for (unsigned i = 0; i < 6; ++i)
        {
          // Fill in b31
          m3(i + 3, j) = b31(i, j);
          // Fill in b32 (up to 6 so we don't need a seperate loop)
          m3(i + 9, j) = b32(i, j);
        }
        // Carry on the for loop
        for (unsigned i = 6; i < 9; ++i)
        {
          // Fill in b32 (6 to 9)
          m3(i + 9, j) = b32(i, j);
        }
      }
    }

    /// Submatrix M4 (B_4 in Bernadou and Boisserie 1997)
    /// This transforms uses the trace functions to determine the normal
    /// derivative degrees of freedom on side at points bi This will be
    template<>
    void BernadouElementBasis<3>::basic_to_local_submatrix_4(
      DenseMatrix<double>& m4) const
    {
      // Constants
      double dFkds0_dot_B2_at_b0(0.0), dFkds0_dot_n0_at_b0(0.0),
        dFkds1_dot_A1_at_b1(0.0), dFkds1_dot_n1_at_b1(0.0);

      // Norms
      double norm2_B2(pow(B2(0), 2) + pow(B2(1), 2)),
        norm2_n0(pow(altitude_vector_1(0), 2) + pow(altitude_vector_1(1), 2));
      double norm2_A1(pow(A1(0), 2) + pow(A1(1), 2)),
        norm2_n1(pow(altitude_vector_2(0), 2) + pow(altitude_vector_2(1), 2));

      // Constant Vectors
      Vector<double> A_1(2), B_2(2), A_2(2), B_1(2);
      A1(A_1);
      A2(A_2);
      B1(B_1);
      B2(B_2);

      // The constant trace vectors at 0.5
      const Vector<double> g1 = g_1(0.5), g2 = g_2(0.5), g3 = g_3(0.5),
                           df1_ds = df_1_ds(0.5), df2_ds = df_2_ds(0.5);

      // Fill in the constants
      for (unsigned i = 0; i < 2; ++i)
      {
        // Jac row 0 at b0
        const double dFkds0_at_b0 =
          A_1[i] - (A_1[i] - B_2[i]) / 4. - (3 * B_1[i] + A_2[i]) / 8.;
        // Jac row 1 at b1
        const double dFkds1_at_b1 =
          B_2[i] + (A_1[i] - B_2[i]) / 4. - (3 * A_2[i] + B_1[i]) / 8.;

        // In B+B this are E1 and E2 respectively
        dFkds0_dot_B2_at_b0 -= dFkds0_at_b0 * B_2[i] / norm2_B2;
        dFkds0_dot_n0_at_b0 += dFkds0_at_b0 * altitude_vector_1(i) / norm2_n0;

        // In B+B this are F1 and F2 respectively
        dFkds1_dot_A1_at_b1 -= dFkds1_at_b1 * A_1[i] / norm2_A1;
        dFkds1_dot_n1_at_b1 += dFkds1_at_b1 * altitude_vector_2(i) / norm2_n1;
      }

      // Loop over the rows i
      for (unsigned i = 0; i < 21; ++i)
      {
        // Evaluate the traces vectors at points di
        m4(i, 0) =
          dFkds0_dot_B2_at_b0 * df2_ds[i] + dFkds0_dot_n0_at_b0 * g2[i];
        m4(i, 1) =
          dFkds1_dot_A1_at_b1 * df1_ds[i] + dFkds1_dot_n1_at_b1 * g1[i];
        m4(i, 2) = -sqrt(2) * g3[i];
      }
    }

    /// Submatrix M4 (B_4 in Bernadou and Boisserie 1997)
    /// This transforms uses the trace functions to determine the normal
    /// derivative degrees of freedom on side at points bi
    template<>
    void BernadouElementBasis<5>::basic_to_local_submatrix_4(
      DenseMatrix<double>& m4) const
    {
      // Constants
      double dFkds0_dot_B2_at_b0(0.0), dFkds0_dot_n0_at_b0(0.0),
        dFkds1_dot_A1_at_b1(0.0), dFkds1_dot_n1_at_b1(0.0);

      // Norms
      double norm2_B2(pow(B2(0), 2) + pow(B2(1), 2)),
        norm2_n0(pow(altitude_vector_1(0), 2) + pow(altitude_vector_1(1), 2));
      double norm2_A1(pow(A1(0), 2) + pow(A1(1), 2)),
        norm2_n1(pow(altitude_vector_2(0), 2) + pow(altitude_vector_2(1), 2));

      // Constant Vectors
      Vector<double> A_1(2), B_2(2), A_2(2), B_1(2), D_1(2), D_2(2);
      A1(A_1);
      A2(A_2);
      B1(B_1);
      B2(B_2);
      D1(D_1);
      D2(D_2);

      // Get the Jacobian at several points
      DenseMatrix<double> j_at_b1(2, 2, 0.0), j_at_b0(2, 2, 0.0);
      Vector<double> s(2);
      s[1] = 0;
      s[0] = 1 / 2.;
      get_basic_jacobian(s, j_at_b1);
      s[0] = 0;
      s[1] = 1 / 2.;
      get_basic_jacobian(s, j_at_b0);


      // The constant trace vectors at 0.5
      const Vector<double> g1 = g_1(0.5), g2 = g_2(0.5), g3 = g_3(0.5),
                           df1_ds = df_1_ds(0.5), df2_ds = df_2_ds(0.5);
      // Fill in the constants
      for (unsigned i = 0; i < 2; ++i)
      {
        //  // Jac row 0 at b0
        //  const double dFkds0_at_b0 = A_1[i] -
        //  (A_1[i]-B_2[i])/4.-(13*B_1[i]+5*A_2[i])/32.
        //    -(D_1[i]+D_2[i])/64.;
        //  // Jac row 1 at b1
        //  const double dFkds1_at_b1 = B_2[i] +
        //  (A_1[i]-B_2[i])/4.-(13*A_2[i]+5*B_1[i])/32.
        //    +(D_1[i]+D_2[i])/64.;

        // In B+B this are E1 and E2 respectively
        //    dFkds0_dot_B2_at_b0 -= dFkds0_at_b0 * B_2[i] / norm2_B2;
        //    dFkds0_dot_n0_at_b0 += dFkds0_at_b0 * altitude_vector_1(i) /
        //    norm2_n0;

        // In B+B this are F1 and F2 respectively
        //    dFkds1_dot_A1_at_b1 -= dFkds1_at_b1 * A_1[i] / norm2_A1;
        //    dFkds1_dot_n1_at_b1 += dFkds1_at_b1 * altitude_vector_2(i) /
        //    norm2_n1;

        dFkds0_dot_B2_at_b0 += j_at_b0(i, 0) * B_2[i] / norm2_B2;
        dFkds0_dot_n0_at_b0 -= j_at_b0(i, 0) * altitude_vector_1(i) / norm2_n0;

        dFkds1_dot_A1_at_b1 += j_at_b1(i, 1) * A_1[i] / norm2_A1;
        dFkds1_dot_n1_at_b1 -= j_at_b1(i, 1) * altitude_vector_2(i) / norm2_n1;
      }

      // Loop over the rows i
      for (unsigned i = 0; i < 21; ++i)
      {
        // Evaluate the traces vectors at points di
        m4(i, 0) =
          dFkds0_dot_B2_at_b0 * df2_ds[i] + dFkds0_dot_n0_at_b0 * g2[i];
        m4(i, 1) =
          dFkds1_dot_A1_at_b1 * df1_ds[i] + dFkds1_dot_n1_at_b1 * g1[i];
        m4(i, 2) = -sqrt(2) * g3[i];
      }
    }

    /// Submatrix M5 (B_5 in Bernadou and Boisserie 1997)
    /// This transforms uses the trace functions to determine the functional
    /// degrees of freedom along edges 1 2 and 3
    template<>
    void BernadouElementBasis<3>::basic_to_local_submatrix_5(
      DenseMatrix<double>& m5) const
    {
      const double nbasis = n_basis_functions();
      // Const vectors
      const Vector<double> trace_at_d1 = f_2(0.75), trace_at_d2 = f_2(0.25);
      const Vector<double> trace_at_d3 = f_1(0.25), trace_at_d4 = f_1(0.75);
      const Vector<double> trace_at_d5 = f_3(0.75), trace_at_d6 = f_3(0.25);
      // Loop over the rows
      for (unsigned i = 0; i < nbasis; ++i)
      {
        // Evaluate the traces vectors at points di
        m5(i, 0) = trace_at_d1[i];
        m5(i, 1) = trace_at_d2[i];
        m5(i, 2) = trace_at_d3[i];
        m5(i, 3) = trace_at_d4[i];
        m5(i, 4) = trace_at_d5[i];
        m5(i, 5) = trace_at_d6[i];
      }
    }

    /// Submatrix M5 (B_5 in Bernadou and Boisserie 1997)
    /// This transforms uses the trace functions to determine the functional
    /// degrees of freedom along edges 1 2 and 3
    template<>
    void BernadouElementBasis<5>::basic_to_local_submatrix_5(
      DenseMatrix<double>& m5) const
    {
      const double nbasis = n_basis_functions();
      // Const vectors
      const Vector<double> trace_at_d1 = f_2(5 / 6.), trace_at_d2 = f_2(2 / 3.);
      const Vector<double> trace_at_d3 = f_2(1 / 3.), trace_at_d4 = f_2(1 / 6.);
      const Vector<double> trace_at_d5 = f_1(1 / 6.), trace_at_d6 = f_1(1 / 3.);
      const Vector<double> trace_at_d7 = f_1(2 / 3.), trace_at_d8 = f_1(5 / 6.);
      const Vector<double> trace_at_d9 = f_3(5 / 6.),
                           trace_at_d10 = f_3(2 / 3.);
      const Vector<double> trace_at_d11 = f_3(1 / 3.),
                           trace_at_d12 = f_3(1 / 6.);
      // Loop over the rows
      for (unsigned i = 0; i < nbasis; ++i)
      {
        // Evaluate the traces vectors at points di
        m5(i, 0) = trace_at_d1[i];
        m5(i, 1) = trace_at_d2[i];
        m5(i, 2) = trace_at_d3[i];
        m5(i, 3) = trace_at_d4[i];
        m5(i, 4) = trace_at_d5[i];
        m5(i, 5) = trace_at_d6[i];
        m5(i, 6) = trace_at_d7[i];
        m5(i, 7) = trace_at_d8[i];
        m5(i, 8) = trace_at_d9[i];
        m5(i, 9) = trace_at_d10[i];
        m5(i, 10) = trace_at_d11[i];
        m5(i, 11) = trace_at_d12[i];
      }
    }

    /// Submatrix M5 (B_6 in Bernadou and Boisserie 1997)
    /// This transforms uses the trace functions to determine the normal
    /// derivative degrees of freeedom at points di This needs specialising for
    /// higher order
    template<>
    void BernadouElementBasis<3>::basic_to_local_submatrix_6(
      DenseMatrix<double>& m6) const
    {
      // Constants
      double dFkds0_dot_B2_at_d0(0.0), dFkds0_dot_n0_at_d0(0.0),
        dFkds0_dot_B2_at_d1(0.0), dFkds0_dot_n0_at_d1(0.0),
        dFkds1_dot_A1_at_d2(0.0), dFkds1_dot_n1_at_d2(0.0),
        dFkds1_dot_A1_at_d3(0.0), dFkds1_dot_n1_at_d3(0.0);

      // Norms
      double norm2_B2(pow(B2(0), 2) + pow(B2(1), 2)),
        norm2_n0(pow(altitude_vector_1(0), 2) + pow(altitude_vector_1(1), 2));
      double norm2_A1(pow(A1(0), 2) + pow(A1(1), 2)),
        norm2_n1(pow(altitude_vector_2(0), 2) + pow(altitude_vector_2(1), 2));

      // Const vectors
      const Vector<double> dtrace_at_d1 = df_2_ds(0.75),
                           dtrace_at_d2 = df_2_ds(0.25);
      const Vector<double> g_trace_at_d1 = g_2(0.75), g_trace_at_d2 = g_2(0.25);
      const Vector<double> dtrace_at_d3 = df_1_ds(0.25),
                           dtrace_at_d4 = df_1_ds(0.75);
      const Vector<double> g_trace_at_d3 = g_1(0.25), g_trace_at_d4 = g_1(0.75);
      const Vector<double> g_trace_at_d5 = g_3(0.75), g_trace_at_d6 = g_3(0.25);

      // Constant Vectors
      Vector<double> A_1(2), B_2(2), A_2(2), B_1(2);
      A1(A_1);
      A2(A_2);
      B1(B_1);
      B2(B_2);

      // Fill in the constants
      for (unsigned i = 0; i < 2; ++i)
      {
        // Jac row 0 at b0
        const double dFkds0_at_d0 = A_1[i] - 9 * (A_1[i] - B_2[i]) / 16. -
                                    3 * (7 * B_1[i] + A_2[i]) / 32.;
        const double dFkds0_at_d1 =
          A_1[i] - (A_1[i] - B_2[i]) / 16. - (5 * B_1[i] + 3 * A_2[i]) / 32.;
        // Jac row 1 at b1
        const double dFkds1_at_d2 =
          B_2[i] + (A_1[i] - B_2[i]) / 16. - (5 * A_2[i] + 3 * B_1[i]) / 32.;
        const double dFkds1_at_d3 = B_2[i] + 9 * (A_1[i] - B_2[i]) / 16. -
                                    3 * (7 * A_2[i] + B_1[i]) / 32.;

        // Fill in more specifically to reduce round--off
        dFkds0_dot_B2_at_d0 -= dFkds0_at_d0 * B_2[i] / norm2_B2;
        dFkds0_dot_n0_at_d0 += dFkds0_at_d0 * altitude_vector_1(i) / norm2_n0;
        dFkds0_dot_B2_at_d1 -= dFkds0_at_d1 * B_2[i] / norm2_B2;
        dFkds0_dot_n0_at_d1 += dFkds0_at_d1 * altitude_vector_1(i) / norm2_n0;

        dFkds1_dot_A1_at_d2 -= dFkds1_at_d2 * A_1[i] / norm2_A1;
        dFkds1_dot_n1_at_d2 += dFkds1_at_d2 * altitude_vector_2(i) / norm2_n1;
        dFkds1_dot_A1_at_d3 -= dFkds1_at_d3 * A_1[i] / norm2_A1;
        dFkds1_dot_n1_at_d3 += dFkds1_at_d3 * altitude_vector_2(i) / norm2_n1;
      }

      // Loop over the rows i
      for (unsigned i = 0; i < 21; ++i)
      {
        // Evaluate the traces vectors at points di
        m6(i, 0) = dFkds0_dot_B2_at_d0 * dtrace_at_d1[i] +
                   dFkds0_dot_n0_at_d0 * g_trace_at_d1[i];

        m6(i, 1) = dFkds0_dot_B2_at_d1 * dtrace_at_d2[i] +
                   dFkds0_dot_n0_at_d1 * g_trace_at_d2[i];

        m6(i, 2) = dFkds1_dot_A1_at_d2 * dtrace_at_d3[i] +
                   dFkds1_dot_n1_at_d2 * g_trace_at_d3[i];

        m6(i, 3) = dFkds1_dot_A1_at_d3 * dtrace_at_d4[i] +
                   dFkds1_dot_n1_at_d3 * g_trace_at_d4[i];

        m6(i, 4) = -sqrt(2) * g_trace_at_d5[i];

        m6(i, 5) = -sqrt(2) * g_trace_at_d6[i];
      }
    }

    /// Submatrix M5 (B_6 in Bernadou and Boisserie 1997)
    /// This transforms uses the trace functions to determine the normal
    /// derivative degrees of freeedom at points di This needs specialising for
    /// higher order
    template<>
    void BernadouElementBasis<5>::basic_to_local_submatrix_6(
      DenseMatrix<double>& m6) const
    {
      // Constants
      double dFkds0_dot_B2_at_d0(0.0), dFkds0_dot_n0_at_d0(0.0),
        dFkds0_dot_B2_at_d1(0.0), dFkds0_dot_n0_at_d1(0.0),
        dFkds0_dot_B2_at_d2(0.0), dFkds0_dot_n0_at_d2(0.0),
        dFkds0_dot_B2_at_d3(0.0), dFkds0_dot_n0_at_d3(0.0),
        dFkds1_dot_A1_at_d4(0.0), dFkds1_dot_n1_at_d4(0.0),
        dFkds1_dot_A1_at_d5(0.0), dFkds1_dot_n1_at_d5(0.0),
        dFkds1_dot_A1_at_d6(0.0), dFkds1_dot_n1_at_d6(0.0),
        dFkds1_dot_A1_at_d7(0.0), dFkds1_dot_n1_at_d7(0.0);

      // Constant Vectors
      Vector<double> A_1(2), B_2(2), A_2(2), B_1(2);
      A1(A_1);
      A2(A_2);
      B1(B_1);
      B2(B_2);

      // Norms
      double norm2_B2(pow(B_2[0], 2) + pow(B_2[1], 2)),
        norm2_n0(pow(altitude_vector_1(0), 2) + pow(altitude_vector_1(1), 2));
      double norm2_A1(pow(A_1[0], 2) + pow(A_1[1], 2)),
        norm2_n1(pow(altitude_vector_2(0), 2) + pow(altitude_vector_2(1), 2));

      // Const vectors
      const Vector<double> df_ds_at_d0 = df_2_ds(5 / 6.),
                           df_ds_at_d1 = df_2_ds(2 / 3.);
      const Vector<double> df_ds_at_d2 = df_2_ds(1 / 3.),
                           df_ds_at_d3 = df_2_ds(1 / 6.);
      const Vector<double> df_ds_at_d4 = df_1_ds(1 / 6.),
                           df_ds_at_d5 = df_1_ds(1 / 3.);
      const Vector<double> df_ds_at_d6 = df_1_ds(2 / 3.),
                           df_ds_at_d7 = df_1_ds(5 / 6.);
      // The g traces
      const Vector<double> g_at_d0 = g_2(5 / 6.), g_at_d1 = g_2(2 / 3.);
      const Vector<double> g_at_d2 = g_2(1 / 3.), g_at_d3 = g_2(1 / 6.);
      const Vector<double> g_at_d4 = g_1(1 / 6.), g_at_d5 = g_1(1 / 3.);
      const Vector<double> g_at_d6 = g_1(2 / 3.), g_at_d7 = g_1(5 / 6.);
      const Vector<double> g_at_d8 = g_3(5 / 6.), g_at_d9 = g_3(2 / 3.);
      const Vector<double> g_at_d10 = g_3(1 / 3.), g_at_d11 = g_3(1 / 6.);

      // Get the Jacobian at several points
      DenseMatrix<double> j_at_d0(2, 2, 0.0), j_at_d1(2, 2, 0.0),
        j_at_d2(2, 2, 0.0), j_at_d3(2, 2, 0.0), j_at_d4(2, 2, 0.0),
        j_at_d5(2, 2, 0.0), j_at_d6(2, 2, 0.0), j_at_d7(2, 2, 0.0);
      // Fill in Jacobian
      Vector<double> s(2, 0.0);
      s[1] = 5 / 6.;
      get_basic_jacobian(s, j_at_d0);
      s[1] = 2 / 3.;
      get_basic_jacobian(s, j_at_d1);
      s[1] = 1 / 3.;
      get_basic_jacobian(s, j_at_d2);
      s[1] = 1 / 6.;
      get_basic_jacobian(s, j_at_d3);

      s[1] = 0.0;
      s[0] = 1 / 6.;
      get_basic_jacobian(s, j_at_d4);
      s[0] = 1 / 3.;
      get_basic_jacobian(s, j_at_d5);
      s[0] = 2 / 3.;
      get_basic_jacobian(s, j_at_d6);
      s[0] = 5 / 6.;
      get_basic_jacobian(s, j_at_d7);

      // Fill in the constants
      for (unsigned i = 0; i < 2; ++i)
      {
        // Fill in the dot products
        dFkds0_dot_B2_at_d0 += j_at_d0(i, 0) * B_2[i] / norm2_B2;
        dFkds0_dot_n0_at_d0 -= j_at_d0(i, 0) * altitude_vector_1(i) / norm2_n0;
        dFkds0_dot_B2_at_d1 += j_at_d1(i, 0) * B_2[i] / norm2_B2;
        dFkds0_dot_n0_at_d1 -= j_at_d1(i, 0) * altitude_vector_1(i) / norm2_n0;

        dFkds0_dot_B2_at_d2 += j_at_d2(i, 0) * B_2[i] / norm2_B2;
        dFkds0_dot_n0_at_d2 -= j_at_d2(i, 0) * altitude_vector_1(i) / norm2_n0;
        dFkds0_dot_B2_at_d3 += j_at_d3(i, 0) * B_2[i] / norm2_B2;
        dFkds0_dot_n0_at_d3 -= j_at_d3(i, 0) * altitude_vector_1(i) / norm2_n0;

        dFkds1_dot_A1_at_d4 += j_at_d4(i, 1) * A_1[i] / norm2_A1;
        dFkds1_dot_n1_at_d4 -= j_at_d4(i, 1) * altitude_vector_2(i) / norm2_n1;
        dFkds1_dot_A1_at_d5 += j_at_d5(i, 1) * A_1[i] / norm2_A1;
        dFkds1_dot_n1_at_d5 -= j_at_d5(i, 1) * altitude_vector_2(i) / norm2_n1;

        dFkds1_dot_A1_at_d6 += j_at_d6(i, 1) * A_1[i] / norm2_A1;
        dFkds1_dot_n1_at_d6 -= j_at_d6(i, 1) * altitude_vector_2(i) / norm2_n1;
        dFkds1_dot_A1_at_d7 += j_at_d7(i, 1) * A_1[i] / norm2_A1;
        dFkds1_dot_n1_at_d7 -= j_at_d7(i, 1) * altitude_vector_2(i) / norm2_n1;
      }

      // Loop over the rows i
      for (unsigned i = 0; i < 21; ++i)
      {
        // Evaluate the traces vectors at points di on face 0
        m6(i, 0) = dFkds0_dot_B2_at_d0 * df_ds_at_d0[i] +
                   dFkds0_dot_n0_at_d0 * g_at_d0[i];

        m6(i, 1) = dFkds0_dot_B2_at_d1 * df_ds_at_d1[i] +
                   dFkds0_dot_n0_at_d1 * g_at_d1[i];

        m6(i, 2) = dFkds0_dot_B2_at_d2 * df_ds_at_d2[i] +
                   dFkds0_dot_n0_at_d2 * g_at_d2[i];

        m6(i, 3) = dFkds0_dot_B2_at_d3 * df_ds_at_d3[i] +
                   dFkds0_dot_n0_at_d3 * g_at_d3[i];

        // Evaluate the traces vectors at points di on face 1
        m6(i, 4) = dFkds1_dot_A1_at_d4 * df_ds_at_d4[i] +
                   dFkds1_dot_n1_at_d4 * g_at_d4[i];

        m6(i, 5) = dFkds1_dot_A1_at_d5 * df_ds_at_d5[i] +
                   dFkds1_dot_n1_at_d5 * g_at_d5[i];

        m6(i, 6) = dFkds1_dot_A1_at_d6 * df_ds_at_d6[i] +
                   dFkds1_dot_n1_at_d6 * g_at_d6[i];

        m6(i, 7) = dFkds1_dot_A1_at_d7 * df_ds_at_d7[i] +
                   dFkds1_dot_n1_at_d7 * g_at_d7[i];

        // Evaluate the traces vectors at points di on face 2
        m6(i, 8) = -sqrt(2) * g_at_d8[i];

        m6(i, 9) = -sqrt(2) * g_at_d9[i];

        m6(i, 10) = -sqrt(2) * g_at_d10[i];

        m6(i, 11) = -sqrt(2) * g_at_d11[i];
      }
    }
    /// Submatrix M5 (B_7 in Bernadou and Boisserie 1997)
    /// The (functional) internal degrees of freedom, simply equal to the values
    /// of function at w(Fk(ei))
    template<unsigned BOUNDARY_ORDER>
    void BernadouElementBasis<BOUNDARY_ORDER>::basic_to_local_submatrix_7(
      DenseMatrix<double>& m7) const
    {
      // Identity Matrix
      for (unsigned i = 0; i < n_internal_dofs(); ++i)
      {
        m7(18 + i, i) = 1.0;
      }
    }

    /// This matrix transforms between the basic dofs (36/55) and the local dofs
    /// (21/28)
    template<unsigned BOUNDARY_ORDER>
    void BernadouElementBasis<BOUNDARY_ORDER>::basic_to_local_matrix(
      DenseMatrix<double>& m) const
    {
      // Handy definitions
      const unsigned nbasis = n_basis_functions();
      const unsigned ninternal = n_internal_dofs();
      const unsigned nmidnodes = n_basic_midside_nodes();

      // The  matrix should be a:
      // DenseMatrix<double> m(nbasis,nbasic,0.0);

      // Initialise the Submatrices

      DenseMatrix<double> m1(nbasis, 3, 0.0), m2(nbasis, 6, 0.0),
        m3(nbasis, 9, 0.0), m4(nbasis, 3, 0.0), m5(nbasis, nmidnodes, 0.0),
        m6(nbasis, nmidnodes, 0.0), m7(nbasis, ninternal, 0.0);

      // Fill it
      basic_to_local_submatrix_1(m1);
      basic_to_local_submatrix_2(m2);
      basic_to_local_submatrix_3(m3);
      basic_to_local_submatrix_4(m4);
      basic_to_local_submatrix_5(m5);
      basic_to_local_submatrix_6(m6);
      basic_to_local_submatrix_7(m7);

      // Now construct the full matrix
      for (unsigned i = 0; i < nbasis; ++i)
      {
        //  fill in m1
        for (unsigned j = 0; j < 3; ++j)
        {
          m(i, j) = m1(i, j);
        }
        //  fill in m2
        for (unsigned j = 0; j < 6; ++j)
        {
          m(i, 3 + j) = m2(i, j);
        }
        //  fill in m3
        for (unsigned j = 0; j < 9; ++j)
        {
          m(i, 9 + j) = m3(i, j);
        }
        //  fill in m4
        for (unsigned j = 0; j < 3; ++j)
        {
          m(i, 18 + j) = m4(i, j);
        }
        //  fill in m5
        for (unsigned j = 0; j < nmidnodes; ++j)
        {
          m(i, 21 + j) = m5(i, j);
        }
        //  fill in m6
        for (unsigned j = 0; j < nmidnodes; ++j)
        {
          m(i, 21 + nmidnodes + j) = m6(i, j);
        }
        //  fill in m7
        for (unsigned j = 0; j < ninternal; ++j)
        {
          m(i, 21 + 2 * nmidnodes + j) = m7(i, j);
        }
      }
    }

    /// Get the shape functions, with s as the basic coordinate
    template<unsigned BOUNDARY_ORDER>
    void BernadouElementBasis<BOUNDARY_ORDER>::shape_basic(
      const Vector<double>& s, Shape& psi, Shape& bpsi) const
    {
      // Handy definitions
      const unsigned nbasis = n_basis_functions();
      const unsigned nbasic = n_basic_basis_functions();
      const unsigned ninternal = n_internal_dofs();

      // Initialise the shape functions
      Shape p7(nbasic);
      DenseMatrix<double> b_matrix(nbasis, nbasic, 0.0);
      DenseMatrix<double> d_matrix(nbasis, nbasis, 0.0);
      DenseMatrix<double> conversion_matrix(nbasis, nbasic, 0.0);
      DenseMatrix<double> gl2basic_matrix(nbasis, nbasic, 0.0);

      // Fill in the matrices
      basic_to_local_matrix(b_matrix);
      local_to_global_matrix(d_matrix);
      full_basic_polynomials(s, p7);

      // Zero the shape functions
      for (unsigned i = 0; i < 3; ++i)
      {
        // Zero nodal dofs
        for (unsigned l = 0; l < 6; ++l)
        {
          psi(i, l) = 0;
        }
      }
      for (unsigned i = 0; i < ninternal; ++i)
      {
        // Zero internal (bubble) dofs
        bpsi(i, 0) = 0;
      }
      // Get the permutation
      unsigned index_shift = 0;
      nodal_index_shift(index_shift);

      // Fill in the shape functions
      for (unsigned k = 0; k < nbasis; ++k)
      {
        for (unsigned l = 0; l < nbasic; ++l)
        {
          // This can be optimised HERE
          // Nodal functional dofs
          for (unsigned i = 0; i < 3; ++i)
          {
            // Permuted i so the global index will work
            unsigned i_global = (i + index_shift) % 3;
            psi(i_global, 0) += d_matrix(i, k) * b_matrix(k, l) * p7[l];
          }
          // Nodal first derivative  dofs
          for (unsigned i = 0; i < 6; ++i)
          {
            // Convert index to node index, inode and dof type index, itype
            const unsigned inode(i / 2), itype(i % 2);
            // Permuted i so the global index will work
            unsigned i_global = (inode + index_shift) % 3;
            psi(i_global, 1 + itype) +=
              d_matrix(i + 3, k) * b_matrix(k, l) * p7[l];
          }
          // Nodal second derivative dofs
          for (unsigned i = 0; i < 9; ++i)
          {
            // Convert index to node index, inode and dof type index, itype
            const unsigned inode(i / 3), itype(i % 3);
            // Permuted i so the global index will work
            unsigned i_global = (inode + index_shift) % 3;
            psi(i_global, 3 + itype) +=
              d_matrix(i + 9, k) * b_matrix(k, l) * p7[l];
          }
          // Now for the bubble functions, only Lagrange type dofs
          for (unsigned inode = 0; inode < ninternal; ++inode)
          {
            // Permuted i so the global index will work
            // Removed the index shift
            unsigned i_global = (inode) /*+ index_shift) % 3*/;
            bpsi(i_global, 0) +=
              d_matrix(inode + 18, k) * b_matrix(k, l) * p7[l];
          }
        }
      }
    }

    /// Get the derivatives of the shape functions, with s as the basic
    /// coordinate
    template<unsigned BOUNDARY_ORDER>
    void BernadouElementBasis<BOUNDARY_ORDER>::d_shape_ds(
      const Vector<double>& s,
      Shape& psi,
      Shape& bpsi,
      DShape& dpsi,
      DShape& dbpsi) const
    {
      // Handy definitions
      const unsigned nbasis = n_basis_functions();
      const unsigned nbasic = n_basic_basis_functions();
      const unsigned ninternal = n_internal_dofs();

      // Initialise the containers
      DenseMatrix<double> b_matrix(nbasis, nbasic, 0.0);
      DenseMatrix<double> d_matrix(nbasis, nbasis, 0.0);
      DenseMatrix<double> conversion_matrix(nbasis, nbasic, 0.0);
      DenseMatrix<double> gl2basic_matrix(nbasis, nbasic, 0.0);
      Shape p7(nbasic);
      DShape dp7(nbasic, 2);

      // Fill in the matrices
      basic_to_local_matrix(b_matrix);
      local_to_global_matrix(d_matrix);
      full_basic_polynomials(s, p7);
      dfull_basic_polynomials(s, dp7);

      // Zero the shape functions
      for (unsigned i = 0; i < 3; ++i)
      {
        // Zero nodal dofs
        for (unsigned l = 0; l < 6; ++l)
        {
          // Shape functions
          psi(i, l) = 0;
          // First Derivatives
          for (unsigned alpha = 0; alpha < 2; ++alpha)
          {
            dpsi(i, l, alpha) = 0;
          }
        }
      }
      for (unsigned i = 0; i < ninternal; ++i)
      {
        // Zero internal (bubble) dofs
        bpsi(i, 0) = 0;
        // First Derivatives
        for (unsigned alpha = 0; alpha < 2; ++alpha)
        {
          dbpsi(i, 0, alpha) = 0;
        }
      }

      /* Could do it like this:

      for(unsigned i=0;i<21;++i)
       for(unsigned j=0;j<36;++j)
         for(unsigned k=0;k<21;++k)
          for(unsigned l=0;l<36;++l)
            conversion_matrix(i,j)+=d_matrix(i,k)*b_matrix(k,l)*a_matrix(l,j);

       If we don't take into account sparsity of the d_matrix - but we can do
      better
      */

      // This way is about six times faster:
      // We reduce the loop over 404 elements of d_matrix to loop over 27
      // Fill in conversion matrix
      // First submatrix, identity i==k
      for (unsigned i = 0; i < 3; ++i)
      {
        for (unsigned l = 0; l < nbasic; ++l)
        {
          conversion_matrix(i, l) += d_matrix(i, i) * b_matrix(i, l);
        }
      }

      // Second submatrices, each is a two by two
      for (unsigned i = 0; i < 2; ++i)
      {
        for (unsigned k = 0; k < 2; ++k)
        {
          for (unsigned l = 0; l < nbasic; ++l)
          {
            conversion_matrix(i + 3, l) +=
              d_matrix(i + 3, k + 3) * b_matrix(k + 3, l);
            conversion_matrix(i + 5, l) +=
              d_matrix(i + 5, k + 5) * b_matrix(k + 5, l);
            conversion_matrix(i + 7, l) +=
              d_matrix(i + 7, k + 7) * b_matrix(k + 7, l);
          }
        }
      }

      // Second submatrices, each is a three by three
      for (unsigned i = 0; i < 3; ++i)
      {
        for (unsigned l = 0; l < nbasic; ++l)
        {
          for (unsigned k = 0; k < 2; ++k)
          {
            conversion_matrix(i + 9, l) +=
              d_matrix(i + 9, k + 9) * b_matrix(k + 9, l);
            conversion_matrix(i + 12, l) +=
              d_matrix(i + 12, k + 11) * b_matrix(k + 11, l);
            conversion_matrix(i + 15, l) +=
              d_matrix(i + 15, k + 13) * b_matrix(k + 13, l);
          }
          // Now do the third column

          conversion_matrix(i + 9, l) += d_matrix(i + 9, 15) * b_matrix(15, l);
          conversion_matrix(i + 12, l) +=
            d_matrix(i + 12, 16) * b_matrix(16, l);
          conversion_matrix(i + 15, l) +=
            d_matrix(i + 15, 17) * b_matrix(17, l);
        }
      }

      // Last submatrix, identity i==k
      for (unsigned i = 0; i < ninternal; ++i)
      {
        for (unsigned l = 0; l < nbasic; ++l)
        {
          conversion_matrix(i + 18, l) +=
            d_matrix(i + 18, i + 18) * b_matrix(i + 18, l);
        }
      }

      // Get the permutation
      unsigned index_shift = 0;
      nodal_index_shift(index_shift);

      // Fill in the shape function
      for (unsigned j = 0; j < nbasic; ++j)
      {
        // Nodal functional dofs
        for (unsigned i = 0; i < 3; ++i)
        {
          // Permuted i so the global index will work
          unsigned i_global = (i + index_shift) % 3;
          psi(i_global, 0) += conversion_matrix(i, j) * p7[j];
          // Local Derivatives
          for (unsigned alpha = 0; alpha < 2; ++alpha)
          {
            dpsi(i_global, 0, alpha) += conversion_matrix(i, j) * dp7(j, alpha);
          }
        }
        // Nodal first derivative  dofs
        for (unsigned i = 0; i < 6; ++i)
        {
          const unsigned inode(i / 2), ideriv(i % 2);
          // Permuted i so the global index will work
          unsigned i_global = (inode + index_shift) % 3;
          psi(i_global, 1 + ideriv) += conversion_matrix(i + 3, j) * p7[j];
          // Local Derivatives
          for (unsigned alpha = 0; alpha < 2; ++alpha)
          {
            dpsi(i_global, 1 + ideriv, alpha) +=
              conversion_matrix(i + 3, j) * dp7(j, alpha);
          }
        }
        // Nodal second derivative  dofs
        for (unsigned i = 0; i < 9; ++i)
        {
          const unsigned inode(i / 3), ideriv(i % 3);
          // Permuted i so the global index will work
          unsigned i_global = (inode + index_shift) % 3;

          psi(i_global, 3 + ideriv) += conversion_matrix(i + 9, j) * p7[j];
          // Local Derivatives
          for (unsigned alpha = 0; alpha < 2; ++alpha)
          {
            dpsi(i_global, 3 + ideriv, alpha) +=
              conversion_matrix(i + 9, j) * dp7(j, alpha);
          }
        }

        // Now for the bubble functions
        for (unsigned i = 0; i < ninternal; ++i)
        {
          // Permuted i so the global index will work
          unsigned i_global = i; // Don't shift these (i + index_shift) % 3;
          bpsi(i_global, 0) += conversion_matrix(i + 18, j) * p7[j];
          // Local Derivatives
          for (unsigned alpha = 0; alpha < 2; ++alpha)
          {
            dbpsi(i_global, 0, alpha) +=
              conversion_matrix(i + 18, j) * dp7(j, alpha);
          }
        }
      }
    }

    /// Get the first derivative shape functions. All of the construction
    /// happens elsewhere this just contains the matrix multiplications (which
    /// is done once
    ///  at the start to avoid doing it six times).
    template<unsigned BOUNDARY_ORDER>
    double BernadouElementBasis<BOUNDARY_ORDER>::d_shape_dx(
      const Vector<double>& s,
      Shape& psi,
      Shape& bpsi,
      DShape& dpsi,
      DShape& dbpsi) const
    {
// check the construction of the elements is complete
#ifdef paranoid
      if (curved_edge == none)
      {
        throw oomphliberror("the element has not been upgraded yet. did \
  you forget to set upe the curved_edge?",
                            oomph_current_function,
                            oomph_exception_location);
      }
#endif
      // Handy definitions
      const unsigned ninternal = n_internal_dofs();

      // Permute the local coordinate
      Vector<double> s_permuted(s);
      permute_shape(s_permuted);

      // Initialise shape
      DShape dpsi_ds(3, 6, 2);

      // Initialise bubble shape
      DShape dbpsi_ds(ninternal, 1, 2);

      // Get the local derivatives
      d_shape_ds(s_permuted, psi, bpsi, dpsi_ds, dbpsi_ds);

      // Now transform
      DenseMatrix<double> jacobian(2, 2, 0.0), inv_jacobian(2, 2, 0.0);

      // Fill in Jacobian and Hessian
      get_basic_jacobian(s_permuted, jacobian);

      // Now invert
      double det = MyShape::invert_two_by_two(jacobian, inv_jacobian);

      // Zero the shape functions
      for (unsigned i = 0; i < 3; ++i)
      {
        // Zero nodal dofs
        for (unsigned l = 0; l < 6; ++l)
        {
          // First Derivatives
          for (unsigned alpha = 0; alpha < 2; ++alpha)
          {
            dpsi(i, l, alpha) = 0;
          }
        }
      }
      for (unsigned i = 0; i < ninternal; ++i)
      {
        // Zero internal (bubble) dofs
        // First Derivatives
        for (unsigned alpha = 0; alpha < 2; ++alpha)
        {
          dbpsi(i, 0, alpha) = 0;
        }
      }

      // Now find the global derivatives at local coordinates
      for (unsigned alpha = 0; alpha < 2; ++alpha)
      {
        for (unsigned beta = 0; beta < 2; ++beta)
        {
          for (unsigned i = 0; i < ninternal; ++i)
          {
            // Transform the first derivatives of the bubble basis functions
            dbpsi(i, 0, beta) +=
              dbpsi_ds(i, 0, alpha) * inv_jacobian(alpha, beta);
          }
          for (unsigned i = 0; i < 3; ++i)
          {
            // Transform the first derivatives of the nodal basis functions
            for (unsigned l = 0; l < 6; ++l)
              dpsi(i, l, beta) +=
                dpsi_ds(i, l, alpha) * inv_jacobian(alpha, beta);
          }
        }
      }

      // Return the det
      return det;
    }

    /// Get the second derivative shape functions with s as the basic
    /// coordinate. All
    ///  of the construction happens elsewhere this just contains the matrix
    ///  multiplication (which is done once at the start to avoid doing it six
    ///  times).
    /// THIS IS A BOTTLENECK
    template<unsigned BOUNDARY_ORDER>
    void BernadouElementBasis<BOUNDARY_ORDER>::d2_shape_ds2(
      const Vector<double>& s,
      Shape& psi,
      Shape& bpsi,
      DShape& dpsi,
      DShape& dbpsi,
      DShape& d2psi,
      DShape& d2bpsi) const
    {
      // Handy definitions
      const unsigned nbasis = n_basis_functions();
      const unsigned nbasic = n_basic_basis_functions();
      const unsigned ninternal = n_internal_dofs();

      // Initialise the containers
      Shape p7(nbasic);
      DShape dp7(nbasic, 2), d2p7(nbasic, 3);

      DenseMatrix<double> b_matrix(nbasis, nbasic, 0.0),
        d_matrix(nbasis, nbasis, 0.0), conversion_matrix(nbasis, nbasic, 0.0);


      // Fill in the matrices
      basic_to_local_matrix(b_matrix);
      local_to_global_matrix(d_matrix);

      // Get the monomials
      full_basic_polynomials(s, p7);
      dfull_basic_polynomials(s, dp7);
      d2full_basic_polynomials(s, d2p7);

      /* Could do it like this:

      for(unsigned i=0;i<nbasis;++i)
       for(unsigned j=0;j<nbasic;++j)
         for(unsigned k=0;k<nbasis;++k)
          for(unsigned l=0;l<nbasic;++l)
            conversion_matrix(i,j)+=d_matrix(i,k)*b_matrix(k,l)*a_matrix(l,j);

       If we don't take into account sparsity of the d_matrix - but we can do
      better
      */

      // This way is about six times faster:
      // We reduce the loop over 404 elements of d_matrix to loop over 27
      // Fill in conversion matrix
      // First submatrix, identity i==k
      for (unsigned i = 0; i < 3; ++i)
      {
        for (unsigned l = 0; l < nbasic; ++l)
        {
          conversion_matrix(i, l) += d_matrix(i, i) * b_matrix(i, l);
        }
      }

      // Second submatrices, each is a two by two
      for (unsigned i = 0; i < 2; ++i)
      {
        for (unsigned k = 0; k < 2; ++k)
        {
          for (unsigned l = 0; l < nbasic; ++l)
          {
            conversion_matrix(i + 3, l) +=
              d_matrix(i + 3, k + 3) * b_matrix(k + 3, l);
            conversion_matrix(i + 5, l) +=
              d_matrix(i + 5, k + 5) * b_matrix(k + 5, l);
            conversion_matrix(i + 7, l) +=
              d_matrix(i + 7, k + 7) * b_matrix(k + 7, l);
          }
        }
      }

      // Second submatrices, each is a three by three
      for (unsigned i = 0; i < 3; ++i)
      {
        for (unsigned l = 0; l < nbasic; ++l)
        {
          for (unsigned k = 0; k < 2; ++k)
          {
            conversion_matrix(i + 9, l) +=
              d_matrix(i + 9, k + 9) * b_matrix(k + 9, l);
            conversion_matrix(i + 12, l) +=
              d_matrix(i + 12, k + 11) * b_matrix(k + 11, l);
            conversion_matrix(i + 15, l) +=
              d_matrix(i + 15, k + 13) * b_matrix(k + 13, l);
          }
          // Now do the third column
          conversion_matrix(i + 9, l) += d_matrix(i + 9, 15) * b_matrix(15, l);
          conversion_matrix(i + 12, l) +=
            d_matrix(i + 12, 16) * b_matrix(16, l);
          conversion_matrix(i + 15, l) +=
            d_matrix(i + 15, 17) * b_matrix(17, l);
        }
      }

      // Last submatrix, identity i==k
      for (unsigned i = 0; i < ninternal; ++i)
      {
        for (unsigned l = 0; l < nbasic; ++l)
        {
          conversion_matrix(i + 18, l) +=
            d_matrix(i + 18, i + 18) * b_matrix(i + 18, l);
        }
      }

      // Zero the shape functions
      for (unsigned i = 0; i < 3; ++i)
      {
        // Zero nodal dofs
        for (unsigned l = 0; l < 6; ++l)
        {
          // Shape functions
          psi(i, l) = 0;
          // First Derivatives
          for (unsigned alpha = 0; alpha < 2; ++alpha)
          {
            dpsi(i, l, alpha) = 0;
            // Second derivatives
            for (unsigned beta = alpha; beta < 2; ++beta)
            {
              d2psi(i, l, alpha + beta) = 0;
            }
          }
        }
      }

      for (unsigned i = 0; i < ninternal; ++i)
      {
        // Zero internal (bubble) dofs
        bpsi(i, 0) = 0;
        // First Derivatives
        for (unsigned alpha = 0; alpha < 2; ++alpha)
        {
          dbpsi(i, 0, alpha) = 0;
          // Second derivatives
          for (unsigned beta = alpha; beta < 2; ++beta)
          {
            d2bpsi(i, 0, alpha + beta) = 0;
          }
        }
      }

      // Get the permutation
      unsigned index_shift = 0;
      nodal_index_shift(index_shift);

      // Fill in the shape functions
      for (unsigned j = 0; j < nbasic; ++j)
      {
        // Nodal functional dofs
        for (unsigned i = 0; i < 3; ++i)
        {
          // Permuted i so the global index will work
          unsigned i_global = (i + index_shift) % 3;
          // Get the shape functions
          psi(i_global, 0) += conversion_matrix(i, j) * p7[j];
          // Local Derivatives
          for (unsigned alpha = 0; alpha < 2; ++alpha)
          {
            dpsi(i_global, 0, alpha) += conversion_matrix(i, j) * dp7(j, alpha);
            // Loop over the second derivatives
            for (unsigned beta = alpha; beta < 2; ++beta)
            {
              d2psi(i_global, 0, alpha + beta) +=
                conversion_matrix(i, j) * d2p7(j, alpha + beta);
            }
          }
        }
        // Nodal first derivative  dofs
        for (unsigned i = 0; i < 6; ++i)
        {
          const unsigned inode(i / 2), ideriv(i % 2);
          // Permuted i so the global dofs are consistent between elements
          unsigned i_global = (inode + index_shift) % 3;

          psi(i_global, 1 + ideriv) += conversion_matrix(i + 3, j) * p7[j];
          // Local Derivatives
          for (unsigned alpha = 0; alpha < 2; ++alpha)
          {
            dpsi(i_global, 1 + ideriv, alpha) +=
              conversion_matrix(i + 3, j) * dp7(j, alpha);
            // Loop over the second derivatives
            for (unsigned beta = alpha; beta < 2; ++beta)
            {
              d2psi(i_global, 1 + ideriv, alpha + beta) +=
                conversion_matrix(i + 3, j) * d2p7(j, alpha + beta);
            }
          }
        }
        // Nodal second derivative  dofs
        for (unsigned i = 0; i < 9; ++i)
        {
          const unsigned inode(i / 3), ideriv(i % 3);
          // Permuted i so the global dofs are consistent between elements
          unsigned i_global = (inode + index_shift) % 3;

          psi(i_global, 3 + ideriv) += conversion_matrix(i + 9, j) * p7[j];
          // Local Derivatives
          for (unsigned alpha = 0; alpha < 2; ++alpha)
          {
            dpsi(i_global, 3 + ideriv, alpha) +=
              conversion_matrix(i + 9, j) * dp7(j, alpha);
            // Loop over the second derivatives
            for (unsigned beta = alpha; beta < 2; ++beta)
            {
              d2psi(i_global, 3 + ideriv, alpha + beta) +=
                conversion_matrix(i + 9, j) * d2p7(j, alpha + beta);
            }
          }
        }
        // Now for the bubble functions
        for (unsigned i = 0; i < ninternal; ++i)
        {
          // Permuted i - to be consistent
          unsigned i_global = i; // (i + index_shift) % 3; DOn't shift these!
          bpsi(i_global, 0) += conversion_matrix(i + 18, j) * p7[j];
          // Local Derivatives
          for (unsigned alpha = 0; alpha < 2; ++alpha)
          {
            dbpsi(i_global, 0, alpha) +=
              conversion_matrix(i + 18, j) * dp7(j, alpha);
            // Loop over the second derivatives
            for (unsigned beta = alpha; beta < 2; ++beta)
            {
              d2bpsi(i_global, 0, alpha + beta) +=
                conversion_matrix(i + 18, j) * d2p7(j, alpha + beta);
            }
          }
        }
      } // End of loop over shape functions
    }

    /// Fill in the Full association matrix for the element: this transforms
    /// from the basic monomials/ shape functions to the physical basic
    /// functions.
    template<unsigned BOUNDARY_ORDER>
    void BernadouElementBasis<BOUNDARY_ORDER>::fill_in_full_association_matrix(
      DenseMatrix<double>& conversion_matrix) const
    {
// check the construction of the elements is complete
#ifdef paranoid
      if (curved_edge == none)
      {
        throw oomphliberror("the element has not been upgraded yet. did \
  you forget to set upe the curved_edge?",
                            oomph_current_function,
                            oomph_exception_location);
      }
#endif
      // Handy definitions
      const unsigned nbasis = n_basis_functions();
      const unsigned nbasic = n_basic_basis_functions();
      const unsigned ninternal = n_internal_dofs();

      // Initialise the containers
      DenseMatrix<double> a_matrix(nbasic, nbasic, 0.0),
        b_matrix(nbasis, nbasic, 0.0), d_matrix(nbasis, nbasis, 0.0);


      // Fill in the matrices
      monomial_to_basic_matrix(a_matrix);
      basic_to_local_matrix(b_matrix);
      local_to_global_matrix(d_matrix);

      /* Could do it like this:

      for(unsigned i=0;i<nbasis;++i)
       for(unsigned j=0;j<nbasic;++j)
         for(unsigned k=0;k<nbasis;++k)
          for(unsigned l=0;l<nbasic;++l)
            conversion_matrix(i,j)+=d_matrix(i,k)*b_matrix(k,l)*a_matrix(l,j);

       If we don't take into account sparsity of the d_matrix - but we can do
      better
      */

      // This way is about six times faster:
      // We reduce the loop over 404 elements of d_matrix to loop over 27
      // Fill in conversion matrix
      // First submatrix, identity i==k
      for (unsigned i = 0; i < 3; ++i)
      {
        for (unsigned j = 0; j < nbasic; ++j)
        {
          for (unsigned l = 0; l < nbasic; ++l)
          {
            conversion_matrix(i, j) +=
              d_matrix(i, i) * b_matrix(i, l) * a_matrix(l, j);
          }
        }
      }

      // Second submatrices, each is a two by two
      for (unsigned i = 0; i < 2; ++i)
      {
        for (unsigned j = 0; j < nbasic; ++j)
        {
          for (unsigned k = 0; k < 2; ++k)
          {
            for (unsigned l = 0; l < nbasic; ++l)
            {
              conversion_matrix(i + 3, j) +=
                d_matrix(i + 3, k + 3) * b_matrix(k + 3, l) * a_matrix(l, j);
              conversion_matrix(i + 5, j) +=
                d_matrix(i + 5, k + 5) * b_matrix(k + 5, l) * a_matrix(l, j);
              conversion_matrix(i + 7, j) +=
                d_matrix(i + 7, k + 7) * b_matrix(k + 7, l) * a_matrix(l, j);
            }
          }
        }
      }

      // Second submatrices, each is a three by three
      for (unsigned i = 0; i < 3; ++i)
      {
        for (unsigned j = 0; j < nbasic; ++j)
        {
          for (unsigned l = 0; l < nbasic; ++l)
          {
            for (unsigned k = 0; k < 2; ++k)
            {
              conversion_matrix(i + 9, j) +=
                d_matrix(i + 9, k + 9) * b_matrix(k + 9, l) * a_matrix(l, j);
              conversion_matrix(i + 12, j) +=
                d_matrix(i + 12, k + 11) * b_matrix(k + 11, l) * a_matrix(l, j);
              conversion_matrix(i + 15, j) +=
                d_matrix(i + 15, k + 13) * b_matrix(k + 13, l) * a_matrix(l, j);
            }
            // Now do the third column

            conversion_matrix(i + 9, j) +=
              d_matrix(i + 9, 15) * b_matrix(15, l) * a_matrix(l, j);
            conversion_matrix(i + 12, j) +=
              d_matrix(i + 12, 16) * b_matrix(16, l) * a_matrix(l, j);
            conversion_matrix(i + 15, j) +=
              d_matrix(i + 15, 17) * b_matrix(17, l) * a_matrix(l, j);
          }
        }
      }

      // Last submatrix, identity i==k
      for (unsigned i = 0; i < ninternal; ++i)
      {
        for (unsigned j = 0; j < nbasic; ++j)
        {
          for (unsigned l = 0; l < nbasic; ++l)
          {
            conversion_matrix(i + 18, j) +=
              d_matrix(i + 18, i + 18) * b_matrix(i + 18, l) * a_matrix(l, j);
          }
        }
      }
    }

    /// Get d2Shape with a precomputed  matrix argument
    template<unsigned BOUNDARY_ORDER>
    void BernadouElementBasis<BOUNDARY_ORDER>::d2_shape_ds2(
      const Vector<double>& s,
      Shape& psi,
      Shape& bpsi,
      DShape& dpsi,
      DShape& dbpsi,
      DShape& d2psi,
      DShape& d2bpsi,
      const DenseMatrix<double>& conversion_matrix) const
    {
      // Handy definitions
      const unsigned nbasic = n_basic_basis_functions();
      const unsigned ninternal = n_internal_dofs();

      // Initialise the containers
      Shape p7(nbasic);
      DShape dp7(nbasic, 2), d2p7(nbasic, 3);

      // Get the monomials
      full_basis_monomials(s, p7);
      dfull_basis_monomials(s, dp7);
      d2full_basis_monomials(s, d2p7);

      // Zero the shape functions
      for (unsigned i = 0; i < 3; ++i)
      {
        // Zero nodal dofs
        for (unsigned l = 0; l < 6; ++l)
        {
          // Shape functions
          psi(i, l) = 0;
          // First Derivatives
          for (unsigned alpha = 0; alpha < 2; ++alpha)
          {
            dpsi(i, l, alpha) = 0;
            // Second derivatives
            for (unsigned beta = alpha; beta < 2; ++beta)
            {
              d2psi(i, l, alpha + beta) = 0;
            }
          }
        }
      }

      for (unsigned i = 0; i < ninternal; ++i)
      {
        // Zero internal (bubble) dofs
        bpsi(i, 0) = 0;
        // First Derivatives
        for (unsigned alpha = 0; alpha < 2; ++alpha)
        {
          dbpsi(i, 0, alpha) = 0;
          // Second derivatives
          for (unsigned beta = alpha; beta < 2; ++beta)
          {
            d2bpsi(i, 0, alpha + beta) = 0;
          }
        }
      }

      // Get the permutation
      unsigned index_shift = 0;
      nodal_index_shift(index_shift);

      // Fill in the shape functions
      for (unsigned j = 0; j < nbasic; ++j)
      {
        // Nodal functional dofs
        for (unsigned i = 0; i < 3; ++i)
        {
          // Permuted i so the global index will work
          unsigned i_global = (i + index_shift) % 3;
          // Get the shape functions
          psi(i_global, 0) += conversion_matrix(i, j) * p7[j];
          // Local Derivatives
          for (unsigned alpha = 0; alpha < 2; ++alpha)
          {
            dpsi(i_global, 0, alpha) += conversion_matrix(i, j) * dp7(j, alpha);
            // Loop over the second derivatives
            for (unsigned beta = alpha; beta < 2; ++beta)
            {
              d2psi(i_global, 0, alpha + beta) +=
                conversion_matrix(i, j) * d2p7(j, alpha + beta);
            }
          }
        }
        // Nodal first derivative  dofs
        for (unsigned i = 0; i < 6; ++i)
        {
          const unsigned inode(i / 2), ideriv(i % 2);
          // Permuted i so the global dofs are consistent between elements
          unsigned i_global = (inode + index_shift) % 3;

          psi(i_global, 1 + ideriv) += conversion_matrix(i + 3, j) * p7[j];
          // Local Derivatives
          for (unsigned alpha = 0; alpha < 2; ++alpha)
          {
            dpsi(i_global, 1 + ideriv, alpha) +=
              conversion_matrix(i + 3, j) * dp7(j, alpha);
            // Loop over the second derivatives
            for (unsigned beta = alpha; beta < 2; ++beta)
            {
              d2psi(i_global, 1 + ideriv, alpha + beta) +=
                conversion_matrix(i + 3, j) * d2p7(j, alpha + beta);
            }
          }
        }
        // Nodal second derivative  dofs
        for (unsigned i = 0; i < 9; ++i)
        {
          const unsigned inode(i / 3), ideriv(i % 3);
          // Permuted i so the global dofs are consistent between elements
          unsigned i_global = (inode + index_shift) % 3;

          psi(i_global, 3 + ideriv) += conversion_matrix(i + 9, j) * p7[j];
          // Local Derivatives
          for (unsigned alpha = 0; alpha < 2; ++alpha)
          {
            dpsi(i_global, 3 + ideriv, alpha) +=
              conversion_matrix(i + 9, j) * dp7(j, alpha);
            // Loop over the second derivatives
            for (unsigned beta = alpha; beta < 2; ++beta)
            {
              d2psi(i_global, 3 + ideriv, alpha + beta) +=
                conversion_matrix(i + 9, j) * d2p7(j, alpha + beta);
            }
          }
        }
        // Now for the bubble functions
        for (unsigned i = 0; i < ninternal; ++i)
        {
          // Permuted i - to be consistent
          // unsigned i_global = (i + index_shift) % 3;
          unsigned i_global = i; // (i + index_shift) % 3; DOn't shift these!
          bpsi(i_global, 0) += conversion_matrix(i + 18, j) * p7[j];
          // Local Derivatives
          for (unsigned alpha = 0; alpha < 2; ++alpha)
          {
            dbpsi(i_global, 0, alpha) +=
              conversion_matrix(i + 18, j) * dp7(j, alpha);
            // Loop over the second derivatives
            for (unsigned beta = alpha; beta < 2; ++beta)
            {
              d2bpsi(i_global, 0, alpha + beta) +=
                conversion_matrix(i + 18, j) * d2p7(j, alpha + beta);
            }
          }
        }
      } // End of loop over shape functions
    }

    /// Get the second derivative shape functions. All of the construction
    /// happens elsewhere this just contains the matrix multiplications (which
    /// is done once
    ///  at the start to avoid doing it six times).
    template<unsigned BOUNDARY_ORDER>
    double BernadouElementBasis<BOUNDARY_ORDER>::d2_shape_dx2(
      const Vector<double>& s,
      Shape& psi,
      Shape& bpsi,
      DShape& dpsi,
      DShape& dbpsi,
      DShape& d2psi,
      DShape& d2bpsi,
      const DenseMatrix<double>& M) const
    {
// check the construction of the elements is complete
#ifdef paranoid
      if (curved_edge == none)
      {
        throw oomphliberror("the element has not been upgraded yet. did \
  you forget to set upe the curved_edge?",
                            oomph_current_function,
                            oomph_exception_location);
      }
#endif
      // Handy definitions
      const unsigned ninternal = n_internal_dofs();
      // Permute the local coordinate
      Vector<double> s_permuted(s);
      permute_shape(s_permuted);

      // Initialise shape
      DShape dpsi_ds(3, 6, 2);
      DShape d2psi_ds2(3, 6, 3);

      // Initialise bubble shape
      DShape dbpsi_ds(ninternal, 1, 2);
      DShape d2bpsi_ds2(ninternal, 1, 3);

      // Get the local derivatives
      d2_shape_ds2(
        s_permuted, psi, bpsi, dpsi_ds, dbpsi_ds, d2psi_ds2, d2bpsi_ds2, M);

      // Now transform
      DenseMatrix<double> jacobian(2, 2, 0.0), inv_jacobian(2, 2, 0.0);
      RankThreeTensor<double> hessian(2, 2, 2, 0.0), d_inv_jac_ds(2, 2, 2, 0.0);

      // Fill in Jacobian and Hessian
      get_basic_jacobian(s_permuted, jacobian);
      get_basic_hessian(s_permuted, hessian);

      // Now invert
      double det = MyShape::invert_two_by_two(jacobian, inv_jacobian);

      // Zero the shape functions
      for (unsigned i = 0; i < 3; ++i)
      {
        // Zero nodal dofs
        for (unsigned l = 0; l < 6; ++l)
        {
          // First Derivatives
          for (unsigned alpha = 0; alpha < 2; ++alpha)
          {
            dpsi(i, l, alpha) = 0;
            // Second derivatives
            for (unsigned beta = alpha; beta < 2; ++beta)
            {
              d2psi(i, l, alpha + beta) = 0;
            }
          }
        }
      }
      for (unsigned i = 0; i < ninternal; ++i)
      {
        // Zero internal (bubble) dofs
        // First Derivatives
        for (unsigned alpha = 0; alpha < 2; ++alpha)
        {
          dbpsi(i, 0, alpha) = 0;
          // Second derivatives
          for (unsigned beta = alpha; beta < 2; ++beta)
          {
            d2bpsi(i, 0, alpha + beta) = 0;
          }
        }
      }

      // Fill in values of the rank three tensor
      // (slightly Bend oomph-lib style guide for more succinct code)
      for (unsigned alpha = 0; alpha < 2; ++alpha)
      {
        for (unsigned beta = 0; beta < 2; ++beta)
        {
          for (unsigned gamma = 0; gamma < 2; ++gamma)
          {
            for (unsigned delta = 0; delta < 2; ++delta)
            {
              for (unsigned zeta = 0; zeta < 2; ++zeta)
              {
                // Fill in the rank three tensor (two contractions so loop over
                // five indices)
                d_inv_jac_ds(alpha, delta, zeta) -= inv_jacobian(alpha, beta) *
                                                    hessian(beta, gamma, zeta) *
                                                    inv_jacobian(gamma, delta);
              }
            }
          }
        }
      }


      // Now find the global derivatives at local coordinates
      for (unsigned alpha = 0; alpha < 2; ++alpha)
      {
        for (unsigned beta = 0; beta < 2; ++beta)
        {
          for (unsigned i = 0; i < ninternal; ++i)
          {
            // Transform the first derivatives of the bubble basis functions
            dbpsi(i, 0, beta) +=
              dbpsi_ds(i, 0, alpha) * inv_jacobian(alpha, beta);
          }
          for (unsigned i = 0; i < 3; ++i)
          {
            // Transform the first derivatives of the nodal basis functions
            for (unsigned l = 0; l < 6; ++l)
            {
              dpsi(i, l, beta) +=
                dpsi_ds(i, l, alpha) * inv_jacobian(alpha, beta);
            }
          }
        }
      }

      // We need six indices to evaluate the summation (bleuagh)
      for (unsigned alpha = 0; alpha < 2; ++alpha)
      {
        for (unsigned beta = alpha; beta < 2; ++beta)
        {
          for (unsigned gamma = 0; gamma < 2; ++gamma)
          {
            for (unsigned delta = 0; delta < 2; ++delta)
            {
              for (unsigned i = 0; i < ninternal; ++i)
              {
                // Fill second derivatives of the bubble basis
                d2bpsi(i, 0, alpha + beta) +=
                  d2bpsi_ds2(i, 0, gamma + delta) * inv_jacobian(gamma, alpha) *
                    inv_jacobian(delta, beta) +
                  dbpsi_ds(i, 0, gamma) * d_inv_jac_ds(gamma, beta, delta) *
                    inv_jacobian(delta, alpha);
              }
              for (unsigned i = 0; i < 3; ++i)
              {
                // Fill in second derivatives of the nodal basis
                for (unsigned l = 0; l < 6; ++l)
                {
                  d2psi(i, l, alpha + beta) +=
                    d2psi_ds2(i, l, gamma + delta) *
                      inv_jacobian(gamma, alpha) * inv_jacobian(delta, beta) +
                    dpsi_ds(i, l, gamma) * d_inv_jac_ds(gamma, beta, delta) *
                      inv_jacobian(delta, alpha);
                }
              }
            }
          }
        }
      }
      // Return the det
      return det;
    }

    /// Get the second derivative shape functions. All of the construction
    /// happens elsewhere this just contains the matrix multiplications (which
    /// is done once
    ///  at the start to avoid doing it six times).
    template<unsigned BOUNDARY_ORDER>
    double BernadouElementBasis<BOUNDARY_ORDER>::d2_shape_dx2(
      const Vector<double>& s,
      Shape& psi,
      Shape& bpsi,
      DShape& dpsi,
      DShape& dbpsi,
      DShape& d2psi,
      DShape& d2bpsi) const
    {
// check the construction of the elements is complete
#ifdef paranoid
      if (curved_edge == none)
      {
        throw oomphliberror("the element has not been upgraded yet. did \
  you forget to set upe the curved_edge?",
                            oomph_current_function,
                            oomph_exception_location);
      }
#endif
      // Handy definitions
      const unsigned ninternal = n_internal_dofs();

      // Permute the local coordinate
      Vector<double> s_permuted(s);
      permute_shape(s_permuted);

      // Initialise shape
      DShape dpsi_ds(3, 6, 2);
      DShape d2psi_ds2(3, 6, 3);

      // Initialise bubble shape
      DShape dbpsi_ds(ninternal, 1, 2);
      DShape d2bpsi_ds2(ninternal, 1, 3);

      // Get the local derivatives
      d2_shape_ds2(
        s_permuted, psi, bpsi, dpsi_ds, dbpsi_ds, d2psi_ds2, d2bpsi_ds2);

      // Now transform
      DenseMatrix<double> jacobian(2, 2, 0.0), inv_jacobian(2, 2, 0.0);
      RankThreeTensor<double> hessian(2, 2, 2, 0.0), d_inv_jac_ds(2, 2, 2, 0.0);

      // Fill in Jacobian and Hessian
      get_basic_jacobian(s_permuted, jacobian);
      get_basic_hessian(s_permuted, hessian);

      // Now invert
      double det = MyShape::invert_two_by_two(jacobian, inv_jacobian);

      // Zero the shape functions
      for (unsigned i = 0; i < 3; ++i)
      {
        // Zero nodal dofs
        for (unsigned l = 0; l < 6; ++l)
        {
          // First Derivatives
          for (unsigned alpha = 0; alpha < 2; ++alpha)
          {
            dpsi(i, l, alpha) = 0;
            // Second derivatives
            for (unsigned beta = alpha; beta < 2; ++beta)
            {
              d2psi(i, l, alpha + beta) = 0;
            }
          }
        }
      }
      for (unsigned i = 0; i < ninternal; ++i)
      {
        // Zero internal (bubble) dofs
        // First Derivatives
        for (unsigned alpha = 0; alpha < 2; ++alpha)
        {
          dbpsi(i, 0, alpha) = 0;
          // Second derivatives
          for (unsigned beta = alpha; beta < 2; ++beta)
          {
            d2bpsi(i, 0, alpha + beta) = 0;
          }
        }
      }

      // Fill in values of the rank three tensor
      // (slightly Bend oomph-lib style guide for more succinct code)
      for (unsigned alpha = 0; alpha < 2; ++alpha)
      {
        for (unsigned beta = 0; beta < 2; ++beta)
        {
          for (unsigned gamma = 0; gamma < 2; ++gamma)
          {
            for (unsigned delta = 0; delta < 2; ++delta)
            {
              for (unsigned zeta = 0; zeta < 2; ++zeta)
              {
                // Fill in the rank three tensor (two contractions so loop over
                // five indices)
                d_inv_jac_ds(alpha, delta, zeta) -= inv_jacobian(alpha, beta) *
                                                    hessian(beta, gamma, zeta) *
                                                    inv_jacobian(gamma, delta);
              }
            }
          }
        }
      }


      // Now find the global derivatives at local coordinates
      for (unsigned alpha = 0; alpha < 2; ++alpha)
      {
        for (unsigned beta = 0; beta < 2; ++beta)
        {
          for (unsigned i = 0; i < ninternal; ++i)
          {
            // Transform the first derivatives of the bubble basis functions
            dbpsi(i, 0, beta) +=
              dbpsi_ds(i, 0, alpha) * inv_jacobian(alpha, beta);
            // Transform the first derivatives of the nodal basis functions
          }
          for (unsigned i = 0; i < 3; ++i)
          {
            for (unsigned l = 0; l < 6; ++l)
            {
              dpsi(i, l, beta) +=
                dpsi_ds(i, l, alpha) * inv_jacobian(alpha, beta);
            }
          }
        }
      }

      // We need six indices to evaluate the summation (bleuagh)
      for (unsigned alpha = 0; alpha < 2; ++alpha)
      {
        for (unsigned beta = alpha; beta < 2; ++beta)
        {
          for (unsigned gamma = 0; gamma < 2; ++gamma)
          {
            for (unsigned delta = 0; delta < 2; ++delta)
            {
              for (unsigned i = 0; i < ninternal; ++i)
              {
                // Fill second derivatives of the bubble basis
                d2bpsi(i, 0, alpha + beta) +=
                  d2bpsi_ds2(i, 0, gamma + delta) * inv_jacobian(gamma, alpha) *
                    inv_jacobian(delta, beta) +
                  dbpsi_ds(i, 0, gamma) * d_inv_jac_ds(gamma, beta, delta) *
                    inv_jacobian(delta, alpha);
              }
              for (unsigned i = 0; i < 3; ++i)
              {
                // Fill in second derivatives of the nodal basis
                for (unsigned l = 0; l < 6; ++l)
                {
                  d2psi(i, l, alpha + beta) +=
                    d2psi_ds2(i, l, gamma + delta) *
                      inv_jacobian(gamma, alpha) * inv_jacobian(delta, beta) +
                    dpsi_ds(i, l, gamma) * d_inv_jac_ds(gamma, beta, delta) *
                      inv_jacobian(delta, alpha);
                }
              }
            }
          }
        }
      }
      // Return the det
      return det;
    }


    // Force build
    // template void
    // BernadouElementBasis<3>::full_basis_monomials<Vector<double> >(const
    // Vector<double>&, Vector<double>&) const;
    // template void
    // BernadouElementBasis<5>::full_basis_monomials<Vector<double> >(const
    // Vector<double>&, Vector<double>&) const;

    /// Locations of the internal dofs (in local coordinates)
    template<>
    const double BernadouElementBasis<3>::Internal_dof_knots[3][2] = {
      {1 / 2., 1 / 4.}, {1 / 4., 1 / 2.}, {1 / 4., 1 / 4.}};

    /// Locations of the internal dofs (in local coordinates)
    template<>
    const double BernadouElementBasis<5>::Internal_dof_knots[10][2] = {
      {1 / 6., 2 / 3.},
      {1 / 6., 1 / 2.},
      {1 / 6., 1 / 3.},
      {1 / 6., 1 / 6.},
      {1 / 3., 1 / 6.},
      {1 / 2., 1 / 6.},
      {2 / 3., 1 / 6.},
      {1 / 2., 1 / 3.},
      {1 / 3., 1 / 2.},
      {1 / 3., 1 / 3.}};
    // Force build
    template class BernadouElementBasis<3>;
    template class BernadouElementBasis<5>;

  } // namespace MyC1CurvedElements

} // end namespace oomph

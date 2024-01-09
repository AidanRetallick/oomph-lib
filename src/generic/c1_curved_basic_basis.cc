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
#include "c1_curved_elements.h"

// The matrix that connverts from the set of monomials to a the shape functions
// defined on the basic element
// this output is directly from the mathematica script that generated the basic
// shape functions.

// We could possibly change this out for a static DenseMatrix which has the LU
// decomposition stored - might give smaller round off and be of comparably
// speed HERE

// We could definitely store this as a static
namespace oomph
{
  namespace MyC1CurvedElements
  {
    /// Precomputed basis polynomials for the 3rd order boundary representation
    template<>
    void BernadouElementBasis<3>::full_basic_polynomials(
      const Vector<double>& s, Shape& phi) const
    {
      // This will be replaced eventually - in favour of the explicit shape
      // functions
      // Get the number of basic basis function
      const unsigned nbasic = n_basic_basis_functions();
      // Initialise the p7 shape functions
      Shape p7(nbasic);
      full_basis_monomials(s, p7);
      // Initialise the basic association matrix
      DenseMatrix<double> a_matrix(nbasic, nbasic, 0.0);
      monomial_to_basic_matrix(a_matrix);
      // Loop over the 36 Basis functions
      for (unsigned i = 0; i < nbasic; ++i)
      {
        // Zero the shape functions
        phi(i) = 0.0;
        // Zero nodal dofs
        for (unsigned j = 0; j < nbasic; ++j)
        {
          phi(i) += a_matrix(i, j) * p7(j);
        }
      }
    }

    /// Precomputed dbasis polynomials for the 3rd order boundary representation
    template<>
    void BernadouElementBasis<3>::dfull_basic_polynomials(
      const Vector<double>& s, DShape& dphi) const
    {
      // This will be replaced eventually - in favour of the explicit shape
      // functions
      // Get the number of basic basis function
      const unsigned nbasic = n_basic_basis_functions();
      // Initialise the p7 shape functions
      DShape dp7(nbasic, 2);
      dfull_basis_monomials(s, dp7);
      // Initialise the basic association matrix
      DenseMatrix<double> a_matrix(nbasic, nbasic, 0.0);
      monomial_to_basic_matrix(a_matrix);
      // Loop over the 36 Basis functions
      for (unsigned i = 0; i < nbasic; ++i)
      {
        for (unsigned alpha = 0; alpha < 2; ++alpha)
        {
          // Zero the shape functions
          dphi(i, alpha) = 0.0;
          // Zero nodal dofs
          for (unsigned j = 0; j < nbasic; ++j)
          {
            dphi(i, alpha) += a_matrix(i, j) * dp7(j, alpha);
          }
        }
      }
    }

    /// Precomputed d2basis polynomials for the 3rd order boundary
    /// representation
    template<>
    void BernadouElementBasis<3>::d2full_basic_polynomials(
      const Vector<double>& s, DShape& d2phi) const
    {
      // This will be replaced eventually - in favour of the explicit shape
      // functions
      // Get the number of basic basis function
      const unsigned nbasic = n_basic_basis_functions();
      // Initialise the p7 shape functions
      DShape d2p7(nbasic, 3);
      d2full_basis_monomials(s, d2p7);
      // Initialise the basic association matrix
      DenseMatrix<double> a_matrix(nbasic, nbasic, 0.0);
      monomial_to_basic_matrix(a_matrix);
      // Loop over the 36 Basis functions
      for (unsigned i = 0; i < nbasic; ++i)
      {
        for (unsigned alpha = 0; alpha < 3; ++alpha)
        {
          // Zero the shape functions
          d2phi(i, alpha) = 0.0;
          // Zero nodal dofs
          for (unsigned j = 0; j < nbasic; ++j)
          {
            d2phi(i, alpha) += a_matrix(i, j) * d2p7(j, alpha);
          }
        }
      }
    }

    /// Precomputed basis polynomials for the 5th order boundary representation
    template<>
    void BernadouElementBasis<5>::full_basic_polynomials(
      const Vector<double>& s_basic, Shape& phi) const
    {
      // For convenience
      const double s = s_basic[0], t = s_basic[1];
      // Precompute
      const double s_p_2 = pow(s, 2), s_p_3 = pow(s, 3), s_p_4 = pow(s, 4),
                   s_p_5 = pow(s, 5), s_p_6 = pow(s, 6), s_p_7 = pow(s, 7);
      const double t_p_2 = pow(t, 2), t_p_3 = pow(t, 3), t_p_4 = pow(t, 4),
                   t_p_5 = pow(t, 5), t_p_6 = pow(t, 6), t_p_7 = pow(t, 7);
      // Now fill in (automatically generated code
      phi[0] =
        (s_p_2 *
         (598360 * s - 8130732 * s_p_2 + 40143300 * s_p_3 - 94869252 * s_p_4 +
          116223876 * s_p_5 - 71313696 * s_p_6 + 17352144 * s_p_7 -
          10246440 * t_p_2 + 126762915 * s * t_p_2 - 538504848 * s_p_2 * t_p_2 +
          994544973 * s_p_3 * t_p_2 - 823436604 * s_p_4 * t_p_2 +
          250880004 * s_p_5 * t_p_2 + 22167939 * t_p_3 - 247491396 * s * t_p_3 +
          924349995 * s_p_2 * t_p_3 - 1297023192 * s_p_3 * t_p_3 +
          592038396 * s_p_4 * t_p_3 - 12586266 * t_p_4 + 92584809 * s * t_p_4 -
          342910908 * s_p_2 * t_p_4 + 327894480 * s_p_3 * t_p_4 +
          5508243 * t_p_5 + 42791328 * s * t_p_5 - 33765336 * s_p_2 * t_p_5 -
          9799704 * t_p_6 - 15545196 * s * t_p_6 + 4956228 * t_p_7)) /
        4000.;

      phi[1] = (t_p_2 *
                (598360 * t - 10246440 * s_p_2 + 126762915 * t * s_p_2 +
                 22167939 * s_p_3 - 247491396 * t * s_p_3 - 12586266 * s_p_4 +
                 92584809 * t * s_p_4 + 5508243 * s_p_5 + 42791328 * t * s_p_5 -
                 9799704 * s_p_6 - 15545196 * t * s_p_6 + 4956228 * s_p_7 -
                 8130732 * t_p_2 - 538504848 * s_p_2 * t_p_2 +
                 924349995 * s_p_3 * t_p_2 - 342910908 * s_p_4 * t_p_2 -
                 33765336 * s_p_5 * t_p_2 + 40143300 * t_p_3 +
                 994544973 * s_p_2 * t_p_3 - 1297023192 * s_p_3 * t_p_3 +
                 327894480 * s_p_4 * t_p_3 - 94869252 * t_p_4 -
                 823436604 * s_p_2 * t_p_4 + 592038396 * s_p_3 * t_p_4 +
                 116223876 * t_p_5 + 250880004 * s_p_2 * t_p_5 -
                 71313696 * t_p_6 + 17352144 * t_p_7)) /
               4000.;

      phi[2] =
        -((-2000 - 4000 * s - 4000 * t - 12000 * s * t - 6000 * s_p_2 -
           24000 * t * s_p_2 + 1978086 * s_p_3 + 3932172 * t * s_p_3 -
           11934864 * s_p_4 - 19937556 * t * s_p_4 + 26368362 * s_p_5 +
           32799168 * t * s_p_5 - 25075656 * s_p_6 - 17352144 * t * s_p_6 +
           8676072 * s_p_7 - 6000 * t_p_2 - 24000 * s * t_p_2 -
           16927218 * s_p_2 * t_p_2 + 53140455 * s_p_3 * t_p_2 -
           15361002 * s_p_4 * t_p_2 - 31939920 * s_p_5 * t_p_2 +
           1978086 * t_p_3 + 3932172 * s * t_p_3 + 53140455 * s_p_2 * t_p_3 -
           197161452 * s_p_3 * t_p_3 + 117457452 * s_p_4 * t_p_3 -
           11934864 * t_p_4 - 19937556 * s * t_p_4 - 15361002 * s_p_2 * t_p_4 +
           117457452 * s_p_3 * t_p_4 + 26368362 * t_p_5 + 32799168 * s * t_p_5 -
           31939920 * s_p_2 * t_p_5 - 25075656 * t_p_6 - 17352144 * s * t_p_6 +
           8676072 * t_p_7) *
          pow(-1 + s + t, 2)) /
        2000.;

      phi[3] =
        -(s_p_2 *
          (6280 * s - 85636 * s_p_2 + 425100 * s_p_3 - 1011996 * s_p_4 +
           1250748 * s_p_5 - 775008 * s_p_6 + 190512 * s_p_7 - 103738 * t_p_2 +
           1284763 * s * t_p_2 - 5469804 * s_p_2 * t_p_2 +
           10136979 * s_p_3 * t_p_2 - 8434692 * s_p_4 * t_p_2 +
           2586492 * s_p_5 * t_p_2 + 222565 * t_p_3 - 2481858 * s * t_p_3 +
           9277785 * s_p_2 * t_p_3 - 13035816 * s_p_3 * t_p_3 +
           5967108 * s_p_4 * t_p_3 - 127188 * t_p_4 + 919827 * s * t_p_4 -
           3415284 * s_p_2 * t_p_4 + 3259440 * s_p_3 * t_p_4 + 59229 * t_p_5 +
           429624 * s * t_p_5 - 331128 * s_p_2 * t_p_5 - 102384 * t_p_6 -
           158436 * s * t_p_6 + 51516 * t_p_7)) /
        400.;

      phi[4] =
        -(t * s_p_2 *
          (6680 - 97796 * s - 9488 * t + 115031 * s * t + 536936 * s_p_2 -
           418950 * t * s_p_2 - 1445940 * s_p_3 + 687627 * t * s_p_3 +
           2043000 * s_p_4 - 601344 * t * s_p_4 - 1450224 * s_p_5 +
           227124 * t * s_p_5 + 406944 * s_p_6 + 55323 * t_p_2 -
           638370 * s * t_p_2 + 1640007 * s_p_2 * t_p_2 -
           1194588 * s_p_3 * t_p_2 + 201204 * s_p_4 * t_p_2 - 122310 * t_p_3 +
           1602207 * s * t_p_3 - 3564648 * s_p_2 * t_p_3 +
           1718496 * s_p_3 * t_p_3 + 18927 * t_p_4 - 1133028 * s * t_p_4 +
           1737936 * s_p_2 * t_p_4 + 102384 * t_p_5 + 158436 * s * t_p_5 -
           51516 * t_p_6)) /
        400.;

      phi[5] = (s * t_p_2 *
                (-6680 + 9488 * s + 97796 * t - 115031 * s * t - 55323 * s_p_2 +
                 638370 * t * s_p_2 + 122310 * s_p_3 - 1602207 * t * s_p_3 -
                 18927 * s_p_4 + 1133028 * t * s_p_4 - 102384 * s_p_5 -
                 158436 * t * s_p_5 + 51516 * s_p_6 - 536936 * t_p_2 +
                 418950 * s * t_p_2 - 1640007 * s_p_2 * t_p_2 +
                 3564648 * s_p_3 * t_p_2 - 1737936 * s_p_4 * t_p_2 +
                 1445940 * t_p_3 - 687627 * s * t_p_3 +
                 1194588 * s_p_2 * t_p_3 - 1718496 * s_p_3 * t_p_3 -
                 2043000 * t_p_4 + 601344 * s * t_p_4 - 201204 * s_p_2 * t_p_4 +
                 1450224 * t_p_5 - 227124 * s * t_p_5 - 406944 * t_p_6)) /
               400.;

      phi[6] =
        -(t_p_2 *
          (6280 * t - 103738 * s_p_2 + 1284763 * t * s_p_2 + 222565 * s_p_3 -
           2481858 * t * s_p_3 - 127188 * s_p_4 + 919827 * t * s_p_4 +
           59229 * s_p_5 + 429624 * t * s_p_5 - 102384 * s_p_6 -
           158436 * t * s_p_6 + 51516 * s_p_7 - 85636 * t_p_2 -
           5469804 * s_p_2 * t_p_2 + 9277785 * s_p_3 * t_p_2 -
           3415284 * s_p_4 * t_p_2 - 331128 * s_p_5 * t_p_2 + 425100 * t_p_3 +
           10136979 * s_p_2 * t_p_3 - 13035816 * s_p_3 * t_p_3 +
           3259440 * s_p_4 * t_p_3 - 1011996 * t_p_4 - 8434692 * s_p_2 * t_p_4 +
           5967108 * s_p_3 * t_p_4 + 1250748 * t_p_5 + 2586492 * s_p_2 * t_p_5 -
           775008 * t_p_6 + 190512 * t_p_7)) /
        400.;

      phi[7] =
        -(s *
          (-200 - 400 * s - 400 * t - 1200 * s * t + 26178 * s_p_2 +
           51156 * t * s_p_2 - 142272 * s_p_3 - 233388 * t * s_p_3 +
           300726 * s_p_4 + 368064 * t * s_p_4 - 279288 * s_p_5 -
           190512 * t * s_p_5 + 95256 * s_p_6 + 29518 * t_p_2 -
           218907 * s * t_p_2 + 355743 * s_p_2 * t_p_2 +
           278154 * s_p_3 * t_p_2 - 563760 * s_p_4 * t_p_2 - 181350 * t_p_3 +
           1104822 * s * t_p_3 - 1916298 * s_p_2 * t_p_3 +
           722196 * s_p_3 * t_p_3 + 448020 * t_p_4 - 1782000 * s * t_p_4 +
           1684800 * s_p_2 * t_p_4 - 495720 * t_p_5 + 939600 * s * t_p_5 +
           203472 * t_p_6) *
          pow(-1 + s + t, 2)) /
        200.;

      phi[8] = -(t *
                 (-200 - 400 * s - 400 * t - 1200 * s * t + 29518 * s_p_2 -
                  218907 * t * s_p_2 - 181350 * s_p_3 + 1104822 * t * s_p_3 +
                  448020 * s_p_4 - 1782000 * t * s_p_4 - 495720 * s_p_5 +
                  939600 * t * s_p_5 + 203472 * s_p_6 + 26178 * t_p_2 +
                  51156 * s * t_p_2 + 355743 * s_p_2 * t_p_2 -
                  1916298 * s_p_3 * t_p_2 + 1684800 * s_p_4 * t_p_2 -
                  142272 * t_p_3 - 233388 * s * t_p_3 + 278154 * s_p_2 * t_p_3 +
                  722196 * s_p_3 * t_p_3 + 300726 * t_p_4 + 368064 * s * t_p_4 -
                  563760 * s_p_2 * t_p_4 - 279288 * t_p_5 - 190512 * s * t_p_5 +
                  95256 * t_p_6) *
                 pow(-1 + s + t, 2)) /
               200.;

      phi[9] =
        (s_p_2 *
         (40 * s - 548 * s_p_2 + 2740 * s_p_3 - 6588 * s_p_4 + 8244 * s_p_5 -
          5184 * s_p_6 + 1296 * s_p_7 - 628 * t_p_2 + 7783 * s * t_p_2 -
          33192 * s_p_2 * t_p_2 + 61677 * s_p_3 * t_p_2 -
          51516 * s_p_4 * t_p_2 + 15876 * s_p_5 * t_p_2 + 1339 * t_p_3 -
          14904 * s * t_p_3 + 55755 * s_p_2 * t_p_3 - 78408 * s_p_3 * t_p_3 +
          35964 * s_p_4 * t_p_3 - 774 * t_p_4 + 5481 * s * t_p_4 -
          20412 * s_p_2 * t_p_4 + 19440 * s_p_3 * t_p_4 + 387 * t_p_5 +
          2592 * s * t_p_5 - 1944 * s_p_2 * t_p_5 - 648 * t_p_6 -
          972 * s * t_p_6 + 324 * t_p_7)) /
        80.;

      phi[10] =
        (t * s_p_2 *
         (40 - 588 * s - 50 * t + 599 * s * t + 3248 * s_p_2 -
          2070 * t * s_p_2 - 8820 * s_p_3 + 3141 * t * s_p_3 + 12600 * s_p_4 -
          2592 * t * s_p_4 - 9072 * s_p_5 + 972 * t * s_p_5 + 2592 * s_p_6 +
          325 * t_p_2 - 3780 * s * t_p_2 + 9621 * s_p_2 * t_p_2 -
          6804 * s_p_3 * t_p_2 + 972 * s_p_4 * t_p_2 - 720 * t_p_3 +
          9621 * s * t_p_3 - 21384 * s_p_2 * t_p_3 + 10368 * s_p_3 * t_p_3 +
          81 * t_p_4 - 6804 * s * t_p_4 + 10368 * s_p_2 * t_p_4 + 648 * t_p_5 +
          972 * s * t_p_5 - 324 * t_p_6)) /
        40.;

      phi[11] =
        (s_p_2 * t_p_2 *
         (10 - 51 * s - 325 * t + 3780 * s * t - 630 * s_p_2 -
          9621 * t * s_p_2 + 2979 * s_p_3 + 6804 * t * s_p_3 - 3888 * s_p_4 -
          972 * t * s_p_4 + 1620 * s_p_5 + 720 * t_p_2 - 9621 * s * t_p_2 +
          21384 * s_p_2 * t_p_2 - 10368 * s_p_3 * t_p_2 - 81 * t_p_3 +
          6804 * s * t_p_3 - 10368 * s_p_2 * t_p_3 - 648 * t_p_4 -
          972 * s * t_p_4 + 324 * t_p_5)) /
        80.;

      phi[12] =
        (s_p_2 * t_p_2 *
         (10 - 325 * s - 51 * t + 3780 * s * t + 720 * s_p_2 -
          9621 * t * s_p_2 - 81 * s_p_3 + 6804 * t * s_p_3 - 648 * s_p_4 -
          972 * t * s_p_4 + 324 * s_p_5 - 630 * t_p_2 - 9621 * s * t_p_2 +
          21384 * s_p_2 * t_p_2 - 10368 * s_p_3 * t_p_2 + 2979 * t_p_3 +
          6804 * s * t_p_3 - 10368 * s_p_2 * t_p_3 - 3888 * t_p_4 -
          972 * s * t_p_4 + 1620 * t_p_5)) /
        80.;

      phi[13] =
        -(s * t_p_2 *
          (-40 + 50 * s + 588 * t - 599 * s * t - 325 * s_p_2 +
           3780 * t * s_p_2 + 720 * s_p_3 - 9621 * t * s_p_3 - 81 * s_p_4 +
           6804 * t * s_p_4 - 648 * s_p_5 - 972 * t * s_p_5 + 324 * s_p_6 -
           3248 * t_p_2 + 2070 * s * t_p_2 - 9621 * s_p_2 * t_p_2 +
           21384 * s_p_3 * t_p_2 - 10368 * s_p_4 * t_p_2 + 8820 * t_p_3 -
           3141 * s * t_p_3 + 6804 * s_p_2 * t_p_3 - 10368 * s_p_3 * t_p_3 -
           12600 * t_p_4 + 2592 * s * t_p_4 - 972 * s_p_2 * t_p_4 +
           9072 * t_p_5 - 972 * s * t_p_5 - 2592 * t_p_6)) /
        40.;

      phi[14] =
        (t_p_2 *
         (40 * t - 628 * s_p_2 + 7783 * t * s_p_2 + 1339 * s_p_3 -
          14904 * t * s_p_3 - 774 * s_p_4 + 5481 * t * s_p_4 + 387 * s_p_5 +
          2592 * t * s_p_5 - 648 * s_p_6 - 972 * t * s_p_6 + 324 * s_p_7 -
          548 * t_p_2 - 33192 * s_p_2 * t_p_2 + 55755 * s_p_3 * t_p_2 -
          20412 * s_p_4 * t_p_2 - 1944 * s_p_5 * t_p_2 + 2740 * t_p_3 +
          61677 * s_p_2 * t_p_3 - 78408 * s_p_3 * t_p_3 +
          19440 * s_p_4 * t_p_3 - 6588 * t_p_4 - 51516 * s_p_2 * t_p_4 +
          35964 * s_p_3 * t_p_4 + 8244 * t_p_5 + 15876 * s_p_2 * t_p_5 -
          5184 * t_p_6 + 1296 * t_p_7)) /
        80.;

      phi[15] =
        -((-1 + 3 * s) * (-1 + 6 * s) * s_p_2 *
          (-20 + 74 * s - 40 * t + 108 * s * t - 90 * s_p_2 - 72 * t * s_p_2 +
           36 * s_p_3 + 279 * t_p_2 - 360 * s * t_p_2 - 234 * t_p_3) *
          pow(-1 + s + t, 2)) /
        40.;

      phi[16] = -(s * t * (-1 + 2 * s + 2 * t) * (-2 + 3 * s + 3 * t) *
                  (-1 + 3 * s + 3 * t) * (-5 + 6 * s + 6 * t) *
                  (-1 + 6 * s + 6 * t) * pow(-1 + s + t, 2)) /
                10.;

      phi[17] =
        -((-1 + 3 * t) * (-1 + 6 * t) * t_p_2 *
          (-20 - 40 * s + 74 * t + 108 * s * t + 279 * s_p_2 - 360 * t * s_p_2 -
           234 * s_p_3 - 90 * t_p_2 - 72 * s * t_p_2 + 36 * t_p_3) *
          pow(-1 + s + t, 2)) /
        40.;

      phi[18] = -16 * s * (-1 + 3 * t) * (-2 + 3 * s + 3 * t) * (-1 + 6 * t) *
                (-5 + 6 * s + 6 * t) * t_p_2 * pow(-1 + s + t, 2);

      phi[19] = -16 * (-1 + 3 * s) * (-1 + 6 * s) * t * (-2 + 3 * s + 3 * t) *
                (-5 + 6 * s + 6 * t) * s_p_2 * pow(-1 + s + t, 2);

      phi[20] = 8 * (-1 + 3 * s) * (-1 + 6 * s) * (-1 + s + t) * (-1 + 3 * t) *
                (-1 + 6 * t) * s_p_2 * t_p_2 * sqrt(2);

      phi[21] = (-1944 * (-1 + 3 * t) * (-1 + 6 * t) * t_p_2 *
                 (16 * t + 32 * s * t - 279 * s_p_2 + 360 * t * s_p_2 +
                  234 * s_p_3 - 40 * t_p_2 - 48 * s * t_p_2 + 24 * t_p_3) *
                 pow(-1 + s + t, 2)) /
                125.;

      phi[22] = (243 * (-1 + 3 * t) * (-1 + 6 * t) * t_p_2 *
                 (10 * t + 20 * s * t - 183 * s_p_2 + 168 * t * s_p_2 +
                  234 * s_p_3 - 22 * t_p_2 - 24 * s * t_p_2 + 12 * t_p_3) *
                 pow(-1 + s + t, 2)) /
                16.;

      phi[23] = (-81 * (-1 + 6 * t) * t_p_2 *
                 (-60 * t - 120 * s * t + 1003 * s_p_2 - 1809 * t * s_p_2 -
                  2898 * s_p_3 + 3078 * t * s_p_3 + 2016 * s_p_4 + 222 * t_p_2 +
                  324 * s * t_p_2 + 648 * s_p_2 * t_p_2 - 270 * t_p_3 -
                  216 * s * t_p_3 + 108 * t_p_4) *
                 pow(-1 + s + t, 2)) /
                16.;

      phi[24] =
        (324 * t_p_2 *
         (480 * t + 960 * s * t - 5929 * s_p_2 + 22458 * t * s_p_2 +
          25473 * s_p_3 - 77868 * t * s_p_3 - 35388 * s_p_4 +
          57888 * t * s_p_4 + 16092 * s_p_5 - 3216 * t_p_2 - 5472 * s * t_p_2 -
          17172 * s_p_2 * t_p_2 + 48600 * s_p_3 * t_p_2 + 7488 * t_p_3 +
          9504 * s * t_p_3 - 2592 * s_p_2 * t_p_3 - 7344 * t_p_4 -
          5184 * s * t_p_4 + 2592 * t_p_5) *
         pow(-1 + s + t, 2)) /
        125.;

      phi[25] =
        (324 * s_p_2 *
         (480 * s + 960 * s * t - 3216 * s_p_2 - 5472 * t * s_p_2 +
          7488 * s_p_3 + 9504 * t * s_p_3 - 7344 * s_p_4 - 5184 * t * s_p_4 +
          2592 * s_p_5 - 5929 * t_p_2 + 22458 * s * t_p_2 -
          17172 * s_p_2 * t_p_2 - 2592 * s_p_3 * t_p_2 + 25473 * t_p_3 -
          77868 * s * t_p_3 + 48600 * s_p_2 * t_p_3 - 35388 * t_p_4 +
          57888 * s * t_p_4 + 16092 * t_p_5) *
         pow(-1 + s + t, 2)) /
        125.;

      phi[26] = (-81 * (-1 + 6 * s) * s_p_2 *
                 (-60 * s - 120 * s * t + 222 * s_p_2 + 324 * t * s_p_2 -
                  270 * s_p_3 - 216 * t * s_p_3 + 108 * s_p_4 + 1003 * t_p_2 -
                  1809 * s * t_p_2 + 648 * s_p_2 * t_p_2 - 2898 * t_p_3 +
                  3078 * s * t_p_3 + 2016 * t_p_4) *
                 pow(-1 + s + t, 2)) /
                16.;

      phi[27] = (243 * (-1 + 3 * s) * (-1 + 6 * s) * s_p_2 *
                 (10 * s + 20 * s * t - 22 * s_p_2 - 24 * t * s_p_2 +
                  12 * s_p_3 - 183 * t_p_2 + 168 * s * t_p_2 + 234 * t_p_3) *
                 pow(-1 + s + t, 2)) /
                16.;

      phi[28] = (-1944 * (-1 + 3 * s) * (-1 + 6 * s) * s_p_2 *
                 (16 * s + 32 * s * t - 40 * s_p_2 - 48 * t * s_p_2 +
                  24 * s_p_3 - 279 * t_p_2 + 360 * s * t_p_2 + 234 * t_p_3) *
                 pow(-1 + s + t, 2)) /
                125.;

      phi[29] =
        (-162 * s_p_2 * t_p_2 *
         (-518 + 7483 * s - 1313 * t + 14000 * s * t - 40700 * s_p_2 -
          21891 * t * s_p_2 + 97059 * s_p_3 - 20016 * t * s_p_3 -
          100800 * s_p_4 + 29700 * t * s_p_4 + 37476 * s_p_5 + 4540 * t_p_2 -
          56751 * s * t_p_2 + 126144 * s_p_2 * t_p_2 - 58968 * s_p_3 * t_p_2 -
          1161 * t_p_3 + 40464 * s * t_p_3 - 63288 * s_p_2 * t_p_3 -
          3168 * t_p_4 - 5292 * s * t_p_4 + 1620 * t_p_5)) /
        125.;

      phi[30] =
        (-81 * s_p_2 * t_p_2 *
         (-274 + 3167 * s + 2021 * t - 22688 * s * t - 11236 * s_p_2 +
          72087 * t * s_p_2 + 15867 * s_p_3 - 82044 * t * s_p_3 - 9144 * s_p_4 +
          30564 * t * s_p_4 + 1620 * s_p_5 - 3382 * t_p_2 + 37767 * s * t_p_2 -
          94392 * s_p_2 * t_p_2 + 58752 * s_p_3 * t_p_2 + 1167 * t_p_3 -
          20268 * s * t_p_3 + 33696 * s_p_2 * t_p_3 + 1008 * t_p_4 +
          2052 * s * t_p_4 - 540 * t_p_5)) /
        32.;

      phi[31] =
        (81 * s_p_2 * t_p_2 *
         (274 - 2021 * s - 3167 * t + 22688 * s * t + 3382 * s_p_2 -
          37767 * t * s_p_2 - 1167 * s_p_3 + 20268 * t * s_p_3 - 1008 * s_p_4 -
          2052 * t * s_p_4 + 540 * s_p_5 + 11236 * t_p_2 - 72087 * s * t_p_2 +
          94392 * s_p_2 * t_p_2 - 33696 * s_p_3 * t_p_2 - 15867 * t_p_3 +
          82044 * s * t_p_3 - 58752 * s_p_2 * t_p_3 + 9144 * t_p_4 -
          30564 * s * t_p_4 - 1620 * t_p_5)) /
        32.;

      phi[32] =
        (-162 * s_p_2 * t_p_2 *
         (-518 - 1313 * s + 7483 * t + 14000 * s * t + 4540 * s_p_2 -
          56751 * t * s_p_2 - 1161 * s_p_3 + 40464 * t * s_p_3 - 3168 * s_p_4 -
          5292 * t * s_p_4 + 1620 * s_p_5 - 40700 * t_p_2 - 21891 * s * t_p_2 +
          126144 * s_p_2 * t_p_2 - 63288 * s_p_3 * t_p_2 + 97059 * t_p_3 -
          20016 * s * t_p_3 - 58968 * s_p_2 * t_p_3 - 100800 * t_p_4 +
          29700 * s * t_p_4 + 37476 * t_p_5)) /
        125.;

      phi[33] = (-648 * s * (-1 + 2 * t) * (-2 + 3 * t) * (-1 + 3 * t) *
                 (-1 + 6 * t) * t_p_2 * pow(-1 + s + t, 2)) /
                25.;

      phi[34] = (81 * s * (-1 + 2 * t) * (-1 + 3 * t) * (-1 + 6 * t) *
                 (-5 + 6 * s + 6 * t) * t_p_2 * pow(-1 + s + t, 2)) /
                4.;

      phi[35] =
        (81 * s * (-1 + 2 * s + 2 * t) * (-2 + 3 * s + 3 * t) * (-1 + 6 * t) *
         (-5 + 6 * s + 6 * t) * t_p_2 * pow(-1 + s + t, 2)) /
        4.;

      phi[36] = (-648 * s * (-1 + 2 * s + 2 * t) * (-2 + 3 * s + 3 * t) *
                 (-1 + 3 * s + 3 * t) * (-5 + 6 * s + 6 * t) * t_p_2 *
                 pow(-1 + s + t, 2)) /
                25.;

      phi[37] = (-648 * t * (-1 + 2 * s + 2 * t) * (-2 + 3 * s + 3 * t) *
                 (-1 + 3 * s + 3 * t) * (-5 + 6 * s + 6 * t) * s_p_2 *
                 pow(-1 + s + t, 2)) /
                25.;

      phi[38] =
        (81 * (-1 + 6 * s) * t * (-1 + 2 * s + 2 * t) * (-2 + 3 * s + 3 * t) *
         (-5 + 6 * s + 6 * t) * s_p_2 * pow(-1 + s + t, 2)) /
        4.;

      phi[39] = (81 * (-1 + 2 * s) * (-1 + 3 * s) * (-1 + 6 * s) * t *
                 (-5 + 6 * s + 6 * t) * s_p_2 * pow(-1 + s + t, 2)) /
                4.;

      phi[40] = (-648 * (-1 + 2 * s) * (-2 + 3 * s) * (-1 + 3 * s) *
                 (-1 + 6 * s) * t * s_p_2 * pow(-1 + s + t, 2)) /
                25.;

      phi[41] = (324 * (-1 + 2 * s) * (-2 + 3 * s) * (-1 + 3 * s) *
                 (-1 + 6 * s) * (-1 + s + t) * s_p_2 * t_p_2 * sqrt(2)) /
                25.;

      phi[42] =
        (81 * (-1 + 2 * s) * (-1 + 3 * s) * (-1 + 6 * s) * (-1 + s + t) *
         (-1 + 6 * t) * s_p_2 * t_p_2 * pow(sqrt(2), -1)) /
        4.;

      phi[43] =
        (81 * (-1 + 6 * s) * (-1 + s + t) * (-1 + 2 * t) * (-1 + 3 * t) *
         (-1 + 6 * t) * s_p_2 * t_p_2 * pow(sqrt(2), -1)) /
        4.;

      phi[44] = (324 * (-1 + s + t) * (-1 + 2 * t) * (-2 + 3 * t) *
                 (-1 + 3 * t) * (-1 + 6 * t) * s_p_2 * t_p_2 * sqrt(2)) /
                25.;

      phi[45] = 2916 * (-1 + 2 * t) * (-1 + 3 * t) * (-1 + 6 * t) * s_p_2 *
                t_p_2 * pow(-1 + s + t, 2);

      phi[46] = -1296 * (-1 + 3 * t) * (-1 + 6 * t) * (-5 + 6 * s + 6 * t) *
                s_p_2 * t_p_2 * pow(-1 + s + t, 2);

      phi[47] = 1296 * (-2 + 3 * s + 3 * t) * (-1 + 6 * t) *
                (-5 + 6 * s + 6 * t) * s_p_2 * t_p_2 * pow(-1 + s + t, 2);

      phi[48] = -2916 * (-1 + 2 * s + 2 * t) * (-2 + 3 * s + 3 * t) *
                (-5 + 6 * s + 6 * t) * s_p_2 * t_p_2 * pow(-1 + s + t, 2);

      phi[49] = 1296 * (-1 + 6 * s) * (-2 + 3 * s + 3 * t) *
                (-5 + 6 * s + 6 * t) * s_p_2 * t_p_2 * pow(-1 + s + t, 2);

      phi[50] = -1296 * (-1 + 3 * s) * (-1 + 6 * s) * (-5 + 6 * s + 6 * t) *
                s_p_2 * t_p_2 * pow(-1 + s + t, 2);

      phi[51] = 2916 * (-1 + 2 * s) * (-1 + 3 * s) * (-1 + 6 * s) * s_p_2 *
                t_p_2 * pow(-1 + s + t, 2);

      phi[52] = 1296 * (-1 + 3 * s) * (-1 + 6 * s) * (-1 + 6 * t) * s_p_2 *
                t_p_2 * pow(-1 + s + t, 2);

      phi[53] = 1296 * (-1 + 6 * s) * (-1 + 3 * t) * (-1 + 6 * t) * s_p_2 *
                t_p_2 * pow(-1 + s + t, 2);

      phi[54] = -729 * (-1 + 6 * s) * (-1 + 6 * t) * (-5 + 6 * s + 6 * t) *
                s_p_2 * t_p_2 * pow(-1 + s + t, 2);
    }

    /// Precomputed dbasis polynomials for the 5th order boundary representation
    template<>
    void BernadouElementBasis<5>::dfull_basic_polynomials(
      const Vector<double>& s_basic, DShape& dphi) const
    {
      // For convenience
      const double s = s_basic[0], t = s_basic[1];

      // Precompute
      const double s_p_2 = pow(s, 2), s_p_3 = pow(s, 3), s_p_4 = pow(s, 4),
                   s_p_5 = pow(s, 5), s_p_6 = pow(s, 6), s_p_7 = pow(s, 7);
      const double t_p_2 = pow(t, 2), t_p_3 = pow(t, 3), t_p_4 = pow(t, 4),
                   t_p_5 = pow(t, 5), t_p_6 = pow(t, 6), t_p_7 = pow(t, 7);
      // Now fill in (automatically generated code)
      dphi(0, 0) =
        (3 * s *
         (598360 * s - 10840976 * s_p_2 + 66905500 * s_p_3 - 189738504 * s_p_4 +
          271189044 * s_p_5 - 190169856 * s_p_6 + 52056432 * s_p_7 -
          6830960 * t_p_2 + 126762915 * s * t_p_2 - 718006464 * s_p_2 * t_p_2 +
          1657574955 * s_p_3 * t_p_2 - 1646873208 * s_p_4 * t_p_2 +
          585386676 * s_p_5 * t_p_2 + 14778626 * t_p_3 - 247491396 * s * t_p_3 +
          1232466660 * s_p_2 * t_p_3 - 2161705320 * s_p_3 * t_p_3 +
          1184076792 * s_p_4 * t_p_3 - 8390844 * t_p_4 + 92584809 * s * t_p_4 -
          457214544 * s_p_2 * t_p_4 + 546490800 * s_p_3 * t_p_4 +
          3672162 * t_p_5 + 42791328 * s * t_p_5 - 45020448 * s_p_2 * t_p_5 -
          6533136 * t_p_6 - 15545196 * s * t_p_6 + 3304152 * t_p_7)) /
        4000.;

      dphi(0, 1) =
        (3 * t * s_p_2 *
         (-6830960 + 84508610 * s + 22167939 * t - 247491396 * s * t -
          359003232 * s_p_2 + 924349995 * t * s_p_2 + 663029982 * s_p_3 -
          1297023192 * t * s_p_3 - 548957736 * s_p_4 + 592038396 * t * s_p_4 +
          167253336 * s_p_5 - 16781688 * t_p_2 + 123446412 * s * t_p_2 -
          457214544 * s_p_2 * t_p_2 + 437192640 * s_p_3 * t_p_2 +
          9180405 * t_p_3 + 71318880 * s * t_p_3 - 56275560 * s_p_2 * t_p_3 -
          19599408 * t_p_4 - 31090392 * s * t_p_4 + 11564532 * t_p_5)) /
        4000.;

      dphi(1, 0) =
        (3 * s * t_p_2 *
         (-6830960 + 22167939 * s + 84508610 * t - 247491396 * s * t -
          16781688 * s_p_2 + 123446412 * t * s_p_2 + 9180405 * s_p_3 +
          71318880 * t * s_p_3 - 19599408 * s_p_4 - 31090392 * t * s_p_4 +
          11564532 * s_p_5 - 359003232 * t_p_2 + 924349995 * s * t_p_2 -
          457214544 * s_p_2 * t_p_2 - 56275560 * s_p_3 * t_p_2 +
          663029982 * t_p_3 - 1297023192 * s * t_p_3 +
          437192640 * s_p_2 * t_p_3 - 548957736 * t_p_4 +
          592038396 * s * t_p_4 + 167253336 * t_p_5)) /
        4000.;

      dphi(1, 1) =
        (3 * t *
         (598360 * t - 6830960 * s_p_2 + 126762915 * t * s_p_2 +
          14778626 * s_p_3 - 247491396 * t * s_p_3 - 8390844 * s_p_4 +
          92584809 * t * s_p_4 + 3672162 * s_p_5 + 42791328 * t * s_p_5 -
          6533136 * s_p_6 - 15545196 * t * s_p_6 + 3304152 * s_p_7 -
          10840976 * t_p_2 - 718006464 * s_p_2 * t_p_2 +
          1232466660 * s_p_3 * t_p_2 - 457214544 * s_p_4 * t_p_2 -
          45020448 * s_p_5 * t_p_2 + 66905500 * t_p_3 +
          1657574955 * s_p_2 * t_p_3 - 2161705320 * s_p_3 * t_p_3 +
          546490800 * s_p_4 * t_p_3 - 189738504 * t_p_4 -
          1646873208 * s_p_2 * t_p_4 + 1184076792 * s_p_3 * t_p_4 +
          271189044 * t_p_5 + 585386676 * s_p_2 * t_p_5 - 190169856 * t_p_6 +
          52056432 * t_p_7)) /
        4000.;

      dphi(2, 0) =
        (-3 * s * (-1 + s + t) *
         (-1986086 * s - 1986086 * s * t + 19209962 * s_p_2 +
          17223876 * t * s_p_2 - 67816998 * s_p_3 - 50593122 * t * s_p_3 +
          111677490 * s_p_4 + 61084368 * t * s_p_4 - 87112584 * s_p_5 -
          26028216 * t * s_p_5 + 26028216 * s_p_6 + 11244812 * t_p_2 -
          71777907 * s * t_p_2 + 82465353 * s_p_2 * t_p_2 +
          77176476 * s_p_3 * t_p_2 - 109230768 * s_p_4 * t_p_2 -
          42779610 * t_p_3 + 321155847 * s * t_p_3 - 505693692 * s_p_2 * t_p_3 +
          181681704 * s_p_3 * t_p_3 + 25730082 * t_p_4 - 335100240 * s * t_p_4 +
          352372356 * s_p_2 * t_p_4 + 43851780 * t_p_5 + 74870892 * s * t_p_5 -
          38645424 * t_p_6)) /
        2000.;

      dphi(2, 1) =
        (3 * t * (-1 + s + t) *
         (1986086 * t + 1986086 * s * t - 11244812 * s_p_2 +
          71777907 * t * s_p_2 + 42779610 * s_p_3 - 321155847 * t * s_p_3 -
          25730082 * s_p_4 + 335100240 * t * s_p_4 - 43851780 * s_p_5 -
          74870892 * t * s_p_5 + 38645424 * s_p_6 - 19209962 * t_p_2 -
          17223876 * s * t_p_2 - 82465353 * s_p_2 * t_p_2 +
          505693692 * s_p_3 * t_p_2 - 352372356 * s_p_4 * t_p_2 +
          67816998 * t_p_3 + 50593122 * s * t_p_3 - 77176476 * s_p_2 * t_p_3 -
          181681704 * s_p_3 * t_p_3 - 111677490 * t_p_4 - 61084368 * s * t_p_4 +
          109230768 * s_p_2 * t_p_4 + 87112584 * t_p_5 + 26028216 * s * t_p_5 -
          26028216 * t_p_6)) /
        2000.;


      dphi(3, 0) =
        -(s * (18840 * s - 342544 * s_p_2 + 2125500 * s_p_3 - 6071976 * s_p_4 +
               8755236 * s_p_5 - 6200064 * s_p_6 + 1714608 * s_p_7 -
               207476 * t_p_2 + 3854289 * s * t_p_2 - 21879216 * s_p_2 * t_p_2 +
               50684895 * s_p_3 * t_p_2 - 50608152 * s_p_4 * t_p_2 +
               18105444 * s_p_5 * t_p_2 + 445130 * t_p_3 - 7445574 * s * t_p_3 +
               37111140 * s_p_2 * t_p_3 - 65179080 * s_p_3 * t_p_3 +
               35802648 * s_p_4 * t_p_3 - 254376 * t_p_4 + 2759481 * s * t_p_4 -
               13661136 * s_p_2 * t_p_4 + 16297200 * s_p_3 * t_p_4 +
               118458 * t_p_5 + 1288872 * s * t_p_5 - 1324512 * s_p_2 * t_p_5 -
               204768 * t_p_6 - 475308 * s * t_p_6 + 103032 * t_p_7)) /
        400.;

      dphi(3, 1) =
        -(t * s_p_2 *
          (-207476 + 2569526 * s + 667695 * t - 7445574 * s * t -
           10939608 * s_p_2 + 27833355 * t * s_p_2 + 20273958 * s_p_3 -
           39107448 * t * s_p_3 - 16869384 * s_p_4 + 17901324 * t * s_p_4 +
           5172984 * s_p_5 - 508752 * t_p_2 + 3679308 * s * t_p_2 -
           13661136 * s_p_2 * t_p_2 + 13037760 * s_p_3 * t_p_2 +
           296145 * t_p_3 + 2148120 * s * t_p_3 - 1655640 * s_p_2 * t_p_3 -
           614304 * t_p_4 - 950616 * s * t_p_4 + 360612 * t_p_5)) /
        400.;

      dphi(4, 0) =
        -(s * t *
          (13360 - 293388 * s - 18976 * t + 345093 * s * t + 2147744 * s_p_2 -
           1675800 * t * s_p_2 - 7229700 * s_p_3 + 3438135 * t * s_p_3 +
           12258000 * s_p_4 - 3608064 * t * s_p_4 - 10151568 * s_p_5 +
           1589868 * t * s_p_5 + 3255552 * s_p_6 + 110646 * t_p_2 -
           1915110 * s * t_p_2 + 6560028 * s_p_2 * t_p_2 -
           5972940 * s_p_3 * t_p_2 + 1207224 * s_p_4 * t_p_2 - 244620 * t_p_3 +
           4806621 * s * t_p_3 - 14258592 * s_p_2 * t_p_3 +
           8592480 * s_p_3 * t_p_3 + 37854 * t_p_4 - 3399084 * s * t_p_4 +
           6951744 * s_p_2 * t_p_4 + 204768 * t_p_5 + 475308 * s * t_p_5 -
           103032 * t_p_6)) /
        400.;

      dphi(4, 1) =
        -(s_p_2 *
          (6680 - 97796 * s - 18976 * t + 230062 * s * t + 536936 * s_p_2 -
           837900 * t * s_p_2 - 1445940 * s_p_3 + 1375254 * t * s_p_3 +
           2043000 * s_p_4 - 1202688 * t * s_p_4 - 1450224 * s_p_5 +
           454248 * t * s_p_5 + 406944 * s_p_6 + 165969 * t_p_2 -
           1915110 * s * t_p_2 + 4920021 * s_p_2 * t_p_2 -
           3583764 * s_p_3 * t_p_2 + 603612 * s_p_4 * t_p_2 - 489240 * t_p_3 +
           6408828 * s * t_p_3 - 14258592 * s_p_2 * t_p_3 +
           6873984 * s_p_3 * t_p_3 + 94635 * t_p_4 - 5665140 * s * t_p_4 +
           8689680 * s_p_2 * t_p_4 + 614304 * t_p_5 + 950616 * s * t_p_5 -
           360612 * t_p_6)) /
        400.;

      dphi(5, 0) =
        -(t_p_2 *
          (6680 - 18976 * s - 97796 * t + 230062 * s * t + 165969 * s_p_2 -
           1915110 * t * s_p_2 - 489240 * s_p_3 + 6408828 * t * s_p_3 +
           94635 * s_p_4 - 5665140 * t * s_p_4 + 614304 * s_p_5 +
           950616 * t * s_p_5 - 360612 * s_p_6 + 536936 * t_p_2 -
           837900 * s * t_p_2 + 4920021 * s_p_2 * t_p_2 -
           14258592 * s_p_3 * t_p_2 + 8689680 * s_p_4 * t_p_2 -
           1445940 * t_p_3 + 1375254 * s * t_p_3 - 3583764 * s_p_2 * t_p_3 +
           6873984 * s_p_3 * t_p_3 + 2043000 * t_p_4 - 1202688 * s * t_p_4 +
           603612 * s_p_2 * t_p_4 - 1450224 * t_p_5 + 454248 * s * t_p_5 +
           406944 * t_p_6)) /
        400.;

      dphi(5, 1) =
        (s * t *
         (-13360 + 18976 * s + 293388 * t - 345093 * s * t - 110646 * s_p_2 +
          1915110 * t * s_p_2 + 244620 * s_p_3 - 4806621 * t * s_p_3 -
          37854 * s_p_4 + 3399084 * t * s_p_4 - 204768 * s_p_5 -
          475308 * t * s_p_5 + 103032 * s_p_6 - 2147744 * t_p_2 +
          1675800 * s * t_p_2 - 6560028 * s_p_2 * t_p_2 +
          14258592 * s_p_3 * t_p_2 - 6951744 * s_p_4 * t_p_2 + 7229700 * t_p_3 -
          3438135 * s * t_p_3 + 5972940 * s_p_2 * t_p_3 -
          8592480 * s_p_3 * t_p_3 - 12258000 * t_p_4 + 3608064 * s * t_p_4 -
          1207224 * s_p_2 * t_p_4 + 10151568 * t_p_5 - 1589868 * s * t_p_5 -
          3255552 * t_p_6)) /
        400.;

      dphi(6, 0) =
        -(s * t_p_2 *
          (-207476 + 667695 * s + 2569526 * t - 7445574 * s * t -
           508752 * s_p_2 + 3679308 * t * s_p_2 + 296145 * s_p_3 +
           2148120 * t * s_p_3 - 614304 * s_p_4 - 950616 * t * s_p_4 +
           360612 * s_p_5 - 10939608 * t_p_2 + 27833355 * s * t_p_2 -
           13661136 * s_p_2 * t_p_2 - 1655640 * s_p_3 * t_p_2 +
           20273958 * t_p_3 - 39107448 * s * t_p_3 + 13037760 * s_p_2 * t_p_3 -
           16869384 * t_p_4 + 17901324 * s * t_p_4 + 5172984 * t_p_5)) /
        400.;

      dphi(6, 1) =
        -(t *
          (18840 * t - 207476 * s_p_2 + 3854289 * t * s_p_2 + 445130 * s_p_3 -
           7445574 * t * s_p_3 - 254376 * s_p_4 + 2759481 * t * s_p_4 +
           118458 * s_p_5 + 1288872 * t * s_p_5 - 204768 * s_p_6 -
           475308 * t * s_p_6 + 103032 * s_p_7 - 342544 * t_p_2 -
           21879216 * s_p_2 * t_p_2 + 37111140 * s_p_3 * t_p_2 -
           13661136 * s_p_4 * t_p_2 - 1324512 * s_p_5 * t_p_2 +
           2125500 * t_p_3 + 50684895 * s_p_2 * t_p_3 -
           65179080 * s_p_3 * t_p_3 + 16297200 * s_p_4 * t_p_3 -
           6071976 * t_p_4 - 50608152 * s_p_2 * t_p_4 +
           35802648 * s_p_3 * t_p_4 + 8755236 * t_p_5 +
           18105444 * s_p_2 * t_p_5 - 6200064 * t_p_6 + 1714608 * t_p_7)) /
        400.;

      dphi(7, 0) =
        -((-1 + s + t) *
          (200 + 200 * s + 200 * t + 400 * s * t - 80134 * s_p_2 -
           79734 * t * s_p_2 + 699978 * s_p_3 + 620244 * t * s_p_3 -
           2357262 * s_p_4 - 1737018 * t * s_p_4 + 3780810 * s_p_5 +
           2043792 * t * s_p_5 - 2901096 * s_p_6 - 857304 * t * s_p_6 +
           857304 * s_p_7 - 29918 * t_p_2 + 523968 * s * t_p_2 -
           1789389 * s_p_2 * t_p_2 - 267453 * s_p_3 * t_p_2 +
           6328044 * s_p_4 * t_p_2 - 5089392 * s_p_5 * t_p_2 + 210868 * t_p_3 -
           3191508 * s * t_p_3 + 11235411 * s_p_2 * t_p_3 -
           11357658 * s_p_3 * t_p_3 + 1514376 * s_p_4 * t_p_3 - 629370 * t_p_4 +
           7117704 * s * t_p_4 - 17931294 * s_p_2 * t_p_4 +
           11312784 * s_p_3 * t_p_4 + 943740 * t_p_5 - 6930360 * s * t_p_5 +
           8812800 * s_p_2 * t_p_5 - 699192 * t_p_6 + 2489616 * s * t_p_6 +
           203472 * t_p_7)) /
        200.;

      dphi(7, 1) =
        (s * t * (-1 + s + t) *
         (30118 - 246625 * s - 331061 * t + 2367072 * s * t + 497916 * s_p_2 -
          5243166 * t * s_p_2 + 272493 * s_p_3 + 3401433 * t * s_p_3 -
          1394010 * s_p_4 + 44226 * t * s_p_4 + 849528 * s_p_5 +
          1349415 * t_p_2 - 7222095 * s * t_p_2 + 11724345 * s_p_2 * t_p_2 -
          5175090 * s_p_3 * t_p_2 - 2583360 * t_p_3 + 8934300 * s * t_p_3 -
          7403400 * s_p_2 * t_p_3 + 2345436 * t_p_4 - 3899016 * s * t_p_4 -
          813888 * t_p_5)) /
        100.;

      dphi(8, 0) =
        -(s * t * (-1 + s + t) *
          (-30118 + 331061 * s + 246625 * t - 2367072 * s * t -
           1349415 * s_p_2 + 7222095 * t * s_p_2 + 2583360 * s_p_3 -
           8934300 * t * s_p_3 - 2345436 * s_p_4 + 3899016 * t * s_p_4 +
           813888 * s_p_5 - 497916 * t_p_2 + 5243166 * s * t_p_2 -
           11724345 * s_p_2 * t_p_2 + 7403400 * s_p_3 * t_p_2 - 272493 * t_p_3 -
           3401433 * s * t_p_3 + 5175090 * s_p_2 * t_p_3 + 1394010 * t_p_4 -
           44226 * s * t_p_4 - 849528 * t_p_5)) /
        100.;

      dphi(8, 1) =
        -((-1 + s + t) *
          (200 + 200 * s + 200 * t + 400 * s * t - 29918 * s_p_2 +
           523968 * t * s_p_2 + 210868 * s_p_3 - 3191508 * t * s_p_3 -
           629370 * s_p_4 + 7117704 * t * s_p_4 + 943740 * s_p_5 -
           6930360 * t * s_p_5 - 699192 * s_p_6 + 2489616 * t * s_p_6 +
           203472 * s_p_7 - 80134 * t_p_2 - 79734 * s * t_p_2 -
           1789389 * s_p_2 * t_p_2 + 11235411 * s_p_3 * t_p_2 -
           17931294 * s_p_4 * t_p_2 + 8812800 * s_p_5 * t_p_2 + 699978 * t_p_3 +
           620244 * s * t_p_3 - 267453 * s_p_2 * t_p_3 -
           11357658 * s_p_3 * t_p_3 + 11312784 * s_p_4 * t_p_3 -
           2357262 * t_p_4 - 1737018 * s * t_p_4 + 6328044 * s_p_2 * t_p_4 +
           1514376 * s_p_3 * t_p_4 + 3780810 * t_p_5 + 2043792 * s * t_p_5 -
           5089392 * s_p_2 * t_p_5 - 2901096 * t_p_6 - 857304 * s * t_p_6 +
           857304 * t_p_7)) /
        200.;

      dphi(9, 0) =
        (s *
         (120 * s - 2192 * s_p_2 + 13700 * s_p_3 - 39528 * s_p_4 +
          57708 * s_p_5 - 41472 * s_p_6 + 11664 * s_p_7 - 1256 * t_p_2 +
          23349 * s * t_p_2 - 132768 * s_p_2 * t_p_2 + 308385 * s_p_3 * t_p_2 -
          309096 * s_p_4 * t_p_2 + 111132 * s_p_5 * t_p_2 + 2678 * t_p_3 -
          44712 * s * t_p_3 + 223020 * s_p_2 * t_p_3 - 392040 * s_p_3 * t_p_3 +
          215784 * s_p_4 * t_p_3 - 1548 * t_p_4 + 16443 * s * t_p_4 -
          81648 * s_p_2 * t_p_4 + 97200 * s_p_3 * t_p_4 + 774 * t_p_5 +
          7776 * s * t_p_5 - 7776 * s_p_2 * t_p_5 - 1296 * t_p_6 -
          2916 * s * t_p_6 + 648 * t_p_7)) /
        80.;


      dphi(9, 1) =
        (t * s_p_2 *
         (-1256 + 15566 * s + 4017 * t - 44712 * s * t - 66384 * s_p_2 +
          167265 * t * s_p_2 + 123354 * s_p_3 - 235224 * t * s_p_3 -
          103032 * s_p_4 + 107892 * t * s_p_4 + 31752 * s_p_5 - 3096 * t_p_2 +
          21924 * s * t_p_2 - 81648 * s_p_2 * t_p_2 + 77760 * s_p_3 * t_p_2 +
          1935 * t_p_3 + 12960 * s * t_p_3 - 9720 * s_p_2 * t_p_3 -
          3888 * t_p_4 - 5832 * s * t_p_4 + 2268 * t_p_5)) /
        80.;

      dphi(10, 0) =
        (s * t *
         (80 - 1764 * s - 100 * t + 1797 * s * t + 12992 * s_p_2 -
          8280 * t * s_p_2 - 44100 * s_p_3 + 15705 * t * s_p_3 + 75600 * s_p_4 -
          15552 * t * s_p_4 - 63504 * s_p_5 + 6804 * t * s_p_5 + 20736 * s_p_6 +
          650 * t_p_2 - 11340 * s * t_p_2 + 38484 * s_p_2 * t_p_2 -
          34020 * s_p_3 * t_p_2 + 5832 * s_p_4 * t_p_2 - 1440 * t_p_3 +
          28863 * s * t_p_3 - 85536 * s_p_2 * t_p_3 + 51840 * s_p_3 * t_p_3 +
          162 * t_p_4 - 20412 * s * t_p_4 + 41472 * s_p_2 * t_p_4 +
          1296 * t_p_5 + 2916 * s * t_p_5 - 648 * t_p_6)) /
        40.;

      dphi(10, 1) =
        (s_p_2 *
         (40 - 588 * s - 100 * t + 1198 * s * t + 3248 * s_p_2 -
          4140 * t * s_p_2 - 8820 * s_p_3 + 6282 * t * s_p_3 + 12600 * s_p_4 -
          5184 * t * s_p_4 - 9072 * s_p_5 + 1944 * t * s_p_5 + 2592 * s_p_6 +
          975 * t_p_2 - 11340 * s * t_p_2 + 28863 * s_p_2 * t_p_2 -
          20412 * s_p_3 * t_p_2 + 2916 * s_p_4 * t_p_2 - 2880 * t_p_3 +
          38484 * s * t_p_3 - 85536 * s_p_2 * t_p_3 + 41472 * s_p_3 * t_p_3 +
          405 * t_p_4 - 34020 * s * t_p_4 + 51840 * s_p_2 * t_p_4 +
          3888 * t_p_5 + 5832 * s * t_p_5 - 2268 * t_p_6)) /
        40.;

      dphi(11, 0) =
        (s * t_p_2 *
         (20 - 153 * s - 650 * t + 11340 * s * t - 2520 * s_p_2 -
          38484 * t * s_p_2 + 14895 * s_p_3 + 34020 * t * s_p_3 -
          23328 * s_p_4 - 5832 * t * s_p_4 + 11340 * s_p_5 + 1440 * t_p_2 -
          28863 * s * t_p_2 + 85536 * s_p_2 * t_p_2 - 51840 * s_p_3 * t_p_2 -
          162 * t_p_3 + 20412 * s * t_p_3 - 41472 * s_p_2 * t_p_3 -
          1296 * t_p_4 - 2916 * s * t_p_4 + 648 * t_p_5)) /
        80.;


      dphi(11, 1) =
        (t * s_p_2 *
         (20 - 102 * s - 975 * t + 11340 * s * t - 1260 * s_p_2 -
          28863 * t * s_p_2 + 5958 * s_p_3 + 20412 * t * s_p_3 - 7776 * s_p_4 -
          2916 * t * s_p_4 + 3240 * s_p_5 + 2880 * t_p_2 - 38484 * s * t_p_2 +
          85536 * s_p_2 * t_p_2 - 41472 * s_p_3 * t_p_2 - 405 * t_p_3 +
          34020 * s * t_p_3 - 51840 * s_p_2 * t_p_3 - 3888 * t_p_4 -
          5832 * s * t_p_4 + 2268 * t_p_5)) /
        80.;

      dphi(12, 0) =
        (s * t_p_2 *
         (20 - 975 * s - 102 * t + 11340 * s * t + 2880 * s_p_2 -
          38484 * t * s_p_2 - 405 * s_p_3 + 34020 * t * s_p_3 - 3888 * s_p_4 -
          5832 * t * s_p_4 + 2268 * s_p_5 - 1260 * t_p_2 - 28863 * s * t_p_2 +
          85536 * s_p_2 * t_p_2 - 51840 * s_p_3 * t_p_2 + 5958 * t_p_3 +
          20412 * s * t_p_3 - 41472 * s_p_2 * t_p_3 - 7776 * t_p_4 -
          2916 * s * t_p_4 + 3240 * t_p_5)) /
        80.;

      dphi(12, 1) =
        (t * s_p_2 *
         (20 - 650 * s - 153 * t + 11340 * s * t + 1440 * s_p_2 -
          28863 * t * s_p_2 - 162 * s_p_3 + 20412 * t * s_p_3 - 1296 * s_p_4 -
          2916 * t * s_p_4 + 648 * s_p_5 - 2520 * t_p_2 - 38484 * s * t_p_2 +
          85536 * s_p_2 * t_p_2 - 41472 * s_p_3 * t_p_2 + 14895 * t_p_3 +
          34020 * s * t_p_3 - 51840 * s_p_2 * t_p_3 - 23328 * t_p_4 -
          5832 * s * t_p_4 + 11340 * t_p_5)) /
        80.;

      dphi(13, 0) =
        (t_p_2 *
         (40 - 100 * s - 588 * t + 1198 * s * t + 975 * s_p_2 -
          11340 * t * s_p_2 - 2880 * s_p_3 + 38484 * t * s_p_3 + 405 * s_p_4 -
          34020 * t * s_p_4 + 3888 * s_p_5 + 5832 * t * s_p_5 - 2268 * s_p_6 +
          3248 * t_p_2 - 4140 * s * t_p_2 + 28863 * s_p_2 * t_p_2 -
          85536 * s_p_3 * t_p_2 + 51840 * s_p_4 * t_p_2 - 8820 * t_p_3 +
          6282 * s * t_p_3 - 20412 * s_p_2 * t_p_3 + 41472 * s_p_3 * t_p_3 +
          12600 * t_p_4 - 5184 * s * t_p_4 + 2916 * s_p_2 * t_p_4 -
          9072 * t_p_5 + 1944 * s * t_p_5 + 2592 * t_p_6)) /
        40.;

      dphi(13, 1) =
        -(s * t *
          (-80 + 100 * s + 1764 * t - 1797 * s * t - 650 * s_p_2 +
           11340 * t * s_p_2 + 1440 * s_p_3 - 28863 * t * s_p_3 - 162 * s_p_4 +
           20412 * t * s_p_4 - 1296 * s_p_5 - 2916 * t * s_p_5 + 648 * s_p_6 -
           12992 * t_p_2 + 8280 * s * t_p_2 - 38484 * s_p_2 * t_p_2 +
           85536 * s_p_3 * t_p_2 - 41472 * s_p_4 * t_p_2 + 44100 * t_p_3 -
           15705 * s * t_p_3 + 34020 * s_p_2 * t_p_3 - 51840 * s_p_3 * t_p_3 -
           75600 * t_p_4 + 15552 * s * t_p_4 - 5832 * s_p_2 * t_p_4 +
           63504 * t_p_5 - 6804 * s * t_p_5 - 20736 * t_p_6)) /
        40.;

      dphi(14, 0) =
        (s * t_p_2 *
         (-1256 + 4017 * s + 15566 * t - 44712 * s * t - 3096 * s_p_2 +
          21924 * t * s_p_2 + 1935 * s_p_3 + 12960 * t * s_p_3 - 3888 * s_p_4 -
          5832 * t * s_p_4 + 2268 * s_p_5 - 66384 * t_p_2 + 167265 * s * t_p_2 -
          81648 * s_p_2 * t_p_2 - 9720 * s_p_3 * t_p_2 + 123354 * t_p_3 -
          235224 * s * t_p_3 + 77760 * s_p_2 * t_p_3 - 103032 * t_p_4 +
          107892 * s * t_p_4 + 31752 * t_p_5)) /
        80.;

      dphi(14, 1) =
        (t * (120 * t - 1256 * s_p_2 + 23349 * t * s_p_2 + 2678 * s_p_3 -
              44712 * t * s_p_3 - 1548 * s_p_4 + 16443 * t * s_p_4 +
              774 * s_p_5 + 7776 * t * s_p_5 - 1296 * s_p_6 - 2916 * t * s_p_6 +
              648 * s_p_7 - 2192 * t_p_2 - 132768 * s_p_2 * t_p_2 +
              223020 * s_p_3 * t_p_2 - 81648 * s_p_4 * t_p_2 -
              7776 * s_p_5 * t_p_2 + 13700 * t_p_3 + 308385 * s_p_2 * t_p_3 -
              392040 * s_p_3 * t_p_3 + 97200 * s_p_4 * t_p_3 - 39528 * t_p_4 -
              309096 * s_p_2 * t_p_4 + 215784 * s_p_3 * t_p_4 + 57708 * t_p_5 +
              111132 * s_p_2 * t_p_5 - 41472 * t_p_6 + 11664 * t_p_7)) /
        80.;

      dphi(15, 0) =
        -(s * (-1 + s + t) *
          (40 - 842 * s + 40 * t - 802 * s * t + 5734 * s_p_2 +
           4932 * t * s_p_2 - 17586 * s_p_3 - 12654 * t * s_p_3 +
           26910 * s_p_4 + 14256 * t * s_p_4 - 20088 * s_p_5 -
           5832 * t * s_p_5 + 5832 * s_p_6 - 638 * t_p_2 + 11133 * s * t_p_2 -
           54459 * s_p_2 * t_p_2 + 94932 * s_p_3 * t_p_2 -
           53136 * s_p_4 * t_p_2 + 1026 * t_p_3 - 15867 * s * t_p_3 +
           60426 * s_p_2 * t_p_3 - 57672 * s_p_3 * t_p_3 - 468 * t_p_4 +
           6318 * s * t_p_4 - 16848 * s_p_2 * t_p_4)) /
        40.;

      dphi(15, 1) =
        (3 * (-1 + 3 * s) * (-1 + 6 * s) * t * (-1 + s + t) * s_p_2 *
         (113 - 267 * s - 303 * t + 357 * s * t + 156 * s_p_2 + 195 * t_p_2)) /
        20.;

      dphi(16, 0) =
        -(t * (-1 + s + t) *
          (10 - 304 * s - 147 * t + 3385 * s * t + 2573 * s_p_2 -
           21195 * t * s_p_2 - 9495 * s_p_3 + 54990 * t * s_p_3 +
           17280 * s_p_4 - 63180 * t * s_p_4 - 15228 * s_p_5 +
           26568 * t * s_p_5 + 5184 * s_p_6 + 812 * t_p_2 - 13905 * s * t_p_2 +
           61290 * s_p_2 * t_p_2 - 100440 * s_p_3 * t_p_2 +
           55080 * s_p_4 * t_p_2 - 2205 * t_p_3 + 26730 * s * t_p_3 -
           74520 * s_p_2 * t_p_3 + 58320 * s_p_3 * t_p_3 + 3150 * t_p_4 -
           24300 * s * t_p_4 + 32400 * s_p_2 * t_p_4 - 2268 * t_p_5 +
           8424 * s * t_p_5 + 648 * t_p_6)) /
        10.;


      dphi(16, 1) =
        -(s * (-1 + s + t) *
          (10 - 147 * s - 304 * t + 3385 * s * t + 812 * s_p_2 -
           13905 * t * s_p_2 - 2205 * s_p_3 + 26730 * t * s_p_3 + 3150 * s_p_4 -
           24300 * t * s_p_4 - 2268 * s_p_5 + 8424 * t * s_p_5 + 648 * s_p_6 +
           2573 * t_p_2 - 21195 * s * t_p_2 + 61290 * s_p_2 * t_p_2 -
           74520 * s_p_3 * t_p_2 + 32400 * s_p_4 * t_p_2 - 9495 * t_p_3 +
           54990 * s * t_p_3 - 100440 * s_p_2 * t_p_3 + 58320 * s_p_3 * t_p_3 +
           17280 * t_p_4 - 63180 * s * t_p_4 + 55080 * s_p_2 * t_p_4 -
           15228 * t_p_5 + 26568 * s * t_p_5 + 5184 * t_p_6)) /
        10.;

      dphi(17, 0) =
        (3 * s * (-1 + s + t) * (-1 + 3 * t) * (-1 + 6 * t) * t_p_2 *
         (113 - 303 * s - 267 * t + 357 * s * t + 195 * s_p_2 + 156 * t_p_2)) /
        20.;

      dphi(17, 1) =
        -(t * (-1 + s + t) *
          (40 + 40 * s - 842 * t - 802 * s * t - 638 * s_p_2 +
           11133 * t * s_p_2 + 1026 * s_p_3 - 15867 * t * s_p_3 - 468 * s_p_4 +
           6318 * t * s_p_4 + 5734 * t_p_2 + 4932 * s * t_p_2 -
           54459 * s_p_2 * t_p_2 + 60426 * s_p_3 * t_p_2 -
           16848 * s_p_4 * t_p_2 - 17586 * t_p_3 - 12654 * s * t_p_3 +
           94932 * s_p_2 * t_p_3 - 57672 * s_p_3 * t_p_3 + 26910 * t_p_4 +
           14256 * s * t_p_4 - 53136 * s_p_2 * t_p_4 - 20088 * t_p_5 -
           5832 * s * t_p_5 + 5832 * t_p_6)) /
        40.;


      dphi(18, 0) =
        -16 * (-1 + s + t) * (-1 + 3 * t) * (-1 + 6 * t) * t_p_2 *
        (-10 + 84 * s + 37 * t - 207 * s * t - 162 * s_p_2 + 198 * t * s_p_2 +
         90 * s_p_3 - 45 * t_p_2 + 126 * s * t_p_2 + 18 * t_p_3);

      dphi(18, 1) =
        -16 * s * t * (-1 + s + t) * (-1 + s + 2 * t) *
        (20 - 54 * s - 351 * t + 837 * s * t + 36 * s_p_2 - 486 * t * s_p_2 +
         1647 * t_p_2 - 3078 * s * t_p_2 + 1296 * s_p_2 * t_p_2 - 2592 * t_p_3 +
         2592 * s * t_p_3 + 1296 * t_p_4);

      dphi(19, 0) =
        -16 * s * t * (-1 + s + t) * (-1 + 2 * s + t) *
        (20 - 351 * s - 54 * t + 837 * s * t + 1647 * s_p_2 - 3078 * t * s_p_2 -
         2592 * s_p_3 + 2592 * t * s_p_3 + 1296 * s_p_4 + 36 * t_p_2 -
         486 * s * t_p_2 + 1296 * s_p_2 * t_p_2);

      dphi(19, 1) =
        -16 * (-1 + 3 * s) * (-1 + 6 * s) * (-1 + s + t) * s_p_2 *
        (-10 + 37 * s + 84 * t - 207 * s * t - 45 * s_p_2 + 126 * t * s_p_2 +
         18 * s_p_3 - 162 * t_p_2 + 198 * s * t_p_2 + 90 * t_p_3);

      dphi(20, 0) = 8 * s * (-1 + 3 * t) * (-1 + 6 * t) *
                    (-2 + 30 * s + 2 * t - 27 * s * t - 108 * s_p_2 +
                     72 * t * s_p_2 + 90 * s_p_3) *
                    t_p_2 * sqrt(2);

      dphi(20, 1) = 8 * (-1 + 3 * s) * (-1 + 6 * s) * t * s_p_2 *
                    (-2 + 2 * s + 30 * t - 27 * s * t - 108 * t_p_2 +
                     72 * s * t_p_2 + 90 * t_p_3) *
                    sqrt(2);

      dphi(21, 0) =
        (-11664 * s * (-1 + s + t) * (-1 + 3 * t) * (-1 + 6 * t) * t_p_2 *
         (93 - 303 * s - 197 * t + 357 * s * t + 195 * s_p_2 + 96 * t_p_2)) /
        125.;

      dphi(21, 1) =
        (-5832 * t * (-1 + s + t) *
         (-16 * t - 16 * s * t + 186 * s_p_2 - 3211 * t * s_p_2 - 342 * s_p_3 +
          5289 * t * s_p_3 + 156 * s_p_4 - 2106 * t * s_p_4 + 272 * t_p_2 +
          256 * s * t_p_2 + 15353 * s_p_2 * t_p_2 - 20142 * s_p_3 * t_p_2 +
          5616 * s_p_4 * t_p_2 - 1488 * t_p_3 - 1232 * s * t_p_3 -
          25644 * s_p_2 * t_p_3 + 19224 * s_p_3 * t_p_3 + 3440 * t_p_4 +
          2208 * s * t_p_4 + 13392 * s_p_2 * t_p_4 - 3504 * t_p_5 -
          1296 * s * t_p_5 + 1296 * t_p_6)) /
        125.;

      dphi(22, 0) =
        (729 * s * (-1 + s + t) * (-1 + 3 * t) * (-1 + 6 * t) * t_p_2 *
         (61 - 239 * s - 107 * t + 229 * s * t + 195 * s_p_2 + 44 * t_p_2)) /
        8.;

      dphi(22, 1) =
        (729 * t * (-1 + s + t) *
         (-10 * t - 10 * s * t + 122 * s_p_2 - 2039 * t * s_p_2 - 278 * s_p_3 +
          4233 * t * s_p_3 + 156 * s_p_4 - 2106 * t * s_p_4 + 166 * t_p_2 +
          156 * s * t_p_2 + 9161 * s_p_2 * t_p_2 - 15534 * s_p_3 * t_p_2 +
          5616 * s_p_4 * t_p_2 - 874 * t_p_3 - 718 * s * t_p_3 -
          13692 * s_p_2 * t_p_3 + 13464 * s_p_3 * t_p_3 + 1918 * t_p_4 +
          1200 * s * t_p_4 + 6192 * s_p_2 * t_p_4 - 1848 * t_p_5 -
          648 * s * t_p_5 + 648 * t_p_6)) /
        16.;

      dphi(23, 0) = (-81 * s * (-1 + s + t) * (-1 + 6 * t) * t_p_2 *
                     (-1003 + 6353 * s + 2632 * t - 12582 * s * t -
                      11277 * s_p_2 + 11727 * t * s_p_2 + 6048 * s_p_3 -
                      1971 * t_p_2 + 5913 * s * t_p_2 + 324 * t_p_3)) /
                    8.;

      dphi(23, 1) =
        (-81 * t * (-1 + s + t) *
         (-180 * t - 180 * s * t + 2006 * s_p_2 - 27133 * t * s_p_2 -
          7802 * s_p_3 + 96471 * t * s_p_3 + 9828 * s_p_4 - 105750 * t * s_p_4 -
          4032 * s_p_5 + 36288 * t * s_p_5 + 2628 * t_p_2 + 2448 * s * t_p_2 +
          80967 * s_p_2 * t_p_2 - 222210 * s_p_3 * t_p_2 +
          134352 * s_p_4 * t_p_2 - 11502 * t_p_3 - 9054 * s * t_p_3 -
          77652 * s_p_2 * t_p_3 + 130248 * s_p_3 * t_p_3 + 21582 * t_p_4 +
          12528 * s * t_p_4 + 19440 * s_p_2 * t_p_4 - 18360 * t_p_5 -
          5832 * s * t_p_5 + 5832 * t_p_6)) /
        16.;

      dphi(24, 0) = (324 * s * (-1 + s + t) * t_p_2 *
                     (11858 - 100135 * s - 53894 * t + 399855 * s * t +
                      268917 * s_p_2 - 762444 * t * s_p_2 - 292788 * s_p_3 +
                      427788 * t * s_p_3 + 112644 * s_p_4 + 62844 * t_p_2 -
                      448092 * s * t_p_2 + 474552 * s_p_2 * t_p_2 -
                      648 * t_p_3 + 135432 * s * t_p_3 - 20736 * t_p_4)) /
                    125.;

      dphi(24, 1) =
        (648 * t * (-1 + s + t) *
         (-720 * t - 720 * s * t + 5929 * s_p_2 - 44105 * t * s_p_2 -
          31402 * s_p_3 + 201435 * t * s_p_3 + 60861 * s_p_4 -
          274410 * t * s_p_4 - 51480 * s_p_5 + 119016 * t * s_p_5 +
          16092 * s_p_6 + 7632 * t_p_2 + 6912 * s * t_p_2 +
          79545 * s_p_2 * t_p_2 - 326214 * s_p_3 * t_p_2 +
          241920 * s_p_4 * t_p_2 - 28368 * t_p_3 - 21456 * s * t_p_3 -
          21276 * s_p_2 * t_p_3 + 139320 * s_p_3 * t_p_3 + 48240 * t_p_4 +
          26784 * s * t_p_4 - 24624 * s_p_2 * t_p_4 - 38448 * t_p_5 -
          11664 * s * t_p_5 + 11664 * t_p_6)) /
        125.;


      dphi(25, 0) =
        (648 * s * (-1 + s + t) *
         (-720 * s - 720 * s * t + 7632 * s_p_2 + 6912 * t * s_p_2 -
          28368 * s_p_3 - 21456 * t * s_p_3 + 48240 * s_p_4 +
          26784 * t * s_p_4 - 38448 * s_p_5 - 11664 * t * s_p_5 +
          11664 * s_p_6 + 5929 * t_p_2 - 44105 * s * t_p_2 +
          79545 * s_p_2 * t_p_2 - 21276 * s_p_3 * t_p_2 -
          24624 * s_p_4 * t_p_2 - 31402 * t_p_3 + 201435 * s * t_p_3 -
          326214 * s_p_2 * t_p_3 + 139320 * s_p_3 * t_p_3 + 60861 * t_p_4 -
          274410 * s * t_p_4 + 241920 * s_p_2 * t_p_4 - 51480 * t_p_5 +
          119016 * s * t_p_5 + 16092 * t_p_6)) /
        125.;

      dphi(25, 1) = (-324 * t * (-1 + s + t) * s_p_2 *
                     (-11858 + 53894 * s + 100135 * t - 399855 * s * t -
                      62844 * s_p_2 + 448092 * t * s_p_2 + 648 * s_p_3 -
                      135432 * t * s_p_3 + 20736 * s_p_4 - 268917 * t_p_2 +
                      762444 * s * t_p_2 - 474552 * s_p_2 * t_p_2 +
                      292788 * t_p_3 - 427788 * s * t_p_3 - 112644 * t_p_4)) /
                    125.;

      dphi(26, 0) =
        (-81 * s * (-1 + s + t) *
         (-180 * s - 180 * s * t + 2628 * s_p_2 + 2448 * t * s_p_2 -
          11502 * s_p_3 - 9054 * t * s_p_3 + 21582 * s_p_4 + 12528 * t * s_p_4 -
          18360 * s_p_5 - 5832 * t * s_p_5 + 5832 * s_p_6 + 2006 * t_p_2 -
          27133 * s * t_p_2 + 80967 * s_p_2 * t_p_2 - 77652 * s_p_3 * t_p_2 +
          19440 * s_p_4 * t_p_2 - 7802 * t_p_3 + 96471 * s * t_p_3 -
          222210 * s_p_2 * t_p_3 + 130248 * s_p_3 * t_p_3 + 9828 * t_p_4 -
          105750 * s * t_p_4 + 134352 * s_p_2 * t_p_4 - 4032 * t_p_5 +
          36288 * s * t_p_5)) /
        16.;


      dphi(26, 1) = (-81 * (-1 + 6 * s) * t * (-1 + s + t) * s_p_2 *
                     (-1003 + 2632 * s + 6353 * t - 12582 * s * t -
                      1971 * s_p_2 + 5913 * t * s_p_2 + 324 * s_p_3 -
                      11277 * t_p_2 + 11727 * s * t_p_2 + 6048 * t_p_3)) /
                    8.;

      dphi(27, 0) =
        (729 * s * (-1 + s + t) *
         (-10 * s - 10 * s * t + 166 * s_p_2 + 156 * t * s_p_2 - 874 * s_p_3 -
          718 * t * s_p_3 + 1918 * s_p_4 + 1200 * t * s_p_4 - 1848 * s_p_5 -
          648 * t * s_p_5 + 648 * s_p_6 + 122 * t_p_2 - 2039 * s * t_p_2 +
          9161 * s_p_2 * t_p_2 - 13692 * s_p_3 * t_p_2 + 6192 * s_p_4 * t_p_2 -
          278 * t_p_3 + 4233 * s * t_p_3 - 15534 * s_p_2 * t_p_3 +
          13464 * s_p_3 * t_p_3 + 156 * t_p_4 - 2106 * s * t_p_4 +
          5616 * s_p_2 * t_p_4)) /
        16.;

      dphi(27, 1) =
        (729 * (-1 + 3 * s) * (-1 + 6 * s) * t * (-1 + s + t) * s_p_2 *
         (61 - 107 * s - 239 * t + 229 * s * t + 44 * s_p_2 + 195 * t_p_2)) /
        8.;

      dphi(28, 0) =
        (-5832 * s * (-1 + s + t) *
         (-16 * s - 16 * s * t + 272 * s_p_2 + 256 * t * s_p_2 - 1488 * s_p_3 -
          1232 * t * s_p_3 + 3440 * s_p_4 + 2208 * t * s_p_4 - 3504 * s_p_5 -
          1296 * t * s_p_5 + 1296 * s_p_6 + 186 * t_p_2 - 3211 * s * t_p_2 +
          15353 * s_p_2 * t_p_2 - 25644 * s_p_3 * t_p_2 +
          13392 * s_p_4 * t_p_2 - 342 * t_p_3 + 5289 * s * t_p_3 -
          20142 * s_p_2 * t_p_3 + 19224 * s_p_3 * t_p_3 + 156 * t_p_4 -
          2106 * s * t_p_4 + 5616 * s_p_2 * t_p_4)) /
        125.;

      dphi(28, 1) =
        (-11664 * (-1 + 3 * s) * (-1 + 6 * s) * t * (-1 + s + t) * s_p_2 *
         (93 - 197 * s - 303 * t + 357 * s * t + 96 * s_p_2 + 195 * t_p_2)) /
        125.;


      dphi(29, 0) =
        (-162 * s * t_p_2 *
         (-1036 + 22449 * s - 2626 * t + 42000 * s * t - 162800 * s_p_2 -
          87564 * t * s_p_2 + 485295 * s_p_3 - 100080 * t * s_p_3 -
          604800 * s_p_4 + 178200 * t * s_p_4 + 262332 * s_p_5 + 9080 * t_p_2 -
          170253 * s * t_p_2 + 504576 * s_p_2 * t_p_2 - 294840 * s_p_3 * t_p_2 -
          2322 * t_p_3 + 121392 * s * t_p_3 - 253152 * s_p_2 * t_p_3 -
          6336 * t_p_4 - 15876 * s * t_p_4 + 3240 * t_p_5)) /
        125.;

      dphi(29, 1) =
        (-162 * t * s_p_2 *
         (-1036 + 14966 * s - 3939 * t + 42000 * s * t - 81400 * s_p_2 -
          65673 * t * s_p_2 + 194118 * s_p_3 - 60048 * t * s_p_3 -
          201600 * s_p_4 + 89100 * t * s_p_4 + 74952 * s_p_5 + 18160 * t_p_2 -
          227004 * s * t_p_2 + 504576 * s_p_2 * t_p_2 - 235872 * s_p_3 * t_p_2 -
          5805 * t_p_3 + 202320 * s * t_p_3 - 316440 * s_p_2 * t_p_3 -
          19008 * t_p_4 - 31752 * s * t_p_4 + 11340 * t_p_5)) /
        125.;

      dphi(30, 0) =
        (-81 * s * t_p_2 *
         (-548 + 9501 * s + 4042 * t - 68064 * s * t - 44944 * s_p_2 +
          288348 * t * s_p_2 + 79335 * s_p_3 - 410220 * t * s_p_3 -
          54864 * s_p_4 + 183384 * t * s_p_4 + 11340 * s_p_5 - 6764 * t_p_2 +
          113301 * s * t_p_2 - 377568 * s_p_2 * t_p_2 + 293760 * s_p_3 * t_p_2 +
          2334 * t_p_3 - 60804 * s * t_p_3 + 134784 * s_p_2 * t_p_3 +
          2016 * t_p_4 + 6156 * s * t_p_4 - 1080 * t_p_5)) /
        32.;

      dphi(30, 1) =
        (-81 * t * s_p_2 *
         (-548 + 6334 * s + 6063 * t - 68064 * s * t - 22472 * s_p_2 +
          216261 * t * s_p_2 + 31734 * s_p_3 - 246132 * t * s_p_3 -
          18288 * s_p_4 + 91692 * t * s_p_4 + 3240 * s_p_5 - 13528 * t_p_2 +
          151068 * s * t_p_2 - 377568 * s_p_2 * t_p_2 + 235008 * s_p_3 * t_p_2 +
          5835 * t_p_3 - 101340 * s * t_p_3 + 168480 * s_p_2 * t_p_3 +
          6048 * t_p_4 + 12312 * s * t_p_4 - 3780 * t_p_5)) /
        32.;

      dphi(31, 0) =
        (81 * s * t_p_2 *
         (548 - 6063 * s - 6334 * t + 68064 * s * t + 13528 * s_p_2 -
          151068 * t * s_p_2 - 5835 * s_p_3 + 101340 * t * s_p_3 -
          6048 * s_p_4 - 12312 * t * s_p_4 + 3780 * s_p_5 + 22472 * t_p_2 -
          216261 * s * t_p_2 + 377568 * s_p_2 * t_p_2 - 168480 * s_p_3 * t_p_2 -
          31734 * t_p_3 + 246132 * s * t_p_3 - 235008 * s_p_2 * t_p_3 +
          18288 * t_p_4 - 91692 * s * t_p_4 - 3240 * t_p_5)) /
        32.;

      dphi(31, 1) =
        (81 * t * s_p_2 *
         (548 - 4042 * s - 9501 * t + 68064 * s * t + 6764 * s_p_2 -
          113301 * t * s_p_2 - 2334 * s_p_3 + 60804 * t * s_p_3 - 2016 * s_p_4 -
          6156 * t * s_p_4 + 1080 * s_p_5 + 44944 * t_p_2 - 288348 * s * t_p_2 +
          377568 * s_p_2 * t_p_2 - 134784 * s_p_3 * t_p_2 - 79335 * t_p_3 +
          410220 * s * t_p_3 - 293760 * s_p_2 * t_p_3 + 54864 * t_p_4 -
          183384 * s * t_p_4 - 11340 * t_p_5)) /
        32.;

      dphi(32, 0) =
        (-162 * s * t_p_2 *
         (-1036 - 3939 * s + 14966 * t + 42000 * s * t + 18160 * s_p_2 -
          227004 * t * s_p_2 - 5805 * s_p_3 + 202320 * t * s_p_3 -
          19008 * s_p_4 - 31752 * t * s_p_4 + 11340 * s_p_5 - 81400 * t_p_2 -
          65673 * s * t_p_2 + 504576 * s_p_2 * t_p_2 - 316440 * s_p_3 * t_p_2 +
          194118 * t_p_3 - 60048 * s * t_p_3 - 235872 * s_p_2 * t_p_3 -
          201600 * t_p_4 + 89100 * s * t_p_4 + 74952 * t_p_5)) /
        125.;

      dphi(32, 1) =
        (-162 * t * s_p_2 *
         (-1036 - 2626 * s + 22449 * t + 42000 * s * t + 9080 * s_p_2 -
          170253 * t * s_p_2 - 2322 * s_p_3 + 121392 * t * s_p_3 -
          6336 * s_p_4 - 15876 * t * s_p_4 + 3240 * s_p_5 - 162800 * t_p_2 -
          87564 * s * t_p_2 + 504576 * s_p_2 * t_p_2 - 253152 * s_p_3 * t_p_2 +
          485295 * t_p_3 - 100080 * s * t_p_3 - 294840 * s_p_2 * t_p_3 -
          604800 * t_p_4 + 178200 * s * t_p_4 + 262332 * t_p_5)) /
        125.;

      dphi(33, 0) = (-648 * (-1 + s + t) * (-1 + 3 * s + t) * (-1 + 2 * t) *
                     (-2 + 3 * t) * (-1 + 3 * t) * (-1 + 6 * t) * t_p_2) /
                    25.;


      dphi(33, 1) = (-648 * s * t * (-1 + s + t) *
                     (-4 + 4 * s + 83 * t - 75 * s * t - 545 * t_p_2 +
                      420 * s * t_p_2 + 1530 * t_p_3 - 900 * s * t_p_3 -
                      1908 * t_p_4 + 648 * s * t_p_4 + 864 * t_p_5)) /
                    25.;

      dphi(34, 0) =
        (81 * (-1 + s + t) * (-1 + 2 * t) * (-1 + 3 * t) * (-1 + 6 * t) *
         t_p_2 * (5 - 27 * s - 11 * t + 30 * s * t + 24 * s_p_2 + 6 * t_p_2)) /
        4.;

      dphi(34, 1) =
        (81 * s * t * (-1 + s + t) *
         (-10 + 22 * s + 203 * t - 405 * s * t - 12 * s_p_2 + 198 * t * s_p_2 -
          1289 * t_p_2 + 2178 * s * t_p_2 - 864 * s_p_2 * t_p_2 + 3456 * t_p_3 -
          4356 * s * t_p_3 + 1080 * s_p_2 * t_p_3 - 4068 * t_p_4 +
          2808 * s * t_p_4 + 1728 * t_p_5)) /
        4.;


      dphi(35, 0) =
        (81 * (-1 + s + t) * (-1 + 6 * t) * t_p_2 *
         (10 - 124 * s - 57 * t + 523 * s * t + 404 * s_p_2 - 1116 * t * s_p_2 -
          504 * s_p_3 + 684 * t * s_p_3 + 216 * s_p_4 + 119 * t_p_2 -
          720 * s * t_p_2 + 756 * s_p_2 * t_p_2 - 108 * t_p_3 +
          324 * s * t_p_3 + 36 * t_p_4)) /
        4.;

      dphi(35, 1) =
        (81 * s * t * (-1 + s + t) *
         (-20 + 114 * s + 361 * t - 1787 * s * t - 238 * s_p_2 +
          3186 * t * s_p_2 + 216 * s_p_3 - 2412 * t * s_p_3 - 72 * s_p_4 +
          648 * t * s_p_4 - 1951 * t_p_2 + 7434 * s * t_p_2 -
          9180 * s_p_2 * t_p_2 + 3672 * s_p_3 * t_p_2 + 4464 * t_p_3 -
          11412 * s * t_p_3 + 7128 * s_p_2 * t_p_3 - 4572 * t_p_4 +
          5832 * s * t_p_4 + 1728 * t_p_5)) /
        4.;

      dphi(36, 0) =
        (-648 * (-1 + s + t) * t_p_2 *
         (-10 + 184 * s + 87 * t - 1237 * s * t - 947 * s_p_2 +
          4611 * t * s_p_2 + 2073 * s_p_3 - 6516 * t * s_p_3 - 2052 * s_p_4 +
          3132 * t * s_p_4 + 756 * s_p_5 - 290 * t_p_2 + 3003 * s * t_p_2 -
          7236 * s_p_2 * t_p_2 + 4968 * s_p_3 * t_p_2 + 465 * t_p_3 -
          3132 * s * t_p_3 + 3672 * s_p_2 * t_p_3 - 360 * t_p_4 +
          1188 * s * t_p_4 + 108 * t_p_5)) /
        25.;


      dphi(36, 1) =
        (-648 * s * t * (-1 + s + t) *
         (-20 + 174 * s + 271 * t - 1817 * s * t - 580 * s_p_2 +
          4398 * t * s_p_2 + 930 * s_p_3 - 4572 * t * s_p_3 - 720 * s_p_4 +
          1728 * t * s_p_4 + 216 * s_p_5 - 1237 * t_p_2 + 6006 * s * t_p_2 -
          9396 * s_p_2 * t_p_2 + 4752 * s_p_3 * t_p_2 + 2538 * t_p_3 -
          7956 * s * t_p_3 + 6048 * s_p_2 * t_p_3 - 2412 * t_p_4 +
          3672 * s * t_p_4 + 864 * t_p_5)) /
        25.;

      dphi(37, 0) =
        (-648 * s * t * (-1 + s + t) *
         (-20 + 271 * s + 174 * t - 1817 * s * t - 1237 * s_p_2 +
          6006 * t * s_p_2 + 2538 * s_p_3 - 7956 * t * s_p_3 - 2412 * s_p_4 +
          3672 * t * s_p_4 + 864 * s_p_5 - 580 * t_p_2 + 4398 * s * t_p_2 -
          9396 * s_p_2 * t_p_2 + 6048 * s_p_3 * t_p_2 + 930 * t_p_3 -
          4572 * s * t_p_3 + 4752 * s_p_2 * t_p_3 - 720 * t_p_4 +
          1728 * s * t_p_4 + 216 * t_p_5)) /
        25.;

      dphi(37, 1) =
        (-648 * (-1 + s + t) * s_p_2 *
         (-10 + 87 * s + 184 * t - 1237 * s * t - 290 * s_p_2 +
          3003 * t * s_p_2 + 465 * s_p_3 - 3132 * t * s_p_3 - 360 * s_p_4 +
          1188 * t * s_p_4 + 108 * s_p_5 - 947 * t_p_2 + 4611 * s * t_p_2 -
          7236 * s_p_2 * t_p_2 + 3672 * s_p_3 * t_p_2 + 2073 * t_p_3 -
          6516 * s * t_p_3 + 4968 * s_p_2 * t_p_3 - 2052 * t_p_4 +
          3132 * s * t_p_4 + 756 * t_p_5)) /
        25.;

      dphi(38, 0) = (81 * s * t * (-1 + s + t) *
                     (-20 + 361 * s + 114 * t - 1787 * s * t - 1951 * s_p_2 +
                      7434 * t * s_p_2 + 4464 * s_p_3 - 11412 * t * s_p_3 -
                      4572 * s_p_4 + 5832 * t * s_p_4 + 1728 * s_p_5 -
                      238 * t_p_2 + 3186 * s * t_p_2 - 9180 * s_p_2 * t_p_2 +
                      7128 * s_p_3 * t_p_2 + 216 * t_p_3 - 2412 * s * t_p_3 +
                      3672 * s_p_2 * t_p_3 - 72 * t_p_4 + 648 * s * t_p_4)) /
                    4.;

      dphi(38, 1) =
        (81 * (-1 + 6 * s) * (-1 + s + t) * s_p_2 *
         (10 - 57 * s - 124 * t + 523 * s * t + 119 * s_p_2 - 720 * t * s_p_2 -
          108 * s_p_3 + 324 * t * s_p_3 + 36 * s_p_4 + 404 * t_p_2 -
          1116 * s * t_p_2 + 756 * s_p_2 * t_p_2 - 504 * t_p_3 +
          684 * s * t_p_3 + 216 * t_p_4)) /
        4.;

      dphi(39, 0) =
        (81 * s * t * (-1 + s + t) *
         (-10 + 203 * s + 22 * t - 405 * s * t - 1289 * s_p_2 +
          2178 * t * s_p_2 + 3456 * s_p_3 - 4356 * t * s_p_3 - 4068 * s_p_4 +
          2808 * t * s_p_4 + 1728 * s_p_5 - 12 * t_p_2 + 198 * s * t_p_2 -
          864 * s_p_2 * t_p_2 + 1080 * s_p_3 * t_p_2)) /
        4.;

      dphi(39, 1) =
        (81 * (-1 + 2 * s) * (-1 + 3 * s) * (-1 + 6 * s) * (-1 + s + t) *
         s_p_2 * (5 - 11 * s - 27 * t + 30 * s * t + 6 * s_p_2 + 24 * t_p_2)) /
        4.;

      dphi(40, 0) = (-648 * s * t * (-1 + s + t) *
                     (-4 + 83 * s + 4 * t - 75 * s * t - 545 * s_p_2 +
                      420 * t * s_p_2 + 1530 * s_p_3 - 900 * t * s_p_3 -
                      1908 * s_p_4 + 648 * t * s_p_4 + 864 * s_p_5)) /
                    25.;

      dphi(40, 1) = (-648 * (-1 + 2 * s) * (-2 + 3 * s) * (-1 + 3 * s) *
                     (-1 + 6 * s) * (-1 + s + t) * (-1 + s + 3 * t) * s_p_2) /
                    25.;


      dphi(41, 0) = (324 * s *
                     (-4 + 81 * s + 4 * t - 75 * s * t - 520 * s_p_2 +
                      420 * t * s_p_2 + 1425 * s_p_3 - 900 * t * s_p_3 -
                      1728 * s_p_4 + 648 * t * s_p_4 + 756 * s_p_5) *
                     t_p_2 * sqrt(2)) /
                    25.;

      dphi(41, 1) =
        (324 * (-1 + 2 * s) * (-2 + 3 * s) * (-1 + 3 * s) * (-1 + 6 * s) * t *
         (-2 + 2 * s + 3 * t) * s_p_2 * sqrt(2)) /
        25.;

      dphi(42, 0) =
        (81 * s * (-1 + 6 * t) *
         (2 - 36 * s - 2 * t + 33 * s * t + 188 * s_p_2 - 144 * t * s_p_2 -
          360 * s_p_3 + 180 * t * s_p_3 + 216 * s_p_4) *
         t_p_2 * pow(sqrt(2), -1)) /
        4.;


      dphi(42, 1) =
        (81 * (-1 + 2 * s) * (-1 + 3 * s) * (-1 + 6 * s) * t * s_p_2 *
         (2 - 2 * s - 21 * t + 18 * s * t + 24 * t_p_2) * pow(sqrt(2), -1)) /
        4.;

      dphi(43, 0) = (81 * s * (-1 + 2 * t) * (-1 + 3 * t) * (-1 + 6 * t) *
                     (2 - 21 * s - 2 * t + 18 * s * t + 24 * s_p_2) * t_p_2 *
                     pow(sqrt(2), -1)) /
                    4.;

      dphi(43, 1) =
        (81 * (-1 + 6 * s) * t * s_p_2 *
         (2 - 2 * s - 36 * t + 33 * s * t + 188 * t_p_2 - 144 * s * t_p_2 -
          360 * t_p_3 + 180 * s * t_p_3 + 216 * t_p_4) *
         pow(sqrt(2), -1)) /
        4.;

      dphi(44, 0) =
        (324 * s * (-1 + 2 * t) * (-2 + 3 * s + 2 * t) * (-2 + 3 * t) *
         (-1 + 3 * t) * (-1 + 6 * t) * t_p_2 * sqrt(2)) /
        25.;

      dphi(44, 1) = (324 * t * s_p_2 *
                     (-4 + 4 * s + 81 * t - 75 * s * t - 520 * t_p_2 +
                      420 * s * t_p_2 + 1425 * t_p_3 - 900 * s * t_p_3 -
                      1728 * t_p_4 + 648 * s * t_p_4 + 756 * t_p_5) *
                     sqrt(2)) /
                    25.;

      dphi(45, 0) = 5832 * s * (-1 + s + t) * (-1 + 2 * s + t) * (-1 + 2 * t) *
                    (-1 + 3 * t) * (-1 + 6 * t) * t_p_2;

      dphi(45, 1) =
        2916 * t * (-1 + s + t) * s_p_2 *
        (2 - 2 * s - 37 * t + 33 * s * t + 199 * t_p_2 - 144 * s * t_p_2 -
         396 * t_p_3 + 180 * s * t_p_3 + 252 * t_p_4);

      dphi(46, 0) = -2592 * s * (-1 + s + t) * (-1 + 3 * t) * (-1 + 6 * t) *
                    t_p_2 *
                    (5 - 19 * s - 11 * t + 21 * s * t + 15 * s_p_2 + 6 * t_p_2);


      dphi(46, 1) =
        -1296 * t * (-1 + s + t) * s_p_2 *
        (10 - 22 * s - 173 * t + 339 * s * t + 12 * s_p_2 - 162 * t * s_p_2 +
         831 * t_p_2 - 1278 * s * t_p_2 + 432 * s_p_2 * t_p_2 - 1404 * t_p_3 +
         1188 * s * t_p_3 + 756 * t_p_4);

      dphi(47, 0) =
        1296 * s * (-1 + s + t) * (-1 + 6 * t) * t_p_2 *
        (-20 + 121 * s + 74 * t - 297 * s * t - 207 * s_p_2 + 252 * t * s_p_2 +
         108 * s_p_3 - 90 * t_p_2 + 180 * s * t_p_2 + 36 * t_p_3);

      dphi(47, 1) =
        1296 * t * (-1 + s + t) * s_p_2 *
        (20 - 74 * s - 301 * t + 963 * s * t + 90 * s_p_2 - 990 * t * s_p_2 -
         36 * s_p_3 + 324 * t * s_p_3 + 1155 * t_p_2 - 2574 * s * t_p_2 +
         1404 * s_p_2 * t_p_2 - 1620 * t_p_3 + 1836 * s * t_p_3 + 756 * t_p_4);

      dphi(48, 0) =
        -2916 * s * (-1 + s + t) * t_p_2 *
        (20 - 181 * s - 114 * t + 761 * s * t + 523 * s_p_2 - 1440 * t * s_p_2 -
         612 * s_p_3 + 828 * t * s_p_3 + 252 * s_p_4 + 238 * t_p_2 -
         1044 * s * t_p_2 + 972 * s_p_2 * t_p_2 - 216 * t_p_3 +
         468 * s * t_p_3 + 72 * t_p_4);


      dphi(48, 1) =
        -2916 * t * (-1 + s + t) * s_p_2 *
        (20 - 114 * s - 181 * t + 761 * s * t + 238 * s_p_2 - 1044 * t * s_p_2 -
         216 * s_p_3 + 468 * t * s_p_3 + 72 * s_p_4 + 523 * t_p_2 -
         1440 * s * t_p_2 + 972 * s_p_2 * t_p_2 - 612 * t_p_3 +
         828 * s * t_p_3 + 252 * t_p_4);

      dphi(49, 0) =
        1296 * s * (-1 + s + t) * t_p_2 *
        (20 - 301 * s - 74 * t + 963 * s * t + 1155 * s_p_2 - 2574 * t * s_p_2 -
         1620 * s_p_3 + 1836 * t * s_p_3 + 756 * s_p_4 + 90 * t_p_2 -
         990 * s * t_p_2 + 1404 * s_p_2 * t_p_2 - 36 * t_p_3 + 324 * s * t_p_3);

      dphi(49, 1) =
        1296 * (-1 + 6 * s) * t * (-1 + s + t) * s_p_2 *
        (-20 + 74 * s + 121 * t - 297 * s * t - 90 * s_p_2 + 180 * t * s_p_2 +
         36 * s_p_3 - 207 * t_p_2 + 252 * s * t_p_2 + 108 * t_p_3);

      dphi(50, 0) =
        -1296 * s * (-1 + s + t) * t_p_2 *
        (10 - 173 * s - 22 * t + 339 * s * t + 831 * s_p_2 - 1278 * t * s_p_2 -
         1404 * s_p_3 + 1188 * t * s_p_3 + 756 * s_p_4 + 12 * t_p_2 -
         162 * s * t_p_2 + 432 * s_p_2 * t_p_2);


      dphi(50, 1) = -2592 * (-1 + 3 * s) * (-1 + 6 * s) * t * (-1 + s + t) *
                    s_p_2 *
                    (5 - 11 * s - 19 * t + 21 * s * t + 6 * s_p_2 + 15 * t_p_2);

      dphi(51, 0) =
        2916 * s * (-1 + s + t) *
        (2 - 37 * s - 2 * t + 33 * s * t + 199 * s_p_2 - 144 * t * s_p_2 -
         396 * s_p_3 + 180 * t * s_p_3 + 252 * s_p_4) *
        t_p_2;

      dphi(51, 1) = 5832 * (-1 + 2 * s) * (-1 + 3 * s) * (-1 + 6 * s) * t *
                    (-1 + s + t) * (-1 + s + 2 * t) * s_p_2;

      dphi(52, 0) = 1296 * s * (-1 + s + t) * (-1 + 6 * t) *
                    (-2 + 31 * s + 2 * t - 27 * s * t - 117 * s_p_2 +
                     72 * t * s_p_2 + 108 * s_p_3) *
                    t_p_2;

      dphi(52, 1) = 2592 * (-1 + 3 * s) * (-1 + 6 * s) * t * (-1 + s + t) *
                    s_p_2 * (1 - s - 11 * t + 9 * s * t + 15 * t_p_2);

      dphi(53, 0) = 2592 * s * (-1 + s + t) * (-1 + 3 * t) * (-1 + 6 * t) *
                    (1 - 11 * s - t + 9 * s * t + 15 * s_p_2) * t_p_2;


      dphi(53, 1) = 1296 * (-1 + 6 * s) * t * (-1 + s + t) * s_p_2 *
                    (-2 + 2 * s + 31 * t - 27 * s * t - 117 * t_p_2 +
                     72 * s * t_p_2 + 108 * t_p_3);

      dphi(54, 0) = -1458 * s * (-1 + s + t) * (-1 + 2 * s + t) * (-1 + 6 * t) *
                    (5 - 54 * s - 6 * t + 54 * s * t + 54 * s_p_2) * t_p_2;

      dphi(54, 1) = -1458 * (-1 + 6 * s) * t * (-1 + s + t) * (-1 + s + 2 * t) *
                    s_p_2 * (5 - 6 * s - 54 * t + 54 * s * t + 54 * t_p_2);
    }

    /// Precomputed d2basis polynomials for the 5th order boundary
    /// representation
    template<>
    void BernadouElementBasis<5>::d2full_basic_polynomials(
      const Vector<double>& s_basic, DShape& d2phi) const
    {
      // For convenience
      const double s = s_basic[0], t = s_basic[1];
      // Precompute
      const double s_p_2 = pow(s, 2), s_p_3 = pow(s, 3), s_p_4 = pow(s, 4),
                   s_p_5 = pow(s, 5), s_p_6 = pow(s, 6), s_p_7 = pow(s, 7);
      const double t_p_2 = pow(t, 2), t_p_3 = pow(t, 3), t_p_4 = pow(t, 4),
                   t_p_5 = pow(t, 5), t_p_6 = pow(t, 6), t_p_7 = pow(t, 7);

      // Now fill in (automatically generated code)
      d2phi(0, 0) =
        (3 *
         (598360 * s - 16261464 * s_p_2 + 133811000 * s_p_3 -
          474346260 * s_p_4 + 813567132 * s_p_5 - 665594496 * s_p_6 +
          208225728 * s_p_7 - 3415480 * t_p_2 + 126762915 * s * t_p_2 -
          1077009696 * s_p_2 * t_p_2 + 3315149910 * s_p_3 * t_p_2 -
          4117183020 * s_p_4 * t_p_2 + 1756160028 * s_p_5 * t_p_2 +
          7389313 * t_p_3 - 247491396 * s * t_p_3 + 1848699990 * s_p_2 * t_p_3 -
          4323410640 * s_p_3 * t_p_3 + 2960191980 * s_p_4 * t_p_3 -
          4195422 * t_p_4 + 92584809 * s * t_p_4 - 685821816 * s_p_2 * t_p_4 +
          1092981600 * s_p_3 * t_p_4 + 1836081 * t_p_5 + 42791328 * s * t_p_5 -
          67530672 * s_p_2 * t_p_5 - 3266568 * t_p_6 - 15545196 * s * t_p_6 +
          1652076 * t_p_7)) /
        2000.;

      d2phi(0, 1) =
        (3 * s * t *
         (-6830960 + 126762915 * s + 22167939 * t - 371237094 * s * t -
          718006464 * s_p_2 + 1848699990 * t * s_p_2 + 1657574955 * s_p_3 -
          3242557980 * t * s_p_3 - 1646873208 * s_p_4 + 1776115188 * t * s_p_4 +
          585386676 * s_p_5 - 16781688 * t_p_2 + 185169618 * s * t_p_2 -
          914429088 * s_p_2 * t_p_2 + 1092981600 * s_p_3 * t_p_2 +
          9180405 * t_p_3 + 106978320 * s * t_p_3 - 112551120 * s_p_2 * t_p_3 -
          19599408 * t_p_4 - 46635588 * s * t_p_4 + 11564532 * t_p_5)) /
        2000.;

      d2phi(0, 2) =
        (3 * pow(s, 2) *
         (-3415480 + 42254305 * s + 22167939 * t - 247491396 * s * t -
          179501616 * s_p_2 + 924349995 * t * s_p_2 + 331514991 * s_p_3 -
          1297023192 * t * s_p_3 - 274478868 * s_p_4 + 592038396 * t * s_p_4 +
          83626668 * s_p_5 - 25172532 * t_p_2 + 185169618 * s * t_p_2 -
          685821816 * s_p_2 * t_p_2 + 655788960 * s_p_3 * t_p_2 +
          18360810 * t_p_3 + 142637760 * s * t_p_3 - 112551120 * s_p_2 * t_p_3 -
          48998520 * t_p_4 - 77725980 * s * t_p_4 + 34693596 * t_p_5)) /
        2000.;

      d2phi(1, 0) =
        (3 * t_p_2 *
         (-3415480 + 22167939 * s + 42254305 * t - 247491396 * s * t -
          25172532 * s_p_2 + 185169618 * t * s_p_2 + 18360810 * s_p_3 +
          142637760 * t * s_p_3 - 48998520 * s_p_4 - 77725980 * t * s_p_4 +
          34693596 * s_p_5 - 179501616 * t_p_2 + 924349995 * s * t_p_2 -
          685821816 * s_p_2 * t_p_2 - 112551120 * s_p_3 * t_p_2 +
          331514991 * t_p_3 - 1297023192 * s * t_p_3 +
          655788960 * s_p_2 * t_p_3 - 274478868 * t_p_4 +
          592038396 * s * t_p_4 + 83626668 * t_p_5)) /
        2000.;

      d2phi(1, 1) =
        (3 * s * t *
         (-6830960 + 22167939 * s + 126762915 * t - 371237094 * s * t -
          16781688 * s_p_2 + 185169618 * t * s_p_2 + 9180405 * s_p_3 +
          106978320 * t * s_p_3 - 19599408 * s_p_4 - 46635588 * t * s_p_4 +
          11564532 * s_p_5 - 718006464 * t_p_2 + 1848699990 * s * t_p_2 -
          914429088 * s_p_2 * t_p_2 - 112551120 * s_p_3 * t_p_2 +
          1657574955 * t_p_3 - 3242557980 * s * t_p_3 +
          1092981600 * s_p_2 * t_p_3 - 1646873208 * t_p_4 +
          1776115188 * s * t_p_4 + 585386676 * t_p_5)) /
        2000.;


      d2phi(1, 2) =
        (3 * (598360 * t - 3415480 * s_p_2 + 126762915 * t * s_p_2 +
              7389313 * s_p_3 - 247491396 * t * s_p_3 - 4195422 * s_p_4 +
              92584809 * t * s_p_4 + 1836081 * s_p_5 + 42791328 * t * s_p_5 -
              3266568 * s_p_6 - 15545196 * t * s_p_6 + 1652076 * s_p_7 -
              16261464 * t_p_2 - 1077009696 * s_p_2 * t_p_2 +
              1848699990 * s_p_3 * t_p_2 - 685821816 * s_p_4 * t_p_2 -
              67530672 * s_p_5 * t_p_2 + 133811000 * t_p_3 +
              3315149910 * s_p_2 * t_p_3 - 4323410640 * s_p_3 * t_p_3 +
              1092981600 * s_p_4 * t_p_3 - 474346260 * t_p_4 -
              4117183020 * s_p_2 * t_p_4 + 2960191980 * s_p_3 * t_p_4 +
              813567132 * t_p_5 + 1756160028 * s_p_2 * t_p_5 -
              665594496 * t_p_6 + 208225728 * t_p_7)) /
        2000.;

      d2phi(2, 0) =
        (-3 *
         (1986086 * s - 31794072 * s_p_2 + 174053920 * s_p_3 -
          448736220 * s_p_4 + 596370222 * s_p_5 - 395992800 * s_p_6 +
          104112864 * s_p_7 - 5622406 * t_p_2 + 81036633 * s * t_p_2 -
          205529076 * s_p_2 * t_p_2 - 90608490 * s_p_3 * t_p_2 +
          618729030 * s_p_4 * t_p_2 - 405776952 * s_p_5 * t_p_2 +
          27012211 * t_p_3 - 435713364 * s * t_p_3 +
          1363972338 * s_p_2 * t_p_3 - 1220397840 * s_p_3 * t_p_3 +
          181127340 * s_p_4 * t_p_3 - 34254846 * t_p_4 + 681986169 * s * t_p_4 -
          1789749432 * s_p_2 * t_p_4 + 1068108120 * s_p_3 * t_p_4 -
          9060849 * t_p_5 - 366119352 * s * t_p_5 + 640864872 * s_p_2 * t_p_5 +
          41248602 * t_p_6 + 36225468 * s * t_p_6 - 19322712 * t_p_7)) /
        1000.;

      d2phi(2, 1) =
        (3 * s * t *
         (11244812 - 81036633 * s - 81036633 * t + 653570046 * s * t +
          137019384 * s_p_2 - 1363972338 * t * s_p_2 + 45304245 * s_p_3 +
          915298380 * t * s_p_3 - 247491612 * s_p_4 - 108676404 * t * s_p_4 +
          135258984 * s_p_5 + 137019384 * t_p_2 - 1363972338 * s * t_p_2 +
          2386332576 * s_p_2 * t_p_2 - 1068108120 * s_p_3 * t_p_2 +
          45304245 * t_p_3 + 915298380 * s * t_p_3 -
          1068108120 * s_p_2 * t_p_3 - 247491612 * t_p_4 -
          108676404 * s * t_p_4 + 135258984 * t_p_5)) /
        1000.;


      d2phi(2, 2) =
        (3 * (-1986086 * t + 5622406 * s_p_2 - 81036633 * t * s_p_2 -
              27012211 * s_p_3 + 435713364 * t * s_p_3 + 34254846 * s_p_4 -
              681986169 * t * s_p_4 + 9060849 * s_p_5 + 366119352 * t * s_p_5 -
              41248602 * s_p_6 - 36225468 * t * s_p_6 + 19322712 * s_p_7 +
              31794072 * t_p_2 + 205529076 * s_p_2 * t_p_2 -
              1363972338 * s_p_3 * t_p_2 + 1789749432 * s_p_4 * t_p_2 -
              640864872 * s_p_5 * t_p_2 - 174053920 * t_p_3 +
              90608490 * s_p_2 * t_p_3 + 1220397840 * s_p_3 * t_p_3 -
              1068108120 * s_p_4 * t_p_3 + 448736220 * t_p_4 -
              618729030 * s_p_2 * t_p_4 - 181127340 * s_p_3 * t_p_4 -
              596370222 * t_p_5 + 405776952 * s_p_2 * t_p_5 +
              395992800 * t_p_6 - 104112864 * t_p_7)) /
        1000.;

      d2phi(3, 0) =
        (-18840 * s + 513816 * s_p_2 - 4251000 * s_p_3 + 15179940 * s_p_4 -
         26265708 * s_p_5 + 21700224 * s_p_6 - 6858432 * s_p_7 +
         103738 * t_p_2 - 3854289 * s * t_p_2 + 32818824 * s_p_2 * t_p_2 -
         101369790 * s_p_3 * t_p_2 + 126520380 * s_p_4 * t_p_2 -
         54316332 * s_p_5 * t_p_2 - 222565 * t_p_3 + 7445574 * s * t_p_3 -
         55666710 * s_p_2 * t_p_3 + 130358160 * s_p_3 * t_p_3 -
         89506620 * s_p_4 * t_p_3 + 127188 * t_p_4 - 2759481 * s * t_p_4 +
         20491704 * s_p_2 * t_p_4 - 32594400 * s_p_3 * t_p_4 - 59229 * t_p_5 -
         1288872 * s * t_p_5 + 1986768 * s_p_2 * t_p_5 + 102384 * t_p_6 +
         475308 * s * t_p_6 - 51516 * t_p_7) /
        200.;

      d2phi(3, 1) =
        -(s * t *
          (-207476 + 3854289 * s + 667695 * t - 11168361 * s * t -
           21879216 * s_p_2 + 55666710 * t * s_p_2 + 50684895 * s_p_3 -
           97768620 * t * s_p_3 - 50608152 * s_p_4 + 53703972 * t * s_p_4 +
           18105444 * s_p_5 - 508752 * t_p_2 + 5518962 * s * t_p_2 -
           27322272 * s_p_2 * t_p_2 + 32594400 * s_p_3 * t_p_2 +
           296145 * t_p_3 + 3222180 * s * t_p_3 - 3311280 * s_p_2 * t_p_3 -
           614304 * t_p_4 - 1425924 * s * t_p_4 + 360612 * t_p_5)) /
        200.;

      d2phi(3, 2) =
        -(s_p_2 *
          (-103738 + 1284763 * s + 667695 * t - 7445574 * s * t -
           5469804 * s_p_2 + 27833355 * t * s_p_2 + 10136979 * s_p_3 -
           39107448 * t * s_p_3 - 8434692 * s_p_4 + 17901324 * t * s_p_4 +
           2586492 * s_p_5 - 763128 * t_p_2 + 5518962 * s * t_p_2 -
           20491704 * s_p_2 * t_p_2 + 19556640 * s_p_3 * t_p_2 +
           592290 * t_p_3 + 4296240 * s * t_p_3 - 3311280 * s_p_2 * t_p_3 -
           1535760 * t_p_4 - 2376540 * s * t_p_4 + 1081836 * t_p_5)) /
        200.;

      d2phi(4, 0) =
        (t * (-6680 + 293388 * s + 9488 * t - 345093 * s * t - 3221616 * s_p_2 +
              2513700 * t * s_p_2 + 14459400 * s_p_3 - 6876270 * t * s_p_3 -
              30645000 * s_p_4 + 9020160 * t * s_p_4 + 30454704 * s_p_5 -
              4769604 * t * s_p_5 - 11394432 * s_p_6 - 55323 * t_p_2 +
              1915110 * s * t_p_2 - 9840042 * s_p_2 * t_p_2 +
              11945880 * s_p_3 * t_p_2 - 3018060 * s_p_4 * t_p_2 +
              122310 * t_p_3 - 4806621 * s * t_p_3 + 21387888 * s_p_2 * t_p_3 -
              17184960 * s_p_3 * t_p_3 - 18927 * t_p_4 + 3399084 * s * t_p_4 -
              10427616 * s_p_2 * t_p_4 - 102384 * t_p_5 - 475308 * s * t_p_5 +
              51516 * t_p_6)) /
        200.;

      d2phi(4, 1) =
        -(s * (6680 - 146694 * s - 18976 * t + 345093 * s * t +
               1073872 * s_p_2 - 1675800 * t * s_p_2 - 3614850 * s_p_3 +
               3438135 * t * s_p_3 + 6129000 * s_p_4 - 3608064 * t * s_p_4 -
               5075784 * s_p_5 + 1589868 * t * s_p_5 + 1627776 * s_p_6 +
               165969 * t_p_2 - 2872665 * s * t_p_2 + 9840042 * s_p_2 * t_p_2 -
               8959410 * s_p_3 * t_p_2 + 1810836 * s_p_4 * t_p_2 -
               489240 * t_p_3 + 9613242 * s * t_p_3 - 28517184 * s_p_2 * t_p_3 +
               17184960 * s_p_3 * t_p_3 + 94635 * t_p_4 - 8497710 * s * t_p_4 +
               17379360 * s_p_2 * t_p_4 + 614304 * t_p_5 + 1425924 * s * t_p_5 -
               360612 * t_p_6)) /
        200.;

      d2phi(4, 2) =
        -(s_p_2 *
          (-9488 + 115031 * s + 165969 * t - 1915110 * s * t - 418950 * s_p_2 +
           4920021 * t * s_p_2 + 687627 * s_p_3 - 3583764 * t * s_p_3 -
           601344 * s_p_4 + 603612 * t * s_p_4 + 227124 * s_p_5 -
           733860 * t_p_2 + 9613242 * s * t_p_2 - 21387888 * s_p_2 * t_p_2 +
           10310976 * s_p_3 * t_p_2 + 189270 * t_p_3 - 11330280 * s * t_p_3 +
           17379360 * s_p_2 * t_p_3 + 1535760 * t_p_4 + 2376540 * s * t_p_4 -
           1081836 * t_p_5)) /
        200.;

      d2phi(5, 0) =
        -(t_p_2 *
          (-9488 + 165969 * s + 115031 * t - 1915110 * s * t - 733860 * s_p_2 +
           9613242 * t * s_p_2 + 189270 * s_p_3 - 11330280 * t * s_p_3 +
           1535760 * s_p_4 + 2376540 * t * s_p_4 - 1081836 * s_p_5 -
           418950 * t_p_2 + 4920021 * s * t_p_2 - 21387888 * s_p_2 * t_p_2 +
           17379360 * s_p_3 * t_p_2 + 687627 * t_p_3 - 3583764 * s * t_p_3 +
           10310976 * s_p_2 * t_p_3 - 601344 * t_p_4 + 603612 * s * t_p_4 +
           227124 * t_p_5)) /
        200.;

      d2phi(5, 1) =
        -(t * (6680 - 18976 * s - 146694 * t + 345093 * s * t + 165969 * s_p_2 -
               2872665 * t * s_p_2 - 489240 * s_p_3 + 9613242 * t * s_p_3 +
               94635 * s_p_4 - 8497710 * t * s_p_4 + 614304 * s_p_5 +
               1425924 * t * s_p_5 - 360612 * s_p_6 + 1073872 * t_p_2 -
               1675800 * s * t_p_2 + 9840042 * s_p_2 * t_p_2 -
               28517184 * s_p_3 * t_p_2 + 17379360 * s_p_4 * t_p_2 -
               3614850 * t_p_3 + 3438135 * s * t_p_3 - 8959410 * s_p_2 * t_p_3 +
               17184960 * s_p_3 * t_p_3 + 6129000 * t_p_4 -
               3608064 * s * t_p_4 + 1810836 * s_p_2 * t_p_4 - 5075784 * t_p_5 +
               1589868 * s * t_p_5 + 1627776 * t_p_6)) /
        200.;

      d2phi(5, 2) =
        (s * (-6680 + 9488 * s + 293388 * t - 345093 * s * t - 55323 * s_p_2 +
              1915110 * t * s_p_2 + 122310 * s_p_3 - 4806621 * t * s_p_3 -
              18927 * s_p_4 + 3399084 * t * s_p_4 - 102384 * s_p_5 -
              475308 * t * s_p_5 + 51516 * s_p_6 - 3221616 * t_p_2 +
              2513700 * s * t_p_2 - 9840042 * s_p_2 * t_p_2 +
              21387888 * s_p_3 * t_p_2 - 10427616 * s_p_4 * t_p_2 +
              14459400 * t_p_3 - 6876270 * s * t_p_3 +
              11945880 * s_p_2 * t_p_3 - 17184960 * s_p_3 * t_p_3 -
              30645000 * t_p_4 + 9020160 * s * t_p_4 - 3018060 * s_p_2 * t_p_4 +
              30454704 * t_p_5 - 4769604 * s * t_p_5 - 11394432 * t_p_6)) /
        200.;

      d2phi(6, 0) =
        -(t_p_2 *
          (-103738 + 667695 * s + 1284763 * t - 7445574 * s * t -
           763128 * s_p_2 + 5518962 * t * s_p_2 + 592290 * s_p_3 +
           4296240 * t * s_p_3 - 1535760 * s_p_4 - 2376540 * t * s_p_4 +
           1081836 * s_p_5 - 5469804 * t_p_2 + 27833355 * s * t_p_2 -
           20491704 * s_p_2 * t_p_2 - 3311280 * s_p_3 * t_p_2 +
           10136979 * t_p_3 - 39107448 * s * t_p_3 + 19556640 * s_p_2 * t_p_3 -
           8434692 * t_p_4 + 17901324 * s * t_p_4 + 2586492 * t_p_5)) /
        200.;

      d2phi(6, 1) =
        -(s * t *
          (-207476 + 667695 * s + 3854289 * t - 11168361 * s * t -
           508752 * s_p_2 + 5518962 * t * s_p_2 + 296145 * s_p_3 +
           3222180 * t * s_p_3 - 614304 * s_p_4 - 1425924 * t * s_p_4 +
           360612 * s_p_5 - 21879216 * t_p_2 + 55666710 * s * t_p_2 -
           27322272 * s_p_2 * t_p_2 - 3311280 * s_p_3 * t_p_2 +
           50684895 * t_p_3 - 97768620 * s * t_p_3 + 32594400 * s_p_2 * t_p_3 -
           50608152 * t_p_4 + 53703972 * s * t_p_4 + 18105444 * t_p_5)) /
        200.;


      d2phi(6, 2) =
        (-18840 * t + 103738 * s_p_2 - 3854289 * t * s_p_2 - 222565 * s_p_3 +
         7445574 * t * s_p_3 + 127188 * s_p_4 - 2759481 * t * s_p_4 -
         59229 * s_p_5 - 1288872 * t * s_p_5 + 102384 * s_p_6 +
         475308 * t * s_p_6 - 51516 * s_p_7 + 513816 * t_p_2 +
         32818824 * s_p_2 * t_p_2 - 55666710 * s_p_3 * t_p_2 +
         20491704 * s_p_4 * t_p_2 + 1986768 * s_p_5 * t_p_2 - 4251000 * t_p_3 -
         101369790 * s_p_2 * t_p_3 + 130358160 * s_p_3 * t_p_3 -
         32594400 * s_p_4 * t_p_3 + 15179940 * t_p_4 +
         126520380 * s_p_2 * t_p_4 - 89506620 * s_p_3 * t_p_4 -
         26265708 * t_p_5 - 54316332 * s_p_2 * t_p_5 + 21700224 * t_p_6 -
         6858432 * t_p_7) /
        200.;

      d2phi(7, 0) =
        (-80334 * s + 1170168 * s_p_2 - 6114480 * s_p_3 + 15345180 * s_p_4 -
         20045718 * s_p_5 + 13154400 * s_p_6 - 3429216 * s_p_7 +
         276743 * t_p_2 - 2233623 * s * t_p_2 + 1352538 * s_p_2 * t_p_2 +
         16665030 * s_p_3 * t_p_2 - 33653070 * s_p_4 * t_p_2 +
         17840088 * s_p_5 * t_p_2 - 1963172 * t_p_3 + 16216308 * s * t_p_3 -
         33488424 * s_p_2 * t_p_3 + 13087980 * s_p_3 * t_p_3 +
         8937540 * s_p_4 * t_p_3 + 5469291 * t_p_4 - 36284409 * s * t_p_4 +
         60902604 * s_p_2 * t_p_4 - 25654320 * s_p_3 * t_p_4 - 7495902 * t_p_5 +
         33674454 * s * t_p_5 - 30188376 * s_p_2 * t_p_5 + 5059584 * t_p_6 -
         11302416 * s * t_p_6 - 1346544 * t_p_7) /
        100.;

      d2phi(7, 1) =
        -(t *
          (30118 - 553486 * s - 361179 * t + 5889516 * s * t + 2233623 * s_p_2 -
           24324462 * t * s_p_2 - 901692 * s_p_3 + 33488424 * t * s_p_3 -
           8332515 * s_p_4 - 9815985 * t * s_p_4 + 13461228 * s_p_5 -
           5362524 * t * s_p_5 - 5946696 * s_p_6 + 1680476 * t_p_2 -
           21877164 * s * t_p_2 + 72568818 * s_p_2 * t_p_2 -
           81203472 * s_p_3 * t_p_2 + 25654320 * s_p_4 * t_p_2 -
           3932775 * t_p_3 + 37479510 * s * t_p_3 - 84186135 * s_p_2 * t_p_3 +
           50313960 * s_p_3 * t_p_3 + 4928796 * t_p_4 - 30357504 * s * t_p_4 +
           33907248 * s_p_2 * t_p_4 - 3159324 * t_p_5 + 9425808 * s * t_p_5 +
           813888 * t_p_6)) /
        100.;

      d2phi(7, 2) =
        (s *
         (-30118 + 276743 * s + 722358 * t - 5889516 * s * t - 744541 * s_p_2 +
          16216308 * t * s_p_2 + 225423 * s_p_3 - 16744212 * t * s_p_3 +
          1666503 * s_p_4 + 3926394 * t * s_p_4 - 2243538 * s_p_5 +
          1787508 * t * s_p_5 + 849528 * s_p_6 - 5041428 * t_p_2 +
          32815746 * s * t_p_2 - 72568818 * s_p_2 * t_p_2 +
          60902604 * s_p_3 * t_p_2 - 15392592 * s_p_4 * t_p_2 +
          15731100 * t_p_3 - 74959020 * s * t_p_3 + 112248180 * s_p_2 * t_p_3 -
          50313960 * s_p_3 * t_p_3 - 24643980 * t_p_4 + 75893760 * s * t_p_4 -
          56512080 * s_p_2 * t_p_4 + 18955944 * t_p_5 - 28277424 * s * t_p_5 -
          5697216 * t_p_6)) /
        100.;

      d2phi(8, 0) =
        (t * (-30118 + 722358 * s + 276743 * t - 5889516 * s * t -
              5041428 * s_p_2 + 32815746 * t * s_p_2 + 15731100 * s_p_3 -
              74959020 * t * s_p_3 - 24643980 * s_p_4 + 75893760 * t * s_p_4 +
              18955944 * s_p_5 - 28277424 * t * s_p_5 - 5697216 * s_p_6 -
              744541 * t_p_2 + 16216308 * s * t_p_2 - 72568818 * s_p_2 * t_p_2 +
              112248180 * s_p_3 * t_p_2 - 56512080 * s_p_4 * t_p_2 +
              225423 * t_p_3 - 16744212 * s * t_p_3 + 60902604 * s_p_2 * t_p_3 -
              50313960 * s_p_3 * t_p_3 + 1666503 * t_p_4 + 3926394 * s * t_p_4 -
              15392592 * s_p_2 * t_p_4 - 2243538 * t_p_5 + 1787508 * s * t_p_5 +
              849528 * t_p_6)) /
        100.;

      d2phi(8, 1) =
        -(s *
          (30118 - 361179 * s - 553486 * t + 5889516 * s * t + 1680476 * s_p_2 -
           21877164 * t * s_p_2 - 3932775 * s_p_3 + 37479510 * t * s_p_3 +
           4928796 * s_p_4 - 30357504 * t * s_p_4 - 3159324 * s_p_5 +
           9425808 * t * s_p_5 + 813888 * s_p_6 + 2233623 * t_p_2 -
           24324462 * s * t_p_2 + 72568818 * s_p_2 * t_p_2 -
           84186135 * s_p_3 * t_p_2 + 33907248 * s_p_4 * t_p_2 -
           901692 * t_p_3 + 33488424 * s * t_p_3 - 81203472 * s_p_2 * t_p_3 +
           50313960 * s_p_3 * t_p_3 - 8332515 * t_p_4 - 9815985 * s * t_p_4 +
           25654320 * s_p_2 * t_p_4 + 13461228 * t_p_5 - 5362524 * s * t_p_5 -
           5946696 * t_p_6)) /
        100.;

      d2phi(8, 2) =
        (-80334 * t + 276743 * s_p_2 - 2233623 * t * s_p_2 - 1963172 * s_p_3 +
         16216308 * t * s_p_3 + 5469291 * s_p_4 - 36284409 * t * s_p_4 -
         7495902 * s_p_5 + 33674454 * t * s_p_5 + 5059584 * s_p_6 -
         11302416 * t * s_p_6 - 1346544 * s_p_7 + 1170168 * t_p_2 +
         1352538 * s_p_2 * t_p_2 - 33488424 * s_p_3 * t_p_2 +
         60902604 * s_p_4 * t_p_2 - 30188376 * s_p_5 * t_p_2 - 6114480 * t_p_3 +
         16665030 * s_p_2 * t_p_3 + 13087980 * s_p_3 * t_p_3 -
         25654320 * s_p_4 * t_p_3 + 15345180 * t_p_4 -
         33653070 * s_p_2 * t_p_4 + 8937540 * s_p_3 * t_p_4 - 20045718 * t_p_5 +
         17840088 * s_p_2 * t_p_5 + 13154400 * t_p_6 - 3429216 * t_p_7) /
        100.;

      d2phi(9, 0) =
        (120 * s - 3288 * s_p_2 + 27400 * s_p_3 - 98820 * s_p_4 +
         173124 * s_p_5 - 145152 * s_p_6 + 46656 * s_p_7 - 628 * t_p_2 +
         23349 * s * t_p_2 - 199152 * s_p_2 * t_p_2 + 616770 * s_p_3 * t_p_2 -
         772740 * s_p_4 * t_p_2 + 333396 * s_p_5 * t_p_2 + 1339 * t_p_3 -
         44712 * s * t_p_3 + 334530 * s_p_2 * t_p_3 - 784080 * s_p_3 * t_p_3 +
         539460 * s_p_4 * t_p_3 - 774 * t_p_4 + 16443 * s * t_p_4 -
         122472 * s_p_2 * t_p_4 + 194400 * s_p_3 * t_p_4 + 387 * t_p_5 +
         7776 * s * t_p_5 - 11664 * s_p_2 * t_p_5 - 648 * t_p_6 -
         2916 * s * t_p_6 + 324 * t_p_7) /
        40.;

      d2phi(9, 1) =
        (s * t *
         (-1256 + 23349 * s + 4017 * t - 67068 * s * t - 132768 * s_p_2 +
          334530 * t * s_p_2 + 308385 * s_p_3 - 588060 * t * s_p_3 -
          309096 * s_p_4 + 323676 * t * s_p_4 + 111132 * s_p_5 - 3096 * t_p_2 +
          32886 * s * t_p_2 - 163296 * s_p_2 * t_p_2 + 194400 * s_p_3 * t_p_2 +
          1935 * t_p_3 + 19440 * s * t_p_3 - 19440 * s_p_2 * t_p_3 -
          3888 * t_p_4 - 8748 * s * t_p_4 + 2268 * t_p_5)) /
        40.;

      d2phi(9, 2) =
        (s_p_2 *
         (-628 + 7783 * s + 4017 * t - 44712 * s * t - 33192 * s_p_2 +
          167265 * t * s_p_2 + 61677 * s_p_3 - 235224 * t * s_p_3 -
          51516 * s_p_4 + 107892 * t * s_p_4 + 15876 * s_p_5 - 4644 * t_p_2 +
          32886 * s * t_p_2 - 122472 * s_p_2 * t_p_2 + 116640 * s_p_3 * t_p_2 +
          3870 * t_p_3 + 25920 * s * t_p_3 - 19440 * s_p_2 * t_p_3 -
          9720 * t_p_4 - 14580 * s * t_p_4 + 6804 * t_p_5)) /
        40.;

      d2phi(10, 0) =
        -(t *
          (-40 + 1764 * s + 50 * t - 1797 * s * t - 19488 * s_p_2 +
           12420 * t * s_p_2 + 88200 * s_p_3 - 31410 * t * s_p_3 -
           189000 * s_p_4 + 38880 * t * s_p_4 + 190512 * s_p_5 -
           20412 * t * s_p_5 - 72576 * s_p_6 - 325 * t_p_2 + 11340 * s * t_p_2 -
           57726 * s_p_2 * t_p_2 + 68040 * s_p_3 * t_p_2 -
           14580 * s_p_4 * t_p_2 + 720 * t_p_3 - 28863 * s * t_p_3 +
           128304 * s_p_2 * t_p_3 - 103680 * s_p_3 * t_p_3 - 81 * t_p_4 +
           20412 * s * t_p_4 - 62208 * s_p_2 * t_p_4 - 648 * t_p_5 -
           2916 * s * t_p_5 + 324 * t_p_6)) /
        20.;

      d2phi(10, 1) =
        (s *
         (40 - 882 * s - 100 * t + 1797 * s * t + 6496 * s_p_2 -
          8280 * t * s_p_2 - 22050 * s_p_3 + 15705 * t * s_p_3 + 37800 * s_p_4 -
          15552 * t * s_p_4 - 31752 * s_p_5 + 6804 * t * s_p_5 + 10368 * s_p_6 +
          975 * t_p_2 - 17010 * s * t_p_2 + 57726 * s_p_2 * t_p_2 -
          51030 * s_p_3 * t_p_2 + 8748 * s_p_4 * t_p_2 - 2880 * t_p_3 +
          57726 * s * t_p_3 - 171072 * s_p_2 * t_p_3 + 103680 * s_p_3 * t_p_3 +
          405 * t_p_4 - 51030 * s * t_p_4 + 103680 * s_p_2 * t_p_4 +
          3888 * t_p_5 + 8748 * s * t_p_5 - 2268 * t_p_6)) /
        20.;

      d2phi(10, 2) =
        (s_p_2 *
         (-50 + 599 * s + 975 * t - 11340 * s * t - 2070 * s_p_2 +
          28863 * t * s_p_2 + 3141 * s_p_3 - 20412 * t * s_p_3 - 2592 * s_p_4 +
          2916 * t * s_p_4 + 972 * s_p_5 - 4320 * t_p_2 + 57726 * s * t_p_2 -
          128304 * s_p_2 * t_p_2 + 62208 * s_p_3 * t_p_2 + 810 * t_p_3 -
          68040 * s * t_p_3 + 103680 * s_p_2 * t_p_3 + 9720 * t_p_4 +
          14580 * s * t_p_4 - 6804 * t_p_5)) /
        20.;

      d2phi(11, 0) =
        (t_p_2 *
         (10 - 153 * s - 325 * t + 11340 * s * t - 3780 * s_p_2 -
          57726 * t * s_p_2 + 29790 * s_p_3 + 68040 * t * s_p_3 -
          58320 * s_p_4 - 14580 * t * s_p_4 + 34020 * s_p_5 + 720 * t_p_2 -
          28863 * s * t_p_2 + 128304 * s_p_2 * t_p_2 - 103680 * s_p_3 * t_p_2 -
          81 * t_p_3 + 20412 * s * t_p_3 - 62208 * s_p_2 * t_p_3 - 648 * t_p_4 -
          2916 * s * t_p_4 + 324 * t_p_5)) /
        40.;

      d2phi(11, 1) =
        (s * t *
         (20 - 153 * s - 975 * t + 17010 * s * t - 2520 * s_p_2 -
          57726 * t * s_p_2 + 14895 * s_p_3 + 51030 * t * s_p_3 -
          23328 * s_p_4 - 8748 * t * s_p_4 + 11340 * s_p_5 + 2880 * t_p_2 -
          57726 * s * t_p_2 + 171072 * s_p_2 * t_p_2 - 103680 * s_p_3 * t_p_2 -
          405 * t_p_3 + 51030 * s * t_p_3 - 103680 * s_p_2 * t_p_3 -
          3888 * t_p_4 - 8748 * s * t_p_4 + 2268 * t_p_5)) /
        40.;

      d2phi(11, 2) =
        (s_p_2 *
         (10 - 51 * s - 975 * t + 11340 * s * t - 630 * s_p_2 -
          28863 * t * s_p_2 + 2979 * s_p_3 + 20412 * t * s_p_3 - 3888 * s_p_4 -
          2916 * t * s_p_4 + 1620 * s_p_5 + 4320 * t_p_2 - 57726 * s * t_p_2 +
          128304 * s_p_2 * t_p_2 - 62208 * s_p_3 * t_p_2 - 810 * t_p_3 +
          68040 * s * t_p_3 - 103680 * s_p_2 * t_p_3 - 9720 * t_p_4 -
          14580 * s * t_p_4 + 6804 * t_p_5)) /
        40.;


      d2phi(12, 0) =
        (t_p_2 *
         (10 - 975 * s - 51 * t + 11340 * s * t + 4320 * s_p_2 -
          57726 * t * s_p_2 - 810 * s_p_3 + 68040 * t * s_p_3 - 9720 * s_p_4 -
          14580 * t * s_p_4 + 6804 * s_p_5 - 630 * t_p_2 - 28863 * s * t_p_2 +
          128304 * s_p_2 * t_p_2 - 103680 * s_p_3 * t_p_2 + 2979 * t_p_3 +
          20412 * s * t_p_3 - 62208 * s_p_2 * t_p_3 - 3888 * t_p_4 -
          2916 * s * t_p_4 + 1620 * t_p_5)) /
        40.;

      d2phi(12, 1) =
        (s * t *
         (20 - 975 * s - 153 * t + 17010 * s * t + 2880 * s_p_2 -
          57726 * t * s_p_2 - 405 * s_p_3 + 51030 * t * s_p_3 - 3888 * s_p_4 -
          8748 * t * s_p_4 + 2268 * s_p_5 - 2520 * t_p_2 - 57726 * s * t_p_2 +
          171072 * s_p_2 * t_p_2 - 103680 * s_p_3 * t_p_2 + 14895 * t_p_3 +
          51030 * s * t_p_3 - 103680 * s_p_2 * t_p_3 - 23328 * t_p_4 -
          8748 * s * t_p_4 + 11340 * t_p_5)) /
        40.;

      d2phi(12, 2) =
        (s_p_2 *
         (10 - 325 * s - 153 * t + 11340 * s * t + 720 * s_p_2 -
          28863 * t * s_p_2 - 81 * s_p_3 + 20412 * t * s_p_3 - 648 * s_p_4 -
          2916 * t * s_p_4 + 324 * s_p_5 - 3780 * t_p_2 - 57726 * s * t_p_2 +
          128304 * s_p_2 * t_p_2 - 62208 * s_p_3 * t_p_2 + 29790 * t_p_3 +
          68040 * s * t_p_3 - 103680 * s_p_2 * t_p_3 - 58320 * t_p_4 -
          14580 * s * t_p_4 + 34020 * t_p_5)) /
        40.;

      d2phi(13, 0) =
        (t_p_2 *
         (-50 + 975 * s + 599 * t - 11340 * s * t - 4320 * s_p_2 +
          57726 * t * s_p_2 + 810 * s_p_3 - 68040 * t * s_p_3 + 9720 * s_p_4 +
          14580 * t * s_p_4 - 6804 * s_p_5 - 2070 * t_p_2 + 28863 * s * t_p_2 -
          128304 * s_p_2 * t_p_2 + 103680 * s_p_3 * t_p_2 + 3141 * t_p_3 -
          20412 * s * t_p_3 + 62208 * s_p_2 * t_p_3 - 2592 * t_p_4 +
          2916 * s * t_p_4 + 972 * t_p_5)) /
        20.;


      d2phi(13, 1) =
        (t *
         (40 - 100 * s - 882 * t + 1797 * s * t + 975 * s_p_2 -
          17010 * t * s_p_2 - 2880 * s_p_3 + 57726 * t * s_p_3 + 405 * s_p_4 -
          51030 * t * s_p_4 + 3888 * s_p_5 + 8748 * t * s_p_5 - 2268 * s_p_6 +
          6496 * t_p_2 - 8280 * s * t_p_2 + 57726 * s_p_2 * t_p_2 -
          171072 * s_p_3 * t_p_2 + 103680 * s_p_4 * t_p_2 - 22050 * t_p_3 +
          15705 * s * t_p_3 - 51030 * s_p_2 * t_p_3 + 103680 * s_p_3 * t_p_3 +
          37800 * t_p_4 - 15552 * s * t_p_4 + 8748 * s_p_2 * t_p_4 -
          31752 * t_p_5 + 6804 * s * t_p_5 + 10368 * t_p_6)) /
        20.;

      d2phi(13, 2) =
        -(s *
          (-40 + 50 * s + 1764 * t - 1797 * s * t - 325 * s_p_2 +
           11340 * t * s_p_2 + 720 * s_p_3 - 28863 * t * s_p_3 - 81 * s_p_4 +
           20412 * t * s_p_4 - 648 * s_p_5 - 2916 * t * s_p_5 + 324 * s_p_6 -
           19488 * t_p_2 + 12420 * s * t_p_2 - 57726 * s_p_2 * t_p_2 +
           128304 * s_p_3 * t_p_2 - 62208 * s_p_4 * t_p_2 + 88200 * t_p_3 -
           31410 * s * t_p_3 + 68040 * s_p_2 * t_p_3 - 103680 * s_p_3 * t_p_3 -
           189000 * t_p_4 + 38880 * s * t_p_4 - 14580 * s_p_2 * t_p_4 +
           190512 * t_p_5 - 20412 * s * t_p_5 - 72576 * t_p_6)) /
        20.;

      d2phi(14, 0) =
        (t_p_2 *
         (-628 + 4017 * s + 7783 * t - 44712 * s * t - 4644 * s_p_2 +
          32886 * t * s_p_2 + 3870 * s_p_3 + 25920 * t * s_p_3 - 9720 * s_p_4 -
          14580 * t * s_p_4 + 6804 * s_p_5 - 33192 * t_p_2 +
          167265 * s * t_p_2 - 122472 * s_p_2 * t_p_2 - 19440 * s_p_3 * t_p_2 +
          61677 * t_p_3 - 235224 * s * t_p_3 + 116640 * s_p_2 * t_p_3 -
          51516 * t_p_4 + 107892 * s * t_p_4 + 15876 * t_p_5)) /
        40.;

      d2phi(14, 1) =
        (s * t *
         (-1256 + 4017 * s + 23349 * t - 67068 * s * t - 3096 * s_p_2 +
          32886 * t * s_p_2 + 1935 * s_p_3 + 19440 * t * s_p_3 - 3888 * s_p_4 -
          8748 * t * s_p_4 + 2268 * s_p_5 - 132768 * t_p_2 +
          334530 * s * t_p_2 - 163296 * s_p_2 * t_p_2 - 19440 * s_p_3 * t_p_2 +
          308385 * t_p_3 - 588060 * s * t_p_3 + 194400 * s_p_2 * t_p_3 -
          309096 * t_p_4 + 323676 * s * t_p_4 + 111132 * t_p_5)) /
        40.;

      d2phi(14, 2) =
        (120 * t - 628 * s_p_2 + 23349 * t * s_p_2 + 1339 * s_p_3 -
         44712 * t * s_p_3 - 774 * s_p_4 + 16443 * t * s_p_4 + 387 * s_p_5 +
         7776 * t * s_p_5 - 648 * s_p_6 - 2916 * t * s_p_6 + 324 * s_p_7 -
         3288 * t_p_2 - 199152 * s_p_2 * t_p_2 + 334530 * s_p_3 * t_p_2 -
         122472 * s_p_4 * t_p_2 - 11664 * s_p_5 * t_p_2 + 27400 * t_p_3 +
         616770 * s_p_2 * t_p_3 - 784080 * s_p_3 * t_p_3 +
         194400 * s_p_4 * t_p_3 - 98820 * t_p_4 - 772740 * s_p_2 * t_p_4 +
         539460 * s_p_3 * t_p_4 + 173124 * t_p_5 + 333396 * s_p_2 * t_p_5 -
         145152 * t_p_6 + 46656 * t_p_7) /
        40.;

      d2phi(15, 0) =
        (20 - 882 * s + 9864 * s_p_2 - 46640 * s_p_3 + 111240 * s_p_4 -
         140994 * s_p_5 + 90720 * s_p_6 - 23328 * s_p_7 - 339 * t_p_2 +
         12573 * s * t_p_2 - 105786 * s_p_2 * t_p_2 + 324090 * s_p_3 * t_p_2 -
         405810 * s_p_4 * t_p_2 + 176904 * s_p_5 * t_p_2 + 832 * t_p_3 -
         28026 * s * t_p_3 + 196128 * s_p_2 * t_p_3 - 426060 * s_p_3 * t_p_3 +
         277020 * s_p_4 * t_p_3 - 747 * t_p_4 + 22653 * s * t_p_4 -
         125388 * s_p_2 * t_p_4 + 149040 * s_p_3 * t_p_4 + 234 * t_p_5 -
         6318 * s * t_p_5 + 25272 * s_p_2 * t_p_5) /
        20.;

      d2phi(15, 1) =
        (3 * s * t *
         (-226 + 4191 * s + 832 * t - 14013 * s * t - 23508 * s_p_2 +
          65376 * t * s_p_2 + 54015 * s_p_3 - 106515 * t * s_p_3 -
          54108 * s_p_4 + 55404 * t * s_p_4 + 19656 * s_p_5 - 996 * t_p_2 +
          15102 * s * t_p_2 - 55728 * s_p_2 * t_p_2 + 49680 * s_p_3 * t_p_2 +
          390 * t_p_3 - 5265 * s * t_p_3 + 14040 * s_p_2 * t_p_3)) /
        20.;

      d2phi(15, 2) = (3 * (-1 + 3 * s) * (-1 + 6 * s) * s_p_2 *
                      (-113 + 380 * s + 832 * t - 1854 * s * t - 423 * s_p_2 +
                       1026 * t * s_p_2 + 156 * s_p_3 - 1494 * t_p_2 +
                       1656 * s * t_p_2 + 780 * t_p_3)) /
                     20.;

      d2phi(16, 0) =
        -(t *
          (157 - 2877 * s - 1918 * t + 27153 * s * t + 18102 * s_p_2 -
           128520 * t * s_p_2 - 53550 * s_p_3 + 270900 * t * s_p_3 +
           81270 * s_p_4 - 262440 * t * s_p_4 - 61236 * s_p_5 +
           95256 * t * s_p_5 + 18144 * s_p_6 + 9051 * t_p_2 -
           96390 * s * t_p_2 + 325080 * s_p_2 * t_p_2 - 437400 * s_p_3 * t_p_2 +
           204120 * s_p_4 * t_p_2 - 21420 * t_p_3 + 162540 * s * t_p_3 -
           349920 * s_p_2 * t_p_3 + 226800 * s_p_3 * t_p_3 + 27090 * t_p_4 -
           131220 * s * t_p_4 + 136080 * s_p_2 * t_p_4 - 17496 * t_p_5 +
           40824 * s * t_p_5 + 4536 * t_p_6)) /
        5.;

      d2phi(16, 1) =
        (10 - 314 * s - 314 * t + 7672 * s * t + 2877 * s_p_2 -
         54306 * t * s_p_2 - 12068 * s_p_3 + 171360 * t * s_p_3 +
         26775 * s_p_4 - 270900 * t * s_p_4 - 32508 * s_p_5 +
         209952 * t * s_p_5 + 20412 * s_p_6 - 63504 * t * s_p_6 - 5184 * s_p_7 +
         2877 * t_p_2 - 54306 * s * t_p_2 + 289170 * s_p_2 * t_p_2 -
         650160 * s_p_3 * t_p_2 + 656100 * s_p_4 * t_p_2 -
         244944 * s_p_5 * t_p_2 - 12068 * t_p_3 + 171360 * s * t_p_3 -
         650160 * s_p_2 * t_p_3 + 933120 * s_p_3 * t_p_3 -
         453600 * s_p_4 * t_p_3 + 26775 * t_p_4 - 270900 * s * t_p_4 +
         656100 * s_p_2 * t_p_4 - 453600 * s_p_3 * t_p_4 - 32508 * t_p_5 +
         209952 * s * t_p_5 - 244944 * s_p_2 * t_p_5 + 20412 * t_p_6 -
         63504 * s * t_p_6 - 5184 * t_p_7) /
        10.;

      d2phi(16, 2) =
        -(s * (157 - 1918 * s - 2877 * t + 27153 * s * t + 9051 * s_p_2 -
               96390 * t * s_p_2 - 21420 * s_p_3 + 162540 * t * s_p_3 +
               27090 * s_p_4 - 131220 * t * s_p_4 - 17496 * s_p_5 +
               40824 * t * s_p_5 + 4536 * s_p_6 + 18102 * t_p_2 -
               128520 * s * t_p_2 + 325080 * s_p_2 * t_p_2 -
               349920 * s_p_3 * t_p_2 + 136080 * s_p_4 * t_p_2 - 53550 * t_p_3 +
               270900 * s * t_p_3 - 437400 * s_p_2 * t_p_3 +
               226800 * s_p_3 * t_p_3 + 81270 * t_p_4 - 262440 * s * t_p_4 +
               204120 * s_p_2 * t_p_4 - 61236 * t_p_5 + 95256 * s * t_p_5 +
               18144 * t_p_6)) /
        5.;

      d2phi(17, 0) = (3 * (-1 + 3 * t) * (-1 + 6 * t) * t_p_2 *
                      (-113 + 832 * s + 380 * t - 1854 * s * t - 1494 * s_p_2 +
                       1656 * t * s_p_2 + 780 * s_p_3 - 423 * t_p_2 +
                       1026 * s * t_p_2 + 156 * t_p_3)) /
                     20.;

      d2phi(17, 1) =
        (3 * s * t *
         (-226 + 832 * s + 4191 * t - 14013 * s * t - 996 * s_p_2 +
          15102 * t * s_p_2 + 390 * s_p_3 - 5265 * t * s_p_3 - 23508 * t_p_2 +
          65376 * s * t_p_2 - 55728 * s_p_2 * t_p_2 + 14040 * s_p_3 * t_p_2 +
          54015 * t_p_3 - 106515 * s * t_p_3 + 49680 * s_p_2 * t_p_3 -
          54108 * t_p_4 + 55404 * s * t_p_4 + 19656 * t_p_5)) /
        20.;

      d2phi(17, 2) =
        (20 - 882 * t - 339 * s_p_2 + 12573 * t * s_p_2 + 832 * s_p_3 -
         28026 * t * s_p_3 - 747 * s_p_4 + 22653 * t * s_p_4 + 234 * s_p_5 -
         6318 * t * s_p_5 + 9864 * t_p_2 - 105786 * s_p_2 * t_p_2 +
         196128 * s_p_3 * t_p_2 - 125388 * s_p_4 * t_p_2 +
         25272 * s_p_5 * t_p_2 - 46640 * t_p_3 + 324090 * s_p_2 * t_p_3 -
         426060 * s_p_3 * t_p_3 + 149040 * s_p_4 * t_p_3 + 111240 * t_p_4 -
         405810 * s_p_2 * t_p_4 + 277020 * s_p_3 * t_p_4 - 140994 * t_p_5 +
         176904 * s_p_2 * t_p_5 + 90720 * t_p_6 - 23328 * t_p_7) /
        20.;

      d2phi(18, 0) =
        -32 * (-1 + 3 * t) * (-1 + 6 * t) * t_p_2 *
        (-47 + 246 * s + 164 * t - 567 * s * t - 378 * s_p_2 + 432 * t * s_p_2 +
         180 * s_p_3 - 189 * t_p_2 + 324 * s * t_p_2 + 72 * t_p_3);

      d2phi(18, 1) =
        -16 * t *
        (20 - 188 * s - 411 * t + 3522 * s * t + 492 * s_p_2 -
         8343 * t * s_p_2 - 504 * s_p_3 + 7668 * t * s_p_3 + 180 * s_p_4 -
         2430 * t * s_p_4 + 2740 * t_p_2 - 20088 * s * t_p_2 +
         39420 * s_p_2 * t_p_2 - 28512 * s_p_3 * t_p_2 + 6480 * s_p_4 * t_p_2 -
         8235 * t_p_3 + 47250 * s * t_p_3 - 65610 * s_p_2 * t_p_3 +
         25920 * s_p_3 * t_p_3 + 12366 * t_p_4 - 48600 * s * t_p_4 +
         34992 * s_p_2 * t_p_4 - 9072 * t_p_5 + 18144 * s * t_p_5 +
         2592 * t_p_6);

      d2phi(18, 2) =
        -32 * s *
        (10 - 47 * s - 411 * t + 1761 * s * t + 82 * s_p_2 - 2781 * t * s_p_2 -
         63 * s_p_3 + 1917 * t * s_p_3 + 18 * s_p_4 - 486 * t * s_p_4 +
         4110 * t_p_2 - 15066 * s * t_p_2 + 19710 * s_p_2 * t_p_2 -
         10692 * s_p_3 * t_p_2 + 1944 * s_p_4 * t_p_2 - 16470 * t_p_3 +
         47250 * s * t_p_3 - 43740 * s_p_2 * t_p_3 + 12960 * s_p_3 * t_p_3 +
         30915 * t_p_4 - 60750 * s * t_p_4 + 29160 * s_p_2 * t_p_4 -
         27216 * t_p_5 + 27216 * s * t_p_5 + 9072 * t_p_6);

      d2phi(19, 0) =
        -32 * t *
        (10 - 411 * s - 47 * t + 1761 * s * t + 4110 * s_p_2 -
         15066 * t * s_p_2 - 16470 * s_p_3 + 47250 * t * s_p_3 + 30915 * s_p_4 -
         60750 * t * s_p_4 - 27216 * s_p_5 + 27216 * t * s_p_5 + 9072 * s_p_6 +
         82 * t_p_2 - 2781 * s * t_p_2 + 19710 * s_p_2 * t_p_2 -
         43740 * s_p_3 * t_p_2 + 29160 * s_p_4 * t_p_2 - 63 * t_p_3 +
         1917 * s * t_p_3 - 10692 * s_p_2 * t_p_3 + 12960 * s_p_3 * t_p_3 +
         18 * t_p_4 - 486 * s * t_p_4 + 1944 * s_p_2 * t_p_4);


      d2phi(19, 1) =
        -16 * s *
        (20 - 411 * s - 188 * t + 3522 * s * t + 2740 * s_p_2 -
         20088 * t * s_p_2 - 8235 * s_p_3 + 47250 * t * s_p_3 + 12366 * s_p_4 -
         48600 * t * s_p_4 - 9072 * s_p_5 + 18144 * t * s_p_5 + 2592 * s_p_6 +
         492 * t_p_2 - 8343 * s * t_p_2 + 39420 * s_p_2 * t_p_2 -
         65610 * s_p_3 * t_p_2 + 34992 * s_p_4 * t_p_2 - 504 * t_p_3 +
         7668 * s * t_p_3 - 28512 * s_p_2 * t_p_3 + 25920 * s_p_3 * t_p_3 +
         180 * t_p_4 - 2430 * s * t_p_4 + 6480 * s_p_2 * t_p_4);

      d2phi(19, 2) =
        -32 * (-1 + 3 * s) * (-1 + 6 * s) * s_p_2 *
        (-47 + 164 * s + 246 * t - 567 * s * t - 189 * s_p_2 + 324 * t * s_p_2 +
         72 * s_p_3 - 378 * t_p_2 + 432 * s * t_p_2 + 180 * t_p_3);

      d2phi(20, 0) = 16 * (-1 + 3 * t) * (-1 + 6 * t) *
                     (-1 + 30 * s + t - 27 * s * t - 162 * s_p_2 +
                      108 * t * s_p_2 + 180 * s_p_3) *
                     t_p_2 * sqrt(2);

      d2phi(20, 1) =
        8 * s * t *
        (-4 + 60 * s + 60 * t - 891 * s * t - 216 * s_p_2 + 3132 * t * s_p_2 +
         180 * s_p_3 - 2430 * t * s_p_3 - 216 * t_p_2 + 3132 * s * t_p_2 -
         10368 * s_p_2 * t_p_2 + 6480 * s_p_3 * t_p_2 + 180 * t_p_3 -
         2430 * s * t_p_3 + 6480 * s_p_2 * t_p_3) *
        sqrt(2);

      d2phi(20, 2) = 16 * (-1 + 3 * s) * (-1 + 6 * s) * s_p_2 *
                     (-1 + s + 30 * t - 27 * s * t - 162 * t_p_2 +
                      108 * s * t_p_2 + 180 * t_p_3) *
                     sqrt(2);

      d2phi(21, 0) = (-11664 * (-1 + 3 * t) * (-1 + 6 * t) * t_p_2 *
                      (-93 + 792 * s + 290 * t - 1714 * s * t - 1494 * s_p_2 +
                       1656 * t * s_p_2 + 780 * s_p_3 - 293 * t_p_2 +
                       906 * s * t_p_2 + 96 * t_p_3)) /
                     125.;

      d2phi(21, 1) =
        (-11664 * s * t *
         (-186 + 792 * s + 3381 * t - 13263 * s * t - 996 * s_p_2 +
          15102 * t * s_p_2 + 390 * s_p_3 - 5265 * t * s_p_3 - 18308 * t_p_2 +
          61176 * s * t_p_2 - 55728 * s_p_2 * t_p_2 + 14040 * s_p_3 * t_p_2 +
          39765 * t_p_3 - 97515 * s * t_p_3 + 49680 * s_p_2 * t_p_3 -
          36828 * t_p_4 + 48924 * s * t_p_4 + 12096 * t_p_5)) /
        125.;

      d2phi(21, 2) =
        (-11664 *
         (16 * t - 93 * s_p_2 + 3381 * t * s_p_2 + 264 * s_p_3 -
          8842 * t * s_p_3 - 249 * s_p_4 + 7551 * t * s_p_4 + 78 * s_p_5 -
          2106 * t * s_p_5 - 432 * t_p_2 - 27462 * s_p_2 * t_p_2 +
          61176 * s_p_3 * t_p_2 - 41796 * s_p_4 * t_p_2 + 8424 * s_p_5 * t_p_2 +
          3520 * t_p_3 + 79530 * s_p_2 * t_p_3 - 130020 * s_p_3 * t_p_3 +
          49680 * s_p_4 * t_p_3 - 12320 * t_p_4 - 92070 * s_p_2 * t_p_4 +
          81540 * s_p_3 * t_p_4 + 20832 * t_p_5 + 36288 * s_p_2 * t_p_5 -
          16800 * t_p_6 + 5184 * t_p_7)) /
        125.;

      d2phi(22, 0) = (729 * (-1 + 3 * t) * (-1 + 6 * t) * t_p_2 *
                      (-61 + 600 * s + 168 * t - 1150 * s * t - 1302 * s_p_2 +
                       1272 * t * s_p_2 + 780 * s_p_3 - 151 * t_p_2 +
                       546 * s * t_p_2 + 44 * t_p_3)) /
                     8.;

      d2phi(22, 1) =
        (729 * s * t *
         (-122 + 600 * s + 2151 * t - 9825 * s * t - 868 * s_p_2 +
          12990 * t * s_p_2 + 390 * s_p_3 - 5265 * t * s_p_3 - 11044 * t_p_2 +
          43392 * s * t_p_2 - 46512 * s_p_2 * t_p_2 + 14040 * s_p_3 * t_p_2 +
          22135 * t_p_3 - 64035 * s * t_p_3 + 38160 * s_p_2 * t_p_3 -
          18684 * t_p_4 + 29484 * s * t_p_4 + 5544 * t_p_5)) /
        8.;

      d2phi(22, 2) =
        (729 * (10 * t - 61 * s_p_2 + 2151 * t * s_p_2 + 200 * s_p_3 -
                6550 * t * s_p_3 - 217 * s_p_4 + 6495 * t * s_p_4 + 78 * s_p_5 -
                2106 * t * s_p_5 - 264 * t_p_2 - 16566 * s_p_2 * t_p_2 +
                43392 * s_p_3 * t_p_2 - 34884 * s_p_4 * t_p_2 +
                8424 * s_p_5 * t_p_2 + 2080 * t_p_3 + 44270 * s_p_2 * t_p_3 -
                85380 * s_p_3 * t_p_3 + 38160 * s_p_4 * t_p_3 - 6980 * t_p_4 -
                46710 * s_p_2 * t_p_4 + 49140 * s_p_3 * t_p_4 + 11298 * t_p_5 +
                16632 * s_p_2 * t_p_5 - 8736 * t_p_6 + 2592 * t_p_7)) /
        8.;

      d2phi(23, 0) = (-81 * (-1 + 6 * t) * t_p_2 *
                      (1003 - 14712 * s - 3635 * t + 43134 * s * t +
                       52890 * s_p_2 - 106758 * t * s_p_2 - 69300 * s_p_3 +
                       71100 * t * s_p_3 + 30240 * s_p_4 + 4603 * t_p_2 -
                       40932 * s * t_p_2 + 52920 * s_p_2 * t_p_2 -
                       2295 * t_p_3 + 12474 * s * t_p_3 + 324 * t_p_4)) /
                     8.;

      d2phi(23, 1) =
        (-81 * s * t *
         (-2006 + 14712 * s + 28959 * t - 197109 * s * t - 35260 * s_p_2 +
          424098 * t * s_p_2 + 34650 * s_p_3 - 365175 * t * s_p_3 -
          12096 * s_p_4 + 108864 * t * s_p_4 - 105652 * t_p_2 +
          599472 * s * t_p_2 - 924624 * s_p_2 * t_p_2 + 426600 * s_p_3 * t_p_2 +
          149565 * t_p_3 - 645165 * s * t_p_3 + 529200 * s_p_2 * t_p_3 -
          84564 * t_p_4 + 224532 * s * t_p_4 + 13608 * t_p_5)) /
        8.;

      d2phi(23, 2) =
        (-81 *
         (180 * t - 1003 * s_p_2 + 28959 * t * s_p_2 + 4904 * s_p_3 -
          131406 * t * s_p_3 - 8815 * s_p_4 + 212049 * t * s_p_4 +
          6930 * s_p_5 - 146070 * t * s_p_5 - 2016 * s_p_6 + 36288 * t * s_p_6 -
          4212 * t_p_2 - 158478 * s_p_2 * t_p_2 + 599472 * s_p_3 * t_p_2 -
          693468 * s_p_4 * t_p_2 + 255960 * s_p_5 * t_p_2 + 28260 * t_p_3 +
          299130 * s_p_2 * t_p_3 - 860220 * s_p_3 * t_p_3 +
          529200 * s_p_4 * t_p_3 - 82710 * t_p_4 - 211410 * s_p_2 * t_p_4 +
          374220 * s_p_3 * t_p_4 + 119826 * t_p_5 + 40824 * s_p_2 * t_p_5 -
          84672 * t_p_6 + 23328 * t_p_7)) /
        8.;

      d2phi(24, 0) =
        (-648 * t_p_2 *
         (5929 - 111993 * s - 32876 * t + 553884 * s * t + 553578 * s_p_2 -
          2146824 * t * s_p_2 - 1123410 * s_p_3 + 2966040 * t * s_p_3 +
          1013580 * s_p_4 - 1351080 * t * s_p_4 - 337932 * s_p_5 +
          58369 * t_p_2 - 910791 * s * t_p_2 + 2527632 * s_p_2 * t_p_2 -
          1804680 * s_p_3 * t_p_2 - 31746 * t_p_3 + 584172 * s * t_p_3 -
          914976 * s_p_2 * t_p_3 - 10044 * t_p_4 - 114696 * s * t_p_4 +
          10368 * t_p_5)) /
        125.;

      d2phi(24, 1) =
        (648 * s * t *
         (-11858 + 111993 * s + 98628 * t - 830826 * s * t - 369052 * s_p_2 +
          2146824 * t * s_p_2 + 561705 * s_p_3 - 2224530 * t * s_p_3 -
          405432 * s_p_4 + 810648 * t * s_p_4 + 112644 * s_p_5 -
          233476 * t_p_2 + 1821582 * s * t_p_2 - 3370176 * s_p_2 * t_p_2 +
          1804680 * s_p_3 * t_p_2 + 158730 * t_p_3 - 1460430 * s * t_p_3 +
          1524960 * s_p_2 * t_p_3 + 60264 * t_p_4 + 344088 * s * t_p_4 -
          72576 * t_p_5)) /
        125.;

      d2phi(24, 2) =
        (648 *
         (1440 * t - 5929 * s_p_2 + 98628 * t * s_p_2 + 37331 * s_p_3 -
          553884 * t * s_p_3 - 92263 * s_p_4 + 1073412 * t * s_p_4 +
          112341 * s_p_5 - 889812 * t * s_p_5 - 67572 * s_p_6 +
          270216 * t * s_p_6 + 16092 * s_p_7 - 25056 * t_p_2 -
          350214 * s_p_2 * t_p_2 + 1821582 * s_p_3 * t_p_2 -
          2527632 * s_p_4 * t_p_2 + 1082808 * s_p_5 * t_p_2 + 144000 * t_p_3 +
          317460 * s_p_2 * t_p_3 - 1947240 * s_p_3 * t_p_3 +
          1524960 * s_p_4 * t_p_3 - 383040 * t_p_4 + 150660 * s_p_2 * t_p_4 +
          573480 * s_p_3 * t_p_4 + 520128 * t_p_5 - 217728 * s_p_2 * t_p_5 -
          350784 * t_p_6 + 93312 * t_p_7)) /
        125.;

      d2phi(25, 0) =
        (648 * (1440 * s - 25056 * s_p_2 + 144000 * s_p_3 - 383040 * s_p_4 +
                520128 * s_p_5 - 350784 * s_p_6 + 93312 * s_p_7 - 5929 * t_p_2 +
                98628 * s * t_p_2 - 350214 * s_p_2 * t_p_2 +
                317460 * s_p_3 * t_p_2 + 150660 * s_p_4 * t_p_2 -
                217728 * s_p_5 * t_p_2 + 37331 * t_p_3 - 553884 * s * t_p_3 +
                1821582 * s_p_2 * t_p_3 - 1947240 * s_p_3 * t_p_3 +
                573480 * s_p_4 * t_p_3 - 92263 * t_p_4 + 1073412 * s * t_p_4 -
                2527632 * s_p_2 * t_p_4 + 1524960 * s_p_3 * t_p_4 +
                112341 * t_p_5 - 889812 * s * t_p_5 + 1082808 * s_p_2 * t_p_5 -
                67572 * t_p_6 + 270216 * s * t_p_6 + 16092 * t_p_7)) /
        125.;

      d2phi(25, 1) =
        (-648 * s * t *
         (11858 - 98628 * s - 111993 * t + 830826 * s * t + 233476 * s_p_2 -
          1821582 * t * s_p_2 - 158730 * s_p_3 + 1460430 * t * s_p_3 -
          60264 * s_p_4 - 344088 * t * s_p_4 + 72576 * s_p_5 + 369052 * t_p_2 -
          2146824 * s * t_p_2 + 3370176 * s_p_2 * t_p_2 -
          1524960 * s_p_3 * t_p_2 - 561705 * t_p_3 + 2224530 * s * t_p_3 -
          1804680 * s_p_2 * t_p_3 + 405432 * t_p_4 - 810648 * s * t_p_4 -
          112644 * t_p_5)) /
        125.;

      d2phi(25, 2) =
        (-648 * s_p_2 *
         (5929 - 32876 * s - 111993 * t + 553884 * s * t + 58369 * s_p_2 -
          910791 * t * s_p_2 - 31746 * s_p_3 + 584172 * t * s_p_3 -
          10044 * s_p_4 - 114696 * t * s_p_4 + 10368 * s_p_5 + 553578 * t_p_2 -
          2146824 * s * t_p_2 + 2527632 * s_p_2 * t_p_2 -
          914976 * s_p_3 * t_p_2 - 1123410 * t_p_3 + 2966040 * s * t_p_3 -
          1804680 * s_p_2 * t_p_3 + 1013580 * t_p_4 - 1351080 * s * t_p_4 -
          337932 * t_p_5)) /
        125.;

      d2phi(26, 0) =
        (-81 *
         (180 * s - 4212 * s_p_2 + 28260 * s_p_3 - 82710 * s_p_4 +
          119826 * s_p_5 - 84672 * s_p_6 + 23328 * s_p_7 - 1003 * t_p_2 +
          28959 * s * t_p_2 - 158478 * s_p_2 * t_p_2 + 299130 * s_p_3 * t_p_2 -
          211410 * s_p_4 * t_p_2 + 40824 * s_p_5 * t_p_2 + 4904 * t_p_3 -
          131406 * s * t_p_3 + 599472 * s_p_2 * t_p_3 - 860220 * s_p_3 * t_p_3 +
          374220 * s_p_4 * t_p_3 - 8815 * t_p_4 + 212049 * s * t_p_4 -
          693468 * s_p_2 * t_p_4 + 529200 * s_p_3 * t_p_4 + 6930 * t_p_5 -
          146070 * s * t_p_5 + 255960 * s_p_2 * t_p_5 - 2016 * t_p_6 +
          36288 * s * t_p_6)) /
        8.;

      d2phi(26, 1) =
        (-81 * s * t *
         (-2006 + 28959 * s + 14712 * t - 197109 * s * t - 105652 * s_p_2 +
          599472 * t * s_p_2 + 149565 * s_p_3 - 645165 * t * s_p_3 -
          84564 * s_p_4 + 224532 * t * s_p_4 + 13608 * s_p_5 - 35260 * t_p_2 +
          424098 * s * t_p_2 - 924624 * s_p_2 * t_p_2 + 529200 * s_p_3 * t_p_2 +
          34650 * t_p_3 - 365175 * s * t_p_3 + 426600 * s_p_2 * t_p_3 -
          12096 * t_p_4 + 108864 * s * t_p_4)) /
        8.;


      d2phi(26, 2) =
        (-81 * (-1 + 6 * s) * s_p_2 *
         (1003 - 3635 * s - 14712 * t + 43134 * s * t + 4603 * s_p_2 -
          40932 * t * s_p_2 - 2295 * s_p_3 + 12474 * t * s_p_3 + 324 * s_p_4 +
          52890 * t_p_2 - 106758 * s * t_p_2 + 52920 * s_p_2 * t_p_2 -
          69300 * t_p_3 + 71100 * s * t_p_3 + 30240 * t_p_4)) /
        8.;

      d2phi(27, 0) =
        (729 *
         (10 * s - 264 * s_p_2 + 2080 * s_p_3 - 6980 * s_p_4 + 11298 * s_p_5 -
          8736 * s_p_6 + 2592 * s_p_7 - 61 * t_p_2 + 2151 * s * t_p_2 -
          16566 * s_p_2 * t_p_2 + 44270 * s_p_3 * t_p_2 -
          46710 * s_p_4 * t_p_2 + 16632 * s_p_5 * t_p_2 + 200 * t_p_3 -
          6550 * s * t_p_3 + 43392 * s_p_2 * t_p_3 - 85380 * s_p_3 * t_p_3 +
          49140 * s_p_4 * t_p_3 - 217 * t_p_4 + 6495 * s * t_p_4 -
          34884 * s_p_2 * t_p_4 + 38160 * s_p_3 * t_p_4 + 78 * t_p_5 -
          2106 * s * t_p_5 + 8424 * s_p_2 * t_p_5)) /
        8.;


      d2phi(27, 1) =
        (729 * s * t *
         (-122 + 2151 * s + 600 * t - 9825 * s * t - 11044 * s_p_2 +
          43392 * t * s_p_2 + 22135 * s_p_3 - 64035 * t * s_p_3 -
          18684 * s_p_4 + 29484 * t * s_p_4 + 5544 * s_p_5 - 868 * t_p_2 +
          12990 * s * t_p_2 - 46512 * s_p_2 * t_p_2 + 38160 * s_p_3 * t_p_2 +
          390 * t_p_3 - 5265 * s * t_p_3 + 14040 * s_p_2 * t_p_3)) /
        8.;

      d2phi(27, 2) = (729 * (-1 + 3 * s) * (-1 + 6 * s) * s_p_2 *
                      (-61 + 168 * s + 600 * t - 1150 * s * t - 151 * s_p_2 +
                       546 * t * s_p_2 + 44 * s_p_3 - 1302 * t_p_2 +
                       1272 * s * t_p_2 + 780 * t_p_3)) /
                     8.;

      d2phi(28, 0) =
        (-11664 *
         (16 * s - 432 * s_p_2 + 3520 * s_p_3 - 12320 * s_p_4 + 20832 * s_p_5 -
          16800 * s_p_6 + 5184 * s_p_7 - 93 * t_p_2 + 3381 * s * t_p_2 -
          27462 * s_p_2 * t_p_2 + 79530 * s_p_3 * t_p_2 -
          92070 * s_p_4 * t_p_2 + 36288 * s_p_5 * t_p_2 + 264 * t_p_3 -
          8842 * s * t_p_3 + 61176 * s_p_2 * t_p_3 - 130020 * s_p_3 * t_p_3 +
          81540 * s_p_4 * t_p_3 - 249 * t_p_4 + 7551 * s * t_p_4 -
          41796 * s_p_2 * t_p_4 + 49680 * s_p_3 * t_p_4 + 78 * t_p_5 -
          2106 * s * t_p_5 + 8424 * s_p_2 * t_p_5)) /
        125.;

      d2phi(28, 1) =
        (-11664 * s * t *
         (-186 + 3381 * s + 792 * t - 13263 * s * t - 18308 * s_p_2 +
          61176 * t * s_p_2 + 39765 * s_p_3 - 97515 * t * s_p_3 -
          36828 * s_p_4 + 48924 * t * s_p_4 + 12096 * s_p_5 - 996 * t_p_2 +
          15102 * s * t_p_2 - 55728 * s_p_2 * t_p_2 + 49680 * s_p_3 * t_p_2 +
          390 * t_p_3 - 5265 * s * t_p_3 + 14040 * s_p_2 * t_p_3)) /
        125.;

      d2phi(28, 2) = (-11664 * (-1 + 3 * s) * (-1 + 6 * s) * s_p_2 *
                      (-93 + 290 * s + 792 * t - 1714 * s * t - 293 * s_p_2 +
                       906 * t * s_p_2 + 96 * s_p_3 - 1494 * t_p_2 +
                       1656 * s * t_p_2 + 780 * t_p_3)) /
                     125.;

      d2phi(29, 0) =
        (-324 * t_p_2 *
         (-518 + 22449 * s - 1313 * t + 42000 * s * t - 244200 * s_p_2 -
          131346 * t * s_p_2 + 970590 * s_p_3 - 200160 * t * s_p_3 -
          1512000 * s_p_4 + 445500 * t * s_p_4 + 786996 * s_p_5 + 4540 * t_p_2 -
          170253 * s * t_p_2 + 756864 * s_p_2 * t_p_2 - 589680 * s_p_3 * t_p_2 -
          1161 * t_p_3 + 121392 * s * t_p_3 - 379728 * s_p_2 * t_p_3 -
          3168 * t_p_4 - 15876 * s * t_p_4 + 1620 * t_p_5)) /
        125.;

      d2phi(29, 1) =
        (-324 * s * t *
         (-1036 + 22449 * s - 3939 * t + 63000 * s * t - 162800 * s_p_2 -
          131346 * t * s_p_2 + 485295 * s_p_3 - 150120 * t * s_p_3 -
          604800 * s_p_4 + 267300 * t * s_p_4 + 262332 * s_p_5 + 18160 * t_p_2 -
          340506 * s * t_p_2 + 1009152 * s_p_2 * t_p_2 -
          589680 * s_p_3 * t_p_2 - 5805 * t_p_3 + 303480 * s * t_p_3 -
          632880 * s_p_2 * t_p_3 - 19008 * t_p_4 - 47628 * s * t_p_4 +
          11340 * t_p_5)) /
        125.;

      d2phi(29, 2) =
        (-324 * s_p_2 *
         (-518 + 7483 * s - 3939 * t + 42000 * s * t - 40700 * s_p_2 -
          65673 * t * s_p_2 + 97059 * s_p_3 - 60048 * t * s_p_3 -
          100800 * s_p_4 + 89100 * t * s_p_4 + 37476 * s_p_5 + 27240 * t_p_2 -
          340506 * s * t_p_2 + 756864 * s_p_2 * t_p_2 - 353808 * s_p_3 * t_p_2 -
          11610 * t_p_3 + 404640 * s * t_p_3 - 632880 * s_p_2 * t_p_3 -
          47520 * t_p_4 - 79380 * s * t_p_4 + 34020 * t_p_5)) /
        125.;

      d2phi(30, 0) =
        (81 * t_p_2 *
         (274 - 9501 * s - 2021 * t + 68064 * s * t + 67416 * s_p_2 -
          432522 * t * s_p_2 - 158670 * s_p_3 + 820440 * t * s_p_3 +
          137160 * s_p_4 - 458460 * t * s_p_4 - 34020 * s_p_5 + 3382 * t_p_2 -
          113301 * s * t_p_2 + 566352 * s_p_2 * t_p_2 - 587520 * s_p_3 * t_p_2 -
          1167 * t_p_3 + 60804 * s * t_p_3 - 202176 * s_p_2 * t_p_3 -
          1008 * t_p_4 - 6156 * s * t_p_4 + 540 * t_p_5)) /
        16.;

      d2phi(30, 1) =
        (-81 * s * t *
         (-548 + 9501 * s + 6063 * t - 102096 * s * t - 44944 * s_p_2 +
          432522 * t * s_p_2 + 79335 * s_p_3 - 615330 * t * s_p_3 -
          54864 * s_p_4 + 275076 * t * s_p_4 + 11340 * s_p_5 - 13528 * t_p_2 +
          226602 * s * t_p_2 - 755136 * s_p_2 * t_p_2 + 587520 * s_p_3 * t_p_2 +
          5835 * t_p_3 - 152010 * s * t_p_3 + 336960 * s_p_2 * t_p_3 +
          6048 * t_p_4 + 18468 * s * t_p_4 - 3780 * t_p_5)) /
        16.;


      d2phi(30, 2) =
        (-81 * s_p_2 *
         (-274 + 3167 * s + 6063 * t - 68064 * s * t - 11236 * s_p_2 +
          216261 * t * s_p_2 + 15867 * s_p_3 - 246132 * t * s_p_3 -
          9144 * s_p_4 + 91692 * t * s_p_4 + 1620 * s_p_5 - 20292 * t_p_2 +
          226602 * s * t_p_2 - 566352 * s_p_2 * t_p_2 + 352512 * s_p_3 * t_p_2 +
          11670 * t_p_3 - 202680 * s * t_p_3 + 336960 * s_p_2 * t_p_3 +
          15120 * t_p_4 + 30780 * s * t_p_4 - 11340 * t_p_5)) /
        16.;

      d2phi(31, 0) =
        (-81 * t_p_2 *
         (-274 + 6063 * s + 3167 * t - 68064 * s * t - 20292 * s_p_2 +
          226602 * t * s_p_2 + 11670 * s_p_3 - 202680 * t * s_p_3 +
          15120 * s_p_4 + 30780 * t * s_p_4 - 11340 * s_p_5 - 11236 * t_p_2 +
          216261 * s * t_p_2 - 566352 * s_p_2 * t_p_2 + 336960 * s_p_3 * t_p_2 +
          15867 * t_p_3 - 246132 * s * t_p_3 + 352512 * s_p_2 * t_p_3 -
          9144 * t_p_4 + 91692 * s * t_p_4 + 1620 * t_p_5)) /
        16.;

      d2phi(31, 1) =
        (81 * s * t *
         (548 - 6063 * s - 9501 * t + 102096 * s * t + 13528 * s_p_2 -
          226602 * t * s_p_2 - 5835 * s_p_3 + 152010 * t * s_p_3 -
          6048 * s_p_4 - 18468 * t * s_p_4 + 3780 * s_p_5 + 44944 * t_p_2 -
          432522 * s * t_p_2 + 755136 * s_p_2 * t_p_2 - 336960 * s_p_3 * t_p_2 -
          79335 * t_p_3 + 615330 * s * t_p_3 - 587520 * s_p_2 * t_p_3 +
          54864 * t_p_4 - 275076 * s * t_p_4 - 11340 * t_p_5)) /
        16.;

      d2phi(31, 2) =
        (81 * s_p_2 *
         (274 - 2021 * s - 9501 * t + 68064 * s * t + 3382 * s_p_2 -
          113301 * t * s_p_2 - 1167 * s_p_3 + 60804 * t * s_p_3 - 1008 * s_p_4 -
          6156 * t * s_p_4 + 540 * s_p_5 + 67416 * t_p_2 - 432522 * s * t_p_2 +
          566352 * s_p_2 * t_p_2 - 202176 * s_p_3 * t_p_2 - 158670 * t_p_3 +
          820440 * s * t_p_3 - 587520 * s_p_2 * t_p_3 + 137160 * t_p_4 -
          458460 * s * t_p_4 - 34020 * t_p_5)) /
        16.;

      d2phi(32, 0) =
        (-324 * t_p_2 *
         (-518 - 3939 * s + 7483 * t + 42000 * s * t + 27240 * s_p_2 -
          340506 * t * s_p_2 - 11610 * s_p_3 + 404640 * t * s_p_3 -
          47520 * s_p_4 - 79380 * t * s_p_4 + 34020 * s_p_5 - 40700 * t_p_2 -
          65673 * s * t_p_2 + 756864 * s_p_2 * t_p_2 - 632880 * s_p_3 * t_p_2 +
          97059 * t_p_3 - 60048 * s * t_p_3 - 353808 * s_p_2 * t_p_3 -
          100800 * t_p_4 + 89100 * s * t_p_4 + 37476 * t_p_5)) /
        125.;

      d2phi(32, 1) =
        (-324 * s * t *
         (-1036 - 3939 * s + 22449 * t + 63000 * s * t + 18160 * s_p_2 -
          340506 * t * s_p_2 - 5805 * s_p_3 + 303480 * t * s_p_3 -
          19008 * s_p_4 - 47628 * t * s_p_4 + 11340 * s_p_5 - 162800 * t_p_2 -
          131346 * s * t_p_2 + 1009152 * s_p_2 * t_p_2 -
          632880 * s_p_3 * t_p_2 + 485295 * t_p_3 - 150120 * s * t_p_3 -
          589680 * s_p_2 * t_p_3 - 604800 * t_p_4 + 267300 * s * t_p_4 +
          262332 * t_p_5)) /
        125.;

      d2phi(32, 2) =
        (-324 * s_p_2 *
         (-518 - 1313 * s + 22449 * t + 42000 * s * t + 4540 * s_p_2 -
          170253 * t * s_p_2 - 1161 * s_p_3 + 121392 * t * s_p_3 -
          3168 * s_p_4 - 15876 * t * s_p_4 + 1620 * s_p_5 - 244200 * t_p_2 -
          131346 * s * t_p_2 + 756864 * s_p_2 * t_p_2 - 379728 * s_p_3 * t_p_2 +
          970590 * t_p_3 - 200160 * s * t_p_3 - 589680 * s_p_2 * t_p_3 -
          1512000 * t_p_4 + 445500 * s * t_p_4 + 786996 * t_p_5)) /
        125.;

      d2phi(33, 0) = (-1296 * (-1 + 2 * t) * (-2 + 3 * s + 2 * t) *
                      (-2 + 3 * t) * (-1 + 3 * t) * (-1 + 6 * t) * t_p_2) /
                     25.;

      d2phi(33, 1) =
        (-648 * t *
         (4 - 16 * s - 87 * t + 324 * s * t + 12 * s_p_2 - 225 * t * s_p_2 +
          628 * t_p_2 - 2080 * s * t_p_2 + 1260 * s_p_2 * t_p_2 - 2075 * t_p_3 +
          5700 * s * t_p_3 - 2700 * s_p_2 * t_p_3 + 3438 * t_p_4 -
          6912 * s * t_p_4 + 1944 * s_p_2 * t_p_4 - 2772 * t_p_5 +
          3024 * s * t_p_5 + 864 * t_p_6)) /
        25.;

      d2phi(33, 2) =
        (-1296 * s *
         (2 - 4 * s - 87 * t + 162 * s * t + 2 * s_p_2 - 75 * t * s_p_2 +
          942 * t_p_2 - 1560 * s * t_p_2 + 630 * s_p_2 * t_p_2 - 4150 * t_p_3 +
          5700 * s * t_p_3 - 1800 * s_p_2 * t_p_3 + 8595 * t_p_4 -
          8640 * s * t_p_4 + 1620 * s_p_2 * t_p_4 - 8316 * t_p_5 +
          4536 * s * t_p_5 + 3024 * t_p_6)) /
        25.;


      d2phi(34, 0) =
        (81 * (-1 + 2 * t) * (-1 + 3 * t) * (-1 + 6 * t) * t_p_2 *
         (16 - 51 * s - 34 * t + 54 * s * t + 36 * s_p_2 + 18 * t_p_2)) /
        2.;

      d2phi(34, 1) =
        (81 * t *
         (10 - 64 * s - 213 * t + 1260 * s * t + 102 * s_p_2 -
          1845 * t * s_p_2 - 48 * s_p_3 + 792 * t * s_p_3 + 1492 * t_p_2 -
          7744 * s * t_p_2 + 9720 * s_p_2 * t_p_2 - 3456 * s_p_3 * t_p_2 -
          4745 * t_p_3 + 19980 * s * t_p_3 - 18900 * s_p_2 * t_p_3 +
          4320 * s_p_3 * t_p_3 + 7524 * t_p_4 - 22464 * s * t_p_4 +
          11664 * s_p_2 * t_p_4 - 5796 * t_p_5 + 9072 * s * t_p_5 +
          1728 * t_p_6)) /
        4.;

      d2phi(34, 2) =
        (81 * s *
         (5 - 16 * s - 213 * t + 630 * s * t + 17 * s_p_2 - 615 * t * s_p_2 -
          6 * s_p_3 + 198 * t * s_p_3 + 2238 * t_p_2 - 5808 * s * t_p_2 +
          4860 * s_p_2 * t_p_2 - 1296 * s_p_3 * t_p_2 - 9490 * t_p_3 +
          19980 * s * t_p_3 - 12600 * s_p_2 * t_p_3 + 2160 * s_p_3 * t_p_3 +
          18810 * t_p_4 - 28080 * s * t_p_4 + 9720 * s_p_2 * t_p_4 -
          17388 * t_p_5 + 13608 * s * t_p_5 + 6048 * t_p_6)) /
        2.;

      d2phi(35, 0) =
        (81 * (-1 + 6 * t) * t_p_2 *
         (67 - 528 * s - 352 * t + 2043 * s * t + 1362 * s_p_2 -
          3456 * t * s_p_2 - 1440 * s_p_3 + 1800 * t * s_p_3 + 540 * s_p_4 +
          681 * t_p_2 - 2592 * s * t_p_2 + 2160 * s_p_2 * t_p_2 - 576 * t_p_3 +
          1080 * s * t_p_3 + 180 * t_p_4)) /
        2.;

      d2phi(35, 1) =
        (81 * t *
         (20 - 268 * s - 381 * t + 4524 * s * t + 1056 * s_p_2 -
          15633 * t * s_p_2 - 1816 * s_p_3 + 23256 * t * s_p_3 + 1440 * s_p_4 -
          15660 * t * s_p_4 - 432 * s_p_5 + 3888 * t * s_p_5 + 2312 * t_p_2 -
          22344 * s * t_p_2 + 59400 * s_p_2 * t_p_2 - 61056 * s_p_3 * t_p_2 +
          21600 * s_p_4 * t_p_2 - 6415 * t_p_3 + 46620 * s * t_p_3 -
          83160 * s_p_2 * t_p_3 + 43200 * s_p_3 * t_p_3 + 9036 * t_p_4 -
          43632 * s * t_p_4 + 38880 * s_p_2 * t_p_4 - 6300 * t_p_5 +
          15120 * s * t_p_5 + 1728 * t_p_6)) /
        4.;

      d2phi(35, 2) =
        (81 * s *
         (10 - 67 * s - 381 * t + 2262 * s * t + 176 * s_p_2 -
          5211 * t * s_p_2 - 227 * s_p_3 + 5814 * t * s_p_3 + 144 * s_p_4 -
          3132 * t * s_p_4 - 36 * s_p_5 + 648 * t * s_p_5 + 3468 * t_p_2 -
          16758 * s * t_p_2 + 29700 * s_p_2 * t_p_2 - 22896 * s_p_3 * t_p_2 +
          6480 * s_p_4 * t_p_2 - 12830 * t_p_3 + 46620 * s * t_p_3 -
          55440 * s_p_2 * t_p_3 + 21600 * s_p_3 * t_p_3 + 22590 * t_p_4 -
          54540 * s * t_p_4 + 32400 * s_p_2 * t_p_4 - 18900 * t_p_5 +
          22680 * s * t_p_5 + 6048 * t_p_6)) /
        2.;

      d2phi(36, 0) =
        (-1296 * t_p_2 *
         (-97 + 1131 * s + 754 * t - 6795 * s * t - 4530 * s_p_2 +
          19800 * t * s_p_2 + 8250 * s_p_3 - 23400 * t * s_p_3 - 7020 * s_p_4 +
          9720 * t * s_p_4 + 2268 * s_p_5 - 2265 * t_p_2 + 14850 * s * t_p_2 -
          28080 * s_p_2 * t_p_2 + 16200 * s_p_3 * t_p_2 + 3300 * t_p_3 -
          14040 * s * t_p_3 + 12960 * s_p_2 * t_p_3 - 2340 * t_p_4 +
          4860 * s * t_p_4 + 648 * t_p_5)) /
        25.;

      d2phi(36, 1) =
        (-648 * t *
         (20 - 388 * s - 291 * t + 4524 * s * t + 2262 * s_p_2 -
          20385 * t * s_p_2 - 6040 * s_p_3 + 39600 * t * s_p_3 + 8250 * s_p_4 -
          35100 * t * s_p_4 - 5616 * s_p_5 + 11664 * t * s_p_5 + 1512 * s_p_6 +
          1508 * t_p_2 - 18120 * s * t_p_2 + 59400 * s_p_2 * t_p_2 -
          74880 * s_p_3 * t_p_2 + 32400 * s_p_4 * t_p_2 - 3775 * t_p_3 +
          33000 * s * t_p_3 - 70200 * s_p_2 * t_p_3 + 43200 * s_p_3 * t_p_3 +
          4950 * t_p_4 - 28080 * s * t_p_4 + 29160 * s_p_2 * t_p_4 -
          3276 * t_p_5 + 9072 * s * t_p_5 + 864 * t_p_6)) /
        25.;


      d2phi(36, 2) =
        (-1296 * s *
         (10 - 97 * s - 291 * t + 2262 * s * t + 377 * s_p_2 -
          6795 * t * s_p_2 - 755 * s_p_3 + 9900 * t * s_p_3 + 825 * s_p_4 -
          7020 * t * s_p_4 - 468 * s_p_5 + 1944 * t * s_p_5 + 108 * s_p_6 +
          2262 * t_p_2 - 13590 * s * t_p_2 + 29700 * s_p_2 * t_p_2 -
          28080 * s_p_3 * t_p_2 + 9720 * s_p_4 * t_p_2 - 7550 * t_p_3 +
          33000 * s * t_p_3 - 46800 * s_p_2 * t_p_3 + 21600 * s_p_3 * t_p_3 +
          12375 * t_p_4 - 35100 * s * t_p_4 + 24300 * s_p_2 * t_p_4 -
          9828 * t_p_5 + 13608 * s * t_p_5 + 3024 * t_p_6)) /
        25.;

      d2phi(37, 0) =
        (-1296 * t *
         (10 - 291 * s - 97 * t + 2262 * s * t + 2262 * s_p_2 -
          13590 * t * s_p_2 - 7550 * s_p_3 + 33000 * t * s_p_3 + 12375 * s_p_4 -
          35100 * t * s_p_4 - 9828 * s_p_5 + 13608 * t * s_p_5 + 3024 * s_p_6 +
          377 * t_p_2 - 6795 * s * t_p_2 + 29700 * s_p_2 * t_p_2 -
          46800 * s_p_3 * t_p_2 + 24300 * s_p_4 * t_p_2 - 755 * t_p_3 +
          9900 * s * t_p_3 - 28080 * s_p_2 * t_p_3 + 21600 * s_p_3 * t_p_3 +
          825 * t_p_4 - 7020 * s * t_p_4 + 9720 * s_p_2 * t_p_4 - 468 * t_p_5 +
          1944 * s * t_p_5 + 108 * t_p_6)) /
        25.;

      d2phi(37, 1) =
        (-648 * s *
         (20 - 291 * s - 388 * t + 4524 * s * t + 1508 * s_p_2 -
          18120 * t * s_p_2 - 3775 * s_p_3 + 33000 * t * s_p_3 + 4950 * s_p_4 -
          28080 * t * s_p_4 - 3276 * s_p_5 + 9072 * t * s_p_5 + 864 * s_p_6 +
          2262 * t_p_2 - 20385 * s * t_p_2 + 59400 * s_p_2 * t_p_2 -
          70200 * s_p_3 * t_p_2 + 29160 * s_p_4 * t_p_2 - 6040 * t_p_3 +
          39600 * s * t_p_3 - 74880 * s_p_2 * t_p_3 + 43200 * s_p_3 * t_p_3 +
          8250 * t_p_4 - 35100 * s * t_p_4 + 32400 * s_p_2 * t_p_4 -
          5616 * t_p_5 + 11664 * s * t_p_5 + 1512 * t_p_6)) /
        25.;


      d2phi(37, 2) =
        (-1296 * s_p_2 *
         (-97 + 754 * s + 1131 * t - 6795 * s * t - 2265 * s_p_2 +
          14850 * t * s_p_2 + 3300 * s_p_3 - 14040 * t * s_p_3 - 2340 * s_p_4 +
          4860 * t * s_p_4 + 648 * s_p_5 - 4530 * t_p_2 + 19800 * s * t_p_2 -
          28080 * s_p_2 * t_p_2 + 12960 * s_p_3 * t_p_2 + 8250 * t_p_3 -
          23400 * s * t_p_3 + 16200 * s_p_2 * t_p_3 - 7020 * t_p_4 +
          9720 * s * t_p_4 + 2268 * t_p_5)) /
        25.;

      d2phi(38, 0) = (81 * t *
                      (10 - 381 * s - 67 * t + 2262 * s * t + 3468 * s_p_2 -
                       16758 * t * s_p_2 - 12830 * s_p_3 + 46620 * t * s_p_3 +
                       22590 * s_p_4 - 54540 * t * s_p_4 - 18900 * s_p_5 +
                       22680 * t * s_p_5 + 6048 * s_p_6 + 176 * t_p_2 -
                       5211 * s * t_p_2 + 29700 * s_p_2 * t_p_2 -
                       55440 * s_p_3 * t_p_2 + 32400 * s_p_4 * t_p_2 -
                       227 * t_p_3 + 5814 * s * t_p_3 - 22896 * s_p_2 * t_p_3 +
                       21600 * s_p_3 * t_p_3 + 144 * t_p_4 - 3132 * s * t_p_4 +
                       6480 * s_p_2 * t_p_4 - 36 * t_p_5 + 648 * s * t_p_5)) /
                     2.;

      d2phi(38, 1) =
        (81 * s *
         (20 - 381 * s - 268 * t + 4524 * s * t + 2312 * s_p_2 -
          22344 * t * s_p_2 - 6415 * s_p_3 + 46620 * t * s_p_3 + 9036 * s_p_4 -
          43632 * t * s_p_4 - 6300 * s_p_5 + 15120 * t * s_p_5 + 1728 * s_p_6 +
          1056 * t_p_2 - 15633 * s * t_p_2 + 59400 * s_p_2 * t_p_2 -
          83160 * s_p_3 * t_p_2 + 38880 * s_p_4 * t_p_2 - 1816 * t_p_3 +
          23256 * s * t_p_3 - 61056 * s_p_2 * t_p_3 + 43200 * s_p_3 * t_p_3 +
          1440 * t_p_4 - 15660 * s * t_p_4 + 21600 * s_p_2 * t_p_4 -
          432 * t_p_5 + 3888 * s * t_p_5)) /
        4.;

      d2phi(38, 2) =
        (81 * (-1 + 6 * s) * s_p_2 *
         (67 - 352 * s - 528 * t + 2043 * s * t + 681 * s_p_2 -
          2592 * t * s_p_2 - 576 * s_p_3 + 1080 * t * s_p_3 + 180 * s_p_4 +
          1362 * t_p_2 - 3456 * s * t_p_2 + 2160 * s_p_2 * t_p_2 -
          1440 * t_p_3 + 1800 * s * t_p_3 + 540 * t_p_4)) /
        2.;

      d2phi(39, 0) =
        (81 * t *
         (5 - 213 * s - 16 * t + 630 * s * t + 2238 * s_p_2 - 5808 * t * s_p_2 -
          9490 * s_p_3 + 19980 * t * s_p_3 + 18810 * s_p_4 - 28080 * t * s_p_4 -
          17388 * s_p_5 + 13608 * t * s_p_5 + 6048 * s_p_6 + 17 * t_p_2 -
          615 * s * t_p_2 + 4860 * s_p_2 * t_p_2 - 12600 * s_p_3 * t_p_2 +
          9720 * s_p_4 * t_p_2 - 6 * t_p_3 + 198 * s * t_p_3 -
          1296 * s_p_2 * t_p_3 + 2160 * s_p_3 * t_p_3)) /
        2.;

      d2phi(39, 1) =
        (81 * s *
         (10 - 213 * s - 64 * t + 1260 * s * t + 1492 * s_p_2 -
          7744 * t * s_p_2 - 4745 * s_p_3 + 19980 * t * s_p_3 + 7524 * s_p_4 -
          22464 * t * s_p_4 - 5796 * s_p_5 + 9072 * t * s_p_5 + 1728 * s_p_6 +
          102 * t_p_2 - 1845 * s * t_p_2 + 9720 * s_p_2 * t_p_2 -
          18900 * s_p_3 * t_p_2 + 11664 * s_p_4 * t_p_2 - 48 * t_p_3 +
          792 * s * t_p_3 - 3456 * s_p_2 * t_p_3 + 4320 * s_p_3 * t_p_3)) /
        4.;

      d2phi(39, 2) =
        (81 * (-1 + 2 * s) * (-1 + 3 * s) * (-1 + 6 * s) * s_p_2 *
         (16 - 34 * s - 51 * t + 54 * s * t + 18 * s_p_2 + 36 * t_p_2)) /
        2.;

      d2phi(40, 0) =
        (-1296 * t *
         (2 - 87 * s - 4 * t + 162 * s * t + 942 * s_p_2 - 1560 * t * s_p_2 -
          4150 * s_p_3 + 5700 * t * s_p_3 + 8595 * s_p_4 - 8640 * t * s_p_4 -
          8316 * s_p_5 + 4536 * t * s_p_5 + 3024 * s_p_6 + 2 * t_p_2 -
          75 * s * t_p_2 + 630 * s_p_2 * t_p_2 - 1800 * s_p_3 * t_p_2 +
          1620 * s_p_4 * t_p_2)) /
        25.;

      d2phi(40, 1) =
        (-648 * s *
         (4 - 87 * s - 16 * t + 324 * s * t + 628 * s_p_2 - 2080 * t * s_p_2 -
          2075 * s_p_3 + 5700 * t * s_p_3 + 3438 * s_p_4 - 6912 * t * s_p_4 -
          2772 * s_p_5 + 3024 * t * s_p_5 + 864 * s_p_6 + 12 * t_p_2 -
          225 * s * t_p_2 + 1260 * s_p_2 * t_p_2 - 2700 * s_p_3 * t_p_2 +
          1944 * s_p_4 * t_p_2)) /
        25.;

      d2phi(40, 2) = (-1296 * (-1 + 2 * s) * (-2 + 3 * s) * (-1 + 3 * s) *
                      (-1 + 6 * s) * (-2 + 2 * s + 3 * t) * s_p_2) /
                     25.;

      d2phi(41, 0) = (648 *
                      (-2 + 81 * s + 2 * t - 75 * s * t - 780 * s_p_2 +
                       630 * t * s_p_2 + 2850 * s_p_3 - 1800 * t * s_p_3 -
                       4320 * s_p_4 + 1620 * t * s_p_4 + 2268 * s_p_5) *
                      t_p_2 * sqrt(2)) /
                     25.;

      d2phi(41, 1) = (324 * s * t *
                      (-8 + 162 * s + 12 * t - 225 * s * t - 1040 * s_p_2 +
                       1260 * t * s_p_2 + 2850 * s_p_3 - 2700 * t * s_p_3 -
                       3456 * s_p_4 + 1944 * t * s_p_4 + 1512 * s_p_5) *
                      sqrt(2)) /
                     25.;

      d2phi(41, 2) = (648 * (-1 + 2 * s) * (-2 + 3 * s) * (-1 + 3 * s) *
                      (-1 + 6 * s) * (-1 + s + 3 * t) * s_p_2 * sqrt(2)) /
                     25.;

      d2phi(42, 0) =
        (81 * (-1 + 6 * t) *
         (1 - 36 * s - t + 33 * s * t + 282 * s_p_2 - 216 * t * s_p_2 -
          720 * s_p_3 + 360 * t * s_p_3 + 540 * s_p_4) *
         t_p_2 * pow(sqrt(2), -1)) /
        2.;

      d2phi(42, 1) =
        (81 * s * t *
         (-4 + 72 * s + 42 * t - 747 * s * t - 376 * s_p_2 + 3816 * t * s_p_2 +
          720 * s_p_3 - 7020 * t * s_p_3 - 432 * s_p_4 + 3888 * t * s_p_4 -
          48 * t_p_2 + 792 * s * t_p_2 - 3456 * s_p_2 * t_p_2 +
          4320 * s_p_3 * t_p_2) *
         pow(sqrt(2), -1)) /
        4.;

      d2phi(42, 2) =
        (81 * (-1 + 2 * s) * (-1 + 3 * s) * (-1 + 6 * s) * s_p_2 *
         (1 - s - 21 * t + 18 * s * t + 36 * t_p_2) * pow(sqrt(2), -1)) /
        2.;


      d2phi(43, 0) = (81 * (-1 + 2 * t) * (-1 + 3 * t) * (-1 + 6 * t) *
                      (1 - 21 * s - t + 18 * s * t + 36 * s_p_2) * t_p_2 *
                      pow(sqrt(2), -1)) /
                     2.;

      d2phi(43, 1) = (81 * s * t *
                      (-4 + 42 * s + 72 * t - 747 * s * t - 48 * s_p_2 +
                       792 * t * s_p_2 - 376 * t_p_2 + 3816 * s * t_p_2 -
                       3456 * s_p_2 * t_p_2 + 720 * t_p_3 - 7020 * s * t_p_3 +
                       4320 * s_p_2 * t_p_3 - 432 * t_p_4 + 3888 * s * t_p_4) *
                      pow(sqrt(2), -1)) /
                     4.;


      d2phi(43, 2) =
        (81 * (-1 + 6 * s) * s_p_2 *
         (1 - s - 36 * t + 33 * s * t + 282 * t_p_2 - 216 * s * t_p_2 -
          720 * t_p_3 + 360 * s * t_p_3 + 540 * t_p_4) *
         pow(sqrt(2), -1)) /
        2.;

      d2phi(44, 0) = (648 * (-1 + 3 * s + t) * (-1 + 2 * t) * (-2 + 3 * t) *
                      (-1 + 3 * t) * (-1 + 6 * t) * t_p_2 * sqrt(2)) /
                     25;

      d2phi(44, 1) = (324 * s * t *
                      (-8 + 12 * s + 162 * t - 225 * s * t - 1040 * t_p_2 +
                       1260 * s * t_p_2 + 2850 * t_p_3 - 2700 * s * t_p_3 -
                       3456 * t_p_4 + 1944 * s * t_p_4 + 1512 * t_p_5) *
                      sqrt(2)) /
                     25;

      d2phi(44, 2) = (648 * s_p_2 *
                      (-2 + 2 * s + 81 * t - 75 * s * t - 780 * t_p_2 +
                       630 * s * t_p_2 + 2850 * t_p_3 - 1800 * s * t_p_3 -
                       4320 * t_p_4 + 1620 * s * t_p_4 + 2268 * t_p_5) *
                      sqrt(2)) /
                     25;

      d2phi(45, 0) = 5832 * (-1 + 2 * t) * (-1 + 3 * t) * (-1 + 6 * t) * t_p_2 *
                     (1 - 6 * s - 2 * t + 6 * s * t + 6 * s_p_2 + t_p_2);

      d2phi(45, 1) =
        5832 * s * t *
        (-2 + 6 * s + 39 * t - 108 * s * t - 4 * s_p_2 + 66 * t * s_p_2 -
         236 * t_p_2 + 564 * s * t_p_2 - 288 * s_p_2 * t_p_2 + 595 * t_p_3 -
         1080 * s * t_p_3 + 360 * s_p_2 * t_p_3 - 648 * t_p_4 +
         648 * s * t_p_4 + 252 * t_p_5);


      d2phi(45, 2) =
        5832 * s_p_2 *
        (-1 + 2 * s + 39 * t - 72 * s * t - s_p_2 + 33 * t * s_p_2 -
         354 * t_p_2 + 564 * s * t_p_2 - 216 * s_p_2 * t_p_2 + 1190 * t_p_3 -
         1440 * s * t_p_3 + 360 * s_p_2 * t_p_3 - 1620 * t_p_4 +
         1080 * s * t_p_4 + 756 * t_p_5);

      d2phi(46, 0) =
        -2592 * (-1 + 3 * t) * (-1 + 6 * t) * t_p_2 *
        (-5 + 48 * s + 16 * t - 102 * s * t - 102 * s_p_2 + 108 * t * s_p_2 +
         60 * s_p_3 - 17 * t_p_2 + 54 * s * t_p_2 + 6 * t_p_3);

      d2phi(46, 1) =
        -2592 * s * t *
        (-10 + 48 * s + 183 * t - 801 * s * t - 68 * s_p_2 + 1026 * t * s_p_2 +
         30 * s_p_3 - 405 * t * s_p_3 - 1004 * t_p_2 + 3672 * s * t_p_2 -
         3744 * s_p_2 * t_p_2 + 1080 * s_p_3 * t_p_2 + 2235 * t_p_3 -
         5805 * s * t_p_3 + 3240 * s_p_2 * t_p_3 - 2160 * t_p_4 +
         2916 * s * t_p_4 + 756 * t_p_5);

      d2phi(46, 2) =
        -2592 * s_p_2 *
        (-5 + 16 * s + 183 * t - 534 * s * t - 17 * s_p_2 + 513 * t * s_p_2 +
         6 * s_p_3 - 162 * t * s_p_3 - 1506 * t_p_2 + 3672 * s * t_p_2 -
         2808 * s_p_2 * t_p_2 + 648 * s_p_3 * t_p_2 + 4470 * t_p_3 -
         7740 * s * t_p_3 + 3240 * s_p_2 * t_p_3 - 5400 * t_p_4 +
         4860 * s * t_p_4 + 2268 * t_p_5);

      d2phi(47, 0) =
        2592 * (-1 + 6 * t) * t_p_2 *
        (10 - 141 * s - 47 * t + 492 * s * t + 492 * s_p_2 - 1134 * t * s_p_2 -
         630 * s_p_3 + 720 * t * s_p_3 + 270 * s_p_4 + 82 * t_p_2 -
         567 * s * t_p_2 + 648 * s_p_2 * t_p_2 - 63 * t_p_3 + 216 * s * t_p_3 +
         18 * t_p_4);

      d2phi(47, 1) =
        2592 * s * t *
        (-20 + 141 * s + 321 * t - 2007 * s * t - 328 * s_p_2 +
         4086 * t * s_p_2 + 315 * s_p_3 - 3375 * t * s_p_3 - 108 * s_p_4 +
         972 * t * s_p_4 - 1456 * t_p_2 + 7038 * s * t_p_2 -
         9936 * s_p_2 * t_p_2 + 4320 * s_p_3 * t_p_2 + 2775 * t_p_3 -
         9045 * s * t_p_3 + 6480 * s_p_2 * t_p_3 - 2376 * t_p_4 +
         3888 * s * t_p_4 + 756 * t_p_5);

      d2phi(47, 2) =
        2592 * s_p_2 *
        (-10 + 47 * s + 321 * t - 1338 * s * t - 82 * s_p_2 + 2043 * t * s_p_2 +
         63 * s_p_3 - 1350 * t * s_p_3 - 18 * s_p_4 + 324 * t * s_p_4 -
         2184 * t_p_2 + 7038 * s * t_p_2 - 7452 * s_p_2 * t_p_2 +
         2592 * s_p_3 * t_p_2 + 5550 * t_p_3 - 12060 * s * t_p_3 +
         6480 * s_p_2 * t_p_3 - 5940 * t_p_4 + 6480 * s * t_p_4 + 2268 * t_p_5);

      d2phi(48, 0) =
        -5832 * t_p_2 *
        (-10 + 201 * s + 67 * t - 1056 * s * t - 1056 * s_p_2 +
         4086 * t * s_p_2 + 2270 * s_p_3 - 5760 * t * s_p_3 - 2160 * s_p_4 +
         2700 * t * s_p_4 + 756 * s_p_5 - 176 * t_p_2 + 2043 * s * t_p_2 -
         5184 * s_p_2 * t_p_2 + 3600 * s_p_3 * t_p_2 + 227 * t_p_3 -
         1728 * s * t_p_3 + 2160 * s_p_2 * t_p_3 - 144 * t_p_4 +
         540 * s * t_p_4 + 36 * t_p_5);

      d2phi(48, 1) =
        -5832 * s * t *
        (-20 + 201 * s + 201 * t - 1584 * s * t - 704 * s_p_2 +
         4086 * t * s_p_2 + 1135 * s_p_3 - 4320 * t * s_p_3 - 864 * s_p_4 +
         1620 * t * s_p_4 + 252 * s_p_5 - 704 * t_p_2 + 4086 * s * t_p_2 -
         6912 * s_p_2 * t_p_2 + 3600 * s_p_3 * t_p_2 + 1135 * t_p_3 -
         4320 * s * t_p_3 + 3600 * s_p_2 * t_p_3 - 864 * t_p_4 +
         1620 * s * t_p_4 + 252 * t_p_5);

      d2phi(48, 2) =
        -5832 * s_p_2 *
        (-10 + 67 * s + 201 * t - 1056 * s * t - 176 * s_p_2 +
         2043 * t * s_p_2 + 227 * s_p_3 - 1728 * t * s_p_3 - 144 * s_p_4 +
         540 * t * s_p_4 + 36 * s_p_5 - 1056 * t_p_2 + 4086 * s * t_p_2 -
         5184 * s_p_2 * t_p_2 + 2160 * s_p_3 * t_p_2 + 2270 * t_p_3 -
         5760 * s * t_p_3 + 3600 * s_p_2 * t_p_3 - 2160 * t_p_4 +
         2700 * s * t_p_4 + 756 * t_p_5);

      d2phi(49, 0) = 2592 * t_p_2 *
                     (-10 + 321 * s + 47 * t - 1338 * s * t - 2184 * s_p_2 +
                      7038 * t * s_p_2 + 5550 * s_p_3 - 12060 * t * s_p_3 -
                      5940 * s_p_4 + 6480 * t * s_p_4 + 2268 * s_p_5 -
                      82 * t_p_2 + 2043 * s * t_p_2 - 7452 * s_p_2 * t_p_2 +
                      6480 * s_p_3 * t_p_2 + 63 * t_p_3 - 1350 * s * t_p_3 +
                      2592 * s_p_2 * t_p_3 - 18 * t_p_4 + 324 * s * t_p_4);

      d2phi(49, 1) = 2592 * s * t *
                     (-20 + 321 * s + 141 * t - 2007 * s * t - 1456 * s_p_2 +
                      7038 * t * s_p_2 + 2775 * s_p_3 - 9045 * t * s_p_3 -
                      2376 * s_p_4 + 3888 * t * s_p_4 + 756 * s_p_5 -
                      328 * t_p_2 + 4086 * s * t_p_2 - 9936 * s_p_2 * t_p_2 +
                      6480 * s_p_3 * t_p_2 + 315 * t_p_3 - 3375 * s * t_p_3 +
                      4320 * s_p_2 * t_p_3 - 108 * t_p_4 + 972 * s * t_p_4);

      d2phi(49, 2) =
        2592 * (-1 + 6 * s) * s_p_2 *
        (10 - 47 * s - 141 * t + 492 * s * t + 82 * s_p_2 - 567 * t * s_p_2 -
         63 * s_p_3 + 216 * t * s_p_3 + 18 * s_p_4 + 492 * t_p_2 -
         1134 * s * t_p_2 + 648 * s_p_2 * t_p_2 - 630 * t_p_3 +
         720 * s * t_p_3 + 270 * t_p_4);

      d2phi(50, 0) =
        -2592 * t_p_2 *
        (-5 + 183 * s + 16 * t - 534 * s * t - 1506 * s_p_2 + 3672 * t * s_p_2 +
         4470 * s_p_3 - 7740 * t * s_p_3 - 5400 * s_p_4 + 4860 * t * s_p_4 +
         2268 * s_p_5 - 17 * t_p_2 + 513 * s * t_p_2 - 2808 * s_p_2 * t_p_2 +
         3240 * s_p_3 * t_p_2 + 6 * t_p_3 - 162 * s * t_p_3 +
         648 * s_p_2 * t_p_3);

      d2phi(50, 1) =
        -2592 * s * t *
        (-10 + 183 * s + 48 * t - 801 * s * t - 1004 * s_p_2 +
         3672 * t * s_p_2 + 2235 * s_p_3 - 5805 * t * s_p_3 - 2160 * s_p_4 +
         2916 * t * s_p_4 + 756 * s_p_5 - 68 * t_p_2 + 1026 * s * t_p_2 -
         3744 * s_p_2 * t_p_2 + 3240 * s_p_3 * t_p_2 + 30 * t_p_3 -
         405 * s * t_p_3 + 1080 * s_p_2 * t_p_3);

      d2phi(50, 2) =
        -2592 * (-1 + 3 * s) * (-1 + 6 * s) * s_p_2 *
        (-5 + 16 * s + 48 * t - 102 * s * t - 17 * s_p_2 + 54 * t * s_p_2 +
         6 * s_p_3 - 102 * t_p_2 + 108 * s * t_p_2 + 60 * t_p_3);

      d2phi(51, 0) =
        5832 * t_p_2 *
        (-1 + 39 * s + 2 * t - 72 * s * t - 354 * s_p_2 + 564 * t * s_p_2 +
         1190 * s_p_3 - 1440 * t * s_p_3 - 1620 * s_p_4 + 1080 * t * s_p_4 +
         756 * s_p_5 - t_p_2 + 33 * s * t_p_2 - 216 * s_p_2 * t_p_2 +
         360 * s_p_3 * t_p_2);


      d2phi(51, 1) =
        5832 * s * t *
        (-2 + 39 * s + 6 * t - 108 * s * t - 236 * s_p_2 + 564 * t * s_p_2 +
         595 * s_p_3 - 1080 * t * s_p_3 - 648 * s_p_4 + 648 * t * s_p_4 +
         252 * s_p_5 - 4 * t_p_2 + 66 * s * t_p_2 - 288 * s_p_2 * t_p_2 +
         360 * s_p_3 * t_p_2);

      d2phi(51, 2) = 5832 * (-1 + 2 * s) * (-1 + 3 * s) * (-1 + 6 * s) * s_p_2 *
                     (1 - 2 * s - 6 * t + 6 * s * t + s_p_2 + 6 * t_p_2);

      d2phi(52, 0) =
        2592 * (-1 + 6 * t) * t_p_2 *
        (1 - 33 * s - 2 * t + 60 * s * t + 222 * s_p_2 - 324 * t * s_p_2 -
         450 * s_p_3 + 360 * t * s_p_3 + 270 * s_p_4 + t_p_2 - 27 * s * t_p_2 +
         108 * s_p_2 * t_p_2);

      d2phi(52, 1) =
        2592 * s * t *
        (-2 + 33 * s + 24 * t - 387 * s * t - 148 * s_p_2 + 1656 * t * s_p_2 +
         225 * s_p_3 - 2295 * t * s_p_3 - 108 * s_p_4 + 972 * t * s_p_4 -
         52 * t_p_2 + 774 * s * t_p_2 - 2736 * s_p_2 * t_p_2 +
         2160 * s_p_3 * t_p_2 + 30 * t_p_3 - 405 * s * t_p_3 +
         1080 * s_p_2 * t_p_3);

      d2phi(52, 2) =
        2592 * (-1 + 3 * s) * (-1 + 6 * s) * s_p_2 *
        (-1 + 2 * s + 24 * t - 42 * s * t - s_p_2 + 18 * t * s_p_2 -
         78 * t_p_2 + 72 * s * t_p_2 + 60 * t_p_3);

      d2phi(53, 0) = 2592 * (-1 + 3 * t) * (-1 + 6 * t) * t_p_2 *
                     (-1 + 24 * s + 2 * t - 42 * s * t - 78 * s_p_2 +
                      72 * t * s_p_2 + 60 * s_p_3 - t_p_2 + 18 * s * t_p_2);

      d2phi(53, 1) = 2592 * s * t *
                     (-2 + 24 * s + 33 * t - 387 * s * t - 52 * s_p_2 +
                      774 * t * s_p_2 + 30 * s_p_3 - 405 * t * s_p_3 -
                      148 * t_p_2 + 1656 * s * t_p_2 - 2736 * s_p_2 * t_p_2 +
                      1080 * s_p_3 * t_p_2 + 225 * t_p_3 - 2295 * s * t_p_3 +
                      2160 * s_p_2 * t_p_3 - 108 * t_p_4 + 972 * s * t_p_4);

      d2phi(53, 2) = 2592 * (-1 + 6 * s) * s_p_2 *
                     (1 - 2 * s - 33 * t + 60 * s * t + s_p_2 - 27 * t * s_p_2 +
                      222 * t_p_2 - 324 * s * t_p_2 + 108 * s_p_2 * t_p_2 -
                      450 * t_p_3 + 360 * s * t_p_3 + 270 * t_p_4);


      d2phi(54, 0) =
        -1458 * (-1 + 6 * t) * t_p_2 *
        (5 - 138 * s - 16 * t + 390 * s * t + 678 * s_p_2 - 1332 * t * s_p_2 -
         1080 * s_p_3 + 1080 * t * s_p_3 + 540 * s_p_4 + 17 * t_p_2 -
         360 * s * t_p_2 + 648 * s_p_2 * t_p_2 - 6 * t_p_3 + 108 * s * t_p_3);

      d2phi(54, 1) = -1458 * s * t *
                     (-10 + 138 * s + 138 * t - 1827 * s * t - 452 * s_p_2 +
                      5400 * t * s_p_2 + 540 * s_p_3 - 5670 * t * s_p_3 -
                      216 * s_p_4 + 1944 * t * s_p_4 - 452 * t_p_2 +
                      5400 * s * t_p_2 - 11520 * s_p_2 * t_p_2 +
                      6480 * s_p_3 * t_p_2 + 540 * t_p_3 - 5670 * s * t_p_3 +
                      6480 * s_p_2 * t_p_3 - 216 * t_p_4 + 1944 * s * t_p_4);


      d2phi(54, 2) =
        -1458 * (-1 + 6 * s) * s_p_2 *
        (5 - 16 * s - 138 * t + 390 * s * t + 17 * s_p_2 - 360 * t * s_p_2 -
         6 * s_p_3 + 108 * t * s_p_3 + 678 * t_p_2 - 1332 * s * t_p_2 +
         648 * s_p_2 * t_p_2 - 1080 * t_p_3 + 1080 * s * t_p_3 + 540 * t_p_4);
    }

  } // namespace MyC1CurvedElements
} // namespace oomph

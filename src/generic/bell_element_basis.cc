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
#include "bell_element_basis.h"

namespace oomph
{
  namespace MyShape
  {
    /// Get the (twice) area of the triangle from the vertices
    double BellElementBasis::get_twice_triangle_area(
      const Vector<double>& v0,
      const Vector<double>& v1,
      const Vector<double>& v2) const
    {
      return v1[0] * v2[1] - v1[1] * v2[0] + v0[1] * v2[0] - v0[0] * v2[1] +
             v0[0] * v1[1] - v0[1] * v1[0];
    }

    /// Get outer normal of side between vertices v0 and v1, assumes
    /// counter-clockwise triangle vertices.
    Vector<double> BellElementBasis::get_outer_normal(
      const Vector<double>& v0, const Vector<double>& v1) const
    {
      Vector<double> normal(2);
      normal[0] = v1[1] - v0[1];
      normal[1] = v0[0] - v1[0];
      return normal;
    }

    /// Basis on a reference element. This follows exactly the notation of M.
    /// Okabe in Comput. Methods Appl. Mech. 117 (1994) 411-421
    void BellElementBasis::d2_basis(const Vector<double>& s,
                                    const Vector<Vector<double>>& v,
                                    Shape& psi,
                                    DShape& dpsi,
                                    DShape& d2psi) const
    {
      // Area of this triangle
      const double A = get_twice_triangle_area(v[0], v[1], v[2]);
      // The area coordinates: area of the triangle made by s and two other
      // vertices
      Vector<double> w(3);
      Vector<Vector<double>> gradw(3, Vector<double>(2));
      for (unsigned i = 0; i < 3; i++)
      {
        // Cyclic over base three
        unsigned j = (i + 1) % 3;
        unsigned k = (j + 1) % 3;

        // Now get the area of the triangle
        w[i] = get_twice_triangle_area(s, v[j], v[k]) / A;
        // Fill in gradient of area coordinates
        gradw[i][0] = (v[j][1] - v[k][1]) / A;
        gradw[i][1] = (v[k][0] - v[j][0]) / A;
      }

      // Normals (not unit) on the reference element
      Vector<Vector<double>> n(3, Vector<double>(2));
      for (unsigned i = 0; i < 3; ++i)
      {
        // Cyclic over base three
        unsigned j = (i + 1) % 3;
        unsigned k = (j + 1) % 3;
        // Normal on side opposite to node i
        n[i] = get_outer_normal(v[j], v[k]);
      }

      // A useful function for simplifying the definitions a little
      DenseMatrix<double> e(3, 3);
      for (unsigned i = 0; i < 3; i++)
      {
        for (unsigned j = 0; j < 3; j++)
        {
          e(i, j) = (gradw[i][0] * n[j][0] + gradw[i][1] * n[j][1]) /
                    (gradw[j][0] * n[j][0] + gradw[j][1] * n[j][1]);
        }
      }

      // Now we define the actual shape functions
      // Loop over the nodes
      for (unsigned i = 0; i < 3; ++i)
      {
        // Cyclic over base three
        unsigned j = (i + 1) % 3;
        unsigned k = (j + 1) % 3;

        // The 1st shape function - related to the zeroth derivative dofs
        psi(i, 0) =
          pow(w[i], 3) * (10 - 15 * w[i] + 6 * w[i] * w[i]) -
          30 * (e(i, j) * w[k] + e(i, k) * w[j]) * w[i] * w[i] * w[j] * w[k];

        // Initialise to the derivative ith respect to alpha
        dpsi(i, 0, 0) = 0.0;
        dpsi(i, 0, 1) = 0.0;

        // The derivatives with respect to wl
        Vector<double> dpsidwl(3, 0.0);
        DenseMatrix<double> d2psidwl(3, 3, 0.0);
        dpsidwl[i] =
          w[i] * w[i] * (30 - 60 * w[i] + 30 * w[i] * w[i]) -
          60 * (e(i, j) * w[k] + e(i, k) * w[j]) * w[i] * w[j] * w[k];
        dpsidwl[j] =
          -30 * (e(i, j) * w[k] + 2. * e(i, k) * w[j]) * w[i] * w[i] * w[k];
        dpsidwl[k] =
          -30 * (2. * e(i, j) * w[k] + e(i, k) * w[j]) * w[i] * w[i] * w[j];

        // Initialise the second partial  derivatives
        d2psi(i, 0, 0) = 0.0;
        d2psi(i, 0, 1) = 0.0;
        d2psi(i, 0, 2) = 0.0;

        // Second derivatives with respect to wl wll.
        d2psidwl(i, i) = 60 * w[i] * (1 - 3 * w[i] + 2 * w[i] * w[i]) -
                         60 * (e(i, j) * w[k] + e(i, k) * w[j]) * w[j] * w[k];
        d2psidwl(j, i) =
          -60 * (e(i, j) * w[k] + 2 * e(i, k) * w[j]) * w[i] * w[k];

        d2psidwl(k, i) =
          -60 * (2 * e(i, j) * w[k] + e(i, k) * w[j]) * w[i] * w[j];

        d2psidwl(i, j) = d2psidwl(j, i);

        d2psidwl(j, j) = -60 * e(i, k) * w[i] * w[i] * w[k];

        d2psidwl(k, j) = -60 * (e(i, j) * w[k] + e(i, k) * w[j]) * w[i] * w[i];

        d2psidwl(i, k) = d2psidwl(k, i);

        d2psidwl(j, k) = d2psidwl(k, j);

        d2psidwl(k, k) = -60 * e(i, j) * w[i] * w[i] * w[j];

        // Now construct the first derivatives
        for (unsigned alpha = 0; alpha < 2; ++alpha)
        {
          for (unsigned l = 0; l < 3; ++l)
          {
            // add on dpsidwi dwi ds
            dpsi(i, 0, alpha) += dpsidwl[l] * gradw[l][alpha];
          }

          // Only need to fill in 3, enumerate in the easiest manner
          for (unsigned beta = 0; beta < alpha + 1; ++beta)
          {
            // Initialise to the derivative ith respect to alpha
            for (unsigned l = 0; l < 3; ++l)
            {
              for (unsigned ll = 0; ll < 3; ++ll)
              {
                // add on d2psid2wi dwi ds
                d2psi(i, 0, alpha + beta) +=
                  d2psidwl(l, ll) * gradw[l][alpha] * gradw[ll][beta];
              }
            }
          }
        }


        // The 2nd and 3rd shape functions - relating to first derivative dofs
        for (unsigned m = 0; m < 2; ++m)
        {
          // mth components of vectors between nodes, defined for convenience
          const double lik(v[i][m] - v[k][m]), lij(v[i][m] - v[j][m]),
            ljk(v[j][m] - v[k][m]), lkj(v[k][m] - v[j][m]);

          psi(i, 1 + m) =
            pow(w[i], 3) * (4 - 3 * w[i]) * (s[m] - v[i][m]) +
            (15 * lik * e(i, j) + 3 * ljk) * pow(w[i] * w[k], 2) * w[j] +
            (15 * lij * e(i, k) + 3 * lkj) * pow(w[i] * w[j], 2) * w[k];

          // The derivatives with respect to wl
          dpsidwl[i] =
            12 * pow(w[i], 2) * (1 - w[i]) * (s[m] - v[i][m]) +
            (15 * lik * e(i, j) + 3 * ljk) * pow(w[k], 2) * w[j] * w[i] * 2. +
            (15 * lij * e(i, k) + 3 * lkj) * pow(w[j], 2) * w[k] * w[i] * 2.;

          dpsidwl[j] =
            (15 * lik * e(i, j) + 3 * ljk) * pow(w[i] * w[k], 2) +
            (15 * lij * e(i, k) + 3 * lkj) * pow(w[i], 2) * w[k] * w[j] * 2.;

          dpsidwl[k] =
            (15 * lik * e(i, j) + 3 * ljk) * pow(w[i], 2) * w[j] * w[k] * 2. +
            (15 * lij * e(i, k) + 3 * lkj) * pow(w[i] * w[j], 2);

          // Initialise to the derivative with respect to alpha
          d2psi(i, 1 + m, 0) = 0.0;
          d2psi(i, 1 + m, 1) = 0.0;
          d2psi(i, 1 + m, 2) = 0.0;

          // Second derivatives with respect to wl wll.
          d2psidwl(i, i) =
            12 * w[i] * (2 - 3 * w[i]) * (s[m] - v[i][m]) +
            (15 * lik * e(i, j) + 3 * ljk) * pow(w[k], 2) * w[j] * 2. +
            (15 * lij * e(i, k) + 3 * lkj) * pow(w[j], 2) * w[k] * 2.;

          d2psidwl(j, i) =
            (15 * lik * e(i, j) + 3 * ljk) * pow(w[k], 2) * w[i] * 2. +
            (15 * lij * e(i, k) + 3 * lkj) * w[j] * w[k] * w[i] * 4.;

          d2psidwl(k, i) =
            (15 * lik * e(i, j) + 3 * ljk) * w[k] * w[j] * w[i] * 4. +
            (15 * lij * e(i, k) + 3 * lkj) * pow(w[j], 2) * w[i] * 2.;

          d2psidwl(i, j) = d2psidwl(j, i);

          d2psidwl(j, j) =
            (15 * lij * e(i, k) + 3 * lkj) * pow(w[i], 2) * w[k] * 2.;

          d2psidwl(k, j) =
            (15 * lik * e(i, j) + 3 * ljk) * pow(w[i], 2) * w[k] * 2. +
            (15 * lij * e(i, k) + 3 * lkj) * pow(w[i], 2) * w[j] * 2.;

          d2psidwl(i, k) = d2psidwl(k, i);

          d2psidwl(j, k) = d2psidwl(k, j);

          d2psidwl(k, k) =
            (15 * lik * e(i, j) + 3 * ljk) * pow(w[i], 2) * w[j] * 2.;

          // The partial derivatives with respect to x_m
          dpsi(i, 1 + m, m) = pow(w[i], 3) * (4 - 3 * w[i]);
          dpsi(i, 1 + m, (m + 1) % 2) = 0.0;

          // The partial derivatives with respect to x_m
          //   d2psi(i,1+m,m)       = 24*pow(w[i],2)*(1-w[i])*gradw[i][0];
          //   d2psi(i,1+m,m+1)     = 24*pow(w[i],2)*(1-w[i])*gradw[i][1];

          // Now construct the first derivatives
          for (unsigned alpha = 0; alpha < 2; ++alpha)
          {
            for (unsigned l = 0; l < 3; ++l)
            {
              // add on dpsidwi dwi ds
              dpsi(i, 1 + m, alpha) += dpsidwl[l] * gradw[l][alpha];
            }

            // Only need to fill in 3, enumerate in the easiest manner
            for (unsigned beta = 0; beta < alpha + 1; ++beta)
            {
              if (alpha == m)
              {
                d2psi(i, 1 + m, alpha + beta) +=
                  12 * pow(w[i], 2) * (1 - w[i]) * gradw[i][beta];
              }
              if (beta == m)
              {
                d2psi(i, 1 + m, alpha + beta) +=
                  12 * pow(w[i], 2) * (1 - w[i]) * gradw[i][alpha];
              }
              for (unsigned l = 0; l < 3; ++l)
              {
                for (unsigned ll = 0; ll < 3; ++ll)
                {
                  // add on d2psid2wi dwi ds
                  d2psi(i, 1 + m, alpha + beta) +=
                    d2psidwl(l, ll) * gradw[l][alpha] * gradw[ll][beta];
                }
              }
            }
          }
        }

        // The 4th and 5th shape functions - relating to the second x and y
        // derivative
        //  dofs
        for (unsigned m = 0; m < 2; ++m)
        {
          // mth components of vectors between nodes, defined for convenience
          const double lik(v[i][m] - v[k][m]), lij(v[i][m] - v[j][m]),
            ljk(v[j][m] - v[k][m]), lkj(v[k][m] - v[j][m]);
          psi(i, 3 + 2 * m) = pow(w[i], 3) * pow(s[m] - v[i][m], 2) / 2. -
                              (5 * lik * e(i, j) + 2 * ljk) * lik *
                                pow(w[i] * w[k], 2) * w[j] / 2. -
                              (5 * lij * e(i, k) + 2 * lkj) * lij *
                                pow(w[i] * w[j], 2) * w[k] / 2.;

          // The derivatives with respect to wl
          Vector<double> dpsidwl(3, 0.0);
          dpsidwl[i] =
            3 * pow(w[i], 2) * pow(s[m] - v[i][m], 2) / 2. -
            (5 * lik * e(i, j) + 2 * ljk) * lik * pow(w[k], 2) * w[j] * w[i] -
            (5 * lij * e(i, k) + 2 * lkj) * lij * pow(w[j], 2) * w[k] * w[i];

          dpsidwl[j] =
            -(5 * lik * e(i, j) + 2 * ljk) * lik * pow(w[i] * w[k], 2) / 2. -
            (5 * lij * e(i, k) + 2 * lkj) * lij * pow(w[i], 2) * w[k] * w[j];

          dpsidwl[k] =
            -(5 * lik * e(i, j) + 2 * ljk) * lik * pow(w[i], 2) * w[j] * w[k] -
            (5 * lij * e(i, k) + 2 * lkj) * lij * pow(w[i] * w[j], 2) / 2.;

          // Second derivatives with respect to wl wll.
          d2psidwl(i, i) =
            6 * w[i] * pow(s[m] - v[i][m], 2) / 2. -
            (5 * lik * e(i, j) + 2 * ljk) * lik * pow(w[k], 2) * w[j] -
            (5 * lij * e(i, k) + 2 * lkj) * lij * pow(w[j], 2) * w[k];

          d2psidwl(j, i) =
            -(5 * lik * e(i, j) + 2 * ljk) * lik * pow(w[k], 2) * w[i] -
            (5 * lij * e(i, k) + 2 * lkj) * lij * w[j] * w[k] * w[i] * 2.;

          d2psidwl(k, i) =
            -(5 * lik * e(i, j) + 2 * ljk) * lik * w[k] * w[j] * w[i] * 2. -
            (5 * lij * e(i, k) + 2 * lkj) * lij * pow(w[j], 2) * w[i];

          d2psidwl(i, j) = d2psidwl(j, i);

          d2psidwl(j, j) =
            -(5 * lij * e(i, k) + 2 * lkj) * lij * pow(w[i], 2) * w[k];

          d2psidwl(k, j) =
            -(5 * lik * e(i, j) + 2 * ljk) * lik * pow(w[i], 2) * w[k] -
            (5 * lij * e(i, k) + 2 * lkj) * lij * pow(w[i], 2) * w[j];

          d2psidwl(i, k) = d2psidwl(k, i);

          d2psidwl(j, k) = d2psidwl(k, j);

          d2psidwl(k, k) =
            -(5 * lik * e(i, j) + 2 * ljk) * lik * pow(w[i], 2) * w[j];

          // The partial derivatives with respect to x_m
          dpsi(i, 3 + 2 * m, m) = pow(w[i], 3) * (s[m] - v[i][m]);
          dpsi(i, 3 + 2 * m, (m + 1) % 2) = 0.0;

          // The second partial derivatives with respect to x_m x_n
          d2psi(i, 3 + 2 * m, 0) = 0.0;
          d2psi(i, 3 + 2 * m, 1) = 0.0;
          d2psi(i, 3 + 2 * m, 2) = 0.0;
          // Direct partial derivative when m==n
          d2psi(i, 3 + 2 * m, m + m) += pow(w[i], 3);

          // Now construct the first derivatives
          for (unsigned alpha = 0; alpha < 2; ++alpha)
          {
            // Initialise to the derivative ith respect to alpha
            for (unsigned l = 0; l < 3; ++l)
            {
              // add on dpsidwi dwi ds
              dpsi(i, 3 + 2 * m, alpha) += dpsidwl[l] * gradw[l][alpha];
            }
            // Only need to fill in 3, enumerate in the easiest manner
            for (unsigned beta = 0; beta < alpha + 1; ++beta)
            {
              if (m == alpha)
                d2psi(i, 3 + 2 * m, alpha + beta) +=
                  3 * pow(w[i], 2) * (s[m] - v[i][m]) * gradw[i][beta];
              if (m == beta)
                d2psi(i, 3 + 2 * m, alpha + beta) +=
                  3 * pow(w[i], 2) * (s[m] - v[i][m]) * gradw[i][alpha];

              // Initialise to the derivative ith respect to alpha
              for (unsigned l = 0; l < 3; ++l)
              {
                for (unsigned ll = 0; ll < 3; ++ll)
                {
                  // add on d2psid2wi dwi ds
                  d2psi(i, 3 + 2 * m, alpha + beta) +=
                    d2psidwl(l, ll) * gradw[l][alpha] * gradw[ll][beta];
                }
              }
            }
          }
        }

        // The  6th shape functions - relating to the cross derivative dofs
        psi(i, 4) = pow(w[i], 3) * (s[0] - v[i][0]) * (s[1] - v[i][1]);
        // The initialise derivatives with respect to wl
        dpsidwl[i] = 3 * pow(w[i], 2) * (s[0] - v[i][0]) * (s[1] - v[i][1]);
        dpsidwl[j] = 0.0;
        dpsidwl[k] = 0.0;

        // The partial derivatives with respect to x_m
        dpsi(i, 4, 0) = pow(w[i], 3) * (s[1] - v[i][1]);
        dpsi(i, 4, 1) = pow(w[i], 3) * (s[0] - v[i][0]);

        // Initialise second partial derivatives with respect to wl
        d2psidwl = DenseMatrix<double>(3, 3, 0.0);
        d2psidwl(i, i) = 6 * w[i] * (s[0] - v[i][0]) * (s[1] - v[i][1]);

        // The second partial derivatives with respect to x_m x_n
        d2psi(i, 4, 0) = 0.0;
        d2psi(i, 4, 1) = pow(w[i], 3);
        d2psi(i, 4, 2) = 0.0;

        d2psi(i, 4, 0 + 0) += 6 * pow(w[i], 2) * (s[1] - v[i][1]) * gradw[i][0];
        d2psi(i, 4, 0 + 1) += 3 * pow(w[i], 2) * (s[1] - v[i][1]) * gradw[i][1];
        d2psi(i, 4, 1 + 1) += 6 * pow(w[i], 2) * (s[0] - v[i][0]) * gradw[i][1];
        d2psi(i, 4, 1 + 0) += 3 * pow(w[i], 2) * (s[0] - v[i][0]) * gradw[i][0];

        // Vectors between nodes, defined for convenience
        Vector<double> lik(2, 0.0), lij(2, 0.0);

        for (unsigned m = 0; m < 2; ++m)
        {
          unsigned n = (m + 1) % 2; // Cyclic in base two
          // Fill in the values
          lik[m] = (v[i][m] - v[k][m]);
          lij[m] = (v[i][m] - v[j][m]);
          lik[n] = (v[i][n] - v[k][n]);
          lij[n] = (v[i][n] - v[j][n]);

          // Vectors compenent between nodes, defined for convenience
          const double lkj(v[k][m] - v[j][m]), ljk(v[j][m] - v[k][m]);

          psi(i, 4) += -(5 * lik[m] * e(i, j) + 2 * ljk) * lik[n] *
                         pow(w[i] * w[k], 2) * w[j] / 2. -
                       (5 * lij[m] * e(i, k) + 2 * lkj) * lij[n] *
                         pow(w[i] * w[j], 2) * w[k] / 2.;

          // Now for the first derivatives
          dpsidwl[i] += -(5 * lik[m] * e(i, j) + 2 * ljk) * lik[n] *
                          pow(w[k], 2) * w[j] * w[i] -
                        (5 * lij[m] * e(i, k) + 2 * lkj) * lij[n] *
                          pow(w[j], 2) * w[k] * w[i];

          dpsidwl[j] += -(5 * lik[m] * e(i, j) + 2 * ljk) * lik[n] *
                          pow(w[i] * w[k], 2) / 2. -
                        (5 * lij[m] * e(i, k) + 2 * lkj) * lij[n] *
                          pow(w[i], 2) * w[k] * w[j];

          dpsidwl[k] += -(5 * lik[m] * e(i, j) + 2 * ljk) * lik[n] *
                          pow(w[i], 2) * w[j] * w[k] -
                        (5 * lij[m] * e(i, k) + 2 * lkj) * lij[n] *
                          pow(w[i] * w[j], 2) / 2.;

          // Second derivatives with respect to wl wll.
          d2psidwl(i, i) +=
            -(5 * lik[m] * e(i, j) + 2 * ljk) * lik[n] * pow(w[k], 2) * w[j] -
            (5 * lij[m] * e(i, k) + 2 * lkj) * lij[n] * pow(w[j], 2) * w[k];

          d2psidwl(j, i) +=
            -(5 * lik[m] * e(i, j) + 2 * ljk) * lik[n] * pow(w[k], 2) * w[i] -
            (5 * lij[m] * e(i, k) + 2 * lkj) * lij[n] * w[j] * w[k] * w[i] * 2.;

          d2psidwl(k, i) +=
            -(5 * lik[m] * e(i, j) + 2 * ljk) * lik[n] * w[k] * w[j] * w[i] *
              2. -
            (5 * lij[m] * e(i, k) + 2 * lkj) * lij[n] * pow(w[j], 2) * w[i];

          d2psidwl(j, j) +=
            -(5 * lij[m] * e(i, k) + 2 * lkj) * lij[n] * pow(w[i], 2) * w[k];

          d2psidwl(k, j) +=
            -(5 * lik[m] * e(i, j) + 2 * ljk) * lik[n] * pow(w[i], 2) * w[k] -
            (5 * lij[m] * e(i, k) + 2 * lkj) * lij[n] * pow(w[i], 2) * w[j];

          d2psidwl(k, k) +=
            -(5 * lik[m] * e(i, j) + 2 * ljk) * lik[n] * pow(w[i], 2) * w[j];
        }

        // The lower triangular terms
        d2psidwl(i, j) = d2psidwl(j, i);

        d2psidwl(i, k) = d2psidwl(k, i);

        d2psidwl(j, k) = d2psidwl(k, j);

        // Now construct the first derivatives
        for (unsigned alpha = 0; alpha < 2; ++alpha)
        {
          // Initialise to the derivative ith respect to alpha
          for (unsigned l = 0; l < 3; ++l)
          {
            // add on dpsidwi dwi ds
            dpsi(i, 4, alpha) += dpsidwl[l] * gradw[l][alpha];
          }
          // Only need to fill in 3, enumerate in the easiest manner
          for (unsigned beta = 0; beta < alpha + 1; ++beta)
          {
            for (unsigned l = 0; l < 3; ++l)
            {
              for (unsigned ll = 0; ll < 3; ++ll)
              {
                // add on d2psid2wi dwi ds
                d2psi(i, 4, alpha + beta) +=
                  d2psidwl(l, ll) * gradw[l][alpha] * gradw[ll][beta];
              }
            }
          }
        }
      }
    }


    /// Basis on a reference element. This follows exactly the notation of M.
    /// Okabe in Comput. Methods Appl. Mech. 117 (1994) 411-421
    double BellElementBasis::d2_basis_eulerian(const Vector<double>& s,
                                               const Vector<Vector<double>>& v,
                                               Shape& psi,
                                               DShape& dpsi,
                                               DShape& d2psi) const
    {
      // Area of this triangle
      const double A = get_twice_triangle_area(v[0], v[1], v[2]);
      // The area coordinates: area of the triangle made by s and two other
      // vertices
      Vector<double> w(3);
      w[0] = s[0];
      w[1] = s[1];
      w[2] = 1 - s[0] - s[1];

      // Get Jacobian
      DenseMatrix<double> J(2, 2);
      J(0, 0) = v[0][0] - v[2][0];
      J(1, 0) = v[0][1] - v[2][1];
      J(0, 1) = v[1][0] - v[2][0];
      J(1, 1) = v[1][1] - v[2][1];

      // Find the x coordinates, J s + x_2
      Vector<double> x(2, 0.0);
      x[0] = J(0, 0) * s[0] + v[2][0] + J(0, 1) * s[1];
      x[1] = J(1, 0) * s[0] + v[2][1] + J(1, 1) * s[1];

      Vector<Vector<double>> gradw(3, Vector<double>(2));
      for (unsigned i = 0; i < 3; i++)
      {
        // Cyclic over base three
        unsigned j = (i + 1) % 3;
        unsigned k = (j + 1) % 3;

        // Fill in gradient of area coordinates
        gradw[i][0] = (v[j][1] - v[k][1]) / A;
        gradw[i][1] = (v[k][0] - v[j][0]) / A;
      }

      // Normals (not unit) on the reference element
      Vector<Vector<double>> n(3, Vector<double>(2));
      for (unsigned i = 0; i < 3; ++i)
      {
        // Cyclic over base three
        unsigned j = (i + 1) % 3;
        unsigned k = (j + 1) % 3;
        // Normal on side opposite to node i
        n[i] = get_outer_normal(v[j], v[k]);
      }

      // A useful function for simplifying the definitions a little
      DenseMatrix<double> e(3, 3);
      for (unsigned i = 0; i < 3; i++)
      {
        for (unsigned j = 0; j < 3; j++)
        {
          e(i, j) = (gradw[i][0] * n[j][0] + gradw[i][1] * n[j][1]) /
                    (gradw[j][0] * n[j][0] + gradw[j][1] * n[j][1]);
        }
      }

      // Now we define the actual shape functions
      // Loop over the nodes
      for (unsigned i = 0; i < 3; ++i)
      {
        // Cyclic over base three
        unsigned j = (i + 1) % 3;
        unsigned k = (j + 1) % 3;

        // The 1st shape function - related to the zeroth derivative dofs
        psi(i, 0) =
          pow(w[i], 3) * (10 - 15 * w[i] + 6 * w[i] * w[i]) -
          30 * (e(i, j) * w[k] + e(i, k) * w[j]) * w[i] * w[i] * w[j] * w[k];

        // Initialise to the derivative ith respect to alpha
        dpsi(i, 0, 0) = 0.0;
        dpsi(i, 0, 1) = 0.0;

        // The derivatives with respect to wl
        Vector<double> dpsidwl(3, 0.0);
        DenseMatrix<double> d2psidwl(3, 3, 0.0);
        dpsidwl[i] =
          w[i] * w[i] * (30 - 60 * w[i] + 30 * w[i] * w[i]) -
          60 * (e(i, j) * w[k] + e(i, k) * w[j]) * w[i] * w[j] * w[k];
        dpsidwl[j] =
          -30 * (e(i, j) * w[k] + 2. * e(i, k) * w[j]) * w[i] * w[i] * w[k];
        dpsidwl[k] =
          -30 * (2. * e(i, j) * w[k] + e(i, k) * w[j]) * w[i] * w[i] * w[j];

        // Initialise the second partial  derivatives
        d2psi(i, 0, 0) = 0.0;
        d2psi(i, 0, 1) = 0.0;
        d2psi(i, 0, 2) = 0.0;

        // Second derivatives with respect to wl wll.
        d2psidwl(i, i) = 60 * w[i] * (1 - 3 * w[i] + 2 * w[i] * w[i]) -
                         60 * (e(i, j) * w[k] + e(i, k) * w[j]) * w[j] * w[k];

        d2psidwl(j, i) =
          -60 * (e(i, j) * w[k] + 2 * e(i, k) * w[j]) * w[i] * w[k];

        d2psidwl(k, i) =
          -60 * (2 * e(i, j) * w[k] + e(i, k) * w[j]) * w[i] * w[j];

        d2psidwl(i, j) = d2psidwl(j, i);

        d2psidwl(j, j) = -60 * e(i, k) * w[i] * w[i] * w[k];

        d2psidwl(k, j) = -60 * (e(i, j) * w[k] + e(i, k) * w[j]) * w[i] * w[i];

        d2psidwl(i, k) = d2psidwl(k, i);

        d2psidwl(j, k) = d2psidwl(k, j);

        d2psidwl(k, k) = -60 * e(i, j) * w[i] * w[i] * w[j];

        // Now construct the first derivatives
        for (unsigned alpha = 0; alpha < 2; ++alpha)
        {
          for (unsigned l = 0; l < 3; ++l)
          {
            // add on dpsidwi dwi ds
            dpsi(i, 0, alpha) += dpsidwl[l] * gradw[l][alpha];
          }

          // Only need to fill in 3, enumerate in the easiest manner
          for (unsigned beta = 0; beta < alpha + 1; ++beta)
          {
            // Initialise to the derivative ith respect to alpha
            for (unsigned l = 0; l < 3; ++l)
            {
              for (unsigned ll = 0; ll < 3; ++ll)
              {
                // add on d2psid2wi dwi ds
                d2psi(i, 0, alpha + beta) +=
                  d2psidwl(l, ll) * gradw[l][alpha] * gradw[ll][beta];
              }
            }
          }
        }

        // The 2nd and 3rd shape functions - relating to first derivative dofs
        for (unsigned m = 0; m < 2; ++m)
        {
          unsigned n = (m + 1) % 2; // Cyclic in base two

          // mth components of vectors between nodes, defined for convenience
          const double lik(v[i][m] - v[k][m]), lij(v[i][m] - v[j][m]),
            ljk(v[j][m] - v[k][m]), lkj(v[k][m] - v[j][m]);

          psi(i, 1 + m) =
            pow(w[i], 3) * (4 - 3 * w[i]) * (x[m] - v[i][m]) +
            (15 * lik * e(i, j) + 3 * ljk) * pow(w[i] * w[k], 2) * w[j] +
            (15 * lij * e(i, k) + 3 * lkj) * pow(w[i] * w[j], 2) * w[k];

          // The derivatives with respect to wl
          dpsidwl[i] =
            12 * pow(w[i], 2) * (1 - w[i]) * (x[m] - v[i][m]) +
            (15 * lik * e(i, j) + 3 * ljk) * pow(w[k], 2) * w[j] * w[i] * 2. +
            (15 * lij * e(i, k) + 3 * lkj) * pow(w[j], 2) * w[k] * w[i] * 2.;

          dpsidwl[j] =
            (15 * lik * e(i, j) + 3 * ljk) * pow(w[i] * w[k], 2) +
            (15 * lij * e(i, k) + 3 * lkj) * pow(w[i], 2) * w[k] * w[j] * 2.;

          dpsidwl[k] =
            (15 * lik * e(i, j) + 3 * ljk) * pow(w[i], 2) * w[j] * w[k] * 2. +
            (15 * lij * e(i, k) + 3 * lkj) * pow(w[i] * w[j], 2);

          // Second derivatives with respect to wl wll.
          d2psidwl(i, i) =
            12 * w[i] * (2 - 3 * w[i]) * (x[m] - v[i][m]) +
            (15 * lik * e(i, j) + 3 * ljk) * pow(w[k], 2) * w[j] * 2. +
            (15 * lij * e(i, k) + 3 * lkj) * pow(w[j], 2) * w[k] * 2.;

          d2psidwl(j, i) =
            (15 * lik * e(i, j) + 3 * ljk) * pow(w[k], 2) * w[i] * 2. +
            (15 * lij * e(i, k) + 3 * lkj) * w[j] * w[k] * w[i] * 4.;

          d2psidwl(k, i) =
            (15 * lik * e(i, j) + 3 * ljk) * w[k] * w[j] * w[i] * 4. +
            (15 * lij * e(i, k) + 3 * lkj) * pow(w[j], 2) * w[i] * 2.;

          d2psidwl(i, j) = d2psidwl(j, i);

          d2psidwl(j, j) =
            (15 * lij * e(i, k) + 3 * lkj) * pow(w[i], 2) * w[k] * 2.;

          d2psidwl(k, j) =
            (15 * lik * e(i, j) + 3 * ljk) * pow(w[i], 2) * w[k] * 2. +
            (15 * lij * e(i, k) + 3 * lkj) * pow(w[i], 2) * w[j] * 2.;

          d2psidwl(i, k) = d2psidwl(k, i);

          d2psidwl(j, k) = d2psidwl(k, j);

          d2psidwl(k, k) =
            (15 * lik * e(i, j) + 3 * ljk) * pow(w[i], 2) * w[j] * 2.;

          // The partial derivatives with respect to x_m
          dpsi(i, 1 + m, m) = pow(w[i], 3) * (4 - 3 * w[i]);
          dpsi(i, 1 + m, n) = 0.0; // n != m

          // The partial derivatives with respect to x_m
          d2psi(i, 1 + m, m + n) = 12 * pow(w[i], 2) * (1 - w[i]) * gradw[i][n];
          d2psi(i, 1 + m, m + m) = 24 * pow(w[i], 2) * (1 - w[i]) * gradw[i][m];
          d2psi(i, 1 + m, n + n) = 0.0;

          // Now construct the first derivatives
          for (unsigned alpha = 0; alpha < 2; ++alpha)
          {
            for (unsigned l = 0; l < 3; ++l)
            {
              // add on dpsidwi dwi ds
              dpsi(i, 1 + m, alpha) += dpsidwl[l] * gradw[l][alpha];
            }

            // Only need to fill in 3, enumerate in the easiest manner
            for (unsigned beta = 0; beta < alpha + 1; ++beta)
            {
              for (unsigned l = 0; l < 3; ++l)
              {
                for (unsigned ll = 0; ll < 3; ++ll)
                {
                  // add on d2psid2wi dwi ds
                  d2psi(i, 1 + m, alpha + beta) +=
                    d2psidwl(l, ll) * gradw[l][alpha] * gradw[ll][beta];
                }
              }
            }
          }
        }

        // The 4th and 6th shape functions - relating to the second x and y
        // derivative
        //  dofs
        for (unsigned m = 0; m < 2; ++m)
        {
          unsigned n = (m + 1) % 2; // Cyclic in base two

          // mth components of vectors between nodes, defined for convenience
          const double lik(v[i][m] - v[k][m]), lij(v[i][m] - v[j][m]),
            ljk(v[j][m] - v[k][m]), lkj(v[k][m] - v[j][m]);

          psi(i, 3 + 2 * m) = pow(w[i], 3) * pow(x[m] - v[i][m], 2) / 2. -
                              (5 * lik * e(i, j) + 2 * ljk) * lik *
                                pow(w[i] * w[k], 2) * w[j] / 2. -
                              (5 * lij * e(i, k) + 2 * lkj) * lij *
                                pow(w[i] * w[j], 2) * w[k] / 2.;

          // The derivatives with respect to wl
          Vector<double> dpsidwl(3, 0.0);
          dpsidwl[i] =
            3 * pow(w[i], 2) * pow(x[m] - v[i][m], 2) / 2. -
            (5 * lik * e(i, j) + 2 * ljk) * lik * pow(w[k], 2) * w[j] * w[i] -
            (5 * lij * e(i, k) + 2 * lkj) * lij * pow(w[j], 2) * w[k] * w[i];

          dpsidwl[j] =
            -(5 * lik * e(i, j) + 2 * ljk) * lik * pow(w[i] * w[k], 2) / 2. -
            (5 * lij * e(i, k) + 2 * lkj) * lij * pow(w[i], 2) * w[k] * w[j];

          dpsidwl[k] =
            -(5 * lik * e(i, j) + 2 * ljk) * lik * pow(w[i], 2) * w[j] * w[k] -
            (5 * lij * e(i, k) + 2 * lkj) * lij * pow(w[i] * w[j], 2) / 2.;

          // Second derivatives with respect to wl wll.
          d2psidwl(i, i) =
            6 * w[i] * pow(x[m] - v[i][m], 2) / 2. -
            (5 * lik * e(i, j) + 2 * ljk) * lik * pow(w[k], 2) * w[j] -
            (5 * lij * e(i, k) + 2 * lkj) * lij * pow(w[j], 2) * w[k];

          d2psidwl(j, i) =
            -(5 * lik * e(i, j) + 2 * ljk) * lik * pow(w[k], 2) * w[i] -
            (5 * lij * e(i, k) + 2 * lkj) * lij * w[j] * w[k] * w[i] * 2.;

          d2psidwl(k, i) =
            -(5 * lik * e(i, j) + 2 * ljk) * lik * w[k] * w[j] * w[i] * 2. -
            (5 * lij * e(i, k) + 2 * lkj) * lij * pow(w[j], 2) * w[i];

          d2psidwl(i, j) = d2psidwl(j, i);

          d2psidwl(j, j) =
            -(5 * lij * e(i, k) + 2 * lkj) * lij * pow(w[i], 2) * w[k];

          d2psidwl(k, j) =
            -(5 * lik * e(i, j) + 2 * ljk) * lik * pow(w[i], 2) * w[k] -
            (5 * lij * e(i, k) + 2 * lkj) * lij * pow(w[i], 2) * w[j];

          d2psidwl(i, k) = d2psidwl(k, i);

          d2psidwl(j, k) = d2psidwl(k, j);

          d2psidwl(k, k) =
            -(5 * lik * e(i, j) + 2 * ljk) * lik * pow(w[i], 2) * w[j];

          // The partial derivatives with respect to x_m
          dpsi(i, 3 + 2 * m, m) = pow(w[i], 3) * (x[m] - v[i][m]);
          dpsi(i, 3 + 2 * m, (m + 1) % 2) = 0.0; // n != m

          // Direct partial second derivative when m==n
          d2psi(i, 3 + 2 * m, m + m) = pow(w[i], 3);
          // Second partial derivatives, d2w_dxn_dl dl_dxm
          d2psi(i, 3 + 2 * m, m + m) +=
            6 * pow(w[i], 2) * (x[m] - v[i][m]) * gradw[i][m];
          d2psi(i, 3 + 2 * m, m + n) =
            3 * pow(w[i], 2) * (x[m] - v[i][m]) * gradw[i][n];
          d2psi(i, 3 + 2 * m, n + n) = 0.0;

          // Now construct the first derivatives
          for (unsigned alpha = 0; alpha < 2; ++alpha)
          {
            // Initialise to the derivative ith respect to alpha
            for (unsigned l = 0; l < 3; ++l)
            {
              // add on dpsidwi dwi ds
              dpsi(i, 3 + 2 * m, alpha) += dpsidwl[l] * gradw[l][alpha];
            }
            // Only need to fill in 3, enumerate in the easiest manner
            for (unsigned beta = 0; beta < alpha + 1; ++beta)
            {
              // Initialise to the derivative ith respect to alpha
              for (unsigned l = 0; l < 3; ++l)
              {
                for (unsigned ll = 0; ll < 3; ++ll)
                {
                  // add on d2psid2wi dwi ds
                  d2psi(i, 3 + 2 * m, alpha + beta) +=
                    d2psidwl(l, ll) * gradw[l][alpha] * gradw[ll][beta];
                }
              }
            }
          }
        }

        // The  5th shape functions - relating to the cross derivative dofs
        psi(i, 4) = pow(w[i], 3) * (x[0] - v[i][0]) * (x[1] - v[i][1]);
        // The initialise derivatives with respect to wl
        dpsidwl[i] = 3 * pow(w[i], 2) * (x[0] - v[i][0]) * (x[1] - v[i][1]);
        dpsidwl[j] = 0.0;
        dpsidwl[k] = 0.0;

        // The partial derivatives with respect to x_m
        dpsi(i, 4, 0) = pow(w[i], 3) * (x[1] - v[i][1]);
        dpsi(i, 4, 1) = pow(w[i], 3) * (x[0] - v[i][0]);

        // Initialise second partial derivatives with respect to wl
        d2psidwl = DenseMatrix<double>(3, 3, 0.0);
        d2psidwl(i, i) = 6 * w[i] * (x[0] - v[i][0]) * (x[1] - v[i][1]);

        // The second partial derivatives with respect to x_m x_n
        d2psi(i, 4, 0) = 0.0;
        d2psi(i, 4, 1) = pow(w[i], 3);
        d2psi(i, 4, 2) = 0.0;

        d2psi(i, 4, 0 + 0) += 6 * pow(w[i], 2) * (x[1] - v[i][1]) * gradw[i][0];
        d2psi(i, 4, 0 + 1) += 3 * pow(w[i], 2) * (x[1] - v[i][1]) * gradw[i][1];
        d2psi(i, 4, 1 + 1) += 6 * pow(w[i], 2) * (x[0] - v[i][0]) * gradw[i][1];
        d2psi(i, 4, 1 + 0) += 3 * pow(w[i], 2) * (x[0] - v[i][0]) * gradw[i][0];

        // Vectors between nodes, defined for convenience
        Vector<double> lik(2, 0.0), lij(2, 0.0);

        for (unsigned m = 0; m < 2; ++m)
        {
          unsigned n = (m + 1) % 2; // Cyclic in base two

          // Fill in the values
          lik[m] = (v[i][m] - v[k][m]);
          lij[m] = (v[i][m] - v[j][m]);
          lik[n] = (v[i][n] - v[k][n]);
          lij[n] = (v[i][n] - v[j][n]);

          // Vectors compenent between nodes, defined for convenience
          const double lkj(v[k][m] - v[j][m]), ljk(v[j][m] - v[k][m]);

          psi(i, 4) += -(5 * lik[m] * e(i, j) + 2 * ljk) * lik[n] *
                         pow(w[i] * w[k], 2) * w[j] / 2. -
                       (5 * lij[m] * e(i, k) + 2 * lkj) * lij[n] *
                         pow(w[i] * w[j], 2) * w[k] / 2.;

          // Now for the first derivatives
          dpsidwl[i] += -(5 * lik[m] * e(i, j) + 2 * ljk) * lik[n] *
                          pow(w[k], 2) * w[j] * w[i] -
                        (5 * lij[m] * e(i, k) + 2 * lkj) * lij[n] *
                          pow(w[j], 2) * w[k] * w[i];

          dpsidwl[j] += -(5 * lik[m] * e(i, j) + 2 * ljk) * lik[n] *
                          pow(w[i] * w[k], 2) / 2. -
                        (5 * lij[m] * e(i, k) + 2 * lkj) * lij[n] *
                          pow(w[i], 2) * w[k] * w[j];

          dpsidwl[k] += -(5 * lik[m] * e(i, j) + 2 * ljk) * lik[n] *
                          pow(w[i], 2) * w[j] * w[k] -
                        (5 * lij[m] * e(i, k) + 2 * lkj) * lij[n] *
                          pow(w[i] * w[j], 2) / 2.;

          // Second derivatives with respect to wl wll.
          d2psidwl(i, i) +=
            -(5 * lik[m] * e(i, j) + 2 * ljk) * lik[n] * pow(w[k], 2) * w[j] -
            (5 * lij[m] * e(i, k) + 2 * lkj) * lij[n] * pow(w[j], 2) * w[k];

          d2psidwl(j, i) +=
            -(5 * lik[m] * e(i, j) + 2 * ljk) * lik[n] * pow(w[k], 2) * w[i] -
            (5 * lij[m] * e(i, k) + 2 * lkj) * lij[n] * w[j] * w[k] * w[i] * 2.;

          d2psidwl(k, i) +=
            -(5 * lik[m] * e(i, j) + 2 * ljk) * lik[n] * w[k] * w[j] * w[i] *
              2. -
            (5 * lij[m] * e(i, k) + 2 * lkj) * lij[n] * pow(w[j], 2) * w[i];

          d2psidwl(j, j) +=
            -(5 * lij[m] * e(i, k) + 2 * lkj) * lij[n] * pow(w[i], 2) * w[k];

          d2psidwl(k, j) +=
            -(5 * lik[m] * e(i, j) + 2 * ljk) * lik[n] * pow(w[i], 2) * w[k] -
            (5 * lij[m] * e(i, k) + 2 * lkj) * lij[n] * pow(w[i], 2) * w[j];

          d2psidwl(k, k) +=
            -(5 * lik[m] * e(i, j) + 2 * ljk) * lik[n] * pow(w[i], 2) * w[j];
        }

        // The lower triangular terms
        d2psidwl(i, j) = d2psidwl(j, i);

        d2psidwl(i, k) = d2psidwl(k, i);

        d2psidwl(j, k) = d2psidwl(k, j);

        // Now construct the first derivatives
        for (unsigned alpha = 0; alpha < 2; ++alpha)
        {
          // Initialise to the derivative with respect to alpha
          for (unsigned l = 0; l < 3; ++l)
          {
            // add on dpsidwi dwi ds
            dpsi(i, 4, alpha) += dpsidwl[l] * gradw[l][alpha];
          }

          // Only need to fill in 3, enumerate in the easiest manner
          for (unsigned beta = 0; beta < alpha + 1; ++beta)
          {
            for (unsigned l = 0; l < 3; ++l)
            {
              for (unsigned ll = 0; ll < 3; ++ll)
              {
                // add on d2psid2wi dwi ds
                d2psi(i, 4, alpha + beta) +=
                  d2psidwl(l, ll) * gradw[l][alpha] * gradw[ll][beta];
              }
            }
          }
        }
      }

      return J(0, 0) * J(1, 1) - J(0, 1) * J(1, 0);
    }
  } // namespace MyShape
} // namespace oomph

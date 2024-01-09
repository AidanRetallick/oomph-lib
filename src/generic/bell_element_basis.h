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
#ifndef OOMPH_MY_BELL_ELEMENTS
#define OOMPH_MY_BELL_ELEMENTS

// oomph-lib headers
#include "../generic/Vector.h"
#include "../generic/shape.h"

namespace oomph
{
  namespace MyShape
  {
    //===========================================================================
    /// Two by two specialisation of function to calculate inverse of a matrix
    //===========================================================================
    inline double invert_two_by_two(const DenseMatrix<double>& jacobian,
                                    DenseMatrix<double>& inverse_jacobian)
    {
      // Calculate the determinant of the matrix
      const double det =
        jacobian(0, 0) * jacobian(1, 1) - jacobian(0, 1) * jacobian(1, 0);

// Report if Matrix is singular or negative
#ifdef PARANOID
      if (fabs(det) < 1e-12)
      {
        std::stringstream error_stream;
        error_stream
          << "The matrix is singular to machine precision : det(M) = " << det
          << ".\n";
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Calculate the inverse of the 2x2 matrix
      inverse_jacobian(0, 0) = jacobian(1, 1) / det;
      inverse_jacobian(0, 1) = -jacobian(0, 1) / det;
      inverse_jacobian(1, 0) = -jacobian(1, 0) / det;
      inverse_jacobian(1, 1) = jacobian(0, 0) / det;

      return det;
    }

    //=============================================================================
    /// Three-by three specialisation of function to calculate inverse of a
    /// matrix
    //=============================================================================
    inline void invert_three_by_three(const DenseMatrix<double>& jacobian,
                                      DenseMatrix<double>& inverse_jacobian)
    {
      // Calculate the determinant of the matrix
      const double det = jacobian(0, 0) * jacobian(1, 1) * jacobian(2, 2) +
                         jacobian(0, 1) * jacobian(1, 2) * jacobian(2, 0) +
                         jacobian(0, 2) * jacobian(1, 0) * jacobian(2, 1) -
                         jacobian(0, 0) * jacobian(1, 2) * jacobian(2, 1) -
                         jacobian(0, 1) * jacobian(1, 0) * jacobian(2, 2) -
                         jacobian(0, 2) * jacobian(1, 1) * jacobian(2, 0);

      // Report if Matrix is singular or negative
#ifdef PARANOID
      if (fabs(det) < 1e-12)
      {
        std::stringstream error_stream;
        error_stream
          << "The matrix is singular to machine precision : det(M) = " << det
          << ".\n";
        throw OomphLibError(
          error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
      }
#endif

      // Calculate the inverse of the 3x3 matrix
      inverse_jacobian(0, 0) =
        (jacobian(1, 1) * jacobian(2, 2) - jacobian(1, 2) * jacobian(2, 1)) /
        det;
      inverse_jacobian(0, 1) =
        -(jacobian(0, 1) * jacobian(2, 2) - jacobian(0, 2) * jacobian(2, 1)) /
        det;
      inverse_jacobian(0, 2) =
        (jacobian(0, 1) * jacobian(1, 2) - jacobian(0, 2) * jacobian(1, 1)) /
        det;
      inverse_jacobian(1, 0) =
        -(jacobian(1, 0) * jacobian(2, 2) - jacobian(1, 2) * jacobian(2, 0)) /
        det;
      inverse_jacobian(1, 1) =
        (jacobian(0, 0) * jacobian(2, 2) - jacobian(0, 2) * jacobian(2, 0)) /
        det;
      inverse_jacobian(1, 2) =
        -(jacobian(0, 0) * jacobian(1, 2) - jacobian(0, 2) * jacobian(1, 0)) /
        det;
      inverse_jacobian(2, 0) =
        (jacobian(1, 0) * jacobian(2, 1) - jacobian(1, 1) * jacobian(2, 0)) /
        det;
      inverse_jacobian(2, 1) =
        -(jacobian(0, 0) * jacobian(2, 1) - jacobian(0, 1) * jacobian(2, 0)) /
        det;
      inverse_jacobian(2, 2) =
        (jacobian(0, 0) * jacobian(1, 1) - jacobian(0, 1) * jacobian(1, 0)) /
        det;
    }

    /// Class to contain the Basis polynomials for the Bell element
    class BellElementBasis
    {
    public:
      /// Constructor
      BellElementBasis() {}

      /// Destructor
      ~BellElementBasis() {}

      // Private copy and assign - so this should cause a compilation error
    private:
      /// Broken copy constructor
      BellElementBasis(BellElementBasis& dummy)
      {
        BrokenCopy::broken_copy("BellElementBasis");
      }
      /// Broken assignment operator
      void operator=(const BellElementBasis&)
      {
        BrokenCopy::broken_assign("BellElementBasis");
      }

    public:
      /// Get the (twice) area of the triangle from the vertices
      double get_twice_triangle_area(const Vector<double>& v0,
                                     const Vector<double>& v1,
                                     const Vector<double>& v2) const;

      /// Get outer normal of side between vertices v0 and v1, assumes
      /// counter-clockwise triangle vertices.
      Vector<double> get_outer_normal(const Vector<double>& v0,
                                      const Vector<double>& v1) const;

      /// Basis for a Bell element. This follows exactly the notation of
      /// M. Okabe in Comput. Methods Appl. Mech. 117 (1994) 411-421
      void d2_basis(const Vector<double>& s,
                    const Vector<Vector<double>>& v,
                    Shape& psi,
                    DShape& dpsi,
                    DShape& d2psi) const;

      /// Explicit basis for a Bell element. This follows exactly the notation
      /// of M. Okabe in Comput. Methods Appl. Mech. 117 (1994) 411-421
      double d2_basis_eulerian(const Vector<double>& s,
                               const Vector<Vector<double>>& v,
                               Shape& psi,
                               DShape& dpsi,
                               DShape& d2psi) const;
    }; // End

  } // namespace MyShape
} // namespace oomph

#endif

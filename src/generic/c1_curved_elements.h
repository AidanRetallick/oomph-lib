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
#ifndef OOMPH_MY_CURVED_C1_ELEMENTS
#define OOMPH_MY_CURVED_C1_ELEMENTS

// oomph-lib headers
#include "Vector.h"
#include "shape.h"
#include "bell_element_basis.h"
#include "my_geom_object.h"

namespace oomph
{
  namespace MyC1CurvedElements
  {
    // HERE replace with class enum in c++11?
    /// enum to enumerate the possible edges that could be curved
    enum Edge
    {
      none = -1,
      zero = 0,
      one = 1,
      two = 2
    };


    // [zdec] I see no point in having this base class?
    class BernadouElementBasisBase
    {
    public:
      BernadouElementBasisBase(){};
      virtual ~BernadouElementBasisBase(){};

      /// Shorthand for a vector of vectors containining the vertices
      typedef Vector<Vector<double>> VertexList;

      /// Get the basis for the  unknowns
      virtual void shape(const Vector<double>& s,
                         Shape& nodal_basis,
                         Shape& bubble_basis) const = 0;

      /// Get the Eulerian first derivatives of the basis functions
      virtual double d_shape_dx(const Vector<double>& s,
                                Shape& psi,
                                Shape& bpsi,
                                DShape& dpsi,
                                DShape& dbpsi) const = 0;

      /// Get the Eulerian second derivatives of the basis functions
      virtual double d2_shape_dx2(const Vector<double>& s,
                                  Shape& psi,
                                  Shape& bpsi,
                                  DShape& dpsi,
                                  DShape& dbpsi,
                                  DShape& d2psi,
                                  DShape& d2bpsi) const = 0;

      /// Get the Eulerian second derivatives of the basis functions, returning
      /// the association matrix used to assemble the basis.
      virtual double d2_shape_dx2(const Vector<double>& s,
                                  Shape& psi,
                                  Shape& bpsi,
                                  DShape& dpsi,
                                  DShape& dbpsi,
                                  DShape& d2psi,
                                  DShape& d2bpsi,
                                  const DenseMatrix<double>& m) const = 0;

      virtual void coordinate_x(const Vector<double>& s,
                                Vector<double>& fk) const = 0;

      virtual void get_jacobian(const Vector<double>& s,
                                DenseMatrix<double>& jacobian) const = 0;

      virtual void upgrade_element(
        const VertexList& verts,
        const double& su,
        const double& so,
        const Edge& curved_edge,
        const CurvilineGeomObject& parametric_curve) = 0;

      /// Check the element
      virtual void self_check(
        const CurvilineGeomObject& parametric_curve) const = 0;

      /// Return the number of basis functions on the physical triangle
      virtual unsigned n_basis_functions() const = 0;

      /// Return the number of bubble basis functions triangle
      virtual unsigned n_internal_dofs() const = 0;

      /// Return the number of basis functions on the basic triangle
      virtual unsigned n_basic_basis_functions() const = 0;

      /// Fill in the association matrix between monomials and unknowns
      virtual void fill_in_full_association_matrix(
        DenseMatrix<double>& conversion_matrix) const = 0;

      /// Get the reference location of the internal dofs for a permuted shape
      virtual void get_internal_dofs_location(
        const unsigned& idof, Vector<double>& s_permute) const = 0;
    };

    /// The BernadouElementBasis class.
    /// These are based on the curved C1 elements of Bernadou and Boisserie 1994
    /// but adapted for use with the Bell elements rather than Argyris. The
    /// structure of these elements is as follows: Fk maps a point on the
    /// reference element to a point on the curved element. It is determined by
    /// the vertex nodal positions (ai) and the prescribed curved boundary
    /// chi(s).
    ///
    /* The Structure */
    /*      @                 @                      @                        */
    /*     /(       Dij      /(          Mij         |\                       */
    /*    /. \      ->      /. \         ->          |. \                     */
    /*   /. . )            /. . )                    |. . \                   */
    /*  @-----@           @-----@                    @-----@                  */
    /*                                                                        */
    /*  Physical (21):      Reference dofs (21):     basic dofs(36): */
    /*   w(ai)               w(ai)                   w(ahati) */
    /*   w,j(ai)             w,tij (ai)              w,j(ahati) */
    /*   w,jk(ai)            w,tij tij(ai)           w,jk(ahati) */
    /*   w(ei)               w,tkj tjk(ai)           w,n(bi) */
    /*                       w(ei)                   w(ei) */
    /*                                               w(di) */
    /*                                               w,n(di) */
    /* where tij are the two tangents at node ai and tjk are the two tangents
     * opposite ai. ai are nodes, bi are midside points and di are mid-midside
     * points*/
    /* ei are internal dofs at points Fk(ei_hat) where ei hat are at (1/4,1/4)
     */
    /* (1/4,1/2) and (1/2,1/4). */
    /// We first map the dofs from global to local dofs such that the
    /// derivatives are now all with respect to local tangents. This is just
    /// done using local rotations of the derivative dofs using a 24x24 matrix
    /// that we refer to as d_matrix. We then map the 21 local dofs onto the 36
    /// basic 'dofs': which aren't really dofs. The mapping is curved so a line
    /// on the basic element will map to a curve on the physical space. For this
    /// reason we need a higher order representation in the basic space: with
    /// the constraint that along the boundaries the traces and normal
    /// derivatives are compatible with the Bell element.
    ///
    /// This means on each boundary the trace must be a fifth order polynomial
    /// of the tangent coordinate, which is defined by the degrees of freedom
    /// shared by the Bell element i.e the dofs: w w,i and w,ij. Similarly the
    /// trace of w,n must be the third order polynomial defined by w,i and w,ij.
    /// The solution is then represented on the curved element as a seventh
    /// order bivariate polynomial that is constrained to have the
    /// aforementioned traces.
    template<unsigned BOUNDARY_ORDER>
    class BernadouElementBasis : public BernadouElementBasisBase
    {
    public:
      /// \short typedef for the edge curve
      typedef void (*ParametricCurveFctPt)(const double& s,
                                           Vector<double>& param_fct);

      /// \short Shorthand for a vector of vectors containining the vertices
      typedef Vector<Vector<double>> VertexList;

      /// Default Constructor
      BernadouElementBasis() : Curved_edge(none) {}

      /// Constructor that takes vertices and start and end parts as arguments.
      void upgrade_element(const VertexList& verts,
                           const double& su,
                           const double& so,
                           const Edge& curved_edge,
                           const CurvilineGeomObject& parametric_curve)
      {
        // Store vertices
        Vertices = verts;
        // Set up the new curved data for the element
        S_ubar = su;
        S_obar = so;
        Curved_edge = curved_edge;
        /// Fill in the function values at vertex 0
        Chi_subar.resize(2);
        D_chi_subar.resize(2);
        D2_chi_subar.resize(2);
        parametric_curve.position(Vector<double>(1, su), Chi_subar);
        parametric_curve.dposition(Vector<double>(1, su), D_chi_subar);
        parametric_curve.d2position(Vector<double>(1, su), D2_chi_subar);
        /// Fill in the function values at vertex 1
        Chi_sobar.resize(2);
        D_chi_sobar.resize(2);
        D2_chi_sobar.resize(2);
        parametric_curve.position(Vector<double>(1, so), Chi_sobar);
        parametric_curve.dposition(Vector<double>(1, so), D_chi_sobar);
        parametric_curve.d2position(Vector<double>(1, so), D2_chi_sobar);
// Check the construction of the elements is complete
#ifdef PARANOID
        self_check(parametric_curve);
#endif
      }


      /// Destructor
      ~BernadouElementBasis() {}

      // Private copy and assign - so this should cause a compilation error
    private:
      /// Broken copy constructor
      BernadouElementBasis(BernadouElementBasis& dummy)
      {
        BrokenCopy::broken_copy("BernadouElementBasis");
      }
      /// Broken assignment operator
      void operator=(const BernadouElementBasis&)
      {
        BrokenCopy::broken_assign("BernadouElementBasis");
      }

    public:
      /// Check the element
      inline void self_check(const CurvilineGeomObject& parametric_curve) const;

      /// Get the physical coordinate
      void coordinate_x(const Vector<double>& s, Vector<double>& fk) const;

      /// Get the vertices const version
      inline VertexList get_vertices() const
      {
        return Vertices;
      }

      /// Get the values of s at start of parametric curve section
      inline const double& get_s_ubar() const
      {
        // If we have upgraded
        if (Curved_edge != none)
        {
          return S_ubar;
        }
        else
        {
          throw OomphLibError("The element has not been upgraded yet. Did \
  you forget to set a Curved_edge?",
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }

      /// Get the values of s at end of parametric curve section
      inline const double& get_s_obar() const
      {
        // If we have upgraded
        if (Curved_edge != none)
        {
          return S_obar;
        }
        else
        {
          throw OomphLibError("The element has not been upgraded yet. Did \
  you forget to set a Curved_edge?",
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }


      /// Get the position at start of parametric curve section
      inline const Vector<double>& get_chi_subar() const
      {
        // If we have upgraded
        if (Curved_edge != none)
        {
          return Chi_subar;
        }
        else
        {
          throw OomphLibError("The element has not been upgraded yet. Did \
  you forget to set a Curved_edge?",
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }


      /// Get the position at start of parametric curve section
      inline const Vector<double>& get_chi_sobar() const
      {
        // If we have upgraded
        if (Curved_edge != none)
        {
          return Chi_sobar;
        }
        else
        {
          throw OomphLibError("The element has not been upgraded yet. Did \
  you forget to set a Curved_edge?",
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }


      /// Get the position at start of parametric curve section
      inline const Vector<double>& get_d_chi_subar() const
      {
        // If we have upgraded
        if (Curved_edge != none)
        {
          return D_chi_subar;
        }
        else
        {
          throw OomphLibError("The element has not been upgraded yet. Did \
  you forget to set a Curved_edge?",
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }


      /// Get the position at start of parametric curve section
      inline const Vector<double>& get_d_chi_sobar() const
      {
        // If we have upgraded
        if (Curved_edge != none)
        {
          return D_chi_sobar;
        }
        else
        {
          throw OomphLibError("The element has not been upgraded yet. Did \
  you forget to set a Curved_edge?",
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }


      /// Get the position at start of parametric curve section
      inline const Vector<double>& get_d2_chi_subar() const
      {
        // If we have upgraded
        if (Curved_edge != none)
        {
          return D2_chi_subar;
        }
        else
        {
          throw OomphLibError("The element has not been upgraded yet. Did \
  you forget to set a Curved_edge?",
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }


      /// Get the position at start of parametric curve section
      inline const Vector<double>& get_d2_chi_sobar() const
      {
        // If we have upgraded
        if (Curved_edge != none)
        {
          return D2_chi_sobar;
        }
        else
        {
          throw OomphLibError("The element has not been upgraded yet. Did \
  you forget to set a Curved_edge?",
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
      }


      /// The approximated (3rd order) polynomial
      void psi_h(const double& s1, Vector<double>& psi_h) const;

      /// Fill in the full association matrix
      void fill_in_full_association_matrix(
        DenseMatrix<double>& conversion_matrix) const;

      /// Return the order of the polynomial on the curved boudnary
      inline unsigned boundary_order() const
      {
        return BOUNDARY_ORDER;
      }

      /// Return the order of the full basis on the basic triangle
      inline unsigned basic_basis_order() const
      {
        return BOUNDARY_ORDER + 4;
      }

      /// Return the number of basis functions on the physical triangle
      inline unsigned n_basic_basis_functions() const
      {
        return (BOUNDARY_ORDER + 5) * (BOUNDARY_ORDER + 6) / 2;
      }

      /// Return the number of basis functions on the physical triangle
      inline unsigned n_basis_functions() const;

      /// Return the number of basis functions on the physical triangle
      inline unsigned n_internal_dofs() const
      {
        return n_basis_functions() - 18;
      }

      /// \short Return the number of linearly dependent midside basic nodes
      /// (not including Argyris dofs that are also eliminated). Each has two
      /// dofs - a normal and a functional dof.
      inline unsigned n_basic_midside_nodes() const
      {
        return (n_basic_basis_functions() - n_internal_dofs() - 21) / 2;
      }

      // Data is private
    private:
      /// The vertices supplied from the element
      VertexList Vertices;

      /// Parametric coordinate at vertex 0 (assuming 2 is always curved edge)
      double S_ubar;

      /// Parametric coordinate at vertex 1 (assuming 2 is always curved edge)
      double S_obar;

      /// The function evaluate at vertex 0 (assuming 2 is always curved edge)
      Vector<double> Chi_subar;

      /// The function evaluate at vertex 1 (assuming 2 is always curved edge)
      Vector<double> Chi_sobar;

      /// The derivative evaluated at vertex 0 (assuming 2 is always curved
      /// edge)
      Vector<double> D_chi_subar;

      /// The derivative evaluated at vertex 1 (assuming 2 is always curved
      /// edge)
      Vector<double> D_chi_sobar;

      /// The 2nd derivative evaluated at vertex 0 (assuming 2 is always curved
      /// edge)
      Vector<double> D2_chi_subar;

      /// The 2nd derivative evaluated at vertex 1 (assuming 2 is always curved
      /// edge)
      Vector<double> D2_chi_sobar;

      /// Whih edge is curved
      Edge Curved_edge;

    protected:
      /*  Protected member functions: */
      /* These functions are used in the construction of shape - but not
       * intended */
      /* for use at the user end */

      /// The mapping F_k - a polynomial degree 3 PRIVATE
      void f_k(const Vector<double>& s, Vector<double>& fk) const;

      /// The basic Jacobian PRIVATE
      void get_basic_jacobian(const Vector<double> s,
                              DenseMatrix<double>& jac) const;

      /// \short The Hessian of the global coordinate (of the vector mapping) -
      /// like a second order Jacobian.
      /*        d^2 x_i
          or:   ----------    (rank 3) with x the global coordinate and s the
         local. d s_i ds_j */
      // PRIVATE
      void get_basic_hessian(const Vector<double>& s,
                             RankThreeTensor<double>& hess) const;


      /// Get the edge permutation, without the index shift
      inline void permute_shape(Vector<double>& s) const;

      /// Get the edge permutation
      inline void get_jacobian_of_permute(DenseMatrix<double>& jac) const;

      /// Get the edge permutation
      inline void nodal_index_shift(unsigned& index_shift) const;

      // These Vectors are used repeatedly in the construction of the shape
      // functions So are inlined.

      /// \short Components of the first of the two tangent vectors at node 0
      /// Vector version  (labelling Ai i in {1,2} anticlockwise)
      inline double A1(const unsigned& i) const
      {
        return Vertices[2][i] - Vertices[0][i];
      }

      /// \short void version filling in the first of the two tangent vectors at
      /// node
      ///  0 Vector version  (labelling Ai i in {1,2} anticlockwise)
      inline void A1(Vector<double>& v) const
      {
        for (unsigned i = 0; i < 2; ++i)
        {
          v[i] = Vertices[2][i] - Vertices[0][i];
        }
      }

      /// \short Components of the second of the two tangent vectors at node 0
      /// Vector version  (labelling Ai i in {1,2} anticlockwise)
      inline double A2(const unsigned& i) const
      {
        return (S_obar - S_ubar) * D_chi_subar[i];
      }

      /// \short void version filling in the second of the two tangent vectors
      /// at node 0 - vector version  (labelling Ai i in {1,2} anticlockwise)
      inline void A2(Vector<double>& v) const
      {
        v[0] = A2(0);
        v[1] = A2(1);
      }

      /// \short Components of the first tangent vector at node 1 and (labelling
      /// Bi i in {1,2} anticlockwise)
      inline double B1(const unsigned& i) const
      {
        return -(S_obar - S_ubar) * D_chi_sobar[i];
      }

      /// \short Fill in first tangent vector at node 1 and (labelling Bi i in
      /// {1,2} anticlockwise)
      inline void B1(Vector<double>& v) const
      {
        v[0] = B1(0);
        v[1] = B1(1);
      }

      /// \short Components of the second tangent vector at node 1 and
      /// (labelling Bi i in {1,2} anticlockwise)
      inline double B2(const unsigned& i) const
      {
        return (Vertices[2][i] - Vertices[1][i]);
      }

      /// Fill in second tangent vector at node 1 and (labelling Bi i in
      /// {1,2} anticlockwise)
      inline void B2(Vector<double>& v) const
      {
        for (unsigned i = 0; i < 2; ++i)
        {
          v[i] = Vertices[2][i] - Vertices[1][i];
        }
      }

      /// The vectors of d2_chi defined at node 0
      inline void D1(Vector<double>& v) const
      {
        v[0] = pow(S_obar - S_ubar, 2) * D2_chi_subar[0];
        v[1] = pow(S_obar - S_ubar, 2) * D2_chi_subar[1];
      }

      /// The vectors of d2_chi defined at node 1
      inline void D2(Vector<double>& v) const
      {
        v[0] = pow(S_obar - S_ubar, 2) * D2_chi_sobar[0];
        v[1] = pow(S_obar - S_ubar, 2) * D2_chi_sobar[1];
      }


      /// Enumerated type that represents the three potential points
      /// along
      ///  a side of the triangle
      enum S_basic_node
      {
        one_quarter = 0,
        one_half = 1,
        three_quarters = 2
      };

      /// These constants are needed in the construction of the Mij matrix to
      /// transform onto the basic element

      /// \short(Signed) area of a parallelogram (or alternatively z component
      /// of cross product of two vectors on x-y plane ). Takes two Vectors as
      /// argument.
      inline double parallelogram_area(const Vector<double>& v0,
                                       const Vector<double>& v1) const
      {
        return v0[0] * v1[1] - v0[1] * v1[0];
      }

      /// \short (Signed) area of a parallelogram (or alternatively z component
      /// of cross product of two vectors on x-y plane ). Takes four compenents
      /// as
      //  arguments.
      inline double parallelogram_area(const double& v0x,
                                       const double& v0y,
                                       const double& v1x,
                                       const double& v1y) const
      {
        return v0x * v1y - v0y * v1x;
      }

      // Define some useful constants
      /// Components of a0 - a1 in direction A1
      inline double a_tilde_1() const
      {
        Vector<double> b2(2, 0), a2(2, 0), a1(2, 0);
        B2(b2);
        A2(a2);
        A1(a1);
        return parallelogram_area(b2, a2) / parallelogram_area(a1, a2) - 1;
      }

      /// Components of a0 - a1 in direction A2
      inline double a_tilde_2() const
      {
        Vector<double> b2(2, 0), a2(2, 0), a1(2, 0);
        B2(b2);
        A2(a2);
        A1(a1);
        return parallelogram_area(a1, b2) / parallelogram_area(a1, a2);
      }

      /// Components of -B1 in directions A1
      inline double a_tildetilde_1() const
      {
        Vector<double> a2(2, 0), b1(2, 0), a1(2, 0);
        A2(a2);
        A1(a1);
        B1(b1);
        return -parallelogram_area(b1, a2) / parallelogram_area(a1, a2);
      }

      /// Components of -B1 in directions A2
      inline double a_tildetilde_2() const
      {
        Vector<double> a2(2, 0), b1(2, 0), a1(2, 0);
        A2(a2);
        A1(a1);
        B1(b1);
        return -parallelogram_area(a1, b1) / parallelogram_area(a1, a2);
      }

      /// Components of D1 in directions A1
      inline double a_utilde_1() const
      {
        Vector<double> a2(2, 0), d1(2, 0), a1(2, 0);
        A2(a2);
        A1(a1);
        D1(d1);
        return parallelogram_area(d1, a2) / parallelogram_area(a1, a2);
      }

      /// Components of D1 in directions A2
      inline double a_utilde_2() const
      {
        Vector<double> a2(2, 0), d1(2, 0), a1(2, 0);
        A2(a2);
        A1(a1);
        D1(d1);
        return parallelogram_area(a1, d1) / parallelogram_area(a1, a2);
      }

      /// Components of a1 - a2 in directions B1
      inline double b_tilde_1() const
      {
        Vector<double> b2(2, 0), b1(2, 0), a1(2, 0);
        B2(b2);
        A1(a1);
        B1(b1);
        return parallelogram_area(a1, b2) / parallelogram_area(b1, b2);
      }

      /// Components of a1 - a2 in directions B2
      inline double b_tilde_2() const
      {
        Vector<double> b2(2, 0), b1(2, 0), a1(2, 0);
        B2(b2);
        A1(a1);
        B1(b1);
        return parallelogram_area(b1, a1) / parallelogram_area(b1, b2) - 1;
      }

      /// Components of A2 in directions B1
      inline double b_tildetilde_1() const
      {
        Vector<double> b2(2, 0), b1(2, 0), a2(2, 0);
        A2(a2);
        B2(b2);
        B1(b1);
        return parallelogram_area(a2, b2) / parallelogram_area(b1, b2);
      }

      /// Components of A2 in directions B2
      inline double b_tildetilde_2() const
      {
        Vector<double> b2(2, 0), b1(2, 0), a2(2, 0);
        A2(a2);
        B2(b2);
        B1(b1);
        return parallelogram_area(b1, a2) / parallelogram_area(b1, b2);
      }

      /// Components of D2 in directions B1
      inline double b_utilde_1() const
      {
        Vector<double> b2(2, 0), b1(2, 0), d2(2, 0);
        D2(d2);
        B2(b2);
        B1(b1);
        return parallelogram_area(d2, b2) / parallelogram_area(b1, b2);
      }

      /// Components of D2 in directions B2
      inline double b_utilde_2() const
      {
        Vector<double> b2(2, 0), b1(2, 0), d2(2, 0);
        D2(d2);
        B2(b2);
        B1(b1);
        return parallelogram_area(b1, d2) / parallelogram_area(b1, b2);
      }

      /// Components of A2 in directions C1
      inline double c_tilde_1() const
      {
        Vector<double> b2(2, 0), a1(2, 0), a2(2, 0);
        A2(a2);
        A1(a1);
        B2(b2);
        return -parallelogram_area(a1, a2) / parallelogram_area(a1, b2);
      }

      /// Components of A2 in directions C2
      inline double c_tilde_2() const
      {
        Vector<double> b2(2, 0), a1(2, 0), a2(2, 0);
        A2(a2);
        A1(a1);
        B2(b2);
        return -parallelogram_area(a2, b2) / parallelogram_area(a1, b2);
      }

      /// Components of B1 in directions C1
      inline double c_tildetilde_1() const
      {
        Vector<double> b2(2, 0), a1(2, 0), b1(2, 0);
        B1(b1);
        A1(a1);
        B2(b2);
        return -parallelogram_area(a1, b1) / parallelogram_area(a1, b2);
      }

      /// Components of B1 in directions C2
      inline double c_tildetilde_2() const
      {
        Vector<double> b2(2, 0), a1(2, 0), b1(2, 0);
        B1(b1);
        A1(a1);
        B2(b2);
        return -parallelogram_area(b1, b2) / parallelogram_area(a1, b2);
      }

      // The altitude vectors: the vector parallel to the inner
      // normal that points from the point ci on the ith side to the node ai.
      /// (Signed) altitude 1.
      inline double altitude_1() const
      {
        return parallelogram_area(A1(0), A1(1), B2(0), B2(1)) /
               sqrt(B2(0) * B2(0) + B2(1) * B2(1));
      }

      /// (Signed) altitude 2
      inline double altitude_2() const
      {
        return parallelogram_area(A1(0), A1(1), B2(0), B2(1)) /
               sqrt(A1(0) * A1(0) + A1(1) * A1(1));
      }

      /// \short Altitude vector 1 (i.e normal vector with length of altitude)
      ///  component i: h1 = (A1 . B2/|B2|) B2/|B2| - A1
      inline double altitude_vector_1(unsigned i) const
      {
        double B2B2(B2(0) * B2(0) + B2(1) * B2(1)),
          B2A1(B2(0) * A1(0) + B2(1) * A1(1));
        return -A1(i) + B2A1 / B2B2 * B2(i);
      }

      /// \short Altitude vector 2 (i.e normal vector with length of altitude)
      /// component i: h2 = (B2 . A1/|A1|) A1/|A1| - B2
      inline double altitude_vector_2(unsigned i) const
      {
        double B2A1(B2(0) * A1(0) + B2(1) * A1(1)),
          A1A1(A1(0) * A1(0) + A1(1) * A1(1));
        return -B2(i) + B2A1 / A1A1 * A1(i);
      }

      /// \short Eccentricity Parameters
      /// eta_1 = 1 -  2 A1.B2 / B2.B2
      inline double eta_1() const
      {
        double B2B2(B2(0) * B2(0) + B2(1) * B2(1)),
          B2A1(B2(0) * A1(0) + B2(1) * A1(1));
        return 1 - 2 * B2A1 / B2B2;
      }

      /// \short Eccentricity Parameters
      /// eta_2 =-1 +  2 A1.B2 / A1.A1
      inline double eta_2() const
      {
        double A1B2(A1(0) * B2(0) + A1(1) * B2(1)),
          A1A1(A1(0) * A1(0) + A1(1) * A1(1));
        return 2 * A1B2 / A1A1 - 1;
      }

      // We need 1d shape functions for the trace of the function and its'
      // normal derivative

      /// \short Get 5th order, 1D hermite shape functions
      /// Dofs: w(0) w(1) w'(0) w'(1) w''(0) w''(1)
      void hermite_shape_1d_5(const double& s, Shape& psi) const;

      /// \short Get 5th order, 1D hermite shape functions at basic nodes
      /// Dofs: w(0) w(1) w'(0) w'(1) w''(0) w''(1)
      void hermite_shape_1d_5(const S_basic_node& s, Shape& psi) const;

      /// \short The local derivative of 5th order, 1D hermite shape functions
      /// Dofs: w(0) w(1) w'(0) w'(1) w''(0) w''(1)
      void d_hermite_shape_1d_5(const double& s, DShape& dpsi) const;

      /// \short The local derivative of 5th order, 1D hermite shape functions
      /// Dofs: w(0) w(1) w'(0) w'(1) w''(0) w''(1) at basic_node
      void d_hermite_shape_1d_5(const S_basic_node& s, DShape& dpsi) const;

      /// \short Get 3rd order, 1D hermite shape functions
      /// Dofs: w(0) w(1) w'(0) w'(1)
      void hermite_shape_1d_3(const double& s, Shape& psi) const;

      /// \short Get 3rd order, 1D hermite shape functions
      /// Dofs: w(0) w(1) w'(0) w'(1) at basic node
      void hermite_shape_1d_3(const S_basic_node& s, Shape& psi) const;

      // Now define the w trace column vectors fi i in {1,2,3}
      // These are effectively column vectors that, when dotted with the local
      // dofs, give the trace along boundary i. They are effectively a set of
      // shape functions for the trace but padded with zeros so they can be used
      // in matrix multiplication.
      /// Padded shape functions for trace on side 1
      Vector<double> f_1(const double& s0) const;
      /// Padded shape functions for trace on side 1  at basic nodes
      Vector<double> f_1(const S_basic_node& s0) const;
      /// Local derivative of padded shape functions for trace on side 1
      Vector<double> df_1_ds(const double& s0) const;
      /// \short Local derivative of padded shape functions for trace on side 1
      /// at basic nodes
      Vector<double> df_1_ds(const S_basic_node& s0) const;

      /// Padded shape functions for trace on side 2
      Vector<double> f_2(const double& s1) const;
      /// Padded shape functions for trace on side 2  at basic nodes
      Vector<double> f_2(const S_basic_node& s1) const;
      /// Local derivative of padded shape functions for trace on side 2
      Vector<double> df_2_ds(const double& s1) const;
      /// \short Local derivative of padded shape functions for trace on side 2
      /// at basic nodes
      Vector<double> df_2_ds(const S_basic_node& s1) const;

      /// Padded shape functions for trace on (curved) side 3
      Vector<double> f_3(const double& s0) const;
      /// Padded shape functions for trace on (curved) side 3 at basic nodes
      Vector<double> f_3(const S_basic_node& s0) const;
      /// Local derivative of padded shape functions for trace on (curved) side
      /// 3
      Vector<double> df_3_ds(const double& s0) const;
      /// Local derivative of padded shape functions for trace on (curved) side
      /// 3 at basic nodes
      Vector<double> df_3_ds(const S_basic_node& s0) const;

      // Now define the  w,n trace column vectors gi i in {1,2,3}
      /// Padded shape functions for normal derivative trace on side 1
      Vector<double> g_1(const double& s0) const;
      /// Padded shape functions for normal derivative trace on side 1
      Vector<double> g_1(const S_basic_node& s0) const;

      /// Padded shape functions for normal derivative trace on side 2
      Vector<double> g_2(const double& s1) const;
      /// Padded shape functions for normal derivative trace on side 2
      Vector<double> g_2(const S_basic_node& s1) const;

      /// Padded shape functions for normal derivative trace on (curved) side 3
      Vector<double> g_3(const double& s0) const;
      /// Padded shape functions for normal derivative trace on (curved) side 3
      Vector<double> g_3(const S_basic_node& s0) const;

      /// Fill in matrix that transforms the global dofs to the local dofs
      void local_to_global_matrix(DenseMatrix<double>& l2g) const;

      /// \short Fill in  matrix that transforms between the basic dofs (36) and
      /// the local dofs (21)
      void basic_to_local_matrix(DenseMatrix<double>& b2l) const;

      /// Fill in submatrix used to construct the full local to basic matrix
      void basic_to_local_submatrix_1(DenseMatrix<double>& M1) const;
      /// Fill in submatrix used to construct the full local to basic matrix
      void basic_to_local_submatrix_2(DenseMatrix<double>& M2) const;
      /// Fill in submatrix used to construct the full local to basic matrix
      void basic_to_local_submatrix_3(DenseMatrix<double>& M3) const;
      /// Fill in submatrix used to construct the full local to basic matrix
      void basic_to_local_submatrix_4(DenseMatrix<double>& M4) const;
      /// Fill in submatrix used to construct the full local to basic matrix
      void basic_to_local_submatrix_5(DenseMatrix<double>& M5) const;
      /// Fill in submatrix used to construct the full local to basic matrix
      void basic_to_local_submatrix_6(DenseMatrix<double>& M6) const;
      /// Fill in submatrix used to construct the full local to basic matrix
      void basic_to_local_submatrix_7(DenseMatrix<double>& M7) const;

      /// \short Fill in the matrix that transforms from the 36(55) basis
      /// monomials of p7(9)  to the 36(55) shape functions on the basic
      /// triangle. It was derived in mathematica and is extremely large.
      void monomial_to_basic_matrix(DenseMatrix<double>& m) const;

      /// \short Fill in the inverse matrix that transforms from the 36(55)
      /// basis monomials of p7(9) (may be accurate/similar time to invert on
      /// the fly) to the 36(55) shape functions on the basic triangle.
      void inverse_monomial_to_basic_matrix(DenseDoubleMatrix& m) const;

      /// \short Get full basis for a generic (BOUNDARY_ORDER+4)th order
      /// bivariate polynomial: i.e 36 basis monomials for generic p7 polynomial
      /// ;   55 basis monomials for generic p9 polynomial.
      void full_basis_monomials(const Vector<double>& s, Shape& pn) const;

      /// Get full basis for the basic element
      void full_basic_polynomials(const Vector<double>& s, Shape& pn) const;

      /// Get first derivatives of full basis for the basic element
      void dfull_basic_polynomials(const Vector<double>& s, DShape& dpn) const;

      /// Get second derivatives of full basis for the basic element
      void d2full_basic_polynomials(const Vector<double>& s,
                                    DShape& d2pn) const;

      /// \short Get first derivatives of the 36(55) basis monomials for generic
      /// p7(9) polynomial
      void dfull_basis_monomials(const Vector<double>& s, DShape& dp7) const;

      /// \short Get second derivatives of the 36(55) basis monomials for
      /// generic
      // p7(9) polynomial
      void d2full_basis_monomials(const Vector<double>& s, DShape& d2p7) const;

      /// Get the basis functions at basic coordinate s
      void shape_basic(const Vector<double>& s, Shape& psi, Shape& bpsi) const;

    public:
      /// Get the basis functions at local coordinate s
      void shape(const Vector<double>& s, Shape& psi, Shape& bpsi) const
      {
// check the construction of the elements is complete
#ifdef PARANOID
        if (Curved_edge == none)
        {
          throw OomphLibError("the element has not been upgraded yet. did \
  you forget to set upe the curved_edge?",
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
        /// Permute the local coordinate
        Vector<double> s_basic(s);
        permute_shape(s_basic);
        shape_basic(s_basic, psi, bpsi);
      }

      /// Get the local to eulerian Jacobian at local coordinate s (not basic)
      void get_jacobian(const Vector<double>& s,
                        DenseMatrix<double>& jacobian) const
      {
// check the construction of the elements is complete
#ifdef PARANOID
        if (Curved_edge == none)
        {
          throw OomphLibError("the element has not been upgraded yet. did \
  you forget to set upe the curved_edge?",
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
        // Permute the local coordinate
        Vector<double> s_basic(s);
        permute_shape(s_basic);
        // Get basic to eulerian jacobian with permuted shape, copy construct
        DenseMatrix<double> jac_b2e(2, 2, 0.0);
        DenseMatrix<double> jac_b2l(2, 2, 0.0);
        get_basic_jacobian(s_basic, jac_b2e);
        get_jacobian_of_permute(jac_b2l);

        jacobian = DenseMatrix<double>(2, 2, 0.0);
        // Now convert to the local to eulerian Jacobian
        for (unsigned i = 0; i < 2; ++i)
        {
          for (unsigned j = 0; j < 2; ++j)
          {
            for (unsigned k = 0; k < 2; ++k)
            {
              // Now take product
              jacobian(i, k) += jac_b2e(i, j) * jac_b2l(j, k);
            }
          }
        }
      }

    protected:
      ///  Get the local first derivatives of the basis functions
      void d_shape_ds(const Vector<double>& s,
                      Shape& psi,
                      Shape& bpsi,
                      DShape& dpsi,
                      DShape& dbpsi) const; // PRIVATE

      // HERE WRITE A PUBLIC dshape_local

      ///  Get the local second derivatives of the basis functions
      void d2_shape_ds2(const Vector<double>& s,
                        Shape& psi,
                        Shape& bpsi,
                        DShape& dpsi,
                        DShape& dbpsi,
                        DShape& d2psi,
                        DShape& d2bpsi) const; // PRIVATE

      /// \short  Get the local second derivatives of the basis functions with
      /// precomputed association matrix
      void d2_shape_ds2(const Vector<double>& s,
                        Shape& psi,
                        Shape& bpsi,
                        DShape& dpsi,
                        DShape& dbpsi,
                        DShape& d2psi,
                        DShape& d2bpsi,
                        const DenseMatrix<double>& m) const; // PRIVATE

      // HERE WRITE A PUBLIC d2shape_local

      // HERE this is a bit dodgy
      /// \short Array to hold the weights and knots (defined in cc file)
      static const double Internal_dof_knots[(BOUNDARY_ORDER == 3 ? 3 : 10)][2];

    public:
      /// Get the Eulerian first derivatives of the basis functions
      double d_shape_dx(const Vector<double>& s,
                        Shape& psi,
                        Shape& bpsi,
                        DShape& dpsi,
                        DShape& dbpsi) const;

      /// Get the Eulerian second derivatives of the basis functions
      double d2_shape_dx2(const Vector<double>& s,
                          Shape& psi,
                          Shape& bpsi,
                          DShape& dpsi,
                          DShape& dbpsi,
                          DShape& d2psi,
                          DShape& d2bpsi) const;

      /// Get the Eulerian second derivatives of the basis functions, returning
      /// the association matrix used to assemble the basis.
      double d2_shape_dx2(const Vector<double>& s,
                          Shape& psi,
                          Shape& bpsi,
                          DShape& dpsi,
                          DShape& dbpsi,
                          DShape& d2psi,
                          DShape& d2bpsi,
                          const DenseMatrix<double>& m) const;

      /// Return the Eulerian coordinate of the ith internal dof.
      void get_internal_dofs_location(const unsigned& idof,
                                      Vector<double>& s_permute) const
      {
// check the construction of the elements is complete
#ifdef PARANOID
        if (Curved_edge == none)
        {
          throw OomphLibError("The element has not been upgraded yet. did \
  you forget to set up the curved_edge?",
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }
#endif
#ifdef RANGE_CHECK
        if (idof > n_internal_dofs() - 1)
        {
          // Output Sensible error message
          char message[100];
          sprintf(message,
                  "Internal dof %i out of bound %i for curved element internal "
                  "`bubble' dofs.",
                  idof,
                  n_internal_dofs());
          throw OomphLibError(
            message, OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
        }
#endif

        // Get the shape
        Vector<double> s_basic(2);
        s_basic[0] = Internal_dof_knots[idof][0]; // HERE RANGE CHECK
        s_basic[1] = Internal_dof_knots[idof][1]; // HERE RANGE CHECK
        permute_shape(s_basic);
        permute_shape(s_basic); // permute once more, to get -1 as it is cyclic
        // Copy over
        s_permute = s_basic;
      }
    };

    /// Return the number of basis functions on the physical triangle
    template<>
    inline unsigned BernadouElementBasis<3>::n_basis_functions() const
    {
      return 21;
    }

    /// Return the number of basis functions on the physical triangle
    template<>
    inline unsigned BernadouElementBasis<5>::n_basis_functions() const
    {
      return 28;
    }

    /// Get the edge permutation
    template<unsigned BOUNDARY_ORDER>
    inline void BernadouElementBasis<BOUNDARY_ORDER>::nodal_index_shift(
      unsigned& index_shift) const
    {
      // Set the index shift
      // If the element has been upgraded
      if (Curved_edge == none)
      {
        // There is no reasonable definition for the shape functions in this
        // case
        throw OomphLibError(
          "The edge has not been set, so an edge permutation cannot be defined. Did \
you forget to set a Curved_edge?",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
      }
      else if (Curved_edge == zero)
      {
        index_shift = 1;
      }
      else if (Curved_edge == one)
      {
        index_shift = 2;
      }
      else // i.e. if (Curved_edge==two)
      {
        index_shift = 0;
      }
    }

    /// Permute shape based on which edge is curved
    template<unsigned BOUNDARY_ORDER>
    inline void BernadouElementBasis<BOUNDARY_ORDER>::permute_shape(
      Vector<double>& s) const
    {
      // Permute the shape coordinate
      // If the element has been upgraded
      if (Curved_edge == none)
      {
        // There is no reasonable definition for the shape functions in this
        // case
        throw OomphLibError(
          "The edge has not been set, so an edge permutation cannot be defined. Did \
you forget to set a Curved_edge?",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
      }
      else if (Curved_edge == zero)
      {
        // We need to permute the local coordinate
        Vector<double> permuted_s(2, 0.0);
        permuted_s[0] = s[1];
        permuted_s[1] = 1 - s[0] - s[1];
        // Copy over
        s = permuted_s;
      }
      else if (Curved_edge == one)
      {
        // We need to permute the local coordinate
        Vector<double> permuted_s(2, 0.0);
        permuted_s[0] = 1 - s[0] - s[1];
        permuted_s[1] = s[0];
        // Copy over
        s = permuted_s;
      }
      else // i.e. if (Curved_edge==two)
      {
        // The local coordinate remains unchanged
        // s = s
      }
    }

    /// Get the jacobian associated with the edge permutation
    template<unsigned BOUNDARY_ORDER>
    inline void BernadouElementBasis<BOUNDARY_ORDER>::get_jacobian_of_permute(
      DenseMatrix<double>& jac) const
    {
      // Permute the shape coordinate
      // If the element has been upgraded
      if (Curved_edge == none)
      {
        // There is no reasonable definition for the shape functions in this
        // case
        throw OomphLibError(
          "The edge has not been set, so an edge permutation cannot be defined. Did \
you forget to set a Curved_edge?",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
      }
      else if (Curved_edge == zero)
      {
        // We need the derivative of the permuted coords wrt the local
        // coordinate
        jac(0, 0) = 0.0;
        jac(0, 1) = 1.0;
        jac(1, 0) = -1.0;
        jac(1, 1) = -1.0;
      }
      else if (Curved_edge == one)
      {
        // We need the derivative of the permuted coords wrt the local
        // coordinate
        jac(0, 0) = -1.0;
        jac(0, 1) = -1.0;
        jac(1, 0) = 1.0;
        jac(1, 1) = 0.0;
      }
      else // i.e. if (Curved_edge==two)
      {
        // The local coordinate remains unchanged
        // We need the derivative of the permuted coords wrt the local
        // coordinate
        jac(0, 0) = 1.0;
        jac(0, 1) = 0.0;
        jac(1, 0) = 0.0;
        jac(1, 1) = 1.0;
      }
    }

    // Inline functions
    /// Self check function:
    /// 1. Checks for inverted elements.
    /// 2. Checks that the specified parametric edge agrees at s_ubar and s_obar
    ///    with the vertices
    /// 3. Checks that the curved edge is not parallel to the adjacent side
    /// 4. Checks that the side lengths do not become close to zero.
    /// 5. Check that the area of the paralellogram traced by the two non curved
    ///    edges and the tangent vectors does not have zero area
    /// 6. Check that the element is not collapsed onto a line.
    /// 7. Extra paranoid check to see if ANY of the denominators needed in the
    ///    construction are zero: these cases should be caught by previous
    ///    checks.
    template<unsigned BOUNDARY_ORDER>
    void BernadouElementBasis<BOUNDARY_ORDER>::self_check(
      const CurvilineGeomObject& parametric_curve_pt) const
    {
      // Tolerance as a static member HERE
      const double tol(1e-15), angle_tol(1e-12);

      // Check that all of the relevant fields have been filled. HERE (HIGHER in
      // complete build of shape function)

      // Check that the element is not inverted HERE
      const bool element_is_inverted = false;
      if (element_is_inverted)
      {
        throw OomphLibError(
          "This element has clockwise vertices which will cause problems in the shape\
definitions.",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
      }

      // Check that the vertices are identical to those provided by the
      // parametric function
      Vector<Vector<double>> local_vertices(3, Vector<double>(2, 0.0));
      Vector<double> vertex_0(2, 0.0), vertex_1(2, 0.0);
      parametric_curve_pt.position(Vector<double>(1, S_ubar), vertex_0);
      parametric_curve_pt.position(Vector<double>(1, S_obar), vertex_1);

      // Magnitude of the difference
      const double diff0 = sqrt(pow(vertex_0[0] - Vertices[0][0], 2) +
                                pow(vertex_0[1] - Vertices[0][1], 2));
      const double diff1 = sqrt(pow(vertex_1[0] - Vertices[1][0], 2) +
                                pow(vertex_0[1] - Vertices[0][1], 2));
      const bool vertices_differ_from_curve = diff0 > tol || diff1 > tol;

      // The parametric curve does not start and end at the vertices.
      if (vertices_differ_from_curve)
      {
        oomph_info << "Difference of " << diff0 << " " << diff1
                   << " between assigned vertices 0"
                   << " and 1 respectively.\n";
        throw OomphLibError(
          "Non zero difference detected between assigned vertices \
and the start and end of the provided Parametric boundary.",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
      }

      // Lengths of the vectors
      const double lA1 = sqrt(A1(0) * A1(0) + A1(1) * A1(1)),
                   lA2 = sqrt(A2(0) * A2(0) + A2(1) * A2(1)),
                   lB1 = sqrt(A1(0) * A1(0) + A1(1) * A1(1)),
                   lB2 = sqrt(A2(0) * A2(0) + A2(1) * A2(1));

      // Check lengths
      if (fabs(lA1) < tol || fabs(lB2) < tol)
      {
        throw OomphLibError(
          "The length of one of the sides became excessively small. Are you sure you \
defined the vertices correctly?",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
      }
      // Check lengths
      if (fabs(lB1) < tol || fabs(lA2) < tol)
      {
        throw OomphLibError(
          "The length of one of the parametric curve tangents became excessively small. \
Are you sure you defined the boundary correctly?",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
      }

      // Check that atilde2 is not close to zero - this happens when A1 is
      // parallel to a1-a2. This should also catch the special case atilde1 = -1
      // which occurs when  A1 is exactly antiparallel to a1-a2 unless the
      // parrallelogram area is zero which is caught below (similar checks for
      // b1 and b2).

      // This can be summarised as a check that the Element isn't flat.
      // This should catch the cases when a_tilde_2 and b_tilde_1 area zero and
      // also the cases when a_tilde_1=-1 and b_tilde_2=-1. The checks below
      // should find the cases when the curve tangent is parallel to the
      // adjacent side.
      const bool element_is_flat_to_tol =
        fabs(parallelogram_area(A1(0), A1(1), B2(0), B2(1)) / (lA1 * lB2)) <
        angle_tol;

      if (element_is_flat_to_tol)
      {
        throw OomphLibError(
          "The cross product between the tangent to side 0 and side 1 became too \
small. This corresponds to a flat element which should not be produced \
by any quality mesh generator. Are you sure you defined the vertices correctly?",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
      }


      // Check that the element won't break due to parallelogram area (angle)
      // being zero  (i.e. the two tangents being at a vertex being linearly
      // dependant) These checks should also catch when c_alpha tilde and double
      // tilde are undefined. Vertex 0
      const bool parallelogram_area_0_is_zero_to_tol =
        fabs(parallelogram_area(A1(0), A1(1), A2(0), A2(1)) / (lA1 * lA2)) <
        angle_tol;

      if (parallelogram_area_0_is_zero_to_tol)
      {
        throw OomphLibError(
          "The cross product between the curved edge tangent and the straight\
adjacent straight side tangent at vertex 0 became too small. These vectors can \
not be linearly dependent. Are you sure you defined the boundary correctly?",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
      }

      // Vertex 1
      const bool parallelogram_area_1_is_zero_to_tol =
        fabs(parallelogram_area(B1(0), B1(1), B2(0), B2(1)) / (lB1 * lB2)) <
        angle_tol;
      if (parallelogram_area_1_is_zero_to_tol)
      {
        throw OomphLibError(
          "The cross product between the curved edge tangent and the straight\
adjacent straight side tangent at vertex 1 became too small. These vectors can \
not be linearly dependent. Are you sure you defined the boundary correctly?",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
      }
      // Now check that A2 is not parallel to B2
      const bool a_tilde_1_is_minus_1_to_tol =
        fabs(parallelogram_area(A2(0), A2(1), B2(0), B2(1)) /
             parallelogram_area(A1(0), A1(1), A2(0), A2(1))) < angle_tol;

      if (a_tilde_1_is_minus_1_to_tol)
      {
        throw OomphLibError(
          "The curve tangent at vertex 0 (A2) and the side 0 tangent (B2) are parallel \
to the prescribed tolerance. These vectors cannot be linearly dependant. Are \
you sure you defined the boundary correctly?\n",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
      }
      // Now check that B1 is not parallel to A1
      const bool b_tilde_2_is_minus_1_to_tol =
        fabs(parallelogram_area(B1(0), B1(1), A1(0), A1(1)) /
             parallelogram_area(B2(0), B2(1), B1(0), B1(1))) < angle_tol;

      if (b_tilde_2_is_minus_1_to_tol)
      {
        throw OomphLibError(
          "The curve tangent at vertex 1 (B1) and the side 1 tangent (A1) are parallel\
 to the prescribed tolerance. These vectors cannot be linearly dependant. Are\
 you sure you defined the boundary correctly?\n",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
      }

      // Now check we have no denominators that get close to zero. These should
      // be caught by previous checks but double check to be certain
      if (fabs(1 + a_tilde_1()) < tol)
      {
        throw OomphLibError("One of the denominators (1+a_tilde_1) \
 in the association matrix turned out to be zero. This should have been caught \
by previous checks and needs further investigation.\n",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

      if (fabs(a_tilde_2()) < tol)
      {
        throw OomphLibError("One of the denominators (a_tilde_2)\
 in the association matrix turned out to be zero. This should have been caught\
 by previous checks and needs further investigation.\n",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

      if (fabs(b_tilde_1()) < tol)
      {
        throw OomphLibError("One of the denominators (b_tilde_1)\
 in the association matrix turned out to be zero. This should have been caught\
 by previous checks and needs further investigation.\n",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

      if (fabs(1 + b_tilde_2()) < tol)
      {
        throw OomphLibError("One of the denominators (1+b_tilde_2)\
 in the association matrix turned out to be zero. This should have been caught\
 by previous checks and needs further investigation.\n",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

      if (fabs(c_tilde_2() * c_tildetilde_1() +
               c_tilde_1() * c_tildetilde_2()) < tol)
      {
        throw OomphLibError(
          "One of the denominators (c_tilde_1 c_tildetilde_2+c_tilde_2 c_tildetilde_1)\
 in the association matrix turned out to be zero. This should have been caught\
 by previous checks and needs further investigation.\n",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
      }
    }
  } // namespace MyC1CurvedElements
} // namespace oomph
#endif

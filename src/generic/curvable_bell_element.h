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
#ifndef OOMPH_CURVABLE_BELL_ELEMENTS_HEADER
#define OOMPH_CURVABLE_BELL_ELEMENTS_HEADER

// [zdec] I was having such a nightmare with this level in the class heirarchy
// that I put the whole thing into Astra 6.0 and argued with it until we both
// agreed. It did a very good job of refactoring and I learnt a lot about smart
// pointers.
//
// I decided to scrap the entire SubparametricTriangleElement class, instead
// using the template free class as a base.
//
// The AI did suggest (by itself!) creating a separation of the "plate-specific"
// stuff from C1 elements in general. And I am allowing it just to see, although
// (whilst the names of the certain methods, e.g. boundary condition helpers,
// are plate-specific) the methods seem quite general for all sorts of boundary
// conditions which are often relative to the normal and are analogous to
// displacement and forcing.
//
// I also moved the management of different fields to the full finite elements
// and keep this a container for the geometry and basis information only. I may
// change my mind on this second point though and go for something in-between,
// where this class has helpers to assign memory based on the number of fields
// but we will see...
//
//   -- Aidan.

// oomph-lib headers
#include <cmath>
#include <memory>
#include "shape.h"
#include "Telements.h"
#include "c1_curved_elements_direct.h"
#include "c1_plate_helper.h"

namespace oomph
{
  //===========================================================================
  /// Template-free base class for curvable Bell Element
  //===========================================================================
  class TemplateFreeCurvableBellElement : public virtual FiniteElement
  {
  public:
    /// Upgrade the element to be curved
    virtual void upgrade_element_to_curved(
      const C1PlateHelper::CurvedEdgeEnumeration& curved_edge,
      const double& s_ubar,
      const double& s_obar,
      C1CurviLine* parametric_edge,
      const unsigned& boundary_order) = 0;

    /// Destroy the element through its template-free interface.
    virtual ~TemplateFreeCurvableBellElement() = default;

    /// Return whether the geometry has a curved edge.
    virtual bool element_is_curved() const = 0;

    /// Return the number of scalar internal Bell basis functions.
    virtual unsigned ninternal_basis() const = 0;

    /// Return the reference position of an internal Bell degree of freedom.
    virtual void get_internal_dofs_location(const unsigned& dof,
                                            Vector<double>& s) const = 0;

    /// Retain the inherited mapping overloads taking shape derivatives.
    using FiniteElement::local_to_eulerian_mapping;

    /// Evaluate the geometry Jacobian (dx_j/ds_i), its inverse and determinant.
    virtual double local_to_eulerian_mapping(
      const Vector<double>& s,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& inverse_jacobian) const = 0;

    /// Evaluate the scalar nodal and internal Bell basis functions.
    /// Callers provide nodal storage for three vertices and six dof types,
    /// and internal storage for ninternal_basis() functions (zero if straight).
    virtual void c1_basis(const Vector<double>& s,
                          Shape& nodal_basis,
                          Shape& bubble_basis) const = 0;

    /// Evaluate the Bell basis and physical first derivatives; return det(J).
    virtual double d_c1_basis_eulerian(const Vector<double>& s,
                                       Shape& nodal_basis,
                                       Shape& bubble_basis,
                                       DShape& dnodal_basis,
                                       DShape& dbubble_basis) const = 0;

    /// Evaluate the Bell basis and physical first/second derivatives;
    /// return det(J). Nodal types are value, x, y, xx, xy, yy derivatives.
    virtual double d2_c1_basis_eulerian(const Vector<double>& s,
                                        Shape& nodal_basis,
                                        Shape& bubble_basis,
                                        DShape& dnodal_basis,
                                        DShape& dbubble_basis,
                                        DShape& d2nodal_basis,
                                        DShape& d2bubble_basis) const = 0;
  };


  //===========================================================================
  /// Boundary-condition interface for C1 plate elements sharing Bell geometry.
  //===========================================================================
  class TemplateFreeC1PlateElement
    : public virtual TemplateFreeCurvableBellElement
  {
  public:
    /// Access function to rotated boundary helper object
    virtual RotatedBoundaryHelper* rotated_boundary_helper_pt() = 0;


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
      TriangleMeshCurviLine* curviline_pt) = 0;


    /// Pin i.e. pin the in-plane and out of plane displacements only.
    /// We also apply implied boundary conditions (e.g. specification of w
    /// also implies dw/dt etc. boundary_values_pt[i] describes boundary
    /// conditions for the three displacement components (i=0,1 in plane (x,y);
    /// i=2: out-of-plane (z)).
    /// curviline_pt provides a pointer to the representation of the curvilinear
    /// boundary in the triangle mesh.
    virtual void pin_specified_boundary(
      const unsigned& b,
      const Vector<BoundaryConditionForC1PlateBending*>& boundary_values_pt,
      TriangleMeshCurviLine* curviline_pt) = 0;


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
    virtual DuplicateNodeConstraintElement* duplicate_constraint_element_factory(
      Node* const& left_node_pt,
      Node* const& right_node_pt,
      C1CurviLine* const& left_boundary_pt,
      C1CurviLine* const& right_boundary_pt,
      Vector<double> const& left_coord,
      Vector<double> const& right_coord) = 0;
  };


  //===========================================================================
  /// Triangle geometry and a scalar Bell basis. NNODE_1D specifies the
  /// available Lagrange node topology. Straight geometry uses only the three
  /// vertices with barycentric interpolateion; curved geometry uses the
  /// Bernadou-style curved map.
  //===========================================================================
  template<unsigned NNODE_1D>
  class CurvableBellElement : public virtual TemplateFreeCurvableBellElement,
                              public TElement<2, NNODE_1D>
  {
  public:
    /// Construct a straight element with the default Bell quadrature rule.
    CurvableBellElement()
      : Curved_edge(C1PlateHelper::CurvedEdgeEnumeration::none),
        Default_integral_pt(new TGauss<2, 13>)
    {
      this->set_integration_scheme(Default_integral_pt.get());
    }

    /// Release the owned basis and default integration rule.
    ~CurvableBellElement() override = default;

    /// Geometry and quadrature ownership cannot be copied.
    CurvableBellElement(const CurvableBellElement&) = delete;

    /// Geometry and quadrature ownership cannot be assigned.
    CurvableBellElement& operator=(const CurvableBellElement&) = delete;

    ///  Boolean function indicating whether element is curved or not
    bool element_is_curved() const override
    {
      return Curved_edge != C1PlateHelper::CurvedEdgeEnumeration::none;
    }


    /// Return the number of basis types at each (vertex) node for the
    /// Bell/Bernadou interpolation
    unsigned nbell_nodal_basis_type() const
    {
      return 6;
    }


    /// Number of internal functions in the scalar basis (zero if straight).
    unsigned ninternal_basis() const override
    {
      return element_is_curved() ?
               Bernadou_element_basis_pt->n_internal_dofs() :
               0;
    }

    /// Borrow the curved basis, or return null for a straight element.
    TemplateFreeBernadouElementBasis* bernadou_element_basis_pt()
    {
      return Bernadou_element_basis_pt.get();
    }

    /// Get the interpolated global coordinates x -- simplex interpolation if
    /// straight edged, uses the Bernadou basis if curved.
    void interpolated_x(const Vector<double>& s,
                        Vector<double>& x) const override
    {
      // Wrapper: call the BernadouElement curved basis if curved
      if (element_is_curved())
      {
        Bernadou_element_basis_pt->coordinate_x(s, x);
      }
      else
      {
        for (unsigned i = 0; i < 2; ++i)
        {
          x[i] = s[0] * this->nodal_position(0, i) +
                 s[1] * this->nodal_position(1, i) +
                 (1.0 - s[0] - s[1]) * this->nodal_position(2, i);
        }
      }
    }

    /// Return the interpolated global coordinate x_i (slow as it just calls
    /// the full interpolated_x and returns the relevant one)
    double interpolated_x(const Vector<double>& s,
                          const unsigned& i) const override
    {
      // Just call interpolated_x and discard the other component (slow)
      Vector<double> x(2, 0.0);
      interpolated_x(s, x);
      return x[i];
    }

    /// Evaluate affine nodal geometry weights; reject curved geometry.
    // [zdec] TODO: It should be possible to implement geometric weights such
    // that the curved Bell basis acts as shape functions. I.e. we can choose
    // X_ij such that:
    //     x_i = sum_j (X_ij * psi_j(s))
    // This would also extend to the dshape and d2shape methods below.
    // Obviously, the pullback F(s) will be interpolated correctly at nodes and
    // internal points, the question is whether F(s) is in the span of the
    // curved Bell basis.
    virtual void shape(const Vector<double>& s, Shape& shape) const override
    {
      if (element_is_curved())
      {
        throw OomphLibError(
          "No shape defined for these elements: use interpolated_x \
to access interpolated eulerian coordinate",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
      }
      // Straight geometry uses only vertex weights.
      else
      {
        simplex_shape(s, shape);
      }
    }

    /// Evaluate affine geometry weights and derivatives; reject curved
    /// geometry.
    virtual void dshape_local(const Vector<double>& s,
                              Shape& shape,
                              DShape& dshape) const override
    {
      if (element_is_curved())
      {
        throw OomphLibError(
          "No shape defined for these elements: use interpolated_x \
to access interpolated eulerian coordinate",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
      }
      // Straight geometry uses only vertex weights.
      else
      {
        dsimplex_shape_local(s, shape, dshape);
      }
    }


    /// Evaluate affine geometry weights and first/second reference derivatives.
    /// Curved geometry has no nodal shape representation: dshape_local rejects
    /// it.
    void d2shape_local(const Vector<double>& s,
                       Shape& psi,
                       DShape& dpsids,
                       DShape& d2psids) const override
    {
      dshape_local(s, psi, dpsids);
      for (unsigned l = 0; l < this->nnode(); ++l)
      {
        for (unsigned i = 0; i < 3; ++i)
        {
          d2psids(l, i) = 0.0;
        }
      }
    }

    /// Retain the inherited mapping overloads taking shape derivatives.
    using TElement<2, NNODE_1D>::local_to_eulerian_mapping;


    /// Local_to_eulerian mapping with local coordinate argument: when upgraded
    /// this uses the Bernadou implementation of the Jacobian
    virtual double local_to_eulerian_mapping(
      const Vector<double>& s,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& inverse_jacobian) const override
    {
      if (element_is_curved())
      {
        // The direct basis returns dx_i/ds_j; oomph uses dx_j/ds_i.
        Bernadou_element_basis_pt->get_jacobian(s, jacobian);
        // Adapt the index convention before inversion.
        double tmp = jacobian(0, 1);
        jacobian(0, 1) = jacobian(1, 0);
        jacobian(1, 0) = tmp;
        // Invert the Jacobian and return the determinant
        return FiniteElement::invert_jacobian<2>(jacobian, inverse_jacobian);
      }
      else
      {
        // Assemble dshape to get the mapping
        Shape psi(this->nnode());
        DShape dpsi(this->nnode(), this->dim());
        dsimplex_shape_local(s, psi, dpsi);
        return TElement<2, NNODE_1D>::local_to_eulerian_mapping(
          dpsi, jacobian, inverse_jacobian);
      }
    }


    /// Return the positive area scaling of the physical geometry map.
    double J_eulerian(const Vector<double>& s) const override
    {
      DenseMatrix<double> jacobian(2, 2), inverse_jacobian(2, 2);
      return std::abs(local_to_eulerian_mapping(s, jacobian, inverse_jacobian));
    }

    /// Return the positive area scaling at a quadrature point.
    double J_eulerian_at_knot(const unsigned& ipt) const override
    {
      Vector<double> s(2);
      for (unsigned i = 0; i < 2; ++i)
      {
        s[i] = this->integral_pt()->knot(ipt, i);
      }
      return J_eulerian(s);
    }

    /// Get the Bell/Bernadou basis for the unknowns
    virtual void c1_basis(const Vector<double>& s,
                          Shape& nodal_basis,
                          Shape& bubble_basis) const override
    {
      if (element_is_curved())
      {
        Bernadou_element_basis_pt->shape(s, nodal_basis, bubble_basis);
      }
      else
      {
        unsigned n_node = this->nvertex_node();
        unsigned n_type = nbell_nodal_basis_type();
        unsigned dim = this->dim();
        unsigned n_deriv = dim;
        unsigned n2_deriv = dim * (dim + 1) / 2;
        Vector<Vector<double>> Verts =
          (Vector<Vector<double>>(n_node, Vector<double>(dim, 0.0)));
        for (unsigned i_vert = 0; i_vert < n_node; ++i_vert)
        {
          for (unsigned j_coord = 0; j_coord < this->dim(); ++j_coord)
          {
            Verts[i_vert][j_coord] = this->nodal_position(i_vert, j_coord);
          }
        }
        DShape dummydshape(n_node, n_type, n_deriv);
        DShape dummyd2shape(n_node, n_type, n2_deriv);
        Bell_basis.d2_basis_eulerian(
          s, Verts, nodal_basis, dummydshape, dummyd2shape);
      }
    };


    /// Get the global (Eulerian) derivative of the basis for the unknowns
    /// Straight-element evaluation uses the Bell basis second-derivative
    /// routine.
    virtual double d_c1_basis_eulerian(const Vector<double>& s,
                                       Shape& nodal_basis,
                                       Shape& bubble_basis,
                                       DShape& dnodal_basis,
                                       DShape& dbubble_basis) const override
    {
      if (element_is_curved())
      {
        return Bernadou_element_basis_pt->d_shape_dx(
          s, nodal_basis, bubble_basis, dnodal_basis, dbubble_basis);
      }
      else
      {
        unsigned n_node = this->nvertex_node();
        unsigned n_type = nbell_nodal_basis_type();
        unsigned dim = this->dim();
        unsigned n2_deriv = dim * (dim + 1) / 2;

        Vector<Vector<double>> Verts = (Vector<Vector<double>>(
          this->nvertex_node(), Vector<double>(this->dim(), 0.0)));
        for (unsigned ivert = 0; ivert < this->nvertex_node(); ++ivert)
        {
          for (unsigned icoord = 0; icoord < this->dim(); ++icoord)
          {
            Verts[ivert][icoord] = this->nodal_position(ivert, icoord);
          }
        }
        DShape dummyd2shape(n_node, n_type, n2_deriv);
        Bell_basis.d2_basis_eulerian(
          s, Verts, nodal_basis, dnodal_basis, dummyd2shape);
        return J_eulerian(s);
      }
    }


    /// Get the global (Eulerian) second derivative of the basis for the
    /// unknowns
    virtual double d2_c1_basis_eulerian(const Vector<double>& s,
                                        Shape& nodal_basis,
                                        Shape& bubble_basis,
                                        DShape& dnodal_basis,
                                        DShape& dbubble_basis,
                                        DShape& d2nodal_basis,
                                        DShape& d2bubble_basis) const override
    {
      if (element_is_curved())
      {
        return Bernadou_element_basis_pt->d2_shape_dx2(s,
                                                       nodal_basis,
                                                       bubble_basis,
                                                       dnodal_basis,
                                                       dbubble_basis,
                                                       d2nodal_basis,
                                                       d2bubble_basis);
      }
      // Use the Bell basis functions if not upgraded
      else
      {
        Vector<Vector<double>> Verts = (Vector<Vector<double>>(
          this->nvertex_node(), Vector<double>(this->dim(), 0.0)));
        for (unsigned ivert = 0; ivert < this->nvertex_node(); ++ivert)
        {
          for (unsigned icoord = 0; icoord < this->dim(); ++icoord)
          {
            Verts[ivert][icoord] = this->nodal_position(ivert, icoord);
          }
        }
        Bell_basis.d2_basis_eulerian(
          s, Verts, nodal_basis, dnodal_basis, d2nodal_basis);
        return J_eulerian(s);
      }
    }


    /// Construct and install a curved basis after successful validation.
    /// The boundary curve is borrowed during construction and is not owned.
    void upgrade_element_to_curved(
      const C1PlateHelper::CurvedEdgeEnumeration& curved_edge,
      const double& s_ubar,
      const double& s_obar,
      C1CurviLine* parametric_edge,
      const unsigned& boundary_order) override
    {
      if (element_is_curved() || parametric_edge == nullptr ||
          curved_edge < C1PlateHelper::zero || curved_edge > C1PlateHelper::two)
      {
        throw OomphLibError("Upgrade requires a straight element, a valid edge "
                            "and a boundary curve.",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }

      std::unique_ptr<TemplateFreeBernadouElementBasis> basis;
      switch (boundary_order)
      {
        case 3:
          basis.reset(new DirectBernadouElementBasis<3>);
          break;
        case 5:
          basis.reset(new DirectBernadouElementBasis<5>);
          break;
        default:
          throw OomphLibError("Only boundary orders 3 and 5 are supported.",
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
      }

      Vector<Vector<double>> vertices(3, Vector<double>(2));
      for (unsigned v = 0; v < 3; ++v)
      {
        for (unsigned i = 0; i < 2; ++i)
        {
          vertices[v][i] =
            this->nodal_position((v + unsigned(curved_edge) + 1) % 3, i);
        }
      }
      basis->upgrade_element(
        vertices, s_ubar, s_obar, curved_edge, *parametric_edge);
      Bernadou_element_basis_pt.swap(basis);
      Curved_edge = curved_edge;
    }

    /// Get the reference location of the internal curved Bell dofs, dof
    void get_internal_dofs_location(const unsigned& dof,
                                    Vector<double>& s) const override
    {
      // If the element has not been upgraded
      if (!element_is_curved())
      {
        throw OomphLibError(
          "There are no internal dofs for these elements as they have not been\
 upgraded to curved elements.",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
      }
      else
      {
        Bernadou_element_basis_pt->get_internal_dofs_location(dof, s);
      }
    };


    /// Function to set the position values of the non-vertex nodes to be their
    /// Lagrangian x and y value.
    ///
    /// The non-vertex nodes are not used in the shape interpolation due to
    /// the subparametric nature of the elements -- however we might expect to
    /// find the global position stored at the Lagrange interpolating nodes
    /// should we look and so this function (re)sets these x values.
    // [zdec] I'm not quite sure why this is even necessary, my guess is it is
    // because the position interpolation is broken before these are initialised
    // the initial assignment is inconsistent. - Aidan.
    void repair_lagrange_node_positions()
    {
      // Loop over the non vertex nodes
      unsigned n_node = this->nnode();
      unsigned n_vertex_node = this->nvertex_node();
      for (unsigned j_node = n_vertex_node; j_node < n_node; j_node++)
      {
        // Store the node pointer
        Node* node_pt = this->node_pt(j_node);
        // Get the local position of the node
        Vector<double> s(2, 0.0);
        this->local_coordinate_of_node(j_node, s);
        // Get the global interpolated position of the node
        Vector<double> interp_x(2, 0.0);
        interpolated_x(s, interp_x);
        // Set the global coordinate of the node using the subparametric
        // Lagrange interpolation
        node_pt->x(0) = interp_x[0];
        node_pt->x(1) = interp_x[1];
      }
    } // End of repair_lagrange_node_positions


  private:
    /// Affine geometry basis: barycentric at vertices, zero elsewhere.
    void simplex_shape(const Vector<double>& s, Shape& psi) const
    {
      // Define the simplex shape functions
      psi[0] = s[0];
      psi[1] = s[1];
      psi[2] = 1.0 - s[0] - s[1];

      // Get the number of vertices and number of total nodes so we can loop
      // over the remaining shape functions
      unsigned n_vertex_node = this->nvertex_node();
      unsigned n_node = this->nnode();
      // Set the rest of shape to zero
      for (unsigned i = n_vertex_node; i < n_node; ++i)
      {
        psi[i] = 0.0;
      }
    }

    /// Evaluate affine geometry weights and their reference derivatives.
    void dsimplex_shape_local(const Vector<double>& s,
                              Shape& psi,
                              DShape& dpsids) const
    {
      // Assign the shape functions
      this->simplex_shape(s, psi);
      // Derivatives
      dpsids(0, 0) = 1.0;
      dpsids(0, 1) = 0.0;
      dpsids(1, 0) = 0.0;
      dpsids(1, 1) = 1.0;
      dpsids(2, 0) = -1.0;
      dpsids(2, 1) = -1.0;

      // Get the number of vertices and number of total nodes so we can loop
      // over the remaining shape functions
      unsigned n_vertex_node = this->nvertex_node();
      unsigned n_node = this->nnode();

      // Set the rest of shape to zero
      for (unsigned i = n_vertex_node; i < n_node; ++i)
      {
        dpsids(i, 0) = 0.0;
        dpsids(i, 1) = 0.0;
      }
    }


    /// Enum to store which edge is curved set to none when element has no
    /// curved edges
    C1PlateHelper::CurvedEdgeEnumeration Curved_edge;

    /// Owned basis and geometry for a curved element; null if straight.
    std::unique_ptr<TemplateFreeBernadouElementBasis> Bernadou_element_basis_pt;

    /// Scalar Bell basis evaluator for straight triangles.
    BellElementBasis Bell_basis;

    /// Owned default quadrature, independent of a subsequently borrowed rule.
    std::unique_ptr<Integral> Default_integral_pt;

  }; // End of CurvableBellElement class
} // namespace oomph
#endif

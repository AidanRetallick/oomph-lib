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
#include "../generic/shape.h"
#include "../generic/Telements.h"
#include "c1_curved_elements.h"
#include "my_geom_object.h"

#ifndef SUBPARAMETRIC_TELEMENTS
#define SUBPARAMETRIC_TELEMENTS
namespace oomph
{
  //==============================================================================
  // SubparametricTElements are subparametric finite elements that provide both
  // basis functions to interpolate the unknowns and shape functions to
  // interpolate the geometry. In general these elements have both nodal basis
  // functions and internal basis functions (corresponding to internal dofs)
  //==============================================================================

  /// Class for subparametric element that have a NNODE_1D-1 order basis for
  /// interpolating unknowns but use simplex shape functions for interpolating
  /// coordinates.
  template<unsigned NNODE_1D>
  class SubparametricTriangleElement : public TElement<2, NNODE_1D>
  {
  public:
    /// Constructor
    SubparametricTriangleElement(){};

    /// Destructor
    ~SubparametricTriangleElement(){};

    // Broken Copy Constructor
    SubparametricTriangleElement(const SubparametricTriangleElement& dummy)
    {
      BrokenCopy::broken_copy(OOMPH_CURRENT_FUNCTION);
    }

    // Broken assignment operator
    void operator=(const SubparametricTriangleElement& dummy)
    {
      BrokenCopy::broken_assign(OOMPH_CURRENT_FUNCTION);
    }

    /// Interpolate global coordinate x using simplex shape functions.
    void interpolated_x(const Vector<double>& s, Vector<double>& x) const
    {
      // Find the number of nodes
      const unsigned n_node = this->nvertex_node();
      // Find the number of positional types
      const unsigned n_position_type = this->nnodal_position_type();
      // Find the dimension stored in the node
      const unsigned nodal_dim = this->nodal_dimension();

      // Assign storage for the local shape function
      Shape psi(this->nnode(), n_position_type);
      // Find the values of shape function
      simplex_shape(s, psi);

      // Loop over the dimensions
      for (unsigned i = 0; i < nodal_dim; i++)
      {
        // Initilialise value of x[i] to zero
        x[i] = 0.0;
        // Loop over the local nodes, vertex nodes are ALWAYS nodes 0 1 2
        for (unsigned l = 0; l < n_node; l++)
        {
          // Loop over the number of dofs
          for (unsigned k = 0; k < n_position_type; k++)
          {
            x[i] += this->nodal_position_gen(l, k, i) * psi(l, k);
          }
        }
      }
    }

    /// TElement<2,2> shape functions, because any higher order geometric scheme
    /// will be redundant for straight sided (Lagrangian) triangular elements.
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

    /// TElement<2,2> dshape functions with respect to local coordinates,
    /// because any higher order geometric scheme will be redundant for straight
    /// sided (Lagrangian) triangular elements.
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

    /// Overloaded shape. Wrapper which ensures we are only retrieving the
    /// simplex shape functions.
    virtual void shape(const Vector<double>& s, Shape& psi) const
    {
      simplex_shape(s, psi);
    }

    /// Overloaded dshape_local. Wrapper which ensures we are only retrieving
    /// the simplex shape functions.
    virtual void dshape_local(const Vector<double>& s,
                              Shape& psi,
                              DShape& dpsids) const
    {
      dsimplex_shape_local(s, psi, dpsids);
    }


    /// Function to set the position values of the non-vertex nodes to be their
    /// Lagrangian x and y value.
    ///
    /// The non-vertex nodes are not used in the shape interpolation due to
    /// the subparametric nature of the elements -- however we might expect to
    /// find the global position stored at the Lagrange interpolating nodes
    /// should we look and so this function (re)sets these x values.
    // [zdec] I'm not quite sure why this is even necessary, my guess is it is
    // because the position interpolation is broken before these are initialised
    // the initial assignment is inconsistent.
    void repair_lagrange_node_positions()
    {
      // Loop over the non vertex nodes
      unsigned n_node = this->nnode();
      unsigned n_vertex_node = this->nvertex_node();
      for(unsigned j_node = n_vertex_node; j_node<n_node; j_node++)
      {
	// Store the node pointer
	Node* node_pt = this->node_pt(j_node);
	// Get the local position of the node
	Vector<double> s(2,0.0);
	this->local_coordinate_of_node(j_node,s);
	// Get the global interpolated position of the node
	Vector<double> interp_x(2,0.0);
	interpolated_x(s, interp_x);
	// Set the global coordinate of the node using the subparametric
	// Lagrange interpolation
	node_pt->x(0) = interp_x[0];
	node_pt->x(1) = interp_x[1];
      }
    } // End of repair_lagrange_node_positions
    
  }; // End of Subparametric TElement class


  //=============================================================================
  /// Curvable Bell element. It inherits simplax shape subparametricity
  /// SubparametricTriangleElement for efficient position interpolation.
  /// It also inherits the NNODE_1D Lagrangian basis.
  /// By default this element is the standard simplex triangle Bell element with
  /// 6dofs per node. The element can be upgraded to curved "triangle" elements
  /// which provide accurate representation of boundary conditions by adding
  /// additional Bernadou basis unknowns to the interior.
  //=============================================================================
  template<unsigned NNODE_1D>
  class CurvableBellElement : public SubparametricTriangleElement<NNODE_1D>
  {
  public:
    /// Constructor that takes the number of fields and a vector of
    /// corresponding bools which is true for each field that should be
    /// interpolated using the Bell/Bernadou basis (by default, one field, Bell
    /// interpolated).
    CurvableBellElement(const unsigned& n_field = 1,
                        const std::vector<bool>& is_bell_interpolated = {true})
      : Curved_edge(MyC1CurvedElements::none),
        Nfield(n_field),
        Field_is_bell_interpolated(is_bell_interpolated),
        First_nodal_type_index_for_field(n_field),
        Index_of_internal_data_for_field(n_field)
    {
      // Use the higher order integration scheme
      TGauss<2, 4>* new_integral_pt = new TGauss<2, 4>;
      this->set_integration_scheme(new_integral_pt);
      // By default, there is no Bernadou basis (straight sided triangle)
      Bernadou_element_basis_pt = 0;
      Association_matrix_pt = 0;

      // Keep track of the number of nodal types
      unsigned n_type_total = 0;

      // Set up each fields data
      for (unsigned i_field = 0; i_field < Nfield; i_field++)
      {
        // Add the (zero) internal data for each field
        Index_of_internal_data_for_field[i_field] =
          this->add_internal_data(new Data(0));

        // Set the index for the first nodal type for each field
        First_nodal_type_index_for_field[i_field] = n_type_total;
        // Add the number of types belonging to i_field to the running total
        n_type_total += nnodal_basis_type_for_field(i_field);
      }
    };


    /// Destructor
    ~CurvableBellElement()
    {
      // Clean Up
      delete Association_matrix_pt;
      delete Bernadou_element_basis_pt;
      // Clean up allocation of integration scheme
      delete this->integral_pt();
    }


    /// Broken copy constructor
    CurvableBellElement(const CurvableBellElement& dummy)
    {
      BrokenCopy::broken_copy(OOMPH_CURRENT_FUNCTION);
    }


    /// Broken copy assignment
    void operator=(const CurvableBellElement& dummy)
    {
      BrokenCopy::broken_assign(OOMPH_CURRENT_FUNCTION);
    }


    /// Alias for enum to enumerate the possible edges that could be curved
    typedef typename MyC1CurvedElements::Edge Edge;


    ///  Boolean function indicating whether element is curved or not
    bool element_is_curved() const
    {
      return Curved_edge != MyC1CurvedElements::none;
    }


    /// Access function for the number of fields
    unsigned nfield() const
    {
      return Nfield;
    }


    /// Return the number of basis types at each (vertex) node for the
    /// Bell/Bernadou interpolation
    unsigned nbell_nodal_basis_type() const
    {
      return 6;
    }


    /// Return a bool telling us whether the ith field interpolated by the
    /// Bell basis, needed when upgrading elements in order to resize internal
    /// data of appropriate fields. Pure virtual as it depends on the equations
    /// requirements.
    virtual bool field_is_bell_interpolated(const unsigned& i_field) const
    {
      return Field_is_bell_interpolated[i_field];
    }


    /// Return the number of nodes for a given field (vertex nodes for Bell
    /// interpolated fields, all nodes for Lagrange interpolated)
    unsigned nnode_for_field(const unsigned& i_field) const
    {
      return (Field_is_bell_interpolated[i_field] ? 3 : this->nnode());
    }


    /// Return the nodal indices for a given field (vertex nodes for Bell
    /// interpolated fields, all nodes for Lagrange interpolated)
    Vector<unsigned> nodal_indices_for_field(const unsigned& i_field) const
    {
      unsigned n_node = nnode_for_field(i_field);
      Vector<unsigned> nodal_indices;
      for (unsigned j_node = 0; j_node < n_node; j_node++)
      {
        nodal_indices.push_back(j_node);
      }
      return nodal_indices;
    }


    /// Return number of nodal dof types for a given field
    unsigned nnodal_basis_type_for_field(const unsigned& i_field) const
    {
      return (Field_is_bell_interpolated[i_field] ? nbell_nodal_basis_type() :
                                                    1);
    };


    /// Return the first nodal type index for field i
    virtual unsigned first_nodal_type_index_for_field(
      const unsigned& i_field) const
    {
      return First_nodal_type_index_for_field[i_field];
    }

    /// Return the nodal value for field i at node j of type k
    virtual double nodal_value_for_field_of_type(const unsigned& i_field,
                                                 const unsigned& j_node,
                                                 const unsigned& k_type) const
    {
      // Get the number of nodal types belonging to other fields before this
      // field.
      unsigned field_first_index = first_nodal_type_index_for_field(i_field);
      return this->raw_nodal_value(j_node, field_first_index + k_type);
    }


    /// Return the index of the internal data for a given field, as each field
    /// is given exactly one piece of internal data to interpolate from, this is
    /// the same as the input field. Virtual as this behaviour may differ in
    /// derived elements.
    virtual unsigned index_of_internal_data_for_field(
      const unsigned& i_field) const
    {
      return i_field;
    }


    /// Return the internal data pt for field i
    virtual Data* internal_data_for_field_pt(const unsigned& i_field) const
    {
      unsigned index = index_of_internal_data_for_field(i_field);
      return this->internal_data_pt(index);
    }


    /// Return number of bubble basis functions for a field. This will be zero
    /// for un-upgraded elements or for Lagrange interpolated fields. Upgrading
    /// introduces additional basis functions depending on the basis
    unsigned ninternal_basis_type_for_field(const unsigned& i_field) const
    {
      return internal_data_for_field_pt(i_field)->nvalue();
    }


    /// Return the internal value for field i of type k
    virtual double internal_value_for_field_of_type(
      const unsigned& i_field, const unsigned& k_type) const
    {
      return internal_data_for_field_pt(i_field)->value(k_type);
    }


    /// Return the t-th history internal value for field i of type k
    virtual double internal_value_for_field_of_type(
      const unsigned& t_time,
      const unsigned& i_field,
      const unsigned& k_type) const
    {
      return internal_data_for_field_pt(i_field)->value(t_time, k_type);
    }


    /// Access function for the Bernadou_element_basis_pt
    MyC1CurvedElements::BernadouElementBasisBase* bernadou_element_basis_pt()
    {
      // [zdec] Should this throw an error if not upgraded or just return null
      // pt?
      return Bernadou_element_basis_pt;
    }

    /// Get the interpolated global coordinates x -- simplex interpolation if
    /// straight edged, uses the Bernadou basis if curved.
    void interpolated_x(const Vector<double>& s, Vector<double>& x) const
    {
      // Wrapper: call the BernadouElement curved basis if curved
      if (element_is_curved())
      {
        Bernadou_element_basis_pt->coordinate_x(s, x);
      }
      else
      {
        SubparametricTriangleElement<NNODE_1D>::interpolated_x(s, x);
      }
    }

    /// Return the interpolated global coordinate x_i (slow as it just calls
    /// the full interpolated_x and returns the relevant one)
    double interpolated_x(const Vector<double>& s, const unsigned& i) const
    {
      // Just call interpolated_x and discard the other component (slow)
      Vector<double> x(2, 0.0);
      interpolated_x(s, x);
      return x[i];
    }

    /// Overloaded shape. Thin wrapper which breaks shape for upgraded elements,
    /// which have a mapping but no shape function defined.
    virtual void shape(const Vector<double>& s, Shape& shape) const
    {
      if (element_is_curved())
      {
        throw OomphLibError(
          "No shape defined for these elements: use interpolated_x \
to access interpolated eulerian coordinate",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
      }
      // Default to TElement shape
      else
      {
        SubparametricTriangleElement<NNODE_1D>::simplex_shape(s, shape);
      }
    }

    /// Overloaded shape. Thin wrapper which breaks shape for upgraded elements,
    /// which have a mapping but no shape function defined.
    virtual void dshape_local(const Vector<double>& s,
                              Shape& shape,
                              DShape& dshape) const
    {
      if (element_is_curved())
      {
        throw OomphLibError(
          "No shape defined for these elements: use interpolated_x \
to access interpolated eulerian coordinate",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
      }
      // Default to TElement shape
      else
      {
        SubparametricTriangleElement<NNODE_1D>::dsimplex_shape_local(
          s, shape, dshape);
      }
    }


    /// Local_to_eulerian mapping with local coordinate argument: when upgraded
    /// this uses the Bernadou implementation of the Jacobian
    virtual double local_to_eulerian_mapping(
      const Vector<double>& s,
      DenseMatrix<double>& jacobian,
      DenseMatrix<double>& inverse_jacobian) const
    {
      if (element_is_curved())
      {
        // Fill in the Jacobian // HERE this provides TRANSVERSE of oomph
        // definition fix lower down if possible // [zdec] what did David mean
        // by this?
        Bernadou_element_basis_pt->get_jacobian(s, jacobian);
        // Temporary fix - just transpose
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
        SubparametricTriangleElement<NNODE_1D>::dsimplex_shape_local(
          s, psi, dpsi);
        return TElement<2, NNODE_1D>::local_to_eulerian_mapping(
          dpsi, jacobian, inverse_jacobian);
      }
    }


    /// Get the Bell/Bernadou basis for the unknowns
    virtual void c1_basis(const Vector<double>& s,
                          Shape& nodal_basis,
                          Shape& bubble_basis) const
    {
      if (element_is_curved())
      {
        Bernadou_element_basis_pt->shape(s, nodal_basis, bubble_basis);
      }
      else
      { // [zdec] inefficient
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


    /// Get the local derivative of the basis for the unknowns
    virtual void d_c1_basis_local(const Vector<double>& s,
                                  Shape& nodal_basis,
                                  Shape& bubble_basis,
                                  DShape& dnodal_basis,
                                  DShape& dbubble_basis) const
    {
      throw OomphLibError("Needs implementing. BLAME AIDAN.",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }


    /// Get the local second derivative of the basis for the unknowns
    virtual void d2_c1_basis_local(const Vector<double>& s,
                                   Shape& nodal_basis,
                                   Shape& bubble_basis,
                                   DShape& dnodal_basis,
                                   DShape& dbubble_basis,
                                   DShape& d2nodal_basis,
                                   DShape& d2bubble_basis) const
    {
      throw OomphLibError("Needs implementing. BLAME AIDAN.",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }


    /// Get the global (Eulerian) derivative of the basis for the unknowns
    // [zdec] Calls the full d2_basis_eulerian (slow and sloppy)
    virtual double d_c1_basis_eulerian(const Vector<double>& s,
                                       Shape& nodal_basis,
                                       Shape& bubble_basis,
                                       DShape& dnodal_basis,
                                       DShape& dbubble_basis) const
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
        return TElement<2, NNODE_1D>::J_eulerian(s);
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
                                        DShape& d2bubble_basis) const
    {
      if (element_is_curved())
      {
        // IF the association matrix_pt is equal to null pointer we don't have
        // one so use the `slow' version
        if (Association_matrix_pt == 0)
        {
          return Bernadou_element_basis_pt->d2_shape_dx2(s,
                                                         nodal_basis,
                                                         bubble_basis,
                                                         dnodal_basis,
                                                         dbubble_basis,
                                                         d2nodal_basis,
                                                         d2bubble_basis);
        }
        // Else use the cached association matrix to compute
        else
        {
          return Bernadou_element_basis_pt->d2_shape_dx2(
            s,
            nodal_basis,
            bubble_basis,
            dnodal_basis,
            dbubble_basis,
            d2nodal_basis,
            d2bubble_basis,
            *Association_matrix_pt);
        }
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
        return TElement<2, NNODE_1D>::J_eulerian(s);
      }
    }


    /// Precompute the association matrix for the curved basis. This is
    /// important for optimisation, as its contruction is expensive. It only
    /// needs constructing once per get_residuals or get_jacobian call.
    void store_association_matrix()
    {
      // If the element has been upgraded
      if (element_is_curved())
      {
        // Create the matrix
        Association_matrix_pt = new DenseMatrix<double>(
          Bernadou_element_basis_pt->n_basis_functions(),
          Bernadou_element_basis_pt->n_basic_basis_functions(),
          0.0);
        // Fill in the matrix
        Bernadou_element_basis_pt->fill_in_full_association_matrix(
          *Association_matrix_pt);
      }
      else
      { /* No association matrix for the Bell Elements. Throw here*/
      }
    };


    /// Delete the association matrix for the curved basis. This is important
    /// for optimisation, as we don't want to be storing a large matrix for
    /// every element.
    void delete_association_matrix()
    {
      // If the element has been upgraded
      if (element_is_curved())
      {
        delete Association_matrix_pt;
        Association_matrix_pt = 0;
      }
      else
      { /* No association matrix for the Bell Elements. Throw here ?*/
      }
    };


    /// Add the a curved element pointer of type BERNADOU_BASIS
    template<typename BERNADOU_BASIS>
    void add_new_curved_basis()
    {
      if (Bernadou_element_basis_pt == 0)
      {
        Bernadou_element_basis_pt = new BERNADOU_BASIS;
      }
      else
      {
        throw OomphLibError(
          "There is already a curved basis for this curvable Bell element.",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
      }
    }


    /// Upgrade the element to be curved
    virtual void upgrade_element_to_curved(const Edge& curved_edge,
                                           const double& s_ubar,
                                           const double& s_obar,
                                           CurvilineGeomObject* parametric_edge,
                                           const unsigned& boundary_order)
    {
      using namespace MyC1CurvedElements;
#ifdef PARANOID
      // Check that we haven't upgraded this element already
      if (Curved_edge != none)
      {
        throw OomphLibError(
          "Cannot upgrade more than a single edge to be curved in C1 Curved Bell \
Elements.",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
      }
#endif
      // Add the curved edge
      Curved_edge = curved_edge;

      Integral* new_integral_pt;
      // Switch for the boundary order
      switch (boundary_order)
      {
        case 3:
	  // [zdec] debug
	  oomph_info << "Upgrading to a 3-er" << std::endl;
          add_new_curved_basis<BernadouElementBasis<3>>();
          new_integral_pt = new TGauss<2, 13>;
          break;
        case 5:
	  // [zdec] debug
	  oomph_info << "Upgrading to a 5-er" << std::endl;
          add_new_curved_basis<BernadouElementBasis<5>>();
          new_integral_pt = new TGauss<2, 16>;
          break;
        default:
          throw OomphLibError("Currently only BernadouElementBasis<3> and "
                              "BernadouElementBasis<5> are implemented.",
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
      }

      // Cannot Null this pointer but we immediatly reset it.
      // Should we have a delete Integral_pt function?
      delete this->integral_pt();
      this->set_integration_scheme(new_integral_pt);

      // Always use same data, just resize according to internal basis
      // requirements.
      unsigned n_field = Nfield;
      unsigned n_bubble = Bernadou_element_basis_pt->n_internal_dofs();
      for (unsigned i_field = 0; i_field < n_field; i_field++)
      {
        if (this->field_is_bell_interpolated(i_field))
        {
          this
            ->internal_data_pt(this->index_of_internal_data_for_field(i_field))
            ->resize(n_bubble);
        }
      }

      // Vector to store nodes
      Vector<Vector<double>> vertices(Vector<Vector<double>>(
        this->nvertex_node(), Vector<double>(this->dim(), 0.0)));
      // Now switch to upgrade
      // The shape functions are designed such that the curved edge is always
      // edge two. So this is where we set that up. Throw an error if an edge is
      // upgraded to none
      if (Curved_edge == none)
      {
        throw OomphLibError(
          "Cannot upgrade edge 'none'. Curved elements must have\
 one side defined by a parametric function.",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
      }
      else
      {
        // Loop over coordinate components then vertices and permute vertices
        for (unsigned i = 0; i < 2; ++i)
        {
          for (unsigned ivert = 0; ivert < this->nvertex_node(); ++ivert)
          {
            vertices[ivert][i] =
              this->node_pt((ivert + unsigned(curved_edge) + 1) % 3)->x(i);
          }
        }
      }
      // Add the vertices to make the shape functions fully functional
      Bernadou_element_basis_pt->upgrade_element(
        vertices, s_ubar, s_obar, curved_edge, *parametric_edge);
    }

    /// Get the reference location of the internal curved Bell dofs, dof
    void get_internal_dofs_location(const unsigned& dof,
                                    Vector<double>& s) const
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


  private:
    /// Enum to store which edge is curved set to none when element has no
    /// curved edges
    MyC1CurvedElements::Edge Curved_edge;

    /// Pointer to Bernadou Element Basis
    MyC1CurvedElements::BernadouElementBasisBase* Bernadou_element_basis_pt;

    /// Basis functions
    MyShape::BellElementBasis Bell_basis;

    /// Pointer to Stored Association matrix
    DenseMatrix<double>* Association_matrix_pt;

    /// Number of fields
    unsigned Nfield;

    /// Vector of bools which dictate whether a field is interpolated
    /// using the Bell/Bernadou basis (use Lagrange if not)
    std::vector<bool> Field_is_bell_interpolated;

    /// Vector of unsigneds which store the first index of a fields nodal types
    /// (should be the sum of all ntypes for fields which come before)
    Vector<unsigned> First_nodal_type_index_for_field;

    /// Indices at which the added internal data is stored
    Vector<unsigned> Index_of_internal_data_for_field;

  }; // End of CurvableBellElement class
} // namespace oomph
#endif

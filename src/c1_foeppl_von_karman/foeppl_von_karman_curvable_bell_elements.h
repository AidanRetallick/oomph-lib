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
#include "../generic/bell_element_basis.h"
#include "../generic/c1_curved_elements.h"
#include "../generic/my_geom_object.h"
#include "../generic/subparametric_Telement.h"

namespace oomph
{
  //============================================================================
  /// FoepplVonKarmanC1CurvableBellElement elements are a subparametric scheme
  /// with  linear Lagrange interpolation for approximating the geometry and
  /// the C1-functions for approximating variables.
  /// These elements use TElement<2,NNODE_1D> with Lagrangian bases to
  /// interpolate the in-plane unknowns u_x, u_y although with sub-parametric
  /// simplex shape interpolation.
  ///
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
    /// Constructor: Call constructors for C1CurvableBellElement and
    /// FoepplVonKarman equations
    FoepplVonKarmanC1CurvableBellElement()
      : CurvableBellElement<NNODE_1D>(Nfield, Field_is_bell_interpolated),
        FoepplVonKarmanEquations(),
	// Rotated_basis_fct_pt(0), // [zdec] old rotation
	Nnodes_to_rotate(0)
    {
      // Use the higher order integration scheme
      TGauss<2, 4>* new_integral_pt = new TGauss<2, 4>;
      this->set_integration_scheme(new_integral_pt);

      // Rotated dof helper
      Rotated_boundary_helper_pt = new RotatedBoundaryHelper(this);
      
      Association_matrix_pt = 0;
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

    // // [zdec] old rotation -- delete
    // /// Function pointer to basis vectors function which sets  basis vectors
    // /// b1 and b2 (which are in general functions of x)
    // typedef void (*BasisVectorsFctPt)(const Vector<double>& x,
    //                                   Vector<double>& b1,
    //                                   Vector<double>& b2,
    //                                   DenseMatrix<double>& db1,
    //                                   DenseMatrix<double>& db2);

    // [zdec] Note for Matthias: this current helper class works for arbitrary
    // boundary nodes on curviline boundaries by storing up to two
    // CurvilinGeomObjects which contain all the necessary parametrisation info.
    // This does not work for straight boundaries where no such object exists
    // however.
    
    //---start of rotation helper class-----------------------------------------
    /// Helper class to contain all the rotation information in the element.
    class RotatedBoundaryHelper
    {
      
    public:

      // [zdec] Do we care that 'Nnode=3' is hard coded?
      /// Constructor: just initialise the member data to their defaults (zeros)
      RotatedBoundaryHelper(FiniteElement* const &parent_element_pt)
	: Parent_element_pt(parent_element_pt),
	  Nboundary_at_node(3,0),
	  Boundary_coordinate_of_node(3,Vector<double>(2,0.0)),
	  Nodal_boundary_parametrisation_pt(3,Vector<CurvilineGeomObject*>(2,0)),
	  Rotation_matrix_at_node(3,DenseMatrix<double>(6,6,0.0))
      {}

      /// Destructor
      ~RotatedBoundaryHelper()
      {}

      /// Acccess for the number of nodes data (read only)
      unsigned nboundary_at_node(const unsigned j_node) const
      {
	return Nboundary_at_node[j_node];
      }
      
      /// Add a new boundary parametrisation to nodes all the nodes in the
      /// vector node_on_boundary
      void add_nodal_boundary_parametrisation(
	const Vector<unsigned> &node_on_boundary,
	const Vector<double> &boundary_coord_of_node,
	CurvilineGeomObject* const &boundary_parametrisation_pt)
      {
	// Loop over all the nodes in node_on_boundary and add the boundary
	// pointer to their vector of boundaries
	unsigned n_node = node_on_boundary.size();
	for(unsigned j=0; j<n_node; j++)
	{
	  // Flag to determine whether we are actuall adding this
	  // parametrisation (we don't if we already have it)
	  bool add_boundary_parametrisation = true;

	  unsigned j_node = node_on_boundary[j];
	  unsigned n_boundary = Nboundary_at_node[j_node];	  
	  
          #ifdef PARANOID
	  // A node can be on a maximum of two boundaries
	  // (if we already have more throw an error)
	  if(n_boundary>1)
	  {
	    std::string error_string = "";
	    std::stringstream error_stringstream;
	    
	    error_stringstream << "Nodes cannot appear on more than two "
			       << "boundaries. An attempt is being made to "
			       << "assign " << j_node
			       << " to a " << n_boundary+1 << "-th boundary."
			       << std::endl;
	    
	    error_stringstream >> error_string;
	    
	    throw OomphLibError(error_string,
				OOMPH_CURRENT_FUNCTION,
				OOMPH_EXCEPTION_LOCATION);
	  }
          #endif
	  
	  // IMPORTANT
	  // If we already have a parametrisation, we need to chack that the
	  // new boundary parametrisation isn't equivalent at this point
	  // (i.e. the parametrisation is not continuous across the two
	  // boundaries)

	  // [zdec] !!! what if t1==t2 but dt1/ds1!=dt2/ds2 ???
	  // Also, need to normalise for different parametrisation velocities

	  // Tolerance to determine whether two tangents are different
	  double diff_tol = 1.0e-8;
	  if(n_boundary==1)
	  {
	    // Get the parametrisation for the first boundary
	    CurvilineGeomObject* old_boundary_parametrisation_pt =
	      Nodal_boundary_parametrisation_pt[j_node][0];
	    Vector<double> old_boundary_coord =
	      {Boundary_coordinate_of_node[j_node][0]};

	    // Get the tangent vector for the first parametrisation
	    Vector<double> old_tangent(2,0.0);
	    old_boundary_parametrisation_pt->dposition(old_boundary_coord,
						       old_tangent);

	    // Get the tangent vector for the second parametrisation
	    Vector<double> new_tangent(2,0.0);
	    boundary_parametrisation_pt->dposition({boundary_coord_of_node[j]},
						   new_tangent);

	    // Compare the two
	    Vector<double> tangent_difference(2,0.0);
	    tangent_difference[0] = new_tangent[0]-old_tangent[0];
	    tangent_difference[1] = new_tangent[1]-old_tangent[1];
	    double tangent_difference_norm =
	      tangent_difference[0]*tangent_difference[0]
	      + tangent_difference[1]*tangent_difference[1];

	    // If the tangents are the same (within tolerance) don't add this
	    // parametrisation
	    if(tangent_difference_norm<diff_tol)
	    {
	      add_boundary_parametrisation=false;
	      oomph_info << "Boundary parametrisations produce tangents ("
			 << old_tangent[0] << ", "
			 << old_tangent[1] << ") and ("
			 << new_tangent[0] << ", "
			 << new_tangent[1] << ") which are equal enough that "
			 << "we consider the boundary continuous."
			 << std::endl;
	    }
	  }
	  
	  if(add_boundary_parametrisation)
	  {
	    // Add the boundary pt
	    Nodal_boundary_parametrisation_pt[j_node][n_boundary] =
	      boundary_parametrisation_pt;
	    
	    // Set the coordinate of node j on this boundary
	    Boundary_coordinate_of_node[j_node][n_boundary] =
	      boundary_coord_of_node[j];
	    
	    // Increment the number of boundaries
	    Nboundary_at_node[j_node]++;
	  
	    // [zdec] DEBUG
	    oomph_info << "NODE " << j_node<< ": " << Parent_element_pt->node_pt(j_node)
		       << " in " << "ELEMENT " << Parent_element_pt
		       << " is on boundary " << boundary_parametrisation_pt << std::endl; 
	    
	    update_rotation_matrices();
	  }	  
	} // end of loop over nodes in node_on_boundary [j]
      } // end of add_nodal_boundary_parametrisation()


      /// Update all rotation matrices (checks if they are needed unless
      /// flag is true)
      void update_rotation_matrices()
      {
	// [zdec] hard coded the three vertex nodes 
	unsigned n_vertex = 3;
	// Loop over each vertex
	for(unsigned j_node=0; j_node<n_vertex; j_node++)
	{
	  // Initialise the two basis vectors and their jacobians
	  Vector<Vector<double>> bi(2, Vector<double>(2, 0.0));
	  Vector<DenseMatrix<double>> dbidx(2, DenseMatrix<double>(2, 2, 0.0));    
	  
	  switch(Nboundary_at_node[j_node])
	  {
	    // Node is not on a boundary, don't rotate
	  case 0:
	    { break; }
	    
	    // Node is on a boundary, rotate to the normal-tangent coordinates
	  case 1:
	    {
	      // Our new coordinate system:
	      //     (l, s)=(normal component, tangent component)
	      // which we define in terms of basis vectors (rescaled)
	      //     ni=dxi/dl / |n|           <-- Jacobian col 1
	      //     ti=dxi/ds / |t|           <-- Jacobian col 2
	      // and their derivatives
	      //     dnidxj=d/dxj(dxi/dl / |n|) <-- Hessian `col' 1
	      //     dtidxj=d/dxj(dxi/ds / |t|) <-- Hessian `col' 2	      

	      // [zdec] we use i and j for brevity
	      // but it should be alpha & beta
	      // Need to write up how the transformation is done

	      // Storage for our basis and derivatives
	      Vector<double> ni(2,0.0);
	      Vector<double> ti(2,0.0);
	      Vector<double> dnids(2,0.0);
	      Vector<double> dtids(2,0.0);

	      // All tensors assumed evaluated on the boundary
	      // Jacobian of inverse mapping
	      DenseMatrix<double> jac_inv(2,2,0.0);
	      // Hessian of mapping [zdec] (not needed because...)
	      Vector<DenseMatrix<double>>
		hess(2,DenseMatrix<double>(2,2,0.0));
	      // Hessian of inverse mapping [zdec] (...this can be found by hand)
	      Vector<DenseMatrix<double>>
		hess_inv(2,DenseMatrix<double>(2,2,0.0));
		
	      // The basis is defined in terms of the boundary parametrisation
	      Vector<double> boundary_coord =
		{Boundary_coordinate_of_node[j_node][0]};
	      CurvilineGeomObject* boundary_pt =
		Nodal_boundary_parametrisation_pt[j_node][0];
	      Vector<double> x(2,0.0);
	      Vector<double> dxids(2,0.0);
	      Vector<double> d2xids2(2,0.0);

	      // Get position (debug)
	      boundary_pt->position(boundary_coord, x);
	      // Get tangent vector
	      boundary_pt->dposition(boundary_coord, dxids);
	      // Get second derivative
	      boundary_pt->d2position(boundary_coord, d2xids2);
	      
	      double mag_t = sqrt(dxids[0]*dxids[0]+dxids[1]*dxids[1]);
	      // ti is the normalised tangent vector
	      ti[0] = dxids[0]/mag_t;
	      ti[1] = dxids[1]/mag_t;
	      // Derivative of (normalised) tangent
	      dtids[0] = d2xids2[0] / std::pow(mag_t,2);
		// - (dxids[0]*d2xids2[0]+dxids[1]*d2xids2[1]) * dxids[0]
		// / std::pow(mag_t,4);
	      dtids[1] = d2xids2[1] / std::pow(mag_t,2);
		// - (dxids[0]*d2xids2[0]+dxids[1]*d2xids2[1]) * dxids[1]
		// / std::pow(mag_t,4);
	      // ni = (t x e_z) implies
	      ni[0] = ti[1];
	      ni[1] =-ti[0];
	      // Same for dnids
	      dnids[0] = dtids[1];
	      dnids[1] =-dtids[0];
	      
	      // Need inverse of mapping to calculate ds/dxi ----------------
	      //   /  dx/dl  dx/ds  \ -1  ___   _1_  /  dy/ds -dx/ds \ .
	      //   \  dy/dl  dy/ds  /     ---   det  \ -dy/dl  dx/dl /
	      //
	      //                          ___  /  dl/dx  dl/dy  \ .  
	      //                          ---  \  ds/dx  ds/dy  /
	      //
	      // Fill out inverse of Jacobian
	      double det = (ni[0]*ti[1] - ni[1]*ti[0]);
	      jac_inv(0,0) = ti[1] / det;
	      jac_inv(0,1) =-ti[0] / det;
	      jac_inv(1,0) =-ni[1] / det; 
	      jac_inv(1,1) = ni[0] / det;

	      // Fill out the Hessian
	      // (unneeded -- can calculate the inverse components by hand)
	      for(unsigned alpha=0; alpha<2; alpha++)
	      {
		// hess[alpha](0,0) = 0.0;
		hess[alpha](0,1) = dnids[alpha];
		hess[alpha](1,0) = dnids[alpha];
		hess[alpha](1,1) = dtids[alpha];
	      }

	      // Fill out inverse of Hessian
	      // H^{-1}abg = J^{-1}ad Hdez J^{-1}eb J^{-1}zg
	      for(unsigned alpha=0; alpha<2; alpha++)
	      {
		for(unsigned beta=0; beta<2; beta++)
		{
		  for(unsigned gamma=0; gamma<2;gamma++)
		  {
		    for(unsigned alpha2=0; alpha2<2; alpha2++)
		    {
		      for(unsigned beta2=0; beta2<2; beta2++)
		      {
			for(unsigned gamma2=0; gamma2<2;gamma2++)
			{
			  hess_inv[alpha](beta,gamma) -= 
			    jac_inv(alpha,alpha2) * hess[alpha2](beta2,gamma2)
			    * jac_inv(beta2,beta) * jac_inv(gamma2,gamma);
			}
		      }
		    }
		  }
		}
	      }

	      // // [zdec] debug
	      // std::ofstream jac_and_hess;
	      // jac_and_hess.open("jac_and_hess_new.csv", std::ios_base::app);
	      // jac_and_hess << "Jacobian inverse:" << std::endl
	      // 		   << bi[0][0] << " " << bi[0][1] << std::endl
	      // 		   << bi[1][0] << " " << bi[1][1] << std::endl
	      // 		   << "Hessian inverse [x]:" << std::endl
	      // 		   << Dbi[0](0,0) << " " << Dbi[0](0,1) << std::endl
	      // 		   << Dbi[0](1,0) << " " << Dbi[0](1,1) << std::endl
	      // 		   << "Hessian inverse [y]:" << std::endl
	      // 		   << Dbi[1](0,0) << " " << Dbi[1](0,1) << std::endl
	      // 		   << Dbi[1](1,0) << " " << Dbi[1](1,1) << std::endl << std::endl;

    

	      // [zdec] debug
	      std::ofstream jac_and_hess;
	      jac_and_hess.open("jac_and_hess_new.csv", std::ios_base::app);
	      jac_and_hess << "Jacobian inverse:" << std::endl
			   << jac_inv(0,0) << " " << jac_inv(0,1) << std::endl
			   << jac_inv(1,0) << " " << jac_inv(1,1) << std::endl
			   << "Hessian inverse [x]:" << std::endl
			   << hess_inv[0](0,0) << " " << hess_inv[0](0,1) << std::endl
			   << hess_inv[0](1,0) << " " << hess_inv[0](1,1) << std::endl
			   << "Hessian inverse [y]:" << std::endl
			   << hess_inv[1](0,0) << " " << hess_inv[1](0,1) << std::endl
			   << hess_inv[1](1,0) << " " << hess_inv[1](1,1) << std::endl << std::endl;
	      jac_and_hess.close();
	      
	      // [zdec] debug
	      std::ofstream debug_stream;
	      debug_stream.open("norm_and_tan.dat", std::ios_base::app);
	      debug_stream << x[0] << " "
			   << x[1] << " "
			   << ni[0] << " "
			   << ni[1] << " "
			   << ti[0] << " "
			   << ti[1] << " "
			   << dnids[0] << " "
			   << dnids[1] << " "
			   << dtids[0] << " "
			   << dtids[1] << " "
			   << d2xids2[0] << " "
			   << d2xids2[1] << std::endl;
	      debug_stream.close();
	      
	      // Fill in the rotation matrix using the new basis
	      fill_in_rotation_matrix_at_node_with_basis(j_node,
							 jac_inv,
							 hess_inv);
	      // [zdec] DEBUG
	      oomph_info << "ELEMENT " << Parent_element_pt
			 << " is updated " << std::endl;
	      break;
	    }

	    // Node is at a corner, rotate to the tangent1-tangent2 coordinates
	  case 2:
	    {
	      // Our new coordinate system:
	      //     (l, s)=(normal component, tangent component)
	      // which we define in terms of basis vectors (rescaled)
	      //     ni=dxi/dl / |n|           <-- Jacobian col 1
	      //     ti=dxi/ds / |t|           <-- Jacobian col 2
	      // and their derivatives
	      //     dnidxj=d/dxj(dxi/dl / |n|) <-- Hessian `col' 1
	      //     dtidxj=d/dxj(dxi/ds / |t|) <-- Hessian `col' 2	      

	      // [zdec] we use i and j for brevity
	      // but it should be alpha & beta
	      // Need to write up how the transformation is done

	      // Storage for our basis and derivatives
	      Vector<double> t1i(2,0.0);
	      Vector<double> t2i(2,0.0);
	      Vector<double> dt1ids1(2,0.0);
	      Vector<double> dt2ids2(2,0.0);

	      // All tensors assumed evaluated on the boundary
	      // Jacobian of inverse mapping
	      DenseMatrix<double> jac_inv(2,2,0.0);
	      // Hessian of mapping
	      Vector<DenseMatrix<double>> hess(2,DenseMatrix<double>(2,2,0.0));
	      // Hessian of inverse mapping
	      Vector<DenseMatrix<double>> hess_inv(2,DenseMatrix<double>(2,2,0.0));
	      
	      // The basis is defined in terms of the two boundary
	      // parametrisation
	      Vector<double> boundary1_coord = {Boundary_coordinate_of_node[j_node][0]};
	      Vector<double> boundary2_coord = {Boundary_coordinate_of_node[j_node][1]};
	      CurvilineGeomObject* boundary1_pt = Nodal_boundary_parametrisation_pt[j_node][0];
	      CurvilineGeomObject* boundary2_pt = Nodal_boundary_parametrisation_pt[j_node][1];
	      Vector<double> x1(2,0.0);
	      Vector<double> x2(2,0.0);
	      Vector<double> dx1ids1(2,0.0);
	      Vector<double> dx2ids2(2,0.0);
	      Vector<double> d2x1ids12(2,0.0);
	      Vector<double> d2x2ids22(2,0.0);

	      // Get positions (debug)
	      boundary1_pt->position(boundary1_coord, x1);
	      boundary2_pt->position(boundary2_coord, x2);
	      // Get tangent vectors
	      boundary1_pt->dposition(boundary1_coord, dx1ids1);
	      boundary2_pt->dposition(boundary2_coord, dx2ids2);
	      // Get second derivatives
	      boundary1_pt->d2position(boundary1_coord, d2x1ids12);
	      boundary2_pt->d2position(boundary2_coord, d2x2ids22);
		
	      double mag_t = 1.0; // sqrt(dxids[0]*dxids[0]+dxids[1]*dxids[1]);
	      // t1i is the first normalised tangent vector
	      t1i[0] = dx1ids1[0]/mag_t;
	      t1i[1] = dx1ids1[1]/mag_t;
	      // Derivative of first tangent (not normalised)
	      dt1ids1[0] = d2x1ids12[0];
	      dt1ids1[1] = d2x1ids12[1];
	      // t2i is the second normalised tangent vector
	      t2i[0] = dx2ids2[0]/mag_t;
	      t2i[1] = dx2ids2[1]/mag_t;
	      // Derivative of first tangent (not normalised)
	      dt2ids2[0] = d2x2ids22[0];
	      dt2ids2[1] = d2x2ids22[1];
	      
	      // Need inverse of mapping to calculate ds/dxi ----------------
	      //   /  dx/dl  dx/ds  \ -1  ___   _1_  /  dy/ds -dx/ds \ .
	      //   \  dy/dl  dy/ds  /     ---   det  \ -dy/dl  dx/dl /
	      //
	      //                          ___  /  dl/dx  dl/dy  \ .  
	      //                          ---  \  ds/dx  ds/dy  /
	      //
	      // Fill out inverse of Jacobian
	      double det = (t1i[0]*t2i[1] - t2i[0]*t1i[1]);
	      jac_inv(0,0) = t2i[1] / det;
	      jac_inv(0,1) =-t2i[0] / det;
	      jac_inv(1,0) =-t1i[1] / det; 
	      jac_inv(1,1) = t1i[0] / det;
	      
	      // Fill out the Hessian
	      // (unneeded -- can calculate the inverse components by hand)
	      for(unsigned alpha=0; alpha<2; alpha++)
	      {
		hess[alpha](0,0) = dt1ids1[alpha];
		// hess[alpha](0,1) = 0.0;
		// hess[alpha](1,0) = 0.0;
		hess[alpha](1,1) = dt2ids2[alpha];
	      }

	      // Fill out inverse of Hessian
	      // H^{-1}abg = J^{-1}ad Hdez J^{-1}eb J^{-1}zg
	      for(unsigned alpha=0; alpha<2; alpha++)
	      {
		for(unsigned beta=0; beta<2; beta++)
		{
		  for(unsigned gamma=0; gamma<2;gamma++)
		  {
		    for(unsigned alpha2=0; alpha2<2; alpha2++)
		    {
		      for(unsigned beta2=0; beta2<2; beta2++)
		      {
			for(unsigned gamma2=0; gamma2<2;gamma2++)
			{
			  hess_inv[alpha](beta,gamma) -=
			  jac_inv(alpha,alpha2) * hess[alpha2](beta2,gamma2)
			  * jac_inv(beta2,beta) * jac_inv(gamma2,gamma);
			} 
		      } 
		    }		  
		  } 
		} 
	      }		
	      
	      // Fill in the rotation matrix using the new basis
	      fill_in_rotation_matrix_at_node_with_basis(j_node,
							 jac_inv,
							 hess_inv);
	      // [zdec] DEBUG
	      oomph_info << "ELEMENT " << Parent_element_pt
			 << " is updated " << std::endl;
	      break;
	    }
	  }
    
	} // end loop over vertices
      } // end of update_rotation_matrices()


      /// Access function to fill out rot_mat using rotation matrix
      void get_rotation_matrix_at_node(const unsigned &j_node,
				       DenseMatrix<double> &rot_mat)
      {
	rot_mat = Rotation_matrix_at_node[j_node];
      }
      
    private:

      /// Helper function to fill in the rotation matrix for a given basis
      void fill_in_rotation_matrix_at_node_with_basis(const unsigned &j_node,
						      const DenseMatrix<double> &jac_inv,
						      const Vector<DenseMatrix<double>> &hess_inv)
      {
	// Rotation matrix, b constructed using submatrices b1, b12, b22
	DenseMatrix<double> b1(2, 2, 0.0), b22(3, 3, 0.0), b12(2, 3, 0.0);
 
	// Fill in the submatrices
	// Loop over the rotated first derivatives
	for (unsigned mu = 0; mu < 2; mu++)
	{
	  // Loop over the unrotated first derivatives
	  for (unsigned alpha = 0; alpha < 2; alpha++)
	  {
	    // Fill in b1 - the Jacobian
	    // Fill in the affine rotation of the first derivatives
	    b1(mu, alpha) = jac_inv(mu, alpha);
	    
	    // Loop over unrotated second derivatives
	    for (unsigned beta = 0; beta < 2; ++beta)
	    {
	      // Avoid double counting the cross derivative
	      if (alpha <= beta)
	      {
		// Define column index
		const unsigned col = alpha + beta;
		
		// Fill in the non-affine part of the rotation of the first
		// derivatives
		b12(mu, col) += hess_inv[mu](alpha, beta);
		// [zdec] debug mixed derivative -- add extra
		if(alpha<beta)
		{
		  //b12(mu, col) -= hess_inv[mu](alpha, beta);
		}
		// Loop over the rotated second derivatives
		for (unsigned nu = 0; nu < 2; nu++)
		{
		  // // Avoid double counting the cross derivative
		  // if (mu <= nu)
		  {		  
		    // Fill in b22 - the Affine part of the Jacobian derivative
		    // Redefine row index for the next submatrix
		    unsigned row_b22 = mu + nu;
		    // Fill in the affine part of the rotation of the second
		    // derivatives [zdec] if( beta>= alpha) ?
		    b22(row_b22, col) += jac_inv(mu, alpha) * jac_inv(nu,beta);
		  }
		}
	      }
	    }
	  }
	}
	// Fill in the submatrices to the full (6x6) matrix
	Rotation_matrix_at_node[j_node](0, 0) = 1.0;
	// Fill in b1 --- the affine contribution to rotation of the
	// first derivatives
	for (unsigned i = 0; i < 2; ++i)
	{
	  for (unsigned j = 0; j < 2; ++j)
	  {
	    Rotation_matrix_at_node[j_node](1 + i, 1 + j) = b1(i, j);
	  }
	}
	// Fill in b21 --- the non-affine (second derivative dependent)
	// rotation of the first derivatives
	for (unsigned i = 0; i < 2; ++i)
	{
	  for (unsigned j = 0; j < 3; ++j)
	  {
	    Rotation_matrix_at_node[j_node](1 + i, 3 + j) = b12(i, j);
	  }
	}
	// Fill in b22 --- the rotation of the second derivatives
	for (unsigned i = 0; i < 3; ++i)
	{
	  for (unsigned j = 0; j < 3; ++j)
	  {
	    Rotation_matrix_at_node[j_node](3 + i, 3 + j) = b22(i, j);
	  }
	}
      } // end fill_in_rotation_matrix_at_node_with_basis

      /// Pointer to the `parent' finite element which this is a helper force
      FiniteElement* Parent_element_pt;
      
      /// A vector of the number of boundaries associated with each node
      /// --- at most this should be 2 (a left and right side parametrisation)
      Vector<unsigned> Nboundary_at_node;

      /// Vector containing
      /// <vector of <boundary parametrised location of node> for each node>
      Vector<Vector<double>> Boundary_coordinate_of_node;
      
      /// Vector containing <vector of <boundary parametrisations> at each node>
      Vector<Vector<CurvilineGeomObject*>> Nodal_boundary_parametrisation_pt;
      
      /// Vector containing <rotation matrix at each node>
      Vector<DenseMatrix<double>> Rotation_matrix_at_node;
      
    };
    //---end of rotation helper class-------------------------------------------

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

    /// Lagrange interpolated shape for in--plane displacements
    void shape_u(const Vector<double>& s, Shape& psi) const;

    /// Lagrange interpolated d_shape for in--plane displacements
    void dshape_u_local(const Vector<double>& s,
                        Shape& psi,
                        DShape& dpsids) const;

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

    //----------------------------------------------------------------------------
    // Interface to FoepplVonKarmanEquations (can this all be (static) data?)

    /// Field indices for u
    virtual Vector<unsigned> u_field_indices() const
    {
      return {0, 1};
    }

    /// Field indices for u
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
      // [zdec] debug
      oomph_info << "We have "
		 << CurvableBellElement<NNODE_1D>::ninternal_basis_type_for_field(
		   w_field_index()) << " in here." << std::endl;

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


  protected:
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

    // End of FoepplVonKarmanEquations interface functions
    //----------------------------------------------------------------------------

  public:
    /// Function to pin all deflection dofs
    void pin_all_deflection_dofs() const;

    /// Function to pin particular in--plane displacement dofs
    void fix_in_plane_displacement_dof(const unsigned& dof_number,
                                       const unsigned& b,
                                       const ScalarFctPt& u);

    /// Function to pin particular out--of--plane displacement dofs
    void fix_out_of_plane_displacement_dof(const unsigned& dof_number,
                                           const unsigned& b,
                                           const ScalarFctPt& w);

    // Is this element curved?
    bool element_is_curved() const
    {
      return CurvableBellElement<NNODE_1D>::element_is_curved();
    }

    /// Access function to rotated boundary helper object
    RotatedBoundaryHelper* rotated_boundary_helper_pt()
    {
      return Rotated_boundary_helper_pt;
    }

  protected:

    /// Pointer to an instance of rotated boundary helper
    RotatedBoundaryHelper* Rotated_boundary_helper_pt;
    
    // // [zdec] Old rotation -- delete
    // /// Get rotation matrices that change the degrees of freedom to the basis
    // /// set by Boundary_parametrisation_pt
    // inline void rotation_matrix_at_node(const unsigned& inode,
    // 					DenseDoubleMatrix& rotation_matrix) const;

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

  public:
    // // [zdec] old rotation
    // /// Set up the rotated degrees of freedom
    // inline void set_up_rotated_dofs(
    //   const unsigned& nnodes_to_rotate,
    //   const Vector<unsigned>& nodes_to_rotate,
    //   const BasisVectorsFctPt& basis_vectors_fct_pt);

    /// Upgrade the Bell element to a curved Bernadou element
    virtual void upgrade_element_to_curved(
      const MyC1CurvedElements::Edge& curved_edge,
      const double& s_ubar,
      const double& s_obar,
      CurvilineGeomObject* parametric_edge,
      const unsigned& boundary_order)
    {
      CurvableBellElement<NNODE_1D>::upgrade_element_to_curved(
        curved_edge, s_ubar, s_obar, parametric_edge, boundary_order);
    }

    /// Required  # of `values' (pinned or dofs)
    /// at node n
    inline unsigned required_nvalue(const unsigned& n) const
    {
      return Initial_Nvalue[n];
    }

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


    // Private Data Members

  private:
    /// Static number of fields (is always 3)
    static const unsigned Nfield;

    /// Static bool vector with the Bell interpolation of the fields
    /// (always only w)
    static const std::vector<bool> Field_is_bell_interpolated;

    /// Static int that holds the number of variables at
    /// nodes: always the same
    static const unsigned Initial_Nvalue[];

    /// Enum to store which edge is curved set to none when element has no
    /// curved edges
    MyC1CurvedElements::Edge Curved_edge;

    /// Pointer to Stored Association matrix
    DenseMatrix<double>* Association_matrix_pt;

    // // [zdec] old rotation
    // /// A Pointer to the function that sets up the rotated basis at point x
    // BasisVectorsFctPt Rotated_basis_fct_pt;
    
    /// Which nodes are we rotating
    Vector<unsigned> Nodes_to_rotate;

    /// Number of nodes to rotate
    unsigned Nnodes_to_rotate;
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


  // //[zdec] old rotation
  // //==============================================================================
  // /// Set up the rotated degrees of freedom: includes a check for the number of
  // /// rotation nodes being greater than three.
  // //==============================================================================
  // template<unsigned NNODE_1D>
  // void FoepplVonKarmanC1CurvableBellElement<NNODE_1D>::set_up_rotated_dofs(
  //   const unsigned& nnodes_to_rotate,
  //   const Vector<unsigned>& nodes_to_rotate,
  //   const BasisVectorsFctPt& basis_vectors_fct_pt)
  // {
  //   // Change the member Nnode_to_rotate
  //   Nnodes_to_rotate = nnodes_to_rotate;
  //   #ifdef PARANOID
  //   // Check that the number of nodes is smaller than 3
  //   if (nnodes_to_rotate > 3)
  //   {
  //     throw OomphLibError(
  // 	"There are only three nodes per element, so we cannot rotate more than three ",
  // 	OOMPH_CURRENT_FUNCTION,
  // 	OOMPH_EXCEPTION_LOCATION);
  //   }
  //   #endif

  //   Nodes_to_rotate = nodes_to_rotate;

  //   // Point to the basis vectors function
  //   Rotated_basis_fct_pt = basis_vectors_fct_pt;
  // }

  
  // // [zdec] old rotation -- delete
  // //==============================================================================
  // /// Rotate the shape functions according to
  // /// w.r.t. global coordinates and return Jacobian of mapping.
  // ///
  // /// Galerkin: Test functions = shape functions
  // //==============================================================================
  // template<unsigned NNODE_1D>
  // void FoepplVonKarmanC1CurvableBellElement<NNODE_1D>::rotation_matrix_at_node(
  //   const unsigned& inode, DenseDoubleMatrix& rotation_matrix) const
  // {
  //   // Initialise x normal and tangent
  //   Vector<double> x(2, 0.0);

  //   // Get the node pointer
  //   Node* nod_pt = this->node_pt(inode);

  //   // Get the position of the vertex
  //   x[0] = nod_pt->x(0);
  //   x[1] = nod_pt->x(1);

  //   // Initialise the two basis vectors
  //   Vector<Vector<double>> bi(2, Vector<double>(2, 0.0));
  //   Vector<DenseMatrix<double>> dbi(2, DenseMatrix<double>(2, 2, 0.0));

  //   (*Rotated_basis_fct_pt)(x,bi[0],bi[1],dbi[0],dbi[1]); 
        
  //   // Rotation matrix, B
  //   DenseMatrix<double> b1(2, 2, 0.0), b22(3, 3, 0.0), b21(3, 2, 0.0);

  //   // Fill in the submatrices
  //   for (unsigned alpha = 0; alpha < 2; ++alpha)
  //   {
  //     for (unsigned beta = 0; beta < 2; ++beta)
  //     {
  //       // Fill in b1 - the Jacobian
  //       // Fill in the rotation of the first derivatives
  //       b1(alpha, beta) = bi[beta][alpha];

  //       // Avoid double counting the cross derivative
  //       if (alpha <= beta)
  //       {
  //         // Define row index
  //         const unsigned row = alpha + beta;
  //         for (unsigned gamma = 0; gamma < 2; ++gamma)
  //         {
  //           // Fill in b21 - the non affine part of the Jacobian derivative
  //           // Define column index
  //           unsigned col_b21 = gamma;
  //           // Fill in the non-affine part of the rotation of the second
  //           // derivatives if( beta>= alpha) ?
  //           b21(row, col_b21) += dbi[gamma](alpha, beta);
  //           for (unsigned delta = 0; delta < 2; ++delta)
  //           {
  //             // Fill in b22 - the Affine part of the Jacobian derivative
  //             // Redefine column index for the next submatrix
  //             unsigned col_b22 = gamma + delta;
  //             // Fill in the affine part of the rotation of the second
  //             // derivatives if( beta>= alpha) ?
  //             b22(row, col_b22) += bi[gamma][alpha] * bi[delta][beta];
  //           }
  //         }
  //       }
  //     }
  //   }


  //   // Fill in the submatrices to the full (6x6) matrix - we need to right
  //   // multiply this matrix so we need the transpose of the Jacobian W dof
  //   // remains the same
  //   rotation_matrix(0, 0) = 1.0;
  //   // Fill in b1
  //   for (unsigned i = 0; i < 2; ++i)
  //   {
  //     for (unsigned j = 0; j < 2; ++j)
  //     {
  //       rotation_matrix(1 + j, 1 + i) = b1(i, j);
  //     }
  //   }
  //   // Fill in b21
  //   for (unsigned i = 0; i < 3; ++i)
  //   {
  //     for (unsigned j = 0; j < 2; ++j)
  //     {
  //       rotation_matrix(1 + j, 3 + i) = b21(i, j);
  //     }
  //   }
  //   // Fill in b22
  //   for (unsigned i = 0; i < 3; ++i)
  //   {
  //     for (unsigned j = 0; j < 3; ++j)
  //     {
  //       rotation_matrix(3 + j, 3 + i) = b22(i, j);
  //     }
  //   }
  // }


  //============================================================================
  /// Fetch the in-plane basis functions and their derivatives w.r.t. global
  /// coordinates at s and return Jacobian of mapping.
  //=============================================================================
  template<unsigned NNODE_1D>
  double FoepplVonKarmanC1CurvableBellElement<
    NNODE_1D>::dbasis_u_eulerian_foeppl_von_karman(const Vector<double>& s,
                                                   Shape& psi_n,
                                                   DShape& dpsi_n_dx) const
  {
    // Initialise and get dpsi w.r.t local coord
    const unsigned n_node = this->nnode();
    DShape dpsids(n_node, this->dim());
    dshape_u_local(s, psi_n, dpsids);
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
    dshape_u_local(s, psi_n, dpsids);
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
    throw OomphLibError("This still needs testing for curved elements.",
                        "void FoepplVonKarmanC1CurvableBellElement<NNODE_1D>::\
shape_and_test_foeppl_von_karman(...)",
                        OOMPH_EXCEPTION_LOCATION);

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
    // Throw if called
    throw OomphLibError("This still needs testing for curved elements.",
                        "void FoepplVonKarmanC1CurvableBellElement<NNODE_1D>::\
dshape_and_dtest_foeppl_von_karman(...)",
                        OOMPH_EXCEPTION_LOCATION); // HERE

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
  /// Rotate the shape functions according to the specified basis on the
  /// boundary. This function does a DenseDoubleMatrix solve to determine
  /// new basis - which could be speeded up by caching the matrices higher
  /// up and performing the LU decomposition only once
  //======================================================================
  template<unsigned NNODE_1D>
  inline void FoepplVonKarmanC1CurvableBellElement<NNODE_1D>::rotate_shape(
    Shape& psi) const
  {
    const unsigned n_dof_types = nw_type_at_each_node();

    // Get the nodes that need rotating
    Vector<unsigned> nodes_to_rotate;
    for(unsigned j_node=0; j_node<3; j_node++)
    {
      if(Rotated_boundary_helper_pt->nboundary_at_node(j_node)>0)
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
      this->Rotated_boundary_helper_pt
      ->get_rotation_matrix_at_node(j_node, rotation_matrix);

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
    for(unsigned j_node=0; j_node<3; j_node++)
    {
      if(Rotated_boundary_helper_pt->nboundary_at_node(j_node)>0)
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
      this->Rotated_boundary_helper_pt
      ->get_rotation_matrix_at_node(j_node, rotation_matrix);

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
    const unsigned n_2ndderiv =
    ((n_dim + 1) * (n_dim)) / 2; 

    // Get the nodes that need rotating
    Vector<unsigned> nodes_to_rotate;
    for(unsigned j_node=0; j_node<3; j_node++)
    {
      if(Rotated_boundary_helper_pt->nboundary_at_node(j_node)>0)
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
      this->Rotated_boundary_helper_pt
	->get_rotation_matrix_at_node(j_node, rotation_matrix);
      
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


  //=============================================================================
  /// Function to pin all deflection dofs
  //=============================================================================
  template<unsigned NNODE_1D>
  void FoepplVonKarmanC1CurvableBellElement<NNODE_1D>::pin_all_deflection_dofs()
    const
  {
    const unsigned w_index = w_field_index();

    // Get nodes
    const unsigned n_node = nw_node();
    const Vector<unsigned> nodes = get_w_node_indices();

    // Curved Bell elements only have deflection dofs at vertices
    for (unsigned j_nodei = 0; j_nodei < n_node; j_nodei++)
    {
      // Get the j_nodei-th node used by i_field
      unsigned j_node = nodes[j_nodei];
      Node* nod_pt = this->node_pt(j_node);

      // Get the number of types at the current node
      unsigned n_type = nw_type_at_each_node();

      // Check if it is on the boundary
      for (unsigned k_type = 0; k_type < n_type; k_type++)
      {
        // Pin and set the value
        nod_pt->pin(2 + k_type);
        nod_pt->set_value(2 + k_type, 0.0);
      }
    }
    
    // Now fix internal dofs
    unsigned n_internal = nw_type_internal();
    for (unsigned k_type = 0; k_type < n_internal; k_type++)
    {
      // Get node
      // Pin and set the value
      this->internal_data_for_field_pt(w_index)->pin(k_type);
    }
  }


  //=============================================================================
  /// Function to pin particular out-of-plane displacement dofs
  //=============================================================================
  template<unsigned NNODE_1D>
  void FoepplVonKarmanC1CurvableBellElement<
    NNODE_1D>::fix_out_of_plane_displacement_dof(const unsigned& j_type,
                                                 const unsigned& b_boundary,
                                                 const ScalarFctPt&
						 specified_deflection_fct_pt)
  {
    const unsigned w_index = w_field_index();
    const unsigned first_nodal_type_index =
    this->first_nodal_type_index_for_field(w_index);
    const unsigned n_vertices = nw_node();
    const unsigned n_type = nw_type_at_each_node();

    #ifdef PARANOID
    // Check that the dof number is a sensible value
    if (j_type >= n_type)
    {
      throw OomphLibError(
        "Foppl von Karman elements only have 6 Hermite deflection degrees\
of freedom at internal points. They are {w ; w,x ; w,y ; w,xx ; w,xy ; w,yy}",
        "FoepplVonKarmanC1CurvableBellElement:fix_out_of_plane_dof()",
        OOMPH_EXCEPTION_LOCATION);
    }
    #endif

    // Bell elements only have deflection dofs at vertices
    for (unsigned n = 0; n < n_vertices; ++n)
    {
      // Get node
      Node* nod_pt = this->node_pt(n);
      // Check if it is on the boundary
      bool is_boundary_node = nod_pt->is_on_boundary(b_boundary);
      if (is_boundary_node)
      {
        // Extract nodal coordinates from node:
        Vector<double> x(2);
        x[0] = nod_pt->x(0);
        x[1] = nod_pt->x(1);
        // Get value
        double value;
        specified_deflection_fct_pt(x, value);
        // Pin and set the value
        nod_pt->pin(first_nodal_type_index + j_type);
        nod_pt->set_value(first_nodal_type_index + j_type, value);
      }
    }
  }


  //=============================================================================
  /// Function to pin particular out-of-plane displacement dofs
  //=============================================================================
  template<unsigned NNODE_1D>
  void FoepplVonKarmanC1CurvableBellElement<
    NNODE_1D>::fix_in_plane_displacement_dof(const unsigned& alpha,
                                             const unsigned& b,
                                             const ScalarFctPt&
					     specified_displacement_fct_pt)
  {
    // Initialise constants that we use in this function
    const unsigned field_index = u_alpha_field_index(alpha);
    const unsigned nodal_type_index =
    this->first_nodal_type_index_for_field(field_index);
    const unsigned n_node = nu_node();
    const unsigned n_type = 2 * nu_type_at_each_node();
    const unsigned dim = this->dim();

    #ifdef PARANOID
    // Check that the dof number is a sensible value
    if (alpha >= n_type)
    {
      throw OomphLibError(
        "Foppl von Karman elements only have 2 in-plane displacement degrees\
of freedom at internal points. They are {ux, uy}",
        "FoepplVonKarmanC1CurvableBellElement:fix_out_of_plane_dof()",
        OOMPH_EXCEPTION_LOCATION);
    }
    #endif

    // Bell elements only have deflection dofs at vertices
    for (unsigned n = 0; n < n_node; ++n)
    {
      // Get node
      Node* nod_pt = this->node_pt(n);
      // Check if it is on the boundary
      bool is_boundary_node = nod_pt->is_on_boundary(b);
      if (is_boundary_node)
      {
        // Extract nodal coordinates from node:
        // Since the element isn't necessarily isoparametric the nodes
        // 'position' is not necessarily correct?
        Vector<double> x(dim), s(dim);
        this->local_coordinate_of_node(n, s);
        interpolated_x(s, x);
        // Fill in value
        double value;
        specified_displacement_fct_pt(x, value);
        // Pin and set the value
        nod_pt->pin(alpha);
        nod_pt->set_value(nodal_type_index, value);
      }
    }
  }
} // namespace oomph


#endif

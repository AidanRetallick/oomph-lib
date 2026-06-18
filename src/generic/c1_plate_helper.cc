#include "c1_plate_helper.h"
#include "subparametric_Telement.h"

namespace oomph
{

 
 /// Static double to define the default tolerance with which the
 /// Newton method in get_zeta(...)
 /// determines the local coordinate zeta associated with a
 /// specified point on the curviline.
 double C1CurviLine::Tol_for_get_zeta=1.0e-12;
 

 //===============================================================================
 /// Document the pin status of the Lagrange multipliers
 //===============================================================================
 void DuplicateNodeConstraintElement::doc_pin_status_of_lagrange_multipliers(
  std::ostream& outfile)
 {
  outfile << Right_node_pt->x(0) << " " << Right_node_pt->x(1) << " : ";
   Data* lagr_data_pt=internal_data_pt(Index_of_lagrange_data);
   if (lagr_data_pt!=0)
    {
     unsigned nval=lagr_data_pt->nvalue();
     for (unsigned i=0;i<nval;i++)
      {
       if (lagr_data_pt->is_pinned(i))
        {
         outfile << " 1 ";
        }
       else
        {
         outfile << " 0 ";
        }
      }
     outfile << std::endl;
    }
  }

  

 //===============================================================================
 /// Function to calculate Jacobian and Hessian of the coordinate mapping
 //===============================================================================
 void DuplicateNodeConstraintElement::get_jac_and_hess_of_coordinate_transform(
  DenseMatrix<double>& jac_of_transform,
  Vector<DenseMatrix<double>>& hess_of_transform)
 {
  //----------------------------------------------------------------------
   // We need the parametrisations either side of the vertex which define
   // the coordinates each node uses for its Hermite dofs.
   Vector<double> left_x(2, 0.0); // [zdec] debug
   Vector<double> right_x(2, 0.0); // [zdec] debug
   Vector<double> left_dxids(2, 0.0);
   Vector<double> left_d2xids2(2, 0.0);
   Vector<double> right_dxids(2, 0.0);
   Vector<double> right_d2xids2(2, 0.0);
   Left_boundary_pt->position(Left_node_coord, left_x); // [zdec] debug 
   Right_boundary_pt->position(Right_node_coord, right_x); // [zdec] debug
   Left_boundary_pt->dposition(Left_node_coord, left_dxids);
   Left_boundary_pt->d2position(Left_node_coord, left_d2xids2);
   Right_boundary_pt->dposition(Right_node_coord, right_dxids);
   Right_boundary_pt->d2position(Right_node_coord, right_d2xids2);

   // Get the speed of each parametrisation
   double left_mag =
    sqrt(left_dxids[0] * left_dxids[0] + left_dxids[1] * left_dxids[1]);
   double right_mag =
    sqrt(right_dxids[0] * right_dxids[0] + right_dxids[1] * right_dxids[1]);

   //----------------------------------------------------------------------
   // Normalise dxids to find the tangent vectors and their
   // derivatives either side of the vertex
   Vector<double> left_ti(2, 0.0);
   Vector<double> left_ni(2, 0.0);
   Vector<double> left_dtids(2, 0.0);
   Vector<double> left_dnids(2, 0.0);
   Vector<double> right_ti(2, 0.0);
   Vector<double> right_ni(2, 0.0);
   Vector<double> right_dtids(2, 0.0);
   Vector<double> right_dnids(2, 0.0);
   for (unsigned alpha = 0; alpha < 2; alpha++)
    {
     // Fill in the tangents either side of the vertex
     left_ti[alpha] = left_dxids[alpha] / left_mag;
     right_ti[alpha] = right_dxids[alpha] / right_mag;
     // Fill in the derivatives of the (normalised) tangents either side of
     // the vertex
     left_dtids[alpha] =
      left_d2xids2[alpha] / std::pow(left_mag, 2) -
      (left_dxids[0] * left_d2xids2[0] + left_dxids[1] * left_d2xids2[1]) *
      left_dxids[alpha] / std::pow(left_mag, 4);
     right_dtids[alpha] = right_d2xids2[alpha] / std::pow(right_mag, 2) -
      (right_dxids[0] * right_d2xids2[0] +
       right_dxids[1] * right_d2xids2[1]) *
      right_dxids[alpha] / std::pow(right_mag, 4);
     // Use these to fill out the corresponding vectors for the normal
     // direction (nx,ny) = (ty,-tx)
    }
   // Use orthogonality to fill in normals and their derivatives
   for (unsigned alpha = 0; alpha < 2; alpha++)
    {
     left_ni[alpha] = pow(-1, alpha) * left_ti[(alpha + 1) % 2];
     right_ni[alpha] = pow(-1, alpha) * right_ti[(alpha + 1) % 2];
     left_dnids[alpha] = pow(-1, alpha) * left_dtids[(alpha + 1) % 2];
     right_dnids[alpha] = pow(-1, alpha) * right_dtids[(alpha + 1) % 2];
    }

   //----------------------------------------------------------------------
   // We need to fill out the Jacobians and Hessians of the boundary
   // coordinates either side of the vertex
   DenseMatrix<double> left_jac(2, 2, 0.0);
   DenseMatrix<double> right_jac(2, 2, 0.0);
   Vector<DenseMatrix<double>> left_hess(2, DenseMatrix<double>(2, 2, 0.0));
   Vector<DenseMatrix<double>> right_hess(2, DenseMatrix<double>(2, 2, 0.0));
   for (unsigned alpha = 0; alpha < 2; alpha++)
    {
     // Fill in Jacobians {{nx,tx},{ny,ty}}
     left_jac(alpha, 0) = left_ni[alpha];
     left_jac(alpha, 1) = left_ti[alpha];
     right_jac(alpha, 0) = right_ni[alpha];
     right_jac(alpha, 1) = right_ti[alpha];
     // Fill in Hessians
     // left_hess[alpha](0,0) = 0.0;
     left_hess[alpha](0, 1) = left_dnids[alpha];
     left_hess[alpha](1, 0) = left_dnids[alpha];
     left_hess[alpha](1, 1) = left_dtids[alpha];
     // right_hess[alpha](0,0) = 0.0;
     right_hess[alpha](0, 1) = right_dnids[alpha];
     right_hess[alpha](1, 0) = right_dnids[alpha];
     right_hess[alpha](1, 1) = right_dtids[alpha];
    }

   //----------------------------------------------------------------------
   // We need the inverse Jacobian and Hessian for the left parametrisation
   DenseMatrix<double> left_jac_inv(2, 2, 0.0);
   Vector<DenseMatrix<double>> left_hess_inv(2,
                                             DenseMatrix<double>(2, 2, 0.0));
   left_jac_inv(0, 0) = left_jac(1, 1);
   left_jac_inv(0, 1) = -left_jac(0, 1);
   left_jac_inv(1, 0) = -left_jac(1, 0);
   left_jac_inv(1, 1) = left_jac(0, 0);
   // Fill out inverse of Hessian
   // H^{-1}abg = J^{-1}ad Hdez J^{-1}eb J^{-1}zg
   for (unsigned alpha = 0; alpha < 2; alpha++)
    {
     for (unsigned beta = 0; beta < 2; beta++)
      {
       for (unsigned gamma = 0; gamma < 2; gamma++)
        {
         for (unsigned alpha2 = 0; alpha2 < 2; alpha2++)
          {
           for (unsigned beta2 = 0; beta2 < 2; beta2++)
            {
             for (unsigned gamma2 = 0; gamma2 < 2; gamma2++)
              {
               left_hess_inv[alpha](beta, gamma) -=
                left_jac_inv(alpha, alpha2) *
                left_hess[alpha2](beta2, gamma2) *
                left_jac_inv(beta2, beta) * left_jac_inv(gamma2, gamma);
              }
            }
          }
        }
      }
    }

   //----------------------------------------------------------------------
   //----------------------------------------------------------------------
   // Use these to calculate the Jacobian of the left->right transform
   //     J = J_{left}^{-1}J_{right}
   // and the Hessian of the left->right transform
   //     H = H_{left}^{-1}J_{right}J_{right} + J_{left}^{-1}H_{right}
   for (unsigned alpha = 0; alpha < 2; alpha++)
    {
     for (unsigned beta = 0; beta < 2; beta++)
      {
       for (unsigned gamma = 0; gamma < 2; gamma++)
        {
         // Add contribution to J
         jac_of_transform(alpha, beta) +=
          left_jac_inv(alpha, gamma) * right_jac(gamma, beta);
         for (unsigned mu = 0; mu < 2; mu++)
          {
           // Add second term contribution to H
           hess_of_transform[alpha](beta, gamma) +=
            left_jac_inv(alpha, mu) * right_hess[mu](beta, gamma);
           for (unsigned nu = 0; nu < 2; nu++)
            {
             // Add first term contribution to H
             hess_of_transform[alpha](beta, gamma) +=
              left_hess_inv[alpha](mu, nu) * right_jac(mu, beta) *
              right_jac(nu, gamma);
            }
          }
        }
      }
    }
  } // End get_jac_and_hess_of_coordinate_transform



 
//==============================================================================
/// Namespace to deal update triangle meshes to deal with C1 elements
//==============================================================================
namespace C1PlateHelper
{
 

 
 /// Boolean to suppress warning about polygonal boundaries
 bool Do_not_warn_about_polygonal_boundaries=false;
  
 /// Output stream to document split elements
 std::ofstream Split_elements_output_stream;
     
 /// Output stream to document duplicated nodes
 std::ofstream Duplicated_node_output_stream;
     
 /// Output stream to document elements whose boundaries have been
 /// upgraded to become curved
 std::ofstream Upgraded_to_curved_edge_element_stream;

 /// Output stream to document nodes at which derivative dofs have been
 /// adjusted to represent derivatives in the normal and tantential direction
 /// We're putputting x,y of the node; the number of the bulk element that uses
 /// the information and the boundary ID of the curvilinear boundary
 std::ofstream Rotated_node_output_stream;
 
 /// Output stream to document elements that contain nodes at which
 /// derivative dofs have been adjusted to represent derivatives in
 /// the normal and tantential direction
 std::ofstream Rotated_element_output_stream;

 /// Boolean vector to keep track of which boundaries are curvilinear
 /// (vs polygonal)
 std::vector<bool> Boundary_is_curvilinear;

 
//==============================================================================
/// Fct to duplicate nodes that span two curved boundaries of the triangle mesh
/// pointed to by bulk_mesh_pt. This makes sure that each node on the boundary 
/// is only  associated with a single (smooth) boundary. Smooth boundaries  of
/// the mesh/domain are available  via the map c1_curviline_pt. The required 
/// continuity of the solution is then imposed by a suitable
/// DuplicateNodeConstraintElement that is created and added to the
/// Mesh pointed to by constraint_mesh_pt.
//==============================================================================
 void duplicate_corner_nodes(Mesh* bulk_mesh_pt, 
                             const std::map<unsigned,C1CurviLine*>& c1_curviline_pt,
                             Mesh* constraint_mesh_pt)
 {

#ifdef PARANOID
#ifdef OOMPH_HAS_MPI
      // Unlikely to work for distributed meshes
      if (bulk_mesh_pt->is_mesh_distributed())
       {
        throw OomphLibError("C1PlateHelper::duplicate_corner_nodes(...) "
                            "is unlikely to work for distributed meshes.\n",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
       }
#endif
#endif

      
  // Storage for collection of nodes that occupy two boundaries together
  // with the boundary IDs (ordered: first < second)
  std::map<Node*,std::pair<unsigned,unsigned>> boundaries_of_boundary_node_pt;

  // Map new_node_pt = newly_created_node_pt[old_node_pt] for replacing
  // nodes in elements
  std::map<Node*,Node*> newly_created_node_pt;

  // Loop over the curvilinear parts of the boundary
  for (const auto& [i_bound, para] : c1_curviline_pt)
  {
   unsigned n_b_node = bulk_mesh_pt->nboundary_node(i_bound);
   for(unsigned i_b_node = 0; i_b_node < n_b_node; i_b_node++)
    {
     // Store the node we are checking
     Node* node_pt = bulk_mesh_pt->boundary_node_pt(i_bound,i_b_node);

     // Pointer to set that contains the boundaries we're on
     std::set<unsigned>* boundaries_pt=0;
     node_pt->get_boundaries_pt(boundaries_pt);
     if (boundaries_pt!=0)
      {
       if (boundaries_pt->size()==2)
        {
         unsigned b_min=UINT_MAX;
         unsigned b_max=0;
         for (unsigned b : (*boundaries_pt))
          {
           if (b<b_min) b_min=b;
           if (b>b_max) b_max=b;           
          }
         // Ordered!
         boundaries_of_boundary_node_pt[node_pt].first=b_min;
         boundaries_of_boundary_node_pt[node_pt].second=b_max;
        }
       else if (boundaries_pt->size()>2)
        {         
         std::stringstream error_message;
         error_message
          << "Sorry; can't currently handle cases where a node is\n"
          << "located on more than two boundaries. Not rocket science but\n"
          << "but somebody needs to implement it"
          << std::endl;
         throw OomphLibError(
          error_message.str(),
          OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
        }
      }
    }
  }
  
  // Here are the nodes that need to be duplicated. We duplicate them
  // on the lower of its two boundaries (this is stored first)
  for (const auto& a : boundaries_of_boundary_node_pt)
   {
    Node* node_to_be_duplicated_pt=a.first;
    unsigned boundary_on_which_node_is_duplicated=a.second.first;
    unsigned boundary_on_which_node_is_left=a.second.second;

    // Doc?
    if (Duplicated_node_output_stream.is_open())
     {
      Duplicated_node_output_stream
       << node_to_be_duplicated_pt->x(0) << " "
       << node_to_be_duplicated_pt->x(1) << " "
       << boundary_on_which_node_is_duplicated << " "
       << boundary_on_which_node_is_left
       << std::endl;
     }
    
    // Find a boundary element that contains the node to be duplicated on
    // the boundary where the node is to be duplicated
    FiniteElement* el_where_node_is_to_be_duplicated_pt=0;

    // Note: in this context we do not need to include the elements 
    // for which the node is the only one on the boundary.
                                                         
    // Store elements and local node number within them
    // that need to have node replaced
    std::set<std::pair<FiniteElement*,unsigned>> node_to_be_replaced_set;
    bool identified_element_that_will_make_the_new_node=false;
    unsigned n_b_el = bulk_mesh_pt->nboundary_element(
     boundary_on_which_node_is_duplicated);
    for (unsigned i_b_el = 0; i_b_el < n_b_el; i_b_el++)
    {
      // Get the element pointer
      FiniteElement* el_pt = bulk_mesh_pt->boundary_element_pt
       (boundary_on_which_node_is_duplicated, i_b_el);


      // oomph_info << "Testing if element " << i_b_el << " of "
      //            << n_b_el << " : " <<  el_pt << " with vertices at \n"
      //            << el_pt->node_pt(0)->x(0) << " "
      //            << el_pt->node_pt(0)->x(1) << "\n"
      //            << el_pt->node_pt(1)->x(0) << " "
      //            << el_pt->node_pt(1)->x(1) << "\n"
      //            << el_pt->node_pt(2)->x(0) << " "
      //            << el_pt->node_pt(2)->x(1) << "\n"
      //            << "contains the to be replaced node "
      //            << node_to_be_duplicated_pt << " on boundary "
      //            << boundary_on_which_node_is_duplicated
      //            << std::endl;
      
      // If the corner node pt is in the element we have found the right
      // element
      int j_local=el_pt->get_node_number(node_to_be_duplicated_pt);
      if (j_local != -1)
       {
        if (identified_element_that_will_make_the_new_node)
         {
          //oomph_info << "    YES: replace but don't build " << std::endl;
          node_to_be_replaced_set.insert(std::make_pair(el_pt,j_local));
         }
        else
         {
          //oomph_info << "    YES: replace and build " << std::endl;
          el_where_node_is_to_be_duplicated_pt = el_pt;
          identified_element_that_will_make_the_new_node=true;
         }
      }
    }

    // Get timestepper of original node
    TimeStepper* time_stepper_pt=node_to_be_duplicated_pt->time_stepper_pt();
     
    // Now we need to create a new node and substitute the element's
    // old corner node for this new one:
    unsigned local_node_number_of_to_be_replaced_node=
     el_where_node_is_to_be_duplicated_pt->
     get_node_number(node_to_be_duplicated_pt);
    
    Node* hierher_old_node_pt=el_where_node_is_to_be_duplicated_pt->
     node_pt(local_node_number_of_to_be_replaced_node);
    
    Node* new_node_pt = el_where_node_is_to_be_duplicated_pt->
     construct_boundary_node(local_node_number_of_to_be_replaced_node,
                             time_stepper_pt);

    // oomph_info << "Dealing with boundaries (dupl,left) "
    //            << boundary_on_which_node_is_duplicated << " "
    //            << boundary_on_which_node_is_left << " "
    //            << " Replacing node " << local_node_number_of_to_be_replaced_node
    //            << " " << hierher_old_node_pt
    //            << " in element " << el_where_node_is_to_be_duplicated_pt
    //            << " with vertices at \n"
    //            << el_where_node_is_to_be_duplicated_pt->node_pt(0)->x(0) << " "
    //            << el_where_node_is_to_be_duplicated_pt->node_pt(0)->x(1) << "\n"
    //            << el_where_node_is_to_be_duplicated_pt->node_pt(1)->x(0) << " "
    //            << el_where_node_is_to_be_duplicated_pt->node_pt(1)->x(1) << "\n"
    //            << el_where_node_is_to_be_duplicated_pt->node_pt(2)->x(0) << " "
    //            << el_where_node_is_to_be_duplicated_pt->node_pt(2)->x(1) << "\n"
    //            << " with new node " << new_node_pt << std::endl;
    
    // Copy the position and other info from the old node into the new node
    new_node_pt->x(0)=node_to_be_duplicated_pt->x(0);
    new_node_pt->x(1)=node_to_be_duplicated_pt->x(1);

    // Then we add this node to the mesh
    bulk_mesh_pt->add_node_pt(new_node_pt);

    // Now overwrite the pointer to the old node with the newly
    // created one in all other boundary elements that shared it
    for (auto [fe_pt, j_local] : node_to_be_replaced_set)
     {
      // oomph_info << "Also replacing node " << j_local << " : "
      //            << fe_pt->node_pt(j_local) << " in element "
      //            << fe_pt
      //            << " with vertices at \n"
      //            << fe_pt->node_pt(0)->x(0) << " "
      //            << fe_pt->node_pt(0)->x(1) << "\n"
      //            << fe_pt->node_pt(1)->x(0) << " "
      //            << fe_pt->node_pt(1)->x(1) << "\n"
      //            << fe_pt->node_pt(2)->x(0) << " "
      //            << fe_pt->node_pt(2)->x(1) << "\n"
      //            << "with the new node " << new_node_pt
      //            << std::endl;
       fe_pt->node_pt(j_local)=new_node_pt;
     }


    
    // Keep track sowe can update the elements
    newly_created_node_pt[node_to_be_duplicated_pt] = new_node_pt;
    
    // // hierher debug
    // {
    //  oomph_info << "hierher BEFORE" << std::endl;
    //  oomph_info << "to be duplicated: " << node_to_be_duplicated_pt << "\n"
    //             << "new             : " << new_node_pt
    //             << std::endl;

    //  std::set<Node*> nod_pt_set={node_to_be_duplicated_pt,new_node_pt};
    //  for (Node* nod_pt : nod_pt_set)
    //   {
    //    // oomph_info <<  "Node " << nod_pt << " at "
    //    //            << nod_pt->x(0) << " "
    //    //            << nod_pt->x(1) << " " 
    //    //            << " is on boundaries: ";
    //    std::set<unsigned>* bnd_set_pt=0;
    //    nod_pt->get_boundaries_pt(bnd_set_pt);       
    //    if (bnd_set_pt!=0)
    //     {
    //      for (const unsigned& b : *bnd_set_pt)
    //       {
    //        //oomph_info << b << " ";
    //        if (nod_pt->boundary_coordinates_have_been_set_up())
    //         {
    //          Vector<double> boundary_zeta(1);
    //          nod_pt->get_coordinates_on_boundary(b, boundary_zeta);
    //          //oomph_info << boundary_zeta[0] << " ; "; 
    //         }
    //        else
    //         {
    //          //oomph_info << " @@@ "  << " ; "; 
    //         }
    //       }
    //     }
    //    //oomph_info << std::endl;
    //   }
    // }



    // Copy across the boundary coordinates
    Vector<double> boundary_zeta_to_be_copied(1);
    bool copy_boundary_coordinate=false;
    if (node_to_be_duplicated_pt->boundary_coordinates_have_been_set_up())
     {
      copy_boundary_coordinate=true;
      node_to_be_duplicated_pt->get_coordinates_on_boundary
       (boundary_on_which_node_is_duplicated, boundary_zeta_to_be_copied);
      // oomph_info << "Getting boundary coordinate on boundary "
      //            << boundary_on_which_node_is_duplicated
      //            << " from " << node_to_be_duplicated_pt << std::endl;
     }
    
    // Then replace the old node for the new one on the boundary
    bulk_mesh_pt->remove_boundary_node(
     boundary_on_which_node_is_duplicated,node_to_be_duplicated_pt);
    bulk_mesh_pt->   add_boundary_node(
     boundary_on_which_node_is_duplicated,new_node_pt);

    // Copy boundary coordinate across
    if (copy_boundary_coordinate)
     {
      // oomph_info << "Adding boundary coordinate on boundary "
      //            << boundary_on_which_node_is_duplicated
      //            << " to " << new_node_pt << std::endl;
      new_node_pt->set_coordinates_on_boundary
       (boundary_on_which_node_is_duplicated, boundary_zeta_to_be_copied);
     }
    
    // // hierher debug
    // {
    //  oomph_info << "hierher AFTER" << std::endl;
    //  oomph_info << "to be duplicated: " << node_to_be_duplicated_pt << "\n"
    //             << "new             : " << new_node_pt
    //             << std::endl;

    //  std::set<Node*> nod_pt_set={node_to_be_duplicated_pt,new_node_pt};
    //  for (Node* nod_pt : nod_pt_set)
    //   {
    //    oomph_info <<  "Node " << nod_pt << " at "
    //               << nod_pt->x(0) << " "
    //               << nod_pt->x(1) << " " 
    //               << " is on boundaries: ";
    //    std::set<unsigned>* bnd_set_pt=0;
    //    nod_pt->get_boundaries_pt(bnd_set_pt);       
    //    if (bnd_set_pt!=0)
    //     {
    //      for (const unsigned& b : *bnd_set_pt)
    //       {
    //        oomph_info << b << " ";
    //        if (nod_pt->boundary_coordinates_have_been_set_up())
    //         {
    //          Vector<double> boundary_zeta(1);
    //          nod_pt->get_coordinates_on_boundary(b, boundary_zeta);
    //          oomph_info << boundary_zeta[0] << " ; "; 
    //         }
    //        else
    //         {
    //          oomph_info << " @@@ "  << " ; "; 
    //         }
    //       }
    //     }
    //    oomph_info << std::endl;
    //   }
    // }

    // The final job is to constrain this duplication using the specialised
    // Lagrange multiplier elements which enforce equality of displacement and
    // its derivatives either side of this corner.
    C1CurviLine* left_parametrisation_pt =
      c1_curviline_pt.at(boundary_on_which_node_is_left);
    C1CurviLine* right_parametrisation_pt =
      c1_curviline_pt.at(boundary_on_which_node_is_duplicated);


    // Get the coordinates on each node on their respective boundaries
    Vector<double> left_boundary_coordinate =
     {left_parametrisation_pt->get_zeta(node_to_be_duplicated_pt->position())};
    Vector<double> right_boundary_coordinate =
     {right_parametrisation_pt->get_zeta(new_node_pt->position())};

    // Create the constraining element using the first bulk element
    // in them mesh. This seems odd/convoluted but is OK, because FvK and
    // KS (both of which inherit from the template-free base class) provide
    // their own specific implementation of the factory that creates the
    // appropriate DuplicateNodeConstraintElement, so we first get a pointer
    // to some (any!) such element...
    TemplateFreeCurvableBellElement* bulk_el_pt=
     dynamic_cast<TemplateFreeCurvableBellElement*>
     (bulk_mesh_pt->element_pt(0));

    // ...and the pointer to the base class then follows on to the specific
    // implementation of the factory function in the specific derived class.
    DuplicateNodeConstraintElement* constraint_element_pt =
     bulk_el_pt->duplicate_constraint_element_factory(node_to_be_duplicated_pt,
                                                      new_node_pt,
                                                      left_parametrisation_pt,
                                                      right_parametrisation_pt,
                                                      left_boundary_coordinate,
                                                      right_boundary_coordinate);
    
    // Add the constraining element to the mesh
    constraint_mesh_pt->add_element_pt(constraint_element_pt);
   }
 }



//==============================================================================
/// A function that upgrades straight sided elements to be curved. This involves
/// Setting up the parametric boundary, F(s) and the first derivative F'(s)
/// We also need to set the edge number of the upgraded element and the positions
/// of the nodes j and k (defined below) and set which edge (k) is to be exterior
///            @ k               .
///           /(                 .
///          /. \                .
///         /._._)               .
///      i @     @ j             .
/// For RESTING or FREE boundaries we need to have a C2 CONTINUOUS boundary
/// representation. That is we need to have a continuous 2nd derivative defined
/// too. This is well discussed in by [Zenisek 1981] (Aplikace matematiky ,
/// Vol. 26 (1981), No. 2, 121--141). This results in the necessity for F''(s)
/// as well. Final optional argument, boundary order can take values 3 and 5 and 
/// represents the order of the polynomial that represents the curved boundary. 
/// Default value of 5 works for all boundary conditions; 3 is faster but only 
/// works for homogeneous clamped boundaries (in a plate context).
/// hierher Aidan: check description of final arg
//=============================================================================
 void upgrade_edge_elements_to_curved_boundaries(
  Mesh* bulk_mesh_pt, 
  const std::map<unsigned,C1CurviLine*>&  c1_curviline_pt,
  const unsigned& boundary_order) 
 {

#ifdef PARANOID
  if ((boundary_order!=3)&&(boundary_order!=5))
   {
    std::stringstream error_message;
    error_message << "You can only upgrade curved boundaries to polynomial\n"
                  << "approximations of order 3 or 5. You specified "
                  << boundary_order << std::endl;
    throw OomphLibError(
     error_message.str(),
     OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
   }
#endif
  
  // Loop over the curvilinear parts of the outer boundary
  for (const auto& [ibound, c1_curve_pt] : c1_curviline_pt)
   {
    // Loop over the bulk elements adjacent to boundary ibound
    const unsigned n_els=bulk_mesh_pt->nboundary_element(ibound);
    for(unsigned e=0; e<n_els; e++)
     {
      // Get pointer to bulk element adjacent to b
      FiniteElement* bulk_el_pt =
       bulk_mesh_pt->boundary_element_pt(ibound,e);

      // Initialise enum for the curved edge
      C1PlateHelper::CurvedEdgeEnumeration edge=
       C1PlateHelper::CurvedEdgeEnumeration::none;
      
      // Loop over all (three) vertex nodes of the element and
      // identify single node that is interior (i.e. not on the current (!)
      // curved boundary)
      unsigned index_of_interior_node = 3;
      unsigned nnode_not_on_curved_boundary = 0;
      const unsigned nnode = 3;
      Vector<Vector<double> > xn(nnode,Vector<double>(2,0.0));
      for(unsigned n=0;n<nnode;++n)
       {
        Node* nod_pt = bulk_el_pt->node_pt(n);
        xn[n][0]=nod_pt->x(0);
        xn[n][1]=nod_pt->x(1);
        
        // Check if it is on the curved boundaries
        bool node_is_on_curved_boundary=false;
        {
         if (nod_pt->is_on_boundary(ibound))  
          {
           node_is_on_curved_boundary=true;
           // oomph_info << "Node " << nod_pt << " at "
           //            << xn[n][0] << " "
           //            << xn[n][1] << " "
           //            << "is on curved boundary "
           //            << ibound << std::endl;
          }
        }
        if (!node_is_on_curved_boundary)
         {
          index_of_interior_node = n;
          nnode_not_on_curved_boundary++;
          // oomph_info << "Node " << nod_pt << " at "
          //            << xn[n][0] << " "
          //            << xn[n][1] << " "
          //            << "is not on curved boundary "
          //            << ibound << std::endl;
         }
       }// end record boundary nodes


    
#ifdef PARANOID

      std::stringstream error_stream;
      error_stream
       << "Problem when upgrading boundary element " << e
       << " " << bulk_el_pt << " on boundary " << ibound
       << "\nIts nodes are at\n"
       << xn[0][0] << " "  << xn[0][1] << "\n"
       << xn[1][0] << " "  << xn[1][1] << "\n"
       << xn[2][0] << " "  << xn[2][1] << "\n"
       << std::endl;
                   
      // Check nnode_on_curved_boundary
      if (nnode_not_on_curved_boundary == 0)
       {
        std::string error_message=
         "No interior nodes. One node per CurvedElement must be interior.\n"+
         error_stream.str();
         
        throw OomphLibError(
         error_message,
         OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
       }
      else if (nnode_not_on_curved_boundary>1)
       {        
        std::stringstream error_stream2;
        error_stream2
         << "Multiple interior nodes: "
         << nnode_not_on_curved_boundary
         << " Only one node per CurvedElement can be interior.\n"
         << error_stream.str();
        throw OomphLibError(
         error_stream2.str(),
         OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
       }
      
#endif
      
      
      // hierher shouldn't these be called zeta (everywhere; sigh)
      // boundary coordinate at the next (cyclic) node after interior
      const double s_ubar =
       c1_curve_pt->get_zeta(xn[(index_of_interior_node+1) % 3]);

      // boundary coordinate at the previous (cyclic) node before interior
      const double s_obar =
       c1_curve_pt->get_zeta(xn[(index_of_interior_node+2) % 3]);
      


      // hierher there's still a problem but we do need to fix it
      // when it happens. Current observation: there's a kink
      // along internal "curved" boundaries

      
      // hierher Aidan check unnecessary (we think!)
// #ifdef PARANOID
//       // Check for inverted elements
//       if (s_ubar>s_obar)
//        {
//         // throw // hierher probably not an issue actually
//         OomphLibError(
//          "Decreasing parametric coordinate. Parametric coordinate must increase as the edge is traversed anti-clockwise.",
//          OOMPH_CURRENT_FUNCTION,
//          OOMPH_EXCEPTION_LOCATION);
//        } // end checks
// #endif
      

      
      // Assign edge case
      edge = static_cast<C1PlateHelper::CurvedEdgeEnumeration>(index_of_interior_node);
  
      
      // Upgrade it
      TemplateFreeCurvableBellElement* curv_el_pt=
       dynamic_cast<TemplateFreeCurvableBellElement*>(bulk_el_pt);
#ifdef PARANOID
      if (curv_el_pt==0)
       {
        throw OomphLibError(
         "Cast to TemplateFreeCurvableBellElement failed",
         OOMPH_CURRENT_FUNCTION,
         OOMPH_EXCEPTION_LOCATION);
       }
#endif

      // By default (this is a generic helper function we use fifth order 
      // polynomials for the boundary representation. third order (the
      // other option) is cheaper but doesn't work for all types of
      // boundary conditions.
      curv_el_pt->upgrade_element_to_curved(edge, s_ubar, s_obar,
                                            c1_curve_pt,
                                            boundary_order);


      // Doc?
      if (Upgraded_to_curved_edge_element_stream.is_open())
       {
        unsigned nplot=5;
        curv_el_pt->output(Upgraded_to_curved_edge_element_stream,nplot);
       }
      
      
     }

   } // end of loop over outer boundaries
  
 } // end_upgrade_elements


 
//======================================================================
/// Function to set up rotated nodes on the boundary: necessary if we want to set
/// up physical boundary conditions on a curved boundary with Hermite type dofs.
/// For example if we know w(n,t) = f(t) (where n and t are the
/// normal and tangent to a boundary) we ALSO know dw/dt and d2w/dt2.
/// NB no rotation is needed if the edges are completely free!
//======================================================================
 void rotate_edge_coordinates(
  Mesh* bulk_mesh_pt, 
  std::map<unsigned,C1CurviLine*> c1_curviline_pt) 
{
 
 // Loop over the bulk elements: Yes, really because we also
 // need to deal with those that only
 // have a single node on the boundary!
 unsigned n_element = bulk_mesh_pt-> nelement();
 for(unsigned e=0; e<n_element; e++)
  {
   // Get pointer to bulk element 
   FiniteElement* el_pt = bulk_mesh_pt->finite_element_pt(e);

   // Does the element have nodes where coordinates have been rotated?
   // hierher Aidan: is this the right way to think about it?
   bool el_rotated=false;
   
   // Loop over the curvilinear parts of the outer boundary
   for (const auto& [b, c1_curve_pt] : c1_curviline_pt)
    {
     // local node numbers of nodes on external boundaries
     Vector<unsigned> boundary_node;
     
     // Boundary coordinates of nodes on the external boundaries
     Vector<double> boundary_coordinate_of_node;
     
     // Loop over vertex nodes (they come first)
     const unsigned nnode=3;
     for (unsigned n=0; n<nnode;++n)
      {
       Node* node_pt=el_pt->node_pt(n);

#ifdef PARANOID
       // Is the node on multiple curved boundaries? If so, die!
       std::set<unsigned>* boundaries_pt=0;
       node_pt->get_boundaries_pt(boundaries_pt);
       if (boundaries_pt!=0)
        {
         if (boundaries_pt->size()>1)
          {
           bool all_boundaries_are_polygonal=true;
           for (unsigned b : *boundaries_pt)
            {
             if (Boundary_is_curvilinear[b])
              {
               all_boundaries_are_polygonal=false;
               break;
              }
            }
           if (!all_boundaries_are_polygonal)
            {
             std::stringstream error_stream;
             error_stream
              << "Node at "
              << node_pt->x(0) << " "
              << node_pt->x(1) << " "
              << "is on " << boundaries_pt->size()
              << " curved boundaries.\n"
              << "This shouldn't happen since it'll be unclear which curved\n"
              << "boundary it's supposed to get its rotated (tangential and\n"
              << "normal) directions from."
              << std::endl;
             throw OomphLibError(
              error_stream.str(),
              OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
            }
          }
        }
#endif
       
       // If on curvilinear boundary b
       if ((node_pt->is_on_boundary(b))&&(Boundary_is_curvilinear[b]))
        {
         boundary_node.push_back(n);
         double coord = c1_curve_pt->get_zeta(el_pt->node_pt(n)->position());
         boundary_coordinate_of_node.push_back(coord);
         // Doc?
         if (Rotated_node_output_stream.is_open())
          {
           Rotated_node_output_stream
            << node_pt->x(0) << " "
            << node_pt->x(1) << " "
            << e << " " 
            << b << " "
            << std::endl;
          }
        }
      }
     
     // If the element has nodes on the curvilinear boundary, setup rotation machinery
     if(!boundary_node.empty())
      {       
       TemplateFreeCurvableBellElement* curv_el_pt=
        dynamic_cast<TemplateFreeCurvableBellElement*>(el_pt);
#ifdef PARANOID
       if (curv_el_pt==0)
        {
         throw OomphLibError(
          "Cast to TemplateFreeCurvableBellElement failed",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
        }
#endif

       // This reinterprets the derivative dofs as derivatives w.r.t.
       // normal and tangent diretions (defined by the c1_curve_pt object,
       // evaluated at the specified boundary coordinate on it), rather than
       // raw x and y derivatives.
       el_rotated=true;
       curv_el_pt->
        rotated_boundary_helper_pt()->
        set_nodal_boundary_parametrisation(boundary_node,
                                           boundary_coordinate_of_node,
                                           c1_curve_pt);
      }
    }

   // Doc?
   if (Rotated_element_output_stream.is_open())
    {
     if (el_rotated)
      {
        unsigned nplot=5;
        el_pt->output(Rotated_element_output_stream,nplot);
      }
    }

  }

 
  } // end rotate_edge_coordinates

 

} // end namespace C1PlateHelper

} // end namespace oomph

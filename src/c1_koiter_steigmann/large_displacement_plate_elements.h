// Header file for the Biharmonic Bell elements
#ifndef OOMPH_NONLINEAR_PLATE_ELEMENTS_HEADER
#define OOMPH_NONLINEAR_PLATE_ELEMENTS_HEADER


#include<sstream>

//OOMPH-LIB headers
#include "../generic/nodes.h"
#include "../generic/oomph_utilities.h"
#include "../generic/Telements.h"
#include "../generic/Subparametric_Telements.h"

namespace oomph
{
//=============================================================
/// A class for all subparametric elements that solve the 2D-
/// Biharmonic equations.
/// \f[
/// \frac{\partial^4 u}{\partial x_i^4} = f(x_j)
/// \f]
/// This contains the generic maths. Shape functions, geometric
/// mapping etc. must get implemented in derived class.
//=============================================================
template <unsigned DIM, unsigned NNODE_1D>
class LargeDisplacementPlateEquations : public virtual BellElement<DIM,NNODE_1D>
{

public:

 /// \short Function pointer to pressure function fct(x,f(x)) --
 /// x is a Vector!
// typedef void (*PressureFctPt)(const Vector<double>& x, double& f);
 typedef void (*StressFctPt)(const Vector<double>& x, const Vector<double>& u,
   const DenseMatrix<double>& strain,  const DenseMatrix<double>& g_tensor,
    DenseMatrix<double>& stress);

 typedef void (*DStressFctPt)(const Vector<double>& x, const Vector<double>& u,
   const DenseMatrix<double>& strain,  const DenseMatrix<double>& g_tensor,
    RankThreeTensor<double>& d_stress_dui, RankFourTensor<double>& d_stress_dstrain);

 /// \short Function pointer to gradient of pressure function  fct(x,g(x)) --
 /// x is a Vector!
 typedef void (*PressureVectorFctPt)(const Vector<double>& x, const
  Vector<double>& u,const DenseMatrix<double>& grad_u, const Vector<double>& n, 
  Vector<double>& pressure_vector);

 /// \short Function pointer to gradient of pressure function  fct(x,g(x)) --
 /// x is a Vector!
 typedef void (*DPressureVectorFctPt)(const Vector<double>& x, const
  Vector<double>& u,const DenseMatrix<double>& grad_u, const Vector<double>& n, 
  DenseMatrix<double>& d_pressure_vector);

 /// \short Function pointer to gradient of pressure function  fct(x,g(x)) --
 /// x is a Vector!
 typedef void (*DPressureVectorDMatrixFctPt)(const Vector<double>& x, const
  Vector<double>& u,const DenseMatrix<double>& grad_u, const Vector<double>& n, 
  RankThreeTensor<double>& d_pressure_vector);

 /// \short Function pointer to the Error Metric we are using
 ///  e.g could be that we are just interested in error on w etc.
 typedef void (*ErrorMetricFctPt)(const Vector<double>& x, const 
  Vector<double>& u, const Vector<double>& u_exact, double& error, double& norm);

 /// \shot Function pointer to the Error Metric we are using if we want multiple
 ///  errors.  e.g could be we want errors seperately on each displacment
 typedef void (*MultipleErrorMetricFctPt)(const Vector<double>& x, const 
  Vector<double>& u, const Vector<double>& u_exact, Vector<double>& error, 
  Vector<double>& norm);

 /// Constructor (must initialise the Pressure_fct_pt to null)
 LargeDisplacementPlateEquations() : Pressure_fct_pt(0),  D_pressure_dr_fct_pt(0), 
  D_pressure_dn_fct_pt(0), D_pressure_d_grad_u_fct_pt(0), Stress_fct_pt(0), 
  D_stress_fct_pt(0), Error_metric_fct_pt(0), Multiple_error_metric_fct_pt(0),
  Number_of_internal_dofs(0), Number_of_internal_dof_types(0),
  Association_matrix_pt(0), Use_finite_difference_jacobian(false), 
  Association_matrix_is_cached(false)
 {
  // Poisson ratio is 0.5 (incompressible) by default.
  Nu_pt = &Default_Nu_Value;
  Eta_u_pt = &Default_Eta_Value;
  Eta_sigma_pt = &Default_Eta_Value;

  // Thickness ratio small (0.001). Might expect thickness independence in this
  // limit depending on the magnitude of the external forcing
  Thickness_pt = &Default_Thickness_Value;
 } 

 /// Broken copy constructor
 LargeDisplacementPlateEquations(const LargeDisplacementPlateEquations& dummy)
  {
   BrokenCopy::broken_copy("LargeDisplacementPlateEquations");
  }

 /// Broken assignment operator
 void operator=(const LargeDisplacementPlateEquations&)
  {
   BrokenCopy::broken_assign("LargeDisplacementPlateEquations");
  }

 /// Enable Finite difference Jacobian
 void enable_finite_difference_jacobian() {Use_finite_difference_jacobian = true;}

 /// Disable Finite difference Jacobian
 void disable_finite_difference_jacobian() {Use_finite_difference_jacobian = false;}

 /// \short Return the index at which the ith displacement l-type (unknown) value
 /// is stored.
 /// In derived multi-physics elements, this function should be overloaded
 /// to reflect the chosen storage scheme. Note that these equations require
 /// that the unknown is always stored at the same index at each node.
 inline unsigned u_index_koiter_model(const unsigned& l, const unsigned& i)
  const 
 {return i*(this->nnodal_position_type()) + l;}

 /// Output with default number of plot points
 void output(std::ostream &outfile)
  {
   const unsigned n_plot=5;
   LargeDisplacementPlateEquations::output(outfile,n_plot);
  }

 /// \short Output FE representation of soln: x,y,u or x,y,z,u at
 /// n_plot^DIM plot points
 void output(std::ostream &outfile, const unsigned &n_plot);

 /// C_style output with default number of plot points
 void output(FILE* file_pt)
  {
   const unsigned n_plot=5;
   LargeDisplacementPlateEquations::output(file_pt,n_plot);
  }

 /// \short C-style output FE representation of soln: x,y,u or x,y,z,u at
 /// n_plot^DIM plot points
 void output(FILE* file_pt, const unsigned &n_plot);

 /// Output exact soln: x,y,u_exact or x,y,z,u_exact at n_plot^DIM plot points
 void output_fct(std::ostream &outfile, const unsigned &n_plot,
                 FiniteElement::SteadyExactSolutionFctPt exact_soln_pt);

 /// \short Output exact soln: x,y,u_exact or x,y,z,u_exact at
 /// n_plot^DIM plot points (dummy time-dependent version to
 /// keep intel compiler happy)
 void output_fct(std::ostream &outfile, const unsigned &n_plot,
                         const double& time,
                         FiniteElement::UnsteadyExactSolutionFctPt
                         exact_soln_pt)
  {
   throw OomphLibError(
    "There is no time-dependent output_fct() for these elements ",
    OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }


 /// Get error against and norm of exact solution
 void compute_error(std::ostream &outfile,
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
                    double& error, double& norm);

 /// Get error against and norm of exact solution
 void compute_error(std::ostream &outfile,
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
                    Vector<double>& error, Vector<double>& norm);

 /// Dummy, time dependent error checker
 inline void compute_error(std::ostream &outfile,
 FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt, const double& time, 
 double& error, double& norm)
  {
  throw OomphLibError(
   "There is no time-dependent compute_error() for these elements",
   OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }

 // Get pointer to association matrix 
 inline DenseMatrix<double> *get_association_matrix_pt()const 
  {return Association_matrix_pt;}

 /// Access function: Pointer to pressure function
 inline PressureVectorFctPt& pressure_fct_pt() {return Pressure_fct_pt;}

 /// Access function: Pointer to pressure function
 inline DPressureVectorFctPt& d_pressure_dn_fct_pt() {return D_pressure_dn_fct_pt;}

 /// Access function: Pointer to pressure function
 inline DPressureVectorFctPt& d_pressure_dr_fct_pt() {return D_pressure_dr_fct_pt;}

 /// Access function: Pointer to pressure function
 inline DPressureVectorDMatrixFctPt& d_pressure_d_grad_u_fct_pt() 
  {return D_pressure_d_grad_u_fct_pt;}

 /// Access function: Pointer to Stress function
 inline StressFctPt& stress_fct_pt() {return Stress_fct_pt;}

 /// Access function: Pointer to Stress function
 inline DStressFctPt& d_stress_fct_pt() {return D_stress_fct_pt;}

 /// Access function: Pointer to error metric function
 inline ErrorMetricFctPt& error_metric_fct_pt() {return Error_metric_fct_pt;}

 /// Access function: Pointer to multiple error metric function
 inline MultipleErrorMetricFctPt& multiple_error_metric_fct_pt() 
  {return Multiple_error_metric_fct_pt;}

 /// Access function: Pointer to pressure function. Const version
 inline PressureVectorFctPt pressure_fct_pt() const {return Pressure_fct_pt;}

 /// Access function: Pointer to pressure function
 inline DPressureVectorFctPt d_pressure_dn_fct_pt() const {return D_pressure_dn_fct_pt;}

 /// Access function: Pointer to pressure function
 inline DPressureVectorFctPt d_pressure_dr_fct_pt() const {return D_pressure_dr_fct_pt;}

 /// Access function: Pointer to pressure function
 inline DPressureVectorDMatrixFctPt d_pressure_d_grad_u_fct_pt() const 
  {return D_pressure_d_grad_u_fct_pt;}

 /// Access function: Pointer to stress function. Const version
 inline StressFctPt stress_fct_pt() const {return Stress_fct_pt;}

 /// Access function: Pointer to stress function. Const version
 inline DStressFctPt d_stress_fct_pt() const {return D_stress_fct_pt;}

 /// Access function: Pointer to error metric function function
 inline ErrorMetricFctPt error_metric_fct_pt() const {return Error_metric_fct_pt;}

 /// Access function: Pointer to multiple error metric function
 inline MultipleErrorMetricFctPt multiple_error_metric_fct_pt() const 
  {return Multiple_error_metric_fct_pt;}

 ///Access function to the Poisson ratio.
 inline const double*& nu_pt() {return Nu_pt;}

 ///Access function to the Poisson ratio.
 inline const double*& thickness_pt() {return Thickness_pt;}

 ///Access function to the Poisson ratio (const version)
 inline const double& get_nu() const {return *Nu_pt;}

 ///Access function to the Poisson ratio (const version)
 inline const double& get_thickness() const {return *Thickness_pt;}

 // Get the kth dof type at internal point l
 virtual double get_u_bubble_dof(const unsigned& l, const unsigned& k) const =0;

 // Get the kth equation at internal point l
 virtual int local_u_bubble_equation(const unsigned& l, const unsigned& k)
   const =0;

 // Precompute the association matrix, pure virtual
 virtual void precompute_association_matrix(DenseMatrix<double>& m)=0;
 // Get the number of basis functions, pure virtual
 virtual double n_basis_functions()=0;
 // Get the number of basic basis functions, pure virtual
 virtual double n_basic_basis_functions()=0;

 /// Get pressure term at (Eulerian) position x. This function is
 /// virtual to allow overloading in multi-physics problems where
 /// the strength of the pressure function might be determined by
 /// another system of equations.
 inline void get_pressure_biharmonic(const unsigned& ipt,
   const Vector<double>& x, const Vector<double>& u, const DenseMatrix<double>&
   grad_u, const Vector<double>& n, Vector<double>& pressure) const
  {
   //If no pressure function has been set, return zero
   if(Pressure_fct_pt==0)
    {
    pressure[0] = 0.0;
    pressure[1] = 0.0;
    pressure[2] = 0.0;
    }
   else
    {
    // Zero the pressure (as a precaution)
    pressure[0] = 0.0;
    pressure[1] = 0.0;
    pressure[2] = 0.0;
    // Get pressure strength
    (*Pressure_fct_pt)(x,u,grad_u,n,pressure);
    }
  }

 /// Get pressure term at (Eulerian) position x. This function is
 /// virtual to allow overloading in multi-physics problems where
 /// the strength of the pressure function might be determined by
 /// another system of equations.
 inline void get_d_pressure_biharmonic_dn(const unsigned& ipt,
   const Vector<double>& x, const Vector<double>& u, const DenseMatrix<double>&
   grad_u, const Vector<double>& n, DenseMatrix<double>& d_pressure_dn) const
  {
   //If no pressure function has been set, return zero
   if(D_pressure_dn_fct_pt==0)
    {
    for(unsigned i=0;i<Number_of_displacements;++i)
     {
     for(unsigned j=0;j<Number_of_displacements;++j)
      { d_pressure_dn(i,j) =0.0;}
     }
    }
   else
    {
    // Zero the pressure (as a precaution)
    for(unsigned i=0;i<Number_of_displacements;++i)
     {
     for(unsigned j=0;j<Number_of_displacements;++j)
      { d_pressure_dn(i,j) =0.0;}
     }
    // Get d_p_dn
    (*D_pressure_dn_fct_pt)(x,u,grad_u,n,d_pressure_dn);
    }
  }

 /// Get pressure term at (Eulerian) position x. This function is
 /// virtual to allow overloading in multi-physics problems where
 /// the strength of the pressure function might be determined by
 /// another system of equations.
 inline void get_d_pressure_biharmonic_dr(const unsigned& ipt,
   const Vector<double>& x, const Vector<double>& u, const DenseMatrix<double>&
   grad_u,const Vector<double>& n,DenseMatrix<double>& d_pressure_dr) const
  {
  //If no pressure function has been set, return zero
  if(D_pressure_dr_fct_pt==0)
   {
   for(unsigned i=0;i<Number_of_displacements;++i)
    {
    for(unsigned j=0;j<Number_of_displacements;++j)
     { d_pressure_dr(i,j) =0.0;}
    }
   }
  else
   {
   // Zero the pressure (as a precaution)
   for(unsigned i=0;i<Number_of_displacements;++i)
    {
    for(unsigned j=0;j<Number_of_displacements;++j)
     { d_pressure_dr(i,j) =0.0;}
    }
   // Get d_p_dn
   (*D_pressure_dr_fct_pt)(x,u,grad_u,n,d_pressure_dr);
   }
  }

 /// Get pressure term at (Eulerian) position x. This function is
 /// virtual to allow overloading in multi-physics problems where
 /// the strength of the pressure function might be determined by
 /// another system of equations.
 inline void get_d_pressure_biharmonic_d_grad_u(const unsigned& ipt,
   const Vector<double>& x, const Vector<double>& u, const DenseMatrix<double>&
   grad_u, const Vector<double>& n, RankThreeTensor<double>& d_pressure_d_grad_u) const
  {
   //If no pressure function has been set, return zero
   if(D_pressure_dn_fct_pt==0)
    {
    for(unsigned i=0;i<Number_of_displacements;++i)
     {
     for(unsigned j=0;j<Number_of_displacements;++j)
      { 
      for(unsigned alpha=0;alpha<2;++alpha)
       {d_pressure_d_grad_u(i,j,alpha) =0.0;}
      }
     }
    }
   else
    {
    // Zero the pressure (as a precaution)
    for(unsigned i=0;i<Number_of_displacements;++i)
     {
     for(unsigned j=0;j<Number_of_displacements;++j)
      { 
      for(unsigned alpha=0;alpha<2;++alpha)
       {d_pressure_d_grad_u(i,j,alpha) =0.0;}
      }
     }
    // Get d_p_dn
    (*D_pressure_d_grad_u_fct_pt)(x,u,grad_u,n,d_pressure_d_grad_u);
    }
  }

 /// Add the element's contribution to its residual vector (wrapper)
 void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   DenseMatrix<double> conversion_matrix (n_basis_functions(),
    n_basic_basis_functions(),0.0);
  // Precompute if not cached
  if(! Association_matrix_is_cached)
   {
   // Precompute the association matrix
   this->precompute_association_matrix(conversion_matrix);
   this->Association_matrix_pt=&conversion_matrix;
   }

  //Call the generic residuals function with flag set to 0
  //using a dummy matrix argument
  fill_in_generic_residual_contribution_biharmonic(
   residuals,GeneralisedElement::Dummy_matrix,0);

  // Reset to zero
  if(! Association_matrix_is_cached)
   {
    this->Association_matrix_pt=0;
   }
  }


 /// Add the element's contribution to its residual vector and
 /// element Jacobian matrix (wrapper)
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                   DenseMatrix<double> &jacobian)
  {
   // Precompute the association matrix
   DenseMatrix<double> conversion_matrix (n_basis_functions(),
    n_basic_basis_functions(),0.0);
   this->precompute_association_matrix(conversion_matrix);
   this->Association_matrix_pt=&conversion_matrix;
   Association_matrix_is_cached = true;
 
   //Call the generic routine with the flag set to 1
   if(! Use_finite_difference_jacobian)
    {
     fill_in_generic_residual_contribution_biharmonic(residuals,jacobian,1);
    }
   else
    {
     // Otherwise call the default
     FiniteElement::fill_in_contribution_to_jacobian(residuals,jacobian);
    }
 
   // Reset to zero
   this->Association_matrix_pt=0;
   Association_matrix_is_cached = false;
  }

  /// Add the element's contribution to its residual vector and
  /// element Jacobian matrix (wrapper)
  void fill_in_contribution_to_jacobian_and_mass_matrix(Vector<double> &residuals,
DenseMatrix<double> &jacobian,DenseMatrix<double> &mass_matrix)
   {
   //Call fill in Jacobian 
   fill_in_contribution_to_jacobian(residuals,jacobian);
   // There is no mass matrix: we will just want J w = 0
 
   // -- COPIED FROM DISPLACMENT FVK EQUATIONS --
   // Dummy diagonal (won't result in global unit matrix but
   // doesn't matter for zero eigenvalue/eigenvector
   unsigned ndof=mass_matrix.nrow();
   for (unsigned i=0;i<ndof;i++)
    {
     mass_matrix(i,i)+=1.0;
    }
   }

 // Return the interpolated unit normal
 virtual void fill_in_metric_tensor(const DenseMatrix<double>& interpolated_drdxi, 
DenseMatrix<double>& g_tensor) = 0;

 // Get the tangent vectors
 virtual void fill_in_tangent_vectors(const DenseMatrix<double>& interpolated_dudxi, 
  DenseMatrix<double>& interpolated_drdxi) = 0 ;

 // Return the determinant of the metric tensor
 virtual double metric_tensor_determinant(const DenseMatrix<double>& 
  interpolated_drdxi) = 0;

 // Return the determinant of the metric tensor
 inline double two_by_two_determinant(const DenseMatrix<double>& g_tensor)
  {
  // Now return determinant
  return g_tensor(0,0)*g_tensor(1,1) - g_tensor(0,1)*g_tensor(1,0);
  }

 // Return the determinant of the metric tensor
 virtual void d_metric_tensor_determinant_du_unknown(const DenseMatrix<double>& 
  g_tensor, const RankFourTensor<double>& d_g_tensor_du_unknown,
  DenseMatrix<double>& d_metric_tensor_determinant_du_unknown) = 0;

 // Return the interpolated unit normal
 virtual void fill_in_unit_normal(const DenseMatrix<double>&
  interpolated_drdxi, Vector<double>& normal) = 0 ;

 // Return the interpolated unit normal
 virtual void fill_in_d_unit_normal_du_unknown(const DenseMatrix<double>&
  interpolated_drdxi, const DenseMatrix<double>& g_tensor,  const RankFourTensor<double>&
  d_g_tensor_du_unknown, RankThreeTensor<double>& d_normal_du_unknown) = 0;

 // Get the (Green Lagrange) strain tensor
 virtual void fill_in_strain_tensor(const DenseMatrix<double>&
  interpolated_dudxi, DenseMatrix<double>& strain) =0 ;

 // Get the (Green Lagrange) strain tensor
 virtual void fill_in_d_g_tensor_du_unknown(const DenseMatrix<double>&
 interpolated_dudxi, RankFourTensor<double>& d_g_tensor_du_unknown) = 0;

 // Get the (Green Lagrange) strain tensor
 virtual void fill_in_d_strain_tensor_du_unknown(const DenseMatrix<double>&
 interpolated_dudxi, RankFourTensor<double>& d_epsilon_tensor_du_unknown) = 0;

 // Fill in the Kirchhoff St Venant stress tensor
 inline void fill_in_kirchhoff_st_venant_stress(const DenseMatrix<double>& strain, 
  DenseMatrix<double>& stress)
  {
  // Zero the stress
  for(unsigned alpha=0; alpha<2; ++alpha)
   {
   for(unsigned beta=0; beta<2; ++beta)
     { stress(alpha,beta) =0.0; }
   }

   // Get Poisson ratio
   double nu=get_nu();
  // Loop over alpha and beta
  for(unsigned alpha=0; alpha<2; ++alpha)
   {
   for(unsigned beta=0; beta<2; ++beta)
    {
     stress(alpha,beta) += (1-nu) * strain(alpha,beta)/(1-nu*nu);
     stress(alpha,alpha) += nu * strain(beta,beta)/(1-nu*nu);
    }
   } 
  }

 // Fill in the Kirchhoff St Venant stress tensor
 inline void fill_in_d_kirchhoff_st_venant_stress_du_unknown(const 
  RankFourTensor<double>& d_strain_du_unknown, RankFourTensor<double>& 
  d_stress_du_unknown)
  {
  // Zero the vector
  for(unsigned alpha=0; alpha<2; ++alpha)
   {
   for(unsigned beta=0; beta<2; ++beta)
    {
    for(unsigned i=0;i<Number_of_displacements;++i)
     {
     d_stress_du_unknown(alpha,beta,i,0) = 0.0;
     d_stress_du_unknown(alpha,beta,i,1) = 0.0;
     d_stress_du_unknown(alpha,beta,i,2) = 0.0;
     }
    }
   }

  // Get Poisson ratio
  double nu=get_nu();
  // Loop over alpha and beta
  for(unsigned alpha=0; alpha<2; ++alpha)
   {
   for(unsigned beta=0; beta<2; ++beta)
    {
    for(unsigned i=0; i<Number_of_displacements; ++i)
     {
     for(unsigned gamma=0; gamma<2; ++gamma)
      {
      d_stress_du_unknown(alpha,beta,i,1+gamma) += (1-nu) * d_strain_du_unknown(alpha,
        beta,i,gamma) / (1-nu*nu);
      d_stress_du_unknown(alpha,alpha,i,1+gamma) += nu * d_strain_du_unknown(beta,
        beta,i,gamma) / (1-nu*nu);
      }
     }
    }
   } 
  }

 // Fill in the curvature tensor
 virtual void fill_in_curvature_tensor(const Vector<double>& unit_normal,
  const DenseMatrix<double>& interpolated_d2rdxi2, DenseMatrix<double>& 
  curvature) = 0;

 // Fill in the curvature tensor
 virtual void fill_in_d_curvature_tensor_du_unknown(const Vector<double>& 
  unit_normal, const DenseMatrix<double>& interpolated_d2rdxi2, const 
  RankThreeTensor<double>&  d_unit_normal_dunknown,  
  RankFourTensor<double>& d_curvature_du_unknown ) = 0;

 // Fill in the moment tensor
 virtual void fill_in_moment_tensor(const DenseMatrix<double>&
  interpolated_d2rxi2,  const Vector<double>& unit_normal, const 
  DenseMatrix<double>& curvature, RankThreeTensor<double>& moment) = 0;

 // Fill in the moment tensor
 virtual void fill_in_d_moment_tensor_du_unknown(const DenseMatrix<double>&
  interpolated_d2rxi2,const Vector<double>& unit_normal,
  const DenseMatrix<double>& curvature, const RankThreeTensor<double>&
  d_unit_normal_du_unknown,  const RankFourTensor<double>& d_curvature_du_unknown,
  RankFiveTensor<double>&  d_moment_du_unknown) = 0;

 // Get the (Green Lagrange) strain tensor
 // \Gamma_{\alpha\beta\gamma} = ( E_{\alpha\beta,gamma} + E_{\gamma\alpha,\beta}
 //   - E_{\beta\gamma,\alpha})
 inline void fill_in_stress_tensor(const Vector<double>& x, const Vector<double>& r, 
  const DenseMatrix<double>& strain, const DenseMatrix<double>& g_tensor, 
  DenseMatrix<double>& stress)
  {
  // Zero the stress
  for(unsigned alpha=0; alpha<2; ++alpha)
   {
   for(unsigned beta=0; beta<2; ++beta)
     { stress(alpha,beta) =0.0; }
   }

  // IF not set use Kirchhoff st venant
  if(Stress_fct_pt==0)
   {
    fill_in_kirchhoff_st_venant_stress(strain, stress);
   }
  // Use the user defined one
  else 
   {
    (*Stress_fct_pt)(x,r,strain,g_tensor,stress);
   }
  }

 // Get the (Green Lagrange) strain tensor
 // \Gamma_{\alpha\beta\gamma} = ( E_{\alpha\beta,gamma} + E_{\gamma\alpha,\beta}
 //   - E_{\beta\gamma,\alpha})
 inline void fill_in_d_stress_tensor_du_unknown(const Vector<double>& x, 
  const Vector<double>& r, const DenseMatrix<double>& strain, 
  const DenseMatrix<double>& g_tensor,  DenseMatrix<double>& stress,
   const RankFourTensor<double>&  d_strain_dunknown, RankFourTensor<double>& 
  d_stress_dunknown)
  {
  // Zero the vector
  for(unsigned alpha=0; alpha<2; ++alpha)
   {
   for(unsigned beta=0; beta<2; ++beta)
    {
    for(unsigned i=0;i<Number_of_displacements;++i)
     {
     // Flatpacked array ui ui,alpha
     for(unsigned k=0;k<3;++k)
      {d_stress_dunknown(alpha,beta,i,k) = 0.0;}
     }
    }
   }

  // IF not set use Kirchhoff st venant
  if(Stress_fct_pt==0 && D_stress_fct_pt==0)
   {
    fill_in_d_kirchhoff_st_venant_stress_du_unknown(d_strain_dunknown, d_stress_dunknown);
   }
  // Use the user defined one
  else if(D_stress_fct_pt != 0)
   {
   // Initialise the tensor (SLOW)
   RankFourTensor<double> d_stress_d_epsilon(2,2,2,2,0.0);
   RankThreeTensor<double> d_stress_du(2,2,Number_of_displacements,0.0);
    // Get the user defined tensor
   (*D_stress_fct_pt)(x,r,strain,g_tensor,d_stress_du,d_stress_d_epsilon);
   for(unsigned i=0; i<Number_of_displacements;++i) 
    {
    for(unsigned alpha=0; alpha<2;++alpha) 
     {
     for(unsigned beta=0; beta<2;++beta) 
      {
      d_stress_dunknown(alpha,beta,i,0) += d_stress_du(alpha,beta,i);
      for(unsigned gamma=0; gamma<2;++gamma) 
       {
       for(unsigned delta=0; delta<2;++delta) 
        {
        for(unsigned mu=0; mu<2;++mu) 
         {
         //Fill in using the rank four tensor (SLOW)
         d_stress_dunknown(alpha,beta,i,1+mu) += d_stress_d_epsilon(alpha,beta,
           gamma,delta)*d_strain_dunknown(gamma,delta,i,mu);
         }
        }
       }
      }
     }
    }
   }
  else 
   {
    throw OomphLibError(
   "You have defined a function for stress but not d_stress_depsilon. A finite\
difference stress is not yet defined so an error has been thrown.",
   OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
   // USE THE USER DEFINED ONE FILL IN
   //     (*D_stress_fct_pt)(strain,stress);
   ///HERE FINITE DIFFERENCE + WARNING
   }
  }

 virtual void fill_in_total_tension(const DenseMatrix<double>& stress, const
RankThreeTensor<double>& christoffel_tensor, const RankThreeTensor<double>&
 moment_tensors, const DenseMatrix<double>& interpolated_dudxi, 
 const DenseMatrix<double>& interpolated_d2rdxi2, DenseMatrix<double>&
 tension_vectors) = 0;

 virtual void d_fill_in_total_tension_du_unknown(const DenseMatrix<double>& stress
   , const RankThreeTensor<double>& christoffel_tensor, const 
   RankThreeTensor<double>& moment_tensors, const DenseMatrix<double>& interpolated_dudxi, 
    const DenseMatrix<double>& interpolated_d2udxi2,
   RankFourTensor<double>& d_stress_du_unknown, const RankFiveTensor<double>& 
   d_christoffel_tensor_du_unknown, const RankFiveTensor<double>& 
   d_moment_tensors_du_unknown,   RankFourTensor<double>&
   d_tension_vectors_du_unknown) =0;

 // Get the Christtoffel tensor of the second kind
 // \Gamma_{\alpha\beta\gamma} = ( E_{\alpha\beta,gamma} + E_{\gamma\alpha,\beta}
 //   - E_{\beta\gamma,\alpha})
 virtual void fill_in_second_christoffel_tensor(const DenseMatrix<double>&
  interpolated_dudxi, const DenseMatrix<double>& interpolated_d2udxi2, 
  RankThreeTensor<double>& gamma_tensor) = 0 ;
 
 // Get the Christtoffel tensor of the second kind
 // \Gamma_{\alpha\beta\gamma} = ( E_{\alpha\beta,gamma} + E_{\gamma\alpha,\beta}
 //   - E_{\beta\gamma,\alpha})
 virtual void fill_in_d_second_christoffel_tensor_dui_unknown(const DenseMatrix<double>&
  interpolated_dudxi, const DenseMatrix<double>& interpolated_d2udxi2, 
  RankFiveTensor<double>& d_gamma_tensor_du_unknown) = 0;

 /// \short Return FE representation of unknown values u(s)
 /// at local coordinate s
 inline void interpolated_u_koiter_plate(const Vector<double> &s,
   Vector<Vector<double> >& interpolated_u) const
  {
   //Find number of position dofs
   const unsigned n_position_type = this->nnodal_position_type();
   // Find the internal dofs
   const unsigned n_b_position_type = this->Number_of_internal_dof_types;
   //Find out how many nodes there are
   const unsigned n_node = this->nnode();
   //Find out how many internal points there are
   const unsigned n_b_node = this->Number_of_internal_dofs;

   //Local c1-shape funtion
   Shape psi(n_node,n_position_type),test(n_node,n_position_type),
    psi_b(n_b_node,n_b_position_type),test_b(n_b_node,n_b_position_type);

   DShape dpsi_dxi(n_node,n_position_type,DIM),dtest_dxi(n_node,n_position_type,DIM),
    dpsi_b_dxi(n_b_node,n_b_position_type,DIM),dtest_b_dxi(n_b_node,n_b_position_type,DIM),
    d2psi_dxi2(n_node,n_position_type,3), d2test_dxi2(n_node,n_position_type,3),
    d2psi_b_dxi2(n_b_node,n_b_position_type,3), d2test_b_dxi2(n_b_node,n_b_position_type,3);

   //Find values of c1-shape function
   d2shape_and_d2test_eulerian_biharmonic(s,psi,psi_b,dpsi_dxi,dpsi_b_dxi,
    d2psi_dxi2,d2psi_b_dxi2,test,test_b,dtest_dxi,dtest_b_dxi,d2test_dxi2,
    d2test_b_dxi2);

   //Interpolated unknown
   //Loop over displacements
   for(unsigned i=0;i<Number_of_displacements;++i)
   { 
    // Loop over nodes  
    for(unsigned l=0;l<n_node;l++)
     {
     // Loop over hermite dofs
     for(unsigned k=0;k<n_position_type;k++)
      {
      // Get the kth nodal value at node l for displacement i
      double u_value_ikl = this->raw_nodal_value(l,u_index_koiter_model(k,i));
      interpolated_u[i][0] += u_value_ikl*psi(l,k);
      interpolated_u[i][1] += u_value_ikl*dpsi_dxi(l,k,0);
      interpolated_u[i][2] += u_value_ikl*dpsi_dxi(l,k,1);
      interpolated_u[i][3] += u_value_ikl*d2psi_dxi2(l,k,0);
      interpolated_u[i][4] += u_value_ikl*d2psi_dxi2(l,k,1);
      interpolated_u[i][5] += u_value_ikl*d2psi_dxi2(l,k,2);
     }
    }
    // Loop over Bubble dofs
    for(unsigned l=0;l<Number_of_internal_dofs;l++)
     {
     // Loop over hermite dofs (only value dofs in internal nodes)
     for(unsigned k=0;k<Number_of_internal_dof_types;k++)
     {
      double u_value = get_u_bubble_dof(l,i);
      interpolated_u[i][0] += u_value * psi_b(l,k);
      interpolated_u[i][1] += u_value*dpsi_b_dxi(l,k,0);
      interpolated_u[i][2] += u_value*dpsi_b_dxi(l,k,1);
      interpolated_u[i][3] += u_value*d2psi_b_dxi2(l,k,0);
      interpolated_u[i][4] += u_value*d2psi_b_dxi2(l,k,1);
      interpolated_u[i][5] += u_value*d2psi_b_dxi2(l,k,2);
     }
    }
   }
  }

 /// \short Self-test: Return 0 for OK
 unsigned self_test();

 /// \short get the coordinate
 virtual void get_coordinate_x(const Vector<double>& s, Vector<double>& x) const=0;

 void pin_all_in_plane_dofs()
  {
   double n_displacements = Number_of_displacements;
   unsigned n_position_type = this->nnodal_position_type();
   // Pin the nodal dofs
   for(unsigned inode=0;inode<this->nnode();++inode)
    {
    Node* nod_pt = this->node_pt(inode);
    for(unsigned k=0;k<n_position_type;++k)
     {
     for(unsigned j=0;j<n_displacements-1;++j)
     {
      // Pin kth value and set to zero
      nod_pt->pin(u_index_koiter_model(k,j));
      nod_pt->set_value(u_index_koiter_model(k,j),0.0);
     }
    }
   }
  // Pin the internal dofs
  const double n_internal_dofs = this->Number_of_internal_dofs;
  for(unsigned i=0;i<n_internal_dofs;++i)
   {
   for(unsigned j=0;j<n_displacements-1;++j)
    {
    Data* internal_data_pt = this->internal_data_pt(1);
    // Pin kth value and set to zero
    internal_data_pt->pin(i+j*n_internal_dofs);
    // HERE need a function that can give us this lookup
    internal_data_pt->set_value(i+j*n_internal_dofs,0.0);
    }
   }  
  }

 void pin_out_of_plane_dofs()
  {
   double n_displacements = Number_of_displacements;
   unsigned n_position_type = this->nnodal_position_type();
   // Pin the nodal dofs
   for(unsigned inode=0;inode<this->nnode();++inode)
    {
     Node* nod_pt = this->node_pt(inode);
     for(unsigned k=0;k<n_position_type;++k)
      {
      for(unsigned j=n_displacements-1;j<n_displacements;++j)
      {
       // Pin kth value and set to zero
       nod_pt->pin(u_index_koiter_model(k,j));
       nod_pt->set_value(u_index_koiter_model(k,j),0.0);
      }
     }
    }
  // Pin the internal dofs
  const double n_internal_dofs = this->Number_of_internal_dofs;
  for(unsigned i=0;i<n_internal_dofs;++i)
   {
   for(unsigned j=n_displacements-1;j<n_displacements;++j)
    {
    Data* internal_data_pt = this->internal_data_pt(1);
    // Pin kth value and set to zero
    internal_data_pt->pin(i+j*n_internal_dofs);
    // HERE need a function that can give us this lookup
    internal_data_pt->set_value(i+j*n_internal_dofs,0.0);
    }
   }  
   }

/// \short HERE eta_u is untested feature so we have disabled it until we
/// have validated it 
//  ///Access function to the z displacement scaling in the displacement.
//  virtual const double*& eta_u_pt() {return Eta_u_pt;}

 ///Access function to the in plane displacment scaling in the displacement.
 virtual const double*& eta_sigma_pt() {return Eta_sigma_pt;}

 ///Access function to the displacement scaling (HERE always 1. - until valid.).
 virtual const double& eta_u()const {return *Eta_u_pt;}

 ///Access function to the in plane displacment scaling.
 virtual const double& eta_sigma()const {return *Eta_sigma_pt;}

protected:

 /// \short Shape/test functions and derivs w.r.t. to global coords at
 /// local coord. s; return  Jacobian of mapping
 virtual double d2shape_and_d2test_eulerian_biharmonic(const Vector<double> &s,
  Shape &psi, Shape& psi_b, DShape &dpsi_dx, DShape &dpsi_b_dx,
  DShape &d2psi_dx2,DShape& d2psi_b_dx2,
  Shape &test, Shape& test_b, DShape &dtest_dx, DShape &dtest_b_dx,
  DShape &d2test_dx2,DShape& d2test_b_dx2) const=0;

 /// \short Shape/test functions and derivs w.r.t. to global coords at
 /// local coord. s; return  Jacobian of mapping
 virtual double dshape_and_dtest_eulerian_biharmonic(const Vector<double> &s,
  Shape &psi, Shape& psi_b, DShape &dpsi_dx, DShape &dpsi_b_dx,
  Shape &test, Shape& test_b, DShape &dtest_dx, DShape &dtest_b_dx) const=0;

 /// \short Shape/test functions at local coordinate s
 virtual void shape_and_test_biharmonic(const Vector<double> &s,
  Shape &psi, Shape& psi_b, Shape &test, Shape& test_b) const=0;

 /// \short Compute element residual Vector only (if flag=and/or element
 /// Jacobian matrix
 virtual void fill_in_generic_residual_contribution_biharmonic(
  Vector<double> &residuals, DenseMatrix<double> &jacobian,
  const unsigned& flag);

 /// Pointer to pressure function:
 PressureVectorFctPt Pressure_fct_pt;

 /// Pointer to pressure function:
 DPressureVectorFctPt D_pressure_dr_fct_pt;

 /// Pointer to pressure function:
 DPressureVectorFctPt D_pressure_dn_fct_pt;

 /// Pointer to pressure function:
 DPressureVectorDMatrixFctPt D_pressure_d_grad_u_fct_pt;

 /// Pointer to stress function
 StressFctPt Stress_fct_pt;

 /// Pointer to epsilon derivative of stress function
 DStressFctPt D_stress_fct_pt;

 /// Pointer to error metric
 ErrorMetricFctPt Error_metric_fct_pt;

 /// Pointer to error metric when we want multiple errors
 MultipleErrorMetricFctPt Multiple_error_metric_fct_pt;

 /// Pointer to Poisson ratio, which this element cannot modify
 const double* Nu_pt;

 /// Pointer to Non dimensional thickness parameter
 const double* Thickness_pt;

 /// \short unsigned that holds the internal 'bubble' dofs the element has -
 // zero for Bell Elements and 3 for C1 curved elements
 unsigned Number_of_internal_dofs;
 
 /// \short unsigned that holds the number if displacments the element has 
 static const unsigned Number_of_displacements;

 /// \short unsigned that holds the number of displacements the element has 
 static const double Default_Nu_Value;

 static const double Default_Eta_Value; //HERE MOVE

 static const double Default_Thickness_Value;

 /// Pointer to displacement scaling
 const double* Eta_u_pt;

 /// Pointer to stress scaling
 const double* Eta_sigma_pt;

 /// \short unsigned that holds the number of types of degree of freedom at each
 // internal point that the element has zero for Bell Elements and 1
 //// for C1 curved elements
 unsigned Number_of_internal_dof_types;

 protected:

 /// Pointer to precomputed matrix that associates shape functions to monomials
 DenseMatrix<double> *Association_matrix_pt;

 /// Flag to use finite difference jacobian
 bool Use_finite_difference_jacobian;

 /// Bool to flag up when assocation matrix is cached
 bool Association_matrix_is_cached;

};

} //end namespace oomph
#endif


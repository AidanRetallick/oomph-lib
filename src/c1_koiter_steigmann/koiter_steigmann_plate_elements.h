// Header file for the Biharmonic Bell elements
#ifndef OOMPH_KOITER_STEIGMANN_PLATE_ELEMENTS_HEADER
#define OOMPH_KOITER_STEIGMANN_PLATE_ELEMENTS_HEADER


#include<sstream>

//OOMPH-LIB headers
#include "large_displacement_plate_elements.h"

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
class KoiterSteigmannPlateEquations : public virtual LargeDisplacementPlateEquations<DIM,NNODE_1D>
{

public:
 /// Constructor (must initialise the Pressure_fct_pt to null)
 KoiterSteigmannPlateEquations() : LargeDisplacementPlateEquations<DIM,NNODE_1D>()  {
  // Default parameters should lead to membrane like behaviour by default:
  // order 1 displacements for small thickness sheet
  // By default the displacements are asumed to be of order 1
 } 

 /// Broken copy constructor
 KoiterSteigmannPlateEquations(const KoiterSteigmannPlateEquations& dummy)
  {
   BrokenCopy::broken_copy("KoiterSteigmannPlateEquations");
  }

 /// Broken assignment operator
 void operator=(const KoiterSteigmannPlateEquations&)
  {
   BrokenCopy::broken_assign("KoiterSteigmannPlateEquations");
  }

// /// \short Output exact soln: x,y,u_exact or x,y,z,u_exact at
// /// n_plot^DIM plot points (dummy time-dependent version to
// /// keep intel compiler happy)
// virtual void output_fct(std::ostream &outfile, const unsigned &n_plot,
//                         const double& time,
//                         FiniteElement::UnsteadyExactSolutionFctPt
//                         exact_soln_pt)
//  {
//   throw OomphLibError(
//    "There is no time-dependent output_fct() for these elements ",
//    OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
//  }

 // Return the interpolated unit normal
 inline void fill_in_metric_tensor(const DenseMatrix<double>& interpolated_drdxi, 
DenseMatrix<double>& g_tensor)
  {
  // Zero the metric tensor 
  for(unsigned alpha=0;alpha<2;++alpha)
   {
    for(unsigned beta=0;beta<2;++beta)
     { g_tensor(alpha,beta)=0.0;  }
   }
  // Fill in metric tensor
  // Loop over plane tangent vectors
  for(unsigned alpha=0; alpha<2; ++alpha)
   {
   // Loop over plane tangent vectors
   for(unsigned beta=0; beta<2; ++beta)
    {
    // LSum over components of tangent vectors
    for(unsigned i=0; i<this->Number_of_displacements; ++i)
     {g_tensor(alpha,beta)+=interpolated_drdxi(i,alpha)
        *interpolated_drdxi(i,beta);}
    }
   }
  }

 // Get the tangent vectors
 inline void fill_in_tangent_vectors(const DenseMatrix<double>& interpolated_dudxi, 
  DenseMatrix<double>& interpolated_drdxi)
 {
 // Zero the tangent vectors 
 for(unsigned i=0;i<this->Number_of_displacements;++i)
  {
   for(unsigned beta=0;beta<2;++beta)
    { interpolated_drdxi(i,beta)=0.0;  }
  }
 // Fill in the tangent vectors
 for(unsigned i=0;i<this->Number_of_displacements;++i)
  {
  for(unsigned beta=0;beta<2;++beta)
   {
   // Scale the displacements
   const double eta_u =this->eta_u();
   // The displacement part
   interpolated_drdxi(i,beta) += eta_u*interpolated_dudxi(i,beta);
   // The coordinate part
   interpolated_drdxi(i,beta) += (i==beta ? 1 : 0);
   }
  }
 } 


 // Return the determinant of the metric tensor
 inline double metric_tensor_determinant(const DenseMatrix<double>& 
  interpolated_drdxi)
  {
  // Intitialise
  DenseMatrix<double> g_tensor(2,2,0.0);
  // Fill in metric tensor
  fill_in_metric_tensor(interpolated_drdxi,g_tensor);
  // Now return determinant
  return g_tensor(0,0)*g_tensor(1,1) - g_tensor(0,1)*g_tensor(1,0);
  }

 // Return the determinant of the metric tensor
 inline void d_metric_tensor_determinant_du_unknown(const DenseMatrix<double>& 
  g_tensor, const RankFourTensor<double>& d_g_tensor_du_unknown,
  DenseMatrix<double>& d_metric_tensor_determinant_du_unknown)
  {
  // // Adjugate metric tensor
  // DenseMatrix<double> adj_g(2,2);
  // adj_g(0,0)= g_tensor(1,1); 
  // adj_g(1,1)= g_tensor(0,0); 
  // adj_g(0,1)=-g_tensor(0,1); 
  // adj_g(1,0)=-g_tensor(1,0);

  // Loop over derivatives
  for(unsigned i=0;i<this->Number_of_displacements;++i)
   {
   d_metric_tensor_determinant_du_unknown(i,0)=0.0;
   d_metric_tensor_determinant_du_unknown(i,1)=0.0;
   for(unsigned gamma=0;gamma<2;++gamma)
    {
    for(unsigned delta=0;delta<2;++delta)
     {
     // Derivative: D(det A)(t) = adj(A) : D(A)(t)
     // d_metric_tensor_determinant_du_unknown[i]+=adj_g(gamma,delta) 
     //   * d_g_tensor_du_unknown(gamma,delta,i) ;
     // But we can write the adjugate like this (as it is symmetric)
     for(unsigned mu=0;mu<2;++mu)
      {
      d_metric_tensor_determinant_du_unknown(i,mu)+=g_tensor((gamma+1)%2,(delta+1)%2) 
         * d_g_tensor_du_unknown(gamma,delta,i,mu)*(gamma!=delta ? -1 : 1);
      } 
     }
    }
   }
  }

 // Return the interpolated unit normal
 inline void fill_in_unit_normal(const DenseMatrix<double>&
  interpolated_drdxi, Vector<double>& normal)
  {
  // Determinant of metric tensor
  const double g =  metric_tensor_determinant(interpolated_drdxi);
  const double sqrt_g = sqrt(g);

  // Zero the normal 
  for(unsigned i=0;i<this->Number_of_displacements;++i)
   { normal[i]=0.0; }

  // Loop over plane tangent vectors
  for(unsigned alpha=0;alpha<2;++alpha)
   {
   // Beta != Alpha
   const unsigned beta = (alpha + 1) % 2;
   Vector<double> g_alpha(3),g_beta(3);

   // Loop over displacement components
   for(unsigned i=0;i<this->Number_of_displacements;++i)
    {
    // Now fill in the tangent vectors
    g_alpha[i] = interpolated_drdxi(i,alpha);
    g_beta[i]  = interpolated_drdxi(i,beta);  
    }
    // Compute the cross product
    Vector<double> tmp(this->Number_of_displacements,0.0);
    VectorHelpers::cross(g_alpha,g_beta,tmp);
    // And fill in the normal 
   for(unsigned i=0;i<this->Number_of_displacements;++i)
    { normal[i] +=  tmp[i]/(2*sqrt_g) * (alpha == 0 ? 1 : -1) ;} 
   } 
  }

 // Return the interpolated unit normal
inline void fill_in_d_unit_normal_du_unknown(const DenseMatrix<double>&
  interpolated_drdxi, const DenseMatrix<double>& g_tensor,  const RankFourTensor<double>&
  d_g_tensor_du_unknown, RankThreeTensor<double>& d_normal_du_unknown)
  {
  // Loop over plane tangent vectors
  for(unsigned i=0;i<this->Number_of_displacements;++i)
   {
    for(unsigned j=0;j<this->Number_of_displacements;++j)
    {
     d_normal_du_unknown(i,j,0) =0.0; 
     d_normal_du_unknown(i,j,1) =0.0; 
    }
   }
    
  // Determinant of metric tensor
  const double g (LargeDisplacementPlateEquations<DIM,NNODE_1D>::
    two_by_two_determinant(g_tensor));
  const double sqrt_g_inv (1/sqrt(g));
  const double sqrt_g_cubed_inv (1/std::pow(g,1.5));

  // Get the metric tensor determinant
  DenseMatrix<double> d_g_du_unknown(this->Number_of_displacements,2,0.0);
  d_metric_tensor_determinant_du_unknown(g_tensor,d_g_tensor_du_unknown,
   d_g_du_unknown);

  // Loop over plane tangent vectors
  double tmp;
  for(unsigned j=0;j<this->Number_of_displacements;++j)
   {
   // Eta parameter
   const double eta_u = this->eta_u();
   for(unsigned alpha=0;alpha<2;++alpha)
    {
    // Beta != Alpha
    const unsigned beta = (alpha + 1) % 2, k = (j+1)%3, i = (j+2)%3 ;

    // Cross product between g_alpha and d_g_beta_du_unknown
    // NB:  dg_beta_duj_unknown[j]  = d2uidx_du_unknown[beta] 
    // Error could be here somewhere 
    tmp = interpolated_drdxi(i,alpha) * eta_u;// * d2uidx_du_unknown[beta];
    d_normal_du_unknown(k,j,beta) +=  tmp * (sqrt_g_inv) * (alpha == 0 ? 1 : -1) ; 

    tmp = interpolated_drdxi(k,beta) * eta_u;// *  d2uidx_du_unknown[alpha];
    d_normal_du_unknown(i,j,alpha) +=  tmp * (sqrt_g_inv) * (alpha == 0 ? 1 : -1) ; 

    // Loop over displacement components
    // And fill in the normal 
    for(unsigned ii=0;ii<this->Number_of_displacements;++ii)
     {
     // The cross product will be epsilon_ijk vj vk
     const unsigned jj=(ii+1)%3,kk=(ii+2)%3;
     // Error could be in 'tmp'
     // Cross product between g_alpha and g_beta
     tmp = interpolated_drdxi(jj,alpha) * interpolated_drdxi(kk,beta) - 
           interpolated_drdxi(kk,alpha) * interpolated_drdxi(jj,beta);
     for(unsigned gamma=0;gamma<2;++gamma)
      {
      d_normal_du_unknown(ii,j,gamma) -= tmp/(4)*sqrt_g_cubed_inv
       * d_g_du_unknown(j,gamma) * (alpha == 0 ? 1 : -1) ; 
      }
     } 
    }
   }
  }

 // Get the (Green Lagrange) strain tensor
 inline void fill_in_strain_tensor(const DenseMatrix<double>&
  interpolated_dudxi, DenseMatrix<double>& strain)
  {
  // Zero the strain
  for(unsigned alpha=0; alpha<2; ++alpha)
   {
   for(unsigned beta=0; beta<2; ++beta)
    { strain(alpha,beta) =0.0; }
   }
  // Loop over alpha and beta
  for(unsigned alpha=0; alpha<2; ++alpha)
   {
   for(unsigned beta=0; beta<2; ++beta)
    {
    // Fill in linear terms
    strain(alpha,beta) +=interpolated_dudxi(alpha,beta)/2.;
    strain(alpha,beta) +=interpolated_dudxi(beta,alpha)/2.;

    // Loop over displacements
    for(unsigned i=0;i<this->Number_of_displacements; ++i)
     {
     // Scale the displacements
     const double eta_u = this->eta_u();
     // Nonlinear terms
     strain(alpha,beta) += eta_u*interpolated_dudxi(i,alpha)
      *interpolated_dudxi(i,beta)/2.;
     }
    }
   } 
  }

// Get the (Green Lagrange) strain tensor
inline void fill_in_d_g_tensor_du_unknown(const DenseMatrix<double>&
 interpolated_dudxi, RankFourTensor<double>& d_g_tensor_du_unknown)
 {
 // Loop over alpha and beta
 for(unsigned alpha=0; alpha<2; ++alpha)
  {
  for(unsigned beta=0; beta<2; ++beta)
   {
   // Zero the vector
   for(unsigned i=0;i<this->Number_of_displacements;++i)
    {
    d_g_tensor_du_unknown(alpha,beta,i,0) = 0.0;
    d_g_tensor_du_unknown(alpha,beta,i,1) = 0.0;
    }
   // Fill in linear terms relating to u_\alpha,\beta
   const double eta_u = this->eta_u();
   d_g_tensor_du_unknown(alpha,beta,beta,alpha)  +=eta_u;
   d_g_tensor_du_unknown(alpha,beta,alpha,beta)  +=eta_u;

   // Loop over displacements
   for(unsigned i=0;i<this->Number_of_displacements; ++i)
    {
    // Scale the displacements
    const double eta_u = this->eta_u();
    // Nonlinear terms
    d_g_tensor_du_unknown(alpha,beta,i,alpha) += eta_u*eta_u*interpolated_dudxi(i
     ,beta);
    d_g_tensor_du_unknown(alpha,beta,i,beta) += eta_u*eta_u*interpolated_dudxi(i
     ,alpha);
    }
   }
  } 
 }

// Get the (Green Lagrange) strain tensor
inline void fill_in_d_strain_tensor_du_unknown(const DenseMatrix<double>&
 interpolated_dudxi, RankFourTensor<double>& d_epsilon_tensor_du_unknown)
 {
 // Loop over alpha and beta
 for(unsigned alpha=0; alpha<2; ++alpha)
  {
  for(unsigned beta=0; beta<2; ++beta)
    {
    // Zero the vector
    for(unsigned i=0;i<this->Number_of_displacements;++i)
     {
      d_epsilon_tensor_du_unknown(alpha,beta,i,0) = 0.0;
      d_epsilon_tensor_du_unknown(alpha,beta,i,1) = 0.0;
     }
   // Fill in linear terms
   d_epsilon_tensor_du_unknown(alpha,beta,beta,alpha) +=1/2.;
   d_epsilon_tensor_du_unknown(alpha,beta,alpha,beta) +=1/2.;

   // Loop over displacements
   for(unsigned i=0;i<this->Number_of_displacements; ++i)
    {
     // Scale the displacements
     const double eta_u = this->eta_u();
     // Nonlinear terms
     d_epsilon_tensor_du_unknown(alpha,beta,i,alpha) += eta_u*interpolated_dudxi(i
      ,beta)/2.;
     d_epsilon_tensor_du_unknown(alpha,beta,i,beta) += eta_u*interpolated_dudxi(i
      ,alpha)/2.;
    }
   }
  } 
 }

 // Fill in the curvature tensor
 inline void fill_in_curvature_tensor(const Vector<double>& unit_normal,
  const DenseMatrix<double>& interpolated_d2rdxi2, DenseMatrix<double>& 
  curvature)
  {
  // Zero the curvature
  for(unsigned alpha=0; alpha<2; ++alpha)
   {
   for(unsigned beta=0; beta<2; ++beta)
    { curvature(alpha,beta) =0.0; }
   }
  // Loop over displacements and then coordinates
  for(unsigned i=0;i<this->Number_of_displacements;++i)
   {
   for(unsigned alpha=0;alpha<2;++alpha)
    {
    for(unsigned beta=0;beta<2;++beta)
     {
      // Scale the displacements
      curvature(alpha,beta) += unit_normal[i]
       *interpolated_d2rdxi2(i,alpha+beta);
     }
    }
   }
  }

 // Fill in the curvature tensor
 inline void fill_in_d_curvature_tensor_du_unknown(const Vector<double>& 
  unit_normal, const DenseMatrix<double>& interpolated_d2rdxi2, const 
  RankThreeTensor<double>&  d_unit_normal_dunknown,  
  RankFourTensor<double>& d_curvature_du_unknown )
  {
  // Zero the vector
  for(unsigned alpha=0; alpha<2; ++alpha)
   {
   for(unsigned beta=0; beta<2; ++beta)
     {
      for(unsigned i=0;i<this->Number_of_displacements;++i)
        {
         // Flatpacked data
         for(unsigned k=0;k<5;++k)
          { d_curvature_du_unknown(alpha,beta,i,k) = 0.0;}
        }
     }
   }
  // Loop over displacements and then coordinates
  for(unsigned i=0;i<this->Number_of_displacements;++i)
   {
   // Scale the displacements
   for(unsigned alpha=0;alpha<2;++alpha)
    {
    for(unsigned beta=0;beta<2;++beta)
     {
      d_curvature_du_unknown(alpha,beta,i,2+alpha+beta) += unit_normal[i];
     for(unsigned j=0;j<this->Number_of_displacements;++j)
      {
       // Note here i is the component (that and we have eta_u_i) and j is
       // unknown
      for(unsigned gamma=0;gamma<2;++gamma)
       d_curvature_du_unknown(alpha,beta,j,gamma) += d_unit_normal_dunknown(i,j,gamma)
        *interpolated_d2rdxi2(i,alpha+beta);
      }
     }
    }
   }
  }

 // Fill in the moment tensor
 inline void fill_in_moment_tensor(const DenseMatrix<double>& interpolated_drxi,
  const Vector<double>& unit_normal,
  const DenseMatrix<double>& curvature, RankThreeTensor<double>& 
  moment)
  {
  // Zero the stress
  for(unsigned i=0; i<this->Number_of_displacements; ++i)
  {
   for(unsigned alpha=0; alpha<2; ++alpha)
    {
    for(unsigned beta=0; beta<2; ++beta)
      { moment(i, alpha,beta) =0.0; }
    }
   }
  // Get thickness
  const double h = this->get_thickness(), nu = this->get_nu();
  // Loop over displacements and then coordinates
  for(unsigned i=0;i<this->Number_of_displacements;++i)
   {
   for(unsigned alpha=0;alpha<2;++alpha)
    {
    for(unsigned beta=0;beta<2;++beta)
     {
      // Get moment i
      moment(i,alpha,beta) += h*h*(1-nu)*unit_normal[i]*curvature(alpha,beta)
            /(12.*(1-nu*nu));
      moment(i,alpha,alpha)+= h*h*nu*unit_normal[i]*curvature(beta,beta)
            /(12.*(1-nu*nu));
     }
    }
   }
  }

 // Fill in the moment tensor
 inline void fill_in_d_moment_tensor_du_unknown(const DenseMatrix<double>&
  interpolated_d2rxi2, const Vector<double>& unit_normal,
  const DenseMatrix<double>& curvature, const RankThreeTensor<double>&
  d_unit_normal_du_unknown,  const RankFourTensor<double>& d_curvature_du_unknown,
  RankFiveTensor<double>&  d_moment_du_unknown)
  {
  // Get thickness
  const double h = this->get_thickness(),nu =this->get_nu();
  // Zero the moment derivative tensor
  for(unsigned i=0;i<this->Number_of_displacements;++i)
   {
   for(unsigned alpha=0;alpha<2;++alpha)
    {
    for(unsigned beta=0;beta<2;++beta)
     { 
     for(unsigned j=0;j<this->Number_of_displacements;++j)
      {
      // Flatpacked first and second derivatives
      for(unsigned k=0;k<5;++k)
       { d_moment_du_unknown(i,alpha,beta,j,k) =0.0;}
      }
     }
    }
   }
  // Loop over displacements and then coordinates
  for(unsigned i=0;i<this->Number_of_displacements;++i)
   {
   for(unsigned alpha=0;alpha<2;++alpha)
    {
    for(unsigned beta=0;beta<2;++beta)
     {
     // Loop over displacements and then coordinates
     for(unsigned j=0;j<this->Number_of_displacements;++j)
      {
      // Get moment i
      // Construct the tensor
      for(unsigned k=0;k<5;++k)
      {
       // First terms
       d_moment_du_unknown(i,alpha,beta,j,k) += (1-nu)* 
          h*h*unit_normal[i]*d_curvature_du_unknown(alpha,beta,j,k)/(12.*(1-nu*nu));
       d_moment_du_unknown(i,alpha,alpha,j,k) += nu* 
          h*h*unit_normal[i]*d_curvature_du_unknown(beta,beta,j,k)/(12.*(1-nu*nu));
      }
      // Second terms
      for(unsigned gamma=0;gamma<2;++gamma)
       {
       d_moment_du_unknown(i,alpha,beta,j,gamma) += (1-nu)* 
          h*h*d_unit_normal_du_unknown(i,j,gamma)*curvature(alpha,beta)/(12.*(1-nu*nu));
       d_moment_du_unknown(i,alpha,alpha,j,gamma) += nu* 
          h*h*d_unit_normal_du_unknown(i,j,gamma)*curvature(beta,beta)/(12.*(1-nu*nu));
       }
      }
     }
    }
   }
  }

 inline void fill_in_total_tension(const DenseMatrix<double>& stress, 
 const RankThreeTensor<double>& christoffel_tensor, const 
 RankThreeTensor<double>&  moment_tensors, const 
 DenseMatrix<double>& interpolated_dudxi,  const DenseMatrix<double>& 
 interpolated_d2rdxi2, DenseMatrix<double>&  tension_vectors)
 {
 // Zero the tension
 for(unsigned i=0;i<this->Number_of_displacements;++i)
  {
  // Loop over inplane components
  for(unsigned gamma=0;gamma<2;++gamma)
   { tension_vectors(i,gamma) =0.0; }
  }

 double delta_ibeta;
 // Loop over displacement components
 for(unsigned i=0;i<this->Number_of_displacements;++i)
  {
  // Loop over inplane components
  for(unsigned gamma=0;gamma<2;++gamma)
   {
   // Loop over inplane components
   for(unsigned beta=0;beta<2;++beta)
    {
    // Kronecker Delta: delta_{i\beta} 
    delta_ibeta = (i==beta? 1.0:0.0);
    // Local copy of eta_u and eta_sigma
    const double eta_u = this->eta_u();
    const double eta_sigma = this->eta_sigma();
    // The diagonal parts of the tangent matrix
    tension_vectors(i,gamma) += (delta_ibeta + eta_u*interpolated_dudxi(i,beta))
     *eta_sigma*stress(gamma,beta);

    for(unsigned alpha=0;alpha<2;++alpha)
     {
     // Shouldn't this be -M_{i \al \be} \Ga^\ga_{\al \be} ?
     tension_vectors(i,gamma) -= eta_u*moment_tensors(i,alpha,beta)*
        christoffel_tensor(gamma,alpha,beta);
     }
    }
   }
  }
 }

 inline void d_fill_in_total_tension_du_unknown(const DenseMatrix<double>& stress
   , const RankThreeTensor<double>& christoffel_tensor, const 
   RankThreeTensor<double>& moment_tensors, const DenseMatrix<double>& interpolated_dudxi, 
   const DenseMatrix<double>& interpolated_d2udxi2,
   RankFourTensor<double>& d_stress_du_unknown, const RankFiveTensor<double>& 
   d_christoffel_tensor_du_unknown, const RankFiveTensor<double>& 
   d_moment_tensors_du_unknown,   RankFourTensor<double>& d_tension_vectors_du_unknown)
 {
  // Zero the derivatives of tension vectors
  for(unsigned i=0;i<this->Number_of_displacements;++i)
   {
   for(unsigned beta=0; beta<2; ++beta)
    {
    for(unsigned j=0;j<this->Number_of_displacements;++j)
     {
     for(unsigned k=0;k<6;++k)
      {
      d_tension_vectors_du_unknown(i,beta,j,k) = 0.0;
      }
     }
    }
   }

  // Loop over displacement components
  for(unsigned i=0;i<this->Number_of_displacements;++i)
   {
   // Because we don't use a d_tangent_dunknown
   //  const double eta_g =(i==2 ? eta_g_z() : eta_g_xy());
   // Loop over inplane components
   double delta_ibeta;
   for(unsigned gamma=0;gamma<2;++gamma)
    {
    // Loop over inplane components
    for(unsigned beta=0;beta<2;++beta)
     {
     // Kronecker Delta: delta_{i\beta} 
     delta_ibeta = (i==beta? 1.0:0.0);
     // Local copy
     const double eta_u = this->eta_u();
     const double eta_sigma = this->eta_sigma();
     // Fill in dtension_du_unknown
     d_tension_vectors_du_unknown(i,gamma,i,1+beta) += eta_sigma*eta_u*stress(gamma,beta);
     for(unsigned j=0;j<this->Number_of_displacements;++j)
      {
      // The diagonal parts of the tangent matrix
      for(unsigned mu=0;mu<3;++mu)
       {
        d_tension_vectors_du_unknown(i,gamma,j,mu) +=
         eta_sigma*d_stress_du_unknown(gamma,beta,j,mu)
          *(delta_ibeta + eta_u* interpolated_dudxi(i,beta));
       }
      }
    for(unsigned alpha=0;alpha<2;++alpha)
     { 
     for(unsigned j=0;j<this->Number_of_displacements;++j)
      {
      for(unsigned k=0;k<5;++k)
       {
       d_tension_vectors_du_unknown(i,gamma,j,1+k) -= eta_u
          *d_moment_tensors_du_unknown(i,alpha,beta,j,k)
          *christoffel_tensor(gamma,alpha,beta);
       d_tension_vectors_du_unknown(i,gamma,j,1+k) -= eta_u
          *moment_tensors(i,alpha,beta)
          *d_christoffel_tensor_du_unknown(gamma,alpha,beta,j,k);
       }
      }
     }
    }
   }
  }
 }

 // Get the Christtoffel tensor of the second kind
 // \Gamma_{\alpha\beta\gamma} = ( E_{\alpha\beta,gamma} + E_{\gamma\alpha,\beta}
 //   - E_{\beta\gamma,\alpha})
 inline void fill_in_second_christoffel_tensor(const DenseMatrix<double>&
  interpolated_dudxi, const DenseMatrix<double>& interpolated_d2udxi2, 
  RankThreeTensor<double>& gamma_tensor)
  {
  // Zero the christoffel
  for(unsigned gamma=0; gamma<2; ++gamma)
   {
   for(unsigned alpha=0; alpha<2; ++alpha)
    {
    for(unsigned beta=0; beta<2; ++beta)
     { 
      gamma_tensor(alpha,beta,gamma)=0.0;
     }
    }
   }

   // An intermediate step
   RankThreeTensor<double> strain_gradient(2,2,2,0.0);

  // Loop over alpha and beta
  for(unsigned gamma=0; gamma<2; ++gamma)
   {
   for(unsigned alpha=0; alpha<2; ++alpha)
    {
    for(unsigned beta=0; beta<2; ++beta)
     {
     // Fill in linear terms
     strain_gradient(alpha,beta,gamma) += interpolated_d2udxi2
       (alpha,beta+gamma)/2. + interpolated_d2udxi2(beta,alpha+gamma)/2.;
     // Loop over displacements
     for(unsigned i=0;i<this->Number_of_displacements; ++i)
      {
      // Scale the displacements
      const double eta_u = this->eta_u();
      // Nonlinear terms
      strain_gradient(alpha,beta,gamma) += eta_u*(interpolated_dudxi(i,alpha)
       *interpolated_d2udxi2(i,beta+gamma))/2.;
      strain_gradient(alpha,beta,gamma) += eta_u*(interpolated_dudxi(i,beta)
       *interpolated_d2udxi2(i,alpha+gamma))/2.;
      }
     }
    } 
   }

   // Fill in Christoffel
   for(unsigned alpha=0; alpha<2; ++alpha)
    {
    for(unsigned beta=0; beta<2; ++beta)
     {
     for(unsigned gamma=0; gamma<2; ++gamma)
      {
      gamma_tensor(alpha,beta,gamma) += strain_gradient(alpha,beta,gamma)+
         strain_gradient(gamma,alpha,beta) - strain_gradient(gamma,beta,alpha); 
      }
     }
    }
  }
 
 // Get the Christtoffel tensor of the second kind
 // \Gamma_{\alpha\beta\gamma} = ( E_{\alpha\beta,gamma} + E_{\gamma\alpha,\beta}
 //   - E_{\beta\gamma,\alpha})
 inline void fill_in_d_second_christoffel_tensor_dui_unknown(const DenseMatrix<double>&
  interpolated_dudxi, const DenseMatrix<double>& interpolated_d2udxi2, 
  RankFiveTensor<double>& d_gamma_tensor_du_unknown)
  {
   // An intermediate step
   RankFiveTensor<double> d_strain_gradient_du_unknown(2,2,2,3,5,0.0);

  // Zero
  for(unsigned gamma=0; gamma<2; ++gamma)
   {
   for(unsigned alpha=0; alpha<2; ++alpha)
    {
    for(unsigned beta=0; beta<2; ++beta)
     {
     for(unsigned i=0; i<this->Number_of_displacements; ++i)
      {
      for(unsigned k=0; k<5; ++k)
       { d_gamma_tensor_du_unknown(alpha,beta,gamma,i,k) = 0.0; }
      }
     }
    }
   }

  // Loop over alpha and beta
  for(unsigned gamma=0; gamma<2; ++gamma)
   {
   for(unsigned alpha=0; alpha<2; ++alpha)
    {
    for(unsigned beta=0; beta<2; ++beta)
     {
     // Fill in linear terms
     d_strain_gradient_du_unknown(alpha,beta,gamma,alpha,2+beta+gamma) += 1/2. ;
     d_strain_gradient_du_unknown(alpha,beta,gamma,beta,2+alpha+gamma) += 1/2. ;
     // Loop over displacements
     for(unsigned i=0;i<this->Number_of_displacements; ++i)
      {
       // Scale the displacements
       const double eta_u =this->eta_u();
       // Nonlinear terms
       d_strain_gradient_du_unknown(alpha,beta,gamma,i,alpha) += eta_u*
        (interpolated_d2udxi2(i,beta+gamma))/2.;
       d_strain_gradient_du_unknown(alpha,beta,gamma,i,beta) += eta_u*
        (interpolated_d2udxi2(i,alpha+gamma))/2.;
       d_strain_gradient_du_unknown(alpha,beta,gamma,i,2+beta+gamma) += eta_u*
         (interpolated_dudxi(i,alpha))/2. ;
       d_strain_gradient_du_unknown(alpha,beta,gamma,i,2+alpha+gamma) += eta_u*
         (interpolated_dudxi(i,beta))/2. ;
      }
     }
    } 
   }
  // Fill in Christoffel
  for(unsigned i=0; i<this->Number_of_displacements; ++i)
   {
   for(unsigned alpha=0; alpha<2; ++alpha)
    {
    for(unsigned beta=0; beta<2; ++beta)
     {
     for(unsigned gamma=0; gamma<2; ++gamma)
      {
      for(unsigned k=0; k<5; ++k)
       {
       d_gamma_tensor_du_unknown(alpha,beta,gamma,i,k) +=
         + d_strain_gradient_du_unknown(alpha,beta,gamma,i,k)
         + d_strain_gradient_du_unknown(gamma,alpha,beta,i,k) 
         - d_strain_gradient_du_unknown(gamma,beta,alpha,i,k); 
       }
      }
     }
    }
   }
  }


 /// \short get the coordinate
 virtual void get_coordinate_x(const Vector<double>& s, Vector<double>& x) const=0;

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

 static const double Default_Eta_Value;
};

} //end namespace oomph
#endif


// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
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
// Config header generated by autoconfig
#ifdef HAVE_CONFIG_H
#include <oomph-lib-config.h>
#endif


// oomph-lib includes
#include "biharmonic_preconditioner.h"


namespace oomph
{
#ifdef OOMPH_HAS_HYPRE

  //=============================================================================
  // defaults settings for the Hypre solver (AMG) when used as the approximate
  // linear solver for the Schur complement (non-compound) linear subsidiary
  // linear systems
  //=============================================================================
  namespace Biharmonic_schur_complement_Hypre_defaults
  {
    /// smoother type - Gauss Seidel: 1
    unsigned AMG_smoother = 1;

    /// amg coarsening strategy: classical Ruge Stueben: 1
    unsigned AMG_coarsening = 1;

    /// number of V cycles: 2
    unsigned N_cycle = 2;

    /// amg strength parameter: 0.25 -- optimal for 2d
    double AMG_strength = 0.25;

    /// jacobi damping -- hierher not used 0.1
    double AMG_jacobi_damping = 0.1;

    /// amg smoother iterations
    unsigned AMG_smoother_iterations = 2;

    /// set the defaults
    void set_defaults(HyprePreconditioner* hypre_prec_pt)
    {
      // use AMG preconditioner
      hypre_prec_pt->hypre_method() = HypreSolver::BoomerAMG;

      // Smoother types
      hypre_prec_pt->amg_simple_smoother() = AMG_smoother;

      // jacobi damping
      //     hypre_prec_pt->amg_damping() = AMG_jacobi_damping;

      // coarsening stategy
      hypre_prec_pt->amg_coarsening() = AMG_coarsening;

      oomph_info << "Current number of v cycles: "
                 << hypre_prec_pt->amg_iterations() << std::endl;

      // number of v-cycles
      hypre_prec_pt->amg_iterations() = N_cycle;

      oomph_info << "Re-assigned number of v cycles: "
                 << hypre_prec_pt->amg_iterations() << std::endl;

      // strength parameter
      hypre_prec_pt->amg_strength() = AMG_strength;

      // hierher new
      oomph_info << "Current number of amg smoother iterations: "
                 << hypre_prec_pt->amg_smoother_iterations() << std::endl;

      hypre_prec_pt->amg_smoother_iterations() = AMG_smoother_iterations;

      oomph_info << "Re-assigned number of amg smoother iterations: "
                 << hypre_prec_pt->amg_smoother_iterations() << std::endl;
    }
  } // namespace Biharmonic_schur_complement_Hypre_defaults
#endif

  //===========================================================================
  /// setup for the biharmonic preconditioner
  //===========================================================================
  void BiharmonicPreconditioner::setup()
  {
    // clean up
    this->clean_up_memory();

    // paranoid check that teh bulk element mesh has been set
#ifdef PARANOID
    if (Bulk_element_mesh_pt == 0)
    {
      std::ostringstream error_message;
      error_message << "The bulk element mesh has not been passed to "
                       "bulk_element_mesh_pt()";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // setup the mesh
    this->set_mesh(0, Bulk_element_mesh_pt);

    // setup the blocks look up schemes
    this->block_setup();

    // determine whether this preconditioner has 4 or 5 block types and set
    // Nblock_types if neccessary
    // unsigned n_row = this->master_nrow();
    // bool nblock_type_check = true;
    // for (unsigned i = 0; i < n_row; i++)
    //  {
    //   if (this->block_number(i) == 4) { nblock_type_check = false; }
    //  }
    // if (nblock_type_check) { Nblock_types = 4; }
    //

    // check the preconditioner type is acceptable
#ifdef PARANOID
    if (Preconditioner_type != 0 && Preconditioner_type != 1 &&
        Preconditioner_type != 2 && Preconditioner_type != 3)
    {
      std::ostringstream error_message;
      error_message << "Preconditioner_type must be equal to 0 (BBD exact), 1 "
                       "(inexact BBD with LU),"
                    << " 2 (inexact BBD with AMG) or 3 (exact BD).";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // create the preconditioners
    bool use_amg = true;
    bool retain_all_blocks = false;
    switch (Preconditioner_type)
    {
        // Exact BBD
      case 0:

        retain_all_blocks = false;
        Sub_preconditioner_1_pt =
          new ExactSubBiharmonicPreconditioner(this, retain_all_blocks);
        Sub_preconditioner_2_pt = new SuperLUPreconditioner;

        oomph_info << "Using exact BBD\n";
        break;

        // Inexact BBD with LU
      case 1:

        use_amg = false;
        Sub_preconditioner_1_pt =
          new InexactSubBiharmonicPreconditioner(this, use_amg);
        Sub_preconditioner_2_pt = new MatrixBasedDiagPreconditioner;
        oomph_info << "Using inexact BBD with LU\n";
        break;


        // Inexact BBD with AMG
      case 2:

        use_amg = true;
        Sub_preconditioner_1_pt =
          new InexactSubBiharmonicPreconditioner(this, use_amg);
        Sub_preconditioner_2_pt = new MatrixBasedDiagPreconditioner;
        oomph_info << "Using inexact BBD with AMG\n";
        break;

        /// Exact BD
      case 3:

        retain_all_blocks = true;
        Sub_preconditioner_1_pt =
          new ExactSubBiharmonicPreconditioner(this, retain_all_blocks);
        Sub_preconditioner_2_pt = new SuperLUPreconditioner;

        oomph_info << "Using exact BD\n";
        break;

      default:

        throw OomphLibError("Wrong type of preconditioner.",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
    }


    // setup sub preconditioner pt 1
    Sub_preconditioner_1_pt->setup(matrix_pt());

    // get the matrix ans setup sub preconditioner pt 2
    CRDoubleMatrix* j_33_pt = new CRDoubleMatrix;
    this->get_block(3, 3, *j_33_pt);
    Sub_preconditioner_2_pt->setup(j_33_pt);
    delete j_33_pt;
    j_33_pt = 0;

    // if the block preconditioner has 5 block types setup the preconditioner
    // for the 5th block diagonal block (Matrix is also diagonal hence a
    // diagonal preconditioner is sufficient in the exact biharmonic
    // preconditioner case as well)
    if (this->nblock_types() == 5)
    {
      // get the matrix for block J_33
      CRDoubleMatrix* j_44_pt = new CRDoubleMatrix;
      this->get_block(4, 4, *j_44_pt);

      // setup the hijacked sub preconditioner
      Hijacked_sub_block_preconditioner_pt = new MatrixBasedDiagPreconditioner;
      Hijacked_sub_block_preconditioner_pt->setup(j_44_pt);
      delete j_44_pt;
      j_44_pt = 0;
    }
  }


  //============================================================================
  /// preconditioner solve for the biharmonic preconditioner
  //============================================================================
  void BiharmonicPreconditioner::preconditioner_solve(const DoubleVector& r,
                                                      DoubleVector& z)
  {
    // zero z
    z.initialise(0.0);

    // solve sub preconditioner 1
    Sub_preconditioner_1_pt->preconditioner_solve(r, z);

    // solve sub preconditioner 2
    DoubleVector block_r;
    get_block_vector(3, r, block_r);
    DoubleVector block_z;
    Sub_preconditioner_2_pt->preconditioner_solve(block_r, block_z);
    return_block_vector(3, block_z, z);

    // solve the hijacked sub block preconditioner if required
    if (this->nblock_types() == 5)
    {
      block_r.clear();
      block_z.clear();
      get_block_vector(4, r, block_r);
      Hijacked_sub_block_preconditioner_pt->preconditioner_solve(block_r,
                                                                 block_z);
      return_block_vector(4, block_z, z);
    }
  }

  //============================================================================
  /// setup for the exact sub biharmonic preconditioner
  //============================================================================
  void ExactSubBiharmonicPreconditioner::setup()
  {
    // clean up memory first
    this->clean_up_memory();

    // setup
    this->block_setup();

    // Number of block types
    unsigned n_block_types = this->nblock_types();

    // check for required number of blocks
#ifdef PARANOID
    if (n_block_types != 3)
    {
      std::ostringstream error_message;
      error_message
        << "This preconditioner requires 3 block types.\n"
        << "It is sub preconditioner for the BiharmonicPreconditioner.\n";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Data type indicating which blocks from the preconditioner matrix we want
    VectorMatrix<BlockSelector> required_blocks(n_block_types, n_block_types);

    // boolean indicating if we want the block or not, stored for readability.
    // Initially this is set to true for all blocks. Later we select which
    // blocks we do not want.
    const bool want_block = true;
    for (unsigned b_i = 0; b_i < n_block_types; b_i++)
    {
      for (unsigned b_j = 0; b_j < n_block_types; b_j++)
      {
        required_blocks[b_i][b_j].select_block(b_i, b_j, want_block);
      }
    }

    // Which blocks do we not want?
    if (!Retain_all_blocks)
    {
      required_blocks[1][2].do_not_want_block();
      required_blocks[2][1].do_not_want_block();
    }

    // Get the preconditioner matrix as defined by required_blocks
    CRDoubleMatrix preconditioner_matrix =
      this->get_concatenated_block(required_blocks);

    // setup the preconditioner
    Sub_preconditioner_pt = new SuperLUPreconditioner;
    Sub_preconditioner_pt->setup(&preconditioner_matrix);

    // preconditioner_matrix will now go out of scope (and is destroyed).
  }


  //============================================================================
  /// preconditioner solve for the exact sub biharmonic preconditioner
  //============================================================================
  void ExactSubBiharmonicPreconditioner::preconditioner_solve(
    const DoubleVector& r, DoubleVector& z)
  {
    // vectors for use within the sub preconditioner
    DoubleVector sub_r;
    DoubleVector sub_z;

    // get the sub r vector
    get_block_ordered_preconditioner_vector(r, sub_r);

    // solve the preconditioner
    Sub_preconditioner_pt->preconditioner_solve(sub_r, sub_z);

    // return the sub z vector to the master z vector
    return_block_ordered_preconditioner_vector(sub_z, z);
  }


  //============================================================================
  /// setup for the inexact sub biharmonic preconditioner
  //============================================================================
  void InexactSubBiharmonicPreconditioner::setup()
  {
    // clean up memory first
    this->clean_up_memory();

    // setup
    this->block_setup();

    // Number of block types
    unsigned n_block_types = this->nblock_types();

    // paranoid check for number of blocks
#ifdef PARANOID
    if (n_block_types != 3)
    {
      std::ostringstream error_message;
      error_message
        << "This preconditioner requires 3 block types.\n"
        << "It is sub preconditioner for the BiharmonicPreconditioner.\n";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // required blocks
    DenseMatrix<bool> required_blocks(n_block_types, n_block_types, false);
    required_blocks(0, 0) = true;
    required_blocks(0, 1) = true;
    required_blocks(1, 0) = true;
    required_blocks(1, 1) = true;
    required_blocks(0, 2) = true;
    required_blocks(2, 0) = true;
    required_blocks(2, 2) = true;

    // Matrix of block matrix pointers
    Matrix_of_block_pointers.resize(n_block_types, n_block_types, 0);

    // get the blocks
    this->get_blocks(required_blocks, Matrix_of_block_pointers);

    // lump the matrix J_11
    Lumped_J_11_preconditioner_pt =
      new MatrixBasedLumpedPreconditioner<CRDoubleMatrix>;
    Lumped_J_11_preconditioner_pt->setup(Matrix_of_block_pointers(1, 1));

    delete Matrix_of_block_pointers(1, 1);
    Matrix_of_block_pointers(1, 1) = 0;

    // lump the matrix J_22
    Lumped_J_22_preconditioner_pt =
      new MatrixBasedLumpedPreconditioner<CRDoubleMatrix>;
    Lumped_J_22_preconditioner_pt->setup(Matrix_of_block_pointers(2, 2));
    delete Matrix_of_block_pointers(2, 2);
    Matrix_of_block_pointers(2, 2) = 0;

    // compute the schur complement
    compute_inexact_schur_complement();

    // create the preconditioner for the S00 Schur complement linear system
    if (Use_amg)
    {
#ifdef OOMPH_HAS_HYPRE
      // Use Hypre Boomer AMG
      S_00_preconditioner_pt = new HyprePreconditioner;
      Biharmonic_schur_complement_Hypre_defaults::set_defaults(
        static_cast<HyprePreconditioner*>(S_00_preconditioner_pt));
#else
      std::ostringstream error_message;
      error_message << "Request AMG solver but oomph-lib does not have HYPRE";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
#endif
    }
    else
    {
      S_00_preconditioner_pt = new SuperLUPreconditioner;
    }

    // setup the preconditioner
    S_00_preconditioner_pt->setup(S_00_pt);

    // clean up
    delete S_00_pt;
    S_00_pt = 0;
  }

  //============================================================================
  /// computes the schur complement for the inexact sub biharmonic
  /// preconditioner
  //============================================================================
  void InexactSubBiharmonicPreconditioner::compute_inexact_schur_complement()
  {
    // if required get pointers to the vector components of J01 and J10
    int* J_01_row_start = 0;
    int* J_01_column_index = 0;
    double* J_01_value = 0;
    int* J_10_row_start = 0;
    int* J_10_column_index = 0;

    // J_01 matrix
    J_01_row_start = Matrix_of_block_pointers(0, 1)->row_start();
    J_01_column_index = Matrix_of_block_pointers(0, 1)->column_index();
    J_01_value = Matrix_of_block_pointers(0, 1)->value();

    // J_10 matrix
    J_10_row_start = Matrix_of_block_pointers(1, 0)->row_start();
    J_10_column_index = Matrix_of_block_pointers(1, 0)->column_index();

    // if required get pointers to the vector components of J01 and J10
    int* J_02_row_start = 0;
    int* J_02_column_index = 0;
    double* J_02_value = 0;
    int* J_20_row_start = 0;
    int* J_20_column_index = 0;

    // J_02 matrix
    J_02_row_start = Matrix_of_block_pointers(0, 2)->row_start();
    J_02_column_index = Matrix_of_block_pointers(0, 2)->column_index();
    J_02_value = Matrix_of_block_pointers(0, 2)->value();

    // J_20 matrix
    J_20_row_start = Matrix_of_block_pointers(2, 0)->row_start();
    J_20_column_index = Matrix_of_block_pointers(2, 0)->column_index();

    // get the inverse lumped vector of J_11 if required
    double* J_11_lumped_and_inverted = 0;
    J_11_lumped_and_inverted =
      Lumped_J_11_preconditioner_pt->inverse_lumped_vector_pt();

    // get the inverse lumped vector of J_22 if required
    double* J_22_lumped_and_inverted = 0;
    J_22_lumped_and_inverted =
      Lumped_J_22_preconditioner_pt->inverse_lumped_vector_pt();

    // size of J00 matrix (and S00 matrix)
    unsigned J_00_nrow = Matrix_of_block_pointers(0, 0)->nrow();

    // vectors for the schur complement
    Vector<int> S_00_row_start(J_00_nrow + 1);
    Vector<int> S_00_column_index;
    Vector<double> S_00_value;

    // number of elements in the x-dimension of the mesh
    unsigned n_element_x =
      dynamic_cast<HermiteQuadMesh<Hijacked<BiharmonicElement<2>>>*>(
        dynamic_cast<BiharmonicPreconditioner*>(
          this->master_block_preconditioner_pt())
          ->bulk_element_mesh_pt())
        ->nelement_in_dim(0);

    // nnz in schur complement (initialised to zero)
    unsigned S_00_nnz = 0;

    // loop over columns of schur complement matrix
    for (unsigned i = 0; i < J_00_nrow; i++)
    {
      // set column_start
      S_00_row_start[i] = S_00_nnz;

      // loop over rows in schur complement matrix
      // the schur complement matrix has 5 well defined bands thus we only
      // perform matrix-matrix multiplication for these bands
      //
      // where the diagonal band is 0 :
      //
      //   band 1 :   -2*n_element_x +/- 5
      //        2 :   -n_element_x +/- 3
      //        3 :   0 +/- 3
      //        4 :   n_element_x +/- 3
      //        5 :   2*n_element_x +/- 5
      //
      // regardless of the type or combination of boundary conditions applied

      // Vector for postion of the bands in S_00
      Vector<std::pair<int, int>> band_position(5);

      // compute the minimum and maximum positions of each band in terms of
      // row number for column j
      // note : static_cast used because max and min don't work on unsigned
      band_position[0].first =
        std::max(0, static_cast<int>(i - n_element_x * 2 - 5));
      band_position[0].second =
        std::max(0,
                 std::min(static_cast<int>(J_00_nrow - 1),
                          static_cast<int>(i - n_element_x * 2 + 5)));
      band_position[1].first =
        std::max(band_position[0].second + 1,
                 std::max(0, static_cast<int>(i - n_element_x - 3)));
      band_position[1].second =
        std::max(0,
                 std::min(static_cast<int>(J_00_nrow - 1),
                          static_cast<int>(i - n_element_x + 3)));
      band_position[2].first = std::max(band_position[1].second + 1,
                                        std::max(0, static_cast<int>(i - 3)));
      band_position[2].second = std::max(
        0, std::min(static_cast<int>(J_00_nrow - 1), static_cast<int>(i + 3)));
      band_position[3].first =
        std::max(band_position[2].second + 1,
                 std::max(0, static_cast<int>(i + n_element_x - 3)));
      band_position[3].second =
        std::max(0,
                 std::min(static_cast<int>(J_00_nrow - 1),
                          static_cast<int>(i + n_element_x + 3)));
      band_position[4].first =
        std::max(band_position[3].second + 1,
                 std::max(0, static_cast<int>(i + n_element_x * 2 - 5)));
      band_position[4].second =
        std::max(0,
                 std::min(static_cast<int>(J_00_nrow - 1),
                          static_cast<int>(i + n_element_x * 2 + 5)));

      // number of bands
      unsigned n_band = 5;

      // loop over the bands
      for (unsigned b = 0; b < n_band; b++)
      {
        // loop over the rows in band b
        for (unsigned j = static_cast<unsigned>(band_position[b].first);
             j <= static_cast<unsigned>(band_position[b].second);
             j++)
        {
          ;

          // temporary value for the computation of S00(i,j)
          double temp_value = Matrix_of_block_pointers(0, 0)->operator()(i, j);

          // iterate through non-zero entries of  column j of A_10
          for (int k = J_01_row_start[i]; k < J_01_row_start[i + 1]; k++)
          {
            if (J_10_column_index[J_10_row_start[J_01_column_index[k]]] <=
                  static_cast<int>(j) &&
                static_cast<int>(j) <=
                  J_10_column_index[J_10_row_start[J_01_column_index[k] + 1] -
                                    1])
            {
              temp_value -= J_01_value[k] *
                            Matrix_of_block_pointers(1, 0)->operator()(
                              J_01_column_index[k], j) *
                            J_11_lumped_and_inverted[J_01_column_index[k]];
            }
          }

          // next compute contribution for A_02*lumped(A_22)'*A_20

          // iterate through non-zero entries of  column j of A_10
          for (int k = J_02_row_start[i]; k < J_02_row_start[i + 1]; k++)
          {
            if (J_20_column_index[J_20_row_start[J_02_column_index[k]]] <=
                  static_cast<int>(j) &&
                static_cast<int>(j) <=
                  J_20_column_index[J_20_row_start[J_02_column_index[k] + 1] -
                                    1])
            {
              temp_value -= J_02_value[k] *
                            Matrix_of_block_pointers(2, 0)->operator()(
                              J_02_column_index[k], j) *
                            J_22_lumped_and_inverted[J_02_column_index[k]];
            }
          }

          // add element to schur complement matrix S00
          if (temp_value != 0.0)
          {
            S_00_nnz++;
            S_00_value.push_back(temp_value);
            S_00_column_index.push_back(j);
          }
        }
      }
    }

    // last entry of s00 column start
    S_00_row_start[J_00_nrow] = S_00_nnz;

    // build the schur complement S00
    S_00_pt = new CRDoubleMatrix(this->block_distribution_pt(0),
                                 J_00_nrow,
                                 S_00_value,
                                 S_00_column_index,
                                 S_00_row_start);

    // replace block J01 with J01*lumped(J11)' (if J11 can be lumped)
    unsigned J_01_nnz = Matrix_of_block_pointers(0, 1)->nnz();
    for (unsigned i = 0; i < J_01_nnz; i++)
    {
      J_01_value[i] *= J_11_lumped_and_inverted[J_01_column_index[i]];
    }

    // replace block J_02 with J_02*lumped(J_22)' (if J22 can be lumped)
    unsigned J_02_nnz = Matrix_of_block_pointers(0, 2)->nnz();
    for (unsigned i = 0; i < J_02_nnz; i++)
    {
      J_02_value[i] *= J_22_lumped_and_inverted[J_02_column_index[i]];
    }
  }


  //============================================================================
  /// preconditioner solve for the inexact sub biharmonic preconditioner
  //============================================================================
  void InexactSubBiharmonicPreconditioner::preconditioner_solve(
    const DoubleVector& r, DoubleVector& z)
  {
    // get the block vectors
    Vector<DoubleVector> block_r(3);
    get_block_vectors(r, block_r);

    // r_0 = r_0 - J_01 * lumped(J_11)'*r_1 - J_02 * lumped(J_22)'*r_2
    // Remember that J_01 has already been premultiplied by lumped(J_11)
    DoubleVector temp;
    Matrix_of_block_pointers(0, 1)->multiply(block_r[1], temp);
    block_r[0] -= temp;
    temp.clear();
    Matrix_of_block_pointers(0, 2)->multiply(block_r[2], temp);
    block_r[0] -= temp;

    // apply the inexact preconditioner
    temp.clear();
    S_00_preconditioner_pt->preconditioner_solve(block_r[0], temp);
    return_block_vector(0, temp, z);

    // solve: lumped(J_11) x_1 = r_1 - J_10 x_0 for x_1
    // remember temp contains r_0 (...or z_0)
    DoubleVector temp2;
    Matrix_of_block_pointers(1, 0)->multiply(temp, temp2);
    block_r[1] -= temp2;
    DoubleVector z_1;
    Lumped_J_11_preconditioner_pt->preconditioner_solve(block_r[1], z_1);
    return_block_vector(1, z_1, z);

    // solve: lumped(J_22) x_2 = r_2 - J_20 x_0 for x_2
    // remember temp contains r_0 (...or z_0)
    temp2.clear();
    Matrix_of_block_pointers(2, 0)->multiply(temp, temp2);
    block_r[2] -= temp2;
    DoubleVector z_2;
    Lumped_J_22_preconditioner_pt->preconditioner_solve(block_r[2], z_2);
    return_block_vector(2, z_2, z);
  }
} // namespace oomph

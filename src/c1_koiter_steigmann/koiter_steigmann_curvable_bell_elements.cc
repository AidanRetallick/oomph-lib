// Non--inline functions for BellBiharmonic elements
#include "koiter_steigmann_curvable_bell_elements.h"

namespace oomph
{
  //=============================================================================
  /// Set the number of fields
  //=============================================================================
  const unsigned KoiterSteigmannC1CurvableBellElement::Nfield = 3;


  //=============================================================================
  /// Set the interpolation of each field
  //=============================================================================
  const std::vector<bool> KoiterSteigmannC1CurvableBellElement::
  Field_is_bell_interpolated = {true, true, true};

  //======================================================================
  /// Set the data for the number of Variables at each node
  //======================================================================
  // template<unsigned DIM,
  //          unsigned NNODE_1D,
  //          unsigned BOUNDARY_ORDER,
  //          template<unsigned DIM_, unsigned NNODE_1D_>
  //          class PLATE_EQUATIONS>
  const unsigned
    KoiterSteigmannC1CurvableBellElement::Initial_Nvalue = 18;


  //====================================================================
  // Force build of templates
  //====================================================================
  // template class LargeDisplacementPlateC1CurvedBellElement<
  //   2,
  //   2,
  //   3,
  //   KoiterSteigmannPlateEquations>;
  // template class LargeDisplacementPlateC1CurvedBellElement<
  //   2,
  //   2,
  //   5,
  //   KoiterSteigmannPlateEquations>;
  // template class LargeDisplacementPlateC1CurvedBellElement<
  //   2,
  //   2,
  //   3,
  //   FoepplVonKarmanCorrectionEquations>;
  // template class LargeDisplacementPlateC1CurvedBellElement<
  //   2,
  //   2,
  //   5,
  //   FoepplVonKarmanCorrectionEquations>;

} // namespace oomph

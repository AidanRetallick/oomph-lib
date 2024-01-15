// Non--inline functions for BellBiharmonic elements
#include "large_displacement_plate_models.h"

namespace oomph
{

//======================================================================
/// Set the data for the number of Variables at each node
//======================================================================
template <unsigned DIM, unsigned NNODE_1D, unsigned BOUNDARY_ORDER, template
<unsigned DIM_, unsigned NNODE_1D_> class PLATE_EQUATIONS >
 const unsigned LargeDisplacementPlateC1CurvedBellElement<DIM,NNODE_1D,BOUNDARY_ORDER,PLATE_EQUATIONS>
  ::Initial_Nvalue = 18;


//====================================================================
// Force build of templates
//====================================================================
template class LargeDisplacementPlateC1CurvedBellElement<2,2,3,KoiterSteigmannPlateEquations>;
template class LargeDisplacementPlateC1CurvedBellElement<2,2,5,KoiterSteigmannPlateEquations>;
template class LargeDisplacementPlateC1CurvedBellElement<2,2,3,FoepplVonKarmanCorrectionEquations>;
template class LargeDisplacementPlateC1CurvedBellElement<2,2,5,FoepplVonKarmanCorrectionEquations>;

}

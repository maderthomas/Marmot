#include "Marmot/MarmotBondSlipLawFactory.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotTesting.h"
#include <memory>

using namespace Marmot;
using namespace Marmot::Testing;

void testLinearElasticBondSlip()
{
  const double props[] = { 200., 5000. };
  auto         law     = std::unique_ptr< MarmotBondSlipLaw >(
    MarmotLibrary::MarmotBondSlipLawFactory::createBondSlipLaw( "LINEARELASTICBONDSLIP", props, 2, 1 ) );

  MarmotBondSlipLaw::State state{ Vector3d::Zero(), 0.0, 0.0, nullptr };
  Matrix3d                 C;
  const Vector3d           slip( 0.1, -0.02, 0.03 );
  law->computeBondStress( state, C, slip, slip, { 0.0, 1.0 } );

  throwExceptionOnFailure( checkIfEqual< double >( state.bondStress, Vector3d( 20., -100., 150. ), 1e-12 ),
                           "linear bond stress" );
  throwExceptionOnFailure( checkIfEqual< double >( C, Matrix3d( Vector3d( 200., 5000., 5000. ).asDiagonal() ), 1e-12 ),
                           "linear bond tangent" );
  throwExceptionOnFailure( checkIfEqual( state.elasticEnergyDensity, 0.5 * slip.dot( state.bondStress ), 1e-12 ),
                           "elastic energy" );
  throwExceptionOnFailure( law->getNumberOfRequiredStateVars() == 0, "no state variables" );
}

void testUnknownBondSlipLawThrows()
{
  const double props[] = { 1., 1. };
  bool         threw   = false;
  try {
    MarmotLibrary::MarmotBondSlipLawFactory::createBondSlipLaw( "NOTABONDSLIPLAW", props, 2, 1 );
  }
  catch ( const std::invalid_argument& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw, "an unknown bond-slip law must be rejected" );
}

int main()
{
  executeTestsAndCollectExceptions( { testLinearElasticBondSlip, testUnknownBondSlipLawThrows } );
  return 0;
}

#include "Marmot/BeamElement.h"
#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

using namespace Marmot;
using namespace Marmot::Testing;

namespace {

  /// A single beam element with its own state, for tests
  struct Beam {
    std::unique_ptr< MarmotElement > element;
    std::vector< double >            coordinates;
    std::vector< double >            properties;
    std::vector< double >            materialProperties;
    std::vector< double >            stateVars;
    int                              nDof;

    Beam( const std::string&    type,
          std::vector< double > coordinates_,
          std::vector< double > properties_,
          const std::string&    material,
          std::vector< double > materialProperties_ )
      : element( MarmotLibrary::MarmotElementFactory::createElement( type, 1 ) ),
        coordinates( std::move( coordinates_ ) ),
        properties( std::move( properties_ ) ),
        materialProperties( std::move( materialProperties_ ) )
    {
      nDof = element->getNDofPerElement();
      element->assignNodeCoordinates( coordinates.data() );
      element->assignProperty( ElementProperties( properties.data(), static_cast< int >( properties.size() ) ) );
      element->initializeYourself();
      element->assignProperty(
        MarmotMaterialSection( material, materialProperties.data(), static_cast< int >( materialProperties.size() ) ) );
      stateVars.assign( element->getNumberOfRequiredStateVars(), 0.0 );
      element->assignStateVars( stateVars.data(), static_cast< int >( stateVars.size() ) );
      element->setInitialConditions( MarmotElement::MarmotMaterialInitialization, nullptr );
    }

    /// residual and tangent for the total dofs U from the committed state; commits if requested
    std::pair< Eigen::VectorXd, Eigen::MatrixXd > kernels( const Eigen::VectorXd& U,
                                                           const Eigen::VectorXd& UOld,
                                                           bool                   commit )
    {
      const auto      backup = stateVars;
      Eigen::VectorXd P      = Eigen::VectorXd::Zero( nDof );
      Eigen::MatrixXd K      = Eigen::MatrixXd::Zero( nDof, nDof );
      Eigen::VectorXd dU     = U - UOld;
      element->computeKernels( U.data(), dU.data(), P.data(), K.data(), 0.0, 1.0 );
      if ( !commit )
        stateVars = backup;
      return { P, K };
    }

    /// central difference tangent from the committed state
    Eigen::MatrixXd numericalTangent( const Eigen::VectorXd& U, const Eigen::VectorXd& UOld, double h )
    {
      Eigen::MatrixXd K( nDof, nDof );
      for ( int j = 0; j < nDof; j++ ) {
        Eigen::VectorXd Up = U, Um = U;
        Up( j ) += h;
        Um( j ) -= h;
        K.col( j ) = ( kernels( Up, UOld, false ).first - kernels( Um, UOld, false ).first ) / ( 2 * h );
      }
      return K;
    }

    /// the linear elastic solution of node 1 clamped, node 2 loaded by f
    Eigen::VectorXd cantilever( const Eigen::VectorXd& f )
    {
      const int       n = nDof / 2;
      const auto      K = kernels( Eigen::VectorXd::Zero( nDof ), Eigen::VectorXd::Zero( nDof ), false ).second;
      Eigen::VectorXd u = K.bottomRightCorner( n, n ).ldlt().solve( f );
      return u;
    }

    Eigen::VectorXd state( const std::string& name, int qp = 0 )
    {
      const auto view = element->getStateView( name, qp );
      return Eigen::Map< const Eigen::VectorXd >( view.stateLocation, view.stateSize );
    }
  };

  const double                E = 210000., nu = 0.3, G = E / ( 2 * ( 1 + nu ) );
  const std::vector< double > linearElastic = { E, nu, 7.85e-9 };

  void check( double value, double expected, double relTol, const std::string& what )
  {
    throwExceptionOnFailure( std::abs( value - expected ) <= relTol * std::abs( expected ),
                             MakeString() << what << ": " << value << " != " << expected );
  }

  void checkTangent( Beam&                  beam,
                     const Eigen::VectorXd& U,
                     const Eigen::VectorXd& UOld,
                     double                 h,
                     double                 relTol,
                     const std::string&     what )
  {
    const auto [P, K]  = beam.kernels( U, UOld, false );
    const auto   numK  = beam.numericalTangent( U, UOld, h );
    const double scale = numK.cwiseAbs().maxCoeff();
    throwExceptionOnFailure( checkIfEqual< double >( K, numK, relTol * scale ),
                             MakeString() << what << ": tangent does not match the numerical one\nK=\n"
                                          << K << "\nnumK=\n"
                                          << numK );
  }

} // namespace

/// the beams are registered, carry displacement and rotation, and report a bar2 shape
void testFactoryNames()
{
  for ( const auto& [name, nDim, nDofPerNode] :
        std::vector< std::tuple< std::string, int, int > >{ { "B23", 2, 3 }, { "B33", 3, 6 } } ) {
    std::unique_ptr< MarmotElement > el( MarmotLibrary::MarmotElementFactory::createElement( name, 1 ) );
    throwExceptionOnFailure( el->getNSpatialDimensions() == nDim && el->getNDofPerElement() == 2 * nDofPerNode &&
                               el->getElementShape() == "bar2" && el->getNNodes() == 2,
                             name + " element description" );
    const auto fields = el->getNodeFields();
    throwExceptionOnFailure( fields.size() == 2 &&
                               fields[0] == std::vector< std::string >{ "displacement", "rotation" },
                             name + " node fields" );
  }
}

/// 2D cantilever (inclined, a single element) under a tip force, a tip moment and an axial force: exact for every
/// profile; the moment at the quadrature points is the exact P (L - x)
void testCantilever2D()
{
  const double          b = 10, h = 20, A = b * h, I = b * h * h * h / 12, L = 500;
  const Eigen::Vector2d t( 0.6, 0.8 ), n( -0.8, 0.6 );
  for ( int profile : { 0, 1, 2 } ) {
    Beam              beam( "B23",
                            { 100., 50., 100. + 300., 50. + 400. },
                            { A, I, double( profile ), 6. },
               "LINEARELASTIC",
               linearElastic );
    const std::string what = "B23 profile " + std::to_string( profile ) + ": ";

    const double    P = 7.5, M = 1200., N = 3000.;
    Eigen::VectorXd f( 3 );
    f << P * n, 0;
    Eigen::VectorXd u = beam.cantilever( f );
    check( u.head( 2 ).dot( n ), P * L * L * L / ( 3 * E * I ), 1e-10, what + "tip deflection (force)" );
    check( u( 2 ), P * L * L / ( 2 * E * I ), 1e-10, what + "tip rotation (force)" );
    throwExceptionOnFailure( std::abs( u.head( 2 ).dot( t ) ) < 1e-12, what + "no axial displacement (force)" );

    Eigen::VectorXd U = Eigen::VectorXd::Zero( 6 );
    U.tail( 3 )       = u;
    beam.kernels( U, Eigen::VectorXd::Zero( 6 ), true );
    const auto xs = beam.element->getCoordinatesAtQuadraturePoints();
    for ( int qp = 0; qp < 3; qp++ ) {
      const double x = ( Eigen::Vector2d( xs[qp][0], xs[qp][1] ) - Eigen::Vector2d( 100., 50. ) ).norm();
      check( beam.state( "bending moment", qp )( 0 ), P * ( L - x ), 1e-10, what + "bending moment" );
      check( beam.state( "curvature", qp )( 0 ), P * ( L - x ) / ( E * I ), 1e-10, what + "curvature" );
    }
    {
      // the linear elastic material does not track its energy; an elastic von Mises material does
      Beam vm( "B23",
               { 100., 50., 100. + 300., 50. + 400. },
               { A, I, double( profile ), 6. },
               "VONMISES",
               { E, nu, 1e12, 0., 0., 1., 7.85e-9 } );
      vm.kernels( U, Eigen::VectorXd::Zero( 6 ), true );
      double energy = 0;
      vm.element->computeInternalEnergy( energy );
      check( energy, 0.5 * P * u.head( 2 ).dot( n ), 1e-10, what + "elastic energy" );
    }

    f << 0, 0, M;
    u = beam.cantilever( f );
    check( u.head( 2 ).dot( n ), M * L * L / ( 2 * E * I ), 1e-10, what + "tip deflection (moment)" );
    check( u( 2 ), M * L / ( E * I ), 1e-10, what + "tip rotation (moment)" );

    f << N * t, 0;
    u = beam.cantilever( f );
    check( u.head( 2 ).dot( t ), N * L / ( E * A ), 1e-10, what + "axial elongation" );
    throwExceptionOnFailure( std::abs( u( 2 ) ) < 1e-14, what + "no rotation (axial)" );
  }
}

/// 3D cantilever in a general orientation: bending about both local axes and torsion
void testCantilever3D()
{
  const double          A = 200, Iy = 1000, Iz = 4000, J = 2500, L = 300;
  const Eigen::Vector3d X1( 10., -20., 5. ), t = Eigen::Vector3d( 1., 2., 2. ) / 3.;
  const Eigen::Vector3d X2 = X1 + L * t;
  const Eigen::Vector3d v( 0., 0., 1. );
  const Eigen::Vector3d e2 = ( v - v.dot( t ) * t ).normalized(), e3 = t.cross( e2 );

  for ( int profile : { 0, 1, 2 } ) {
    Beam              beam( "B33",
                            { X1( 0 ), X1( 1 ), X1( 2 ), X2( 0 ), X2( 1 ), X2( 2 ) },
                            { A, Iy, Iz, J, v( 0 ), v( 1 ), v( 2 ), double( profile ), 4. },
               "LINEARELASTIC",
               linearElastic );
    const std::string what = "B33 profile " + std::to_string( profile ) + ": ";
    const double      P = 5., Tq = 2000.;

    Eigen::VectorXd f( 6 ), u;
    f << P * e2, Eigen::Vector3d::Zero();
    u = beam.cantilever( f );
    check( u.head( 3 ).dot( e2 ), P * L * L * L / ( 3 * E * Iz ), 1e-10, what + "deflection along e2" );
    check( u.tail( 3 ).dot( e3 ), P * L * L / ( 2 * E * Iz ), 1e-10, what + "rotation about e3" );
    throwExceptionOnFailure( std::abs( u.head( 3 ).dot( e3 ) ) + std::abs( u.tail( 3 ).dot( e2 ) ) < 1e-12,
                             what + "no coupling e2 -> e3" );

    f << P * e3, Eigen::Vector3d::Zero();
    u = beam.cantilever( f );
    check( u.head( 3 ).dot( e3 ), P * L * L * L / ( 3 * E * Iy ), 1e-10, what + "deflection along e3" );
    check( u.tail( 3 ).dot( e2 ), -P * L * L / ( 2 * E * Iy ), 1e-10, what + "rotation about e2" );

    f << Eigen::Vector3d::Zero(), Tq * t;
    u = beam.cantilever( f );
    check( u.tail( 3 ).dot( t ), Tq * L / ( G * J ), 1e-10, what + "twist" );
    throwExceptionOnFailure( u.head( 3 ).norm() < 1e-12, what + "no displacement under torsion" );

    Eigen::VectorXd U = Eigen::VectorXd::Zero( 12 );
    U.tail( 6 )       = u;
    beam.kernels( U, Eigen::VectorXd::Zero( 12 ), true );
    check( beam.state( "torque", 1 )( 0 ), Tq, 1e-10, what + "torque" );
    check( beam.state( "twist rate", 1 )( 0 ), Tq / ( G * J ), 1e-10, what + "twist rate" );
  }
}

/// (infinitesimal) rigid body motions are stress free
void testRigidBodyMotions()
{
  {
    Beam            beam( "B23", { 1., 2., 41., 32. }, { 20., 50. }, "LINEARELASTIC", linearElastic );
    const double    w = 1e-3;
    Eigen::VectorXd U( 6 );
    U << 0.3 - w * 2., 0.1 + w * 1., w, 0.3 - w * 32., 0.1 + w * 41., w;
    const auto P = beam.kernels( U, Eigen::VectorXd::Zero( 6 ), false ).first;
    throwExceptionOnFailure( P.norm() < 1e-10, MakeString() << "B23 rigid body motion: |P| = " << P.norm() );
  }
  {
    Beam                  beam( "B33",
                                { 1., 2., 3., 41., 32., -7. },
                                { 20., 50., 60., 80., 0., 0., 1. },
               "LINEARELASTIC",
               linearElastic );
    const Eigen::Vector3d c( 0.3, -0.2, 0.1 ), w( 1e-3, -2e-3, 0.5e-3 );
    Eigen::VectorXd       U( 12 );
    U << c + w.cross( Eigen::Vector3d( 1., 2., 3. ) ), w, c + w.cross( Eigen::Vector3d( 41., 32., -7. ) ), w;
    const auto P = beam.kernels( U, Eigen::VectorXd::Zero( 12 ), false ).first;
    throwExceptionOnFailure( P.norm() < 1e-10, MakeString() << "B33 rigid body motion: |P| = " << P.norm() );
  }
}

/// consistent nodal loads of a uniform line load: total force b A L and fixed-end moments q L^2 / 12
void testBodyForce()
{
  const double    A = 20, L = 50;
  Beam            beam( "B23", { 0., 0., 30., 40. }, { A, 50. }, "LINEARELASTIC", linearElastic );
  Eigen::VectorXd P    = Eigen::VectorXd::Zero( 6 );
  const double    b[2] = { 0., -2. };
  beam.element->computeBodyForce( P.data(), nullptr, b, nullptr, 0, 1 );
  check( P( 1 ) + P( 4 ), -2. * A * L, 1e-12, "body force: total" );
  const double q = -2. * A * 0.6; // transverse component (normal (-0.8, 0.6))
  check( P( 2 ), q * L * L / 12, 1e-12, "body force: moment at node 1" );
  check( P( 5 ), -q * L * L / 12, 1e-12, "body force: moment at node 2" );
}

/// fiber section with von Mises plasticity: consistent tangents in plastic states (bending, axial force, torsion)
void testPlasticTangents()
{
  //                                     E        nu   fy    H     dfy  delta rho
  const std::vector< double > vm = { 210000., 0.3, 500., 2000., 0., 1., 7.85e-9 };
  {
    Beam            beam( "B23", { 0., 0., 60., 80. }, { 200., 6666.67, 1., 6. }, "VONMISES", vm );
    Eigen::VectorXd UOld = Eigen::VectorXd::Zero( 6 ), U( 6 );
    U << 0., 0., 0.02, 0.15, 0.05, 0.06;
    beam.kernels( U, UOld, true );
    UOld = U;
    U << 0., 0., 0.025, 0.17, 0.06, 0.07;
    checkTangent( beam, U, UOld, 1e-8, 1e-5, "B23 von Mises" );
  }
  {
    Beam            beam( "B33",
                          { 0., 0., 0., 60., 80., 0. },
                          { 200., 1666.67, 6666.67, 4500., 0., 0., 1., 1., 4. },
               "VONMISES",
               vm );
    Eigen::VectorXd UOld = Eigen::VectorXd::Zero( 12 ), U( 12 );
    U << 0., 0., 0., 0., 0., 0., 0.15, 0.05, 0.1, 0.05, -0.02, 0.06;
    beam.kernels( U, UOld, true );
    UOld = U;
    U << 0., 0., 0., 0., 0., 0., 0.17, 0.06, 0.12, 0.06, -0.025, 0.07;
    checkTangent( beam, U, UOld, 1e-8, 1e-5, "B33 von Mises (bending, axial, torsion)" );
    throwExceptionOnFailure( beam.state( "dissipation", 0 )( 0 ) > 0, "B33 von Mises: plastic dissipation" );
  }
}

/// fully plastic rectangle (perfect von Mises plasticity): the plastic moment of the fibers, and the axial yield force
void testPlasticLimits()
{
  const double                b = 10, h = 20, A = b * h, I = b * h * h * h / 12, L = 100, fy = 500;
  const int                   n  = 20;
  const std::vector< double > vm = { E, nu, fy, 0., 0., 1., 7.85e-9 };

  Beam beam( "B23", { 0., 0., L, 0. }, { A, I, 1., double( n ) }, "VONMISES", vm );
  // pure bending: kappa = 2 Theta / L, up to 50 x the first yield curvature 2 fy / (E h)
  const double    kappaMax = 50 * 2 * fy / ( E * h );
  Eigen::VectorXd UOld     = Eigen::VectorXd::Zero( 6 );
  for ( int i = 1; i <= 50; i++ ) {
    const double    Theta = 0.5 * i / 50. * kappaMax * L;
    Eigen::VectorXd U( 6 );
    U << 0, 0, -Theta, 0, 0, Theta;
    beam.kernels( U, UOld, true );
    UOld = U;
  }
  // the fiber coordinates are scaled by n / sqrt(n^2 - 1) to reproduce I exactly
  const double Mp = fy * b * h * h / 4 * n / std::sqrt( n * n - 1. );
  for ( int qp = 0; qp < 3; qp++ ) {
    check( beam.state( "bending moment", qp )( 0 ), Mp, 1e-6, "plastic moment" );
    throwExceptionOnFailure( std::abs( beam.state( "normal force", qp )( 0 ) ) < 1e-6 * fy * A, "pure bending" );
  }

  Beam bar( "B23", { 0., 0., L, 0. }, { A, I }, "VONMISES", vm );
  UOld = Eigen::VectorXd::Zero( 6 );
  for ( int i = 1; i <= 10; i++ ) {
    Eigen::VectorXd U = Eigen::VectorXd::Zero( 6 );
    U( 3 )            = i * 0.05;
    bar.kernels( U, UOld, true );
    UOld = U;
  }
  check( bar.state( "normal force", 1 )( 0 ), fy * A, 1e-8, "axial yield force" );
}

/// invalid input is rejected
void testInvalidInput()
{
  bool thrown = false;
  try {
    Beam beam( "B33", { 0., 0., 0., 0., 0., 10. }, { 20., 50., 60., 80., 0., 0., 1. }, "LINEARELASTIC", linearElastic );
  }
  catch ( const std::invalid_argument& ) {
    thrown = true;
  }
  throwExceptionOnFailure( thrown, "orientation vector parallel to the axis must be rejected" );
}

int main()
{
  executeTestsAndCollectExceptions( {
    testFactoryNames,
    testCantilever2D,
    testCantilever3D,
    testRigidBodyMotions,
    testBodyForce,
    testPlasticTangents,
    testPlasticLimits,
    testInvalidInput,
  } );
  return 0;
}

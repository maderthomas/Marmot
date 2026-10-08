#include "Marmot/MarmotMeshfreeKernelFunctionBSpline1stOrderBoxed.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <cmath>

namespace Marmot::Meshfree {

  MarmotMeshfreeKernelFunctionBSpline1stOrderBoxed::MarmotMeshfreeKernelFunctionBSpline1stOrderBoxed(
    double* centerCoord,
    int     dim,
    double  supportRadius )
    : _centerCoord( centerCoord ), _supportRadius( supportRadius ), _dim( dim )
  {
  }

  double MarmotMeshfreeKernelFunctionBSpline1stOrderBoxed::computeKernelFunction( const double* coord ) const
  {

    double res = 1.0;

    for ( int i = 0; i < _dim; i++ ) {
      res *= computeBSpline1stOrder( coord[i] - _centerCoord[i] );
    }

    return res;
  }

  void MarmotMeshfreeKernelFunctionBSpline1stOrderBoxed::computeKernelFunctionGradient( const double* coord,
                                                                                        double*       grad ) const
  {
    for ( int i = 0; i < _dim; i++ ) {
      grad[i] = 0;
    }

    for ( int i = 0; i < _dim; i++ ) {
      double res = 1.0;
      for ( int j = 0; j < _dim; j++ ) {
        if ( i == j )
          res *= computeBSpline1stOrderGradient( coord[j] - _centerCoord[j] );
        else
          res *= computeBSpline1stOrder( coord[j] - _centerCoord[j] );
      }
      grad[i] = res;
    }
  }

  double MarmotMeshfreeKernelFunctionBSpline1stOrderBoxed::computeBSpline1stOrder( double coord_minus_center ) const
  {
    const double z = std::abs( coord_minus_center ) / _supportRadius;
    if ( z < 1.0 ) {
      return 1.0 - z;
    }

    return 0;
  }

  double MarmotMeshfreeKernelFunctionBSpline1stOrderBoxed::computeBSpline1stOrderGradient(
    double coord_minus_center ) const
  {
    const double z         = std::abs( coord_minus_center ) / _supportRadius;
    const double dz_dcoord = coord_minus_center > 0 ? 1.0 / _supportRadius : -1.0 / _supportRadius;

    return ( z < 1.0 && z > 0.0 ) ? dz_dcoord * -1.0 : 0;
  }

  const double* MarmotMeshfreeKernelFunctionBSpline1stOrderBoxed::getCenterCoordinates() const
  {
    return _centerCoord;
  }

  void MarmotMeshfreeKernelFunctionBSpline1stOrderBoxed::moveTo( const double* coordinate )
  {
    for ( int i = 0; i < _dim; i++ ) {
      _centerCoord[i] = coordinate[i];
    }
  }

  bool MarmotMeshfreeKernelFunctionBSpline1stOrderBoxed::isInSupport( const double* coord ) const
  {
    return computeKernelFunction( coord ) > 0;
  }

  void MarmotMeshfreeKernelFunctionBSpline1stOrderBoxed::getBoundingBox( double* min, double* max ) const
  {
    for ( int i = 0; i < _dim; i++ ) {
      min[i] = _centerCoord[i] - _supportRadius;
      max[i] = _centerCoord[i] + _supportRadius;
    }
  }

}; // namespace Marmot::Meshfree

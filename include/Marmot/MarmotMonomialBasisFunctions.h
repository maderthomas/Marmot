/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck,
 * 2020 - today
 *
 * festigkeitslehre@uibk.ac.at
 *
 * Matthias Neuner matthias.neuner@uibk.ac.at
 *
 * This file is part of the MAteRialMOdellingToolbox (marmot).
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of marmot.
 * ---------------------------------------------------------------------
 */

#pragma once

#include <Eigen/Core>
#include <cmath>

namespace Marmot::Math {

  inline int computeSizeOfMonomialBasisVector( int order, int dim )
  {
    // compute the size of the H vector based on the completeness order and the dimension
    // for Eq. (3.67) in the book by Belytschko, Chen, Hillman.

    int size = 0;
    for ( int i = 0; i <= order; i++ ) {
      if ( dim == 1 )
        size += 1;
      else
        size += computeSizeOfMonomialBasisVector( order - i, dim - 1 );
    }
    return size;
  }

  inline int _computeMonomialBasisRecursion( int                    order,
                                             const Eigen::VectorXd& x,
                                             Eigen::VectorXd&       res,
                                             int                    idxEnd,
                                             int                    dim )
  {
    for ( int i = 0; i <= order; i++ ) {

      const int idxStart = idxEnd;
      if ( dim > 1 )
        idxEnd = _computeMonomialBasisRecursion( order - i, x, res, idxEnd, dim - 1 );
      else {
        idxEnd++;
      }

      for ( int idx = idxStart; idx < idxEnd; idx++ )
        res( idx ) *= std::pow( x[dim - 1], i );
    }
    return idxEnd;
  }

  inline int _computeMonomialBasisGradientRecursion( int                    order,
                                                     const Eigen::VectorXd& x,
                                                     Eigen::MatrixXd&       res,
                                                     int                    idxEnd,
                                                     int                    dim )
  {
    // column j of res is grad_j of the monomial: the derivative factor in dimension j
    // times the PLAIN powers of all other dimensions. Each recursion level multiplies
    // every column by its factor (res starts at ones), so mixed monomials like x*y get
    // grad = (y, x), not the bare derivative factors (1, 1).
    const int nDimTotal = static_cast< int >( x.size() );
    for ( int i = 0; i <= order; i++ ) {

      const int idxStart = idxEnd;

      if ( dim > 1 )
        idxEnd = _computeMonomialBasisGradientRecursion( order - i, x, res, idxEnd, dim - 1 );
      else {
        idxEnd++;
      }
      const double f  = std::pow( x[dim - 1], i );
      const double df = i > 0 ? i * std::pow( x[dim - 1], i - 1 ) : 0.0;
      for ( int idx = idxStart; idx < idxEnd; idx++ )
        for ( int j = 0; j < nDimTotal; j++ )
          res( idx, j ) *= ( j == dim - 1 ) ? df : f;
    }
    return idxEnd;
  }

  inline void computeMonomialBasis( int order, const Eigen::VectorXd& x, Eigen::VectorXd& res )
  {
    res.setOnes();
    _computeMonomialBasisRecursion( order, x, res, 0, x.size() );
  }

  inline void computeMonomialBasisGradient( int order, const Eigen::VectorXd& x, Eigen::MatrixXd& res )
  {
    res.setOnes();
    _computeMonomialBasisGradientRecursion( order, x, res, 0, x.size() );
  }

} // namespace Marmot::Math

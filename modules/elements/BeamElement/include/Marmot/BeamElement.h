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
#include "Marmot/MarmotElement.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/MarmotStateVarVectorManager.h"
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <numbers>
#include <string>
#include <vector>

namespace Marmot::Elements {

  /**
   * @class Marmot::Elements::BeamElement
   * @brief Two-node Euler-Bernoulli beam (frame) element in two or three spatial dimensions, geometrically linear,
   * with a fiber section for any hypoelastic material of Marmot.
   *
   * @tparam nDim number of spatial dimensions: 2 (nodal dofs \f$u_x, u_y, \theta_z\f$) or 3 (nodal dofs
   *              \f$u_x, u_y, u_z, \theta_x, \theta_y, \theta_z\f$)
   *
   * Node fields: @b displacement and @b rotation (2D: one rotation about the out-of-plane axis; 3D: the rotation
   * vector). The dofs are ordered node by node, displacement before rotation.
   *
   * @par Kinematics
   * In the local frame \f$( \boldsymbol{e}_1, \boldsymbol{e}_2, \boldsymbol{e}_3 )\f$ (\f$\boldsymbol{e}_1\f$ along
   * the axis from node 1 to node 2, local coordinates \f$( x, y, z )\f$ with \f$y, z\f$ in the cross section) the
   * axial displacement \f$u\f$ and the twist \f$\varphi\f$ are interpolated linearly, the transverse displacements
   * \f$v, w\f$ by cubic Hermite polynomials together with the nodal rotations \f$\theta_z = v'\f$ and
   * \f$\theta_y = - w'\f$ (Euler-Bernoulli: no shear deformation, hence no shear locking). The generalized strains are
   * the axial strain \f$\varepsilon_0 = u'\f$, the curvatures \f$\kappa_y = \theta_y' = -w''\f$,
   * \f$\kappa_z = \theta_z' = v''\f$ and the rate of twist \f$\chi = \varphi'\f$; in 2D only \f$\varepsilon_0\f$ and
   * \f$\kappa = \kappa_z\f$. The element reproduces the exact nodal displacements and rotations of a linear elastic
   * beam loaded at its nodes (e.g., a cantilever under a tip force or moment with a single element).
   *
   * Geometrically linear: rotations must remain small (a rigid body rotation of finite size strains the element).
   *
   * @par Fiber section
   * The section is integrated with fibers at \f$( y_f, z_f )\f$ with areas \f$A_f\f$. A fiber has the axial strain
   * \f$\varepsilon_{11} = \varepsilon_0 + z_f \kappa_y - y_f \kappa_z\f$ and, in 3D, the torsional shear strains
   * \f$\gamma_{12} = - s z_f \chi\f$, \f$\gamma_{13} = s y_f \chi\f$; its material is evaluated in the stress state
   * \f$\sigma_{22} = \sigma_{33} = \sigma_{23} = 0\f$ (Newton iteration on the free strains, condensed tangent). The
   * section forces \f$( N, M_y, M_z, T )\f$ are the work conjugates of \f$( \varepsilon_0, \kappa_y, \kappa_z, \chi
   * )\f$. The fiber coordinates are scaled such that the fibers reproduce the given area and second moments of area
   * exactly, and \f$s = \sqrt{ J / \sum_f ( y_f^2 + z_f^2 ) A_f }\f$ such that the torsional stiffness is exactly
   * \f$G J\f$ for an elastic material. Hence, for a linear elastic material, \f$N = E A \varepsilon_0\f$,
   * \f$M_y = E I_y \kappa_y\f$, \f$M_z = E I_z \kappa_z\f$, \f$T = G J \chi\f$ exactly, for any profile. Profiles:
   *
   * - 0, @b generic: the fewest fibers reproducing \f$A, I\f$: 2 fibers at \f$y = \pm \sqrt{I/A}\f$ (2D) or 4 fibers
   *   at \f$( \pm \sqrt{I_z/A}, \pm \sqrt{I_y/A} )\f$ (3D). Exact for elastic materials; for inelastic materials it
   *   is an idealized sandwich section (full plastification at first yield).
   * - 1, @b rectangle: a grid of \f$n\f$ (2D: layers through the height) or \f$n \times n\f$ (3D) fibers.
   * - 2, @b circle: \f$\max( 1, n/2 )\f$ rings of equal area, each with \f$2n\f$ (rounded up to a multiple of 4)
   *   sectors.
   *
   * Element properties:
   * - 2D: \f$[ A, I, \text{profile}, n ]\f$ with \f$I = \int y^2 \mathrm{d}A\f$ (bending in the plane).
   * - 3D: \f$[ A, I_y, I_z, J, v_x, v_y, v_z, \text{profile}, n ]\f$ with \f$I_y = \int z^2 \mathrm{d}A\f$,
   *   \f$I_z = \int y^2 \mathrm{d}A\f$, the torsion constant \f$J\f$, and the orientation vector \f$\boldsymbol{v}\f$
   *   whose part normal to the axis is the local \f$y\f$ axis \f$\boldsymbol{e}_2\f$; \f$\boldsymbol{e}_3 =
   *   \boldsymbol{e}_1 \times \boldsymbol{e}_2\f$.
   *
   * profile and \f$n\f$ are optional (defaults: generic, \f$n = 8\f$). The section integration along the axis uses 3
   * Gauss points (exact for an elastic beam). The material of all fibers is the same; it is a hypoelastic material
   * (MarmotMaterialHypoElastic) of Marmot, which sees the characteristic length \f$L/3\f$.
   *
   * Quadrature point states: 2D: @b normal force, @b bending moment, @b axial strain, @b curvature; 3D: @b normal
   * force, @b bending moments (\f$M_y, M_z\f$), @b torque, @b axial strain, @b curvatures (\f$\kappa_y, \kappa_z\f$),
   * @b twist rate; both: @b section forces and @b section strains (all generalized forces or strains), @b elastic
   * energy, @b dissipation (integrated over the quadrature point's volume). The state of fiber \f$k\f$ (stress
   * \f$\sigma_{11}, \sigma_{12}, \sigma_{13}\f$, and the material states) is accessible as "fiber k stress" or "fiber
   * k <material state name>".
   *
   * Not implemented: distributed (surface) loads, inertia.
   */
  template < int nDim >
  class BeamElement : public MarmotElement {

    static_assert( nDim == 2 || nDim == 3, "BeamElement: nDim must be 2 or 3" );

  public:
    /// The cross section profile, which determines the fiber layout
    enum Profile {
      Generic   = 0,
      Rectangle = 1,
      Circle    = 2,
    };

    static constexpr int nNodes          = 2;
    static constexpr int nRotations      = nDim == 2 ? 1 : 3;
    static constexpr int nDofPerNode     = nDim + nRotations;
    static constexpr int sizeLoadVector  = nNodes * nDofPerNode;
    static constexpr int nGeneralized    = nDim == 2 ? 2 : 4; ///< generalized strains / section forces
    static constexpr int nProperties     = nDim == 2 ? 2 : 7; ///< required element properties
    static constexpr int nFiberStateVars = 5;                 ///< fiber stress (3), energy and dissipation density

    using VectorDim  = Eigen::Matrix< double, nDim, 1 >;
    using MatrixDim  = Eigen::Matrix< double, nDim, nDim >;
    using RhsSized   = Eigen::Matrix< double, sizeLoadVector, 1 >;
    using KSized     = Eigen::Matrix< double, sizeLoadVector, sizeLoadVector >;
    using BSized     = Eigen::Matrix< double, nGeneralized, sizeLoadVector >;
    using VectorGen  = Eigen::Matrix< double, nGeneralized, 1 >;
    using MatrixGen  = Eigen::Matrix< double, nGeneralized, nGeneralized >;
    using FiberS     = Eigen::Matrix< double, 3, nGeneralized >; ///< generalized strains -> fiber strains
    using InterpSize = Eigen::Matrix< double, nDim, sizeLoadVector >;

    /// A fiber of the cross section
    struct Fiber {
      double y;    ///< local y coordinate
      double z;    ///< local z coordinate (0 in 2D)
      double area; ///< area
      FiberS S;    ///< generalized strains -> fiber strains [eps11, gamma12, gamma13]
    };

    /// A quadrature point along the axis
    struct QuadraturePoint {
      const double xi;        ///< parent coordinate in [-1, 1]
      const double weight;    ///< Gauss weight
      BSized       B;         ///< generalized strains of the global dofs
      double       dV = 0.0;  ///< length x weight / 2
      double*      stateVars; ///< the state of this point
      QuadraturePoint( double xi, double weight ) : xi( xi ), weight( weight ), stateVars( nullptr ) {}
    };

    /// Element label
    const int elLabel;

    /// Element properties, see the class documentation
    Eigen::Map< const Eigen::VectorXd > elementProperties;
    /// Reference nodal coordinates, [dim, node]
    Eigen::Map< const Eigen::Matrix< double, nDim, nNodes > > coordinates;

    /// Length of the element
    double length = 0.0;
    /// Rows: the local axes e1, e2(, e3) in global coordinates
    MatrixDim R = MatrixDim::Identity();
    /// Global -> local dof transformation
    KSized T = KSized::Identity();
    /// The fibers of the cross section
    std::vector< Fiber > fibers;
    /// The torsion scaling factor s of the fiber shear strains
    double torsionScaling = 1.0;
    /// The material, shared by all fibers (a hypoelastic material carries no state itself)
    std::unique_ptr< MarmotMaterialHypoElastic > material;
    /// The quadrature points along the axis
    std::vector< QuadraturePoint > qps;

    /// The names of the per-point states and their lengths, see the class documentation
    static const std::vector< std::pair< std::string, int > >& qpStateEntries()
    {
      static const std::vector< std::pair< std::string, int > >
        entries = nDim == 2 ? std::vector< std::pair< std::string, int > >{ { "normal force", 1 },
                                                                            { "bending moment", 1 },
                                                                            { "axial strain", 1 },
                                                                            { "curvature", 1 },
                                                                            { "elastic energy", 1 },
                                                                            { "dissipation", 1 } }
                            : std::vector< std::pair< std::string, int > >{ { "normal force", 1 },
                                                                            { "bending moments", 2 },
                                                                            { "torque", 1 },
                                                                            { "axial strain", 1 },
                                                                            { "curvatures", 2 },
                                                                            { "twist rate", 1 },
                                                                            { "elastic energy", 1 },
                                                                            { "dissipation", 1 } };
      return entries;
    }

    /// offsets in the state of a point
    static constexpr int idxSectionForces  = 0;
    static constexpr int idxSectionStrains = nGeneralized;
    static constexpr int idxElasticEnergy  = 2 * nGeneralized;
    static constexpr int idxDissipation    = 2 * nGeneralized + 1;
    static constexpr int nQpStateVarsOwn   = 2 * nGeneralized + 2;

    explicit BeamElement( int elementID )
      : elLabel( elementID ), elementProperties( nullptr, 0 ), coordinates( nullptr )
    {
      for ( const auto& qpInfo : FiniteElement::Quadrature::Spatial1D::gaussPointList3 )
        qps.emplace_back( qpInfo.xi( 0 ), qpInfo.weight );
    }

    int nMaterialStateVars() const { return material ? material->getNumberOfRequiredStateVars() : 0; }

    int nFiberStride() const { return nFiberStateVars + nMaterialStateVars(); }

    int nQpStateVars() const { return nQpStateVarsOwn + static_cast< int >( fibers.size() ) * nFiberStride(); }

    int getNumberOfRequiredStateVars() override { return static_cast< int >( qps.size() ) * nQpStateVars(); }

    std::vector< std::vector< std::string > > getNodeFields() override
    {
      return std::vector< std::vector< std::string > >( nNodes, { "displacement", "rotation" } );
    }

    std::vector< int > getDofIndicesPermutationPattern() override
    {
      std::vector< int > pattern( sizeLoadVector );
      for ( int i = 0; i < sizeLoadVector; i++ )
        pattern[i] = i;
      return pattern;
    }

    int getNNodes() override { return nNodes; }

    int getNSpatialDimensions() override { return nDim; }

    int getNDofPerElement() override { return sizeLoadVector; }

    std::string getElementShape() override { return "bar2"; }

    void assignStateVars( double* stateVars, int nStateVars ) override
    {
      if ( nStateVars < getNumberOfRequiredStateVars() )
        throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": too few state variables" );
      for ( size_t i = 0; i < qps.size(); i++ )
        qps[i].stateVars = stateVars + i * nQpStateVars();
    }

    void assignProperty( const ElementProperties& property ) override
    {
      if ( property.nElementProperties < nProperties )
        throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": a beam element in " << nDim
                                                  << "D requires " << nProperties << " section properties" );
      new ( &elementProperties )
        Eigen::Map< const Eigen::VectorXd >( property.elementProperties, property.nElementProperties );
    }

    void assignProperty( const MarmotMaterialSection& section ) override
    {
      material = std::unique_ptr< MarmotMaterialHypoElastic >(
        MarmotLibrary::MarmotMaterialHypoElasticFactory::createMaterial( section.materialName,
                                                                         section.materialProperties,
                                                                         section.nMaterialProperties,
                                                                         elLabel ) );
      material->setCharacteristicElementLength( length / static_cast< double >( qps.size() ) );
    }

    void assignNodeCoordinates( const double* coords ) override
    {
      new ( &coordinates ) Eigen::Map< const Eigen::Matrix< double, nDim, nNodes > >( coords );
    }

    /// The section properties [A, Iy, Iz, J] (2D: Iy = J = 0, Iz = I)
    std::array< double, 4 > sectionProperties() const
    {
      if constexpr ( nDim == 2 )
        return { elementProperties[0], 0.0, elementProperties[1], 0.0 };
      else
        return { elementProperties[0], elementProperties[1], elementProperties[2], elementProperties[3] };
    }

    int profile() const
    {
      return elementProperties.size() > nProperties ? static_cast< int >( elementProperties[nProperties] ) : Generic;
    }

    int nFibersPerDirection() const
    {
      return elementProperties.size() > nProperties + 1 ? static_cast< int >( elementProperties[nProperties + 1] ) : 8;
    }

    /// The fibers of the profile, scaled to the area and the second moments of area
    void makeFibers()
    {
      const auto [A, Iy, Iz, J] = sectionProperties();
      if ( A <= 0 || Iz <= 0 || ( nDim == 3 && ( Iy <= 0 || J <= 0 ) ) )
        throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": element " << elLabel
                                                  << ": the section properties must be positive" );

      std::vector< std::array< double, 3 > > raw; // y, z, relative area
      const int                              n = nFibersPerDirection();
      switch ( profile() ) {
      case Generic: {
        if constexpr ( nDim == 2 )
          raw = { { -1, 0, 1 }, { 1, 0, 1 } };
        else
          raw = { { -1, -1, 1 }, { 1, -1, 1 }, { 1, 1, 1 }, { -1, 1, 1 } };
        break;
      }
      case Rectangle: {
        if ( n < 2 )
          throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": a rectangle needs n >= 2 fibers" );
        for ( int i = 0; i < n; i++ )
          for ( int j = 0; j < ( nDim == 2 ? 1 : n ); j++ )
            raw.push_back( { ( i + 0.5 ) / n - 0.5, nDim == 2 ? 0.0 : ( j + 0.5 ) / n - 0.5, 1.0 } );
        break;
      }
      case Circle: {
        if ( n < 2 )
          throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": a circle needs n >= 2" );
        const int nRings   = std::max( 1, n / 2 );
        const int nSectors = 4 * ( ( 2 * n + 3 ) / 4 );
        for ( int r = 0; r < nRings; r++ ) {
          // rings of equal area, fibers at the ring's centroidal radius
          const double ri = std::sqrt( static_cast< double >( r ) / nRings );
          const double ro = std::sqrt( static_cast< double >( r + 1 ) / nRings );
          const double rc = 2. / 3. * ( ro * ro * ro - ri * ri * ri ) / ( ro * ro - ri * ri );
          for ( int s = 0; s < nSectors; s++ ) {
            const double phi = ( s + 0.5 ) * 2 * std::numbers::pi / nSectors;
            raw.push_back( { rc * std::cos( phi ), nDim == 2 ? 0.0 : rc * std::sin( phi ), 1.0 } );
          }
        }
        break;
      }
      default: throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": unknown profile " << profile() );
      }

      double sumA = 0, sumYY = 0, sumZZ = 0;
      for ( const auto& f : raw ) {
        sumA += f[2];
        sumYY += f[0] * f[0] * f[2];
        sumZZ += f[1] * f[1] * f[2];
      }
      // scale the areas to A and the coordinates to Iz = sum y^2 A_f and Iy = sum z^2 A_f
      const double scaleY = std::sqrt( Iz / A * sumA / sumYY );
      const double scaleZ = nDim == 2 ? 0.0 : std::sqrt( Iy / A * sumA / sumZZ );

      fibers.clear();
      double polar = 0;
      for ( const auto& f : raw ) {
        Fiber fiber{ f[0] * scaleY, f[1] * scaleZ, f[2] / sumA * A, FiberS::Zero() };
        polar += ( fiber.y * fiber.y + fiber.z * fiber.z ) * fiber.area;
        fibers.push_back( fiber );
      }
      torsionScaling = nDim == 2 ? 0.0 : std::sqrt( J / polar );

      for ( auto& f : fibers ) {
        if constexpr ( nDim == 2 )
          f.S << 1, -f.y, //
            0, 0,         //
            0, 0;
        else
          f.S << 1, f.z, -f.y, 0,           //
            0, 0, 0, -torsionScaling * f.z, //
            0, 0, 0, torsionScaling * f.y;
      }
    }

    /// The local displacements [u, v(, w)] at the parent coordinate xi in terms of the global dofs
    InterpSize displacementInterpolation( double xi ) const
    {
      const double s  = 0.5 * ( xi + 1 );
      const double L  = length;
      const double H1 = 1 - 3 * s * s + 2 * s * s * s, H2 = L * ( s - 2 * s * s + s * s * s );
      const double H3 = 3 * s * s - 2 * s * s * s, H4 = L * ( -s * s + s * s * s );

      Eigen::Matrix< double, nDim, sizeLoadVector > N  = Eigen::Matrix< double, nDim, sizeLoadVector >::Zero();
      constexpr int                                 n2 = nDofPerNode;
      N( 0, 0 )                                        = 1 - s;
      N( 0, n2 )                                       = s;
      if constexpr ( nDim == 2 ) {
        N( 1, 1 ) = H1, N( 1, 2 ) = H2, N( 1, n2 + 1 ) = H3, N( 1, n2 + 2 ) = H4;
      }
      else {
        N( 1, 1 ) = H1, N( 1, 5 ) = H2, N( 1, n2 + 1 ) = H3, N( 1, n2 + 5 ) = H4;   // v with theta_z
        N( 2, 2 ) = H1, N( 2, 4 ) = -H2, N( 2, n2 + 2 ) = H3, N( 2, n2 + 4 ) = -H4; // w with theta_y = -w'
      }
      return N * T;
    }

    /// The generalized strains at the parent coordinate xi in terms of the global dofs
    BSized strainOperator( double xi ) const
    {
      const double s = 0.5 * ( xi + 1 );
      const double L = length;
      // second derivatives of the Hermite polynomials w.r.t. x
      const double d1 = ( -6 + 12 * s ) / ( L * L ), d2 = ( -4 + 6 * s ) / L;
      const double d3 = ( 6 - 12 * s ) / ( L * L ), d4 = ( -2 + 6 * s ) / L;

      BSized        B  = BSized::Zero();
      constexpr int n2 = nDofPerNode;
      B( 0, 0 )        = -1 / L;
      B( 0, n2 )       = 1 / L;
      if constexpr ( nDim == 2 ) {
        B( 1, 1 ) = d1, B( 1, 2 ) = d2, B( 1, n2 + 1 ) = d3, B( 1, n2 + 2 ) = d4; // kappa = v''
      }
      else {
        B( 1, 2 ) = -d1, B( 1, 4 ) = d2, B( 1, n2 + 2 ) = -d3, B( 1, n2 + 4 ) = d4; // kappa_y = -w''
        B( 2, 1 ) = d1, B( 2, 5 ) = d2, B( 2, n2 + 1 ) = d3, B( 2, n2 + 5 ) = d4;   // kappa_z = v''
        B( 3, 3 ) = -1 / L, B( 3, n2 + 3 ) = 1 / L;                                 // chi = phi'
      }
      return B * T;
    }

    void initializeYourself() override
    {
      const VectorDim axis = coordinates.col( 1 ) - coordinates.col( 0 );
      length               = axis.norm();
      if ( length <= 0 )
        throw std::invalid_argument( MakeString()
                                     << __PRETTY_FUNCTION__ << ": element " << elLabel << " has zero length" );
      const VectorDim e1 = axis / length;
      if constexpr ( nDim == 2 ) {
        R.row( 0 ) = e1.transpose();
        R.row( 1 ) << -e1( 1 ), e1( 0 );
      }
      else {
        const Eigen::Vector3d v  = elementProperties.template segment< 3 >( 4 );
        const Eigen::Vector3d vn = v - v.dot( e1 ) * e1;
        if ( vn.norm() < 1e-8 * v.norm() || v.norm() == 0 )
          throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": element " << elLabel
                                                    << ": the orientation vector must not be parallel to the axis" );
        const Eigen::Vector3d e2 = vn.normalized();
        R.row( 0 )               = e1.transpose();
        R.row( 1 )               = e2.transpose();
        R.row( 2 )               = e1.cross( e2 ).transpose();
      }

      T.setZero();
      for ( int a = 0; a < nNodes; a++ ) {
        T.template block< nDim, nDim >( a * nDofPerNode, a * nDofPerNode ) = R;
        if constexpr ( nDim == 2 )
          T( a * nDofPerNode + 2, a * nDofPerNode + 2 ) = 1;
        else
          T.template block< 3, 3 >( a * nDofPerNode + 3, a * nDofPerNode + 3 ) = R;
      }

      makeFibers();

      for ( auto& qp : qps ) {
        qp.B  = strainOperator( qp.xi );
        qp.dV = 0.5 * length * qp.weight;
      }
      if ( material )
        material->setCharacteristicElementLength( length / static_cast< double >( qps.size() ) );
    }

    double* fiberState( const QuadraturePoint& qp, size_t f ) const
    {
      return qp.stateVars + nQpStateVarsOwn + f * nFiberStride();
    }

    void setInitialConditions( StateTypes state, const double* ) override
    {
      switch ( state ) {
      case MarmotElement::MarmotMaterialInitialization: {
        for ( auto& qp : qps )
          for ( size_t f = 0; f < fibers.size(); f++ ) {
            double* fs = fiberState( qp, f );
            std::fill( fs, fs + nFiberStateVars, 0.0 );
            material->initializeYourself( fs + nFiberStateVars, nMaterialStateVars() );
          }
        break;
      }
      default: throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": invalid initial condition" );
      }
    }

    /**
     * The stress of a fiber for the strain increment [d eps11, d gamma12, d gamma13], with sigma22 = sigma33 = sigma23
     * = 0 found by a Newton iteration on the free strains, and the condensed tangent.
     *
     * @param fs fiber state: [sigma11, sigma12, sigma13, energy density, dissipation density, material state]
     */
    void computeFiberStress( double* fs, Eigen::Matrix3d& C, const Eigen::Vector3d& dEps, double time, double dT ) const
    {
      using MHE                               = MarmotMaterialHypoElastic;
      static constexpr std::array< int, 3 > p = { 0, 3, 4 }; // prescribed: eps11, gamma12, gamma13
      static constexpr std::array< int, 3 > f = { 1, 2, 5 }; // free: sigma22 = sigma33 = sigma23 = 0

      const int                   nMat = nMaterialStateVars();
      double*                     sv   = fs + nFiberStateVars;
      const std::vector< double > svOld( sv, sv + nMat );
      Marmot::Vector6d            dE = Marmot::Vector6d::Zero();
      Marmot::Matrix6d            C6;
      MHE::state3D                state;
      for ( int i = 0; i < 3; i++ )
        dE( p[i] ) = dEps( i );

      const double sigmaScale = std::max( 1.0, Eigen::Map< const Eigen::Vector3d >( fs ).cwiseAbs().maxCoeff() );
      for ( int iteration = 0;; iteration++ ) {
        std::copy( svOld.begin(), svOld.end(), sv );
        state.stress.setZero();
        for ( int i = 0; i < 3; i++ )
          state.stress( p[i] ) = fs[i];
        state.elasticEnergyDensity = fs[3];
        state.dissipation          = fs[4];
        state.stateVars            = sv;

        C6.setZero();
        material->computeStress( state, C6, dE, { time, dT } );

        Eigen::Vector3d r, sP;
        Eigen::Matrix3d Cff, Cfp, Cpf, Cpp;
        for ( int i = 0; i < 3; i++ ) {
          r( i )  = state.stress( f[i] );
          sP( i ) = state.stress( p[i] );
          for ( int j = 0; j < 3; j++ ) {
            Cff( i, j ) = C6( f[i], f[j] );
            Cfp( i, j ) = C6( f[i], p[j] );
            Cpf( i, j ) = C6( p[i], f[j] );
            Cpp( i, j ) = C6( p[i], p[j] );
          }
        }
        const auto   CffLU    = Cff.fullPivLu();
        const double residual = r.cwiseAbs().maxCoeff();
        if ( residual <= 1e-12 * std::max( sigmaScale, sP.cwiseAbs().maxCoeff() ) ||
             ( iteration > 7 && residual <= 1e-9 * std::max( sigmaScale, sP.cwiseAbs().maxCoeff() ) ) ) {
          for ( int i = 0; i < 3; i++ )
            fs[i] = sP( i );
          fs[3] = state.elasticEnergyDensity;
          fs[4] = state.dissipation;
          C     = Cpp - Cpf * CffLU.solve( Cfp );
          return;
        }
        if ( iteration >= 15 || !CffLU.isInvertible() ) {
          MarmotJournal::warningToMSG( "BeamElement: fiber stress iteration requires cutback" );
          throw Marmot::StressUpdateFailed( MakeString() << __PRETTY_FUNCTION__ << ": element " << elLabel
                                                         << ": the fiber stress iteration did not converge" );
        }
        const Eigen::Vector3d ddE = CffLU.solve( r );
        for ( int i = 0; i < 3; i++ )
          dE( f[i] ) -= ddE( i );
      }
    }

    /// Section forces and tangent for the increment of the generalized strains; updates the state of the point
    void computeSection( QuadraturePoint& qp, VectorGen& s, MatrixGen& D, const VectorGen& dE, double time, double dT )
    {
      s.setZero();
      D.setZero();
      double energy = 0, dissipation = 0;
      for ( size_t k = 0; k < fibers.size(); k++ ) {
        const Fiber&    fiber = fibers[k];
        double*         fs    = fiberState( qp, k );
        Eigen::Matrix3d C;
        computeFiberStress( fs, C, fiber.S * dE, time, dT );
        const Eigen::Map< const Eigen::Vector3d > sigma( fs );
        s += fiber.S.transpose() * sigma * fiber.area;
        D += fiber.S.transpose() * C * fiber.S * fiber.area;
        energy += fs[3] * fiber.area;
        dissipation += fs[4] * fiber.area;
      }
      qp.stateVars[idxElasticEnergy] = energy * qp.dV;
      qp.stateVars[idxDissipation]   = dissipation * qp.dV;
    }

    void computeKernels( const double* QTotal_, const double* dQ_, double* Pe_, double* Ke_, double time, double dT )
      override
    {
      Eigen::Map< const RhsSized > QTotal( QTotal_ );
      Eigen::Map< const RhsSized > dQ( dQ_ );
      Eigen::Map< RhsSized >       Pe( Pe_ );
      Eigen::Map< KSized >         Ke( Ke_ );

      for ( auto& qp : qps ) {
        VectorGen s;
        MatrixGen D;
        computeSection( qp, s, D, qp.B * dQ, time, dT );

        Eigen::Map< VectorGen >( qp.stateVars + idxSectionForces )  = s;
        Eigen::Map< VectorGen >( qp.stateVars + idxSectionStrains ) = qp.B * QTotal;

        Pe += qp.B.transpose() * s * qp.dV;
        Ke += qp.B.transpose() * D * qp.B * qp.dV;
      }
    }

    void computeKernelsExplicit( const double* QTotal, const double* dQ, double* Pe, double time, double dT ) override
    {
      KSized K = KSized::Zero();
      computeKernels( QTotal, dQ, Pe, K.data(), time, dT );
    }

    void computeDistributedLoad( DistributedLoadTypes,
                                 double*,
                                 double*,
                                 int,
                                 const double*,
                                 const double*,
                                 double,
                                 double ) override
    {
      throw std::invalid_argument( MakeString()
                                   << __PRETTY_FUNCTION__ << ": distributed loads are not supported by beam elements" );
    }

    /// Consistent nodal forces and moments of a body force (per volume), i.e., a line load b A along the axis
    void computeBodyForce( double* P_, double*, const double* load, const double*, double, double ) override
    {
      Eigen::Map< RhsSized >        P( P_ );
      Eigen::Map< const VectorDim > b( load );
      const VectorDim               q = R * b * sectionProperties()[0]; // local line load
      for ( const auto& qp : qps )
        P += displacementInterpolation( qp.xi ).transpose() * q * qp.dV;
    }

    void computeInternalEnergy( double& internalEnergy ) override
    {
      for ( const auto& qp : qps )
        internalEnergy += qp.stateVars[idxElasticEnergy];
    }

    StateView getStateView( const std::string& stateName, int qpNumber ) override
    {
      const auto& qp = qps[qpNumber];
      if ( stateName == "section forces" )
        return { qp.stateVars + idxSectionForces, nGeneralized };
      if ( stateName == "section strains" )
        return { qp.stateVars + idxSectionStrains, nGeneralized };
      int offset = 0;
      for ( const auto& [name, len] : qpStateEntries() ) {
        if ( name == stateName )
          return { qp.stateVars + offset, len };
        offset += len;
      }
      // "fiber <k> <state>"
      if ( stateName.rfind( "fiber ", 0 ) == 0 ) {
        const size_t space = stateName.find( ' ', 6 );
        if ( space != std::string::npos ) {
          const size_t k = std::stoul( stateName.substr( 6, space - 6 ) );
          if ( k >= fibers.size() )
            throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": fiber " << k << " does not exist" );
          double*           fs   = fiberState( qp, k );
          const std::string name = stateName.substr( space + 1 );
          if ( name == "stress" )
            return { fs, 3 };
          return material->getStateView( name, fs + nFiberStateVars );
        }
      }
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": unknown state " << stateName );
    }

    std::vector< double > getCoordinatesAtCenter() override
    {
      const VectorDim x = 0.5 * ( coordinates.col( 0 ) + coordinates.col( 1 ) );
      return std::vector< double >( x.data(), x.data() + nDim );
    }

    std::vector< std::vector< double > > getCoordinatesAtQuadraturePoints() override
    {
      std::vector< std::vector< double > > result;
      for ( const auto& qp : qps ) {
        const double    s = 0.5 * ( qp.xi + 1 );
        const VectorDim x = ( 1 - s ) * coordinates.col( 0 ) + s * coordinates.col( 1 );
        result.emplace_back( x.data(), x.data() + nDim );
      }
      return result;
    }

    int getNumberOfQuadraturePoints() override { return static_cast< int >( qps.size() ); }
  };

} // namespace Marmot::Elements

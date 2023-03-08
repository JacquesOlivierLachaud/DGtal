/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
/**
 * @file tangency-explorer.cpp
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5807), University of Savoie, France
 *
 * @date 2022/06/15
 *
 * This file is part of the DGtal library.
 */
#include <iostream>
#include <string>
#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/helpers/Shortcuts.h>
#include <DGtal/helpers/ShortcutsGeometry.h>
#include <DGtal/shapes/SurfaceMesh.h>
#include "DGtal/geometry/volumes/TangencyComputer.h"
#include "DGtal/geometry/tools/QuickHull.h"
#include "DGtal/geometry/surfaces/estimation/PlaneProbingTetrahedronEstimator.h"
#include "SymmetricConvexExpander.h"
#include "SymmetricSeparator.h"

#include "polyscope/pick.h"
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/point_cloud.h>

#include "ConfigExamples.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace DGtal;
using namespace Z3i;

// Using standard 3D digital space.
typedef Shortcuts<Z3i::KSpace>         SH3;
typedef ShortcutsGeometry<Z3i::KSpace> SHG3;
// The following typedefs are useful
typedef SurfaceMesh< RealPoint, RealVector >  SurfMesh;
typedef SurfMesh::Face                        Face;
typedef SurfMesh::Vertex                      Vertex;
typedef DGtal::ConvexHullIntegralKernel< 3 > Kernel3D;
typedef DGtal::QuickHull< Kernel3D >         QuickHull3D;
typedef QuickHull3D::IndexRange              IndexRange;

//Polyscope global
polyscope::SurfaceMesh *psMesh;
polyscope::SurfaceMesh *psDualMesh;
polyscope::SurfaceMesh *psTriangle;
polyscope::SurfaceMesh *psDisks;
polyscope::SurfaceMesh *psSeps;
polyscope::PointCloud*  psCloud;
polyscope::PointCloud*  psCloudCvx;
SurfMesh surfmesh;
SurfMesh dual_surfmesh;
float gridstep   = 1.0;
int   vertex_idx = -1;
int   face_idx   = -1;
int   edge_idx   = -1;
float Time = 0.0;
int   MinNbPoints = 10;
bool  is_selected = false;
Point selected_kpoint; // valid if selection == true
bool  PSym = true;
bool  enforceFC = false;
bool  filterTB = true;
bool  fastSymSep  = false;
bool  onlySurfels = false;
float Thickness = 0.0;
float maxThickness = 2.0;

std::vector< Point >     digital_points;
std::vector< Point >     cell_points; // cell as Khalimsky coordinates.
KSpace K;
DGtal::DigitalConvexity< KSpace > dconv;
DGtal::TangencyComputer< KSpace > TC;
DGtal::LatticeSetByIntervals< Space > LS;
typedef DGtal::TangencyComputer< KSpace >::Index Index;
typedef std::vector< Index >  Indices;
typedef double                Scalar;
typedef std::vector< Scalar > Scalars;

struct UnorderedPointSetPredicate
{
  typedef DGtal::Z3i::Point Point;
  const std::unordered_set< Point >* myS;
  explicit UnorderedPointSetPredicate( const std::unordered_set< Point >& S )
    : myS( &S ) {}
  bool operator()( const Point& p ) const
  { return myS->count( p ) != 0; }
};

std::unordered_set< Point > unorderedSet;
std::unordered_set< Point > immInterior;
std::unordered_set< Point > immExterior;


struct BinaryImagePredicate
{
  typedef DGtal::Z3i::Point Point;
  typedef DGtal::Z3i::Domain Domain;
  typedef DGtal::ImageContainerBySTLVector<Domain, bool> BinaryImage;
  DGtal::CountedPtr<BinaryImage> myImage;
  Domain myDomain;
  explicit BinaryImagePredicate( DGtal::CountedPtr<BinaryImage> PtrI )
    : myImage( PtrI ), myDomain( PtrI->domain() ) {}
  bool operator()( const Point& p ) const
  { return myDomain.isInside( p ) && (*myImage)( p ); }
};

BinaryImagePredicate* ptrImagePredicate;

// ----------------------------------------------------------------------
// utilities pointel
Point pointelRealPoint2Point( RealPoint p )
{
  RealPoint sp = RealPoint( round( p[ 0 ] / gridstep + 0.5 ),
                            round( p[ 1 ] / gridstep + 0.5 ),
                            round( p[ 2 ] / gridstep + 0.5 ) );
  return Point( sp[ 0 ], sp[ 1 ], sp[ 2 ] );
}
RealPoint pointelPoint2RealPoint( Point q )
{
  return RealPoint( gridstep * ( q[ 0 ] - 0.5 ),
                    gridstep * ( q[ 1 ] - 0.5 ),
                    gridstep * ( q[ 2 ] - 0.5 ) );
}
RealPoint pointelPoint2RealPoint( RealPoint q )
{
  return RealPoint( gridstep * ( q[ 0 ] - 0.5 ),
                    gridstep * ( q[ 1 ] - 0.5 ),
                    gridstep * ( q[ 2 ] - 0.5 ) );
}
void embedPointels( const std::vector< Point >& vq, std::vector< RealPoint >& vp )
{
  vp.resize( vq.size() );
  for ( auto i = 0; i < vp.size(); ++i )
    vp[ i ] = pointelPoint2RealPoint( vq[ i ] );
}
void embedPointels( const std::vector< RealPoint >& vq, std::vector< RealPoint >& vp )
{
  vp.resize( vq.size() );
  for ( auto i = 0; i < vp.size(); ++i )
    vp[ i ] = pointelPoint2RealPoint( vq[ i ] );
}
void digitizePointels( const std::vector< RealPoint >& vp, std::vector< Point >& vq )
{
  vq.resize( vp.size() );
  for ( auto i = 0; i < vq.size(); ++i )
    vq[ i ] = pointelRealPoint2Point( vp[ i ] );
}

// ----------------------------------------------------------------------
// utilities voxel
Point voxelRealPoint2Point( RealPoint p )
{
  RealPoint sp = RealPoint( round( p[ 0 ] / gridstep ),
                            round( p[ 1 ] / gridstep ),
                            round( p[ 2 ] / gridstep ) );
  return Point( sp[ 0 ], sp[ 1 ], sp[ 2 ] );
}
RealPoint voxelPoint2RealPoint( Point q )
{
  return RealPoint( gridstep * ( q[ 0 ] ),
                    gridstep * ( q[ 1 ] ),
                    gridstep * ( q[ 2 ] ) );
}
void embedVoxels( const std::vector< Point >& vq, std::vector< RealPoint >& vp )
{
  vp.resize( vq.size() );
  for ( auto i = 0; i < vp.size(); ++i )
    vp[ i ] = voxelPoint2RealPoint( vq[ i ] );
}
void digitizeVoxels( const std::vector< RealPoint >& vp, std::vector< Point >& vq )
{
  vq.resize( vp.size() );
  for ( auto i = 0; i < vq.size(); ++i )
    vq[ i ] = voxelRealPoint2Point( vp[ i ] );
}

// @return the indices of all the points of X different from p that are cotangent to p.
Indices tangentCone( const Point& p )
{
  return TC.getCotangentPoints( p );
}

Scalars distances( const Point& p, const Indices& idx )
{
  Scalars D( idx.size() );
  for ( Index i = 0; i < idx.size(); i++ )
    D[ i ] = ( TC.point( i ) - p ).norm();
  return D;
}

///////////////////////////////////////////////////////////////////////////////
std::vector< Point >
findCorners( const std::unordered_set< Point >& S,
             const std::vector< Vector >& In,
             const std::vector< Vector >& Out )
{
  std::vector< Point > C;
  for ( auto&& p : S )
    {
      bool corner = true;
      for ( auto&& n : In )
        if ( ! S.count( p+n ) ) { corner = false; break; }
      if ( ! corner ) continue;
      for ( auto&& n : Out )
        if ( S.count( p+n ) ) { corner = false; break; }
      if ( corner ) C.push_back( p );
    }
  return C;
}

void computeQuadrant( int q,
                      std::vector< Vector >& In,
                      std::vector< Vector >& Out )
{
  In.clear();
  Out.clear();
  In.push_back( Vector( q & 0x1 ? 1 : -1, 0, 0 ) );
  In.push_back( Vector( 0, q & 0x2 ? 1 : -1, 0 ) );
  In.push_back( Vector( 0, 0, q & 0x4 ? 1 : -1 ) );
  Out.push_back( Vector( q & 0x1 ? 1 : -1, q & 0x2 ? 1 : -1, q & 0x4 ? 1 : -1 ) );
  Vector D = ( In[ 1 ] - In[ 0 ] ).crossProduct( In[ 2 ] - In[ 0 ] );
  if ( D.dot( Out[ 0 ] ) < 0.0 ) std::swap( In[ 1 ], In[ 2 ] );
  In.push_back( In[ 0 ]+In[ 1 ] );
  In.push_back( In[ 0 ]+In[ 2 ] );
  In.push_back( In[ 1 ]+In[ 2 ] );
}


void computePlanes()
{
  // - mode specifies the candidate set, it is one of { ProbingMode::H, ProbingMode::R, ProbingMode::R1 }.
  using Estimator
    = DGtal::PlaneProbingTetrahedronEstimator< UnorderedPointSetPredicate,
                                               ProbingMode::R >;
  trace.beginBlock( "Compute planes" );
  std::vector< RealPoint > positions;
  std::vector< Point >     vertices;
  std::vector< std::vector<SH3::SurfaceMesh::Vertex> > faces;
  std::vector< Vector >  In;
  std::vector< Vector > Out;
  std::unordered_set< Point > S( digital_points.cbegin(), digital_points.cend() );
  UnorderedPointSetPredicate predS( S );
  Index i = 0;
  for ( int q = 0; q < 8; q++ ) // for each quadrant
    {
      computeQuadrant( q, In, Out );
      std::vector< Point > corners = findCorners( S, In, Out );
      std::cout << "Found " << corners.size() << " in Q" << q << std::endl;
      std::array<Point, 3> m = { In[ 0 ], In[ 1 ], In[ 2 ] };
      for ( auto&& p : corners )
        {
          Estimator estimator( p, m, predS );
          // if ( estimator.hexagonState()
          //      != Estimator::Neighborhood::HexagonState::Planar)
          //   continue;
          std::vector<SH3::SurfaceMesh::Vertex> triangle { i, i+1, i+2 };
          auto v = estimator.vertices();
          faces.push_back( triangle );
          vertices.push_back( v[ 0 ] );
          vertices.push_back( v[ 1 ] );
          vertices.push_back( v[ 2 ] );
          while (estimator.advance().first) {
            auto state = estimator.hexagonState();
            if (state == Estimator::Neighborhood::HexagonState::Planar) {
              auto v = estimator.vertices();
              if ( S.count( v[ 0 ] ) && S.count( v[ 1 ] ) && S.count( v[ 2 ] ) )
                {
                  std::vector< Point > X { v[ 0 ], v[ 1 ], v[ 2 ] };
                  auto P = dconv.makePolytope( X );
                  if ( dconv.isFullySubconvex( P, LS ) )
                      // // && TC.arePointsCotangent( v[ 0 ], v[ 1 ], v[ 2 ] ) ) 
                    {
                      vertices[ i   ] = v[ 0 ];
                      vertices[ i+1 ] = v[ 1 ];
                      vertices[ i+2 ] = v[ 2 ];
                    }
                }
            }          
          }
          i += 3;
          // }
        }
    }
  Time = trace.endBlock();
  embedPointels( vertices, positions );
  psTriangle = polyscope::registerSurfaceMesh("Triangle", positions, faces);
}

///////////////////////////////////////////////////////////////////////////////
void computeTangentCone()
{
  if ( digital_points.empty() ) return;
  if ( vertex_idx < 0 || vertex_idx >= digital_points.size() ) return;
  const auto p = digital_points[ vertex_idx ];
  trace.beginBlock( "Compute tangent cone" );
  auto local_X_idx = TC.getCotangentPoints( p );
  std::cout << "#cone=" << local_X_idx.size() << std::endl;
  local_X_idx.push_back( vertex_idx );
  std::vector< Point > local_X;
  std::vector< RealPoint > emb_local_X;
  for ( auto idx : local_X_idx )
    local_X.push_back( TC.point( idx ) );
  std::vector< double > values( local_X.size(), 0.0 );
  values.back() = 1.0;
  embedPointels( local_X, emb_local_X );
  psCloud = polyscope::registerPointCloud( "Tangent cone", emb_local_X );
  psCloud->setPointRadius( gridstep / 100.0 );
  psCloud->addScalarQuantity( "Classification", values );
  QuickHull3D hull;
  hull.setInput( local_X, false );
  hull.computeConvexHull();
  std::vector< Point > positions;
  std::vector< RealPoint > emb_positions;
  hull.getVertexPositions( positions );
  Time = trace.endBlock();
  embedPointels( positions, emb_positions );
  psCloudCvx = polyscope::registerPointCloud( "Tangent cone vertices", emb_positions );
  psCloudCvx->setPointRadius( gridstep / 50.0  );
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// SYMMETRIC CONVEX SET
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void computeSymmetricConvexSet()
{
  if ( digital_points.empty() ) return;
  if ( ! is_selected ) return;
  trace.beginBlock( "Compute symmetric convex set" );
  Point lo = Point::diagonal(-10000);
  Point up = Point::diagonal( 10000);
  UnorderedPointSetPredicate predS( unorderedSet );
  SymmetricConvexExpander< KSpace, UnorderedPointSetPredicate > SCE
    ( predS, selected_kpoint, lo, up );
  while ( SCE.advance( enforceFC ) )
    if ( PSym && ! SCE.myPerfectSymmetry 
         && SCE.current().second >= SCE.myPerfectSymmetryRadius )
      break;
  std::cout << "#symcvx=" << SCE.myPoints.size() << std::endl;
  Time = trace.endBlock();
  // Display as point cloud.
  std::vector< Point > positions( SCE.myPoints.cbegin(), SCE.myPoints.cend() );
  std::vector< RealPoint > emb_positions;
  embedPointels( positions, emb_positions );
  psCloud = polyscope::registerPointCloud( "Symmetric convex set", emb_positions );
  psCloud->setPointRadius( gridstep / 100.0 );
  // Display as disk
  PerfectSymmetricSet< Space > symset;
  bool ok = symset.init( SCE.kCenter(), SCE.myPerfectSymmetryRadius,
                         SCE.myPoints.cbegin(), SCE.myPoints.cend(), true );
  trace.info() << "Symmetric set is " << ( ok ? "OK" : "ERROR" ) << std::endl;
  std::vector< RealPoint > disks_positions;
  std::vector< RealPoint > embedded_disks_positions;
  std::vector<std::vector<std::size_t> > disks_faces;
  int index = 0;
  symset.addEllipse( index, disks_positions, disks_faces, 18 );
  embedPointels( disks_positions, embedded_disks_positions );
  psDisks = polyscope::registerSurfaceMesh( "Tangent ellipses",
                                            embedded_disks_positions, disks_faces );
  
}

void computeAllSymmetricConvexSet()
{
  if ( digital_points.empty() ) return;
  PerfectSymmetricTangentBundle< Space > TB;
  
  trace.beginBlock( "Compute all symmetric convex set" );
  Point lo = Point::diagonal(-10000);
  Point up = Point::diagonal( 10000);
  UnorderedPointSetPredicate predS( unorderedSet );
  for ( auto i = 0; i < cell_points.size(); i++ )
    {
      SymmetricConvexExpander< KSpace, UnorderedPointSetPredicate > SCE
        ( predS, cell_points[ i ], lo, up );
      while ( SCE.advance( enforceFC ) )
        if ( PSym && ! SCE.myPerfectSymmetry 
             && SCE.current().second >= SCE.myPerfectSymmetryRadius )
          break;
      PerfectSymmetricSet< Space > symset;
      if ( SCE.myPoints.size() >= MinNbPoints )
        {
          bool ok = symset.init( SCE.kCenter(), SCE.myPerfectSymmetryRadius,
                                 SCE.myPoints.cbegin(), SCE.myPoints.cend(), true );
          // Add to tangent bundle.
          if ( ok ) TB.add( symset, filterTB );
        }
      if ( i % 1000 == 0 )
        std::cout << i << "/" << cell_points.size()
                  << " #TB=" << TB.tangent_sets.size() << std::endl;
    }
  std::cout << "#TB=" << TB.tangent_sets.size() << std::endl;
  Time = trace.endBlock();
  // Display as disk
  std::vector< double >    eccentricities;
  std::vector< RealPoint > disks_positions;
  std::vector< RealPoint > embedded_disks_positions;
  std::vector<std::vector<std::size_t> > disks_faces;
  TB.addEllipses( disks_positions, disks_faces, eccentricities, 0, 18 );
  embedPointels( disks_positions, embedded_disks_positions );
  psDisks = polyscope::registerSurfaceMesh( "Tangent disks",
                                            embedded_disks_positions, disks_faces );
  psDisks->addFaceScalarQuantity("eccentricities", eccentricities);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// SYMMETRIC SEPARATOR
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void computeSymmetricSeparator()
{
  if ( digital_points.empty() ) return;
  if ( ! is_selected ) return;
  trace.beginBlock( "Compute symmetric convex set" );
  //SymmetricSeparator< KSpace, BinaryImagePredicate > SS( *ptrImagePredicate );
  UnorderedPointSetPredicate iinterior( immInterior );
  UnorderedPointSetPredicate iexterior( immExterior );
  SymmetricSeparator< KSpace,
                      UnorderedPointSetPredicate,
                      UnorderedPointSetPredicate > SS( iinterior, iexterior );
  bool ok = SS.init( selected_kpoint );
  while ( ok ) ok = fastSymSep ? SS.advanceFast() : SS.advance();
  std::cout << "#symsep=" << SS.myPoints.size() << std::endl;
  Time = trace.endBlock();
  // Display as point cloud.
  std::vector< Point > positions( SS.myPoints.cbegin(), SS.myPoints.cend() );
  std::vector< RealPoint > emb_positions;
  embedVoxels( positions, emb_positions );
  psCloud = polyscope::registerPointCloud( "Symmetric separator set", emb_positions );
  psCloud->setPointRadius( gridstep / 100.0 );
  //std::cout << "Sym=" << SS.kCenter() << std::endl;
  //std::cout << "In: ";
  //for ( auto p : positions ) std::cout << " " << p;
  for ( auto& p : positions ) p = SS.symmetric( p );
  //std::cout << std::endl << "Out:";
  //for ( auto p : positions ) std::cout << " " << p;
  //std::cout << std::endl;
  emb_positions.clear();
  embedVoxels( positions, emb_positions );
  psCloud = polyscope::registerPointCloud( "Symmetric separator out set", emb_positions );
  psCloud->setPointRadius( gridstep / 100.0 );
  // Display as point cloud.
  positions = std::vector< Point >( SS.myActive.cbegin(), SS.myActive.cend() );
  emb_positions.clear();
  embedVoxels( positions, emb_positions );
  psCloud = polyscope::registerPointCloud( "Active points", emb_positions );
  psCloud->setPointRadius( gridstep / 100.0 );

  // QuickHull3D hull;
  // hull.setInput( positions, false );
  // hull.computeConvexHull();
  // std::cout << hull << std::endl;
  // std::vector< Point > hull_positions;
  // hull.getVertexPositions( hull_positions );
  // std::vector< IndexRange > hull_facet_vertices;
  // bool ok2 = hull.getFacetVertices( hull_facet_vertices );
  // if ( ! ok2 ) trace.error() << "Bad facet computation" << std::endl;
  // polyscope::registerSurfaceMesh("Separator hull", hull_positions, hull_facet_vertices );

  SymmetricSeparatorSet< Z3i::Space > SSS( SS.myKCenter, SS.myPoints.cbegin(), SS.myPoints.cend() );
  std::vector< RealPoint >  sandwich_positions;
  std::vector< IndexRange > sandwich_facet_vertices;
  Index idx = 0;
  if ( SSS.isValid() )
    {
      SSS.addSeparator( idx, sandwich_positions, sandwich_facet_vertices );
      polyscope::registerSurfaceMesh("Sandwich", sandwich_positions, sandwich_facet_vertices );
    }
  bool ok4 = SSS.check();
  SSS.selfDisplay( std::cout );
  std::cout << ( ok4 ? " OK" : " ERROR" ) << std::endl;
  Thickness = SSS.thickness();
  // std::vector< Point > X( SS.myPoints.cbegin(), SS.myPoints.cend() );
  // for ( auto p : SS.myPoints ) X.push_back( SS.symmetric( p ) );
  // QuickHull3D hull;
  // hull.setInput( X, false );
  // hull.computeConvexHull();
  // std::vector< Point > hull_positions;
  // hull.getVertexPositions( hull_positions );
  // std::vector< IndexRange > hull_facet_vertices;
  // bool ok2 = hull.getFacetVertices( hull_facet_vertices );
  // if ( ! ok2 ) trace.error() << "Bad facet computation" << std::endl;
  // std::cout << hull << std::endl;
  // // Look for smallest diameter
  // typedef QuickHull3D::HalfSpace HalfSpace;
  // std::vector< HalfSpace > facet_halfspaces;
  // bool ok3 = hull.getFacetHalfSpaces( facet_halfspaces );
  // if ( ! ok3 ) {
  //   trace.error() << "Bad facet half-spaces computation" << std::endl;
  //   return;
  // }
  // int best_i = -1;
  // Scalar best_d = std::numeric_limits< Scalar >::infinity();
  // for ( auto i = 0; i < facet_halfspaces.size(); i++ )
  //   {
  //     Point  s = SS.symmetric( hull_positions[ hull_facet_vertices[ i ][ 0 ] ] );
  //     auto  FN = facet_halfspaces[ i ].internalNormal();
  //     auto  Fc = facet_halfspaces[ i ].internalIntercept();
  //     Scalar d = fabs( FN.dot( s ) - Fc );
  //     d       /= FN.norm();
  //     if  ( d < best_d ) { best_i = i; best_d = d; }
  //   }
  // const auto& fv = hull_facet_vertices[ best_i ];
  // const auto  nb = fv.size();
  // std::vector< RealPoint >  sandwich_positions( 2*nb );
  // std::vector< IndexRange > sandwich_facet_vertices( 2 );
  // sandwich_facet_vertices[ 0 ].resize( nb );
  // sandwich_facet_vertices[ 1 ].resize( nb );
  // for ( auto i = 0; i < nb; i++ )
  //   {
  //     const Point p = hull_positions[ fv[ i ] ];
  //     const Point q = SS.symmetric( p );
  //     sandwich_positions     [ i ]      = voxelPoint2RealPoint( p );
  //     sandwich_positions     [ nb + i ] = voxelPoint2RealPoint( q );
  //     sandwich_facet_vertices[ 0 ][ i ] = i;
  //     sandwich_facet_vertices[ 1 ][ i ] = nb + i;
  //   }
  // polyscope::registerSurfaceMesh("Sandwich", sandwich_positions, sandwich_facet_vertices );

}

void computeAllSymmetricSeparatorSet()
{
  if ( digital_points.empty() ) return;
  SymmetricSeparatorBundle< Space > TB;
  
  trace.beginBlock( "Compute all symmetric separators" );
  UnorderedPointSetPredicate iinterior( immInterior );
  UnorderedPointSetPredicate iexterior( immExterior );
  SymmetricSeparator< KSpace,
                      UnorderedPointSetPredicate,
                      UnorderedPointSetPredicate > SS( iinterior, iexterior );

  for ( auto i = 0; i < cell_points.size(); i++ )
    {
      if ( onlySurfels && K.uDim( K.uCell( cell_points[ i ] ) ) != 2 ) continue;
      // Compute Symmetric separator
      bool ok = SS.init( cell_points[ i ] );
      while ( ok ) ok = fastSymSep ? SS.advanceFast() : SS.advance();
      // Compute its proxy
      SymmetricSeparatorSet< Z3i::Space > SSS( SS.myKCenter,
                                               SS.myPoints.cbegin(), SS.myPoints.cend() );
      if ( SSS.isValid() && SSS.thickness() < maxThickness )
        TB.add( SSS, filterTB );
      if ( i % 1000 == 0 )
        trace.info() << i << "/" << cell_points.size()
                  << " #TB=" << TB.tangent_sets.size() << std::endl;
    }
  trace.info()<< "#TB=" << TB.tangent_sets.size() << std::endl;
  Time = trace.endBlock();

  // Display as separators
  // std::vector< double >    eccentricities;
  std::vector< RealPoint > seps_positions;
  std::vector<std::vector<std::size_t> > seps_faces;
  TB.addSeparators( seps_positions, seps_faces, 0 );
  psSeps = polyscope::registerSurfaceMesh( "Tangent seps",
                                            seps_positions, seps_faces );
  // psSeps->addFaceScalarQuantity("eccentricities", eccentricities);
}



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// TRIANGLES
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

struct TriangleContext
{
  Point   p;
  Indices cone_P; // the vector of indices of points in the tangent cone to p.
  Scalars d_P;    // the associated vector of distance to point p.
  Scalar  max_d;  // the maximum of d_P
  Scalar  min_d;  // the minimum distance between all pair of points in {P,A,B,C}
  Scalar  best;   // the current best D=d_P(A)+d_P(B)+d_P(C)+d_A(B)+d_A(C)+d_B(C)
  Scalar  cur_qA; 
  Scalar  cur_qAB; 
};

Scalar bestTriangleAB( TriangleContext& ctx,
                       const Index A, const Index B, Index& C )
{
  if ( ( ctx.best - ctx.cur_qAB ) > 4.0 * ctx.max_d ) return 0.0;
  Index cur_C;
  Scalar q = 0.0;
  const Point a = TC.point( ctx.cone_P[ A ] );
  const Point b = TC.point( ctx.cone_P[ B ] );
  for ( cur_C = B + 1; cur_C < ctx.cone_P.size(); cur_C++ )
    {
      const Scalar d_PC = ctx.d_P[ cur_C ];
      if ( d_PC < ctx.min_d ) continue;
      const Point     c = TC.point( ctx.cone_P[ cur_C ] );
      const Scalar d_AC = ( c - a ).norm();
      if ( d_AC < ctx.min_d ) continue;
      const Scalar d_BC = ( c - b ).norm();
      if ( d_BC < ctx.min_d ) continue;
      const Scalar cur_q = ctx.cur_qAB + d_PC + d_AC + d_BC
        + sqrt( ( ( b - a ).crossProduct( c - a ) ).norm() );
      if ( cur_q > q )
        {
          if ( ! TC.arePointsCotangent( a, c ) ) continue;
          if ( ! TC.arePointsCotangent( b, c ) ) continue;
          C = cur_C; q = cur_q;
        }
    }
  return q;
}

Scalar bestTriangleA( TriangleContext& ctx,
                      const Index A,
                      Index& B, Index& C )
{
  Index cur_B, cur_C;
  Scalar q = 0.0;
  const Point a = TC.point( ctx.cone_P[ A ] );
  for ( cur_B = A + 1; cur_B < ctx.cone_P.size(); cur_B++ )
    {
      const Scalar d_PB = ctx.d_P[ cur_B ]; 
      if ( d_PB < ctx.min_d ) continue;
      const Point     b = TC.point( ctx.cone_P[ cur_B ] );
      const Scalar d_AB = ( b - a ).norm();
      if ( d_AB < ctx.min_d ) continue;
      if ( ! TC.arePointsCotangent( a, b ) ) continue;
      ctx.cur_qAB = ctx.cur_qA + d_PB + d_AB;
      const Scalar result = bestTriangleAB( ctx, A, cur_B, cur_C );
      if ( result > q )
        {
          B = cur_B; C = cur_C; q = result;
        }
    }
  return q;
}

Scalar bestTriangle( TriangleContext& ctx,
                     Index& A, Index& B, Index& C )
{
  Index cur_A, cur_B, cur_C;
  for ( cur_A = 0; cur_A < ctx.cone_P.size(); cur_A++ )
    {
      std::cout << " " << cur_A;
      std::cout.flush();
      const Scalar d_PA = ctx.d_P[ cur_A ];
      if ( d_PA < ctx.min_d ) continue;
      ctx.cur_qA = d_PA;
      const Scalar result = bestTriangleA( ctx, cur_A, cur_B, cur_C );
      if ( result > ctx.best )
        {
          A = cur_A; B = cur_B; C = cur_C;
          ctx.best = result;
        }
    }
  return ctx.best;
}

void computeGreatTriangle()
{
  if ( digital_points.empty() ) return;
  if ( vertex_idx < 0 || vertex_idx >= digital_points.size() ) return;
  trace.beginBlock( "Compute best triangle" );
  TriangleContext ctx;
  ctx.p      = digital_points[ vertex_idx ];
  ctx.cone_P = tangentCone( ctx.p );
  std::cout << "#cone=" << ctx.cone_P.size() << std::endl;
  ctx.cone_P.push_back( vertex_idx );
  ctx.d_P    = distances( ctx.p, ctx.cone_P );
  ctx.max_d  = *( std::max_element( ctx.d_P.cbegin(), ctx.d_P.cend() ) );
  ctx.min_d  = ctx.max_d / 3.0;
  ctx.best   = 0.0;
  Index A, B, C;
  Scalar d = bestTriangle( ctx, A, B, C );
  Time = trace.endBlock();
  std::cout << "Best triangle " << A << " " << B << " " << C
            << " d=" << ctx.best << std::endl;
  if ( d == 0 ) return;
  std::vector< std::vector<SH3::SurfaceMesh::Vertex> > faces;
  std::vector<SH3::SurfaceMesh::Vertex> triangle { 0, 1, 2 };
  faces.push_back( triangle );
  std::vector<RealPoint> positions;
  std::vector< Point > vertices
    { TC.point( ctx.cone_P[ A ] ),
      TC.point( ctx.cone_P[ B ] ),
      TC.point( ctx.cone_P[ C ] ) };
  embedPointels( vertices, positions );
  psTriangle = polyscope::registerSurfaceMesh("Triangle", positions, faces);
}

void myCallback()
{
  // Select a vertex with the mouse
  if (polyscope::pick::haveSelection()) {
    bool goodSelection = false;
    auto selection = polyscope::pick::getSelection();
    auto selectedSurface = static_cast<polyscope::SurfaceMesh*>(selection.first);
    int idx = selection.second;

    // Only authorize selection on the input surface and the reconstruction
    auto surf = polyscope::getSurfaceMesh("Input surface");
    goodSelection = goodSelection || (selectedSurface == surf);
    const auto nv = selectedSurface->nVertices(); 
    const auto nf = selectedSurface->nFaces(); 
    const auto ne = selectedSurface->nEdges(); 
    // Validate that it its a face index
    if ( goodSelection )
      {
        if ( idx < nv )
          {
            vertex_idx = idx;
            is_selected  = true;
            selected_kpoint = digital_points[ vertex_idx ] * 2;
            std::ostringstream otext;
            otext << "Selected vertex = " << vertex_idx
                  << " pos=" << selected_kpoint; 
            ImGui::Text( "%s", otext.str().c_str() );
          }
        else if ( idx - nv < nf )
          {
            face_idx  = idx - nv;
            is_selected  = true;
            selected_kpoint = Point::zero;
            for ( auto i : selectedSurface->faces[ face_idx ] )
              selected_kpoint += digital_points[ i ];
            selected_kpoint /= 2;
            std::ostringstream otext;
            otext << "Selected face = " << face_idx
                  << " pos=" << selected_kpoint; 
            ImGui::Text( "%s", otext.str().c_str() );
          }
        else if ( idx - nv - nf < ne )
          {
            edge_idx  = idx - nv - nf;
            is_selected  = false; // true; // ne fonctionne pas
            selected_kpoint = Point::zero;
            for ( auto i : selectedSurface->edgeIndices[ edge_idx ] )
              selected_kpoint += digital_points[ i ];
            selected_kpoint /= 2;
            std::ostringstream otext;
            otext << "Selected edge = " << edge_idx
                  << " pos=" << selected_kpoint; 
            ImGui::Text( "%s", otext.str().c_str() );
          }
        else
          {
            vertex_idx = -1;
            face_idx   = -1;
            edge_idx   = -1;
            is_selected= false;
          }
      }
  }
  if (ImGui::Button("Compute tangent cone"))
    computeTangentCone();
  if (ImGui::Button("Compute symmetric separator"))
    computeSymmetricSeparator();
  ImGui::Checkbox("Fast symmetric separator", &fastSymSep );
  ImGui::Text( "Thickness = %f", Thickness );
  if (ImGui::Button("Compute all symmetric separators"))
    computeAllSymmetricSeparatorSet();
  ImGui::Checkbox("Only surfels", &onlySurfels );  
  ImGui::SliderFloat("#max thickness for symmetric set", &maxThickness, 0, 10.0);
  if (ImGui::Button("Compute symmetric convex set"))
    computeSymmetricConvexSet();
  ImGui::Checkbox("Perfect symmetry", &PSym );
  ImGui::Checkbox("Full convexity", &enforceFC );
  if (ImGui::Button("Compute all symmetric convex set"))
    computeAllSymmetricConvexSet();
  ImGui::SliderInt("#min points for symmetric set", &MinNbPoints, 5, 1000);
  ImGui::Checkbox("filter tangent bundle", &filterTB );
  if (ImGui::Button("Compute great triangle"))
    computeGreatTriangle();
  if (ImGui::Button("Compute planes"))
    computePlanes();
  ImGui::Text( "Computation time = %f ms", Time );
}

int main( int argc, char* argv[] )
{
  auto params = SH3::defaultParameters() | SHG3::defaultParameters() |  SHG3::parametersGeometryEstimation();
  params("surfaceComponents", "All");
  std::string filename = examplesPath + std::string("/samples/bunny-32.vol");
  if ( argc > 1 ) filename = std::string( argv[ 1 ] );
  auto binary_image    = SH3::makeBinaryImage(filename, params );
  K                    = SH3::getKSpace( binary_image, params );
  auto surface         = SH3::makeDigitalSurface( binary_image, K, params );
  auto surfels         = SH3::getSurfelRange( surface, params );
  auto primalSurface   = SH3::makePrimalSurfaceMesh(surface);
  SH3::Surfel2Index s2i;
  auto dualSurface     = SH3::makeDualPolygonalSurface( s2i, surface );  

  ptrImagePredicate = new BinaryImagePredicate( binary_image );
  
  //Need to convert the faces
  std::vector<std::vector<SH3::SurfaceMesh::Vertex>> faces;
  std::vector<RealPoint> positions;
  for(auto face= 0 ; face < primalSurface->nbFaces(); ++face)
    faces.push_back(primalSurface->incidentVertices( face ));
  
  //Recasting to vector of vertices
  positions = primalSurface->positions();
  
  surfmesh = SurfMesh(positions.begin(),
                      positions.end(),
                      faces.begin(),
                      faces.end());
  std::vector<std::vector<SH3::SurfaceMesh::Vertex>> dual_faces;
  std::vector<RealPoint> dual_positions;
  for(auto face= 0 ; face < dualSurface->nbFaces(); ++face)
    dual_faces.push_back( dualSurface->verticesAroundFace( face ));
    
    //Recasting to vector of vertices
  for ( auto vtx = 0; vtx < dualSurface->nbVertices(); ++vtx )
    dual_positions.push_back( dualSurface->position( vtx ) );
    
  dual_surfmesh = SurfMesh(dual_positions.begin(),
                           dual_positions.end(),
                           dual_faces.begin(),
                           dual_faces.end());
  std::cout << surfmesh << std::endl;
  std::cout << "number of non-manifold Edges = "
            << surfmesh.computeNonManifoldEdges().size() << std::endl;
  std::cout << dual_surfmesh << std::endl;
  std::cout << "number of non-manifold Edges = "
            << dual_surfmesh.computeNonManifoldEdges().size() << std::endl;
  // Make digital surface
  digitizePointels( positions, digital_points );
  trace.info() << "Surface has " << digital_points.size() << " pointels." << std::endl;
  dconv = DGtal::DigitalConvexity< KSpace >( K );
  TC    = DGtal::TangencyComputer< KSpace >( K );
  TC.init( digital_points.cbegin(), digital_points.cend() );
  LS    = DGtal::LatticeSetByIntervals< Space >
    ( digital_points.cbegin(), digital_points.cend(), 0 ).starOfPoints();
  trace.info() << "#cell_cover = " << TC.cellCover().nbCells() << std::endl;
  trace.info() << "#lattice_cover = " << LS.size() << std::endl;

  // Compute all cell points
  typedef Z3i::Cell Cell;
  std::set< Point > all_cells;
  for ( auto&& s : surfels )
    {
      const Cell us  = K.unsigns( s );
      const auto k   = K.uOrthDir( us );
      const Cell l0  = K.uIncident( us, (k+1)%3, false );
      const Cell l1  = K.uIncident( us, (k+1)%3, true );
      const Cell p00 = K.uIncident( l0, (k+2)%3, false );
      const Cell p01 = K.uIncident( l0, (k+2)%3, true );
      const Cell p10 = K.uIncident( l1, (k+2)%3, false );
      const Cell p11 = K.uIncident( l1, (k+2)%3, true );
      all_cells.insert( K.uKCoords( us ) );
      all_cells.insert( K.uKCoords( l0 ) );
      all_cells.insert( K.uKCoords( l1 ) );
      all_cells.insert( K.uKCoords( p00 ) );
      all_cells.insert( K.uKCoords( p01 ) );
      all_cells.insert( K.uKCoords( p10 ) );
      all_cells.insert( K.uKCoords( p11 ) );
    }
  cell_points = std::vector< Point >( all_cells.cbegin(), all_cells.cend() );
  all_cells.clear();
  std::cout << "Number of primal cells = " << cell_points.size() << std::endl;
  
  // Make predicate for pointel-based digital surface
  unorderedSet = std::unordered_set< Point >( digital_points.cbegin(),
                                              digital_points.cend() );

  // Prepare predicate for voxel based digital surface analysis
  for ( auto&& s : surfels )
    {
      const auto k       = K.sOrthDir( s );
      const auto direct  = K.sDirect( s, k );
      const auto int_vox  = K.sIncident( s, k, direct );
      const auto ext_vox = K.sIncident( s, k, ! direct );
      immInterior.insert( K.sCoords( int_vox ) );
      immExterior.insert( K.sCoords( ext_vox ) );
    }
  
  // Initialize polyscope
  polyscope::init();

  psMesh = polyscope::registerSurfaceMesh("Input surface", positions, faces);
  psDualMesh = polyscope::registerSurfaceMesh("Input dual surface", dual_positions, dual_faces);

  // Set the callback function
  polyscope::state::userCallback = myCallback;
  polyscope::show();
  return EXIT_SUCCESS;
}

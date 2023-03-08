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
 **/

#pragma once

/**
 * @file SymmetricConvexExpander.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5807), University of Savoie, France
 *
 * @date 2022/09/02
 *
 * This file is part of the DGtal library.
 */

#if defined(SymmetricConvexExpander_RECURSES)
#error Recursive header files inclusion detected in SymmetricConvexExpander.h
#else // defined(SymmetricConvexExpander_RECURSES)
/** Prevents recursive inclusion of headers. */
#define SymmetricConvexExpander_RECURSES

#if !defined SymmetricConvexExpander_h
/** Prevents repeated inclusion of headers. */
#define SymmetricConvexExpander_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/math/linalg/SimpleMatrix.h"
#include "DGtal/math/linalg/EigenDecomposition.h"
#include "DGtal/topology/CCellularGridSpaceND.h"
#include "DGtal/geometry/volumes/DigitalConvexity.h"
#include "DGtal/geometry/surfaces/ParallelStrip.h"
#include "DGtal/geometry/surfaces/COBAGenericStandardPlaneComputer.h"
#include "DGtal/shapes/PolygonalSurface.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // class SymmetricConvexExpander
  /**
   * Description of class 'SymmetricConvexExpander' <p>
   *
   * \brief Aim: SymmetricConvexExpander computes symmetric fully
   * convex subsets of a given digital set.
   *
   * @tparam TKSpace an arbitrary model of CCellularGridSpaceND.
   * @tparam TPointPredicate an arbitrary model of predicate Point -> bool
   **/

  template < typename TKSpace,
             typename TPointPredicate >
  class SymmetricConvexExpander
  {
    BOOST_CONCEPT_ASSERT(( concepts::CCellularGridSpaceND< TKSpace > ));
 
  public:
    typedef SymmetricConvexExpander< TKSpace, TPointPredicate > Self;
    typedef TKSpace                         KSpace;
    typedef TPointPredicate                 PointPredicate;
    typedef typename KSpace::Integer        Integer;
    typedef typename KSpace::Point          Point;
    typedef typename KSpace::Vector         Vector;
    typedef typename KSpace::Space          Space;
    typedef std::size_t                     Size;
    typedef DGtal::BoundedLatticePolytope < Space > LatticePolytope;
    typedef DGtal::BoundedRationalPolytope< Space > RationalPolytope;
    typedef std::vector<Point>              PointRange;
    typedef std::vector<Vector>             VectorRange;
    typedef std::unordered_set<Point>       PointSet;
    
    static const Dimension dimension = KSpace::dimension;

    typedef std::pair< Point, Integer >     Node;

    // Inversion order since priority queue output max element.
    struct NodeComparator {
      /// Default constructor. 
      NodeComparator() = default;
      // p < q iff p.second < q.second
      bool operator()( const Node& p, const Node& q ) const
      {
        return p.second > q.second;
      }
    };

    typedef std::priority_queue< Node, std::vector<Node>, NodeComparator > NodeQueue;

    /// Constructor from predicate and symmetry center point.
    SymmetricConvexExpander( const PointPredicate& predicate,
                             const Point& kcenter,
                             const Point& lo, 
                             const Point& hi )
      : myPredicate( &predicate ), myConvexity( lo, hi )
    {
      init( kcenter );
    }

    /// @return the center with doubled coordinates.
    Point kCenter() const { return myKCenter; }
    
    bool predicate( const Point& p ) const
    {
      return (*myPredicate)( p );
    }
    
    void init( const Point& kcenter )
    {
      myKCenter = kcenter;
      myPoints.clear();
      myQ = NodeQueue();
      myM.clear();
      myPerfectSymmetry       = true;
      myPerfectSymmetryRadius = 0;
      // The starting points depend on the parity of the coordinates of the center.
      // There are from 1 to 2^d starting points.
      PointRange points;
      const auto x = myKCenter[ 0 ];
      if ( x % 2 == 0 )
        points.push_back( Point::base( 0, x / 2 ) );
      else
        {
          points.push_back( Point::base( 0, (x-1) / 2 ) );
          points.push_back( Point::base( 0, (x+1) / 2 ) );
        }
      for ( Dimension k = 1; k < dimension; k++ )
        {
          const auto n = points.size();
          const auto y = myKCenter[ k ];
          if ( y % 2 == 0 )
            {
              for ( auto i = 0; i < n; i++ )
                points[ i ][ k ] = y / 2;
            }
          else
            {
              points.resize( 2*n );
              const auto z  = (y-1)/2;
              const auto z1 = z + 1;
              for ( auto i = 0; i < n; i++ )
                {
                  points[ i ][ k ] = z;
                  Point q = points[ i ];
                  q[ k ]  = z1;
                  points[ i+n ]    = q;
                }                  
            }
        }
      // Keep only the points that satisfy the predicate.
      for ( auto&& p : points )
        {
          const Point sp = symmetric( p );
          if ( ! myM.count( p )
               && predicate( p ) && predicate( sp ) )
            {
              Node n( p, (2*p - myKCenter).squaredNorm() );
              myQ.push( n );
              myM.insert( p );
              myM.insert( sp );
            }
        }
    }

    /// Advance of one symmetric point.
    bool advance( bool enforce_full_convexity )
    {
      while ( ! finished() )
        {
          const auto p = current().first;
          const auto d = current().second; // current ring distance
          const auto sp = symmetric( p );
          myPoints.insert( p );
          myPoints.insert( sp );
          PointRange X( myPoints.cbegin(), myPoints.cend() );
          if ( enforce_full_convexity && ! myConvexity.isFullyConvex( X ) )
            {
              myPoints.erase( p );
              myPoints.erase( sp );
              ignore();
            }
          else
            {
              expand();
              return true;
            }
        }
      return false;
    }
    
    /**
       @return a const reference on the current visited vertex. The
       node is a pair <Vertex,Data> where the second term is the
       topological distance to the start vertex or set.

       NB: valid only if not 'finished()'.
     */
    const Node& current() const
    {
      return myQ.top();
    }

    /**
       Goes to the next vertex but ignores the current vertex for
       determining the future visited vertices. Otherwise said, no
       future visited vertex will have this vertex as a father.

       NB: valid only if not 'finished()'.
     */
    void ignore()
    {
      myQ.pop();
    }

    /**
       Goes to the next vertex and take into account the current
       vertex for determining the future visited vertices.
       NB: valid only if not 'finished()'.
     */
    void expand()
    {
      const Point p = current().first;
      myQ.pop();
      myPoints.insert( p );
      myPoints.insert( symmetric( p ) );
      const auto next_points = next( p );
      for ( auto&& q : next_points )
        {
          if ( ! myM.count( q ) )
            {
              const auto  sq   = symmetric( q );
              const bool  q_in = predicate( q );
              const bool sq_in = predicate( sq );
              const Integer d2 = (2*q - myKCenter).squaredNorm();
              if ( q_in && sq_in )
                {
                  Node n( q, d2 );
                  //Node sn( sq, n.second ); //(2*q - myKCenter).squaredNorm() );
                  myQ.push( n );
                  //myQ.push( sn );
                  myM.insert(  q );
                  myM.insert( sq );
                  // if ( myPerfectSymmetry )
                  //   myPerfectSymmetryRadius = std::max( myPerfectSymmetryRadius,
                  //                                       n.second );
                }
              else if ( ( q_in && ! sq_in ) || ( ! q_in && sq_in ) )
                {
                  if ( myPerfectSymmetry ) {
                    //std::cout << (2*q - myKCenter) << " " << d2 << std::endl;
                    myPerfectSymmetryRadius = d2;
                    myPerfectSymmetry = false;
                  } else {
                    myPerfectSymmetryRadius = std::min( d2, myPerfectSymmetryRadius );
                  }
                }
            }
        }
    }

    /**
       @return 'true' if all possible elements have been visited.
     */
    bool finished() const
    {
      return myQ.empty();
    }
    
    
    Point symmetric( const Point& p ) const
    {
      return myKCenter - p;
    }

    PointRange next( const Point& p ) const
    {
      PointRange N;
      Point d = 2*p - myKCenter;
      for ( Dimension i = 0; i < dimension; i++ )
        {
          if ( d[ i ] >= 0 ) N.push_back( p + Point::base( i, 1 ) );
          if ( d[ i ] <= 0 ) N.push_back( p - Point::base( i, 1 ) );
        }
      return N;
    }

    /// The predicate that every point must satisfy
    const PointPredicate* myPredicate;

    /// The digital convexity object
    DigitalConvexity< KSpace > myConvexity;
    
    /// Symmetry center (with doubled coordinates to represent half-integers).
    Point myKCenter;

    /// Symmetric range of lattice points, sorted.
    PointSet myPoints;

    /// Queue of points (only one point is inserted in the queue for a
    /// symmetric pair of point).
    NodeQueue  myQ;
    /// Marked points, i.e. points already in the queue or in the object.
    PointSet   myM;

    /// True iff the set and its local complement are symmetric.
    bool myPerfectSymmetry;
    /// Upper bound on the max distance of perfect symmetry.
    Integer myPerfectSymmetryRadius;
    
  };

template < typename TSpace >
struct PerfectSymmetricSet {
  typedef TSpace                            Space;
  typedef typename Space::Point             Point;
  typedef typename Space::RealPoint         RealPoint;
  typedef typename RealPoint::Coordinate    Scalar;
  typedef typename Space::Integer           Integer;
  typedef DGtal::int64_t                    InternalInteger;
  typedef DGtal::COBAGenericStandardPlaneComputer< Space, InternalInteger > PlaneComputer;
  typedef typename PlaneComputer::Primitive Plane;
  typedef std::vector< std::size_t >        Indices;
  typedef typename Space::RealVector        RealVector;
  typedef DGtal::SimpleMatrix< double, 3, 3 > RealTensor;
  
  Point   myKCenter; ///< center with doubled coordinates
  RealPoint myCenter; ///< center with real coordinates
  Scalar myRadius; ///< radius of set rounded above.
  Integer mySquaredRadius; ///< squared radius of set
  Plane   myPlane; ///< the standard plane in which lies the set
  RealTensor myPCA;
  
  PerfectSymmetricSet() = default;

  /// Initialise the set with a center, its squared radius and the range of points.
  /// @return 'true' iff all the range of points belong to a standard plane.
  template <typename PointIterator>
  bool init( Point kcenter, Integer sqradius, PointIterator b, PointIterator e,
             bool ellipsoid = false )
  {
    myKCenter = kcenter;
    myCenter  = RealPoint( myKCenter[ 0 ], myKCenter[ 1 ], myKCenter[ 2 ] );
    myCenter /= 2.0;
    mySquaredRadius = sqradius;
    myRadius = 0.5*sqrt( double( sqradius ) );
    // Determine the standard plane that includes the given range of points.
    PlaneComputer PC;
    PC.init( InternalInteger( 2*myRadius ) );
    bool ok = true;
    for ( auto it = b ; it != e; it++ )
      ok = ok && PC.extend( *it );
    myPlane = PC.primitive();
    // Compute PCA
    if ( ! ellipsoid ) return ok;
    int n = 0;
    for ( auto it = b ; it != e; it++, n++ )
      {
        auto v = (2 * (*it)) - kcenter;
        for ( auto i = 0; i < 3; i++ )
          for ( auto j = 0; j < 3; j++ )
            myPCA( i, j ) += double( v[ i ] ) * double( v[ j ] );
      }
    myPCA /= 4.0 * double( n-1 );
    return ok;
  }
  
  bool isInside( const PerfectSymmetricSet& other ) const
  {
    if ( mySquaredRadius > other.mySquaredRadius ) return false;
    if ( ! other.myPlane( myCenter ) ) return false;
    double d = (other.myCenter - myCenter).norm();
    return ( d + myRadius ) < other.myRadius;
  }

  // return the eccentricity, i.e. 0
  double addDisk( int& index, 
                  std::vector< RealPoint >& positions,
                  std::vector< Indices >&  faces,
                  const int nb ) const
  {
    RealPoint x( 1.0, 0.0, 0.0 );
    RealPoint y( 0.0, 1.0, 0.0 );
    RealPoint u = x.crossProduct( myPlane.normal() );
    if ( u.norm() < 0.1 ) u = y.crossProduct( myPlane.normal() );
    RealPoint v = myPlane.normal().crossProduct( u );
    u /= u.norm();
    v /= v.norm();
    Indices face( nb );
    for ( int i = 0; i < nb; i++ )
      {
        const double angle = double(i) * 2.0 * M_PI / double( nb );
        const RealPoint p  = myCenter + myRadius * ( cos(angle) * u + sin(angle) * v );
        positions.push_back( p );
        face[ i ] = index + i;
      }
    faces.push_back( face );
    index += nb;
    return 0.0;
  }

  // @return the eccentricity (>= 0.0)
  double addEllipse( int& index, 
                     std::vector< RealPoint >& positions,
                     std::vector< Indices >&  faces,
                     const int nb ) const
  {
    // Diagonalize PCA
    RealTensor V; // unit eigen vectors in columns V.column( 0 ) = normal vector
    RealVector E; //  eigen values (smallest to highest).
    DGtal::EigenDecomposition< 3, double >::getEigenDecomposition( myPCA, V, E );
    RealVector u = V.column( 2 ); // great ellipse axis
    RealVector v = V.column( 1 ); // small ellipse axis
    u *= 2.0*sqrt( E[ 2 ] );
    v *= 2.0*sqrt( E[ 1 ] );
    Indices face( nb );
    for ( int i = 0; i < nb; i++ )
      {
        const double angle = double(i) * 2.0 * M_PI / double( nb );
        const RealPoint p  = myCenter + ( cos(angle) * u + sin(angle) * v );
        positions.push_back( p );
        face[ i ] = index + i;
      }
    faces.push_back( face );
    index += nb;
    return sqrt( std::max( 0.0, 1.0 - E[ 1 ] / E[ 2 ] ) );
  }
};

//  Simplified version of PerfectSymmetricSet that computes only its PCA.
template < typename TSpace >
struct PCAPerfectSymmetricSet {
  typedef TSpace                            Space;
  typedef typename Space::Point             Point;
  typedef typename Space::RealPoint         RealPoint;
  typedef typename RealPoint::Coordinate    Scalar;
  typedef typename Space::Integer           Integer;
  typedef DGtal::int64_t                    InternalInteger;
  typedef std::vector< std::size_t >        Indices;
  typedef typename Space::RealVector        RealVector;
  typedef DGtal::SimpleMatrix< double, 3, 3 > RealTensor;
  
  Point   myKCenter; ///< center with doubled coordinates
  RealPoint myCenter; ///< center with real coordinates
  RealTensor myPCA; ///< PCA tensor
  RealVector myA;   ///< the greatest half-axis
  RealVector myB;   ///< the smallest half-axis
  RealVector myN;   ///< the normal direction (as a unit vector)
  
  PCAPerfectSymmetricSet() = default;

  /// Initialise the set with a center and the range of points.
  template <typename PointIterator>
  void init( Point kcenter, PointIterator b, PointIterator e )
  {
    myKCenter = kcenter;
    myCenter  = RealPoint( myKCenter[ 0 ], myKCenter[ 1 ], myKCenter[ 2 ] );
    myCenter /= 2.0;
    // Compute PCA
    int n = 0;
    for ( auto it = b ; it != e; it++, n++ )
      {
        auto v = (2 * (*it)) - kcenter;
        for ( auto i = 0; i < 3; i++ )
          for ( auto j = 0; j < 3; j++ )
            myPCA( i, j ) += double( v[ i ] ) * double( v[ j ] );
      }
    myPCA /= 4.0 * double( n-1 );

    // Diagonalize PCA
    RealTensor V; // unit eigen vectors in columns V.column( 0 ) = normal vector
    RealVector E; //  eigen values (smallest to highest).
    DGtal::EigenDecomposition< 3, double >::getEigenDecomposition( myPCA, V, E );
    myA = V.column( 2 ) * ( 2.0*sqrt( E[ 2 ] ) ); // great ellipse axis
    myB = V.column( 1 ) * ( 2.0*sqrt( E[ 1 ] ) ); // small ellipse axis
    myN = myA.crossProduct( myB ).getNormalized(); // normal direction (not oriented).
  }
  
};


template < typename TSpace >
struct PerfectSymmetricTangentBundle {
  typedef TSpace                        Space;
  typedef PerfectSymmetricSet< Space > SymmetricSet;
  typedef typename Space::Point     Point;
  typedef typename Space::RealPoint RealPoint;
  typedef typename RealPoint::Coordinate Scalar;
  typedef std::vector< std::size_t > Indices;

  std::vector< SymmetricSet > tangent_sets;
  
  PerfectSymmetricTangentBundle() = default;

  void add( const SymmetricSet& sset, bool filter )
  {
    if ( filter )
      {
        for ( const auto& s : tangent_sets )
          if ( sset.isInside( s ) ) return;
        std::size_t i = 0;
        while ( i < tangent_sets.size() )
          {
            if ( tangent_sets[ i ].isInside( sset ) )
              {
                std::swap( tangent_sets[ i ], tangent_sets.back() );
                tangent_sets.pop_back();
              }
            else i++;
          }
      }
    tangent_sets.push_back( sset );
  }

  void addDisks( std::vector< RealPoint >& positions,
                 std::vector< Indices >&  faces,
                 int index,
                 const int nb ) const
  {
    for ( const auto& s : tangent_sets )
      s.addDisk( index, positions, faces, nb );
  }
  void addEllipses( std::vector< RealPoint >& positions,
                    std::vector< Indices >&   faces,
                    std::vector< Scalar >&    eccentricities,
                    int index,
                    const int nb ) const
  {
    for ( const auto& s : tangent_sets )
      {
        double e = s.addEllipse( index, positions, faces, nb );
        eccentricities.push_back( e );
      }
  }
};



  // Represents the reconstruction of the dual surface with fully convex tangent planes.
template < typename TKSpace >
struct DualReconstruction {
  typedef TKSpace KSpace;
  typedef std::size_t                 Size;
  typedef std::size_t                 Index;
  typedef typename KSpace::Point      Point;
  typedef typename KSpace::Vector     Vector;
  typedef typename KSpace::Space      Space;
  typedef typename Space::RealPoint   RealPoint;
  typedef typename Space::RealVector  RealVector;
  typedef typename Vector::Component  Integer;
  typedef typename RealVector::Component Scalar;
  typedef typename KSpace::Cell       Cell;
  typedef typename KSpace::SCell      SCell;
  typedef std::pair< Vector, Scalar > Plane;
  typedef DGtal::SimpleMatrix< double, 3, 3 > RealTensor;

  typedef ::DGtal::PolygonalSurface<RealPoint> PolygonalSurface;
  typedef typename PolygonalSurface::Vertex Vertex;

  KSpace                   K; ///< Khalimsky space in which everything lies.
  CountedPtr< PolygonalSurface > dualSurface; ///< the dual surface
  std::vector< SCell >     surfels;    ///< the  range of surfels
  std::map< Point,Index >  surfel2idx; ///< the map surfel kcoords -> its index
  std::vector< Scalar >    intercepts; ///< the array surfel idx -> height of dual reconstruction.
  std::vector< Index >     plane_indices; ///< the array surfel idx -> plane idx
  std::vector< Plane >     planes;     ///< the array plane idx -> plane of dual reconstruction.
  std::vector< std::vector<Index> > dual_faces;
  std::vector< RealPoint >          dual_positions;
  double                   gridstep;
  IntegerComputer< Integer > ic;
  std::vector< RealVector > dual_vertex_normals;
  std::vector< RealVector > dual_face_normals;
  
  DualReconstruction( const KSpace& aK,
                      CountedPtr< PolygonalSurface > aDualSurface,
                      const std::vector< SCell >&    theSurfels,
                      const std::map< Point,Index >& s2i,
                      double h = 1.0 )
    : K( aK ), dualSurface( aDualSurface ), surfels( theSurfels ), surfel2idx( s2i ),
      gridstep( h )
  {
    intercepts       = std::vector< double >( surfels.size(), 0.0 );
    for ( Index i = 0; i < surfels.size(); i++ )
      {
        SCell     s = surfels[ i ];
        Dimension k = K.sOrthDir( s );
        bool orient = K.sDirect( s, k ); // toward interior
        Point invox = K.sCoords( K.sIncident( s, k, orient   ) );
        Point exvox = K.sCoords( K.sIncident( s, k, ! orient ) );
        Vector    N = exvox - invox;
        plane_indices.push_back( i );
        planes.push_back( std::make_pair( N, N.dot( invox ) ) );
        dual_vertex_normals.push_back( RealVector( N[ 0 ], N[ 1 ], N[ 2 ] ) );
      }
    for(auto face= 0 ; face < dualSurface->nbFaces(); ++face)
      {
        auto V = dualSurface->verticesAroundFace( face );
        dual_faces.push_back( V );
        RealVector N;
        for ( auto v : V ) N += dual_vertex_normals[ v ];
        N /= N.norm();
        dual_face_normals.push_back( N );
      }
    
    //Recasting to vector of vertices
    for ( auto vtx = 0; vtx < dualSurface->nbVertices(); ++vtx )
      dual_positions.push_back( dualSurface->position( vtx ) );
    std::cout << "#surfels=" << surfels.size() << std::endl;
    std::cout << "#intercepts=" << intercepts.size() << std::endl;
    std::cout << "#dual_positions=" << dual_positions.size() << std::endl;
  }

  void updateDualPositions()
  {
    std::cout << "Updating positions" << std::endl;

    for ( auto i = 0; i < surfels.size(); i++ )
      {
        auto int_vox = K.interiorVoxel( surfels[ i ] );
        auto ext_vox = K.exteriorVoxel( surfels[ i ] );
        auto int_p   = voxelPoint2RealPoint( int_vox ); //K.uCoords( int_vox ) );
        auto ext_p   = voxelPoint2RealPoint( ext_vox ); //K.uCoords( ext_vox ) );
        const double    s = intercepts[ i ] + 0.001; 
        const RealPoint q = ( 1.0 - s ) * int_p + s * ext_p;
        dual_positions[ i ] = q;
      }
  }

  void reduce( Vector& N )
  {
    auto g1 = ic.gcd( N[ 0 ], N[ 1 ] );
    auto g2 = ic.gcd( N[ 1 ], N[ 2 ] );
    auto g  = ic.gcd( g1, g2 );
    N      /= g;
  }
  /// @param X the vertices of the current face
  /// @param cells the cells intersected by the face.
  void updateReconstructionFromCells( const std::vector< Point >& X,
                                      const std::vector< Point >& cells )
  {
    // Compute plane
    Vector N = ( X[ 1 ] - X[ 0 ] ).crossProduct( X[ 2 ] - X[ 0 ] );
    reduce( N );
    const auto   N1 = N.norm( Vector::NormType::L_1 );
    const auto    a = N.dot( X[ 0 ] );
    Index plane_idx_pos = planes.size();
    Index plane_idx_neg = planes.size()+1;
    int          nb = 0;
    // Extract 1-cells which are dual to surfels
    for ( auto&& kp : cells )
      {
        // Look for dimension 1 cells.
        const Cell c = K.uCell( kp );
        if ( K.uDim( c ) != 1 ) continue;
        // Compute dual surfel
        const Dimension t = *K.uDirs( c );
        const Cell p0 = K.uIncident( c, t, false );
        const Cell p1 = K.uIncident( c, t, true );
        const Point dual_kp = K.uCoords( p0 ) + K.uCoords( p1 ) + Point::diagonal(1);
        const auto it = surfel2idx.find( dual_kp );
        if ( it == surfel2idx.cend() ) continue;
        // Compute and update intercept
        const Size  idx     = it->second;
        const SCell surfel  = surfels[ idx ];
        const Point int_vox = K.interiorVoxel( surfel );
        const Point ext_vox = K.exteriorVoxel( surfel );
        const auto  int_val = N.dot( int_vox );
        const auto  ext_val = N.dot( ext_vox );
        // std::cout << " int_val=" << int_val << " a=" << a << " ext_val=" << ext_val;
        // if ( ( int_val <= a && ext_val <= a ) || ( int_val >= a && ext_val >= a ) )
        // {
        //   if ( ( int_val < a && ext_val < a ) || ( int_val > a && ext_val > a ) )
        //     trace.warning() << "Bad intersection" << std::endl;
        //   continue;
        // }
        if ( int_val == ext_val && a != int_val ) continue;
        const double s     = (a == int_val)
          ? 0.0
          : (double)( a - int_val ) / (double) (ext_val - int_val );
        const double old_s = intercepts[ idx ];
        // const auto     dot = N.dot( planes[ plane_indices[ idx ] ].first );
        // if ( old_s < s ) std::cout  << " s=" << old_s << " -> " << s << std::endl;
        if ( ( old_s < ( s - 0.0001 ) )
             || ( ( fabs( old_s - s ) <= 0.0001 )
                  && ( planes[ plane_indices[ idx ] ].first.norm( Vector::NormType::L_1 ) ) < N1 ) )
          { // this face intersects the dual linel above.
            intercepts[ idx ] = std::max( old_s, s );
            nb += 1;
            plane_indices[ idx ] = ( ext_val > int_val ) ? plane_idx_pos : plane_idx_neg;
          }
      }
    if ( nb > 0 ) {
      planes.push_back( std::make_pair( N, a ) );
      planes.push_back( std::make_pair( -N, -a ) );
    }
  }

  RealPoint voxelPoint2RealPoint( Point q ) const
  {
    return RealPoint( gridstep * ( q[ 0 ] ),
                      gridstep * ( q[ 1 ] ),
                      gridstep * ( q[ 2 ] ) );
  }
  Point pointelRealPoint2Point( RealPoint p ) const
  {
    RealPoint sp = RealPoint( round( p[ 0 ] / gridstep + 0.5 ),
                              round( p[ 1 ] / gridstep + 0.5 ),
                              round( p[ 2 ] / gridstep + 0.5 ) );
    return Point( sp[ 0 ], sp[ 1 ], sp[ 2 ] );
  }
  RealPoint pointelPoint2RealPoint( RealPoint q ) const
  {
    return RealPoint( gridstep * ( q[ 0 ] - 0.5 ),
                      gridstep * ( q[ 1 ] - 0.5 ),
                      gridstep * ( q[ 2 ] - 0.5 ) );
  }

  bool isSamePlane( Index p_idx1, Index p_idx2 ) const
  {
    return ( planes[ p_idx1 ].first  == planes[ p_idx2 ].first )
      &&   ( planes[ p_idx1 ].second == planes[ p_idx2 ].second );
  }

  void cleanUpPlaneFaces()
  {
    std::set< Index > all_plane_indices;
    for ( auto idx : plane_indices ) all_plane_indices.insert( idx );
    std::vector< Plane > planes_copy = planes;
    planes.resize( all_plane_indices.size() );
    Index nidx = 0;
    std::map< Index, Index > renumber;
    for ( auto idx : all_plane_indices )
      {
        renumber[ idx ] = nidx;
        planes[ nidx ]  = planes_copy[ idx ];
        nidx           += 1;
      }
    for ( int i = 0; i < plane_indices.size(); i++ )
      plane_indices[ i ] = renumber[ plane_indices[ i ] ];
  }
  
  void updatePlaneFaces()
  {
    Size merged = 0;
    for ( Vertex v = 0; v < dualSurface->nbVertices(); ++v )
      {
        std::vector< Vertex > N;
        auto it = std::back_inserter( N );
        dualSurface->writeNeighbors( it, v );
        for ( auto w : N )
          {
            if ( v < w )
              {
                auto pv = plane_indices[ v ];
                auto pw = plane_indices[ w ];
                if ( ( pv != pw ) && isSamePlane( pv, pw ) )
                  {
                    plane_indices[ w ] = pv;
                    merged += 1;
                  }
              }
          }
      }
    std::cout << "Merged " << merged << " planar faces" << std::endl;
  }

  std::vector< RealVector > getNormalVectors() const
  {
    std::vector< RealVector > normals( surfels.size() );
    for ( Index i = 0; i < surfels.size(); i++ )
      {
        const auto  N = planes[ plane_indices[ i ] ].first;
        RealVector rN = RealVector( N[ 0 ], N[ 1 ], N[ 2 ] );
        normals[ i ]  = rN.getNormalized(); 
      }
    return normals;
  }
  std::vector< Vector > getExactNormalVectors() const
  {
    std::vector< Vector > exact_normals( surfels.size() );
    for ( Index i = 0; i < surfels.size(); i++ )
      {
        const auto  N = planes[ plane_indices[ i ] ].first;
        //RealVector rN = RealVector( N[ 0 ], N[ 1 ], N[ 2 ] );
        exact_normals[ i ]  = N;
      }
    return exact_normals;
  }
  std::vector< Scalar > getExactIntercepts() const
  {
    std::vector< Scalar > exact_intercepts( surfels.size() );
    for ( Index i = 0; i < surfels.size(); i++ )
      {
        const auto  alpha = planes[ plane_indices[ i ] ].second;
        exact_intercepts[ i ]  = alpha;
      }
    return exact_intercepts;
  }

  std::vector< RealPoint > getLowerPoints() const
  {
    std::vector< RealPoint > L;
    for ( Index i = 0; i < surfels.size(); i++ )
      {
        if ( intercepts[ i ] == 0.0 )
          {
            auto int_vox = K.interiorVoxel( surfels[ i ] );
            auto int_p   = voxelPoint2RealPoint( int_vox );
            L.push_back( int_p );
          }
      }
    return L;
  }

  RealPoint pointelPosition( RealPoint pp,
                             std::vector< Index > surfel_indices, int & type ) const
  {
    std::set< Index > local_plane_indices;
    for ( auto idx : surfel_indices ) local_plane_indices.insert( plane_indices[ idx ] );
    std::vector< RealVector > A;
    std::vector< Scalar >     B;
    for ( auto pidx : local_plane_indices )
      {
        Vector N = planes[ pidx ].first;
        RealVector rN( N[ 0 ], N[ 1 ], N[ 2 ] );
        Scalar l = rN.norm(); 
        A.push_back( rN / l );
        B.push_back( Scalar( planes[ pidx ].second ) / l );
      }
    type = A.size();
    // we project p onto the first plane, since all of them are equivalent or parallel.
    RealPoint q;
    for ( auto i = 0; i < A.size(); i++ )
      {
        Scalar t = ( pp.dot( A[ i ] ) - B[ i ] ) / A[ i ].squaredNorm();
        q += ( pp - t * A[ i ] );
      }
    return q / A.size();
  }

  RealPoint pointelPosition2( RealPoint pp,
                             std::vector< Index > surfel_indices, int & type ) const
  {
    typedef DGtal::SimpleMatrix< double, 3, 3 > RealTensor;
    auto p = pp;// pointelRealPoint2Point( pp );
    RealPoint pointel( p[ 0 ], p[ 1 ], p[ 2 ] );
    std::set< Index > local_plane_indices;
    for ( auto idx : surfel_indices ) local_plane_indices.insert( plane_indices[ idx ] );
    std::vector< RealVector > A;
    std::vector< Scalar >     B;
    for ( auto pidx : local_plane_indices )
      {
        Vector N = planes[ pidx ].first;
        RealVector rN( N[ 0 ], N[ 1 ], N[ 2 ] );
        Scalar l = rN.norm(); 
        A.push_back( rN / l );
        B.push_back( Scalar( planes[ pidx ].second ) / l );
      }
    type = 3;
    if ( A.size() >= 3 )
      {  // Trying full rank. We perform linear regression
        RealTensor M;
        RealVector C;
        for ( int k = 0; k < A.size(); k++ )
          for ( int j = 0; j < 3; j++ ) {
            for ( int i = 0; i < 3; i++ ) {
              M( i, j ) += A[ k ][ i ] * A[ k ][ j ];
            }
            C[ j ] += A[ k ][ j ] * B[ k ];
          }
        const auto d  = M.determinant();
        if ( d != 0 ) {
          const auto P = M.inverse() * C;
          auto newp = RealPoint( P[ 0 ], P[ 1 ], P[ 2 ] );
          return ( newp - p ).squaredNorm() <= 1.0 ? newp : p;
        }
      }
    type = 2;
    if ( A.size() >= 2 )
      { // Trying to find two independent vectors
        RealVector D = RealVector::zero;
        for ( int i = 0; i < A.size(); i++ )
          for ( int j = i+1; j < A.size(); j++ )
            {
              D = A[ i ].crossProduct( A[ j ] );
              if ( D != RealVector::zero )
                { // we project p onto the line that is the intersection of the planes
                  RealTensor M;
                  RealVector C;
                  for ( int k = 0; k < 3; k++ ) {
                    M( 0, k ) = A[ i ][ k ];
                    M( 1, k ) = A[ j ][ k ];
                    M( 2, k ) = D[ k ];
                  }
                  C = RealVector( B[ 0 ], B[ 1 ], pointel.dot( D ) );
                  const auto d  = M.determinant();
                  if ( d == 0 ) std::cerr << "defficient rank" << std::endl;
                  else
                    {
                      const auto P = M.inverse() * C;
                      auto newp = RealPoint( P[ 0 ], P[ 1 ], P[ 2 ] );
                      return ( newp - p ).squaredNorm() <= 1.0 ? newp : p;
                    }
                }
            }
      }
    type = 1;
    // we project p onto the first plane, since all of them are equivalent or parallel.
    Scalar t = ( pointel.dot( A[ 0 ] ) - B[ 0 ] ) / A[0].squaredNorm();
    return ( pointel - t * A[ 0 ] );
  }
  
   
  std::pair< std::vector< RealVector >, std::vector< Scalar > >
  getDualPlanarities() const
  {
    const int nb = dualSurface->nbFaces();
    std::vector< RealVector > N( nb );
    std::vector< Scalar >     P( nb );
    for ( Index f = 0; f < nb; f++ )
      std::tie( N[ f ], P[ f ] ) = dualPlanarity( f );
    return std::make_pair( N, P );
  }
  
  std::pair< RealVector, Scalar > dualPlanarity( Index face ) const
  {
    auto V = dualSurface->verticesAroundFace( face );
    auto N = dual_face_normals[ face ];
    Scalar      best_error = std::numeric_limits<Scalar>::infinity();
    RealVector best_normal = RealVector::zero;
    Dimension      best_i  = 3;
    // Check all axis
    for ( Dimension i = 0; i < 3; i++ )
      {
        Dimension  j = (i+1)%3;
        Dimension  k = (i+2)%3;
        // Retrieve data
        std::vector< RealPoint > Input;
        for ( auto v : V ) Input.push_back( dual_positions[ v ] );
        RealTensor M;
        RealVector B;
        Scalar error = 0.0;
        for ( int v = 0; v < V.size(); v++ )
          {
            const RealVector& D = Input[ v ];
            M( i, i ) += 1.0;           // sum 1
            M( j, j ) += D[ j ]*D[ j ]; // sum y_i^2
            M( k, k ) += D[ k ]*D[ k ]; // sum z_i^2
            M( j, i ) -= D[ j ];        // sum -y_i
            M( i, j ) -= D[ j ];        // sum -y_i
            M( k, i ) -= D[ k ];        // sum -z_i
            M( i, k ) -= D[ k ];        // sum -z_i
            M( j, k ) += D[ j ]*D[ k ]; // sum y_i z_i
            M( k, j ) += D[ j ]*D[ k ]; // sum y_i z_i
            B[ i ]    += D[ i ];
            B[ j ]    += -D[ j ]*D[ i ];
            B[ k ]    += -D[ k ]*D[ i ];
            error     += D[ i ]*D[ i ];
          }
        if ( fabs( M.determinant() ) > 1e-8 ) {
          auto R = M.inverse() * B;
          RealVector rN = R;
          rN[ i ] = 1.0;
          if ( N.dot( rN ) < 0.0 ) rN *= -1;  
          error  += R.dot( M * R ) - 2.0 * R.dot( B );
          if ( error < best_error ) {
            best_i = i;
            best_error = error;
            best_normal = rN.getNormalized();
          }
        }
      }
    return best_error == std::numeric_limits<Scalar>::infinity()
      ? std::make_pair( N, -1.0 )
      : std::make_pair( best_normal, best_error );
  }

    
  std::vector<Index>
  segmentDualSurface( double epsilon = 0.00001 )
  {
    const Index Invalid = (Index) -1;
    std::vector< RealVector > N; //< normals
    std::vector< Scalar >     P; //< planarities
    std::tie( N, P ) = getDualPlanarities();
    std::vector<Index> labels( dualSurface->nbFaces(), Invalid );
    typedef std::pair< Scalar, Index > Node;
    std::priority_queue< Node, std::vector< Node >, std::greater< Node > > Q; //< queue
    std::queue< Index > QI; //< queue of invalid faces
    // Create queue
    for ( Index f = 0; f < dualSurface->nbFaces(); f++ )
      if ( P[ f ] == -1.0 ) QI.push( f );
      else Q.push( std::make_pair( P[ f ], f ) );
    std::cout << "Faces #valid=" << Q.size() << " #invalid=" << QI.size() << std::endl;
    // Process queue
    Index label = 0;
    while ( ! Q.empty() )
      {
        Index f = Q.top().second; Q.pop();
        if ( labels[ f ]  != Invalid ) continue; // already labelled
        std::queue< Index > localQ;
        std::set< Index >   M;
        localQ.push( f );
        M.insert( f );
        while ( ! localQ.empty() )
          {
            Index g = localQ.front(); localQ.pop();
            labels[ g ] = label;
            auto arcs = dualSurface->arcsAroundFace( g );
            for ( auto&& a : arcs )
              {
                auto ng = dualSurface->faceAroundArc( dualSurface->opposite( a ) );
                if ( ng == dualSurface->INVALID_FACE ) continue; // boundary
                if ( M.count( ng ) ) continue;          // already processed
                if ( ( P[ ng ] == -1 ) || ( ( N[ f ] - N[ ng ] ).norm() < epsilon ) )
                  {
                    localQ.push( ng );
                    M.insert( ng );
                  }
              }
          }
        label += 1;
      }
    return labels;
  }

};

  
} // namespace DGtal


//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined SymmetricConvexExpander_h

#endif // !defined SymmetricConvexExpander_RECURSES

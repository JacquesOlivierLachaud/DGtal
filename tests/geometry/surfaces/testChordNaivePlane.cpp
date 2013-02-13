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

/**
 * @file testChordNaivePlane.cpp
 * @ingroup Tests
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2012/03/05
 *
 * Functions for testing class ChordNaivePlane.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <cstdlib>
#include <iostream>
#include <QtGui/qapplication.h>
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/DrawWithDisplay3DModifier.h"
#include "DGtal/io/Color.h"
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/kernel/CPointPredicate.h"
#include "DGtal/geometry/surfaces/ChordNaivePlane.h"
#include "DGtal/geometry/surfaces/ChordGenericNaivePlane.h"
///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;

///////////////////////////////////////////////////////////////////////////////
// Functions for testing class ChordNaivePlane.
///////////////////////////////////////////////////////////////////////////////

template <typename Integer>
Integer getRandomInteger( Integer first, Integer after_last )
{
  Integer r = (Integer) random();
  return ( r % (after_last - first) ) + first;
}

/**
 * Checks the naive plane d <= ax+by+cz <= d + max(|a|,|b|,|c|)-1
 */
template <typename Integer, typename NaivePlane>
bool
checkPlane( Integer a, Integer b, Integer c, Integer d, 
            int diameter, unsigned int nbtries )
{
  typedef typename NaivePlane::Point Point;
  typedef typename Point::Component PointInteger;
  IntegerComputer<Integer> ic;
  Integer absA = ic.abs( a );
  Integer absB = ic.abs( b );
  Integer absC = ic.abs( c );
  Integer x, y, z;
  Dimension axis;
  if ( ( absA >= absB ) && ( absA >= absC ) )
    axis = 0;
  else if ( ( absB >= absA ) && ( absB >= absC ) )
    axis = 1;
  else
    axis = 2;
  Point p;
  NaivePlane plane;
  plane.init( axis, 1, 1 );
  // Checks that points within the naive plane are correctly recognized.
  unsigned int nb = 0;
  unsigned int nbok = 0;
  while ( nb != nbtries )
    {
      p[ 0 ] = getRandomInteger<PointInteger>( -diameter+1, diameter ); 
      p[ 1 ] = getRandomInteger<PointInteger>( -diameter+1, diameter ); 
      p[ 2 ] = getRandomInteger<PointInteger>( -diameter+1, diameter );
      x = (Integer) p[ 0 ];
      y = (Integer) p[ 1 ];
      z = (Integer) p[ 2 ];
      switch ( axis ) {
      case 0: p[ 0 ] = NumberTraits<Integer>::castToInt64_t( ic.ceilDiv( d - b * y - c * z, a ) ); break;
      case 1: p[ 1 ] = NumberTraits<Integer>::castToInt64_t( ic.ceilDiv( d - a * x - c * z, b ) ); break;
      case 2: p[ 2 ] = NumberTraits<Integer>::castToInt64_t( ic.ceilDiv( d - a * x - b * y, c ) ); break;
      } 
      bool ok_ext = plane.isExtendable( p ); // should be ok
      bool ok = plane.extend( p ); // should be ok
      ++nb, nbok += ok_ext ? 1 : 0;
      ++nb, nbok += ok ? 1 : 0;
      if ( ! ok )
        {
          std::cerr << "[ERROR] p=" << p << " NOT IN plane=" << plane << std::endl;
          for ( typename NaivePlane::ConstIterator it = plane.begin(), itE = plane.end();
                it != itE; ++it )
            std::cerr << " " << *it;
          std::cerr << endl;
          std::cerr << "d <= a*x+b*y+c*z <= d+max(a,b,c)"
                    << d << " <= " << a << "*" << p[0] 
                    << "+" << b << "*" << p[1] 
                    << "+" << c << "*" << p[2]
                    << " = " << (a*p[0]+b*p[1]+c*p[2])
                    << std::endl;
          break;
        }
      if ( ! ok_ext )
        {
          std::cerr << "[ERROR] p=" << p << " was NOT extendable IN plane=" << plane << std::endl;
          break;
        }
      // else
      //   std::cerr << "[OK] p=" << p << " IN plane=" << plane << std::endl;
    }

  // Checks that points outside the naive plane are correctly recognized as outliers.
  while ( nb != (nbtries * 11 ) / 10 )
    {
      p[ 0 ] = getRandomInteger<PointInteger>( -diameter+1, diameter ); 
      p[ 1 ] = getRandomInteger<PointInteger>( -diameter+1, diameter ); 
      p[ 2 ] = getRandomInteger<PointInteger>( -diameter+1, diameter );
      x = (Integer) p[ 0 ];
      y = (Integer) p[ 1 ];
      z = (Integer) p[ 2 ];
      switch ( axis ) {
      case 0: p[ 0 ] = NumberTraits<Integer>::castToInt64_t( ic.ceilDiv( d - b * y - c * z, a ) ); break;
      case 1: p[ 1 ] = NumberTraits<Integer>::castToInt64_t( ic.ceilDiv( d - a * x - c * z, b ) ); break;
      case 2: p[ 2 ] = NumberTraits<Integer>::castToInt64_t( ic.ceilDiv( d - a * x - b * y, c ) ); break;
      } 
      PointInteger tmp = getRandomInteger<PointInteger>( 2, 5 ) 
        * (2*getRandomInteger<PointInteger>( 0, 2 ) - 1 );
      p[ axis ] += tmp;
      bool ok_ext = ! plane.isExtendable( p ); // should *not* be ok
      bool ok = ! plane.extend( p ); // should *not* be ok
      ++nb, nbok += ok ? 1 : 0;
      ++nb, nbok += ok_ext ? 1 : 0;
      if ( ! ok )
        {
          std::cerr << "[ERROR] p=" << p << " IN plane=" << plane << std::endl;
          break;
        }
      if ( ! ok_ext )
        {
          std::cerr << "[ERROR] p=" << p << " was extendable IN plane=" << plane << std::endl;
          break;
        }
      // else
      //   std::cerr << "[OK] p=" << p << " IN plane=" << plane << std::endl;
    }
  return nb == nbok;
}

/**
 * Checks the naive plane d <= ax+by+cz <= d + max(|a|,|b|,|c|)-1
 */
template <typename Integer, typename NaivePlane>
bool
checkPlaneGroupExtension( Integer a, Integer b, Integer c, Integer d, 
                          int diameter, unsigned int nbtries )
{
  typedef typename NaivePlane::Point Point;
  typedef typename Point::Component PointInteger;
  IntegerComputer<Integer> ic;
  Integer absA = ic.abs( a );
  Integer absB = ic.abs( b );
  Integer absC = ic.abs( c );
  Integer x, y, z;
  Dimension axis;
  if ( ( absA >= absB ) && ( absA >= absC ) )
    axis = 0;
  else if ( ( absB >= absA ) && ( absB >= absC ) )
    axis = 1;
  else
    axis = 2;
  Point p;
  NaivePlane plane;
  plane.init( axis, 1, 1 );
  // Checks that points within the naive plane are correctly recognized.
  unsigned int nb = 0;
  unsigned int nbok = 0;
  while ( nb < nbtries )
    {
      std::vector<Point> points( 5 );
      for ( unsigned int i = 0; i < 5; ++i )
        {
          Point & p = points[ i ];
          p[ 0 ] = getRandomInteger<PointInteger>( -diameter+1, diameter ); 
          p[ 1 ] = getRandomInteger<PointInteger>( -diameter+1, diameter ); 
          p[ 2 ] = getRandomInteger<PointInteger>( -diameter+1, diameter );
          x = (Integer) p[ 0 ];
          y = (Integer) p[ 1 ];
          z = (Integer) p[ 2 ];
          switch ( axis ) {
          case 0: p[ 0 ] = NumberTraits<Integer>::castToInt64_t
              ( ic.ceilDiv( d - b * y - c * z, a ) ); break;
          case 1: p[ 1 ] = NumberTraits<Integer>::castToInt64_t
              ( ic.ceilDiv( d - a * x - c * z, b ) ); break;
          case 2: p[ 2 ] = NumberTraits<Integer>::castToInt64_t
              ( ic.ceilDiv( d - a * x - b * y, c ) ); break;
          } 
        }
      bool ok_ext = plane.isExtendable( points.begin(), points.end() ); // should be ok
      bool ok = plane.extend( points.begin(), points.end() ); // should be ok
      ++nb, nbok += ok_ext ? 1 : 0;
      ++nb, nbok += ok ? 1 : 0;
      if ( ! ok )
        {
          std::cerr << "[ERROR] p=" << points[ 0 ] << " NOT IN plane=" << plane << std::endl;
          for ( typename NaivePlane::ConstIterator it = plane.begin(), itE = plane.end();
                it != itE; ++it )
            std::cerr << " " << *it;
          std::cerr << endl;
          std::cerr << "d <= a*x+b*y+c*z <= d+max(a,b,c)"
                    << d << " <= " << a << "*" << p[0] 
                    << "+" << b << "*" << p[1] 
                    << "+" << c << "*" << p[2]
                    << " = " << (a*p[0]+b*p[1]+c*p[2])
                    << std::endl;
          break;
        }
      if ( ! ok_ext )
        {
          std::cerr << "[ERROR] p=" << p << " was NOT extendable IN plane=" << plane << std::endl;
          break;
        }
      // else
      //   std::cerr << "[OK] p=" << p << " IN plane=" << plane << std::endl;
    }

  // Checks that points outside the naive plane are correctly recognized as outliers.
  while ( nb < (nbtries * 11 ) / 10 )
    {
      p[ 0 ] = getRandomInteger<PointInteger>( -diameter+1, diameter ); 
      p[ 1 ] = getRandomInteger<PointInteger>( -diameter+1, diameter ); 
      p[ 2 ] = getRandomInteger<PointInteger>( -diameter+1, diameter );
      x = (Integer) p[ 0 ];
      y = (Integer) p[ 1 ];
      z = (Integer) p[ 2 ];
      switch ( axis ) {
      case 0: p[ 0 ] = NumberTraits<Integer>::castToInt64_t( ic.ceilDiv( d - b * y - c * z, a ) ); break;
      case 1: p[ 1 ] = NumberTraits<Integer>::castToInt64_t( ic.ceilDiv( d - a * x - c * z, b ) ); break;
      case 2: p[ 2 ] = NumberTraits<Integer>::castToInt64_t( ic.ceilDiv( d - a * x - b * y, c ) ); break;
      } 
      PointInteger tmp = getRandomInteger<PointInteger>( 2, 5 ) 
        * (2*getRandomInteger<PointInteger>( 0, 2 ) - 1 );
      p[ axis ] += tmp;
      bool ok_ext = ! plane.isExtendable( p ); // should *not* be ok
      bool ok = ! plane.extend( p ); // should *not* be ok
      ++nb, nbok += ok ? 1 : 0;
      ++nb, nbok += ok_ext ? 1 : 0;
      if ( ! ok )
        {
          std::cerr << "[ERROR] p=" << p << " IN plane=" << plane << std::endl;
          break;
        }
      if ( ! ok_ext )
        {
          std::cerr << "[ERROR] p=" << p << " was extendable IN plane=" << plane << std::endl;
          break;
        }
      // else
      //   std::cerr << "[OK] p=" << p << " IN plane=" << plane << std::endl;
    }
  return nb == nbok;
}



/**
 * Checks the naive plane d <= ax+by+cz <= d + max(|a|,|b|,|c|)-1
 */
template <typename Integer, typename GenericNaivePlane>
bool
checkGenericPlane( Integer a, Integer b, Integer c, Integer d, 
                   int diameter, unsigned int nbtries )
{
  typedef typename GenericNaivePlane::Point Point;
  typedef typename Point::Component PointInteger;
  IntegerComputer<Integer> ic;
  Integer absA = ic.abs( a );
  Integer absB = ic.abs( b );
  Integer absC = ic.abs( c );
  Integer x, y, z;
  Dimension axis;
  if ( ( absA >= absB ) && ( absA >= absC ) )
    axis = 0;
  else if ( ( absB >= absA ) && ( absB >= absC ) )
    axis = 1;
  else
    axis = 2;
  Point p;
  GenericNaivePlane plane;
  plane.init( 1, 1 );
  // Checks that points within the naive plane are correctly recognized.
  unsigned int nb = 0;
  unsigned int nbok = 0;
  while ( nb != nbtries )
    {
      p[ 0 ] = getRandomInteger<PointInteger>( -diameter+1, diameter ); 
      p[ 1 ] = getRandomInteger<PointInteger>( -diameter+1, diameter ); 
      p[ 2 ] = getRandomInteger<PointInteger>( -diameter+1, diameter );
      x = (Integer) p[ 0 ];
      y = (Integer) p[ 1 ];
      z = (Integer) p[ 2 ];
      switch ( axis ) {
      case 0: p[ 0 ] = NumberTraits<Integer>::castToInt64_t( ic.ceilDiv( d - b * y - c * z, a ) ); break;
      case 1: p[ 1 ] = NumberTraits<Integer>::castToInt64_t( ic.ceilDiv( d - a * x - c * z, b ) ); break;
      case 2: p[ 2 ] = NumberTraits<Integer>::castToInt64_t( ic.ceilDiv( d - a * x - b * y, c ) ); break;
      } 
      bool ok_ext = plane.isExtendable( p ); // should be ok
      bool ok = plane.extend( p ); // should be ok
      ++nb, nbok += ok_ext ? 1 : 0;
      ++nb, nbok += ok ? 1 : 0;
      if ( ! ok )
        {
          std::cerr << "[ERROR] p=" << p << " NOT IN plane=" << plane << std::endl;
          for ( typename GenericNaivePlane::ConstIterator it = plane.begin(), itE = plane.end();
                it != itE; ++it )
            std::cerr << " " << *it;
          std::cerr << endl;
          std::cerr << "d <= a*x+b*y+c*z <= d+max(a,b,c)"
                    << d << " <= " << a << "*" << p[0] 
                    << "+" << b << "*" << p[1] 
                    << "+" << c << "*" << p[2]
                    << " = " << (a*p[0]+b*p[1]+c*p[2])
                    << std::endl;
          break;
        }
      if ( ! ok_ext )
        {
          std::cerr << "[ERROR] p=" << p << " was NOT extendable IN plane=" << plane << std::endl;
          break;
        }
      // else
      //   std::cerr << "[OK] p=" << p << " IN plane=" << plane << std::endl;
    }

  // Checks that points outside the naive plane are correctly recognized as outliers.
  while ( nb != (nbtries * 11 ) / 10 )
    {
      p[ 0 ] = getRandomInteger<PointInteger>( -diameter+1, diameter ); 
      p[ 1 ] = getRandomInteger<PointInteger>( -diameter+1, diameter ); 
      p[ 2 ] = getRandomInteger<PointInteger>( -diameter+1, diameter );
      x = (Integer) p[ 0 ];
      y = (Integer) p[ 1 ];
      z = (Integer) p[ 2 ];
      switch ( axis ) {
      case 0: p[ 0 ] = NumberTraits<Integer>::castToInt64_t( ic.ceilDiv( d - b * y - c * z, a ) ); break;
      case 1: p[ 1 ] = NumberTraits<Integer>::castToInt64_t( ic.ceilDiv( d - a * x - c * z, b ) ); break;
      case 2: p[ 2 ] = NumberTraits<Integer>::castToInt64_t( ic.ceilDiv( d - a * x - b * y, c ) ); break;
      } 
      PointInteger tmp = getRandomInteger<PointInteger>( 2, 5 ) 
        * (2*getRandomInteger<PointInteger>( 0, 2 ) - 1 );
      p[ axis ] += tmp;
      bool ok_ext = ! plane.isExtendable( p ); // should *not* be ok
      bool ok = ! plane.extend( p ); // should *not* be ok
      ++nb, nbok += ok ? 1 : 0;
      ++nb, nbok += ok_ext ? 1 : 0;
      if ( ! ok )
        {
          std::cerr << "[ERROR] p=" << p << " IN plane=" << plane << std::endl;
          break;
        }
      if ( ! ok_ext )
        {
          std::cerr << "[ERROR] p=" << p << " was extendable IN plane=" << plane << std::endl;
          break;
        }
      // else
      //   std::cerr << "[OK] p=" << p << " IN plane=" << plane << std::endl;
    }
  std::cerr << "plane = " << plane << std::endl;
  return nb == nbok;
}


template <typename Integer, typename NaivePlane>
bool
checkPlanes( unsigned int nbplanes, int diameter, unsigned int nbtries )
{
  //using namespace Z3i;
  //typedef ChordNaivePlane<Z3, Integer> NaivePlane;
  unsigned int nb = 0;
  unsigned int nbok = 0;
  for ( unsigned int nbp = 0; nbp < nbplanes; ++nbp )
    {
      Integer a = getRandomInteger<Integer>( (Integer) 0, (Integer) diameter / 2 ); 
      Integer b = getRandomInteger<Integer>( (Integer) 0, (Integer) diameter / 2 ); 
      Integer c = getRandomInteger<Integer>( (Integer) 0, (Integer) diameter / 2 ); 
      Integer d = getRandomInteger<Integer>( (Integer) 0, (Integer) diameter / 2 ); 
      if ( ( a != 0 ) || ( b != 0 ) || ( c != 0 ) )
        {
          ++nb, nbok += checkPlane<Integer, NaivePlane>( a, b, c, d, diameter, nbtries ) ? 1 : 0;
          if ( nb != nbok )
            {
              std::cerr << "[ERROR] (Simple extension) for plane " << a << " * x + " 
                        << b << " * y + " << c << " * z = " << d << std::endl;
              break;
            }
          ++nb, nbok += checkPlaneGroupExtension<Integer, NaivePlane>( a, b, c, d, diameter, nbtries ) ? 1 : 0;
          if ( nb != nbok )
            {
              std::cerr << "[ERROR] (Group extension) for plane " << a << " * x + " 
                        << b << " * y + " << c << " * z = " << d << std::endl;
              break;
            }
        }
    }
  return nb == nbok;
}

/**
 * Checks the naive plane d <= ax+by+cz <= d + max(|a|,|b|,|c|)-1
 */
template <typename Integer, typename NaivePlane>
bool
checkWidth( Integer a, Integer b, Integer c, Integer d, 
            int diameter, unsigned int nbtries )
{
  typedef typename NaivePlane::Point Point;
  typedef typename NaivePlane::InternalScalar InternalScalar;
  IntegerComputer<Integer> ic;
  Integer absA = ic.abs( a );
  Integer absB = ic.abs( b );
  Integer absC = ic.abs( c );
  Integer x, y, z;
  Dimension axis;
  if ( ( absA >= absB ) && ( absA >= absC ) )
    axis = 0;
  else if ( ( absB >= absA ) && ( absB >= absC ) )
    axis = 1;
  else
    axis = 2;
  // Checks that points within the naive plane are correctly recognized.
  unsigned int nb = 0;
  unsigned int nbok = 0;
  std::vector<Point> points( nbtries );
  for ( unsigned int i = 0; i != nbtries; ++i )
    {
      Point & p = points[ i ];
      p[ 0 ] = getRandomInteger<Integer>( -diameter+1, diameter ); 
      p[ 1 ] = getRandomInteger<Integer>( -diameter+1, diameter ); 
      p[ 2 ] = getRandomInteger<Integer>( -diameter+1, diameter );
      x = (Integer) p[ 0 ];
      y = (Integer) p[ 1 ];
      z = (Integer) p[ 2 ];
      switch ( axis ) {
      case 0: p[ 0 ] = NumberTraits<Integer>::castToInt64_t( ic.ceilDiv( d - b * y - c * z, a ) ); break;
      case 1: p[ 1 ] = NumberTraits<Integer>::castToInt64_t( ic.ceilDiv( d - a * x - c * z, b ) ); break;
      case 2: p[ 2 ] = NumberTraits<Integer>::castToInt64_t( ic.ceilDiv( d - a * x - b * y, c ) ); break;
      } 
    }
  trace.beginBlock( "Computing axis width." );
  trace.info() << "- plane is " 
               << d << " <= " << a << "*x"
               << "+" << b << "*y"
               << "+" << c << "*z"
               << " <= d + max(|a|,|b|,|c|)"
               << std::endl;
  trace.info() << "- " << points.size() << " points tested in diameter " << diameter
               << std::endl;
  double min = -1.0;
  for ( unsigned int i = 0; i < 3; ++i )
    {
      std::pair<InternalScalar, InternalScalar> width 
        = NaivePlane::computeAxisWidth( i, points.begin(), points.end() );
      double wn = NumberTraits<InternalScalar>::castToDouble( width.first );
      double wd = NumberTraits<InternalScalar>::castToDouble( width.second );
      trace.info() << "  (" << i << ") width=" << (wn/wd) << std::endl;
      if ( min < 0.0 ) min = wn/wd;
      else if ( wn/wd < min ) min = wn/wd;
    }
  ++nb, nbok += (min < 1.0 ) ? 1 : 0;
  trace.info() << "(" << nbok << "/" << nb << ") min width = " << min
               << " < 1.0" << std::endl;
  ++nb, nbok += (0.9 < min ) ? 1 : 0;
  trace.info() << "(" << nbok << "/" << nb << ") min width = " << min
               << " > 0.9" << std::endl;
  trace.endBlock();
  return nb == nbok;
}

template <typename Integer, typename NaivePlane>
bool
checkWidths( unsigned int nbplanes, int diameter, unsigned int nbtries )
{
  //using namespace Z3i;
  //typedef ChordNaivePlane<Z3, Integer> NaivePlane;
  unsigned int nb = 0;
  unsigned int nbok = 0;
  for ( unsigned int nbp = 0; nbp < nbplanes; ++nbp )
    {
      Integer a = getRandomInteger<Integer>( (Integer) 0, (Integer) diameter / 2 ); 
      Integer b = getRandomInteger<Integer>( (Integer) 0, (Integer) diameter / 2 ); 
      Integer c = getRandomInteger<Integer>( (Integer) 0, (Integer) diameter / 2 ); 
      Integer d = getRandomInteger<Integer>( (Integer) 0, (Integer) diameter / 2 ); 
      if ( ( a != 0 ) || ( b != 0 ) || ( c != 0 ) )
        {
          ++nb, nbok += checkWidth<Integer, NaivePlane>( a, b, c, d, diameter, nbtries ) ? 1 : 0;
          if ( nb != nbok )
            {
              std::cerr << "[ERROR] (checkWidth) for plane " << a << " * x + " 
                        << b << " * y + " << c << " * z = " << d << std::endl;
              break;
            }
        }
    }
  return nb == nbok;
}


/**
 * Example of a test. To be completed.
 *
 */
bool testChordNaivePlane()
{
  unsigned int nbok = 0;
  unsigned int nb = 0;
  typedef DGtal::int64_t Integer;
  typedef DGtal::Z3i::Point Point;
  typedef ChordNaivePlane<Point, Integer> NaivePlane;
  typedef ChordGenericNaivePlane<Point, Integer> GenericNaivePlane;

  BOOST_CONCEPT_ASSERT(( CPointPredicate< NaivePlane > ));
  BOOST_CONCEPT_ASSERT(( boost::ForwardContainer< NaivePlane > ));
  BOOST_CONCEPT_ASSERT(( CPointPredicate< GenericNaivePlane > ));
  BOOST_CONCEPT_ASSERT(( boost::ForwardContainer< GenericNaivePlane > ));

  trace.beginBlock ( "Testing block: ChordNaivePlane instantiation." );
  NaivePlane plane;
  Point pt0( 0, 0, 0 );
  plane.init( 2, 1, 1 );
  bool pt0_inside = plane.extend( pt0 );
  trace.info() << "(" << nbok << "/" << nb << ") Plane=" << plane
               << std::endl;
  Point pt1( Point( 2, 0, 0 ) );
  bool pt1_inside = plane.extend( pt1 );
  ++nb, nbok += pt1_inside == true ? 1 : 0;
  trace.info() << "(" << nbok << "/" << nb << ") add " << pt1 
               << " Plane=" << plane << std::endl;
  Point pt2( Point( 0, 2, 2 ) );
  bool pt2_inside = plane.extend( pt2 );
  ++nb, nbok += pt2_inside == true ? 1 : 0;
  trace.info() << "(" << nbok << "/" << nb << ") add " << pt2 
               << " Plane=" << plane << std::endl;

  Point pt3( Point( 1, 1, 1 ) );
  bool pt3_inside = plane.extend( pt3 );
  ++nb, nbok += pt3_inside == true ? 1 : 0;
  trace.info() << "(" << nbok << "/" << nb << ") add " << pt3
               << " Plane=" << plane << std::endl;

  Point pt4( Point( -10, -10, 10 ) );
  bool pt4_inside = plane.extend( pt4 );
  ++nb, nbok += pt4_inside == false ? 1 : 0;
  trace.info() << "(" << nbok << "/" << nb << ") impossible add " << pt4
               << " Plane=" << plane << std::endl;

  Point pt5 = pt2 + Point( 1, 0, 1 );
  bool pt5_inside = plane.extend( pt5 );
  ++nb, nbok += pt5_inside == true ? 1 : 0;
  trace.info() << "(" << nbok << "/" << nb << ") add " << pt5
               << " Plane=" << plane << std::endl;

  Point pt6 = pt5 + Point( 6, 0, 2 );
  bool pt6_inside = plane.extend( pt6 );
  ++nb, nbok += pt6_inside == true ? 1 : 0;
  trace.info() << "(" << nbok << "/" << nb << ") add " << pt6
               << " Plane=" << plane << std::endl;

  NaivePlane plane2;
  plane2.init( 2, 1, 1 );
  plane2.extend( Point( 10, 0, 0 ) );
  plane2.extend( Point( 0, 8, 0 ) );
  plane2.extend( Point( 0, 0, 6 ) );
  trace.info() << "(" << nbok << "/" << nb << ") "
               << " Plane2=" << plane2 << std::endl;

  ++nb, nbok += checkPlane<Integer,NaivePlane>( 11, 5, 19, 20, 100, 100 ) ? 1 : 0;
  trace.info() << "(" << nbok << "/" << nb 
               << ") checkPlane<Integer,NaivePlane>( 11, 5, 19, 20, 100, 100 )"
               << std::endl;

  ++nb, nbok += checkGenericPlane<Integer,GenericNaivePlane>( 11, 5, 19, 20, 100, 100 ) ? 1 : 0;
  trace.info() << "(" << nbok << "/" << nb 
               << ") checkGenericPlane<Integer,GenericNaivePlane>( 11, 5, 19, 20, 100, 100 )"
               << std::endl;
  ++nb, nbok += checkGenericPlane<Integer,GenericNaivePlane>( 17, 33, 7, 10, 100, 100 ) ? 1 : 0;
  trace.info() << "(" << nbok << "/" << nb 
               << ") checkGenericPlane<Integer,GenericNaivePlane>( 17, 33, 7, 10, 100, 100 )"
               << std::endl;
  ++nb, nbok += checkPlane<Integer,NaivePlane>( 15, 8, 13, 15, 100, 100 ) ? 1 : 0;
  trace.info() << "(" << nbok << "/" << nb 
                << ") checkPlane<Integer,NaivePlane>( 15, 8, 13, 15, 100, 100 )"
                << std::endl;
  ++nb, nbok += checkGenericPlane<Integer,GenericNaivePlane>( 15, 8, 13, 15, 100, 100 ) ? 1 : 0;
  trace.info() << "(" << nbok << "/" << nb 
               << ") checkGenericPlane<Integer,GenericNaivePlane>( 15, 8, 13, 15, 100, 100 )"
               << std::endl;
  trace.endBlock();

  {
    trace.beginBlock ( "Testing block: ChordNaivePlane vertical instantiation." );
    NaivePlane plane;
    Point pt0( 0, 0, 0 );
    plane.init( 2, 5, 2 );
    bool pt0_inside = plane.extend( pt0 );
    ++nb, nbok += pt0_inside == true ? 1 : 0;
    trace.info() << "(" << nbok << "/" << nb << ") Plane=" << plane
                 << std::endl;
    Point pt1( 3, 2, 2 );
    bool pt1_inside = plane.extend( pt1 );
    ++nb, nbok += pt1_inside == true ? 1 : 0;
    trace.info() << "(" << nbok << "/" << nb << ") Plane=" << plane
                 << std::endl;
    Point pt2( 0, 0, 1 );
    bool pt2_inside = plane.extend( pt2 );
    ++nb, nbok += pt2_inside == true ? 1 : 0;
    trace.info() << "(" << nbok << "/" << nb << ") Plane=" << plane
                 << std::endl;
    Point pt3 = pt1 + Point( 0, 0, 1 );
    bool pt3_inside = plane.extend( pt3 );
    ++nb, nbok += pt3_inside == true ? 1 : 0;
    trace.info() << "(" << nbok << "/" << nb << ") Plane=" << plane
                 << std::endl;
    Point pt4 = pt3 + Point( 0, 0, 1 );
    bool pt4_inside = plane.extend( pt4 );
    ++nb, nbok += pt4_inside == true ? 1 : 0;
    trace.info() << "(" << nbok << "/" << nb << ") Plane=" << plane
                 << std::endl;
    trace.endBlock();
  }

  {
    trace.beginBlock ( "Testing block: ChordNaivePlane vertical instantiation 2." );
    NaivePlane plane;
    plane.init( 1, 1, 1 );
    Point pt0( -6, -3, 5 );
    bool pt0_inside = plane.extend( pt0 );
    ++nb, nbok += pt0_inside == true ? 1 : 0;
    trace.info() << "(" << nbok << "/" << nb << ") Plane=" << plane
                 << std::endl;
    Point pt1( 4, 4, -5 );
    bool pt1_inside = plane.extend( pt1 );
    ++nb, nbok += pt1_inside == true ? 1 : 0;
    trace.info() << "(" << nbok << "/" << nb << ") Plane=" << plane
                 << std::endl;
    Point pt2( -5, -2, 4 );
    bool pt2_inside = plane.extend( pt2 );
    ++nb, nbok += pt2_inside == true ? 1 : 0;
    trace.info() << "(" << nbok << "/" << nb << ") Plane=" << plane
                 << std::endl;
    trace.endBlock();
  }

  return nbok == nb;
}

template <typename NaivePlane>
bool 
checkManyPlanes( unsigned int diameter,
                 unsigned int nbplanes, 
                 unsigned int nbpoints )
{
  unsigned int nbok = 0;
  unsigned int nb = 0;
  typedef typename NaivePlane::InternalScalar Scalar;
  stringstream ss (stringstream::out);
  ss << "Testing block: Diameter is " << diameter << ". Check " << nbplanes << " planes with " << nbpoints << " points each.";
  trace.beginBlock ( ss.str() );
  ++nb, nbok += checkPlanes<Scalar,NaivePlane>( nbplanes, diameter, nbpoints ) ? 1 : 0;
  trace.info() << "(" << nbok << "/" << nb 
               << ") checkPlanes<Scalar,NaivePlane>()"
               << std::endl;
  trace.endBlock();
  return nbok == nb;
}

///////////////////////////////////////////////////////////////////////////////
// Standard services - public :

int main( int argc, char** argv )
{
  using namespace Z3i;

  // QApplication application(argc,argv);
  // Viewer3D viewer;
  // viewer.setWindowTitle("simpleViewer");
  // viewer.show();
  // Point pt0( -6, -3, 5 );
  // Point pt1( 4, 4, -5 );
  // Point pt2( -5, -2, 4 );
  // Domain domain( pt0, pt1 );
  // viewer << SetMode3D( pt0.className(), "Paving" );
  // viewer << domain;
  // viewer << CustomColors3D(Color(250, 200,0, 100),Color(250, 200,0, 50));
  // viewer << pt0 << pt1 << pt2;
  // viewer << Display3D::updateDisplay;
  // application.exec();

  // Max diameter is ~20 for int32_t, ~500 for int64_t, any with BigInteger.
  trace.beginBlock ( "Testing class ChordNaivePlane" );
  bool res = true 
    && testChordNaivePlane()
    && checkManyPlanes<ChordNaivePlane<Z3i::Point, DGtal::int32_t> >( 4, 100, 200 )
    && checkManyPlanes<ChordNaivePlane<Z3i::Point, DGtal::int32_t> >( 8, 100, 200 )
    && checkManyPlanes<ChordNaivePlane<Z3i::Point, DGtal::int32_t> >( 20, 100, 200 )
    && checkManyPlanes<ChordNaivePlane<Z3i::Point, DGtal::int32_t> >( 100, 100, 200 )
    && checkManyPlanes<ChordNaivePlane<Z3i::Point, DGtal::int64_t> >( 2000, 100, 200 )
    && checkWidths<DGtal::int64_t, ChordNaivePlane<Z3i::Point, DGtal::int64_t> >( 1000, 1000000, 1000 );
    // && checkManyPlanes<ChordNaivePlane<Z3, DGtal::int64_t> >( 500, 100, 200 )
    // && checkManyPlanes<ChordNaivePlane<Z3, DGtal::BigInteger> >( 10000, 10, 200 );
  trace.emphase() << ( res ? "Passed." : "Error." ) << endl;
  trace.endBlock();
  return res ? 0 : 1;
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

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
 * @file qglViewer.cpp
 * @ingroup Examples
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2011/01/04
 *
 * An example file named qglViewer.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <queue>
#include <QImageReader>
#include <QtGui/qapplication.h>
#include "DGtal/topology/Object.h"
#include "DGtal/shapes/Shapes.h"
#include "DGtal/helpers/StdDefs.h"

///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;

///////////////////////////////////////////////////////////////////////////////

template <typename Object, typename Map>
void 
generateSimplenessTable( const typename Object::DigitalTopology & dt,
			 Map & map )
{
  typedef typename Object::DigitalTopology DigitalTopology;
  typedef typename Object::DigitalSet DigitalSet;
  typedef typename Object::Point Point;
  typedef typename DigitalSet::Domain Domain;

  DigitalSet shapeSet;
  Object shape( dt, shape_set );
  Point p1 = Point::diagonal( -1 );
  Point p2 = Point::diagonal(  1 );
  Point c = Point::diagonal( 0 );
  Domain domain( p1, p2 );
  unsigned int k = 0;
  for ( Domain::ConstIterator it = domain.begin(); it != domain.end(); ++it )
    if ( *it != c ) ++k;
  ASSERT( ( k < 32 )
	  && "[generateSimplenessTable] number of configurations is too high." );
  unsigned int nbCfg = 1 << k;
  for ( unsigned int cfg = 0; cfg < nbCfg; ++cfg )
    {
      shape.pointSet().clear();
      shape.pointSet().insert( c );
      unsigned int mask = 1;
      for ( Domain::ConstIterator it = domain.begin(); it != domain.end(); ++it )
	{
	  if ( cfg & mask ) shape.pointSet().insert( *it );
	  mask <<= 1;
	}
      bool simple = shape.isSimple( c );
      std::cerr << "Simple[ " << cfg << " ] = " << simple << std::endl;
      map[ cfg ] = simple;
    }
}


int main( int argc, char** argv )
{

  using namespace Z2i;
  trace.beginBlock ( "Generate 2d table for 4-8 topology" );
  std::vector<bool> table4_8( 256 );
  generateSimplenessTable< Object4_8 >( dt4_8, table4_8 );
  trace.endBlock();

  trace.beginBlock ( "Generate 2d table for 8-4 topology" );
  std::vector<bool> table8-4( 256 );
  generateSimplenessTable< Object8_4 >( dt8_4, table8_4 );
  trace.endBlock();

  return 0;
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


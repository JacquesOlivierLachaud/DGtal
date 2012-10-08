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
 * @file logoDGtal.cpp
 * @ingroup Examples
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2012/03/05
 * 
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <QtGui/qapplication.h>
#include "DGtal/base/Common.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/helpers/StdDefs.h"
///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;
using namespace DGtal::Z3i;

///////////////////////////////////////////////////////////////////////////////

int main( int argc, char** argv )
{
  QApplication application(argc,argv);
  trace.beginBlock ( "Generate DGtal logo (without drop shadow)" );

  //! [logoDGtal-main]
  uint8_t dg[ 7 ][ 12 ] = { 
    { 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0 },
    { 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0 },
    { 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0 },
    { 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1 },
    { 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1 },
    { 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1 },
    { 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0 } };
  uint8_t tal[ 7 ][ 12 ] = { 
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 },
    { 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 },
    { 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0 },
    { 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0 },
    { 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0 } };
 
  Point p1( 0, 0, 0 );
  Point p2( 12, 14, 7 );
  Domain domain( p1, p2 );
  Color greyblue1( 190, 190, 210 );
  Color greyblue2( 160, 160, 180 );
  Color greyblue3( 100, 100, 120 );
  Color greyblue4( 60, 60, 80 );
  Color green1( 60, 255, 80 );
 
  Viewer3D viewer;
  viewer.show();
  viewer << CustomColors3D( greyblue2, greyblue2 );
  for ( unsigned int x = 0; x < 12; ++x )
    for ( unsigned int y = 0; y < 7; ++y )
      if ( dg[ y ][ x ] != 0 ) 
        viewer << Point( 0, 11-x, 6-y ) << Point( 1, 11-x, 6-y );
  viewer << CustomColors3D( greyblue3, greyblue3 );
  for ( unsigned int x = 0; x < 12; ++x )
    for ( unsigned int y = 0; y < 7; ++y )
      if ( tal[ y ][ x ] != 0 ) 
        viewer << Point( x+2, 0, 6-y ) << Point( x+2, 1, 6-y );
  KSpace K;
  K.init( p1, p2, true );
  Cell surfel = K.uCell( Point( 1, 1, 0 ) );
  for ( unsigned int x = 0; x < 13; ++x )
    for ( unsigned int y = 0; y < 12; ++y )
      {
        Color greyblue4( 60, 60, 80, 255 - 5*x - 5*y );
        Color green1( 60, 255, 80, 255 - 5*x - 5*y );
        viewer << ( ( (x+y) % 2 ) == 0 
                    ? CustomColors3D( greyblue4, greyblue4 )
                    : CustomColors3D( green1, green1 ) );
        viewer << K.uCell( Point( x, y, 0 ), surfel );
      }
  surfel = K.uCell( Point( 1, 0, 1 ) );
  for ( unsigned int x = 0; x < 13; ++x )
    for ( unsigned int z = 0; z < 7; ++z )
      {
        Color greyblue4( 60, 60, 80, 255 - 20*z - 10*x );
        Color green1( 60, 255, 80, 255 - 20*z - 10*x );
        viewer << ( ( (x+z) % 2 ) == 0 
                    ? CustomColors3D( greyblue4, greyblue4 )
                    : CustomColors3D( green1, green1 ) );
        viewer << K.uCell( Point( x, 12, z ), surfel );
      }
  surfel = K.uCell( Point( 0, 1, 1 ) );
  for ( unsigned int y = 0; y < 12; ++y )
    for ( unsigned int z = 0; z < 7; ++z )
      {
        Color greyblue4( 60, 60, 80, 255 - 20*z - 10*y );
        Color green1( 60, 255, 80, 255 - 20*z - 10*y );
        viewer << ( ( (y+z) % 2 ) == 0 
                    ? CustomColors3D( greyblue4, greyblue4 )
                    : CustomColors3D( green1, green1 ) );
        viewer << K.uCell( Point( 13, y, z ), surfel );
      }
  viewer << Viewer3D::updateDisplay;

  //! [logoDGtal-main]
  trace.endBlock();
  return application.exec();
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

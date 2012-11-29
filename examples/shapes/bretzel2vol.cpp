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
 * @file bretzel2vol.cpp
 * @ingroup Examples
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5807), University of Savoie, France
 *
 * @date 2012/02/12
 *
 * Tools for creating a 3D bretzel in vol file. More generally, it shows how to build a volume from a 3D parametric curve with some thickness.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <cmath>
#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/shapes/Shapes.h"
#include "DGtal/images/imagesSetsUtils/ImageFromSet.h"
#include "DGtal/io/writers/VolWriter.h"

///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;

typedef Z3i::Space Space3D;
typedef Z3i::Point Point3D;
typedef Space3D::RealPoint RealPoint3D;
typedef RealPoint3D::Component Component;

// in [-1:1]^2
RealPoint3D parametricBretzel( Component t )
{
  // plot [-0.35:6.63] 0.75*sin(t), sin(1.5*t)
  return RealPoint3D( sin( t ), sin( 1.5*t ), 0.1 * cos( t ) + 0.05 * sin( t/2.0 ) );
}

Point3D digitize( RealPoint3D p, Component h )
{
  return Point3D( p[ 0 ] / h, p[ 1 ] / h, p[ 2 ] / h );
}

///////////////////////////////////////////////////////////////////////////////
// Standard services - public :

int main( int argc, char** argv )
{
  typedef Z3i::Domain Domain3D;
  typedef Z3i::DigitalSet DigitalSet3D;
  //Default image selector = STLVector
  typedef ImageSelector<Domain3D, unsigned char>::Type Image3D;

  if ( argc == 1 )
    {
      std::cout << "Usage: bretzel2vol <N> [<E>] [<file.vol>]" << std::endl
                << "       Creates a bretzel in .vol 3D file." << std::endl
                << "       - N : resolution for the bretzel (=> =~ 2Nx2Nx2E bretzel)" << std::endl
                << "       - E : inner radius or thickness of the bretzel." << std::endl;
      return 0;
    }
  unsigned int N = ( argc >= 2 ) ? atoi( argv[ 1 ] ) : 50;
  unsigned int E = ( argc >= 3 ) ? atoi( argv[ 2 ] ) : 6;
  std::string filename =( argc >= 4 ) ? string( argv[ 3 ] ) : "bretzel.vol";
  double h = 1.0 / (double) N;
  Point3D low( -N-E-2, -N-E-2, -E-(N+1)/7 );
  Point3D up( N+E+2, N+E+2, E+(N+1)/7 );
  Domain3D domain( low, up );

  DigitalSet3D aSet( domain );

  for ( Component t = -0.31; t < 6.59; t += h )
    {
      Point3D p2 = digitize( parametricBretzel( t ), h );
      Shapes<Domain3D>::addNorm2Ball( aSet, p2, E );
    }
  Image3D vol = ImageFromSet<Image3D>::create
    ( aSet, 1, 1, aSet.begin(), aSet.end() );
  VolWriter<Image3D>::exportVol( filename.c_str(), vol);
 
  return 0;
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

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
 * @file convex-and-concave-parts.cpp
 * @ingroup Examples
 * @author Tristan Roussillon (\c tristan.roussillon@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 * @date 2011/01/24
 *
 * An example file named convex-and-concave-parts.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include "DGtal/base/Common.h"
#include "DGtal/io/boards/Board2D.h"
#include "DGtal/io/Color.h"
#include "DGtal/shapes/Shapes.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/geometry/curves/ArithmeticalDSS.h"
#include "DGtal/geometry/curves/FreemanChain.h"
#include "DGtal/geometry/curves/SaturatedSegmentation.h"
///////////////////////////////////////////////////////////////////////////////
#include "ConfigExamples.h"

using namespace std;
using namespace DGtal;
using namespace Z2i;

///////////////////////////////////////////////////////////////////////////////

typedef FreemanChain<int> Range; 

///////////////////////////////////////////////////////////////////////////////
enum GeometricType {
  Unknown, Convex, Concave, ConvexToConcave, ConcaveToConvex
};

/**
 * Drawing a segmentation
 */
template <typename Iterator, typename Board>
void drawCCP(const Iterator& itb, const Iterator& ite, Board& aBoard)
{
  
  typedef typename Iterator::SegmentComputer::ConstIterator PointIterator; 

  aBoard << SetMode( "ArithmeticalDSS", "BoundingBox" );
  string aStyleName = "ArithmeticalDSS/BoundingBox";

  typedef std::vector< std::pair< PointIterator, GeometricType > > GeometryContainer;
  typedef typename GeometryContainer::const_iterator GeometryConstIterator;
  GeometryContainer geometry;
  PointIterator last;
  for (Iterator i = itb; i != ite; ++i) {
     
    //choose pen color
    CustomPenColor* aPenColor;

    if ( !(i.intersectNext() && i.intersectPrevious()) ) {

      aPenColor = new CustomPenColor( Color::Black );

    } else {

      //end points

      PointIterator begin = i->begin();  --begin; 
      PointIterator end = i->end();

      //parameters
      int mu = i->getMu();
      int omega = i->getOmega();

      //configurations
      PointIterator itp;
      if ( (i->getRemainder(begin)<=mu-1)&&
           (i->getRemainder(end)<=mu-1) ) {                //concave
        aPenColor = new CustomPenColor( Color::Green);
        for ( itp = i->begin(); itp != i->end(); ++itp )
          {
            if ( i->getRemainder( itp ) == mu+omega-1 )
              {
                geometry.push_back( std::make_pair( itp, Concave ) );
                break;
              }
          }
      } else if ( (i->getRemainder(begin)>=mu+omega)&&
            (i->getRemainder(end)>=mu+omega) ) {           //convex
        aPenColor = new CustomPenColor( Color::Blue );
        for ( itp = i->begin(); itp != i->end(); ++itp )
          {
            if ( i->getRemainder( itp ) == mu )
              {
                geometry.push_back( std::make_pair( itp, Convex ) );
                break;
              }
          }
      } else if ( (i->getRemainder(begin)>=mu+omega)&&
            (i->getRemainder(end)<=mu-1) ) {               //convex to concave
        aPenColor = new CustomPenColor( Color::Yellow );
        for ( itp = i->begin(); itp != i->end(); ++itp )
          {
            if ( i->getRemainder( itp ) == mu )
              {
                geometry.push_back( std::make_pair( itp, ConvexToConcave ) );
                break;
              }
          }
        for ( ; itp != i->end(); ++itp )
          {
            if ( i->getRemainder( itp ) == mu+omega-1 )
              {
                geometry.push_back( std::make_pair( itp, Concave ) );
                break;
              }
          }
      } else if ( (i->getRemainder(begin)<=mu-1)&&
            (i->getRemainder(end)>=mu+omega) ) {           //concave to convex
        aPenColor = new CustomPenColor( Color::Yellow );
        for ( itp = i->begin(); itp != i->end(); ++itp )
          {
            if ( i->getRemainder( itp ) == mu+omega-1 )
              {
                geometry.push_back( std::make_pair( itp, ConcaveToConvex ) );
                break;
              }
          }
        for ( ; itp != i->end(); ++itp )
          {
            if ( i->getRemainder( itp ) == mu )
              {
                geometry.push_back( std::make_pair( itp, Convex ) );
                break;
              }
          }
      } else {                                                    //pb
        aPenColor = new CustomPenColor( Color::Red );
      }

    }
    trace.info() << "MS" << std::endl;
    // draw each segment
    aBoard << CustomStyle( aStyleName, aPenColor )
           << *i; 
    Iterator inext = boost::next( i );
    if ( inext == ite ) last = i->end();
  } 

  // draw convexity for vertices.
  PointIterator itp = itb->begin();
  string className = (*itp).className();
  string ptStyleName = className + "/Grid";
  aBoard << SetMode( className, "Grid" );
  GeometryContainer allGeometry;
  GeometricType tfill = Unknown;
  for ( GeometryConstIterator itGeom = geometry.begin();
        itGeom != geometry.end(); ++itGeom )
    {
      for ( ; ( itp != itGeom->first ) && ( itp != last ); ++itp )
        {
          if ( allGeometry.empty() || ( allGeometry.back().first != itp ) )
            allGeometry.push_back( std::make_pair( itp, tfill ) );
          else //if ( ( tfill == ConvexToConcave ) || ( tfill == ConcaveToConvex ) )
            allGeometry.back() = std::make_pair( itp, tfill );
          trace.info() << "p=" << *itp << " obj=" << *(itGeom->first);
          if ( tfill == Convex ) trace.info() << " Convex.";
          if ( tfill == Concave ) trace.info() << " Concave.";
          if ( tfill == ConvexToConcave ) trace.info() << " ConvexToConcave.";
          if ( tfill == ConcaveToConvex ) trace.info() << " ConcaveToConvex.";
          trace.info() << std::endl;
        }
      tfill = itGeom->second;
    }
  for ( GeometryConstIterator itGeom = allGeometry.begin(), itGeomEnd = allGeometry.end(); 
        itGeom != itGeomEnd; ++itGeom ) 
    {
      //choose pen color
      CustomColors* aPenColor;
      if ( itGeom->second == Concave )
        aPenColor = new CustomColors( Color::Green, Color::Green );
      else if ( itGeom->second == Convex )
        aPenColor = new CustomColors( Color::Blue, Color::Blue );
      else if ( itGeom->second == ConvexToConcave )
        aPenColor = new CustomColors( Color::Yellow, Color::Yellow );
      else if ( itGeom->second == ConcaveToConvex )
        aPenColor = new CustomColors( Color::Red, Color::Red );
      else
        aPenColor = new CustomColors( Color::Black, Color::Black );
      aBoard << CustomStyle( ptStyleName, aPenColor ) 
             << *(itGeom->first);
    }
}

/**
 * saturated segmentation of a range
 */
template <typename Iterator, typename Board>
void segmentationIntoMaximalDSSs(const Iterator& itb, const Iterator& ite, 
                                 Board& aBoard)
{
  typedef typename IteratorCirculatorTraits<Iterator>::Value::Coordinate Coordinate; 
  typedef ArithmeticalDSS<Iterator,Coordinate,4> RecognitionAlgorithm;
  typedef SaturatedSegmentation<RecognitionAlgorithm> Segmentation;

  RecognitionAlgorithm algo;
  Segmentation s(itb,ite,algo);
  
  typename Segmentation::SegmentComputerIterator i = s.begin();
  typename Segmentation::SegmentComputerIterator end = s.end();

  drawCCP<typename Segmentation::SegmentComputerIterator, Board>
  (i,end,aBoard); 

}


int main( int argc, char** argv )
{

  trace.beginBlock ( "Example convex-and-concave-parts" );

  trace.info() << "Args:";
  for ( int i = 0; i < argc; ++i )
    trace.info() << " " << argv[ i ];
  trace.info() << endl;


  string codes; 
  if (argc >= 2) codes = argv[1];
  else codes = "030030330303303030300001010101101011010000030330303303030300001010110101011010000033"; 

  stringstream ss(stringstream::in | stringstream::out);
  ss << "0 0 " << codes << endl;
  Range theContour( ss );
  
  trace.info() << "Processing of " << ss.str() << endl;

  //Maximal Segments
  Board2D aBoard;
  aBoard
   << SetMode( "PointVector", "Grid" )
   << theContour;

  segmentationIntoMaximalDSSs(theContour.begin(), theContour.end(), aBoard);

  aBoard.saveSVG("convex-and-concave-parts-2.svg");
  aBoard.saveEPS("convex-and-concave-parts-2.eps");

  trace.endBlock();

  return 0;
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

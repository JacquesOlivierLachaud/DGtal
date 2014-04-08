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
 * @file testDynDichoSet.cpp
 * @ingroup Tests
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2014/04/08
 *
 * Functions for testing class DynDichoSet.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <set>
#include "DGtal/base/Common.h"
#include "DGtal/base/DynDichoSet.h"
///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;

///////////////////////////////////////////////////////////////////////////////
// Functions for testing class DynDichoSet.
///////////////////////////////////////////////////////////////////////////////
/**
 * Example of a test. To be completed.
 *
 */
template <typename Set>
bool testDynDichoSet( const std::string & name, int N, int R )
{
  unsigned int nbok = 0;
  unsigned int nb = 0;
  std::string str = "Testing set " + name; 
  trace.beginBlock ( str.c_str() );
  trace.info() << "- " << N << " insertions of random values in 0,..., " 
               << (R-1) << std::endl; 
  Set aSet;
  for ( int i = 0; i < N; ++i )
    aSet.insert( rand() % R );
  trace.info() << "- set size = " << aSet.size() << std::endl;
  trace.endBlock();

  // nbok += true ? 1 : 0; 
  // nb++;
  // trace.info() << "(" << nbok << "/" << nb << ") "
  //              << "true == true" << std::endl;
  
  return nbok == nb;
}

///////////////////////////////////////////////////////////////////////////////
// Standard services - public :

int main( int argc, char** argv )
{
  trace.beginBlock ( "Testing class DynDichoSet" );
  trace.info() << "Args:";
  for ( int i = 0; i < argc; ++i )
    trace.info() << " " << argv[ i ];
  trace.info() << endl;
  typedef std::set<int> StdSet;
  typedef DGtal::DynDichoSet<int> DDSet;
  bool res = 
    testDynDichoSet<StdSet>( "std::set<int>", 100000, 1000000 )
    && testDynDichoSet<DDSet>( "DGtal::DynDichoSet<int>", 100000, 1000000 )
    && testDynDichoSet<StdSet>( "std::set<int>", 1000000, 1000000 )
    && testDynDichoSet<DDSet>( "DGtal::DynDichoSet<int>", 1000000, 1000000 ); 
  trace.emphase() << ( res ? "Passed." : "Error." ) << endl;
  trace.endBlock();
  return res ? 0 : 1;
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

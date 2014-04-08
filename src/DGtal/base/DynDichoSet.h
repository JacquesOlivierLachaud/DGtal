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
 * @file DynDichoSet.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2014/04/08
 *
 * Header file for module DynDichoSet.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(DynDichoSet_RECURSES)
#error Recursive header files inclusion detected in DynDichoSet.h
#else // defined(DynDichoSet_RECURSES)
/** Prevents recursive inclusion of headers. */
#define DynDichoSet_RECURSES

#if !defined DynDichoSet_h
/** Prevents repeated inclusion of headers. */
#define DynDichoSet_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <vector>
#include "DGtal/base/Common.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // template class DynDichoSet
  /**
   * Description of template class 'DynDichoSet' <p>
   * \brief Aim:
   */
  template <typename T>
  class DynDichoSet
  {
  public:
    typedef DynDichoSet<T> Self;
    typedef std::vector<T> Container;
    typedef typename Container::value_type value_type;
    typedef typename Container::size_type size_type;
    typedef typename Container::iterator iterator;              //< To be modified.
    typedef typename Container::const_iterator const_iterator;  //< To be modified.
    typedef const_iterator ConstIterator;
    typedef iterator Iterator;
    typedef size_type Size;

    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Destructor.
     */
    ~DynDichoSet();

    /**
     * Constructor.
     */
    DynDichoSet();

    /**
     * Copy constructor.
     * @param other the object to clone.
     */
    DynDichoSet ( const DynDichoSet & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     */
    DynDichoSet & operator= ( const DynDichoSet & other );

    Size size() const;
    ConstIterator end() const;
    void insert( const T& val );
    ConstIterator find( const T& val ) const;

    // ----------------------- Interface --------------------------------------
  public:

    /**
     * Writes/Displays the object on an output stream.
     * @param out the output stream where the object is written.
     */
    void selfDisplay ( std::ostream & out ) const;

    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool isValid() const;

    // ------------------------- Protected Datas ------------------------------
  private:
    // ------------------------- Private Datas --------------------------------
  private:
    void mergeInplace( Size j );
    void mergeAfter( Size j );

    // ------------------------- Hidden services ------------------------------
  protected:
    ConstIterator dichotomy( Size i, Size j, const T& v ) const;

    // ------------------------- Internals ------------------------------------
  private:
    Size myNb;
    std::vector<T> myData;
  }; // end of class DynDichoSet


  /**
   * Overloads 'operator<<' for displaying objects of class 'DynDichoSet'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'DynDichoSet' to write.
   * @return the output stream after the writing.
   */
  template <typename T>
  std::ostream&
  operator<< ( std::ostream & out, const DynDichoSet<T> & object );

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "DGtal/base/DynDichoSet.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined DynDichoSet_h

#undef DynDichoSet_RECURSES
#endif // else defined(DynDichoSet_RECURSES)

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
 * @file ChordGenericNaivePlane.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2012/09/20
 *
 * Header file for module ChordGenericNaivePlane.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(ChordGenericNaivePlane_RECURSES)
#error Recursive header files inclusion detected in ChordGenericNaivePlane.h
#else // defined(ChordGenericNaivePlane_RECURSES)
/** Prevents recursive inclusion of headers. */
#define ChordGenericNaivePlane_RECURSES

#if !defined ChordGenericNaivePlane_h
/** Prevents repeated inclusion of headers. */
#define ChordGenericNaivePlane_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <set>
#include "DGtal/base/Common.h"
#include "DGtal/kernel/CInteger.h"
#include "DGtal/kernel/CSpace.h"
#include "DGtal/kernel/SpaceND.h"
#include "DGtal/kernel/PointVector.h"
#include "DGtal/arithmetic/IntegerComputer.h"
#include "DGtal/geometry/surfaces/ChordNaivePlane.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // template class ChordGenericNaivePlane
  /**
   * Description of template class 'ChordGenericNaivePlane' <p> \brief
   * Aim: A class that recognizes pieces of digital planes of given
   * axis width. When the width is 1, it corresponds to naive
   * planes. Contrary to ChordNaivePlane, the axis is \b not specified
   * at initialization of the object. This class uses three instances
   * of ChordNaivePlane, one per axis.
   *
   * As a (3D) geometric primitive, it obeys to a subset of the
   * concept CSegmentComputer. It is copy constructible,
   * assignable. It is iterable (inner type ConstIterator, begin(),
   * end()). You may clear() it. It has methods \ref extend(), extend(
   * InputIterator, InputIterator) and \ref isExtendable(),
   * isExtendable(InputIterator, InputIterator).  The object stores
   * all the distinct points \c p such that 'extend( \c p )' was
   * successful. It is thus a model of boost::ForwardContainer (non
   * mutable).
   *
   * It is also a model of CPointPredicate (returns 'true' iff a point
   * is within the current bounds).
   *
   * Note on complexity: See ChordNaivePlane. Although it uses three
   * instances of ChordNaivePlane, the recognition is \b not three
   * times slower. Indeed, recognition stops quickly on bad axes.
   *
   * @tparam TPoint specifies the type of input points (digital or not). 
   *
   * @tparam TInternalScalar specifies the type of scalar used in
   * internal computations, generally a more precise type than
   * TPoint::Component. For instance, for digital points, the type
   * should be able to hold integers of order (2*D^3) if D is the
   * diameter of the set of digital points.
   *
   @code
   typedef SpaceND<3,int> Z3;
   typedef ChordGenericNaivePlane< Z3::Point, int64_t > NaivePlane;
   NaivePlane plane;
   plane.init( 1, 1 ); // width is 1/1 => naive 
   plane.extend( Point( 10, 0, 0 ) ); // return 'true'
   plane.extend( Point( 0, 8, 0 ) );  // return 'true'
   plane.extend( Point( 0, 0, 6 ) );  // return 'true'
   plane.extend( Point( 5, 5, 5 ) );  // return 'false'
   // There is no naive plane going through the 3 first points and the last one.
   @endcode
   *
   * Model of boost::DefaultConstructible, boost::CopyConstructible,
   * boost::Assignable, boost::ForwardContainer, CPointPredicate.
   */
  template < typename TPoint,
             typename TInternalScalar >
  class ChordGenericNaivePlane
  {
    // BOOST_CONCEPT_ASSERT(( CPoint< TPoint > ));
    BOOST_CONCEPT_ASSERT(( CSignedNumber< TInternalScalar > ));
    BOOST_STATIC_ASSERT(( TPoint::dimension == 3 ));

    // ----------------------- public types ------------------------------
  public:
    typedef TPoint Point;
    typedef TInternalScalar InternalScalar;
    typedef Point Vector;
    typedef typename Vector::Component Component;
    typedef typename Point::Coordinate Coordinate;
    typedef InternalScalar InternalVector[ 3 ];

    typedef std::set< Point > PointSet;
    typedef typename PointSet::size_type Size;
    typedef typename PointSet::const_iterator ConstIterator;
    typedef typename PointSet::iterator Iterator;

    // ----------------------- std public types ------------------------------
  public:
    typedef typename PointSet::const_iterator const_iterator;
    typedef typename PointSet::const_pointer const_pointer;
    typedef typename PointSet::const_reference const_reference;
    typedef typename PointSet::value_type value_type;
    typedef typename PointSet::difference_type difference_type;
    typedef typename PointSet::size_type size_type;

    // ----------------------- internal types ------------------------------
  private:
    typedef ChordNaivePlane< Point, InternalScalar > ChordComputer;
    typedef std::vector<Dimension>::iterator AxisIterator;
    typedef std::vector<Dimension>::const_iterator AxisConstIterator;
    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Destructor.
     */
    ~ChordGenericNaivePlane();

    /**
     * Constructor. The object is not valid and should be initialized.
     * @see init
     */
    ChordGenericNaivePlane();

    /**
     * Copy constructor.
     * @param other the object to clone.
     */
    ChordGenericNaivePlane ( const ChordGenericNaivePlane & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     */
    ChordGenericNaivePlane & operator= ( const ChordGenericNaivePlane & other );

    /**
       @return an active axis (or the active axis when there is only one).
    */
    Dimension active() const;

    /**
     * Clear the object, free memory. The plane keeps its main axis,
     * diameter and width, but contains no point.
     */
    void clear();

    /**
     * All these parameters cannot be changed during the process.
     * After this call, the object is in a consistent state and can
     * accept new points for recognition. Calls clear so that the
     * object is ready to be extended.
     *
     * @param widthNumerator the maximal axis-width (x,y,or z) for the
     * plane is defined as the rational number \a widthNumerator / \a
     * widthDenominator (default is 1/1, i.e. naive plane).
     *
     * @param widthDenominator the maximal axis-width (x,y,or z) for
     * the plane is defined as the rational number \a widthNumerator /
     * \a widthDenominator (default is 1/1, i.e. naive plane).
     */
    void init( InternalScalar widthNumerator = NumberTraits< InternalScalar >::ONE, 
               InternalScalar widthDenominator = NumberTraits< InternalScalar >::ONE );

    //-------------------- model of ForwardContainer -----------------------------
  public:

    /**
     * @return the number of distinct points in the current naive plane.
     */
    Size size() const;

    /**
     * @return 'true' if and only if this object contains no point.
     */
    bool empty() const;

    /**
     * @return a const iterator pointing on the first point stored in the current naive plane.
     */
    ConstIterator begin() const;

    /**
     * @return a const iterator pointing after the last point stored in the current naive plane.
     */
    ConstIterator end() const;

    /**
     * NB: std version.
     * @return the maximal allowed number of points in the current naive plane.
     * @see maxSize
     */
    Size max_size() const;

    /**
     * same as max_size
     * @return the maximal allowed number of points in the current naive plane.
     */
    Size maxSize() const;


    //-------------------- model of CPointPredicate -----------------------------
  public:

    /**
     * Checks if the point \a p is in the current digital
     * plane. Therefore, a ChordGenericNaivePlane is a model of
     * CPointPredicate.
     *
     * @param p any 3D point.
     *
     * @return 'true' if it is in the current plane, false otherwise.
     */
    bool operator()( const Point & p ) const;

    //-------------------- model of CIncrementalPrimitiveComputer -----------------------------
  public:

    /**
     * Adds the point \a p to this plane if it is within the current
     * bounds. The plane parameters are not updated.
     *
     * @param p any 3D point (in the specified diameter).
     *
     * @return 'true' if \a p is in the plane, 'false' otherwise (the
     * object is then in its original state).
     */
    bool extendAsIs( const Point & p );

    /**
     * Adds the point \a p and checks if we have still a digital plane
     * of specified width. The plane parameters may be updated so as
     * to include the new point.
     *
     * @param p any 3D point (in the specified diameter).
     *
     * @return 'true' if it is still a plane, 'false' otherwise (the
     * object is then in its original state).
     */
    bool extend( const Point & p );

    /**
     * Checks if we have still a digital plane of specified width when
     * adding point \a p. The object is left unchanged whatever the
     * returned value. The invariant is 'this->isExtendable( p ) ==
     * true <=> this->extend( p ) == true'.
     *
     * @param p any 3D point (in the specified diameter).
     *
     * @return 'true' if this is still a plane, 'false' otherwise.
     */
    bool isExtendable( const Point & p ) const;

    //-------------------- model of CAdditivePrimitiveComputer -----------------------------
  public:

    /**
     * Adds the range of points [\a it, \a itE) and checks if we have
     * still a digital plane of specified width. The plane parameters
     * may be updated so as to include all the new points. All points
     * pointed by iterators should be in the diameter of this object.
     *
     * @tparam TInputIterator any model of InputIterator on Point.
     * @param it an iterator on the first element of the range of 3D points.
     * @param itE an iterator after the last element of the range of 3D points.
     *
     * @return 'true' if it is still a plane, 'false' otherwise (the
     * object is then in its original state).
     */
    template <typename TInputIterator>
    bool extend( TInputIterator it, TInputIterator itE );

    /**
     * Checks if we have still a digital plane of specified width when
     * adding the range of points [\a it, \a itE). The object is left
     * unchanged whatever the returned value.  All points pointed by
     * iterators should be in the diameter of this object. The
     * invariant is 'this->isExtendable( it, itE ) == true <=>
     * this->extend( it, itE ) == true'.
     *
     * @tparam TInputIterator any model of InputIterator on Point.
     * @param it an iterator on the first element of the range of 3D points.
     * @param itE an iterator after the last element of the range of 3D points.
     *
     * @return 'true' if this is still a plane, 'false' otherwise.
     */
    template <typename TInputIterator>
    bool isExtendable( TInputIterator it, TInputIterator itE ) const;

    //-------------------- Parameters services -----------------------------
  public:

    /**
     * @tparam Vector3D any type T such that T.operator[](int i)
     * returns a reference to a double. i ranges in 0,1,2.
     *
     * @param (updates) the current normal vector 
     */
    template <typename Vector3D>
    void getNormal( Vector3D & normal ) const;

    /**
     * @tparam Vector3D any type T such that T.operator[](int i)
     * returns a reference to a double. i ranges in 0,1,2.
     *
     * @param (updates) the current unit normal vector 
     */
    template <typename Vector3D>
    void getUnitNormal( Vector3D & normal ) const;

    /**
     * If n is the unit normal to the current plane, then n.x >= min
     * and n.x <= max are the two half-planes defining it.
     *
     * @param min the lower bound (corresponding to the unit vector).
     * @param max the upper bound (corresponding to the unit vector).
     */
    void getBounds( double & min, double & max ) const;

    /**
     * @pre ! empty()
     * @return the current minimal point of the plane, i.e. the one
     * with the smallest scalar product with the current normal
     * vector. Note that other points may also have a minimum value.
     */
    const Point & minimalPoint() const;

    /**
     * @pre ! empty()
     * @return the current maximal point of the plane, i.e. the one
     * with the highest scalar product with the current normal
     * vector. Note that other points may also have a maximum value.
     */
    const Point & maximalPoint() const;

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

    // ------------------------- Private Datas --------------------------------
  private:
    std::vector<Dimension> myAxes; /**< The list of active plane axes. Starts with {0,1,2}. At least one. */
    ChordComputer myComputers[ 3 ]; /**< The three COBA plane computers. */
    mutable std::vector<Dimension> _axesToErase; /**< Useful when erasing axes. */
    // ------------------------- Hidden services ------------------------------
  protected:


    // ------------------------- Internals ------------------------------------
  private:
  }; // end of class ChordGenericNaivePlane


  /**
   * Overloads 'operator<<' for displaying objects of class 'ChordGenericNaivePlane'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'ChordGenericNaivePlane' to write.
   * @return the output stream after the writing.
   */
  template <typename TPoint, typename TInternalScalar>
  std::ostream&
  operator<< ( std::ostream & out, const ChordGenericNaivePlane<TPoint, TInternalScalar> & object );

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "DGtal/geometry/surfaces/ChordGenericNaivePlane.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined ChordGenericNaivePlane_h

#undef ChordGenericNaivePlane_RECURSES
#endif // else defined(ChordGenericNaivePlane_RECURSES)

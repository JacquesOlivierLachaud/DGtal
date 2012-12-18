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
 * @file MinCircumcircleCurvatureEstimator.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2012/12/18
 *
 * Header file for module MinCircumcircleCurvatureEstimator.cpp
 *
 * This file is part of the DGtal library.
 *
 * @see testL1LengthEstimator.cpp, testLengthEstimators.cpp
 */

#if defined(MinCircumcircleCurvatureEstimator_RECURSES)
#error Recursive header files inclusion detected in MinCircumcircleCurvatureEstimator.h
#else // defined(MinCircumcircleCurvatureEstimator_RECURSES)
/** Prevents recursive inclusion of headers. */
#define MinCircumcircleCurvatureEstimator_RECURSES

#if !defined MinCircumcircleCurvatureEstimator_h
/** Prevents repeated inclusion of headers. */
#define MinCircumcircleCurvatureEstimator_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include "DGtal/base/Common.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // template class MinCircumcircleCurvatureEstimator
  /**
   * Description of template class 'MinCircumcircleCurvatureEstimator'
   * <p> \brief Aim: an estimator that uses the notion of global
   * curvature (a mapping r(x,y,z) -> RadiusCircumcircle(x, y, z ),
   * where x,y,z are points of the curve).
   *
   * @todo For now limited to 4-connected arithmetic DSS.
   * 
   * Model of @href CLocalCurveGeometricEstimator.
   *
   * @tparam TConstIterator a model of iterator on points. Type
   * IteratorCirculatorTraits<ConstIterator>::Value::Coordinate must
   * be defined.
   */
  template <typename TConstIterator>
  class MinCircumcircleCurvatureEstimator
  {
    // ----------------------- Standard services ------------------------------
  public:


    ///@todo CONCEPT CHECK sur ConstIterator
    typedef TConstIterator ConstIterator;

    typedef double Quantity;
  

    /**
     * Default Constructor.
     */
    MinCircumcircleCurvatureEstimator();
    
    /**
     * Destructor.
     */
    ~MinCircumcircleCurvatureEstimator();

  
    // ----------------------- Interface --------------------------------------
  public:
    
    /** 
     * Initialize the measure computation.
     * 
     * @param h grid size (must be >0).
     * @param itb begin iterator
     * @param ite end iterator
     */
    void init( const double h, const ConstIterator& itb, 
	       const ConstIterator& ite );
    

    /** 
     * Computation of the curvature at point pointed by \a it.
     * @pre init() method must be called before.
     * @param it an iterator on a point of the curve.
     * 
     * @return the curvature estimation at this point.
     */
    Quantity eval( const ConstIterator& it ) const;

    /** 
     * Computation of the curvature at all points in [itb, ite).
     *
     * @pre init() method must be called before.
     * @param itb first point to analyse.
     * @param ite after the last point to analyse.
     * @param ito an output iterator where curvature estimations are written.
     *
     * @return the last value of the output iterator.
     */
    template <typename OutputIterator>
    OutputIterator eval( ConstIterator itb,
			 ConstIterator ite,
			 OutputIterator ito ) const;

 
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
    
    ///Grid size.
    double myH;

    ///Copy of the range.
    ConstIterator myBeginIt;
    ConstIterator myEndIt;

    ///Boolean to make sure that init() has been called before eval().
    bool myIsInitBefore;

    typedef typename IteratorCirculatorTraits<ConstIterator>::Value::Coordinate Coordinate; 
    typedef ArithmeticalDSS<ConstIterator,Coordinate,4> DSSAlgorithm;
    typedef SaturatedSegmentation<DSSAlgorithm> Segmentation;

    /// Memorizes the segmentation into maximal DSS.
    Segmentation* mySegmentation;
      
    
  private:

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    MinCircumcircleCurvatureEstimator ( const MinCircumcircleCurvatureEstimator & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    MinCircumcircleCurvatureEstimator & operator= ( const MinCircumcircleCurvatureEstimator & other );

    // ------------------------- Internals ------------------------------------
  private:

  }; // end of class MinCircumcircleCurvatureEstimator


  /**
   * Overloads 'operator<<' for displaying objects of class 'MinCircumcircleCurvatureEstimator'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'MinCircumcircleCurvatureEstimator' to write.
   * @return the output stream after the writing.
   */
  template <typename T>
  std::ostream&
  operator<< ( std::ostream & out, const MinCircumcircleCurvatureEstimator<T> & object );

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "DGtal/geometry/curves/estimation/MinCircumcircleCurvatureEstimator.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined MinCircumcircleCurvatureEstimator_h

#undef MinCircumcircleCurvatureEstimator_RECURSES
#endif // else defined(MinCircumcircleCurvatureEstimator_RECURSES)

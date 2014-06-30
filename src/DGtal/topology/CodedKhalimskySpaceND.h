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
 * @file CodedKhalimskySpaceND.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5807), University of Savoie, France
 *
 * @date 2014/06/27
 *
 * Header file for module CodedKhalimskySpaceND.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(CodedKhalimskySpaceND_RECURSES)
#error Recursive header files inclusion detected in CodedKhalimskySpaceND.h
#else // defined(CodedKhalimskySpaceND_RECURSES)
/** Prevents recursive inclusion of headers. */
#define CodedKhalimskySpaceND_RECURSES

#if !defined CodedKhalimskySpaceND_h
/** Prevents repeated inclusion of headers. */
#define CodedKhalimskySpaceND_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <set>
#include <map>
#include <algorithm>
#include "DGtal/base/Common.h"
#include "DGtal/kernel/CInteger.h"
#include "DGtal/kernel/CUnsignedNumber.h"
#include "DGtal/kernel/PointVector.h"
#include "DGtal/kernel/SpaceND.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{


  /**
     @brief This class is useful for looping on all "interesting" coordinates of a
     cell. For instance, surfels in Z3 have two interesting coordinates (the
     ones spanned by the surfel).
     @code
     KSpace::Cell p;
     KnSpace::DirIterator q;
     for ( q = ks.uDirs( p ); q != 0; ++q )
     {
     KSpace::Dimension dir = *q;
     ...
     }
     @endcode
  */
  template < DGtal::Dimension dim, 
             typename TInteger = DGtal::int32_t>
  class CodedCellDirectionIterator
  {
    //Integer must be signed to characterize a ring.
    BOOST_CONCEPT_ASSERT(( CInteger<TInteger> ) );

  public:
    /// Self type.
    typedef CodedCellDirectionIterator<dim, TInteger> Self;
    /// Arithmetic ring induced by (+,-,*) and Integer numbers.
    typedef TInteger Integer;

  public:
    /**
     * Constructor from dirs.
     * @param cell any unsigned cell
     */
    CodedCellDirectionIterator( Integer dirs );

    /**
     * @return the current direction.
     */
    Dimension operator*() const;

    /**
     * Pre-increment. Go to next direction.
     */
    CodedCellDirectionIterator& operator++();

    /**
     * Fast comparison with unsigned integer (unused
     * parameter). Comparison is 'false' at the end of the iteration.
     *
     * @return 'true' if the iterator is finished.
     */
    bool operator!=( const Integer ) const;

    /**
     * @return 'true' if the iteration is ended.
     */
    bool end() const;

    /**
     * Slow comparison with other iterator. Useful to check for end of loop.
     * @param other any direction iterator.
     */
    bool operator!=( const CodedCellDirectionIterator & other ) const;

    /**
     * Slow comparison with other iterator.
     * @param other any direction iterator.
     */
    bool operator==( const CodedCellDirectionIterator & other ) const;

  private:
    /** the directions to visit. */
    Integer myDirs;
    /** the current direction. */
    Dimension myDir;

  private:
    /** Look for next valid coordinate. */
    void find();
  };


  /////////////////////////////////////////////////////////////////////////////
  // template class CodedKhalimskySpaceND
  /**
   * Description of template class 'CodedKhalimskySpaceND' <p>
   *
   * \brief Aim: This class is a model of CCellularGridSpaceND. It
   * represents the cubical grid as a cell complex, whose cells are
   * defined as a coded array of integers. The topology of the cells
   * is defined by the parity of the coordinates (even: closed, odd:
   * open). Both topology and coordinates are coded as one integer of
   * type TCode.
   *
   * This space is finite. The user should choose between a closed
   * (default) cell space or an open cell space at initialisation (see
   * method CodedKhalimskySpaceND::init).
   *
   * @tparam dim the dimension of the digital space.
   * @tparam TInteger the Integer class used to specify the arithmetic computations on coordinates (default type = int32).
   * @tparam TCode the unsigned integer class used to code all coordinates, topology and orientation of a cell (default type = DGtal::uint64_t).
   *
   * @note Essentially a backport from
   * [ImaGene](https://gforge.liris.cnrs.fr/projects/imagene).
  */
  template < Dimension dim,
             typename TInteger = DGtal::int32_t,
             typename TCode = DGtal::uint64_t>
  class CodedKhalimskySpaceND
  {
    //Integer must be signed to characterize a ring.
    BOOST_CONCEPT_ASSERT(( CInteger<TInteger> ) );
    BOOST_CONCEPT_ASSERT(( CUnsignedNumber<TCode> ) );

  public:
    /// Self type.
    typedef CodedKhalimskySpaceND<dim, TInteger, TCode> Self;
    /// Arithmetic ring induced by (+,-,*) and Integer numbers.
    typedef TInteger Integer;
    /// The code storing integer coordinates, topology and orientation
    typedef TCode Code;

    ///Type used to represent sizes in the digital space.
    typedef typename NumberTraits<Integer>::UnsignedVersion Size;

    /**
     * A bit field describes a subpart of the code of a cell. There is
     * one field per coordinate, one field to store the topology, and a field for
     * the orientation for signed cells.
     * Codes a field mask, its inverse, the shift, and the number of masked bits
     * within an instance of Code.
     */
    struct BitField {
      /// Default constructor.
      BitField() {}
      /**
       * Constructor provided for convenience. If your mask is like 0x00f0, the constructor
       * should be called with (0x0f00, 8, 4).
       * @param _mask the mask for the bit field
       * @param _shift the index of the first non-zero bit
       * @param _nb_bits the number of non-zero bits
       */
      BitField( Code _mask, Dimension _shift, Dimension _nb_bits ) 
        : mask( _mask ), inv_mask( ~mask ), shift( _shift ), nb_bits( _nb_bits ) {}
      /// The consecutive bits of the word where the field is defined.
      Code mask;
      /// The bits of the word not in the field.
      Code inv_mask;
      /// The shift in the word to reach the field.
      Dimension shift;
      /// The number of consecutive bits of the field.
      Dimension nb_bits;

      /**
       * Given a cell code \a c, returns the part of the code
       * corresponding to the given field.
       * @param c any cell or signed cell code 
       * @return the code \a c masked by the mask of \a bf.
       */
      Code select( Code c ) const
      { return c & mask; }

      /**
       * Given a cell code \a c, returns the part of the code
       * without the part corresponding to the given field.
       * @param c any cell or signed cell code 
       * @return the code \a c masked by the inverted mask of \a bf.
       */
      Code unselect( Code c ) const
      { return c & inv_mask; }

      /**
       * Given a cell code \a c, returns the part of the code
       * with the part corresponding to the given field flipped.
       * @param c any cell or signed cell code 
       * @return the code \a c flipped by the  mask of \a bf.
       */
      Code inverse( Code c ) const
      { return c ^ mask; }

      /**
       * Given a cell code \a c, flips the part corresponding to this field.
       * @param[in,out] c any cell or signed cell code 
       */
      void flips( Code& c ) const
      { c ^= mask; }

      /**
       * @param c any cell or signed cell code 
       * @return the code \a c with all bits of the bit field set to 1.
       */
      Code one( Code c ) const
      { return c | mask; }

      /**
       * @param c any cell or signed cell code 
       * @return the code \a c with all bits of the bit field set to 0.
       */
      Code zero( Code c ) const
      { return c & inv_mask; }

 
      /**
       * Given a cell code \a c and a value \a v, modifies the
       * code to set the given value in the given field.
       * @param c any cell or signed cell code 
       * @param v an integer value (within bounds).
       * @return the code \a c updated.
       */
      Code codeWithValue( Code c, Integer v ) const
      { return ( c & inv_mask ) | codeValue( v ); }

      /**
       * Given a cell code \a c and a value \a v, modifies the
       * code to set the given value in the given field.
       * @param[in,out] c any cell or signed cell code 
       * @param v an integer value (within bounds).
       */
      void setValue( Code& c, Integer v ) const
      { c &= inv_mask; c |= codeValue( v ); }

      /**
       * Given a cell code \a c and a coded value \a d
       * returns the code with the given coded value in the given field.
       */
      Code codeWithCodedValue( Code c, Code d ) const
      { return ( c & inv_mask ) | d; }

      /**
       * Given a cell code \a c, a field \a bf and a coded value [FCODE], 
       * modifies the code to set the given coded value in the given field.
       */
      void changeCode( Code& c, Code d ) const
      { c = ( c & inv_mask ) | d; }

      /**
       * Given a value \a v, returns the coded value.
       */
      Code codeValue( Integer v ) const
      { return ((Code)v) << shift; }

      /**
       * Given a field \a bf and a cell code \a c, returns the value of this 
       * field.
       */
      Integer value( Code c ) const
      { return ( c & mask ) >> shift; }

      /**
       * Given a bit index \a b in the field, returns the mask for
       * checking the bit of this field.
       */
      Code bitSelector( Dimension b ) const
      { return NumberTraits<Code>::ONE << ( b + shift ); }

      /**
       * Given a cell code \a c and a bit index \a b in the 
       * field, returns true/false depending if this bit is set in the specified
       * field of the code.
       */
      bool bitCheck( Code c, Dimension b ) const
      { return c & bitSelector( b ); }

    };
    /// a Vector of BitField.
    typedef std::vector<BitField> BitFields;


    // Cells
    /// type for common cells of cellular space.
    typedef Code Cell;
    /// type for common signed cells of cellular space.
    typedef Code SCell;
    /// type for an array of codes
    typedef std::vector<Code> Codes;
    /// type for signed d-1-cells of cellular space.
    typedef SCell Surfel;
    /// type for sign of signed cells
    typedef bool Sign;
    /// type for iterating over the directions of a cell.
    typedef CodedCellDirectionIterator< dim, Integer > DirIterator;

    //Points and Vectors
    typedef PointVector< dim, Integer > Point;
    typedef PointVector< dim, Integer > Vector;

    typedef SpaceND<dim, Integer> Space;
    typedef CodedKhalimskySpaceND<dim, Integer> KhalimskySpace;

#if defined ( WIN32 )
    // static constants
    static const Dimension dimension = dim;
    static const Dimension DIM = dim;
    static const Sign POS = true;
    static const Sign NEG = false;
#else
    // static constants
    static const Dimension dimension = dim;
    static const Dimension DIM;
    static const Sign POS;
    static const Sign NEG;
#endif //WIN32

    /**
     * Type used for representing a sequence of (signed or unsigned) cells
     * @tparam CellType either Cell or SCell.
     */
    template <typename CellType>
    struct AnyCellCollection : public std::vector<CellType> {
      typedef CellType Value;
      typedef typename std::vector<CellType> Container;
      typedef typename std::vector<CellType>::iterator Iterator;
      typedef typename std::vector<CellType>::const_iterator ConstIterator;
    };

    /// Type used for Neighborhoods, Incident cells, Faces and Cofaces
    typedef AnyCellCollection<Cell> Cells;
    /// Type used for Neighborhoods, Incident signed cells, signed Faces and Cofaces
    typedef AnyCellCollection<SCell> SCells;

    // Sets, Maps
    /// Preferred type for defining a set of Cell(s).
    typedef std::set<Cell> CellSet;
    /// Preferred type for defining a set of SCell(s).
    typedef std::set<SCell> SCellSet;
    /// Preferred type for defining a set of surfels (always signed cells).
    typedef std::set<SCell> SurfelSet;
    /// Template rebinding for defining the type that is a mapping
    /// Cell -> Value.
    template <typename Value> struct CellMap {
      typedef std::map<SCell,Value> Type;
    };
    /// Template rebinding for defining the type that is a mapping
    /// SCell -> Value.
    template <typename Value> struct SCellMap {
      typedef std::map<SCell,Value> Type;
    };
    /// Template rebinding for defining the type that is a mapping
    /// SCell -> Value.
    template <typename Value> struct SurfelMap {
      typedef std::map<SCell,Value> Type;
    };
    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Destructor.
     */
    ~CodedKhalimskySpaceND();

    /**
     * Default onstructor.
     * The object is invalid.
     */
    CodedKhalimskySpaceND();

    /**
     * Specifies the upper and lower bounds for the maximal cells in
     * this space.
     *
     * @param lower the lowest point in this space (digital coords)
     * @param upper the upper point in this space (digital coords)
     * @param closed 'true' if this space is closed, 'false' if open.
     *
     * @return true if the initialization was valid (ie, such bounds
     * are representable with these integers).
     */
    bool init( const Point & lower,
               const Point & upper,
               bool closed );

    // ------------------------- Basic services ------------------------------
  public:

    /**
       @param k a coordinate (from 0 to 'dim()-1').
       @return the width of the space in the [k]-dimension.
    */
    Size size( Dimension k ) const;
    /**
       @param k a coordinate (from 0 to 'dim()-1').
       @return the minimal coordinate in the [k]-dimension.
     */
    Integer min( Dimension k ) const;
    /**
       @param k a coordinate (from 0 to 'dim()-1').
       @return the maximal coordinate in the [k]-dimension.
     */
    Integer max( Dimension k ) const;
    /**
       @return the lower bound for digital points in this space.
    */
    const Point & lowerBound() const;
    /**
       @return the upper bound for digital points in this space.
    */
    const Point & upperBound() const;
    /**
       @return the lower bound for cells in this space.
    */
    Cell lowerCell() const;
    /**
       @return the upper bound for cells in this space.
    */
    Cell upperCell() const;

    /**
       @return 'true' iff the space is closed.
    */
    bool isSpaceClosed() const;

    // ----------------------- Cell creation services --------------------------
  public:

    /**
     * From the Khalimsky coordinates of a cell, builds the
     * corresponding unsigned cell.
     *
     * @param kp an integer point (Khalimsky coordinates of cell).
     * @return the unsigned cell.
     */
    Cell uCell( Point kp ) const;

    /**
     * From the digital coordinates of a point in Zn and a cell type,
     * builds the corresponding cell.
     *
     * @param p an integer point (digital coordinates of cell).
     * @param c another cell defining the topology.
     *
     * @return the cell having the topology of [c] and the given
     * digital coordinates [p].
     */
    Cell uCell( Point p, Cell c ) const;

    /**
     * From the Khalimsky coordinates of a cell and a sign, builds the
     * corresponding unsigned cell.
     *
     * @param kp an integer point (Khalimsky coordinates of cell).
     * @param sign the sign of the cell (either POS or NEG).
     * @return the unsigned cell.
     */
    SCell sCell( Point kp, Sign sign = POS ) const;

    /**
     * From the digital coordinates of a point in Zn and a signed cell
     * type, builds the corresponding signed cell.
     *
     * @param p an integer point (digital coordinates of cell).
     * @param c another cell defining the topology and sign.
     *
     * @return the cell having the topology and sign of [c] and the given
     * digital coordinates [p].
     */
    SCell sCell( Point p, SCell c ) const;

    /**
     * From the digital coordinates of a point in Zn, creates the spel
     * (cell of maximal dimension) with these coordinates.
     *
     * @param p an integer point (digital coordinates of cell).
     *
     * @return the spel having the given digital coordinates [p].
     */
    Cell uSpel( const Point & p ) const;

    /**
     * From the digital coordinates of a point in Zn, creates the spel
     * (cell of maximal dimension) with these coordinates.
     *
     * @param p an integer point (digital coordinates of cell).
     * @param sign the sign of the cell (either POS or NEG).
     *
     * @return the signed spel having the given digital coordinates [p].
     */
    SCell sSpel( const Point & p, Sign sign = POS ) const;

    /**
     * From the digital coordinates of a point in Zn, creates the pointel
     * (cell of dimension 0) with these coordinates.
     *
     * @param p an integer point (digital coordinates of cell).
     *
     * @return the pointel having the given digital coordinates [p].
     */
    Cell uPointel( const Point & p ) const;

    /**
     * From the digital coordinates of a point in Zn, creates the pointel
     * (cell of dimension 0) with these coordinates.
     *
     * @param p an integer point (digital coordinates of cell).
     * @param sign the sign of the cell (either POS or NEG).
     *
     * @return the signed pointel having the given digital coordinates [p].
     */
    SCell sPointel( const Point & p, Sign sign = POS ) const;


    // ----------------------- Read accessors to cells ------------------------
  public:
    /**
     * @param c any unsigned cell.
     * @param k any valid dimension.
     * @return its Khalimsky coordinate along [k].
     */
    Integer uKCoord( Cell c, Dimension k ) const;

    /**
     * @param c any unsigned cell.
     * @param k any valid dimension.
     * @return its digital coordinate  along [k].
     */
    Integer uCoord( Cell c, Dimension k ) const;

    /**
     * @param c any unsigned cell.
     * @return its Khalimsky coordinates.
     */
    Point uKCoords( Cell c ) const;

    /**
     * @param c any unsigned cell.
     * @return its digital coordinates.
     */
    Point uCoords( Cell c ) const;

    /**
     * @param c any signed cell.
     * @param k any valid dimension.
     * @return its Khalimsky coordinate along [k].
     */
    Integer sKCoord( SCell c, Dimension k ) const;

    /**
     * @param c any signed cell.
     * @param k any valid dimension.
     * @return its digital coordinate  along [k].
     */
    Integer sCoord( SCell c, Dimension k ) const;

    /**
     * @param c any signed cell.
     * @return its Khalimsky coordinates.
     */
    Point sKCoords( SCell c ) const;

    /**
     * @param c any signed cell.
     * @return its digital coordinates.
     */
    Point sCoords( SCell c ) const;

    /**
     * @param c any signed cell.
     * @return its sign.
     */
    Sign sSign( SCell c ) const;

    // ----------------------- Write accessors to cells ------------------------
  public:

    /**
     * Sets the [k]-th Khalimsky coordinate of [c] to [i].
     * @param c any unsigned cell.
     * @param k any valid dimension.
     * @param i an integer coordinate within the space.
     */
    void uSetKCoord( Cell & c, Dimension k, Integer i ) const;

    /**
     * Sets the [k]-th Khalimsky coordinate of [c] to [i].
     * @param c any signed cell.
     * @param k any valid dimension.
     * @param i an integer coordinate within the space.
     */
    void sSetKCoord( SCell & c, Dimension k, Integer i ) const;

    /**
     * Sets the [k]-th digital coordinate of [c] to [i].
     * @param c any unsigned cell.
     * @param k any valid dimension.
     * @param i an integer coordinate within the space.
     */
    void uSetCoord( Cell & c, Dimension k, Integer i ) const;

    /**
     * Sets the [k]-th digital coordinate of [c] to [i].
     * @param c any signed cell.
     * @param k any valid dimension.
     * @param i an integer coordinate within the space.
     */
    void sSetCoord( SCell & c, Dimension k, Integer i ) const;

    /**
     * Sets the Khalimsky coordinates of [c] to [kp].
     * @param c any unsigned cell.
     * @param kp the new Khalimsky coordinates for [c].
     */
    void uSetKCoords( Cell & c, const Point & kp ) const;

    /**
     * Sets the Khalimsky coordinates of [c] to [kp].
     * @param c any signed cell.
     * @param kp the new Khalimsky coordinates for [c].
     */
    void sSetKCoords( SCell & c, const Point & kp ) const;

    /**
     * Sets the digital coordinates of [c] to [kp].
     * @param c any unsigned cell.
     * @param kp the new Khalimsky coordinates for [c].
     */
    void uSetCoords( Cell & c, const Point & kp ) const;

    /**
     * Sets the digital coordinates of [c] to [kp].
     * @param c any signed cell.
     * @param kp the new Khalimsky coordinates for [c].
     */
    void sSetCoords( SCell & c, const Point & kp ) const;

    /**
     * Sets the sign of the cell.
     * @param c (modified) any signed cell.
     * @param s any sign.
     */
    void sSetSign( SCell & c, Sign s ) const;

    // -------------------- Conversion signed/unsigned ------------------------
  public:
    /**
     * Creates a signed cell from an unsigned one and a given sign.
     * @param p any unsigned cell.
     * @param s a sign.
     * @return the signed version of the cell [p] with sign [s].
     */
    SCell signs( Cell p, Sign s ) const;

    /**
     * Creates an unsigned cell from a signed one.
     * @param p any signed cell.
     * @return the unsigned version of the cell [p].
     */
    Cell unsigns( SCell p ) const;

    /**
     * Creates the signed cell with the inverse sign of [p].
     * @param p any signed cell.
     * @return the cell [p] with opposite sign.
     */
    SCell sOpp( SCell p ) const;

    // ------------------------- Cell topology services -----------------------
  public:
    /**
     * @param p any unsigned cell.
     * @return the topology word of [p].
     */
    Integer uTopology( Cell p ) const;

    /**
     * @param p any signed cell.
     * @return the topology word of [p].
     */
    Integer sTopology( SCell p ) const;

    /**
     * @param p any unsigned cell.
     * @return the dimension of the cell [p].
     */
    Dimension uDim( Cell p ) const;

    /**
     * @param p any signed cell.
     * @return the dimension of the cell [p].
     */
    Dimension sDim( SCell p ) const;

    /**
     * @param b any unsigned cell.
     * @return 'true' if [b] is a surfel (spans all but one coordinate).
     */
    bool uIsSurfel( Cell b ) const;

    /**
     * @param b any signed cell.
     * @return 'true' if [b] is a surfel (spans all but one coordinate).
     */
    bool sIsSurfel( SCell b ) const;

    /**
       @param p any cell.
       @param k any direction.
       @return 'true' if [p] is open along the direction [k].
    */
    bool uIsOpen( Cell p, Dimension k ) const;

    /**
       @param p any signed cell.
       @param k any direction.
       @return 'true' if [p] is open along the direction [k].
    */
    bool sIsOpen( SCell p, Dimension k ) const;

    // -------------------- Iterator services for cells ------------------------
  public:

    /**
       Given an unsigned cell [p], returns an iterator to iterate over
       each coordinate the cell spans. (A spel spans all coordinates; a
       surfel all but one, etc). Example:

       @code
       KSpace::Cell p;
       ...
       for ( KSpace::DirIterator q = ks.uDirs( p ); q != 0; ++q )
       {
         Dimension dir = *q;
       ...
       }
       @endcode

       @param p any unsigned cell.

       @return an iterator that points on the first coordinate spanned
       by the cell.
    */
    DirIterator uDirs( Cell p ) const;

    /**
       Given a signed cell [p], returns an iterator to iterate over
       each coordinate the cell spans. (A spel spans all coordinates; a
       surfel all but one, etc). Example:

       @code
       KSpace::SCell p;
       ...
       for ( KSpace::DirIterator q = ks.uDirs( p ); q != 0; ++q )
       {
         Dimension dir = *q;
       ...
       }
       @endcode

       @param p any signed cell.

       @return an iterator that points on the first coordinate spanned
       by the cell.
    */
    DirIterator sDirs( SCell p ) const;

    /**
       Given an unsigned cell [p], returns an iterator to iterate over each
       coordinate the cell does not span. (A spel spans all coordinates;
       a surfel all but one, etc). Example:

       @code
       KSpace::Cell p;
       ...
       for ( KSpace::DirIterator q = ks.uOrthDirs( p ); q != 0; ++q )
       {
         Dimension dir = *q;
       ...
       }
       @endcode

       @param p any unsigned cell.

       @return an iterator that points on the first coordinate spanned
       by the cell.
    */
    DirIterator uOrthDirs( Cell p ) const;

    /**
       Given a signed cell [p], returns an iterator to iterate over each
       coordinate the cell does not span. (A spel spans all coordinates;
       a surfel all but one, etc). Example:

       @code
       KSpace::SCell p;
       ...
       for ( KSpace::DirIterator q = ks.uOrthDirs( p ); q != 0; ++q )
       {
         Dimension dir = *q;
       ...
       }
       @endcode

       @param p any signed cell.

       @return an iterator that points on the first coordinate spanned
       by the cell.
    */
    DirIterator sOrthDirs( SCell p ) const;

    /**
       Given an unsigned surfel [s], returns its orthogonal direction (ie,
       the coordinate where the surfel is closed).

       @param s an unsigned surfel
       @return the orthogonal direction of [s]
    */
    Dimension uOrthDir( Cell s ) const;

    /**
       Given a signed surfel [s], returns its orthogonal direction (ie,
       the coordinate where the surfel is closed).

       @param s a signed surfel
       @return the orthogonal direction of [s]
    */
    Dimension sOrthDir( SCell s ) const;

    // -------------------- Unsigned cell geometry services --------------------
  public:

    /**
       @return the first cell of the space with the same type as [p].
    */
    Cell uFirst( Cell p ) const;

    /**
       @return the last cell of the space with the same type as [p].
    */
    Cell uLast( Cell p ) const;

    /**
       NB: you can go out of the space.
       @param p any cell.
       @param k the coordinate that is changed.

       @return the same element as [p] except for the incremented
       coordinate [k].
    */
    Cell uGetIncr( Cell p, Dimension k ) const;

    /**
       Useful to check if you are going out of the space.
       @param p any cell.
       @param k the tested coordinate.

       @return true if [p] cannot have its [k]-coordinate augmented
       without leaving the space.
    */
    bool uIsMax( Cell p, Dimension k ) const;


    /**
       Useful to check if you are going out of the space.
       @param p any cell.
       @param k the tested coordinate.

       @return true if [p] has its [k]-coordinate within the allowed bounds.
    */
    bool uIsInside( Cell p, Dimension k ) const;


    /**
       Useful to check if you are going out of the space.
       @param p any cell.
       @param k the concerned coordinate.

       @return the cell similar to [p] but with the maximum allowed
       [k]-coordinate.
    */
    Cell uGetMax( Cell p, Dimension k ) const;

    /**
       NB: you can go out of the space.
       @param p any cell.
       @param k the coordinate that is changed.

       @return the same element as [p] except for an decremented
       coordinate [k].
    */
    Cell uGetDecr( Cell p, Dimension k ) const;

    /**
       Useful to check if you are going out of the space.
       @param p any cell.
       @param k the tested coordinate.

       @return true if [p] cannot have its [k]-coordinate decreased
       without leaving the space.
    */
    bool uIsMin( Cell p, Dimension k ) const;

    /**
       Useful to check if you are going out of the space.
       @param p any cell.
       @param k the concerned coordinate.

       @return the cell similar to [p] but with the minimum allowed
       [k]-coordinate.
    */
    Cell uGetMin( Cell p, Dimension k ) const;


    /**
       NB: you can go out of the space.
       @param p any cell.
       @param k the coordinate that is changed.
       @param x the increment.

       @return the same element as [p] except for a coordinate [k]
       incremented with x.
    */
    Cell uGetAdd( Cell p, Dimension k, Integer x ) const;

    /**
       NB: you can go out of the space.
       @param p any cell.
       @param k the coordinate that is changed.
       @param x the decrement.

       @return the same element as [p] except for a coordinate [k]
       decremented with x.
    */
    Cell uGetSub( Cell p, Dimension k, Integer x ) const;

    /**
       Useful to check if you are going out of the space.
       @param p any cell.
       @param k the coordinate that is tested.
       @return the number of increment to do to reach the maximum value.
    */
    Integer uDistanceToMax( Cell p, Dimension k ) const;

    /**
       Useful to check if you are going out of the space.
       @param p any cell.
       @param k the coordinate that is tested.

       @return the number of decrement to do to reach the minimum
       value.
    */
    Integer uDistanceToMin( Cell p, Dimension k ) const;

    /**
       Add the vector [vec] to [p].
       NB: you can go out of the space.
       @param p any cell.
       @param vec any pointel.
       @return the unsigned code of the cell [p] translated by [coord].
    */
    Cell uTranslation( Cell p, const Vector & vec ) const;

    /**
       Return the projection of [p] along the [k]th direction toward
       [bound]. Otherwise said, p[ k ] == bound[ k ] afterwards.

       @param p any cell.
       @param bound the element acting as bound (same topology as p).
       @param k the concerned coordinate.
       @return the projection.
    */
    Cell uProjection( Cell p, Cell bound, Dimension k ) const;

    /**
       Projects [p] along the [k]th direction toward
       [bound]. Otherwise said, p[ k ] == bound[ k ] afterwards.

       @param [in,out] p any cell.
       @param [in] bound the element acting as bound (same topology as p).
       @param [in] k the concerned coordinate.
    */
    void uProject( Cell& p, Cell bound, Dimension k ) const;

    /**
       Increment the cell [p] to its next position (as classically done in
       a scanning). Example:

       \code
       KSpace K;
       Cell first, last; // lower and upper bounds
       Cell p = first;
       do
       { // ... whatever [p] is the current cell
       }
       while ( K.uNext( p, first, last ) );
       \endcode

       @param p any cell.
       @param lower the lower bound.
       @param upper the upper bound.

       @return true if p is still within the bounds, false if the
       scanning is finished.
    */
    bool uNext( Cell& p, Cell lower, Cell upper ) const;

    // -------------------- Signed cell geometry services --------------------
  public:

    /**
       @return the first cell of the space with the same type as [p].
    */
    SCell sFirst( SCell p ) const;

    /**
       @return the last cell of the space with the same type as [p].
    */
    SCell sLast( SCell p ) const;

    /**
       NB: you can go out of the space.
       @param p any cell.
       @param k the coordinate that is changed.

       @return the same element as [p] except for the incremented
       coordinate [k].
    */
    SCell sGetIncr( SCell p, Dimension k ) const;

    /**
       Useful to check if you are going out of the space.
       @param p any cell.
       @param k the tested coordinate.

       @return true if [p] cannot have its [k]-coordinate augmented
       without leaving the space.
    */
    bool sIsMax( SCell p, Dimension k ) const;

    /**
       Useful to check if you are going out of the space.
       @param p any cell.
       @param k the tested coordinate.

       @return true if [p] has its [k]-coordinate within the allowed bounds.
    */
    bool sIsInside( SCell p, Dimension k ) const;

    /**
       Useful to check if you are going out of the space.
       @param p any cell.
       @param k the concerned coordinate.

       @return the cell similar to [p] but with the maximum allowed
       [k]-coordinate.
    */
    SCell sGetMax( SCell p, Dimension k ) const;

    /**
       NB: you can go out of the space.
       @param p any cell.
       @param k the coordinate that is changed.

       @return the same element as [p] except for an decremented
       coordinate [k].
    */
    SCell sGetDecr( SCell p, Dimension k ) const;

    /**
       Useful to check if you are going out of the space.
       @param p any cell.
       @param k the tested coordinate.

       @return true if [p] cannot have its [k]-coordinate decreased
       without leaving the space.
    */
    bool sIsMin( SCell p, Dimension k ) const;

    /**
       Useful to check if you are going out of the space.
       @param p any cell.
       @param k the concerned coordinate.

       @return the cell similar to [p] but with the minimum allowed
       [k]-coordinate.
    */
    SCell sGetMin( SCell p, Dimension k ) const;

    /**
       NB: you can go out of the space.
       @param p any cell.
       @param k the coordinate that is changed.
       @param x the increment.

       @return the same element as [p] except for a coordinate [k]
       incremented with x.
    */
    SCell sGetAdd( SCell p, Dimension k, Integer x ) const;

    /**
       NB: you can go out of the space.
       @param p any cell.
       @param k the coordinate that is changed.
       @param x the decrement.

       @return the same element as [p] except for a coordinate [k]
       decremented with x.
    */
    SCell sGetSub( SCell p, Dimension k, Integer x ) const;

    /**
       Useful to check if you are going out of the space.
       @param p any cell.
       @param k the coordinate that is tested.
       @return the number of increment to do to reach the maximum value.
    */
    Integer sDistanceToMax( SCell p, Dimension k ) const;

    /**
       Useful to check if you are going out of the space.
       @param p any cell.
       @param k the coordinate that is tested.

       @return the number of decrement to do to reach the minimum
       value.
    */
    Integer sDistanceToMin( SCell p, Dimension k ) const;

    /**
       Add the vector [vec] to [p].
       NB: you can go out of the space.
       @param p any cell.
       @param vec any pointel.
       @return the signed code of the cell [p] translated by [coord].
    */
    SCell sTranslation( SCell p, const Vector & vec ) const;

    /**
       Return the projection of [p] along the [k]th direction toward
       [bound]. Otherwise said, p[ k ] == bound[ k ] afterwards.

       @param p any cell.
       @param bound the element acting as bound (same topology as p).
       @param k the concerned coordinate.
       @return the projection.
    */
    SCell sProjection( SCell p, SCell bound, Dimension k ) const;

    /**
       Projects [p] along the [k]th direction toward
       [bound]. Otherwise said, p[ k ] == bound[ k ] afterwards.

       @param p any cell.
       @param bound the element acting as bound (same topology as p).
       @param k the concerned coordinate.
    */
    void sProject( SCell & p, SCell bound, Dimension k ) const;

    /**
       Increment the cell [p] to its next position (as classically done in
       a scanning). Example:

       \code
       KSpace K;
       Cell first, last; // lower and upper bounds
       Cell p = first;
       do
       { // ... whatever [p] is the current cell
       }
       while ( K.uNext( p, first, last ) );
       \endcode

       @param p any cell.
       @param lower the lower bound.
       @param upper the upper bound.

       @return true if p is still within the bounds, false if the
       scanning is finished.
    */
    bool sNext( SCell & p, SCell lower, SCell upper ) const;

    // ----------------------- Neighborhood services --------------------------
  public:

    /**
       Computes the 1-neighborhood of the cell [c] and returns
       it. It is the set of cells with same topology that are adjacent
       to [c] and which are within the bounds of this space.

       @param cell the unsigned cell of interest.
       @return the cells of the 1-neighborhood of [cell].
    */
    Cells uNeighborhood( Cell cell ) const;

    /**
       Computes the 1-neighborhood of the cell [c] and returns
       it. It is the set of cells with same topology that are adjacent
       to [c] and which are within the bounds of this space.

       @param cell the signed cell of interest.
       @return the cells of the 1-neighborhood of [cell].
    */
    SCells sNeighborhood( SCell cell ) const;

    /**
       Computes the proper 1-neighborhood of the cell [c] and returns
       it. It is the set of cells with same topology that are adjacent
       to [c], different from [c] and which are within the bounds of
       this space.

       @param cell the unsigned cell of interest.
       @return the cells of the proper 1-neighborhood of [cell].
    */
    Cells uProperNeighborhood( Cell cell ) const;

    /**
       Computes the proper 1-neighborhood of the cell [c] and returns
       it. It is the set of cells with same topology that are adjacent
       to [c], different from [c] and which are within the bounds of
       this space.

       @param cell the signed cell of interest.
       @return the cells of the proper 1-neighborhood of [cell].
    */
    SCells sProperNeighborhood( SCell cell ) const;

    /**
       NB: you can go out of the space.
       @param p any cell.
       @param k the coordinate that is changed.
       @param up if 'true' the orientation is forward along axis
       [k], otherwise backward.

       @return the adjacent element to [p] along axis [k] in the given
       direction and orientation.

       @note It is an alias to 'up ? uGetIncr( p, k ) : uGetDecr( p, k )'.
    */
    Cell uAdjacent( Cell p, Dimension k, bool up ) const;
    /**
       NB: you can go out of the space.
       @param p any cell.
       @param k the coordinate that is changed.
       @param up if 'true' the orientation is forward along axis
       [k], otherwise backward.

       @return the adjacent element to [p] along axis [k] in the given
       direction and orientation.

       @note It is an alias to 'up ? sGetIncr( p, k ) : sGetDecr( p, k )'.
    */
    SCell sAdjacent( SCell p, Dimension k, bool up ) const;

    // ----------------------- Incidence services --------------------------
  public:

    /**
       @param c any unsigned cell.
       @param k any coordinate.

       @param up if 'true' the orientation is forward along axis
       [k], otherwise backward.

       @return the forward or backward unsigned cell incident to [c]
       along axis [k], depending on [forward].

       @note It may be a lower incident cell if [c] is open along axis
       [k], else an upper incident cell.

       @note The cell should have an incident cell in this
       direction/orientation.
    */
    Cell uIncident( Cell c, Dimension k, bool up ) const;

    /**
       @param c any signed cell.
       @param k any coordinate.

       @param up if 'true' the orientation is forward along axis
       [k], otherwise backward.

       @return the forward or backward signed cell incident to [c]
       along axis [k], depending on [forward]. It is worthy to note
       that the forward and backward cell have opposite
       sign. Furthermore, the sign of these cells is defined so as to
       satisfy a boundary operator.

       @note It may be a lower incident cell if [c] is open along axis
       [k], else an upper incident cell.

       @note The cell should have an incident cell in this
       direction/orientation.
    */
    SCell sIncident( SCell c, Dimension k, bool up ) const;

    /**
       @param c any unsigned cell.
       @return the cells directly low incident to c in this space.
    */
    Cells uLowerIncident( Cell c ) const;

    /**
       @param c any unsigned cell.
       @return the cells directly up incident to c in this space.
    */
    Cells uUpperIncident( Cell c ) const;

    /**
       @param c any signed cell.
       @return the signed cells directly low incident to c in this space.
       @note it is the lower boundary of c expressed as a list of signed cells.
    */
    SCells sLowerIncident( SCell c ) const;

    /**
       @param c any signed cell.
       @return the signed cells directly up incident to c in this space.
       @note it is the upper boundary of c expressed as a list of signed cells.
    */
    SCells sUpperIncident( SCell c ) const;

    /**
       @param c any unsigned cell.
       @return the proper faces of [c] (chain of lower incidence).
    */
    Cells uFaces( Cell c ) const;

    /**
       @param c any unsigned cell.
       @return the proper cofaces of [c] (chain of upper incidence).
    */
    Cells uCoFaces( Cell c ) const;

    /**
       Return 'true' if the direct orientation of [p] along [k] is in
       the positive coordinate direction. The direct orientation in a
       direction allows to go from positive incident cells to positive
       incident cells.  This means that
       @code
       K.sSign( K.sIncident( p, k, K.sDirect( p, k ) ) ) == K.POS
       @endcode
       is always true.

       @param p any signed cell.
       @param k any coordinate.

       @return the direct orientation of [p] along [k] (true is
       upward, false is backward).
    */
    bool sDirect( SCell p, Dimension k ) const;

    /**
       @param p any signed cell.
       @param k any coordinate.

       @return the direct incident cell of [p] along [k] (the incident
       cell along [k] whose sign is positive).
    */
    SCell sDirectIncident( SCell p, Dimension k ) const;

    /**
       @param p any signed cell.
       @param k any coordinate.

       @return the indirect incident cell of [p] along [k] (the incident
       cell along [k] whose sign is negative).
    */
    SCell sIndirectIncident( SCell p, Dimension k ) const;


    // ----------------------- Interface --------------------------------------
  public:

    /**
     * Writes/Displays the object on an output stream.
     * @param out the output stream where the object is written.
     */
    void selfDisplay ( std::ostream & out ) const;

    /**
     * Displays a cell in a readable format.
     * @param[in,out] out the output stream where the object is written.
     * @param c the cell to display.
     */
    void printCell( std::ostream & out, Cell c ) const;

    /**
     * Displays a signed cell in a readable format.
     * @param[in,out] out the output stream where the object is written.
     * @param c the signed cell to display.
     */
    void printSCell( std::ostream & out, SCell c ) const;

    /**
     * Displays a cell in a readable format and outputs endl.
     * @param[in,out] out the output stream where the object is written.
     * @param c the cell to display.
     */
    void printlnCell( std::ostream & out, Cell c ) const;

    /**
     * Displays a signed cell in a readable format  and outputs endl.
     * @param[in,out] out the output stream where the object is written.
     * @param c the signed cell to display.
     */
    void printlnSCell( std::ostream & out, SCell c ) const;

    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool isValid() const;

    // ------------------------- Protected Datas ------------------------------
  private:
    // ------------------------- Private Datas --------------------------------
  private:

    /// Lowest digital point in the space
    Point myLower;
    /// Uppermost digital point in the space
    Point myUpper;
    /// Lowest cell in the space (in Khalimsky Coordinates)
    Point myLowerCell;
    /// Uppermost cell in the space (in Khalimsky Coordinates)
    Point myUpperCell;
    /// Lowest valid cell in the space (in Khalimsky Coordinates)
    /// (greater than myLowerCell if ! myIsClosed)
    Point myLowerValidCell;
    /// Uppermost valid cell in the space (in Khalimsky Coordinates)
    /// (smaller than myUpperCell if ! myIsClosed)
    Point myUpperValidCell;
    /// If 'true' the space is closed (low dimensional cells on
    /// boundary), otherwise it is open (high dimensional cells on
    /// boundary).
    bool myIsClosed;
    /// The lowest codes for cell (array of size 1 << dimension)
    Codes myLowerCellCodes;
    /// The uppermost codes for cell (array of size 1 << dimension)
    Codes myUpperCellCodes;

    /// Size per dimension. It is also one more than the maximum
    /// coordinate allowed.
    Point mySizes;
    /// Size of the space, ie its number of spels.
    Code myNbSpels;
    /// Number of bits used to represent a cell in this space.
    Dimension myNbBits;
    /// Number of types of cells.
    Dimension myNbCellTypes;
    /// It is the coded value of the maxima in all dimensions as a code
    Code myCodedAllCoordMax;
    /**
     * For each coordinate, it is the coded value of the maximum.
     * Ex: if the y-coordinate is shift by 8 and has maximum 10, then
     * myCodedCoordMaxs[ 1 ] = &0x00000a00
     * NB: the coded value of the minimum is always 0x00000000.
     */
    Codes myCodedCoordMaxs;
    /**
     * For each coordinate, it is the coded value of the positive increment.
     * Ex: if the y-coordinate is shift by 8, then
     * myCodedCoordIncrs[ 1 ] = &0x00000100
     */
    Codes myCodedCoordIncrs;
    /**
     * For each possible cell type, it is the coded value of a positive 
     * increment in each coordinate set to 1 of the cell type. 
     * Ex: if the y-coordinate is shift by 8, then
     * myUCodedCellTypeIncrs[ 3 ] = &0x00000101
     */
    Codes myCodedCellTypeIncrs;
    /**
     * This array stores the different coded vertices of any cell-type. 
     * For instance, for the code 5 (101), then there are four possible
     * vertices: 0x00000000, 0x00000001, 0x00010000, 0x00010001
     * if x-coordinate and y-coordinate take each 8 bits.
     * Note that you should shift the code by '<< dim()' to access this array.
     * To get the number of vertices, just use 'KnUtils::countSetBits'.
     * @see KnUtils::countSetBits
     */
    Codes myCodedCoordsVtx;

    /// Field for selecting the coordinates field of a coded Cell or a SCell
    /// looks like: 0x00008fff for a 256x128 image
    BitField myBFAllCoords;
    /// Field array to reach the k-th coordinate within a coded Cell or SCell.
    BitFields myBFCoords;
    /// Field to reach the sign in a coded SCell (no meaning for Cell).
    BitField myBFSign;
    /// Field to reach all directions and the sign in a coded SCell (no meaning for Cell).
    BitField myBFAllDirsAndSign;
    /// Field to select all directions in a coded Cell or SCell.
    BitField myBFAllDirs;
    /// Field array to select one direction in a coded Cell or SCell.
    BitFields myBFDirs;
    /// LUT direction code -> orthogonal direction, valid for d-1
    /// cells. Invalid value is d.
    std::vector<Dimension> myOrthDir;
    /// LUT direction code -> tangent direction, valid for
    /// 1-cells. Invalid value is d.
    std::vector<Dimension> myTgtDir;

    // ------------------------- Incidence attributes ---------------------------
  private:
    /**
     * LUT to store if two cells are low-incident or not (size (1<<2*dim())).
     * Given two cells [p1] and [p2] of cell type [t1] and [t2], then
     * m_low_incident[ t1 << (1<<dim()) + t2 ] stores a boolean.
     * If 'false' then they cannot be low incident.
     * Otherwise, you should check the table 'm_uid_coded_coords_vtx' to get
     * the different possible coordinate increments where the bit of the code
     * is 1. 
     */
    std::vector<bool> myLowIncident;

    /**
     * LUT to store if two cells are up-incident or not (size (1<<2*dim())).
     * Given two cells [p1] and [p2] of cell type [t1] and [t2], then
     * m_up_incident[ t1 << (1<<dim()) + t2 ] stores a boolean.
     * If 'false' then they cannot be up incident.
     * Otherwise, you should check the table 'm_uid_coded_coords_vtx' to get
     * the different possible coordinate increments where the bit of the code
     * is 1. 
     */
    std::vector<bool> myUpIncident;

    /**
     * Array of size (1<<dim())*dim(). For a given unsigned cell-type [c] and a
     * given direction [k], stores the sign of the permutation of coordinates to 
     * get back [c]. 'true' is +1, 'false' is -1. 
     * Useful to compute sign of incidence.
     * Exemples: <pre>
     * if [c]=zyx and [k]=y then [c]-[k]=zx and sign(zyx,y.zx) = -1 (1 perm)
     * if [c]=zyx and [k]=x then [c]-[k]=zy and sign(zyx,x.zy) = +1 (2 perms)
     * if [c]=zyx and [k]=x then [c]-[k]=zy and sign(zyx,x.zy) = +1 (0 perms)
     * if [c]=tzyx and [k]=y then [c]-[k]=tzx and sign(tzyx,y.tzx) = -1 (2 perms)
     * To get the index in the array: m_permutation[ ( k << n ) + c ]
     * </pre>
     */
    std::vector<bool> myPermutation;

    /**
     * Array of size (1<<(dim()+1)*dim(). For a given signed cell-type [c] and
     * a given direction [k], stores the signed cell-type of the first incident.
     * It is a LUT for computing 1-incident cells.
     * To get the index in the array: m_permutation[ ( k << (n+1) ) + c ]
     * @see mySIncident2;
     */
    Codes mySIncident1;

    /**
     * Array of size (1<<(dim()+1)*dim(). For a given signed cell-type [c] and
     * a given direction [k], stores the signed cell-type of the second incident.
     * It is a LUT for computing 1-incident cells.
     * To get the index in the array: m_permutation[ ( k << (n+1) ) + c ]
     * @see mySIncident1
     */
    Codes mySIncident2;

    /**
     * Array of size (1<<(dim()+1)^2. For two given signed cell-type [c1] and
     * [c2], stores the low incidence value (c1:c2). Note that this array is a
     * very sparse matrix, since this value may be different from 0 only if 
     * c2 is one-dimension less than c1.
     * It is a LUT for computing incidence matrices.
     * To get the index in the array: m_permutation[ ( c1 << (n+1) ) + c2 ]
     */
    std::vector<DGtal::int8_t> mySLowIncidenceMatrix;

    /**
     * Array of size (1<<(dim()+1)*dim(). For a given sign cell-type
     * [c] and a given direction [k], stores 'true' at index `(k <<
     * (dim+1)) + c` if the direct orientation means a positive
     * coordinate displacement, 'false' otherwise. The direct
     * orientation is the one where the 1-incident cell is positive.
     */
    std::vector<bool> mySDirectOrientation;

    /// Stores all the binomials (n,k) with n = dimension.
    std::vector<Integer> myBinomials;

    // ------------------------- Vector attributes ------------------------------
    /**
     * An array of vector which represents the different vector basis 
     * for each type of cell. Any vector formed with linear combinations of
     * vector of the basis has a null scalar product with any vector of the
     * orthogonal basis.
     */
    std::vector<Vector> myBasis;
    /**
     * An array of vector which represents the different orthogonal vector basis 
     * for each type of cell.  Any vector formed with linear combinations of
     * vector of the orthogonal basis has a null scalar product with any vector
     * of the basis.
     */
    std::vector<Vector> myOrthBasis;


    // ------------------------- Hidden services ------------------------------
  protected:


  private:



    // ------------------------- Internals ------------------------------------
  private:
    /**
       Used by uFaces for computing incident faces.
    */
    void uAddFaces( Cells& faces, Cell c, Dimension axis ) const;

    /**
       Used by uCoFaces for computing incident cofaces.
    */
    void uAddCoFaces( Cells& cofaces, Cell c, Dimension axis ) const;

    /// Internal methods to compute incidence LUT.
    void computeIncidenceLUT();
    /// Internal methods to compute permutation LUT.
    void computePermutationLUT();
    /// Internal methods to compute direct orientation LUT.
    void computeDirectOrientationLUT();
    /// Internal methods to compute binomials LUT.
    void computeBinomials();

  }; // end of class CodedKhalimskySpaceND


  /**
   * Overloads 'operator<<' for displaying objects of class 'CodedKhalimskySpaceND'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'CodedKhalimskySpaceND' to write.
   * @return the output stream after the writing.
   */
  template < DGtal::Dimension dim, 
             typename TInteger,
             typename TCode >
  std::ostream&
  operator<< ( std::ostream & out,
               const CodedKhalimskySpaceND<dim, TInteger, TCode> & object );

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "DGtal/topology/CodedKhalimskySpaceND.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined CodedKhalimskySpaceND_h

#undef CodedKhalimskySpaceND_RECURSES
#endif // else defined(CodedKhalimskySpaceND_RECURSES)

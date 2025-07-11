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
 * @file ShortcutsGeometry.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2018/06/25
 *
 * Header file for module ShortcutsGeometry.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(ShortcutsGeometry_RECURSES)
#error Recursive header files inclusion detected in ShortcutsGeometry.h
#else // defined(ShortcutsGeometry_RECURSES)
/** Prevents recursive inclusion of headers. */
#define ShortcutsGeometry_RECURSES

#if !defined ShortcutsGeometry_h
/** Prevents repeated inclusion of headers. */
#define ShortcutsGeometry_h

//////////////////////////////////////////////////////////////////////////////
#include "DGtal/helpers/Shortcuts.h"
#include "DGtal/kernel/BasicPointPredicates.h"
#include "DGtal/geometry/volumes/distance/LpMetric.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "DGtal/geometry/volumes/distance/VoronoiMap.h"
#include "DGtal/geometry/volumes/distance/DistanceTransformation.h"
#include "DGtal/geometry/surfaces/estimation/TrueDigitalSurfaceLocalEstimator.h"
#include "DGtal/geometry/surfaces/estimation/VoronoiCovarianceMeasureOnDigitalSurface.h"
#include "DGtal/geometry/surfaces/estimation/VCMDigitalSurfaceLocalEstimator.h"
#include "DGtal/geometry/surfaces/estimation/IIGeometricFunctors.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantVolumeEstimator.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantCovarianceEstimator.h"
#include "DGtal/geometry/meshes/CorrectedNormalCurrentComputer.h"

#include "DGtal/dec/DiscreteExteriorCalculusFactory.h"
#include "DGtal/dec/ATSolver2D.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  namespace sgf = ::DGtal::functors::ShapeGeometricFunctors;

  /////////////////////////////////////////////////////////////////////////////
  // template class ShortcutsGeometry
  /**
   * Description of template class 'ShortcutsGeometry' <p> \brief Aim: This
   * class is used to simplify shape and surface creation. With it,
   * you can create new shapes and surface in a few lines. The
   * drawback is that you use specific types or objects, which could
   * lead to faster code or more compact data structures.
   *
   * @tparam TKSpace any cellular grid space, a model of
   * concepts::CCellularGridSpaceND like KhalimskySpaceND.
   */
  template  < typename TKSpace >
    class ShortcutsGeometry : public Shortcuts< TKSpace >
    {
      BOOST_CONCEPT_ASSERT(( concepts::CCellularGridSpaceND< TKSpace > ));
    public:
      typedef Shortcuts< TKSpace >                     Base;
      typedef ShortcutsGeometry< TKSpace >             Self;
      using Base::parametersKSpace;
      using Base::getKSpace;
      using Base::parametersDigitizedImplicitShape3D;

      // ----------------------- Usual space types --------------------------------------
    public:

      /// Digital cellular space
      typedef TKSpace                                  KSpace;
      /// Digital space
      typedef typename KSpace::Space                   Space;
      /// Integer numbers
      typedef typename Space::Integer                  Integer;
      /// Point with integer coordinates.
      typedef typename Space::Point                    Point;
      /// Vector with integer coordinates.
      typedef typename Space::Vector                   Vector;
      /// Vector with floating-point coordinates.
      typedef typename Space::RealVector               RealVector;
      /// Point with floating-point coordinates.
      typedef typename Space::RealPoint                RealPoint;
      /// Floating-point numbers.
      typedef typename RealVector::Component           Scalar;
      /// An (hyper-)rectangular domain.
      typedef HyperRectDomain<Space>                   Domain;
      /// The type for 8-bits gray-scale elements.
      typedef unsigned char                            GrayScale;

      // ----------------------- Shortcut types --------------------------------------
    public:
      /// defines a multi-variate polynomial : RealPoint -> Scalar
      typedef MPolynomial< Space::dimension, Scalar >      ScalarPolynomial;
      /// defines an implicit shape of the space, which is the
      /// zero-level set of a ScalarPolynomial.
      typedef ImplicitPolynomial3Shape<Space>              ImplicitShape3D;
      /// defines the digitization of an implicit shape.
      typedef GaussDigitizer< Space, ImplicitShape3D >     DigitizedImplicitShape3D;
      /// defines a black and white image with (hyper-)rectangular domain.
      typedef ImageContainerBySTLVector<Domain, bool>      BinaryImage;
      /// defines a grey-level image with (hyper-)rectangular domain.
      typedef ImageContainerBySTLVector<Domain, GrayScale> GrayScaleImage;
      /// defines a float image with (hyper-)rectangular domain.
      typedef ImageContainerBySTLVector<Domain, float>     FloatImage;
      /// defines a double image with (hyper-)rectangular domain.
      typedef ImageContainerBySTLVector<Domain, double>    DoubleImage;
      /// defines a set of surfels
      typedef typename KSpace::SurfelSet                   SurfelSet;
      /// defines a light container that represents a connected digital
      /// surface over a binary image.
      typedef LightImplicitDigitalSurface< KSpace, BinaryImage >  LightSurfaceContainer;
      /// defines a connected digital surface over a binary image.
      typedef ::DGtal::DigitalSurface< LightSurfaceContainer >    LightDigitalSurface;
      /// defines a heavy container that represents any digital surface.
      typedef SetOfSurfels< KSpace, SurfelSet >                   ExplicitSurfaceContainer;
      /// defines an arbitrary digital surface over a binary image.
      typedef ::DGtal::DigitalSurface< ExplicitSurfaceContainer > DigitalSurface;
      /// defines a connected or not indexed digital surface.
      typedef IndexedDigitalSurface< ExplicitSurfaceContainer >   IdxDigitalSurface;
      typedef typename LightDigitalSurface::Surfel                Surfel;
      typedef typename LightDigitalSurface::Cell                  Cell;
      typedef typename LightDigitalSurface::SCell                 SCell;
      typedef typename LightDigitalSurface::Vertex                Vertex;
      typedef typename LightDigitalSurface::Arc                   Arc;
      typedef typename LightDigitalSurface::Face                  Face;
      typedef typename LightDigitalSurface::ArcRange              ArcRange;
      typedef typename IdxDigitalSurface::Vertex                  IdxSurfel;
      typedef typename IdxDigitalSurface::Vertex                  IdxVertex;
      typedef typename IdxDigitalSurface::Arc                     IdxArc;
      typedef typename IdxDigitalSurface::ArcRange                IdxArcRange;
      typedef std::set< IdxSurfel >                               IdxSurfelSet;
      typedef std::vector< Surfel >                               SurfelRange;
      typedef std::vector< Cell >                                 CellRange;
      typedef std::vector< IdxSurfel >                            IdxSurfelRange;
      typedef std::vector< Scalar >                               Scalars;
      typedef std::vector< RealVector >                           RealVectors;
      typedef std::vector< RealPoint >                            RealPoints;

      typedef ::DGtal::Statistic<Scalar>                          ScalarStatistic;

      typedef sgf::ShapePositionFunctor<ImplicitShape3D>          PositionFunctor;
      typedef sgf::ShapeNormalVectorFunctor<ImplicitShape3D>      NormalFunctor;
      typedef sgf::ShapeMeanCurvatureFunctor<ImplicitShape3D>     MeanCurvatureFunctor;
      typedef sgf::ShapeGaussianCurvatureFunctor<ImplicitShape3D> GaussianCurvatureFunctor;
      typedef sgf::ShapeFirstPrincipalCurvatureFunctor<ImplicitShape3D> FirstPrincipalCurvatureFunctor;
      typedef sgf::ShapeSecondPrincipalCurvatureFunctor<ImplicitShape3D> SecondPrincipalCurvatureFunctor;
      typedef sgf::ShapeFirstPrincipalDirectionFunctor<ImplicitShape3D> FirstPrincipalDirectionFunctor;
      typedef sgf::ShapeSecondPrincipalDirectionFunctor<ImplicitShape3D> SecondPrincipalDirectionFunctor;
      typedef sgf::ShapePrincipalCurvaturesAndDirectionsFunctor<ImplicitShape3D> PrincipalCurvaturesAndDirectionsFunctor;

      typedef typename functors::IIPrincipalCurvaturesAndDirectionsFunctor<Space>::Quantity   CurvatureTensorQuantity;
      typedef std::vector< CurvatureTensorQuantity >              CurvatureTensorQuantities;

      typedef CorrectedNormalCurrentComputer<RealPoint, RealVector> CNCComputer;

      typedef TrueDigitalSurfaceLocalEstimator
        < KSpace, ImplicitShape3D, PositionFunctor >                TruePositionEstimator;
      typedef TrueDigitalSurfaceLocalEstimator
        < KSpace, ImplicitShape3D, NormalFunctor >                  TrueNormalEstimator;
      typedef TrueDigitalSurfaceLocalEstimator
        < KSpace, ImplicitShape3D, MeanCurvatureFunctor >           TrueMeanCurvatureEstimator;
      typedef TrueDigitalSurfaceLocalEstimator
        < KSpace, ImplicitShape3D, GaussianCurvatureFunctor >       TrueGaussianCurvatureEstimator;
      typedef TrueDigitalSurfaceLocalEstimator
        < KSpace, ImplicitShape3D, FirstPrincipalCurvatureFunctor > TrueFirstPrincipalCurvatureEstimator;
      typedef TrueDigitalSurfaceLocalEstimator
        < KSpace, ImplicitShape3D, SecondPrincipalCurvatureFunctor > TrueSecondPrincipalCurvatureEstimator;
      typedef TrueDigitalSurfaceLocalEstimator
        < KSpace, ImplicitShape3D, FirstPrincipalDirectionFunctor > TrueFirstPrincipalDirectionEstimator;
      typedef TrueDigitalSurfaceLocalEstimator
        < KSpace, ImplicitShape3D, SecondPrincipalDirectionFunctor > TrueSecondPrincipalDirectionEstimator;
      typedef TrueDigitalSurfaceLocalEstimator
        < KSpace, ImplicitShape3D, PrincipalCurvaturesAndDirectionsFunctor > TruePrincipalCurvaturesAndDirectionsEstimator;

      typedef ::DGtal::Mesh<RealPoint>                            Mesh;
      typedef ::DGtal::TriangulatedSurface<RealPoint>             TriangulatedSurface;
      typedef ::DGtal::PolygonalSurface<RealPoint>                PolygonalSurface;
      typedef std::map<Surfel, IdxSurfel>                         Surfel2Index;
      typedef std::map<Cell,   IdxVertex>                         Cell2Index;

      typedef DigitalSetByAssociativeContainer<Domain, std::unordered_set<typename Domain::Point>> DigitalSet;
      typedef functors::NotPointPredicate<DigitalSet> VoronoiPointPredicate;

      // ----------------------- Static services --------------------------------------
    public:

      // ----------------------- Exact geometry services ------------------------------
      /// @name Exact geometry services
      /// @{
    public:

      /// @return the parameters used in ShortcutsGeometry.
      static Parameters defaultParameters()
      {
        return parametersShapeGeometry()
          | parametersGeometryEstimation()
          | parametersATApproximation()
          | parametersVoronoiMap();
      }

      /// @return the parameters and their default values which are used
      /// to approximate the geometry of continuous shape.
      ///   - projectionMaxIter [    20]: the maximum number of iteration for the projection.
      ///   - projectionAccuracy[0.0001]: the zero-proximity stop criterion during projection.
      ///   - projectionGamma   [   0.5]: the damping coefficient of the projection.
      ///   - gridstep [  1.0]: the gridstep that defines the digitization (often called h).
      static Parameters parametersShapeGeometry()
      {
        return Parameters
          ( "projectionMaxIter",  20 )
          ( "projectionAccuracy", 0.0001 )
          ( "projectionGamma",    0.5 )
          ( "gridstep",           1.0 );
      }

      /// Given a space \a K, an implicit \a shape, a sequence of \a
      /// surfels, and a gridstep \a h, returns the closest positions on
      /// the surface at the specified surfels, in the same order.
      ///
      /// @note The surfel centroids are iteratively projected onto the
      /// implicit surface through a damped Newton process.
      ///
      /// @param[in] shape the implicit shape.
      /// @param[in] K the Khalimsky space whose domain encompasses the digital shape.
      /// @param[in] surfels the sequence of surfels that we project onto the shape's surface
      /// @param[in] params the parameters:
      ///   - projectionMaxIter [    20]: the maximum number of iteration for the projection.
      ///   - projectionAccuracy[0.0001]: the zero-proximity stop criterion during projection.
      ///   - projectionGamma   [   0.5]: the damping coefficient of the projection.
      ///   - gridstep [  1.0]: the gridstep that defines the digitization (often called h).
      ///
      /// @return the vector containing the true normals, in the same
      /// order as \a surfels.
      static RealPoints
        getPositions
        ( CountedPtr<ImplicitShape3D> shape,
          const KSpace&               K,
          const SurfelRange&          surfels,
          const Parameters&           params = parametersShapeGeometry() )
      {
        RealVectors         n_true_estimations;
        TruePositionEstimator true_estimator;
        int     maxIter = params[ "projectionMaxIter"  ].as<int>();
        double accuracy = params[ "projectionAccuracy" ].as<double>();
        double    gamma = params[ "projectionGamma"    ].as<double>();
        Scalar gridstep = params[ "gridstep"           ].as<Scalar>();
        true_estimator.attach( *shape );
        true_estimator.setParams( K, PositionFunctor(), maxIter, accuracy, gamma );
        true_estimator.init( gridstep, surfels.begin(), surfels.end() );
        true_estimator.eval( surfels.begin(), surfels.end(),
                             std::back_inserter( n_true_estimations ) );
        return n_true_estimations;
      }

      /// Given an implicit \a shape and a sequence of \a points,
      /// returns the closest positions on the surface at the specified
      /// points, in the same order.
      ///
      /// @param[in] shape the implicit shape.
      /// @param[in] points the sequence of points that we project onto the shape's surface.
      /// @param[in] params the parameters:
      ///   - projectionMaxIter [    20]: the maximum number of iteration for the projection.
      ///   - projectionAccuracy[0.0001]: the zero-proximity stop criterion during projection.
      ///   - projectionGamma   [   0.5]: the damping coefficient of the projection.
      ///   - gridstep [  1.0]: the gridstep that defines the digitization (often called h).
      ///
      /// @return the vector containing the projected points.
      static RealPoints
        getPositions
        ( CountedPtr<ImplicitShape3D> shape,
          const RealPoints&           points,
          const Parameters&           params = parametersShapeGeometry() )
      {
        RealPoints proj_points( points.size() );
        int     maxIter = params[ "projectionMaxIter"  ].as<int>();
        double accuracy = params[ "projectionAccuracy" ].as<double>();
        double    gamma = params[ "projectionGamma"    ].as<double>();
        for ( unsigned int i = 0; i < points.size(); ++i )
          proj_points[ i ] = shape->nearestPoint( points[ i ], accuracy,
                                                  maxIter, gamma );
        return proj_points;
      }

      /// Given a space \a K, an implicit \a shape, a sequence of \a
      /// surfels, and a gridstep \a h, returns the normal vectors at the
      /// specified surfels, in the same order.
      ///
      /// @note that the normal vector is approximated by projecting the
      /// surfel centroid onto the implicit 3D shape.
      ///
      /// @param[in] shape the implicit shape.
      /// @param[in] K the Khalimsky space whose domain encompasses the digital shape.
      /// @param[in] surfels the sequence of surfels at which we compute the normals
      /// @param[in] params the parameters:
      ///   - projectionMaxIter [    20]: the maximum number of iteration for the projection.
      ///   - projectionAccuracy[0.0001]: the zero-proximity stop criterion during projection.
      ///   - projectionGamma   [   0.5]: the damping coefficient of the projection.
      ///   - gridstep [  1.0]: the gridstep that defines the digitization (often called h).
      ///
      /// @return the vector containing the true normals, in the same
      /// order as \a surfels.
      static RealVectors
        getNormalVectors
        ( CountedPtr<ImplicitShape3D> shape,
          const KSpace&               K,
          const SurfelRange&          surfels,
          const Parameters&           params = parametersShapeGeometry() )
      {
        RealVectors         n_true_estimations;
        TrueNormalEstimator true_estimator;
        int     maxIter = params[ "projectionMaxIter"  ].as<int>();
        double accuracy = params[ "projectionAccuracy" ].as<double>();
        double    gamma = params[ "projectionGamma"    ].as<double>();
        Scalar gridstep = params[ "gridstep"           ].as<Scalar>();
        true_estimator.attach( *shape );
        true_estimator.setParams( K, NormalFunctor(), maxIter, accuracy, gamma );
        true_estimator.init( gridstep, surfels.begin(), surfels.end() );
        true_estimator.eval( surfels.begin(), surfels.end(),
                             std::back_inserter( n_true_estimations ) );
        return n_true_estimations;
      }

      /// Given a space \a K, an implicit \a shape, a sequence of \a
      /// surfels, and a gridstep \a h, returns the mean curvatures at the
      /// specified surfels, in the same order.
      ///
      /// @note that the mean curvature is approximated by projecting the
      /// surfel centroid onto the implicit 3D shape.
      ///
      /// @param[in] shape the implicit shape.
      /// @param[in] K the Khalimsky space whose domain encompasses the digital shape.
      /// @param[in] surfels the sequence of surfels at which we compute the normals
      /// @param[in] params the parameters:
      ///   - projectionMaxIter [    20]: the maximum number of iteration for the projection.
      ///   - projectionAccuracy[0.0001]: the zero-proximity stop criterion during projection.
      ///   - projectionGamma   [   0.5]: the damping coefficient of the projection.
      ///   - gridstep [  1.0]: the gridstep that defines the digitization (often called h).
      ///
      /// @return the vector containing the mean curvatures, in the same
      /// order as \a surfels.
      static Scalars
        getMeanCurvatures
        ( CountedPtr<ImplicitShape3D> shape,
          const KSpace&               K,
          const SurfelRange&          surfels,
          const Parameters&           params = parametersShapeGeometry() )
      {
        Scalars                n_true_estimations;
        TrueMeanCurvatureEstimator true_estimator;
        int     maxIter = params[ "projectionMaxIter"  ].as<int>();
        double accuracy = params[ "projectionAccuracy" ].as<double>();
        double    gamma = params[ "projectionGamma"    ].as<double>();
        Scalar gridstep = params[ "gridstep"           ].as<Scalar>();
        true_estimator.attach( *shape );
        true_estimator.setParams( K, MeanCurvatureFunctor(), maxIter, accuracy, gamma );
        true_estimator.init( gridstep, surfels.begin(), surfels.end() );
        true_estimator.eval( surfels.begin(), surfels.end(),
                             std::back_inserter( n_true_estimations ) );
        return n_true_estimations;
      }

    /// Given a SurfaceMesh, compute mean curvature at each face using 
    /// CorrectedNormalCurrent method.
    ///
    /// @warning In this code, only triangles with barycenters strictly inside the sphere are 
    /// considered.
    ///
    /// @param mesh The surface mesh
    /// @param faces The faces to compute curvature at
    /// @param params
    ///   - unit_u: Whether the computed normals should be normalized or not
    ///   - r-radius        [   3.0]: the constant for kernel radius parameter r in r(h)=r h^alpha (VCM,II,Trivial).
    ///   - alpha           [  0.33]: the parameter alpha in r(h)=r h^alpha (VCM, II, CNC)."
    ///   - gridstep        [   1.0]: the digitization gridstep (often denoted by h).
    /// @return The curvatures for each face of the mesh, in the same order `faces` 
    static Scalars
      getCNCMeanCurvatures
      ( CountedPtr<typename Base::SurfaceMesh>  mesh, 
        const typename Base::SurfaceMesh::Faces faces, 
        const Parameters&                       params = parametersShapeGeometry() )
      {
        using Face = typename Base::SurfaceMesh::Face;

        bool unit_u = params["unit_u"].as<int>();
        double radius = params["r-radius"].as<double>();
        double alpha  = params["alpha"].as<double>();
        double h      = params["gridstep"].as<double>();
        if ( alpha != 1.0 ) radius *= pow( h, alpha-1.0 );

        CNCComputer computer(*mesh, unit_u);
        
        const auto& mu0 = computer.computeMu0();
        const auto& mu1 = computer.computeMu1();
        
        Scalars curvatures(faces.size());
        for (size_t i = 0; i < faces.size(); ++i)
        {
          const auto center = mesh->faceCentroid(faces[i]);
          const auto area = mu0.measure(center, radius, faces[i]);
          const auto lmu1 = mu1.measure(center, radius, faces[i]);
          curvatures[i] = CNCComputer::meanCurvature(area, lmu1);
        }

        return curvatures;
      }

      /// Given a SurfaceMesh, compute mean curvature at each face using
      /// CorrectedNormalCurrent method.
      ///
      /// This overloads compute curvature for each face of the mesh.
      ///
      /// @warning In this code, only triangles with barycenters
      /// strictly inside the sphere are considered.
      ///
      /// @param mesh The surface mesh
      /// @param params
      ///   - unit_u: Whether the computed normals should be normalized or not
      ///   - r-radius        [   3.0]: the constant for kernel radius parameter
      ///   r in r(h)=r h^alpha (VCM,II,Trivial).
      ///   - alpha           [  0.33]: the parameter alpha in r(h)=r h^alpha
      ///   (VCM, II, CNC)."
      ///   - gridstep        [   1.0]: the digitization gridstep (often denoted
      ///   by h).
      /// @return The curvatures for each face of mesh, in the the order given by
      /// the mesh
      static Scalars getCNCMeanCurvatures(
      CountedPtr<typename Base::SurfaceMesh> mesh,
      const Parameters & params = parametersShapeGeometry() )
      {
        std::vector<typename Base::SurfaceMesh::Face> allFaces(mesh->nbFaces());
        std::iota(allFaces.begin(), allFaces.end(), 0);

        return getCNCMeanCurvatures(mesh, allFaces, params);
      }

      /// Given a SurfaceMesh, compute mean curvature at each face using
      /// CorrectedNormalCurrent method.
      ///
      /// @warning In this code, only triangles with barycenters strictly inside
      /// the sphere are considered.
      ///
      /// @tparam Any digital object convertible to surface mesh via
      /// Shortcuts::makePrimalSurfaceMesh
      /// @param digitalObject A digital object
      /// @param params
      ///   - unit_u: Whether the computed normals should be normalized or not
      ///   - r-radius        [   3.0]: the constant for kernel radius parameter
      ///   r in r(h)=r h^alpha (VCM,II,Trivial).
      ///   - alpha           [  0.33]: the parameter alpha in r(h)=r h^alpha
      ///   (VCM, II, CNC)."
      ///   - gridstep        [   1.0]: the digitization gridstep (often denoted
      ///   by h).
      /// @return The curvatures for each face of the triangulated surface object
      template <typename T>
      static Scalars getCNCMeanCurvatures(
      T & digitalObject, const Parameters & params = parametersShapeGeometry() )
      {
        CountedPtr<typename Base::SurfaceMesh> mesh = Base::makePrimalSurfaceMesh(digitalObject);
        return getCNCMeanCurvatures(mesh, params);
      }

      /// Given a space \a K, an implicit \a shape, a sequence of \a
      /// surfels, and a gridstep \a h, returns the gaussian curvatures at the
      /// specified surfels, in the same order.
      ///
      /// @note that the gaussian curvature is approximated by projecting the
      /// surfel centroid onto the implicit 3D shape.
      ///
      /// @param[in] shape the implicit shape.
      /// @param[in] K the Khalimsky space whose domain encompasses the digital shape.
      /// @param[in] surfels the sequence of surfels at which we compute the normals
      ///
      /// @param[in] params the parameters:
      ///   - projectionMaxIter [    20]: the maximum number of iteration for the projection.
      ///   - projectionAccuracy[0.0001]: the zero-proximity stop criterion during projection.
      ///   - projectionGamma   [   0.5]: the damping coefficient of the projection.
      ///   - gridstep [  1.0]: the gridstep that defines the digitization (often called h).
      ///
      /// @return the vector containing the gaussian curvatures, in the same
      /// order as \a surfels.
      static Scalars
        getGaussianCurvatures
        ( CountedPtr<ImplicitShape3D> shape,
          const KSpace&               K,
          const SurfelRange&          surfels,
          const Parameters&           params = parametersShapeGeometry() )
      {
        Scalars                n_true_estimations;
        TrueGaussianCurvatureEstimator true_estimator;
        int     maxIter = params[ "projectionMaxIter"  ].as<int>();
        double accuracy = params[ "projectionAccuracy" ].as<double>();
        double    gamma = params[ "projectionGamma"    ].as<double>();
        Scalar gridstep = params[ "gridstep"           ].as<Scalar>();
        true_estimator.attach( *shape );
        true_estimator.setParams( K, GaussianCurvatureFunctor(), maxIter, accuracy, gamma );
        true_estimator.init( gridstep, surfels.begin(), surfels.end() );
        true_estimator.eval( surfels.begin(), surfels.end(),
                             std::back_inserter( n_true_estimations ) );
        return n_true_estimations;
      }

      /// Given a SurfaceMesh, compute gaussian curvature at each face using
      /// CorrectedNormalCurrent method.
      ///
      /// @warning  In this code, only triangles with barycenters strictly
      /// inside the sphere are considered.
      ///
      /// @param mesh The surface mesh
      /// @param faces The faces to compute curvature at
      /// @param params
      ///   - unit_u: Whether the computed normals should be normalized or not
      ///   - r-radius        [   3.0]: the constant for kernel radius parameter
      ///   r in r(h)=r h^alpha (VCM,II,Trivial).
      ///   - alpha           [  0.33]: the parameter alpha in r(h)=r h^alpha
      ///   (VCM, II, CNC)."
      ///   - gridstep        [   1.0]: the digitization gridstep (often denoted
      ///   by h).
      /// @return The curvatures for each face of the mesh, in the same order
      /// `faces`
      static Scalars getCNCGaussianCurvatures(
      CountedPtr<typename Base::SurfaceMesh> mesh,
      const typename Base::SurfaceMesh::Faces & faces,
      const Parameters & params = parametersShapeGeometry() )
      {
        using Face = typename Base::SurfaceMesh::Face;

        bool unit_u = params["unit_u"].as<int>();
        double radius = params["r-radius"].as<double>();
        double alpha  = params["alpha"].as<double>();
        double h      = params["gridstep"].as<double>();
        if ( alpha != 1.0 ) radius *= pow( h, alpha-1.0 );

        CNCComputer computer(*mesh, unit_u);
        
        const auto& mu0 = computer.computeMu0();
        const auto& mu2 = computer.computeMu2();

        Scalars curvatures(faces.size());
        for (size_t i = 0; i < faces.size(); ++i)
        {
          const auto center = mesh->faceCentroid(faces[i]);
          const auto area = mu0.measure(center, radius, faces[i]);
          const auto lmu2 = mu2.measure(center, radius, faces[i]);
          curvatures[i] = CNCComputer::GaussianCurvature(area, lmu2);
        }
        
        return curvatures;
      }

      /// Given a SurfaceMesh, compute gaussian curvature at each face using
      /// CorrectedNormalCurrent method.
      ///
      /// This overloads compute curvature for each face of the mesh.
      ///
      /// @warning  In this code, only triangles with barycenters strictly
      /// inside the sphere are considered.
      ///
      /// @param mesh The surface mesh
      /// @param params
      ///   - unit_u: Whether the computed normals should be normalized or not
      ///   - r-radius        [   3.0]: the constant for kernel radius parameter
      ///   r in r(h)=r h^alpha (VCM,II,Trivial).
      ///   - alpha           [  0.33]: the parameter alpha in r(h)=r h^alpha
      ///   (VCM, II, CNC)."
      ///   - gridstep        [   1.0]: the digitization gridstep (often denoted
      ///   by h).
      /// @return The curvatures for each face of the mesh, in thn the order
      /// given by the mesh
      static Scalars getCNCGaussianCurvatures(
      CountedPtr<typename Base::SurfaceMesh> mesh,
      const Parameters & params = parametersShapeGeometry() )
      {
        std::vector<typename Base::SurfaceMesh::Face> allFaces(mesh->nbFaces());
        std::iota(allFaces.begin(), allFaces.end(), 0);
        
        return getCNCGaussianCurvatures(mesh, allFaces, params);
      }

      /// Given a SurfaceMesh, compute mean curvature at each face using
      /// CorrectedNormalCurrent method.
      ///
      /// @warning  In this code, only triangles with barycenters strictly
      /// inside the sphere are considered.
      ///
      /// @tparam T Any digital object convertible to surface mesh via
      /// Shortcuts::makePrimalSurfaceMesh
      /// @param digitalObject A digital object
      /// @param params
      ///   - unit_u: Whether the computed normals should be normalized or not
      ///   - r-radius        [   3.0]: the constant for kernel radius parameter
      ///   r in r(h)=r h^alpha (VCM,II,Trivial).
      ///   - alpha           [  0.33]: the parameter alpha in r(h)=r h^alpha
      ///   (VCM, II, CNC)."
      ///   - gridstep        [   1.0]: the digitization gridstep (often denoted
      ///   by h).
      /// @return The curvatures for each face of the triangulated surface object
      template <typename T>
      static Scalars getCNCGaussianCurvatures(
      T & digitalObject, const Parameters & params = parametersShapeGeometry() )
      {
        CountedPtr<typename Base::SurfaceMesh> mesh = Base::makePrimalSurfaceMesh(digitalObject);
        return getCNCGaussianCurvatures(mesh, params);
      }

      /// Given a space \a K, an implicit \a shape, a sequence of \a principal curvatures 
      /// at the specified surfels, in the same order.
      ///
      /// @note that the first principal curvature is approximated by projecting the
      /// surfel centroid onto the implicit 3D shape.
      ///
      /// @param[in] shape the implicit shape.
      /// @param[in] K the Khalimsky space whose domain encompasses the digital shape.
      /// @param[in] surfels the sequence of surfels at which we compute the normals
      ///
      /// @param[in] params the parameters:
      ///   - projectionMaxIter [    20]: the maximum number of iteration for the projection.
      ///   - projectionAccuracy[0.0001]: the zero-proximity stop criterion during projection.
      ///   - projectionGamma   [   0.5]: the damping coefficient of the projection.
      ///   - gridstep [  1.0]: the gridstep that defines the digitization (often called h).
      ///
      /// @return the vector containing the first principal curvatures, in the same
      /// order as \a surfels.
      static Scalars
      getFirstPrincipalCurvatures
        ( CountedPtr<ImplicitShape3D> shape,
          const KSpace&               K,
          const SurfelRange&          surfels,
          const Parameters&           params = parametersShapeGeometry() )
      {
        Scalars n_true_estimations;
        TrueFirstPrincipalCurvatureEstimator true_estimator;
        int     maxIter = params[ "projectionMaxIter"  ].as<int>();
        double accuracy = params[ "projectionAccuracy" ].as<double>();
        double    gamma = params[ "projectionGamma"    ].as<double>();
        Scalar gridstep = params[ "gridstep"           ].as<Scalar>();
        true_estimator.attach( *shape );
        true_estimator.setParams( K, FirstPrincipalCurvatureFunctor(),
				  maxIter, accuracy, gamma );
        true_estimator.init( gridstep, surfels.begin(), surfels.end() );
        true_estimator.eval( surfels.begin(), surfels.end(),
                             std::back_inserter( n_true_estimations ) );
        return n_true_estimations;
      }

      /// Given a space \a K, an implicit \a shape, a sequence of \a
      /// surfels, and a gridstep \a h, returns the second (greatest)
      /// principal curvatures at the specified surfels, in the same
      /// order.
      ///
      /// @note that the second principal curvature is approximated by projecting the
      /// surfel centroid onto the implicit 3D shape.
      ///
      /// @param[in] shape the implicit shape.
      /// @param[in] K the Khalimsky space whose domain encompasses the digital shape.
      /// @param[in] surfels the sequence of surfels at which we compute the normals
      ///
      /// @param[in] params the parameters:
      ///   - projectionMaxIter [    20]: the maximum number of iteration for the projection.
      ///   - projectionAccuracy[0.0001]: the zero-proximity stop criterion during projection.
      ///   - projectionGamma   [   0.5]: the damping coefficient of the projection.
      ///   - gridstep [  1.0]: the gridstep that defines the digitization (often called h).
      ///
      /// @return the vector containing the second principal curvatures, in the same
      /// order as \a surfels.
      static Scalars
      getSecondPrincipalCurvatures
        ( CountedPtr<ImplicitShape3D> shape,
          const KSpace&               K,
          const SurfelRange&          surfels,
          const Parameters&           params = parametersShapeGeometry() )
      {
        Scalars n_true_estimations;
        TrueSecondPrincipalCurvatureEstimator true_estimator;
        int     maxIter = params[ "projectionMaxIter"  ].as<int>();
        double accuracy = params[ "projectionAccuracy" ].as<double>();
        double    gamma = params[ "projectionGamma"    ].as<double>();
        Scalar gridstep = params[ "gridstep"           ].as<Scalar>();
        true_estimator.attach( *shape );
        true_estimator.setParams( K, SecondPrincipalCurvatureFunctor(),
				  maxIter, accuracy, gamma );
        true_estimator.init( gridstep, surfels.begin(), surfels.end() );
        true_estimator.eval( surfels.begin(), surfels.end(),
                             std::back_inserter( n_true_estimations ) );
        return n_true_estimations;
      }

      /// Given a space \a K, an implicit \a shape, a sequence of \a
      /// surfels, and a gridstep \a h, returns the first principal
      /// directions (corresponding to the smallest principal
      /// curvature) at the specified surfels, in the same order.
      ///
      /// @note that the first principal direction is approximated by projecting the
      /// surfel centroid onto the implicit 3D shape.
      ///
      /// @param[in] shape the implicit shape.
      /// @param[in] K the Khalimsky space whose domain encompasses the digital shape.
      /// @param[in] surfels the sequence of surfels at which we compute the normals
      ///
      /// @param[in] params the parameters:
      ///   - projectionMaxIter [    20]: the maximum number of iteration for the projection.
      ///   - projectionAccuracy[0.0001]: the zero-proximity stop criterion during projection.
      ///   - projectionGamma   [   0.5]: the damping coefficient of the projection.
      ///   - gridstep [  1.0]: the gridstep that defines the digitization (often called h).
      ///
      /// @return the vector containing the first principal directions, in the same
      /// order as \a surfels.
      static RealVectors
      getFirstPrincipalDirections
        ( CountedPtr<ImplicitShape3D> shape,
          const KSpace&               K,
          const SurfelRange&          surfels,
          const Parameters&           params = parametersShapeGeometry() )
      {
        RealVectors n_true_estimations;
        TrueFirstPrincipalDirectionEstimator true_estimator;
        int     maxIter = params[ "projectionMaxIter"  ].as<int>();
        double accuracy = params[ "projectionAccuracy" ].as<double>();
        double    gamma = params[ "projectionGamma"    ].as<double>();
        Scalar gridstep = params[ "gridstep"           ].as<Scalar>();
        true_estimator.attach( *shape );
        true_estimator.setParams( K, FirstPrincipalDirectionFunctor(),
				  maxIter, accuracy, gamma );
        true_estimator.init( gridstep, surfels.begin(), surfels.end() );
        true_estimator.eval( surfels.begin(), surfels.end(),
                             std::back_inserter( n_true_estimations ) );
        return n_true_estimations;
      }

      /// Given a space \a K, an implicit \a shape, a sequence of \a
      /// surfels, and a gridstep \a h, returns the second principal
      /// directions (corresponding to the greatest principal
      /// curvature) at the specified surfels, in the same order.
      ///
      /// @note that the second principal direction is approximated by projecting the
      /// surfel centroid onto the implicit 3D shape.
      ///
      /// @param[in] shape the implicit shape.
      /// @param[in] K the Khalimsky space whose domain encompasses the digital shape.
      /// @param[in] surfels the sequence of surfels at which we compute the normals
      ///
      /// @param[in] params the parameters:
      ///   - projectionMaxIter [    20]: the maximum number of iteration for the projection.
      ///   - projectionAccuracy[0.0001]: the zero-proximity stop criterion during projection.
      ///   - projectionGamma   [   0.5]: the damping coefficient of the projection.
      ///   - gridstep [  1.0]: the gridstep that defines the digitization (often called h).
      ///
      /// @return the vector containing the second principal directions, in the same
      /// order as \a surfels.
      static RealVectors
      getSecondPrincipalDirections
        ( CountedPtr<ImplicitShape3D> shape,
          const KSpace&               K,
          const SurfelRange&          surfels,
          const Parameters&           params = parametersShapeGeometry() )
      {
        RealVectors n_true_estimations;
        TrueSecondPrincipalDirectionEstimator true_estimator;
        int     maxIter = params[ "projectionMaxIter"  ].as<int>();
        double accuracy = params[ "projectionAccuracy" ].as<double>();
        double    gamma = params[ "projectionGamma"    ].as<double>();
        Scalar gridstep = params[ "gridstep"           ].as<Scalar>();
        true_estimator.attach( *shape );
        true_estimator.setParams( K, SecondPrincipalDirectionFunctor(),
				  maxIter, accuracy, gamma );
        true_estimator.init( gridstep, surfels.begin(), surfels.end() );
        true_estimator.eval( surfels.begin(), surfels.end(),
                             std::back_inserter( n_true_estimations ) );
        return n_true_estimations;
      }

      /// Given a space \a K, an implicit \a shape, a sequence of \a
      /// surfels, and a gridstep \a h, returns the principal
      /// curvatures and principal directions as a tuple (k1, k2, d1, d2) at the
      /// specified surfels, in the same order.
      ///
      /// @note that the second principal direction is approximated by projecting the
      /// surfel centroid onto the implicit 3D shape.
      ///
      /// @param[in] shape the implicit shape.
      /// @param[in] K the Khalimsky space whose domain encompasses the digital shape.
      /// @param[in] surfels the sequence of surfels at which we compute the normals
      ///
      /// @param[in] params the parameters:
      ///   - projectionMaxIter [    20]: the maximum number of iteration for the projection.
      ///   - projectionAccuracy[0.0001]: the zero-proximity stop criterion during projection.
      ///   - projectionGamma   [   0.5]: the damping coefficient of the projection.
      ///   - gridstep [  1.0]: the gridstep that defines the digitization (often called h).
      ///
      /// @return the vector containing the principal curvatures and
      /// principal directions as a tuple (k1, k2, d1, d2), in the
      /// same order as \a surfels.
      static CurvatureTensorQuantities
      getPrincipalCurvaturesAndDirections
        ( CountedPtr<ImplicitShape3D> shape,
          const KSpace&               K,
          const SurfelRange&          surfels,
          const Parameters&           params = parametersShapeGeometry() )
      {
        CurvatureTensorQuantities n_true_estimations;
        TruePrincipalCurvaturesAndDirectionsEstimator true_estimator;
        int     maxIter = params[ "projectionMaxIter"  ].as<int>();
        double accuracy = params[ "projectionAccuracy" ].as<double>();
        double    gamma = params[ "projectionGamma"    ].as<double>();
        Scalar gridstep = params[ "gridstep"           ].as<Scalar>();
        true_estimator.attach( *shape );
        true_estimator.setParams( K, PrincipalCurvaturesAndDirectionsFunctor(),
				  maxIter, accuracy, gamma );
        true_estimator.init( gridstep, surfels.begin(), surfels.end() );
        true_estimator.eval( surfels.begin(), surfels.end(),
                             std::back_inserter( n_true_estimations ) );
        return n_true_estimations;
      }

      /// Given a SurfaceMesh, compute principal curvature at each face using
      /// CorrectedNormalCurrent method.
      ///
      /// @note If no normals are provided for the faces, the normals will be
      /// computed (and set) using vertex normals if they exist and positions
      /// otherwise.
      ///
      /// @warning  In this code, only triangles with barycenters strictly
      /// inside the sphere are considered.
      ///
      /// @param[in,out] mesh The surface mesh. The mesh will be modified if no
      /// face normals are provided.
      /// @param[in] faces The faces to compute curvature at
      /// @param[in] params
      ///   - unit_u: Whether the computed normals should be normalized or not
      ///   - r-radius        [   3.0]: the constant for kernel radius parameter
      ///   r in r(h)=r h^alpha (VCM,II,Trivial).
      ///   - alpha           [  0.33]: the parameter alpha in r(h)=r h^alpha
      ///   (VCM, II, CNC)."
      ///   - gridstep        [   1.0]: the digitization gridstep (often denoted
      ///   by h).
      /// @return The principal curvatures for each face of the mesh, in the same
      /// order as faces. The result is a 4-element tuples: [first curvatures,
      /// second curvatures, first directions, second directions].
      static std::tuple<Scalars, Scalars, RealVectors, RealVectors>
        getCNCPrincipalCurvaturesAndDirections
        ( CountedPtr<typename Base::SurfaceMesh>   mesh, 
          const typename Base::SurfaceMesh::Faces& faces, 
          const Parameters&                        params = parametersShapeGeometry() )
        {
          using Face = typename Base::SurfaceMesh::Face;

          bool unit_u = params["unit_u"].as<int>();
          double radius = params["r-radius"].as<double>();
          double alpha  = params["alpha"].as<double>();
          double h      = params["gridstep"].as<double>();
          if ( alpha != 1.0 ) radius *= pow( h, alpha-1.0 );

          CNCComputer computer(*mesh, unit_u);
          
          const auto& mu0  = computer.computeMu0();
          const auto& muxy = computer.computeMuXY();
          
          if (mesh->faceNormals().size() == 0)
          {
            // Try to use vertex normals if any
            if (mesh->vertexNormals().size() == 0)
              mesh->computeFaceNormalsFromPositions();
            else 
              mesh->computeFaceNormalsFromVertexNormals();
          }
          
          const auto& normals = mesh->faceNormals();

          Scalars k1(faces.size()), k2(faces.size());
          RealVectors d1(faces.size()), d2(faces.size());

          for (size_t i = 0; i < faces.size(); ++i)
          {
            const auto center = mesh->faceCentroid(faces[i]);
            const auto area  = mu0 .measure(center, radius, faces[i]);
            const auto lmuxy = muxy.measure(center, radius, faces[i]);
            std::tie(k1[i], k2[i], d1[i], d2[i]) = 
                CNCComputer::principalCurvatures(area, lmuxy, normals[faces[i]]);
          }

          return std::make_tuple(k1, k2, d1, d2);
        }

        /// Given a SurfaceMesh, compute principal curvature at each face using
        /// CorrectedNormalCurrent method.
        ///
        /// This overloads compute curvature for each face of the mesh.
        ///
        /// @note If no normals are provided for the faces, the normals will be
        /// computed (and set) using vertex normals if they exist and positions
        /// otherwise.
        ///
        /// @warning  In this code, only triangles with barycenters strictly
        /// inside the sphere are considered.
        ///
        /// @param[in,out] mesh The surface mesh. The mesh will be modified if
        /// no face normals are provided.
        /// @param[in] params
        ///   - unit_u: Whether the computed normals should be normalized or not
        ///   - r-radius        [   3.0]: the constant for kernel radius
        ///   parameter r in r(h)=r h^alpha (VCM,II,Trivial).
        ///   - alpha           [  0.33]: the parameter alpha in r(h)=r h^alpha
        ///   (VCM, II, CNC)."
        ///   - gridstep        [   1.0]: the digitization gridstep (often
        ///   denoted by h).
        /// @return The principal curvatures for each face of the mesh, in the
        /// same order as mesh faces. The result is a 4-element tuples: [first
        /// curvatures, second curvatures, first directions, second directions].
        static std::tuple<Scalars, Scalars, RealVectors, RealVectors>
        getCNCPrincipalCurvaturesAndDirections(
        CountedPtr<typename Base::SurfaceMesh> mesh,
        const Parameters & params = parametersShapeGeometry() )
        {
          std::vector<typename Base::SurfaceMesh::Face> allFaces(mesh->nbFaces());
          std::iota(allFaces.begin(), allFaces.end(), 0);

          return getCNCPrincipalCurvaturesAndDirections(mesh, allFaces, params);
        }

        /// Given a SurfaceMesh, compute principal curvature at each face using
        /// CorrectedNormalCurrent method.
        ///
        /// @warning In this code, only triangles with barycenters strictly
        /// inside the sphere are considered.
        ///
        /// @tparam T Any digital object convertible to surface mesh via
        /// Shortcuts::makePrimalSurfaceMesh
        /// @param digitalObject A digital object
        /// @param params
        ///   - unit_u: Whether the computed normals should be normalized or not
        ///   - r-radius        [   3.0]: the constant for kernel radius
        ///   parameter r in r(h)=r h^alpha (VCM,II,Trivial).
        ///   - alpha           [  0.33]: the parameter alpha in r(h)=r h^alpha
        ///   (VCM, II, CNC)."
        ///   - gridstep        [   1.0]: the digitization gridstep (often
        ///   denoted by h).
        /// @return The curvatures for each face of the triangulated surface
        /// object
        template <typename T>
        static std::tuple<Scalars, Scalars, RealVectors, RealVectors>
        getCNCPrincipalCurvaturesAndDirections(
        T & digitalObject,
        const Parameters & params = parametersShapeGeometry() )
        {
          CountedPtr<typename Base::SurfaceMesh> mesh = Base::makePrimalSurfaceMesh(digitalObject);
          return getCNCPrincipalCurvaturesAndDirections(mesh, params);
        }
         /// @}

      // --------------------------- geometry estimation ------------------------------
      /// @name Geometry estimation services
      /// @{
    public:

      /// @return the parameters and their default values which are used
      /// to estimate the geometry of a digital surface.
      ///   - verbose         [     1]: verbose trace mode 0: silent, 1: verbose.
      ///   - t-ring          [   3.0]: the radius used when computing convolved trivial normals (it is a graph distance, not related to the grid step).
      ///   - R-radius        [  10.0]: the constant for distance parameter R in R(h)=R h^alpha (VCM).
      ///   - r-radius        [   3.0]: the constant for kernel radius parameter r in r(h)=r h^alpha (VCM,II,Trivial).
      ///   - kernel          [ "hat"]: the kernel integration function chi_r, either "hat" or "ball". )
      ///   - alpha           [  0.33]: the parameter alpha in r(h)=r h^alpha (VCM, II)."
      ///   - surfelEmbedding [     0]: the surfel -> point embedding for VCM estimator: 0: Pointels, 1: InnerSpel, 2: OuterSpel.
      ///   - unit_u [0]: Use unit normals for (CNC) curvature computations.
      static Parameters parametersGeometryEstimation()
      {
        return Parameters
          ( "verbose",           1 )
          ( "t-ring",          3.0 )
          ( "kernel",        "hat" )
          ( "R-radius",       10.0 )
          ( "r-radius",        3.0 )
          ( "alpha",          0.33 )
          ( "surfelEmbedding",   0 )
          ( "unit_u"         ,   0 );
      }

      /// Given a digital space \a K and a vector of \a surfels,
      /// returns the trivial normals at the specified surfels, in the
      /// same order.
      ///
      /// @param[in] K the Khalimsky space whose domain encompasses the digital shape.
      /// @param[in] surfels the sequence of surfels at which we compute the normals
      ///
      /// @return the vector containing the trivial normal vectors, in the
      /// same order as \a surfels.
      static RealVectors
        getTrivialNormalVectors( const KSpace&      K,
                                 const SurfelRange& surfels )
      {
        std::vector< RealVector > result;
        for ( auto s : surfels )
          {
            Dimension  k = K.sOrthDir( s );
            bool  direct = K.sDirect( s, k );
            RealVector t = RealVector::zero;
            t[ k ]       = direct ? -1.0 : 1.0;
            result.push_back( t );
          }
        return result;
      }

      /// Given a digital surface \a surface, a sequence of \a surfels,
      /// and some parameters \a params, returns the convolved trivial
      /// normal vector estimations at the specified surfels, in the
      /// same order.
      ///
      /// @tparam TAnyDigitalSurface either kind of DigitalSurface, like ShortcutsGeometry::LightDigitalSurface or ShortcutsGeometry::DigitalSurface.
      ///
      /// @param[in] surface the digital surface
      /// @param[in] surfels the sequence of surfels at which we compute the normals
      /// @param[in] params the parameters:
      ///   - verbose         [     1]: verbose trace mode 0: silent, 1: verbose.
      ///   - t-ring          [   3.0]: the radius used when computing convolved trivial normals (it is a graph distance, not related to the grid step).
      ///
      /// @return the vector containing the estimated normals, in the
      /// same order as \a surfels.
      template <typename TAnyDigitalSurface>
        static RealVectors
        getCTrivialNormalVectors
        ( CountedPtr<TAnyDigitalSurface> surface,
          const SurfelRange&             surfels,
          const Parameters&              params = parametersGeometryEstimation() )
        {
          int    verbose = params[ "verbose"  ].as<int>();
          Scalar       t = params[ "t-ring"   ].as<double>();
          typedef typename TAnyDigitalSurface::DigitalSurfaceContainer  SurfaceContainer;
          typedef LpMetric<Space>                                       Metric;
          typedef functors::HatFunction<Scalar>                         Functor;
          typedef functors::ElementaryConvolutionNormalVectorEstimator
            < Surfel, CanonicSCellEmbedder<KSpace> >                    SurfelFunctor;
          typedef LocalEstimatorFromSurfelFunctorAdapter
            < SurfaceContainer, Metric, SurfelFunctor, Functor>         NormalEstimator;
          if ( verbose > 0 )
            trace.info() << "- CTrivial normal t-ring=" << t << " (discrete)" << std::endl;
          const Functor fct( 1.0, t );
          const KSpace &  K = surface->container().space();
          Metric    aMetric( 2.0 );
          CanonicSCellEmbedder<KSpace> canonic_embedder( K );
          std::vector< RealVector >    n_estimations;
          SurfelFunctor                surfelFct( canonic_embedder, 1.0 );
          NormalEstimator              estimator;
          estimator.attach( *surface);
          estimator.setParams( aMetric, surfelFct, fct, t );
          estimator.init( 1.0, surfels.begin(), surfels.end());
          estimator.eval( surfels.begin(), surfels.end(),
                          std::back_inserter( n_estimations ) );
          std::transform( n_estimations.cbegin(), n_estimations.cend(), n_estimations.begin(),
                          [] ( RealVector v ) { return -v; } );
          return n_estimations;
        }

      /// Given a digital surface \a surface, a sequence of \a surfels,
      /// and some parameters \a params, returns the normal Voronoi
      /// Covariance Measure (VCM) estimation at the specified surfels,
      /// in the same order.
      ///
      /// @tparam TAnyDigitalSurface either kind of DigitalSurface, like ShortcutsGeometry::LightDigitalSurface or ShortcutsGeometry::DigitalSurface.
      ///
      /// @param[in] surface the digital surface
      /// @param[in] surfels the sequence of surfels at which we compute the normals
      /// @param[in] params the parameters:
      ///   - verbose         [     1]: verbose trace mode 0: silent, 1: verbose.
      ///   - t-ring          [   3.0]: the radius used when computing convolved trivial normals (it is a graph distance, not related to the grid step).
      ///   - R-radius        [  10.0]: the constant for distance parameter R in R(h)=R h^alpha (VCM).
      ///   - r-radius        [   3.0]: the constant for kernel radius parameter r in r(h)=r h^alpha (VCM,II,Trivial).
      ///   - kernel          [ "hat"]: the kernel integration function chi_r, either "hat" or "ball". )
      ///   - alpha           [  0.33]: the parameter alpha in r(h)=r h^alpha (VCM, II)."
      ///   - surfelEmbedding [     0]: the surfel -> point embedding for VCM estimator: 0: Pointels, 1: InnerSpel, 2: OuterSpel.
      ///   - gridstep [  1.0]: the gridstep that defines the digitization (often called h).
      ///
      /// @return the vector containing the estimated normals, in the
      /// same order as \a surfels.
      template <typename TAnyDigitalSurface>
        static RealVectors
        getVCMNormalVectors
        ( CountedPtr<TAnyDigitalSurface> surface,
          const SurfelRange&             surfels,
          const Parameters&              params = parametersGeometryEstimation() )
        {
          typedef ExactPredicateLpSeparableMetric<Space,2> Metric;
          typedef typename TAnyDigitalSurface::DigitalSurfaceContainer SurfaceContainer;
          RealVectors n_estimations;
          int        verbose = params[ "verbose"   ].as<int>();
          std::string kernel = params[ "kernel"    ].as<std::string>();
          Scalar      h      = params[ "gridstep"  ].as<Scalar>();
          Scalar      R      = params[ "R-radius"  ].as<Scalar>();
          Scalar      r      = params[ "r-radius"  ].as<Scalar>();
          Scalar      t      = params[ "t-ring"    ].as<Scalar>();
          Scalar      alpha  = params[ "alpha"     ].as<Scalar>();
          int      embedding = params[ "embedding" ].as<int>();
          // Adjust parameters according to gridstep if specified.
          if ( alpha != 1.0 ) R *= pow( h, alpha-1.0 );
          if ( alpha != 1.0 ) r *= pow( h, alpha-1.0 );
          Surfel2PointEmbedding embType = embedding == 0 ? Pointels :
            embedding == 1 ? InnerSpel : OuterSpel;
          if ( verbose > 0 )
            {
              trace.info() << "- VCM normal kernel=" << kernel << " emb=" << embedding
                           << " alpha=" << alpha << std::endl;
              trace.info() << "- VCM normal r=" << (r*h)  << " (continuous) "
                           << r << " (discrete)" << std::endl;
              trace.info() << "- VCM normal R=" << (R*h)  << " (continuous) "
                           << R << " (discrete)" << std::endl;
              trace.info() << "- VCM normal t=" << t << " (discrete)" << std::endl;
            }
          if ( kernel == "hat" )
            {
              typedef functors::HatPointFunction<Point,Scalar>             KernelFunction;
              typedef VoronoiCovarianceMeasureOnDigitalSurface
                < SurfaceContainer, Metric, KernelFunction >               VCMOnSurface;
              typedef functors::VCMNormalVectorFunctor<VCMOnSurface>       NormalVFunctor;
              typedef VCMDigitalSurfaceLocalEstimator
                < SurfaceContainer, Metric, KernelFunction, NormalVFunctor> VCMNormalEstimator;
              KernelFunction chi_r( 1.0, r );
              VCMNormalEstimator estimator;
              estimator.attach( *surface );
              estimator.setParams( embType, R, r, chi_r, t, Metric(), verbose > 0 );
              estimator.init( h, surfels.begin(), surfels.end() );
              estimator.eval( surfels.begin(), surfels.end(),
                              std::back_inserter( n_estimations ) );
            }
          else if ( kernel == "ball" )
            {
              typedef functors::BallConstantPointFunction<Point,Scalar>    KernelFunction;
              typedef VoronoiCovarianceMeasureOnDigitalSurface
                < SurfaceContainer, Metric, KernelFunction >               VCMOnSurface;
              typedef functors::VCMNormalVectorFunctor<VCMOnSurface>       NormalVFunctor;
              typedef VCMDigitalSurfaceLocalEstimator
                < SurfaceContainer, Metric, KernelFunction, NormalVFunctor> VCMNormalEstimator;
              KernelFunction chi_r( 1.0, r );
              VCMNormalEstimator estimator;
              estimator.attach( *surface );
              estimator.setParams( embType, R, r, chi_r, t, Metric(), verbose > 0 );
              estimator.init( h, surfels.begin(), surfels.end() );
              estimator.eval( surfels.begin(), surfels.end(),
                              std::back_inserter( n_estimations ) );
            }
          else
            {
              trace.warning() << "[ShortcutsGeometry::getVCMNormalVectors] Unknown kernel: "
                              << kernel << std::endl;
            }
          return n_estimations;
        }

      /// Given a digital shape \a bimage, a sequence of \a surfels,
      /// and some parameters \a params, returns the normal Integral
      /// Invariant (II) estimation at the specified surfels, in the
      /// same order.
      ///
      /// @param[in] bimage the characteristic function of the shape as a binary image (inside is true, outside is false).
      /// @param[in] surfels the sequence of surfels at which we compute the normals
      /// @param[in] params the parameters:
      ///   - verbose         [     1]: verbose trace mode 0: silent, 1: verbose.
      ///   - r-radius        [   3.0]: the constant for kernel radius parameter r in r(h)=r h^alpha (VCM,II,Trivial).
      ///   - alpha           [  0.33]: the parameter alpha in r(h)=r h^alpha (VCM, II)."
      ///   - gridstep        [   1.0]: the digitization gridstep (often denoted by h).
      ///
      /// @return the vector containing the estimated normals, in the
      /// same order as \a surfels.
      ///
      /// @note Be careful, normals are reoriented with respect to
      /// Trivial normals. If you wish a more robust orientation, use
      /// getCTrivialNormalVectors.
      ///
      /// @note It is better to have surfels in a specific order, as
      /// given for instance by a depth-first traversal (see @ref getSurfelRange)
      static RealVectors
        getIINormalVectors( CountedPtr<BinaryImage> bimage,
                            const SurfelRange&      surfels,
                            const Parameters&       params
                            = parametersGeometryEstimation()
                            | parametersKSpace() )
      {
        auto K =  getKSpace( bimage, params );
        return getIINormalVectors( *bimage, K, surfels, params );
      }

      /// Given a digitized implicit shape \a dshape, a sequence of \a surfels,
      /// and some parameters \a params, returns the normal Integral
      /// Invariant (II) estimation at the specified surfels, in the
      /// same order.
      ///
      /// @param[in] dshape the digitized implicit shape, which is an
      /// implicitly defined characteristic function.
      ///
      /// @param[in] surfels the sequence of surfels at which we compute the normals
      /// @param[in] params the parameters:
      ///   - verbose         [     1]: verbose trace mode 0: silent, 1: verbose.
      ///   - r-radius        [   3.0]: the constant for kernel radius parameter r in r(h)=r h^alpha (VCM,II,Trivial).
      ///   - alpha           [  0.33]: the parameter alpha in r(h)=r h^alpha (VCM, II)."
      ///   - gridstep        [   1.0]: the digitization gridstep (often denoted by h).
      ///   - minAABB         [ -10.0]: the min value of the AABB bounding box (domain)
      ///   - maxAABB         [  10.0]: the max value of the AABB bounding box (domain)
      ///   - offset          [   5.0]: the digital dilation of the digital space,
      ///                       useful when you process shapes and that you add noise.
      ///   - closed          [     1]: specifies if the Khalimsky space is closed (!=0) or not (==0)
      ///
      /// @return the vector containing the estimated normals, in the
      /// same order as \a surfels.
      ///
      /// @note Be careful, normals are reoriented with respect to
      /// Trivial normals. If you wish a more robust orientation, use
      /// getCTrivialNormalVectors.
      ///
      /// @note It is better to have surfels in a specific order, as
      /// given for instance by a depth-first traversal (see @ref getSurfelRange)
      static RealVectors
        getIINormalVectors( CountedPtr< DigitizedImplicitShape3D > dshape,
                            const SurfelRange&      surfels,
                            const Parameters&       params
                            = parametersGeometryEstimation()
                            | parametersKSpace()
                            | parametersDigitizedImplicitShape3D() )
      {
        auto K =  getKSpace( params );
        return getIINormalVectors( *dshape, K, surfels, params );
      }

      /// Given an arbitrary PointPredicate \a shape: Point -> boolean, a Khalimsky
      /// space \a K, a sequence of \a surfels, and some parameters \a
      /// params, returns the normal Integral Invariant (II) estimation
      /// at the specified surfels, in the same order.
      ///
      /// @tparam TPointPredicate any type of map Point -> boolean.
      /// @param[in] shape a function Point -> boolean telling if you are inside the shape.
      /// @param[in] K the Khalimsky space where the shape and surfels live.
      /// @param[in] surfels the sequence of surfels at which we compute the normals
      /// @param[in] params the parameters:
      ///   - verbose         [     1]: verbose trace mode 0: silent, 1: verbose.
      ///   - r-radius        [   3.0]: the constant for kernel radius parameter r in r(h)=r h^alpha (VCM,II,Trivial).
      ///   - alpha           [  0.33]: the parameter alpha in r(h)=r h^alpha (VCM, II)."
      ///   - gridstep        [   1.0]: the digitization gridstep (often denoted by h).
      ///
      /// @return the vector containing the estimated normals, in the
      /// same order as \a surfels.
      ///
      /// @note Be careful, normals are reoriented with respect to
      /// Trivial normals. If you wish a more robust orientation, use
      /// getCTrivialNormalVectors.
      ///
      /// @note It is better to have surfels in a specific order, as
      /// given for instance by a depth-first traversal (see @ref getSurfelRange)
      template <typename TPointPredicate>
        static RealVectors
        getIINormalVectors( const TPointPredicate&  shape,
                            const KSpace&           K,
                            const SurfelRange&      surfels,
                            const Parameters&       params
                            = parametersGeometryEstimation()
                            | parametersKSpace() )
        {
          typedef functors::IINormalDirectionFunctor<Space> IINormalFunctor;
          typedef IntegralInvariantCovarianceEstimator
            <KSpace, TPointPredicate, IINormalFunctor>          IINormalEstimator;

          RealVectors n_estimations;
          int        verbose = params[ "verbose"   ].as<int>();
          Scalar     h       = params[ "gridstep"  ].as<Scalar>();
          Scalar     r       = params[ "r-radius"  ].as<Scalar>();
          Scalar     alpha   = params[ "alpha"     ].as<Scalar>();
          if ( alpha != 1.0 ) r *= pow( h, alpha-1.0 );
          if ( verbose > 0 )
            {
              trace.info() << "- II normal alpha=" << alpha << std::endl;
              trace.info() << "- II normal r=" << (r*h)  << " (continuous) "
                           << r << " (discrete)" << std::endl;
            }
          IINormalFunctor     functor;
          functor.init( h, r*h );
          IINormalEstimator   ii_estimator( functor );
          ii_estimator.attach( K, shape );
          ii_estimator.setParams( r );
          ii_estimator.init( h, surfels.begin(), surfels.end() );
          ii_estimator.eval( surfels.begin(), surfels.end(),
                             std::back_inserter( n_estimations ) );
          const RealVectors n_trivial = getTrivialNormalVectors( K, surfels );
          orientVectors( n_estimations, n_trivial );
          return n_estimations;
        }


      /// Given a digital shape \a bimage, a sequence of \a surfels,
      /// and some parameters \a vm, returns the mean curvature Integral
      /// Invariant (II) estimation at the specified surfels, in the
      /// same order.
      ///
      /// @param[in] bimage the characteristic function of the shape as a binary image (inside is true, outside is false).
      /// @param[in] surfels the sequence of surfels at which we compute the mean curvatures
      /// @param[in] params the parameters:
      ///   - verbose         [     1]: verbose trace mode 0: silent, 1: verbose.
      ///   - r-radius        [   3.0]: the constant for kernel radius parameter r in r(h)=r h^alpha (VCM,II,Trivial).
      ///   - alpha           [  0.33]: the parameter alpha in r(h)=r h^alpha (VCM, II)."
      ///   - gridstep        [   1.0]: the digitization gridstep (often denoted by h).
      ///
      /// @return the vector containing the estimated mean curvatures, in the
      /// same order as \a surfels.
      ///
      /// @note The function is faster when surfels are in a specific order, as
      /// given for instance by a depth-first traversal (see @ref getSurfelRange)
      static Scalars
        getIIMeanCurvatures( CountedPtr<BinaryImage> bimage,
                             const SurfelRange&      surfels,
                             const Parameters&       params
                             = parametersGeometryEstimation()
                             | parametersKSpace() )
      {
        auto K =  getKSpace( bimage, params );
        return getIIMeanCurvatures( *bimage, K, surfels, params );
      }

      /// Given a digitized implicit shape \a dshape, a sequence of \a surfels,
      /// and some parameters \a params, returns the mean curvature Integral
      /// Invariant (II) estimation at the specified surfels, in the
      /// same order.
      ///
      /// @param[in] dshape the digitized implicit shape, which is an
      /// implicitly defined characteristic function.
      ///
      /// @param[in] surfels the sequence of surfels at which we compute the normals
      /// @param[in] params the parameters:
      ///   - verbose         [     1]: verbose trace mode 0: silent, 1: verbose.
      ///   - r-radius        [   3.0]: the constant for kernel radius parameter r in r(h)=r h^alpha (VCM,II,Trivial).
      ///   - alpha           [  0.33]: the parameter alpha in r(h)=r h^alpha (VCM, II)."
      ///   - gridstep        [   1.0]: the digitization gridstep (often denoted by h).
      ///   - minAABB         [ -10.0]: the min value of the AABB bounding box (domain)
      ///   - maxAABB         [  10.0]: the max value of the AABB bounding box (domain)
      ///   - offset          [   5.0]: the digital dilation of the digital space,
      ///                       useful when you process shapes and that you add noise.
      ///   - closed          [     1]: specifies if the Khalimsky space is closed (!=0) or not (==0)
      ///
      /// @return the vector containing the estimated mean curvatures, in the
      /// same order as \a surfels.
      ///
      /// @note It is better to have surfels in a specific order, as
      /// given for instance by a depth-first traversal (see @ref getSurfelRange)
      static Scalars
        getIIMeanCurvatures( CountedPtr< DigitizedImplicitShape3D > dshape,
                             const SurfelRange&      surfels,
                             const Parameters&       params
                             = parametersGeometryEstimation()
                             | parametersKSpace()
                             | parametersDigitizedImplicitShape3D() )
      {
        auto K =  getKSpace( params );
        return getIIMeanCurvatures( *dshape, K, surfels, params );
      }


      /// Given an arbitrary PointPredicate \a shape: Point -> boolean, a Khalimsky
      /// space \a K, a sequence of \a surfels, and some parameters \a
      /// params, returns the mean curvature Integral
      /// Invariant (II) estimation at the specified surfels, in the
      /// same order.
      ///
      /// @tparam TPointPredicate any type of map Point -> boolean.
      /// @param[in] shape a function Point -> boolean telling if you are inside the shape.
      /// @param[in] K the Khalimsky space where the shape and surfels live.
      /// @param[in] surfels the sequence of surfels at which we compute the mean curvatures
      /// @param[in] params the parameters:
      ///   - verbose         [     1]: verbose trace mode 0: silent, 1: verbose.
      ///   - r-radius        [   3.0]: the constant for kernel radius parameter r in r(h)=r h^alpha (VCM,II,Trivial).
      ///   - alpha           [  0.33]: the parameter alpha in r(h)=r h^alpha (VCM, II)."
      ///   - gridstep        [   1.0]: the digitization gridstep (often denoted by h).
      ///
      /// @return the vector containing the estimated mean curvatures, in the
      /// same order as \a surfels.
      ///
      /// @note The function is faster when surfels are in a specific order, as
      /// given for instance by a depth-first traversal (see @ref getSurfelRange)
      template <typename TPointPredicate>
        static Scalars
        getIIMeanCurvatures( const TPointPredicate&  shape,
                             const KSpace&           K,
                             const SurfelRange&      surfels,
                             const Parameters&       params
                             = parametersGeometryEstimation()
                             | parametersKSpace() )
        {
          typedef functors::IIMeanCurvature3DFunctor<Space> IIMeanCurvFunctor;
          typedef IntegralInvariantVolumeEstimator
            <KSpace, TPointPredicate, IIMeanCurvFunctor>    IIMeanCurvEstimator;

          Scalars  mc_estimations;
          int      verbose = params[ "verbose"   ].as<int>();
          Scalar   h       = params[ "gridstep"  ].as<Scalar>();
          Scalar   r       = params[ "r-radius"  ].as<Scalar>();
          Scalar   alpha   = params[ "alpha"     ].as<Scalar>();
          if ( alpha != 1.0 ) r *= pow( h, alpha-1.0 );
          if ( verbose > 0 )
            {
              trace.info() << "- II mean curvature alpha=" << alpha << std::endl;
              trace.info() << "- II mean curvature r=" << (r*h)  << " (continuous) "
                           << r << " (discrete)" << std::endl;
            }
          IIMeanCurvFunctor   functor;
          functor.init( h, r*h );
          IIMeanCurvEstimator ii_estimator( functor );
          ii_estimator.attach( K, shape );
          ii_estimator.setParams( r );
          ii_estimator.init( h, surfels.begin(), surfels.end() );
          ii_estimator.eval( surfels.begin(), surfels.end(),
                             std::back_inserter( mc_estimations ) );
          return mc_estimations;
        }

      /// Given a digital shape \a bimage, a sequence of \a surfels,
      /// and some parameters \a vm, returns the Gaussian curvature Integral
      /// Invariant (II) estimation at the specified surfels, in the
      /// same order.
      ///
      /// @param[in] bimage the characteristic function of the shape as a binary image (inside is true, outside is false).
      /// @param[in] surfels the sequence of surfels at which we compute the Gaussian curvatures
      /// @param[in] params the parameters:
      ///   - verbose         [     1]: verbose trace mode 0: silent, 1: verbose.
      ///   - r-radius        [   3.0]: the constant for kernel radius parameter r in r(h)=r h^alpha (VCM,II,Trivial).
      ///   - alpha           [  0.33]: the parameter alpha in r(h)=r h^alpha (VCM, II)."
      ///   - gridstep        [   1.0]: the digitization gridstep (often denoted by h).
      ///
      /// @return the vector containing the estimated Gaussian curvatures, in the
      /// same order as \a surfels.
      ///
      /// @note The function is faster when surfels are in a specific order, as
      /// given for instance by a depth-first traversal (see @ref getSurfelRange)
      static Scalars
        getIIGaussianCurvatures( CountedPtr<BinaryImage> bimage,
                                 const SurfelRange&      surfels,
                                 const Parameters&       params
                                 = parametersGeometryEstimation()
                                 | parametersKSpace() )
      {
        auto K =  getKSpace( bimage, params );
        return getIIGaussianCurvatures( *bimage, K, surfels, params );
      }

      /// Given a digitized implicit shape \a dshape, a sequence of \a surfels,
      /// and some parameters \a params, returns the Gaussian curvature Integral
      /// Invariant (II) estimation at the specified surfels, in the
      /// same order.
      ///
      /// @param[in] dshape the digitized implicit shape, which is an
      /// implicitly defined characteristic function.
      ///
      /// @param[in] surfels the sequence of surfels at which we compute the normals
      /// @param[in] params the parameters:
      ///   - verbose         [     1]: verbose trace mode 0: silent, 1: verbose.
      ///   - r-radius        [   3.0]: the constant for kernel radius parameter r in r(h)=r h^alpha (VCM,II,Trivial).
      ///   - alpha           [  0.33]: the parameter alpha in r(h)=r h^alpha (VCM, II)."
      ///   - gridstep        [   1.0]: the digitization gridstep (often denoted by h).
      ///   - minAABB         [ -10.0]: the min value of the AABB bounding box (domain)
      ///   - maxAABB         [  10.0]: the max value of the AABB bounding box (domain)
      ///   - offset          [   5.0]: the digital dilation of the digital space,
      ///                       useful when you process shapes and that you add noise.
      ///   - closed          [     1]: specifies if the Khalimsky space is closed (!=0) or not (==0)
      ///
      /// @return the vector containing the estimated Gaussian curvatures, in the
      /// same order as \a surfels.
      ///
      /// @note It is better to have surfels in a specific order, as
      /// given for instance by a depth-first traversal (see @ref getSurfelRange)
      static Scalars
        getIIGaussianCurvatures( CountedPtr< DigitizedImplicitShape3D > dshape,
                                 const SurfelRange&      surfels,
                                 const Parameters&       params
                                 = parametersGeometryEstimation()
                                 | parametersKSpace()
                                 | parametersDigitizedImplicitShape3D() )
      {
        auto K =  getKSpace( params );
        return getIIGaussianCurvatures( *dshape, K, surfels, params );
      }


      /// Given an arbitrary PointPredicate \a shape: Point -> boolean, a Khalimsky
      /// space \a K, a sequence of \a surfels, and some parameters \a
      /// params, returns the Gaussian curvature Integral
      /// Invariant (II) estimation at the specified surfels, in the
      /// same order.
      ///
      /// @tparam TPointPredicate any type of map Point -> boolean.
      /// @param[in] shape a function Point -> boolean telling if you are inside the shape.
      /// @param[in] K the Khalimsky space where the shape and surfels live.
      /// @param[in] surfels the sequence of surfels at which we compute the Gaussian curvatures
      /// @param[in] params the parameters:
      ///   - verbose         [     1]: verbose trace mode 0: silent, 1: verbose.
      ///   - r-radius        [   3.0]: the constant for kernel radius parameter r in r(h)=r h^alpha (VCM,II,Trivial).
      ///   - alpha           [  0.33]: the parameter alpha in r(h)=r h^alpha (VCM, II)."
      ///   - gridstep        [   1.0]: the digitization gridstep (often denoted by h).
      ///
      /// @return the vector containing the estimated Gaussian curvatures, in the
      /// same order as \a surfels.
      ///
      /// @note The function is faster when surfels are in a specific order, as
      /// given for instance by a depth-first traversal (see @ref getSurfelRange)
      template <typename TPointPredicate>
        static Scalars
        getIIGaussianCurvatures( const TPointPredicate&  shape,
                                 const KSpace&           K,
                                 const SurfelRange&      surfels,
                                 const Parameters&       params
                                 = parametersGeometryEstimation()
                                 | parametersKSpace() )
        {
          typedef functors::IIGaussianCurvature3DFunctor<Space> IIGaussianCurvFunctor;
          typedef IntegralInvariantCovarianceEstimator
            <KSpace, TPointPredicate, IIGaussianCurvFunctor>    IIGaussianCurvEstimator;

          Scalars  mc_estimations;
          int      verbose = params[ "verbose"   ].as<int>();
          Scalar   h       = params[ "gridstep"  ].as<Scalar>();
          Scalar   r       = params[ "r-radius"  ].as<Scalar>();
          Scalar   alpha   = params[ "alpha"     ].as<Scalar>();
          if ( alpha != 1.0 ) r *= pow( h, alpha-1.0 );
          if ( verbose > 0 )
            {
              trace.info() << "- II Gaussian curvature alpha=" << alpha << std::endl;
              trace.info() << "- II Gaussian curvature r=" << (r*h)  << " (continuous) "
                           << r << " (discrete)" << std::endl;
            }
          IIGaussianCurvFunctor   functor;
          functor.init( h, r*h );
          IIGaussianCurvEstimator ii_estimator( functor );
          ii_estimator.attach( K, shape );
          ii_estimator.setParams( r );
          ii_estimator.init( h, surfels.begin(), surfels.end() );
          ii_estimator.eval( surfels.begin(), surfels.end(),
                             std::back_inserter( mc_estimations ) );
          return mc_estimations;
        }

      /// Given a digital shape \a bimage, a sequence of \a surfels,
      /// and some parameters \a vm, returns the principal curvatures and
      /// directions using an Integral
      /// Invariant (II) estimation at the specified surfels, in the
      /// same order.
      ///
      /// @param[in] bimage the characteristic function of the shape as a binary image (inside is true, outside is false).
      /// @param[in] surfels the sequence of surfels at which we compute the Gaussian curvatures
      /// @param[in] params the parameters:
      ///   - verbose         [     1]: verbose trace mode 0: silent, 1: verbose.
      ///   - r-radius        [   3.0]: the constant for kernel radius parameter r in r(h)=r h^alpha (VCM,II,Trivial).
      ///   - alpha           [  0.33]: the parameter alpha in r(h)=r h^alpha (VCM, II)."
      ///   - gridstep        [   1.0]: the digitization gridstep (often denoted by h).
      ///
      /// @return the vector containing the estimated Gaussian curvatures, in the
      /// same order as \a surfels.
      ///
      /// @note The function is faster when surfels are in a specific order, as
      /// given for instance by a depth-first traversal (see @ref getSurfelRange)
      static CurvatureTensorQuantities
      getIIPrincipalCurvaturesAndDirections( CountedPtr<BinaryImage> bimage,
                                          const SurfelRange&      surfels,
                                          const Parameters&       params
                                            = parametersGeometryEstimation()
                                            | parametersKSpace() )
      {
        auto K =  getKSpace( bimage, params );
        return getIIPrincipalCurvaturesAndDirections( *bimage, K, surfels, params );
      }

      /// Given a digital shape \a dshape, a sequence of \a surfels,
      /// and some parameters \a vm, returns the principal curvatures and
      /// directions using an Integral
      /// Invariant (II) estimation at the specified surfels, in the
      /// same order.
      ///
      /// @param[in] dshape the digitized implicit shape, which is an
      /// implicitly defined characteristic function.
      ///
      /// @param[in] surfels the sequence of surfels at which we compute the normals
      /// @param[in] params the parameters:
      ///   - verbose         [     1]: verbose trace mode 0: silent, 1: verbose.
      ///   - r-radius        [   3.0]: the constant for kernel radius parameter r in r(h)=r h^alpha (VCM,II,Trivial).
      ///   - alpha           [  0.33]: the parameter alpha in r(h)=r h^alpha (VCM, II)."
      ///   - gridstep        [   1.0]: the digitization gridstep (often denoted by h).
      ///   - minAABB         [ -10.0]: the min value of the AABB bounding box (domain)
      ///   - maxAABB         [  10.0]: the max value of the AABB bounding box (domain)
      ///   - offset          [   5.0]: the digital dilation of the digital space,
      ///                       useful when you process shapes adding some noise.
      ///   - closed          [     1]: specifies if the Khalimsky space is closed (!=0) or not (==0)
      ///
      /// @return the vector containing the estimated principal curvatures and directions, in the
      /// same order as \a surfels.
      ///
      /// @note It is better to have surfels in a specific order, as
      /// given for instance by a depth-first traversal (see @ref getSurfelRange)
      static CurvatureTensorQuantities
      getIIPrincipalCurvaturesAndDirections( CountedPtr< DigitizedImplicitShape3D > dshape,
                                          const SurfelRange&      surfels,
                                          const Parameters&       params
                                            = parametersGeometryEstimation()
                                            | parametersKSpace()
                                            | parametersDigitizedImplicitShape3D() )
      {
        auto K =  getKSpace( params );
        return getIIPrincipalCurvaturesAndDirections( *dshape, K, surfels, params );
      }


      /// Given an arbitrary PointPredicate \a shape: Point -> boolean, a Khalimsky
      /// space \a K, a sequence of \a surfels, and some parameters \a
      /// params, returns the principal curvatures and directions using Integral
      /// Invariant (II) estimation at the specified surfels, in the
      /// same order.
      ///
      /// @tparam TPointPredicate any type of map Point -> boolean.
      /// @param[in] shape a function Point -> boolean telling if you are inside the shape.
      /// @param[in] K the Khalimsky space where the shape and surfels live.
      /// @param[in] surfels the sequence of surfels at which we compute the Gaussian curvatures
      /// @param[in] params the parameters:
      ///   - verbose         [     1]: verbose trace mode 0: silent, 1: verbose.
      ///   - r-radius        [   3.0]: the constant for kernel radius parameter r in r(h)=r h^alpha (VCM,II,Trivial).;
      ///   - alpha           [  0.33]: the parameter alpha in r(h)=r h^alpha (VCM, II)."
      ///   - gridstep        [   1.0]: the digitization gridstep (often denoted by h).
      ///
      /// @return the vector containing the estimated principal curvatures and directions,
      ///  in the same order as \a surfels.
      ///
      /// @note The function is faster when surfels are in a specific order, as
      /// given for instance by a depth-first traversal (see @ref getSurfelRange)
      template <typename TPointPredicate>
      static CurvatureTensorQuantities
      getIIPrincipalCurvaturesAndDirections( const TPointPredicate&  shape,
                                          const KSpace&           K,
                                          const SurfelRange&      surfels,
                                          const Parameters&       params
                                            = parametersGeometryEstimation()
                                            | parametersKSpace() )
      {
        typedef functors::IIPrincipalCurvaturesAndDirectionsFunctor<Space> IICurvFunctor;
        typedef IntegralInvariantCovarianceEstimator<KSpace, TPointPredicate, IICurvFunctor>    IICurvEstimator;

        CurvatureTensorQuantities  mc_estimations;
        int      verbose = params[ "verbose"   ].as<int>();
        Scalar   h       = params[ "gridstep"  ].as<Scalar>();
        Scalar   r       = params[ "r-radius"  ].as<Scalar>();
        Scalar   alpha   = params[ "alpha"     ].as<Scalar>();
        if ( alpha != 1.0 ) r *= pow( h, alpha-1.0 );
        if ( verbose > 0 )
        {
          trace.info() << "- II principal curvatures and directions alpha=" << alpha << std::endl;
          trace.info() << "- II principal curvatures and directions r=" << (r*h)  << " (continuous) "
          << r << " (discrete)" << std::endl;
        }
        IICurvFunctor   functor;
        functor.init( h, r*h );
        IICurvEstimator ii_estimator( functor );
        ii_estimator.attach( K, shape );
        ii_estimator.setParams( r );
        ii_estimator.init( h, surfels.begin(), surfels.end() );
        ii_estimator.eval( surfels.begin(), surfels.end(),
                          std::back_inserter( mc_estimations ) );
        return mc_estimations;
      }

      /// @}

      // --------------------------- AT approximation ------------------------------
      /// @name AT approximation services
      /// @{
    public:

      /// @return the parameters and their default values which are used
      /// to compute Ambrosio-Tortorelli piecewise-smooth approximation of a function.
      ///   - at-enabled      [  1     ]: 1 if AT is enabled, 0 otherwise.
      ///   - at-alpha        [  0.1   ]: parameter alpha in AT (data fit)
      ///   - at-lambda       [  0.025 ]: parameter lambda in AT (1/length of discontinuities)
      ///   - at-epsilon      [  0.5   ]: (last value of) parameter epsilon in AT (width of discontinuities)
      ///   - at-epsilon-start[  2.0   ]: first value for parameter epsilon in Gamma-convergence optimization (sequence of AT optimization with decreasing epsilon)
      ///   - at-epsilon-ratio[  2.0   ]: ratio between two consecutive epsilon value in Gamma-convergence optimization (sequence of AT optimization with decreasing epsilon)
      ///   - at-max-iter     [ 10     ]: maximum number of alternate minization in AT optimization
      ///   - at-diff-v-max   [  0.0001]: stopping criterion that measures the loo-norm of the evolution of \a v between two iterations
      ///   - at-v-policy     ["Maximum"]: the policy when outputing feature vector v onto cells: "Average"|"Minimum"|"Maximum"
      ///
      static Parameters parametersATApproximation()
      {
        return Parameters
          ( "at-enabled",        1 )
          ( "at-alpha",          0.1 )
          ( "at-lambda",         0.025 )
          ( "at-epsilon",        0.25 )
          ( "at-epsilon-start",  2.0 )
          ( "at-epsilon-ratio",  2.0 )
          ( "at-max-iter",      10 )
          ( "at-diff-v-max",     0.0001 )
          ( "at-v-policy",   "Maximum" );
      }

      /// Given any digital \a surface, a surfel range \a surfels, and an input vector field \a input,
      /// returns a piece-smooth approximation of \a input using Ambrosio-Tortorelli functional.
      ///
      /// @see \ref moduleGenericAT
      ///
      /// @tparam TAnyDigitalSurface either kind of DigitalSurface, like ShortcutsGeometry::LightDigitalSurface or ShortcutsGeometry::DigitalSurface.
      /// @tparam VectorFieldInput the type of vector field for input values (RandomAccess container)
      ///
      /// @param[in] surface the digital surface
      /// @param[in] surfels the sequence of surfels at which we compute the normals
      /// @param[in] params the parameters:
      ///   - at-alpha        [  0.1   ]: parameter alpha in AT (data fit)
      ///   - at-lambda       [  0.025 ]: parameter lambda in AT (1/length of discontinuities)
      ///   - at-epsilon      [  0.5   ]: (last value of) parameter epsilon in AT (width of discontinuities)
      ///   - at-epsilon-start[  2.0   ]: first value for parameter epsilon in Gamma-convergence optimization (sequence of AT optimization with decreasing epsilon)
      ///   - at-epsilon-ratio[  2.0   ]: ratio between two consecutive epsilon value in Gamma-convergence optimization (sequence of AT optimization with decreasing epsilon)
      ///   - at-max-iter     [ 10     ]: maximum number of alternate minization in AT optimization
      ///   - at-diff-v-max   [  0.0001]: stopping criterion that measures the loo-norm of the evolution of \a v between two iterations
      /// @param[in] input the input vector field (a vector of vector values)
      ///
      /// @return the piecewise-smooth approximation of \a input.
      ///
      template <typename TAnyDigitalSurface,
                typename VectorFieldInput>
      static
      VectorFieldInput
      getATVectorFieldApproximation( CountedPtr<TAnyDigitalSurface> surface,
                                     const SurfelRange&             surfels,
                                     const VectorFieldInput&        input,
                                     const Parameters&              params
                                     = parametersATApproximation() | parametersGeometryEstimation() )
      {
        (void)surface; //param not used. FIXME: JOL

        int      verbose   = params[ "verbose"          ].as<int>();
        Scalar   alpha_at  = params[ "at-alpha"         ].as<Scalar>();
        Scalar   lambda_at = params[ "at-lambda"        ].as<Scalar>();
        Scalar   epsilon1  = params[ "at-epsilon-start" ].as<Scalar>();
        Scalar   epsilon2  = params[ "at-epsilon"       ].as<Scalar>();
        Scalar   epsilonr  = params[ "at-epsilon-ratio" ].as<Scalar>();
        int      max_iter  = params[ "at-max-iter"      ].as<int>();
        Scalar   diff_v_max= params[ "at-diff-v-max"    ].as<Scalar>();
        typedef DiscreteExteriorCalculusFactory<EigenLinearAlgebraBackend> CalculusFactory;
        const auto calculus = CalculusFactory::createFromNSCells<2>( surfels.cbegin(), surfels.cend() );
        ATSolver2D< KSpace > at_solver( calculus, verbose );
        at_solver.initInputVectorFieldU2( input, surfels.cbegin(), surfels.cend() );
        at_solver.setUp( alpha_at, lambda_at );
        at_solver.solveGammaConvergence( epsilon1, epsilon2, epsilonr, false, diff_v_max, max_iter );
        auto output = input;
        at_solver.getOutputVectorFieldU2( output, surfels.cbegin(), surfels.cend() );
        return output;
      }

      /// Given any digital \a surface, a surfel range \a surfels, and an input vector field \a input,
      /// returns a piece-smooth approximation of \a input using Ambrosio-Tortorelli functional.
      /// Given a range of pointels, linels or 2-cells [itB,itE), it
      /// also outputs the feature vector \a features, corresponding to
      /// 0-form \a v in AT (the average of \a v for linels/surfels).
      ///
      /// @see \ref moduleGenericAT
      ///
      /// @tparam TAnyDigitalSurface either kind of DigitalSurface, like ShortcutsGeometry::LightDigitalSurface or ShortcutsGeometry::DigitalSurface.
      /// @tparam VectorFieldInput the type of vector field for input values (RandomAccess container)
      /// @tparam CellRangeConstIterator the type of iterator for traversing a range of cells
      ///
      /// @param[out] features the vector of scalar feature values (a scalar field where 1 means continuity and 0 discontinuity in the reconstruction), evaluated in the range[itB,itE).
      /// @param[in] itB the start of the range of cells.
      /// @param[in] itE past the end of the range of cells.
      /// @param[in] surface the digital surface
      /// @param[in] surfels the sequence of surfels at which we compute the normals
      /// @param[in] params the parameters:
      ///   - at-alpha        [  0.1   ]: parameter alpha in AT (data fit)
      ///   - at-lambda       [  0.025 ]: parameter lambda in AT (1/length of discontinuities)
      ///   - at-epsilon      [  0.5   ]: (last value of) parameter epsilon in AT (width of discontinuities)
      ///   - at-epsilon-start[  2.0   ]: first value for parameter epsilon in Gamma-convergence optimization (sequence of AT optimization with decreasing epsilon)
      ///   - at-epsilon-ratio[  2.0   ]: ratio between two consecutive epsilon value in Gamma-convergence optimization (sequence of AT optimization with decreasing epsilon)
      ///   - at-max-iter     [ 10     ]: maximum number of alternate minization in AT optimization
      ///   - at-diff-v-max   [  0.0001]: stopping criterion that measures the loo-norm of the evolution of \a v between two iterations
      ///   - at-v-policy     ["Maximum"]: the policy when outputing feature vector v onto cells: "Average"|"Minimum"|"Maximum"
      /// @param[in] input the input vector field (a vector of vector values)
      ///
      /// @return the piecewise-smooth approximation of \a input.
      ///
      template <typename TAnyDigitalSurface,
                typename VectorFieldInput,
                typename CellRangeConstIterator>
      static
      VectorFieldInput
      getATVectorFieldApproximation( Scalars&                       features,
                                     CellRangeConstIterator         itB,
                                     CellRangeConstIterator         itE,
                                     CountedPtr<TAnyDigitalSurface> surface,
                                     const SurfelRange&             surfels,
                                     const VectorFieldInput&        input,
                                     const Parameters&              params
                                     = parametersATApproximation() | parametersGeometryEstimation() )
      {
        (void)surface; //param not used FIXME: JOL
        
        int      verbose   = params[ "verbose"          ].as<int>();
        Scalar   alpha_at  = params[ "at-alpha"         ].as<Scalar>();
        Scalar   lambda_at = params[ "at-lambda"        ].as<Scalar>();
        Scalar   epsilon1  = params[ "at-epsilon-start" ].as<Scalar>();
        Scalar   epsilon2  = params[ "at-epsilon"       ].as<Scalar>();
        Scalar   epsilonr  = params[ "at-epsilon-ratio" ].as<Scalar>();
        int      max_iter  = params[ "at-max-iter"      ].as<int>();
        Scalar   diff_v_max= params[ "at-diff-v-max"    ].as<Scalar>();
        std::string policy = params[ "at-v-policy"      ].as<std::string>();
        typedef DiscreteExteriorCalculusFactory<EigenLinearAlgebraBackend> CalculusFactory;
        const auto calculus = CalculusFactory::createFromNSCells<2>( surfels.cbegin(), surfels.cend() );
        ATSolver2D< KSpace > at_solver( calculus, verbose );
        at_solver.initInputVectorFieldU2( input, surfels.cbegin(), surfels.cend() );
        at_solver.setUp( alpha_at, lambda_at );
        at_solver.solveGammaConvergence( epsilon1, epsilon2, epsilonr, false, diff_v_max, max_iter );
        auto output = input;
        at_solver.getOutputVectorFieldU2( output, surfels.cbegin(), surfels.cend() );
        auto p = ( policy == "Average" ) ? at_solver.Average
          :      ( policy == "Minimum" ) ? at_solver.Minimum
          :                                at_solver.Maximum;
        at_solver.getOutputScalarFieldV0( features, itB, itE, p );
        return output;
      }

      /// Given any digital \a surface, a surfel range \a surfels, and
      /// an input scalar field \a input, returns a piece-smooth
      /// approximation of \a input using Ambrosio-Tortorelli
      /// functional.
      ///
      /// @see \ref moduleGenericAT
      ///
      /// @tparam TAnyDigitalSurface either kind of DigitalSurface, like ShortcutsGeometry::LightDigitalSurface or ShortcutsGeometry::DigitalSurface.
      ///
      /// @param[in] surface the digital surface
      /// @param[in] surfels the sequence of surfels at which we compute the normals
      /// @param[in] params the parameters:
      ///   - at-alpha        [  0.1   ]: parameter alpha in AT (data fit)
      ///   - at-lambda       [  0.025 ]: parameter lambda in AT (1/length of discontinuities)
      ///   - at-epsilon      [  0.5   ]: (last value of) parameter epsilon in AT (width of discontinuities)
      ///   - at-epsilon-start[  2.0   ]: first value for parameter epsilon in Gamma-convergence optimization (sequence of AT optimization with decreasing epsilon)
      ///   - at-epsilon-ratio[  2.0   ]: ratio between two consecutive epsilon value in Gamma-convergence optimization (sequence of AT optimization with decreasing epsilon)
      ///   - at-max-iter     [ 10     ]: maximum number of alternate minization in AT optimization
      ///   - at-diff-v-max   [  0.0001]: stopping criterion that measures the loo-norm of the evolution of \a v between two iterations
      /// @param[in] input the input scalar field (a vector of scalar values)
      ///
      /// @return the piecewise-smooth approximation of \a input.
      ///
      template <typename TAnyDigitalSurface>
      static
      Scalars
      getATScalarFieldApproximation( CountedPtr<TAnyDigitalSurface> surface,
                                     const SurfelRange&             surfels,
                                     const Scalars&                 input,
                                     const Parameters&              params
                                     = parametersATApproximation() | parametersGeometryEstimation() )
      {
        (void)surface; //param not used FIXME: JOL

        int      verbose   = params[ "verbose"          ].as<int>();
        Scalar   alpha_at  = params[ "at-alpha"         ].as<Scalar>();
        Scalar   lambda_at = params[ "at-lambda"        ].as<Scalar>();
        Scalar   epsilon1  = params[ "at-epsilon-start" ].as<Scalar>();
        Scalar   epsilon2  = params[ "at-epsilon"       ].as<Scalar>();
        Scalar   epsilonr  = params[ "at-epsilon-ratio" ].as<Scalar>();
        int      max_iter  = params[ "at-max-iter"      ].as<int>();
        Scalar   diff_v_max= params[ "at-diff-v-max"    ].as<Scalar>();
        typedef DiscreteExteriorCalculusFactory<EigenLinearAlgebraBackend> CalculusFactory;
        const auto calculus = CalculusFactory::createFromNSCells<2>( surfels.cbegin(), surfels.cend() );
        ATSolver2D< KSpace > at_solver( calculus, verbose );
        at_solver.initInputScalarFieldU2( input, surfels.cbegin(), surfels.cend() );
        at_solver.setUp( alpha_at, lambda_at );
        at_solver.solveGammaConvergence( epsilon1, epsilon2, epsilonr, false, diff_v_max, max_iter );
        auto output = input;
        at_solver.getOutputScalarFieldU2( output, surfels.cbegin(), surfels.cend() );
        return output;
      }

      /// Given any digital \a surface, a surfel range \a surfels, and
      /// an input scalar field \a input, returns a piece-smooth
      /// approximation of \a input using Ambrosio-Tortorelli
      /// functional.  Given a range of pointels, linels or 2-cells
      /// [itB,itE), it also outputs the feature vector \a features,
      /// corresponding to 0-form \a v in AT (the average of \a v for
      /// linels/surfels).
      ///
      /// @see \ref moduleGenericAT
      ///
      /// @tparam TAnyDigitalSurface either kind of DigitalSurface, like ShortcutsGeometry::LightDigitalSurface or ShortcutsGeometry::DigitalSurface.
      /// @tparam CellRangeConstIterator the type of iterator for traversing a range of cells
      ///
      /// @param[out] features the vector of scalar feature values (a
      /// scalar field where 1 means continuity and 0 discontinuity in
      /// the reconstruction), evaluated in the range[itB,itE).
      ///
      /// @param[in] itB the start of the range of cells.
      /// @param[in] itE past the end of the range of cells.
      /// @param[in] surface the digital surface
      /// @param[in] surfels the sequence of surfels at which we compute the normals
      /// @param[in] params the parameters:
      ///   - at-alpha        [  0.1   ]: parameter alpha in AT (data fit)
      ///   - at-lambda       [  0.025 ]: parameter lambda in AT (1/length of discontinuities)
      ///   - at-epsilon      [  0.5   ]: (last value of) parameter epsilon in AT (width of discontinuities)
      ///   - at-epsilon-start[  2.0   ]: first value for parameter epsilon in Gamma-convergence optimization (sequence of AT optimization with decreasing epsilon)
      ///   - at-epsilon-ratio[  2.0   ]: ratio between two consecutive epsilon value in Gamma-convergence optimization (sequence of AT optimization with decreasing epsilon)
      ///   - at-max-iter     [ 10     ]: maximum number of alternate minization in AT optimization
      ///   - at-diff-v-max   [  0.0001]: stopping criterion that measures the loo-norm of the evolution of \a v between two iterations
      ///   - at-v-policy     ["Maximum"]: the policy when outputing feature vector v onto cells: "Average"|"Minimum"|"Maximum"
      /// @param[in] input the input scalar field (a vector of scalar values)
      ///
      /// @return the piecewise-smooth approximation of \a input.
      ///
      template <typename TAnyDigitalSurface,
                typename CellRangeConstIterator>
      static
      Scalars
      getATScalarFieldApproximation( Scalars&                       features,
                                     CellRangeConstIterator         itB,
                                     CellRangeConstIterator         itE,
                                     CountedPtr<TAnyDigitalSurface> surface,
                                     const SurfelRange&             surfels,
                                     const Scalars&                 input,
                                     const Parameters&              params
                                     = parametersATApproximation() | parametersGeometryEstimation() )
      {
        (void)surface; //param not used FIXME: JOL
        
        int      verbose   = params[ "verbose"          ].as<int>();
        Scalar   alpha_at  = params[ "at-alpha"         ].as<Scalar>();
        Scalar   lambda_at = params[ "at-lambda"        ].as<Scalar>();
        Scalar   epsilon1  = params[ "at-epsilon-start" ].as<Scalar>();
        Scalar   epsilon2  = params[ "at-epsilon"       ].as<Scalar>();
        Scalar   epsilonr  = params[ "at-epsilon-ratio" ].as<Scalar>();
        int      max_iter  = params[ "at-max-iter"      ].as<int>();
        Scalar   diff_v_max= params[ "at-diff-v-max"    ].as<Scalar>();
        std::string policy = params[ "at-v-policy"      ].as<std::string>();
        typedef DiscreteExteriorCalculusFactory<EigenLinearAlgebraBackend> CalculusFactory;
        const auto calculus = CalculusFactory::createFromNSCells<2>( surfels.cbegin(), surfels.cend() );
        ATSolver2D< KSpace > at_solver( calculus, verbose );
        at_solver.initInputScalarFieldU2( input, surfels.cbegin(), surfels.cend() );
        at_solver.setUp( alpha_at, lambda_at );
        at_solver.solveGammaConvergence( epsilon1, epsilon2, epsilonr, false, diff_v_max, max_iter );
        auto output = input;
        at_solver.getOutputScalarFieldU2( output, surfels.cbegin(), surfels.cend() );
        auto p = ( policy == "Average" ) ? at_solver.Average
          :      ( policy == "Minimum" ) ? at_solver.Minimum
          :                                at_solver.Maximum;
        at_solver.getOutputScalarFieldV0( features, itB, itE, p );
        return output;
      }

      /// @}

      // ------------------------- Error measures services -------------------------
      /// @name Error measure services
      /// @{
    public:

      /// Orient \a v so that it points in the same direction as \a
      /// ref_v (scalar product is then non-negative afterwards).
      ///
      /// @param[in,out] v the vectors to reorient.
      /// @param[in]    ref_v the vectors having the reference orientation.
      static void
        orientVectors( RealVectors&       v,
                       const RealVectors& ref_v )
      {
        std::transform( ref_v.cbegin(), ref_v.cend(), v.cbegin(), v.begin(),
                        [] ( RealVector rw, RealVector w )
                        { return rw.dot( w ) >= 0.0 ? w : -w; } );
      }

      /// Computes the statistic of a vector of scalars
      ///
      /// @param[in] v a vector of scalars
      /// @return its statistic.
      static ScalarStatistic
        getStatistic( const Scalars& v )
      {
        ScalarStatistic stat;
        stat.addValues( v.begin(), v.end() );
        stat.terminate();
        return stat;
      }

      /// Computes the statistic that measures the angle differences
      /// between the two arrays of unit vectors.
      ///
      /// @param[in] v1 the first array of unit vectors (normals)
      /// @param[in] v2 the second array of unit vectors (normals)
      /// @return the vector of angle differences.
      static Scalars
        getVectorsAngleDeviation( const RealVectors& v1,
                                  const RealVectors& v2 )
      {
        Scalars v( v1.size() );
        if ( v1.size() == v2.size() )
          {
            auto outIt = v.begin();
            for ( auto it1 = v1.cbegin(), it2 = v2.cbegin(), itE1 = v1.cend();
                  it1 != itE1; ++it1, ++it2 )
              {
                Scalar angle_error = acos( (*it1).dot( *it2 ) );
                *outIt++ = angle_error;
              }
          }
        else
          {
            trace.warning() << "[ShortcutsGeometry::getVectorsAngleDeviation]"
                            << " v1.size()=" << v1.size() << " should be equal to "
                            << " v2.size()=" << v2.size() << std::endl;
          }
        return v;
      }

      /// Computes the absolute difference between each element of the two vectors.
      /// @param[in] v1 any vector of values.
      /// @param[in] v2 any vector of values.
      /// @return the vector composed of elemenst |v1[i]-v2[i]|.
      static Scalars
        getScalarsAbsoluteDifference( const Scalars & v1,
                                      const Scalars & v2 )
      {
        Scalars result( v1.size() );
        std::transform( v2.cbegin(), v2.cend(), v1.cbegin(), result.begin(),
                        [] ( Scalar val1, Scalar val2 )
                        { return fabs( val1 - val2 ); } );
        return result;
      }

      /// Computes the l2-norm of v1-v2, ie the square root of the
      /// mean-squared error of the two vectors.
      ///
      /// @param[in] v1 any vector of values.
      /// @param[in] v2 any vector of values.
      /// @return the normL2 of v1-v2, ie. sqrt( 1/n sum_i (v1[i]-v2[i])^2 ).
      static Scalar
        getScalarsNormL2( const Scalars & v1,
                          const Scalars & v2 )
      {
        Scalar sum = 0;
        for ( unsigned int i = 0; i < v1.size(); i++ )
          sum += ( v1[ i ] - v2[ i ] ) * ( v1[ i ] - v2[ i ] );
        return sqrt( sum / v1.size() );
      }

      /// Computes the l1-norm of v1-v2, ie the average of the absolute
      /// differences of the two vectors.
      ///
      /// @param[in] v1 any vector of values.
      /// @param[in] v2 any vector of values.
      /// @return the normL1 of v1-v2, ie. 1/n sum_i |v1[i]-v2[i]|.
      static Scalar
        getScalarsNormL1( const Scalars & v1,
                          const Scalars & v2 )
      {
        Scalar sum = 0;
        for ( unsigned int i = 0; i < v1.size(); i++ )
          sum += fabs( v1[ i ] - v2[ i ] );
        return sum / v1.size();
      }

      /// Computes the loo-norm of v1-v2, ie the maximum of the absolute
      /// differences of the two vectors.
      ///
      /// @param[in] v1 any vector of values.
      /// @param[in] v2 any vector of values.
      /// @return the normLoo of v1-v2, ie. max_i |v1[i]-v2[i]|.
      static Scalar
        getScalarsNormLoo( const Scalars & v1,
                           const Scalars & v2 )
      {
        Scalar loo = 0;
        for ( unsigned int i = 0; i < v1.size(); i++ )
          loo = std::max( loo, fabs( v1[ i ] - v2[ i ] ) );
        return loo;
      }
      /// @}
      
      // ----------------------- VoronoiMap services ------------------------------
    public:
      /// @name VoronoiMap services
      /// @{

      /// @return the parameters and their default values which are used
      /// in VoronoiMap and DistanceTransformation
      //    - toroidal-x [false]: If the domain is toroidal in the first  dimension
      //    - toroidal-y [false]: If the domain is toroidal in the second dimension
      //    - toroidal-z [false]: If the domain is toroidal in the third  dimension
      static Parameters parametersVoronoiMap() {
        return Parameters
          // Toricity might be moved elsewhere as this is quite a general parameter
          ( "toroidal-x" , false )
          ( "toroidal-y" , false )
          ( "toroidal-z" , false );
      }



      /// @brief Computes the VoronoiMap of a domain, where sites are given through a range.
      ///
      /// @tparam p The exponent in the Lp metric
      /// @tparam PointRange An iterable of points (std::vector, DGtal::DigitalSet*, ...)
      ///
      /// @param domain The associated space to compute VoronoiMap on
      /// @param sites The list of sites
      /// @param params the parameters
      //    - toroidal-x [false]: If the domain is toroidal in the first  dimension
      //    - toroidal-y [false]: If the domain is toroidal in the second dimension
      //    - toroidal-z [false]: If the domain is toroidal in the third  dimension
      /// 
      /// @return The VoronoiMap within a domain with prescribed sites
      template<uint32_t p, typename PointRange>
      static VoronoiMap<Space, VoronoiPointPredicate, ExactPredicateLpSeparableMetric<Space, p>>
        getVoronoiMap(Domain domain, 
                      const PointRange& sites,
                      const Parameters& params = parametersVoronoiMap())
      {
        using Metric = ExactPredicateLpSeparableMetric<Space, p>;
        using Map = VoronoiMap<Space, VoronoiPointPredicate, Metric>;
        DigitalSet set(domain); set.insert(sites.begin(), sites.end());
        VoronoiPointPredicate predicate(set);
        Metric metric;

        typename Map::PeriodicitySpec specs = {false, false, false};
        if (params["toroidal-x"].as<int>()) specs[0] = true;
        if (params["toroidal-y"].as<int>()) specs[1] = true;
        if (params["toroidal-z"].as<int>()) specs[2] = true;


        // Do not return a pointer here for two reasons:
        //  - The distance transform will not be passed anywhere else
        //  - The operator() is less accessible with pointers.
        return Map(domain, predicate, metric, specs);
      }

      /// @brief Computes the VoronoiMap of a domain, where sites are given through a range.
      ///
      /// @tparam p The exponent in the Lp metric
      /// @tparam PointRange An iterable of points (std::vector, DGtal::DigitalSet*, ...)
      ///
      /// @param domain The associated space to compute VoronoiMap on
      /// @param sites The list of sites
      /// @param params the parameters
      //    - toroidal-x [false]: If the domain is toroidal in the first  dimension
      //    - toroidal-y [false]: If the domain is toroidal in the second dimension
      //    - toroidal-z [false]: If the domain is toroidal in the third  dimension
      /// 
      /// @return The VoronoiMap within a domain with prescribed sites
      template<uint32_t p, typename PointRange>
      static VoronoiMap<Space, VoronoiPointPredicate, ExactPredicateLpSeparableMetric<Space, p>>
        getVoronoiMap(CountedPtr<Domain> domain, 
                      const PointRange& sites,
                      const Parameters& params = parametersVoronoiMap())
      {
        using Metric = ExactPredicateLpSeparableMetric<Space, p>;
        using Map = VoronoiMap<Space, VoronoiPointPredicate, Metric>;
        DigitalSet set(*domain); set.insert(sites.begin(), sites.end());
        VoronoiPointPredicate predicate(set);
        Metric metric;

        typename Map::PeriodicitySpec specs = {false, false, false};
        if (params["toroidal-x"].as<int>()) specs[0] = true;
        if (params["toroidal-y"].as<int>()) specs[1] = true;
        if (params["toroidal-z"].as<int>()) specs[2] = true;

        // Do not return a pointer here for two reasons:
        //  - The distance transform will not be passed anywhere else
        //  - The operator() is less accessible with pointers.
        return Map(*domain, predicate, metric, specs);
      }

      /// @brief Computes the VoronoiMap of a domain, where sites are given through a range.
      ///
      /// @note: This overloads return a distance transformation, ie. where operator() returns
      /// the distance to the closest site.
      ///
      /// @tparam p The exponent in the Lp metric
      /// @tparam PointRange An iterable of points (std::vector, DGtal::DigitalSet*, ...)
      ///
      /// @param domain The associated space to compute VoronoiMap on
      /// @param sites The list of sites
      /// @param params the parameters
      //    - toroidal-x [false]: If the domain is toroidal in the first  dimension
      //    - toroidal-y [false]: If the domain is toroidal in the second dimension
      //    - toroidal-z [false]: If the domain is toroidal in the third  dimension
      /// 
      /// @return The DistanceTransformation within a domain with prescribed sites
      template<uint32_t p, typename PointRange>
      static DistanceTransformation<Space, VoronoiPointPredicate, ExactPredicateLpSeparableMetric<Space, p>>
        getDistanceTransformation(Domain domain, 
                                   const PointRange& sites,
                                   const Parameters& params = parametersVoronoiMap())
      {
        using Metric = ExactPredicateLpSeparableMetric<Space, p>;
        using Map = DistanceTransformation<Space, VoronoiPointPredicate, Metric>;
        DigitalSet set(domain); set.insert(sites.begin(), sites.end());
        VoronoiPointPredicate predicate(set);
        Metric metric;

        typename Map::PeriodicitySpec specs = {false, false, false};
        if (params["toroidal-x"].as<int>()) specs[0] = true;
        if (params["toroidal-y"].as<int>()) specs[1] = true;
        if (params["toroidal-z"].as<int>()) specs[2] = true;

        // Do not return a pointer here for two reasons:
        //  - The distance transform will not be passed anywhere else
        //  - The operator() is less accessible with pointers.
        return Map(domain, predicate, metric, specs);
      }

      /// @brief Computes the VoronoiMap of a domain, where sites are given through a range.
      ///
      /// @note: This overloads return a distance transformation, ie. where operator() returns
      /// the distance to the closest site.
      ///
      /// @tparam p The exponent in the Lp metric
      /// @tparam PointRange An iterable of points (std::vector, DGtal::DigitalSet*, ...)
      ///
      /// @param domain The associated space to compute VoronoiMap on
      /// @param sites The list of sites
      /// @param params the parameters
      //    - toroidal-x [false]: If the domain is toroidal in the first  dimension
      //    - toroidal-y [false]: If the domain is toroidal in the second dimension
      //    - toroidal-z [false]: If the domain is toroidal in the third  dimension
      /// 
      /// @return The DistanceTransformation within a domain with prescribed sites
      template<uint32_t p, typename PointRange>
      static DistanceTransformation<Space, VoronoiPointPredicate, ExactPredicateLpSeparableMetric<Space, p>>
        getDistanceTransformation(CountedPtr<Domain> domain, 
                              const PointRange& sites,
                              const Parameters& params = parametersVoronoiMap())
      {
        using Metric = ExactPredicateLpSeparableMetric<Space, p>;
        using Map = DistanceTransformation<Space, VoronoiPointPredicate, Metric>;
        DigitalSet set(*domain); set.insert(sites.begin(), sites.end());
        VoronoiPointPredicate predicate(set);
        Metric metric;

        typename Map::PeriodicitySpec specs = {false, false, false};
        if (params["toroidal-x"].as<int>()) specs[0] = true;
        if (params["toroidal-y"].as<int>()) specs[1] = true;
        if (params["toroidal-z"].as<int>()) specs[2] = true;
        
        // Do not return a pointer here for two reasons:
        //  - The distance transform will not be passed anywhere else
        //  - The operator() is less accessible with pointers.
        return Map(*domain, predicate, metric, specs);
      }

      /// @brief Computes the direction of the closest site of a range of points
      /// 
      /// @tparam p The exponent in the Lp metric
      /// @tparam PointRange The range of point 
      /// @tparam PointRangeSites The range of sites
      ///
      /// @param points The one to compute the closest site of
      /// @param sites The list of sites
      /// @param params Parameters
      //    - toroidal-x [false]: If the domain is toroidal in the first  dimension
      //    - toroidal-y [false]: If the domain is toroidal in the second dimension
      //    - toroidal-z [false]: If the domain is toroidal in the third  dimension
      /// 
      /// @return A vector of direction to the closest in the same order as 'points'.
      template<uint32_t p, typename PointRangeSites, typename PointRange>
      static std::vector<Vector> getDirectionToClosestSite(
        const PointRange& points, 
        const PointRangeSites& sites,
        const Parameters& params = parametersVoronoiMap())
      {
        using Metric = ExactPredicateLpSeparableMetric<Space, p>;
        using Map = VoronoiMap<Space, VoronoiPointPredicate, Metric>;

        // Compute domain of points
        Point pmin = *points.begin();
        Point pmax = pmin;

        size_t pCount = 0;
        for (auto it = points.begin(); it != points.end(); ++it) 
        {
          pCount ++;
          for (size_t i = 0; i < Space::dimension; ++i)
          {
            pmin[i] = std::min(pmin[i], (*it)[i] - 1);
            pmax[i] = std::max(pmax[i], (*it)[i] + 1);
          }
        }

        for (auto it = sites.begin(); it != sites.end(); ++it) 
        {
          for (size_t i = 0; i < Space::dimension; ++i)
          {
            pmin[i] = std::min(pmin[i], (*it)[i] - 1);
            pmax[i] = std::max(pmax[i], (*it)[i] + 1);
          }
        }

        Domain domain(pmin, pmax);

        DigitalSet set(domain); set.insert(sites.begin(), sites.end());
        VoronoiPointPredicate predicate(set);
        Metric metric;


        typename Map::PeriodicitySpec specs = {false, false, false};
        if (params["toroidal-x"].as<int>()) specs[0] = true;
        if (params["toroidal-y"].as<int>()) specs[1] = true;
        if (params["toroidal-z"].as<int>()) specs[2] = true;

        auto map = Map(domain, predicate, metric, specs);

        std::vector<Vector> directions(pCount);
        size_t i = 0;
        for (auto it = points.begin(); it != points.end(); ++it)
        {
          directions[i++] = map(*it);
        }
        return directions;
      }

      /// @brief Computes the distance of the closest site of a range of points
      /// 
      /// @tparam p The exponent in the Lp metric
      /// @tparam PointRange The range of point 
      /// @tparam PointRangeSites The range of sites
      ///
      /// @param points The one to compute the closest site of
      /// @param sites The list of sites
      /// @param params Parameters
      //    - toroidal-x [false]: If the domain is toroidal in the first  dimension
      //    - toroidal-y [false]: If the domain is toroidal in the second dimension
      //    - toroidal-z [false]: If the domain is toroidal in the third  dimension
      /// 
      /// @return A vector of distances to the closest in the same order as 'points'.
      template<uint32_t p, typename PointRangeSites, typename PointRange>
      static std::vector<typename ExactPredicateLpSeparableMetric<Space, p>::Value> getDistanceToClosestSite(
        const PointRange& points, 
        const PointRangeSites& sites,
        const Parameters& params = parametersVoronoiMap())
      {
        using Metric = ExactPredicateLpSeparableMetric<Space, p>;
        using Map = DistanceTransformation<Space, VoronoiPointPredicate, Metric>;

        // Compute domain of points
        Point pmin = *points.begin();
        Point pmax = pmin;

        size_t pCount = 0;
        for (auto it = points.begin(); it != points.end(); ++it) 
        {
          pCount ++;
          for (size_t i = 0; i < Space::dimension; ++i)
          {
            pmin[i] = std::min(pmin[i], (*it)[i] - 1);
            pmax[i] = std::max(pmax[i], (*it)[i] + 1);
          }
        }

        for (auto it = sites.begin(); it != sites.end(); ++it) 
        {
          for (size_t i = 0; i < Space::dimension; ++i)
          {
            pmin[i] = std::min(pmin[i], (*it)[i] - 1);
            pmax[i] = std::max(pmax[i], (*it)[i] + 1);
          }
        }

        Domain domain(pmin, pmax);

        DigitalSet set(domain); set.insert(sites.begin(), sites.end());
        VoronoiPointPredicate predicate(set);
        Metric metric;

        typename Map::PeriodicitySpec specs = {false, false, false};
        if (params["toroidal-x"].as<int>()) specs[0] = true;
        if (params["toroidal-y"].as<int>()) specs[1] = true;
        if (params["toroidal-z"].as<int>()) specs[2] = true;

        auto map = Map(domain, predicate, metric, specs);

        std::vector<typename Metric::Value> directions(pCount);
        size_t i = 0;
        for (auto it = points.begin(); it != points.end(); ++it)
        {
          directions[i++] = map(*it);
        }
        return directions;
      }


      /// @}

      // ----------------------- Standard services ------------------------------
      /// @name Standard services
      /// @{
    public:

      /**
       * Default constructor.
       */
      ShortcutsGeometry() = delete;

      /**
       * Destructor.
       */
      ~ShortcutsGeometry() = delete;

      /**
       * Copy constructor.
       * @param other the object to clone.
       */
      ShortcutsGeometry ( const ShortcutsGeometry & other ) = delete;

      /**
       * Move constructor.
       * @param other the object to move.
       */
      ShortcutsGeometry ( ShortcutsGeometry && other ) = delete;

      /**
       * Copy assignment operator.
       * @param other the object to copy.
       * @return a reference on 'this'.
       */
      ShortcutsGeometry & operator= ( const ShortcutsGeometry & other ) = delete;

      /**
       * Move assignment operator.
       * @param other the object to move.
       * @return a reference on 'this'.
       */
      ShortcutsGeometry & operator= ( ShortcutsGeometry && other ) = delete;

      /// @}

      // ----------------------- Interface --------------------------------------
    public:

      // ------------------------- Protected Datas ------------------------------
    protected:

      // ------------------------- Private Datas --------------------------------
    private:

      // ------------------------- Hidden services ------------------------------
    protected:

      // ------------------------- Internals ------------------------------------
    private:

    }; // end of class ShortcutsGeometry


} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined ShortcutsGeometry_h

#undef ShortcutsGeometry_RECURSES
#endif // else defined(ShortcutsGeometry_RECURSES)

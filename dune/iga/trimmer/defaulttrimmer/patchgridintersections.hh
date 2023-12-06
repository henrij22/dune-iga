// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#pragma once

#include <variant>
/** \file
 * @brief The TrimmedPatchGridLeafIntersection and TrimmedLevelIntersection classes
 */

namespace Dune::IGANEW::DefaultTrim {

  // External forward declarations
  template <class Grid>
  struct HostGridAccess;

  /** @brief An intersection with a leaf neighbor element
   * \ingroup PatchGrid
   * Mesh entities of codimension 0 ("elements") allow to visit all neighbors, where
   * a neighbor is an entity of codimension 0 which has a common entity of codimension 1
   * These neighbors are accessed via a IntersectionIterator. This allows the implement
   * non-matching meshes. The number of neighbors may be different from the number
   * of an element!
   */
  template <class GridImp>
  class TrimmedLeafIntersection {
    friend  typename GridImp::Traits::LeafIntersectionIterator;


    friend struct HostGridAccess<typename std::remove_const<GridImp>::type>;

    constexpr static int dim   = GridImp::dimension;
    constexpr static int mydim = GridImp::dimension - 1;

    constexpr static int dimworld = GridImp::dimensionworld;
    using Trimmer             = typename GridImp::Trimmer;



    using HostLeafIntersection = typename GridImp::Trimmer::TrimmerTraits::HostLeafIntersection;


   public:
    // The type used to store coordinates
    typedef typename GridImp::ctype ctype;
    using LocalCoordinate = FieldVector<ctype, mydim>;

    typedef typename GridImp::Trimmer::TrimmerTraits::template Codim<1>::LocalParameterSpaceGeometry Geometry;
    typedef typename GridImp::Trimmer::TrimmerTraits::template Codim<1>::LocalGeometry LocalGeometry;
    // typedef typename GridImp::template Codim<0>::Entity Entity;
    typedef FieldVector<ctype, dim> NormalVector;

    using ParameterSpaceGridEntity= typename GridImp::Trimmer::TrimmerTraits::template Codim<0>::ParameterSpaceGridEntity;

    using IntersectionGeometry = GeometryKernel::NURBSPatch<mydim,dim,ctype>;

    TrimmedLeafIntersection() {}

    TrimmedLeafIntersection(const GridImp* parameterSpaceGrid, const HostLeafIntersection& hostIntersection)
        : patchGrid_(parameterSpaceGrid) ,hostIntersection_{hostIntersection}{}

    TrimmedLeafIntersection(const GridImp* parameterSpaceGrid, HostLeafIntersection&& hostIntersection)
        : patchGrid_(parameterSpaceGrid) ,hostIntersection_{hostIntersection}{}
    HostLeafIntersection  hostIntersection_;
    bool operator==(const TrimmedLeafIntersection& other) const {
      // DUNE_THROW(NotImplemented, "equals not implemented");
      return hostIntersection_ == other.hostIntersection_;
    }

    //! returns the inside entity
    ParameterSpaceGridEntity inside() const {
      DUNE_THROW(NotImplemented, "inside not implemented");
      return ParameterSpaceGridEntity(patchGrid_,hostIntersection_.inside(),{});
    }

    //! return Entity on the outside of this intersection
    //! (that is the neighboring Entity)
    ParameterSpaceGridEntity outside() const {
      DUNE_THROW(NotImplemented, "outside not implemented");
      return ParameterSpaceGridEntity(patchGrid_,hostIntersection_.outside(),{});

    }

    //! return true if intersection is with boundary.
    [[nodiscard]] bool boundary() const {
      DUNE_THROW(NotImplemented, "boundary not implemented");
    }

    /** @brief Return unit outer normal (length == 1)
     *
     *   The returned vector is the normal at the center() of the
     *     intersection's geometry.
     *       It is scaled to have unit length. */
    NormalVector centerUnitOuterNormal() const {
      DUNE_THROW(NotImplemented, "centerUnitOuterNormal not implemented");
      //! @todo compute jacobian create normal by cross-product, the cross-product has to take into account the normal points away from the inside element
    }

    //! return true if across the edge an neighbor on this level exists
    bool neighbor() const {
      DUNE_THROW(NotImplemented, "neighbor not implemented");
      return {};
    }

    //! return the boundary segment index
    size_t boundarySegmentIndex() const {
      DUNE_THROW(NotImplemented, "boundarySegmentIndex not implemented");
      return {};
    }

    //! Return true if this is a conforming intersection
    bool conforming() const {
      DUNE_THROW(NotImplemented, "conforming not implemented");
      return {};
    }

    //! Geometry type of an intersection
    GeometryType type() const {
      DUNE_THROW(NotImplemented, "type not implemented");
      return {};
    }

    //! @todo this function should return how this inersection resides in the inside host element.
    //! Therefore, this function should provide this intersection and maps to the 0..1 space
    LocalGeometry geometryInInside() const {
      DUNE_THROW(NotImplemented, "geometryInInside not implemented");

      return std::variant_alternative_t<1,typename LocalGeometry::Variant>{};
    }

    //! Same as above
    LocalGeometry geometryInOutside() const {
      DUNE_THROW(NotImplemented, "geometryInOutside not implemented");

      return std::variant_alternative_t<1,typename LocalGeometry::Variant>{};
    }

    //! geometry of the intersection this geometry should map into the knotspan domain
    using TrimmedParameterSpaceGeometry= typename GridImp::Trimmer::TrimmerTraits::template Codim<1>::TrimmedParameterSpaceGeometry;

    TrimmedParameterSpaceGeometry geometry() const {
      // @todo trim this will be wrong as soon as the intersection geometry has a special geoemtry

      return {};
    }

    //! local number of codim 1 entity in self where intersection is contained in
    int indexInInside() const {
      DUNE_THROW(NotImplemented, "indexInInside not implemented");
    }

    //! local number of codim 1 entity in neighbor where intersection is contained
    int indexInOutside() const {
      DUNE_THROW(NotImplemented, "indexInOutside not implemented");
    }

    //! return outer normal
    FieldVector<ctype, GridImp::dimensionworld> outerNormal(const LocalCoordinate& local) const {
      DUNE_THROW(NotImplemented, "outerNormal not implemented");

      return {};
    }

    //! return outer normal multiplied by the integration element
    FieldVector<ctype, GridImp::dimensionworld> integrationOuterNormal(const LocalCoordinate& local) const {
      DUNE_THROW(NotImplemented, "integrationOuterNormal not implemented");

      return {};
    }

    //! return unit outer normal
    FieldVector<ctype, GridImp::dimensionworld> unitOuterNormal(const LocalCoordinate& local) const {
      DUNE_THROW(NotImplemented, "unitOuterNormal not implemented");

      return {};
    }

   private:
    //**********************************************************
    //  private methods
    //**********************************************************

    const GridImp* patchGrid_{nullptr};
    IntersectionGeometry geo;
  };

  template <class GridImp>
  class TrimmedLevelIntersection {
    friend  typename GridImp::Traits::LevelIntersectionIterator;

    friend struct HostGridAccess<typename std::remove_const<GridImp>::type>;

    constexpr static int dim   = GridImp::dimension;
    constexpr static int mydim = GridImp::dimension - 1;

    constexpr static int dimworld = GridImp::dimensionworld;

    using Trimmer = typename GridImp::Trimmer;

    using ParameterSpaceGridEntity= typename GridImp::Trimmer::TrimmerTraits::template Codim<0>::ParameterSpaceGridEntity;
    typedef typename GridImp::Trimmer::TrimmerTraits::template Codim<1>::LocalParameterSpaceGeometry Geometry;
    typedef typename GridImp::Trimmer::TrimmerTraits::template Codim<1>::LocalGeometry LocalGeometry;

    // typedef typename GridImp::GridFamily::LevelIntersection HostLevelIntersection;

    using HostLevelIntersection = typename GridImp::Trimmer::TrimmerTraits::HostLevelIntersection;


    using MatrixHelper = typename MultiLinearGeometryTraits<double>::MatrixHelper;

   public:
    // The type used to store coordinates
    typedef typename GridImp::ctype ctype;
    using LocalCoordinate = FieldVector<ctype, mydim>;

    // typedef typename GridImp::template Codim<1>::Geometry Geometry;
    // typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;
    // typedef typename GridImp::template Codim<0>::Entity Entity;
    typedef FieldVector<ctype, dimworld> NormalVector;

    TrimmedLevelIntersection() = default;

    TrimmedLevelIntersection(const GridImp* identityGrid, const HostLevelIntersection& hostIntersection)
        : patchGrid_(identityGrid),hostIntersection_{hostIntersection}{}

    TrimmedLevelIntersection(const GridImp* identityGrid, HostLevelIntersection&& hostIntersection)
        : patchGrid_(identityGrid),hostIntersection_{hostIntersection}{}

    HostLevelIntersection hostIntersection_;
    [[nodiscard]] bool operator==(const TrimmedLevelIntersection& other) const {
      // DUNE_THROW(NotImplemented, "equals not implemented");

      return hostIntersection_ == other.hostIntersection_;
    }

    //! return Entity on the inside of this intersection
    //! (that is the Entity where we started this Iterator)
    [[nodiscard]] ParameterSpaceGridEntity inside() const {
      DUNE_THROW(NotImplemented, "inside not implemented");

      return ParameterSpaceGridEntity(patchGrid_,hostIntersection_.inside(),{});
    }

    //! return Entity on the outside of this intersection
    //! (that is the neighboring Entity)
    [[nodiscard]] ParameterSpaceGridEntity outside() const {
      DUNE_THROW(NotImplemented, "outside not implemented");

      return ParameterSpaceGridEntity(patchGrid_,hostIntersection_.outside(),{});
    }

    /** @brief return true if intersection is with boundary.
     */
    [[nodiscard]] bool boundary() const {
      DUNE_THROW(NotImplemented, "boundary not implemented");

      return {};
    }

    /** @brief Return unit outer normal (length == 1)
     *
     *   The returned vector is the normal at the center() of the
     *     intersection's geometry.
     *       It is scaled to have unit length. */
    NormalVector centerUnitOuterNormal() const {
      DUNE_THROW(NotImplemented, "centerUnitOuterNormal not implemented");

      return {};
    }

    //! return true if across the edge an neighbor on this level exists
    [[nodiscard]] bool neighbor() const {
      DUNE_THROW(NotImplemented, "neighbor not implemented");

      return {};
    }

    //! return the boundary segment index
    [[nodiscard]] size_t boundarySegmentIndex() const {
      DUNE_THROW(NotImplemented, "boundarySegmentIndex not implemented");

      return {};
    }

    //! Return true if this is a conforming intersection
    [[nodiscard]] bool conforming() const {
      DUNE_THROW(NotImplemented, "conforming not implemented");

      return {};
    }

    //! Geometry type of an intersection
    [[nodiscard]] GeometryType type() const {
      DUNE_THROW(NotImplemented, "type not implemented");

      return {};
    }

    //! intersection of codimension 1 of this neighbor with element where
    //! iteration started.
    //! Here returned element is in LOCAL coordinates of the element
    //! where iteration started.
    [[nodiscard]] LocalGeometry geometryInInside() const {
      DUNE_THROW(NotImplemented, "geometryInInside not implemented");

      return std::variant_alternative_t<1,typename LocalGeometry::Variant>{};
    }

    //! intersection of codimension 1 of this neighbor with element where iteration started.
    //! Here returned element is in LOCAL coordinates of neighbor
    [[nodiscard]] LocalGeometry geometryInOutside() const {
      DUNE_THROW(NotImplemented, "geometryInOutside not implemented");

      return std::variant_alternative_t<1,typename LocalGeometry::Variant>{};
    }

    //! intersection of codimension 1 of this neighbor with element where iteration started.
    //! Here returned element is in GLOBAL coordinates of the element where iteration started.
    using TrimmedParameterSpaceGeometry= typename GridImp::Trimmer::TrimmerTraits::template Codim<1>::TrimmedParameterSpaceGeometry;
    [[nodiscard]] TrimmedParameterSpaceGeometry geometry() const {
      DUNE_THROW(NotImplemented, "geometry not implemented");

      return {};
    }

    //! local number of codim 1 entity in self where intersection is contained in
    [[nodiscard]] int indexInInside() const {
      DUNE_THROW(NotImplemented, "indexInInside not implemented");

      return {};
    }

    //! local number of codim 1 entity in neighbor where intersection is contained
    [[nodiscard]] int indexInOutside() const {
      DUNE_THROW(NotImplemented, "indexInOutside not implemented");

      return {};
    }

    //! return outer normal
    [[nodiscard]] FieldVector<ctype, dimworld> outerNormal(const LocalCoordinate& local) const {
      DUNE_THROW(NotImplemented, "outerNormal not implemented");

      return {};
    }

    //! return outer normal multiplied by the integration element
    [[nodiscard]] FieldVector<ctype, dimworld> integrationOuterNormal(const LocalCoordinate& local) const {
      DUNE_THROW(NotImplemented, "integrationOuterNormal not implemented");

      return {};
    }

    //! return unit outer normal
    [[nodiscard]] FieldVector<ctype, dimworld> unitOuterNormal(const LocalCoordinate& local) const {
      DUNE_THROW(NotImplemented, "unitOuterNormal not implemented");

      return {};
    }

   private:
    const GridImp* patchGrid_;

    // HostLevelIntersection parameterSpaceIntersection;
  };

}  // namespace Dune::IGANEW

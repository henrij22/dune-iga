// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#pragma once
#include "patchgridintersections.hh"

/** \file
 * @brief The TrimmedPatchGridLeafIntersection and TrimmedLevelIntersection classes
 */

// @todo same as gridintersection, one variant that forwards everything to hostEntity and one trimmed one
// @todo Also make this for leaf and level the same

namespace Dune::IGANEW::DefaultTrim {

// External forward declarations
template <class Grid>
struct HostGridAccess;

namespace Impl {

  enum class IntersectionType
  {
    Level,
    Leaf
  };

  namespace IntersectionTraits {

    template <class GridImp, IntersectionType type_>
    using HostIntersection = std::conditional_t<type_ == IntersectionType::Leaf,
                                                typename GridImp::Trimmer::TrimmerTraits::HostLeafIntersection,
                                                typename GridImp::Trimmer::TrimmerTraits::HostLevelIntersection>;

    template <class GridImp>
    using Geometry = typename GridImp::Trimmer::TrimmerTraits::template Codim<1>::LocalParameterSpaceGeometry;

    template <class GridImp>
    using LocalGeometry = typename GridImp::Trimmer::TrimmerTraits::template Codim<1>::LocalGeometry;

    template <class GridImp>
    using ParameterSpaceGridEntity =
        typename GridImp::Trimmer::TrimmerTraits::template Codim<0>::ParameterSpaceGridEntity;

  } // namespace IntersectionTraits

  template <class GridImp, IntersectionType type_>
  struct TimmedIntersectionImpl
  {
    constexpr static int dim      = GridImp::dimension;
    constexpr static int mydim    = GridImp::dimension - 1;
    constexpr static int dimworld = GridImp::dimension;
    using ctype                   = typename GridImp::ctype;
    using IdType                  = typename GridImp::GridFamily::TrimmerTraits::GlobalIdSetId;
    using NormalVector            = FieldVector<ctype, dim>;
    using LocalCoordinate         = FieldVector<ctype, mydim>;
    using EdgeInfo                = typename GridImp::Trimmer::TrimmerTraits::template Codim<1>::EntityInfo;

    using IntersectionGeometry     = GeometryKernel::NURBSPatch<mydim, dim, ctype>;
    using HostLeafIntersection     = IntersectionTraits::HostIntersection<GridImp, type_>;
    using ParameterSpaceGridEntity = IntersectionTraits::ParameterSpaceGridEntity<GridImp>;
    using LocalGeometry            = IntersectionTraits::LocalGeometry<GridImp>;
    using Geometry                 = IntersectionTraits::Geometry<GridImp>;

    TimmedIntersectionImpl() = default;

    TimmedIntersectionImpl(const GridImp* patchGrid, IdType& insideElementId, const EdgeInfo& edgeinfo)
        : patchGrid_(patchGrid),
          edgeInfo_(edgeinfo),
          insideElementId_(insideElementId),
          geo_(geometryForIdx(
              patchGrid_->trimmer().entityContainer_.idToElementInfoMap[insideElementId_].indexInLvlStorage)) {
      assert(geo_.domain()[0].isUnitDomain());
    }

    ParameterSpaceGridEntity inside() const {
      return patchGrid_->trimmer().entityContainer_.template entity<0>(insideElementId_, level());
    }

    ParameterSpaceGridEntity outside() const {
      DUNE_THROW(GridError, "No Outside Entities for trimmed Intersections");
    }

    int level() const {
      if constexpr (type_ == IntersectionType::Leaf)
        return patchGrid_->maxLevel();
      else
        return edgeInfo_.lvl;
    }

    bool boundary() const {
      return true;
    }

    NormalVector centerUnitOuterNormal() const {}

    bool neighbor() const {
      return false;
    }

    size_t boundarySegmentIndex() const {
      return 0;
    }

    bool conforming() const {
      return false;
    }

    GeometryType type() const {
      return GeometryTypes::none(mydim);
    }

    LocalGeometry geometryInInside() const {
      DUNE_THROW(NotImplemented, "");
    }

    LocalGeometry geometryInOutside() const {
      DUNE_THROW(GridError, "No Outside Entities for trimmed Intersections");
    }

    Geometry geometry() const {
      return geo_;
    }

    int indexInInside() const {
      DUNE_THROW(NotImplemented, "");
    }

    int indexInOutside() const {
      DUNE_THROW(NotImplemented, "");

    }

    NormalVector outerNormal(const LocalCoordinate& local) const {
      DUNE_THROW(NotImplemented, "");
    }

    NormalVector integrationOuterNormal(const LocalCoordinate& local) const {
      DUNE_THROW(NotImplemented, "");
    }

    NormalVector unitOuterNormal(const LocalCoordinate& local) const {
      DUNE_THROW(NotImplemented, "");
    }

  private:
    const GridImp* patchGrid_{};
    EdgeInfo edgeInfo_;
    IdType insideElementId_;
    IntersectionGeometry geo_;
  };

  template <class GridImp, IntersectionType type_>
  struct TimmedHostIntersectionImpl
  {
    constexpr static int dim      = GridImp::dimension;
    constexpr static int mydim    = GridImp::dimension - 1;
    constexpr static int dimworld = GridImp::dimension;
    using ctype                   = typename GridImp::ctype;
    using IdType                  = typename GridImp::GridFamily::TrimmerTraits::GlobalIdSetId;
    using NormalVector            = FieldVector<ctype, dim>;
    using LocalCoordinate         = FieldVector<ctype, mydim>;
    using EdgeInfo                = typename GridImp::Trimmer::TrimmerTraits::template Codim<1>::EntityInfo;

    using IntersectionGeometry     = GeometryKernel::NURBSPatch<mydim, dim, ctype>;
    using HostLeafIntersection     = IntersectionTraits::HostIntersection<GridImp, type_>;
    using ParameterSpaceGridEntity = IntersectionTraits::ParameterSpaceGridEntity<GridImp>;
    using LocalGeometry            = IntersectionTraits::LocalGeometry<GridImp>;
    using Geometry                 = IntersectionTraits::Geometry<GridImp>;

    int level() const {
      if constexpr (type_ == IntersectionType::Leaf)
        return patchGrid_->maxLevel();
      else
        return hostIntersection_.inside().level();
    }
    IdType insideElementId() const {
      auto hostId = patchGrid_->trimmer().parameterSpaceGrid().globalIdSet().id(hostIntersection_.inside());
      return {.entityIdType = IdType::EntityIdType::host, .id = hostId};
    }

    TimmedHostIntersectionImpl() = default;

    TimmedHostIntersectionImpl(const GridImp* patchGrid, const HostLeafIntersection& hostIntersection,
                               const EdgeInfo& edgeinfo)
        : patchGrid_(patchGrid),
          hostIntersection_(hostIntersection),
          edgeInfo_(edgeinfo),
          geo_(geometryForIdx(
              patchGrid_->trimmer().entityContainer_.idToElementInfoMap[insideElementId()].indexInLvlStorage)) {
      assert(geo_.domain()[0].isUnitDomain());
    }

    ParameterSpaceGridEntity inside() const {
      return patchGrid_->trimmer().entityContainer_.template entity<0>(insideElementId(), level());
    }

    ParameterSpaceGridEntity outside() const {
      auto hostId      = patchGrid_->trimmer().parameterSpaceGrid().globalIdSet().id(hostIntersection_.outside());
      IdType elementId = {.entityIdType = IdType::EntityIdType::host, .id = hostId};
      return patchGrid_->trimmer().entityContainer_.template entity<0>(elementId, level());
    }

    bool boundary() const {
      return hostIntersection_.boundary();
    }

    NormalVector centerUnitOuterNormal() const {
      return hostIntersection_.centerUnitOuterNormal();
    }

    bool neighbor() const {
      return hostIntersection_.neighbor();
    }

    // return the boundary segment index
    size_t boundarySegmentIndex() const {
      return 0;
    }

    // Return true if this is a conforming intersection
    bool conforming() const {
      return hostIntersection_.conforming();
    }

    // This is up to debate, if its not GeometryTypes::none(mydim);
    GeometryType type() const {
      return hostIntersection_.type();
    }

    // todo
    LocalGeometry geometryInInside() const {
      return hostIntersection_.geometryInInside();
    }

    // todo
    LocalGeometry geometryInOutside() const {
      return hostIntersection_.geometryInOutside();
    }

    // todo
    Geometry geometry() const {
      return hostIntersection_.geometry();
    }

    int indexInInside() const {
      return hostIntersection_.indexInInside();
    }

    int indexInOutside() const {
      return hostIntersection_.indexInOutside();
    }

    // todo
    NormalVector outerNormal(const LocalCoordinate& local) const {
      return hostIntersection_.outerNormal(local);
    }

    // todo
    NormalVector integrationOuterNormal(const LocalCoordinate& local) const {
      return hostIntersection_.integrationOuterNormal(local);
    }

    // todo
    NormalVector unitOuterNormal(const LocalCoordinate& local) const {
      return hostIntersection_.integrationOuterNormal(local);
    }

  private:
    const GridImp* patchGrid_{};
    HostLeafIntersection hostIntersection_;
    EdgeInfo edgeInfo_;
    IntersectionGeometry geo_;
  };

  template <class GridImp, IntersectionType type_>
  struct HostIntersectionImpl
  {
    constexpr static int dim      = GridImp::dimension;
    constexpr static int mydim    = GridImp::dimension - 1;
    constexpr static int dimworld = GridImp::dimension;
    using ctype                   = typename GridImp::ctype;
    using IdType                  = typename GridImp::GridFamily::TrimmerTraits::GlobalIdSetId;
    using NormalVector            = FieldVector<ctype, dim>;
    using LocalCoordinate         = FieldVector<ctype, mydim>;

    using IntersectionGeometry     = GeometryKernel::NURBSPatch<mydim, dim, ctype>;
    using HostLeafIntersection     = IntersectionTraits::HostIntersection<GridImp, type_>;
    using ParameterSpaceGridEntity = IntersectionTraits::ParameterSpaceGridEntity<GridImp>;
    using LocalGeometry            = IntersectionTraits::LocalGeometry<GridImp>;
    using Geometry                 = IntersectionTraits::Geometry<GridImp>;

    HostIntersectionImpl() = default;

    HostIntersectionImpl(const GridImp* patchGrid, const HostLeafIntersection& hostIntersection)
        : patchGrid_(patchGrid),
          hostIntersection_(hostIntersection) {}

    int level() const {
      if constexpr (type_ == IntersectionType::Leaf)
        return patchGrid_->maxLevel();
      else
        return hostIntersection_.inside().level();
    }

    ParameterSpaceGridEntity inside() const {
      auto hostId      = patchGrid_->trimmer().parameterSpaceGrid().globalIdSet().id(hostIntersection_.inside());
      IdType elementId = {.entityIdType = IdType::EntityIdType::host, .id = hostId};
      return patchGrid_->trimmer().entityContainer_.template entity<0>(elementId, level());
    }

    ParameterSpaceGridEntity outside() const {
      auto hostId      = patchGrid_->trimmer().parameterSpaceGrid().globalIdSet().id(hostIntersection_.outside());
      IdType elementId = {.entityIdType = IdType::EntityIdType::host, .id = hostId};
      return patchGrid_->trimmer().entityContainer_.template entity<0>(elementId, level());
    }

    bool boundary() const {
      return hostIntersection_.boundary();
    }

    NormalVector centerUnitOuterNormal() const {
      return hostIntersection_.centerUnitOuterNormal();
    }

    bool neighbor() const {
      return hostIntersection_.neighbor();
    }

    // return the boundary segment index
    size_t boundarySegmentIndex() const {
      return 0;
      // This is not implmented in SubGrid
    }

    // Return true if this is a conforming intersection
    bool conforming() const {
      return hostIntersection_.conforming();
    }

    GeometryType type() const {
      return hostIntersection_.type();
    }

    LocalGeometry geometryInInside() const {
      return hostIntersection_.geometryInInside();
    }

    LocalGeometry geometryInOutside() const {
      return hostIntersection_.geometryInOutside();
    }

    Geometry geometry() const {
      return hostIntersection_.geometry();
    }

    int indexInInside() const {
      return hostIntersection_.indexInInside();
    }

    int indexInOutside() const {
      return hostIntersection_.indexInOutside();
    }

    NormalVector outerNormal(const LocalCoordinate& local) const {
      return hostIntersection_.outerNormal(local);
    }

    NormalVector integrationOuterNormal(const LocalCoordinate& local) const {
      return hostIntersection_.integrationOuterNormal(local);
    }

    NormalVector unitOuterNormal(const LocalCoordinate& local) const {
      return hostIntersection_.integrationOuterNormal(local);
    }

  private:
    const GridImp* patchGrid_{};
    IntersectionGeometry geo;
    HostLeafIntersection hostIntersection_;
  };

  template <class GridImp, IntersectionType type_>
  struct IntersectionVariant
  {
    template <class Implementation>
    explicit IntersectionVariant(const Implementation& impl)
        : impl_(impl) {}
    IntersectionVariant() = default;

    auto visit(auto&& lambda) const {
      return std::visit(lambda, impl_);
    }

    auto inside() const {
      return visit([](const auto& impl) { return impl.inside(); });
    }

    auto outside() const {
      return visit([](const auto& impl) { return impl.outside(); });
    }

    bool boundary() const {
      return visit([](const auto& impl) { return impl.boundary(); });
    }

    auto centerUnitOuterNormal() const {
      return visit([](const auto& impl) { return impl.centerUnitOuterNormal(); });
    }

    bool neighbor() const {
      return visit([](const auto& impl) { return impl.neighbor(); });
    }

    size_t boundarySegmentIndex() const {
      return visit([](const auto& impl) { return impl.boundarySegmentIndex(); });
    }

    bool conforming() const {
      return visit([](const auto& impl) { return impl.conforming(); });
    }

    auto type() const {
      return visit([](const auto& impl) { return impl.type(); });
    }

    auto geometryInInside() const {
      return visit([](const auto& impl) { return impl.geometryInInside(); });
    }

    auto geometryInOutside() const {
      return visit([](const auto& impl) { return impl.geometryInOutside(); });
    }

    auto geometry() const {
      return visit([](const auto& impl) { return impl.geometry(); });
    }

    int indexInInside() const {
      return visit([](const auto& impl) { return impl.indexInInside(); });
    }

    int indexInOutside() const {
      return visit([](const auto& impl) { return impl.indexInOutside(); });
    }

    auto outerNormal(const auto& local) const {
      return visit([&](const auto& impl) { return impl.outerNormal(local); });
    }

    auto integrationOuterNormal(const auto& local) const {
      return visit([&](const auto& impl) { return impl.integrationOuterNormal(local); });
    }

    auto unitOuterNormal(const auto& local) const {
      return visit([&](const auto& impl) { return impl.unitOuterNormal(local); });
    }

  private:
    std::variant<TimmedIntersectionImpl<GridImp, type_>, TimmedHostIntersectionImpl<GridImp, type_>,
                 HostIntersectionImpl<GridImp, type_>>
        impl_{};
  };

} // namespace Impl

/** @brief An intersection with a leaf neighbor element
 * \ingroup PatchGrid
 * Mesh entities of codimension 0 ("elements") allow to visit all neighbors, where
 * a neighbor is an entity of codimension 0 which has a common entity of codimension 1
 * These neighbors are accessed via a IntersectionIterator. This allows the implement
 * non-matching meshes. The number of neighbors may be different from the number
 * of an element!
 */
template <class GridImp>
class TrimmedLeafIntersection
{
  friend typename GridImp::Traits::LeafIntersectionIterator;

  friend struct HostGridAccess<std::remove_const_t<GridImp>>;

  constexpr static int dim      = GridImp::dimension;
  constexpr static int mydim    = GridImp::dimension - 1;
  constexpr static int dimworld = GridImp::dimension;

  using Trimmer = typename GridImp::Trimmer;

  using HostLeafIntersection = typename GridImp::Trimmer::TrimmerTraits::HostLeafIntersection;
  using EdgeInfo             = typename GridImp::Trimmer::TrimmerTraits::template Codim<1>::EntityInfo;

public:
  // The type used to store coordinates
  typedef typename GridImp::ctype ctype;
  using LocalCoordinate = FieldVector<ctype, mydim>;

  typedef typename GridImp::Trimmer::TrimmerTraits::template Codim<1>::LocalParameterSpaceGeometry Geometry;
  typedef typename GridImp::Trimmer::TrimmerTraits::template Codim<1>::LocalGeometry LocalGeometry;
  typedef FieldVector<ctype, dim> NormalVector;

  using ParameterSpaceGridEntity =
      typename GridImp::Trimmer::TrimmerTraits::template Codim<0>::ParameterSpaceGridEntity;

  using IntersectionGeometry = GeometryKernel::NURBSPatch<mydim, dim, ctype>;
  using IntersectionImpl     = Impl::IntersectionVariant<GridImp, Impl::IntersectionType::Leaf>;

  TrimmedLeafIntersection() = default;

  TrimmedLeafIntersection(const GridImp* patchGrid, const HostLeafIntersection& hostIntersection)
      : underlying_{Impl::HostIntersectionImpl<GridImp, Impl::IntersectionType::Leaf>(patchGrid, hostIntersection)} {}

  // Trimmed Host Intersection
  TrimmedLeafIntersection(const GridImp* patchGrid, const HostLeafIntersection& hostIntersection,
                          const EdgeInfo& edgeInfo)
      : underlying_{Impl::TimmedHostIntersectionImpl<GridImp, Impl::IntersectionType::Leaf>(patchGrid, hostIntersection,
                                                                                            edgeInfo)} {}

  // TrimmedLeafIntersection(const GridImp* parameterSpaceGrid, HostLeafIntersection&& hostIntersection)
  //     : patchGrid_(parameterSpaceGrid),
  //       hostIntersection_{hostIntersection} {}

  bool operator==(const TrimmedLeafIntersection& other) const {
    DUNE_THROW(NotImplemented, "");
  }
  using IdType = typename GridImp::GridFamily::TrimmerTraits::GlobalIdSetId;

  // returns the inside entity
  ParameterSpaceGridEntity inside() const {
    return underlying_.inside();
  }

  // return Entity on the outside of this intersection
  // (that is the neighboring Entity)
  ParameterSpaceGridEntity outside() const {
    return underlying_.outside();
  }

  // return true if intersection is with boundary.
  [[nodiscard]] bool boundary() const {
    return underlying_.boundary();
  }

  /** @brief Return unit outer normal (length == 1)
   *
   *   The returned vector is the normal at the center() of the
   *     intersection's geometry.
   *       It is scaled to have unit length. */
  NormalVector centerUnitOuterNormal() const {
    return underlying_.centerUnitOuterNormal();
  }

  // return true if across the edge an neighbor on this level exists
  bool neighbor() const {
    return underlying_.neighbor();
  }

  // return the boundary segment index
  size_t boundarySegmentIndex() const {
    return underlying_.boundarySegmentIndex();
  }

  // Return true if this is a conforming intersection
  bool conforming() const {
    return underlying_.conforming();
  }

  // Geometry type of an intersection
  GeometryType type() const {
    return underlying_.type();
  }

  LocalGeometry geometryInInside() const {
    return underlying_.geometryInInside();
  }

  LocalGeometry geometryInOutside() const {
    return underlying_.geometryInOutside();
  }

  Geometry geometry() const {
    return underlying_.geometry();
  }

  // local number of codim 1 entity in self where intersection is contained in
  int indexInInside() const {
    return underlying_.indexInInside();
  }

  // local number of codim 1 entity in neighbor where intersection is contained
  int indexInOutside() const {
    return underlying_.indexInOutside();
  }

  // return outer normal
  FieldVector<ctype, dim> outerNormal(const LocalCoordinate& local) const {
    return underlying_.outerNormal(local);
  }

  // return outer normal multiplied by the integration element
  FieldVector<ctype, dim> integrationOuterNormal(const LocalCoordinate& local) const {
    return underlying_.integrationOuterNormal(local);
  }

  // return unit outer normal
  FieldVector<ctype, dim> unitOuterNormal(const LocalCoordinate& local) const {
    return underlying_.unitOuterNormal(local);
  }

private:
  IntersectionImpl underlying_{};
};

template <class GridImp>
class TrimmedLevelIntersection
{
  friend typename GridImp::Traits::LevelIntersectionIterator;

  friend struct HostGridAccess<std::remove_const_t<GridImp>>;

  constexpr static int dim   = GridImp::dimension;
  constexpr static int mydim = GridImp::dimension - 1;

  constexpr static int dimworld = dim;

  using Trimmer = typename GridImp::Trimmer;

  using ParameterSpaceGridEntity =
      typename GridImp::Trimmer::TrimmerTraits::template Codim<0>::ParameterSpaceGridEntity;
  typedef typename GridImp::Trimmer::TrimmerTraits::template Codim<1>::LocalParameterSpaceGeometry Geometry;
  typedef typename GridImp::Trimmer::TrimmerTraits::template Codim<1>::LocalGeometry LocalGeometry;

  using HostLevelIntersection = typename GridImp::Trimmer::TrimmerTraits::HostLevelIntersection;

  using MatrixHelper = MultiLinearGeometryTraits<double>::MatrixHelper;

public:
  using IntersectionImpl = Impl::IntersectionVariant<GridImp, Impl::IntersectionType::Level>;

  typedef typename GridImp::ctype ctype;
  using LocalCoordinate = FieldVector<ctype, mydim>;

  typedef FieldVector<ctype, dimworld> NormalVector;

  TrimmedLevelIntersection() = default;

  TrimmedLevelIntersection(const GridImp* parameterSpaceGrid, const HostLevelIntersection& hostIntersection)
      : underlying_{
            Impl::HostIntersectionImpl<GridImp, Impl::IntersectionType::Level>(parameterSpaceGrid, hostIntersection)} {}

  // TrimmedLevelIntersection(const GridImp* patchGrid, HostLevelIntersection&& hostIntersection)
  //     : patchGrid_(patchGrid),
  //       hostIntersection_{hostIntersection} {}

  bool operator==(const TrimmedLevelIntersection& other) const {
    DUNE_THROW(NotImplemented, "");
  }
  using IdType = typename GridImp::GridFamily::TrimmerTraits::GlobalIdSetId;

  // returns the inside entity
  ParameterSpaceGridEntity inside() const {
    return underlying_.inside();
  }

  // return Entity on the outside of this intersection
  // (that is the neighboring Entity)
  ParameterSpaceGridEntity outside() const {
    return underlying_.outside();
  }

  // return true if intersection is with boundary.
  [[nodiscard]] bool boundary() const {
    return underlying_.boundary();
  }

  /** @brief Return unit outer normal (length == 1)
   *
   *   The returned vector is the normal at the center() of the
   *     intersection's geometry.
   *       It is scaled to have unit length. */
  NormalVector centerUnitOuterNormal() const {
    return underlying_.centerUnitOuterNormal();
  }

  // return true if across the edge an neighbor on this level exists
  bool neighbor() const {
    return underlying_.neighbor();
  }

  // return the boundary segment index
  size_t boundarySegmentIndex() const {
    return underlying_.boundarySegmentIndex();
  }

  // Return true if this is a conforming intersection
  bool conforming() const {
    return underlying_.conforming();
  }

  // Geometry type of an intersection
  GeometryType type() const {
    return underlying_.type();
  }

  LocalGeometry geometryInInside() const {
    return underlying_.geometryInInside();
  }

  LocalGeometry geometryInOutside() const {
    return underlying_.geometryInOutside();
  }

  Geometry geometry() const {
    return underlying_.geometry();
  }

  // local number of codim 1 entity in self where intersection is contained in
  int indexInInside() const {
    return underlying_.indexInInside();
  }

  // local number of codim 1 entity in neighbor where intersection is contained
  int indexInOutside() const {
    return underlying_.indexInOutside();
  }

  // return outer normal
  FieldVector<ctype, dim> outerNormal(const LocalCoordinate& local) const {
    return underlying_.outerNormal(local);
  }

  // return outer normal multiplied by the integration element
  FieldVector<ctype, dim> integrationOuterNormal(const LocalCoordinate& local) const {
    return underlying_.integrationOuterNormal(local);
  }

  // return unit outer normal
  FieldVector<ctype, dim> unitOuterNormal(const LocalCoordinate& local) const {
    return underlying_.unitOuterNormal(local);
  }

private:
  IntersectionImpl underlying_{};
};

} // namespace Dune::IGANEW::DefaultTrim

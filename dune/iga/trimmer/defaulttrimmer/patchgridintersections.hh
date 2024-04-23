// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#pragma once

/** \file
 * @brief The TrimmedPatchGridLeafIntersection and TrimmedLevelIntersection classes
 */

// @todo same as gridintersection, one variant that forwards everything to hostEntity and one trimmed one
// @todo Also make this for leaf and level the same


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
class TrimmedLeafIntersection
{
  friend typename GridImp::Traits::LeafIntersectionIterator;

  friend struct HostGridAccess<std::remove_const_t<GridImp>>;

  constexpr static int dim   = GridImp::dimension;
  constexpr static int mydim = GridImp::dimension - 1;

  constexpr static int dimworld = GridImp::dimension;
  using Trimmer                 = typename GridImp::Trimmer;

  using HostLeafIntersection = typename GridImp::Trimmer::TrimmerTraits::HostLeafIntersection;

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

  TrimmedLeafIntersection() = default;

  TrimmedLeafIntersection(const GridImp* parameterSpaceGrid, const HostLeafIntersection& hostIntersection)
      : patchGrid_(parameterSpaceGrid),
        hostIntersection_{hostIntersection} {
  }

  TrimmedLeafIntersection(const GridImp* parameterSpaceGrid, HostLeafIntersection&& hostIntersection)
      : patchGrid_(parameterSpaceGrid),
        hostIntersection_{hostIntersection} {
  }
  HostLeafIntersection hostIntersection_;
  bool operator==(const TrimmedLeafIntersection& other) const {
    // DUNE_THROW(NotImplemented, "equals not implemented");
    return hostIntersection_ == other.hostIntersection_;
  }
  using IdType = typename GridImp::GridFamily::TrimmerTraits::GlobalIdSetId;

  //  returns the inside entity
  ParameterSpaceGridEntity inside() const {
    // DUNE_THROW(NotImplemented, "inside not implemented");
    auto hostId      = patchGrid_->trimmer().parameterSpaceGrid().globalIdSet().id(hostIntersection_.inside());
    IdType elementId = {.entityIdType = IdType::EntityIdType::host, .id = hostId};
    return patchGrid_->trimmer().entityContainer_.template entity<0>(elementId, patchGrid_->maxLevel());
  }

  //  return Entity on the outside of this intersection
  //  (that is the neighboring Entity)
  ParameterSpaceGridEntity outside() const {
    // DUNE_THROW(NotImplemented, "outside not implemented");
    auto hostId      = patchGrid_->trimmer().parameterSpaceGrid().globalIdSet().id(hostIntersection_.outside());
    IdType elementId = {.entityIdType = IdType::EntityIdType::host, .id = hostId};
    return patchGrid_->trimmer().entityContainer_.template entity<0>(elementId, patchGrid_->maxLevel());
  }

  //  return true if intersection is with boundary.
  [[nodiscard]] bool boundary() const {
    return hostIntersection_.boundary();
  }

  /** @brief Return unit outer normal (length == 1)
   *
   *   The returned vector is the normal at the center() of the
   *     intersection's geometry.
   *       It is scaled to have unit length. */
  NormalVector centerUnitOuterNormal() const {
    DUNE_THROW(NotImplemented, "centerUnitOuterNormal not implemented");
    //  @todo compute jacobian create normal by cross-product, the cross-product has to take into account the normal
    //  points away from the inside element
  }

  //  return true if across the edge an neighbor on this level exists
  bool neighbor() const {
    return hostIntersection_.neighbor();
  }

  //  return the boundary segment index
  size_t boundarySegmentIndex() const {
    return 0;
    // This is not implmented in SubGrid
  }

  //  Return true if this is a conforming intersection
  bool conforming() const {
    return hostIntersection_.conforming();
  }

  //  Geometry type of an intersection
  GeometryType type() const {
    return hostIntersection_.type();
  }

  //  @todo this function should return how this inersection resides in the inside host element.
  //  Therefore, this function should provide this intersection and maps to the 0..1 space
  LocalGeometry geometryInInside() const {
    return hostIntersection_.geometryInInside();
  }

  //  Same as above
  LocalGeometry geometryInOutside() const {
    return hostIntersection_.geometryInOutside();
  }

  //  geometry of the intersection this geometry should map into the knotspan domain
  using TrimmedParameterSpaceGeometry =
      typename GridImp::Trimmer::TrimmerTraits::template Codim<1>::TrimmedParameterSpaceGeometry;
  using LocalParameterSpaceGeometry =
      typename GridImp::Trimmer::TrimmerTraits::template Codim<1>::LocalParameterSpaceGeometry;

  LocalParameterSpaceGeometry geometry() const {
    return hostIntersection_.geometry();
  }

  //  local number of codim 1 entity in self where intersection is contained in
  int indexInInside() const {
    return hostIntersection_.indexInInside();
  }

  //  local number of codim 1 entity in neighbor where intersection is contained
  int indexInOutside() const {
    return hostIntersection_.indexInOutside();
  }

  //  return outer normal
  FieldVector<ctype, dim> outerNormal(const LocalCoordinate& local) const {
    return hostIntersection_.outerNormal(local);
  }

  //  return outer normal multiplied by the integration element
  FieldVector<ctype, dim> integrationOuterNormal(const LocalCoordinate& local) const {
    return hostIntersection_.integrationOuterNormal(local);
  }

  //  return unit outer normal
  FieldVector<ctype, dim> unitOuterNormal(const LocalCoordinate& local) const {
    return hostIntersection_.unitOuterNormal(local);
  }

private:
  const GridImp* patchGrid_{nullptr};
  IntersectionGeometry geo;
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
  typedef typename GridImp::ctype ctype;
  using LocalCoordinate = FieldVector<ctype, mydim>;

  typedef FieldVector<ctype, dimworld> NormalVector;

  TrimmedLevelIntersection() = default;

  TrimmedLevelIntersection(const GridImp* patchGrid, const HostLevelIntersection& hostIntersection)
      : patchGrid_(patchGrid),
        hostIntersection_{hostIntersection} {
  }

  TrimmedLevelIntersection(const GridImp* patchGrid, HostLevelIntersection&& hostIntersection)
      : patchGrid_(patchGrid),
        hostIntersection_{hostIntersection} {
  }

  [[nodiscard]] bool operator==(const TrimmedLevelIntersection& other) const {
    return hostIntersection_ == other.hostIntersection_;
  }

  using IdType = typename GridImp::GridFamily::TrimmerTraits::GlobalIdSetId;
  //  return Entity on the inside of this intersection
  //  (that is the Entity where we started this Iterator)
  [[nodiscard]] ParameterSpaceGridEntity inside() const {
    // DUNE_THROW(NotImplemented, "inside not implemented");
    auto hostId      = patchGrid_->trimmer().parameterSpaceGrid().globalIdSet().id(hostIntersection_.inside());
    IdType elementId = {.entityIdType = IdType::EntityIdType::host, .id = hostId};
    //  @todo Don't contruct this on the fly?
    return patchGrid_->trimmer().entityContainer_.template entity<0>(elementId, hostIntersection_.inside().level());
  }

  //  return Entity on the outside of this intersection
  //  (that is the neighboring Entity)
  [[nodiscard]] ParameterSpaceGridEntity outside() const {
    // DUNE_THROW(NotImplemented, "outside not implemented");
    auto hostId      = patchGrid_->trimmer().parameterSpaceGrid().globalIdSet().id(hostIntersection_.outside());
    IdType elementId = {.entityIdType = IdType::EntityIdType::host, .id = hostId};
    //  @todo Don't contruct this on the fly?trim
    return patchGrid_->trimmer().entityContainer_.template entity<0>(elementId, hostIntersection_.outside().level());
  }

  /** @brief return true if intersection is with boundary.
   */
  [[nodiscard]] bool boundary() const {
    return hostIntersection_.boundary();

    DUNE_THROW(NotImplemented, "boundary not implemented");

    return {};
  }

  /** @brief Return unit outer normal (length == 1)
   *
   *   The returned vector is the normal at the center() of the
   *     intersection's geometry.
   *       It is scaled to have unit length. */
  NormalVector centerUnitOuterNormal() const {
    return hostIntersection_.centerUnitOuterNormal();
  }

  //  return true if across the edge an neighbor on this level exists
  [[nodiscard]] bool neighbor() const {
    return hostIntersection_.neighbor();
  }

  //  return the boundary segment index
  [[nodiscard]] size_t boundarySegmentIndex() const {
    return 0;
    // This is not implmented in SubGrid
    return hostIntersection_.boundarySegmentIndex();
  }

  //  Return true if this is a conforming intersection, within one patch we are always conforming
  [[nodiscard]] bool conforming() const {
    return true;
  }

  //  Geometry type of an intersection
  [[nodiscard]] GeometryType type() const {
    return GeometryTypes::line;
  }

  //  intersection of codimension 1 of this neighbor with element where
  //  iteration started.
  //  Here returned element is in LOCAL coordinates of the element
  //  where iteration started.
  [[nodiscard]] LocalGeometry geometryInInside() const {
    return hostIntersection_.geometryInInside();
  }

  //  intersection of codimension 1 of this neighbor with element where iteration started.
  //  Here returned element is in LOCAL coordinates of neighbor
  [[nodiscard]] LocalGeometry geometryInOutside() const {
    return hostIntersection_.geometryInOutside();
  }

  //  intersection of codimension 1 of this neighbor with element where iteration started.
  //  Here returned element is in GLOBAL coordinates of the element where iteration started.
  using TrimmedParameterSpaceGeometry =
      typename GridImp::Trimmer::TrimmerTraits::template Codim<1>::TrimmedParameterSpaceGeometry;
  using LocalParameterSpaceGeometry =
      typename GridImp::Trimmer::TrimmerTraits::template Codim<1>::LocalParameterSpaceGeometry;

  LocalParameterSpaceGeometry geometry() const {
    return hostIntersection_.geometry();
  }

  //  local number of codim 1 entity in self where intersection is contained in
  [[nodiscard]] int indexInInside() const {
    return hostIntersection_.indexInInside();
  }

  //  local number of codim 1 entity in neighbor where intersection is contained
  [[nodiscard]] int indexInOutside() const {
    return hostIntersection_.indexInOutside();
  }

  //  return outer normal
  [[nodiscard]] FieldVector<ctype, dimworld> outerNormal(const LocalCoordinate& local) const {
    return hostIntersection_.outerNormal(local);
  }

  //  return outer normal multiplied by the integration element
  [[nodiscard]] FieldVector<ctype, dimworld> integrationOuterNormal(const LocalCoordinate& local) const {
    return hostIntersection_.integrationOuterNormal(local);
  }

  //  return unit outer normal
  [[nodiscard]] FieldVector<ctype, dimworld> unitOuterNormal(const LocalCoordinate& local) const {
    return hostIntersection_.integrationOuterNormal(local);
  }

private:
  const GridImp* patchGrid_;
  HostLevelIntersection hostIntersection_;
};

} // namespace Dune::IGANEW::DefaultTrim

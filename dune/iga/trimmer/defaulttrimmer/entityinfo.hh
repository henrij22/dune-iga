// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/common/reservedvector.hh>

namespace Dune::IGANEW::DefaultTrim {

template <typename HostIdType>
struct IdType;

/**
 *
 * @tparam Traits The trimmer traits
 * @tparam codim
 */
template <typename Traits, int codim>
struct EntityInfoImpl
{
};

template <typename Traits>
struct EntityInfoImpl<Traits, 2>
{
  static constexpr int codimension = 2;
  using HostIdType                 = typename Traits::ParameterSpaceGrid::GlobalIdSet::IdType;
  using EntitySeedType = typename Traits::ParameterSpaceGrid::template Codim<codimension>::Entity::EntitySeed;
  using TrimInfo       = typename Traits::ElementTrimData::VertexInfo;

  unsigned int indexInLvlStorage{std::numeric_limits<unsigned int>::max()};
  int lvl{};
  bool trimmed{false};
  IdType<HostIdType> id;
  EntitySeedType hostSeed{};

  std::optional<TrimInfo> trimInfo{};

  bool isTrimmed() const {
    return trimmed;
  }

  bool isValid() const {
    return indexInLvlStorage != std::numeric_limits<unsigned int>::max();
  }
};

template <typename Traits>
struct EntityInfoImpl<Traits, 1>
{
  static constexpr int codimension = 1;
  using HostIdType                 = typename Traits::ParameterSpaceGrid::GlobalIdSet::IdType;
  using EntitySeedType = typename Traits::ParameterSpaceGrid::template Codim<codimension>::Entity::EntitySeed;
  using TrimmedEntityGeometry =
      typename Traits::template Codim<codimension>::TrimmedParameterSpaceGeometry::PatchGeometry;

  using TrimInfo = typename Traits::ElementTrimData::EdgeInfo;

  unsigned int indexInLvlStorage{std::numeric_limits<unsigned int>::max()};
  int lvl{};
  bool trimmed{false};
  IdType<HostIdType> id;
  EntitySeedType hostSeed{};

  struct GeometryMap
  {
    unsigned int indexOfInsideElementinLvl;
    TrimmedEntityGeometry geometry;
  };

  /**
   * \brief for a given inside index return the edge geometry in the inside element
   * @param idx lvl index of inside element
   */
  TrimmedEntityGeometry geometryForIdx(unsigned int idx) const {
    if (auto it = iteratorForGeometryForIdx(idx); it != edgeGeometries.end())
      return it->geometry;
    DUNE_THROW(IOError, "geometryForIdx couldn't find geometry for idx");
  }

  /**
   * \brief for a given inside index return the edge geometry in the outside element
   * @param idx lvl index of inside element
   */
  TrimmedEntityGeometry otherGeometryForIdx(unsigned int idx) const {
    assert(edgeGeometries.size() == 2);
    if (auto it = iteratorForGeometryForIdx(idx); it != edgeGeometries.end()) {
      auto itD = std::ranges::distance(edgeGeometries.begin(), it);
      return itD == 1 ? edgeGeometries[0].geometry : edgeGeometries[1].geometry;
    }
    DUNE_THROW(IOError, "otherGeometryForIdx couldn't find geometry for idx");
  }

  std::vector<GeometryMap> edgeGeometries{};
  std::optional<TrimInfo> trimInfo{};

  bool isTrimmed() const {
    return trimmed;
  }
  bool isTrimmedHost() const {
    return hostSeed.isValid() and isTrimmed();
  }
  bool isValid() const {
    return indexInLvlStorage != std::numeric_limits<unsigned int>::max();
  }

private:
  auto iteratorForGeometryForIdx(unsigned int idx) const {
    return std::ranges::find_if(edgeGeometries,
                                [&](const auto& geoMap) { return geoMap.indexOfInsideElementinLvl == idx; });
  }
};

template <typename Traits>
struct EntityInfoImpl<Traits, 0>
{
  static constexpr int codimension = 0;
  using HostIdType                 = typename Traits::ParameterSpaceGrid::GlobalIdSet::IdType;
  using EntitySeedType = typename Traits::ParameterSpaceGrid::template Codim<codimension>::Entity::EntitySeed;

  unsigned int indexInLvlStorage{std::numeric_limits<unsigned int>::max()};
  unsigned int unTrimmedIndexInLvl{std::numeric_limits<unsigned int>::max()};
  unsigned int trimmedIndexInLvl{std::numeric_limits<unsigned int>::max()};
  unsigned int hostIndexInLvl{std::numeric_limits<unsigned int>::max()};
  int lvl{};
  bool trimmed{false};
  IdType<HostIdType> id{};
  EntitySeedType hostSeed{};

  auto isTrimmed() const {
    return trimmed;
  }
  bool isValid() const {
    return indexInLvlStorage != std::numeric_limits<unsigned int>::max();
  }

  std::optional<IdType<HostIdType>> fatherId;
  ReservedVector<IdType<HostIdType>, 4> decendantIds;
};

} // namespace Dune::IGANEW::DefaultTrim
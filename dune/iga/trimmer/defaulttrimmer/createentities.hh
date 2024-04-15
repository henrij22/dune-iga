// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/iga/trimmer/defaulttrimmer/trimmingutils/cliputils.hh>

namespace Dune::IGANEW::DefaultTrim {

template <int dim, int dimworld, typename ScalarType>
auto TrimmerImpl<dim, dimworld, ScalarType>::makeElementID(const HostEntity<0>& ele) -> GlobalIdType {
  auto& globalIdSetParameterSpace = parameterSpaceGrid_->globalIdSet();
  auto hostId                     = globalIdSetParameterSpace.id(ele);
  return {.entityIdType = GlobalIdType::EntityIdType::host, .id = hostId};
}

template <int dim, int dimworld, typename ScalarType>
auto TrimmerImpl<dim, dimworld, ScalarType>::idForTrimmedHostEdge(typename TrimmerTraits::PersistentIndexType hostEdgeId, const typename ElementTrimData::EdgeInfo& trimmedEdge) -> GlobalIdType {
  auto& globalIdSetParameterSpace = parameterSpaceGrid_->globalIdSet();

  // We have to check weather the edgeId is already there

  GlobalIdType edgeId = {.entityIdType = GlobalIdType::EntityIdType::newId,
                         .id           = globalIdSet_->newFreeIndex(),
                         .hostId       = std::make_optional(hostEdgeId)};

  return edgeId;
}

template <int dim, int dimworld, typename ScalarType>
void TrimmerImpl<dim, dimworld, ScalarType>::createAndSaveElementInfo(const std::tuple<unsigned int, unsigned int, int>& indices,
                                                                      const HostEntity<0>& ele, bool trimmed) {
  auto& globalIdSetParameterSpace = parameterSpaceGrid_->globalIdSet();

  std::optional<GlobalIdType> fatherId{};
  if (ele.hasFather())
    fatherId = {.entityIdType = GlobalIdType::EntityIdType::host, .id = globalIdSetParameterSpace.id(ele.father())};

  GlobalIdType elementId = makeElementID(ele);

  const unsigned int unTrimmedElementIndex = std::get<0>(indices);
  const unsigned int trimmedElementIndex   = std::get<1>(indices);
  const int newLevel              = std::get<2>(indices);
  EntityInfo<0> elementInfo       = {
            .indexInLvlStorage   = trimmedElementIndex + unTrimmedElementIndex,
            .unTrimmedIndexInLvl = unTrimmedElementIndex,
            .lvl                 = newLevel,
            .stemFromTrim        = trimmed,
            .id                  = elementId,
            .hostSeed            = ele.seed(),
            .fatherId            = fatherId,
  };

  entityContainer_.idToElementInfoMap.insert({elementId, elementInfo});

  // If we have a father we have to add us as his son, this can be faster, we can store in decendantIds,
  // the indexInLvlStorage and lvl, which would provide faster access
  if (fatherId.has_value()) {
    entityContainer_.idToElementInfoMap.at(fatherId.value()).decendantIds.push_back(elementId);
    entityContainer_.template entity<0>(fatherId.value(), newLevel - 1).entityInfo_.decendantIds.push_back(elementId);
  }
  std::get<0>(entityContainer_.entityImps_.back()).emplace_back(grid_, ele, elementInfo);
}

template <int dim, int dimworld, typename ScalarType>
void TrimmerImpl<dim, dimworld, ScalarType>::collectElementEdges(int level, const HostEntity<0>& ele,
                                                                 const ElementTrimData& eleTrimData) {
  auto& indexSet = parameterSpaceGrid_->levelGridView(level).indexSet();
  GlobalIdType elementId = makeElementID(ele);

  auto& elementEdgeIndices        = entityContainer_.globalEdgesIdOfElementsMap_[elementId];
  auto& globalIdSetParameterSpace = parameterSpaceGrid_->globalIdSet();

  // Define some lambdas (can be moved to seperate functions later)
  auto addHostEdge = [&](const unsigned int localEdgeIndex) {
    auto hostEdgeId = globalIdSetParameterSpace.subId(ele, localEdgeIndex, 1);
    GlobalIdType edgeId   = {.entityIdType = GlobalIdType::EntityIdType::host, .id = hostEdgeId};
    elementEdgeIndices.emplace_back(edgeId);

    if (entityContainer_.idToEdgeInfoMap.contains(edgeId))
      return;

    auto edge = ele.template subEntity<1>(localEdgeIndex);
    EntityInfo<1> edgeInfo{
      .indexInLvlStorage = indexSet.index(edge), .lvl = level, .stemFromTrim = false, .id = edgeId, .hostSeed = edge.seed()};
    entityContainer_.idToEdgeInfoMap.insert({edgeId, edgeInfo});

    // store vertex ids of edges, for subIndex method of indexSet
    auto& edgeVertexIndices = entityContainer_.globalVertexIdOfEdgesMap_[edgeId];

    if (edgeVertexIndices.size() < 2) {
      const auto& cube = ReferenceElements<ctype, mydimension>::cube();
      for (auto vertexLocalIndexWRTElement : cube.subEntities(localEdgeIndex, 1, 2)) {
        auto hostVertexId = globalIdSetParameterSpace.subId(ele, vertexLocalIndexWRTElement, 2);
        GlobalIdType vertexId   = {.entityIdType = GlobalIdType::EntityIdType::host, .id = hostVertexId};
        edgeVertexIndices.emplace_back(vertexId);
      }
    }
  };

  auto addTrimmedHostEdge = [&](int localEdgeIndex, const typename ElementTrimData::EdgeInfo& edgeOfTrimmedElement) {
    // To get neighborhoud information we save the original hostId
    auto hostEdgeId = globalIdSetParameterSpace.subId(ele, localEdgeIndex, 1);
    GlobalIdType edgeId = idForTrimmedHostEdge(hostEdgeId, edgeOfTrimmedElement);

    // \todo We have to somehow check if this edge is already present in the idSet
    // Check for hostId, then for coordinates
    // We have to get a unique .id

    // GlobalIdType edgeId = {.entityIdType = GlobalIdType::EntityIdType::newId,
    //                        .id           = globalIdSet_->newFreeIndex(),
    //                        .hostId       = std::make_optional(hostEdgeId)};

  };

  if (eleTrimData.flag() == ElementTrimFlag::full) {
    for (auto localEdgeIndex : Dune::range(ele.subEntities(1))) {
      addHostEdge(localEdgeIndex);
    }
  } else /* trimmed */ {
    for (auto edgeOfTrimmedElement : eleTrimData.edges()) {
      if (edgeOfTrimmedElement.isHost and not edgeOfTrimmedElement.isTrimmed) {
        // This is a host Edge with full length (do some as above)
        int localEdgeIndex = Util::edgeIndexMapping[edgeOfTrimmedElement.idx];
        addHostEdge(localEdgeIndex);
      } else if (edgeOfTrimmedElement.isHost and edgeOfTrimmedElement.isTrimmed) {
        // This is a host Edge which is partially trimmed
        auto localEdgeIndex    = Util::edgeIndexMapping[edgeOfTrimmedElement.idx];
        addTrimmedHostEdge(localEdgeIndex, edgeOfTrimmedElement);
      }
    }
  }
}

template <int dim, int dimworld, typename ScalarType>
void TrimmerImpl<dim, dimworld, ScalarType>::collectElementVertices(int level, const HostEntity<0>& ele,
                                                                    const ElementTrimData& eleTrimData) {
  GlobalIdType elementId = makeElementID(ele);
  auto& indexSet = parameterSpaceGrid_->levelGridView(level).indexSet();

  auto& elementVertexIndices      = entityContainer_.globalVerticesIdOfElementsMap[elementId];
  auto& globalIdSetParameterSpace = parameterSpaceGrid_->globalIdSet();

  if (eleTrimData.flag() == ElementTrimFlag::full) {
    // Save eleVertices in local order
    for (auto localVertexIndex : Dune::range(ele.subEntities(2))) {
      auto hostVertexId = globalIdSetParameterSpace.subId(ele, localVertexIndex, 2);
      GlobalIdType vertexId   = {.entityIdType = GlobalIdType::EntityIdType::host, .id = hostVertexId};
      elementVertexIndices.emplace_back(vertexId);

      auto vertex = ele.template subEntity<2>(localVertexIndex);
      EntityInfo<2> vertexInfo{
        .indexInLvlStorage = indexSet.index(vertex), .lvl = level, .stemFromTrim = false, .id = vertexId, .hostSeed = vertex.seed()};
      entityContainer_.idToVertexInfoMap.back().insert({vertexId, vertexInfo});
    }

  } else /* trimmed */ {

  }
}

template <int dim, int dimworld, typename ScalarType>
void TrimmerImpl<dim, dimworld, ScalarType>::createSubEntities(int level) {
  using EdgeEntity = typename TrimmerTraits::template Codim<1>::ParameterSpaceGridEntity;
  using VertexEntity = typename TrimmerTraits::template Codim<2>::ParameterSpaceGridEntity;

  auto& globalIdSetParameterSpace = parameterSpaceGrid_->globalIdSet();
  auto& indexSet = parameterSpaceGrid_->levelGridView(level).indexSet();
  auto& vertexMap = entityContainer_.idToVertexInfoMap.back();
  auto& edgeMap = entityContainer_.idToEdgeInfoMap;

  // we are resizing the containers, so we can add the entities in the indexSet order, the index is obtained from
  // the infos (also this is a bit more efficient as pushing back all the time)

  auto& vertexContainer  = std::get<2>(entityContainer_.entityImps_.back());
  vertexContainer.resize(vertexMap.size());
  for (auto& [vertexId, vertexInfo] : vertexMap) {
    auto vertex = parameterSpaceGrid_->entity(vertexInfo.hostSeed);
    vertexContainer[vertexInfo.indexInLvlStorage] = VertexEntity{grid_, vertex, vertexInfo};
  }

  auto& edgeContainer    = std::get<1>(entityContainer_.entityImps_.back());
  edgeContainer.resize(entityContainer_.sizeOfInfos(1, level));
  for (auto& [edgeId, edgeInfo] : edgeMap) {
    if (edgeInfo.lvl != level)
      continue;
    auto edge = parameterSpaceGrid_->entity(edgeInfo.hostSeed);
    edgeContainer[edgeInfo.indexInLvlStorage] = EdgeEntity{grid_, edge, edgeInfo};
  }
}
} // namespace Dune::IGANEW::DefaultTrim
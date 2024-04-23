// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#pragma once

#include <dune/grid/common/intersection.hh>
#include <dune/iga/hierarchicpatch/patchgridentity.hh>
#include <dune/iga/hierarchicpatch/patchgridintersections.hh>

/** \file
 * @brief The PatchGridLeafIntersectionIterator and PatchGridLevelIntersectionIterator classes
 */

namespace Dune::IGANEW::DefaultTrim {

namespace Impl {

  enum class PositionToken
  {
    Begin,
    End
  };

  enum IntersectionIteratorType
  {
    Level,
    Leaf
  };

  namespace IntersectionTraits {
    template <class GridImp, IntersectionIteratorType type_>
    using Intersection = std::conditional_t<type_ == Leaf, typename GridImp::Traits::LeafIntersection,
                                            typename GridImp::Traits::LevelIntersection>;

    template <class GridImp, IntersectionIteratorType type_>
    using ParameterSpaceIntersection =
        std::conditional_t<type_ == Leaf, typename GridImp::Trimmer::TrimmerTraits::ParameterSpaceLeafIntersection,
                           typename GridImp::Trimmer::TrimmerTraits::ParameterSpaceLevelIntersection>;

    template <class GridImp, IntersectionIteratorType type_>
    using ParameterSpaceIntersectionIterator =
        std::conditional_t<type_ == Leaf, typename GridImp::ParameterSpaceGrid::LeafGridView::IntersectionIterator,
                           typename GridImp::ParameterSpaceGrid::LevelGridView::IntersectionIterator>;

    template <class GridImp>
    using TrimInfo = typename GridImp::Trimmer::ElementTrimData;

  } // namespace IntersectionTraits

  template <class GridImp, IntersectionIteratorType type_>
  struct TrimmedIntersectionIterator
  {
    using Intersection               = IntersectionTraits::Intersection<GridImp, type_>;
    using ParameterSpaceIntersection = IntersectionTraits::ParameterSpaceIntersection<GridImp, type_>;
    using TrimInfo                   = IntersectionTraits::TrimInfo<GridImp>;

    TrimmedIntersectionIterator() = default;

    TrimmedIntersectionIterator(const GridImp* parameterSpaceGrid, const TrimInfo& trimInfo, PositionToken pos)
        : parameterSpaceGrid_(parameterSpaceGrid),
          trimData_(trimInfo),
          maxIterator_(trimInfo.size(1)),
          iterator_(pos == PositionToken::Begin ? 0 : maxIterator_) {}

    bool equals(const TrimmedIntersectionIterator& other) const {
      return trimData_ == other.trimData_ and iterator_ == other.iterator_;
    }

    void increment() {
      ++iterator_;
    }
    Intersection dereference() const {
      return Intersection{};
    }

  private:
    const GridImp* parameterSpaceGrid_{};
    std::optional<TrimInfo> trimData_{};

    unsigned int maxIterator_{};
    unsigned int iterator_{};
  };

  template <class GridImp, IntersectionIteratorType type_>
  struct HostIntersectionIterator
  {
    using Intersection                       = IntersectionTraits::Intersection<GridImp, type_>;
    using ParameterSpaceIntersection         = IntersectionTraits::ParameterSpaceIntersection<GridImp, type_>;
    using ParameterSpaceIntersectionIterator = IntersectionTraits::ParameterSpaceIntersectionIterator<GridImp, type_>;

    HostIntersectionIterator() = default;

    HostIntersectionIterator(const GridImp* parameterSpaceGrid, const ParameterSpaceIntersectionIterator& hostIterator)
        : parameterSpaceGrid_(parameterSpaceGrid),
          hostIterator_(hostIterator) {}

    bool equals(const HostIntersectionIterator& other) const {
      return hostIterator_ == other.hostIterator_;
    }

    void increment() {
      ++hostIterator_;
    }
    Intersection dereference() const {
      auto parameterspaceIntersection = ParameterSpaceIntersection(parameterSpaceGrid_, *hostIterator_);
      auto realIntersection = typename Intersection::Implementation(parameterSpaceGrid_, parameterspaceIntersection);
      return Intersection(realIntersection);
    }

  private:
    const GridImp* parameterSpaceGrid_{};
    ParameterSpaceIntersectionIterator hostIterator_{};
  };

  template <class GridImp, IntersectionIteratorType type_>
  struct IntersectionIteratorVariant
  {
    auto visit(auto&& lambda) const {
      return std::visit(lambda, impl_);
    }
    auto visit(auto&& lambda) {
      return std::visit(lambda, impl_);
    }

    template <class Implementation>
    explicit IntersectionIteratorVariant(const Implementation& impl)
        : impl_(impl) {}
    IntersectionIteratorVariant() = default;

    void increment() {
      visit([](auto& impl) { impl.increment(); });
    }
    auto dereference() const {
      return visit([](const auto& impl) { return impl.dereference(); });
    }

    bool equals(const IntersectionIteratorVariant& other) const {
      return visit([&]<typename T>(const T& impl) { return impl.equals(std::get<T>(other.impl_)); });
    }

  private:
    std::variant<TrimmedIntersectionIterator<GridImp, type_>, HostIntersectionIterator<GridImp, type_>> impl_{};
  };
} // namespace Impl

/** @brief Iterator over all element neighbors
 * @ingroup PatchGrid
 * Mesh entities of codimension 0 ("elements") allow to visit all neighbors, where
 * a neighbor is an entity of codimension 0 which has a common entity of codimension 1
 * These neighbors are accessed via a IntersectionIterator. This allows the implement
 * non-matching meshes. The number of neighbors may be different from the number
 * of an element!
 */
template <class GridImp>
class PatchGridLeafIntersectionIterator
{
  constexpr static int dim      = GridImp::dimension;
  constexpr static int dimworld = GridImp::dimensionworld;
  typedef typename GridImp::ctype ctype;

  using ParameterSpaceIntersectionIterator = typename GridImp::ParameterSpaceGrid::LeafGridView::IntersectionIterator;
  using ParameterSpaceLeafIntersection     = typename GridImp::Trimmer::TrimmerTraits::ParameterSpaceLeafIntersection;
  using LeafIntersection                   = typename GridImp::Traits::LeafIntersection;
  using TrimInfo                           = typename GridImp::Trimmer::ElementTrimData;

public:
  using Intersection  = Dune::Intersection<const GridImp, PatchGridLeafIntersection<GridImp>>;
  using IteratorImpl  = Impl::IntersectionIteratorVariant<GridImp, Impl::Leaf>;
  using PositionToken = Impl::PositionToken;

  PatchGridLeafIntersectionIterator() = default;

  PatchGridLeafIntersectionIterator(const GridImp* parameterSpaceGrid,
                                    const ParameterSpaceIntersectionIterator& hostIterator)
      : underlyingIterator_{Impl::HostIntersectionIterator<GridImp, Impl::Leaf>(parameterSpaceGrid, hostIterator)} {}

  PatchGridLeafIntersectionIterator(const GridImp* parameterSpaceGrid, const TrimInfo& trimInfo, PositionToken position)
      : underlyingIterator_{
            Impl::TrimmedIntersectionIterator<GridImp, Impl::Leaf>(parameterSpaceGrid, trimInfo, position)} {}

  // equality
  bool equals(const PatchGridLeafIntersectionIterator& other) const {
    return underlyingIterator_.equals(other.underlyingIterator_);
  }

  // prefix increment
  void increment() {
    underlyingIterator_.increment();
  }

  // @brief dereferencing
  LeafIntersection dereference() const {
    return underlyingIterator_.dereference();
  }

private:
  IteratorImpl underlyingIterator_{};
};

template <class GridImp>
class PatchGridLevelIntersectionIterator
{
  constexpr static int dim      = GridImp::dimension;
  constexpr static int dimworld = GridImp::dimensionworld;
  typedef typename GridImp::ctype ctype;

  using ParameterSpaceIntersectionIterator = typename GridImp::ParameterSpaceGrid::LevelGridView::IntersectionIterator;
  using ParameterSpaceLevelIntersection    = typename GridImp::Trimmer::TrimmerTraits::ParameterSpaceLevelIntersection;
  using LevelIntersection                  = typename GridImp::Traits::LevelIntersection;
  using TrimInfo                           = typename GridImp::Trimmer::ElementTrimData;

public:
  using Intersection = Dune::Intersection<const GridImp, PatchGridLevelIntersection<GridImp>>;
  using IteratorImpl = Impl::IntersectionIteratorVariant<GridImp, Impl::Level>;

  using PositionToken = Impl::PositionToken;

  PatchGridLevelIntersectionIterator() = default;

  PatchGridLevelIntersectionIterator(const GridImp* parameterSpaceGrid,
                                     const ParameterSpaceIntersectionIterator& hostIterator)
      : underlyingIterator_{Impl::HostIntersectionIterator<GridImp, Impl::Level>(parameterSpaceGrid, hostIterator)} {}

  PatchGridLevelIntersectionIterator(const GridImp* parameterSpaceGrid, const TrimInfo& trimInfo,
                                     PositionToken position)
      : underlyingIterator_{
            Impl::TrimmedIntersectionIterator<GridImp, Impl::Level>(parameterSpaceGrid, trimInfo, position)} {}

  // equality
  bool equals(const PatchGridLevelIntersectionIterator& other) const {
    return underlyingIterator_.equals(other.underlyingIterator_);
  }

  // prefix increment
  void increment() {
    underlyingIterator_.increment();
  }

  // @brief dereferencing
  LevelIntersection dereference() const {
    return underlyingIterator_.dereference();
  }

private:
  IteratorImpl underlyingIterator_{};
};

} // namespace Dune::IGANEW::DefaultTrim

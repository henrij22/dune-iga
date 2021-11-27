//
// Created by lex on 23.11.21.
//

#pragma once

namespace Dune::IGA {
  namespace Impl {
    static constexpr int noNeighbor = -1;
  }

  template <typename NURBSIntersection>
  class NURBSGridInterSectionIterator;

  template <std::integral auto mydim, typename GridViewImp>
  class NURBSintersection {
  public:
    using Iterator      = NURBSGridInterSectionIterator<NURBSintersection>;
    using Entity        = typename GridViewImp::template Codim<0>::Entity;
    using Geometry      = typename GridViewImp::template Codim<1>::Geometry;
    using LocalGeometry = typename GridViewImp::template Codim<1>::LocalGeometry;

    NURBSintersection(int innerLocalIndex, int outerLocalIndex, unsigned int innerDirectIndex, unsigned int outerDirectIndex,
                      const GridViewImp& gridView)
        : gridView_{&gridView},
          innerDirectIndex_{innerDirectIndex},
          outerDirectIndex_{outerDirectIndex},
          innerLocalIndex_{innerLocalIndex},
          outerLocalIndex_{outerLocalIndex} {}
    /** coordinate type */
    typedef double ctype;

    /** \brief Dimension of the cube element */
    static constexpr std::integral auto mydimension = mydim;

    /** \brief Dimension of the world space that the cube element is embedded in*/
    static constexpr std::integral auto dimworld = GridViewImp::dimensionworld;

    /** \brief Type used for a vector of element coordinates */
    using LocalCoordinate = typename GridViewImp::NurbsGridLinearAlgebraTraits::template FixedVectorType<mydim>;

    /** \brief Type used for a vector of world coordinates */
    using GlobalCoordinate = typename GridViewImp::NurbsGridLinearAlgebraTraits::template FixedVectorType<dimworld>;
    [[nodiscard]] bool boundary() const { return outerDirectIndex_ != Impl::noNeighbor; }
    [[nodiscard]] bool neighbor() const { return outerDirectIndex_ != Impl::noNeighbor; }
    [[nodiscard]] bool conforming() const { return true; }
    [[nodiscard]] GeometryType type() const { return GeometryTypes::cube(mydim); }

    Entity inside() const { return gridView_->template getEntity<0>(innerDirectIndex_); }
    Entity outside() const {
      assert(neighbor() && "Outer Element does not exist.");
      return gridView_->template getEntity<0>(outerDirectIndex_);
    }
    [[nodiscard]] int indexInInside() const { return innerLocalIndex_; }
    [[nodiscard]] int indexInOutside() const { return outerLocalIndex_; }
    LocalGeometry geometryInInside() const { return LocalGeometry(innerLocalIndex_); }
    LocalGeometry geometryInOutside() const { return LocalGeometry(outerLocalIndex_); }
    Geometry geometry() const { return inside().template subEntity<1>(innerLocalIndex_).geometry(); }
    GlobalCoordinate outerNormal(const LocalCoordinate& xi) const {
      if constexpr (mydim == 0)
        return this->geometry().jacobianTransposed(xi)[0];
      else if constexpr (mydim == 1) {  // edges
        if constexpr (dimworld == 2) {  // edges in R2
          const auto innerJacobianTransposed = this->geometry().jacobianTransposed(xi)[0];
          switch (innerLocalIndex_) {
            case 0:
            case 3: return GlobalCoordinate({-innerJacobianTransposed[1], innerJacobianTransposed[0]});
            case 1:
            case 2: return GlobalCoordinate({innerJacobianTransposed[1], -innerJacobianTransposed[0]});
            default: __builtin_unreachable();
          }
        } else if constexpr (dimworld == 3)  // edges in R3
        {
          const auto xiInElementCoords        = geometryInInside().global(xi);
          const auto insideJacobianTransposed = inside().geometry().jacobianTransposed(xiInElementCoords);
          auto normal                         = cross(insideJacobianTransposed[0], insideJacobianTransposed[1]);
          switch (innerLocalIndex_) {
            case 0: return cross(normal, insideJacobianTransposed[1]);
            case 1: return cross(insideJacobianTransposed[1], normal);
            case 2: return cross(insideJacobianTransposed[0], normal);
            case 3: return cross(normal, insideJacobianTransposed[0]);
            default: __builtin_unreachable();
          }
        }
      } else if constexpr (mydim == 2 && dimworld == 3) {  // surfaces in R3
        const auto innerJacobianTransposed = this->geometry().jacobianTransposed(xi);
        switch (innerLocalIndex_) {
          case 0:
          case 3:
          case 4: return cross(innerJacobianTransposed[1], innerJacobianTransposed[0]);
          case 1:
          case 2:
          case 5: return cross(innerJacobianTransposed[0], innerJacobianTransposed[1]);
          default: __builtin_unreachable();
        }
      }
      __builtin_unreachable();
    }
    GlobalCoordinate unitOuterNormal(const LocalCoordinate& xi) const {
      auto N = this->outerNormal(xi);
      return N / N.two_norm();
    }
    GlobalCoordinate integrationOuterNormal(const LocalCoordinate& xi) const {
      return this->unitOuterNormal() * this->geometry().integrationElement(xi);
    }
    [[nodiscard]] std::size_t boundarySegmentIndex() const {
      throw std::logic_error("boundarySegmentIndex Not implemented");
      return 0;
    }

    auto operator<=>(const NURBSintersection&) const = default;

  private:
    const GridViewImp* gridView_;
    unsigned int innerDirectIndex_;
    int innerLocalIndex_;
    unsigned int outerDirectIndex_;
    int outerLocalIndex_;
  };
}  // namespace Dune::IGA
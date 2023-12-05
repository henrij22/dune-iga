// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file trimmer.hh
 * @brief Definition of the identity trimmer class.
 * @author Alexander Müller <mueller@ibb.uni-stuttgart.de>
 * @date 2023
 */

#pragma once

#include "patchgridindexsets.hh"
#include "patchgridleafiterator.hh"

#include <dune/geometry/referenceelements.hh>

#include <dune/grid/concepts.hh>
#include <dune/grid/yaspgrid.hh>

#include "dune/iga/hierarchicpatch/patchgridfwd.hh"
#include <dune/iga/hierarchicpatch/patchgridentityseed.hh>
#include <dune/iga/hierarchicpatch/patchgridgeometry.hh>
#include <dune/iga/hierarchicpatch/patchgridintersectioniterator.hh>
#include <dune/iga/hierarchicpatch/patchgridleveliterator.hh>
#include <dune/iga/hierarchicpatch/patchgridlocalgeometry.hh>
#include <dune/iga/hierarchicpatch/patchgridview.hh>

namespace Dune::IGANEW {

  namespace GeometryKernel {
    template <int dim_, int dimworld_, typename ScalarType>
    class NURBSPatch;
  }

  namespace IdentityTrim {

    /**
     * @brief Parameter struct representing parameters for the trimming operation.
     */
    struct Parameter {};

    /**
     * @brief ElementTrimData struct representing trim data for an element.
     * @tparam mydim_ Dimension of the element.
     * @tparam ScalarType Scalar type for the coordinates.
     */
    template <int mydim_, typename ScalarType>
    struct ElementTrimData {};

    /**
     * @brief ElementTrimDataContainer struct representing a container for element trim data.
     * @tparam ParameterSpaceGrid Type of the parameter space grid.
     */
    template <typename ParameterSpaceGrid>
    struct ElementTrimDataContainer {};

    /**
     * @brief PatchTrimData struct representing trim data for a patch.
     * @tparam dim Dimension of the patch.
     * @tparam ScalarType Scalar type for the coordinates.
     */
    template <int dim, typename ScalarType>
    struct PatchTrimData {};
    template <int dim, int dimworld, typename ScalarType>
    class Trimmer;
    template <int dim, int dimworld, typename ScalarType>
    struct PatchGridFamily {
      using ctype = ScalarType;
      using Grid  = PatchGrid<dim, dimworld, PatchGridFamily, ScalarType>;
      using Trimmer        = Trimmer<dim,dimworld,ScalarType>;

      using GlobalIdSetType = PatchGridGlobalIdSet<const Grid>;
      using LocalIdSetType  = PatchGridLocalIdSet<const Grid>;
      using LevelIndexSet   = PatchGridLevelIndexSet<const Grid>;
      using LeafIndexSet    = PatchGridLeafIndexSet<const Grid>;
      template <int codim, PartitionIteratorType pitype>
      using LeafIterator = PatchGridLeafIterator<codim, pitype, const Grid>;
      // using GlobalIdSetIdType = typename GlobalIdSetType<const Grid>;

      struct TrimmerTraits {
        using ParameterSpaceGrid
            = YaspGrid<dim, TensorProductCoordinates<ScalarType, dim>>;  ///< Type of the Parametric grid
      };
      // clang-format off
      typedef GridTraits<
        dim, dimworld, Grid,
      PatchGridGeometry,
      PatchGridEntity,
      PatchGridLevelIterator,
      PatchGridLeafIntersection,
      PatchGridLevelIntersection,
      PatchGridLeafIntersectionIterator,
      PatchGridLevelIntersectionIterator,
      PatchGridHierarchicIterator,
      PatchGridLeafIterator,
      LevelIndexSet,
      LeafIndexSet,
      GlobalIdSetType,
      typename TrimmerTraits::ParameterSpaceGrid::Traits::GlobalIdSet::IdType,
      LocalIdSetType,
      typename TrimmerTraits::ParameterSpaceGrid::Traits::LocalIdSet::IdType,
      Communication<No_Comm>,
      PatchGridLevelGridViewTraits,
      PatchGridLeafGridViewTraits,
      PatchGridEntitySeed,
      PatchGridLocalGeometry,
          typename TrimmerTraits::ParameterSpaceGrid::Traits::LevelIndexSet::IndexType,
          typename TrimmerTraits::ParameterSpaceGrid::Traits::LevelIndexSet::Types,
          typename TrimmerTraits::ParameterSpaceGrid::Traits::LeafIndexSet::IndexType,
          typename TrimmerTraits::ParameterSpaceGrid::Traits::LeafIndexSet::Types>
          Traits;
      // clang-format on
      // using GlobalIdSetType =  PatchGridGlobalIdSet<const Grid>;

      // template <int codim, PartitionIteratorType pitype>
      // using LeafIterator = PatchGridLeafIterator<codim,pitype,GridImp>;
    };

    /**
     * @brief Trimmer struct representing a trimmer with identity trimming (no trimming at all).
     * @ingroup Trimmer
     * @tparam dim Dimension of the patch.
     * @tparam ScalarType Scalar type for the coordinates.
     */
    template <int dim, int dimworld, typename ScalarType>
    class Trimmer {
     public:
      using GridFamily = PatchGridFamily<dim, dimworld, ScalarType>;  ///< Scalar type for the coordinates.
      using Traits = typename GridFamily::Traits;
      using TrimmerTraits = typename GridFamily::TrimmerTraits;

      using GridImp                       = typename Traits::Grid;
      static constexpr int mydimension    = GridImp::dimension;     ///< Dimension of the patch.
      static constexpr int dimensionworld = GridImp::dimensionworld;  ///< Dimension of the world.
      using ctype                         = ScalarType;

      template <int codim>
      static constexpr bool isLocalGeometryLinear
          = true;  ///< boolean for the linearity of the local geometry, for the untrimmed case this is always true
      static constexpr bool isAlwaysTrivial = true;  ///< Boolean indicating if the trimming is always trivial, no
                                                     ///< trimming or simple deletion of element.

      /**
       * @brief Default constructor for Trimmer.
       */
      Trimmer() = default;

      using ParameterSpaceGrid = typename TrimmerTraits::ParameterSpaceGrid;  ///< Type of the Parametric grid
      template <int mydim>
      using ReferenceElementType =
          typename Dune::Geo::ReferenceElements<ctype, mydimension>::ReferenceElement;  ///< Reference element type.

      using ElementTrimData = ElementTrimData<ParameterSpaceGrid::dimension,
                                              typename ParameterSpaceGrid::ctype>;  ///< Element trim data type.
      using PatchTrimData   = PatchTrimData<ParameterSpaceGrid::dimension,
                                          typename ParameterSpaceGrid::ctype>;  ///< Patch trim data type.
      using ElementTrimDataContainer
          = ElementTrimDataContainer<ParameterSpaceGrid>;  ///< Container for element trim data.

      template <int codim>
      using LocalParameterSpaceGeometry = typename ParameterSpaceGrid::template Codim<codim>::Geometry;

      template <int codim, PartitionIteratorType pitype>
      using ParameterSpaceLeafIterator =
          typename ParameterSpaceGrid::template Codim<codim>::template Partition<pitype>::LeafIterator;

      using LocalIdSetType = PatchGridLocalIdSet<GridImp>;
      using LeafIndexSet   = PatchGridLeafIndexSet<GridImp>;
      using LevelIndexSet  = PatchGridLevelIndexSet<GridImp>;

      using EntityContainerType = Empty;

      template <int codim>
      using LocalGeometry = typename ParameterSpaceGrid::template Codim<codim>::LocalGeometry;

      template <int codim>
      using ParameterSpaceGridEntity = typename ParameterSpaceGrid::template Codim<codim>::Entity;

      using ParameterType = Parameter;  ///< Type for trimming parameters.

      /**
       * @brief Get the reference element for a given entity.
       * @tparam EntityType Type of the entity.
       * @param entity The entity for which the reference element is requested.
       * @return Reference element for the entity.
       */
      template </* Dune::Concept::Entity */ typename EntityType>
      static auto referenceElement(const EntityType& entity) {
        return Dune::referenceElement<ctype, EntityType::mydimension>(entity.type());
      }

      /**
       * @brief Create the parameter space grid based on the patch and trim data.
       * @tparam dimworld Dimension of the world.
       * @param patchData NURBS patch data.
       * @param trimData Optional patch trim data.
       */
      void createParameterSpaceGrid(const GeometryKernel::NURBSPatch<mydimension, dimensionworld, ctype>& patch,
                                    const std::optional<PatchTrimData>& trimData) {
        parameterSpaceGrid_ = std::make_unique<ParameterSpaceGrid>(patch.uniqueKnotVector());
      }

      /**
       * @brief Constructor for Trimmer with patch and trim data.
       * @tparam dimworld Dimension of the world.
       * @param patch NURBS patch data.
       * @param trimData Optional patch trim data.
       */
      Trimmer(const GeometryKernel::NURBSPatch<mydimension, dimensionworld, ctype>& patch,
              const std::optional<PatchTrimData>& trimData) {
        createParameterSpaceGrid(patch, trimData);
      }

      /**
       * @brief Pass parameters to the trimmer.
       * @param par The parameters.
       */
      void setup(const ParameterType&) {}

      /**
       * @brief Get a const reference to the parameter space grid.
       * @return Const reference to the parameter space grid.
       */
      const ParameterSpaceGrid& parameterSpaceGrid() const { return *parameterSpaceGrid_; }

      /**
       * @brief Get a reference to the parameter space grid.
       * @return Reference to the parameter space grid.
       */
      ParameterSpaceGrid& parameterSpaceGrid() { return *parameterSpaceGrid_; }

      /**
       * @brief Refine the grid globally.
       * @param ref Number of refinement levels.
       */
      auto globalRefine(int ref) {
        // fill up container
        // patchTrimData,trimDatas_;

        // @todo Trim move the refine here from the grid
        ;
      }

     private:
      std::unique_ptr<ParameterSpaceGrid> parameterSpaceGrid_;  ///< The parameter space grid.
    };

  }  // namespace IdentityTrim
}  // namespace Dune::IGANEW

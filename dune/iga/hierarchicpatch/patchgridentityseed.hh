// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#pragma once

/**
 * \file
 * @brief The PatchGridEntitySeed class
 */

namespace Dune::IGANEW {

  /**
   * @brief The EntitySeed class provides the minimal information needed to restore an Entity using the grid.
   * \ingroup PatchGrid
   *
   */
  template <int codim, class GridImp>
  class PatchGridEntitySeed {
   protected:
    // Entity type of the hostgrid
    typedef typename GridImp::ParameterSpaceGrid::Traits::template Codim<codim>::Entity ParameterSpaceGridEntity;

    // EntitySeed type of the hostgrid
    typedef
        typename GridImp::ParameterSpaceGrid::Traits::template Codim<codim>::EntitySeed ParameterSpaceGridEntitySeed;

   public:
    constexpr static int codimension = codim;

    /**
     * @brief Construct an empty (i.e. isValid() == false) seed.
     */
    PatchGridEntitySeed() = default;

    /**
     * @brief Create EntitySeed from hostgrid Entity
     *
     * We call hostEntity.seed() directly in the constructor
     * of PatchGridEntitySeed to allow for return value optimization.
     */
    PatchGridEntitySeed(const ParameterSpaceGridEntity& hostEntity) : hostEntitySeed_(hostEntity.seed()) {}

    /**
     * @brief Get stored ParameterSpaceGridEntitySeed
     */
    const ParameterSpaceGridEntitySeed& hostEntitySeed() const { return hostEntitySeed_; }

    /**
     * @brief Check whether it is safe to create an Entity from this Seed
     */
    bool isValid() const { return hostEntitySeed_.isValid(); }

   private:
    ParameterSpaceGridEntitySeed hostEntitySeed_;
  };

}  // namespace Dune::IGANEW

// #define DUNE_IDENTITY_GRID_ENTITY_SEED_HH
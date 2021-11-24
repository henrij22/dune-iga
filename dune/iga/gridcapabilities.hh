//
// Created by lex on 16.11.21.
//

#pragma once

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/test/checkindexset.hh>
#include <dune/grid/test/gridcheck.hh>

namespace Dune::Capabilities {
  template <std::integral auto dim, std::integral auto dimworld, int codim>
  requires(codim == 0 || (codim == 1 && dim < 3) || (codim == 2 && dim == 3)
           || codim == dim) struct hasEntity<Dune::IGA::NURBSGrid<dim, dimworld>, codim> {
    static const bool v = true;
  };
}  // namespace Dune::Capabilities

template <std::integral auto dim, std::integral auto dimworld>
struct Dune::EnableBoundarySegmentIndexCheck<Dune::IGA::NURBSGrid<dim, dimworld>> : public std::false_type {};
template <std::integral auto dim, std::integral auto dimworld>
struct EnableLevelIntersectionIteratorCheck<Dune::IGA::NURBSGrid<dim, dimworld>> {
  static const bool v = true;
};

// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#define DUNE_CHECK_BOUNDS
#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <cfenv>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>

#include <dune/iga/hierarchicpatch/patchgridfactory.hh>
#include <dune/iga/patchgrid.hh>
#include <dune/iga/trimmer/defaulttrimmer/trimmer.hh>

using namespace Dune::IGANEW;

auto testIbraReader() {
  Dune::TestSuite t("", Dune::TestSuite::ThrowPolicy::ThrowOnRequired);

  using PatchGrid   = PatchGrid<2, 2, DefaultTrim::PatchGridFamily>;
  using GridFactory = Dune::GridFactory<PatchGrid>;

  auto gridFactory = GridFactory();
  gridFactory.insertTrimParameters(GridFactory::TrimParameterType{100});

  const std::vector testCases{std::tuple<std::string, int, int>{"auxiliaryfiles/element_trim_xb.ibra", 0, 3},
  {"auxiliaryfiles/element_trim.ibra", 0, 3}, {"auxiliaryfiles/trim_2edges.ibra", 0, 3}, {"auxiliaryfiles/trim_multi.ibra", 0, 3},
  {"auxiliaryfiles/surface-hole.ibra", 1, 3}, {"auxiliaryfiles/surface-hole-skew.ibra", 1, 3},
    {"auxiliaryfiles/surface-hole-square.ibra", 1, 3}};

  for (auto& [name, min, max] : testCases) {
    for (int i = min; i <= max; i++) {
      gridFactory.insertJson(name, true, {i, i});
      t.checkNoThrow([&]{auto grid = gridFactory.createGrid();});
    }
  }

  return t;
}

// auto testIbraReaderWithSquareHole() {
//   Dune::TestSuite t("", Dune::TestSuite::ThrowPolicy::ThrowOnRequired);
//
//   using PatchGrid   = PatchGrid<2, 2, DefaultTrim::PatchGridFamily>;
//   using GridFactory = Dune::GridFactory<PatchGrid>;
//
//   auto gridFactory = GridFactory();
//   gridFactory.insertTrimParameters(GridFactory::TrimParameterType{100});
//   gridFactory.insertJson("auxiliaryfiles/surface-hole-square.ibra", true, {2, 2});
//
//   t.checkNoThrow([&]{auto grid = gridFactory.createGrid();});
//
//   return t;
// }

int main(int argc, char** argv) try {
  feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

  // Initialize MPI, if necessary
  Dune::MPIHelper::instance(argc, argv);
  Dune::TestSuite t("", Dune::TestSuite::ThrowPolicy::ThrowOnRequired);

  t.subTest(testIbraReader());

  // auto success = t.report();

  return t.exit();
} catch (Dune::Exception& e) {
  std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
}

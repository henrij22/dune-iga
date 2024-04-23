// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later
#define DUNE_CHECK_BOUNDS
#ifdef HAVE_CONFIG_H
  #include "config.h"
#endif

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/iga/trimmer/defaulttrimmer/trimmingutils/indextransformations.hh>

using namespace Dune;

auto testTransformations() {
  TestSuite t;

  constexpr std::array vertexIndexMapping = {0u, 1u, 3u, 2u};
  constexpr std::array edgeIndexMapping   = {2u, 1u, 3u, 0u};

  using Transformations = Dune::IGANEW::DefaultTrim::Transformations;

  for (const auto i : Dune::range(4u)) {
    t.check(Transformations::mapToDune(2, i) == vertexIndexMapping[i]);
    t.check(Transformations::mapToDune(1, i) == edgeIndexMapping[i]);
  }

  constexpr std::array vertexIndexBackMapping = {0u, 1u, 2u, 3u};
  constexpr std::array edgeIndexBackMapping   = {3u, 1u, 0u, 2u};

  for (const auto i : Dune::range(4u)) {
    t.check(Transformations::mapToTrimmer(2, i) == vertexIndexBackMapping[i]);
    t.check(Transformations::mapToTrimmer(1, i) == edgeIndexBackMapping[i]);
  }

  return t;
}

#include <cfenv>
int main(int argc, char** argv) try {
  feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

  // Initialize MPI, if necessary
  Dune::MPIHelper::instance(argc, argv);
  Dune::TestSuite t("", Dune::TestSuite::ThrowPolicy::ThrowOnRequired);

  t.report();

  return t.exit();
} catch (Dune::Exception& e) {
  std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
}
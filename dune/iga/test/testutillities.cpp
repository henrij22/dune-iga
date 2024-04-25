// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later
#define DUNE_CHECK_BOUNDS
#ifdef HAVE_CONFIG_H
  #include "config.h"
#endif

#include <dune/common/exceptions.hh>
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
    t.check(Transformations::mapToDune(2, i) == vertexIndexMapping[i])
        << "Map to Dune for Vertex not correct, mapping says its " << Transformations::mapToDune(2, i)
        << " but should be " << vertexIndexMapping[i];
    t.check(Transformations::mapToDune(1, i) == edgeIndexMapping[i])
        << "Map to Dune for Edge not correct, mapping says its " << Transformations::mapToDune(1, i)
        << " but should be " << edgeIndexMapping[i];
  }

  constexpr std::array vertexIndexBackMapping = {0u, 1u, 3u, 2u};
  constexpr std::array edgeIndexBackMapping   = {3u, 1u, 0u, 2u};

  for (const auto i : Dune::range(4u)) {
    t.check(Transformations::mapToTrimmer(2, i) == vertexIndexBackMapping[i])
        << "Map to Trimmer for Vertex not correct, mapping says its " << Transformations::mapToTrimmer(2, i)
        << " but should be " << vertexIndexBackMapping[i];
    ;
    t.check(Transformations::mapToTrimmer(1, i) == edgeIndexBackMapping[i])
        << "Map to Trimmer for Edge not correct, mapping says its " << Transformations::mapToTrimmer(1, i)
        << " but should be " << edgeIndexBackMapping[i];
    ;
  }

  return t;
}

#include <cfenv>
int main(int argc, char** argv) try {
  feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

  // Initialize MPI, if necessary
  Dune::MPIHelper::instance(argc, argv);
  Dune::TestSuite t("", Dune::TestSuite::ThrowPolicy::ThrowOnRequired);

  t.subTest(testTransformations());

  t.report();

  return t.exit();
} catch (Dune::Exception& e) {
  std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
}
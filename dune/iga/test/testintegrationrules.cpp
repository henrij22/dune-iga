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
#include <dune/iga/patchgrid.hh>
#include <dune/iga/trimmer/defaulttrimmer/integrationrules/simplexintegrationrulegenerator.hh>
#include <dune/iga/trimmer/defaulttrimmer/trimmer.hh>

using namespace Dune;
using namespace Dune::IGANEW;

auto testElement1(Dune::TestSuite& t, int refLevel) {
  constexpr int gridDim  = 2;
  constexpr int dimworld = 2;

  using PatchGrid   = PatchGrid<gridDim, dimworld, DefaultTrim::PatchGridFamily>;
  using GridFactory = Dune::GridFactory<PatchGrid>;

  using IntegrationRuleGenerator = DefaultTrim::SimplexIntegrationRuleGenerator<const PatchGrid>;

  auto gridFactory = GridFactory();
  gridFactory.insertTrimParameters(GridFactory::TrimParameterType{100});
  gridFactory.insertJson("auxiliaryfiles/element_trim.ibra", true, {refLevel, refLevel});

  const auto grid = gridFactory.createGrid();
  auto gv         = grid->leafGridView();

  double area{0};
  for (const auto& ele : elements(gv)) {
    auto qR  = IntegrationRuleGenerator::createIntegrationRule(ele, grid.get(), 2 * gridDim);
    auto geo = ele.geometry();

    for (auto [gp, w] : qR) {
      area += geo.integrationElement(gp) * w;
    }
  }
  std::cout << "Area: " << area << std::endl;
}

#include <cfenv>
int main(int argc, char** argv) try {
  // feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

  // Initialize MPI, if necessary
  Dune::MPIHelper::instance(argc, argv);
  Dune::TestSuite t("", Dune::TestSuite::ThrowPolicy::ThrowOnRequired);

  testElement1(t, 0);
  testElement1(t, 1);

  t.report();

  return t.exit();
} catch (Dune::Exception& e) {
  std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
}
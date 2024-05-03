// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#define DUNE_CHECK_BOUNDS
#define CHECK_RESERVEDVECTOR
#ifdef HAVE_CONFIG_H
  #include "config.h"
#endif

#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/test/basistest.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/iga/hierarchicpatch/patchgrid.hh>
#include <dune/iga/nurbsbasis.hh>
#include <dune/iga/trimmer/defaulttrimmer/trimmer.hh>
#include <dune/iga/trimmer/identitytrimmer/trimmer.hh>

using namespace Dune;

auto testNurbsBasis() {

  // Create test case
  using PatchGrid   = IGANEW::PatchGrid<2, 2, IGANEW::DefaultTrim::PatchGridFamily>;
  using GridFactory = Dune::GridFactory<PatchGrid>;

  auto igaGridFactory = GridFactory();
  igaGridFactory.insertJson("auxiliaryfiles/element_trim.ibra", true, {1, 1});
  auto grid = igaGridFactory.createGrid();

  auto gridView = grid->leafGridView();

  Dune::TestSuite test(TestSuite::ThrowPolicy::AlwaysThrow);



  using GridView = decltype(gridView);
  Dune::Functions::NurbsBasis<GridView> basis(gridView, gridView.impl().patchData());

  // Test open knot vectors
  std::cout << "  Testing B-spline basis with open knot vectors" << std::endl;

  {
    using namespace Functions::BasisFactory;
    // Check basis created via its constructor
    Functions::NurbsBasis<GridView> basis2(gridView, nurbs());
    test.subTest(checkBasis(basis2, EnableContinuityCheck(), EnableContinuityCheck()));
    Dune::Functions::forEachBoundaryDOF(basis2, [](auto&& localIndex) {});
  }

  // {
  //   // Check basis created via its constructor
  //   Functions::NurbsBasis<GridView> basis2(gridView);
  //   test.subTest(checkBasis(basis2, EnableContinuityCheck(), EnableContinuityCheck()));
  // }
  //
  // {
  //   // Check basis created via makeBasis
  //   using namespace Functions::BasisFactory;
  //   auto basis2 = makeBasis(gridView, nurbs());
  //   test.subTest(checkBasis(basis2, EnableContinuityCheck(), EnableContinuityCheck()));
  // }
  //
  // {
  //   // Check whether a B-Spline basis can be combined with other bases.
  //   using namespace Functions::BasisFactory;
  //   auto basis2 = makeBasis(gridView, power<2>(nurbs()));
  //   test.subTest(checkBasis(basis2, EnableContinuityCheck(), EnableContinuityCheck()));
  // }
  //
  // {
  //   grid->degreeElevateOnAllLevels({1, 1});
  //   auto gridViewNew = grid->leafGridView();
  //   // Check lower order basis created via its constructor
  //   using namespace Functions::BasisFactory;
  //   Functions::NurbsBasis<GridView> basis2(gridViewNew, nurbs(degreeElevate(1, 1)));
  //   test.subTest(checkBasis(basis2, EnableContinuityCheck(), EnableContinuityCheck()));
  // }
  //
  // {
  //   grid->degreeElevateOnAllLevels({0, 1});
  //   auto gridViewNew = grid->leafGridView();
  //   // Check lower order basis created via its constructor
  //   using namespace Functions::BasisFactory;
  //   Functions::NurbsBasis<GridView> basis2(gridViewNew, nurbs(degreeElevate(1, 0)));
  //   test.subTest(checkBasis(basis2, EnableContinuityCheck(), EnableContinuityCheck()));
  // }

  return test;
}

int main(int argc, char** argv) try {
  // Initialize MPI, if necessary
  MPIHelper::instance(argc, argv);

  TestSuite t("");
  t.subTest(testNurbsBasis());

  return t.exit();
} catch (Dune::Exception& e) {
  std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
}
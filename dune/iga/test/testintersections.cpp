// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later
#define DUNE_CHECK_BOUNDS
#ifdef HAVE_CONFIG_H
  #include "config.h"
#endif

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/fvector.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/common/tupleutility.hh>
#include <dune/common/tuplevector.hh>
#include <dune/iga/patchgrid.hh>
#include <dune/iga/trimmer/defaulttrimmer/trimmer.hh>

using namespace Dune::IGANEW;
using namespace Dune;

std::vector<GeometryType> geometryTypes(bool trimmed, int refLevel) {
  assert(refLevel == 0);
  if (not trimmed)
    return {GeometryTypes::cube(1), GeometryTypes::cube(1), GeometryTypes::cube(1), GeometryTypes::cube(1)};
  else
    return {GeometryTypes::cube(1), GeometryTypes::none(1), GeometryTypes::none(1), GeometryTypes::none(1), GeometryTypes::cube(1)};

}

std::vector<Dune::FieldVector<double, 2>> unitOuterNormals(bool trimmed, int refLevel) {
  assert(refLevel == 0);
  if (not trimmed)
    return {
        {-1,  0},
        { 1,  0},
        { 0, -1},
        { 0,  1}
    };
  DUNE_THROW(IOError, "");
}

template <typename PatchGrid>
auto testIntersections(auto& grid, bool trimmed, int refLevel) {
  Dune::TestSuite t("", Dune::TestSuite::ThrowPolicy::AlwaysThrow);

  using LeafGridView  = typename PatchGrid::LeafGridView;
  using LevelGridView = typename PatchGrid::LevelGridView;

  using EntitySeed = typename PatchGrid::template Codim<0>::EntitySeed;
  using NormalVector = Dune::FieldVector<double, 2>;

  // NumIntersection, UnitOuterNormal, centerUnitOuterNormal, outerNormal, geometryType, insideEntitySeed
  using ResultTuple = std::tuple<int, std::vector<NormalVector>, std::vector<NormalVector>, std::vector<NormalVector>,
                                 std::vector<GeometryType>, std::vector<EntitySeed>>;

  ResultTuple resLevel{};
  ResultTuple resLeaf{};

  auto gvs = Dune::TupleVector<LevelGridView, LeafGridView>(grid.levelGridView(grid.maxLevel()), grid.leafGridView());
  Hybrid::forEach(gvs, [&]<typename GV>(const GV& gridView) {
    ResultTuple resTuple{};

    int globalIntersectionCount{0};
    for (const auto& ele : elements(gridView)) {
      int eleIntersectionCount{};
      for (const auto& intersection : intersections(gridView, ele)) {
        ++eleIntersectionCount;
        ++globalIntersectionCount;

        std::get<1>(resTuple).push_back(intersection.unitOuterNormal({0.5}));
        std::get<2>(resTuple).push_back(intersection.centerUnitOuterNormal());
        std::get<3>(resTuple).push_back(intersection.outerNormal({0.5}));
        std::get<4>(resTuple).push_back(intersection.type());
        std::get<5>(resTuple).push_back(intersection.inside().seed());
      }
      t.check(eleIntersectionCount == ele.subEntities(1))
          << "There should be " << ele.subEntities(1) << " intersections, but there are " << eleIntersectionCount;
    }

    std::get<0>(resTuple) = globalIntersectionCount;
    // Save appropriate tuple
    if constexpr (std::is_same_v<GV, LeafGridView>)
      resLeaf = resTuple;
    else
      resLevel = resTuple;
  });

  // Cehck level == leaf
  t.check(std::get<0>(resLevel) == std::get<0>(resLeaf));

  const auto expectedOuterNormals = unitOuterNormals(trimmed, refLevel);
  auto expectedGeometryTypes      = geometryTypes(trimmed, refLevel);
  for (const auto i : Dune::range(std::get<0>(resLevel))) {
    // Check outerNormals
    auto lvlUnitOuterNormal  = std::get<1>(resLevel)[i];
    auto leafUnitOuterNormal = std::get<1>(resLeaf)[i];

    auto lvlCenterOuterNormal  = std::get<2>(resLevel)[i];
    auto leafCenterOuterNormal = std::get<2>(resLeaf)[i];

    auto lvlOuterNormal  = std::get<3>(resLevel)[i];
    auto leafOuterNormal = std::get<3>(resLeaf)[i];

    t.check(FloatCmp::eq(lvlUnitOuterNormal, leafUnitOuterNormal));
    t.check(FloatCmp::eq(lvlCenterOuterNormal, leafCenterOuterNormal));
    t.check(FloatCmp::eq(lvlOuterNormal, leafOuterNormal));

    t.check(FloatCmp::eq(lvlUnitOuterNormal, expectedOuterNormals[i]));
    t.check(FloatCmp::eq(lvlCenterOuterNormal, expectedOuterNormals[i]));

    // Check unit length
    t.check(FloatCmp::eq(lvlUnitOuterNormal.two_norm(), 1.0));

    // Check linear dependence of unitOuterNormal und outerNormal
    auto testJ = Dune::FieldMatrix<double, 2>{lvlUnitOuterNormal, lvlOuterNormal};
    t.check(FloatCmp::eq(testJ.determinant(), 0.0));

    // Check geometryTypes
    auto levelGeoType = std::get<4>(resLevel)[i];
    auto leafGeoType  = std::get<4>(resLeaf)[i];

    t.check(levelGeoType == leafGeoType);
    t.check(levelGeoType == expectedGeometryTypes[i]);

    // Seed
    auto levelEntityInside = grid.entity(std::get<5>(resLevel)[i]);
    auto leafEntityInside = grid.entity(std::get<5>(resLeaf)[i]);

    t.check(levelEntityInside == leafEntityInside);
  }

  return t;
}

auto makeTestCase(Dune::TestSuite& t, bool trimmed, int refLevel) {
  constexpr int gridDim  = 2;
  constexpr int dimworld = 2;

  using PatchGrid   = PatchGrid<gridDim, dimworld, DefaultTrim::PatchGridFamily>;
  using GridFactory = Dune::GridFactory<PatchGrid>;

  auto gridFactory = GridFactory();
  gridFactory.insertTrimParameters(GridFactory::TrimParameterType{100});
  gridFactory.insertJson("auxiliaryfiles/element_trim.ibra", trimmed, {refLevel, refLevel});

  const auto grid = gridFactory.createGrid();
  t.subTest(testIntersections<PatchGrid>(*grid, trimmed, refLevel));
}

#include <cfenv>
int main(int argc, char** argv) try {
  // feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

  /******************
   * These are unit tests for the intersection implementation, we are mostly testing geometry, and not correct indices,
   * etc, as this is tested by the gridtests
   ******************/

  // Initialize MPI, if necessary
  Dune::MPIHelper::instance(argc, argv);
  Dune::TestSuite t("", Dune::TestSuite::ThrowPolicy::ThrowOnRequired);

  makeTestCase(t, false, 0);
  // makeTestCase(t, true, 0);

  t.report();

  return t.exit();
} catch (Dune::Exception& e) {
  std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
}
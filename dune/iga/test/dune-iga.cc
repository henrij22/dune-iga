// SPDX-FileCopyrightText: 2022 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later

// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/test/checkentitylifetime.hh>
#include <dune/grid/test/checkiterators.hh>
#include <dune/grid/test/checkjacobians.hh>
#include <dune/functions/functionspacebases/test/basistest.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/iga/nurbsgrid.hh>
#include <dune/iga/nurbspatchgeometry.h>
#include <dune/iga/nurbstrimutils.hh>
#include <dune/iga/ibraReader.hh>

#define TEST_ALL

// Uncomment for grid and geometry checks -- makes copilation time way longer
#ifdef TEST_ALL

#include <dune/grid/test/checkgeometry.hh>
#include <dune/grid/test/checkindexset.hh>
#include <dune/grid/test/gridcheck.hh>
#include <dune/iga/gridcapabilities.hh>
#include <dune/grid/io/file/printgrid.hh>
#include <dune/iga/nurbstrimfunctionality.hh>

#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>

#include <dune/functions/functionspacebases/flatmultiindex.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/test/checkgeometry.hh>
#include <dune/grid/io/file/printgrid.hh>
#include <dune/iga/dunelinearalgebratraits.hh>
#include <dune/iga/nurbspatch.hh>
#include <dune/iga/nurbsbasis.hh>
#include <dune/iga/igaDataCollector.h>
#include <dune/vtk/vtkwriter.hh>

#endif

#include "plotFunctionality.h"



using namespace Dune;
using namespace Dune::IGA;

#if 0
template <typename T, int worldDim, int Items>
struct Compare {
  constexpr bool operator()(const std::array<FieldVector<double, worldDim>, Items>& lhs,
                            const std::array<FieldVector<double, worldDim>, Items>& rhs) const {
    return std::ranges::lexicographical_compare(std::ranges::join_view(lhs), std::ranges::join_view(rhs));
  };
};

auto checkUniqueEdges(const auto& gridView) {
  TestSuite t;

  constexpr int gridDimensionworld = std::remove_cvref_t<decltype(gridView)>::dimensionworld;
  constexpr int gridDimension      = std::remove_cvref_t<decltype(gridView)>::dimension;
  std::set<std::array<FieldVector<double, gridDimensionworld>, 2>, Compare<double, gridDimensionworld, 2>>
      edgeVertexPairSet;
  for (int eleIndex = 0; auto&& element : elements(gridView)) {
    edgeVertexPairSet.clear();
    for (int edgeIndex = 0; edgeIndex < element.subEntities(gridDimension - 1); ++edgeIndex) {
      auto edge = element.template subEntity<gridDimension - 1>(edgeIndex);
      std::array<FieldVector<double, gridDimensionworld>, 2> pair;
      for (int c = 0; c < edge.geometry().corners(); ++c)
        pair[c] = edge.geometry().corner(c);
      bool inserted = edgeVertexPairSet.insert(pair).second;
      t.require(inserted) << "Duplicate edge detected in Element " << eleIndex << " Edges: " << pair[0] << ", "
                          << pair[1];
    }
    ++eleIndex;
  }
  return t;
}

auto checkUniqueSurfaces(const auto& gridView) {
  TestSuite t;

  constexpr int gridDimensionworld = std::remove_cvref_t<decltype(gridView)>::dimensionworld;
  std::set<std::array<FieldVector<double, gridDimensionworld>, 4>, Compare<double, gridDimensionworld, 4>>
      edgeVertexPairSet;
  for (int eleIndex = 0; auto&& element : elements(gridView)) {
    edgeVertexPairSet.clear();
    for (int edgeIndex = 0; edgeIndex < element.subEntities(2); ++edgeIndex) {
      auto edge = element.template subEntity<2>(edgeIndex);
      std::array<FieldVector<double, gridDimensionworld>, 4> tuple;
      for (int c = 0; c < edge.geometry().corners(); ++c)
        tuple[c] = edge.geometry().corner(c);

      t.require(edgeVertexPairSet.insert(tuple).second)
          << "Duplicate surface detected in Element " << eleIndex << " Surfaces: " << tuple[0] << ", " << tuple[1]
          << ", " << tuple[2] << ", " << tuple[3];
    }
    ++eleIndex;
  }
  return t;
}

auto testNURBSGridCurve() {
  ////////////////////////////////////////////////////////////////
  //  Second test
  //  A B-Spline curve of dimWorld 3
  ////////////////////////////////////////////////////////////////

  // parameters
  int subSampling = 80;

  const auto dim      = 1;
  const auto dimworld = 3;

  const std::array<int, dim> order                     = {3};
  const std::array<std::vector<double>, dim> knotSpans = {{{0, 0, 0, 0, 1, 3, 4, 4, 4, 5, 5, 5, 5}}};
  //  const std::array<std::vector<double>, dim> knotSpans = {{{ 0, 0, 1,1}}};
  using ControlPoint = Dune::IGA::NURBSPatchData<dim, dimworld>::ControlPointType;

  const std::vector<ControlPoint> endPoints
      = {{.p = {1, 3, 4}, .w = 1}, {.p = {2, 2, 2}, .w = 3}, {.p = {3, 4, 5}, .w = 1},
         {.p = {5, 1, 7}, .w = 2}, {.p = {4, 7, 2}, .w = 1}, {.p = {8, 6, 2}, .w = 1},
         {.p = {2, 9, 9}, .w = 7}, {.p = {1, 4, 3}, .w = 1}, {.p = {1, 7, 1}, .w = 5}};

  std::array<int, dim> dimsize = {static_cast<int>(endPoints.size())};
  auto controlNet              = Dune::IGA::NURBSPatchData<dim, dimworld>::ControlPointNetType(dimsize, endPoints);

  Dune::IGA::NURBSPatchData<dim, dimworld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = order;
  patchData.controlPoints = controlNet;
  auto additionalKnots    = std::vector<double>(2);
  additionalKnots[0]      = 0.5;
  additionalKnots[1]      = 3.5;
  patchData               = knotRefinement<dim>(patchData, additionalKnots, 0);
  patchData               = degreeElevate(patchData, 0, 1);

  TestSuite t;
  IGA::NURBSGrid<dim, dimworld> grid(patchData);
  grid.globalRefine(3);
  auto gridView        = grid.leafGridView();
  const auto& indexSet = gridView.indexSet();

  for (int eleIndex = 0;
       auto& ele :
       elements(gridView))  // This test also exists in grid check, but it is more convenient to debug it here
  {
    const int numCorners = ele.subEntities(dim);
    for (int c = 0; c < numCorners; ++c) {
      auto vertex = ele.template subEntity<dim>(c).geometry();
      auto elegeo = ele.geometry();
      t.check(Dune::FloatCmp::eq((elegeo.corner(c) - vertex.corner(0)).two_norm(), 0.0))
          << "Corner " << c << " " << elegeo.corner(c) << " Alt: " << vertex.corner(0);
    }
    ++eleIndex;
  }

  Dune::RefinementIntervals refinementIntervals1(subSampling);
  SubsamplingVTKWriter<decltype(gridView)> vtkWriter(gridView, refinementIntervals1);
  //  vtkWriter.write("NURBSGridTest-CurveNewFineResample");
  vtkWriter.write("NURBSGridTest-CurveNewFineResample-R");
  //  vtkWriter.write("NURBSGridTest-CurveNewFineResample_knotRefine");
  gridcheck(grid);
  return t;
}

void testNURBSGridSurface() {
  TestSuite t;
  int subSampling = 10;

  //////////////////////////////////////////////////////////////
  // Create a 2d NURBS surface in 3d
  //////////////////////////////////////////////////////////////

  const auto dim                   = 2;
  const auto dimworld              = 3;
  const std::array<int, dim> order = {2, 2};

  const std::array<std::vector<double>, dim> knotSpans = {{{0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1}}};
  //  const std::vector<std::vector<FieldVector<double, dimworld> > > controlPointsold
  //      = {{{0, 0, 1}, {1, 0, 1}, {2, 0, 2}}, {{0, 1, 0}, {1, 1, 0}, {2, 1, 0}}, {{0, 2, 1}, {1, 2, 2}, {2, 2, 2}}};
  using ControlPoint = Dune::IGA::NURBSPatchData<dim, dimworld>::ControlPointType;

  const std::vector<std::vector<ControlPoint>> endPoints
      = {{{.p = {0, 0, 1}, .w = 2}, {.p = {1, 0, 1}, .w = 2}, {.p = {2, 0, 2}, .w = 1}},
         {{.p = {0, 1, 0}, .w = 1}, {.p = {1, 1, 0}, .w = 4}, {.p = {2, 1, 0}, .w = 1}},
         {{.p = {0, 2, 1}, .w = 1}, {.p = {1, 2, 2}, .w = 2}, {.p = {2, 2, 2}, .w = 4}}};

  std::array<int, dim> dimsize = {static_cast<int>(endPoints.size()), static_cast<int>(endPoints[0].size())};

  auto controlNet = Dune::IGA::NURBSPatchData<dim, dimworld>::ControlPointNetType(dimsize, endPoints);

  Dune::IGA::NURBSPatchData<dim, dimworld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = order;
  patchData.controlPoints = controlNet;
  IGA::NURBSGrid<dim, dimworld> grid(patchData);
  auto gridView        = grid.leafGridView();
  const auto& indexSet = gridView.indexSet();

  Dune::RefinementIntervals refinementIntervals1(subSampling);
  SubsamplingVTKWriter<decltype(gridView)> vtkWriter(gridView, refinementIntervals1);
  vtkWriter.write("NURBSGridTest-Surface");

  t.subTest(checkUniqueEdges(gridView));
}

auto test3DGrid() {
  constexpr std::size_t dim        = 3;
  constexpr std::size_t dimworld   = 3;
  const std::array<int, dim> order = {2, 2, 2};
  // quarter cylindrical hyperSurface
  const double lx = 2;
  const double ly = 1;
  const double lz = 1;

  using ControlPoint = Dune::IGA::NURBSPatchData<dim, dimworld>::ControlPointType;
  Dune::IGA::NURBSPatchData<dim, dimworld> nurbsPatchData;
  nurbsPatchData.knotSpans = {{{0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1}}};

  std::vector<std::vector<std::vector<ControlPoint>>> controlp;
  std::array<int, dim> dimSize = {3, 3, 3};
  for (int i = 0; i < dimSize[0]; ++i) {
    controlp.emplace_back();
    for (int j = 0; j < dimSize[1]; ++j) {
      controlp[i].emplace_back();
      for (int k = 0; k < dimSize[2]; ++k) {
        controlp[i][j].push_back(
            {.p = {i * i * lx / (dimSize[0] - 1) + 1, 2 * i * j * k * ly / (dimSize[1] - 1) + (k + 1) + (j + 1),
                   k * k * lz / (dimSize[2] - 1)},
             .w = 1});
      }
    }
  }

  nurbsPatchData.controlPoints = MultiDimensionNet(dimSize, controlp);
  nurbsPatchData.degree        = order;

  auto additionalKnots = std::vector<double>(2);
  additionalKnots[0]   = 0.1;
  additionalKnots[1]   = 0.3;
  //  additionalKnots[2] = 0.6;
  //  additionalKnots[1] = 3.5;
  //  nurbsPatchData = knotRefinement<dim>(nurbsPatchData, additionalKnots, 2);
  //  nurbsPatchData = degreeElevate(nurbsPatchData,0,1);
  //  nurbsPatchData = degreeElevate(nurbsPatchData, 1, 2);
  //  nurbsPatchData = degreeElevate(nurbsPatchData,2,1);
  IGA::NURBSGrid<3, 3> grid(nurbsPatchData);
  //  grid.globalRefine(1);
  //  gridcheck(grid);
  //  grid.globalRefineInDirection(0,1);
  //  gridcheck(grid);
  //  grid.globalRefineInDirection(1,2);
  //  gridcheck(grid);
  //  grid.globalRefineInDirection(2, 3);
  //  gridcheck(grid);

  auto gridView = grid.leafGridView();
  TestSuite t;
  t.subTest(checkUniqueEdges(gridView));
  t.subTest(checkUniqueSurfaces(gridView));

  const int subSampling = 10;
  Dune::RefinementIntervals refinementIntervals1(subSampling);
  SubsamplingVTKWriter vtkWriter(gridView, refinementIntervals1);
  vtkWriter.write("NURBSGridTest-Solid2");

  Dune::GeometryChecker<decltype(grid)> geometryChecker;
  geometryChecker.checkGeometry(gridView);
  Dune::checkIndexSet(grid, gridView, std::cout);

  checkEntityLifetime(gridView, gridView.size(0));

  for (auto&& elegeo : elements(gridView) | std::views::transform([](const auto& ele) { return ele.geometry(); }))
    checkJacobians(elegeo);

  checkIterators(gridView);

  return t;
}

void testFactory() {
  constexpr auto dim               = 2UL;
  constexpr auto dimworld          = 3UL;
  const std::array<int, dim> order = {2, 1};
  const double invsqr2             = 1.0 / std::sqrt(2.0);
  // quarter cylindrical surface
  const double l   = 10;
  const double rad = 5;
  //  const std::array<std::vector<double>, dim> knotSpans = {{{0, 0, 0,0.5, 1, 1, 1}, {0, 0, 1, 1}}};

  Dune::IGA::NURBSPatchData<dim, dimworld> nurbsPatchData;
  nurbsPatchData.knotSpans = {{{0, 0, 0, 1, 1, 1}, {0, 0, 1, 1}}};

  using ControlPoint = Dune::IGA::NURBSPatchData<dim, dimworld>::ControlPointType;
  const std::vector<std::vector<ControlPoint>> endPoints
      = {{{.p = {0, 0, rad}, .w = 1}, {.p = {0, l, rad}, .w = 1}},
         {{.p = {rad, 0, rad}, .w = invsqr2}, {.p = {rad, l, rad}, .w = invsqr2}},
         //          {{.p = {rad*2, 0,   0}, .w =       1},  {.p = {rad*2, l*2,   0}, .w = 1     }},
         {{.p = {rad, 0, 0}, .w = 1}, {.p = {rad, l, 0}, .w = 1}}};

  std::array<int, dim> dimsize = {static_cast<int>(endPoints.size()), static_cast<int>(endPoints[0].size())};
  nurbsPatchData.controlPoints = Dune::IGA::NURBSPatchData<dim, dimworld>::ControlPointNetType(dimsize, endPoints);
  nurbsPatchData.degree        = order;
};

auto testTorusGeometry() {
  const double R       = 2.0;
  const double r       = 1.0;
  auto circle          = makeCircularArc(r);
  auto nurbsPatchData  = makeSurfaceOfRevolution(circle, {R, 0, 0}, {0, 1, 0}, 360.0);
  nurbsPatchData       = degreeElevate(nurbsPatchData, 0, 1);
  nurbsPatchData       = degreeElevate(nurbsPatchData, 1, 2);
  auto additionalKnots = std::vector<double>(1);
  additionalKnots[0]   = 0.1;
  nurbsPatchData       = knotRefinement<2>(nurbsPatchData, additionalKnots, 1);

  IGA::NURBSGrid<2, 3> grid(nurbsPatchData);
  grid.globalRefine(1);
  grid.globalRefineInDirection(1, 1);
  grid.globalRefineInDirection(0, 2);

  auto gridView = grid.leafGridView();

  const int subSampling = 2;
  Dune::RefinementIntervals refinementIntervals1(subSampling);
  SubsamplingVTKWriter vtkWriter(gridView, refinementIntervals1);
  vtkWriter.write("NURBSGridTest-SurfaceRevolutionFLAT");

  TestSuite test;
  Dune::GeometryChecker<IGA::NURBSGrid<2UL, 3UL>> geometryChecker;
  //  geometryChecker.checkGeometry(gridView);
  Dune::checkIndexSet(grid, gridView, std::cout);

  double area = 0.0;
  for (auto&& ele : elements(gridView)) {
    area += ele.geometry().volume();
  }
  const double pi                        = std::numbers::pi_v<double>;
  const double referenceTorusSurfaceArea = 4.0 * pi * pi * r * R;
  test.check(area - referenceTorusSurfaceArea < 1e-4, "The integrated area of the torus hyperSurface is wrong!");

  double gaussBonnet = 0.0;
  for (auto& ele : elements(gridView)) {
    const auto rule
        = Dune::QuadratureRules<double, 2>::rule(ele.type(), 2 * (*std::ranges::max_element(nurbsPatchData.degree)));
    for (auto& gp : rule) {
      const auto Kinc = ele.geometry().impl().gaussianCurvature(gp.position());
      const auto Kmax = 1 / (r * (R + r));
      const auto Kmin = -1 / (r * (R - r));
      test.check(Kinc < Kmax && Kinc > Kmin, "The Gaussian curvature should be within bounds");
      gaussBonnet += Kinc * gp.weight() * ele.geometry().integrationElement(gp.position());
    }
  }

  test.check(std::abs(gaussBonnet) < 1e-5,
             "Gauss-Bonnet theorem dictates a vanishing integrated Gaussian curvature for the torus!");
  checkEntityLifetime(gridView, gridView.size(0));

  for (auto&& elegeo : elements(gridView) | std::views::transform([](const auto& ele) { return ele.geometry(); }))
    checkJacobians(elegeo);

  checkIterators(gridView);

  gridcheck(grid);
  test.subTest(checkUniqueEdges(gridView));
  return test;
}

auto testNURBSSurface() {
  // parameters
  int subSampling = 5;

  //////////////////////////////////////////////////////////////
  // Create a 2d NURBS surface in 3d
  //////////////////////////////////////////////////////////////

  const auto dim                   = 2UL;
  const auto dimworld              = 3UL;
  const std::array<int, dim> order = {2, 2};

  const std::array<std::vector<double>, dim> knotSpans = {{{0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1}}};
  using ControlPoint                                   = Dune::IGA::NURBSPatchData<dim, dimworld>::ControlPointType;

  const std::vector<std::vector<ControlPoint>> endPoints
      = {{{.p = {0, 0, 1}, .w = 2}, {.p = {1, 0, 1}, .w = 2}, {.p = {2, 0, 2}, .w = 1}},
         {{.p = {0, 1, 0}, .w = 1}, {.p = {1, 1, 0}, .w = 4}, {.p = {2, 1, 0}, .w = 1}},
         {{.p = {0, 2, 1}, .w = 1}, {.p = {1, 2, 2}, .w = 2}, {.p = {2, 2, 2}, .w = 4}}};

  std::array<int, dim> dimsize = {static_cast<int>(endPoints.size()), static_cast<int>(endPoints[0].size())};
  //

  //  auto weightNet  = MultiDimensionNet<dim, double>(dimsize, weight);
  auto controlNet = Dune::IGA::NURBSPatchData<dim, dimworld>::ControlPointNetType(dimsize, endPoints);

  IGA::NURBSPatch<dim, dimworld> patch(knotSpans, controlNet, order);

  TestSuite testSuite;
  testSuite.check(patch.size(0) == 1);
  testSuite.check(patch.size(1) == 4);
  testSuite.check(patch.size(2) == 9);

  return testSuite;
}

auto testNURBSCurve() {
  // parameters
  unsigned int subSampling = 5;

  ////////////////////////////////////////////////////////////////
  //  Create a B-spline curve in 3d
  ////////////////////////////////////////////////////////////////

  const int dim      = 1;
  const int dimworld = 3;

  const std::array<int, dim> order                     = {2};
  const std::array<std::vector<double>, dim> knotSpans = {{{0, 0, 0, 1, 1, 2, 3, 4, 4, 5, 5, 5}}};

  using ControlPoint = Dune::IGA::NURBSPatchData<dim, dimworld>::ControlPointType;
  const std::vector<ControlPoint> endPoints
      = {{.p = {1, 3, 4}, .w = 2}, {.p = {2, 2, 2}, .w = 2}, {.p = {3, 4, 5}, .w = 1},
         {.p = {5, 1, 7}, .w = 1}, {.p = {4, 7, 2}, .w = 4}, {.p = {8, 6, 2}, .w = 2},
         {.p = {2, 9, 9}, .w = 1}, {.p = {1, 4, 3}, .w = 2}, {.p = {1, 7, 1}, .w = 4}};

  std::array<int, dim> dimsize = {static_cast<int>(endPoints.size())};
  auto controlNet              = Dune::IGA::NURBSPatchData<dim, dimworld>::ControlPointNetType(dimsize, endPoints);

  IGA::NURBSPatch<dim, dimworld> patch(knotSpans, controlNet, order);

  TestSuite testSuite;
  testSuite.check(patch.size(0) == 5);
  testSuite.check(patch.size(1) == endPoints.size());
  return testSuite;
}

void testNurbsGridCylinder() {
  ////////////////////////////////////////////////////////////////
  //  First test
  //  A B-Spline surface of dimWorld 3
  ////////////////////////////////////////////////////////////////

  // parameters
  int subSampling = 1;

  //////////////////////////////////////////////////////////////
  // Create a 2d B-spline grid in 3d
  //////////////////////////////////////////////////////////////
  constexpr auto dim               = 2UL;
  constexpr auto dimworld          = 3UL;
  const std::array<int, dim> order = {2, 1};
  const double invsqr2             = 1.0 / std::sqrt(2.0);
  // quarter cylindrical surface
  const double l   = 10;
  const double rad = 5;
  //  const std::array<std::vector<double>, dim> knotSpans = {{{0, 0, 0,0.5, 1, 1, 1}, {0, 0, 1, 1}}};
  const std::array<std::vector<double>, dim> knotSpans = {{{0, 0, 0, 1, 1, 1}, {0, 0, 1, 1}}};

  using ControlPoint = Dune::IGA::NURBSPatchData<dim, dimworld>::ControlPointType;
  const std::vector<std::vector<ControlPoint>> endPoints
      = {{{.p = {0, 0, rad}, .w = 1}, {.p = {0, l, rad}, .w = 1}},
         {{.p = {rad, 0, rad}, .w = invsqr2}, {.p = {rad, l, rad}, .w = invsqr2}},
         //          {{.p = {rad*2, 0,   0}, .w =       1},  {.p = {rad*2, l*2,   0}, .w = 1     }},
         {{.p = {rad, 0, 0}, .w = 1}, {.p = {rad, l, 0}, .w = 1}}};

  std::array<int, dim> dimsize = {static_cast<int>(endPoints.size()), static_cast<int>(endPoints[0].size())};
  auto controlNet              = Dune::IGA::NURBSPatchData<dim, dimworld>::ControlPointNetType(dimsize, endPoints);

  Dune::IGA::NURBSPatchData<dim, dimworld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = order;
  patchData.controlPoints = controlNet;
  IGA::NURBSGrid<dim, dimworld> grid(patchData);
  grid.globalRefine(5);
  auto gridView        = grid.leafGridView();
  const auto& indexSet = gridView.indexSet();

  TestSuite testSuite;

  IGA::NURBSPatch<dim, dimworld> nurbsPatch(knotSpans, controlNet, order);

  testSuite.check(nurbsPatch.size(0) == 1);
  testSuite.check(nurbsPatch.size(1) == 4);
  testSuite.check(nurbsPatch.size(2) == 4);

  //! Test code for VTKWriter, please uncomment to inspect the remaining errors
  Dune::RefinementIntervals refinementIntervals1(subSampling);
  VTKWriter vtkWriter(gridView);

  vtkWriter.write("ZylRefine");

  testSuite.subTest(checkUniqueEdges(gridView));
}

// template <int dim>
auto testNurbsBasis() {
  ////////////////////////////////////////////////////////////////
  //  First test
  //  A B-Spline surface of dimWorld 3
  ////////////////////////////////////////////////////////////////

  // parameters
  int subSampling = 1;

  //////////////////////////////////////////////////////////////
  // Create a 2d B-spline grid in 3d
  //////////////////////////////////////////////////////////////
  constexpr auto dim               = 2UL;
  constexpr auto dimworld          = 3UL;
  const std::array<int, dim> order = {2, 1};
  const double invsqr2             = 1.0 / std::sqrt(2.0);
  // quarter cylindrical surface
  const double l   = 10;
  const double rad = 5;
  //  const std::array<std::vector<double>, dim> knotSpans = {{{0, 0, 0,0.5, 1, 1, 1}, {0, 0, 1, 1}}};

  using ControlPoint = Dune::IGA::NURBSPatchData<dim, dimworld>::ControlPointType;
  Dune::IGA::NURBSPatchData<dim, dimworld> nurbsPatchData;
  nurbsPatchData.knotSpans = {{{0, 0, 0, 1, 1, 1}, {0, 0, 1, 1}}};

  nurbsPatchData.controlPoints
      = {{{.p = {0, 0, rad}, .w = 1}, {.p = {0, l, rad}, .w = 1}},
         {{.p = {rad, 0, rad}, .w = invsqr2}, {.p = {rad, l, rad}, .w = invsqr2}},
         //          {{.p = {rad*2, 0,   0}, .w =       1},  {.p = {rad*2, l*2,   0}, .w = 1     }},
         {{.p = {rad, 0, 0}, .w = 1}, {.p = {rad, l, 0}, .w = 1}}};
  nurbsPatchData.degree = order;

  IGA::NURBSGrid<dim, dimworld> grid(nurbsPatchData);
  //  grid.globalRefine(1);
  grid.globalRefineInDirection(0, 2);
  //  grid.globalRefineInDirection(1, 3);
  auto gridView        = grid.leafGridView();
  const auto& indexSet = gridView.indexSet();

  TestSuite test;

  //! Test code for VTKWriter, please uncomment to inspect the remaining errors
  Dune::RefinementIntervals refinementIntervals1(subSampling);
  SubsamplingVTKWriter<decltype(gridView)> vtkWriter(gridView, refinementIntervals1);

  vtkWriter.write("ZylRefine");
  using GridView = decltype(gridView);
  Dune::Functions::NurbsBasis<GridView> basis(gridView, gridView.impl().getPatchData());

  // Test open knot vectors
  std::cout << "  Testing B-spline basis with open knot vectors" << std::endl;

  {
    // Check basis created via its constructor
    Functions::NurbsBasis<GridView> basis2(gridView, gridView.impl().getPatchData());
    test.subTest(checkBasis(basis2, EnableContinuityCheck(), EnableContinuityCheck()));
  }

  {
    // Check basis created via its constructor
    Functions::NurbsBasis<GridView> basis2(gridView);
    test.subTest(checkBasis(basis2, EnableContinuityCheck(), EnableContinuityCheck()));
  }

  {
    // Check basis created via makeBasis
    using namespace Functions::BasisFactory;
    auto basis2 = makeBasis(gridView, nurbs<dim>(gridView.impl().getPatchData()));
    test.subTest(checkBasis(basis2, EnableContinuityCheck(), EnableContinuityCheck()));
  }

  {
    // Check whether a B-Spline basis can be combined with other bases.
    using namespace Functions::BasisFactory;
    auto basis2 = makeBasis(gridView, power<2>(gridView.impl().getPreBasis()));
    test.subTest(checkBasis(basis2, EnableContinuityCheck(), EnableContinuityCheck()));
  }
  return test;
}

auto testBsplineBasisFunctions() {
  std::vector<double> knots = {0, 0, 0, 0.5, 0.5, 2, 2, 3, 3, 3};
  int degree                = 2;
  TestSuite test;

  auto N = Dune::IGA::BsplineBasis1D<double>::basisFunctions(0.3, knots, degree);
  using Dune::FloatCmp::eq;
  test.check(eq(N[0], 0.16), "P=2,N0,u=0.3");
  test.check(eq(N[1], 0.48), "P=2,N1,u=0.3");
  test.check(eq(N[2], 0.36), "P=2,N2,u=0.3");

  N = Dune::IGA::BsplineBasis1D<double>::basisFunctions(0.5, knots, degree);  // try knot span boundary

  test.check(eq(N[0], 1.0), "P=2,N0,u=0.5");
  test.check(eq(N[1], 0.0), "P=2,N1,u=0.5");
  test.check(eq(N[2], 0.0), "P=2,N2,u=0.5");

  N = Dune::IGA::BsplineBasis1D<double>::basisFunctions(0.0, knots, degree);  // try left end

  test.check(eq(N[0], 1.0), "P=2,N0,u=0.0");
  test.check(eq(N[1], 0.0), "P=2,N1,u=0.0");
  test.check(eq(N[2], 0.0), "P=2,N2,u=0.0");

  N = Dune::IGA::BsplineBasis1D<double>::basisFunctions(3.0, knots, degree);  // try right end

  test.check(eq(N[0], 0.0), "P=2,N0,u=3");
  test.check(eq(N[1], 0.0), "P=2,N1,u=3");
  test.check(eq(N[2], 1.0), "P=2,N2,u=3");

  knots  = {0, 0, 0, 0.5, 1, 1, 1};
  degree = 2;

  N = Dune::IGA::BsplineBasis1D<double>::basisFunctions(0.1, knots, degree);
  test.check(eq(N[0], 0.64), "P=2,N0,u=0.1");
  test.check(eq(N[1], 0.34), "P=2,N1,u=0.1");
  test.check(eq(N[2], 0.02), "P=2,N2,u=0.1");

  N = Dune::IGA::BsplineBasis1D<double>::basisFunctions(0.01, knots, degree);
  test.check(eq(N[0], 0.9604), "P=2,N0,u=0.01");
  test.check(eq(N[1], 0.0394), "P=2,N1,u=0.01");
  test.check(eq(N[2], 0.0002), "P=2,N2,u=0.01");

  knots  = {0, 0, 0, 0, 0.5, 1, 1, 2, 2, 2, 2};
  degree = 3;

  auto N2 = Dune::IGA::BsplineBasis1D<double>::basisFunctions(1.45, knots, degree);
  test.check(eq(N2[0], 0.1109166666666667), "P=3,N0,u=1.45");
  test.check(eq(N2[1], 0.4638333333333333), "P=3,N1,u=1.45");
  test.check(eq(N2[2], 0.334125), "P=3,N2,u=1.45");
  test.check(eq(N2[3], 0.091125), "P=3,N3,u=1.45");

  auto dN = Dune::IGA::BsplineBasis1D<double>::basisFunctionDerivatives(1.45, knots, degree, 3);
  // check values
  test.check(eq(dN[0][0], 0.1109166666666667), "P=3,dN00,u=1.45");
  test.check(eq(dN[0][1], 0.4638333333333333), "P=3,dN01,u=1.45");
  test.check(eq(dN[0][2], 0.334125), "P=3,dN02,u=1.45");
  test.check(eq(dN[0][3], 0.091125), "P=3,dN03,u=1.45");

  // check first derivatives
  test.check(eq(dN[1][0], -0.605), "P=3,dN10,u=1.45");
  test.check(eq(dN[1][1], -0.88), "P=3,dN11,u=1.45");
  test.check(eq(dN[1][2], 0.8775), "P=3,dN12,u=1.45");
  test.check(eq(dN[1][3], 0.6075), "P=3,dN13,u=1.45");

  // check second derivatives
  test.check(eq(dN[2][0], 2.2), "P=3,dN20,u=1.45");
  test.check(eq(dN[2][1], -2.8), "P=3,dN21,u=1.45");
  test.check(eq(dN[2][2], -2.1), "P=3,dN22,u=1.45");
  test.check(eq(dN[2][3], 2.7), "P=3,dN23,u=1.45");

  // check third derivatives
  test.check(eq(dN[3][0], -4.0), "P=3,dN30,u=1.45");
  test.check(eq(dN[3][1], 16.0), "P=3,dN31,u=1.45");
  test.check(eq(dN[3][2], -18.0), "P=3,dN32,u=1.45");
  test.check(eq(dN[3][3], 6.0), "P=3,dN33,u=1.45");
  // https://godbolt.org/z/Ta3fzW553
  auto Nf                           = Dune::IGA::BsplineBasis1D<double>(knots, degree);
  std::vector<double> NAtEvalPoints = {1,
                                       0.1714677640603567,
                                       0.001371742112482855,
                                       0.0740740740740741,
                                       0.00274348422496571,
                                       0.4682213077274805,
                                       0.1975308641975309,
                                       0.05852766346593506,
                                       0.007315957933241894,
                                       0};
  for (int i = 0; i < NAtEvalPoints.size(); ++i) {
    test.check(eq(Nf(i / (NAtEvalPoints.size() - 1.0) * 2.0)[0], NAtEvalPoints[i]));
  }

  std::array<double, 2> xieta{0.2, 0.25};
  std::array<std::vector<double>, 2> knots2 = {{{0, 0, 0, 0.5, 0.5, 2, 2, 3, 3, 3}, {0, 0, 0, 2, 2, 2}}};
  std::array<int, 2> degree2{2, 2};
  const std::vector<std::vector<double>> weights2
      = {{{1, 2, 3, 4, 5, 6, 7}, {8, 9, 10, 11, 12, 13, 14}, {15, 16, 17, 18, 19, 20, 21}}};
  std::array<int, 2> dimsize = {static_cast<int>(weights2.size()), static_cast<int>(weights2[0].size())};
  MultiDimensionNet<2UL, double> weightNet(dimsize, weights2);

  auto N_Nurbs = Dune::IGA::Nurbs<2>::basisFunctions(xieta, knots2, degree2, weightNet).directGetAll();

  test.check(N_Nurbs.size() == (degree2[0] + 1) * (degree2[1] + 1));

  test.check(eq(N_Nurbs[0], 0.04023722627737226), "Nurbs2d P=2,N0");  // check ansatzfunctions in domain
  test.check(eq(N_Nurbs[1], 0.4291970802919708), "Nurbs2d P=2,N1");
  test.check(eq(N_Nurbs[2], 0.2682481751824818), "Nurbs2d P=2,N2");
  test.check(eq(N_Nurbs[3], 0.02299270072992701), "Nurbs2d P=2,N3");
  test.check(eq(N_Nurbs[4], 0.137956204379562), "Nurbs2d P=2,N4");
  test.check(eq(N_Nurbs[5], 0.08175182481751825), "Nurbs2d P=2,N5");
  test.check(eq(N_Nurbs[6], 0.002463503649635036), "Nurbs2d P=2,N6");
  test.check(eq(N_Nurbs[7], 0.01094890510948905), "Nurbs2d P=2,N7");
  test.check(eq(N_Nurbs[8], 0.006204379562043796), "Nurbs2d P=2,N8");
  test.check(eq(std::accumulate(N_Nurbs.begin(), N_Nurbs.end(), 0.0), 1.0), "partition of unity in domain");

  xieta   = {0, 0.1};
  N_Nurbs = Dune::IGA::Nurbs<2>::basisFunctions(xieta, knots2, degree2, weightNet).directGetAll();
  test.check(eq(N_Nurbs[0], 0.8204545454545455), "Nurbs P=2,N0");  // check ansatzfunctions on boundaries
  test.check(eq(N_Nurbs[1], 0.0), "Nurbs P=2,N1");
  test.check(eq(N_Nurbs[2], 0.0), "Nurbs P=2,N2");
  test.check(eq(N_Nurbs[3], 0.1727272727272728), "Nurbs P=2,N3");
  test.check(eq(N_Nurbs[4], 0.0), "Nurbs P=2,N4");
  test.check(eq(N_Nurbs[5], 0.0), "Nurbs P=2,N5");
  test.check(eq(N_Nurbs[6], 0.00681818181818182), "Nurbs P=2,N6");
  test.check(eq(N_Nurbs[7], 0.0), "Nurbs P=2,N7");
  test.check(eq(N_Nurbs[8], 0.0), "Nurbs P=2,N8");

  test.check(eq(std::accumulate(N_Nurbs.begin(), N_Nurbs.end(), 0.0), 1.0), "partition of unity on boundary");
  xieta         = {0, 0.1};
  auto dN_Nurbs = Dune::IGA::Nurbs<2>::basisFunctionDerivatives(xieta, knots2, degree2, weightNet, 5);

  test.check(eq(dN_Nurbs.get({0, 0}).get({0, 0}), 0.8204545454545455),
             "Nurbs P=2,N0");  // check ansatzfunctions on boundaries
  test.check(eq(dN_Nurbs.get({0, 0}).get({1, 0}), 0.0), "Nurbs P=2,N1");
  test.check(eq(dN_Nurbs.get({0, 0}).get({2, 0}), 0.0), "Nurbs P=2,N2");
  test.check(eq(dN_Nurbs.get({0, 0}).get({0, 1}), 0.1727272727272728), "Nurbs P=2,N3");
  test.check(eq(dN_Nurbs.get({0, 0}).get({1, 1}), 0.0), "Nurbs P=2,N4");
  test.check(eq(dN_Nurbs.get({0, 0}).get({2, 1}), 0.0), "Nurbs P=2,N5");
  test.check(eq(dN_Nurbs.get({0, 0}).get({0, 2}), 0.00681818181818182), "Nurbs P=2,N6");
  test.check(eq(dN_Nurbs.get({0, 0}).get({1, 2}), 0.0), "Nurbs P=2,N7");
  test.check(eq(dN_Nurbs.get({0, 0}).get({2, 2}), 0.0), "Nurbs P=2,N8");

  // first derivative in u dir
  test.check(eq(dN_Nurbs.get({1, 0}).get({0, 0}), -24.166115702479335),
             "Nurbs P=2,N0");  // check ansatzfunctions on boundaries
  test.check(eq(dN_Nurbs.get({1, 0}).get({1, 0}), 26.25454545454545), "Nurbs P=2,N1");
  test.check(eq(dN_Nurbs.get({1, 0}).get({2, 0}), 0.0), "Nurbs P=2,N2");
  test.check(eq(dN_Nurbs.get({1, 0}).get({0, 1}), -5.0876033057851231), "Nurbs P=2,N3");
  test.check(eq(dN_Nurbs.get({1, 0}).get({1, 1}), 3.1090909090909089), "Nurbs P=2,N4");
  test.check(eq(dN_Nurbs.get({1, 0}).get({2, 1}), 0.0), "Nurbs P=2,N5");
  test.check(eq(dN_Nurbs.get({1, 0}).get({0, 2}), -0.20082644628099172), "Nurbs P=2,N6");
  test.check(eq(dN_Nurbs.get({1, 0}).get({1, 2}), 0.090909090909090925), "Nurbs P=2,N7");
  test.check(eq(dN_Nurbs.get({1, 0}).get({2, 2}), 0.0), "Nurbs P=2,N8");

  // second derivative in u dir
  test.check(eq(dN_Nurbs.get({2, 0}).get({0, 0}), 1236.8386175807659),
             "Nurbs P=2,N0");  // check ansatzfunctions on boundaries
  test.check(eq(dN_Nurbs.get({2, 0}).get({1, 0}), -1441.6132231404956), "Nurbs P=2,N1");
  test.check(eq(dN_Nurbs.get({2, 0}).get({2, 0}), 98.454545454545439), "Nurbs P=2,N2");
  test.check(eq(dN_Nurbs.get({2, 0}).get({0, 1}), 260.38707738542439), "Nurbs P=2,N3");
  test.check(eq(dN_Nurbs.get({2, 0}).get({1, 1}), -170.71735537190079), "Nurbs P=2,N4");
  test.check(eq(dN_Nurbs.get({2, 0}).get({2, 1}), 11.054545454545453), "Nurbs P=2,N5");
  test.check(eq(dN_Nurbs.get({2, 0}).get({0, 2}), 10.278437265214125), "Nurbs P=2,N6");
  test.check(eq(dN_Nurbs.get({2, 0}).get({1, 2}), -4.9917355371900829), "Nurbs P=2,N7");
  test.check(eq(dN_Nurbs.get({2, 0}).get({2, 2}), 0.30909090909090914), "Nurbs P=2,N8");

  // third derivative in u dir
  test.check(eq(dN_Nurbs.get({3, 0}).get({0, 0}), -94449.494433440283),
             "Nurbs P=2,N0");  // check ansatzfunctions on boundaries
  test.check(eq(dN_Nurbs.get({3, 0}).get({1, 0}), 110086.82794891056), "Nurbs P=2,N1");
  test.check(eq(dN_Nurbs.get({3, 0}).get({2, 0}), -7518.3471074380141), "Nurbs P=2,N2");
  test.check(eq(dN_Nurbs.get({3, 0}).get({0, 1}), -19884.104091250589), "Nurbs P=2,N3");
  test.check(eq(dN_Nurbs.get({3, 0}).get({1, 1}), 13036.598046581514), "Nurbs P=2,N4");
  test.check(eq(dN_Nurbs.get({3, 0}).get({2, 1}), -844.16528925619821), "Nurbs P=2,N5");
  test.check(eq(dN_Nurbs.get({3, 0}).get({0, 2}), -784.89884570726042), "Nurbs P=2,N6");
  test.check(eq(dN_Nurbs.get({3, 0}).get({1, 2}), 381.18707738542452), "Nurbs P=2,N7");
  test.check(eq(dN_Nurbs.get({3, 0}).get({2, 2}), -23.603305785123968), "Nurbs P=2,N8");

  // fourth derivative in u dir
  test.check(eq(dN_Nurbs.get({4, 0}).get({0, 0}), 9616675.7968593743),
             "Nurbs P=2,N0");  // check ansatzfunctions on boundaries
  test.check(eq(dN_Nurbs.get({4, 0}).get({1, 0}), -11208840.663889075), "Nurbs P=2,N1");
  test.check(eq(dN_Nurbs.get({4, 0}).get({2, 0}), 765504.432757325), "Nurbs P=2,N2");
  test.check(eq(dN_Nurbs.get({4, 0}).get({0, 1}), 2024563.3256546052), "Nurbs P=2,N3");
  test.check(eq(dN_Nurbs.get({4, 0}).get({1, 1}), -1327362.7101973905), "Nurbs P=2,N4");
  test.check(eq(dN_Nurbs.get({4, 0}).get({2, 1}), 85951.374906085635), "Nurbs P=2,N5");
  test.check(eq(dN_Nurbs.get({4, 0}).get({0, 2}), 79916.973381102871), "Nurbs P=2,N6");
  test.check(eq(dN_Nurbs.get({4, 0}).get({1, 2}), -38811.775151970498), "Nurbs P=2,N7");
  test.check(eq(dN_Nurbs.get({4, 0}).get({2, 2}), 2403.2456799398947), "Nurbs P=2,N8");

  // fifth derivative in u dir
  test.check(eq(dN_Nurbs.get({5, 0}).get({0, 0}), -1223940555.9639204),
             "Nurbs P=2,N0");  // check ansatzfunctions on boundaries
  test.check(eq(dN_Nurbs.get({5, 0}).get({1, 0}), 1426579720.8586094), "Nurbs P=2,N1");
  test.check(eq(dN_Nurbs.get({5, 0}).get({2, 0}), -97427836.896386802), "Nurbs P=2,N2");
  test.check(eq(dN_Nurbs.get({5, 0}).get({0, 1}), -257671695.99240425), "Nurbs P=2,N3");
  test.check(eq(dN_Nurbs.get({5, 0}).get({1, 1}), 168937072.20694059), "Nurbs P=2,N4");
  test.check(eq(dN_Nurbs.get({5, 0}).get({2, 1}), -10939265.897138171), "Nurbs P=2,N5");
  test.check(eq(dN_Nurbs.get({5, 0}).get({0, 2}), -10171251.15759491), "Nurbs P=2,N6");
  test.check(eq(dN_Nurbs.get({5, 0}).get({1, 2}), 4939680.4738871539), "Nurbs P=2,N7");
  test.check(eq(dN_Nurbs.get({5, 0}).get({2, 2}), -305867.63199235022), "Nurbs P=2,N8");

  // first derivative in v dir
  test.check(eq(dN_Nurbs.get({0, 1}).get({0, 0}), -1.609504132231405),
             "Nurbs P=2,N0");  // check ansatzfunctions on boundaries
  test.check(eq(dN_Nurbs.get({0, 1}).get({1, 0}), 0.0), "Nurbs P=2,N1");
  test.check(eq(dN_Nurbs.get({0, 1}).get({2, 0}), 0.0), "Nurbs P=2,N2");
  test.check(eq(dN_Nurbs.get({0, 1}).get({0, 1}), 1.479338842975206), "Nurbs P=2,N3");
  test.check(eq(dN_Nurbs.get({0, 1}).get({1, 1}), 0.0), "Nurbs P=2,N4");
  test.check(eq(dN_Nurbs.get({0, 1}).get({2, 1}), 0.0), "Nurbs P=2,N5");
  test.check(eq(dN_Nurbs.get({0, 1}).get({0, 2}), 0.1301652892561984), "Nurbs P=2,N6");
  test.check(eq(dN_Nurbs.get({0, 1}).get({1, 2}), 0.0), "Nurbs P=2,N7");
  test.check(eq(dN_Nurbs.get({0, 1}).get({2, 2}), 0.0), "Nurbs P=2,N8");

  // second derivative in v dir
  test.check(eq(dN_Nurbs.get({0, 2}).get({0, 0}), 3.380916604057099),
             "Nurbs P=2,N0");  // check ansatzfunctions on boundaries
  test.check(eq(dN_Nurbs.get({0, 2}).get({1, 0}), 0.0), "Nurbs P=2,N1");
  test.check(eq(dN_Nurbs.get({0, 2}).get({2, 0}), 0.0), "Nurbs P=2,N2");
  test.check(eq(dN_Nurbs.get({0, 2}).get({0, 1}), -4.507888805409466), "Nurbs P=2,N3");
  test.check(eq(dN_Nurbs.get({0, 2}).get({1, 1}), 0.0), "Nurbs P=2,N4");
  test.check(eq(dN_Nurbs.get({0, 2}).get({2, 1}), 0.0), "Nurbs P=2,N5");
  test.check(eq(dN_Nurbs.get({0, 2}).get({0, 2}), 1.126972201352367), "Nurbs P=2,N6");
  test.check(eq(dN_Nurbs.get({0, 2}).get({1, 2}), 0.0), "Nurbs P=2,N7");
  test.check(eq(dN_Nurbs.get({0, 2}).get({2, 2}), 0.0), "Nurbs P=2,N8");

  // third derivative in v dir
  test.check(eq(dN_Nurbs.get({0, 3}).get({0, 0}), -9.220681647428448),
             "Nurbs P=2,N0");  // check ansatzfunctions on boundaries
  test.check(eq(dN_Nurbs.get({0, 3}).get({1, 0}), 0.0), "Nurbs P=2,N1");
  test.check(eq(dN_Nurbs.get({0, 3}).get({2, 0}), 0.0), "Nurbs P=2,N2");
  test.check(eq(dN_Nurbs.get({0, 3}).get({0, 1}), 12.29424219657127), "Nurbs P=2,N3");
  test.check(eq(dN_Nurbs.get({0, 3}).get({1, 1}), 0.0), "Nurbs P=2,N4");
  test.check(eq(dN_Nurbs.get({0, 3}).get({2, 1}), 0.0), "Nurbs P=2,N5");
  test.check(eq(dN_Nurbs.get({0, 3}).get({0, 2}), -3.073560549142818), "Nurbs P=2,N6");
  test.check(eq(dN_Nurbs.get({0, 3}).get({1, 2}), 0.0), "Nurbs P=2,N7");
  test.check(eq(dN_Nurbs.get({0, 3}).get({2, 2}), 0.0), "Nurbs P=2,N8");

  // fourth derivative in v dir
  test.check(eq(dN_Nurbs.get({0, 4}).get({0, 0}), 33.52975144519434),
             "Nurbs P=2,N0");  // check ansatzfunctions on boundaries
  test.check(eq(dN_Nurbs.get({0, 4}).get({1, 0}), 0.0), "Nurbs P=2,N1");
  test.check(eq(dN_Nurbs.get({0, 4}).get({2, 0}), 0.0), "Nurbs P=2,N2");
  test.check(eq(dN_Nurbs.get({0, 4}).get({0, 1}), -44.70633526025914), "Nurbs P=2,N3");
  test.check(eq(dN_Nurbs.get({0, 4}).get({1, 1}), 0.0), "Nurbs P=2,N4");
  test.check(eq(dN_Nurbs.get({0, 4}).get({2, 1}), 0.0), "Nurbs P=2,N5");
  test.check(eq(dN_Nurbs.get({0, 4}).get({0, 2}), 11.17658381506479), "Nurbs P=2,N6");
  test.check(eq(dN_Nurbs.get({0, 4}).get({1, 2}), 0.0), "Nurbs P=2,N7");
  test.check(eq(dN_Nurbs.get({0, 4}).get({2, 2}), 0.0), "Nurbs P=2,N8");

  // fifth derivative in v dir
  test.check(eq(dN_Nurbs.get({0, 5}).get({0, 0}), -152.4079611145197),
             "Nurbs P=2,N0");  // check ansatzfunctions on boundaries
  test.check(eq(dN_Nurbs.get({0, 5}).get({1, 0}), 0.0), "Nurbs P=2,N1");
  test.check(eq(dN_Nurbs.get({0, 5}).get({2, 0}), 0.0), "Nurbs P=2,N2");
  test.check(eq(dN_Nurbs.get({0, 5}).get({0, 1}), 203.2106148193597), "Nurbs P=2,N3");
  test.check(eq(dN_Nurbs.get({0, 5}).get({1, 1}), 0.0), "Nurbs P=2,N4");
  test.check(eq(dN_Nurbs.get({0, 5}).get({2, 1}), 0.0), "Nurbs P=2,N5");
  test.check(eq(dN_Nurbs.get({0, 5}).get({0, 2}), -50.80265370483993), "Nurbs P=2,N6");
  test.check(eq(dN_Nurbs.get({0, 5}).get({1, 2}), 0.0), "Nurbs P=2,N7");
  test.check(eq(dN_Nurbs.get({0, 5}).get({2, 2}), 0.0), "Nurbs P=2,N8");

  // 1,1 mixed derivative
  test.check(eq(dN_Nurbs.get({1, 1}).get({0, 0}), 66.39293764087149),
             "Nurbs P=2,N0");  // check ansatzfunctions on boundaries
  test.check(eq(dN_Nurbs.get({1, 1}).get({1, 0}), -51.50413223140495), "Nurbs P=2,N1");
  test.check(eq(dN_Nurbs.get({1, 1}).get({2, 0}), 0.0), "Nurbs P=2,N2");
  test.check(eq(dN_Nurbs.get({1, 1}).get({0, 1}), -39.5762584522915), "Nurbs P=2,N3");
  test.check(eq(dN_Nurbs.get({1, 1}).get({1, 1}), 26.62809917355371), "Nurbs P=2,N4");
  test.check(eq(dN_Nurbs.get({1, 1}).get({2, 1}), 0.0), "Nurbs P=2,N5");
  test.check(eq(dN_Nurbs.get({1, 1}).get({0, 2}), -3.676183320811419), "Nurbs P=2,N6");
  test.check(eq(dN_Nurbs.get({1, 1}).get({1, 2}), 1.735537190082645), "Nurbs P=2,N7");
  test.check(eq(dN_Nurbs.get({1, 1}).get({2, 2}), 0.0), "Nurbs P=2,N8");

  // 2,1 mixed derivative
  test.check(eq(dN_Nurbs.get({2, 1}).get({0, 0}), -4511.311932245063),
             "Nurbs P=2,N0");  // check ansatzfunctions on boundaries
  test.check(eq(dN_Nurbs.get({2, 1}).get({1, 0}), 4043.131480090156), "Nurbs P=2,N1");
  test.check(eq(dN_Nurbs.get({2, 1}).get({2, 0}), -193.1404958677686), "Nurbs P=2,N2");
  test.check(eq(dN_Nurbs.get({2, 1}).get({0, 1}), 1791.166723584455), "Nurbs P=2,N3");
  test.check(eq(dN_Nurbs.get({2, 1}).get({1, 1}), -1318.232907588279), "Nurbs P=2,N4");
  test.check(eq(dN_Nurbs.get({2, 1}).get({2, 1}), 94.67768595041321), "Nurbs P=2,N5");
  test.check(eq(dN_Nurbs.get({2, 1}).get({0, 2}), 178.898026091114), "Nurbs P=2,N6");
  test.check(eq(dN_Nurbs.get({2, 1}).get({1, 2}), -91.08940646130728), "Nurbs P=2,N7");
  test.check(eq(dN_Nurbs.get({2, 1}).get({2, 2}), 5.900826446280992), "Nurbs P=2,N8");

  // std::cout<<std::setprecision(16)<<std::endl;
  // std::cout<<dN_Nurbs.get({3,1}).get({0,0})<<std::endl;
  // std::cout<<dN_Nurbs.get({3,1}).get({1,0})<<std::endl;
  // std::cout<<dN_Nurbs.get({3,1}).get({2,0})<<std::endl;
  // std::cout<<dN_Nurbs.get({3,1}).get({0,1})<<std::endl;
  // std::cout<<dN_Nurbs.get({3,1}).get({1,1})<<std::endl;
  // std::cout<<dN_Nurbs.get({3,1}).get({2,1})<<std::endl;
  // std::cout<<dN_Nurbs.get({3,1}).get({0,2})<<std::endl;
  // std::cout<<dN_Nurbs.get({3,1}).get({1,2})<<std::endl;
  // std::cout<<dN_Nurbs.get({3,1}).get({2,2})<<std::endl;

  // 3,1 mixed derivative
  test.check(eq(dN_Nurbs.get({3, 1}).get({0, 0}), 430363.3606745686),
             "Nurbs P=2,N0");  // check ansatzfunctions on boundaries
  test.check(eq(dN_Nurbs.get({3, 1}).get({1, 0}), -408827.1566149851), "Nurbs P=2,N1");
  test.check(eq(dN_Nurbs.get({3, 1}).get({2, 0}), 21583.77160030052), "Nurbs P=2,N2");
  test.check(eq(dN_Nurbs.get({3, 1}).get({0, 1}), -118703.546081676), "Nurbs P=2,N3");
  test.check(eq(dN_Nurbs.get({3, 1}).get({1, 1}), 88813.60562803084), "Nurbs P=2,N4");
  test.check(eq(dN_Nurbs.get({3, 1}).get({2, 1}), -6462.50939143501), "Nurbs P=2,N5");
  test.check(eq(dN_Nurbs.get({3, 1}).get({0, 2}), -12947.75940540574), "Nurbs P=2,N6");
  test.check(eq(dN_Nurbs.get({3, 1}).get({1, 2}), 6609.384604876715), "Nurbs P=2,N7");
  test.check(eq(dN_Nurbs.get({3, 1}).get({2, 2}), -429.1510142749812), "Nurbs P=2,N8");
  return test;
}

#include <dune/grid/test/checkindexset.hh>
#include <dune/iga/gridcapabilities.hh>
void gridCheck() {
  TestSuite test;

  auto circle = makeCircularArc();
  circle.controlPoints.directGet(0).p[1] += 3;
  const auto patch = makeSurfaceOfRevolution(circle, {2.0, 0, 0}, {0, 1, 0}, 360.0);

  IGA::NURBSGrid<2UL, 3UL> grid(patch);
  grid.globalRefine(2);
  //  grid.globalRefineInDirection(0,2);
  auto gridView = grid.leafGridView();

  //  Dune::checkIndexSet(grid,gridView_,std::cout);
}

template <typename C>
bool checkIfFinite(const C& v) {
  for (auto& vi : v) {
    if constexpr (std::is_arithmetic_v<std::remove_cvref_t<decltype(vi)>>) {
      if (not std::isfinite(vi)) return false;
    } else {
      for (auto& vii : vi)
        if constexpr (std::is_arithmetic_v<std::remove_cvref_t<decltype(vii)>>) {
          if (not std::isfinite(vii)) return false;
        } else {
          for (auto& viii : vii)
            if (not std::isfinite(viii)) return false;
        }
    }
  }
  return true;
};

template <typename C>
bool checkIfZero(const C& v) {
  for (auto& vi : v) {
    if constexpr (std::is_arithmetic_v<std::remove_cvref_t<decltype(vi)>>) {
      if (Dune::FloatCmp::ne(vi, 0.0)) return false;
    } else {
      for (auto& vii : vi)
        if constexpr (std::is_arithmetic_v<std::remove_cvref_t<decltype(vii)>>) {
          if (Dune::FloatCmp::ne(vii, 0.0)) return false;
        } else {
          for (auto& viii : vii)
            if (Dune::FloatCmp::ne(viii, 0.0)) return false;
        }
    }
  }
  return true;
};

auto testPlate() {
  constexpr int gridDim                = 2;
  constexpr auto dimworld              = 2;
  const std::array<int, gridDim> order = {2, 2};
  TestSuite t;

  const std::array<std::vector<double>, gridDim> knotSpans = {{{0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1}}};

  using ControlPoint = Dune::IGA::NURBSPatchData<gridDim, dimworld>::ControlPointType;

  const std::vector<std::vector<ControlPoint>> endPoints
      = {{{.p = {0, 0}, .w = 1}, {.p = {0.5, 0}, .w = 1}, {.p = {1, 0}, .w = 1}},
         {{.p = {0, 0.5}, .w = 1}, {.p = {0.5, 0.5}, .w = 1}, {.p = {1, 0.5}, .w = 1}},
         {{.p = {0, 1}, .w = 1}, {.p = {0.5, 1}, .w = 1}, {.p = {1, 1}, .w = 1}}};

  std::array<int, gridDim> dimsize = {(int)(endPoints.size()), (int)(endPoints[0].size())};

  auto controlNet = Dune::IGA::NURBSPatchData<gridDim, dimworld>::ControlPointNetType(dimsize, endPoints);
  using Grid      = Dune::IGA::NURBSGrid<gridDim, dimworld>;

  Dune::IGA::NURBSPatchData<gridDim, dimworld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = order;
  patchData.controlPoints = controlNet;
  auto grid               = std::make_shared<Grid>(patchData);
  grid->globalRefine(0);
  auto gridView = grid->leafGridView();
  Dune::GeometryChecker<decltype(grid)> geometryChecker;
  geometryChecker.checkGeometry(gridView);
  Dune::checkIndexSet(*grid, gridView, std::cout);

  checkEntityLifetime(gridView, gridView.size(0));

  for (auto&& elegeo : elements(gridView) | std::views::transform([](const auto& ele) { return ele.geometry(); }))
    checkJacobians(elegeo);

  checkIterators(gridView);
  auto basis     = Dune::Functions::BasisFactory::makeBasis(gridView, gridView.impl().getPreBasis());
  auto localView = basis.localView();
  std::vector<Dune::FieldVector<double, 1>> N;
  std::vector<Dune::FieldVector<double, 1>> ddNi02;
  std::vector<Dune::FieldVector<double, 1>> ddNi20;
  std::vector<Dune::FieldVector<double, 1>> ddNi11;
  std::vector<Dune::FieldMatrix<double, 1, 2>> dN;
  for (auto& element : elements(gridView)) {
    localView.bind(element);
    const auto& fe     = localView.tree().finiteElement();
    auto& localBasis   = fe.localBasis();
    auto referenceEle  = referenceElement(element);
    auto geo           = element.geometry();
    auto localGeometry = referenceEle.template geometry<0>(0);
    for (int i = 0; i < localGeometry.corners(); ++i) {
      auto local  = localGeometry.corner(i);
      auto global = geo.global(local);
      auto J      = geo.jacobianTransposed(local);
      auto J2     = geo.impl().secondDerivativeOfPosition(local);

      t.check(checkIfFinite(global));
      t.check(checkIfFinite(J));
      t.check(checkIfZero(J2));  // J2 should be zero since the grid is linearly parametrized

      localBasis.evaluateFunction(local, N);
      localBasis.evaluateJacobian(local, dN);
      localBasis.partial({2, 0}, local, ddNi20);
      localBasis.partial({1, 1}, local, ddNi11);
      localBasis.partial({0, 2}, local, ddNi02);

      t.check(checkIfFinite(N));
      t.check(checkIfFinite(dN));
      t.check(checkIfFinite(ddNi20));
      t.check(checkIfFinite(ddNi11));
      t.check(checkIfFinite(ddNi02));
    }
  }
  t.subTest(checkUniqueEdges(gridView));
  return t;
}



auto testPatchGeometryCurve() {
  TestSuite t;

  const auto dim      = 1;
  const auto dimworld = 2;

  const std::array<int, dim> order                     = {3};
  const std::array<std::vector<double>, dim> knotSpans = {{{0, 0, 0, 0, 1, 1, 1, 1}}};

  using ControlPoint = Dune::IGA::NURBSPatchData<dim, dimworld>::ControlPointType;

  const std::vector<ControlPoint> controlPoints
      = {{.p = {-4, -4}, .w = 1},
         {.p = {-3, 2.8}, .w = 2.5},
         {.p = {2, -4}, .w = 1},
         {.p = {4, 4}, .w = 1} };

  std::array<int, dim> dimsize = {static_cast<int>(controlPoints.size())};
  auto controlNet              = Dune::IGA::NURBSPatchData<dim, dimworld>::ControlPointNetType(dimsize, controlPoints);

  Dune::IGA::NURBSPatchData<dim, dimworld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = order;
  patchData.controlPoints = controlNet;

  // Make Geometry
  NURBSPatchGeometry<dim, dimworld> geometry(std::make_shared<Dune::IGA::NURBSPatchData<dim, dimworld>>(patchData));

  auto p0 = geometry.global({0.0});
  t.check(Dune::FloatCmp::eq(p0, {-4, -4}));

  auto p1 = geometry.global({0.5});
  t.check(Dune::FloatCmp::eq(p1, {-1.32, 0.72}));

  auto p2 = geometry.global({1});
  t.check(Dune::FloatCmp::eq(p2, {4, 4}));

  auto p3 = geometry.global({0.25});
  t.check(Dune::FloatCmp::eq(p3, {-2.7607655502392343, 0.4688995215311005}));

  // Test Operator ()
  auto p4 = geometry({0.4});
  t.check(Dune::FloatCmp::eq(p4, {-1.9854368932038828, 0.7669902912621357}));

  // Check derivative
  auto jc0 = geometry.jacobianTransposed({0});
  t.check(Dune::FloatCmp::eq(jc0[0], {7.5, 51}));

  auto jc1 = geometry.jacobianTransposed({0.5});
  t.check(Dune::FloatCmp::eq(jc1[0], {7.4496, -0.9216}));

  // Check local function
  auto u0 = geometry.local({-4, -4});
  t.check(Dune::FloatCmp::eq(u0, {0}));

  auto u1 = geometry.local({-1.32, 0.72});
  t.check(Dune::FloatCmp::eq(u1, {0.5}));

  // geomdl reports 13.230641820866644 for the length of the curve. The volume function approaches this value,
  // if you use a lot of gauß-points
  auto len = geometry.volume();

  // Check corners
  t.check(geometry.corners() == 2);
  std::array<FieldVector<double, 2>, 2> expectedCorners{{
      FieldVector<double, 2>{-4, -4},
      FieldVector<double, 2>{4, 4}
  }};
  for (int i = 0; i < 2; ++i)
    t.check(Dune::FloatCmp::eq(geometry.corner(i), expectedCorners[i]));


  return t;
}

auto testPatchGeometrySurface() {
  TestSuite t;

  const auto dim                   = 2;
  const auto dimworld              = 3;
  const std::array<int, dim> order = {2, 2};

  const std::array<std::vector<double>, dim> knotSpans = {{{0, 0, 0, 1, 1, 1}, {0, 0, 0, 1, 1, 1}}};

  using ControlPoint = Dune::IGA::NURBSPatchData<dim, dimworld>::ControlPointType;

  const std::vector<std::vector<ControlPoint>> controlPoints
      = {{{.p = {0, 0, 1}, .w = 1}, {.p = {1, 0, 1}, .w = 1}, {.p = {2, 0, 2}, .w = 1}},
         {{.p = {0, 1, 0}, .w = 1}, {.p = {1, 1, 0}, .w = 1}, {.p = {2, 1, 0}, .w = 1}},
         {{.p = {0, 2, 1}, .w = 1}, {.p = {1, 2, 2}, .w = 1}, {.p = {2, 2, 2}, .w = 1}}};

  std::array<int, dim> dimsize = {static_cast<int>(controlPoints.size()), static_cast<int>(controlPoints[0].size())};

  auto controlNet = Dune::IGA::NURBSPatchData<dim, dimworld>::ControlPointNetType(dimsize, controlPoints);

  Dune::IGA::NURBSPatchData<dim, dimworld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = order;
  patchData.controlPoints = controlNet;

  // Make Geometry
  NURBSPatchGeometry<dim, dimworld> geometry(std::make_shared<Dune::IGA::NURBSPatchData<dim, dimworld>>(patchData));

  auto p1 = geometry.global(FieldVector<double, 2>{0.5, 0.5});
  t.check(Dune::FloatCmp::eq(p1, {1.0, 1.0, 0.75}));

  auto p2 = geometry.global(FieldVector<double, 2>{0, 0});
  t.check(Dune::FloatCmp::eq(p2, {0, 0, 1}));

  auto p3 = geometry.global(FieldVector<double, 2>{0, 1});
  t.check(Dune::FloatCmp::eq(p3, {2, 0, 2}));

  // Check derivative
  auto jc1 = geometry.jacobianTransposed({0.5, 0.5});
  t.check(Dune::FloatCmp::eq(jc1[0], {0.0, 2.0, 0.5}));
  t.check(Dune::FloatCmp::eq(jc1[1], {2.0, 0.0, 0.5}));

  // Check local function
  auto u1 = geometry.local({1.0, 1.0, 0.75});
  t.check(Dune::FloatCmp::eq(u1, {0.5, 0.5}));

  auto u2 = geometry.local({2, 0, 2});
  t.check(Dune::FloatCmp::eq(u2, {0, 1}));

  // Check corners
  t.check(geometry.corners() == 4);

  std::array<FieldVector<double, 3>, 4> expectedCorners{{
      FieldVector<double, 3>{0, 0, 1},
      FieldVector<double, 3>{0, 2, 1},
      FieldVector<double, 3>{2, 0, 2},
      FieldVector<double, 3>{2, 2, 2}
  }};
  for (int i = 0; i < 4; ++i)
    t.check(Dune::FloatCmp::eq(geometry.corner(i), expectedCorners[i]));

  // Check domain
  t.check(Dune::FloatCmp::eq(geometry.domain()[0][0], 0.0));
  t.check(Dune::FloatCmp::eq(geometry.domain()[0][1], 1.0));

  // Check domain midpoint
  t.check(Dune::FloatCmp::eq(geometry.domainMidPoint()[0], 0.5));

  return t;
}


auto testIbraReader()
{
  TestSuite t;

  std::shared_ptr<NURBSGrid<2,2>> grid = IbraReader<2, 2>::read("auxiliaryFiles/element.ibra");

  // Check n_ele = 1, n_vert = 4
  t.check(grid->size(0) == 1);
  t.check(grid->size(2) == 4);

  grid->globalRefine(1);

  // Check n_ele = 4, n_vert = 9 after refinement
  t.check(grid->size(0) == 4);
  t.check(grid->size(2) == 9);

  // Test degree (maybe test degree elevate)
  t.check(grid->leafGridView().impl().getPatchData().degree[0] == 1);
  t.check(grid->leafGridView().impl().getPatchData().degree[1] == 1);

  // Enumerate elements and check position of centers
  auto gV = grid->leafGridView();
  std::vector<FieldVector<double, 2>> expectedElementCenters{{0.25, 0.25}, {0.75, 0.25}, {0.25, 0.75}, {0.75, 0.75}};
  const auto& indexSet = grid->leafGridView().indexSet();
  for (auto& ele : elements(gV))
    t.check(ele.geometry().center() == expectedElementCenters[indexSet.index(ele)]);

  // Test shell structure, no trim functionality right now
  std::shared_ptr<NURBSGrid<2,3>> grid3D = IbraReader<2, 3>::read("auxiliaryFiles/schale.ibra");
  grid3D->globalRefine(2);
  VTKWriter<NURBSGrid<2,3>::GridView> vtkWriter(grid3D->leafGridView());
  vtkWriter.write("grid");

  return t;
}

std::array<int, 3> getAmountOfTrimFlags(const auto& gridView) {
  int trimmedCounter {0};
  int emptyCounter {0};
  int fullCounter {0};
  for (auto& ele : elements(gridView))
    switch (ele.impl().getTrimFlag()) {
      case ElementTrimFlag::full:
        fullCounter++;
        break;
      case ElementTrimFlag::empty:
        emptyCounter++;
        break;
      case ElementTrimFlag::trimmed:
        trimmedCounter++;
        break;
    }

  return {trimmedCounter, emptyCounter, fullCounter};
}

auto testTrimImpactWithRefinement() {
  TestSuite t;

  // O refinement, 1 trimmed
  std::shared_ptr<NURBSGrid<2,2>> grid = IbraReader<2, 2>::read("auxiliaryFiles/element_trim.ibra");

  Plot::plotParametricGridAndPhysicalGrid(grid, "0");
  Plot::plotEveryReconstructedGrid(grid, "0");

  // 1 refinement: 3 trimmed, 0 empty, 1 full
  grid->globalRefine(1);

  Plot::plotParametricGridAndPhysicalGrid(grid, "1");
  Plot::plotEveryReconstructedGrid(grid, "1");

//  // 2 refinement: 6 trimmed, 2 empty, 8 full
//  grid->globalRefine(1);
//
//  Plot::plotParametricGridAndPhysicalGrid(grid, "2");
//  Plot::plotEveryReconstructedGrid(grid, "2");



  return t;
}


auto testMultiParametrisation() {
  TestSuite t;

  // 0 refinement, 1 trimmed
  std::shared_ptr<NURBSGrid<2,2>> grid = IbraReader<2, 2>::read("auxiliaryFiles/element_trim_2p.ibra");
  auto trimFlagCounter1 = getAmountOfTrimFlags(grid->leafGridView());
  t.check(trimFlagCounter1[0] == 1);

  Plot::plotParametricGridAndPhysicalGrid(grid, "0");
  Plot::plotEveryReconstructedGrid(grid, "0");

  // 1 refinement 3 trimmed, 1 full
  grid->globalRefine(1);

  auto trimFlagCounter2 = getAmountOfTrimFlags(grid->leafGridView());
  t.check(trimFlagCounter2[0] == 3);
  t.check(trimFlagCounter2[1] == 0);
  t.check(trimFlagCounter2[2] == 1);

  Plot::plotParametricGridAndPhysicalGrid(grid, "1");
  Plot::plotEveryReconstructedGrid(grid, "1");

  return t;
}

auto testNURBSSurfaceTrim() {
  TestSuite t;

  std::shared_ptr<NURBSGrid<2,2>> grid = IbraReader<2, 2>::read("auxiliaryFiles/nurbs_1.ibra");

  grid->globalRefine(4);
  t.check(grid->size(0) == 256);

  // 2nd surface trimmed circle
  std::shared_ptr<NURBSGrid<2,2>> grid2 = IbraReader<2, 2>::read("auxiliaryFiles/circle_trim.ibra");

  Plot::plotParametricGridAndPhysicalGrid(grid2, "0");
  Plot::plotEveryReconstructedGrid(grid2, "0");

  grid2->globalRefine(1);
  Plot::plotParametricGridAndPhysicalGrid(grid2, "1");
  Plot::plotEveryReconstructedGrid(grid2, "1");

  return t;
}

auto testPipeGeometry() {
  TestSuite t;

  std::shared_ptr<NURBSGrid<2,2>> grid = IbraReader<2, 2>::read("auxiliaryFiles/pipe_trim.ibra");
  grid->globalRefine(1);
  Plot::plotParametricGridAndPhysicalGrid(grid, "_pipe");
  Plot::plotEveryReconstructedGrid(grid, "_pipe");

  return t;
}


auto furtherExamples() {
  TestSuite t;

  std::shared_ptr<NURBSGrid<2,2>> grid = IbraReader<2, 2>::read("auxiliaryFiles/element_trim4.ibra");
  grid->globalRefine(1);

  Plot::plotParametricGridAndPhysicalGrid(grid, "_ex1");
  Plot::plotEveryReconstructedGrid(grid, "_ex1");
  Plot::saveEveryReconstructedGrid(grid, "_ex1");

  std::shared_ptr<NURBSGrid<2,2>> grid2 = IbraReader<2, 2>::read("auxiliaryFiles/element_trim_Xa.ibra");
  grid2->globalRefine(1);

  Plot::plotParametricGridAndPhysicalGrid(grid2, "_ex2");
  Plot::plotEveryReconstructedGrid(grid2, "_ex2");
  Plot::saveEveryReconstructedGrid(grid2, "_ex2");

  return t;
}

auto testHoleGeometry() {
  TestSuite t;

  std::shared_ptr<NURBSGrid<2,2>> grid = IbraReader<2, 2>::read("auxiliaryFiles/element_hole_circle.ibra");
  grid->globalRefine(3);

  Plot::plotParametricGridAndPhysicalGrid(grid, "_hole");
  Dune::printGrid(*grid, MPIHelper::instance(), "plot_hole/grid");
  Plot::plotEveryReconstructedGrid(grid, "_hole");

  VTKWriter vtkWriter(grid->leafGridView());
  vtkWriter.write("plot_hole/grid");

  return t;
}

auto testIntegrationPoints() {
  TestSuite t;

  // O refinement, 1 trimmed
  std::shared_ptr<NURBSGrid<2,2>> grid = IbraReader<2, 2>::read("auxiliaryFiles/element_trim.ibra");
  double area = 0;
  int order = 1;

  std::vector<Dune::QuadraturePoint<double, 2>> ipVec;
  for(auto& ele : elements(grid->leafGridView())) {
    ele.impl().getIntegrationPoints(ipVec, order);
    auto geo = ele.geometry();
    for (auto& ip : ipVec)
      area += geo.integrationElement(ip.position()) * ip.weight();
  }
  t.check(Dune::FloatCmp::eq(area, 0.737416, 1e-4));

  return t;
}

auto testGeometry() {
  TestSuite t;

  // O refinement, 1 trimmed
  std::shared_ptr<NURBSGrid<2,2>> grid = IbraReader<2, 2>::read("auxiliaryFiles/element_trim.ibra");
  for (auto& ele : elements(grid->leafGridView()))
    t.check(Dune::FloatCmp::eq(ele.geometry().volume(), 0.737416, 1e-4));

  return t;
}

auto testMapsInTrimmedPatch() {
  TestSuite t;

  // O refinement, 1 trimmed
  std::shared_ptr<NURBSGrid<2,2>> grid = IbraReader<2, 2>::read("auxiliaryFiles/element_trim.ibra");
  auto& patch = grid->getPatch();

  t.check(patch.getDirectIndex<0>(0) == 0);
  t.check(patch.getRealIndex<0>(0) == 0);

  // 1 refinement: 3 trimmed, 0 empty, 1 full
  grid->globalRefine(1);

  auto& patch_1_1             = grid->getPatch();
  auto [full, trimmed, empty] = patch_1_1.getAmountOfElementTrimTypes();
  t.check(full == 1);
  t.check(empty == 0);
  t.check(trimmed == 3);

  // As n_f + n_t = n, there has to be a 1 to 1 mapping of the indices
  for (int i = 0; i < 4; ++i) {
    t.check(patch_1_1.getDirectIndex<0>(i) == i);
    t.check(patch_1_1.getRealIndex<0>(i) == i);
  }

  // Load next example Grid
  std::shared_ptr<NURBSGrid<2,2>> grid2 = IbraReader<2, 2>::read("auxiliaryFiles/element_trim_Xb.ibra");
  grid2->globalRefine(1);
  auto& patch_2_1 = grid2->getPatch();

  // After one refinement we should have 3 trimmed one empty element
  auto [full2, trimmed2, empty2] = patch_2_1.getAmountOfElementTrimTypes();

  t.check(full2 == 0);
  t.check(empty2 == 1);
  t.check(trimmed2 == 3);

  // The second element is the first element with a direct Index and so forth
  for (int i = 1; i < 4; ++i) {
    t.check(patch_2_1.getRealIndex<0>(i) == i - 1);
    t.check(patch_2_1.getDirectIndex<0>(i - 1) == i);
  }

  return t;
}
#endif

#ifdef TEST_ALL

auto testEntityFunctionality2() {
  TestSuite t;

  std::shared_ptr<NURBSGrid<2,2>> grid = IbraReader<2, 2>::read("auxiliaryFiles/element_trim.ibra");
  grid->globalRefine(1);

  int counter = 0;
  for (auto& ele : elements(grid->leafGridView())) {
    counter++;
    t.check(ele.impl().getTrimFlag() != ElementTrimFlag::empty);
  }
  t.check(counter == 4);
  t.check(grid->size(0) == 4);

  std::vector<Dune::FieldVector<double, 2>> expectedElementCenters {{0.25, 0.25}, {0.75, 0.25}, {0.25, 0.75}, {0.75, 0.75}};
  std::vector<Dune::FieldVector<double, 2>> expectedEdgeCenters { {0.5, 0.25}, {1, 0.25}, {0.75, 0}, {0.75, 0.5}};
  std::vector<Dune::FieldVector<double, 2>> expectedElementCorners { {0.5, 0}, {1, 0}, {0.5, 0.5}, {1, 0.5}};

  Dune::FieldVector<double, 2> shift {1.7573970089500115, 0.0};
  std::ranges::for_each(expectedElementCenters, [&shift](auto& c) {c += shift;});
  std::ranges::for_each(expectedEdgeCenters, [&shift](auto& c) {c += shift;});
  std::ranges::for_each(expectedElementCorners, [&shift](auto& c) {c += shift;});

  for (int i = 0; auto& ele : elements(grid->leafGridView())) {
    t.check(Dune::FloatCmp::eq(expectedElementCenters.at(i), ele.geometry().center()));

    if (i == 1) {
      for (int j = 0; auto& intersection : intersections(grid->leafGridView(), ele)) {
        t.check(Dune::FloatCmp::eq(expectedEdgeCenters.at(j), intersection.geometry().center()));
        ++j;
      }
      for (int j = 0; j < 4; ++j)
        t.check(Dune::FloatCmp::eq(expectedElementCorners.at(j), ele.subEntity<2>(j).geometry().center()));
    }
    ++i;
  }

  // Test Dune Stuff
//  Dune::GeometryChecker<typename decltype(grid)::element_type> geometryChecker;
//  geometryChecker.checkGeometry(grid->leafGridView());
//  Dune::checkIndexSet(*grid, grid->leafGridView(), std::cout);
//  gridcheck(*grid);
//
//  Dune::printGrid(*grid, Dune::MPIHelper::instance());

  // Grid 2
  std::shared_ptr<NURBSGrid<2,2>> grid3 = IbraReader<2, 2>::read("auxiliaryFiles/element_hole_circle.ibra");
  grid3->globalRefine(2);

  // Test Dune Stuff
  Dune::GeometryChecker<typename decltype(grid3)::element_type> geometryChecker2;
  geometryChecker2.checkGeometry(grid3->leafGridView());
  Dune::checkIndexSet(*grid3, grid3->leafGridView(), std::cout);
  gridcheck(*grid3);

  return t;
}



auto testEntityFunctionality() {
  TestSuite t;

  std::shared_ptr<NURBSGrid<2,2>> grid = IbraReader<2, 2>::read("auxiliaryFiles/Element_trim_Xb.ibra");
  grid->globalRefine(1);

    int counter = 0;
    for (auto& ele : elements(grid->leafGridView())) {
      counter++;
      t.check(ele.impl().getTrimFlag() != ElementTrimFlag::empty);
    }
    t.check(counter == 3);
    t.check(grid->size(0) == 3);

    // Reference geometry
    Plot::plotGridView(grid->leafGridView(), "plot/test1");

    // Test element geometry (unfortunately the whole thing is shifted to the right and the top
    Dune::FieldVector<double, 2> shift {1.6142394486009448, 2.632929598653166};

    std::vector<Dune::FieldVector<double, 2>> expectedElementCenters {{0.75, 0.25}, {0.25, 0.75}, {0.75, 0.75}};
    std::vector<Dune::FieldVector<double, 2>> expectedEdgeCenters { {0.5, 0.25}, {1, 0.25}, {0.75, 0}, {0.75, 0.5}};
    std::vector<Dune::FieldVector<double, 2>> expectedElementCorners { {0.5, 0}, {1, 0}, {0.5, 0.5}, {1, 0.5}};

    std::ranges::for_each(expectedElementCenters, [&shift](auto& c) {c += shift;});
    std::ranges::for_each(expectedEdgeCenters, [&shift](auto& c) {c += shift;});
    std::ranges::for_each(expectedElementCorners, [&shift](auto& c) {c += shift;});

    for (int i = 0; auto& ele : elements(grid->leafGridView())) {
      t.check(Dune::FloatCmp::eq(expectedElementCenters.at(i), ele.geometry().center()));

      if (i == 0) {
        for (int j = 0; auto& intersection : intersections(grid->leafGridView(), ele)) {
          t.check(Dune::FloatCmp::eq(expectedEdgeCenters.at(j), intersection.geometry().center()));
          ++j;
        }
        for (int j = 0; j < 4; ++j)
          t.check(Dune::FloatCmp::eq(expectedElementCorners.at(j), ele.subEntity<2>(j).geometry().center()));
      }
      ++i;
    }

    // TODO Write tests for subentities


    // Test Size functions
    auto gV = grid->leafGridView();

    t.check(gV.size(0) == 3);
    t.check(gV.size(Dune::GeometryTypes::none(2)) == 3);
    t.check(gV.size(Dune::GeometryTypes::cube(2)) == 0);

//   Test Dune Stuff
    Dune::GeometryChecker<typename decltype(grid)::element_type> geometryChecker;
    geometryChecker.checkGeometry(grid->leafGridView());
    Dune::checkIndexSet(*grid, grid->leafGridView(), std::cout);

//   If this yields the correct boundaries, then we are happy
    Dune::printGrid(*grid, Dune::MPIHelper::instance());

  return t;
}

auto testDataCollectorAndVtkWriter() {
  TestSuite t;

  std::shared_ptr<NURBSGrid<2,2>> grid = IbraReader<2, 2>::read("auxiliaryFiles/element_hole_circle.ibra");
  grid->globalRefine(3);

const auto gv = grid->leafGridView();
  Dune::Vtk::DiscontinuousIgaDataCollector dataCollector1(gv);

  Dune::VtkUnstructuredGridWriter writer2(dataCollector1, Vtk::FormatTypes::ASCII);
  auto lambdaf= [](auto x){ return Dune::FieldVector<double,2>({std::sin(x[0]),std::cos(3*x[0])+std::sin(4*x[1])});};
  auto lambaGV = Dune::Functions::makeAnalyticGridViewFunction(lambdaf, gv);

  writer2.addPointData(lambaGV, Dune::VTK::FieldInfo("displacement", Dune::VTK::FieldInfo::Type::vector, 2));
  writer2.write("TestFile");

  return t;
}


auto checkDuneGeometryAndGrid() {
  TestSuite t;

  std::shared_ptr<NURBSGrid<2,2>> grid1 = IbraReader<2, 2>::read("auxiliaryFiles/element_trim.ibra", false);
  grid1->globalRefine(1);
  gridcheck(*grid1);
  printGrid(*grid1, MPIHelper::instance());


//  std::shared_ptr<NURBSGrid<2,2>> grid = IbraReader<2, 2>::read("auxiliaryFiles/element_trim_Xb.ibra", true);
//  grid->globalRefine(1);
//
////  Dune::GeometryChecker<typename decltype(grid)::element_type> geometryChecker;
////  geometryChecker.checkGeometry(grid->leafGridView());
////  Dune::checkIndexSet(*grid, grid->leafGridView(), std::cout);
//
//  gridcheck(*grid);


  return t;
}
#endif

auto testTrimFunctionality() {
  TestSuite t;

  constexpr int gridDim {2};
  constexpr int dimWorld {2};
  const std::array<std::vector<double>, gridDim> knotSpans{{{0, 0, 1, 1}, {0, 0, 1, 1}}};

  using ControlPoint = Dune::IGA::NURBSPatchData<gridDim, dimWorld>::ControlPointType;
  const double Lx    = 1;
  const double Ly    = 1;
  const std::vector<std::vector<ControlPoint>> controlPoints{{{.p = {0, 0}, .w = 1}, {.p = {0, Ly}, .w = 1}},
                                                             {{.p = {Lx, 0}, .w = 1}, {.p = {Lx, Ly}, .w = 1}}};

  std::array<int, gridDim> dimSize{{2, 2}};

  auto controlNet{Dune::IGA::NURBSPatchData<gridDim, dimWorld>::ControlPointNetType(dimSize, controlPoints)};

  using Grid = Dune::IGA::NURBSGrid<gridDim, dimWorld>;
  Dune::IGA::NURBSPatchData<gridDim, dimWorld> patchData;

  patchData.knotSpans     = knotSpans;
  patchData.degree        = {1, 1};
  patchData.controlPoints = controlNet;

  Grid grid{patchData};

  Clipper2Lib::PathsD clip;
  clip.push_back(Clipper2Lib::MakePathD("0,0, 0,0.5, 0.4,0.6, 0.5,0.65, 0.7, 0.4, 0.9, 0.45, 1,0.5, 1,0"));

  using namespace Dune::IGA::Impl::Trim;
  ClippingResult res;
  for (auto& ele : elements(grid.leafGridView())) {
    res = clipElement(ele, clip).second.value();
  }

  return t;
}


auto testTrimmedElementGrid() {
  TestSuite t;

  std::shared_ptr<NURBSGrid<2,2>> grid = IbraReader<2, 2>::read("auxiliaryFiles/element_trim.ibra");

  for (int i = 0; i < 3; ++i) {
    grid->globalRefine(1);
    for (int j = 0; const auto& ele : elements(grid->leafGridView())) {
      auto repr = grid->getPatch().getTrimmedElementRepresentation(j++);
      if (repr->isTrimmed()) {
        bool hasOverlap = repr->checkGridForOverlappingElements();
        t.check(!hasOverlap, "Grid Overlap at ref: " + std::to_string(i+1) + ", ele: " + std::to_string(j-1));
        if (hasOverlap) {
          auto gV = repr->gridView();
          Plot::plotGridView(gV, "overlapCheck/grid_" + std::to_string(i+1) + "_" + std::to_string(j-1));
        }
      }
    }
  }

  return t;
}


auto testNurbsBasis2() {
  TestSuite t;

//  std::shared_ptr<NURBSGrid<2,2>> grid = IbraReader<2, 2>::read("auxiliaryFiles/element_trim_Xb.ibra", true);
//  grid->globalRefine(2);
//
//  auto gridView = grid->leafGridView();
//  using GridView = decltype(gridView);
//  Dune::Functions::NurbsBasis<GridView> basis(gridView, gridView.impl().getPatchData());
//
//  std::cout << "Grid 1\n";
//  t.subTest(checkBasis(basis, EnableContinuityCheck(), EnableContinuityCheck()));
//
//  // Degree 2

  std::shared_ptr<NURBSGrid<2,2>> grid2 = IbraReader<2, 2>::read("auxiliaryFiles/Element_trim_Xb.ibra", true, {1, 1});
  grid2->globalRefine(1);
  printGrid(*grid2, MPIHelper::instance());

  auto gridView2 = grid2->leafGridView();
  using GridView = decltype(gridView2);
  Dune::Functions::NurbsBasis<GridView> basis2(gridView2);

  std::cout << "Grid 2\n";
  t.subTest(checkBasis(basis2, EnableContinuityCheck(), EnableContinuityCheck()));

  // Hole geometry

//  std::shared_ptr<NURBSGrid<2,2>> grid3 = IbraReader<2, 2>::read("auxiliaryFiles/element_hole_circle.ibra", true, {1, 1});
//  grid3->globalRefine(3);
//
//  auto gridView3 = grid3->leafGridView();
//  Dune::Functions::NurbsBasis<decltype(gridView3)> basis3(gridView3);
//
//  std::cout << "Hole\n";
//  t.subTest(checkBasis(basis3, EnableContinuityCheck(), EnableContinuityCheck()));

  return t;
}


int main(int argc, char** argv) try {
  // Initialize MPI, if necessary
  MPIHelper::instance(argc, argv);
  TestSuite t;

  //t.subTest(testNurbsBasis2());
  //t.subTest(checkDuneGeometryAndGrid());

  t.subTest(testTrimFunctionality());
//  t.subTest(testTrimmedElementGrid());

  //t.subTest(testIntegrationPoints());
  //t.subTest(testGeometry());
  //t.subTest(testHoleGeometry());

  //t.subTest(testPatchGeometryCurve());
  //t.subTest(testPatchGeometrySurface());
//
//  t.subTest(testMapsInTrimmedPatch());
 t.subTest(testEntityFunctionality());
  t.subTest(testEntityFunctionality2());
//
//
//  //t.subTest(testIbraReader());
//  t.subTest(testTrimImpactWithRefinement());
  // t.subTest(testEntityFunctionality());
  t.subTest(testDataCollectorAndVtkWriter());

  //t.subTest(testMultiParametrisation());
  //t.subTest(testNURBSSurfaceTrim());

  //t.subTest(testPipeGeometry());
  //t.subTest(furtherExamples());

  //t.subTest(testNurbsBasis());

#if 0
  t.subTest(test3DGrid());
  t.subTest(testNURBSGridCurve());
  t.subTest(testPlate());
  testNurbsGridCylinder();
  t.subTest(testTorusGeometry());



  gridCheck();
  t.subTest(testBsplineBasisFunctions());

#endif
  t.report();

  return 0;
} catch (Dune::Exception& e) {
  std::cerr << "Dune reported error: " << e << std::endl;
}

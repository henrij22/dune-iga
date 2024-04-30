// SPDX-FileCopyrightText: 2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later
#pragma once

#include <mapbox/earcut.hpp>

#include <dune/geometry/quadraturerules.hh>
#include <dune/iga/geometrykernel/geohelper.hh>

namespace Dune::IGANEW::DefaultTrim {

template <typename GridImp>
struct SimplexIntegrationRuleGenerator
{
  using PatchElement = typename GridImp::Traits::template Codim<0>::Entity;

  static constexpr int dim = GridImp::dimension;

  struct Parameters
  {
    int maxBoundaryDivisions{5};
    double targetTolerance{1e-4}; // not used atm
  };

  static auto createIntegrationRule(const PatchElement& element, GridImp* grid, int quadratureOrder,
                                    const Parameters& parameters  = Parameters{},
                                    const QuadratureType::Enum qt = QuadratureType::GaussLegendre) {
    auto gv = grid->levelGridView(element.level());

    std::vector<Point> vertices{};
    std::vector<Element> elements{};

    for (const auto& intersection : intersections(gv, element)) {
      auto geometryInInside = intersection.geometryInInside();

      std::vector<Point> pointsT = splitBoundary(geometryInInside, parameters);
      vertices.insert(vertices.end(), pointsT.begin(), pointsT.end());
    }

    auto indices = triangulate(vertices);
    for (auto it = indices.begin(); it < indices.end(); it += 3)
      elements.emplace_back(GeometryTypes::triangle, std::vector<FieldVector<double, dim>>{
                                                         vertices[*it], vertices[*(it + 1)], vertices[*(it + 2)]});

    return makeQuadratureRule(elements, quadratureOrder, qt);
  }

private:
  using Point   = FieldVector<double, dim>;
  using Element = MultiLinearGeometry<double, dim, dim>;
  using Index   = std::uint64_t;

  static std::vector<Point> splitBoundary(const auto& localGeometry, const Parameters& parameters) {
    std::vector<Point> points;
    points.reserve(parameters.maxBoundaryDivisions);

    // todo at the moment only boundaryDivisions
    for (auto local : Utilities::linspace(0.0, 1.0, parameters.maxBoundaryDivisions))
      points.push_back(localGeometry.global({local}));

    return points;
  }
  static auto triangulate(const std::vector<Point>& points) -> std::vector<Index> {
    std::vector<std::vector<Point>> polygonInput;
    polygonInput.push_back(points);

    return mapbox::earcut<Index>(polygonInput);
  }
  static auto makeQuadratureRule(const std::vector<Element>& elements, int quadratureOrder,
                                 const QuadratureType::Enum qt) {
    assert(qt == QuadratureType::GaussLegendre);

    Dune::QuadratureRule<double, dim> vector{};

    for (auto& subElement : elements) {
      const auto& rule = Dune::QuadratureRules<double, dim>::rule(subElement.type(), quadratureOrder, qt);
      for (auto ip : rule) {
        auto globalInSpan = subElement.global(ip.position());
        vector.emplace_back(globalInSpan, ip.weight() * subElement.integrationElement(ip.position()));
      }
    }
    return vector;
  }
};

} // namespace Dune::IGANEW::DefaultTrim

// Add support for Dune::FieldVector in Earcut
namespace mapbox::util {

template <typename T>
struct nth<0, Dune::FieldVector<T, 2>>
{
  inline static auto get(const Dune::FieldVector<T, 2>& t) {
    return t[0];
  };
};

template <typename T>
struct nth<1, Dune::FieldVector<T, 2>>
{
  inline static auto get(const Dune::FieldVector<T, 2>& t) {
    return t[1];
  };
};
} // namespace mapbox::util

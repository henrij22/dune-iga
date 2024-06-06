// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/iga/integrationrules/simplexintegrationrulegenerator.hh>

namespace Dune::IGA {
namespace DefaultTrim {
  template <typename GridImp>
  struct DefaultIntegrationRuleGenerator
  {
    using Generator = SimplexIntegrationRuleGenerator<GridImp>;
    static auto integrationRule() {
      return [](const auto& element, int order, QuadratureType::Enum qt = QuadratureType::GaussLegendre) {
        const auto parameters =
            typename Generator::Parameters{.boundaryDivisions = Preferences::getInstance().boundaryDivisions(),
                                           .targetAccuracy    = Preferences::getInstance().targetAccuracy()};
        return Generator::createIntegrationRule(element, order, parameters, qt);
      };
    }
  };
}



template <typename GridImp>
struct IntegrationRuleHolder
{
  using PatchElement       = typename GridImp::Traits::template Codim<0>::Entity;
  static constexpr int dim = GridImp::dimension;

  using FunctionType = std::function<QuadratureRule<double, dim>(const PatchElement&, int, QuadratureType::Enum)>;

  IntegrationRuleHolder()
      : generator_(DefaultTrim::DefaultIntegrationRuleGenerator<GridImp>::integrationRule()) {}

  void integrationRule(FunctionType generator) {
    generator_ = generator;
  }

  FunctionType integrationRule() const {
    return generator_;
  }

private:
  FunctionType generator_;
};
} // namespace Dune::IGA
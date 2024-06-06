// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "integrationrules/simplexintegrationrulegenerator.hh"

#include <mutex>

namespace Dune::IGA::DefaultTrim {
class Preferences
{
public:
  // Get the singleton instance
  static Preferences& getInstance() {
    static Preferences instance; // Guaranteed to be destroyed and initialized on first use
    return instance;
  }

  // Disable copy constructor and assignment operator
  Preferences(const Preferences&)            = delete;
  Preferences& operator=(const Preferences&) = delete;

  int boundaryDivisions() const {
    return boundaryDivisions_;
  }

  void boundaryDivisions(int _boundaryDivisions) {
    std::lock_guard lock(mtx);
    this->boundaryDivisions_ = _boundaryDivisions;
  }

  double targetAccuracy() const {
    return targetAccuracy_;
  }
  void targetAccuracy(double _targetAccuracy) {
    std::lock_guard lock(mtx);
    this->targetAccuracy_ = _targetAccuracy;
  }

  bool reportTrimmedElementGeometryTypeAsNone() const {
    return reportTrimmedElementGeometryTypeAsNone_;
  }

  void reportTrimmedElementGeometryTypeAsNone(bool _reportTrimmedElementGeometryTypeAsNone) {
    std::lock_guard lock(mtx);
    this->reportTrimmedElementGeometryTypeAsNone_ = _reportTrimmedElementGeometryTypeAsNone;
  }

  bool reconstructTrimmedLocalGeometry() const {
    return reconstructTrimmedLocalGeometry_;
  }

  void reconstructTrimmedLocalGeometry(bool _reconstructTrimmedLocalGeometry) {
    std::lock_guard lock(mtx);
    this->reconstructTrimmedLocalGeometry_ = _reconstructTrimmedLocalGeometry;
  }

private:
  // Private constructor for Singleton pattern
  Preferences()
      : boundaryDivisions_(5),
        targetAccuracy_{1} {}

  std::mutex mtx; // Mutex for thread safety
  int boundaryDivisions_;
  double targetAccuracy_;
  bool reportTrimmedElementGeometryTypeAsNone_{true};
  bool reconstructTrimmedLocalGeometry_{true};
};

template <typename GridImp>
struct DefaultIntegrationRuleGenerator
{
  using Generator = SimplexIntegrationRuleGenerator<GridImp>;
  static auto integrationRule() {
    return [](const auto& element, int order,
                                QuadratureType::Enum qt =  QuadratureType::GaussLegendre) {
      const auto parameters = typename Generator::Parameters{
        .boundaryDivisions = Preferences::getInstance().boundaryDivisions(),
        .targetAccuracy    = Preferences::getInstance().targetAccuracy()};
      return Generator::createIntegrationRule(element, order, parameters, qt);
    };
  }
};

template <typename GridImp>
struct IntegrationRuleHolder
{
  using PatchElement       = typename GridImp::Traits::template Codim<0>::Entity;
  static constexpr int dim = GridImp::dimension;

  using FunctionType = std::function<QuadratureRule<double, dim>(const PatchElement&, int, QuadratureType::Enum)>;

  IntegrationRuleHolder()
      : generator_(DefaultIntegrationRuleGenerator<GridImp>::integrationRule()) {}

  void integrationRule(FunctionType generator) {
    generator_ = generator;
  }

  FunctionType integrationRule() const {
    return generator_;
  }

private:
  FunctionType generator_;
};
} // namespace Dune::IGA::DefaultTrim
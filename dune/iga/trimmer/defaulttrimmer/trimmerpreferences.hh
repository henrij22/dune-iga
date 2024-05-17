// SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <mutex>

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

  bool reportTrimmedIntersections() const {
    return reportTrimmedIntersections_;
  }

  void reportTrimmedIntersections(bool _reportTrimmedIntersections_) {
    std::lock_guard lock(mtx);
    this->reportTrimmedIntersections_ = _reportTrimmedIntersections_;
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
  bool reportTrimmedIntersections_{true};
};
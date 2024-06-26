// SPDX-FileCopyrightText: 2022-2024 The dune-iga developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "dune/python/iga/gridenums.hh"
#include <dune/iga/parameterspace/default/trimmerpreferences.hh>
#include <dune/python/pybind11/pybind11.h>

PYBIND11_MODULE(_iga, m) {
  pybind11::enum_<Dune::Python::IGA::Reader> reader(m, "reader");
  reader.value("json", Dune::Python::IGA::Reader::json);

  m.def(
      "registerParameterSpacePreferences",
      [](int boundaryDivisions = 5, double targetAccuracy = 1) {
        Dune::IGA::DefaultParameterSpace::Preferences::getInstance().targetAccuracy(targetAccuracy);
        Dune::IGA::DefaultParameterSpace::Preferences::getInstance().boundaryDivisions(boundaryDivisions);
      },
      pybind11::arg("boundaryDivisions") = 5, pybind11::arg("targetAccuracy") = 1.0);
}

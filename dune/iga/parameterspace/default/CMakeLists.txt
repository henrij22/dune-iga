# SPDX-FileCopyrightText: 2022-2024 The dune-iga developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

add_subdirectory(trimmingutils)

set(HEADERS
    createentities.hh
    createlevel.hh
    elementtrimdata.hh
    entitycontainer.hh
    entityinfo.hh
    idset.hh
    patchgridentityseed.hh
    patchgridhierarchiciterator.hh
    patchgridindexsets.hh
    patchgridintersectioniterator.hh
    patchgridintersections.hh
    patchgridleafiterator.hh
    patchgridleveliterator.hh
    patchtrimdata.hh
    referenceelement.hh
    trimelement.hh
    parameterspaceentity.hh
    trimmedlocalgeometry.hh
    parameterspace.hh
    trimmerpreferences.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/iga/parameterspace/default)

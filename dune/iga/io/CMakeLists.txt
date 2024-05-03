# SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

add_subdirectory(vtk)

# install headers

set(HEADERS
        griddrawer.hh
        ibrageometry.hh
        ibrajsonreader.hh
        ibrareader.hh
)

install(FILES ${HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/iga/io
)

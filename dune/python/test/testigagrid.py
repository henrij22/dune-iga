from dune.iga import (
    IGAGrid,
    ControlPointNet,
    ControlPoint,
    NurbsPatchData
)

from dune.iga import reader as readeriga
from dune.iga.basis import defaultGlobalBasis, Power, Lagrange, Nurbs
from dune.common import FieldVector
from dune.grid import gridFunction


if __name__ == "__main__":
    inputParameter = dict(
                file_path="../../iga/test/auxiliaryfiles/element.ibra",
                reader=readeriga.json,
                pre_knot_refine=(1, 1),
            )
    gridView = IGAGrid(inputParameter, dimgrid=2, dimworld=2)
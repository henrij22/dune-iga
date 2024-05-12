from dune.iga import IGAGrid, IGAGridType, ControlPointNet, ControlPoint, NurbsPatchData

from dune.iga import reader as readeriga

if __name__ == "__main__":

    inputParameter = dict(
        file_path="../../iga/test/auxiliaryfiles/element.ibra",
        reader=readeriga.json,
        pre_knot_refine=(1, 1),
    )
    gridView = IGAGrid(inputParameter, dimgrid=2, dimworld=2, gridType=IGAGridType.Identity)

    gridView = IGAGrid(inputParameter, dimgrid=2, dimworld=2, gridType=IGAGridType.Default)


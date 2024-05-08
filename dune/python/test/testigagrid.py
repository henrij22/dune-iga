from dune.iga import reader as readeriga
from dune.iga import IGAGrid

if __name__ == "__main__":
    inputParameter = dict(
                file_path="../../iga/test/auxiliaryfiles/element.ibra",
                reader=readeriga.json,
                post_knot_refine=(1, 1),
            )
    gridView = IGAGrid(inputParameter, dimgrid=2, dimworld=2)
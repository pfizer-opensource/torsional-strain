from floe.api import OEMolIStreamCube
from torsion.utils.process_sd_data import print_torsion


class PrintTorsion(OEMolIStreamCube):
    """Print the torsion of interest in the molecules streamed through this cube."""

    title = "Torsion Printer"
    classification = [["Input"]]
    tags = ["dataset", "read"]

    def __iter__(self):
        parent_iter = super(DisplaySDData, self).__iter__()
        for mol in parent_iter:
            print_torsion(mol)
            oechem.OESetSDData(mol, "number of atoms", str(mol.NumAtoms()))
            yield mol

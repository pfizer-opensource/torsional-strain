from floe.api import OEMolComputeCube

from openeye import oechem

class FilterBrI(OEMolComputeCube):
    """
        Remove compounds with Bromine or Iodine
    """

    title = "Remove compounds with Br or I"
    classification = [["Cheminformatics"]]

    def process(self, mol, port):
        """
            Search for any Bromine atom and remove it.
        """
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == oechem.OEElemNo_Br or atom.GetAtomicNum() == oechem.OEElemNo_I:
                self.failure.emit(mol)
        self.success.emit(mol)


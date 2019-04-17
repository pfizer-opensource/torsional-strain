import numpy as np
from openeye import oechem

def get_bond_lengths(mc_mol):
    bond_lengths = []

    if mc_mol.__class__ == oechem.OEGraphMol:
        for bond in mc_mol.GetBonds():
            atom_pos0 = np.array(mc_mol.GetCoords()[bond.GetBgnIdx()])
            atom_pos1 = np.array(mc_mol.GetCoords()[bond.GetEndIdx()])
            bond_lengths.append(np.linalg.norm(atom_pos1 - atom_pos0))
    elif mc_mol.__class__ == oechem.OEMol:
        for conf in mc_mol.GetConfs():
            conf_bond_lengths = []
            for bond in conf.GetBonds():
                atom_pos0 = np.array(conf.GetCoords()[bond.GetBgnIdx()])
                atom_pos1 = np.array(conf.GetCoords()[bond.GetEndIdx()])
                conf_bond_lengths.append(np.linalg.norm(atom_pos1 - atom_pos0))
            bond_lengths.append(conf_bond_lengths)
    else:
        raise ValueError("Unrecognized molecule type.")
    return np.array(bond_lengths)

def has_undesirable_elements(mol):
    '''
    returns True if molecule contains any element other than
    H, C, N, O, F, S, Cl, or P

    @param mol:
    @type mol: OEGraphMol
    @return: bool
    '''
    atomsHC = oechem.OEOrAtom(oechem.OEIsHydrogen(), oechem.OEIsCarbon())
    atomsNO = oechem.OEOrAtom(oechem.OEIsNitrogen(), oechem.OEIsOxygen())
    atomsFS = oechem.OEOrAtom(oechem.OEHasAtomicNum(9), oechem.OEIsSulfur())
    atomsHCNO = oechem.OEOrAtom(atomsHC, atomsNO)
    atomsHCNOFS = oechem.OEOrAtom(atomsHCNO, atomsFS)
    atomsHCNOFSCl = oechem.OEOrAtom(atomsHCNOFS, oechem.OEHasAtomicNum(17))
    atomsHCNOFSClP = oechem.OEOrAtom(atomsHCNOFSCl, oechem.OEIsPhosphorus())

    undesirable_atom = mol.GetAtom(oechem.OENotAtom(atomsHCNOFSClP))
    if undesirable_atom is not None:
        return True

    return False

def is_undesirable_molecule(mol):
    if has_undesirable_elements(mol):
        return True
    if oechem.OECount(mol, oechem.IsRotor()) == 0:
        return True

    return False

def get_modified_molecule_inchi(mol):
    title = mol.GetTitle().replace(',', '')
    inchi = oechem.OECreateInChIKey(mol)
    return inchi + title


from openeye import oechem


class hasDoubleBondO(oechem.OEUnaryAtomPred):
    def __call__(self, atom):
        for bond in atom.GetBonds():
            if bond.GetOrder() == 2 and bond.GetNbr(atom).IsOxygen():
                return True
        return False


def isAmideRotor(bond):
    if bond.GetOrder() != 1:
        return False
    atomB = bond.GetBgn()
    atomE = bond.GetEnd()
    pred = hasDoubleBondO()
    if atomB.IsCarbon() and atomE.IsNitrogen() and pred(atomB):
        return True
    if atomB.IsNitrogen() and atomE.IsCarbon() and pred(atomE):
        return True
    return False


def isMethylRotor(bond):
    if bond.GetOrder() != 1:
        return False
    atomB = bond.GetBgn()
    atomE = bond.GetEnd()

    if atomB.IsHydrogen() or atomE.IsHydrogen():
        return False

    def isMethylCarbon(atom):
        return (
            atom.GetAtomicNum() == oechem.OEElemNo_C
            and atom.GetHvyDegree() == 1
            and atom.GetTotalHCount() == 3
        )

    return isMethylCarbon(atomB) or isMethylCarbon(atomE)


def isEtherRotor(bond):
    if bond.GetOrder() != 1:
        return False
    atomB = bond.GetBgn()
    atomE = bond.GetEnd()

    isEtherOxygen = oechem.OEMatchAtom("[OX2][C,c]")
    return (atomB.IsCarbon() and isEtherOxygen(atomE)) or (
        atomE.IsCarbon() and isEtherOxygen(atomB)
    )


def isRotatableBond(bond):
    inRing = oechem.OEBondIsInRing()

    return (not inRing(bond)) and (
        isAmideRotor(bond) or isMethylRotor(bond) or isEtherRotor(bond)
    )


# Torsion Library
torsion_library = [
    "[C,N,c:1][NX3:2][C:3](=[O])[C,N,c,O:4] 0 180",  # amides are flipped cis and trans
    "[#1:1][NX3H:2][C:3](=[O])[C,N,c,O:4] 0",  # primary amides are NOT flipped
    "[*:1][C,c:2][OX2:3][*:4] 0 180",  # hydroxyls and ethers are rotated 180 degrees
    "[H:1][CH3:2]-!@[!#1:3][*:4] 0.1 180",  # methyls are rotated 180 degrees
    "[H:1][CH3:2]-!@[!#1:3]=[*:4] 0.1 180",
]


class distance_predicate(oechem.OEUnaryBondPred):
    def __init__(self, atom1_idx, atom2_idx):
        oechem.OEUnaryBondPred.__init__(self)
        self.atom1_idx = atom1_idx
        self.atom2_idx = atom2_idx

    def __call__(self, bond):
        atomB = bond.GetBgn()
        atomE = bond.GetEnd()
        mol = bond.GetParent()
        atom1 = mol.GetAtom(oechem.OEHasAtomIdx(self.atom1_idx))
        atom2 = mol.GetAtom(oechem.OEHasAtomIdx(self.atom2_idx))
        return (
            max(
                oechem.OEGetPathLength(atomB, atom1),
                oechem.OEGetPathLength(atomE, atom1),
                oechem.OEGetPathLength(atomB, atom2),
                oechem.OEGetPathLength(atomE, atom2),
            )
            <= 3
        )

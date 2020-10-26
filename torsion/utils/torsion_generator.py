from __future__ import print_function

"""
Script to extract torsion fragments from parent molecules.

Torsion fragments are extracted based on a set of rules

"""
import sys, os
import subprocess, tempfile

from openeye.oechem import *
from openeye.oemedchem import *

# CORINA_CMD = "corina -d wh -d preserve -t n"

# SCRATCH_DIR = tempfile.mkdtemp()

TORSION_ATOMS_FRAG_KEY = "TORSION_ATOMS_FRAGMENT"
NUM_ADJACENT_TORSIONS_KEY = "NUM_ADJACENT_TORSIONS"
ADJACENT_TORSION_ATOMS_KEY = "ADJACENT_TORSION_ATOMS"
FRAGMENT_TO_PARENT_ATOMS_KEY = "FRAGMENT_TO_PARENT_ATOMS"
TORSION_ATOMS_PARENT_MOL_KEY = "TORSION_ATOMS_ParentMol"


class TorsionGenerator:
    @staticmethod
    def IsOrtho(atom, torsion):
        """

        @param atom:
        @param torsion:
        @type atom: OEAtombase
        @type torsion: OETorsion
        @return: bool
        """
        numRingAtoms = 0
        for a, b in zip(
            OEShortestPath(atom, torsion.b), OEShortestPath(atom, torsion.c)
        ):
            if a.IsInRing():
                numRingAtoms += 1
            if b.IsInRing():
                numRingAtoms += 1

        return numRingAtoms <= 4

    @staticmethod
    def GetNbrs(atomSet):
        """
        Returns atoms connected to atomSet atoms,
        excluding those from atomSet

        @param atomSet: OEAtomBondSet
        @type atomSet: OEAtomBondSet
        @return: None
        """
        allNbrs = []
        for atom in atomSet.GetAtoms():
            nbrs = atom.GetAtoms()
            for nbr in nbrs:
                if not atomSet.HasAtom(nbr):
                    allNbrs.append(nbr)

        return allNbrs

    @staticmethod
    def GetSameRingAtoms(mol, atomSet):
        """

        @type mol: OEGraphMol
        @type atomSet: OEAtomBondSet
        @return list[OEAtombase]
        """
        OEFindRingAtomsAndBonds(mol)
        numRings, ringIdxPerAtom = OEDetermineRingSystems(mol)
        toKeepRings = {}
        for atom in atomSet.GetAtoms():
            if ringIdxPerAtom[atom.GetIdx()] > 0:
                toKeepRings[ringIdxPerAtom[atom.GetIdx()]] = True

        ringAtoms = []
        for i in range(0, len(ringIdxPerAtom)):
            if ringIdxPerAtom[i] in toKeepRings:
                ringAtoms.append(mol.GetAtom(OEHasAtomIdx(i)))

        return ringAtoms

    @staticmethod
    def GetFuncGroups(mol):
        """

        :param mol:
        :return:
        """
        funcGrps = []
        for funcGrp in OEGetFuncGroupFragments(mol):
            if OECount(funcGrp, OEIsHeavy()) > 5:
                continue
            if OECount(funcGrp, OEIsHetero()) == 0:
                continue
            if OECount(funcGrp, OEAtomIsInRing()) > 0:
                continue

            funcGrps.append(OEAtomBondSet(funcGrp))

        return funcGrps

    @staticmethod
    def GetTorsions(mol):
        """
        Goes through each rotatable bond in the molecule
        and extracts torsion atoms (a-b-c-d)

        Core torsion atoms are extended by one bond
        If core or extended atoms are part of a ring,
        then entire ring is kept

        Keep ortho substitution

        Keep functional groups that have at least one atom overlap
        with the core/extended torsion atoms

        Functional group inclusion criteria:
        - <= 5 heavy atoms
        - must contain at least one hetero atom
        - non-ring

        Add methyl cap if bond involving hetero atom is broken

        @param mol: OEGraphMol
        @type mol: OEGraphMol
        @return: list[OEGraphMol]
        """
        OEAssignHybridization(mol)
        funcGrps = TorsionGenerator.GetFuncGroups(mol)
        includedTorsions = OEAtomBondSet()
        torsionMols = []
        for atom in mol.GetAtoms():
            atom.SetData("idx", atom.GetIdx() + 1)

        for torsion in OEGetTorsions(mol, IsRotor()):
            if (
                torsion.a.IsHydrogen()
                or torsion.b.IsHydrogen()
                or torsion.c.IsHydrogen()
                or torsion.d.IsHydrogen()
            ):
                continue

            torsion_bond = mol.GetBond(torsion.b, torsion.c)
            if includedTorsions.HasBond(torsion_bond):
                continue
            # if includedTorsions.HasAtom(torsion.b) and \
            #     includedTorsions.HasAtom(torsion.c):
            #     continue

            # revert map idx to zero in original mol
            for atom in mol.GetAtoms():
                atom.SetMapIdx(0)

            # includedTorsions.AddAtom(torsion.b)
            # includedTorsions.AddAtom(torsion.c)
            includedTorsions.AddBond(torsion_bond)

            torsionSet = OEAtomBondSet(mol.GetBonds())
            torsionSet.AddAtoms([torsion.a, torsion.b, torsion.c, torsion.d])
            for atom in torsionSet.GetAtoms():
                atom.SetMapIdx(1)

            # extend core torsion atoms by one bond
            nbrs = TorsionGenerator.GetNbrs(torsionSet)
            torsionSet.AddAtoms(nbrs)

            # include ring atoms
            ringAtoms = TorsionGenerator.GetSameRingAtoms(mol, torsionSet)
            torsionSet.AddAtoms(ringAtoms)

            for atom in torsionSet.GetAtoms():
                if not atom.GetMapIdx() == 1:
                    atom.SetMapIdx(2)

            # add functional groups that overlap with torsion set
            TorsionGenerator.AddFuncGroupAtoms(funcGrps, torsionSet)

            # add relevant ring atoms (ortho substituents and ring H)
            TorsionGenerator.AddRelevantRingAtoms(mol, torsion, torsionSet)

            # special treatment for C=O
            for atom in torsionSet.GetAtoms(
                OEAndAtom(OEIsOxygen(), OEIsAtomHybridization(OEHybridization_sp2))
            ):
                for nbr in atom.GetAtoms():
                    if torsionSet.HasAtom(nbr):
                        for nbr2 in nbr.GetAtoms(OEIsHeavy()):
                            if not torsionSet.HasAtom(nbr2):
                                nbr2.SetMapIdx(2)
                                torsionSet.AddAtom(nbr2)

            # mark bridging atom and cap if needed
            BRIDGE_ATOM_IDX = 4
            TorsionGenerator.MarkBridgingAtoms(BRIDGE_ATOM_IDX, mol, torsionSet)

            A_IDX = 11
            B_IDX = 12
            C_IDX = 13
            D_IDX = 14
            torsion.a.SetMapIdx(A_IDX)
            torsion.b.SetMapIdx(B_IDX)
            torsion.c.SetMapIdx(C_IDX)
            torsion.d.SetMapIdx(D_IDX)

            torsionMol = OEGraphMol()
            OESubsetMol(torsionMol, mol, torsionSet, True)
            torsionMol.Sweep()
            torsionMols.append(torsionMol)

            # change bridge atom to Carbon
            for atom in torsionMol.GetAtoms(OEHasMapIdx(BRIDGE_ATOM_IDX)):
                atom.SetAtomicNum(OEElemNo_C)
                # adjust implicit atom count
                explicit_valence = atom.GetExplicitValence()
                if explicit_valence < 4:
                    atom.SetImplicitHCount(4 - explicit_valence)

            TorsionGenerator.SetSDData(A_IDX, B_IDX, C_IDX, D_IDX, torsion, torsionMol)

            # set map idx to zero in torsion mol
            for atom in torsionMol.GetAtoms():
                atom.SetMapIdx(0)

        return torsionMols

    @staticmethod
    def SetSDData(A_IDX, B_IDX, C_IDX, D_IDX, torsion, torsionMol):
        taIdx = torsionMol.GetAtom(OEHasMapIdx(A_IDX)).GetIdx() + 1
        tbIdx = torsionMol.GetAtom(OEHasMapIdx(B_IDX)).GetIdx() + 1
        tcIdx = torsionMol.GetAtom(OEHasMapIdx(C_IDX)).GetIdx() + 1
        tdIdx = torsionMol.GetAtom(OEHasMapIdx(D_IDX)).GetIdx() + 1
        apStr = "cs1:0:1;1%{}:1%{}:1%{}:1%{}".format(taIdx, tbIdx, tcIdx, tdIdx)
        OESetSDData(torsionMol, "TORSION_ATOMPROP", apStr)
        fragTorAtoms = "{} {} {} {}".format(taIdx, tbIdx, tcIdx, tdIdx)
        OESetSDData(torsionMol, TORSION_ATOMS_FRAG_KEY, fragTorAtoms)
        parentTorAtoms = "{} {} {} {}".format(
            torsion.a.GetIdx() + 1,
            torsion.b.GetIdx() + 1,
            torsion.c.GetIdx() + 1,
            torsion.d.GetIdx() + 1,
        )
        OESetSDData(torsionMol, TORSION_ATOMS_PARENT_MOL_KEY, parentTorAtoms)

        atom_map = ""
        for atom in torsionMol.GetAtoms():
            atom_map += str(atom.GetIdx() + 1) + "_" + str(atom.GetData("idx")) + "-"
        atom_map = atom_map[:-1]
        OESetSDData(torsionMol, FRAGMENT_TO_PARENT_ATOMS_KEY, atom_map)

    @staticmethod
    def MarkBridgingAtoms(BRIDGE_ATOM_IDX, mol, torsionSet):
        NorOorS = OEOrAtom(OEOrAtom(OEIsNitrogen(), OEIsOxygen()), OEIsSulfur())
        for atom in mol.GetAtoms(OEAndAtom(OEHasMapIdx(2), NorOorS)):
            for nbr in atom.GetAtoms(OEIsHeavy()):
                if not torsionSet.HasAtom(nbr):
                    if nbr.GetMapIdx() == 0:
                        torsionSet.AddAtom(nbr)
                        if nbr.GetHvyDegree() == 1:
                            nbr.SetMapIdx(3)
                            continue

                        nbr.SetMapIdx(BRIDGE_ATOM_IDX)

    @staticmethod
    def AddRelevantRingAtoms(mol, torsion, torsionSet):
        atom1or2 = OEOrAtom(OEHasMapIdx(1), OEHasMapIdx(2))
        ringNbrs = []
        for atom in mol.GetAtoms(OEAndAtom(OEAtomIsInRing(), atom1or2)):
            for nbr in atom.GetAtoms(
                OEAndAtom(OENotAtom(atom1or2), OENotAtom(OEAtomIsInRing()))
            ):
                if nbr.IsHydrogen():
                    ringNbrs.append(nbr)
                    continue

                if TorsionGenerator.IsOrtho(nbr, torsion):
                    ringNbrs.append(nbr)
        for nbr in ringNbrs:
            if not torsionSet.HasAtom(nbr):
                nbr.SetMapIdx(2)
                torsionSet.AddAtom(nbr)

    @staticmethod
    def AddFuncGroupAtoms(funcGrps, torsionSet):
        addGrps = []
        for funcGrp in funcGrps:
            for atom in funcGrp.GetAtoms():
                if torsionSet.HasAtom(atom):
                    addGrps.append(funcGrp)
                    break
        for grp in addGrps:
            for atom in grp.GetAtoms():
                if not torsionSet.HasAtom(atom):
                    atom.SetMapIdx(2)
                    torsionSet.AddAtom(atom)

    @staticmethod
    def GetMinPathLength(refTorsion, adjTorsion):
        """
        Returns path length between the two torsions

        @param refTorsion: OETorsion
        @param adjTorsion: OETorsion
        @return: int
        """
        minPathLen = 1000
        for refAtom in [refTorsion.b, refTorsion.c]:
            for torAtom in [adjTorsion.b, adjTorsion.c]:
                pathLen = OEGetPathLength(refAtom, torAtom)
                if pathLen < minPathLen:
                    minPathLen = pathLen

        return minPathLen

    @staticmethod
    def GetAdjacentTorsions(mol, refTorsion):
        """
        Returns all torsions that are 0 or 1 path length away from
        the reference torsion

        @param mol: OEGraphMol
        @param refTorsion: OETorsion
        @return: int
        """
        adjTorsions = []
        PATH_LENGTH_THRESHOLD = 1
        torset = {str(refTorsion.b.GetIdx()) + "_" + str(refTorsion.c.GetIdx()): True}
        torset[str(refTorsion.c.GetIdx()) + "_" + str(refTorsion.b.GetIdx())] = True
        pred = OEAndBond(OEHasOrder(1), OENotBond(OEBondIsInRing()))
        for adjTorsion in OEGetTorsions(mol, pred):
            # skip nitrile
            order_ab = adjTorsion.a.GetBond(adjTorsion.b).GetOrder()
            order_cd = adjTorsion.c.GetBond(adjTorsion.d).GetOrder()
            if order_ab == 3 or order_cd == 3:
                continue

            # skip torsions involving terminal -N-H
            if adjTorsion.a.IsHydrogen() and adjTorsion.b.IsNitrogen():
                continue
            if adjTorsion.d.IsHydrogen() and adjTorsion.c.IsNitrogen():
                continue

            key1 = str(adjTorsion.b.GetIdx()) + "_" + str(adjTorsion.c.GetIdx())
            key2 = str(adjTorsion.c.GetIdx()) + "_" + str(adjTorsion.b.GetIdx())
            if key1 in torset or key2 in torset:
                continue

            pathLen = TorsionGenerator.GetMinPathLength(refTorsion, adjTorsion)
            if pathLen <= PATH_LENGTH_THRESHOLD:
                adjTorsions.append(adjTorsion)
                torset[key1] = True
                torset[key2] = True

        return adjTorsions


def get_molecule_torsion_fragments(mol):
    # generate torsion fragments from the input molecule
    torgen = TorsionGenerator()
    tormols = torgen.GetTorsions(mol)

    ## process torsion fragments using corina
    ## add missing hydrogens and neutralize
    # fragfile = tempfile.NamedTemporaryFile(suffix=".sdf").name
    # ofs = oemolostream(fragfile)
    # for tormol in tormols:
    #    if OECount(tormol, OEIsHeavy()) > 25:
    #        continue
    #    OEWriteMolecule(ofs, tormol)
    # ofs.close()
    #
    # corinafile = tempfile.NamedTemporaryFile(dir=SCRATCH_DIR, suffix=".sdf").name
    # corinaProg = "{} {} {}".format(CORINA_CMD, fragfile, corinafile)
    ## print("Running corina: ", corinafile)
    # subprocess.call(corinaProg, shell=True)
    #
    # if os.path.exists(corinafile):
    #    # retrieve corina output molecules
    #    ifs = oemolistream(corinafile)
    #    frag_mols = []
    #    for mol in ifs.GetOEGraphMols():
    #        frag_mols.append(OEGraphMol(mol))
    #
    #    ifs.close()
    #
    #    return frag_mols
    # else:
    #    return tormols
    return tormols


def gen_torsion_fragments(mol):
    return get_molecule_torsion_fragments(mol)


def get_fragment_to_parent_atom_mapping(parent_mol, frag_mol):
    try:
        mapping_data = OEGetSDData(frag_mol, FRAGMENT_TO_PARENT_ATOMS_KEY)
        idx_map = dict(
            map(int, idx_pair.split("_")) for idx_pair in mapping_data.split("-")
        )
        atom_map = {}
        for frag_idx, parent_idx in idx_map.items():
            frag_atom = frag_mol.GetAtom(OEHasAtomIdx(frag_idx - 1))
            parent_atom = parent_mol.GetAtom(OEHasAtomIdx(parent_idx - 1))
            if frag_atom is not None and parent_atom is not None:
                atom_map[frag_atom] = parent_atom

        return atom_map
    except Exception as e:
        return {}


def get_modified_inchi_key(mol, atoms):
    """
    Generates InChIKey for the input molecule.
    Passed atoms of the molecule are mutated (atomic number changed)
    based on the mapping defined in the function.

    @param mol:
    @param atoms:
    @type mol: OEGraphMol
    @type atoms: list[OEAtombase]
    @return: str
    """
    copy_mol = OEGraphMol(mol)
    atom_map = {
        OEElemNo_C: OEElemNo_Pb,
        OEElemNo_N: OEElemNo_Bi,
        OEElemNo_O: OEElemNo_Po,
        OEElemNo_F: OEElemNo_At,
        OEElemNo_S: OEElemNo_Te,
        OEElemNo_Cl: OEElemNo_I,
        OEElemNo_P: OEElemNo_Sb,
        OEElemNo_Br: 117,
        OEElemNo_I: 118,
    }
    for ref_atom in atoms:
        copy_atom = copy_mol.GetAtom(OEHasAtomIdx(ref_atom.GetIdx()))
        if copy_atom is None:
            raise Exception("Null atom found")
        copy_atom.SetAtomicNum(atom_map[copy_atom.GetAtomicNum()])

    return OECreateInChIKey(copy_mol)


def CreateInchiKeyPlus(mol):
    """

    @param mol: molecule with torsion atom index
    @type mol: OEGraphMol
    @return: str
    """
    inchiKey = OECreateInChIKey(mol)
    try:
        atomIndices = list(
            map(int, OEGetSDData(mol, TORSION_ATOMS_FRAG_KEY).strip().split())
        )
        a = mol.GetAtom(OEHasAtomIdx(atomIndices[0] - 1))
        b = mol.GetAtom(OEHasAtomIdx(atomIndices[1] - 1))
        c = mol.GetAtom(OEHasAtomIdx(atomIndices[2] - 1))
        d = mol.GetAtom(OEHasAtomIdx(atomIndices[3] - 1))

        ad = a.GetAtomicNum() * d.GetAtomicNum()
        bc = b.GetAtomicNum() * c.GetAtomicNum()
        adAro = int(a.IsAromatic()) + int(d.IsAromatic())
        bcAro = int(b.IsAromatic()) + int(c.IsAromatic())

        count1 = len(list(OEGetSubtree(b, c)))
        count2 = len(list(OEGetSubtree(c, b)))
        count = count1 * count2

        inchiKey = inchiKey + str(ad) + str(bc) + str(adAro) + str(bcAro) + str(count)
    except Exception as e:
        print("CreateInchiKeyPlus ", e)

    return inchiKey


def CreateTorsionInchiKey(mol):
    """

    @param mol: molecule with torsion atom index
    @type mol: OEGraphMol
    @return: str|None
    """
    try:
        tormol = OEGraphMol(mol)
        atomIndices = list(
            map(int, OEGetSDData(tormol, TORSION_ATOMS_FRAG_KEY).strip().split())
        )
        a = tormol.GetAtom(OEHasAtomIdx(atomIndices[0] - 1))
        b = tormol.GetAtom(OEHasAtomIdx(atomIndices[1] - 1))
        c = tormol.GetAtom(OEHasAtomIdx(atomIndices[2] - 1))
        d = tormol.GetAtom(OEHasAtomIdx(atomIndices[3] - 1))

        replace_map = {
            OEElemNo_C: OEElemNo_Pb,
            OEElemNo_N: OEElemNo_Bi,
            OEElemNo_O: OEElemNo_Po,
            OEElemNo_F: OEElemNo_At,
            OEElemNo_S: OEElemNo_Te,
            OEElemNo_Cl: OEElemNo_I,
            OEElemNo_P: OEElemNo_Sb,
        }
        # bond = tormol.GetBond(b, c)
        # tormol.DeleteBond(bond)

        # for hatom in tormol.GetAtoms(OEIsHydrogen()):
        #    hatom.SetData("originalH", True)

        # fill the valence. otherwise, inchi generation will throw warnings
        # OEAssignImplicitHydrogens(tormol)
        # OEAddExplicitHydrogens(tormol)
        # for hatom in tormol.GetAtoms(OEIsHydrogen()):
        #    if not hatom.HasData("originalH"):
        #        hatom.SetAtomicNum(OEElemNo_Li)

        # inchiKey = OECreateInChIKey(mol) + OECreateInChIKey(tormol)
        b.SetAtomicNum(replace_map[b.GetAtomicNum()])
        c.SetAtomicNum(replace_map[c.GetAtomicNum()])
        inchiKey = OECreateInChIKey(mol) + OECreateInChIKey(tormol)
        return inchiKey
    except Exception as e:
        print("CreateTorsionInchiKey ", e)
        return None

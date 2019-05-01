import logging
import collections
import numpy as np
from openeye import oechem
from torsion.utils import write_energy_profile_to_sddata
from torsion.utils import get_sd_data
from torsion.utils import get_torsion_oeatom_list
from torsion.core import normalize_coordinates, get_dihedral
from torsion.utils import get_modified_inchi_key
from torsion.utils import get_molecule_torsion_fragments
from torsion.utils import get_fragment_to_parent_atom_mapping
from torsion.utils import get_modified_molecule_inchi
from torsion.utils import generate_energy_profile_sd_data_1d, ENERGY_PROFILE_TAG

STRAIN_TAG = 'STRAIN'

def construct_dihedral_energy_profile(torsion_conformers, num_points=24):
    angle_list = np.array([360*i/num_points for i in range(num_points)])

    num_confs = 0
    profile = np.full(num_points, np.nan)
    for mol in torsion_conformers:
        if not mol:
            continue
        num_confs += 1
        conf = mol.GetActive()
        conf_title = get_sd_data(conf, "CONFORMER_LABEL")
        tor_atoms = get_sd_data(mol, "TORSION_ATOMS_ParentMol").split()
        parent_name = conf_title[:-3]
        dih_label = '_'.join(str(x) for x in tor_atoms)
        fragment_label = parent_name + '_' + dih_label
        angle_idx = int(conf_title[-2:])
        
        profile[angle_idx] = np.float(get_sd_data(conf, 'PSI4_ENERGY')) 
        logging.debug("angle_idx: %d", angle_idx)
        logging.debug("Psi4 Energy: %f", float(get_sd_data(conf, 'PSI4_ENERGY')))

    # check for angles where no energies are available
    for angle in angle_list[np.all(np.isnan(profile))]:
        logging.warning("Warning: No energies found for angle {:.1f} for fragment: {}".format(angle, fragment_label))

    # calculate relative energies
    min_energy = np.nanmin(profile)
    profile -= min_energy
    profile[np.isnan(profile)] = -1     # set nans to -1
    torsional_strain = np.column_stack(( angle_list,
                                         profile
                                        ))
    
    # combine conformers
    output_conformers = oechem.OEMol(torsion_conformers[0])
    output_conformers.DeleteConfs()
    title = fragment_label
    output_conformers.SetTitle(title)

    # setup normalization
    torsion_tag = 'TORSION_ATOMS_FRAGMENT'
    torsion_atoms_in_fragment = get_sd_data(mol, torsion_tag).split()
    print(torsion_atoms_in_fragment)
    dihedral_atom_indices = [int(x)-1 for x in torsion_atoms_in_fragment]
    dih, _ = get_dihedral(output_conformers, dihedral_atom_indices)

    for old_conf in torsion_conformers:
        if old_conf:
            new_conf = output_conformers.NewConf(old_conf)
            normalize_coordinates(new_conf, dih)
            oechem.OEClearSDData(new_conf)
            for dp in oechem.OEGetSDDataPairs(old_conf.GetActive()):
                if dp.GetTag() not in ['OEConfTitle', 'CONFORMER_LABEL']:
                    oechem.OESetSDData(new_conf, dp.GetTag(), dp.GetValue())
            torsion_angle = get_sd_data(old_conf, 'TORSION_ANGLE')
            title = fragment_label +\
                    ': Angle ' + torsion_angle
            new_conf.SetTitle(title)

    write_energy_profile_to_sddata(output_conformers, torsional_strain.copy())
    
    # Calculate all possible torsion inchi keys for this fragment
    torsion_inchi_list = []
    inchi_key = oechem.OECreateInChIKey(output_conformers)
    _, b, c, _ = get_torsion_oeatom_list(output_conformers)
    for a in b.GetAtoms(oechem.OEIsHeavy()):
        for d in c.GetAtoms(oechem.OEIsHeavy()):
            if a.GetIdx() == c.GetIdx() or d.GetIdx() == b.GetIdx():
                continue

            torsion_inchi = inchi_key + get_modified_inchi_key(output_conformers, [a, b, c, d])
            torsion_inchi_list.append(torsion_inchi)

    return output_conformers, torsional_strain, torsion_inchi_list


def cal_molecule_torsion_strain(mol, profiles_map):
    '''

    @type mol: oechem.OEGraphMol|oechem.OEMol
    :param mol:
    :param profiles_map:
    :return:
    '''
    if type(mol) is oechem.OEMol:
        graph_mol = oechem.OEGraphMol(mol.GetActive())
        data = extract_molecule_torsion_data(graph_mol)
    else:
        data = extract_molecule_torsion_data(mol)

    if data is not None:
        _, tor_map = data

        if type(mol) is oechem.OEGraphMol:
            for tor_inchi, tor_data_list in tor_map.items():
                if tor_inchi in profiles_map:
                    for tor_data in tor_data_list:
                        _, b_idx, c_idx, _, angle = tor_data
                        bond = mol.GetBond(mol.GetAtom(oechem.OEHasAtomIdx(b_idx)),
                                           mol.GetAtom(oechem.OEHasAtomIdx(c_idx)))
                        if bond is not None:
                            strain_energy = profiles_map[tor_inchi](angle)
                            if strain_energy < 0:
                                strain_energy = 0
                            if bond.HasData(STRAIN_TAG) and bond.GetData(STRAIN_TAG) > strain_energy:
                                bond.SetData(STRAIN_TAG, strain_energy)

            total_strain = 0.0
            for bond in mol.GetBonds():
                if bond.HasData(STRAIN_TAG):
                    total_strain += bond.GetData(STRAIN_TAG)

            mol.SetData(STRAIN_TAG, total_strain)

        elif type(mol) is oechem.OEMol:
            for conf in mol.GetConfs():
                bondIdx2energy = {}; bondIdx2profile = {};
                bondIdx2toratoms = {}; bondIdx2angles = {}
                for tor_inchi, tor_data_list in tor_map.items():
                    if tor_inchi in profiles_map:
                        for tor_data in tor_data_list:
                            a_idx, b_idx, c_idx, d_idx, _ = tor_data
                            b_atm = conf.GetAtom(oechem.OEHasAtomIdx(b_idx))
                            c_atm = conf.GetAtom(oechem.OEHasAtomIdx(c_idx))
                            bond = conf.GetBond(b_atm, c_atm)
                            if bond is not None:
                                a_atm = conf.GetAtom(oechem.OEHasAtomIdx(a_idx))
                                d_atm = conf.GetAtom(oechem.OEHasAtomIdx(d_idx))
                                angle = oechem.OEGetTorsion(conf, a_atm, b_atm, c_atm, d_atm)*oechem.Rad2Deg
                                strain_energy = float(profiles_map[tor_inchi](angle))
                                if strain_energy < 0:
                                    strain_energy = 0

                                x = range(-165,181,15)
                                y = []
                                for a in x:
                                    y.append(float(profiles_map[tor_inchi](a)))

                                energy_profile = generate_energy_profile_sd_data_1d(list(zip(x, y)))
                                bondIdx = bond.GetIdx()
                                if bondIdx not in bondIdx2energy:
                                    bondIdx2energy[bondIdx] = strain_energy
                                    bondIdx2profile[bondIdx] = energy_profile
                                    bondIdx2toratoms[bondIdx] = [a_idx+1, b_idx+1, c_idx+1, d_idx+1]
                                    bondIdx2angles[bondIdx] = angle
                                elif bondIdx2energy[bondIdx] > strain_energy:
                                    bondIdx2energy[bondIdx] = strain_energy
                                    bondIdx2profile[bondIdx] = energy_profile
                                    bondIdx2toratoms[bondIdx] = [a_idx+1, b_idx+1, c_idx+1, d_idx+1]
                                    bondIdx2angles[bondIdx] = angle

                # sd property place holder
                oechem.OESetSDData(conf, 'QM_STRAIN', '0.0')
                oechem.OESetSDData(conf, 'NUM_QM_TORSION_PROFILES', '0')
                oechem.OESetSDData(conf, 'NUM_MISSING_QM_TORSIONS', '-1')

                total_strain = 0.0
                tor_idx = 1
                for tor_count, bond in enumerate(conf.GetBonds(oechem.OEIsRotor())):
                    bidx = bond.GetIdx()
                    if bidx in bondIdx2energy:
                        total_strain += bondIdx2energy[bidx]

                        tmp = ':1%'.join(list(map(str, bondIdx2toratoms[bidx])))
                        tor_atomprop = 'cs1:0:1;1%' + tmp
                        oechem.OESetSDData(conf, 'QM_TORSION_ATOMS_%d_'%tor_idx + 'ATOMPROP', tor_atomprop)
                        oechem.OESetSDData(conf, 'QM_TORSION_%d_'%tor_idx + ENERGY_PROFILE_TAG, bondIdx2profile[bidx])
                        oechem.OESetSDData(conf, 'TORSION_ANGLE_%d'%tor_idx, '%.1f'%bondIdx2angles[bidx])
                        oechem.OESetSDData(conf, 'QM_STRAIN_TORSION_%d'%tor_idx, '%.1f'%bondIdx2energy[bidx])
                        tor_idx += 1

                oechem.OESetSDData(conf, 'QM_STRAIN', '%.1f'%total_strain)
                oechem.OESetSDData(conf, 'NUM_QM_TORSION_PROFILES', '%d'%(tor_count+1))
                num_missing_torsions = (tor_count + 1) - (tor_idx - 1)
                oechem.OESetSDData(conf, 'NUM_MISSING_QM_TORSIONS', '%d'%num_missing_torsions)
                conf.SetData(STRAIN_TAG, total_strain)


def get_dihedral_inchi_key(mol):
    try:
        a, b, c, d = get_torsion_oeatom_list(mol)

        ad = a.GetAtomicNum() * d.GetAtomicNum()
        bc = b.GetAtomicNum() * c.GetAtomicNum()
        adAro = int(a.IsAromatic()) + int(d.IsAromatic())
        bcAro = int(b.IsAromatic()) + int(c.IsAromatic())

        count1 = len(list(oechem.OEGetSubtree(b, c)))
        count2 = len(list(oechem.OEGetSubtree(c, b)))
        count = count1 * count2

        inchiKey = oechem.OECreateInChIKey(mol)
        inchiKey = inchiKey + str(ad) + str(bc) + str(adAro) + str(bcAro) + str(count)
        return inchiKey
    except Exception as e:
        logging.warning(e)
        return None


def get_generic_dihedral_inchi_key(mol):
    '''
    generates dihedral inchi key after mutating only the central two atoms
    '''
    try:
        _, b, c, _ = get_torsion_oeatom_list(mol)
        modified_inchi = get_modified_inchi_key(mol, [b, c])
        inchiKey = oechem.OECreateInChIKey(mol) + modified_inchi
        return inchiKey
    except Exception as e:
        logging.warning(e)
        return None

    
def get_specific_dihedral_inchi_key(mol):
    '''
    generates unique dihedral inchi key by mutating all four dihedral atoms 
    '''
    try:
        a, b, c, d = get_torsion_oeatom_list(mol)
        modified_inchi = get_modified_inchi_key(mol, [a, b, c, d])
        inchiKey = oechem.OECreateInChIKey(mol) + modified_inchi
        return inchiKey
    except Exception as e:
        logging.warning(e)
        return None

    
def extract_molecule_torsion_data(parent_mol, frag_mols = None):
    '''
    extract dihedral angle associated with each torsion motif in the input molecule
    Torsion motifs are represented using generic modified inchi (central two atoms)
    and specific modified inchi (4 torsion atoms)

    @param parent_mol:
    @type parent_mol: oechem.OEGraphMol
    @return: tuple(str, dict[str, list[float]])
    '''
    if frag_mols is None:
        frag_mols = get_molecule_torsion_fragments(parent_mol)

    torsion_data = collections.defaultdict(list)
    for frag_mol in frag_mols:
        inchi_key = oechem.OECreateInChIKey(frag_mol)
        atom_map = get_fragment_to_parent_atom_mapping(parent_mol, frag_mol)

        try:
            _, b, c, _ = get_torsion_oeatom_list(frag_mol)

            for a in b.GetAtoms(oechem.OEIsHeavy()):
                for d in c.GetAtoms(oechem.OEIsHeavy()):
                    if a.GetIdx() == c.GetIdx() or d.GetIdx() == b.GetIdx():
                        continue

                    ap = atom_map[a]
                    bp = atom_map[b]
                    cp = atom_map[c]
                    dp = atom_map[d]

                    if a.GetAtomicNum() == ap.GetAtomicNum() \
                            and b.GetAtomicNum() == bp.GetAtomicNum() \
                            and c.GetAtomicNum() == cp.GetAtomicNum() \
                            and d.GetAtomicNum() == dp.GetAtomicNum():
                        angle = oechem.OEGetTorsion(parent_mol, ap, bp, cp, dp)*oechem.Rad2Deg
                        torsion_inchi = inchi_key + get_modified_inchi_key(
                                                            frag_mol, [a, b, c, d])

                        torsion_data[torsion_inchi].append((ap.GetIdx(), bp.GetIdx(), cp.GetIdx(), dp.GetIdx(), angle))

        except Exception as e:
            logging.warning(e)
            continue

    parent_inchi = get_modified_molecule_inchi(parent_mol)

    return (parent_inchi, torsion_data)

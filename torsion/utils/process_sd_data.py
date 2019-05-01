from openeye import oechem
import numpy as np
from scipy.interpolate import interp1d

ENERGY_PROFILE_TAG = "energy_profile_XYPLOT"
TORSION_ATOMS_FRAGMENT_TAG = "TORSION_ATOMS_FRAGMENT"
SPECIFIC_INCHI_TAG = "specific_inchi"


def has_sd_data(mol, tag):
    if oechem.OEHasSDData(mol, tag):
        return True
    if oechem.OEHasSDData(mol.GetActive(), tag):
        return True
    return False


def get_sd_data(mol, tag):
    try:
        if oechem.OEHasSDData(mol, tag):
            return oechem.OEGetSDData(mol, tag)
        if oechem.OEHasSDData(mol.GetActive(), tag):
            return oechem.OEGetSDData(mol.GetActive(), tag)
    except AttributeError as e:
        print(e)
        return ""


def get_torsion_oeatom_list(mol, tag="TORSION_ATOMS_FRAGMENT"):
    if has_sd_data(mol, tag):
        torsion_atoms = get_sd_data(mol, tag)
        try:
            torsion_atoms_idx = list(map(int, torsion_atoms.split()))
            torsion_oeatoms = map(
                lambda idx: mol.GetAtom(oechem.OEHasAtomIdx(idx - 1)), torsion_atoms_idx
            )
            return list(torsion_oeatoms)
        except Exception as e:
            print(e)
            return None


def get_torsion_oebond(mol, tag="TORSION_ATOMS_FRAGMENT"):
    torsion_atoms = get_torsion_oeatom_list(mol, tag)
    try:
        return mol.GetBond(torsion_atoms[1], torsion_atoms[2])
    except Exception as e:
        print(e)
        return None


def delete_sd_data(mol, tag, locator_tag):
    if oechem.OEHasSDData(mol, locator_tag):
        return oechem.OEDeleteSDData(mol, tag)
    elif oechem.OEHasSDData(mol.GetActive(), locator_tag):
        return oechem.OEDeleteSDData(mol.GetActive(), tag)
    return False


def dump_sd_data(mol):
    print("Data Attached at the molecule level:")
    for dp in oechem.OEGetSDDataPairs(mol):
        print(dp.GetTag(), ":", dp.GetValue())
    if type(mol) == oechem.OEMol:
        print("\n\n" + 10 * "-" + "Data Attached to Conformers:")
        for conf_id, conf in enumerate(mol.GetConfs()):
            print("Data attached to conformer {}:".format(conf_id))
            for dp in oechem.OEGetSDDataPairs(conf):
                print(dp.GetTag(), ":", dp.GetValue())
    print()


def print_torsion(mol):
    print("The Dihedral of interest is: " + get_sd_data(mol, "TORSION_ATOMS_FRAGMENT"))


def generate_energy_profile_sddata(energy_profile):
    """Generates sd data string to store an energy profile
    
    Inputs
        energy_profile: Numpy array with 2 rows
    
    Outputs
        sddata: string to write to SDF file
    """

    # Test arguments
    if energy_profile.shape[1] != 2:
        raise ValueError("Energy profile does not have 2 columns.")

    # Map angle onto -180 to 180
    energy_profile[:, 0] = np.mod(energy_profile[:, 0], 360)
    energy_profile[energy_profile[:, 0] > 180, 0] = (
        energy_profile[energy_profile[:, 0] > 180, 0] - 360
    )

    # Sort
    idx = np.argsort(energy_profile[:, 0])
    energy_profile = energy_profile[idx, :]

    # Get relative energies
    energy_profile[:, 1] = energy_profile[:, 1] - np.min(energy_profile[:, 1])

    # add a row for energy at -180
    relE_at_180 = energy_profile[-1, 1]
    energy_profile = np.vstack((np.array([-180, relE_at_180]), energy_profile))

    ymax = np.max(energy_profile[:, 1]) * 1.1

    # Generate sddata string
    sddata = (
        "Torsion scan profile (data=["
        + ",".join(["{:.2f}@{:.2f}".format(row[1], row[0]) for row in energy_profile])
        + "], "
        + "xlab=[Torsion angle (deg)], xScale=[-180.0, 180.0], "
        + "ylab=[Energy (kcal/mol)], yScale=[0,{:.2f}], type=[spline])".format(ymax)
    )

    return sddata


def write_energy_profile_to_sddata(mol, energy_profile):
    """
        Writes energy profile to a OEGraphMol object.
        This will fail for multi-conformer OEMol objects.
    """
    data = generate_energy_profile_sddata(energy_profile)
    return oechem.OESetSDData(mol, ENERGY_PROFILE_TAG, data)


def save_sddata(mol, data_kv_pairs):
    """
        Writes data from the input dictionary as SD property
        If the input is a multi-conformer molecule, the data
        will be written to each conformer

    :param mol: OEMol|OEGraphMol
    :param data_kv_pairs: dict[str, str]
    :return: None
    """
    if type(mol) == oechem.OEGraphMol:
        for k, v in data_kv_pairs.items():
            oechem.OESetSDData(mol, k, str(v))
    elif type(mol) == oechem.OEMol:
        for conf in mol.GetConfs():
            for k, v in data_kv_pairs.items():
                oechem.OESetSDData(conf, k, str(v))
    else:
        print("Unknown object type encountered in save_sddata")


def sanitize_fragment(mol):
    approved_tags = [
        "TORSION_ATOMPROP",
        "TORSION_ATOMS_FRAGMENT",
        "TORSION_ATOMS_ParentMol",
        "COUNT",
    ]

    for dp in oechem.OEGetSDDataPairs(mol):
        if dp.GetTag() in approved_tags:
            continue
        oechem.OEDeleteSDData(mol, dp.GetTag())


def is_amide_torsion(mol):
    torsion_atoms = get_torsion_oeatom_list(mol, "TORSION_ATOMS_FRAGMENT")
    try:
        if (
            torsion_atoms[1].IsCarbon()
            and oechem.OEHasDoubleBondO(torsion_atoms[1])
            and torsion_atoms[2].IsNitrogen()
            and torsion_atoms[2].GetExplicitHCount() > 0
        ):
            return True
        if (
            torsion_atoms[1].IsNitrogen()
            and torsion_atoms[1].GetExplicitHCount() > 0
            and torsion_atoms[2].IsCarbon()
            and oechem.OEHasDoubleBondO(torsion_atoms[2])
        ):
            return True
    except IndexError as e:
        print(e)
        return False

    return False


def has_protonated_double_bonded_ring_nitrogen(mol):
    for atom in mol.GetAtoms(
        oechem.OEAndAtom(
            oechem.OEAndAtom(oechem.OEIsNitrogen(), oechem.OEHasFormalCharge(1)),
            oechem.OEAtomIsInRing(),
        )
    ):
        for bond in atom.GetBonds():
            if bond.GetOrder() == 2:
                return True

    return False


def has_nitrogen_with_negative_formal_charge(mol):
    for atom in mol.GetAtoms(oechem.OEIsNitrogen()):
        if atom.GetFormalCharge() < 0:
            return True

    return False


def get_profile_xy_data(mol, tag=ENERGY_PROFILE_TAG):
    if not has_sd_data(mol, tag):
        return None

    profile = get_sd_data(mol, tag)
    return extract_numeric_data_from_profile_str(profile)


def extract_numeric_data_from_profile_str(profile):
    try:
        profile = profile.split("data=[")[1]
        profile = profile.split("]")[0]
        xyData = [map(float, item.split("@")) for item in profile.split(",")]
        y, x = zip(*xyData)

        return np.array(x), np.array(y)
    except Exception as e:
        print(e)


def get_profile_interp1d(mol):
    try:
        x, y = get_profile_xy_data(mol)
        f = interp1d(np.array(x), np.array(y), kind="cubic", bounds_error=False)
        return f
    except Exception as e:
        print(e)


def copy_QM_calc_parameters(srcMol, dstMol):
    for tag in [
        "PSI4_OPT_METHOD",
        "PSI4_OPT_BASIS",
        "PSI4_SPE_METHOD",
        "PSI4_SPE_BASIS",
        "scf_type",
        "fail_on_maxiter",
        "guess_basis",
        "use_soscf",
        "dft_radial_points",
        "dft_spherical_points",
        "num_processors",
    ]:
        value = ""
        if has_sd_data(srcMol, tag):
            value = get_sd_data(srcMol, tag)
        oechem.OESetSDData(dstMol, tag, value)


def generate_energy_profile_sd_data_1d(data):
    angles, energies = zip(*data)

    angles = list(angles)
    energies = list(energies)
    angles.insert(0, -180)
    energies.insert(0, energies[-1])

    min_energy = min(energies)
    rel_energies = [energy - min_energy for energy in energies]

    # Generate sddata string
    sddata = (
        "Torsion scan profile (data=["
        + ",".join(
            [
                "{:.2f}@{:.2f}".format(energy, angle)
                for angle, energy in zip(angles, rel_energies)
            ]
        )
        + "], "
        + "xlab=[Torsion angle (deg)], xScale=[-180.0, 180.0], "
        + "ylab=[Energy (kcal/mol)], yScale=[0,{:.2f}], type=[spline])".format(
            max(rel_energies)
        )
    )

    return sddata


def extract_energy_profile_sd_data1d(sddata):
    """Convert profile string attached as sd data to a vector of energies
    
    Parameters
    ----------
    sddata : str
    String that contains the profile as [energy@angle,energy@angle,...] within
    the first set of square angle brackets in the string.
    
    Returns
    ----------
    energies: np.array (1D)
    angles: np.array (1D)
    
    Example
    ----------
    energies, angles = extract_energy_profile_sd_data1d(sddata)

    """
    energies = np.array(
        [
            float(x[0])
            for x in [
                row.split("@") for row in sddata.split("[")[1].split("]")[0].split(",")
            ]
        ]
    )
    angles = np.array(
        [
            float(x[1])
            for x in [
                row.split("@") for row in sddata.split("[")[1].split("]")[0].split(",")
            ]
        ]
    )

    return angles, energies

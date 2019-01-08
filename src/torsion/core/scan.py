import psi4
import logging
import numpy as np
from openeye import oechem
from torsion.utils import save_sddata, get_sd_data
from torsion.psi4wrapper import get_conf_from_psi4_mol, psi4_calculation


def get_dihedral(mol, dihedralAtomIndices):
    """Get the dihedral corresponding to the indices in dihedralAtomIndices.

    Note that the indices are zero-indexed. The torsion of interest is the bond
    between atoms with indices 1 and 2.
    """
    #dihedralAtomIndices = [int(x) for x in get_sd_data(mol, 'TORSION_ATOMS_FRAGMENT').split()]
    dih = oechem.OEAtomBondSet()
    tor = oechem.OEAtomBondSet()
    for i in range(3):
        srcIdx = dihedralAtomIndices[i]
        destIdx = dihedralAtomIndices[i+1]
        src = mol.GetAtom(oechem.OEHasAtomIdx(srcIdx))
        dest = mol.GetAtom(oechem.OEHasAtomIdx(destIdx))
        dih.AddAtom(src)
        bond = mol.GetBond(src,dest)
        dih.AddBond(bond)
        if i == 1:
            tor.AddBond(bond)
    dih.AddAtom(dest)

    return dih, tor


def normalize_coordinates(mol, dih):
    """Reset the coordinates of an OEGraphMol object with respect to a dihedral.

    The first three atoms of the dihedral will be set to the following
    coordinates:
    |Atom| Position |
    |----|----------|
    |  0 | (0,0,0)  |
    |  1 | (x1,0,0) |
    |  2 | (x2,y2,0)|

    This defines the coordinate system required to determine the coordinates of
    all the other atoms in the fragment.

    Arguments:
        mol: An OEGraphMol object.
        dih: An OEAtomBondSet with atoms belonging to mol.

    Returns:
        success: A bool indicating whether the transformation was successful.
    """
    dih_atoms = [atom for atom in dih.GetAtoms()]

    # Translate the fragment until the first atom of the dihedral is at the origin
    origin = mol.GetCoords()[dih_atoms[0].GetIdx()]
    shift = oechem.OEDoubleArray([-x for x in origin])
    oechem.OETranslate(mol, shift)

    # Get coordinates of atoms 2 and 3 in the dihedral
    u = np.array(mol.GetCoords()[dih_atoms[1].GetIdx()])
    v = np.array(mol.GetCoords()[dih_atoms[2].GetIdx()])

    # Get three orthogonal unit vectors
    u_hat = u/np.linalg.norm(u)

    v_hat = v - ((np.dot(u, v)*u)/np.dot(u, u))  # Gram-Schmidt Orthogonalization
    v_hat = v_hat/np.linalg.norm(v_hat)

    w_hat = np.cross(u_hat, v_hat)

    # Assemble unit vectors into rotation matrix
    R = np.stack((u_hat, v_hat, w_hat), axis=0)
    rotation_matrix = oechem.OEDoubleArray(R.flatten())
    oechem.OERotate(mol, rotation_matrix)


def calculate_energy(mol, dih,
                     spe_method='SCF',
                     spe_basis='6-31G',
                     geom_opt_technique='None',
                     opt_method='SCF',
                     opt_basis='6-31G',
                     geom_maxiter=100,
                     only_selected_conf=False,
                     molden_output=False,
                     **psi4opts):
    """Calculates the energy for a single conformer at a single dihedral angle.
    """

    # Argument validation
    if geom_opt_technique not in ['None', 'QM', 'MM']:
        geom_opt_technique = 'None'

    parent_torsion_tag = 'TORSION_ATOMS_ParentMol'
    torsion_atoms_in_parent = get_sd_data(mol, parent_torsion_tag).split()
    dih_name = mol.GetTitle()+ '_' + '_'.join(torsion_atoms_in_parent)

    if only_selected_conf:
        conf_selection_tag = 'SELECTED_CONFORMER'
        if not mol.HasData(conf_selection_tag):
            raise ValueError("Could not find 'SELECTED_CONFORMER' Tag in %s.", dih_name)
        key_conf_id = mol.GetIntData(conf_selection_tag)

    for conf in mol.GetConfs():
        if only_selected_conf:
            if conf.GetIdx() != key_conf_id:
                continue

        conf_name = get_sd_data(conf, 'CONFORMER_LABEL')
        if only_selected_conf:
            logging.debug("Only running psi4 calculation for %s" % conf_name)
        else:
            logging.debug("Running psi4 calculation for %s" % conf_name)
        angle_deg = conf.GetDoubleData('TORSION_ANGLE')

        normalize_coordinates(conf, dih)
        try :
            energy, psi_mol, wavefcn = psi4_calculation(conf, dih, conf_name, angle_deg,
                                            spe_method,
                                            spe_basis,
                                            geom_opt_technique,
                                            opt_method,
                                            opt_basis,
                                            geom_maxiter,
                                            **psi4opts)
        except ValueError as e:
            logging.error(e)
            raise ValueError('Failed to run calculation for conformer %s angle %d.',
                    conf_name, angle_deg)

        PSI4_OPT_METHOD_KEY = 'PSI4_OPT_METHOD'
        PSI4_OPT_BASIS_KEY = 'PSI4_OPT_BASIS'

        prev_opt_method = oechem.OEGetSDData(conf, PSI4_OPT_METHOD_KEY)
        if len(prev_opt_method) > 0:
            prev_opt_method += '_'
        prev_opt_basis = oechem.OEGetSDData(conf, PSI4_OPT_BASIS_KEY)
        if len(prev_opt_basis) > 0:
            prev_opt_basis += '_'

        complete_opt_method = prev_opt_method + str(opt_method)
        complete_opt_basis = prev_opt_basis + str(opt_basis)

        get_conf_from_psi4_mol(mol, psi_mol, conf)
        logging.debug("Completed psi4 calculation for %s with energy %f" % (conf_name, energy))
        conf.SetEnergy(energy)
        oechem.OESetSDData(conf, 'PSI4_ENERGY', str(energy))
        conf.SetDoubleData('PSI4_ENERGY', energy)
        if geom_opt_technique == 'QM':
            oechem.OESetSDData(conf, PSI4_OPT_METHOD_KEY, complete_opt_method)
            oechem.OESetSDData(conf, PSI4_OPT_BASIS_KEY, complete_opt_basis)
        oechem.OESetSDData(conf, 'PSI4_SPE_METHOD', str(spe_method))
        oechem.OESetSDData(conf, 'PSI4_SPE_BASIS', str(spe_basis))

        # write wavefuction
        if molden_output:
            wf_filename = '{}_{}.molden'.format('wavefcn', conf_name)
            psi4.molden(wavefcn, wf_filename)
            try:
                with open(wf_filename, 'r') as fptr:
                    molden_data = fptr.read()
                    conf.SetData("MOLDEN_DATA", molden_data)
            except IOError as e:
                print("Unable to save wave function data ", wf_filename, " because of error: ", e)

    if psi4opts:
        save_sddata(mol, psi4opts)



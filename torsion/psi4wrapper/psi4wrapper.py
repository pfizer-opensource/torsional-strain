import sys
import time
import logging
from openeye import oechem
import psi4


def get_psi4_mol_from_conf(conf):
    """Use an OEConfBase to construct a psi4 molecule"""
    xyz = oechem.OEFloatArray(3)
    mol_str = str()
    total_charge = 0
    for atom in conf.GetAtoms():
        total_charge += atom.GetFormalCharge()
        conf.GetCoords(atom, xyz)
        # scaled_xyz = [c / psi4.constants.bohr2angstroms for c in xyz]
        mol_str += '{} {} {} {}\n'.format(oechem.OEGetAtomicSymbol(atom.GetAtomicNum()),
                                            # scaled_xyz[0], scaled_xyz[1], scaled_xyz[2])
                                            xyz[0], xyz[1], xyz[2])

    psi4_mol = psi4.geometry(mol_str)
    psi4_mol.set_molecular_charge(total_charge)
    #psi_mol.set_name('{}, conf {}'.format(conf.GetTitle(), conf.GetIdx()))
    logging.debug("\n\n\nNew Psi4 molecule created with geometry: \n" +\
            psi4_mol.create_psi4_string_from_molecule())

    return psi4_mol


def get_conf_from_psi4_mol(mol, psi_mol, conf=None):
    """
    Take a psi4 molecule and return a new conf from a molecule with coordinates
    copies
    """
    coords = psi_mol.geometry().to_array()
    coords = [[c * psi4.constants.bohr2angstroms for c in xyz] for xyz in coords]
    if conf:
        new_conf = conf
    else:
        new_conf = mol.NewConf()
    for atom in new_conf.GetAtoms():
        new_conf.SetCoords(atom, coords[atom.GetIdx()])

    new_conf.SetDimension(3)

    return new_conf


def setup_psi4(dft_radial_points=75,
              dft_spherical_points=302,
              num_processors=1,
              guess_basis='false',
              use_soscf='false',
              scf_type='DIRECT',
              fail_on_maxiter='false',
              g_convergence=None,
              max_disp_g_convergence=None):
    """Setup the control parameters for the psi4 job"""
    opts = dict()
    # opts['guess'] = 'sad'  # defaults to auto
    opts['scf_type'] = scf_type
    opts['fail_on_maxiter'] = fail_on_maxiter
    opts['basis_guess'] = guess_basis
    opts['soscf'] = use_soscf
    opts['dft_radial_points'] = dft_radial_points
    opts['dft_spherical_points'] = dft_spherical_points
    if g_convergence is not None:
        opts['g_convergence'] = g_convergence
    if max_disp_g_convergence is not None:
        opts['max_disp_g_convergence'] = max_disp_g_convergence

    psi4.set_options(opts)
    psi4.set_num_threads(num_processors)


def cleanup_psi4():
    """Cleanup the control parameters and temporary files between independent
    psi4 jobs"""
    # cleanup psi4 from previous calculations
    psi4.core.clean()
    psi4.core.clean_options()
    psi4.core.clean_variables()
    # psi4.core.opt_clean()          # causes SEG FAULTS


def run_psi4(calculation, mol, method, basis, **psi4opts):
    """Execute current psi4 setup with passed molecule.  Handle exceptions and
    stderr/stdout"""
    logging.debug('Beginning Psi4 {} calculation with {}/{}'.format(calculation,
                                                                method,
                                                                basis))

    # capture stderr & stdout and write to log
    # std_out = sys.stdout  # this is probably already handled by the cube
    # std_err = sys.stderr
    # psi_out = io.StringIO()
    # psi_err = io.StringIO()
    # sys.stdout = psi_out
    # sys.stderr = psi_err

    # execute with exception handling
    wave_fcn = None
    try:
        setup_psi4(**psi4opts)

        # logging.warn(str(os.path.abspath(psidatadir)))
        psi4.set_options({'basis': basis})
        if calculation == 'energy':
            logging.info('Psi4 call {}: {}/{}'.format(calculation, method, basis))
            ret_val, wave_fcn = psi4.energy(method, molecule=mol, return_wfn=True)
        elif calculation == 'optimize':
            logging.info('Psi4 call {}: {}/{}'.format(calculation, method, basis))
            ret_val, wave_fcn = psi4.optimize(method, molecule=mol, return_wfn=True)
        else:
            raise ValueError("Unrecognized calculation type.")

    except psi4.ValidationError as e:
        logging.error('Failed to read setup: {}'.format(e))
        ret_val, wave_fcn = None, None
    except psi4.ConvergenceError as e:
        logging.error('Failed to converge: {}'.format(e))
        logging.error('error msg {}'.format(e))
        ret_val, wave_fcn = None, None
    except psi4.Dftd3Error as e:
        logging.error('Dftd3Error {}'.format(e))
        ret_val, wave_fcn = None, None
    except psi4.PsiException as e:
        logging.error('Failed to properly execute psi4 {} ({}/{})'.format(calculation, method, basis))
        logging.error('Psi4 Error: {}'.format(e))
        ret_val, wave_fcn = None, None
    except Exception as e:
        logging.error('Unexpected error in psi4: {}'.format(e))
        logging.error("Unexpected error in psi4: %s", sys.exc_info()[0])
        ret_val, wave_fcn = None, None
    else:
        logging.info('Successful psi4 execution {} ({}/{})'.format(calculation, method, basis))

    cleanup_psi4()

    # return and reset stderr, stdout
    # psi4.core.flush_outfile()  # this is probably already handled by the cube
    # logging.info('errorlen {} outlen {}'.format(len(psi_err.getvalue()),len(psi_out.getvalue())))
    # logging.error(psi_err.getvalue())
    # logging.info(psi_out.getvalue())
    # sys.stdout = std_out
    # sys.stderr = std_err
    return ret_val, wave_fcn



def psi4_calculation(conf, dih, conf_name, angle_deg,
                     spe_method,
                     spe_basis,
                     geom_opt_technique,
                     opt_method,
                     opt_basis,
                     geom_maxiter=100,
                     **psi4opts):

    dih_atoms = [x for x in dih.GetAtoms()]
    dih_string = ' '.join(['{}'.format(x.GetIdx()+1) for x in dih_atoms])

    psi4_log_filename = 'psi4_' + conf_name + '.dat'
    psi4.core.set_output_file(psi4_log_filename, False)

    psi_mol = get_psi4_mol_from_conf(conf)

    wf_filename = None
    if geom_opt_technique == 'QM':
        logging.info('Start geometry optimization on conformer %s angle %d',
                                                        conf_name, angle_deg)
        geom_opt_start = time.time()
        logging.debug("Frozen Dihedral: '%s'", dih_string)
        psi4.set_options({'frozen_dihedral': dih_string,
                          'geom_maxiter': geom_maxiter,
                          'dynamic_level': 1})
        psi_mol.update_geometry()
        calc_value, wavefcn = run_psi4('optimize', psi_mol,
                                       opt_method, opt_basis, **psi4opts)
        if calc_value is None:
            logging.warning('Failed to optimize geometry for conformer %s angle %d.',
                                                        conf_name, angle_deg)
            # get logfile to help debug failures
            try:
                psi4.core.flush_outfile()
                with open(psi4_log_filename, 'r') as fptr:
                    psi4_log_data = fptr.read()
                    conf.SetData("PSI4_LOG", psi4_log_data)
            except IOError as e:
                print("Unable to retrieve psi4 log data ", psi4_log_filename, " because of error: ", e)

            raise ValueError('Failed to optimize geometry for conformer %s angle %.1f.',
                                                        conf_name, angle_deg)
        logging.info('Optimized geometry for conformer %s angle %d in %.1f s.',
                        conf_name, angle_deg, time.time()-geom_opt_start)
        calc_value *= psi4.constants.hartree2kcalmol

        # write wavefuction
        wf_filename = '{}_{}.molden'.format('wavefcn_opt', conf_name)
        psi4.molden(wavefcn, wf_filename)

    elif geom_opt_technique == 'MM':
        #TODO: Implement molecular mechanics-based optimization
        pass
    else:
        pass

    if (geom_opt_technique == 'None') or (opt_method != spe_method) or (opt_basis != spe_basis):
        logging.info('Calculate single point energy of conformer %s angle %d',
                                                        conf_name, angle_deg)
        spe_start = time.time()
        calc_value, wavefcn = run_psi4('energy', psi_mol,
                                        spe_method, spe_basis, **psi4opts)
        if calc_value is None:
            logging.error('Failed to calculate single point energy for conformer %s angle %d.',
                          conf_name, angle_deg)
            raise ValueError('Failed to calculate single point energy for conformer %s angle %d.',
                             conf_name, angle_deg)
        logging.info('Calculated single point energy for conformer %s angle %d in %.1f s.',
                     conf_name, angle_deg, time.time()-spe_start)
        calc_value *= psi4.constants.hartree2kcalmol

    # attached psi4 log to each file
    try:
        # extract minimal geometry optimization data for success
        psi4.core.flush_outfile()
        with open(psi4_log_filename, 'r') as fptr:
            psi4_log_str = str()
            for line in fptr.readlines():
                if line.find('~') > 0 or line.find('Exception') > 0 or line.find('Optimization') > 0:
                    psi4_log_str += line
            conf.SetData("PSI4_LOG", psi4_log_str)
    except IOError as e:
        print("Unable to save psi4 log data ", psi4_log_filename, " because of error: ", e)

    return calc_value, psi_mol, wavefcn

import datetime
import tempfile
import os
import socket
import traceback
from sh import which

from openeye import oechem

import psi4

from floe import constants
from floe.api import parameter
from floe.api import OEMolComputeCube, MoleculeOutputPort, ParallelOEMolComputeCube, ParallelMixin

from torsion.utils import get_sd_data, write_energy_profile_to_sddata, save_sddata
from torsion.core import get_dihedral, calculate_energy


class SerialPsi4EnergyCalculation(OEMolComputeCube):
    """Calculate the energy using psi4 for a single conformer.
    """
    title = 'PSI4 Energy Calculation (Serial)'
    description = """Flexible 3rd party QM engine.  Allows a variety of calculations, methods and basis sets using
                     the trusted DFT code from Dave Sherill's group at Georgia Tech."""
    classification = [['Energetics', 'DFT', 'PSI4']]
    tags = [tag for tag_list in classification for tag in tag_list]

    success = MoleculeOutputPort('success')
    failure = MoleculeOutputPort('failure')
    system_failure = MoleculeOutputPort('system_failure')

    spe_method = parameter.StringParameter(
        'spe_method',
        title='Method for QM Single-point Energy Calculation',
        default='B3LYP',
        choices=['SCF', 'B3LYP', 'B3LYP-D', 'B2PLYP', 'B97', 'B97-D', 'PBE-D', 'CCSD', 'SAPT0',
                 'MP2', 'MP4', 'HF', 'hf3c', 'hf-d3bj', 'pbeh3c'],
        description='QM method for the final single-point calculation')

    spe_basis = parameter.StringParameter(
        'spe_basis',
        title='Basis Set for QM Single-point Energy Calculation',
        default='6-31G**',
        choices=[None, 'cc-pVDZ', 'cc-pVTZ', 'cc-pVQZ', 'aug-cc-pVDZ', 'aug-cc-pVTZ', 'minix', 'sto-3g',
                 '3-21G', '6-31G', '6-31G*', '6-31G**', '6-31+G*', '6-31+G**', '6-31++G*', '6-31++G**', '6-311+G**'],
        description='The basis set for the final single-point calculation')

    geom_opt_technique = parameter.StringParameter(
        'geom_opt_technique',
        title='Type of Geometry Optimization',
        default='QM',
        choices=['None', 'QM'])

    opt_method = parameter.StringParameter(
        'opt_method',
        title='Method for QM Geometry Optimization',
        default='B3LYP',
        choices=['SCF', 'B3LYP', 'B3LYP-D', 'B97', 'B97-D', 'PBE-D', 'CCSD', 'SAPT0',
                 'MP2', 'MP4', 'HF', 'hf3c', 'hf-d3bj', 'pbeh3c'],
        description='QM method for the initial optimization')

    opt_basis = parameter.StringParameter(
        'opt_basis',
        title='Basis Set for QM Geometry Optimization',
        default='6-31G*',
        choices=[None, 'cc-pVDZ', 'cc-pVTZ', 'cc-pVQZ', 'aug-cc-pVDZ', 'aug-cc-pVTZ', 'minix', 'sto-3g',
                 '3-21G', '6-31G', '6-31G*', '6-31G**', '6-31+G*', '6-31+G**', '6-31++G*', '6-31++G**', '6-311+G**'],
        description='The basis set for the initial optimization.')

    geom_maxiter = parameter.IntegerParameter(
        'geom_maxiter',
        title='Maximum number of geometry optimization steps.',
        default=200,
        min_value=1,
        max_value=1000,
        description="""Maximum number of geometry optimization steps.""")

    dft_radial_points = parameter.IntegerParameter(
        'dft_radial_points',
        title='DFT Grid Radial Points',
        default=50,
        min_value=25,
        max_value=500,
        level=constants.ADVANCED,
        description="""Parameter controls the radial density of the DFT grid.  DFT calculation speeds are very
                        sensitive to this parameter.""")

    dft_spherical_points = parameter.IntegerParameter(
        'dft_spherical_points',
        title='DFT GRid Spherical Points',
        default=194,
        min_value=25,
        max_value=500,
        description="""Parameter controls the spherical density of the DFT grid.  DFT calculation speeds are sensitive
                        to this parameter.""",
        level=constants.ADVANCED)

    num_processors = parameter.IntegerParameter(
        'num_processors',
        title='Number of Processors',
        default=1,
        min_value=0,
        max_value=32,
        level=constants.ADVANCED,
        help_text="""Number of processors 1-32 (0 indicates to use all the processors on the machine).  This is the
                      Number of processors for each Orion worker instance.  Unless your Orion worker instances have
                      been setup to avoid sharing resources, the should remain set at 1.""")

    guess_basis = parameter.BooleanParameter(
        'guess_basis',
        title='Basis set guess',
        level=constants.ADVANCED,
        required=False,
        default=False,
        help_text="""Psi4 advanced parameter:  Accelerate convergence by performing a preliminary scf with this small 
        basis set followed by projection into the full target basis. A value of TRUE turns on projection using the 
        3-21G small basis set.
        http://www.psicode.org/psi4manual/master/autodir_options_c/scf__basis_guess.html""")

    use_soscf = parameter.BooleanParameter(
        'use_soscf',
        title='Use Second-Order SCF',
        level=constants.ADVANCED,
        required=False,
        default=False,
        help_text="""Psi4 advance parameter: Do use second-order SCF convergence methods?
        http://www.psicode.org/psi4manual/master/autodir_options_c/scf__soscf.html""")

    scf_type = parameter.StringParameter(
        'scf_type',
        title='SCF Type',
        level=constants.ADVANCED,
        required=False,
        default='DIRECT',
        choices=['DIRECT', 'DF', 'PK', 'OUT_OF_CORE', 'PS', 'INDEPENDENT', 'GTFOCK'],
        help_text="""Psi4 parameter: SCF Type.
        """
    )

    only_selected_conformer = parameter.BooleanParameter(
        'only_selected_conformer',
        title='Calculate Energy only for a selected conformer',
        default=False,
        help_text="""If this is set, energy is calculated only for a single 
        conformer of each OEMol passed as input to this cube. The conformer for
        which the energy is calculated is determined by the 'SELECTED_CONFORMER'
        integer data tag on the OEMol.""")

    molden_output = parameter.BooleanParameter(
        'molden_output',
        title='Attach electronic wave function from molden file as SD data',
        default=False,
        help_text="""If this is set, electronic wave function from
        molden file will be attached as SD data.""")

    g_convergence = parameter.StringParameter(
        'g_convergence',
        title='Psi4 g_convergence parameter',
        default='QCHEM',
        choices=['QCHEM', 'MOLPRO', 'GAU', 'GAU_LOOSE', 'GAU_TIGHT', 'INTERFRAG_TIGHT', 'GAU_VERYTIGHT', 'TURBOMOLE', 'CFOUR', 'NWCHEM_LOOSE'],
        help_text="""Allows selection of a psi4 convergence criteria.  See: http://www.psicode.org/psi4manual/master/autodoc_glossary_options_c.html#term-g-convergence-optking""")

    max_disp_g_convergence = parameter.DecimalParameter(
        'max_disp_g_convergence',
        title='Psi4 max_disp_g_convergence parameter',
        default=1.2e-3,
        help_text="""Psi4 Maximum displacement convergence criteria.  NOTE: For loose optimization, try 5.0e-2.""")

    def begin(self):
        psi_path = which('psi4')
        psi_path, tail = os.path.split(psi_path)
        psi_path, tail = os.path.split(psi_path)
        os.environ['PSI'] = psi_path
        conda_path = '/'.join(psi_path.split('/')[:-1])
        os.environ['PSIPATH'] = psi_path

        os.environ['PSI_SCRATCH'] = tempfile.gettempdir()
        os.environ['PSI_SCRATCH_LOCAL'] = tempfile.mkdtemp()

        psi4.set_memory('1 GB')

        self.basis_guess_str = 'false'
        if self.args.guess_basis and self.args.spe_basis not in ['sto-3g', '3-21G', 'minix'] and self.args.opt_basis not in ['3-21G', 'sto-3g', 'minix']:
            self.basis_guess_str = 'true'

        self.use_soscf_str = 'false'
        if self.args.use_soscf and self.args.spe_basis not in ['sto-3g', '3-21G', 'minix'] and self.args.opt_basis not in ['3-21G', 'sto-3g', 'minix']:
            self.use_soscf_str = 'true'

        self.psi4opts = {
                'scf_type' : self.args.scf_type,
                'fail_on_maxiter' : 'false',
                'guess_basis' : self.basis_guess_str,
                'use_soscf' : self.use_soscf_str,
                'dft_radial_points': self.args.dft_radial_points,
                'dft_spherical_points': self.args.dft_spherical_points,
                'num_processors': self.args.num_processors,
                'g_convergence': self.args.g_convergence,
                'max_disp_g_convergence': self.args.max_disp_g_convergence
                }

    def end(self):
        if os.path.exists(os.environ['PSI_SCRATCH_LOCAL']):
            os.rmdir(os.environ['PSI_SCRATCH_LOCAL'])

    def process(self, mol, port):
        parent_torsion_tag = 'TORSION_ATOMS_ParentMol'
        torsion_atoms_in_parent = get_sd_data(mol, parent_torsion_tag).split()
        dih_name = mol.GetTitle()+ '_' + '_'.join(torsion_atoms_in_parent)

        torsion_tag = 'TORSION_ATOMS_FRAGMENT'
        torsion_atoms_in_fragment = get_sd_data(mol, torsion_tag).split()
        dihedral_atom_indices = [int(x)-1 for x in torsion_atoms_in_fragment]
        if dihedral_atom_indices is None:
            self.log.warning('Unable to find labelled torsion in %s' % dih_name)
            self.failure.emit(mol)
            return

        try:
            if self.args.only_selected_conformer:
                conf_selection_tag = 'SELECTED_CONFORMER'
                key_conf_id = mol.GetIntData(conf_selection_tag)
                for conf in mol.GetConfs():
                    if conf.GetIdx() != key_conf_id:
                        continue
                conf_name = get_sd_data(conf, 'CONFORMER_LABEL')
            else:
                conf_name = get_sd_data(mol, 'CONFORMER_LABEL')
            time_stamp = "{:%Y-%m-%d %H:%M:%S}".format(datetime.datetime.now())
            hostname = socket.gethostname()
            self.log.info("Starting psi4 calculation for %s on %s at %s" % (conf_name, hostname, time_stamp))

            if self.args.only_selected_conformer:
                oechem.OESetSDData(conf, '%s start time' % self.name, time_stamp)
            else:
                oechem.OESetSDData(mol, '%s start time' % self.name, time_stamp)


            dih, _ = get_dihedral(mol, dihedral_atom_indices)
            calculate_energy(mol, dih,
                             spe_method=self.args.spe_method,
                             spe_basis=self.args.spe_basis,
                             geom_opt_technique=self.args.geom_opt_technique,
                             opt_method=self.args.opt_method,
                             opt_basis=self.args.opt_basis,
                             geom_maxiter=self.args.geom_maxiter,
                             only_selected_conf=self.args.only_selected_conformer,
                             molden_output=self.args.molden_output,
                             **self.psi4opts)

            if self.args.only_selected_conformer:
                conf_selection_tag = 'SELECTED_CONFORMER'
                key_conf_id = mol.GetIntData(conf_selection_tag)
                for conf in mol.GetConfs():
                    if conf.GetIdx() != key_conf_id:
                        continue
                conf_name = get_sd_data(conf, 'CONFORMER_LABEL')
            else:
                conf_name = get_sd_data(mol, 'CONFORMER_LABEL')
            time_stamp = "{:%Y-%m-%d %H:%M:%S}".format(datetime.datetime.now())
            hostname = socket.gethostname()
            self.log.info("Completed psi4 calculation for %s on %s at %s" % (conf_name, hostname, time_stamp))

            if self.args.only_selected_conformer:
                oechem.OESetSDData(conf, '%s end time' % self.name, time_stamp)
            else:
                oechem.OESetSDData(mol, '%s end time' % self.name, time_stamp)

            self.success.emit(mol)
        except Exception as e:
            print(e)
#            traceback.print_stack()
            self.log.error("Error with {} {}".format(mol.GetTitle(), e))
            self.failure.emit(mol)


class ParallelPsi4EnergyCalculation(ParallelMixin, SerialPsi4EnergyCalculation):
    title = 'Calculate Psi4 Energy (Parallel)'
    classification = [["Energetics", "Psi4", "DFT"]]

    parameter_overrides = {
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_timeout": {"default": 43200.0},  # Default 12 hour limit (units are seconds)
        "item_count": {"default": 1},  # 1 molecule at a time
        "max_failures": {"default": 1}, # only 1 failure permitted

        #"spe_method": {"hidden": False},
        #"spe_basis": {"hidden": False},
        #"opt_method": {"hidden": False},
        #"opt_basis": {"hidden": False},
        #"geom_maxiter": {"hidden": False},
        #"molden_output": {"hidden": False},

        #"dft_radial_points": {"hidden": True},
        #"dft_spherical_points": {"hidden": True},
        #"guess_basis": {"hidden": True},
        #"geom_opt_technique": {"hidden": True},
        #"g_convergence": {"hidden": True},
        #"max_disp_g_convergence": {"hidden": True},
        #"num_processors": {"hidden": True},
        #"only_selected_conformer": {"hidden": True},
        #"scf_type": {"hidden": True},
        #"use_soscf": {"hidden": True},
    }

    def process_failed(self, data, port, last_error):
        print("Parallel cube failed to process {} from {} with error: {}".format(data, port, last_error))
        self.system_failure.emit(data)


class HiddenParamParallelPsi4EnergyCalculation(ParallelMixin, SerialPsi4EnergyCalculation):
    """
    Class to mimic Psi4 cube but with no exposed parameters
    """
    title = 'Hidden parameter - Calculate Psi4 Energy (Parallel)'
    classification = [["Energetics", "Psi4", "DFT"]]

    parameter_overrides = {
        "prefetch_count": {"default": 1},  # 1 molecule at a time
        "item_timeout": {"default": 43200.0},  # Default 12 hour limit (units are seconds)
        "item_count": {"default": 1},  # 1 molecule at a time
        "max_failures": {"default": 1}, # only 1 failure permitted
        "num_processors": {"hidden": True},
        "dft_radial_points": {"hidden": True},
        "dft_spherical_points": {"hidden": True},
        "g_convergence": {"hidden": True},
        "max_disp_g_convergence": {"hidden": True},
        "guess_basis": {"hidden": True},
        "use_soscf": {"hidden": True},
        "scf_type": {"hidden": True},
        "spe_method": {"hidden": True},
        "spe_basis": {"hidden": True},
        "opt_method": {"hidden": True},
        "opt_basis": {"hidden": True},
        "geom_maxiter": {"hidden": True},
        "only_selected_conformer": {"hidden": True},
        "geom_opt_technique": {"hidden": True},
        "molden_output": {"hidden": True},
    }

    def process_failed(self, data, port, last_error):
        print("Parallel cube failed to process {} from {} with error: {}".format(data, port, last_error))
        self.system_failure.emit(data)



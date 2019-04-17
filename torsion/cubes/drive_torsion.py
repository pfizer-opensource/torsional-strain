import os
import datetime
from floe import constants
from floe.api import parameter
from floe.api import MoleculeOutputPort
from floe.api import OEMolComputeCube, ParallelOEMolComputeCube, Cube

from openeye import oechem
from torsion.core import get_dihedral
from torsion.utils import get_sd_data

class DriveTorsion(Cube):
    """
        Drive the primary torsion.
    """

    cube_type = constants.CUBE_COMPUTE

    num_points = parameter.IntegerParameter(
        'num_points',
        title='Number of points at which to sample the dihedral.',
        default=24,
        min_value=4,
        max_value=36,
        description="""The number of evenly spaced torsion angles to sample in
        order to determine the enthalpy surface.""")

    to_energy_calc = MoleculeOutputPort('to_energy_calc')

    conf_selection_tag = 'SELECTED_CONFORMER'

    def process(self, mol, port):
        """
            The input to this cube will be an OEMol with one or more conformers
            with "CONFORMER_LABEL" SD Data of the form 'XY-1234567_1_2_3_4_00_00'
        """
        num_confs = mol.NumConfs()

        last_conf = mol.GetActive()
        last_conf_name = oechem.OEGetSDData(last_conf, 'CONFORMER_LABEL')
        self.log.info("Processing conformer {} on {} at {:%Y-%m-%d %H:%M:%S}".format(last_conf_name, os.environ['HOSTNAME'], datetime.datetime.now()))

        if num_confs == self.args.num_points:
            self.success.emit(mol)
            self.log.info("Completed scan for {} on {} at {:%Y-%m-%d %H:%M:%S}".format(mol.GetTitle(), os.environ['HOSTNAME'], datetime.datetime.now()))
            return

        if num_confs == 1 and not mol.HasData(self.conf_selection_tag):
            self.log.info("Conformer {} is a fresh starting conformer on {} at {:%Y-%m-%d %H:%M:%S}".format(mol.GetTitle(), os.environ['HOSTNAME'], datetime.datetime.now()))
            mol.SetIntData(self.conf_selection_tag, last_conf.GetIdx())
            last_conf.SetDoubleData('TORSION_ANGLE', 0.0)
            oechem.OESetSDData(last_conf, 'TORSION_ANGLE', '0.0')
            self.log.info("Sending conformer {} to energy calculation from {} at {:%Y-%m-%d %H:%M:%S}".format(last_conf_name, os.environ['HOSTNAME'], datetime.datetime.now()))
            self.to_energy_calc.emit(mol)
            return

        try:
            torsion_tag = 'TORSION_ATOMS_FRAGMENT'
            torsion_atoms_in_fragment = get_sd_data(mol, torsion_tag).split()
            dihedral_atom_indices = [int(x)-1 for x in torsion_atoms_in_fragment]

            dih, _ = get_dihedral(mol, dihedral_atom_indices)
            dih_atoms = [x for x in dih.GetAtoms()]

            # if the last energy calculation failed
            if not oechem.OEHasSDData(last_conf, 'PSI4_ENERGY'):
                self.log.info("Conformer {} found to have NO ENERGY on {} at {:%Y-%m-%d %H:%M:%S}".format(last_conf_name, os.environ['HOSTNAME'], datetime.datetime.now()))
                mol.PopActive()
                last_conf = mol.GetActive()

            new_conf = mol.NewConf(last_conf)
            mol.PushActive(new_conf)
            conf_no = num_confs
            conformer_label = last_conf_name[:-3] +\
                             '_{:02d}'.format(conf_no)
            oechem.OESetSDData(new_conf, "CONFORMER_LABEL", conformer_label)

            angle = num_confs*2*oechem.Pi/self.args.num_points
            angle_deg = oechem.Rad2Deg*angle
            new_conf.SetDoubleData('TORSION_ANGLE', angle_deg)
            oechem.OESetSDData(new_conf, 'TORSION_ANGLE', '{:.1f}'.format(angle_deg))

            if not oechem.OESetTorsion(new_conf,
                    dih_atoms[0], dih_atoms[1], dih_atoms[2], dih_atoms[3],
                    angle):
                self.log.error("Could not rotate conformer {} by {:.1f} on {} at {:%Y-%m-%d %H:%M:%S}".format(last_conf_name, angle_deg, os.environ['HOSTNAME'], datetime.datetime.now()))

            mol.SetIntData(self.conf_selection_tag, new_conf.GetIdx())
            self.log.info("Sending conformer {} to energy calculation from {} at {:%Y-%m-%d %H:%M:%S}".format(conformer_label, os.environ['HOSTNAME'], datetime.datetime.now()))
            self.to_energy_calc.emit(mol)

        except Exception as e:
            self.log.error("COuld not drive torsion in  conformer {} on {} at {:%Y-%m-%d %H:%M:%S}: {}".format(last_conf_name, os.environ['HOSTNAME'], datetime.datetime.now(), e))
            self.failure.emit(mol)

class SerialDriveTorsion(DriveTorsion, OEMolComputeCube):
    """
        Drive the primary torsion.
    """
    cube_type = constants.CUBE_COMPUTE
    title = "Drive the Primary Torsion"
    classification = [["Cheminformatics"]]

class ParallelDriveTorsion(DriveTorsion, ParallelOEMolComputeCube):
    """
        Drive the primary torsion.
    """
    cube_type = constants.CUBE_COMPUTE_PARALLEL
    title = "Drive the Primary Torsion (Parallel)"
    classification = [["Cheminformatics"]]

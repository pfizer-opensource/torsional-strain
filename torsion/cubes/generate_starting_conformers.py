from openeye import oechem

from datarecord import OEMolRecord
from floe.api import ParallelMixin, IntegerParameter, DecimalParameter, BooleanParameter
from cuberecord.cubes import OEMolRecordCube, InOutMolFieldMixin

from torsion.conf import gen_starting_confs, split_confs, torsion_library
from torsion.utils import get_sd_data


class GenerateStartingConfs(OEMolRecordCube, InOutMolFieldMixin):
    """ Generate starting conformations for torsion scan
    """
    title = "Generate Multiple Starting Conformations"
    classification = [["Cheminformatics"]]

    max_confs = IntegerParameter(
        'max_confs',
        title='Maximum Conformers',
        required=True,
        default=20,
        help_text='Maximum number of starting conformations to use in QM torsion driving experiment. '+
        'NOTE: If you do not want to generate starting conformations, set this to 1 (not recommended).')

    rms_cutoff =DecimalParameter(
        'rms_cutoff',
        title='RMSD Cutoff',
        required=True,
        default=0.0,
        help_text="""Minimum RMSD between any two conformers in the final ensemble.
        NOTE: A value of 0.0 permits methyl rotors to be sampled.""")

    energy_window = DecimalParameter(
        'energy_window',
        title='Energy Window',
        required=True,
        default=25,
        help_text="""Acceptable energy difference between alternate starting
        conformers.""")

    split_output = BooleanParameter(
            'split_output',
            title='Split Output',
            required=False,
            default=False,
            help_text="""Emit each starting conformer as its own OEMol.
            """)

    def begin(self):
        self.torsion_library = torsion_library

    def process(self, record, port):
        if record.has_value(self.args.in_mol_field):
            mol = record.get_value(self.args.in_mol_field)
            fragmentLabel = mol.GetTitle() + '_' + '_'.join(get_sd_data(mol, 'TORSION_ATOMS_ParentMol').split())
            try:
                starting_conformers = gen_starting_confs(mol, self.torsion_library,
                                                         max_one_bond_away=True,
                                                         num_conformers=self.args.max_confs,
                                                         rms_cutoff=self.args.rms_cutoff,
                                                         energy_window=self.args.energy_window)
                if self.args.split_output:
                    for pose in split_confs(starting_conformers):
                        pose_record = OEMolRecord()
                        pose_record.set_mol(pose)
                        self.success.emit(pose_record)
                else:
                    multi_conf_record = OEMolRecord()
                    multi_conf_record.set_mol(starting_conformers)
                    self.success.emit(multi_conf_record)
                    
                self.log.info('%d starting conformers generated for fragment %s.' % ( starting_conformers.NumConfs(), fragmentLabel))

            except Exception as e:
                self.log.error("Could not generate conformers for fragment %s: %s" % (fragmentLabel, e))
                self.failure.emit(record)
            
        else:
            self.log.error("Could not find molecules in OEMolRecord!")
            self.failure.emit(record)

            
class ParallelGenerateStartingConfs(ParallelMixin, GenerateStartingConfs):
    """
        Sample select conformations of amide, methyl, and hydroxyl rotors upto
        two bonds away from the dihedral of interest.
    """
    title = "Generate Multiple Starting Conformations (Parallel)"

    parameter_overrides = {
        "prefetch_count": {"default": 10},  # 10 molecules at a time
        "item_timeout": {"default": 100.0},  # (units are seconds)
        "item_count": {"default": 1},  # 1 molecule at a time
        "max_failures": {"default": 1}, # only 1 failure permitted
    }


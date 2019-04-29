from openeye import oechem

from datarecord import OEMolRecord
from floe.api import ParallelMixin, IntegerParameter, BooleanParameter
from cuberecord.cubes import OEMolRecordCube, InOutMolFieldMixin

from torsion.core import get_dihedral
from torsion.conf import gen_torsional_confs, split_confs, torsion_library, get_best_conf
from torsion.utils import get_sd_data


class GenerateTorsionalConfs(OEMolRecordCube, InOutMolFieldMixin):
    """ Generate conformers by rotating the primary torsion.
    """
    title = "Generate Multiple Torsional Conformations"
    classification = [["Cheminformatics"]]

    num_points = IntegerParameter(
        'num_points',
        title='Number of torsional conformers to generate.',
        default=24,
        min_value=1,
        max_value=36,
        description="""The number of evenly spaced torsion angles to sample 
        when generating torsional conformers.""")

    split_confs= BooleanParameter(
        'split_confs',
        title='Emit each conformer separately',
        default=True,
        description="""Whether conformers should be emitted separately or as part of a single molecule.""")

    best_conf= BooleanParameter(
        'best_conf',
        title='For each torsion select single best conformer',
        default=True,
        description="""Whether single best conformer should be emitted for each dihedral angle.""")

    def begin(self):
        self.torsion_library = torsion_library

    def process(self, record, port):
        if record.has_value(self.args.in_mol_field):
            mol = record.get_value(self.args.in_mol_field)
            fragmentLabel = mol.GetTitle() + '_' + \
                             '_'.join(get_sd_data(mol, 'TORSION_ATOMS_ParentMol').split())

            torsion_tag = 'TORSION_ATOMS_FRAGMENT'
            torsion_atoms_in_fragment = get_sd_data(mol, torsion_tag).split()
            dihedral_atom_indices = [int(x)-1 for x in torsion_atoms_in_fragment]

            try:
                dih, _  = get_dihedral(mol, dihedral_atom_indices)
                if self.args.best_conf:
                    torsional_conformers = get_best_conf(mol, dih, self.args.num_points)
                else:
                    torsional_conformers = gen_torsional_confs(mol, dih,
                                                               self.args.num_points,
                                                               include_input=False)

                if self.args.split_confs:
                    for pose in split_confs(torsional_conformers):
                        record = OEMolRecord()
                        record.set_mol(pose)
                        self.success.emit(record)
                else:
                    record = OEMolRecord()
                    record.set_mol(torsional_conformers)
                    self.success.emit(torsional_conformers)

                self.log.info('%d torsional conformers generated for fragment %s.' % (
                                    torsional_conformers.NumConfs(), fragmentLabel))

            except Exception as e:
                self.log.error("Could not generate conformers for fragment %s: %s" % (fragmentLabel, e))
                self.failure.emit(record)
            
        else:
            self.log.error("Could not find molecules in OEMolRecord!")
            self.failure.emit(record)
            


class ParallelGenerateTorsionalConfs(ParallelMixin, GenerateTorsionalConfs):
    """
        Generate conformers by rotating the primary torsion.
    """
    title = "Generate Multiple Torsional Conformations (Parallel)"

    parameter_overrides = {
        "prefetch_count": {"default": 10},  # 10 molecules at a time
        "item_timeout": {"default": 100.0},  # (units are seconds)
        "item_count": {"default": 1},  # 1 molecule at a time
        "max_failures": {"default": 1}, # only 1 failure permitted
    }


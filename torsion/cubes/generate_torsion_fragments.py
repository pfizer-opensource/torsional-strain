from openeye import oechem

from datarecord import OEMolRecord
from floe.api import ParallelMixin
from cuberecord.cubes import OEMolRecordCube, InOutMolFieldMixin

from torsion.utils import gen_torsion_fragments


class GenerateFragments(OEMolRecordCube, InOutMolFieldMixin):
    """ Generate torsion fragments
    """
    title = "Generate Torsion Fragments"
    classification = [["Cheminformatics"]]

    def process(self, record, port):
        if record.has_value(self.args.in_mol_field):
            mol = record.get_value(self.args.in_mol_field)
            try:
                fragments = gen_torsion_fragments(mol)
                if len(fragments) == 0:
                    raise('Unable to generate torsion fragments')

                for fragment in fragments:
                    fragment_record = OEMolRecord()
                    fragment_record.set_mol(fragment)
                    self.success.emit(fragment_record)
                    
                self.log.info('%d torsion fragments generated for molecule %s.' % (
                                    len(fragments), mol.GetTitle()))
            except Exception as e:
                self.log.error("Could not generate torsion fragments %s" % mol.GetTitle())
                self.failure.emit(mol)
        else:
            self.failure.emit(record)


class ParallelGenerateFragments(GenerateFragments, ParallelMixin):
    """ Generate torsion fragments
    """

    parameter_overrides = {
        "prefetch_count": {"default": 10},  # 10 molecules at a time
        "item_timeout": {"default": 100.0},  # (units are seconds)
        "item_count": {"default": 1},  # 1 molecule at a time
        "max_failures": {"default": 1}, # only 1 failure permitted
    }


from floe import constants
from floe.api import parameter
from floe.api import OEMolComputeCube, ParallelOEMolComputeCube, Cube

from openeye import oechem
from torsion.utils import get_sd_data, gen_torsion_fragments


class GenerateFragments(Cube):
    """
        Generate torsion fragments
    """

    cube_type = constants.CUBE_COMPUTE

    def begin(self):
        pass

    def process(self, mol, port):
        try:
            fragments = gen_torsion_fragments(mol)
            if len(fragments) == 0:
                raise('Unable to generate torsion fragments')

            for fragment in fragments:
                self.success.emit(fragment)

            self.log.info('%d torsion fragments generated for molecule %s.' % (
                                len(fragments), mol.GetTitle()))

        except Exception as e:
            self.log.error("Could not generate torsion fragments %s" % mol.GetTitle())
            self.failure.emit(mol)

class SerialGenerateFragments(GenerateFragments, OEMolComputeCube):
    """
        Generate torsion fragments
    """
    cube_type = constants.CUBE_COMPUTE
    title = "Generate Torsion Fragments"
    classification = [["Cheminformatics"]]

class ParallelGenerateFragments(GenerateFragments, ParallelOEMolComputeCube):
    """
        Generate torsion fragments
    """
    cube_type = constants.CUBE_COMPUTE_PARALLEL
    title = "Generate Torsion Fragments (Parallel)"
    classification = [["Cheminformatics"]]

    parameter_overrides = {
        "prefetch_count": {"default": 10},  # 10 molecules at a time
        "item_timeout": {"default": 100.0},  # (units are seconds)
        "item_count": {"default": 1},  # 1 molecule at a time
        "max_failures": {"default": 1}, # only 1 failure permitted
    }


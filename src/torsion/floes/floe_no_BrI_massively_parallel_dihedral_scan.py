import logging
import socket
print('host:', socket.gethostname())
logging.basicConfig(level=logging.DEBUG, 
        format='%(asctime)s' + socket.gethostname() + '%(levelname)-8s  %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S')
from torsion.cubes import ParallelPsi4EnergyCalculation, HiddenParamParallelPsi4EnergyCalculation
from torsion.cubes import ParallelGenerateStartingConfs, ParallelGenerateTorsionalConfs
from torsion.cubes import FilterBrI
from floe.api import WorkFloe
from floe.api import ParallelCubeGroup
from floe.api import OEMolOStreamCube
from floe.api import OEMolIStreamCube
from floe import constants

# Declare Floe, add metadata for UI
job = WorkFloe('Filtered Massively Parallel 2-stage Dihedral Scan of Torsions')
job.description = """
Filtered Calculate the Energy Profile associated with a Dihedral Scan in a massively parallel manner.
"""
job.classification = [["Torsion"]]

# Declare Cubes
ifs = OEMolIStreamCube('ifs')
ifs.promote_parameter('data_in', promoted_name='ifs')
ifs.parameter_overrides["download_format"] = {"hidden": True}
ifs.parameter_overrides["limit"] = {"hidden":True}

element_filter = FilterBrI('element_filter')

# Conformer generation
confgenCube = ParallelGenerateStartingConfs('confgenCube')
confgenCube.promote_parameter('max_confs', required=True)
confgenCube.promote_parameter('rms_cutoff', required=True)
confgenCube.parameter_overrides["energy_window"] = {"hidden": True}
confgenCube.parameter_overrides["split_output"] = {"hidden": True}

confgen_failure = OEMolOStreamCube('confgen_failure')
confgen_failure.promote_parameter('data_out', promoted_name='confgen_failure',
                                  title='Conformer Generation Failures',
                                  description='Conformer Generation Failures',
                                  default='confgen_failures')
confgen_failure.promote_parameter('buffered', default=False)
confgen_failure.parameter_overrides["buffered"] = {"hidden": True}

# Search torsions & select best conformer for each
torsgenCube = ParallelGenerateTorsionalConfs('torsional_conf_gen')
torsgenCube.promote_parameter('num_points',
                              required=True,
                              title='Number of torsional conformers to generate.',
                              description="""The number of evenly spaced torsion angles to sample 
                              when generating torsional conformers.""")
torsgenCube.parameter_overrides["split_confs"] = {"hidden": True}
torsgenCube.parameter_overrides["best_conf"] = {"hidden": True}

# Conformer/torsion cube group
group = ParallelCubeGroup(cubes=[confgenCube, torsgenCube])
job.add_group(group)

# Fast QM optimization with hf3c to get approximate geometry
hf3cCube = HiddenParamParallelPsi4EnergyCalculation('hf3c_precalc')
hf3cCube.set_parameters(spe_method='hf3c', spe_basis='minix', opt_method='hf3c', opt_basis='minix', geom_maxiter=200)


# Full QM optimization & SP for energy
psi4EnergyCube = ParallelPsi4EnergyCalculation('parallel_psi4_energy_calculation')
psi4EnergyCube.set_parameters(geom_maxiter=200)
psi4EnergyCube.promote_parameter('spe_method', required=True)
psi4EnergyCube.promote_parameter('spe_basis', required=True)
psi4EnergyCube.promote_parameter('molden_output', required=True)
psi4EnergyCube.promote_parameter('opt_method', required=True)
psi4EnergyCube.promote_parameter('opt_basis', required=True)

# hf3c failure handling
failhf3c = OEMolOStreamCube('failhf3c')
failhf3c.promote_parameter('data_out', promoted_name='hf3cfail', title='hf3cfail',
                          description="hf3c Failures",
                          default='hf3c_failures')
failhf3c.promote_parameter('buffered', default=False)
failhf3c.parameter_overrides["buffered"] = {"hidden": True}


# torsion driving failures
failfs1 = OEMolOStreamCube('failfs1')
failfs1.promote_parameter('data_out', promoted_name='torfail', title='torfail',
                          description="Torsional Conformer Generation Failures",
                          default='torgen_failures')
failfs1.promote_parameter('buffered', default=False)
failfs1.parameter_overrides["buffered"] = {"hidden": True}


# final molecular output
ofs = OEMolOStreamCube('ofs')
ofs.promote_parameter('data_out', promoted_name='ofs',
                      description='Floe output',
                      default='output',
                      title='Successes')
ofs.promote_parameter('buffered', default=False)
ofs.parameter_overrides["buffered"] = {"hidden": True}

# psi4 cube failures
failfs2 = OEMolOStreamCube('failfs2')
failfs2.promote_parameter('data_out', promoted_name='psifail', title='psifail',
                          description="Psi4 Failures",
                          default='psi4_failures')
failfs2.promote_parameter('buffered', default=False)
failfs2.parameter_overrides["buffered"] = {"hidden": True}

# psi4 system failures
failfs3 = OEMolOStreamCube('failfs3')
failfs3.promote_parameter('data_out', promoted_name='sysfail', title='sysfail',
                          description="Psi4 System Failures",
                          default='psi4_sysfailures')
failfs3.promote_parameter('buffered', default=False)
failfs3.parameter_overrides["buffered"] = {"hidden": True}

# Bromine/Iodine filtering failures
failBrI = OEMolOStreamCube('failBrI')
failBrI.promote_parameter('data_out', title='fail_BrI', default='fail_BrI')
failBrI.promote_parameter('buffered', default=False)
failBrI.parameter_overrides["buffered"] = {"hidden": True}

# Add Cubes to Floe
cubes = [ifs, element_filter, confgenCube, torsgenCube, hf3cCube, psi4EnergyCube, ofs, 
         failhf3c, failfs1, failfs2, failfs3, failBrI, confgen_failure]
[job.add_cube(c) for c in cubes]

# Connect ports
ifs.success.connect(element_filter.intake)
element_filter.success.connect(confgenCube.intake)
element_filter.failure.connect(failBrI.intake)
confgenCube.success.connect(torsgenCube.intake)
confgenCube.failure.connect(confgen_failure.intake)
torsgenCube.success.connect(hf3cCube.intake)
torsgenCube.failure.connect(failfs1.intake)
hf3cCube.success.connect(psi4EnergyCube.intake)
hf3cCube.failure.connect(failhf3c.intake)
hf3cCube.system_failure.connect(psi4EnergyCube.intake)
psi4EnergyCube.success.connect(ofs.intake)
psi4EnergyCube.failure.connect(failfs2.intake)
psi4EnergyCube.system_failure.connect(failfs3.intake)

# If called from command line, run the floe
if __name__ == "__main__":
    job.run()

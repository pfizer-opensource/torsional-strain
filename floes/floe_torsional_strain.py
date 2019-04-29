import os
import logging
logging.basicConfig(level=logging.DEBUG,
        format='%(asctime)s' + os.getenv('HOSTNAME') + '%(levelname)-8s  %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S')
from floe.api import WorkFloe
from cuberecord import DatasetReaderCube, DatasetWriterCube
from torsion.cubes import ParallelPsi4EnergyCalculation
from torsion.cubes import GenerateFragments, GenerateStartingConfs, GenerateTorsionalConfs

# Declare Floe, add metadata for UI
job = WorkFloe('Torsional Strain Energy')
job.description = """
Calculate torsional strain energy.
"""
job.classification = [["Torsional Strain"]]

# Declare cubes
ifs = DatasetReaderCube('ifs')
fraggenCube = GenerateFragments('fragment_generation')
fraggen_failure = DatasetWriterCube('fraggen_failure')
confgenCube = GenerateStartingConfs('starting_conf_gen')
confgen_failure = DatasetWriterCube('confgen_failure')
torsgenCube = GenerateTorsionalConfs('torsional_conf_gen')

# geometry optimization using smaller basis set (e.g. minix)
psi4EnergyCube1 = ParallelPsi4EnergyCalculation('parallel_psi4_energy_calculation1')
psi4EnergyCube1.title = 'Psi4_Cube1'
psi4EnergyCube1.set_parameters(opt_method = 'hf3c')
psi4EnergyCube1.set_parameters(spe_method = 'hf3c')
psi4EnergyCube1.set_parameters(opt_basis = 'minix')
psi4EnergyCube1.set_parameters(spe_basis = 'minix')

# geometry optimization using larger basis set (e.g. 6-31G*)
psi4EnergyCube2 = ParallelPsi4EnergyCalculation('parallel_psi4_energy_calculation2')
psi4EnergyCube2.title = 'Psi4_Cube2'
psi4EnergyCube2.opt_method = 'B3LYP'
psi4EnergyCube2.spe_method = 'B3LYP'
psi4EnergyCube2.opt_basis = '6-31G*'
psi4EnergyCube2.spe_basis = '6-31G**'

sysfail = DatasetWriterCube('sysfail')
failfs1 = DatasetWriterCube('failfs1')
failfs2 = DatasetWriterCube('failfs2')
failfs3 = DatasetWriterCube('failfs3')
ofs = DatasetWriterCube('ofs')
torfrags_ofs = DatasetWriterCube('torfrags_ofs')
startconfs_ofs = DatasetWriterCube('startconfs_ofs')
torconfs_ofs = DatasetWriterCube('torconfs_ofs')
[job.add_cube(n) for n in [ifs, fraggenCube, confgenCube, torsgenCube, psi4EnergyCube1,
                           psi4EnergyCube2, ofs, failfs1, failfs2, failfs3, fraggen_failure,
                           confgen_failure, sysfail, torfrags_ofs, startconfs_ofs, torconfs_ofs]]

# Promote parameters
ifs.promote_parameter('data_in', promoted_name='ifs')
fraggen_failure.promote_parameter('data_out', promoted_name='fraggen_failure',
                                  title='Fragment Generation Failures',
                                  description='Fragment Generation Failures',
                                  default='fraggen_failures')
confgenCube.promote_parameter('max_confs',
                              title='Maximum Alternate Starting Conformers',
                              description='Maximum number of starting conformations to use in QM torsion driving experiment. '+
                                          'NOTE: If you do not want to generate starting conformations, set this to 1.')
confgen_failure.promote_parameter('data_out', promoted_name='confgen_failure',
                                  title='Conformer Generation Failures',
                                  description='Conformer Generation Failures',
                                  default='confgen_failures')
torsgenCube.promote_parameter('num_points',
                              title='Number of torsional conformers to generate.',
                              description="""The number of evenly spaced torsion angles to sample 
                              when generating torsional conformers.""")
sysfail.promote_parameter('data_out', promoted_name='sysfail', title='System Failures')
failfs1.promote_parameter('data_out', promoted_name='torfail', title='torfail',
                          description="Torsional Conformer Generation Failures",
                          default='torgen_failures')
ofs.promote_parameter('data_out', promoted_name='ofs',
                      description='Floe output',
                      default='output',
                      title='Successes')
failfs2.promote_parameter('data_out', promoted_name='psifail', title='psifail',
                          description="Psi4 Failures",
                          default='psi4_failures')
failfs3.promote_parameter('data_out', promoted_name='sysfail', title='sysfail',
                          description="Psi4 System Failures",
                          default='psi4_sysfailures')
torfrags_ofs.promote_parameter('data_out', promoted_name='torfrags_ofs',
                      description='Torsional Fragments output',
                      default='torfrags',
                      title='Torsional Fragments')
startconfs_ofs.promote_parameter('data_out', promoted_name='startconfs_ofs',
                      description='Starting Conformers output',
                      default='startconfs',
                      title='Starting Conformers')
torconfs_ofs.promote_parameter('data_out', promoted_name='torconfs_ofs',
                      description='Torsional Conformers output',
                      default='torconfs',
                      title='Torsional Conformers')


# Connect ports
ifs.success.connect(fraggenCube.intake)
fraggenCube.success.connect(confgenCube.intake)
fraggenCube.success.connect(torfrags_ofs.intake)
fraggenCube.failure.connect(fraggen_failure.intake)

confgenCube.success.connect(torsgenCube.intake)
confgenCube.success.connect(startconfs_ofs.intake)
confgenCube.failure.connect(confgen_failure.intake)

torsgenCube.success.connect(psi4EnergyCube1.intake)
torsgenCube.success.connect(torconfs_ofs.intake)
torsgenCube.failure.connect(failfs1.intake)

psi4EnergyCube1.success.connect(psi4EnergyCube2.intake)
psi4EnergyCube1.failure.connect(psi4EnergyCube2.intake)
psi4EnergyCube1.system_failure.connect(psi4EnergyCube2.intake)
psi4EnergyCube1.system_failure.connect(sysfail.intake)

psi4EnergyCube2.success.connect(ofs.intake)
psi4EnergyCube2.failure.connect(failfs2.intake)
psi4EnergyCube2.system_failure.connect(failfs3.intake)

# If called from command line, run the floe
if __name__ == "__main__":
    job.run()

import os
import logging
logging.basicConfig(level=logging.DEBUG,
        format='%(asctime)s' + os.getenv('HOSTNAME') + '%(levelname)-8s  %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S')
from torsion.cubes import ParallelPsi4EnergyCalculation, SerialGenerateStartingConfs, SerialGenerateTorsionalConfs
from floe.api import WorkFloe
from floe.api import OEMolOStreamCube
from floe.api import OEMolIStreamCube

# Declare Floe, add metadata for UI
job = WorkFloe('Massively Parallel Multi Step Dihedral Scan of Torsions')
job.description = """
Calculate the Energy Profile associated with a Dihedral Scan in a massively parallel manner.
"""
job.classification = [["Torsion"]]

# Declare Cubes
ifs = OEMolIStreamCube('ifs')
ifs.promote_parameter('data_in', promoted_name='ifs')

confgenCube = SerialGenerateStartingConfs('starting_conf_gen')
confgenCube.promote_parameter('max_confs',
                              title='Maximum Alternate Starting Conformers',
                              description='Maximum number of starting conformations to use in QM torsion driving experiment. '+
                                          'NOTE: If you do not want to generate starting conformations, set this to 1.')

confgen_failure = OEMolOStreamCube('confgen_failure')
confgen_failure.promote_parameter('data_out', promoted_name='confgen_failure',
                                  title='Conformer Generation Failures',
                                  description='Conformer Generation Failures',
                                  default='confgen_failures')

torsgenCube = SerialGenerateTorsionalConfs('torsional_conf_gen')
torsgenCube.promote_parameter('num_points',
                              title='Number of torsional conformers to generate.',
                              description="""The number of evenly spaced torsion angles to sample 
                              when generating torsional conformers.""")

# geometry optimization using smaller basis set (e.g. 6-31G)
psi4EnergyCube1 = ParallelPsi4EnergyCalculation('parallel_psi4_energy_calculation1')
psi4EnergyCube1.title = 'Psi4_Cube1'
psi4EnergyCube1.opt_basis = '6-31G'
psi4EnergyCube1.spe_basis = '6-31G'

# geometry optimization using larger basis set (e.g. 6-31G*)
psi4EnergyCube2 = ParallelPsi4EnergyCalculation('parallel_psi4_energy_calculation2')
psi4EnergyCube2.title = 'Psi4_Cube2'
psi4EnergyCube2.opt_basis = '6-31G*'
psi4EnergyCube2.spe_basis = '6-31G**'

sysfail = OEMolOStreamCube('sysfail')
sysfail.promote_parameter('data_out', promoted_name='sysfail', title='System Failures')

failfs1 = OEMolOStreamCube('failfs1')
failfs1.promote_parameter('data_out', promoted_name='torfail', title='torfail',
                          description="Torsional Conformer Generation Failures",
                          default='torgen_failures')

ofs = OEMolOStreamCube('ofs')
ofs.promote_parameter('data_out', promoted_name='ofs',
                      description='Floe output',
                      default='output',
                      title='Successes')
failfs2 = OEMolOStreamCube('failfs2')
failfs2.promote_parameter('data_out', promoted_name='psifail', title='psifail',
                          description="Psi4 Failures",
                          default='psi4_failures')

failfs3 = OEMolOStreamCube('failfs3')
failfs3.promote_parameter('data_out', promoted_name='sysfail', title='sysfail',
                          description="Psi4 System Failures",
                          default='psi4_sysfailures')

# Add Cubes to Floe
[job.add_cube(n) for n in [ifs, confgenCube, torsgenCube, psi4EnergyCube1,
                           psi4EnergyCube2, ofs, failfs1, failfs2, failfs3,
                           confgen_failure, sysfail]]

# Connect ports
ifs.success.connect(confgenCube.intake)
confgenCube.success.connect(torsgenCube.intake)
confgenCube.failure.connect(confgen_failure.intake)
torsgenCube.success.connect(psi4EnergyCube1.intake)
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

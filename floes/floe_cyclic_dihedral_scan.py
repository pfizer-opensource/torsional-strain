import os
import logging
logging.basicConfig(level=logging.DEBUG,
        format='%(asctime)s' + os.getenv('HOSTNAME') + '%(levelname)-8s  %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S')
from torsion.cubes import ParallelPsi4EnergyCalculation, SerialGenerateStartingConfs, SerialGenerateTorsionalConfs, SerialDriveTorsion
from floe.api import WorkFloe
from floe.api import OEMolOStreamCube
from floe.api import OEMolIStreamCube

# Declare Floe, add metadata for UI
job = WorkFloe('Cyclic Dihedral Scan of Torsions')
job.description = """
Calculate the Energy Profile associated with a Dihedral Scan in a cyclic manner.
"""
job.classification = [["Torsion"]]

# Declare Cubes
ifs = OEMolIStreamCube('ifs')
ifs.promote_parameter('data_in', promoted_name='ifs')

confgenCube = SerialGenerateStartingConfs('starting_conf_gen')
confgenCube.promote_parameter('max_confs',
                              title='Maximum Alternate Starting Conformers',
                              description="""Maximum number of starting
                              conformations to use in QM torsion driving
                              experiment. NOTE: If you do not want to generate
                              starting conformations, set this to 1.""")


# We will use the following cube to split the starting conformers into
# individual OEMols with the appropriate conformer labels
torsgenCube = SerialGenerateTorsionalConfs('torsional_conf_gen')
torsgenCube.set_parameters(num_points=1)    # Don't Drive Torsions here.

torsionDriverCube = SerialDriveTorsion('torsion_driver')
torsionDriverCube.promote_parameter('num_points',
                              title="""Number of points at which to sample
                              the dihedral.""",
                              description="""The number of evenly spaced
                              torsion angles to sample in order to determine
                              the enthalpy surface.""")

psi4EnergyCube = ParallelPsi4EnergyCalculation('parallel_psi4_energy_calculation')
psi4EnergyCube.promote_parameter('spe_method')
psi4EnergyCube.promote_parameter('spe_basis')
psi4EnergyCube.promote_parameter('geom_opt_technique')
psi4EnergyCube.promote_parameter('opt_method')
psi4EnergyCube.promote_parameter('opt_basis')
psi4EnergyCube.promote_parameter('geom_maxiter')
psi4EnergyCube.set_parameters(only_selected_conformer=True)

ofs = OEMolOStreamCube('ofs')
ofs.promote_parameter('data_out', promoted_name='ofs')

# Add Cubes to Floe
[job.add_cube(n) for n in [ifs, confgenCube, torsgenCube, torsionDriverCube, psi4EnergyCube, ofs]]

# Connect ports
ifs.success.connect(confgenCube.intake)
confgenCube.success.connect(torsgenCube.intake)
torsgenCube.success.connect(torsionDriverCube.intake)
torsionDriverCube.to_energy_calc.connect(psi4EnergyCube.intake)
psi4EnergyCube.success.connect(torsionDriverCube.intake)
psi4EnergyCube.failure.connect(torsionDriverCube.intake)
psi4EnergyCube.system_failure.connect(torsionDriverCube.intake)
torsionDriverCube.success.connect(ofs.intake)

# If called from command line, run the floe
if __name__ == "__main__":
    job.run()

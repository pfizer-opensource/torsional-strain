import logging
import socket
print('host:', socket.gethostname())
logging.basicConfig(level=logging.DEBUG, 
        format='%(asctime)s' + socket.gethostname() + '%(levelname)-8s  %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S')
from torsion.cubes import ParallelGenerateStartingConfs, ParallelGenerateTorsionalConfs
from floe.api import WorkFloe
from floe.api import ParallelCubeGroup
from floe.api import OEMolOStreamCube
from floe.api import OEMolIStreamCube

# Declare Floe, add metadata for UI
job = WorkFloe('Generate Starting Conformers')
job.description = """
Generate Starting and Torsional Conformers
"""
job.classification = [["Torsion"]]

# Declare Cubes
ifs = OEMolIStreamCube('ifs')
ifs.promote_parameter('data_in', promoted_name='ifs')
ifs.parameter_overrides["download_format"] = {"hidden": True}
ifs.parameter_overrides["limit"] = {"hidden":True}

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

torsgen_failure = OEMolOStreamCube('torsgen_failure')
torsgen_failure.promote_parameter('data_out', promoted_name='torsgen_failure',
                                  title='Torsional Conformer Generation Failures',
                                  description='Torsional Conformer Generation Failures',
                                  default='torsgen_failures')
torsgen_failure.promote_parameter('buffered', default=False)
torsgen_failure.parameter_overrides["buffered"] = {"hidden": True}

# Conformer/torsion cube group
group = ParallelCubeGroup(cubes=[confgenCube, torsgenCube])
job.add_group(group)

# final molecular output
ofs = OEMolOStreamCube('ofs')
ofs.promote_parameter('data_out', promoted_name='ofs',
                      description='Floe output',
                      default='output',
                      title='Successes')
ofs.promote_parameter('buffered', default=False)
ofs.parameter_overrides["buffered"] = {"hidden": True}

# Add Cubes to Floe
cubes = [ifs, confgenCube, torsgenCube, ofs, confgen_failure, torsgen_failure]
[job.add_cube(c) for c in cubes]

# Connect ports
ifs.success.connect(confgenCube.intake)
confgenCube.success.connect(torsgenCube.intake)
confgenCube.failure.connect(confgen_failure.intake)
torsgenCube.success.connect(ofs.intake)
torsgenCube.failure.connect(torsgen_failure.intake)

import logging
import socket
print('host:', socket.gethostname())
logging.basicConfig(level=logging.DEBUG,
        format='%(asctime)s' + socket.gethostname() + '%(levelname)-8s  %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S')
from torsion.cubes import ParallelGenerateFragments
from floe.api import WorkFloe
from floe.api import ParallelCubeGroup
from floe.api import OEMolOStreamCube
from floe.api import OEMolIStreamCube

# Declare Floe, add metadata for UI
job = WorkFloe('Generate Torsional Fragments')
job.description = """
Generate Torsional Fragments
"""
job.classification = [["Torsion"]]

# Declare Cubes
ifs = OEMolIStreamCube('ifs')
ifs.promote_parameter('data_in', promoted_name='ifs')
ifs.parameter_overrides["download_format"] = {"hidden": True}
ifs.parameter_overrides["limit"] = {"hidden":True}

# Fragment generation
fraggenCube = ParallelGenerateFragments('fraggenCube')

fraggen_failure = OEMolOStreamCube('fraggen_failure')
fraggen_failure.promote_parameter('data_out', promoted_name='fraggen_failure',
                                  title='Fragment Generation Failures',
                                  description='Fragment Generation Failures',
                                  default='fraggen_failures')
fraggen_failure.promote_parameter('buffered', default=False)
fraggen_failure.parameter_overrides["buffered"] = {"hidden": True}


# final molecular output
ofs = OEMolOStreamCube('ofs')
ofs.promote_parameter('data_out', promoted_name='ofs',
                      description='Floe output',
                      default='output',
                      title='Successes')
ofs.promote_parameter('buffered', default=False)
ofs.parameter_overrides["buffered"] = {"hidden": True}

# Add Cubes to Floe
cubes = [ifs, fraggenCube, ofs, fraggen_failure]
[job.add_cube(c) for c in cubes]

# Connect ports
ifs.success.connect(fraggenCube.intake)
fraggenCube.success.connect(ofs.intake)
fraggenCube.failure.connect(fraggen_failure.intake)

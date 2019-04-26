import logging
import socket
print('host:', socket.gethostname())
logging.basicConfig(level=logging.DEBUG,
        format='%(asctime)s' + socket.gethostname() + '%(levelname)-8s  %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S')

from floe.api import WorkFloe
from cuberecord.cubes import DatasetReaderCube, DatasetWriterCube
from torsion.cubes import ParallelGenerateFragments

# Declare Floe, add metadata for UI
job = WorkFloe('Generate Torsional Fragments')
job.description = """
Generate Torsional Fragments
"""
job.classification = [["Torsion"]]

# Declare cubes
ifs = DatasetReaderCube('ifs')
fraggenCube = ParallelGenerateFragments('fraggenCube')
fraggen_failure = DatasetWriterCube('fraggen_failure')
ofs = DatasetWriterCube('ofs')
cubes = [ifs, fraggenCube, ofs, fraggen_failure]


# Promote parameters
ifs.promote_parameter('data_in', promoted_name='data_in')
ofs.promote_parameter('data_out', promoted_name='ofs',
                      description='Floe output',
                      default='output',
                      title='Successes')
fraggen_failure.promote_parameter('data_out', promoted_name='fraggen_failure',
                                  title='Fragment Generation Failures',
                                  description='Fragment Generation Failures',
                                  default='fraggen_failures')

# Add Cubes to Floe
[job.add_cube(c) for c in cubes]

# Connect ports
ifs.success.connect(fraggenCube.intake)
fraggenCube.success.connect(ofs.intake)
fraggenCube.failure.connect(fraggen_failure.intake)

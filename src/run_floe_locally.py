import sys
from torsion.floes import MassivelyParallelDihedralScanFloe
from torsion.floes import CyclicDihedralScanFloe
from torsion.floes import MassivelyParallelMultiStepDihedralScanFloe

prog_name = sys.argv[0]

mode = 'multistep'
if '--floe-mode' in sys.argv:
    arg_index = [i for i, x in enumerate(sys.argv) if x == '--floe-mode']
    if len(arg_index) > 1:
        raise ValueError("Error: --floe-mode flag found more than once.")
    else:
        arg_index = arg_index[0]
        mode = sys.argv[arg_index + 1]
    if mode not in ['serial', 'parallel', 'multistep']:
        raise ValueError('Invalid value for mode. Only "serial", "parallel", or "multistep" are allowed.')

    del sys.argv[arg_index]    # remove '--mode' flag
    del sys.argv[arg_index]    # remove mode argument

#print("Mode: {}".format(mode))
#print(sys.argv)

if mode == 'parallel':
    print("Running MassivelyParallelDihedralScanFloe.")
    MassivelyParallelDihedralScanFloe.run()
elif mode == 'multistep':
    print("Running MassivelyParallelMultiStepDihedralScanFloe.")
    MassivelyParallelMultiStepDihedralScanFloe.run()
elif mode == 'serial':
    print("Running CyclicDihedralScanFloe.")
    CyclicDihedralScanFloe.run()

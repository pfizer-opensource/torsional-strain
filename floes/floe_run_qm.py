import logging
import socket
print('host:', socket.gethostname())
logging.basicConfig(level=logging.DEBUG, 
        format='%(asctime)s' + socket.gethostname() + '%(levelname)-8s  %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S')
from torsion.cubes import ParallelPsi4EnergyCalculation, HiddenParamParallelPsi4EnergyCalculation
from floe.api import WorkFloe
from floe.api import OEMolOStreamCube
from floe.api import OEMolIStreamCube

# Declare Floe, add metadata for UI
job = WorkFloe('Parallel 4-stage QM Energy Calculations')
job.description = """
Rerun the QM energy calculations for specific (potentially problematic) fragments 
that failed Dihedral Torsional Scan floes.
"""
job.classification = [["Torsion"]]

# Declare Cubes
ifs = OEMolIStreamCube('ifs')
ifs.promote_parameter('data_in', promoted_name='ifs')
ifs.parameter_overrides["download_format"] = {"hidden": True}
ifs.parameter_overrides["limit"] = {"hidden":True}

# Fast QM optimization with hf3c to get approximate geometry
hf3cCube = HiddenParamParallelPsi4EnergyCalculation('hf3c_precalc')
hf3cCube.set_parameters(spe_method='hf3c', spe_basis='minix', opt_method='hf3c', opt_basis='minix', geom_maxiter=200)

# Fast QM optimization with 3-21G to get approximate geometry
b3lyp321gCube = HiddenParamParallelPsi4EnergyCalculation('b3lyp321g_precalc')
b3lyp321gCube.set_parameters(spe_method='B3LYP', spe_basis='3-21G', opt_method='B3LYP', opt_basis='3-21G', geom_maxiter=200)

# Fast QM optimization with 6-11G to get approximate geometry
b3lyp631gCube = HiddenParamParallelPsi4EnergyCalculation('b3lyp631g_precalc')
b3lyp631gCube.set_parameters(spe_method='B3LYP', spe_basis='6-31G', opt_method='B3LYP', opt_basis='6-31G', geom_maxiter=200)

# Full QM optimization & SP for energy
psi4EnergyCube = ParallelPsi4EnergyCalculation('parallel_psi4_energy_calculation')
psi4EnergyCube.set_parameters(geom_maxiter=200)
psi4EnergyCube.promote_parameter('spe_method', required=True)
psi4EnergyCube.promote_parameter('spe_basis', required=True)
psi4EnergyCube.promote_parameter('molden_output', required=True)
psi4EnergyCube.promote_parameter('opt_method', required=True)
psi4EnergyCube.promote_parameter('opt_basis', required=True)

# hf3c failure handling
failFastQM = OEMolOStreamCube('fail_fast_QM')
failFastQM.promote_parameter('data_out', promoted_name='fast_qm_fail', title='fast_qm_fail',
                             description="fast QM Failures",
                             default='fast_QM_failures')
failFastQM.promote_parameter('buffered', default=False)
failFastQM.parameter_overrides["buffered"] = {"hidden": True}

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

# Add Cubes to Floe
cubes = [ifs, hf3cCube, psi4EnergyCube,
         ofs, failFastQM, failfs2, failfs3]
[job.add_cube(c) for c in cubes]

# Connect ports
ifs.success.connect(hf3cCube.intake)
hf3cCube.success.connect(psi4EnergyCube.intake)
hf3cCube.failure.connect(failFastQM.intake)
hf3cCube.system_failure.connect(psi4EnergyCube.intake)

psi4EnergyCube.success.connect(ofs.intake)
psi4EnergyCube.failure.connect(failfs2.intake)
psi4EnergyCube.system_failure.connect(failfs3.intake)

# If called from command line, run the floe
if __name__ == "__main__":
    job.run()

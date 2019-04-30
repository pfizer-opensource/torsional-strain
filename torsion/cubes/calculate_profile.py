import numpy as np
from scipy.interpolate import interp1d
from openeye import oechem

from torsion.utils import get_sd_data
from torsion.analysis import construct_dihedral_energy_profile, \
                             cal_molecule_torsion_strain

from cuberecord.cubes import OEMolRecordCube, InOutMolFieldMixin, InitMolRecordMixin
from cuberecord.ports import RecordOutputPort
from datarecord import OEMolRecord


class ProfileAssembler(OEMolRecordCube, InOutMolFieldMixin, InitMolRecordMixin):
    """Assemble torsional conformers and compute torsional energy profiles.
    """
    
    title = 'Torsional Profile Assembler'
    description = """
        Assembles QM geometry optimization+single point energy 
        calculations on individual torsional conformers into 
        torsional energy profiles.
    """
    classification = [['Energetics', 'Torsion']]
    tags = [tag for tag_list in classification for tag in tag_list]
    
    qm_conf_output = RecordOutputPort('qm_conf_output')
    
    def begin(self):
        self.fragment_library = {}
        self.input_molecules = {}
        for mol in self.get_init_molecules():
            self.log.warn(f"Cube initialized with {mol.GetTitle()}.")
            self.input_molecules[mol.GetTitle()] = oechem.OEMol(mol)
    
    
    def process(self, record, port):
        if record.has_value(self.args.in_mol_field):
            mol = record.get_value(self.args.in_mol_field)
        else:
            self.log.error("Could not find molecules in OEMolRecord")
            self.failure.emit(record)
            return
        
        conf = mol.GetActive()
        conf_title = get_sd_data(conf, "CONFORMER_LABEL")
        tor_atoms = get_sd_data(mol, "TORSION_ATOMS_ParentMol").split()
        parent_name = conf_title[:-3]
        dih_label = '_'.join(str(x) for x in tor_atoms)
        conf_idx = int(conf_title[-2:])
        
        if parent_name in self.fragment_library:
            if dih_label in self.fragment_library[parent_name]:
                self.fragment_library[parent_name][dih_label][conf_idx] = mol.CreateCopy()
            else:
                self.fragment_library[parent_name][dih_label] = [None]*24
                self.fragment_library[parent_name][dih_label][conf_idx] = mol.CreateCopy()
        else:
            self.fragment_library[parent_name] = {}
            self.fragment_library[parent_name][dih_label] = [None]*24
            self.fragment_library[parent_name][dih_label][conf_idx] = mol.CreateCopy()
            
        
    def end(self):
        for mol_name in self.fragment_library.keys():
            self.log.warn(f"Processing results for {mol_name}")
            profile_map = {}
            for dih_label, conformer_list in self.fragment_library[mol_name].items():
                results = construct_dihedral_energy_profile(conformer_list, 24)
                output_conformers, torsional_profile, torsion_inchi_list = results
                
                conf_output_record = OEMolRecord()
                conf_output_record.set_mol(output_conformers)
                self.qm_conf_output.emit(conf_output_record)
                
                xp = torsional_profile[:,0]
                yp = torsional_profile[:,1]
                xp1 = np.concatenate((xp[:-1]-360, xp, xp[1:]+360))
                yp1 = np.concatenate((yp[:-1], yp, yp[1:]))
                for inchi in torsion_inchi_list:
                    profile_map[inchi] = interp1d(xp1, yp1, kind='cubic', bounds_error=False)
                
            input_mol = self.input_molecules[mol_name]
            output_mol = input_mol.CreateCopy()
            cal_molecule_torsion_strain(output_mol, profile_map)
            record = OEMolRecord()
            record.set_mol(output_mol)
            self.success.emit(record)
            
            
            
            

                        
#!/usr/bin/env python3

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdmolfiles, Descriptors, rdMolTransforms, rdMolDescriptors, rdmolops
from rdkit.Chem import TorsionFingerprints, rdMolTransforms, PropertyMol, rdDistGeom, rdmolops, rdMolDescriptors, Descriptors
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions, GetStereoisomerCount
from rdkit import RDConfig
import numpy as np 
import pandas as pd
import multiprocessing as mp
import os, sys, glob, subprocess
from rdkit.Chem.AtomPairs import Utils
import itertools, copy, math
from itertools import combinations, product, chain
from rdkit.Geometry import Point3D
import itertools
from math import pow, sqrt
from collections import Counter
from sklearn.preprocessing import StandardScaler, MinMaxScaler

from functools import reduce
from openbabel import pybel

#ob_log_handler = pybel.ob.OBMessageHandler()
#ob_log_handler.SetOutputLevel(0) # Set to 0 to suppress all but critical messages
pybel.ob.obErrorLog.StopLogging()

class Assign_Atrop:
    """
    This class will recognize bisphosphine atropisomers with a biaryl system 
    and will return atrop_ls consisting of each conformer (in the 'CONFORMERS' dir)
    and the associated type ('A' or 'B')

    Can add in functionality to rename tag of conformer/structure based on atropisomer
    type. This differentiation is needed due to MacroModel not recognizing the chirality
    of atropisomers during the conformational search
    """
    def __init__(self):
        
        self.conf_dir = 'CONFORMERS'
        

    def pybel_conv(self, fn, in_format, out_format):
        """Converts a file from one format to another using Pybel."""
        mol = next(pybel.readfile(in_format, f'{fn}.{in_format}'))
        
        output = pybel.Outputfile(out_format, f'{fn}-TEMP.{out_format}', overwrite=True)
        output.write(mol)
        output.close()
        
        return None

    def get_mol(self, xyz_file):
        """Converts .xyz file to .sdf and returns RDKit molecule object."""
        mol2_file = xyz_file.replace('.xyz', '-TEMP.sdf')
        self.pybel_conv(xyz_file.replace('.xyz', ''), 'xyz', 'sdf')
        
        if os.path.exists(mol2_file):
            ligs = Chem.SDMolSupplier(mol2_file)
            mol = ligs[0]
            os.remove(mol2_file)
        else:
            mol=None
        
        return mol

    def rotate(self, mol, coords, axis, theta):
        """Rotates atoms in a molecule around a specified axis by theta radians."""
        Omol = copy.deepcopy(mol)
        for coord in coords:
            x, y, z = coord[1], coord[2], coord[3]
            if axis == 'X':
                Y = y * math.cos(theta) - z * math.sin(theta)
                Z = y * math.sin(theta) + z * math.cos(theta)
                mol.GetConformer().SetAtomPosition(coord[0], Point3D(x, Y, Z))
            elif axis == 'Y':
                X = x * math.cos(theta) + z * math.sin(theta)
                Z = z * math.cos(theta) - x * math.sin(theta)
                mol.GetConformer().SetAtomPosition(coord[0], Point3D(X, y, Z))
            elif axis == 'Z':
                X = x * math.cos(theta) - y * math.sin(theta)
                Y = x * math.sin(theta) + y * math.cos(theta)
                mol.GetConformer().SetAtomPosition(coord[0], Point3D(X, Y, z))

        return mol


    def calculate_angle(self, a, b, c):
        """Calculates angle in degrees between three 3D coordinates."""     
        ba = np.array(a) - np.array(b)
        bc = np.array(c) - np.array(b)

        cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
        angle = np.arccos(cosine_angle) 
        return np.degrees(angle)


    def adjust(self, mol, coords, adjust_by):
        """Shifts all atom coordinates in the molecule by adjust_by vector."""
        for coord in coords:
            x = coord[1] - adjust_by[0]
            y = coord[2] - adjust_by[1]
            z = coord[3] - adjust_by[2]
            mol.GetConformer().SetAtomPosition(coord[0], Point3D(x, y, z))
        return mol



    def get3Dcoords(self, mol, conf):
        """Extracts 3D coordinates of atoms in a molecule's conformer."""
        coords = []
        for atom in mol.GetAtoms():
            crds = list(conf.GetAtomPosition(atom.GetIdx()))
            coords.append([atom.GetIdx(),crds[0],crds[1],crds[2]])
        return coords


    def update_coords(self, mol):
        """Updates and returns the 3D coordinates from the molecule's conformer."""
        emol = copy.deepcopy(mol)
        confe = emol.GetConformer()
        return self.get3Dcoords(emol, confe)



      
    def get_adjust_atoms(self, mol, ligand_P_idxs, atrop_bond):
        valid = False
        bid_path = list(rdmolops.GetShortestPath(mol, ligand_P_idxs[0], ligand_P_idxs[1]))
        for atom in mol.GetAtomWithIdx(bid_path[-2]).GetNeighbors():
            if atom.GetSymbol() != 'H' and atom.GetIdx() not in bid_path:
                atom_nei = atom.GetIdx()
                valid = True
        
        
        if not valid:
            return None
            
        
        sorted_path = [
            bid_path[1], 
            bid_path[int((len(bid_path)/2)-1)], 
            bid_path[int((len(bid_path)/2))], 
            bid_path[-2], 
            atom_nei
        ]

        return sorted_path



    def update_all(self, mol, coords, angle_atoms, axis):
        """Performs rotation and returns updated coordinates."""
        a, b, c = angle_atoms
        angle = np.radians(self.calculate_angle(a, b, c))
        return mol, self.rotate(mol, coords, axis, angle), self.update_coords(mol)


    
    def adjust_mol(self, mol, sorted_path):
        typ=None
        x0y0z0 = [0,0,0]
        
        x1y0z0 = [1,0,0]
        x0y1z0 = [0,1,0]
        x0y0z1 = [0,0,1]
        
        
        coords = self.update_coords(mol)
        
        
        coord_path = []
        for atm in sorted_path:
            coord = [x for x in coords if x[0] == atm]
            coord_path.append(coord)
        
        
        adjust_by = [(coord_path[2][0][1]-0), (coord_path[2][0][2]-0), (coord_path[2][0][3]-0)] 
        mol = self.adjust(mol, coords, adjust_by)
        coords = self.update_coords(mol)
        
        
        coord_path = []
        for atm in sorted_path:
            coord = [x for x in coords if x[0] == atm]
            coord_path.append(coord)
        
        
        Omol = copy.deepcopy(mol)
        
        angle_atoms = [x0y0z1, x0y0z0, [0, coord_path[3][0][2], coord_path[3][0][3]]]
        Omol, mol, coords = self.update_all(Omol, coords, angle_atoms, 'X')
        
        coord_path = []
        for atm in sorted_path:
            coord = [x for x in coords if x[0] == atm]
            coord_path.append(coord)
        
        
        
        b3 = [coord_path[2][0][1], coord_path[2][0][2], coord_path[2][0][3]]
        b4 = [coord_path[3][0][1], coord_path[3][0][2], coord_path[3][0][3]]
        
        
       
        Omol = copy.deepcopy(mol)
        
        angle_atoms = [x1y0z0, x0y0z0, [coord_path[3][0][1], 0, coord_path[3][0][3]]]
        Omol, mol, coords = self.update_all(Omol, coords, angle_atoms, 'Y')
        
        coord_path = []
        for atm in sorted_path:
            coord = [x for x in coords if x[0] == atm]
            coord_path.append(coord)
        
        a1, a2, a3, a4, a5 = coord_path[0][0], coord_path[1][0], coord_path[2][0], coord_path[3][0], coord_path[4][0]
        
        if a1[3] < 0 and a5[2] < 0:
            typ = 'A'
        elif a1[3] > 0 and a5[2] > 0:
            typ = 'A'
        elif a1[3] > 0 and a5[2] < 0:
            typ = 'B'
        
        elif a1[3] < 0 and a5[2] > 0:
            typ = 'B'
        
        if a5[1] < 0:
            if typ == 'A':
                typ == 'B'
            elif typ == 'B':
                typ == 'A'
        
        
        return typ
            
            
    
        
    def get_inter_dists(self, conf, inter_idx, inter_dist):
        for idx_pair in inter_idx:
            at1Coord = np.array(conf.GetAtomPosition(idx_pair[0]))
            at2Coord = np.array(conf.GetAtomPosition(idx_pair[1]))
            at_dist = np.linalg.norm(at2Coord - at1Coord)
            inter_dist.append(at_dist)
        return inter_dist


    def assign_atropisomer_type(self, fn):
        """Assigns atropisomer type 'A' or 'B' for a conformer file."""
        biaryl_smi = 'C1(C2=CC=CC=C2)=CC=CC=C1'
        benz_smi = 'C1=CC=CC=C1'
        mol = self.get_mol(fn)
        if not mol or not mol.HasSubstructMatch(Chem.MolFromSmiles(biaryl_smi)):
            return None

        biaryl_match = mol.GetSubstructMatches(Chem.MolFromSmiles(biaryl_smi))
        sing_ring = mol.GetSubstructMatches(Chem.MolFromSmiles(benz_smi))

        bi_ls = list(biaryl_match[0])
        matched_rings = [list(ring) for ring in sing_ring if all(x in bi_ls for x in ring)]

        atrop_bond = []
        for bond in mol.GetBonds():
            at1, at2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if at1 in bi_ls and at2 in bi_ls:
                atrop_bond.extend([at1, at2])
            for ring in matched_rings:
                if all(x in ring for x in [at1, at2]):
                    if at1 in atrop_bond: atrop_bond.remove(at1)
                    if at2 in atrop_bond: atrop_bond.remove(at2)

        P_atoms = sorted([atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == 'P'])
        sorted_path = self.get_adjust_atoms(mol, P_atoms, atrop_bond)
        return self.adjust_mol(mol, sorted_path) if sorted_path else None



    def detect_atropisomers(self):
        """Detects and returns list of ligands with assigned atropisomer types."""
        ### coordinating_atoms are the two atom symbols that 
        ### would be bound in a bidentate fashion 
        atrop_ls = []
        updated_fns = []
        bad_list = []
        fns = [fn for fn in glob.glob(f'{self.conf_dir}/*.xyz') if 'errocene' not in fn]
        for fn in fns:
            try:
                
                biaryl_smi = 'C1(C2=CC=CC=C2)=CC=CC=C1'
                benz_smi = 'C1=CC=CC=C1'
                biaryl_mol = Chem.MolFromSmiles(biaryl_smi)
                benz_mol = Chem.MolFromSmiles(benz_smi)
                
                
                mol = self.get_mol(fn)
                
                Omol = copy.deepcopy(mol)
                
                if not mol.HasSubstructMatch(biaryl_mol):
                    continue
                
                biaryl_match = mol.GetSubstructMatches(biaryl_mol)
                
                
                sing_ring = mol.GetSubstructMatches(benz_mol)
                matched_rings = []
                bi_ls = list(biaryl_match[0])
                for ring in sing_ring:
                    if all(x in bi_ls for x in list(ring)):
                        matched_rings.append(list(ring))
                
                atrop_bond = []
                for bond in mol.GetBonds():
                    at1 = (bond.GetBeginAtom()).GetIdx()
                    at2 = (bond.GetEndAtom()).GetIdx()
                    if at1 in bi_ls and at2 in bi_ls:
                        atrop_bond.append(at1)
                        atrop_bond.append(at2)
                        
                    for ring in matched_rings:
                        if all(x in list(ring) for x in [at1,at2]):
                            atrop_bond.remove(at1)
                            atrop_bond.remove(at2)
                
                
                frag_mol = Chem.FragmentOnBonds(mol, atrop_bond)
                frag_mol = Chem.DeleteSubstructs(frag_mol, Chem.MolFromSmiles('[Pd+2]C([CH]c1ccccc1)[CH]c1ccccc1'))
                frag_smi = Chem.MolToSmiles(frag_mol)
                
                num_mols = len(frag_smi.split('.'))
                if num_mols != 2:
                    continue
                
                if not ('P' in str(frag_smi.split('.')[0])) or not ('P' in str(frag_smi.split('.')[1])):
                    continue
                
                
                tag = os.path.basename(fn.replace('.xyz', ''))
                
                
                updated_fns.append(tag)
            except:
                bad_list.append(fn)
            
        
        
        updated_fns = list(set(updated_fns))
        for fn in updated_fns:
            typ = self.assign_atropisomer_type(f'{self.conf_dir}/{fn}.xyz')
            atrop_ls.append([fn, typ])
        
        return atrop_ls



#atrop_j = Assign_Atrop()
#atrop_ls = atrop_j.detect_atropisomers()
#print(atrop_ls)

### NEW CLASS

class curate_data:
    """
    A data processing class for cleaning, transforming, and analyzing ligand feature data.
    Performs Boltzmann averaging, high-energy pruning, structural parsing, and dataset curation.
    """

    def __init__(
        self,
        exclude_irrelavent_features=True,
        boltz_99=True,
        remove_high_energy_confs=False,
        drop_low_data=False,
        remove_rmsd_nrg_duplicates=False,
        export_csv=False,
        outfile=None
    ):
        self.csv_tags = [
            'ENERGIES', 'ORBITALS-TZVP', 'NBO-TZVP',
            'COSMO_general', 'COSMO_mecn', 'GEOMETRIES',
            'STERIMOL2VEC', 'VOL2VEC', 'BURIED-VOLUME'
        ]

        self.sdf_dir = 'SDFs'
        self.conf_dir = 'CONFORMERS'
        self.data_dir = 'CSV_DATA'

        # Constants
        self.kb = 8.3144621  # J / K / mol
        self.KCAL_TO_KJ = 4.184
        self.J_TO_KCAL = 4.184 * 1000.0
        self.T = 298.15
        self.AU_TO_KCAL = 627.5

        # Feature exclusions
        self.excluded_features = [
            'c.p.', 'trace', 'Edisp', 'E8', 'E6', 'log10',
            'E_vdw', 'cp_mix', 'ring_cor', 'H_HB', 'H_bond(accept)',
            'H_bond(donor)', 'H_int', 'H_MF', 'H_vdw'
        ]

        self.exclude_irrelavent_features = exclude_irrelavent_features
        self.boltz_99 = boltz_99
        self.remove_high_energy_confs = remove_high_energy_confs
        self.export_csv = export_csv
        self.drop_low_data = drop_low_data
        self.remove_rmsd_nrg_duplicates = remove_rmsd_nrg_duplicates
        self.outfile = outfile



    def high_e_filter(self, energy, e_thr=10.0, quiet=True):
        """Returns indices of conformers with relative energy above threshold."""
        mine = min(energy)
        high_e=[]
        for i in range(len(energy)):
            if energy[i] > mine + e_thr:
                if not quiet:
                    print("HIGH E Conformer!",i)
                high_e.append(i)
        return(high_e)
    
    ## Remove high NRG structures based on threshold 
    
    def rm_HiNRG(self, df):
        """Removes conformers that exceed a defined energy threshold within each ligand group."""
        df_noHiE = pd.DataFrame()
        
        for lig in list(set(list(df['Ligand'].values))):
            nrow = df[['Conformer','def2-TZVP-GAS']].loc[df['Ligand'] == lig]
            high_e_idxs = self.high_e_filter(nrow['def2-TZVP-GAS'].values, e_thr=4.0)
            remove_confs = [nrow['Conformer'].values[x] for x in high_e_idxs]
            nrow = df.loc[(df['Ligand'] == lig) & (~df['Conformer'].isin(remove_confs))]
            df_noHiE = pd.concat([nrow, df_noHiE])
        
        df_noHiE.reset_index()
        
        return df_noHiE
    
    def get_boltz_avg(self, prop, energy):
        """Computes Boltzmann-weighted average of a property based on conformer energies."""
        gas_const = self.kb / 1000  # kJ / mol K
        temperature = 298.15
        boltz_sum = 0.0
        e_zero = 0
        
        for j, e in enumerate(energy):
            if j == 0: e_zero = e
            e_rel = e - e_zero
            boltz_sum += math.exp(-e_rel/(gas_const*self.T))
        weights = []
        for e in energy:
            if j == 0: e_zero = e
            e_rel = e - e_zero
            weight = math.exp(-e_rel/(gas_const*self.T)) / boltz_sum
            weights.append(weight)
        boltz_avg = 0.0
        for i,p in enumerate(prop):
            boltz_avg += p * weights[i]
        
        return boltz_avg

    def boltz_avg_data(self, dataset, energy_type):
        """Applies Boltzmann averaging to all properties grouped by ligand."""
        ligands, props, boltz = [], [], []
        weighted_df = pd.DataFrame(columns=['Ligand','Property','Boltz_avg'])
        
        for ligand in list(set(list(dataset['Ligand'].values))):
            
            nrow = dataset.loc[dataset['Ligand'] == ligand]
            
            for key in nrow.keys():
                if key in ['Ligand','Conformer']:
                    continue
                
                boltz_avg = self.get_boltz_avg(nrow[key],nrow[energy_type])
                ligands.append(ligand)
                props.append(key)
                boltz.append(boltz_avg)
                
        
        weighted_df['Ligand'] = ligands
        weighted_df['Property'] = props
        weighted_df['Boltz_avg'] = boltz
        
        
        return weighted_df


    def scale_data(self, df, min_r, max_r):
        """Scales numerical data in DataFrame using min-max scaling."""
        scaled_df = df.copy()
        numeric_list = ['int', 'float', 'number'] 
        # select columns based on the above list
        numeric_data = scaled_df.select_dtypes(include=numeric_list)
        scaler = MinMaxScaler(feature_range=(min_r,max_r))
        scaled_df[numeric_data.columns] = scaler.fit_transform(scaled_df[numeric_data.columns])
        
        return scaled_df



    def get_bridge_ls(self, df):
        remove = np.array(
            [1,2,3,4,5,8,9,10,11,12,13,14,15,16,17,
            18,19,20,21,22,23,24,25,26,27,28,29,30,31],
            dtype=int
        )
        
        remove = remove - 1
        remove = np.flip(remove)
        
        bridge_ls=[]
        ligands = list(set(list(df['Ligand'].values)))
        for lig in ligands:
            if 'A_cs' in lig:
                fn = f'{lig.replace("A_cs","_cs")}_out.sdf'
            elif 'B_cs' in lig:
                fn = f'{lig.replace("B_cs","_cs")}_out.sdf'
            else:
                fn = f'{lig}_out.sdf'
            filename = f'{self.sdf_dir}/{fn}'
            if not os.path.exists(filename):
                filename = filename.replace('_out.sdf', '_1_out.sdf')

            if not os.path.exists(filename):
                continue

            suppl = Chem.SDMolSupplier(filename, sanitize=True, removeHs=False)
            mol = suppl[0]
            
            if mol is None:
                continue

            newmol = Chem.EditableMol(mol)
            for atom in remove:
                newmol.RemoveAtom(int(atom))
            newmol = newmol.GetMol()
            Chem.SanitizeMol(newmol)
            newmol = Chem.RemoveHs(newmol)
            Chem.Compute2DCoords(newmol)
            
                
            P_idx = [atom.GetIdx() for atom in newmol.GetAtoms() if atom.GetSymbol() == 'P']
            path = Chem.GetShortestPath(newmol,P_idx[0],P_idx[1])
            if len(path) == 0:
                n_bridge_atoms=None
            else:
                n_bridge_atoms=len(path)-2
            
            bridge_type = []
            for idx in path[1:-1]:
                atom = newmol.GetAtomWithIdx(idx)
                if str(atom.GetHybridization()) != 'SP':
                    
                    bridge_type.append(atom.GetSymbol()+'('+str(atom.GetHybridization())+')')
                else:
                    bridge_type.append(atom.GetSymbol()+'(SP2)')
            
            bridge = Counter(bridge_type)
            
            if bridge == None:
                bridge_str=None
            else:
                label = []
                for key in sorted (bridge.keys()) :
                    
                    label.append(f'{bridge[key]}{key}')
                    
                bridge_str= ','.join(label)
            
            
            #newrow['bridge_counter'] = bridge
            bridge_ls.append([lig, bridge_str])
        
        return bridge_ls


    def get_fps(self, fn):
        """Returns RDKit fingerprint from an SDF file if readable."""
        if not os.path.exists(fn): 
            return None
         
        mol = Chem.SDMolSupplier(fn, sanitize=True,removeHs=True)[0]
        if not mol: 
            return None
         
        fp = Chem.RDKFingerprint(mol)
         
        return fp
         
         
    def check_dup_2D(self, fn1, fn2):
        """Compares two ligands via 2D fingerprint similarity (Tanimoto)."""
        fp1, fp2 = self.get_fps(fn1), self.get_fps(fn2)
        fn1 = os.path.basename(fn1.replace('_out.sdf',''))
        
        fn2 = os.path.basename(fn2.replace('_out.sdf',''))
        
        if '99' in fn1 or '99' in fn2:
            fna, fnb = fn1, fn2
            fn1, fn2 = fnb, fna
        
        t = [fn1 if fp1 else None]
            
            
        if fp1 and fp2:
            tp = [[fn1, fn2], Chem.DataStructs.TanimotoSimilarity(fp1,fp2)]
        
        
        elif fp1:
            tp = [[fn1, None], 0]
        
        elif fp2:
            tp = [[fn2, None], 0]
        
        else:
            tp = [[None, None], -1]
            
        
        
        return tp
    

            
    def get_unique_ligands(self):
        """Identifies and groups structurally redundant ligands by 2D similarity."""
        x_ls = list(range(1, 198))
        struct_nums = [x for x in x_ls if x % 2 != 0]
        duplicates_2D = [(f'{self.sdf_dir}/{x}_cs_out.sdf',f'{self.sdf_dir}/{x+1}_cs_out.sdf') for x in struct_nums]
        
        for fn in glob.iglob(f'{self.sdf_dir}/ferrocene*_allyl_cs_out.sdf'):
            fn2 = fn.replace('_cs_', '_2_cs_')
            duplicates_2D.append((fn, fn2))
        
        for fn in glob.iglob(f'{self.sdf_dir}/*_1_out.sdf'):
            fn2 = fn.replace('_1_out', '_2_out')
            duplicates_2D.append((fn, fn2))
        
        
        tan_pairs = []
        for fn in duplicates_2D:
            fp1, fp2 = self.get_fps(fn[0]), self.get_fps(fn[1])
            tan_pair = self.check_dup_2D(fn[0], fn[1])
            tan_pairs.append(tan_pair)
        
        
        same_2Ds = [(x[0][0], x[0][1]) for x in tan_pairs if float(x[1]) == 1]
        diff_2Ds = [(x[0][0], x[0][1]) for x in tan_pairs if 0 <= float(x[1]) < 1]
        error_2Ds = [(x[0][0], x[0][1]) for x in tan_pairs if x[1] == -1]
        
        
        all_valid = same_2Ds + diff_2Ds
        unique_ligand_sets = [x for x in all_valid if x[0] is not None and x[1] is not None]
        
        
        return unique_ligand_sets


    def get_rmsd(self, sdf_filename, i, j):
        suppl = Chem.SDMolSupplier(sdf_filename, sanitize=True,removeHs=True)
        # "checking rmsd for",i,j
        mol_i = suppl[i]
        mol_j = suppl[j]
        try:
            rmsd = Chem.GetBestRMS(mol_i, mol_j)
        except:
            rmsd = 100.0
            #print("no match with confs",i,j)
        return rmsd


    def check_rmsd_duplicates(self, sdf_filename, e_thr=0.31375, i_thr=np.array([181.,181.,181.]), r_thr=0.2):
        confs = []
        energy = []

        name_tag = os.path.basename(sdf_filename.split('_out.sdf')[0])
        
        lines = open(f'{self.data_dir}/ENERGIES.csv', 'r').readlines()
        for i, line in enumerate(lines):
            if i != 0:
                x = line.split(',')[0]
                if name_tag == x:
                    conf = line.split(',')[1]
                    e = line.split(',')[2]
                    energy.append(float(e))
                    confs.append(conf)

        num_confs = len(energy)
        
        d_energy = np.zeros((num_confs, num_confs), dtype=float)
        
        for i in range(num_confs):
            energy_i = energy[i]
            
            for j in range(num_confs):
                energy_j = energy[j]
                
                if i >= j: 
                    continue 
                d_energy[i,j] = d_energy[j,i] = abs(energy_i - energy_j)
                
        
        #check if rmsd file exists, if not create new one 
        csv_name = sdf_filename.split('_out.sdf')[0] + '.csv'
        if len(glob.glob(csv_name)) == 0:
            csv_exists = False
            rmsd_df  = pd.DataFrame(columns = ['mol1','mol2','rmsd'])
        else:
            csv_exists = True
            rmsd_df = pd.read_csv(csv_name,header=0)
            
        
        #check duplicates
        duplicate = []
        for i in range(num_confs):
            for j in range(num_confs):
                if i >= j: 
                    continue 
                if d_energy[i,j] <= e_thr: 
                    #within E threshold 
                    

                    #within rmsd threshold
                    if csv_exists:
                        rmsd = rmsd_df.loc[rmsd_df['mol1'] == float(i)].loc[rmsd_df['mol2']== float(j),['rmsd']].values

                    else:
                        try:        
                            rmsd = self.get_rmsd(sdf_filename,i,j)
                            df = {'mol1': i, 'mol2': j, 'rmsd': rmsd}
                            rmsd_df = rmsd_df.append(df, ignore_index = True)
                        except:
                            continue
  
                    
                    if rmsd <= r_thr:
                        
                        #duplicate.append(j)
                        duplicate.append([confs[j],name_tag])
        
        if not csv_exists:
            rmsd_df.to_csv(csv_name,index=False)

        #duplicate = sorted(list(set(duplicate)))
        return duplicate




    def export_cregen_xyz(self, sdf_fns, energy, num_confs):
        first_mol_flag = False
        second_mol_flag = False
        for n,fn in enumerate(sdf_fns):
            file = sdf_fns[n][0]
            conf_ls = sdf_fns[n][1]
            suppl = Chem.SDMolSupplier(file,sanitize=True,removeHs=False)
            
        min_e = sys.float_info.max
        min_e_idx = 400
        for i in range(num_confs):
            
            if min_e > energy[i]:
                min_e = energy[i]
                min_e_idx = i
    
        
        xyz_data = []
        min_conf = []
        
        n_num=0
        for n,fn in enumerate(sdf_fns):
            file = sdf_fns[n][0]
            conf_ls = sdf_fns[n][1]
            suppl = Chem.SDMolSupplier(file,sanitize=True,removeHs=False)
            for i,mol in enumerate(suppl):
                data = Chem.MolToXYZBlock(mol)
                data = data.split('\n')
                data = data[0:-1]
                
                
                data[1] = f'Energy = {energy[n_num]/627.5} !{file.split(".sdf")[0]}_{str(i)}'
    
                if n_num == min_e_idx:
                    min_conf = data
      
                xyz_data.extend(data)
                
                n_num+=1
        
        
        # get new xyz filename
        combined_xyz_name = []
        sdf_filenames = [x[0] for x in sdf_fns]
        for file in sdf_filenames:
            combined_xyz_name.append(file.split('_out.sdf')[0])
        combined_xyz_name = '_'.join(combined_xyz_name)
        
        best_xyz_name = combined_xyz_name + '_best.xyz'
        combined_xyz_name = combined_xyz_name + '.xyz'
        
        
        with open(combined_xyz_name, 'w') as file:
            [file.write(line+'\n') for line in xyz_data]
        with open(best_xyz_name, 'w') as file: 
            [file.write(line+'\n') for line in min_conf]
            
        return best_xyz_name, combined_xyz_name


    def run_cregen_99(self):
        fns = []
        
        for xyz in glob.iglob(f'{self.conf_dir}/*.xyz'):
            fn = os.path.basename(xyz)
            fn = fn.split('_CONF_')[0]
            fns.append(x[:-2])
        
        fns = list(set(fns))
        
        energy_list = []
        
        
        for fn in fns:
            num_confs = 0
            nrgs = []
            
            fns_pair = [x for x in self.unique_ligand_sets if fn in x[0]]
            fns_ls=[]
            for f in fns_pair:
                
                conf_ls = []
                lines = open(f'{self.data_dir}/ENERGIES.csv', 'r').readlines()
                for i, line in enumerate(lines):
                    if i != 0:
                        x = line.split(',')[0]
                        e = float(line.split(',')[2])
                        conf = line.split(',')[1]
                        conf = int(conf.replace('CONF_',''))
                        if fn in x:
                            nrgs.append(e)
                            
                            num_confs += 1
                            conf_ls.append(conf)
                fns_ls.append([f'{f}_out.sdf',conf_ls])
            
            best_xyz_name, combined_xyz_name = self.export_cregen_xyz(fns_ls, nrgs, num_confs)
            # Runs crest conformation filters 
            os.system(f'crest {best_xyz_name} -cregen {combined_xyz_name} -notopo -ethr 0.31375 -rthr 0.2 -ewin 4.0')
    

    # 99% Ensemble Parsing
    def cregen_99_confs(self, df):
        conf_list = []
        df99 = pd.DataFrame()
        
        for file in sorted(glob.glob('cregen_final_filter/files/*.out')):
            with open(file) as f:
                lines = f.readlines()
            filtered_confs= []
            conf_start, conf_end = -1,-1
            
            for i in range(len(lines)):
                if 'Erel/kcal' in lines[i]:
                    conf_start = i+1
                if 'T /K' in lines[i]:
                    conf_end = i
            
            ensemble_count = 0
            confs = []
            for line in lines[conf_start:conf_end]:
                if len(line.split()) > 5:
                    ensemble_count += float(line.split()[4])
                    confs.append(line.split()[-1])
                if ensemble_count > 0.99: 
                    break
            
            conf_list.append(confs)
        
        flat_confs_ls = list(chain(*conf_list))
        flat_confs_ls = list(set(flat_confs_ls))
        
        
        for mol_name in flat_confs_ls:
            x = mol_name.split('_out_')
            ligand = x[0]
            y = str(int(x[1]) + 1)
            conf = f'CONF_{y.zfill(4)}'
            nrow = df.loc[(df['Ligand'] == ligand) & (df['Conformer'] == conf)]
            df99 = pd.concat([nrow, df99])
        
        df99.reset_index()

        return df99

    def generate_dataset(self):
        """
        Main function to generate processed dataset:
        - Loads CSV data
        - Computes features
        - Applies Boltzmann filtering if enabled
        - Removes irrelevant features

        Returns:
            pd.DataFrame: Cleaned and feature-enhanced ligand dataset.
        """
        self.unique_ligand_sets = self.get_unique_ligands()

        df_list = [pd.read_csv(f'{self.data_dir}/{fn}.csv') for fn in self.csv_tags]
        combined_df = pd.concat(df_list).groupby(['Ligand', 'Conformer'], as_index=False).first()
        combined_df = combined_df.drop_duplicates().dropna(how='all')
        combined_df = combined_df.dropna(thresh=combined_df.shape[1] - 3)

        df = combined_df.copy()

        # Add derived features
        df['HOMO-LUMO Gap'] = df['LUMO'] - df['HOMO']
        df['Area/Volume Ratio'] = df['Area'] / df['Volume']
        df['Pd-P difference'] = abs(df['5-6'] - df['5-7'])
        df['Avg. Pd-P distance'] = (df['5-6'] + df['5-7']) / 2
        df['Atom67 difference'] = abs(df['Atom6'] - df['Atom7'])
        df['Vbur_67 difference'] = abs(df['Vbur_6'] - df['Vbur_7'])

        e_types = ["def2-SVP", "def2-TZVP-GAS", "def2-TZVP-COSMO"]
        for e in e_types:
            df[e] = abs(df[e] / self.KCAL_TO_KJ)

        df['dG-solv'] = df['def2-TZVP-GAS'] + df['c.p.'] + df['Gsolv']
        df['dG'] = df['def2-TZVP-COSMO'] + df['c.p.']

        print(df)
        if self.remove_rmsd_nrg_duplicates:
            rm_rmsd_dupl = []
            for sdf in glob.glob(f'{self.sdf_dir}/*.sdf'):
                dupl = self.check_rmsd_duplicates(sdf)
                #name_tag = os.path.basename(sdf.split('_out.sdf')[0])

                #rm_rmsd_dupl.append([name_tag, dupl])

            removal_set = set(map(tuple, dupl))  # Convert to set of tuples for fast lookup

            mask = df.apply(lambda row: (row['Conformer'], row['Ligand']) in removal_set, axis=1)
            df = df[~mask].reset_index(drop=True)

        print(df)

        if self.boltz_99:
            df = self.cregen_99_confs(df)

        Gsolvs = ['G(toluene)', 'G(MeOH)', 'G(dmf)', 'G(water)', 'G(MeCN)']
        for solv in Gsolvs:
            dGsolv = f'dGs({solv})'
            df[dGsolv] = df[solv] + df['c.p.'] + df['Gsolv']

        if self.exclude_irrelavent_features:
            df = df.drop(columns=self.excluded_features, errors='ignore')

        return df


df_j = curate_data(
    exclude_irrelavent_features=True, 
    boltz_99=True,  
    remove_high_energy_confs=False, 
    drop_low_data=False, 
    remove_rmsd_nrg_duplicates=True, 
    export_csv=False, 
    outfile=None
)

df = df_j.generate_dataset()
#print(unique_ligand_sets)
#df.to_csv('EXAMPLE.csv', index=False)












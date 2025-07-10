

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem, rdmolfiles, rdPartialCharges, rdmolops
import os, sys, argparse
import subprocess as sp
from math import floor
import math, glob
import numpy as np



class extend_ligand_metal:
    def __init__(self):
        self.m = ''
        
        
        
    def mae_to_mol2(self, mae_file):
        
        name = os.path.basename(mae_file)
        mol2_file = name.replace('.mae', '-TEMP.mol2')
        conversion = sp.Popen([f'mol2convert -imae {mae_file} -omol2 {mol2_file}'], shell=True)
        conversion.communicate()
        
        
                
                
        inputfile = open(mol2_file, 'r').readlines()
        writefile = open(mol2_file,'w')
        
        for line in inputfile:
            if 'PD' in line or 'Pd' in line:
                if 'Any' in line:
                    
                    writefile.write(line.replace('Any', 'Pd '))
                    
                elif line.split(' ')[5] == 'H':
                             
                    writefile.write(line.replace(' H ', ' Pd'))
                
                else:
                    writefile.write(line)
            
            else:
                writefile.write(line)
        writefile.close()
        
        mol = Chem.MolFromMol2File(mol2_file, removeHs=False)
        #os.system(f'rm -r {mol2_file}')
        print(mol2_file)
        return mol
    
    
    def extend_vect(self, len_extend, start_point, end_point):
        
        x1 = end_point[0]
        y1 = end_point[1]
        z1 = end_point[2]
        x0 = start_point[0]
        y0 = start_point[1]
        z0 = start_point[2]
        
        L1 = math.sqrt((x1-x0)**2+(y1-y0)**2+(z1-z0)**2)
        L2 = len_extend
        
        x = (x1-x0)*(L2/L1) + x0
        y = (y1-y0)*(L2/L1) + y0
        z = (z1-z0)*(L2/L1) + z0
        
        new_point = [x,y,z]
        
        return new_point
        
    
    def extend_metal_pos(self, fn, num, len_extend=3.0):
        print(fn)
        mol = self.mae_to_mol2(fn)
        P_idxs = [x.GetIdx() for x in mol.GetAtoms() if x.GetSymbol() == 'P']
        Rh_idxs = [x.GetIdx() for x in mol.GetAtoms() if x.GetSymbol() == 'Pd']
        
            
            
        inputfile = open(fn, 'r').readlines()
        
        """
        P_idxs = []
        Rh_idxs = []
        for line in inputfile:
            try:
                if line.split(' ')[10] == 'RH':
                    Rh_idxs.append(line.split('')[0])
                
                elif line.split(' ')[10] == 'P':
                    P_idxs.append(line.split('')[0])
            except:
                continue
        """     
                
        all_P_env_coords = []
        dots = []
        for i, line in enumerate(inputfile):
            if ':::' in line:
                dots.append(i)
        
        for i, line in enumerate(inputfile):
            
            if dots[2] < i < dots[3]:
                if float(line.split()[0]) == float(P_idxs[0]+1):
                    P1, P1i = [float(line.split()[2]),float(line.split()[3]),float(line.split()[4])], i
                elif float(line.split()[0]) == float(P_idxs[1]+1):
                    P2, P2i = [float(line.split()[2]),float(line.split()[3]),float(line.split()[4])], i
                elif float(line.split()[0]) == float(Rh_idxs[0]+1):
                    Rh, Rhi = [float(line.split()[2]),float(line.split()[3]),float(line.split()[4])], i
        
        P1_P2_mid = [
            (P1[0]+P2[0])/2,
            (P1[1]+P2[1])/2,
            (P1[2]+P2[2])/2
        ]
        
        
        
        
        
        Rh_pos = self.extend_vect(len_extend, P1_P2_mid, Rh)
        x = "{:.6f}".format(Rh_pos[0])
        y = "{:.6f}".format(Rh_pos[1])
        z = "{:.6f}".format(Rh_pos[2])
        
        
        new_fn = fn.replace('-2.mae', f'-{num}.mae')
        
        wf = open(new_fn, 'w')
        for i, line in enumerate(inputfile):
            if i == Rhi:
                xo = line.split()[2]
                zo = line.split()[4]
                posx = line.find(xo)
                posz = line.find(zo)
                posz = posz + len(zo)
                line = f'{line[:posx]}{x} {y} {z}{line[posz:]}'
                
                
            wf.write(line)
            
            
        wf.close()
        
        return new_fn
    
    

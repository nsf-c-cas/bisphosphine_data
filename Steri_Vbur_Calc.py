import argparse
import os
import glob
import math
import numpy as np
import pandas as pd
from concurrent import futures
from rdkit import Chem
from rdkit.Chem import rdmolfiles, rdmolops
import dbstep.Dbstep as db

# Constants
GAS_CONST = 1.987204e-3  # kcal/mol·K
TEMP = 298.15
E_WIN = 10.0
MAX_WORKERS = 36
INTERVAL = np.arange(0.0, 12.5, 0.5)
SCAN_INT = '0.0:12.0:0.5'
ATOM1 = "5"
ATOM2 = "6,7"
EXCLUDE = ",".join(str(i) for i in range(1, 32))

def parse_args():
    parser = argparse.ArgumentParser(description="Steric scan job runner.")
    parser.add_argument('--calc', required=True, choices=['Sterimol', 'Vbur', 'V2V'],
                        help="Calculation type: 'Sterimol', 'Vbur', or 'V2V'")
    return parser.parse_args()

def load_energy_data(fns):
    energy_dict = {fn: [] for fn in fns}
    with open('CSV_DATA/ENERGIES.csv', 'r') as f:
        for line in f.readlines()[1:]:
            mol, conf, e = line.strip().split(',')[:3]
            if mol in energy_dict:
                energy_dict[mol].append([conf, e])
    return [[k, v] for k, v in energy_dict.items()]

def dbstep_scan(file, calc_type):
    if calc_type == 'Sterimol':
        mol = db.dbstep(file, measure="classic", atom1=ATOM1, atom2=ATOM2,
                        exclude=EXCLUDE, sterimol=True, commandline=True)
        return mol.Bmin, mol.Bmax, mol.L
    elif calc_type == 'Vbur':
        results = []
        for atom in ['5', '6', '7']:
            mol = db.dbstep(file, atom1=atom, exclude=EXCLUDE, volume=True, commandline=True)
            results.append(mol.bur_vol)
        return tuple(results)
    elif calc_type == 'V2V':
        mol = db.dbstep(file, atom1=ATOM1, atom2=ATOM2, exclude=EXCLUDE,
                        sterimol=True, volume=True, scan=SCAN_INT, commandline=True)
        return mol.Bmin, mol.Bmax, mol.bur_vol, mol.bur_shell

def run_jobs(job_list, calc_type):
    results = []
    with futures.ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures_list = [executor.submit(dbstep_scan, job[0], calc_type) for job in job_list]
        for fut in futures_list:
            results.append(fut.result())
    return results

def setup_jobs(filenames, energy_list):
    job_list, n_confs = [], []
    for m, filename in enumerate(filenames):
        n_confs.append(0)
        confs = [f"CONF_{line.split('_CONF_')[1][:4]}" for line in open(filename) if '_CONF_' in line]
        suppl = Chem.SDMolSupplier(filename, sanitize=False, removeHs=False)
        mol_energies = [x[1] for x in energy_list if x[0] in os.path.basename(filename)][0]

        e_zero, boltz_sum = 0, 0
        for n, mol in enumerate(suppl):
            mol_energy = float([x[1] for x in mol_energies if x[0] == confs[n]][0])
            if n == 0: e_zero = mol_energy
            e_rel = mol_energy - e_zero
            if e_rel < E_WIN:
                boltz_sum += math.exp(-e_rel / (GAS_CONST * TEMP))

        for n, mol in enumerate(suppl):
            mol_energy = float([x[1] for x in mol_energies if x[0] == confs[n]][0])
            if n == 0: e_zero = mol_energy
            e_rel = mol_energy - e_zero
            if e_rel < E_WIN:
                boltz_fac = math.exp(-e_rel * 1000 / (GAS_CONST * TEMP)) / boltz_sum
                smiles = Chem.MolToSmiles(mol)
                cn = confs[n].replace('CONF_', '')
                xyz_file = filename.split('.')[0] + f'_{cn}.xyz'
                rdmolfiles.MolToXYZFile(mol, xyz_file)
                job_list.append([xyz_file, m, n_confs[m], smiles, e_rel, boltz_fac])
                n_confs[m] += 1
    return job_list, len(filenames), n_confs

def build_dataframe(calc, job_list, results, filenames, n_confs):
    if calc == 'Sterimol':
        columns = ['molecule', 'smiles', 'e_rel', 'boltz', 'Bmin', 'Bmax', 'L']
    elif calc == 'Vbur':
        columns = ['molecule', 'smiles', 'e_rel', 'boltz', 'Vbur_5', 'Vbur_6', 'Vbur_7']
    else:
        columns = ['molecule', 'smiles', 'e_rel', 'boltz'] + \
                  [f'Bmin_{r}' for r in INTERVAL] + [f'Bmax_{r}' for r in INTERVAL] + \
                  [f'Vbur_{r}' for r in INTERVAL] + [f'Vshell_{r}' for r in INTERVAL]

    master_df = pd.DataFrame(columns=columns)
    wshell_df = pd.DataFrame(columns=columns) if calc == 'V2V' else None

    entry = 0
    for i, name in enumerate(filenames):
        for j in range(n_confs[i]):
            for n, job in enumerate(job_list):
                mol, conf = job[1], job[2]
                if mol == i and conf == j:
                    smiles, e_rel, boltz = job[3:]
                    data = [name + "_" + str(j), smiles, e_rel, boltz]

                    if calc == 'Sterimol':
                        data += list(results[n])
                    elif calc == 'Vbur':
                        data += list(results[n])
                    elif calc == 'V2V':
                        bmin, bmax, vol, shell = results[n]
                        data += bmin + bmax + vol + shell
                        for k in range(len(INTERVAL)):
                            # Accumulate Boltzmann-weighted contributions
                            wshell_df.at[i, columns[4 + k]] = wshell_df.get(columns[4 + k], 0) + boltz * bmin[k]
                            wshell_df.at[i, columns[4 + len(INTERVAL) + k]] = wshell_df.get(columns[4 + len(INTERVAL) + k], 0) + boltz * bmax[k]
                            wshell_df.at[i, columns[4 + 2 * len(INTERVAL) + k]] = wshell_df.get(columns[4 + 2 * len(INTERVAL) + k], 0) + boltz * vol[k]
                            wshell_df.at[i, columns[4 + 3 * len(INTERVAL) + k]] = wshell_df.get(columns[4 + 3 * len(INTERVAL) + k], 0) + boltz * shell[k]

                    master_df.loc[entry] = data
                    entry += 1

    return master_df, wshell_df

def main():
    args = parse_args()
    calc = args.calc

    filenames = glob.glob('SDFs/*.sdf')
    fns = [os.path.basename(fn.split('_out')[0]) for fn in filenames]
    energy_list = load_energy_data(fns)

    job_list, n_mols, n_confs = setup_jobs(filenames, energy_list)
    results = run_jobs(job_list, calc)
    master_df, wshell_df = build_dataframe(calc, job_list, results, filenames, n_confs)

    if calc == 'Sterimol':
        master_df.to_csv('raw_overall_data.csv', index=False)
    elif calc == 'Vbur':
        master_df.to_csv('raw_overall_vbur.csv', index=False)
    elif calc == 'V2V':
        wshell_df.drop(columns=['e_rel', 'boltz'], inplace=True, errors='ignore')
        wshell_df.to_csv('weighted_shell_scan.csv', index=False)
        master_df.to_csv('raw_data.csv', index=False)

if __name__ == "__main__":
    main()



### USAGE
"""
python script.py --calc Sterimol
python script.py --calc Vbur
python script.py --calc V2V

"""






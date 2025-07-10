#!/usr/bin/env python

import os, subprocess, glob, math, re, argparse, sys



from check_atrop import *
from extend_metal_position import *
from submit_cs_job import *




class merge_structures:
    def __init__(self, settings1, settings2, job_settings, merge_run, template):
        
        
        
        for k,v in settings1.items():
            setattr(self, k, v)
        
        for k,v in settings2.items():
            setattr(self, k, v)
        
        for k,v in job_settings.items():
            setattr(self, k, v)
            
        self.cwd = os.getcwd()
        self.old_cwd = self.cwd.split('RUN-')[0]
        self.q2mm_path = os.path.join(self.old_cwd, self.q2mm)
        self.template = os.path.join(self.old_cwd, template)
        self.key = 'f_m_ct { '
        self.job_settings = job_settings
        self.merge_run = merge_run
        
        self.merged_TS_path = self.merged_path
        if not os.path.exists(self.merged_TS_path):
            os.mkdir(self.merged_TS_path)



    def cs_com(self, com_fn, chig_line):

        wf = open(com_fn,'w')
        wf.write(f'{com_fn.replace("_cs.com", "_mini.mae")}\n')
        wf.write(f'{com_fn.replace(".com", ".mae")}\n')
        wf.write(' DEBG      55    179      0      0     0.0000     0.0000     0.0000     0.0000\n')
        wf.write(' SEED   40000      0      0      0     0.0000     0.0000     0.0000     0.0000\n')
        wf.write(' MMOD       0      1      0      0     0.0000     0.0000     0.0000     0.0000\n')
        wf.write(' FFLD       2      1      0      0     1.0000     0.0000     0.0000     0.0000\n')
        wf.write(' BDCO       0      0      0      0    41.5692 99999.0000     0.0000     0.0000\n')
        wf.write(' READ       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n')
        wf.write(f'{chig_line}')
        wf.write(' CRMS       0      0      0      0     0.0000     0.5000     0.0000     0.0000\n')
        wf.write(f' MCMM{self.cs_steps.rjust(8)}      0      0      0     0.0000     0.0000     0.0000     0.0000\n')
        wf.write(' MCNV       1      5      0      0     0.0000     0.0000     0.0000     0.0000\n')
        wf.write(' MCSS       2      0      0      0    21.0000     0.0000     0.0000     0.0000\n')
        wf.write(' MCOP       1      0      0      0     0.0000     0.0000     0.0000     0.0000\n')
        wf.write(' DEMX       0    166      0      0    21.0000    42.0000     0.0000     0.0000\n')
        wf.write(' MSYM       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n')
        wf.write(' AUTO       0      2     -1      1     0.0000    -1.0000     0.0000     0.0000\n')
        wf.write(' CONV       2      0      0      0     0.0500     0.0000     0.0000     0.0000\n')
        wf.write(' MINI       9      1    500      0     0.0000     0.0000     0.0000     0.0000\n')
        wf.write(' END        0      0      0      0     0.0000     0.0000     0.0000     0.0000\n')
        wf.close()
        
    def mini_com(self, com_fn, chig_line):
        wf = open(com_fn.replace('_cs.com', '_mini.com'),'w')
        wf.write(f'{com_fn.replace("_cs.com", ".mae")}\n')
        wf.write(f'{com_fn.replace("_cs.com", "_mini.mae")}\n')
        wf.write(' DEBG      55    179      0      0     0.0000     0.0000     0.0000     0.0000\n')
        wf.write(' SEED   40000      0      0      0     0.0000     0.0000     0.0000     0.0000\n')
        wf.write(' FFLD       2      1      0      0     1.0000     0.0000     0.0000     0.0000\n')
        wf.write(' EXNB       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n')
        wf.write(' BDCO       0      0      0      0    89.4427 99999.0000     0.0000     0.0000\n')
        wf.write(' CRMS       0      0      0      0     4.1840     0.2500    60.0000     0.0000\n')
        wf.write(' BGIN       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n')
        wf.write(' READ       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n')
        wf.write(f'{chig_line}')
        wf.write(' CONV       2      0      0      0     0.0500     0.0000     0.0000     0.0000\n')
        wf.write(' MINI       1      0   2500      0     0.0000     0.0000     0.0000     0.0000\n')
        wf.write(' END        0      0      0      0     0.0000     0.0000     0.0000     0.0000\n')
        wf.close()
    
    def re_com(self, com_fn, chig_line):
        wf = open(com_fn.replace('.com', '_re.com'),'w')
        wf.write(f'{com_fn.replace(".com", ".mae")}\n')
        wf.write(f'{com_fn.replace(".com", "_re.mae")}\n')
        wf.write(' DEBG      55    179      0      0     0.0000     0.0000     0.0000     0.0000\n')
        wf.write(' SEED   40000      0      0      0     0.0000     0.0000     0.0000     0.0000\n')
        wf.write(' FFLD       2      1      0      0     1.0000     0.0000     0.0000     0.0000\n')
        wf.write(' EXNB       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n')
        wf.write(' BGIN       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n')
        wf.write(' READ       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n')
        wf.write(f'{chig_line}')
        wf.write(' COMP       0      0      0      0     0.0000     0.0000     2.0000     0.0000\n')
        wf.write(' MINI       9      0   2500      0     0.0000     0.0000     0.0000     0.0000\n')
        wf.write(' END        0      0      0      0     0.0000     0.0000     0.0000     0.0000\n')
        wf.close()
    
    
    def write_coms(self, fn):
        chig_line=''
        com_fn = fn.replace('.mae', '_cs.com')
        chig=False
        out_fn = fn.replace('.mae', '_cs.mae')
        if os.path.exists(com_fn):
            chig=True
            
            lines = open(com_fn, 'r').readlines()
            for line in lines:
                if 'CHIG' in line:
                    chig_line+=line
                    
                elif 'FXTA' in line:
                    chig_line+=line
        
        
            
        conversion = subprocess.Popen(
        [f'$SCHRODINGER/run {self.q2mm_path}/screen/setup_com_from_mae.py {fn} {com_fn} {out_fn} -n {self.cs_steps}'], shell=True)
        conversion.communicate()
    
        
        
        
        if not chig:
            lines = open(com_fn, 'r').readlines()
            
            for line in lines:
                if 'CHIG' in line:
                    chig_line+=line
            
        
        self.mini_com(com_fn, chig_line)
        self.re_com(com_fn, chig_line)
        self.cs_com(com_fn, chig_line)
    
    
    
    def com_writer(self):
        j1 = find_atrop(self.merged_TS_path, self.q2mm_path)
        j1.run()
        
        
        
        
        for fn in glob.iglob(f'{self.merged_TS_path}/*.mae'):
            if '_cs' not in fn and '_mini' not in fn:
                self.write_coms(fn)
            
            
                
    def write_job(self):
        j1 = get_jobs(self.job_settings, self.MACROMODEL_JOB_SCRIPT, self.q2mm_path)
        
        j1.run_all(self.merged_TS_path, submit_job=self.submit_job)    

    
    
    def split_merged(self, merge_tag):
        out_fn = f'{merge_tag}-X.mae'
        lines = open(out_fn, 'r').readlines()
        num = 0
        fnt = []
        fn_d = {}
        
        for i, line in enumerate(lines):
            if self.key in line:
                fnt.append(i)
        
        for num in range(int(len(fnt))):
            fnf = '{\n s_m_m2io_version\n :::\n 2.0.0\n}\n\n' 
            
            new_lines = []
            for i, line in enumerate(lines):
                try:
                    if i < fnt[num+1]:
                        new_lines.append(line)
                except:
                    if i >= fnt[num]:
                        new_lines.append(line)
            
            
            fn_split = f'{self.merged_TS_path}/{merge_tag}_{str(num+1)}.mae'
            print(fn_split)
            if num+1 == int(len(fnt)):
                with open(fn_split, 'w') as wf:
                    wf.write(fnf)
                    for line in new_lines:
                        wf.write(line)
            
            else:
                with open(fn_split, 'w') as wf:
                    for line in new_lines:
                        wf.write(line)
            
            
            
        
            
        os.remove(out_fn)         
    
    
    
    
    def run_merge(self, structs):
        run = True
        ligand = os.path.join(self.old_cwd, f'{self.ligand_dir}/{structs[0]}-2.mae')
        substrate = os.path.join(self.old_cwd, f'{self.substrate_dir}/{structs[1]}-2.mae')
        
        if not os.path.exists(substrate):
            run=False
            
        
        if not os.path.exists(ligand):
            run=False
            
        merge_tag = f'TS-{structs[0]}_{structs[1]}'
        check_fn = os.path.join(self.merged_TS_path, f'{merge_tag}_1.mae')
        if os.path.exists(check_fn):
            run=False
        if run:
            
            j2 = extend_ligand_metal()
            ligand_ext = j2.extend_metal_pos(ligand, f'9{structs[2]}')
            conversion = subprocess.Popen(
            [f'$SCHRODINGER/run {self.q2mm_path}/screen/merge.py -g {self.template} -g {substrate} -g {ligand_ext} -o {merge_tag}-X.mae'], shell=True)
            conversion.communicate()
            
            self.split_merged(merge_tag)
            
        
        return None
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
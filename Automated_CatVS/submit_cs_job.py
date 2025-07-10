#!/usr/bin/env python3

import argparse
import subprocess as sp
import sys, math, os, glob, time, re
from multiprocessing import Pool


import structure_property_settings




class get_jobs:
    def __init__(self, job_settings, MACROMODEL_JOB_SCRIPT, q2mm):
        
        
        settings={
            'MAIL' : 'ae',
            'QUEUE' : 'long',
            'RETRY' : 'n'
        }
        
        self.settings = settings
        for k, v in self.settings.items():
            setattr(self, k, v)
        
        
        for k, v in job_settings.items():
            setattr(self, k, v)
            
        self.temp = '{mytemp}'
        self.user = '{USER}'
        self.jid = '{JOB_ID}'
        self.mm_job_script = MACROMODEL_JOB_SCRIPT
        self.q2mm = q2mm
        
                
    def split_files(self, filenames, n):
        for x in range(0, len(filenames), n):
            every_chunk = filenames[x: n+x]
            if len(every_chunk) < n:
                every_chunk = every_chunk + \
                    [None for y in range(n-len(every_chunk))]
            yield every_chunk
            
    
    def setup_mm_scripts(self, filenames, job_name):
        cs_match = re.match('cs', job_name)
        re_match = re.match('re', job_name)
        with open(job_name, 'w') as f:
            f.write(self.mm_job_script.format(
                    self.email, self.MAIL, self.QUEUE, 
                    self.RETRY, 'mytemp', self.user, 
                    self.jid, self.temp, self.q2mm, 
                    self.temp, self.temp, self.temp, self.temp
                    )
                )
            
            
            
            
            for fn in filenames: 
                f.write(f'bmin -WAIT {fn}\n')
                    
                    
        
        print(f'WROTE: {job_name}')
        if self.submit_job:
            os.system(f'qsub {job_name}')
            
    
    
    
    def get_mm_scripts(self, filenames):
        
        if int(self.number_of_mm_jobs) == 1:
            job_name = f'MM_job_1.sh'
            self.setup_mm_scripts(filenames, job_name)
            
        else:
            num_per_job = math.ceil(len(filenames)/int(self.number_of_mm_jobs)) 
            split_files = list(self.split_files(filenames, num_per_job))
            for i in range(int(self.number_of_mm_jobs)):
                job_name = f'MM_job_{str(i+1)}.sh'
                self.setup_mm_scripts(split_files[i], job_name)


   
                
    
    def run_all(self, merged_dr, submit_job=False):
        self.submit_job = submit_job
        
        fns = list(glob.iglob(f'{merged_dr}/*.mae'))
        filenames = []
        for fn in fns:  
            if '_cs' not in fn:
                fn1 = fn.replace('.mae', '_mini')
                fn2 = fn.replace('.mae', '_cs')
                fn3 = fn.replace('.mae', '_cs_re')
                filenames.append(fn1) 
                filenames.append(fn2) 
                filenames.append(fn3)  
        
        
        
        
        self.get_mm_scripts(filenames)
        






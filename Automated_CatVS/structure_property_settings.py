#!/usr/bin/env python3

MACROMODEL_JOB_SCRIPT = \
"""#!/bin/bash

#$ -M {}
#$ -m {}
#$ -q {}
#$ -r {}
#$ -pe smp 1

module load schrodinger/2018u3


"""


job_settings = {
    'number_of_mm_jobs' : '5',
    'email' : 'bstenforcrcroot@gmail.com',
    'MACROMODEL_JOB_SCRIPT' : MACROMODEL_JOB_SCRIPT,
    'cs_steps' : '15000',                              # Number of steps for conformational search (15000)
    'submit_job' : True
}
### NOTE there is no path to mm3.fld and atom.typ listed as these
### cannot have a different name (Schrodinger only recognizes those file names)


# To use different forcefields, a directory must be created with a seperate mm3.fld and atom.typ 
# file for each unique force field. you must submit the job script for conformational searches
# in the specified directory containing the unique parameter files (see 'merged_TS_path' in the dict below to change)


in_path_settings={
    'ligand_dir' : 'Ligands',                          # Predefined folders containing all of ligands to be screened
    'substrate_dir' : 'Substrates',                    # Predefined folders containing all of ligands to be screened
    'q2mm' : 'q2mm',                                   # Path to the Q2MM/CatVS folder (traditionally just 'q2mm')
}
    
    
out_path_settings={
    'merged_path' : 'MERGED_TS'                        # Path used for a unique run (-FF1 as a tag to differentiate the FFs)
}
    





## The following commands are covered in more detail in the q2mm/screen/README.md file
mae_commands="""
 b_cs_first_match_only
 b_cs_use_substructure
"""




# Order of numbers tracks to above commands. Below shows b_cs_first_match_only as 1 and b_cs_use_substructure as 0.
# Please see more info in the READ_ME files in the q2mm folder 
ligand_true_false="""  
  1
  0
"""

substrate_true_false="""  
  0
  0
"""


structure_settings = {
    'mae_commands' : mae_commands,
    'substrate_true_false' : substrate_true_false,
    'ligand_true_false' : ligand_true_false,
    'substrate_pattern' : 'O=C-N-C=C',             # A pattern which can match with all substrates to be merged onto a template
    'ligand_pattern' : 'P-[Rh]-P',
}




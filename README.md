
# Bisphosphine Conformational Dependance Overview

Available code for the investigation into conformational dependace of features for Pd[allyl] bisphosphine complexes. The reported observation can be found at INSERT DOI/CITATION.




## Installation

Download the available code and data from this github branch. 

For visualization: utilize the existing data file containing features of all conformers (Bisphos_all_features.csv) along with the 'CSV_DATA' , 'SDFs', and 'CONFORMERS' directories. The 'Violin_Plots.py' file offers three options for plotting and exporting figures.

If you want to further filter the conformations, inlcude energy features, or further curate/add to the data, use the 'Compile_Data.py' file. 

If you plan on creating new ligand complexes, you must  use our in house Q2MM package (https://github.com/ericchansen/q2mm). The original automated workflow for the Q2MM process used for this work can be found at 'q2mm/smiles_to_catvs/scripts'. This workflow must be used to do the initial coordination to palladium. Once this is done, the code in the 'Automated_CatVS' directory can be used alond with the corrresponding reaction template (template.mae) and allylic substrate (allyl.mae). Requires access to MacroModel.

BASF code?


    
## Usage/Examples

#### Steri_Vbur_Calc.py

This code allows  for the calculation of steric parameters from the program developed in the Paton lab, DBSTEP (https://github.com/patonlab/DBSTEP). Install DBSTEP before use. 

Examples for getting Sterimol, Sterimol2Vec, Percent Buried Volume and Vol2Vec parameter sets:

'''
python3 script.py --calc Sterimol
python3 script.py --calc Vbur
python3 script.py --calc V2V
'''


#### Compile_Data.py

The Compile_Data.py code uses the inital .xyz conformers to further filter (via various methods) the conformations and extract relevent descriptors from the raw data in 'CSV_DATA'. Funtions in this file are also used in the below Violin_Plots.py file and must be used together. 

The code is split into two classes. The first class('Assign_Atrop') is a lengthy segment on identifying atropisomers within a conformational ensamble. The second class ('curate_data') concatinates and curates the data.

#### Violin_Plots.py

Three plotting options...

For below plotting, all ligands are used in a single plot, one for each specified property.
'''python
properties = ['Vbur_5','Atom5','Atom6']
j.plot_BB_all_ligs(properties, dpi)
'''
Example output for 'Bmax_3.5':


For below plotting, define exact properties and ligands. No backbone key will be given here.
'''python
properties = ['Vbur_5','Atom5','Atom6']
ligands = ['100_cs', '101_cs']
outfig = 'Sterics_NBO_2_ligands'
j.plot_defined_ligs_props(properties, ligands, dpi, outfig)
'''

Below plotting is for the feature range bar graph, can inlcude/exclude ferrocene structures.
'''python
outfile = 'Feature_ranges_all_ligands'
j.plot_feat_range(dpi, outfile, ferro=True, noferro=False)
'''

#### Automated_CatVS

A few python files in here...

'extend_metal_position.py' uses a .mae file of the ligand coordinated to palladium and extendds the metal opposide of the phosphine linker. This helps ensure phosphine groups will not be trapped due to coordination and is only needed in some cases.

'merge.py' utilizes the preexisting merge code in the Q2MM package ('q2mm/screen/merge.py') along with:
- The reaction template ('template.mae') 
- The allylic substrate ('allyl.mae') 
- The MM3 custom parameter file (mm3.fld)
- The atom type file (atom.typ)

All provided in the directory. This code will merge the species, giving the Pd[allyl] complex and write the command files (.com) for MacroModel. Three command files are written for an initial minimization, Monte Carlo search, and redundant conformer elimination. 

'submit_cs_job.py' takes these command files and writes them seuentially to a job file (.sh). The number of jobs is user defined and can be set in the 'structure_property_settings.py' file.


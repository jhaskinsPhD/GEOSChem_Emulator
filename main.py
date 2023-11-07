import os 
import sys  

sys.path.append('/uufs/chpc.utah.edu/common/home/u6044586/python_scripts/modules/GEOSChem_Emulator')
import GEOSChem_Emulator as gcem

# Define the path to the KPP file of the version of GEOS-Chem you want to emulate: 
kppfile='/uufs/chpc.utah.edu/common/home/u6044586/python_scripts/modules/GEOSChem_Emulator/Input_Files/fullchem.eqn'

# Define the path to the gckPP Rates file(s) of the version of GEOS-Chem you want to emulate: 
rate_files='/uufs/chpc.utah.edu/common/home/u6044586/python_scripts/modules/GEOSChem_Emulator/Input_Files/gckpp_Rates.F90'

# Define the version of GEOS-Chem you want to emulate: 
GC_version='13.3.4';

# Define how the j-value matching should be done to translate GEOS-Chem 
# photolysis rates to their corresponding MCM photolysis rate. Allowed options are (case insensitive):  
jmap_type='CrossSec_Match'

# Decide if you want to include heterogeneous reactions in your mechanism or not: 
include_het_rxns=False    
    
# Set the path where your output files to use in F0AM will be stored: 
#output_dir='/uufs/chpc.utah.edu/common/home/u6044586/python_scripts/modules/GEOSChem_Emulator/Outpu_Files/10_31_Testing/'
output_dir='/uufs/chpc.utah.edu/common/home/u6044586/MATLAB/F0AM_v4.4.2/Chem/GEOSChem/v13-3-4/NEW/'

# Decide whether or not to overwrite files in the output_dir if they exist or to append a version #... 
overwrite=True

# Call the function that does it all! 
fjx= gcem.make_GC_mechanism(kppfile, rate_files, GC_version, jmap_type, include_het_rxns, output_dir, overwrite=overwrite)




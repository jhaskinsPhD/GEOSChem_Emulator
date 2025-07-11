import os 
import sys  

sys.path.append('/uufs/chpc.utah.edu/common/home/u6044586/python_scripts/modules/GEOSChem_Emulator_CLASS')
import GEOSChem_Emulator_dictv2_class as gcem

def list_directories(path):
    ls= [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]
    return  ls

def create_directory(path, directory_name):
    # Construct the full path for the new directory
    new_directory_path = os.path.join(path, directory_name)
    
    # Check if the directory already exists
    if not os.path.exists(new_directory_path):
        # Create the directory
        os.makedirs(new_directory_path)

    return new_directory_path

# Define the version of GEOS-Chem you want to emulate: 
GC_version='14.5.2';

# Define path where GEOS-Chem input KPP files are stored: 
inpt_pth='/uufs/chpc.utah.edu/common/home/u6044586/python_scripts/modules/GEOSChem_Emulator_CLASS/Input_Files/GC_Mech_Files'

###################################################################################
# DO FULL-CHEM FILES: 
####################################################################################
# Define the path to the KPP file of the version of GEOS-Chem you want to emulate: 
kppfile=os.path.join(inpt_pth,'v'+GC_version,'fullchem.eqn')

# Define the path to the gckPP Rates file(s) of the version of GEOS-Chem you want to emulate: 
rate_files=[os.path.join(inpt_pth,'v'+GC_version,'gckpp_Rates.F90'),
            os.path.join(inpt_pth,'v'+GC_version,'fullchem_RateLawFuncs.F90'),
            os.path.join(inpt_pth,'v'+GC_version,'rateLawUtilFuncs.F90')]

# Define how the j-value matching should be done to translate GEOS-Chem 
# photolysis rates to their corresponding MCM photolysis rate. Allowed options are (case insensitive):  
jmap_type='CrossSec_Match'

# Decide if you want to include heterogeneous reactions in your mechanism or not: 
include_het_rxns=False

# Set the path where your output files to use in F0AM will be stored: 
foam_pth='/uufs/chpc.utah.edu/common/home/u6044586/MATLAB/F0AM_v4.4.2/Chem/GEOSChem'
output_dir= create_directory(foam_pth,'v'+GC_version)

# Decide whether or not to overwrite files in the output_dir if they exist or to append a version #... 
overwrite=True

# Call the function that does it all! 
foam_dict,kpp_dict,inds2drop,fixed_vars, fjx_df, jmap_dict = gcem.make_GC_mechanism(kppfile, 
                                                                                    rate_files, 
                                                                                    GC_version, 
                                                                                    jmap_type, 
                                                                                    include_het_rxns, 
                                                                                    output_dir, 
                                                                                    overwrite=overwrite)


###################################################################################
# DO Hg Mechanism Files: 
####################################################################################
# Define the path to the KPP file of the version of GEOS-Chem you want to emulate: 
hg_GC_version='14.5.2-Hg'; include_het_rxns= True
hg_kppfile=os.path.join(inpt_pth,'v'+hg_GC_version,'Hg.eqn')

# Define the path to the gckPP Rates file(s) of the version of GEOS-Chem you want to emulate: 
hg_rate_files=[os.path.join(inpt_pth,'v'+hg_GC_version,'Hg_RateLawFuncs.F90')]

hg_output_dir= create_directory(foam_pth,'v'+hg_GC_version)
hg_foam_dict,hg_kpp_dict,hg_inds2drop,hg_fixed_vars, hg_fjx_df, hg_jmap_dict = gcem.make_GC_mechanism(hg_kppfile, 
                                                                                    hg_rate_files,
                                                                                    hg_GC_version, 
                                                                                    jmap_type, 
                                                                                    include_het_rxns, 
                                                                                    hg_output_dir, 
                                                                                    overwrite=True)

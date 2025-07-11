#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions utilized to parse a KPP file & return its contents. 

@author: Prof. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD
"""
import os
import sys
import yaml
import re
import pandas as pd
from collections import defaultdict

sys.path.append('/uufs/chpc.utah.edu/common/home/u6044586/python_scripts/modules/GEOSChem_Emulator_CLASS/')
import utils_class as utils
from mechanism_class import *
from prep_for_foam_class import *
import photolysis_class as hv 
# %%---------------------------------------------------------------------------
# L2-function used to peice together reactions that may have been split over multiple
# lines in the KPP file before line is passed to parse_rxns_and_rates() within read_kpp(). 
#-------------------------------------------------------------------------------

def _fix_multiline_rxns(line, line_stripped, count, kppfile):
    """ Function used to parse reactions that are listed on multiple lines in 
    the KPP file to connect them into a single line with the entire reaction/rate on one line. 
    These lines that should get sent to this function will look like this in the KPP file: 

        ln#     'Some text in KPP file line'
        41      'MVKOHOO + HO2 = 0.360MCO3 + 0.360GLYC +' 
        42      '  0.335MVKHP + 0.050MGLY + 0.050CH2O :  GCARR(2.12E-13, 0.0E+00, 1300.0); {2019/11/06; Bates2019; KHB}' 

    After this function completes, the outputted "line" should look like this so the rreaction parser will work! 
         'O3 + MO2 = CH2O + HO2 + O2 :   GCARR(2.90E-16, 0.0E+00, -1000.0);   {2014/02/03; Eastham2014; SDE}' 
                                     ^split #1                            ^ split #2  

    INPUTS: 
    -------

        (1)  line - String of the current KPP file line you're parsing for info. This line 
               will contain only half of the reaction and will end with a '+' because 
               that is what triggers this function being called! 

        (2)  line_stripped- String that is the same as "line", but has no spaces in it. 

        (3) count- Int containing the current line number in the file. Is the index of the 
               line that contains the first half of the reaction's stuff. 

        (4) kpp_file- String containing full path to the KPP input file that's being parsed. 
            Needed here so we can read the NEXT lines with info about this continued reaction. 


    OUTPUTS: 
    --------

        (1) line- Updated string that has all the info about this reaction in a single string. 
              It is crafted by concatonating continuing lines together, but does not contain any new line chars. 

        (2) line_stripped- Updated string containing the same info as updated "line", but has all spaces removed. 

        (3) count- Updated int containing line number of the last line you read in. Is the index of the line 
                in the KPP file where the split reaction ENDS! Updated so the next iteration doesn't 
                try to parse a line starting in the middle of this split reaction... 

    REQUIREMENTS: 
    -------------
        LIBRARIES:           None

        CUSTOM FUNCTIONS:    None

    USAGE: 
    ------  
        CALLED WITHIN:  read_kpp() via a call to make_GEOSCHEM_AllRxns_file() 

        OUTPUT USAGE: Updated output "line" next passed to parse_rxns_and_rates() within 
                        read_kpp() called within make_GEOSCHEM_AllRxns_file(). 

    AUTHOR: 
    -------
        Prof. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD

    """

    # Read more lines til you get the entire reaction, rate & comments in a single line!
    while line_stripped.endswith('+') == True or ':' not in line_stripped:
        nline = kppfile.readline()   # Read next line from file
        count = count+1                 # Update line counter variable
        # append this new line to our current line str.
        line = line + nline.strip()
        # Update line_stripped to not contain spaces!
        line_stripped = line.strip().replace(' ', '')

    # Now that you got the entire reaction on this line, make sure to replace 
    # any new line characters with a space!
    line = line.replace('\n', ' ')
    line_stripped = line_stripped.replace('\n', ' ')

    # Return the count, new line, and line_stripped vars to continue on.
    return line, line_stripped, count


# %%---------------------------------------------------------------------------
# L2-function used in read_kpp() to see if KPP line is blank or not. 
# -----------------------------------------------------------------------------

def _is_blank(line):
    """Function to determine if a line is blank in the KPP file or not after 
    all spaces have been removed from it."""
    
    tf = True if len(utils.str_multi_replace(line, [' ', '\n','\t'], rep_all=''))==0 else False
    return tf

# %%---------------------------------------------------------------------------
# L2-function used to parse the KPP lines with info on variable definitions in 
# the KPP file & update a dict with that info.  
#------------------------------------------------------------------------------

def _parse_var_defs(line:str, tracer_info:dict, fixed_vars:list, is_fixed:bool,
                    spell_defs_out:bool=True):
    """Function to parse variable defintion lines in the KPP File header & dump 
    var names/ information on vars into a dictionary. These lines look like: 
        
        'A3O2       = IGNORE; {CH3CH2CH2OO; Primary RO2 from C3H8}'  
        
    F'n should work as  expect as long as the varname only appears before the 
    first '=' sign in the string and var info appears within the outermost '{ }' 
    pair in the part of the string that comes after the first '='. 

    INPUTS: 
    -------
        (1) line        - STR of the current KPP file line you're parsing for info. 

        (2) tracer_info - DICT holding var names (as keys) and var info as values 
                          read in from previous lines to update with info from input line. 
                        
        (3) tracer_type - DICT with 2 keys 'varying' and 'fixed' pointing to lists 
                          of tracers of that type to update with info from input line
                        
        (4) is_fixed    - BOOL indicating if the variable is a fixed var or not. 
                          Typically assigned to value of file_locaiton['is_deffix']
                          in read_kpp(). 

    OUTPUTS:
    -------- 
       (1) tracer_info - Updated DICT with this line's new defined variable/info in it. 
        
       (2) fixed_vars - List of all fixed variables. 
 
    REQUIREMENTS: 
     -------------
        LIBRARIES:           import re 

        CUSTOM FUNCTIONS:    _simplifiy_tracer_dict()

    USAGE: 
    -----  
        CALLED WITHIN:  read_kpp() from within make_GEOSCHEM_AllRxns_file(). 

        OUTPUT USAGE:   Output dict "var_info" passed as output of 
                         read_kpp() and make_GEOSCHEM_AllRxns_file()

    AUTHOR: 
    -------
        Prof. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD
        """
    # Variable definition lines  in KPP files look like this: 
    #   "LIMO       = IGNORE; {C10H16; Limonene}  "
    #   "LIMO2      = IGNORE; {RO2 from LIMO}     "
    #   "APINP      = IGNORE;                     " 
    
    # Split the line at the first '=' to extract the variable name
    parts = line.split('=', 1); varname = parts[0].strip()
    
    # If there's no value after '=', assign an empty string to 'rest'
    rest = parts[1].strip() if len(parts) > 1 else ''
    
    # Check if there's additional information about the var on the rest of the 
    # lineby using regex to extract any content from within the outermost '{}' 
    match = re.search(r'{(.*)}', rest) # rest = "IGNORE; {C10H16; Limonene}  " 
    
    if match:
        # If you find something between two outer {} chars, then that's the var_info. 
        # so, remove any trailing/leading sapces fomr it. 
        var_info = match.group(1).strip()
    else:
        # If you didn't find those chars, then just set var info to an empy str. 
        var_info = ''
        
    # Add the variable name and info to the dictionary
    tracer_info[varname] = var_info
    
    # If this tracer has a fixed value & isn't in list, append to list of fixed vars.  
    if is_fixed and varname not in fixed_vars:
        fixed_vars.append(varname)
        
    return tracer_info, fixed_vars

# %% --------------------------------------------------------------------------
# Helper fn's used exclusively within read_kpp() & each other, to ID & update 
# what section of the KPP file is currently being read in. No external libraries 
# required for these, but do depend on one another. 
#------------------------------------------------------------------------------

# L-2 Function called within within read_kpp()
def _print_loc(file_loc, line): 
    """Helper function calledto print the file line & which file location keys 
    are currently set to true for debugging purposes"""
    
    # Get a list of the keys that are currently set to "TRUE" . 
    true_keys = '.'.join([key for key, value in file_loc.items() if value])
    
    # Print Line & which Key is set to tru: 
    print(f'{true_keys} ...  Line:"{line}"')
    
    return 

# L-3 function called within _update_loc_in_KPP()  
def _toggle_file_loc(file_loc:dict, turn_on:str=''):
    """Helper function called to turn one bool in the file_loc dictionary to True 
    & all the others to False since you can ONLY be in one section of the file at a time."""
    
    # Check that key to turn on is actually in the dictionary first
    if turn_on not in list(file_loc.keys()): 
        raise ValueError(f"Key '{turn_on}' not found in thefile_location dictionary."+\
                         f"Valid keys include: {','.join([k for k in list(file_loc.keys())])}")
            
    for key in list(file_loc.keys()):
        if key==turn_on: 
            file_loc[key]=True
        else: 
            file_loc[key]=False 
            
    return file_loc
    
# L-2 Function called within read_kpp() to track location in KPP file
def _update_loc_in_KPP(line:str, file_loc:dict,before_parsing:bool, ): 
    """Function to update bools in file_loc dict that tell you what section of a 
    KPP file you are on / entering.""" 
    
    # Strip all white space, & convert to lower case before comparing to loc strings
    line_mod = line.strip().replace(' ', '').lower()
    
    #--------------------------------------------------------------------------
    # Define lists of words that ALL must appear in a commented out line of the 
    # KPP file that would indicate we are entering a new section of the file:
    #--------------------------------------------------------------------------
    # "// %%%%% Reactions extracted from sulfate_mod.F90 (MSL, BMY)             %%%%%"
    sulf_mod_keywords= ['//', 'sulfate_mod', 'reactions']
    # "// %%%%% Gas-phase chemistry reactions                                   %%%%%"
    gas_rxns_keywords= ['//', 'gas', 'phase',  'reactions']
    # "// %%%%% Heterogeneous chemistry reactions                               %%%%%"
    het_rxns_keywords= ['//', 'heterogeneous', 'reactions']
    # "// %%%%% Photolysis reactions                                            %%%%%"
    photo_rx_keywords= ['//', 'photolysis', 'reactions']
                                                        
    #--------------------------------------------------------------------------
    # Look for indicaters in the current line that we're entering a new section. 
    # Turn the old section bool to FALSE and the new section bool to TRUE. 
    # If debug is set to true, it will also print the (un-modified line) that 
    # and what keys are true so you can see what is triggering the toggle switch.
    #--------------------------------------------------------------------------
    if '#defvar' in line_mod: 
        # When "#defvar" appears on the current line, you will no longer be in 
        # the KPP file header section, but will be in the variable definition 
        # section on the next line. 
        if before_parsing==True: 
            file_loc['in_header'] = False # Turn OFF before parsing line. 
        else: 
            file_loc['in_defvar'] = True # Turn ON after parsing line. 
        
    elif '#deffix' in line_mod: 
        # When "#defix" appears on the current line, you will no longer be in 
        # the variable definition section, but will enter the section where
        # *fixed* tracers names are declared section on the next line. 
        if before_parsing==True: 
            file_loc['in_defvar'] = False # Turn OFF before parsing line. 
        else: 
            file_loc['in_deffix'] = True # Turn ON after parsing line. 
        
    elif '#equations' in line_mod: 
        # You are no longer in section defining fixed vars if #equations appears 
        # in the current line! You will now be entering a sulf/gas/het/photo rxn 
        # section (but for max flexibility if the order the rxn sections are in 
        # changes in later iterations, we're not tunring anything on yet... 
        file_loc['in_deffix'] = False
        
    elif all([words in line_mod for words in sulf_mod_keywords]):
        # You will be in the sulfate_mod reaction section on the next line if
        # all these strings appear in the current line.
        file_loc=_toggle_file_loc(file_loc, turn_on='in_sulf_rxns')
        
    elif all([words in line_mod for words in gas_rxns_keywords]):
        # You will be in the gas-phase reaction section on the next line if 
        # all these strings appear in the current line.
        file_loc=_toggle_file_loc(file_loc, turn_on='in_gas_rxns')
        
    elif all([words in line_mod for words in het_rxns_keywords]):
        # You will be in the heterogenous reaction section on the next line if 
        #all these strings appear in the current line.
        file_loc=_toggle_file_loc(file_loc, turn_on='in_het_rxns')
        
    elif all([words in line_mod for words in photo_rx_keywords]):
        # You will be in the photolysis reaction section on the next line if 
        # all these strings appear in the current line.
        file_loc=_toggle_file_loc(file_loc, turn_on='in_photo_rxns')
    
    return file_loc
                             
# L2-Function Called in read_kpp()
def _simplify_tracer_dict(input_dict): 
    """Function to not allow circular refs to items in tracer_dict"""
    resolved_dict = {}
    for key, value in input_dict.items():
        # Split into formula and description if ';' exists
        if ';' in value:
            parts = value.split(';', 1)
            formula = parts[0].strip()
            description = parts[1].strip()
        else:
            formula = ''
            description = value.strip()

        pattern = r'\bfrom\s+(\w+)\b'
        for match in re.finditer(pattern, description):
            ref_key_candidate = match.group(1)
            # ONLY proceed if ref_key_candidate is actually a key in input_dict
            if ref_key_candidate in input_dict:
                ref_value = input_dict[ref_key_candidate]
                ref_parts = ref_value.split(';', 1)
                ref_formula = ref_parts[0].strip() if ref_parts else ''
                ref_name = ref_parts[1].strip() if len(ref_parts) > 1 else ref_key_candidate

                original_phrase = match.group(0)
                replacement = f'{original_phrase} (e.g. {ref_name}/{ref_formula})'
                description = description.replace(original_phrase, replacement, 1)

        # Reconstruct full string
        if formula:
            resolved_entry = f"{formula}; {description}".strip('; ')
        else:
            resolved_entry = description

        resolved_dict[key] = resolved_entry

    return resolved_dict
# %% ##########################################################################
# read_kpp() is the main wrapper function to read in a KPP file and output its 
# contents. Called within make_GEOSCHEM_AllRxns_file(). All helper functions 
# exclusively called by read_kpp() are defined herein. But it is also uses
# dependent on some functions defined in utils.py which are also used in 
# other portions of the GC Emulator & are NOT exclusive to read_kpp().
###############################################################################

def read_kpp(kppfile,debug:bool=False, return_df:bool=False, out2yaml:bool=False, 
             out2excel:bool=False, output_dir:str='',output_fname:str='kpp_mech',
             overwrite:bool=False, spell_defs_out:bool=True):
    """Function to parse a KPP file and extract info about tracers, reactions,
    & rates in a pandas dataframe. Output can also be provided in other ways using 
    keywords. 

    INPUTS:  
    -------  
      (1) kppfile - STR containing the full path and name of a KPP file.  
                    
      (2) debug - (OPTIONAL) BOOL for printing additional output for debugging. 
                  Default is FALSE. 
                  
      (3) return_df - (OPTIONAL) BOOL indicating if output kpp_dict, should be 
                  returned as a pandas dataframe rather than dictionary (default). 
                  Default is  set to FALSE, as dictionary is easier for parsing later, 
                  but as df, output is much more read-able for humans...  

      (3) out2yaml - (OPTIONAL) BOOL indicating if outputs should be saved to a 
                     yaml file. Default is FALSE.  
       
      (4) out2excel - (OPTIONAL) BOOL indicating if outputs shoudl be saved to 
                      an excel .xlsx sheet. Default is FALSE. 

       ------ below only used if either out2yaml or out2excel is true ---------
          
      (5) output_dir - (OPTIONAL) STR containing the path where the output files
                        should be saved.
      
      (6) output_fname - (OPTIONAL) STR with name of the output file to write. 
                         Default is 'kpp_mech.{ext}'
      
      (7) overwrite - (OPTIONAL) BOOL indicating whether or not output files 
                        should overwrite any exisiting files at output_dir
                        with output_fname or not. 
       
    OUTPUTS:  
    -------- 
       (1) tracer_info - DICT with keys as tracer names, values = long name/ info
       
       (2) tracer_type - DICT with keys 'varying' and 'fixed' pointing to lists 
                         of tracers that are variables or fixed. 
       
       (3) kpp_dict - DICT with info on all reactions from KPP file.  


    AUTHOR: 
    -------
       Prof. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD
    """
    ###########################################################################
    # Initialize all vars to hold outputs/ keep track of where we are!
    ###########################################################################

    # Get # of lines to expect so we know when we've reached the end of the file!
    num_lines = sum(1 for line in open(kppfile))
    
    # Initialize variable used to count the line number being parsed that's 
    # updated in the while loop when we start to read the file. 
    count = 0  
    
    # Initialize dict to hold the tracer names (keys) & info (value) declared in 
    # KPP file and list of all fixed tracers
    tracer_info={}; fixed_vars= []

    # Define boolean dictionary with key/values that will tell us what part of 
    # the KPP file we are currently parsing.
    file_location = dict({
        'in_header': True,     # True oat file start & before vars defined. 
        'in_defvar': False,    # True if parsing the *variable* tracer name definitions
        'in_deffix': False,    # True if parsing the *fixed* tracer name definitions 
        'in_sulf_rxns': False, # True if parsing reactions from sulfate_mod.F (het)
        'in_gas_rxns': False,  # True if parsing gas phase reactions (not photo/het)
        'in_het_rxns': False,  # True if parsing heterogenous reactions
        'in_photo_rxns': False # True if parsing photolysis reactions
    })

    # Initialize output as empty dict- Contents assigned in _parse_rxns_and_rates(). 
    kpp_mech=kpp_mechanism(); 
    
    ##############################################################################
    # Open the KPP file in "read" mode & loop over all lines until end is reached!
    ##############################################################################
    file = open(kppfile, 'r')
    while True:  # loop over each line of the KPP file.
        line = file.readline()  # Read next line from file

        # Strip all white space from current line.
        line_stripped = line.strip().replace(' ', '')
        
        # Update line counter variable now that you read in this line.
        count = count+1

        # ---------------------------------------------------------------------
        # BEFORE parsing the line, set file location booleans to FALSE if the 
        # current line indicates you're entering a new section!
        # ---------------------------------------------------------------------
        file_location= _update_loc_in_KPP(line, file_location,before_parsing=True)
        
        if debug is True and _is_blank(line)==False: _print_loc(file_location, line)
        
        # Break if reached end of file, don't attempt to parse line.
        if count > num_lines:
            break
        # ---------------------------------------------------------------------
        # If you're in the part of the file where variable definitions are made
        # & tracer names are listed then pass to ==> parse_var_defs()
        # ---------------------------------------------------------------------
        if (((file_location['in_defvar'] == True) or (file_location['in_deffix'] == True))
            and (_is_blank(line) == False)):

            # Pass (un-stripped) line to parse_var_defs() & return updated version of 
            # the output dictionaries, tracer_info that contains varnames, defitions, 
            # the outlut list of fixed tracer names.
            
            tracer_info, fixed_vars= _parse_var_defs(line, tracer_info,fixed_vars, 
                                        is_fixed=file_location['in_deffix'])
            
        # ---------------------------------------------------------------------
        # If you're in part of the file with chemical equations & rates then 
        # pass it to ==> parse_rxns_and_rates() to add info to the output kpp_dict. 
        # ---------------------------------------------------------------------
        if any(cond == True for cond in [file_location['in_sulf_rxns'], 
                                         file_location['in_gas_rxns'], 
                                         file_location['in_photo_rxns'],
                                         file_location['in_het_rxns']]):

            # Don't parse commented out or blank lines.
            if (('//' != line[0:2]) and (_is_blank(line) == False)):

                # Decide if you've been given a line that contains a partial reaction 
                # & fix it before parsing as rxn.
                if ':' not in line:
                    # These lines won't contain a ':' and look like this:
                    #   ln #     Content of Line
                    # -------------------------------------------------------------
                    #   ln 42   'MVKOHOO + HO2 = 0.360MCO3 + 0.360GLYC +'   **(Would trigger when parsing this line)**
                    #   ln 43   '  0.335MVKHP + 0.050MGLY + 0.050CH2O :  GCARR(2.12E-13, 0.0E+00, 1300.0); {2019/11/06; Bates2019; KHB}'

                    line, line_stripped, count = _fix_multiline_rxns(line, line_stripped, count, file)
                    
                # Figure out reaction type based on location in file: 
                is_gas= file_location['in_gas_rxns']
                is_het=True if any ([file_location['in_het_rxns'],file_location['in_sulf_rxns']]) else False 
                is_photo= file_location['in_photo_rxns'] 
                
                # (After fixing if necessary), create a kpp_reaction object from this line: 
                rxn_obj=kpp_reaction(line=line, is_gas=is_gas,is_het=is_het, 
                                     is_photo=is_photo, rxn_arrow='=', 
                                     rcts2ignore=['hv','M'], fixed_vars=fixed_vars)
                
                # And append this reaction oject to the mechanism object: 
                kpp_mech.add_rxn(rxn_obj=rxn_obj)

            
        # ----------------------------------------------------------------------
        # AFTER parsing the current line, set file location booleans to TRUE if 
        # the current line indicates the NEXT line will be in a new section!
        # ----------------------------------------------------------------------
        file_location= _update_loc_in_KPP(line, file_location, before_parsing=False)
        
    # While loop will eventually break once we found the end of the KPP file.
    
    file.close() # Close the KPP file.
    
    # if out2yaml==True: 
    #     # Check that the filename is available & return the full path to write the file: 
    #     yaml_file=utils.check_filename(filename=output_fname, default_name= 'kpp_mech', ext='.yml', 
    #                        savepath=output_dir, overwrite=overwrite, return_full=True, quiet=True)
    #     # Write the dictionary to a YAML file
    #     with open(yaml_file, 'w') as f:   
    #         yaml.dump(kpp_dict, f)
    #     print(f'Output of read_kpp() saved at: \n\t{yaml_file}')
    
    # if return_df== True or out2excel== True: 
    #     # First convert from a record based dict to an index based on so df 
    #     # formatting makes sense... 
    #     kpp_dict = utils._convert_to_category_dict(kpp_dict) 
        
    #     # Convert dict to dataframe with all this information... 
    #     kpp_df = pd.DataFrame(kpp_dict)
        
    # if out2excel==True: 
    #     # Check that the filename is available & return the full path to write the file: 
    #     excel_file=utils.check_filename(filename=output_fname, default_name= 'KPP_Mechanism', ext='.xlsx', 
    #                        savepath=output_dir, overwrite=overwrite, return_full=True, quiet=True)
    #     # Write the dataframe to an excel file: 
    #     kpp_df.to_excel(excel_file)
    #     print(f'Output of read_kpp() saved at: \n\t{excel_file}')
    
        
    # if spell_defs_out is True: 
    #     tracer_info=dict(sorted(tracer_info.items(), key=lambda item: len(item[0]), reverse=True))
    #     tracer_info=_simplify_tracer_dict(tracer_info)
    
    # if return_df is False: 
    #     return tracer_info, fixed_vars, kpp_dict 
    # else: 
    #     return tracer_info, fixed_vars, kpp_df
    
    return tracer_info, fixed_vars, kpp_mech





###############################################################################
##  Standalone debugging of read_kpp_mod.py: 
###############################################################################
outdir='/uufs/chpc.utah.edu/common/home/u6044586/MATLAB/F0AM_v4.4.2/Chem/GEOSChem/v14.5.2'
kppfile='/uufs/chpc.utah.edu/common/home/u6044586/python_scripts/modules/GEOSChem_Emulator/Input_Files/GC_Mech_Files/v14.5.2/fullchem.eqn'

# tracer_info, fixed_vars, kpp_mech = read_kpp(kppfile, output_dir=outdir,
#                                       return_df=True, debug=False)


GC_version='14.5.2';
inpt_pth='/uufs/chpc.utah.edu/common/home/u6044586/python_scripts/modules/GEOSChem_Emulator/Input_Files/GC_Mech_Files'

# Define the path to the gckPP Rates file(s) of the version of GEOS-Chem you want to emulate: 
rate_files=[os.path.join(inpt_pth,'v'+GC_version,'gckpp_Rates.F90'),
            os.path.join(inpt_pth,'v'+GC_version,'fullchem_RateLawFuncs.F90'),
            os.path.join(inpt_pth,'v'+GC_version,'rateLawUtilFuncs.F90')]


# Store all the user input in a nice verified dictionary (passed to many functions)
user_inp={'kppfile':kppfile, 
          'rate_files':rate_files,
          'GC_version': GC_version,
          'jmap_type': 'CrossSec_Match',
          'include_het_rxns':False,
          'input_paths': utils.get_my_relative_paths(), 
          'output_dir':outdir,
          'verbose':False}

# Read in the right FJX file for this GC-Version and build the j-mapping
# dictionary & FJX dataframes to use in mapping KPP J's to F0AM J's.
jmap_dict, fjx_df, user_inp = hv.create_jmap_dict(user_inp, tracer_info)


# %%

lines=['NO3 + hv = NO2 + O :     PHOTOL(12);     {2014/02/03; Eastham2014; SDE}',
        'PPN = RCO3 + NO2 :       GCJPLEQ_acabab(9.00d-29, 14000.0d0, 9.00d-28, 8.9d0, 7.7d-12, 0.2d0, 0.6d0);',
        'ClO + ClO = Cl2 + O2 :   GCARR_ac(1.00d-12, -1590.0d0);       {2014/02/03; Eastham2014; SDE}',
        'RCO3 + NO2 {+M} = PPN :        GCJPLPR_abab(9.00d-28, 8.9d0, 7.7d-12, 0.2d0, 0.6d0);       {JPL Eval 17}',
        'ClOO {+M} = Cl + O2 {+M} :     GCJPLEQ_acabab(6.60d-25, 2502.0d0, 2.20d-33, 3.1d+00, 1.8d-10, 0.0d0, 0.6d0);   {JPL 15-10; XW}',
        'O + O2 {+M} = O3 {+M} :     GCARR_ab(6.00d-34, 2.4d0)*NUMDEN;        {2014/02/03; Eastham2014; SDE}'
        ]

rs=[]
is_photo=[True, False, False, False, False, False] 
for i in range(0, len(lines)): 
    rxn_obj=kpp_reaction(line=lines[i], is_gas=True,is_het=False, 
                          is_photo=is_photo[i], rxn_arrow='=', 
                          rcts2ignore=['hv','M'],prds2ignore=['hv','M'],
                          fixed_vars=['H2','O2','N2'])
    
    
    f=foam_reaction(rxn_obj, jmap_dict=jmap_dict, ignore=['hv','M'])
    f.disp(rxn_obj=rxn_obj) 
    print('\n\n')

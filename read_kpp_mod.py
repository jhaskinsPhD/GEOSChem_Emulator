#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions utilized to parse a KPP file & return its contents. 

@author: Prof. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD
"""

import sys
import yaml
import re
import pandas as pd
from collections import defaultdict

sys.path.append('/uufs/chpc.utah.edu/common/home/u6044586/python_scripts/modules/GEOSChem_Emulator/')
import utils


# %%---------------------------------------------------------------------------
# Functions called within read_kpp() to parse lines chemical reactions & 
# extract contents about a reaction into an output dictionary. 
#------------------------------------------------------------------------------

# L4-function used in parse_reaction_str() to seperate stoichiometry from a 
# tracer name in a single input string. 

def _sep_stoich_n_tracer(pair:str): 
    """Helper function to take a single input string & seperate it based on where the 
    first alpha numeirc character appears.Anything before first letter 
    is assumed to be stoichimetry. Anthing after assumed to be tracer name. 
    If first char is letter, stoichiometery assumed to be 1.
    Called in parse_reaction_str() """ 
    
    ind=0; # Initialize index counting where we're at in input string. 
    
    # First char is always the sign. 
    sign=pair[0]; pair=pair[1:]
    while True:  # Loop through string from left to right til you find a letter. 
    
        # If this char is a letter OR you reached the end of the string ... 
        if pair[ind].isalpha() or  ind==len(pair)-1: 
            
            # Then assume the tracer name begins here & fills rest of str. 
            tracer=pair[ind:]
            
            if len(pair[:ind].replace(' ',''))==0: 
                # If nothing is before the compound name, then the stoichiometry 
                # must be 1. 
                stoich=float(1) 
            else: 
                # Otherwise, the stoichiometry is what comes before the tracer name. 
                stoich=float(pair[:ind])     
                
            break # You got cmpd & stoich so break out of while loop, move to 
                   # the next name/stoich pair. 
                       
        ind=ind+1 # Update index for next ieration if still doing while loop. 
    
    #Update stoich to be netative if it is... 
    stoich=-1*stoich if sign=='-' else stoich
              
    return stoich, tracer 
        
# L4-function used in parse_indv_reaction_str() to get list of stuff between any
# set of input delimiter. 

def _get_grps_between_delims(s:str):
    """Helper function to seperate a string/ its sign along any input delimiters in the input
    list called in parse_reaction_str()"""
    delimiters=['+','-']
    
    # Join the delimiters list into a regex pattern string
    delimiters_str = '|'.join(map(re.escape, ['+','-']))  # Escape delimiters for regex
    
    # Using regex to split but capture the delimiters
    pattern = f'(?<=[{re.escape(delimiters_str)}])|(?=[{re.escape(delimiters_str)}])'
    parts = re.split(pattern, s)
    
    # Assume delimter is '+' if no delim preceeds the group. 
    current_delimiter = '+'
    
    result=[]
    for part in parts:
        if part in delimiters:
            current_delimiter = part
        elif part:
            result.append(current_delimiter+part)
    
    return result
    
# L3-function used in _parse_rxns_and_rates() to get list of products/reactants 
# & stoichiometry of a input reaction string. 

def parse_indv_reaction_str(rxn_str:str, rxn_arrow:str='=',
                       rcts2ignore:list=[], prds2ignore:list=[]):
    """Function to parse a chemical reaction string and seperate the stoichiometry 
       from the tracer names for reactants/products seperately. NOTE: Function 
       won't work properly if tracer names can begin with numbers. 
    
    Inputs: 
    -------
        rxn_str - A STRING you'd like to parse to seperate into reactant & product 
                  name/stoichioemtry pairs like: 
                     'TOLU + OH = TRO2 + 1.920CH2O + 0.260GLYX + 0.215MGLY + OH'
                     
        rxn_arrow - A STR that indciates what char sequence seperates the reactants 
                     from products in the input rxn strings. Default is set to '=' 
                     as it is in KPP files, but may be ->', '-->', '<->' in other 
                     mechanism types. 
        
        pair_delims - A list of STRS that are allowed to spereate tracer name/
                       stoichiometry pairs from one another in rxn.This list should 
                       NOT include the delimiter for the reaction arrow. NOTE: 
                       Some mechanisms allow '-' instead of just '+' / why it
                       can be a list. 
                       
        rcts2ignore - A list of STRS to ignore in output lists if they are parsed 
                      as tracernames in the reactants. Can be useful if you don't 
                      want to parse 'O3+hv=O1D' as having a reactant named 'hv' 
                      with stoich=1. 
                      
        prds2ignore - A list of STRS to ignore in output lists if they are parsed 
                         as tracer names in the products. Can be useful if you don't 
                         want to actually track 'H2O' or 'O2' explicitily.           
                      
    Outputs: 
    -------
        Dictionary containing input reaction, & lists  of reactants/products & stoichiometry.
        
            out ['reaction']   = 'TOLU + OH = TRO2 + 1.920CH2O + 0.260GLYX + 0.215MGLY + OH'
                ['rct_cmpds']  = ['TOLU','OH']
                ['rct_stoich'] = [1.0,1.0]
                ['prd_cmpds']  = ['TRO2' 'CH2O', 'GLYX', 'MGLY', 'OH']
                ['prd_stoich'] = [1.0, 1.920 ,0.260, 0.215, 1.0]      

    REQUIREMENTS: 
    -------------
        LIBRARIES:          None 

        CUSTOM FUNCTIONS:    
                           _get_grps_between_delims() 
                           _sep_stoich_n_tracer
                            _convert_fortran_nums() 

    USAGE: 
    -----  
        Called Within: parse_rxns_and_rates() from  read_kpp() within
                      make_GEOSCHEM_AllRxns_file(). 
                      
        Output Usage: Contents used to fill the contents of kpp_dict for 
                    ['rct_cmpds','rct_stoich', 'prd_cmpds', prd_stoich']

    AUTHOR: 
    -------
        Prof. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD
    
    """
    # Ensure the input string has no spaces: 
    rxn_str=rxn_str.replace(' ', '') # rxn_str='TOLU+OH=TRO2+1.920CH2O+0.260GLYX+0.215MGLY+OH'
    
    # If rxn string contains the expected delimiter: 
    if rxn_arrow in rxn_str: 
        # Split reaction string on the reaction arrow delimiter first: 
        rct_half,prd_half =rxn_str.split(rxn_arrow)
        # rct_half='TOLU+OH'   &  prd_half='=TRO2+1.920CH2O+0.260GLYX+0.215MGLY+OH'
    else: 
        raise ValueError('Reaction arrow delimiter "{rxn_arrow}" could not be '+\
                         'found in the following reaction: \n\t{rxn_str}')
        
    # Seperate the two halves of the rxn along any pairs of delimiters inputted
    # giving us potential name/stoich pairs (e.g. Extract a list of whatever is 
    # between '+' or '-' in 'A+B-C+D' style strings: 
        
    #  Get list of reactant name &  stoich pairs & if positive/negative! 
    rct_pairs= _get_grps_between_delims(rct_half) # ['TOLU','OH']
    
    #  Get list of product  name &  stoich pairs & if positive/negative! 
    prd_pairs= _get_grps_between_delims(prd_half)# [ 'TRO2', '1.920CH2O','0.260GLYX','0.215MGLY','OH'] 
                
    # Define output dictionary: 
    out={'reaction':rxn_str, 'rct_cmpds':[], 'rct_stoich':[], 'prd_cmpds':[], 'prd_stoich':[]} 
    
    # Loop over each tracer/stoich string in the reactant half of the string: 
    for pair in rct_pairs: 
        # Seperate the stoich from tracer in this individual pair
        stoich_i, tracer_i= _sep_stoich_n_tracer(pair) 
        # But only add this tracer & its stoich if its NOT in the list of rcts to ignore... 
        if tracer_i not in rcts2ignore: 
            out['rct_cmpds'].append(tracer_i)
            out['rct_stoich'].append(stoich_i)    
        
    #Loop over each tracer/stoich string in the reactant half of the string: 
    for pair in prd_pairs: 
        # Seperate the stoich from tracer in this individual pair
        stoich_i, tracer_i= _sep_stoich_n_tracer(pair) 
        # But only add this tracer & its stoich if its NOT in the list of prds to ignore... 
        if tracer_i not in prds2ignore:
            out['prd_cmpds'].append(tracer_i)
            out['prd_stoich'].append(stoich_i)   

        
    return out 

# L2-function used in read_kpp() to extract info from a line in a KPP file & return 
# updated dict with info about this specific reaction line.  Called from read_kpp()
# within make_GEOSCHEM_AllRxns_file().

def _parse_rxns_and_rates(line: str, kpp_dict: dict, file_location: dict):
    """Function used to parse reaction lines in KPP file and extract their info 
    into dictionary, kpp_dict. 

    INPUTS:
    -------

        (1) line - STR of line with reaction, rate, and (maybe) a comment. These lines look like this:     
                      Reaction               Rate                              Comment
                "O3 + NO = NO2 + O2 :  GCARR_ac(3.00d-12, -1500.0d0); {2014/02/03; Eastham2014; SDE}"

        (2)  kpp_dict - (Nested) DICT. Top level key is reaction index (int) which 
                         points to a DICT for each entry that has the following 
                         key/value pairs: 
                             
                KEYS           :           VALUES
                ----------------------------------------------------------------------
                'reaction'            : STR with reaction 
                'rate'           : STR with full rate (func & arg) for this reaction 
                'comments'       : STR with comments about this reaction (or blank). 
                'is_gas'         : BOOL of whether rxn is gas phase or not. 
                'is_het'         : BOOL of whether rxn is heterogeneous or not. 
                'is_photo'       : BOOL of whether rxn is a photolysis rxn or not. 
                'is_3body'       : BOOL of whether rxn includes 3rd body reactant, M, or not. 
                'rct_cmpds'      : LIST of all reactants in this reaction.  
                'rct_stoich'     : LIST of stoiciometry of reactants in reaction (same len as rct_cmpds)
                'prd_cmpds'      : LIST of all products in this reaction. 
                'prd_stoich'     : LIST of stoiciometry of reactants in reaction (same len as prd_cmpds)
                'KPP_line'       : STR with original line from KPP file
                'rate_function'  : STR with rate function (only- no args or constants) used. 

        (3) file_location - DICT indicating current location in KPP file ... 
        
     OUTPUTS: 
     -------   
         (1)  kpp_dict    - Updated DICT with items from this line added to it. 


    REQUIREMENTS: 
    -------------
        LIBRARIES:           import re 

        CUSTOM FUNCTIONS:    
                            utils._str_multi_replace()
                            utils._convert_fortran_numbers()
                            utils._convert_fortran_operators
                            parse_indv_reaction_str()


    USAGE: 
    -----  
        Called Within: read_kpp() within make_GEOSCHEM_AllRxns_file(). 
        Output Usage: Output dict 'kpp_dict' used in read_kpp() to sep info for dif rxn types into own dicts! 

    AUTHOR: 
    -------
        Prof. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD

    """
    # Clean up line to maintain in output: 
    clean_line= utils.str_multi_replace(line, ['\n', '\t'], rep_all='').strip()
    
    #####################################################################################
    # Split Line in to reaction, rate, rate function, and comment & Clean up for MATLAB!
    #####################################################################################
    # Split on colon to get rxn only
    rxn_parts = line.split(':')  
    
    # Split on ';' to get rates, remove trailing/tailing spaces
    rate = rxn_parts[1].split(';')[0].strip() 
        
    # Define Booleans that tell us what type of rxn it is: 
    rxn_with_M=True if '{+M}' in rxn_parts[0] else False  
    rxn_is_gas= file_location['in_gas_rxns']
    rxn_is_photo= file_location['in_photo_rxns']
    rxn_is_het=True if any ([file_location['in_het_rxns'],file_location['in_sulf_rxns']]) else False 

    # If this part of the line after the reaction has a comment, then split that 
    # on'{' and clean it up.
    if '{' in rxn_parts[1]:
        replace_me = dict({'{': '', '}': '', ';': ',', '\n': '','\t':''})
        comment = utils.str_multi_replace(rxn_parts[1].split('{')[1], replace_me)
    else:
        comment = ''  # Otherwise make comment blank!

    # Update comment to reflect whether or not the rxn was a 3 body rxn with M or not.
    if rxn_with_M == True: comment= '(3-Body rxn with +M) '+comment
    
    # Now remove spaces, and '{+M'}, '{' , '}' characters from the reaction
    rxn = utils.str_multi_replace(rxn_parts[0], ['{+M}', '{', '}', ' '], rep_all='').strip()
 
    # Split entire reaction into Reactants & Products & their stoichiometry:
    sdict= parse_indv_reaction_str(rxn, rxn_arrow='=', 
                           rcts2ignore=['hv'], prds2ignore=[])
    
    # Clean up rates so the #s no longer have FORTRAN style numbers / operators.
    # NOTE: This does NOT convert multiplicaiton/division to MATALB style, just 
    # makes the syntax readable / in modern notation. 
    rate=utils._convert_fortran_numbers(rate)
    rate=utils._convert_fortran_operators(rate)

    # Get function name only from the rate constant defs/ The lines look like this: 
    # 'GCARR_ac(2.60e-13, 1300)' but some maybe sums of calls to multiple functions 
    # like so: GCARR_ac(2.60e-13, 1300) + GCARR_ac(2.60e-13, 1300) and some may 
    # contain exp(), log(), or even log10()... So use cray regex: 
    pattern = r'\b((?!exp\b|log\b|log10\b)[A-Za-z_]+[a-zA-Z0-9_]*)\('

    # Use a list comprehension to do regex search  & return get just the 
    # function(s) names referenced in the rated defintion in a list. 
    functions = [match_i[:] for match_i in re.findall(pattern, rate)]     
        
    # ------------------------------------------------------------------------- 
    # Now save everything into dict where a generic index is the key: 
    # ------------------------------------------------------------------------- 
    
    # Figure out what indexes are already taken / filled out: 
    inds=list(kpp_dict.keys())
    
    # If empty, start at 0, otherwise, add 1 to last value. 
    rxn_ind=0 if len(inds) ==0 else inds[-1]+1

    # Add & populate nested dict with info about this specific reaction to kpp_dict. 
    kpp_dict[rxn_ind]= dict({
            'reaction': rxn,         # STR with reaction 
            'rate' : rate,           # STR with rate (func & arg) for this rxn
            'comments': comment,     # STR with comments about this rxn (or blank). 
            'is_gas': rxn_is_gas,    # BOOL of whether rxn is gas phase or not. 
            'is_het': rxn_is_het,    # BOOL of whether rxn is heterogeneous or not. 
            'is_photo': rxn_is_photo,# BOOL of whether rxn is a photolysis rxn or not. 
            'is_3body': rxn_with_M,  # BOOL of whether rxn includes 3rd body reactant, M, or not. 
            'rct_cmpds':  sdict['rct_cmpds'],  # LIST of all reactants in this reaction.  
            'rct_stoich': sdict['rct_stoich'], # LIST of rct stoich in rxn (same len as rct_cmpds)
            'prd_cmpds':  sdict['prd_cmpds'],  # LIST of all products in this reaction. 
            'prd_stoich': sdict['prd_stoich'], # LIST of prd stoich in rxn (same len as prd_cmpds)
            'KPP_line': clean_line,     # STR with original line from KPP file
            'rate_function': functions, # LIST of STR(s) with rate function(s) used in rxn. 
            })
                
    return kpp_dict


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
    kpp_dict={}; 
    
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

                # (After fixing if necessary), parse this reaction and add all info about this rxn to lists:
                kpp_dict = _parse_rxns_and_rates(line, kpp_dict, file_location)
            
        # ----------------------------------------------------------------------
        # AFTER parsing the current line, set file location booleans to TRUE if 
        # the current line indicates the NEXT line will be in a new section!
        # ----------------------------------------------------------------------
        file_location= _update_loc_in_KPP(line, file_location, before_parsing=False)
        
    # While loop will eventually break once we found the end of the KPP file.
    
    file.close() # Close the KPP file.
    
    if out2yaml==True: 
        # Check that the filename is available & return the full path to write the file: 
        yaml_file=utils.check_filename(filename=output_fname, default_name= 'kpp_mech', ext='.yml', 
                           savepath=output_dir, overwrite=overwrite, return_full=True, quiet=True)
        # Write the dictionary to a YAML file
        with open(yaml_file, 'w') as f:   
            yaml.dump(kpp_dict, f)
        print(f'Output of read_kpp() saved at: \n\t{yaml_file}')
    
    if return_df== True or out2excel== True: 
        # First convert from a record based dict to an index based on so df 
        # formatting makes sense... 
        kpp_dict = utils._convert_to_category_dict(kpp_dict) 
        
        # Convert dict to dataframe with all this information... 
        kpp_df = pd.DataFrame(kpp_dict)
        
    if out2excel==True: 
        # Check that the filename is available & return the full path to write the file: 
        excel_file=utils.check_filename(filename=output_fname, default_name= 'KPP_Mechanism', ext='.xlsx', 
                           savepath=output_dir, overwrite=overwrite, return_full=True, quiet=True)
        # Write the dataframe to an excel file: 
        kpp_df.to_excel(excel_file)
        print(f'Output of read_kpp() saved at: \n\t{excel_file}')
    
        
    if spell_defs_out is True: 
        tracer_info=dict(sorted(tracer_info.items(), key=lambda item: len(item[0]), reverse=True))
        tracer_info=_simplify_tracer_dict(tracer_info)
    
    if return_df is False: 
        return tracer_info, fixed_vars, kpp_dict 
    else: 
        return tracer_info, fixed_vars, kpp_df





###############################################################################
##  Standalone debugging of read_kpp_mod.py: 
###############################################################################
# outdir='/uufs/chpc.utah.edu/common/home/u6044586/MATLAB/F0AM_v4.4.2/Chem/GEOSChem/v14.5.2'
# kppfile='/uufs/chpc.utah.edu/common/home/u6044586/python_scripts/modules/GEOSChem_Emulator/Input_Files/GC_Mech_Files/v14.5.2/fullchem.eqn'

# tracer_info, fixed_vars, kpp_dict = read_kpp(kppfile, output_dir=outdir,
#                                      return_df=True, debug=True)

###############################################################################
##  OLD VERSION FROM GEOSCHEM_EMULATOR_dictv.py
###############################################################################
    
# # %%---------------------------------------------------------------------------
# # -------------- (A.2) parse_kpp_main() & Sub-Functions -----------------------
# # -----------------------------------------------------------------------------
# # A.2.3.2.2 L4-function used in parse_reaction_str() called from parse_rxns_and_rates()
# #within parse_rxns_and_rates() from within parse_kpp_main() from within make_GEOSCHEM_AllRxns_file():
# def _sep_stoich_n_tracer(pair:str): 
#     """Helper function to take a single input string & seperate it based on where the 
#     first alpha numeirc character appears.Anything before first letter 
#     is assumed to be stoichimetry. Anthing after assumed to be tracer name. 
#     If first char is letter, stoichiometery assumed to be 1.
#     Called in parse_reaction_str() """ 
    
#     ind=0; # Initialize index counting where we're at in input string. 
    
#     while True:  # Loop through string from left to right til you find a letter. 
    
#         # If this char is a letter OR you reached the end of the string ... 
#         if pair[ind].isalpha() or  ind==len(pair)-1: 
            
#             # Then assume the tracer name begins here & fills rest of str. 
#             tracer=pair[ind:]
            
#             if len(pair[:ind].replace(' ',''))==0: 
#                 # If nothing is before the compound name, then the stoichiometry 
#                 # must be 1. 
#                 stoich=float(1) 
#             else: 
#                 # Otherwise, the stoichiometry is what comes before the tracer name. 
#                 stoich=float(pair[:ind])     
                
#             break # You got cmpd & stoich so break out of while loop, move to 
#                    # the next name/stoich pair. 
                       
#         ind=ind+1 # Update index for next ieration if still doing while loop. 
    
#     return stoich, tracer 
        
# # A.2.3.2.1 L4-function used in parse_reaction_str() called from parse_rxns_and_rates()
# #within parse_rxns_and_rates() from within parse_kpp_main() from within make_GEOSCHEM_AllRxns_file():
# def _get_grps_between_delims(s:str, delimiters:list):
#     """Helper function to seperate a string along any input delimiters in the input
#     list called in parse_reaction_str()"""
    
#     # Join the delimiters list into a string (big regex OR statement). 
#     delimiters_str = ''.join(delimiters)
    
#     # Escape any special regex characters
#     pattern = rf'[^{re.escape(delimiters_str)}]+'
    
#     # Find all sequences of characters not in delimiters
#     groups = re.findall(pattern, s)
    
#     # Don't return any empty groups: 
#     return [g for g in groups if g]

# # A.2.3.2 L3-function used in parse_rxns_and_rates() called within parse_rxns_and_rates() from 
# # within parse_kpp_main() from within make_GEOSCHEM_AllRxns_file():
# def parse_reaction_str(rxn_str:str, rxn_arrow:str='=', pair_delims:list=['+','-'],
#                        rcts2ignore:list=[], prds2ignore:list=[]):
#     """Function to parse a chemical reaction string and seperate the stoichiometry 
#        from the tracer names for reactants/products seperately. NOTE: Function 
#        won't work properly if tracer names can begin with numbers. 
    
#     Inputs: 
#     -------
#         rxn_str - A STRING you'd like to parse to seperate into reactant & product 
#                   name/stoichioemtry pairs like: 
#                      'TOLU + OH = TRO2 + 1.920CH2O + 0.260GLYX + 0.215MGLY + OH'
                     
#         rxn_arrow - A STR that indciates what char sequence seperates the reactants 
#                      from products in the input rxn strings. Default is set to '=' 
#                      as it is in KPP files, but may be ->', '-->', '<->' in other 
#                      mechanism types. 
        
#         pair_delims - A list of STRS that are allowed to spereate tracer name/
#                        stoichiometry pairs from one another in rxn.This list should 
#                        NOT include the delimiter for the reaction arrow. NOTE: 
#                        Some mechanisms allow '-' instead of just '+' / why it
#                        can be a list. 
                       
#         rcts2ignore - A list of STRS to ignore in output lists if they are parsed 
#                       as tracernames in the reactants. Can be useful if you don't 
#                       want to parse 'O3+hv=O1D' as having a reactant named 'hv' 
#                       with stoich=1. 
                      
#         prds2ignore - A list of STRS to ignore in output lists if they are parsed 
#                          as tracer names in the products. Can be useful if you don't 
#                          want to actually track 'H2O' or 'O2' explicitily.           
                      
#     Outputs: 
#     -------
#         Dictionary containing input reaction, & lists  of reactants/products & stoichiometry. 
#             out ['reaction']   = 'TOLU + OH = TRO2 + 1.920CH2O + 0.260GLYX + 0.215MGLY + OH'
#                 ['rct_cmpds']  = ['TOLU','OH']
#                 ['rct_stoich'] = [1.0,1.0]
#                 ['prd_cmpds']  = ['TRO2' 'CH2O', 'GLYX', 'MGLY', 'OH']
#                 ['prd_stoich'] = [1.0, 1.920 ,0.260, 0.215, 1.0]
           
#     Author: 
#     -------
#         Dr. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jdhask
#     """
#     # Ensure the input string has no spaces: 
#     rxn_str=rxn_str.replace(' ', '') # rxn_str='TOLU+OH=TRO2+1.920CH2O+0.260GLYX+0.215MGLY+OH'
    
#     # If rxn string contains the expected delimiter: 
#     if rxn_arrow in rxn_str: 
#         # Split reaction string on the reaction arrow delimiter first: 
#         rct_half,prd_half =rxn_str.split(rxn_arrow)
#         # rct_half='TOLU+OH'   &  prd_half='=TRO2+1.920CH2O+0.260GLYX+0.215MGLY+OH'
#     else: 
#         raise ValueError('Reaction arrow delimiter "{rxn_arrow}" could not be '+\
#                          'found in the following reaction: \n\t{rxn_str}')
        
#     # Seperate the two halves of the rxn along any pairs of delimiters inputted
#     # giving us potential name/stoich pairs (e.g. Extract a list of whatever is 
#     # between '+' or '-' in 'A+B-C+D' style strings: 
        
#     #  Get list of reactant name &  stoich pairs: 
#     rct_pairs= _get_grps_between_delims(rct_half, pair_delims) # ['TOLU','OH']
    
#     #  Get list of product  name &  stoich pairs: 
#     prd_pairs= _get_grps_between_delims(prd_half, pair_delims)# [ 'TRO2', '1.920CH2O','0.260GLYX','0.215MGLY','OH'] 
                
#     # Define output dictionary: 
#     out={'reaction':rxn_str, 'rct_cmpds':[], 'rct_stoich':[], 'prd_cmpds':[], 'prd_stoich':[]} 
    
#     # Loop over each tracer/stoich string in the reactant half of the string: 
#     for pair in rct_pairs: 
#         # Seperate the stoich from tracer in this individual pair
#         stoich_i, tracer_i= _sep_stoich_n_tracer(pair) 
#         # But only add this tracer & its stoich if its NOT in the list of rcts to ignore... 
#         if tracer_i not in rcts2ignore: 
#             out['rct_cmpds'].append(tracer_i)
#             out['rct_stoich'].append(stoich_i)    
        
#     #Loop over each tracer/stoich string in the reactant half of the string: 
#     for pair in prd_pairs: 
#         # Seperate the stoich from tracer in this individual pair
#         stoich_i, tracer_i= _sep_stoich_n_tracer(pair) 
#         # But only add this tracer & its stoich if its NOT in the list of prds to ignore... 
#         if tracer_i not in prds2ignore:
#             out['prd_cmpds'].append(tracer_i)
#             out['prd_stoich'].append(stoich_i)   
    
#     return out 

# # A.2.3.1 L3-function used within parse_rxns_and_rates() within parse_kpp_main() within make_GEOSCHEM_AllRxns_file()
# def convert_fortran_nums(in_string):
#     """ Function to take an input string and convert Fortran-style numeric literals 
#         and (select) common matehmatical function names within it into modern 
#         scientific notation syntax. Used to convert numbers in KPP files to something you 
#         can write to a MATLAB/python/CSV/excel file. Used to convert strings with 
#         reaction rates like:
            
#                      In KPP File:                        as Output 
#             'GCARR_ac(1.00d-14, -490.0d0)'   to   'GCARR_ac(1.00e-14, -490.0)'

#         Can handle: 
#         - Fortran exponents like 'd', 'D', 'e', 'E' (e.g., '1.23d-4', '5.67D+8')
#         - Removes the '_dp' suffix used in Fortran for double precision #s
#         - Converts function names EXP, LOG, LOG10 to lowercase. 

#     INPUTS: 
#     ------  
#         (1) in_string - single STR containing (maybe) some numerical FORTRAN notation 

#     OUTPUTS: 
#     --------
#         (2) out_string - updated STR o with FORTRAN math converted to (modern) math.

#     USAGE: 
#     ------ 
#         Called Within: parse_rxns_and_rates() from within parse_KPP_main() from 
#                        within make_GEOSCHEM_AllRxns_file(). Can be used as standalone. 

#     REQUIREMENTS: 
#     -------------
#         LIBRARIES:           import re

#         CUSTOM FUNCTIONS:    nested/inner/helper function: _fnum_replace() 

#     AUTHOR: 
#     -------
#         Prof. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD
#     """
#     ###########################################################################
#     #  Part I: Transform Fortran Style NUMBERS into MATLAB style numbers... 
#     ###########################################################################

#     # Define a regular expression pattern to match FORTRAN-style numbers: 
#     fnum_pattern = re.compile(r'([-+]?\d*\.\d+|[-+]?\d+)([dDeE]([+-]?\d+))?(_dp)?')
    
#     #  ^^This will match numbers that may start with an optional + or - sign.
#     # The number can either be a decimal or an integer. It can have an exponent 
#     # part, which starts with either'd', 'D', 'e', or 'E', followed by an optional
#     # '+' or '-', and then digits (for example, d-4 or E+10).Finally, it can 
#     # optionally end with the suffix '_dp', which is used in Fortran to indicate 
#     # double precision numbers. After finding these, string subs are only done 
#     # when the number matched has a dif format (eg. doesn't mess with integers).
#     # Example matches &(subsequent) substitions below: 
#     #      '42'        -->  '42'       (intergers are same).
#     #      '-5.67d-8'  -->  '-5.67e-8'
#     #      '3.29E+00'  -->  '3.29'
#     #      '-5.0E-3_dp'-->  '-5.0e-3'

#     def _fnum_replace(match):
#         """
#         Nested/Helper/Local function to convert a matched Fortran number to 
#         more modern scientific notation format. 
        
#         Inputs:  match - a re.Match object (not tuple/ it should have methods 
#                                             like .group() within it). 
#         Outputs - a replacement string that has: 
#                         - exponents noted ONLY with 'e' (not 'd', 'D', or 'E"). 
#                         - no '_dp' suffix 
#                         - no exponent if the exponent was zero (e.g.'E+00','d0', 'e0')
#                         - an exponent sign ONLY if negative (e.g. no e+4).  
#         """
        
#         # Extract the relevant groups from the (individual) match:
#         number = match.group(1)        # base # (e.g., '1.23', '-4.56')
#         exponent = match.group(3)      # Sign and number pf exponent only (e.g., '-2', '+3', or None)
        
#         if exponent is None or exponent == '':
#             # Case 1: There is no explicit 'e'/'E'/'d'/'D'in the matched pattern 
#             #         Occurs for plain old numbers like '-3.14' or '+42'.
#             #         So, just return the number itself! 
#                 return number
#         else:
#             # Case 2: The number has an explicit exponent (e.g. "1.23e-4", "5.67E+8")
#             #          Convert exponent str to integer is (e.g. -4, +8, 0, etc.). 
#             exponent_value = int(exponent)
#             if exponent_value == 0:
#                 # Case 2a: If the value of exponent is 0, just remove it / only return number. 
#                 #          (e.g. '2.14d0' --> '2.14' or '42eE+00' --> '42'). 
#                 return number
#             else:
#                 # Case 2b: If the exponent is not 0/''/None & does exist, then 
#                 #          return it using the 'e' notation between the # & exponent 
#                 #          VALUE. Only retains sign of exponent sign if negative. 
#                 #          (e.g. '3.14d-4'--> '3.14e-4' but '3.14e+4'--> '3.14e4'
#                 return f"{number}e{exponent_value}"

#     # Use the re.sub() function from the regex module to: 
#     #    (1) Find all substrings in the input string that match the fortran number pattern. 
#     #    (2) For each match, call the internal _fnum_replace() function, passing 
#     #        in the match object to generate a new replacement substring (bye Fortran syntax!)
#     #    (3) Replace each matched substring with the new substring returned by _fnum_replace(). 
#     #    (4) Assign the result to a single output string after all substring replacements. 
#     out_string = fnum_pattern.sub(_fnum_replace, in_string)

#     ###########################################################################
#     #  Part II: Transform Fortran Style Math FUNCTIONS into lowercase functions: 
#     ###########################################################################
#     # Loop over some common math function names / convert them to lowercase 
#     for fname, funct_pattern in [('EXP', r'EXP\((.*?)\)'),   # Pattern to match 'EXP(' ... ')'
#                                ('LOG', r'LOG\((.*?)\)'),      # Pattern to match 'LOG(' ... ')'
#                                ('LOG10', r'LOG10\((.*?)\)')]: # Pattern to match 'LOG10(' ... ')'
#     # Use the re.sub() function from the regex module to: 
#     #    (1) Find all substrings in the input string that match the fortran function pattern. 
#     #    (2) For each match, call the lambda function, passing in the match object to 
#     #        generate a new replacement substring. It converts the function name to 
#     #        lowercase using 'fname.lower()' while maintaining the original arguments 
#     #        to the function with 'm.group(1)'.
#     #    (3) Replace each matched substring with the new substring returned by the lambda function. 
#     #    (4) Assign the result to a single output string after all substring replacements. 
#     # This results in replacing, for example, 'EXP(1.23)' with 'exp(1.23)', or 'LOG10(4.56)' with 'log10(4.56)'

#         out_string = re.sub(funct_pattern,lambda m: f"{fname.lower()}({m.group(1)})",out_string)
    
#     return out_string

# # A.2.3 L2-function used in parse_kpp_main() called within make_GEOSCHEM_AllRxns_file():
# def _is_blank(line):
#     """Function to determine if a line is blank or not after all spaces have been removed."""
#     tf = True if len(line.replace(' ', '')) == 0 else False
#     return tf

# # A.2.3 L2-function used in parse_kpp_main() called within make_GEOSCHEM_AllRxns_file():
# def parse_rxns_and_rates(line: str, kpp_dict: dict, file_location: dict):
#     """Function used to parse reaction lines in KPP file and extract their info into dictionary, kpp_dict. 

#     INPUTS:
#     -------

#         (1) line - STR of line with reaction, rate, and (maybe) a comment. These lines look like this:     
#                       Reaction               Rate                              Comment
#                 "O3 + NO = NO2 + O2 :  GCARR_ac(3.00d-12, -1500.0d0); {2014/02/03; Eastham2014; SDE}"

#         (2)  kpp_dict - DICT with the following Key/Value Pairs.
#                 KEYS           :           VALUES
#                 ----------------------------------------------------------------------
#                 KPP_line       : List of original line from KPP file
#                 rxn            : List of strings with full reaction only... 
#                 reactants      : *Nested* list containing all reactants per reaction.  
#                 rct_stoich     : *Nested* list of reactants stoichiometry per reaction.  
#                 products       : *Nested* list containing all products per reaction. 
#                 prd_stoich     : *Nested* list of products stoichiometry per reaction. 
#                 is_gas         : Boolean List of whether rxn is gas phase or not. 
#                 is_het         : Boolean List of whether rxn is heterogeneous or not. 
#                 is_photo       : Boolean list of whether rxn is a photolysis rxn or not. 
#                 is_3body       : Boolean list of whether rxn includes 3rd body reactant, M, or not. 
#                 rate           : List of strings with full rate for each reaction 
#                 rate_function  : List of strings with rate function for each reaction
#                 comments       : List of comments for each reaction 

#                 NOTE: All lists in the dict, kpp_dict, should be the same length upon input and output to 
#                     this function, as they all are all "indexed "by the reaction they refer to. 

#         (3) file_location - Dictionary with the following Key/Value Pairs: 

#                 KEYS           :           VALUES
#                 ----------------------------------------------------------------------
#                 in_header      : Boolean. True if currently parsing the header, otherwise False.   
#                 in_defvar      : Boolean. True if currently parsing the fortran variable defintion section
#                 in_deffix      : Boolean. True if currently parsing the "#deffix" section  
#                 in_sulf_rxns   : Boolean. True if currently parsing sulfate mod het reactions  
#                 in_gas_rxns    : Boolean. True if currently parsing Gas Phase reactions
#                 in_het_rxns    : Boolean. True if currently parsing Heterogeneous reactions
#                 in_photo_rxns  : Boolean. True if currently parsing Photolysis reactions

#      OUTPUTS: 
#      -------   
#          (1)  kpp_dict       - Updated DICT with the same Key/Value Pairs as above, but outputted 
#                              lists in kpp_dict will be +1 in len since they will now include info 
#                              about the reaction contained on the input line. This function checks 
#                              for len consistency @ end. 

#         (2) file_location - Updated DICT with the same Key/Value Pairs as above. 

#     REQUIREMENTS: 
#     -------------
#         LIBRARIES:           None

#         CUSTOM FUNCTIONS:    
#                             _str_multi_replace()
#                             sep_stoich_vs_tracer()
#                             fortran2matlab_nums()
                            

#     USAGE: 
#     -----  
#         Called Within: parse_kpp_main() within make_GEOSCHEM_AllRxns_file(). 
#         Output Usage: Output dict 'kpp_dict' used in parse_kpp_main() to sep info for dif rxn types into own dicts! 

#     AUTHOR: 
#     -------
#         Prof. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD

#     """
#     # Clean up line to maintain in output: 
#     clean_line= utils.str_multi_replace(line, ['\n', '\t'], rep_all='').strip()
    
#     #####################################################################################
#     # Split Line in to reaction, rate, rate function, and comment & Clean up for MATLAB!
#     #####################################################################################
#     rxn_parts = line.split(':')  # Split on colon to get rxn only
#     rate = rxn_parts[1].split(';')[0].strip() # Split on ';' to get rates, remove trailing/tailing spaces
        
#     # Define Booleans that tell us what type of rxn it is: 
#     rxn_with_M=True if '{+M}' in rxn_parts[0] else False  
#     rxn_is_gas= file_location['in_gas_rxns']
#     rxn_is_photo= file_location['in_photo_rxns']
#     rxn_is_het=True if any ([file_location['in_het_rxns'],file_location['in_sulf_rxns']]) else False 

#     # If this part of the line after the reaction has a comment, then split that on'{' and clean it up.
#     if '{' in rxn_parts[1]:
#         replace_me = dict({'{': '', '}': '', ';': ',', '\n': '','\t':''})
#         comment = utils.str_multi_replace(rxn_parts[1].split('{')[1], replace_me)
#     else:
#         comment = ''  # Otherwise make comment blank!

#     # Update comment to reflect whether or not the rxn was a 3 body rxn with M or not.
#     if rxn_with_M == True: comment= '(3-Body rxn with +M) '+comment
    
#     # Now remove spaces, and '{+M'}, '{' , '}' characters from the reaction
#     rxn = utils.str_multi_replace(rxn_parts[0], ['{+M}', '{', '}', ' '], rep_all='').strip()

#     # Split entire reaction into Reactants & Products & their stoichiometry:
#     stoich_dict= parse_reaction_str(rxn, rxn_arrow='=', pair_delims=['+','-'],
#                            rcts2ignore=['hv'], prds2ignore=[])
    
#     # Clean up rates so the #s no longer have FORTRAN style numbers (does NOT
#     # change operators), only #s and case ofcommon math functions: 
#     rate=convert_fortran_nums(rate)

#     # Get function name only from the rate constant defs/ The lines look like this: 
#     # 'GCARR_ac(2.60e-13, 1300)' but some maybe sums of calls to multiple functions 
#     # like so: GCARR_ac(2.60e-13, 1300) + GCARR_ac(2.60e-13, 1300) and some may 
#     # contain exp(), log(), or even log10()... So use cray regex: 
#     pattern = r'\b((?!exp\b|log\b|log10\b)[A-Za-z_]+[a-zA-Z0-9_]*)\('

#     # Use a list comprehension to apply the regex to each string in rates
#     functions = [match_i[:] for match_i in re.findall(pattern, rate)]     
        
#     # ------------------------------------------------------------------------- 
#     # Now save everything into dict where rxn is the key: 
#     # ------------------------------------------------------------------------- 
#     # Update with new index on each call. 
#     inds=list(kpp_dict.keys())
#     new_ind=0 if len(inds) ==0 else inds[-1]+1
    
#     kpp_dict[new_ind]= dict({# List of original KPP file lines.
#                          'KPP_line': clean_line,   
#                          # STR with parsed reaction: 
#                          'reaction': rxn,   
#                          # List of compounds in reactants: 
#                          'rct_cmpds': stoich_dict['rct_cmpds'],
#                          # Stoichiometry of compounds in reactants (same order!!)  
#                          'rct_stoich': stoich_dict['rct_stoich'],
#                          # List of compounds in products: 
#                          'prd_cmpds': stoich_dict['prd_cmpds'],
#                          # Stoichiometry of compounds in products (same order!!)  
#                          'prd_stoich': stoich_dict['prd_stoich'],
#                          #STR with rate constant for reaction 
#                          'rate': rate,
#                          # STR with rate function 
#                          'rate_function': functions,
#                          # STR comment (if rxn has one) 
#                          'comments': comment,    
#                          # BOOL  of whether rxn is gas phase or not. 
#                          'is_gas': rxn_is_gas,
#                          # Boolean of whether rxn is heterogeneous or not.
#                          'is_het': rxn_is_het,
#                          # BOOL  of whether rxn is a photolysis rxn or not.
#                          'is_photo': rxn_is_photo,
#                          # BOOL  of whether rxn is a 3-bdy rxn w/ M or not.
#                          'is_3body': rxn_with_M,
#                          })

#     return kpp_dict

# # A.2.2 L2-function used in parse_kpp_main() called within make_GEOSCHEM_AllRxns_file():
# def fix_multiline_rxns(line, line_stripped, count, kppfile):
#     """ Function used to parse reactions that are listed on multiple lines in 
#     the KPP file to connect them into a single line with the entire reaction/rate on one line. 
#     These lines that should get sent to this function will look like this in the KPP file: 

#         ln#     'Some text in KPP file line'
#         41      'MVKOHOO + HO2 = 0.360MCO3 + 0.360GLYC +' 
#         42      '  0.335MVKHP + 0.050MGLY + 0.050CH2O :  GCARR(2.12E-13, 0.0E+00, 1300.0); {2019/11/06; Bates2019; KHB}' 

#     After this function completes, the outputted "line" should look like this so the rreaction parser will work! 
#          'O3 + MO2 = CH2O + HO2 + O2 :   GCARR(2.90E-16, 0.0E+00, -1000.0);   {2014/02/03; Eastham2014; SDE}' 
#                                      ^split #1                            ^ split #2  

#     INPUTS: 
#     -------

#         (1)  line - String of the current KPP file line you're parsing for info. This line 
#                will contain only half of the reaction and will end with a '+' because 
#                that is what triggers this function being called! 

#         (2)  line_stripped- String that is the same as "line", but has no spaces in it. 

#         (3) count- Int containing the current line number in the file. Is the index of the 
#                line that contains the first half of the reaction's stuff. 

#         (4) kpp_file- String containing full path to the KPP input file that's being parsed. 
#             Needed here so we can read the NEXT lines with info about this continued reaction. 


#     OUTPUTS: 
#     --------

#         (1) line- Updated string that has all the info about this reaction in a single string. 
#               It is crafted by concatonating continuing lines together, but does not contain any new line chars. 

#         (2) line_stripped- Updated string containing the same info as updated "line", but has all spaces removed. 

#         (3) count- Updated int containing line number of the last line you read in. Is the index of the line 
#                 in the KPP file where the split reaction ENDS! Updated so the next iteration doesn't 
#                 try to parse a line starting in the middle of this split reaction... 

#     REQUIREMENTS: 
#     -------------
#         LIBRARIES:           None

#         CUSTOM FUNCTIONS:    None

#     USAGE: 
#     ------  
#         CALLED WITHIN:  parse_kpp_main() via a call to make_GEOSCHEM_AllRxns_file() 

#         OUTPUT USAGE: Updated output "line" next passed to parse_rxns_and_rates() within 
#                         parse_kpp_main() called within make_GEOSCHEM_AllRxns_file(). 

#     AUTHOR: 
#     -------
#         Prof. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD

#     CHANGE LOG: 
#     ----------
#         09/28/2023    JDH Created 
#         11/15/2024: JDH added "':' not in line_stripped: " to while statement 
#                     to work with Alfie's custom version of 14.5.0... 

#     """

#     # Read more lines til you get the entire reaction, rate & comments in a single line!
#     while line_stripped.endswith('+') == True or ':' not in line_stripped:
#         nline = kppfile.readline()   # Read next line from file
#         count = count+1                 # Update line counter variable
#         # append this new line to our current line str.
#         line = line + nline.strip()
#         # Update line_stripped to not contain spaces!
#         line_stripped = line.strip().replace(' ', '')

#     # Now that you got the entire reaction on this line, make sure to replace any new line characters with a space!
#     line = line.replace('\n', ' ')
#     line_stripped = line_stripped.replace('\n', ' ')

#     # Return the count, new line, and line_stripped vars to continue on.
#     return line, line_stripped, count

# # A.2.1 L2-function used in parse_kpp_main() called within make_GEOSCHEM_AllRxns_file():
# def parse_var_defs(line: str, tracers: list = [], long_name: list = [], split_on: str = 'ignore'):
#     """Function to parse variable defintion lines in the KPP File header. These lines look like: 
#         'A3O2       = IGNORE; {CH3CH2CH2OO; Primary RO2 from C3H8}'  

#     INPUTS: 
#     -------
#         (1) line      - STR of the current KPP file line you're parsing for info. 

#         (2) tracers   - LIST containing names of all tracers that have already been defined 
#                          by parsing previous lines. Only contains tracer abbreviation! 

#         (3) long_name - LIST containing extra comment info about the tracer being defined that 
#                          has been compiled by parsing previous lines. 

#         (4) split_on  - (OPTIONAL) STR we want to split the input line on... 
#                         It is case insensitive & space insensitive. Default is set to 
#                         'ignore' to work with most recent GC KPP files. 

#     OUTPUTS:
#     -------- 
#         tracers   - Updated LIST with this line's new defined variable in it.

#         long_name - Updated LIST with this line's tracer info in it. 

#     REQUIREMENTS: 
#      -------------
#         LIBRARIES:           None

#         CUSTOM FUNCTIONS:    _str_multi_replace()

#     USAGE: 
#     -----  
#         CALLED WITHIN:  parse_kpp_main() via a call within make_GEOSCHEM_AllRxns_file(). 

#         OUTPUT USAGE:   Zipped into dict "var_info" passed as output of make_GEOSCHEM_AllRxns_file()

#     AUTHOR: 
#     -------
#         Prof. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD

#     CHANGE LOG: 
#     ----------
#         09/28/2023    JDH Created 
#         10/30/2023    JDH updated documentation. 

#     """

#     if split_on in line.lower():
#         # Figure out index of where to split line between def and comments. Case insensitive.
#         ind = line.lower().index(split_on)

#         # Remove spaces and "=" from line, to get only the tracer name being defined!
#         tracer = utils.str_multi_replace(line[:ind],  [' ', '='], rep_all='')

#         # Keep tracer name of this compound defined in a list.
#         tracers.append(tracer)

#         # Pull out the rest of the comment about this tracer using the index.
#         tracer_info = line[ind+len(split_on):]

#         # Remove any remaining weird characters we don't care about (";", "{", "}").
#         if tracer_info[0] == ';':
#             tracer_info = tracer_info[1:]
#         tracer_info = utils.str_multi_replace(
#             tracer_info,  ['{', '}'], rep_all='').strip()

#         # Now add the "long name" slash comment about this tracer to a list!
#         long_name.append(tracer_info)

#     return tracers, long_name

# # A.2 L1-function called within make_GEOSCHEM_AllRxns_file()
# def parse_kpp_main(kppfile, out2yaml:bool=True, out2excel:bool=False, debug:bool=False):
#     """Function to parse a KPP file and extract info about tracers, reactions,
#     & rates in a pandas dataframe. Output can also be provided in other ways using 
#     keywords. 

#     INPUTS:  
#     -------  
#        (1) kppfile - STR containing the full path and name of a KPP file.  

#        (2) out2yaml - BOOL indicating if outputs should be saved to a yaml file. 
       
#        (3) out2excel - BOOL indicating if outputs shoudl be saved to an excel .xlsx sheet. 
              
#        (4) debug- BOOL for printing additional output for debugging. Default is FALSE. 
       
#     OUTPUTS:  
#     -------- 
#        (1) var_info - Dictionary with all tracer names / long names read in from KPP. 
       
#        (2) kpp_dict - Dictonary with reactions, rates, stoichiometry, comments, 
#                       & reaction type for each reaction in KPP file plus more.
#                       Is saved in yaml file if output2yaml is set to True. 
       
#        (3) kpp_df - Same info in kpp_dict, but as a dataframe for easy viewing.
#                     Is saved in .xlsx file if output2excel is set to True.
       

#     REQUIREMENTS: 
#     -------------
#         LIBRARIES:           None

#         CUSTOM FUNCTIONS:   parse_var_defs()
#                             _is_blank()
#                             fix_multiline_rxns()
#                             parse_rxns_and_rates()
#                             utils.check_filename()
                            

#     USAGE: 
#     -----  
#         CALLED WITHIN: make_GEOSCHEM_AllRxns_file()

#         OUTPUT USAGE: Info in kpp_dict translated to F0AM... 

#     AUTHOR: 
#     -------
#        Prof. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD
#     """
#     ###########################################################################
#     # Initialize all vars to hold outputs/ keep track of where we are!
#     ###########################################################################

#     # Get # of lines to expect so we know when we've reached the end of the file!
#     num_lines = sum(1 for line in open(kppfile))
#     tracers = []
#     long_name = []  # Make empty lists to hold tracers & tracer info!
#     count = 0  # Initialize line counter variable.

#     # Define boolean dictionary with key/values that will tell us what part of the
#     # KPP file we are currently parsing.
#     file_location = dict({
#         'in_header': True, # Boolean. True on initiation & if currently parsing header.
#         'in_defvar': False,# Boolean. True if currently parsing the Fortran variable definition section
#         'in_deffix': False,  # Boolean. True if currently parsing the "#deffix" section
#         'in_sulf_rxns': False,  # Boolean. True if currently parsing sulfate mod het reactions
#         'in_gas_rxns': False,  # Boolean. True if currently parsing gas phase reactions
#         'in_het_rxns': False,  # Boolean. True if currently parsing heterogenous reactions
#         'in_photo_rxns': False  # Boolean. True if currently parsing photolysis reactions
#     })


#     ##############################################################################
#     # Open the KPP file in "read" mode & Loop over all lines until end is reached!
#     ##############################################################################
#     file = open(kppfile, 'r')
#     kpp_dict={} 
#     while True:  # loop over each line of the KPP file.
    
#         line = file.readline()  # Read next line from file
        
#         # Strip all white space from current line.
#         line_stripped = line.strip().replace(' ', '')
        
#         # Update line counter variable now that you read in this line.
#         count = count+1

#         # --------------------------------------------------------------------------------------------------------------------------
#         # BEFORE parsing the line, set file location booleans to FALSE if the current line indicates you're entering a new section!
#         # --------------------------------------------------------------------------------------------------------------------------
#         # You are no longer in the header section once #defvar appears on the line!
#         file_location['in_header'] = False if '#defvar' in line_stripped.lower() else file_location['in_header']
#         # You are no longer in the variable defintion section if #deffix appears in the current line!
#         file_location['in_defvar'] = False if '#deffix' in line_stripped.lower() else file_location['in_defvar']
#         # You are no longer in section defining fixed vars if #equations appears in the current line!
#         file_location['in_deffix'] = False if '#equations' in line_stripped.lower() else file_location['in_deffix']
#         # You are no longer in the sulfate_mod reaction section on the next line if 
#         # all of these strs: '//', 'gas', 'phase', 'chemistry', & 'reactions' appear in the current line.
#         file_location['in_sulf_rxns'] = False if all([item_i in line_stripped.lower() \
#                                         for item_i in ['//', 'gas', 'phase', 'chemistry', 'reactions']]) \
#                                         else file_location['in_sulf_rxns']
#         # You are no longer in the gas phase reaction section on the next line if 
#         # all of these strs: '//', 'heterogeneous','chemistry', & 'reactions' appear in the current line!
#         file_location['in_gas_rxns'] = False if all([item_i in line_stripped.lower() \
#                                        for item_i in ['//', 'heterogeneous', 'chemistry', 'reactions']]) \
#                                        else file_location['in_gas_rxns']
#         # You are no longer in the photolysis reaction section n the next line if 
#         # all of these strs: '//', 'photolysis', & 'reactions' appear in the current line!
#         file_location['in_het_rxns'] = False if all([item_i in line_stripped.lower() \
#                                        for item_i in ['//', 'photolysis', 'reactions']]) \
#                                        else file_location['in_het_rxns']
        
#         # Get a list of the keys that are currently set to "TRUE" (useful for debugging).. 
#         true_keys = '.'.join([key for key, value in file_location.items() if value])
#         if debug == True: 
#             print(f'Line: "{line}"\n\tFile Location: {true_keys} \n')

#         # Break if reached end of file, don't attempt to parse line.
#         if count > num_lines:
#             break
#         # ------------------------------------------------------------------------------------------
#         # If you're in the part of the file where variable definitions are made  //
#         # tracer names are listed then pass to ==> parse_var_defs()
#         # ------------------------------------------------------------------------------------------
#         if file_location['in_defvar'] == True or file_location['in_deffix'] == True:

#             # Pass to parse_var_defs() & return updated lists with tracer name and long name.
#             tracers, long_name = parse_var_defs(line, tracers=tracers, long_name=long_name, split_on='ignore')

#         # -----------------------------------------------------------------------------------------
#         # If you're in part of the file with chemical equations & rates ==> parse_rxns_and_rates()
#         # -----------------------------------------------------------------------------------------
#         if any(cond == True for cond in [file_location['in_sulf_rxns'], 
#                                          file_location['in_gas_rxns'], 
#                                          file_location['in_photo_rxns'],
#                                          file_location['in_het_rxns']]):

#             # Don't parse commented out or blank lines.
#             if (('//' != line[0:2]) and (_is_blank(line) == False)):

#                 # Decide if you've been given a line that contains a partial reaction 
#                 # & fix it before parsing as rxn.
#                 if ':' not in line:
        
#                     # These lines won't contain a ':' and look like this:
#                     #   ln #     Content of Line
#                     # -------------------------------------------------------------
#                     #   ln 42   'MVKOHOO + HO2 = 0.360MCO3 + 0.360GLYC +'   **(Would trigger when parsing this line)**
#                     #   ln 43   '  0.335MVKHP + 0.050MGLY + 0.050CH2O :  GCARR(2.12E-13, 0.0E+00, 1300.0); {2019/11/06; Bates2019; KHB}'

#                     line, line_stripped, count = fix_multiline_rxns(line, line_stripped, count, file)

#                 # (After fixing if necessary), parse this reaction and add all info about this rxn to lists:
#                 kpp_dict = parse_rxns_and_rates(line, kpp_dict, file_location)
                

#         # --------------------------------------------------------------------------------------------------------------------------
#         # AFTER parsing the line, set file location booleans to TRUE if the current line indicates you're entering a new section!
#         # --------------------------------------------------------------------------------------------------------------------------
#         # You will be in the variable definition section on the next line if '#defvar' appears in the current line.
#         file_location['in_defvar'] = True if '#defvar' in line_stripped.lower() else file_location['in_defvar']
#         # You will be in the deffix section on the next line if #deffix appears in the current line.
#         file_location['in_deffix'] = True if '#deffix' in line_stripped.lower() else file_location['in_deffix']
#         # You will be in the sulfate_mod reaction section on the next line if
#         # all these strings appear in the current line.
#         file_location['in_sulf_rxns'] = True if all([item_i in line_stripped.lower() \
#                                         for item_i in ['//', 'sulfate_mod', 'reactions']]) \
#                                         else file_location['in_sulf_rxns']
#         # You will be in the gas-phase reaction section on the next line if 
#         # all these strings appear in the current line.
#         file_location['in_gas_rxns'] = True if all([item_i in line_stripped.lower() \
#                                             for item_i in ['//', 'gas', 'phase', 'chemistry', 'reactions']]) \
#                                             else file_location['in_gas_rxns']
#         # You will be in the heterogenous reaction section on the next line if 
#         #all these strings appear in the current line.
#         file_location['in_het_rxns'] = True if all([item_i in line_stripped.lower() \
#                                             for item_i in ['//', 'heterogeneous', 'chemistry', 'reactions']]) \
#                                             else file_location['in_het_rxns']
#         # You will be in the photolysis reaction section on the next line if 
#         # all these strings appear in the current line.
#         file_location['in_photo_rxns'] = True if all([item_i in line_stripped.lower(
#         ) for item_i in ['//', 'photolysis', 'reactions']]) else file_location['in_photo_rxns']

#     # While loop will break once we found the end of the KPP file and have parsed all lines.
#     file.close()  # Close the KPP file.

#     # Build var_info, a dictionary connecting tracer names & long names.
#     var_info = dict(zip(tracers, long_name))
    
#     # --------------------------------------------------------------------------
#     # Seperate out reactions from kpp_dict into dif dicts for each type of rxn.
#     # --------------------------------------------------------------------------
#     # First, just make sure kpp_dict has right keys & consistent list lengths!
#     #_check_kpp_dict(kpp_dict, check_lens=True, ignore_extras=False)

#     # Make dataframe with all this information... 
#     kpp_df = pd.DataFrame(kpp_dict)
    
#     if out2yaml==True: 
#         yaml_file=utils.check_filename(filename='kpp_mech2', default_name= 'kpp_mech2', ext='.yml', 
#                            savepath=master_outdir, overwrite=True, return_full=True, quiet=True)
#         # Write the dictionary to a YAML file
#         with open(yaml_file, 'w') as f:   
#             yaml.dump(kpp_dict, f)
#         print(f'Output of () saved as dictionary at: \n\t{yaml_file}')
    
#     if out2excel==True: 
#         excel_file=utils.check_filename(filename='kpp_mech2', default_name= 'kpp_mech2', ext='.xlsx', 
#                            savepath=master_outdir, overwrite=True, return_full=True, quiet=True)
#         # Write the dataframe to an excel file: 
#         kpp_df.to_excel(excel_file)
#         print(f'Output of parse_kpp() saved as dataframe at: \n\t{excel_file}')
            
#     return var_info, kpp_df, kpp_dict 




#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  8 00:36:37 2025

@author: u6044586
"""
import sys
import re 
import numpy as np 

sys.path.append('/uufs/chpc.utah.edu/common/home/u6044586/python_scripts/modules/GEOSChem_Emulator/')
import utils_class as utils

###############################################################################
#                  Functions for comparing reactions/ mechanisms: 
################################################################################
def not_cmt(rxn_obj): 
    tf= False if rxn_obj.comment_out is True else True 
    return tf

def update_attr(rxn_obj0, attr, new_value, verbose:bool=True):
    """Update a specified attribute of the Reaction object with a new value."""
    
    # Make a copy, don't mess with the original ... 
    rxn_obj=rxn_obj0.copy() 
    if hasattr(rxn_obj, attr):
        setattr(rxn_obj, attr, new_value)
        #if verbose: 
           # print(f"{attr} updated to {new_value}.")
    #else:
       # if verbose:
            #print(f"Attribute {attr} does not exist in the kpp_reaction object.")
            
    return rxn_obj

def calc_rxn_similarity(rxn1, rxn2): 
    ###########################################################################
    #                            Check/Calc reactant similarity: 
    ###########################################################################
    # Check if the 2 rxns have the same reactant compounds: 
    same_rcts =True if set(rxn1.unq_rct_cmpds)==set(rxn2.unq_rct_cmpds) else False
    
    # Calculate the degree of similaritiy of the reactant compounds' stoichiometry: 
    comp_rct_stoich=[]
    for cmpd in list(set(rxn1.unq_rct_cmpds+rxn2.unq_rct_cmpds)):
        
        # If this compound is a reactant in both rxns AND it's stoich is same in both: 
        if ((cmpd in list(rxn1.unq_rct_stoich.keys())) and 
            (cmpd in list(rxn2.unq_rct_stoich.keys()))and 
            (rxn1.unq_rct_stoich[cmpd]==rxn2.unq_rct_stoich[cmpd])): 
                # then append a 1 to the list comparing their stoichiometry: 
                comp_rct_stoich.append(1) 
        else: 
            # Otherwise append a 0 (for either dif stoich or only in 1 of the rxns). 
            comp_rct_stoich.append(0) 
    
    # Calculate similarity score for reactants as percentage of the total unique
    # reactants occuring in either reaction where stoich is same in other reaction. 
    #      Note: Cannot be 100% if the reactions don't have identical reactants. 
    rct_similarity=(sum(comp_rct_stoich)/len(comp_rct_stoich))*100 
    
    ###########################################################################
    #                            Repeat for Products: 
    ###########################################################################
    same_prds =True if set(rxn1.unq_prd_cmpds)==set(rxn2.unq_prd_cmpds) else False

    # Calculate the degree of similaritiy of the product compounds' stoichiometry: 
    comp_prd_stoich=[]
    for cmpd in list(set(rxn1.unq_prd_cmpds+rxn2.unq_prd_cmpds)):
        
        # If this compound is a reactant in both rxns AND it's stoich is same in both: 
        if ((cmpd in list(rxn1.unq_prd_stoich.keys())) and 
            (cmpd in list(rxn2.unq_prd_stoich.keys()))and 
            (rxn1.unq_prd_stoich[cmpd]==rxn2.unq_prd_stoich[cmpd]) ): 
            # then append a 1 to the list comparing their stoichiometry: 
                comp_prd_stoich.append(1) 
        else: 
            # Otherwise append a 0 (for either dif stoich or only in 1 of the rxns). 
            comp_prd_stoich.append(0) 
    
    # Calculate similartiy score for products as percentage of the total unique # of
    # products occuring in either reaction where stoich is same in other reaction. 
    #      Note: Cannot be 100% if the reactions don't have identical products. 
    prd_similarity=(sum(comp_prd_stoich)/len(comp_prd_stoich))*100     
    
    return same_rcts, rct_similarity, same_prds, prd_similarity


def are_reactions_equal(rxn1, rxn2):
    """Checks if two reactions are equal based on reactants and products with their stoichiometry."""
    
    same_rcts, rct_similarity, same_prds, prd_similarity=calc_rxn_similarity(rxn1, rxn2)
    
    if ((same_rcts==True) and (rct_similarity==100) and 
        (same_prds==True) and (prd_similarity==100) and (rxn1.rate==rxn2.rate)): 
        return True
    else:
        return False
        
        
###############################################################################
#               Functions for defining "KPP Reaction"  class: 
###############################################################################
def sep_stoich_n_tracer(pair:str): 
    """Helper function to take a single input string & seperate it based on where 
    the first alpha numeirc character appears. Anything before first letter 
    is assumed to be stoichimetry. Anthing after assumed to be tracer name. 
    If first char is letter, stoichiometery assumed to be 1.""" 
    
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
        
def get_grps_between_delims(s:str):
    """Helper function used to seperate one half of a rxn along any pairs 
    of delimiters inputted, giving us a list of individual tracer name/stoich 
    pairs. Tracer name/stoichiomety seperated later... (e.g. This just extracts 
    a list of whatever is between '+' or '-' in 'A+B-C+D' style strings). 
    """
    delimiters=['+','-']
    
    # Join the delimiters list into a regex pattern string
    delim_str = '|'.join(map(re.escape, delimiters))  # Escape delimiters for regex
    
    # Using regex to split but capture the delimiters
    pattern = f'(?<=[{re.escape(delim_str)}])|(?=[{re.escape(delim_str)}])'
    parts = re.split(pattern, s)
    
    # Assume delimter is '+' if no delim preceeds the group (e.g. assume all 
    # stoichiometry is positive unless it is explicityly noted as being negative). 
    current_delimiter = '+'
    
    result=[]
    for part in parts:
        if part in delimiters:
            current_delimiter = part
        elif part:
            result.append(current_delimiter+part)
    
    return result

def parse_indv_reaction_str(rxn_str:str, rxn_arrow:str='=',ignore:list=['M','hv'], 
                             is_photo:bool=False):
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
                       
        ignore - A list of STRS to ignore in output lists if they are parsed 
                      as tracernames in the reactants/products. Can be useful if 
                      you don't want to parse 'O3+hv=O1D' as having a reactant 
                      named 'hv' with stoich=1. 
                      
                      
    Outputs: 
    -------
        Dictionary containing input reaction, & lists of reactants/products &
        stoichiometry.NOTE: For duplicate entries of products/reactants, the 
        default behavior is to accumulate the stoichiometry b/c keys must be unique. 
        So, 'OH + OH = H2O2' results in  out['reactants']= {'OH':2}
        
        Example: 
            out['reaction']   = 'TOLU + OH = TRO2 + 1.920CH2O + 0.260GLYX + 0.215MGLY + OH'
            
                ['reactants']  = {'TOLU': 1.0,
                                    'OH': 1.0}
                                  
                ['products']  = {'TRO2' : 1.0,
                                  'CH2O': 1.920, 
                                  'GLYX': 0.260, 
                                  'MGLY': 0.215, 
                                  'OH'  : 1.0}
   

    REQUIREMENTS: 
    -------------
        LIBRARIES:          None 

        CUSTOM FUNCTIONS:    
                           get_grps_between_delims() 
                           sep_stoich_n_tracer

    USAGE: 
    -----  
        Called Within: parse_rxns_and_rates()
                      
        Output Usage: Contents used to fill the contents of reaction class for 
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
    rct_pairs= get_grps_between_delims(rct_half) # ['TOLU','OH']
    
    #  Get list of product  name &  stoich pairs & if positive/negative! 
    prd_pairs= get_grps_between_delims(prd_half)# [ 'TRO2', '1.920CH2O','0.260GLYX','0.215MGLY','OH'] 
                

    # Define output dictionary: 
    out={'rct_cmpds':{},'rct_stoich':{},'prd_cmpds':{},'prd_stoich':{}}
        
    # Loop over each tracer/stoich string in the reactant half of the string: 
    r_ct=0                  
    for pair in rct_pairs: 
        # Seperate the stoich from tracer in this individual pair
        stoich_i, tracer_i= sep_stoich_n_tracer(pair) 
        
        # But only add this tracer & its stoich if its NOT in the list of cmpds to ignore... 
        if tracer_i not in ignore: 
            out['rct_cmpds'][r_ct]= tracer_i
            out['rct_stoich'][r_ct]= stoich_i
            r_ct=r_ct+1
        
    # Loop over each tracer/stoich string in the reactant half of the string: 
    p_ct=0
    for pair in prd_pairs: 
        # Seperate the stoich from tracer in this individual product stoich/tracer pair
        stoich_i, tracer_i= sep_stoich_n_tracer(pair) 
                
        # Only add this tracer & its stoich if its NOT in the list of cmpds to ignore... 
        if tracer_i not in ignore: 
            out['prd_cmpds'][p_ct]= tracer_i
            out['prd_stoich'][p_ct]= stoich_i
            p_ct=p_ct+1
                    
    return out
        
class kpp_reaction(object):
    """Generator for a single reaction. Holds a reaction's reactants, products, 
    and rate priding a way to retrieve these items easily. 

    Attributes
    ----------
    'reaction'       : STR with reaction 
    'rxn_cmts'       : STR with comments about this reaction (or blank). 
    'rate'           : STR with full rate (func & arg) for this reaction 
    'rate_cmts'       : STR with comments about this  (or blank). 
    'is_gas'         : BOOL of whether rxn is gas phase or not. 
    'is_het'         : BOOL of whether rxn is heterogeneous or not. 
    'is_photo'       : BOOL of whether rxn is a photolysis rxn or not. 
    'is_3body'       : BOOL of whether rxn increactionludes 3rd body reactant, M, or not. 
    'rct_cmpds'      : LIST of all reactants in this reaction.  
    'rct_stoich'     : Dict where keys are th tracer names of products (in prd_cmpds) and values 
                            are the stoichiomentry of that specific product. 
    'prd_cmpds'      : LIST of all products in this reaction. 
    'prd_stoich'     : Dict where keys are th tracer names of products (in prd_cmpds) and values 
                            are the stoichiomentry of that specific product. 
    'original_line'       : STR with original line from KPP file
    'rate_function'  : STR with rate function (only- no args or constants) used.
    'comment_out'        :    False 
                
    """
    def __init__(self, line:str, is_gas:bool, is_het:bool, is_photo:bool,
                 rxn_arrow:str='=',ignore:list=['hv','M'],
                 fixed_vars=None):
        """Initialize info about a rxn read from a KPP line for a KPP_reaction object
        
        INPUTS:
        -------
    
            (1) line - STR of line with reaction, rate, and (maybe) a comment. These lines look like this:   
                
                          Reaction               Rate                              Comment
                    "O3 + NO = NO2 + O2 :  GCARR_ac(3.00d-12, -1500.0d0); {2014/02/03; Eastham2014; SDE}"
                

        REQUIREMENTS: 
        -------------
            LIBRARIES:           import re 
    
            CUSTOM FUNCTIONS:    
                                utils._str_multi_replace()
                                utils._convert_fortran_numbers()
                                utils._convert_fortran_operators
                                parse_indv_reaction_str()
    
        AUTHOR: 
        -------
            Prof. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD
    
        """
        # Clean up line to maintain in output: 
        clean_line= utils.str_multi_replace(line, ['\n', '\t'], rep_all='').strip()
        
        #####################################################################################
        # Split Line in to reaction, rate, rate function, and comment & Clean up!
        #####################################################################################
        comment=''
        
        # Split on colon to get rxn only
        rxn_parts = line.split(':')  
        
        # Split on ';' to get rates, remove trailing/tailing spaces
        rate = rxn_parts[1].split(';')[0].strip() 
            
        if '{+M}' in rxn_parts[0]:
            # Define Boolean that tell us if its a 3 body rxn or not.
            rxn_with_M=True 
            # Remove the '{' and '}' from around +M in 3 body rxn strings: 
            rxn_parts[0]=rxn_parts[0].replace('{+M}', '+M')
        else: 
            rxn_with_M=False; 
            
        # If this part of the line after the reaction has a comment, then split that 
        # on'{' and clean it up.
        if '{' in rxn_parts[1]:
            replace_me = dict({'{': '', '}': '', ';': ',', '\n': '','\t':''})
            comment = utils.str_multi_replace(rxn_parts[1].split('{')[1], replace_me)
    
        # Now remove spaces, and '{' , '}' characters from the reaction
        rxn = utils.str_multi_replace(rxn_parts[0], ['{+M}', '{', '}', ' '], rep_all='').strip()
     
        # Split entire reaction into Reactants & Products & their stoichiometry:
        sdict= parse_indv_reaction_str(rxn, rxn_arrow='=', ignore=['hv', 'M'], is_photo=is_photo)
        
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
        # Now save everything as an attribute of this rxn object: 
        # -------------------------------------------------------------------------
        
        # Add & populate self attributes with info about this specific reaction to kpp_dict. 
        self.reaction=rxn   # STR with reaction 
        self.rxn_cmts=comment # STR with comments about this rxn (or blank). 
        self.rate=rate # STR with rate (func & arg) for this rxn
        self.rate_cmts='' # STR with comments about this rate (or blank). 
        self.is_3body=rxn_with_M # BOOL of whether rxn includes 3rd body reactant, M, or not. 
        self.is_gas=is_gas # BOOL of whether rxn is gas phase or not. 
        self.is_het=is_het # BOOL of whether rxn is heterogeneous or not. 
        self.is_photo=is_photo # BOOL of whether rxn is a photolysis rxn or not. 

        self.KPP_line=clean_line # STR with original line from KPP file
        self.rate_function=functions # LIST of STR(s) with rate function(s) used in rxn.
        self.fixed_vars=fixed_vars # List of STR(s) that are "fixed" varables (in mechanism). 
        
        # DICTIONARIES indexed by reactant/product index containing NON-UNIQUE
        # entries for reactants/products if orginally written that way in mech: 
        # e.g.                               for a rxn like 'ClO + ClO = Cl2 + O2' 
        self.indv_rct_cmpds=sdict['rct_cmpds']   #  would be = {0: 'ClO',  1: 'ClO'} 
        self.indv_rct_stoich=sdict['rct_stoich'] #  would be = {0:  1.0,   1: 1.0} 
        self.indv_prd_cmpds=sdict['prd_cmpds']   #  would be = {0: 'Cl2',  1: 'O2'} 
        self.indv_prd_stoich=sdict['prd_stoich'] #  would be = {0:  1.0,   1: 1.0} 
        
        #################################################################################
        # Now get a version that accumulates the stoichiometry for non-unique rct/prds: 
        #################################################################################
        r={}; 
        for num,name in self.indv_rct_cmpds.items(): 
            # If this tracer is NOT in the list of reactant keys already, then add 
            # it pointing to this stoich! 
            if name not in list(r.keys()): 
                r[name]= self.indv_rct_stoich[num] 
            else: # Otherwise add it to what's already there: 
                r[name]=r[name]+self.indv_rct_stoich[num]
        
        p={};
        for num,name in self.indv_prd_cmpds.items():
            # If this tracer is NOT in the list of product keys already, then add 
            # it pointing to this stoich! 
            if name not in list(p.keys()): 
                p[name]= self.indv_prd_stoich[num] 
            else: # Otherwise add it to what's already there: 
                p[name]=p[name]+self.indv_prd_stoich[num] 
        
        self.unq_rct_cmpds=list(r.keys())  # LIST of all UNIQUE reactants in this reaction. 
        self.unq_rct_stoich=r   # DICT with TOTAL stoich of rcts as values, rct_cmpds as keys. 
        self.unq_prd_cmpds=list(p.keys()) # LIST of all UNIQUE products in this reaction. 
        self.unq_prd_stoich=p   # DICT with TOTAL stoich of prds as values, prd_cmpds as keys. 
        
        self.comment_out=False # No commented out rxns exist when readi ninf rom KPP file. 
        
    def copy(self):
        """Create a deep copy of this Mechanism object."""
        import copy 
        return copy.deepcopy(self) 

###############################################################################
#            Functions for defining "F0AM_Reaction class"
#                & Converting to "KPP_Reaction" class. 
###############################################################################

def format_cmt(cmt:str): 
    """Remove leading/trailing spaces from a comment and ensure it does NOT 
    have a leading '%' char already (added elsewhere). Also remove new line characters. """ 
    cmt_stripped=cmt.strip() 
    
    if len(cmt_stripped) >0: 
        if cmt_stripped[0]=='%' and len(cmt_stripped)>1: 
            cmt=cmt_stripped[1:].strip()
        
        cmt_stripped=cmt_stripped.replace('\n', '')
    
    return cmt_stripped 
    
def add_met_to_rates(string:str): 
    # Define a regex pattern that can sperate rate function names from rate function
    # args. This splits the function from the args based on the OUTERMOST matching parenthesis.  
    pattern = re.compile(r'([A-Za-z_][A-Za-z0-9_]*)\s*\((.*)\)')
    
    def _rep(match,string ):  
        """ Helper function for adding 'Met' arg into rate defs"""
        mch_str=match.group(0) # Thing that got matched 
        name = match.group(1) # function name 
        args = match.group(2) # function args 
        
        if name not in ['exp', 'log', 'log10']: # Don't add 'Met' inside math functions! 
            # If ')' occurs BEFORE '(' in the matched "args", you didn't seperate 
            # the function / args correctly. Can happen with strings like: 
            #    'GCARR_ac(3.1e-12) + GCARR_adc(1e-1)' which would match 
            #      name='GACC_ac'   and    args='3.1e-12) + GCARR_adc(1e-1'
            
            if '(' in args and ')' in args and args.index('(',0)> args.index(')',0): 
                # So, re-do the check with a new pattern ending on the FIRST
                # parenthesis instead seperate them correctly.  
                pattern2 = re.compile(r'([A-Za-z_][A-Za-z0-9_]*)\s*\(([^)]*)\)')
                rep_dict2= {match2.group(0): _rep(match2,string) for match2 in re.finditer(pattern2, string)}
                
                # now, replace all the recursive mataches accordingly... 
                result=mch_str;
                for key,value in rep_dict2.items(): 
                    if key in result: 
                        result=result.replace(key,value)
                return result 
            # seperate argument string into arg list: 
            args_list = [arg.strip() for arg in args.split(',') if arg.strip()]
            
            # Add Met as first in arg list: 
            args_list.insert(0, 'Met')
            new_args = ', '.join(args_list)
            result = f"{name}({new_args})"
            
            return result 
        else: 
            # Don't add it to math functions... 
            return string
        
    # Replace all matches in string accordingly & return. 
    rep_dict= {match.group(0): _rep(match, string) for match in re.finditer(pattern, string)}
    
    out_string=string;
    for key,value in rep_dict.items(): 
        if key in out_string: 
            out_string=out_string.replace(key,value)

    return out_string  

def format_rxnrate(rxn_obj, jmap_dict:dict): 
    
    #----------------------------------------------------------------------
    # Ensure any FORTRAN operators in rate functs are in MATLAB syntax
    #----------------------------------------------------------------------
    # NOTE: In standard Fortran, the arithmetic operators like "*"  & 
    # "/" must be immediately adjacent to operands with no intervening 
    # spaces. So,the regex patterns used here are look behind/ahead
    # assertions to ensure the '*' or '/' is not preceded or followed by 
    # whitespace so that we ONLY replace mathematical operators this way  
    
    # Replace division operators '/' with './' for MATLAB
    rate2use = re.sub(r'(?<=\S)/(?=\S)', './', rxn_obj.rate ) 
    
    # Replace multiplication operators '*' with '.*' for MATLAB
    rate2use = re.sub(r'(?<=\S)\*(?=\S)', '.*', rate2use )
    
    #-----------------------------------------------------------------------
    # If its NOT Photolysis,  add 'Met' as an argument to any rate functions:
    #-----------------------------------------------------------------------
   
    if rxn_obj.is_photo==False:

            # If it IS a function, then add 'Met' as an argument for all rate 
            # constant functions so we can indeed use Temp/ M/RH as vars in
            # the rate constant functions 
            rate2use=add_met_to_rates(rate2use)
            rate_cmts=rxn_obj.rate_cmts 
            rate_funct=rxn_obj.rate_function 
            
            # Fixed vars not in met, but rxn rates require them for units to work 
            for fvar in rxn_obj.fixed_vars:
                if fvar in rxn_obj.unq_rct_cmpds: 
                    rate2use=f'({rate2use}).*({fvar}_conc)'
                    rate_cmts=rate_cmts+ f' Must be rate*[{fvar}] for units to work'     
    else: 
        #-----------------------------------------------------------------------
        # If it IS a Photolysis Rxn, we need to map the KPP rate to the F0AM rate.  
        #------------------------------------------------------------------------
        #   Coming out of KPP      -->   Map to Named J Value using jmap_dict
        #   rate = 'PHOTOL(IND)'   -->            rate= 'jO3'
        
        # First, make sure we have a mapped named for PHOTOL(IND) in the jmap dict: 
        if rate2use not in list(jmap_dict['Inds_to_JNames'].keys()):
            raise ValueError("In function, prep4f0am(), we could not find a J-Value defined in F0AM \n" +
                             "within 'jmap_dict' corresponding to: '"+rate2use +
                             "'\n\nThis likely means you need to update 'jmap_dict' to account for new F0AM/FJX mapping j-value changes!")
        else:  # If we do have a match for it,...
        
            # We want to use the name of the rate like jO3 WITHOUT the 
            # surrounding quote characters as the rate to use in F0AM:
            rate2use= jmap_dict['Inds_to_JNames'][rate2use].replace("'", '')
            # And make the rate comment the info about where that J comes from: 
            rate_cmts= jmap_dict['JNames_to_Vals'][rate2use]['Info']
            
            # Set the "rate function" as the Jname WITH the surrouding quote chars: 
            rate_funct=jmap_dict['Inds_to_JNames']
            
    return rate2use, rate_cmts, rate_funct

class foam_reaction(object): 
    """Generator to convert a kpp_rxn class to a foam_rxn class: 

    Attributes
    ----------
    'Rnames'         : STR with rxn "Rnames{i} = '{rxn}'; % {rxn_cmts}"  
    'K_Line'         : STR with rxn rate "k(:,i) = {rate}; % {rate_cmts}"  
    'rate_function'  : LIST with all rate functions referenced in this K_Line
    'G_Line'         : STR listing (indv) reactants: "Gstr{i,1} = 'OH';  Gstr{i,2} = 'HO2';" 
    'F_line'         : STR adjusting stoichiometry for not-fixed/ignored vars: 
                            "fOH(i) = fOH(i)-2.0; fH2O2(i) = fH2O2(i)+1.0;" 

    'comment_out'    : BOOL of whether the rxn should be commented out in final mech or not.     
    'is_het'         : BOOL of whether rxn is heterogeneous or not. 
    'is_photo'       : BOOL of whether rxn is a photolysis rxn or not. 
  
    """
    def __init__(self, rxn_obj, jmap_dict:dict, ignore:list=['hv','M']):
        
        # Decide if rxn should be commented out or not using val assigned in input reaction. 
        self.comment_out=rxn_obj.comment_out
        
        #######################################################################
        # (1) Format Reaction Line:  "Rnames{i}='HNO3+SALAAL=NIT'; % [Potential Comment]" 
        #######################################################################
        # Ensure rxn comment does not have a leading '%'/ new line chars/ spaces around it.  
        rxn_cmt=format_cmt(rxn_obj.rxn_cmts) 
        
        # Format reaction line for F0AM:    
        if len(rxn_cmt)>0: 
           self.Rxn_Line="Rnames{i} = '"+rxn_obj.reaction+f"'; % {rxn_cmt}"
        else: 
           self.Rxn_Line="Rnames{i} = '"+rxn_obj.reaction+ "';"
           
        #######################################################################
        # (2) Format K-Line:  "k(:,i) = GCARR_ac(Met, 4.80e-11, 250.0, Met); % [Potential Comment]" 
        #######################################################################
        # Get info for making k-Line in F0AM from the input KPP.rate object.
        # The wrapper function "format_rxnrate()" does the following: 
        #    (1) Converts FORTRAN math operators to MATLAB syntax
        #    (2) Adds 'Met' as an arg to all input rate functions
        #    (3) Multiplies rates where fixed vars are reactants by fixed var concentrations 
        #    (4) Maps KPP Photolysis rates "PHOTOL(IND)" to named J-rates "jO3" using jmap_dict
        rate2use, rate_cmts, rate_funct=format_rxnrate(rxn_obj, jmap_dict=jmap_dict)
        
        # Ensure rate comment does not have a leading '%'/ new line chars/ spaces around it. 
        rate_cmts=format_cmt(rate_cmts) 
        
        # Format rate line for F0AM:  
        if len(rate_cmts)>0: 
            self.K_Line=f"k(:,i)={rate2use}; % {rate_cmts}"  
        else: 
            self.K_Line=f"k(:,i)={rate2use};"
            
        self.rate_function=rate_funct
        
        #######################################################################
        # (3) Format G-Line:  "Gstr{i,1} = 'OH';  Gstr{i,2} = 'HO2';" 
        #######################################################################
        g_strs=[] # Make empty list to hold a string for each reactant in this rxn. 

        # Get all reactants that AREN't fixed variables in the rxn. 
        mod_rcts=[cmpd for n,cmpd in rxn_obj.indv_rct_cmpds.items() if cmpd not in rxn_obj.fixed_vars]
        
        for n, rct_i in enumerate(mod_rcts):
            # Craft a g-str for it keeping track of index properly (rct 1 or 2?).
            g_strs.append("Gstr{i,"+str(int(n+1))+"}='"+rct_i+"'; ")
            
        # Join all together with a space to get a single G-Line for this rxn.
        self.G_Line =' '.join(g_strs)
        
        #######################################################################
        # (4) Format F-Line:  "fOH(i) = fOH(i)-2.0; fH2O2(i) = fH2O2(i)+1.0;" 
        ####################################################################### 
        fstrs=[] # Initialize a list the hold individual f-strings for rcts/prds in rxn: 
            
        # Loop over unique compounds in reactants & write f str to adjust their stoich: 
        for cmpd, stoich in rxn_obj.unq_rct_stoich.items() :
            # Postivie reactant stoichiometry leads to LOSS of Reactants, so operator 
            # to adjust rct stoichiometry in f_line should be NEGATIVE when stoich is positive. 
            operator='-' if stoich > 0 else '+'
            if cmpd not in rxn_obj.fixed_vars+ignore: # Don't adjust concs of fixed vars! 
                fstrs.append(f"f{cmpd}(i)=f{cmpd}(i){operator}{np.abs(stoich)};")   
            
        # Loop over unique compounds in products & write str to adjust their stoich: 
        for cmpd, stoich in rxn_obj.unq_prd_stoich.items():
        # Postivie product stoichiometry leads to PRODUCTION of PRODUCTS, so operator 
        # to adjust prd stoichiometry in f_line should be POSITIVE when stoich is positive. 
            operator='+' if stoich > 0 else '-'
            if cmpd not in rxn_obj.fixed_vars+ignore: # Don't adjust concs of fixed vars! 
                fstrs.append(f"f{cmpd}(i)=f{cmpd}(i){operator}{np.abs(stoich)};")       
        
        # Join all the strs togeter with a space for the f-line for this rxn: 
        self.F_Line  =' '.join(fstrs) 
        
        
        
    def disp(self, rxn_obj):
        print('Original Rxn: ',rxn_obj.reaction)
        print('Original Rate: ',rxn_obj.rate)
        print('Fixed_Vars: ',  rxn_obj.fixed_vars)


        print(self.Rxn_Line)
        print(self.K_Line)
        print(self.rate_function)
        print(self.G_Line)
        print(self.F_Line)
        return 
           


###############################################################################
#    Functions for defining a general "Mechanism" class: 
###############################################################################
def get_species_list(rxn_list):
    """Get a list of unique species for non-commented out rxns in an input 
    list of reaction objects"""
    species_set = set()
    # Loop over reaction list 
    for rxn_i in rxn_list:
        # Only look at rates required for non-commented out reactions
        if rxn_i.comment_out is not True: 
            # Only add rcts/prds to the species set if you don't have em already! 
            species_set.update(rxn_i.unq_rct_cmpds)
            species_set.update(rxn_i.unq_prd_cmpds)
    
    sp_list= sorted(list(species_set)) 
    
    # Return the species list sored alphabetically 
    return sp_list

def get_nRxns(rxn_list):
    """Calculate the number of (non-commented out) reactions in an input 
    list of reaction objects"""
    n_rxns=0;
    n_rxns= sum([1 for rxn_i in rxn_list if rxn_i.comment_out is not True])
    return n_rxns 

def get_reqJs(rxn_list): 
    """Return a list of all required named Photolysis rates and all required 
    rate constant functions used/referenced in an input list of reaction objects"""
    
    req_Js=set();
    
    # Loop over all reactions in the mechanism: 
    for rxn_i in rxn_list:
        # Only look at rates required for non-commented out reactions
        if rxn_i.comment_out is not True and rxn_i.is_photo is True: 
            # Store photolysis rate names in req_Js
            req_Js.update(rxn_i.rate_function)
         
    return list(req_Js)

def get_reqKs(rxn_list): 
    """Return a list of all required rate constant functions used/referenced 
    in an input list of reaction objects"""
    
    req_Ks=set()
    
    # Loop over all reactions in the mechanism: 
    for rxn_i in rxn_list:
        
        # Only look at rates required for non-commented out reactions
        if rxn_i.comment_out is not True and rxn_i.is_photo is False: 
            req_Ks.update(rxn_i.rate_function)
                
    return list(req_Ks)
            
    
class mechanism(object): 
    """Generator for a number of reactions. Holds information on all reactions 
    within a mechanism easily
    
    Attributes
    ----------
    rxn_list: LIST
        List containing invidivual reaction objects that make up the mechanism. 
    
    species_list: LIST
        List containing unique species involed in any reaction in the mechanism. 
        
    n_rxns: INT
        Number of total reactions in the mechanism. 
        
    n_species: INT
        Number of total species in the mechanism. 
        
    required_Ks: LIST
    
    required_Js: LIST 
    """ 
    def __init__(self, rxn_list=None):
        
        # Intialize vars' we'll set in looping over input rxn list: 
        self.rxn_list = []; 
        
        if rxn_list: 
            # Loop over input reaction list and a dd all reactions to the mech rxn_list
            for rxn_i in rxn_list:
                self.add_rxn(rxn_i)
                   
        # Get list of unqique species (len) participating in non-commented out rxns: 
        self.species_list = get_species_list(self.rxn_list)
        self.n_species = len(self.species_list) 
        
        # Get number of non-commented out rxn in mechanism: 
        self.n_rxns =  get_nRxns(self.rxn_list)

        # Get lists of unique Js/Ks referened in rate functions of non-commented out
        # reactions that are required for the mechanism to compile 
        self.req_Js=get_reqJs(self.rxn_list)
        self.req_Ks=get_reqKs(self.rxn_list)
        
        # Indicate if the mechanism has ANY commented out reactions: 
        self.includes_cmts=any([r.comment_out for r in self.rxn_list])
     
    def copy(self):
        """Method to create a deep copy of this Mechanism object."""
        import copy 
        return copy.deepcopy(self) 
    
    def add_rxn(self, rxn_obj):
        """Method to add a reaction to a mechanism. Input is a reaction object. """
        if not any(are_reactions_equal(rxn_obj, r) for r in self.rxn_list):
            self.rxn_list.append(rxn_obj)
            
            # After adding thix rxn, update all the other attributes of the
            # mechanism accordingly. but, won't change if rxn_obj.comment_out is True. 
            self.species_list = get_species_list(self.rxn_list) 
            self.n_species = len(self.species_list) 
            self.n_rxns = get_nRxns(self.rxn_list)
            self.req_Js=get_reqJs(self.rxn_list)
            self.req_Ks=get_reqKs(self.rxn_list)
            self.includes_cmts=any([r.comment_out for r in self.rxn_list])
        else:
            print(f"Duplicate reaction not added: {rxn_obj.reaction}")

    def remove_rxn(self, rxn_obj):
        """Method to remove a reaction from a mechanism. Input is a reaction object."""
        
        rm=False # initalize var to let us know if anything was removed: 
            
        # Loop over all reactions currently in the rxn list: 
        for existing_rxn in self.rxn_list:
            
            # Test to see if this rxn in our mech is equal to the input one or not: 
            if are_reactions_equal(rxn_obj, existing_rxn):
                
                # If so, then remove the rxn from our list & print it: 
                self.rxn_list.remove(existing_rxn)
                #print(f"Removed reaction: {existing_rxn.reaction}")
                rm=True
                break
            
        # If you removed a rxn, update all the other attributes of the
        # mechanism accordingly. But, won't change if existing_rxn.comment_out was True.
        if rm ==True: 
            self.species_list = get_species_list(self.rxn_list) 
            self.n_species = len(self.species_list) 
            self.n_rxns = get_nRxns(self.rxn_list)
            self.req_Js=get_reqJs(self.rxn_list)
            self.req_Ks=get_reqKs(self.rxn_list)
            self.includes_cmts=any([r.comment_out for r in self.rxn_list])
        else: 
            # Otherwise, let em' know you never found a match: 
            print(f" Could NOT remove: {rxn_obj.reaction}...  not found in mechanism's list of reactions.")

    def comment_rxn(self, rxn_obj):
        """Method to comment out a reaction from a mechanism. Input is a reaction object."""
        
        changed=False # initalize var to let us know if anything was commented: 
            
        # Loop over all reactions currently in the rxn list: 
        for existing_rxn in self.rxn_list:
            
            # Test to see if this rxn in our mech is equal to the input one or not & f
            if are_reactions_equal(rxn_obj, existing_rxn): 
                
                # If the rxn was NOT already commented out... then: 
                if not_cmt(existing_rxn):
                    # Remove the non-commented out verson from the list 
                    self.remove_rxn(existing_rxn)
                    # Add back in a commented out version of this reaction: 
                    self.add_rxn(update_attr(existing_rxn, 'comment_out', True))
                    changed= True 
                # Beack cuze either changed it or it didn't need to be changed: 
                break
            #else:
                # Otherwise, let em' know you never found a match to change 
                #print(f"{rxn_obj.reaction} not found in mechanism's list of reactions.")
            
            # If you changed a rxn, update all the other attributes of the
            # mechanism accordingly. 
            if changed ==True: 
                self.species_list = get_species_list(self.rxn_list) 
                self.n_species = len(self.species_list) 
                self.n_rxns = get_nRxns(self.rxn_list)
                self.req_Js=get_reqJs(self.rxn_list)
                self.req_Ks=get_reqKs(self.rxn_list)
        
    def select_by_rcts(self, rct_list:list):
          """Return a sub-selection of the mechanism with only rxns involving specific reactants"""
          rxn_list=[]
          for rxn_i in self.rxn_list:
              if any([rct_i in rxn_i.unq_rct_cmpds for rct_i in rct_list]):
                  rxn_list.append(rxn_i)
          rct_mech=mechanism(rxn_list=rxn_list)
          return rct_mech  
      
    def select_by_prds(self, prd_list:list):
          """Return a sub-selection of the mechanism with only rxns involving specific products"""
          rxn_list=[]
          for rxn_i in self.rxn_list:
              if any([prd_i in rxn_i.unq_prd_cmpds for prd_i in prd_list]):
                  rxn_list.append(rxn_i)
          prd_mech=mechanism(rxn_list=rxn_list)
          return prd_mech 
      
    def select_by_cmpds(self, cmpd_list:list):
          """Return a sub-selection of the mechanism with only rxns involving specific compounds"""
          rxn_list=[]
          for rxn_i in self.rxn_list:
              if any([cmpd_i in rxn_i.unq_prd_cmpds+rxn_i.unq_rct_cmpds for cmpd_i in cmpd_list]):
                  rxn_list.append(rxn_i)
          prd_mech=mechanism(rxn_list=rxn_list)
          return prd_mech 
      
    def select_photolysis_rxns(self):
        """Return a sub-selection of the mechanism with phototlysis rxns only"""
        return mechanism(rxn_list=[rxn for rxn in self.rxn_list if rxn.is_photo])
    
    def select_het_rxns(self):
        """Return a sub-selection of the mechanism with heterogeneous rxns only"""
        return mechanism(rxn_list=[rxn for rxn in self.rxn_list if rxn.is_het])
    
    def select_nhet_rxns(self):
        """Return a sub-selection of the mechanism with non-het rxns only"""
        return mechanism(rxn_list=[rxn for rxn in self.rxn_list if not rxn.is_het])
    
    def select_gas_rxns(self):
        """Return a sub-selection of the mechanism with gas phase rxns only"""
        return mechanism(rxn_list=[rxn for rxn in self.rxn_list if rxn.is_gas])
    
    def select_3body_rxns(self):
        """Return a sub-selection of the mechanism with 3-body rxns with M only"""
        return mechanism(rxn_list=[rxn for rxn in self.rxn_list if rxn.is_3body])
    
    def select_without_cmts(self):
        """Return a mechanism with only reactions that aren't commented out included:""" 
        return mechanism(rxn_list=[r for r in self.rxn_list if not_cmt(r)])
    
    def disp_rxn_list(self):
        [print(r.reaction) for r in self.rxn_list]
        return 
    
    def get_unq_prds(self):
        """Return a list of unique products in the mechanism: """
        prods=[]
        for rxn_i in self.rxn_list:
            [prods.append(prd_i) for prd_i in rxn_i.unq_prd_cmpds if prd_i not in prods]
        
        return prods
    
    

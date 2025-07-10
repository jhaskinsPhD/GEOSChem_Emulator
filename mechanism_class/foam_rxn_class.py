#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module containing classes & functions that define and
help manipulate pulling info about a given reaction string and mechanism. 

@author: u6044586
"""
import re 
import numpy as np

###############################################################################
#    Functions for defining "Reaction"  Class: 
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

def get_Rnames_line(rxn_obj:reaction): 
    if len(rxn_cmt) > 0: 
        Rnames="Rnames{i}='{rxn_obj.rxn_str}'; %{rxn_obj.rxn_cmt}' 
    else: 
        Rnames="Rnames{i}='{rxn_obj.rxn_str}'; %{rxn_obj.rxn_cmt}' 

def parse_indv_reaction_str(rxn_str:str, rxn_arrow:str='=',rcts2ignore:list=['hv'], 
                            prds2ignore:list=[]):
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
    out={'reaction':rxn_str, 'reactants':{}, 'products':{}} 
    
    # Loop over each tracer/stoich string in the reactant half of the string: 
    is_photo= None; is_3body= None
    for pair in rct_pairs: 
        # Seperate the stoich from tracer in this individual pair
        stoich_i, tracer_i= sep_stoich_n_tracer(pair) 
        
        if tracer_i =='hv': is_photo=True 
        if tracer_i =='M': is_3body=True 
        
        # But only add this tracer & its stoich if its NOT in the list of rcts to ignore... 
        if tracer_i not in rcts2ignore: 
            # If this tracer isn't in the list of reactant keys already, then add it pointing to this stoich! 
            if tracer_i not in list(out['reactants'].keys()): 
                out['reactants'][tracer_i]= stoich_i
                
            # If it is already there, we need to add up the stoichiometry... 
            else:
                old_stoich = out['reactants'][tracer_i]
                out['reactants'][tracer_i]= old_stoich+ stoich_i
                     
        
    # Loop over each tracer/stoich string in the reactant half of the string: 
    for pair in prd_pairs: 
        # Seperate the stoich from tracer in this individual product stoich/tracer pair
        stoich_i, tracer_i= sep_stoich_n_tracer(pair) 
                
        # Only add this tracer & its stoich if its NOT in the list of prds to ignore... 
        if tracer_i not in prds2ignore: 
            
            # If this tracer is NOT in the list of product keys already, then add 
            # it pointing to this stoich! 
            if tracer_i not in list(out['products'].keys()): 
                out['products'][tracer_i]= stoich_i
                
            # If it is already there, we need to add up the stoichiometry... 
            else:
                old_stoich = out['products'][tracer_i]
                out['products'][tracer_i]= old_stoich + stoich_i
                    
    return out, is_photo, is_3body
 
class kpp_reaction(object):
    """Generator for a single reaction. Holds a reaction's reactants, products, 
    and rate priding a way to retrieve these items easily. 

    Attributes
    ----------
    'reaction'       : STR with reaction 
    'rate'           : STR with full rate (func & arg) for this reaction 
    'comments'       : STR with comments about this reaction (or blank). 
    'is_gas'         : BOOL of whether rxn is gas phase or not. 
    'is_het'         : BOOL of whether rxn is heterogeneous or not. 
    'is_photo'       : BOOL of whether rxn is a photolysis rxn or not. 
    'is_3body'       : BOOL of whether rxn includes 3rd body reactant, M, or not. 
    'rct_cmpds'      : LIST of all reactants in this reaction.  
    'rct_stoich'     : Dict where keys are th tracer names of products (in prd_cmpds) and values 
                            are the stoichiomentry of that specific product. 
    'prd_cmpds'      : LIST of all products in this reaction. 
    'prd_stoich'     : Dict where keys are th tracer names of products (in prd_cmpds) and values 
                            are the stoichiomentry of that specific product. 
    'KPP_line'       : STR with original line from KPP file
    'rate_function'  : STR with rate function (only- no args or constants) used. 
                
    """

    def __init__(self, kppline, rxn_arrow:str='=',rcts2ignore:list=['hv'], prds2ignore:list=[]):
        
        out, is_photo, is_3body= parse_indv_reaction_str(rxn_str, 
                                     rxn_arrow=rxn_arrow,
                                     rcts2ignore=rcts2ignore, 
                                     prds2ignore=prds2ignore)
        
        
        
                
                
                
        
        
        self.rxn_str = rxn_str
        self.rate_str = rate_str
        self.rct_cmpds = list(out['reactants'].keys())
        self.rct_stoich = out['reactants']
        self.prd_cmpds = list(out['products'].keys())
        self.prd_stoich = out['products']
        
        if is_photo==True: 
            self.rxn_type='photolysis'
        elif is_3body ==True: 
            self.rxn_type='3-body'
        else: 
            self.rxn_type=None
            
###############################################################################
#    Functions for defining "Mechanism" Class: 
###############################################################################
def get_species_in_mech(rxn_list):
    """Get a list of unique species in a list of reaction objects"""
    species_set = set()
    for rxn_i in rxn_list:
        species_set.update(rxn_i.rct_cmpds)
        species_set.update(rxn_i.prd_cmpds)
    return sorted(list(species_set))

def calc_rxn_similarity(rxn1, rxn2): 
    ###########################################################################
    #                            Check/Calc reactant similarity: 
    ###########################################################################
    # Check if the 2 rxns have the same reactant compounds: 
    same_rcts =True if set(rxn1.rct_cmpds)==set(rxn2.rct_cmpds) else False
    
    # Calculate the degree of similaritiy of the reactant compounds' stoichiometry: 
    comp_rct_stoich=[]
    for cmpd in list(set(rxn1.rct_cmpds+rxn2.rct_cmpds)):
        
        # If this compound is a reactant in both rxns AND it's stoich is same in both: 
        if ((cmpd in list(rxn1.rct_stoich.keys())) and 
            (cmpd in list(rxn2.rct_stoich.keys()))and 
            (rxn1.rct_stoich[cmpd]==rxn2.rct_stoich[cmpd])): 
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
    same_prds =True if set(rxn1.prd_cmpds)==set(rxn2.prd_cmpds) else False

    # Calculate the degree of similaritiy of the product compounds' stoichiometry: 
    comp_prd_stoich=[]
    for cmpd in list(set(rxn1.prd_cmpds+rxn2.prd_cmpds)):
        
        # If this compound is a reactant in both rxns AND it's stoich is same in both: 
        if ((cmpd in list(rxn1.prd_stoich.keys())) and 
            (cmpd in list(rxn2.prd_stoich.keys()))and 
            (rxn1.prd_stoich[cmpd]==rxn2.prd_stoich[cmpd]) ): 
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
        (same_prds==True) and (prd_similarity==100) and (rxn1.rate_str==rxn2.rate_str)): 
        return True
    else:
        return False
    
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
    """ 
    def __init__(self, rxn_list=None):
        self.rxn_list = []
        if rxn_list:
            for rxn in rxn_list:
                self.add_reaction(rxn)
        self.species_list = get_species_in_mech(self.rxn_list)
        self.n_rxns = len(self.rxn_list)
        self.n_species = len(self.species_list)

    def add_reaction(self, rxn_obj):
        """Function to add a reaction to a mechanism. Input is a reaction object. """
        if isinstance(rxn_obj, reaction):
            if not any(are_reactions_equal(rxn_obj, r) for r in self.rxn_list):
                self.rxn_list.append(rxn_obj)
                self.species_list = get_species_in_mech(self.rxn_list)
                self.n_rxns = len(self.rxn_list)
                self.n_species = len(self.species_list)
            else:
                print(f"Duplicate reaction not added: {rxn_obj.rxn_str}")

    def remove_reaction(self, rxn_obj):
        """Function to remove a reaction from a mechanism. Input is a reaction object."""

        if isinstance(rxn_obj, reaction):
            for existing_rxn in self.rxn_list:
                if are_reactions_equal(rxn_obj, existing_rxn):
                    self.rxn_list.remove(existing_rxn)
                    print(f"Removed reaction: {existing_rxn.rxn_str}")
                    break
            else:
                print(f"{rxn_obj.rxn_str} not found in mechanism's list of reactions.")
            
            self.species_list = get_species_in_mech(self.rxn_list)
            self.n_rxns = len(self.rxn_list)
            self.n_species = len(self.species_list)
            

            
            
    


            
            
            


        
        
    


        

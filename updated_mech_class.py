#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 10 18:21:48 2025

@author: u6044586
"""
import re 

##############################################################################
#               Functions for defining "KPP Reaction"  class: 
###############################################################################
def sep_stoich_n_tracer(inp_str:str): 
    """ Seperates the stoichiometry from a tracer name in an input string. 

   **Strengths & Assumptions: **
       - Assumes the 1st char in input string is the stoichiometry sign (`+` or `-`) 
       - Assumes 1st alphabetic character it encountres is start of tracer name
             & anything between sign/tracer name is stoichiometry. 
       - Works well for strings like '+1.920CH2O', '-0.26GLYX', or '+CH2O'

   **Potential Pitfalls / Limitations:**
   
       - Will not correctly parse tracer names containing special characters 
       such as underscores, hyphens, or other non-alphanumeric characters.
       - Will interpret the characters before the first alphabetic as stoichiometry, 
       which may be incorrect if the string format varies unexpectedly.
       - If `pair` has no alphabetic character, the method will error or produce 
       unintended results.
       - Assumes that `pair` has been pre-validated to contain a sign at 
       the start and valid characters

   Example Usage:
       stoich, tracer = sep_stoich_n_tracer('+CH4')       # stoich=1.0,   tracer='CH4'
       stoich, tracer = sep_stoich_n_tracer('+0.260GLYX') # stoich=0.260, tracer='GLYX'
       stoich, tracer = sep_stoich_n_tracer('-1.5DUMMY')  # stoich=-1.5,  tracer='DUMMY'

   """
   
    #---------------------------------------------------------------------------
    # Check that input is formatted as expected (safe guards against bad parsing!)
    #---------------------------------------------------------------------------
    if not isinstance(inp_str, str):
         raise TypeError(f"Input must be a string, but got {type(inp_str)}.")
         
    inp_str=inp_str.replace(' ' ,'') # remove spaces. 
    
    if len(inp_str) <= 1:
        raise ValueError("Input string {inp_str } is too short. Must be >len(1) including sign, (optional stoich), & tracer.")

    # Validate assumption that string always has a sign character at start
    sign = inp_str[0]; str2sep = inp_str[1:] # tracer stoich/name is rest. 
    if sign not in ['+', '-']:
        raise ValueError(f"Input must start with '+' or '-', but got '{sign}' in '{inp_str}'.")
    
    # Check for invalid chars that would mess up splitting tracer/stoich using .isalpha()... 
    if not str2sep.replace('.', '').isalnum():
        raise ValueError(f"Input '{str2sep}' contains invalid characters. Allowed: digits, letters &  '.'")
    
    # Find the first alphabetic character in the string (start of tracer name). 
    tracer_start = None
    for idx, ch in enumerate(str2sep):
        if ch.isalpha():
            tracer_start = idx
            break
        
    # Make sure we found the start to a tracer name: 
    if tracer_start is None:
        raise ValueError(f"No alphabetic characters found in '{inp_str}'. Cannot determine tracer name.")

    # Extract the potential stoich string and tracer name 
    stoich_str = str2sep[0:tracer_start]
    tracer = str2sep[tracer_start:]

    # Validate stoichiometry is a number & set to 1 if none was included in front of tracer. 
    if len(stoich_str) == 0:
        stoich = float(1.0)
    else:
        try:
            stoich = float(stoich_str)
        except ValueError:
            raise ValueError(f"Cannot convert '{stoich_str}' to float in '{str2sep}'.")
            
    # Update stoich to be negative if it's incoming sign indicated it should be. 
    stoich=-1*stoich if sign=='-' else stoich
              
    return stoich, tracer 
       

  
def split_halfrxn(rxn_half:str, allow_neg_stoich:bool=False):
    """
    Splits a reaction half into tokens with signs, e.g.,
    'TRO2+1.920CH2O-GLYX' -> ['+TRO2', '+1.920CH2O', '-GLYX']
    
    Args: 
    -----
         rxn_half - <STR> 1/2 of a chemical reaction string that's been split 
                          at the "reaction arrow". What's before = reactants, 
                          whats after = products. 
        
        allow_neg_stoich: <BOOL>  Physically, rxns can NEVER have "negative stochiometry" , 
                          but some model mechanisms allow this. Default is FALSE.
                          (e.g. assume all stoich is '+'). Setting to TRUE if '-' 
                          not there not an issue unless tracer names include '-'. 
    
    Returns: 
    --------
        
        rxnhalf_parts - <LIST> Where each item is a stochiometry/tracer names pair 
                         that includes a preceeding '+' or '-' indicating the 
                         sign of the stoichiometry as written in the input string.
                         
    AUTHOR: 
    -------
    Prof. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD
    """
    
    # Make sure the input string doesn't have any spaces: 
    rxn_half=rxn_half.replace(' ' ,'') 
    
    # Only split the input string at a '-' char if we know it can have neg 
    # stoichimetry( otherwise, the '-' might be part of the tracer name). 
    allowed_delims=['+','-']  if allow_neg_stoich is True else ['+'] 
        
    # Add regex escapes to each delimiter in the input list & combine into a
    # single regex "OR" statement... resulting in delim_str = '\\+|\\-' 
    delim_str = '|'.join(map(re.escape, allowed_delims))  
    
    # Split the string just before/after any delimiter. 
    #  IN:  rxn_half =     'TRO2 + 1.920CH2O + 0.215MGLY + OH - DUM'  
    # OUT:  tokens= ['TRO2', '+','1.920CH2O','+','0.215MGLY','+','OH','-','DUM'] 
    pattern = f'(?<=[{re.escape(delim_str)}])|(?=[{re.escape(delim_str)}])'
    tokens = re.split(pattern, rxn_half)
    
    # If there was no delimiter preceeding the FIRST item in the list, tokens, then 
    # assume it was '+' (e.g. all stoichiometry is positive unless noted!)
    if tokens and tokens[0] not in allowed_delims:
       tokens.insert(0, '+')
       
    rxnhalf_parts=[] # Initialize the output list

    # Loop over every item in the resulting tokens list. 
    for t in tokens:
    
        # If this item is one of our allowed delimiters ('+' or '-'), keep it. 
        if t in allowed_delims:
            prev_sign = t
        # If not, its a stoich/tracer pair, so combine it with previous sign & 
        # add that combined str to the output list: 
        elif t: 
            rxnhalf_parts.append(prev_sign + t)
    
    # Check to make sure you got it right. Everything but the first item in 
    # results should exist in the spaceless input str - otherwise you messed up.  
    # and no strings with just two delims should be in there at all. 
    if (not all([combo in rxn_half for combo in rxnhalf_parts[1:]]) and 
        all([mistake!=rxnhalf_parts for mistake in ['+', '-', '++' , '--', '+-']])):
        print('Uh oh!') 
    
    return rxnhalf_parts

def parse_indv_reaction_str(rxn_str:str, rxn_arrow:str='=',allow_neg_stoich:bool=False,
                            ignore:list=['M','hv']):
    """Function to parse a chemical reaction string and seperate the stoichiometry 
       from the tracer names for reactants/products seperately. N
       OTE: Function 
       won't work properly if tracer names can begin with numbers. 
    
    Inputs: 
    -------
        rxn_str - A STRING representing a chemical reaction that you'd like to 
                  seperate into reactant & product name/stoichioemtry pairs: 
                (e.g. 'TOLU + OH = TRO2 + 1.920CH2O + 0.260GLYX + 0.215MGLY + OH')
                     
        rxn_arrow - A STR that indciates what character sequence seperates the 
                    reactants from products in the input rxn strings. Default is 
                    set to '='  but may be ':', '->', '-->', '<->' in some mechanisms. 
        
        allow_neg_stoich: <BOOL> Indicating if a minus  sign is allowed to delimit 
                          pairs of stochiometry/compounds within the rxn half or not. 
                          Physically, rxns can NEVER have "negative stochiometry" , 
                          but some model mechanisms allow this. Default is FALSE
                          (e.g. assume all stoich is '+'). 
                       
        ignore - A list of STRS to ignore in output lists if they are parsed 
                      as tracernames in the reactants/products. Can be useful if 
                      you don't want to parse 'O3+hv=O1D' as having a reactant 
                      named 'hv' with stoich=1. 
                      
                      
    Outputs: 
    -------
        2 Dictionares, "products" and "reactants" each with keys "cmpds" and "stoich" 
        pointing to a nested dict that has keys with correspondng indexes and 
        values as strings of compounds names or floats with that compounds' stoichiomtry. 
        
        Example: 
           'TOLU + OH = TRO2 + 1.920CH2O + 0.260GLYX + 0.215MGLY + OH - DUM'
            
            reactants['cmpds']   = {0: 'TOLU', 1:  'OH'}
            reactants[['stoich'] = {0: 1.0,    1:  1.0}
            products['cmpds']    = {0: 'TRO2', 1:  'CH2O', 2: 'GLYX', 3:'MGLY', 4: 'OH', 5: 'DUM' }
            products['stoich']   = {0: 1.0,    1:  1.920,  2: 0.260,  3: 0.215, 4: 1.0,  5: -1.0 }              
  
    """
    # Ensure the input string has no spaces: 
    rxn_str=rxn_str.replace(' ', '')
    
    # Count occurrences of rxn_arrow
    count = rxn_str.count(rxn_arrow)
    if count == 0:
        raise ValueError(f"CAN'T PARSE REACTION. Reaction arrow delimiter '{rxn_arrow}' could NOT be found in the reaction string:\n\t{rxn_str}")
    elif count > 1:
        raise ValueError(f"CAN'T PARSE REACTION. Multiple reaction arrows '{rxn_arrow}' were found in reaction string:\n\t{rxn_str}")

    # Since only 1 occurrence is allowed, safe to split reaction string on the reaction arrow 
    rct_half,prd_half =rxn_str.split(rxn_arrow)
    
    # Get list of stoichiometry/tracer names in reactants (with signs of stoich!) 
    rct_pairs= split_halfrxn(rct_half, allow_neg_stoich=allow_neg_stoich ) 
    
    # Get list of stoichiometry/tracer names in reactants (with signs of stoich!) 
    #     [ '+TRO2', +'1.920CH2O','+0.260GLYX',-OH'] 
    prd_pairs= split_halfrxn(prd_half,allow_neg_stoich=allow_neg_stoich )
                
    # Prepare output dictionaries
    reactants={'cmpds':{},'stoich':{}}
    products={'cmpds':{},'stoich':{}}
        
    # Seperate stoichiometry from tracers in reactants: 
    r_ct=0                  
    for pair in rct_pairs: 
        # Seperate the stoich from tracer in this individual pair
        stoich_i, tracer_i= sep_stoich_n_tracer(pair) 
        
        # But only add this tracer & its stoich if its NOT in the list of cmpds to ignore... 
        if tracer_i not in ignore: 
            reactants['cmpds'][r_ct]= tracer_i
            reactants['stoich'][r_ct]= stoich_i
            r_ct += 1 
         
    p_ct=0 # Do the same fro the products: 
    for pair in prd_pairs: 
        stoich_i, tracer_i= sep_stoich_n_tracer(pair) 
        if tracer_i not in ignore: 
            products['cmpds'][p_ct]= tracer_i
            products['stoich'][p_ct]= stoich_i
            p_ct += 1
                    
    return reactants, products


def test_rxn_parsing():
    """Debugging Function to see if things in this script are working as expected""" 
    arrows=['=','->','-->','=>', '<->', '<-->', ':']
    arr=arrows[1]
    
    test_strings = [
        # Simple Valid reactions with dif .arrows (should pass) 
        f"X + YY {arr} XYZ",                   #letters only. 
        f"H2 + O2 {arr}  H2O",                 # addding #s between tracers. 
        f"MTPA + 2OH {arr} 0.075LIMO2 + 0.67APINO2 +0.255BPINO2", # adding stoicihoemtry: 
        f"O3 + hv {arr}  O2 + O",              # with hv 
        f"SO2  + SALAAL + O3  {arr} SO4 - SALAAL",   # negative stoichiometry. 
        
        # Edge cases:  ??? 
        "TOLU +OH -> TRO2 - CH4",            # Dash in arrow & in stoichiometry. 
        "TOLU +  .3OH -->  TRO2 - CH4 ",     # Dash in arrow & in stoichiometry. 
        f"2A+-3CH2O{arr}TR2 +  3CH4   ",     # multiple signs
        f"X + --Y {arr} Z",                  # double minus
        f"A +B {arr} C ++ D",                # double plus
        
        # Sign and illegitimate characters (should error)
        f"TOLU +OH {arr} 1.2.TRO2 - CH4",          # Too many .s 
        f"TOLU + O3 +{arr} TRO2 + CH4",            # += is malformed; should cause error
        f"TOLU + O3 {arr} 0.123TRO2 + CH4++",     # trailing sign plus
        f"TOLU + O3 {arr} TRO2 + CH4$"  ,         # invalid character '$'
        f"4TOLU + H_OH {arr} TRO2 + CH4",         # invalid character '_'
        "I2O2 {+M} -> 0.008IO + 0.996I + {+M}",  # invalid character '{' 
    
        # Multiple arrows (should error) 
        f"TOLU + 2O{arr}TR2+O2{arr}DUM"]
    


    for i, test_str in enumerate(test_strings, 1):
        try:
            reactants, products = parse_indv_reaction_str(test_str, rxn_arrow=arr, allow_neg_stoich=True)
            rcts=' , '.join(['<'+str(reactants['stoich'][i])+'> "'+reactants['cmpds'][i]+'"' for i in list(reactants['cmpds'].keys())]) 
            prds= ' , '.join(['<'+str(products['stoich'][i])+'> "'+products['cmpds'][i]+'"' for i in list(products['cmpds'].keys())]) 
            print(f"Test {i}: PASS     '{test_str}'")
            print(f"  Reactants: {rcts}")
            print(f"  Products: {prds}\n")
        except Exception as e:
            print(f"Test {i}: FAIL     '{test_str}'")
            print(f"  Error: {e}\n")
                
    return 

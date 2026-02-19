#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 11:44:28 2023

@author: u6044586
"""
import os 
import sys 
import inspect
import re 
import numpy as np 
from collections import defaultdict

###############################################################################
# Functions used to convert math in fortran files to MATLAB syntax. 
# Called from within parse_rxns_and_rates()from within read_kpp(). Here because 
# they are also called from make_import_GC_rates() 
###############################################################################

def _convert_fortran_operators(in_string): 
    """Function to take an input string and convert (select) FORTRAN style 
     mathematical function names within it into modern function names. 
     
     Specifically: 
        - Converts function names EXP, LOG, LOG10 to lowercase.
        - Converts ".FALSE. to 'false' and '.TRUE.' to 'true' 
        
    REQUIREMENTS: import re 
    
    USAGE: Called within parse_rxns_and_rates() from within read_kpp() 
           to clean up any lines with rate functions containing these math operators  
    """

    # Define a list of tuple pairs holding the function name & a regex pattern to match it: 
    operators=[('EXP', r'EXP\((.*?)\)'),      # Pattern to match 'EXP(' ... ')'
               ('LOG', r'LOG\((.*?)\)'),      # Pattern to match 'LOG(' ... ')'
               ('LOG10', r'LOG10\((.*?)\)')   # Pattern to match 'LOG10(' ... ')'
               ] 
        
    # Initialize "out_string" as copy of input (what's returned if no subs needed).
    out_string=in_string
        
    # Loop over some common math function names & convert them to lowercase 
    # while maintaining function arguments entirely. 
    for fname, funct_pattern in operators: 
        # Use the re.sub() function from the regex module to do all this at once: 
        #   
        #    (1) Find all substrings in the input string that match the 
        #        fortran function pattern. 
        #
        #    (2) For each match, call the lambda function, passing in the match 
        #        object to generate a new replacement substring. It converts the 
        #        function name to lowercase using 'fname.lower()' while maintaining 
        #        the original arguments to the function with 'm.group(1)'.
        #
        #    (3) Replace each matched substring with the new substring returned 
        #        by the lambda function. 
        #
        #    (4) Assign the result to a single output string after all substring
        #        replacements (on this loop) have been completed. 
        #
        # This results in replacing  'EXP(1.23)'   with 'exp(1.23)'     OR ... 
        #                            'LOG10(4.56)' with 'log10(4.56)'

        out_string = re.sub(funct_pattern,lambda m: f"{fname.lower()}({m.group(1)})",out_string)
    
    # Now just do some simple string replacement for the true/false operators... 
    out_string=out_string.replace('.FALSE.', 'false')
    out_string=out_string.replace('.TRUE.', 'true')
        
    return out_string
 
        
def _convert_fortran_numbers(in_string:str ):
    """ Function to take an input string and convert Fortran-style numeric literals 
        into modern scientific notation syntax. Used to convert numbers in KPP 
        files to something you can write to a MATLAB/python/CSV/excel file. 
        Used to convert strings with reaction rates like:
            
                     In KPP File:                        as Output 
            'GCARR_ac(1.00d-14, -490.0d0)'   to   'GCARR_ac(1.00e-14, -490.0)'

        Can handle: 
        - Fortran exponents like 'd', 'D', 'e', 'E' (e.g., '1.23d-4', '5.67D+8')
        - Removes the '_dp' suffix used in Fortran for double precision #s

    INPUTS: 
    ------  
        (1) in_string - single STR containing (maybe) some numerical FORTRAN notation 

    OUTPUTS: 
    --------
        (2) out_string - updated STR with FORTRAN numbers converted to (modern) math.

    USAGE: 
    ------ 
        Called Within: parse_rxns_and_rates() from within read_kpp() from 
                       within make_GEOSCHEM_AllRxns_file(). Can be used as standalone. 

    REQUIREMENTS: 
    -------------
        LIBRARIES:           import re

        CUSTOM FUNCTIONS:    nested/inner/helper function: _fnum_replace() 

    AUTHOR: 
    -------
        Prof. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD
    """
    ###########################################################################
    #  Part I: Transform Fortran Style NUMBERS into MATLAB style numbers... 
    ###########################################################################

    # Define a regular expression pattern to match FORTRAN-style numbers: 
    fnum_pattern = re.compile(r'([-+]?\d*\.\d+|[-+]?\d+)([dDeE]([+-]?\d+))?(_dp)?')
    
    #  ^^This will match numbers that may start with an optional + or - sign.
    # The number can either be a decimal or an integer. It can have an exponent 
    # part, which starts with either: 'd', 'D', 'e', or 'E', followed by an optional
    # '+' or '-', and then digits (for example, d-4 or E+10). Finally, it can 
    # optionally end with the suffix '_dp', which is used in Fortran to indicate 
    # double precision numbers. After finding these, string subs are only done 
    # when the number matched has a dif format (eg. doesn't mess with integers).
    # Example matches &(subsequent) substitions below: 
    #      '42'        -->  '42'       (intergers are same).
    #      '10.3e5'    -->  10.3e5'    (no update needed). 
    #      '-5.67d-8'  -->  '-5.67e-8' 
    #      '3.29D+00'  -->  '3.29'  # (when exp evals to 1, removes exp)
    #      '-5.0E-3_dp'-->  '-5.0e-3'

    def _fnum_replace(match):
        """
        Nested/Helper/Local function to convert a matched Fortran number to 
        more modern scientific notation format. 
        
        Inputs:  match - a re.Match object (not tuple/ it should have methods 
                                            like .group() within it). 
        Outputs - a replacement string that has: 
                        - exponents noted ONLY with 'e' (not 'd', 'D', or 'E"). 
                        - no '_dp' suffix 
                        - no exponent if the exponent was zero (e.g.'E+00','d0', 'e0')
                        - an exponent sign ONLY if negative (e.g. no e+4).  
        """
        
        # Extract the relevant groups from the (individual) match:
        number = match.group(1)        # base # (e.g., '1.23', '-4.56')
        exponent = match.group(3)      # Sign and number pf exponent only (e.g., '-2', '+3', or None)
        
        if exponent is None or exponent == '':
            # Case 1: There is no explicit 'e'/'E'/'d'/'D'in the matched pattern 
            #         Occurs for plain old numbers like '-3.14' or '+42'.
            #         So, just return the number itself! 
                return number
        else:
            # Case 2: The number has an explicit exponent (e.g. "1.23e-4", "5.67E+8")
            #          Convert exponent str to integer is (e.g. -4, +8, 0, etc.). 
            exponent_value = int(exponent)
            if exponent_value == 0:
                # Case 2a: If the value of exponent is 0, just remove it / only return number. 
                #          (e.g. '2.14d0' --> '2.14' or '42eE+00' --> '42'). 
                return number
            else:
                # Case 2b: If the exponent is not 0/''/None & does exist, then 
                #          return it using the 'e' notation between the # & exponent 
                #          VALUE. Only retains sign of exponent sign if negative. 
                #          (e.g. '3.14d-4'--> '3.14e-4' but '3.14e+4'--> '3.14e4'
                return f"{number}e{exponent_value}"

    # Use the re.sub() function from the regex module to: 
    #    (1) Find all substrings in the input string that match the fortran number pattern. 
    #    (2) For each match, call the internal _fnum_replace() function, passing 
    #        in the match object to generate a new replacement substring (bye Fortran syntax!)
    #    (3) Replace each matched substring with the new substring returned by _fnum_replace(). 
    #    (4) Assign the result to a single output string after all substring replacements. 
    
    out_string = fnum_pattern.sub(_fnum_replace, in_string)

    return out_string
  
    
 ###############################################################################
 # Functions used to convert from categorical dict to record index dict... 
 ############################################################################### 


def _convert_to_index_dict(category_dict):
    """
    Convert categorical dict into a record based index dictionary. 
        
   (INPUT) 'category_dict':  (Useful for grouping info abou all rxns by category)
       
       - KEYS = string CATEGORIES, VALUES = a LIST of items.
       - All lists under each key have equal lengths.
       - Elements at a particular index in each list corresponds to 
         a matching record in any other list from a dif cateogory key. 
        
         #               INPUT DICT. STRUCTURE: 
         #   category key   = list[index_0 ,  ... ,  index_n ]
         dict ['rxns']      = [  'A+B=C+D' ,  ... , 'D+E=F+G' ] 
         dict ['rct_cmpds'] = [ ['A','B']  ,  ... , ['D,'E']  ] 
                                                     
         
   (OUTPUT) 'index_dict' : (useful for getting all info about a single rxn) 
   
        - KEYS= integer INDEX, VALUES= Sub-Dictionary. 
        - Sub-dictionaries have same string CATEGORY KEYS as input dict, with
          that point to individual elements from the lists in the input dict.
        - The index key mirrors the position of elements in the intput dict lists
        
        #                    OUTPUT DICT. STRUCTURE 
        # list index key = dict{ category keys:  values in list at index_key}
        dict [0]   = {'rxns: 'A+B=C+D' , 'rct_cmpds': ['A','B'] } 
        dict [...] = {'rxns:    ...    , 'rct_cmpds': ...       } 
        dict [n]   = {'rxns: 'D+E=F+G' , 'rct_cmpds': ['D','E'] } 
    """
    # initialize output dictionary
    index_dict = {}
    
    # Figure out how many indvidual records there are (e.g. len of lists) 
    list_len = len(next(iter(category_dict.values())))
    
    # Loop over each index in all lists: 
    for indx in range(list_len):
        # Inialize the sub-dictionary: 
        index_dict[indx] = {}
        # Loop over each categocial key: 
        for cat_key in list(category_dict.keys()):
            # Assign value in output dict for this cateogry at this ind. 
            index_dict[indx][cat_key] = category_dict[cat_key][indx]
    
    return index_dict


def _convert_to_category_dict(index_dict):
    """
     Convert a record based index dictionary into a categorical dictionary .
          
    (INPUT) 'index_dict' : (useful for getting all info about a single rxn) 
    
         - KEYS= integer INDEX, VALUES= Sub-Dictionary. 
         - Sub-dictionaries have same string CATEGORY KEYS as input dict, with
           that point to individual elements from the lists in the input dict.
         - The index key mirrors the position of elements in the intput dict lists
         
         #                    INPUT DICT. STRUCTURE 
         # list index key = dict{ category keys:  values in list at index_key}
         dict [0]   = {'rxns: 'A+B=C+D' , 'rct_cmpds': ['A','B'] } 
         dict [...] = {'rxns:    ...    , 'rct_cmpds': ...       } 
         dict [n]   = {'rxns: 'D+E=F+G' , 'rct_cmpds': ['D','E'] } 
                       
    (OUTPUT) 'category_dict':  (Useful for grouping info about all rxns by category)
        
        - KEYS = string CATEGORIES, VALUES = a LIST of items.
        - All lists under each key have equal lengths.
        - Elements at a particular index in each list corresponds to 
          a matching record in any other list from a dif cateogory key. 
         
          #               OUTPUT DICT. STRUCTURE: 
          #   category key   = list[index_0 ,  ... ,  index_n ]
          dict ['rxns']      = [  'A+B=C+D' ,  ... , 'D+E=F+G' ] 
          dict ['rct_cmpds'] = [ ['A','B']  ,  ... , ['D,'E']  ] 
     """
    # Get all category keys: 
    all_cat_keys = next(iter(index_dict.values())).keys()
    
    # Inialize output dict with category keys pointing to empty lists: 
    category_dict = {cat_key: [] for cat_key in all_cat_keys}

    # Loop over each index (sorted so stuff is all in same order!!!) 
    for indx in sorted(index_dict.keys()):
        # At this index, loop over all category keys 
        for cat_key in all_cat_keys:
            # Append the contents for this category to the appropriate list. 
            category_dict[cat_key].append(index_dict[indx][cat_key])
    
    return category_dict


###############################################################################
# GENERAL UTILITY FUNCTIONS used in make_GCAllRxns.py & make_import_GC_rates.py
###############################################################################
 
def _edit_line(line:str, placeholder:str, val2insert:str): 
 
    nline = line.replace(placeholder, val2insert) if placeholder in line.strip() else line 
        
    return nline 


def str_multi_replace(string:str, baddies, rep_all:str=''): 
    """ Replace multiple items in a string. 
    
    INPUTS:
    -------
        (1) string  -  String you want to replace items in. 
        
        (2) baddies - LIST of chars to replace with rep_all 
                                    --OR--
                      DICT where keys are items to replace and values are replacements. 
                 
        (3) rep_all - (OPTIONAL) String to replace items in list baddies with. Not relevant if baddies is dict. 
        
    REQUIREMENTS: 
    -------------
        LIBRARIES:           import numpy as np 
                             import sys 
        CUSTOM FUNCTIONS:    None
                    
    AUTHOR: 
    -------
    Prof. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD
    
    CHANGE_LOG: 
    ----------
    1/14/2022    JDH Created for pyMCM 
    9/28/2023    JDH Copy/Pasted to GEOSChem_Emulator, updated documentation. 
    """
    # Check user input type of "baddies" to make sure its a list or dict! 
    if all([type(baddies)!=type(ex) for ex in [list([]),np.zeros(1), dict({})]]): 
        raise TypeError("Input to _str_multi_replace() for 'baddies' must be either: \n"+
                        "    (1) a LIST or np.ARRAY of strings you'd like to replace \n"+
                        "                        OR \n" + 
                        "    (2) a DICT where keys are strings to replace and values are replacements.").with_traceback(sys.exc_info()[2])
    
    if type(baddies)!=type(dict({})):  # If baddies is a list: 
        for item in baddies: 
            string=string.replace(item, rep_all)
    else:  # If baddies is a dict! 
        for key in list(baddies.keys()):
            string=string.replace(key, baddies[key])
        
    return string 


def enforce_list_inds(inds,list_ofs):
    """Return items in a list from a list of their indices.
    Mostly used when taking indices from unq2List(rxn_list) and 
    applying them in K,GStr, F lists.
    
    Author: 
    -------
        Dr. Jessica D. Haskins    GitHub: @jhaskinsPhD  Email: jhaskins@alum.mit.edu
    
    Change Log: 
    ----------
    1/14/2022    JDH Created 
    """ 
    
    out=[];# Empty list to hold things for output. 
    
    # Loop over all items turn into arrays so we can index them and append to out
    for item in list(list_ofs): 
        this_list=[]; 
        this_list= [item[i] for i in list(inds) if item[i] not in this_list]
        out.append(this_list)

    return out


def flatten_nested_list(list_in, drop_dupes=False): 
    """Function to take a nested list and flatten into a 1d list.
    
    Author: 
    -------
        Dr. Jessica D. Haskins    GitHub: @jhaskinsPhD  Email: jhaskins@alum.mit.edu
    
    Change Log: 
    ----------
    1/14/2022    JDH Created 
    """ 
    if len(list_in) > 0: 
        is_nested=True if type(list_in[0]) ==type(['v','d']) else False 
        if is_nested is True: 
            list_out=[]; 
            if drop_dupes is True: 
                [list_out.append(r) for list_i in list_in for r in list_i if r not in list_out]
            else: 
                [list_out.append(r) for list_i in list_in for r in list_i]
        else: 
            list_out=list_in 
    else: 
        list_out=[]
            
    return list_out


def drop_dupes(listy): 
    """Function to return a list with no dupes.
    
    Author: 
    -------
        Dr. Jessica D. Haskins    GitHub: @jhaskinsPhD  Email: jhaskins@alum.mit.edu
    
    Change Log: 
    ----------
    1/14/2022    JDH Created 
    
    """ 
    dupes=list_dupes(listy, thresh=0) # List of where eveything is. 
    
    list_out=[]
    for unq in dupes.keys(): 
        first= dupes[unq][0]
        list_out.append(listy[first])
    
    return list_out


def list_dupes(seq, thresh=1): 
    """Finds the location of ALL DUPLICATES in a list in a single pass.
    returns a dictionary with keys= all duplicate values and values = index in seq.
    
    Author: 
    -------
        Dr. Jessica D. Haskins    GitHub: @jhaskinsPhD  Email: jhaskins@alum.mit.edu
    
    Change Log: 
    ----------
    1/14/2022    JDH Created 
    """
    out=dict({}); tally = defaultdict(list)
    
    # Nothing is a list. 
    if all([ type(item)!= type(list()) for item in seq]) : 
        [tally[item].append(i) for i,item in enumerate(seq)]
    else: 
        # If the values in seq are a list, then join them together with a comma 
        # because dict keys can't be lists they must be strings.
        for i,item in enumerate(seq):
            skip=False
            if type(item)==type(list()) and len(item) > 1: 
                item=','.join([str(it) for it in item])
            elif type(item)==type(list()) and len(item)==1: item=item[0] 
            elif type(item)==type(list()) and len(item)<1: skip=True 
            if skip is False: tally[item].append(i)
    
    out= {key:locs for key,locs in tally.items() if len(locs)>thresh}

    return  out


def find_in_list(val2find, inlist, replace_with='', drop_all=False, drop_first=False,
                 match_i=None, match_if_contains=False, partial_match=False): 
    """ Function to find a value within a list, and return all indcies of that match 
    Options include returning a list with all matches to val2find dropped,
    or a list with just the first occurance dropped, or a list with that value replaced with 
    replace_with. 
    
    If you're looking for a match within a list (e.g. if val2find is a list!),
    then use match_i to set the index of where the val2find should be located in the inlist.
    
    Match if contains will return a match to any item in inlist which contains a partial match to val2find.
    
    Author: 
    -------
        Dr. Jessica D. Haskins    GitHub: @jhaskinsPhD  Email: jhaskins@alum.mit.edu
    
    Change Log: 
    ----------
    1/14/2022    JDH Created 
    
    """
    outlist=inlist; indx=[];  # Set the output to the input list. Pop values later if found. 

    # Returns dictionary of all items in inlist where keys are the items and values are the index! 
    inds2all= list_dupes(inlist, thresh=0) 
    
    # If you're trying to find a list within a list!... 
    # then join the target list in the same way list_dupes(in_list).. with commas.
    if type(val2find)==type(list()): 
        if len(val2find) >1: # only do it if you list has more than one item in it. 
                val2find=','.join([str(it) for it in val2find])
        else: #Otherwise find the one value that's in that list! 
            if len(val2find)==1: val2find=val2find[0] 
    
    # If the user wants an exact match and the key of it is in the inds2all. then get its index. 
    if val2find in inds2all.keys() and match_if_contains is False: 
            indx=inds2all[val2find] # List of indices of match to this val
    elif partial_match is True: 
         indx=[inds2all[ky] for ky in list(inds2all.keys()) if val2find in ky] # See if any of the keys contain your val. 
         if len(indx) > 0:  indx=flatten_nested_list(indx)
        
    if match_if_contains is True: 
        indx=[inds2all[key] for key in inds2all if val2find in key.split(',')]
        if len(indx) >0:  indx=flatten_nested_list(indx, drop_dupes=True)
    
    # Replace val2find with replace_with... 
    if replace_with != '' and  (len(indx)>0):
        for ind in indx: outlist[ind]=replace_with 
    
    # Remove ALL occurances of the match from the input list. 
    if drop_all is True and  (len(indx)>0):
        indx.sort(reverse=True) # Sort in decending order (biggest first!!!) 
        if (match_i is None):
            [outlist.pop(i) for i in indx ]
        else: 
            [outlist.pop(i[match_i]) for i in indx]
        
    # Remove only the first occurance of the match from the input list. 
    # or the match at match_i from the input list... 
    if (drop_first is True) and (len(indx)>0):
        indx.sort(); # Either drop first or put match_i in as the "first" to drop... 
        if match_i is not None: indx.pop(match_i); indx.insert(match_i)
        outlist.pop(indx[0])
                
    if ((drop_all is True) or (drop_first is True)) and (len(outlist) == len(inlist)): 
        sys.exit("ERROR in find_in_list(): Couldn't find val2find in inlist.")
        
    return indx, outlist


def join_list_for_MATLAB(join_str:str, ls0:list,min_len:int =75, add_semicolon:bool =False, 
                         preface:str ='', comment:bool=False, comment_skip:list=list([]),
                         insert:str='',insert_skip:list=list([]),adj_ln1_width:int=0): 
    """Function to take a list and join it together with a delimiter as str.join(ls) does 
    except, you actually insert a MATLAB style line break where the lines start to get 
    long. Doesn't append delimiter to beginning or end.
    
    INPUTS: 
    -------
        (1) join_str         - STR delimiter you want to use to join the list: ',', '+'  , etc.
        
        (2) ls               - LIST of items you want to join with 'join_str'. 
        
        (3) min_len          - INTEGER that is the minimum line length a string must be before inserting a line break. 
        
        (4) add_semicolon    - BOOL indicating  wehther to add a semicolon at the end of the string or not. 
        
        (5) comment          - BOOL indicating whether the entire list should be commented our or not. 
        
            (6) comment_skip - LIST of integers corresponding to line #s you DON'T want to comment out 
                               NOTE: This arguement is only relevent if COMMENT=TRUE.
            
        (7) insert           - STR of anything you'd like to insert before every new line (e.g. tab spaceing?)
        
            (8) insert_skip  - LIST of integers corresponding to line #s you DON'T want to insert things at. 
                               NOTE: This arguement is only relevent if INSERT=TRUE.
        
    OUTPUT: 
    -------  
        (1) line - A single STR that has all items in list joined with delimiters and 
                    MATLAB style line breaks where lines are long. 
    REQUIREMENTS: 
    -------------
        LIBRARIES:           None
            
        CUSTOM FUNCTIONS:    None

    AUTHOR: 
    -------
        Dr. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD
                                 
    """
    # Parse inputs options: 
    end_char=';' if add_semicolon is True  else '' 
        
    # Make a copy of the input list to join so you're not changing the original! 
    ls2join=ls0.copy()
    
    # Initialize vars: 
    ln_ls=list([]); lines=list([]); ln_num=1; 
      
    if len(preface) > 0: 
        # Add insert before preface if asked: 
        if 1 not in insert_skip: 
            preface=insert+preface
            
        # Ensure preface is a comment if its supposed to be: 
        if comment is True and 1 not in comment_skip and preface.startswith('%') is False: 
            preface='%'+preface
            
        # Now that we've added the adjustments to the preface as needed, figure out 
        # what the new first line len should be so it doesn't spill over min_len 
        # when the length of the preface is included: 
        ln1_len= min_len-len(preface)
        ln1_has_insert=True
    elif adj_ln1_width>0: 
        ln1_len= np.max([1,min_len-adj_ln1_width])
        ln1_has_insert=False
    else: 
        ln1_len=min_len
        ln1_has_insert=False

    # Loop over all items in the input list to join: 
    for i, item in enumerate(ls2join): 
        
        # Make sure the item is a string... Join only takes string lists! 
        if type(item) != str: item=str(item) 
        
        # Don't add a blank to the list of stuff to join...  
        if len(item) > 0: 
            
            if len(ln_ls)==0: # Indicates we're on a new line: 
                
                # Check to see if the user wants to insert anything before this new line: 
                if insert != '' and ln_num not in insert_skip: 
                    # Don't add insert on line 1 if the preface already has it... 
                    if (ln_num>1) or ((ln_num==1) and (ln1_has_insert is False)): 
                        item=insert+item
                    
                # Make it a comment if they wanted:         
                if comment is True and ln_num not in comment_skip: 
                    # Don't add comment on line 1 if the preface already has it...
                    if (ln_num>1) or ((ln_num==1) and (ln1_has_insert is False)): 
                        if item.startswith('%') is False: item='%'+item
                              
            # Add the (not empty) string to the list of things to join on THIS line.
            ln_ls.append(item) 
            
            # Join everything in ln_ls with this new addition together with "join_str" 
            ln= join_str.join(ln_ls) 
                        
            # If the length of the new string is greater than your min length, you need a new line! 
            # or if you're on the very LAST item in the list, then you need to add it to the list of stuff to write. 
            if (len(ln) > min_len) or (i== len(ls2join)-1) or ((ln_num==1) and (len(ln) > ln1_len)):
        
                lines.append(ln) # Append this full line to the list of lines to join at the end! 
               
                ln_num= ln_num+1 # Update the line number to be the NEXT line... 
                ln_ls=[];        # Reset list holding items to put on a single line 
                
    # Join each line together with the MATLAB line break and new line character 
    to_join= join_str+'...\n' 
    out_str=preface+to_join.join(lines)+end_char
    
    return out_str


def _parse_filename(filename, quiet:bool=True): 
    """Function to take some input string of a filename and seperate it into its 
    path, extension, and actual name. Helper function for check_filename().
    
    INPUTS: 
    -------
       (1) filename - STR with full path and filename  
       
       (2) quiet    - BOOL of whether or not to print output. Default is True. 
    
    OUTPUTS: 
    --------
       (1) opath - STR containing ONLY the path to the file  
                   NOTE: This outputdoes not have a '/' after it, so use os.path.join() 
                         to combine the path with filename! 
                         
       (2) oname - STR containing ONLY the filename
       
       (3) oext - STR contaiing ONLY the file's extension. (e.g. '.m', '.txt', ''.py')
        
    REQUIREMENTS: 
    -------------
        LIBRARIES:           import os
        CUSTOM FUNCTIONS:    None
    
    AUTHOR: 
    -------
        Dr. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD
    
    CHANGE LOG: 
    ----------
    10/03/2023    JDH Created 
     
    """
    # Split filename into name/ extension: 
    [path_and_name,oext] = os.path.splitext(filename)
    
    # Split the file minus its extnsion (if it has one) into the savepath & name 
    opath, oname= os.path.split(path_and_name)
    
    # Print intput/ output if quiet is False: 
    if quiet is False: 
        print("Original Filename = '{}'".format(filename))
        print("    path='{}'\n    name='{}'\n    ext='{}'".format(opath,oname,oext))
        
    return opath, oname, oext
    

def check_filename(filename:str='', default_name:str= '', ext:str='', 
                   savepath:str='', overwrite: bool=False, return_full:bool=False,
                   quiet: bool =True): 
    """ Highly flexible function to check if a file path is valid, if that file 
    exists already, with options to rename the file following some naming 
    convention with an updated version # if you don't want to overwrite the existing file and 
    an option to either return the path/name.extension as a full string or as separate outputs. 
    
    INPUTS: 
    -------
    
            (1) filename     - STR with filename. Valid forms for this can inlude path/extension or not: 
                            
                '/path/to/my_file.txt' - name with path and extension  
                'my_file.txt'          -  name with extension 
                'my_file'              - name without path or extension
            
            (2) default_name - (OPTIONAL) STR with a "default" filename. Accepted forms for this can 
                                inlude path/ extension or not, as above. 
               
            (3)  ext         - (OPTIONAL) STR of the file's extension. NOTE: If this arg is 
                                not passed, either 'filename ' or 'default_name' must contain 
                                the extension of the file. 
                            
            (4) savepath    - (OPTIONAL) STR with path of the file.NOTE: If this arg is 
                                not passed, either 'filename ' or 'default_name' must contain 
                                the path of the file. 
                                
            (5) overwrite    - (OPTIONAL) BOOLEAN indicating whether you'd like  to overwrite 
                                the file if it exists. Default is FALSE which will create a new 
                                file with a version # string attached to it like my_file_v3.txt
            
            (6) return_full  - (OPTIONAL) BOOLEAN  indicating whether you'd like to return
                                a STR with the absolute path to the file+ext (if True) or 
                                if you'd like to return a list containing the [path, file+ext]. 
                                Default is to return a list seperating path from file... 
                                
            (7) quiet        - (OPTIONAL) BOOLEAN indicating whether or not to print 
                                "pass" + filename or not. Default is set to FALSE (suppress "pass" output).
    
    OUTPUT: 
    ------- 
        EITHER:                                
               [savepath, filename+ext] - LIST with (updated) strings if return_full=FALSE (DEFAULT) 
        OR: 
                savepath+filename+ext   - STR of absolute path to file if return_full=TRUE
    REQUIREMENTS: 
    -------------
        LIBRARIES:           import os 
                             import sys 
                             import re 
                             import inspect 
            
        CUSTOM FUNCTIONS:    _parse_filename()     (defined above).       
        
    AUTHOR: 
    -------
        Prof. Jessica D. Haskins    GitHub: @jhaskinsPhD  Email: jessica.haskins@utah.edu
    
    CHANGE LOG: 
    ----------
    1/18/2022    JDH Created 
    7/19/2022    JDH add "quiet" option. 
    10/12/2023   JDH added user input parsing errors, updated version # creation using regex, 
                    & updated default "savepath" to be set to directory of the script calling this function
                    + drastically updated documnetation. Now depends on regex & inspect libraries. 
    
    """
    ###############################################################################################
    # (A) Parse user inputs for erros & figure out what combo of name/path/extension options were used
    ###############################################################################################
    using_default=False # Set boolean to indicate if we're using the default name (vs the filename)... initally False!
    opath=None
    
    # (A1) Determine if the user passed input to savepath or not & set opath if they did to 'savepath' 
    if len(savepath) ==0:
        opath_defined=False 
    else:
        opath=savepath; opath_defined=True
            
    # (A2) If the user didn't pass "filename", but did pass a default_name, set filename to default_name:  
    if len(filename)==0 and len(default_name) >0: 
        filename=default_name; using_default=True; 
    
    # (A3) If they didn't specify either a filename or a default name, throw an error
    elif len(filename)==0 and len(default_name)==0: 
        raise ValueError("INPUT ERROR (A.3) in check_filename(): \n"+
                         "   You must specify either 'filename' or a 'default_name' as input arguements.").with_traceback(sys.exc_info()[2])
    
    # Create string to define whether you're using 'default_name' or 'filename' as input... 
    inp_type='default_filename'if using_default is True else 'filename' 
    
    # Parse the filename to see if it has a path/name/extension attached to it! 
    ipath, name_only, iext = _parse_filename(filename)
    
    # Set boolean to indicate if the filename contains a path or not.. 
    path_in_filename=True  if len(ipath)>0 else False 
    
    # (A.4) Check parsed extension, make sure it matches input to "ext", and that some extension is defined! 
    if len(ext) > 0: 
        # (A.4.1) If user passed input for ext, remove spaces and make sure it has the 
        #     period in front of it in the "ext" var for the file's extension...
        ext='.' +ext.strip() if (('.' not in ext) or (ext.strip()[0]!='.')) else ext.strip()
        
        # (A.4.2) Throw an error if the parsed extension doesn't match the input extension! 
        if len(iext) > 0 and ext != iext: 
            raise ValueError("INPUT ERROR (A.4.2) in check_filename():\n"+
                             "    You passed inconsistent extensions in the input for '"+inp_type+"',  and 'ext'.\n"+
                             "        '"+inp_type+"' = '"+filename+"\n"+
                             "        'ext' = '"+ext+"'\n"
                             "    Check_filename() thinks: '"+iext +"' != '"+ext+"'.").with_traceback(sys.exc_info()[2])
    
    # (A.4.3)If no extension was passed, used iext as the extension!
    if len(ext)==0 and len(iext) >0: ext=iext
    
    # (A.4.4) Throw an error if no file extension can be found: 
    elif len(ext)==0 and len(iext)==0: 
        raise ValueError("INPUT ERROR (A.4.4) in check_filename():\n"+
                         "    You must include the extension of the file you wish to 'check' for one of the following input arguements: \n"+
                         "       (1)'filename'\n"+
                         "       (2) 'default_name'\n"+
                         "       (3) 'ext'").with_traceback(sys.exc_info()[2])
    
    ###########################################################################
    # (B) Figure out the file's path, name, and extension! 
    ###########################################################################
    
    # (B.1) If no savepath is defined and there IS a path in the filename, use that as opath... 
    if opath_defined is False and path_in_filename is True: 
        opath=ipath; opath_defined=True        
    
    # (B.2) If there is NOT a path in the filename/default_name set opath to the path 
    #       of the file that called this function. 
    if opath_defined is False and path_in_filename is False: 
            
        # Get the filename of the script that called this function: 
        parent_filename=inspect.stack()[1].filename
                
        # Set opath to the path of the parent script & tell the rest of the code we found opath! 
        opath, dum= os.path.split(parent_filename); opath_defined=True
                
        # Print warning that we arbitrarily set the savepath...
        raise Warning("No 'savepath' could be found in the input args. It will be set it to the following path,\n"+
                      "    corresponding to that of the script that called check_filename(): \n\n"+
                      "savepath = '"+opath+"'")
    
    if opath_defined is False: 
        raise ValueError('DEBUGGING ISSUE: The path to the file in check_filename() was not ever assigned... ')

    # Check that the output path actually exists! 
    if not os.path.exists(opath): 
        raise NotADirectoryError("In check_filename(), the following parsed directory of the file to 'check' does not exist: \n\n"+ 
                                 "   path_checked ='"+opath+"'\n\n"+
                                 "Your inputs this was parsed from were: \n"
                                 "   "+inp_type+" = '"+filename+"'\n"+ 
                                 "   savepath = '"+savepath+"'\n\n"+
                                 "If the path_checked looks wrong compared to these, its an issue within check_filename()...").with_traceback(sys.exc_info()[2])
                                 
    # Define the absolute path of the file you're checking... 
    fullname=os.path.join(opath,name_only+ext)
    
    ###########################################################################
    # (C) Do file overwriting / renaming if asked... 
    ###########################################################################
    n=1; #initialize counter variable! 
    
    # If/ While what we have set "fullname" to exits already ...  
    while os.path.isfile(fullname):

        ##########  C.1 Delete file if it exists & "overwrite" is True ########
        if overwrite is True:
            os.remove(fullname); just_file=name_only+ext
            
        ##########  C.2 Otherwise append a version # to the file ##############
        else:
            # Use a regex string to pull out the last version # of file
            # This regex looks for 1+ # following "_v" ... Pulls out ['39'] from 'my_file_v39.txt'
            last_ver=re.findall(r'j([a-zA-Z0-9]+)',name_only)
            
            # C.2.1 If the file exists but does NOT have a version #, then add one! 
            if len(last_ver) == 0: 
                just_file=name_only+'_v'+str(n)+ext
            
            # C.2.2 If the file exists and DOES have a version number already... 
            else: 
                # Update version # to be +1 and convert to a string. 
                new_ver='_v'+str(np.int32(last_ver[0]) +n)
                                 
                # Then get the old name, without the version #, and update "just_file"
                just_file=name_only.split('_v'+last_ver[0])[0]+new_ver
            
            # After changing the filename accordingly, update "fullname" to check if it exists! 
            fullname=os.path.join(opath,just_file)
            
        n=n+1 # Update counter variable! 
            
        # You now havea valid (unique) file path, name, extension combo!
        if os.path.isfile(fullname) is False: break 
    
    ###########################################################################
    # (D) Prep Output
    ###########################################################################
    
    # Print "PASS" info once we get here if asked... 
    if quiet is False: print('Check_filename() found no errors for:'+ fullname)
    
    # Decide what to return to the user... 
    if return_full is True: 
        return fullname
    
    else: # Have to add '/' after opath because pullit it out with os doesn't keep it & 
          # os.path.join() adds it for us when we make fullname...  
        return [opath+'/', just_file]


def get_my_relative_paths(): 
    """Function used to identify where you've installed the GEOS-Chem Emulator on 
    your computer so we can tell it how to find all the relative inputs it needs!""" 
    
    # Directory where the GEOS-Chem Emulator is stored: 
    GCEM_path=os.path.dirname(os.path.abspath(__file__))
        
    # Where we expect to find Input Files: 
    input_path=os.path.join(GCEM_path, 'Input_Files')
        
    # Where the dif photolysis FJX files are located: 
    photolysis_paths=dict({'FJX_Files':                     os.path.join(input_path, 'Photolysis_Files/FJX_Files'),
                           'FJX_input_by_GC-version.xlsx':  os.path.join(input_path, 'Photolysis_Files/FJX_input_by_GC-version.xlsx'),
                           'FJX_cross_sect_to_jvalues.xlsx':os.path.join(input_path, 'Photolysis_Files/FJX_cross_sect_to_jvalues.xlsx'),
                           'jmap_type_help': os.path.join(input_path, 'Photolysis_Files/jmap_type_README.txt')
                           })
    
    # Where the templates for writing the F0AM files are located: 
    template_paths=dict({'GasRxns':      os.path.join(input_path, 'Template_Files/GEOSChem_GasRxns_template.txt'), 
                         'HetRxns':      os.path.join(input_path, 'Template_Files/GEOSChem_HetRxns_template.txt'), 
                         'K':      os.path.join(input_path, 'Template_Files/GEOSChem_K_template.txt'), 
                         'J':      os.path.join(input_path, 'Template_Files/GEOSChem_J_template.txt'), 
                         'import_GC_rates': os.path.join(input_path, 'Template_Files/import_GC_rates_template.txt')
                         }) 
            
    # Dump into ouput dictionary for all: 
    my_paths= dict({'Photolysis_Files': photolysis_paths,
                    'Template_Files': template_paths}) 
        
    return my_paths





                    


    



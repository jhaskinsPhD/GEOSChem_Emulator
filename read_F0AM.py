#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  7 14:25:05 2025

@author: u6044586
"""
import sys 
import csv 
import numpy as np 

sys.path.append('/uufs/chpc.utah.edu/common/home/u6044586/python_scripts/modules/GEOSChem_Emulator/')
from foam_rxn_class import * 

def check_for_duplicates(inp_list):
    unique = set(); duplicates = set()
    
    for item in inp_list:
        if item in unique:
            duplicates.add(item)
        else:
            unique.add(item)

    return list(unique), list(duplicates)

def compare_lists(dict1:dict, dict2:dict ):
    # Unpack inputs from dictionarys. Keys= list names, lists= items. Expects 1 key only. 
    list1_name, list1= dict1.items();
    list2_name, list2= dict2.items()
    
    # Turn lists into sets for easy comparisons: 
    set1 = set(list1); set2 = set(list2)

    # Figure out diffs in the sets: 
    unq_to_list1 = set1 - set2
    unq_to_list2 = set2 - set1
    common_items = set1 & set2

    # Put contents into output dict: 
    out_dict={'Unq_to_{list1_name}': in_list1_not_in_list2,
              'Unq_to_{list2_name}': in_list2_not_in_list1,
              'In_Both':common_items}
    
    return out_dict


def is_MATLAB_blank(line) :
    # Assum line is not blank: 
    is_blank=False 
    
    # Strip all white space from the line & replace the new line character... 
    line_stripped=line.replace(' ', '').replace('\n','') 
    
    # If its now empy or only has a '%' comment char on it , its blank. 
    if len(line_stripped)==0 or len(line_stripped.replace('%','')) ==0: 
        is_blank=True
    
    return is_blank 
        
def is_MATLAB_comment(line) :
    # Assume line is not a comment: 
    is_comment=False 
    
    # Strip all white space from the line & replace the new line char... 
    line_stripped=line.replace(' ', '').replace('\n','') 
    
    if len(line_stripped)>0 and line_stripped[0] =='%':
        # If it begins with a '%' char in matlab, then it is a comment: 
        is_comment=True
    
    return is_comment 


def read_MATLAB_cell_list(file, line, count,drop_dupes:bool=True): 

    # Make dict where line number is key, and value is line parsed when making list: 
    orig_list_lines={count:line} 
    
    # Split the line on the char marking the beggining of a MATLAB cell list: 
    ln_parts=line.split('{')
    
    ########################################################################################
    # Check to make sure the first line containing the cell list is formatted how you think:
    ########################################################################################
    if len(ln_parts)==1: 
        # If line doesn't actually contain a '{' char indicating start of list thow error. 
        raise ValueError("The following line passed as input to read_MATLAB_cell_list() "+\
                         "does not appear to contain the character '{' expected to "+\
                         "indicate that a list is beginning on this line:"+\
                         f"\t {line}").with_traceback(sys.exc_info()[2])
    elif len(ln_parts) !=2: 
        # If line contains MORE THAN 1 '{' character, throw a warning saying idk whht to do: 
        raise Warning("The following line passed as input to read_MATLAB_cell_list() "+\
                      "contains more than 1 '{' character. Not sure what "+\
                      "to do with the rest of it. "+\
                       f"\t {line}").with_traceback(sys.exc_info()[2])
        # Re-join subsequent items with the '{' delimiter assuming its part of the item:  
        list_str= '{'.join(ln_parts[1:])
                         
    else: 
        # Just pull out the part of the line that has the list items in it 
        list_str=ln_parts[1]

    
    # ONLY keep items that come BEFORE a comment character if it appears on the line:
    if '%' in list_str: 
        list_str=list_str.split('%')[0]
        
    # Remove leading/trailing spaces or new line characters: 
    list_str=list_str.strip().replace('\n','')
        
    # Read more lines til you get to the line at the end of the list (does 
    # nothing if the list finishes on this line!). 
    while list_str.endswith('...') == True or '}' not in list_str:
        
        nline = file.readline()   # Read next line from file
        count = count+1 # Update line counter variable
        
        orig_list_lines[count]=nline # Keep line in dict of orig lines parsed: 
        
        # ONLY keep items that come BEFORE a comment character if it appears on the line: 
        if '%' in nline: 
            nline=nline.split('%')[0]
                    
        # Remove any leading / trailing spaces or new line chars from the line:  
        nline=nline.strip().replace('\n','')
    
        # Append this new line to our current string with the full list: 
        list_str = list_str + nline.strip()
    
     
    # Now that you have the full list in a single string, remove all the
    # "'" for strs "..." for line continues and '}' for close of list chars: 
    list_str= list_str.strip().replace('\n','').replace('...','').replace("'",'').replace('}', '')
    
    # Now split the list on the ';' or ',' delimiters 
    if ';' in list_str: 
        my_list=list_str.split(';') 
    elif ',' in list_str: 
        my_list=list_str.split(',') 
        
    # Remove any leading/trailing spaces around the items in the list (& empties)
    my_list=[item.strip() for item in my_list if len(item.strip())>0] 
    
    # Make dict with line number start/stop containing original lines: 
    
    # Return the parsed list, new line number, dict with line numbers/orginal list of list
    return my_list, count, orig_list_lines
    










def read_foam_mech(filepath): 
    #------------------------------------------------------------------------------
    # Initialize Variables Used in While Loop:  
    #------------------------------------------------------------------------------
    # Get total number of of lines so we know when we've reached the end of the file
    num_lines = sum(1 for line in open(filepath, 'r'))
    
    # Initalize line number and reaction number indices: 
    count = 0;  rxn_indx=-1 
        
    # Initalize output dicts that will hold all info on mechanism:
    all_info={};mech_info={}; 
    
    # Initalize lists that will hold commented out lines for specific parts of mechanism: 
    header_lines=[]; after_sp_list=[]; after_ro2_list=[];before_first_rxn=[]; 
    
    # Set bools that tells us what part of the file we're currently parseing: 
    file_loc={'in_header':True,
              'past_sp_list':False,
              'past_ro2_list':False,
              'past_addspecies':False} 
        
    # %%   
    # Open the file & start parsing:    
    file = open(filepath, 'r')
    
    while True:  # Loop over each line of the input F0AM mechanism file: 

        line0 = file.readline() # Read next line from file
        line=line0              # Keep original line, & make copy we're gonna manipulate. 
        count=count+1          # Adjust line counting variabl
        
        # Figure out if the line we're parsing is a commented out line or is a blank line. 
        is_cmt= is_MATLAB_comment(line)   
        
        ###########################################################################
        # (1) Parsing F0AM lines betwee:      HEADER  & SPECIES LIST
        ###########################################################################
        if file_loc['past_sp_list'] is False: 
                
            if is_cmt==True: 
                #------------------------------------------------------------------
                # Parse the Header
                #------------------------------------------------------------------
                # Pull out declared # of species in the mechanism frpm the header : 
                if all([word in line for word in ['#', 'of', 'species', '=']]):
                    number_of_what, equals, this_value = line.rpartition("=")
                    all_info['nSp_in_header']={'value':int(this_value.strip()),
                                               'orig_line': line}                
                    
                # Pull out declared # of reactions in the mechanism from the header : 
                elif all([word in line for word in ['#', 'of', 'reactions', '=']]):   
                    number_of_what, equals, this_value = line.rpartition("=")
                    all_info['nRxns_in_header']={'value':int(this_value.strip()),
                                                 'orig_line': line}
                else: 
                    # Store all comments that come before the species list here:
                    header_lines.append(line) 
                    
            else: # if its not a comment & we're not past the species list yet: 
                #------------------------------------------------------------------
                # Parse the Species list: 
                #------------------------------------------------------------------
                if 'SpeciesToAdd' in line: 
    
                    # Read in all next lines til end of species list is found: 
                    sp_list, count, sp_list_lines= read_MATLAB_cell_list(file, line, count)
             
                    # Update bool that tells us about our file locaiton... 
                    file_loc['past_sp_list']=True
                    
        ###########################################################################
        # (2) Parsing F0AM lines between:      SPECIES LIST & RO2 LIST
        ###########################################################################
        elif file_loc['past_sp_list'] is True and file_loc['past_ro2_list'] is False:
            
            if is_cmt==True: 
                # Store comments occuring between the species list &  RO2 lists here: 
                after_sp_list.append(line) 
            else:
                #------------------------------------------------------------------
                # Parse the RO2 list: 
                #------------------------------------------------------------------
                if 'RO2ToAdd' in line:   # Parse the RO2 list: 
                    
                    # Read in all lines til end of the RO2 list is found: 
                    ro2_list, count,ro2_list_lines= read_MATLAB_cell_list(file, line, count)
                            
                    # Update bool that tells us we're past the line containing the RO2 list: 
                    file_loc['past_ro2_list']=True
        
        ###########################################################################
        # (3) Parsing F0AM lines between:       RO2 LIST , 'AddSpecies', &  1st Rxn 
        ###########################################################################
        elif file_loc['past_ro2_list'] is True and file_loc['past_addspecies'] is False:
            # Store comments occuring between the RO2 list & before "Add Speices" here:
            if is_cmt==True:  
                after_ro2_list.append(line)  
            
            elif 'AddSpecies' in line.strip(): 
                # Update bool that tells us we're past the line containing declaring 'AddSpeices' 
                file_loc['past_addspecies']=True
        
        ###########################################################################
        # (4) Parsing F0AM lines:            Containing all reactions
        ###########################################################################
        elif file_loc['past_addspecies'] is True: 
            
            if is_cmt is False: 
    
                # Seperate any in-comments from the uncommented out part of the line. 
                if '%' in line:
                    # Split the line into the non-comment & comment parts: 
                    line_parts=line.split('%') 
                    line=line_parts[0];  inline_cmt='%'.join(line_parts[1:])
                    
                    # If the comment has nothing after replacing spaces/ new line chars, 
                    # then leave it as blank: 
                    if len(inline_cmt.replace(' ','').replace('\n','')) ==0: inline_cmt=''
                else: 
                    inline_cmt='' 
                
                line_stripped=line.replace(' ', '') # Remove any spaces from the line. 
                
                # We'll assume the 'Rnames{i}='... denotes beginning of a new rxn: 
                # and will update the 'rxn_indx' used to store info on the rxn, k_line, 
                # g_line, and f_lines about this reaction in the output mech_info dictionary. 
                if 'Rnames{i}=' in line_stripped: 
                    # Update the rxn_index: 
                    rxn_indx=rxn_indx+1
                    mech_info[rxn_indx]=dict({'rxn_line':line,
                                              'rxn_cmts':inline_cmt})
                elif 'k(:,i)=' in line_stripped:  
                    mech_info[rxn_indx]['k_line']=line
                    mech_info[rxn_indx]['k_cmts']=inline_cmt
    
                elif 'Gstr{i,1}=' in line_stripped:  
                    mech_info[rxn_indx]['g_line']=line
                    mech_info[rxn_indx]['g_cmts']=inline_cmt
                elif '(i)=f' in line_stripped:  
                    mech_info[rxn_indx]['f_line']=line
                    mech_info[rxn_indx]['f_cmts']=inline_cmt
                    
            else: 
                # If it is a commented line in the middle of all this info about rxns:  
                if rxn_indx==-1: 
                    # Store comments occuring after "AddSpecies" and before first 
                    # uncommented 'Rnames{i}=' line here: 
                    before_first_rxn.append(line)
                else:
                    # If its truely a random comment line in this mech, store it in the 
                    # previous rxns dict as a comment to come afterwards: 
                    if 'cmts_after' not in list (mech_info[rxn_indx]): 
                        mech_info[rxn_indx]['cmts_after']= line
                    else: 
                        mech_info[rxn_indx]['cmts_after']=mech_info[rxn_indx]['cmts_after']+'/n'+line
                
        # If you reach the end of the file the break out of the while loop: 
        if count > num_lines:
          break
      
    
      
    
    # Add all the other info the the output: 
    all_info['1_header_lines']=header_lines 
    all_info['2_sp_list']=sp_list
    all_info['3_cmts_after_sp_list']=after_sp_list
    all_info['4_ro2_list']=ro2_list
    all_info['5_cmts_after_ro2_list']=after_ro2_list
    all_info['6_cmts_before_first_rxn']=before_first_rxn
       
        #  Return stuff we parsed: 
    return all_info, mech_info










filepath='/uufs/chpc.utah.edu/common/home/u6044586/MATLAB/F0AM_v4.4.2/Chem/GEOSChem/v14.5.2/GEOSCHEM_GasRxns_v1452.m'

all_info, mech_info= read_foam_mech(filepath) 


# %% 
rxn_list=[]

for k,d in mech_info.items(): 
    rxn_i= d['rxn_line'].replace('Rnames{i}=','').replace("'",'').replace(';','').replace('\n','')
    rate_i=d['k_line'].replace('k(:,i)=','').replace("'",'').replace(';','').replace('\n','')
    
    rxn_obj=reaction(rxn_str=rxn_i, rate_str=rate_i, rxn_arrow='=', rcts2ignore=['hv'])
    
    rxn_list.append(rxn_obj) 
    
my_mech=mechanism(rxn_list=rxn_list)   
    
    




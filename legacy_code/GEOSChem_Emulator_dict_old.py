#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 12:23:52 2023

@AUTHOR: u6044586
"""


import os
import sys
import re
import numpy as np
import pandas as pd
from datetime import date

sys.path.append('/uufs/chpc.utah.edu/common/home/u6044586/python_scripts/modules/GEOSChem_Emulator/')
from read_kpp_mod import read_kpp 
import photolysis as hv 
import utils



# %%---------------------------------------------------------------------------
# -------------- (A.4) make_header_for_F0AM_mech() & Sub-Functions ------------
# -----------------------------------------------------------------------------
# A.4 L1-function called within make_GEOSCHEM_AllRxns_file()
def eliminate_dead_ends(kpp_dict:dict, verbose:bool=True):
    """ID Reactions from stuff that's ONLY formed via a Heterogenous reaction
       Basically cut dead end reactions that we don't need...

       Custom Functions: 
           flatten_nested_list  --> in pyMCM/F0AM_Tools/utils.py
           enforce_list_inds    --> in pyMCM/F0AM_Tools/utils.py
           find_in_list         --> in pyMCM/F0AM_Tools/utils.py
           drop_dupes           --> in pyMCM/F0AM_Tools/utils.py
           """
    # Make a dict with info only on het/not-het rxns: 
    het_info={k:v for k,v in kpp_dict.items() if v['is_het']==True}
    nhet_info={k:v for k,v in kpp_dict.items() if v['is_het']==False}
    
    # Extract (nested) lists of all reactants/products in from any non-het rxns ...
    nhet_inds=[rx_ind for rx_ind,rx_dict in nhet_info.items()] 
    nhet_rxns=[rx_dict['reaction'] for rx_ind,rx_dict in nhet_info.items()] 
    nhet_rct_cmpds=[rx_dict['rct_cmpds'] for rx_ind,rx_dict in nhet_info.items()] 
    nhet_prd_cmpds=[rx_dict['prd_cmpds'] for rx_ind,rx_dict in nhet_info.items()] 

    # Extract a (nested) list of all reactants/products in from any het rxns ...
    het_prd_cmpds=[rx_dict['prd_cmpds'] for rx_ind,rx_dict in het_info.items()] 

    # Start with a a list of all UNIQUE products formed in any het rxns 
    het_prds = utils.flatten_nested_list(het_prd_cmpds, drop_dupes=True) 
    
    #Iinitialize vars controlling the while loop: 
    growing = True
    start = True
    while growing:
        # Initialize everything on first pass through to nested lists of rcts/prds 
        # in the nhet mech & empty lists of stuff to hold what we'll drop on exit. 
        if start == True:  
            rct_a = nhet_rct_cmpds
            prd_a = nhet_prd_cmpds
            inds2drop = []
            sp2drop = []
            start = False
        else:
            # Pretend you dropped all the dead-end reactions you found on last interation 
            # Now, what cmpds are ONLY formed in het rxns are reactants in nhet rxns? 
            # ((Since we dropped some prds from nhet mech, there may be new species 
            # that are now ONLY formed in the het rxns / why we're in a while loop.))
            inds2keep = [indx for indx in range(0, len(nhet_rxns)) if indx not in inds2drop]
            [rct_a, prd_a] = utils.enforce_list_inds(inds2keep, [nhet_rct_cmpds, nhet_prd_cmpds])

        # Get a unique list of products that are still being formed  in the non-het rxns.
        prds = utils.flatten_nested_list(prd_a, drop_dupes=True)

        # Get current size of list of rxn inds2drop before we drop any more on this iteration 
        sz0 = len(inds2drop)  

        # Find all species formed in het rxns that are NOT formed from any non-het rxns
        dead=[]
        [dead.append(sp) for sp in het_prds if (sp not in prds) and (sp not in dead)]

        # Loop over all species ID'd as only forming through a HET rxn.
        for cmpd in dead:
            
            # Find all reactions in the non-het rxns that this species is a REACTANT in.
            inds, args = utils.find_in_list(cmpd, nhet_rct_cmpds, match_if_contains=True)

            # Add the index of this rxn (in non-het mech) to list we'll drop.
            [inds2drop.append(i) for i in inds if i not in inds2drop]

            # Add the species to the list of compounds you know you need to drop.
            if cmpd not in sp2drop:
                sp2drop.append(cmpd)

        sz1 = len(inds2drop) # Now take size of list of rxn inds to drop ... 
        if sz0 == sz1: break  # if list stopped growing, stop while loop!

    # Get a list of the top level key in the kpp_dict for the rxns we will DROP 
    # from the entire mechuanism (those in the nhet mech that depend on het stuff) : 
    KPP_inds2drop = [nhet_inds[indx] for indx in inds2drop]
    rxns2drop=[v['reaction'] for k,v in kpp_dict.items() if k in KPP_inds2drop]
    
    # Print out info on what rxns are being commented out if verbose is True: 
    if verbose is True and len(rxns2drop) > 0:
        print('Commenting out these gas/ photolysis rxns from mech because they ' +
              'involve species only formed in het reactions, which are not included:\n')
        [print('    '+r) for r in rxns2drop]
        print('\n')
    
    return  KPP_inds2drop
# %%---------------------------------------------------------------------------
# -------------- (A.3) KPP_Dict_to_F0AM_IO() & Sub-Functions ------------------
# -----------------------------------------------------------------------------

# A.3.3 L2-function called within KPP_Dict_to_F0AM_IO within make_GEOSCHEM_AllRxns_file()
def _append_to_ndict_list(ddict_in: dict, key1: str, key2: str, values, allow_dupes: bool = False):
    """Function to append a list of values to a key in a 2 level nested dictionary that 
    points to a list. 

    INPUTS: 
    -------
       (1) ddict_in    - (Outer/Full) Nested DICT you want to add values within 

       (2) key1        - STR of the first key to access the value of the inner nested dict

       (3) key2        - STR of the second key to access the value of the nested dict

       (4) values      - LIST of values you want to append to the inner dict. 

       (5) allow_dupes - BOOL indicating whether or not the ouput should be allowed 
                       to contain duplicates or not. 
    OUTPUTS: 
    --------
       (1)  ddict     - Revised input DICT with new values appended to inner nested dict list. 

    REQUIREMENTS: 
    -------------
        LIBRARIES:           None

        CUSTOM FUNCTIONS:    None

    AUTHOR: 
    -------
        Prof. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD

    CHANGE LOG: 
    ----------
        10/12/2023 JDH Created
    """
    # Make a copy of the input so we don't change it.
    ddict = ddict_in.copy()

    # Make sure the input is a list so the fuction works.
    if type(values) != list:
        values = [values]

    # If key 2 doesn't exist, add the values as a list underneath it.
    if key2 not in list(ddict[key1].keys()):
        ddict[key1][key2] = values
    else:
        # Don't worry about dupes, just append the values to the nested dict.
        if allow_dupes is True:
            ddict[key1][key2] = ddict[key1][key2]+values
        else:
            # Loop over values and only add those that aren't in there already.
            for v in values:
                if v not in ddict[key1][key2]:
                    ddict[key1][key2] = ddict[key1][key2]+[v]

    return ddict



# %%---------------------------------------------------------------------------
# -------------- (A.3) Prep/Write to F0AM  & Sub-Functions---------------------
# -----------------------------------------------------------------------------
# A.3.3 L2-function called within KPP_Dict_to_F0AM_IO within make_GEOSCHEM_AllRxns_file()
def _append_to_ndict_list(ddict_in: dict, key1: str, key2: str, values, allow_dupes: bool = False):
    """Function to append a list of values to a key in a 2 level nested dictionary that 
    points to a list. 

    INPUTS: 
    -------
       (1) ddict_in    - (Outer/Full) Nested DICT you want to add values within 

       (2) key1        - STR of the first key to access the value of the inner nested dict

       (3) key2        - STR of the second key to access the value of the nested dict

       (4) values      - LIST of values you want to append to the inner dict. 

       (5) allow_dupes - BOOL indicating whether or not the ouput should be allowed 
                       to contain duplicates or not. 
    OUTPUTS: 
    --------
       (1)  ddict     - Revised input DICT with new values appended to inner nested dict list. 

    REQUIREMENTS: 
    -------------
        LIBRARIES:           None

        CUSTOM FUNCTIONS:    None

    AUTHOR: 
    -------
        Prof. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD

    CHANGE LOG: 
    ----------
        10/12/2023 JDH Created
    """
    # Make a copy of the input so we don't change it.
    ddict = ddict_in.copy()

    # Make sure the input is a list so the fuction works.
    if type(values) != list:
        values = [values]

    # If key 2 doesn't exist, add the values as a list underneath it.
    if key2 not in list(ddict[key1].keys()):
        ddict[key1][key2] = values
    else:
        # Don't worry about dupes, just append the values to the nested dict.
        if allow_dupes is True:
            ddict[key1][key2] = ddict[key1][key2]+values
        else:
            # Loop over values and only add those that aren't in there already.
            for v in values:
                if v not in ddict[key1][key2]:
                    ddict[key1][key2] = ddict[key1][key2]+[v]

    return ddict

def _format_f_line(data, fixed_vars): 
    """Helper function to prepare/format the F-line for F0AM""" 
    #Initialize a dictionary with keys as the (unique) species involved in 
    # the rxn and their values as net stoichimetry of the rxn. Allows us to 
    # combine stoichiometry before writing f-line for reactions where cmpds 
    # appear in both the prds & rcts or twice in prds/rcts individually :  
        
    # (1) CMPD appears multiple times in prd/rcts:     'MO2+MO2=2CH2O +20HO2'
    #     Here, the rxn loses 2 mol MO2, despite it appearing twice... 
    #     The F0AM F-line should be: "fMO2(i)=fMO2(i)-2; ..."
    #      rather than "fMO2(i)=fMO2(i)-1; fMO2(i)=fMO2(i)-1; ..."
    # (2) Has CMPD in RCTS & PRDS:     'R7P+OH=0.500OH+0.500RCHO'
    #     Here, the rxn loses 1 mol OH, but regenerates 0.5 mols of OH.
    #     The F0AM F-line should be:   
    #       fR7P(i)=fR7P(i)-1; fOH(i)=fOH(i)-0.5; fRCHO(i)=fRCHO(i)+0.5;
                                 
    net_dict = {};fstrs=[]
    
    def accumulate_stoich(cmpd:str, stoich, net_dict:dict):
        # Helper function to accumulate stoichiometry for input compound 
        # in dictionary or assign it to 0 if no value exists. 
        if cmpd not in net_dict:
            net_dict[cmpd] = 0.0
        net_dict[cmpd] += stoich
        return net_dict
    
    # Accumulate reactant stoichiometry  
    for cmpd, stoich in zip(data['rct_cmpds'], data['rct_stoich']):
        net_dict=accumulate_stoich(cmpd, -stoich, net_dict)
        
    # Accumulate product stoich 
    for cmpd, stoich in  zip(data['prd_cmpds'], data['prd_stoich']):
        net_dict=accumulate_stoich(cmpd, stoich, net_dict)
        
    # Store the results
    for cmpd, stoich in net_dict.items():
        if stoich!=0: 
            operator='+' if stoich > 0 else '-'
            if cmpd not in fixed_vars: # Don't adjust concs of fixed vars! 
                fstrs.append(f"f{cmpd}(i)=f{cmpd}(i){operator}{np.abs(stoich)};")
    data['F_Line'] =' '.join(fstrs)

    return data 

def _add_met_to_rates(string): 
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
            
            # Add Met to arg list: 
            args_list.append('Met')
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
    
def prep4f0am(kpp_dict, inds2drop, fixed_vars, fjx_df, jmap_dict): 
    not_num=[]; num= []; 
    
    # Create F0AM Dicitonary to hold all output: 
    foam_dict={'nhet':{'lines':{},
                       'req_ks':[],
                       'req_js':[],
                       'sp_list':[],
                       'rate_dict':{}, 
                       'comment_out':[]},
               'het':{'lines':{},
                       'req_ks':[],
                       'req_js':[],
                       'sp_list':[], 
                       'rate_dict':{}}}
    
    # Loop over all reactions we have info about from KPP & convert to lines for F0AM: 
    for indx in list(kpp_dict.keys()):
            
        # Extract dictionary about this specific reaction parsed from KPP File: 
        data= kpp_dict[indx]

        #######################################################################
        # (1) Format Reaction Line:  "Rnames{i}='HNO3+SALAAL=NIT'; % Comment?" 
        #######################################################################
        # Format reaction line for F0AM: 
        data['Rxn_Line']="Rnames{i}='"+data['reaction']+ "';"

        # If the KPP file had a comment about that reaction, add it back in... 
        if len(data['comments'])>0: 
           data['Rxn_Line']=  data['Rxn_Line'] + f" % {data['comments']}" 

        #######################################################################
        # (2) Format K-Line:  "k(:,i) = GCARR_ac(4.80e-11, 250.0, Met);" 
        #######################################################################
        #----------------------------------------------------------------------
        # Ensure any FORTRAN operators in rate functs are in MATLAB syntax
        #----------------------------------------------------------------------
        # NOTE: In standard Fortran, the arithmetic operators like "*"  & 
        # "/" must be immediately adjacent to operands with no intervening 
        # spaces. So,the regex patterns used here are look behind/ahead
        # assertions to ensure the '*' or '/' is not preceded or followed by 
        # whitespace so that we ONLY replace mathematical operators this way  
        
        # Replace division operators '/' with './' for MATLAB
        rate = re.sub(r'(?<=\S)/(?=\S)', './', data['rate']) 
        
        # Replace multiplication operators '*' with '.*' for MATLAB
        rate = re.sub(r'(?<=\S)\*(?=\S)', '.*', rate)
        
        #-----------------------------------------------------------------------
        # If its NOT Photolysis,  add 'Met' as an argument to any rate functions:
        #-----------------------------------------------------------------------
        if data['is_photo']==False:
            
            # Figure out if the rate is a number or a function by trying to 
            # convert the string holding it into a floating point number... 
            try:
                new_r = float(rate)
            except ValueError:
                new_r=rate 
            
            # If the rate  does indeed convert to a number without an error, then 
            # we'll write the rate constant assignment directly into the 
            # GEOSChem_AllRxns.m file (NOT as a string)...  just adding the 
            # 'k(:,i)=' and ';' part of the F0AM lines tp the rate constant. 
            if isinstance(new_r, float):  
               data['K_Line']=f"k(:,i) = {rate};" 
            else:
                
                # Since its a func, we'll make it a named rate in the GEOSChem_AllRxns.m 
                #(e.g. written as a string) that is later defined in the 
                # GEOSChem_K.m & import_rates_GC.m files... 
                data['K_Line']=f"k(:,i) = '{rate}';" 
                
                # If it IS a function, then add 'Met' as an argument for all rate 
                # constant functions so we can indeed use Temp/ M/RH as vars in
                # the rate constant functions & add it to key/value pair: 
                if data['is_het']==False: 
                    if "'"+rate+"'" not in list(foam_dict['nhet']['rate_dict'].keys()): 
                        foam_dict['nhet']['rate_dict']["'"+rate+"'"]=_add_met_to_rates(rate)
                else: 
                    if "'"+rate+"'" not in list(foam_dict['het']['rate_dict'].keys()): 
                        foam_dict['het']['rate_dict']["'"+rate+"'"]=_add_met_to_rates(rate)
                    

            # Add all reaction rate function names to output dicts too. 
            data['Function']=data['rate_function']
        else: 
            #-----------------------------------------------------------------
            # If it IS Photolysis, map PHOTOL(IND) --> Actual named J value. 
            #------------------------------------------------------------------
            # Use the KPP rate parsed as the "key" to look up the j-value match 
            # in "jmap_dict". Both keys & rates should be formatted as 'PHOTOL(IND)'.
            
            # If this photolysis rate does NOT have a correponding match to a
            # j-value we KNOW is already defined in F0AM within 'jdict', pop an error.
            if rate not in list(jmap_dict['Inds_to_Names'].keys()):
                raise ValueError("In function, prep4f0am(), we could not find a J-Value defined in F0AM \n" +
                                 "within 'jmap_dict' corresponding to: '"+rate +
                                 "'\n\nThis likely means you need to update 'jmap_dict' to account for new F0AM/FJX mapping j-value changes!")
            else:  
                # Otherwise, use the mapped value in J-Dict to define rate in F0AM. 
               
                # Pull out comment on where the value used in PHOTOL(IND) comes from: 
                jind = fjx_df.index[fjx_df['PHOTOL(IND)'] == rate.replace(' ', '')]
                j_info = fjx_df.loc[jind[0], 'Info']
                
                # Add the info on the rxn / prepre the rate for F0AM as we did abmove. 
                data['K_Line']=f"k(:,i)={jmap_dict['Inds_to_Names'][rate]}; {j_info}"
                
                # Likewise, keep up with the Named Js too. 
                data['Function']=jmap_dict['Inds_to_Names'][rate]
        
        #######################################################################
        # (3) Format G-Line:  "Gstr{i,1} = 'OH';  Gstr{i,2} = 'HO2';" 
        #######################################################################
        g_strs=[] # Make empty list to hold a string for each reactant in this rxn. 

        # Loop over all reactants  that AREN't Fixed variables in the rxn. 
        mod_rcts=[cmpd for cmpd in data['rct_cmpds'] if cmpd not in  fixed_vars]
        for n, rct_i in enumerate(mod_rcts):
            # Craft a g-str for it keeping track of index properly (rct 1 or 2?).
            g_strs.append("Gstr{i,"+str(int(n+1))+"}='"+rct_i+"'; ")
            
        # Join all together with a space to get a single G-Line for this rxn.
        data['G_Line'] =' '.join(g_strs)
        
        #######################################################################
        # (4) Format F-Line:  "fOH(i) = fOH(i)-2.0; fH2O2(i) = fH2O2(i)+1.0;" 
        ####################################################################### 
        # Accumulate all the stoichiomentry onto single f-line expression: 
        data=_format_f_line(data,fixed_vars)
        
        #######################################################################
        # (5) WRAP THE LINES TOGETHER 
        ####################################################################### 
        # Decide whether or not the rxn should be commented out: 
        char2add='%' if data['is_het']==False and indx in inds2drop else '' 
            
        data['F0AM_Lines']='\n'.join([char2add+'i=i+1;',
                                 char2add+data['Rxn_Line'],
                                 char2add+data['K_Line'],
                                 char2add+data['G_Line'],
                                 char2add+data['F_Line']])
                                
        #######################################################################
        # (6) ORGANIZE OUTPUT & WRITE TO OUTPUT DICT:  
        ####################################################################### 
        # Define what will be keys of the subdicts for the het/non-het dictionaries. 
        fdict_keys=['F0AM_Lines','Rxn_Line' , 'K_Line', 'G_Line', 'F_Line']
        
        # ---------------------------------------------------------------------
        # If the reaction is a NON-HET RNX, put info in thet non-het output dict:
        # ---------------------------------------------------------------------
        if data['is_het']==False: 
            
            # Decide what index # this is (e.g. where info will go in all lists): 
            nhet_ind=0 if len(foam_dict['nhet']['lines']) ==0 else list(foam_dict['nhet']['lines'].keys())[-1]+1
            
            # Add info about this rnx index to the output dicts for all keys: 
            foam_dict['nhet']['lines'][nhet_ind]={k:data[k] for k in fdict_keys}
            
            if indx in inds2drop: 
                # If we ARE commenting this rxn out, then add its index to the list we're 
                # commenting out so we can group writing info about these xns together in mech: 
                foam_dict=_append_to_ndict_list(foam_dict, 'nhet', 'comment_out', nhet_ind, allow_dupes= False)
                
            else: 
                # If we are NOT commenting out this reaction b/c it DOESN't include 
                # het only reactants, then also add the species involved & required 
                # K's/Js to the species list / list of required Ks & Js. 
                
                # Add species involved in this reaction to non-het mech species list if we don't have them already:   
                foam_dict=_append_to_ndict_list(foam_dict, 'nhet', 'sp_list', data['rct_cmpds']+data['prd_cmpds'], allow_dupes= False)
                
                # Add required k's or js from this rxn to master list for non-het mech if we don't have them already: 
                if data['is_photo'] ==False: 
                    foam_dict=_append_to_ndict_list(foam_dict, 'nhet', 'req_ks', data['Function'], allow_dupes= False)
                else: 
                    foam_dict=_append_to_ndict_list(foam_dict, 'nhet', 'req_js', data['Function'], allow_dupes= False)
            
        # ---------------------------------------------------------------------
        # OTHERWISE.. if the reaction is a HET RNX, put info in the het output dict:
        # ---------------------------------------------------------------------
        else: 
                                         
            # Add F0AM lines to dict about this reaction: 
            het_ind=0 if len(foam_dict['het']['lines']) ==0 else list(foam_dict['het']['lines'].keys())[-1]+1
            foam_dict['het']['lines'][het_ind]={k:data[k] for k in fdict_keys}
            
            # Add species involved in this reaction to het mech species list if we don't have them already:   
            foam_dict=_append_to_ndict_list(foam_dict, 'het', 'sp_list', data['rct_cmpds']+data['prd_cmpds'], allow_dupes= False)
            
            # Add required k's or js from this rxn to master list for non-het mech if we don't have them already: 
            if data['is_photo'] ==False: 
                foam_dict=_append_to_ndict_list(foam_dict, 'het', 'req_ks', data['Function'], allow_dupes= False)
            else: 
                foam_dict=_append_to_ndict_list(foam_dict, 'het', 'req_js', data['Function'], allow_dupes= False)
    
    #######################################################################
    # (7 ) Renove het only species from  non-het sp list & get # of sp/rxns: 
    ###################################################################### 
    # Remove any species declared in the nhet mech from the het mech: 
    foam_dict['het']['sp_list']= [sp for sp in foam_dict['het']['sp_list'] if sp not in foam_dict['nhet']['sp_list']]

    # Add count # of sp and # of rxns in each & add to output dict: 
    for k in ['nhet', 'het']:
        foam_dict[k]['n_sp']=len(foam_dict[k]['sp_list'])
        foam_dict[k]['n_rxns']=len(foam_dict[k]['lines'])
        
    return foam_dict


# %%###########################################################################
# -----------------(A) CREATE THE GEOSCHEM_AllRxns.m File ---------------------
###############################################################################

# A.5 L1-function called within make_GEOSCHEM_AllRxns_file()
def edit_GCRxns_template(user_inp:dict, foam_dict:dict,  template_file:str,
                         has_excluded_rxns:bool,
                         het_species_list=[]):
    """Function to read in the GEOSChem_Rxns.m templates (gas/het) and edit them 
    with info specific to this mechanism as needed. Returns list of lines to write. 

    INPUTS: 
    -------
       (1) user_inp -DICT with, at minimum the following keys: 
           
              'kppfile'       - STR with full path to the KPP input file that's being parsed. 

              'GC_version'    - STR formatted as '13.3.4' corresponding to what GEOS-Chem 
                             version you are seeking to create a F0AM file for. 
                             
              'jmap_type'     - STR indicating how J-values were mapped for mechanism. 

       (2) foam_dict     - DICT containing all info about mechanism for F0AM 
                           (NOTE: Input should be subdict held in under key 'het' or nhet'). 
                                   
       (3) tempalte _file  -  STR With path to template file to read in / edit  
       
       (4) has_excluded_rxns     - BOOL indcating if ther eare excluded rxns to consider/ comment out
       
       (5) het_species_list - List of Het only species (used only when doing nhe list)... 

       

    OUTPUTS: 
    --------
       (1) lines2write    - List of lines that should be written to 'GEOSCHEM_Rxns.m'

    REQUIREMENTS: 
    -------------
        LIBRARIES:          from datetime import date  

        CUSTOM FUNCTIONS:    utils.join_list_for_MATLAB()
                             utils._edit_line()

    AUTHOR: 
    -------
        Prof. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD
    """
    # Unpack vars from user_inp dict referenced herein: 
    kppfile= user_inp['kppfile'] 
    GC_version= user_inp['GC_version'] 
    jmap_type= user_inp['jmap_type'] 
     
    # Create a formatted list of all the required rate functions (K's) that must 
    # be imported & named J's. This info goes IN THE HEADER (commented out). Format 
    # it so both are in in a commented out comma delimited list with line breaks 
    # where  necessary (for long MATLAB lists that continue on multiple lines). 
    hdr_k_list= utils.join_list_for_MATLAB(',', [k+'()' for k in foam_dict['req_ks']], # Add '()' after func name... 
                                           insert='    ', adj_ln1_width=33,
                                           insert_skip=[1],comment=True, 
                                           comment_skip=[1]) 
    
    hdr_j_list= utils.join_list_for_MATLAB(',', foam_dict['req_js'],
                                           insert='    ',  adj_ln1_width=33,
                                           insert_skip=[1],comment=True, 
                                           comment_skip=[1]) 
   
    # Create nicely formatted list of all species referenced in the mechanism
    # used in the "SpeciesToAdd" line... 
    sp_list=utils.join_list_for_MATLAB(';', ["'"+sp+"'" for sp in sorted(foam_dict['sp_list'])], 
                                       comment=False, add_semicolon=False, adj_ln1_width=33)
   
    # Also create list of heterogeneous (only) species: 
    if len(het_species_list) >0: 
        het_species_list=utils.join_list_for_MATLAB(';', ["'"+sp+"'" for sp in sorted(het_species_list)], 
                                           comment=True, add_semicolon=False, 
                                           insert='       ',adj_ln1_width=33,
                                           insert_skip=[1],comment_skip=[1])
    else: 
        het_species_list= '{NONE FOUND};'
        
    ###########################################################################
    # Read in GC_Rates Template and customize the header used ...
    ###########################################################################
    # Initialize list of lines we'll write to the output file...
    lines2write = list([])

    # Figure out how many lines there are in the "Template" file so you know when 
    # you reach the end.
    num_lines = sum(1 for line in open(template_file))
    
    # Initilize a list to hold comments about excluded reactions & a bool that'll 
    # tell us if we're parsing that part of the file yet or not. 
    excluded_rxn_cmts = [];   marker_found = False

    # Using readlines(), open the template file and read in line by line.
    inF = open(template_file, 'r')
    lines = inF.readlines()
    count = 0  # Initialize line counter variable.
    
    for line in lines:
         # On most lines we don't need to modify anything, but there are several
         # "Trigger" strings where we need to insert info about the mechanism.
         # pass line to func that looks for trigger & modifies line accordingly if found. 
        
        # Add to HEADER what version of GEOS-Chem this mechanism is based on to header: 
        line = utils._edit_line(line, '**GC_VERSION**', GC_version) 
        
        # Add to HEADER what date this fiel was generated on in the function header
        line = utils._edit_line(line, '**THE_DATE**', str(date.today())) 
        
        # Add the path to the KPP file that was used to generate this mechanism
        line = utils._edit_line(line, '**KPP_FILE_PATH**', kppfile)  
        
        # Add to HEADER the number of species declared in the mechanism: 
        line = utils._edit_line(line, '**N_SPECIES**', str(int(foam_dict['n_sp']))) 
        
        # Add to HEADER the number of reactions in the mechanism: 
        line = utils._edit_line(line, '**N_RXNS**', str(int(foam_dict['n_rxns']))) 
        
        # Add to HEADER  how j-values were mapped: 
        line = utils._edit_line(line, '**JMAP_TYPE**', jmap_type) 
        
        # Add to HEADER which named J-values are required for mechanism: 
        line = utils._edit_line(line, '**HDR_LIST_OF_REQUIRED_J_NAMES**', hdr_j_list) 
        
        # Add to HEADER which functions are required for rate constants: 
        line = utils._edit_line(line, '**HDR_LIST_OF_REQUIRED_RATE_FUNCTIONS**', hdr_k_list) 
        
        # Add the species list 
        line = utils._edit_line(line, '**SPECIES_LIST**', sp_list) 
        
        # Check if the marker telling us the rest of the lines are abput the 
        # excluded reactions has (already) been found on subsequent lines or not: 
        if marker_found:
            # Add the (het) species list where appropriate
            line = utils._edit_line(line, '**HET_SPECIES_LIST**', het_species_list) 
            
            # Add all subsequent lines about the excluded rxns to the list: 
            excluded_rxn_cmts.append(line.strip()+'\n')
            
        # If not, check to see if the marker is in this line or not: 
        elif '**STOP_IF_NO_EXCLUDED_RXNS**' in line:
             marker_found = True
        else: 
            # If neither of these, then we're still parsing some other file part, 
            # so after modifying the lines (if needed), save the line to output list!
            lines2write.append(line)

        count += 1  # Update counter.

        if count > num_lines:
            break  # Exit while loop if you reach end of file.

    # Close the template file now that you got all the good stuff out.
    inF.close()
    
    # Decide if we should actually include the comments about the excluded 
    # reactions in the returned header lines or not & add them in if so: 
    if has_excluded_rxns is True: 
        lines2write=lines2write+excluded_rxn_cmts
    
    return lines2write


def make_GEOSCHEM_Rxns_file(user_inp:dict, foam_dict:dict, output_dir: str = '', 
                            overwrite: bool = False):
    
    """Function to convert a GEOS-Chem KPP File to a F0AM Compliant Mechanism. 
    Prelim outputs for het/non-het mechanism contained in foam_dict. If the user 
    does NOT want to include het rxns, then we find stuff only created in the 
    preliminary mech from a het rxn andremove that chemistry & rebuild the prelim 
    mech! Next, we build a custom header for the GEOSChem_AllRxns.m file and writes
    the final mechanism with its header to the output file. This function OUTPUTS 
    two things: a list of all unique reaction rate definations that must be defined
    in the GEOSChem_K.m file and a dictionary with keys about each unique rate 
    call and what rxns use that rate so that we can build the reaction rates 
    file needed in F0AM.     
    INPUTS: 
    -------
    (1) user_inp - DICT containing at minimum the following keys: 
        'kppfile'          - STR with full path to the KPP input file that's being parsed. 
    
        'GC_version'       - STR formatted as '13.3.4' corresponding to what GEOS-Chem 
                            version you are seeking to create a F0AM file for. 
    
        'jmap_type'        - STR indicating HOW the j-value mapping from GEOS-Chem to 
                            MCM j-values should be done. Default is 'CrossSecMatch'.
    
        'include_het_rxns' - BOOL indicating whether or not heterogenous reactions should 
                                be included in the output mechanisms or not  
    
        'Template_Paths' - STR with absolute paths to dif template files. 
    
    (2) foam_dict - Prelim dict with info on all rxns to write... 

    (3) output_dir     - (OPTIONAL) STR with full path where the outputfile should be written. 
                                Default is set to: "Your_Path_To/GC_Emulator/Output_Files/"

    (4) overwrite      - (OPTIONAL) BOOL of whether or not to overwite the output file or 
                                create a new one if it exists already. Default is False. 
    
    REQUIREMENTS: 
    -------------
        LIBRARIES:           import numpy as np 

        CUSTOM FUNCTIONS:    
                             make_header_for_Rxns()
                             utils.check_filename()
    AUTHOR: 
    -------
        Dr. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD
    """
    ###########################################################################
    #   Write the  GEOSCHem_GasRxns.m File (for all Non-Heterogeneous rxns): 
    ###########################################################################
    # Decide if the non-het/ gas phase mechanism has any excluded rxns or not 
    # we need to comment out in the GEOSCHEM_GasRxns.m file: 
    has_excluded_rxns = True if len(foam_dict['nhet']['comment_out']) > 0 else False 
    
    # Customize the lines that need to be written to "GEOSChem_GasRxns.m"
    # before the first reaction starts (Header, Species List, RO2 list & call to 
    # "Add Species").
    nhet_header= edit_GCRxns_template(user_inp,
                                      foam_dict['nhet'], 
                                      has_excluded_rxns=has_excluded_rxns, 
                                      template_file=user_inp['Template_Paths']['GasRxns'], 
                                      het_species_list= foam_dict['het']['sp_list'])
                                                                     
    # Decide what to call the output file & make sure the path to it exists...
    gas_rxns_file = utils.check_filename(default_name='GEOSCHEM_GasRxns.m', 
                                         ext='.m',savepath=output_dir, 
                                         overwrite=overwrite, return_full=True)
    # Open the Output file to Write
    out_gasF = open(gas_rxns_file, "w")

    # Write the (customized) header lines to the file: 
    [out_gasF.write(line) for line in nhet_header]
    
    # Skip a line
    out_gasF.write('\n') 
    
    # Write the commented out reactions (involving Het Rxns ONLY) first: 
    if has_excluded_rxns: 
        for ind in foam_dict['nhet']['comment_out']: 
            d=foam_dict['nhet']['lines'][ind]
            out_gasF.write(d['F0AM_Lines']+'\n\n') 
            
        # Write section break seperating the commented out rxns & the normal ones: 
        out_gasF.write('% '+'='*68+'\n') 
        out_gasF.write( '% GAS-PHASE / PHOTOTLYSIS ONLY REACTIONS \n')   
        out_gasF.write('% '+'='*68+'\n') 
       
    # Now begin writing the non-het normal rxn lines: 
    for i,d in foam_dict['nhet']['lines'].items(): 
        # Only write this rxn if this rxn has not already been written above
        # in the commented out section: 
        if i not in foam_dict['nhet']['comment_out']:
            out_gasF.write(d['F0AM_Lines']+'\n\n')

    out_gasF.close()  # Close the output file
    print('GEOSCHEM_GasRxns.m file for F0AM saved at: \n\t'+gas_rxns_file)
    
    if user_inp['include_het_rxns'] is True:
        ###########################################################################
        #   Write the  GEOSCHem_HetRxns.m File (for all Het rxns): 
        ###########################################################################
        het_header=edit_GCRxns_template(user_inp, foam_dict['het'], 
                                           has_excluded_rxns=has_excluded_rxns, 
                                           template_file=user_inp['Template_Files']['HetRxns'])
            
        # Decide what to call the output file & make sure the path to it exists...
        het_rxns_file = utils.check_filename(default_name='GEOSCHEM_HetRxns.m', 
                                             ext='.m',savepath=output_dir, 
                                             overwrite=overwrite, return_full=True)
    
        # Open the Output file to Write
        out_hetF = open(het_rxns_file, "w")
    
        # Write the header lines to the file: 
        [out_hetF.write(line) for line in het_header]
        
        # Skip a line
        out_hetF.write('\n') 
        
        # Now begin writing the non-het rxn lines: 
        for i,d in foam_dict['het']['lines'].items(): 
            out_hetF.write(d['F0AM_Lines']+'\n\n')   
            
        out_hetF.close()  # Close the output file
        print('GEOSCHEM_HetRxns.m file for F0AM saved at: \n\t'+het_rxns_file)
        
    return 


# %%###########################################################################
# -----------------(B) CREATE THE GEOSCHEM_K.m File ----------------------------
###############################################################################
def make_GEOSCHEM_K_file(kppfile:str,  GC_version:str, template_file: str, foam_dict: dict, 
                          output_fname: str = 'GEOSCHEM_K.m', output_dir: str = '', overwrite: bool = False):
    """ 
    Function that creates the F0AM 'GEOSCHEM_K.m' file where all rate functions 
    used in the F0AM 'GEOSCHEM_AllRxns.m' are defined. 

    This File is created by reading in the "template" file for GEOSCHEM_K.m, 
    replacing all lines with info specific to this mechanism, and defines rates 
    with "names" that are just STRINGS of the call to the actual function and 
    values that point to a handle of that function that's been converted from 
    the KPP's rate def in fortran, converted to MATLAB. The call to the corresponding
    MATLAB function also passes the input "Met", the MATLAB structure with Temp, M,
    etc. info so that function can access those vars. 

    INPUTS:
    -------

        kppfile       - STR with full path to the KPP input file mehch was built from. 

        GC_version    - STR formatted as '13.3.4' corresponding to what GEOS-Chem 
                        version you are seeking to create a F0AM file for. 

        template_file - STRING containing the full path to the template file for 
                        GEOSChem_K_template.txt

        rxn_rate_dict - DICT with reaction rate names/ functions used in mechansm
    
        output_fname  - (OPTIONAL) STR with name of the outputfile to write. 
                        Default is 'GEOSCHEM_K.m'

        output_dir    - (OPTIONAL) STR with full path where the outputfile should 
                         be written. Default is to "GC_Emulator/Output_Files/"

        overwrite     - (OPTIONAL) BOOL of whether or not to overwite the output 
                        file or create a new one if it exists already. Default is False. 


    OUTPUTS:
    --------
        A F0AM compliant rate defintion file saved at output_dir+output_fname. 

    REQUIREMENTS: 
    -------------
        LIBRARIES:           from datetime import date  

        CUSTOM FUNCTIONS:    utils.join_list_for_MATLAB()
                             utils.check_filename()
                             utils._edit_line() 

    USAGE: 
    -----  
        CALLED WITHIN: 

        OUTPUT USAGE: 


    AUTHOR: 
    -------
            Dr. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD


    """
    
    rxn_rate_dict=foam_dict['nhet']['rate_dict']
    
    
    # Names of (nique) rate functions required in mechanism: 
    ks_required=foam_dict['nhet']['req_ks']
    
    # Create a formatted list of all the required rate functions (K's) that must 
    # be imported & named J's. This info goes IN THE HEADER (commented out). Format 
    # its in a commented out comma delimited list with line breaks 
    # where necessary (for long MATLAB lists that continue on multiple lines). 
    
    k_functs= [k+'()' for k in ks_required] # Add '()' after func name... 
    hdr_k_list= utils.join_list_for_MATLAB(',', k_functs,
                                           insert='    ', 
                                           insert_skip=[1],comment=True, 
                                           comment_skip=[1]) 
    
    # Now, create the correctly formatted list w/ same info that is actually used 
    # to IMPORT THESE RATE FUNCTION HANDLES in the actual MATLAB code 
    # unlike above, this is NOT commented/ in header of the file:  
    ks_as_strs=["'"+k+"'" for k in ks_required] 
    import_rate_list= utils.join_list_for_MATLAB(',', ks_as_strs, 
                                          adj_ln1_width=38,
                                          insert='   ', 
                                          insert_skip=[1])
    
    rate_list= utils.join_list_for_MATLAB(',', ks_required, insert='   ', 
                                          insert_skip=[1])
    
    # Create (all) the lines that define rates with names in this file like these:
    # Normal F0AM Example |   What ours will look like:
    # --------------------|-------------------------------------
    # i=i+1;              |    i=i+1;
    # Knames{i} = 'KDEC'; |    Knames{i} = 'GCARR_ac(1300,5)';
    # krx(:,i) = 1e6;     |    krx(:,i) = GCARR_ac(Met,1300,5);
    
    rate_def_ls = list([])
    
    # Loop over the string assigned as the Kname in the GEOSCHem_Rxns., files (keys)
    all_k_names=list(rxn_rate_dict.keys())
    for k_i, kname_i in enumerate(all_k_names):
        
        iline = "i=i+1;\n"  
        k_dec = "Knames{i} = "+kname_i+";\n"
        k_def= "krx(:,i) = "+rxn_rate_dict[kname_i]+";\n"
    
        # Add this indv rate defintion to the list of rate_defs
        rate_def_ls.append(iline+k_dec+k_def)
    
    # Join all rate_defs into a single string....
    rate_defs = '\n'.join(rate_def_ls)
   
    ###########################################################################
    # Read in GC_Rates Template and sub in specific info for our versio nof the file.
    ###########################################################################
    # Initialize list of lines we'll write to the  output file...
    all_k_lines=list()

    # Figure out how many lines there are in the "Template" file so you know when you reach the end.
    num_lines = sum(1 for line in open(template_file))
    
    # Using readlines(), open the template file and read in line by line.
    inF = open(template_file, 'r')
    lines = inF.readlines()
    count = 0  # Initialize line counter variable.
    
    for line in lines:
         # On most lines we don't need to modify anything, but there are several
         # "Trigger" strings where we need to insert info about the mechanism.
         # pass line to func that looks for trigger & modifies line accordingly if found. 
        
        # Add to HEADER what version of GEOS-Chem this mechanism is based on to header: 
        line = utils._edit_line(line, '**GC_VERSION**', GC_version) 
        
        # Add to HEADER what date this fiel was generated on in the function header
        line = utils._edit_line(line, '**THE_DATE**', str(date.today())) 
        
        # Add the path to the KPP file that was used to generate this mechanism
        line = utils._edit_line(line, '**KPP_FILE_PATH**', kppfile)  
       
        # Add to HEADER which functions are required for rate constants: 
        line = utils._edit_line(line, '**HDR_LIST_OF_REQUIRED_RATE_FUNCTIONS**', hdr_k_list) 
    
        # Update actual MATLAB CODE to import only the handles of the rate
        # functions used in the mechanism in the rate dictionary 
        line = utils._edit_line(line, '**LIST_OF_RATE_FUNCTS_TO_IMPORT**', import_rate_list) 
    
        # Add the line that unpacks the handles from the rate dictionary for use
        line = utils._edit_line(line, '**LIST_OF_IMPORTED_FUNCTION_HANDLES**', rate_list) 
        
        # Add number of Ks required to create correct len arrays in F0AM. 
        line = utils._edit_line(line, '**NUMBER_OF_RATE_CONSTANTS**', str(int(len(all_k_names)))) 
        
        # Add to HEADER which functions are required for rate constants: 
        line = utils._edit_line(line, '**IMPORT_RATES_CALL**', hdr_k_list) 
        
        # Add all lines defining rates to the file in the rate def section: 
        line = utils._edit_line(line, '**INSERT_RATE_DEFINITIONS**',rate_defs ) 
       
        # Add (now possibly modified) line to list of lines to write: 
        all_k_lines.append(line)
        
        count += 1  # Update counter.

        if count > num_lines:
            break  # Exit while loop if you reach end of file.

    # Close the template file now that you got all the good stuff out.
    inF.close()
   
    # Decide what to call the output file & make sure the path to it exists...
    k_file = utils.check_filename(filename=output_fname, default_name='GEOSChem_K.m', 
                                         ext='.m',savepath=output_dir, 
                                         overwrite=overwrite, return_full=True)
    
    # Open the Output file to Write
    out_kF = open(k_file, "w")
    
    # Write the (customized) lines to the file: 
    [out_kF.write(line) for line in all_k_lines]
    
    # Close the new KK file and print where it is: 
    out_kF.close()
    print('GEOSCHEM_K.m file for F0AM saved at: \n\t'+k_file)
    
    return

def make_GEOSCHEM_J_file(fjx_version:str,  GC_version:str, template_file: str, 
                         foam_dict: dict, fjx_df, 
                         output_fname: str = 'GEOSCHEM_J.m',
                         output_dir: str = '', overwrite: bool = False):  
    
    
    # Create a formatted list of all the required named J's. This info goes 
    #IN THE HEADER (commented out). Format it so its in a commented out comma 
    # delimited list with line breaks where  necessary (for long MATLAB lists 
    # that continue on multiple lines). 
    hdr_j_list= utils.join_list_for_MATLAB(',', foam_dict['nhet']['req_js'],
                                           insert='    ',  adj_ln1_width=33,
                                           insert_skip=[1],comment=True, 
                                           comment_skip=[1]) 
   
    # Create lists to hold the lines assiging Js 
    # with direct MCM analogoues & others (seperately) 
    mcm_js=[];not_mcm_js=[];
    
    # Loop over all the Js required for the Non-Heterogenous Rxns: 
    for j_name in foam_dict['nhet']['req_js']:
        # Find all the info on this J-nickname in the FJD_df: 
        row=fjx_df.loc[fjx_df['J-NickName']==j_name.replace("'",'')]
        
        # Pull out the line to assign this named value to a known value: 
        assign=row['GEOSChem_Jline']
        
        if len(assign)>0: # Only write it if it's not empty... 
            line ='J.'+assign.values[0]+'\n' # a dd a space after this line. 
    
            # Keep all MCM analogous in their own list 
            if row['Pure_MCM_Analogue'].values[0]==True: 
                mcm_js.append(line)
            else: # seperate from the non-MCM analogues: 
                not_mcm_js.append(line)
                
    # Figure out how many lines there are in the "Template" file so you know when 
    # you reach the end.
    num_lines = sum(1 for line in open(template_file))
    
    # Using readlines(), open the template file and read in line by line.
    inF = open(template_file, 'r')
    lines = inF.readlines()
    count = 0  # Initialize line counter variable.
    
    lines_to_write=[] 
    for line in lines:
        
        # Add to HEADER what version of GEOS-Chem this mechanism is based on to header: 
        line = utils._edit_line(line, '**GC_VERSION**', GC_version) 
        
        # Add to HEADER what date this fiel was generated on in the function header
        line = utils._edit_line(line, '**THE_DATE**', str(date.today())) 
        
        # Add to HEADER what version of GEOS-Chem this mechanism is based on to header: 
        line = utils._edit_line(line, '**FJX_FILEPATH**', fjx_version) 

        # Add lines with the MCM direct analogue J-Values here: 
        line = utils._edit_line(line, '**HDR_LIST_OF_REQUIRED_J_NAMES**', hdr_j_list) 

        # Add lines with the MCM direct analogue J-Values here: 
        line = utils._edit_line(line, '**INSERT_DIRECTLY_MAPPED_JS**', '\n'.join(mcm_js)) 
        
        # Add lines with the non-MCM direct analogue J-Values here: 
        line = utils._edit_line(line, '**INSERT_OTHER_JS**', '\n'.join(not_mcm_js)) 
        
        # Append possibly modified line to list to write to output file: 
        lines_to_write.append(line) 
        
        count += 1  # Update counter.

        if count > num_lines:
            break  # Exit while loop if you reach end of file.
            
        
    # Close the template file: 
    inF.close() 
    
    # Decide what to call the output file & make sure the path to it exists...
    j_file = utils.check_filename(filename=output_fname, default_name='GEOSCHEM_J.m', 
                                         ext='.m',savepath=output_dir, 
                                         overwrite=overwrite, return_full=True)

    # Open the Output file to Write
    out_JF = open(j_file, "w")

    # Write the header lines to the file: 
    [out_JF.write(line) for line in lines_to_write]
    
    out_JF.close()  # Close the output file
    print('GEOSCHEM_HetRxns.m file for F0AM saved at: \n\t'+j_file)
    
    return 
                
                
        
        
# %%###########################################################################
# -----------------(C) CREATE THE import_GC_Rates.m File -----------------------
###############################################################################

# C.4 L1-function called inside make_import_GC_rates_file()
def make_GCrates_wrapper_funct(template_file, all_rfuncts):
    """Function that reads in the template for "import_GC_rates()" and modifies 
    a few lines in that top level "wrapper" function that's used to call each 
    FORTRAN--> MATRLAB converted function. Returns a list of lines to write to a file 
    with those cutomizations done.

    INPUTS: 
    -------

        (1) template_file - STR containing the full path to the template file for import_GC_Rates()

        (2) all_rfuncts   - A LIST of all rate functions used in the mechanism that will be defined! 

    OUTPUTS: 
    --------
        (1) master_funct_lines - LIST containing lines to write to the MATLAB Rates file. 
                            Will contain the custom Import_GC_Rates function and a section header 
                            where the FORTRAN--> MATLAB converted functions should live BELOW it. 
    REQUIREMENTS: 
    -------------
        LIBRARIES:           None

        CUSTOM FUNCTIONS:    join_list_for_MATLAB (defined in GEOSCHEM_Emulator/utils.py)

    USAGE: 
    -----  
        CALLED WITHIN: make_import_GC_rates_file()

        OUTPUT USAGE: Used to make wrapper function for output of make_import_GC_rates_file()
                     that is written as output to the import_GC_Rates.m file outputted from that fuct.

    AUTHOR: 
    -------
        Dr. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD

    CHANGE LOG: 
    ----------
    10/3/2023    JDH Created 

    """
    # Initialize Output list.
    master_funct_lines = list([])

    # Figure out how many lines there are in the "Template" file so you know when you reach the end.
    num_lines = sum(1 for line in open(template_file))

    # Using readlines(), open the template file and read in line by line.
    file = open(template_file, 'r')
    lines = file.readlines()

    count = 0  # Initialize line counter variable.

    # Loop over each line...
    for line in lines:
        # On most lines we don't need to modify anything, but there are several
        # "Trigger" strings where we need to insert info about the specific
        # rates that are being defined in this file where we need to make modifications.

        # Update the template's Commented Out Example to contain an example for importing "all" the rates used in this mech.
        if 'CMT_INSERT_FUNCTION_NAMES' in line.strip():
            cmt_list_names = utils.join_list_for_MATLAB(
                ',', all_rfuncts, insert='    %       ')
            cmt_list_names = cmt_list_names[12:]  # remove first '    %       '
            line = line.replace('CMT_INSERT_FUNCTION_NAMES', cmt_list_names)

        # Update the actual function where rate_keys are defined to include each
        # of the actual functions referenced in this file with their names.
        if 'INSERT_INDV_FUNCTION_NAMES' in line.strip():
            rfuncts_as_strs = ["'"+rate+"'" for rate in all_rfuncts]
            list_names = utils.join_list_for_MATLAB(
                ',', rfuncts_as_strs, insert='                ')
            list_names = list_names[16:]  # Remove first spaces...
            line = line.replace('INSERT_INDV_FUNCTION_NAMES', list_names)

        # Update the master function where rate handles are defined as keys in dict!
        if 'INSERT_INDV_FUNCTION_HANDLES' in line.strip():
            handles = ['@'+funct for funct in all_rfuncts]
            list_handles = utils.join_list_for_MATLAB(
                ',', handles, insert='                ')
            list_handles = list_handles[16:]  # Remove first spaces...
            line = line.replace('INSERT_INDV_FUNCTION_HANDLES', list_handles)

        # After modifying lines (if needed) in master function, save each to output list!
        master_funct_lines.append(line)

        count += 1  # Update counter.

        if count > num_lines:
            break  # Exit while loop if you reach end of file.

    # While loop will break once we found the end of the template file and have parsed all lines.
    file.close()  # Close the template file.

    return master_funct_lines

# C.3 L1-function called inside make_import_GC_rates_file()
def redo_function_header(line: str, spl_str: str = 'FUNCTION'):
    """Function used to reformat FORTRANT Function definition lines into MATALB 
    compliant function definitions. Works for two distict cases now: 

                           FORTRAN:                                                 MATLAB
    Case #1     'FUNCTION ARRPLUS_ade(a0,d0,e0) RESULT(k)'            --> 'function [k] = ARRPLUS_ade(a0, d0, e0)'
    Case #2     'REAL(kind=dp) FUNCTION FALL( A0,B0,C0,A1,B1,C1,CF)'  --> 'function [FALL] = FALL(A0,B0,C0,A1,B1,C1,CF) '

    NOTE: Case #1 is triggered by the word "RESULT" appearing in line, all other inputs assume Case #2. 

    INPUTS:
    -------
        (1) line    - STR of the function line you'd like to reformat

        (2) spl_str - STR that the line should be "split" on. Default is 'FUNCTION'

    OUTPUTS: 
    --------
        (1) matlab_dec - Updated STR of input line in MATLAB compliant formatting. 

    REQUIREMENTS: 
    -------------
        LIBRARIES:           None
        CUSTOM FUNCTIONS:    _str_multi_replace (defined in GEOSChem_Emulator/utils.py)

    USAGE: 
    -----  
        CALLED WITHIN: make_import_GC_rates_file()

        OUTPUT USAGE: Used to write function headers in output of make_import_GC_rates_file: 
                     'import_GCrates.m''

    AUTHOR: 
    -------
        Prof. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD

    CHANGE LOG: 
    ----------
        10/03/23: JDH Created 
    """

    if 'RESULT' in line:
       # Turn FORTRAN:
        #   'FUNCTION ARRPLUS_ade( a0, d0, e0 ) RESULT( k )'
       # Into MATLAB:
        #   'function [k] = ARRPLUS_ade(a0, d0, e0 )

        # Split the input line at the word "FUNCTION" to just get the rest:
        # e.g: ['', 'ARRPLUS_ade( a0, d0, e0 ) RESULT( k )']
        l = line.split(spl_str)[1].strip()

        # Get the "result" output var ("k") isolated from the word "RESULT" and "(", and ")"
        l = l.split('RESULT')  # e.g: ['ARRPLUS_ade( a0, d0, e0 )' , '( k )']
        output_var = utils.str_multi_replace(
            l[1].strip(), ['(', ')'], rep_all='').strip()  # e.g: 'k'

        # Now remove spaces from the function title & arguements:
        funct_title_args = l[0].strip()  # e.g:  ARRPLUS_ade(a0,d0,e0)'

        # Add "Met" as an input arg so we can always retrieve the Met Vars in MATLAB as inputs to this function:
        title_ls = funct_title_args.split(
            '(')  # e.g: ['ARRPLUS_ade', 'a0,d0,e0)']
        # e.g: 'ARRPLUS_ade(Met,a0,d0,e0)']
        new_title_args = title_ls[0]+'(Met,'+title_ls[1]

        # Add it all back together how MATLAB expects:
        matlab_dec = 'function ['+output_var+'] = ' + \
            new_title_args.replace(' ', '')

    else:
        # Turn FORTRAN:
        #   'REAL(kind=dp) FUNCTION FALL( A0, B0, C0, A1, B1, C1, CF )'
        # Into MATLAB:
        #   'function [FALL] = FALL(Met,A0,B0,C0,A1,B1,C1,CF)'

        # Split the input line at the word "FUNCTION" to just get the rest & remove spaces
        l = line.split(spl_str)[1].strip()  # e.g: 'FALL(A0,B0,C0,A1,B1,C1,CF)'

        # Split the function on the "(" character:
        title_ls = l.split('(')  # e.g. ['FALL', 'A0,B0,C0,A1,B1,C1,CF)']

        # Isolate just the function name (to get the output var) & remove spaces:
        output_var = title_ls[0].strip()  # e.g: 'FALL'

        # Insert "Met" as an input arg so we can always retrieve the Met Vars in MATLAB as inputs to this function:
        new_title_args = title_ls[0]+'( Met,'+title_ls[1]

        # Add it all back together how MATLAB expects:
        matlab_dec = 'function ['+output_var+'] = ' + \
            new_title_args.replace(' ', '')

    return matlab_dec

# C.2 L1-function called inside make_import_GC_rates_file()
def insert_after_indent(line, insert):
    """Function to insert some txt AFTER the FIRST character of 'line' that isn't a space 
    and BEFORE anything else that appears on that line. Basically, will preserve 
    spaces used for indentation and can be used to insert a MATLAB '%' to make somthing a comment...

    INPUTS: 
    -------
        (1) line   - STR of whatever line you want to insert somthing in. 

        (2) insert - STR of whatever you'd like to insert 

    OUTPUTS: 
    --------
        (1) line - Updated STR of the line with the insertion .... 

    REQUIREMENTS: 
    -------------
        LIBRARIES:            None.
        CUSTOM FUNCTIONS:     None. 

    USAGE: 
    -----  
        CALLED WITHIN: make_import_GC_rates_file() to insert a '%' in front of the 
            first non space character, to perserve MATLAB indentations. 

    AUTHOR: 
    -------
        Prof. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD

    CHANGE LOG: 
    ----------
        10/03/2023 JDH Created.

    """

    # Figure out the index of everything that isn't a space inthe input line:
    not_spaces = [i for i, ll in enumerate(line) if ll != ' ']

    # If what you want to insert is not equal to the very first NON-Space character:
    if line[not_spaces[0]] != insert:
        # Then re-shuffle your line to be: (1) whatever came before the spaces,
        # plus (2) what you want to insert, and then plus (3) whatever came after it...
        line = line[0:not_spaces[0]] + insert + line[not_spaces[0]:]
    return line

# C.1 L1-function called inside make_import_GC_rates_file()
def is_comment(line):
    """Utility Function to determine if a line is supposed to be a comment in MATLAB or FORTRAN

    USED WITHIN: make_import_GC_rates_file()

    AUTHOR: 
    -------
        Prof. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD

    CHANGE LOG: 
    ----------
        10/03/2023 JDH Created.
    """
    not_spaces = [i for i, ll in enumerate(line) if ll != ' ']
    if len(not_spaces) > 0:
        tf = False if line[not_spaces[0]] not in ['%', '!'] else True
    else:
        tf = False

    return tf

# C.1 Function used in to create the import_GC_rates.m output file:
def make_import_GC_rates_file(rate_files: list, all_rfuncts: list, template_file: str,
                              include_het_rxns: bool = False, output_fname: str = 'import_GC_rates',
                              output_dir: str = '', overwrite: bool = False):
    """ Function to convert all functioned defined in a GEOS-Chem gckpp_Rates.F90 file 
    into a MATLAB compliant set of functions contained in a single MATLAB wrapper 
    function that is called in 'GEOSCHEM_K.m' and returns the handles to each 
    individual MATLAB function within it given the name of the function you want to return. 
    This function essientially reads in a gckpp file and converts it into a single place 
    where all rate functions are deinfed in MATLAB. 

    INPUTS: 
    -------
        (1) rate_files       - LIST with the full path and names of all rate files you need to convert. 
                               NOTE: Accepts files like gckpp_Rates.F90 and RateLawFuncs.F90, HetRates.F90, etc. 

        (2) all_rfuncts      - LIST of all unique rate functions that are called in your mechanism 
                               in the KPP.eqn file

        (3) template_file    - STR with full path and name of the import_GC_Rates_template.txt file 
                                that should be stored under: 
                                'Your/Path/To/GEOSChem_Emulator/Input_Files/MATLAB_Templates/importGC_Rates_template.txt'

        (4) include_het_rxns - BOOL indicating whether or not heterogenous reactions should 
                                be included in the output mechanisms or not                   

        (5) output_fname     - (OPTIONAL) STR with name of the output file to write. Default is 'import_GC_rates.m'

        (6) output_dir     - (OPTIONAL) STR with full path where the outputfile should be written. 
                                Default is set to: "Your_Path_To/GC_Emulator/Output_Files/"

        (7) overwrite        - (OPTIONAL) BOOL of whether or not to overwite the output file or 
                                create a new one if it exists already. Default is False. 

    OUTPUTS: 
    --------
       (1)  A MATLAB compliant .m file to define & import all rate functions required 
            saved at output_dir+outfile_name.

    REQUIREMENTS: 
    -------------
        LIBRARIES:           import os
                             import sys

        CUSTOM FUNCTIONS:    is_comment()
                             insert_after_indent()
                             redo_function_header()
                             utils._str_multi_replace()
                             utils.check_filename()
                             make_GCrates_wrapper_funct()

    USAGE: 
    -----  
        CALLED WITHIN: 

        OUTPUT USAGE: 

    AUTHOR: 
    -------
        Prof. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD

    CHANGE LOG: 
    ----------
        01/14/2022 JDH Created.
        10/31/2023 JDH modified documentation

    """
    # -------------------------------------------------------------------------------------
    # ------                     Prep Input and Output Files  ------------
    # -------------------------------------------------------------------------------------

    # File name changes ...  gckpp_Rates for v <13.3.0 or fullchem_RateLaw_Funcs for v >=13.3.0
    mechnm = os.path.basename(rate_files[0]).split('.')[0]

    # -------------------------------------------------------------------------------------
    # ------     Convert all the Fortran functions to matlab functions!  ------------
    # -------------------------------------------------------------------------------------
    key = ''
    all_lines = list()
    fdict = dict({})
    for file in rate_files:
        count = 0
        pass_go = False  # Initialize counter

        # Count # of lines in file the original file so we know when we reached the end.
        num_lines = sum(1 for line in open(
            file, encoding="UTF-8", errors='ignore'))

        # Open the Input file to Read
        inF = open(file, 'r', encoding="UTF-8", errors='ignore')

        while True:  # Loop over each line of the KPP rate file.

            line = inF.readline()  # Read next line from file
            line = line.replace('\n', '')  # Remove new line characters.
            valid_starts = ['RATE-LAW FUNCTIONS FOR GAS-PHASE REACTIONS',
                            'BEGIN RATE LAW FUNCTIONS', 'ARRHENIUS FUNCTIONS']
            valid_ends = ['End INLINED Rate Law Functions', 'END MODULE '+mechnm+'_RateLawFuncs',
                          'RATE-LAW FUNCTIONS FOR HETEROGENEOUS REACTIONS']

            if any([item in line.upper() for item in valid_starts]):
                # Boolean of whether to start writing file or not.
                pass_go = True

            if pass_go is True:
                if len(line) > 0:
                    if '::' in line and 'PARAMETER' in line:  # keep parameter definitions.
                        line = insert_after_indent(line, line.split('::')[
                                                   1].strip() + '&&baddie&&')
                        line = line.split('&&baddie&&')[0]
                    # Comment out all var definitions. Don't need them in MATLAB.
                    if '::' in line and is_comment(line) is False:
                        line = insert_after_indent(line, '%')

                    # Replace "END FUNCTION" with 'end' so it won't think the "end"line is a line with a function to fix...
                    if 'END FUNCTION' in line and is_comment(line) is False:
                        line = 'end'

                    # Replace "END SUBROUTINE" with 'end' so it won't think the "end"line is a line with a function to fix...
                    if 'END SUBROUTINE' in line and is_comment(line) is False:
                        line = 'end'

                    # Make sure we comment out "intrinsic variables.. b/c I have no idea what that means & neither does MATLAB
                    if 'INTRINSIC' in line and is_comment(line) is False:
                        line = insert_after_indent(line, '%')

                    # Move things around at function declarations to make MATLAB happy.
                    if 'FUNCTION' in line and is_comment(line) is False:
                        line = redo_function_header(line)

                    # Only after fixing functions can we make sure REAL declarations are commented out.
                    if 'REAL' in line and is_comment(line) is False:
                        line = '%'+line

                    # Now also comment out declariations of Character vars...
                    if 'CHARACTER' in line and is_comment(line) is False:
                        line = '%'+line

                    if 'WRITE' in line.upper() and is_comment(line) is False:
                        line = "disp('" + "'".join(line.split("'")[-2:])+")"

                    # Check for use of state global vars in Fortran, which use a '%' but in middle of the line, not at beginning
                    if (is_comment(line) is False) and (' %' in line):
                        line = line.replace(' %', '')
                    if (is_comment(line) is False) and ('%' in line):
                        line = line.replace('%', '')

                    # If its a one line "if" statement then need to add an "end" for MATLAB.
                    if 'IF' in line and 'THEN' not in line and 'END IF' not in line and \
                            '&' not in line and 'ENDIF' not in line and 'DIFF' not in line:
                        ln_ls = line.split(')')
                        test = ln_ls[0]
                        if len(ln_ls) == 1:
                            line = test+'); end '
                        elif len(ln_ls) > 1:
                            contents = ')'.join(ln_ls[1:])
                            sub = False
                            for cmt in ['!', '%']:
                                if cmt in contents:
                                    content = contents.split('!')
                                    line = test+');' + \
                                        content[0]+'; end %' + \
                                        ' '.join(content[1:])
                                    sub = True
                            if sub == False:
                                line = test+');'+contents+'; end '
                        else:
                            print(line, len(ln_ls))

                    # Move things around at function declarations to make MATLAB happy.
                    if 'SUBROUTINE' in line and is_comment(line) is False:
                        line = redo_function_header(line, spl_str='SUBROUTINE')

                    # Otherwise, FORTRAN --> MATLAB is a series of operator substitutions...
                    fortran_2_matlab = dict({'!': '%',
                                             '**': '.^',
                                             '*': '.*',
                                             '/': './',
                                             'EXP': 'exp',
                                             'temp': 'TEMP',
                                             'MAX': 'max',
                                             'LOG10': 'log10',
                                             '&': '...',
                                             'ELSE IF': 'elseif', 'ELSEIF': 'elseif',
                                             'ELSE': 'else',
                                             'END MODULE': '% END MODULE',
                                             'ENDIF': 'end',  'ENDif': 'end', 'END IF': 'end',
                                             'IF': 'if',
                                             'DBLE': 'double',
                                             'e+': 'e',
                                             'E-': 'e-',
                                             'E+': 'e',
                                             'RETURN': '',
                                             'SQRT': 'sqrt',
                                             # Not a great solution just yet ..
                                             'CALL': "",
                                             'ENDDO': 'end',
                                             '.TRUE.': 'true', '.FALSE.': 'false',
                                             '.eq.': '==',  '.EQ.': '==',
                                             '.not.': '~',  '.NOT.': '~',
                                             '.lt.': '<',   '.LT.': '<',
                                             '.gt.': '>',   '.GT.': '>',
                                             '.ne.': '~=',  '.NE.': '~=',
                                             '.ge.': '>=',  '.GE.': '>=',
                                             '.le.': '<=',  '.LE.': '<=',
                                             '.and.': '&&', '.AND.': '&&',
                                             '.or.': '||', '.OR.': '||'})

                    line = utils.str_multi_replace(line, fortran_2_matlab)

                    # Remove things like double precision and 'Then' statements - unnecessary.
                    line = utils.str_multi_replace(
                        line, ['_dp', 'THEN'], rep_all='')

                    # Fix anything my substitutitons messed up.
                    line = utils.str_multi_replace(line, dict({'../': './', '..*': '.*', '..^': '.^',
                                                               ' .* ': '.*', ' ./ ': './',
                                                              '); ) ;': '));', '); ;': ');'}))

                    # Make sure all the lines are suppressing their content in MATLAB by adding a ';' where appropriate.
                    if 'end' not in line and 'function' not in line and is_comment(line) is False and line[-3:] != '...':
                        if '%' not in line:
                            if len(line.replace(' ', '')) != 0:
                                line = line+';'  # Make sure stuff is on line if place a ;
                        else:
                            cmt_i = line.index('%')
                            line = line[0:cmt_i]+';'+line[cmt_i:]

                    # Collect all lines from every function in a dictionary ....
                    if not any([True if e in line else False for e in valid_ends]):
                        if line.strip() not in ['%', ';', '']:  # Don't store empties!
                            if 'function [' in line:
                                key = line.split('=')[1].split(
                                    '(')[0].replace(' ', '')
                                fdict[key] = [line]
                            elif line != 'end' and key != '':
                                if line.strip()[0] == '%':
                                    line = '    % '+line.strip()[1:].strip()
                                all_lines.append(line)
                            elif line == 'end' and key != '':
                                fdict[key] = fdict[key]+all_lines+[line]
                                key = ''
                                all_lines = list()

            count = count+1  # Update counter var.

            # Break if reached end of file or end of the functions for rate constants.
            if count > num_lines or any([True if e in line else False for e in valid_ends]):
                break
        inF.close()  # Close the input file

    # Make sure the output file path exists... Set default and and ext incase users don't!
    rate_file = utils.check_filename(filename=output_fname, default_name='import_GC_rates', ext='.m',
                                     savepath=output_dir, overwrite=overwrite, return_full=True)

    outF = open(rate_file, "w")  # Open the Output file to Write

    # -------------------------------------------------------------------------------------
    # ------Write the Wrapper Function for All rates to the Output Rates File! ------------
    # -------------------------------------------------------------------------------------

    # Create the wrapper function that needs to be written to the "rates" file
    wrapper_funct_lines = make_GCrates_wrapper_funct(
        template_file, all_rfuncts)

    # Write this "Wrapper" function to the top of the output rates file.
    for ln in wrapper_funct_lines:
        if '\n' in ln:
            outF.write(ln)
        else:
            outF.write(ln+'\n')

    # Write a line to retrieve all global vars inside each function.
    get_globals = '    % Pass input struct, "Met" to function "extract_MetVars()" to get vars this function may depend on...\n' + \
        '    [TEMP, PRESS, NUMDEN, H2O, TEMP_OVER_K300, K300_OVER_TEMP,INV_TEMP,SR_TEMP,RELHUM]= extract_MetVars(Met);\n'
    written=[] 
    for key in fdict:  # Loop over all functions
        if key in all_rfuncts or key.replace('_sp', '') in all_rfuncts:
            # What its called in our rate dictionary: 
            rate_key=key.replace('_sp', '') if key.replace('_sp', '') in all_rfuncts else key 
            written.append(rate_key)
            
            
            glbs_tf = False
            # list of lines we need to write for each function.
            ln_list = fdict[key]
            for i, line in enumerate(ln_list):
                outF.write(line+'\n')
                if glbs_tf is False and ')' in line:
                    outF.write(get_globals+'\n')
                    glbs_tf = True

    outF.close()  # Close the output file
    print('MATLAB compliant version of gckpp_Rates.F90 called in "GEOSCHEM_K.m" to define all rate functions saved at: \n\t'+rate_file)
    
    print('MISSING:') 
    missing=[print(f'\t{f}()') for f in all_rfuncts if f not in written] 
    return missing


# %%###########################################################################
# -----------------(D) Wrapper Function to Make ALL Files -----------------------
###############################################################################
# D.1 L1-Function used in make_GC_mechanism() to check user inputs to make_GC_mechanism()



def make_GC_mechanism(kppfile, rate_files, GC_version: str, jmap_type: str, include_het_rxns: bool,
                      output_dir: str, overwrite: bool = False, verbose:bool=True):
    """Main function to write a F0AM version of a GEOS-Chem mechanism.

    INPUTS:

        (1) kppfile - string containing the path to a GEOS-Chem KPP file (path_to/fullchem.eqn')
        (2) rate_files- list of strings OR string contaiing the path to a GEOS-Chem a 'gckpp_Rates.F90' file. 
        (3) GC_version - Version of GEOS-Chem being Emulated (so chooses right J-files!). 
        (4) jmap_type - STRING indicating HOW you want J-values to be mapped. 
        (5) include_het_rxns-  bool indicating whether or not you want heterogenous rxns to be included in the resulting mech. 
        (6) output_dir -str containing the path where the output F0AM files should be located. 
        (7) overwrite - bool indicating whether or not output files should overwrite any exisiting files there or not. 
        (8) verbose - bool. Set to true to print off more info about things. 

    """
    
    # Check the user's inputs to make sure all files we need exist and are located 
    # where we expect them. Dump it all into a nice dictionary passed around later. 
    user_inp= check_user_inputs(kppfile=kppfile, 
                                 rate_files=rate_files, 
                                 GC_version=GC_version,
                                 jmap_type=jmap_type, 
                                 include_het_rxns=include_het_rxns,
                                 output_dir=output_dir,
                                 verbose=verbose)
    
    # Read in the right FJX file for this GC-Version and build the j-mapping
    # dictionary & FJX dataframes to use in mapping KPP J's to F0AM J's.
    jmap_dict, fjx_df, user_inp = hv.create_jmap_dict(user_inp)
    
    # Parse the KPP file. Get list of rxns and rates we need to write a F0AM file.
    print('Parsing KPP file...') 
    tracer_info, fixed_vars, kpp_dict = read_kpp(kppfile, output_dir=output_dir)
    
    # If NOT including heterogeneous reactions in output mech, we should be able to 
    # eliminate some reactions from the mechanism (commenting them out). (e.g. Why
    # photolyze ClNO2 if there is none ever made b/c formation is het only?).
    inds2drop=[] # Inialize list of inds to drop from kpp_dict... 
    
    if include_het_rxns is False:
        inds2drop=eliminate_dead_ends(kpp_dict,include_het_rxns=include_het_rxns)
    
    # Prep reactions & rates for F0AM (e.g. turn kpp_dict info into format F0AM 
    # The output "foam dict" has 2 top level keys: 'het' & 'nhet' w/ same info 
    # about stuff that will go in the het/nonhet mechanism files (if het is requested): 
    print('Preparing for F0AM...')
    foam_dict=prep4f0am(kpp_dict, inds2drop, fixed_vars, fjx_df, jmap_dict)

    # Write the GEOSChem_GasRxns.m file (&, if asked, the GEOSChem_HetRxns.m file). 
    make_GEOSCHEM_Rxns_file(user_inp=user_inp,
                            foam_dict=foam_dict, 
                            output_dir=output_dir,
                            overwrite=overwrite)
    
    # Write the GEOSChem_K.m file  
    make_GEOSCHEM_K_file(user_inp=user_inp,
                          foam_dict=foam_dict,
                          output_fname='GEOSCHEM_K.m',
                          output_dir=output_dir,
                          overwrite=overwrite)

    # # Create import_GC_rates.m file where all rate functions are defined which is
    # # referenced and used as is without needing modifications in the GEOSChem_K.m file.
    if type(rate_files) == str:
        rate_files = [rate_files]

    all_rfuncts=foam_dict['nhet']['req_ks']
    fdict = make_import_GC_rates_file(rate_files=rate_files,
                                      all_rfuncts=all_rfuncts,
                                      template_file=template_paths['import_GC_rates'],
                                      include_het_rxns=include_het_rxns,
                                      output_fname='import_GC_rates.m',
                                      output_dir=output_dir,
                                      overwrite=overwrite)

    make_GEOSCHEM_J_file(fjx_version=fjx_version,  
                         GC_version=GC_version,
                         template_file=template_paths['J'],
                         foam_dict=foam_dict,
                         fjx_df=fjx_df,
                         output_fname='GEOSCHEM_J.m',
                         output_dir=output_dir,
                         overwrite=overwrite)
                            
    
    return foam_dict,kpp_dict,inds2drop,fixed_vars, fjx_df, jmap_dict 




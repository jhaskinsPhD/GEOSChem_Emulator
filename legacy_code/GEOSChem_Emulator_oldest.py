#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 12:23:52 2023

@AUTHOR: u6044586
"""

import utils
import os
import sys
import yaml
import re
import numpy as np
import pandas as pd
from datetime import date

sys.path.append(
    '/uufs/chpc.utah.edu/common/home/u6044586/python_scripts/modules/GEOSChem_Emulator/')

sys.path.append(
    '/uufs/chpc.utah.edu/common/home/u6044586/python_scripts/modules/pyF0AM/')

import pyF0AM as fio

ignore_list=[ 'hv']

# %%---------------------------------------------------------------------------
# -------------- (A.4) make_header_for_F0AM_mech() & Sub-Functions ------------
# -----------------------------------------------------------------------------
def make_header_for_GEOSChem_J(photo_functs):
    """Function to make a descriptive header for the GEOSChem_J.m file. 

       (4) photo_functs  - LIST of unique J-values that must be defined in F0AM 
                            for the resulting mechanism to work outputted from KPP_Dict_to_F0AM_IO()  

    """
    cmt_space1 = '%                    '
    cmt_space2 = '%              '
    j_info = '%    j_values={' + ','.join([f if np.mod(
        i, 15) != 0 or i == 0 else cmt_space2+f for i, f in enumerate(photo_functs)])+'};'
    jzero = 'Jzero=0.*J3; % Set photolysis=0 for stratospheric only rxns \n'

    return

# A.5 L1-function called within make_GEOSCHEM_AllRxns_file()
def make_header_for_AllRxns(kppfile, GC_version, r_functs, photo_functs, jmap_type):
    """Function to make a descriptive header for the generated mechanism file. 

    INPUTS: 
    -------
       (1) kppfile       - STR with full path to the KPP input file that's being parsed. 

       (2) GC_version    - STR formatted as '13.3.4' corresponding to what GEOS-Chem 
                             version you are seeking to create a F0AM file for. 

       (3) r_functs      - LIST of all names of pre-defined K-Rates used in the mechanism. 

       (4) photo_functs - LIST of all names of pre-defined J-Rates used in the mechanism. 

    OUTPUTS: 
    --------
       (1) mech_title    -

       (2) mech_comments -


    REQUIREMENTS: 
    -------------
        LIBRARIES:           import numpy as np 

        CUSTOM FUNCTIONS:    pylb.join_list_for_MATLAB

    USAGE: 
    -----  
        CALLED WITHIN: make_GEOSCHEM_AllRxns_file() 

        OUTPUT USAGE: Passed to fio.write_mech_to_file to write mech. 

    AUTHOR: 
    -------
        Prof. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD

    CHANGE LOG: 
    ----------
        10/03/23 JDH Created 
        10/30/23 JDH Modified documentation 
    """
    r_functs = np.unique(r_functs).tolist()

    # Join the list of J-values/ Rate Functions used in this mech to put into the mechanism's header.
    cmt_space1 = '%                    '
    cmt_space2 = '%              '
    rate_info = '%    rate_functions={' + ','.join([f + '()' if np.mod(i, 5) != 0 \
                                                    or i == 0 else cmt_space1+f + '()' \
                                                        for i, f in enumerate(r_functs)])+'};'
    j_info = '%    j_values={' + ','.join([f if np.mod(i, 15) != 0 \
                                           or i == 0 else cmt_space2+f \
                                               for i, f in enumerate(photo_functs)])+'};'

    # Make Title Containing details + list of functions this rxn file requires to run.
    mech_title = ''.join(['% F0AM Compliant GEOS-Chem Mechanism for version'+GC_version + ' generated on '+str(date.today())+' from KPP file: ',
                          '%    '+kppfile + ' using jmap_type =' + jmap_type])

    mech_comments = ''.join(['% List of reaction rates that must be defined in "GEOSChem_K.m" for this mech to work:\n',
                             rate_info,
                            '\n\n',
                             '% List of photolysis rates that must be defined in "GEOSChem_J.m" for this mech to work:\n',
                             j_info])

    return mech_title, mech_comments

# A.4 L1-function called within make_GEOSCHEM_AllRxns_file()
def eliminate_dead_ends(het_info, nhet_info):
    """ID Reactions from stuff that's ONLY formed via a Heterogenous reaction
       Basically cut dead end reactions that we don't need...

       Custom Functions: 
           flatten_nested_list  --> in pyMCM/F0AM_Tools/utils.py
           enforce_list_inds    --> in pyMCM/F0AM_Tools/utils.py
           find_in_list         --> in pyMCM/F0AM_Tools/utils.py
           drop_dupes           --> in pyMCM/F0AM_Tools/utils.py
           """

    # Unpack stuff from het / non-het mechanism: 
    [het_sp_list, het_rxns, fhet_s, het_gs, het_ks, het_rct_cmpds, 
     het_prd_cmpds, het_rct_ylds, het_prd_ylds] = het_info
    
    [nhet_sp_list, nhet_rxns, nfhet_s, nhet_gs, nhet_ks, 
     nhet_rct_cmpds, nhet_prd_cmpds, nhet_rct_ylds, nhet_prd_ylds] = nhet_info
    
    # A list of all products formed in any het rxns ...
    het_prds = utils.flatten_nested_list(het_prd_cmpds, drop_dupes=True)

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

        # Get a unique list of reactants and products that are still left in the non-het rxns.
        prds = utils.flatten_nested_list(prd_a, drop_dupes=True)

        # Get size of list of rxn inds to drop on this iteration of the while loop. 
        sz0 = len(inds2drop)  

        # Find all species formed in het rxns that are NOT formed from any non-het rxns
        dead=[]
        [dead.append(sp) for sp in het_prds if (sp not in prds) and (sp not in ignore_list+dead)]

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

    # Get a list of the indexes of the rxns we will KEEP in the nhet mech: 
    inds2keep = [indx for indx in range(0, len(nhet_rxns)) if indx not in inds2drop]
    rxns2drop = [nhet_rxns[indx] for indx in inds2drop]

    # Keep the reactions that you don't want to drop! 
    nhet_info_NEW = utils.enforce_list_inds(inds2keep, nhet_info[1:]) # not species list (not same len).
    new_rcts=nhet_info_NEW[4];new_prds=nhet_info_NEW[5];
    
    # Make a new species list using stuff in the NEW reactants and products list.
    new_sp = utils.drop_dupes(
        utils.flatten_nested_list(new_rcts, drop_dupes=True)+\
        utils.flatten_nested_list(new_prds, drop_dupes=True))
        
    # Remove any species that aren't supposed to be declared in mech (in ignore list)... 
    [new_sp.remove(item) for item in ignore_list+['hv'] if item in new_sp]
    new_sp.sort()
    
    # Package it all up into output info dump (same format of nested lists as input).  
    nhet_info_NEW=[new_sp]+nhet_info_NEW
    
    return nhet_info_NEW, rxns2drop, inds2drop

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


# A.3.1 L2-function called within KPP_Dict_to_F0AM_IO within make_GEOSCHEM_AllRxns_file()
def convert_FJX_PhotoRxns_to_MCM(rxn_df, jmap_dict, fjx_df):
    """Function used to format photolysis reactions and convert the rates for 
    photolysis reactions into lists we can use to build an MCM compliant mechanism. 
    Specifically, this is where the actual mapping of FJX/KPP photolysis rates into 
    their corresponding MCM photolysis rate is done. 

    INPUTS: 
    ------- 
       (1) rxn_df -             Pandas df of info about all reactions generated by 
                                the function "parse_KPP". 

       (2) jmap_dict  -         Mapping Dictionary with keys corresponding to FJX phototlysis rates and 
                                values corresponding to F0AM corresponding rate values.

       (3) fjx_df     -         Pandas df outputted from create_jmap_dict() 

    OUTPUTS: 
    --------
        (1) rxn_df -             Updated df with J rates having their wFJX values converted to 
                                 their F0AM counterpart and formatted with k(:,i) in front of them 
                                 and with a ';' behind! 


    REQUIREMENTS: 
    -------------
        LIBRARIES:           import sys 
                             import numpy as np 

        CUSTOM FUNCTIONS:    

    USAGE: 
    -----  
        Called Within: KPP_Dict_to_F0AM_IO() within make_GEOSCHEM_AllRxns_file(). 

        Output Usage: Directly passed as output of KPP_Dict_To_F0AM_IO() 

    AUTHOR: 
    -------
        Dr. Jessica D. Haskins (jessica.haskins@utah.edu)   GitHub: @jhaskinsPhD 

    CHANGE LOG: 
    ----------
        08/10/2022: JDH Created 
        10/02/2023: JDH modified to use inputs to 'photo_dict' and pass J-Value mapping work to 
                    updated function, get_all_MCM_Js_needed(). 
    """

    # Extract info about photolysis rxns only: 
    photo_df =rxn_df[rxn_df['is_photo']==1]
        
    # Check that all rate functions in user input photo_dict['rate_function'] are
    # indeed "PHOTOL" (aka that all rates used in what we think are the
    # photolysis reactions are actually photolysis rates. Otherwise pop an error!
    unq_rates = np.unique(photo_df['rate_function'])
    
    if len(unq_rates) != 1 or unq_rates[0] != 'PHOTOL':
        other_rates = '\n'.join(
            ['    '+rate for rate in unq_rates if rate != 'PHOTOL'])
        raise ValueError("Not all rate functions in the 'photo_dict['rate_functions]' list passed to the function, \n" +
                         "convert_FJX_PhotoRxns_to_MCM(), are 'PHOTOL' as expected.+ \n" +
                         "The following rate functions were also found: "+other_rates +
                         "\n\n  This likely means something weird is happening with: \n" +
                         "    (1) How you detect a rxn is a photolysis rxn in the function, parse_KPP(). \n" +
                         "                   OR \n   " +
                         "    (2) That something weird is happening in the function, _enforce_list_dict_inds(), \n"
                         "        which is used to sub-select photolysis info in the dict, 'kpp_dict', within the function, parse_KPP()").with_traceback(sys.exc_info()[2])

    # Loop over all rates used in photolysis reactions:
    for i in photo_df.index:
        
        # Use the KPP rate parsed as the "key" to look up the j-value match in "jdict".
        # Both keys & rates should be formatted as 'PHOTOL(IND)'.
        rate= photo_df['rate'][i]

        # If this photolysis rate does NOT have a correponding match to an MCM j-value in 'jdict, pop an error
        if rate not in list(jmap_dict.keys()):
            raise ValueError("In function, convert_FJX_rates_to_MCM(), we could not find a MCM J-Value \n" +
                             "match in 'jdict' corresponding to: '"+rate +
                             "'\n\nThis likely means you need to update 'jdict' to account for new MCM/FJX mapping j-value changes!")

        else:  # Otherwise, use the mapped value in J-Dict to define rate in MCM!

            # Define the j-value rate that should be used for this reaction!
            jind = fjx_df.index[fjx_df['PHOTOL(IND)'] == rate.replace(' ', '')]
            j_info = fjx_df.loc[jind[0], 'Info']
            rxn_df.at[i,'rates4foam']='k(:,i) = '+jmap_dict[rate]+';'+j_info
            rxn_df.at[i,'rate_function']=jmap_dict[rate].replace("'",'')
            
    return rxn_df

# %%---------------------------------------------------------------------------
# -------------- (A.2) parse_kpp_main() & Sub-Functions -----------------------
# -----------------------------------------------------------------------------
# def prep_df4f0am(rxn_df): 
    
#     rxns= rxn_df['rxn'].to_list()
#     cmts=rxn_df['comments'].to_list()
#     rates= rxn_df['rate'].to_list()

#     rxns4foam=["Rnames{i}='"+rxns[i]+ "'; % "+cmts[i] if len(cmts[i])>0 else "Rnames{i}='"+rxns[i]+ "';" for i in range(0,len(rxns))]
#     rates4foam=["k(:,i) = "+rates[i]+ ";" for i in range(0,len(rates))]
#     foam_df=pd.DataFrame({'Reactions':rxns4foam, 'Rates':rates4foam}
    
#     return foam_df
# %%
def prep4f0am(rxn_df, fjx_df, jmap_dict): 
    
    # Define a local function to create the reaction string for F0AM by combining 
    # the KPP parsed rxn & appending the KPP comments (if any). 
    def _fmt_foam_rxn(row):
        if len(row['comments']) > 0:
            return "Rnames{i} = " +f"'{row['rxn']}'; % {row['comments']}"
        else:
            return "Rnames{i} = " +f"'{row['rxn']}';"
        
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
    # Define a local function to create the reaction rate string for F0AM by 
    # using the KPP parsed rate & doing the photolysis mapping.  
    def _fmt_foam_rate(row):
        rate= row['rate']
        
        ###########################################################################
        # Pt 1: Transform Fortran Style Math OPERATORS into MATLAB operators: 
        ###########################################################################
        # Ensure any FORTRAN operators in rate degs are converted to MATLAB syntax ... 
        # In standard Fortran, the arithmetic operators like * (multiplication) and 
        # / (division) must be immediately adjacent to operands with no intervening 
        # spaces. So, here, the regex patterns used are look behind/ahead assertions 
        # to ensure the '*' or '/' is not preceded or followed by whitespace so we 
        # ONLY replace mathematical operators... 
        
        # Replace division operators '/' with './' for MATLAB
        rate = re.sub(r'(?<=\S)/(?=\S)', './', rate) 
        
        # Replace multiplication operators '*' with '.*' for MATLAB
        rate = re.sub(r'(?<=\S)\*(?=\S)', '.*', rate)
        
        if row['is_photo']==False:
            # Add 'Met' as an argument for all rate constant functions: 
            rate= _add_met_to_rates(rate)
                
            # For all non-photolysis rates, just format it for f0Am by adding k str & ';'. 
            rate4foam=f"k(:,i) = {rate};"
            func4foam=row[['rate_function']].to_list()
        else: 
            ###########################################################################
            #  MAP Photolysis Functions to Named Photolysis Rates: 
            ###########################################################################
            # Use the KPP rate parsed as the "key" to look up the j-value match in "jdict".
            # Both keys & rates should be formatted as 'PHOTOL(IND)'.
            
            # If this photolysis rate does NOT have a correponding match to a
            # j-value we KNOW is already defined in F0AM within 'jdict', pop an error.
            if rate not in list(jmap_dict.keys()):
                raise ValueError("In function, prep4f0am(), we could not find a J-Value defined in F0AM \n" +
                                 "within 'jdict' corresponding to: '"+rate +
                                 "'\n\nThis likely means you need to update 'jdict' to account for new F0AM/FJX mapping j-value changes!")
            else:  # Otherwise, use the mapped value in J-Dict to define rate in F0AM. 
                # Define the j-value rate that should be used for this reaction!
                jind = fjx_df.index[fjx_df['PHOTOL(IND)'] == rate.replace(' ', '')]
                j_info = fjx_df.loc[jind[0], 'Info']
                rate4foam=f"k(:,i) = {jmap_dict[rate]}; {j_info}"
                func4foam=jmap_dict[rate]
          
        return [rate4foam, func4foam]
    
    
    rinfo=rxn_df.apply(_fmt_foam_rate, axis=1)
    rates=[]; functs=[]; 
    for ind in rinfo.index: 
        rates.append(rinfo[ind][0])
        if type(rinfo[ind][1]) ==list: 
            functs.append(rinfo[ind][1][0])
        else: 
            functs.append(rinfo[ind][1])
            
    
    
    # Apply functions to columns
    foam_df = pd.DataFrame({
        'Reactions': rxn_df.apply(_fmt_foam_rxn, axis=1),
        'Rates': rates,
        'Rate_Functions': functs,
        'is_het': rxn_df['is_het'],
        'is_photo': rxn_df['is_photo'],
        'is_gas':rxn_df['is_gas']})

    return foam_df
# %%
# A.2.3.4 L3-function used within parse_rxns_and_rates() within parse_kpp_main() within make_GEOSCHEM_AllRxns_file()
def _check_list_lens(list_dict):
    """ Function to check the lengths of all input lists contained in the values of list_dict. 
    Will only pop an error if any list is not the same len as the others. 

    INPUTS: 
    -------
        list_dict - a Dictionary with list names as key and corresponding lists as values. 

    OUTPUTS: 
    -------
        NONE unless the lists are different lengths, in which case an error will arise. 

    REQUIREMENTS: 
    -------------
            LIBRARIES:           import sys
            CUSTOM FUNCTIONS:    None.    

    USAGE:    
    ------
        Used to ensure that all lists in dictionary "kpp_dict" that contain info 
        about reactions in a KPP file are all the same length to guard against 
        having dif lengths that might mess up indexing later... 

        CALLED WITHIN: parse_rxns_and_rates() within parse_kpp_main() within make_GEOSCHEM_AllRxns_file()

    AUTHOR: 
    -------
        Prof. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD

    CHANGE LOG: 
    ----------
    09/28/2023    JDH Created 
    10/31/2023    JDH Moved to make_GCAllRxns.py and updated documentation. 

    """

    # Get list of lengths of each list in list_dict (keys=varname, value= list)
    list_lens = [len(list_dict[list_i]) for list_i in list(list_dict.keys())]

    # If you have more than 1 unique value of a list length, then throw an error
    if len(np.unique(list_lens)) > 1:

        # Pull out the name of the function & line number where this error was found!
        # Name of parent function where error check initiated.
        err_function = sys._getframe(1).f_code.co_name
        # Line # of parent function error check initiated
        err_line = str(sys._getframe(1).f_lineno)

        # Make a string that will list all the lengths of the lists using their dict keys as names...
        err_display = '\n'.join(['    Len('+key+') = '+str(list_lens[i])
                                for i, key in enumerate(list(list_dict.keys()))])

        raise ValueError('Lists are not the same length in '+err_function+'() on Line # '
                         + err_line+': \n'+err_display).with_traceback(sys.exc_info()[2])

    return

# A.2.3.3 L3-function used within parse_rxns_and_rates() within parse_kpp_main() within make_GEOSCHEM_AllRxns_file()
def _check_kpp_dict(kpp_dict: dict, check_lens: bool = False, ignore_extras: bool = False):
    """Function to check if kpp_dict has all the right key/value pairs functions expect!
    Option to check the list lengths in kpp_dict as well... Returns nothing if "kpp_dict" is good!

    INPUTS: 
    -------
        (1) kpp_dict         - DICT with list names as key and corresponding lists as value. 

        (2) check_lens    - BOOL indicating wehther or not to check that the lists in the 
                            DICT, kpp_dict, stored as values are indeed the same length. 

        (3) ignore_extras - BOOL of whether you'd like to check kpp_dict but ignore any 
                            extra keys you weren't expecting... Useful if kpp_dict gets modified.'

    OUTPUT: 
        (1) NONE if kpp_dict is all good and lists are same length, otherwise pops an error. 

    REQUIREMENTS: 
    -------------
         LIBRARIES:           import sys

         CUSTOM FUNCTIONS:    check_list_lens() (defined above)

    USAGE: 
    -----  
    Used to ensure that "kpp_dict" has all the keys we expect everywhere its used 
    in this function. Optionally will check the lengths of lists in "kpp_dict". 
    Nice to ensure that I haven't messed up the keys in one function or another.

         CALLED WITHIN: parse_rxns_and_rates() within parse_kpp_main() 
                        within make_GEOSCHEM_AllRxns_file()

    AUTHOR: 
    -------
        Prof. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD
    """

    # Define a list of the expected keys in kpp_dict...
    expected_keys = ['KPP_line','rxn', 'reactants', 'rct_stoich', 'products', 'prd_stoich',
                     'is_gas',  'is_het', 'is_photo','is_3body', 'rate', 'rate_function',
                     'comments']

    # Figure out what keys in kpp_dict are extra or missing...
    missing = [x for x in expected_keys if x not in list(kpp_dict.keys())]
    extra = [x for x in list(kpp_dict.keys()) if x not in expected_keys]

    # If you have Missing and/or extra keys in your kpp_dict... pop an error!
    if any([True if len(listy) != 0 else False for listy in [missing, extra]]):

        # Pull out the name of the function & line number where this error was found!
        # Name of parent function where error check initiated.
        err_function = sys._getframe(1).f_code.co_name
        # Line # of parent function error check initiated
        err_line = str(sys._getframe(1).f_lineno)

        if len(missing) != 0 and len(extra) == 0:  # Only missing some keys...
            raise ValueError('In '+err_function+'() on Line # ' + err_line +
                             ', the dictionary, kpp_dict, is missing the following required keys: \n    ' +
                             ','.join(missing)+'\n').with_traceback(sys.exc_info()[2])
        # Only have extra keys...
        elif len(extra) != 0 and len(missing) == 0 and ignore_extras is False:
            raise ValueError('In '+err_function+'() on Line # ' + err_line +
                             ', the dictionary, kpp_dict, contains the following extra keys: \n    ' +
                             ','.join(extra)+'\n').with_traceback(sys.exc_info()[2])
        # Have missing keys and extra keys!
        elif len(missing) != 0 and len(extra) != 0 and ignore_extras is False:
            raise ValueError('In '+err_function+'() on Line # ' + err_line +
                             ', the dictionary, kpp_dict, is missing the following required keys: \n    ' +
                             ','.join(missing)+'\n     and has the following extra keys: \n    ' +
                             ','.join(extra)+'\n').with_traceback(sys.exc_info()[2])

    # Optionally, check if all the lengths of lists in kpp_dict are consistent!
    if check_lens is True:
        _check_list_lens(kpp_dict)

    return

# A.2.3.2 L3-function used within parse_rxns_and_rates() within parse_kpp_main() within make_GEOSCHEM_AllRxns_file()
def convert_fortran_nums(in_string):
    """ Function to take an input string and convert Fortran-style numeric literals 
        and (select) common matehmatical function names within it into modern 
        scientific notation syntax. Used to convert numbers in KPP files to something you 
        can write to a MATLAB/python/CSV/excel file. Used to convert strings with 
        reaction rates like:
            
                     In KPP File:                        as Output 
            'GCARR_ac(1.00d-14, -490.0d0)'   to   'GCARR_ac(1.00e-14, -490.0)'

        Can handle: 
        - Fortran exponents like 'd', 'D', 'e', 'E' (e.g., '1.23d-4', '5.67D+8')
        - Removes the '_dp' suffix used in Fortran for double precision #s
        - Converts function names EXP, LOG, LOG10 to lowercase. 

    INPUTS: 
    ------  
        (1) in_string - single STR containing (maybe) some numerical FORTRAN notation 

    OUTPUTS: 
    --------
        (2) out_string - updated STR o with FORTRAN math converted to (modern) math.

    USAGE: 
    ------ 
        Called Within: parse_rxns_and_rates() from within parse_KPP_main() from 
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
    # part, which starts with either'd', 'D', 'e', or 'E', followed by an optional
    # '+' or '-', and then digits (for example, d-4 or E+10).Finally, it can 
    # optionally end with the suffix '_dp', which is used in Fortran to indicate 
    # double precision numbers. After finding these, string subs are only done 
    # when the number matched has a dif format (eg. doesn't mess with integers).
    # Example matches &(subsequent) substitions below: 
    #      '42'        -->  '42'       (intergers are same).
    #      '-5.67d-8'  -->  '-5.67e-8'
    #      '3.29E+00'  -->  '3.29'
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

    ###########################################################################
    #  Part II: Transform Fortran Style Math FUNCTIONS into lowercase functions: 
    ###########################################################################
    # Loop over some common math function names / convert them to lowercase 
    for fname, funct_pattern in [('EXP', r'EXP\((.*?)\)'),   # Pattern to match 'EXP(' ... ')'
                               ('LOG', r'LOG\((.*?)\)'),      # Pattern to match 'LOG(' ... ')'
                               ('LOG10', r'LOG10\((.*?)\)')]: # Pattern to match 'LOG10(' ... ')'
    # Use the re.sub() function from the regex module to: 
    #    (1) Find all substrings in the input string that match the fortran function pattern. 
    #    (2) For each match, call the lambda function, passing in the match object to 
    #        generate a new replacement substring. It converts the function name to 
    #        lowercase using 'fname.lower()' while maintaining the original arguments 
    #        to the function with 'm.group(1)'.
    #    (3) Replace each matched substring with the new substring returned by the lambda function. 
    #    (4) Assign the result to a single output string after all substring replacements. 
    # This results in replacing, for example, 'EXP(1.23)' with 'exp(1.23)', or 'LOG10(4.56)' with 'log10(4.56)'

        out_string = re.sub(funct_pattern,lambda m: f"{fname.lower()}({m.group(1)})",out_string)
    
    return out_string

# A.2.3.1 L3-function used within parse_rxns_and_rates() within parse_kpp_main() within make_GEOSCHEM_AllRxns_file()
def _sep_stoich_vs_tracer(rxn_str: str, seps: list = ['+', '=', '->', '-->', '<->']):
    """Function used to seperate all stoichiometry from compounds in a string into 2 lists
        using custom string seperators to accomodate dif reaction writing styles. 
        NOTE: You do NOT have to pass a whole rxn- you can use recusively after splitting 
        prods from rcts, but function won't work properly if chemical compounds names
        can begin with numbers. 

    INPUTS:
    -------
        rxn_str -  String with reaction like: 'TOLU + OH = TRO2 + 1.920CH2O + 0.260GLYX + 0.215MGLY + OH' 

        seps    -  (OPTIONAL) list of strings to use to seperate reaction. Default is to 
                    use all of these to seperate rxns: ['+','=', '->', '-->', '<->']

    OUTPUTS: 
    --------
        cmpds  - A string list containing all the compounds in the reaction passed. 
                 For example above, would look like: ['TOLU', 'OH', 'TRO2', 'CH2O', 'GLYX', 'MGLY', 'OH']

        stoich- A list of floating point  numbers containing the stoichiometric coefficents 
                that correspond to the compounds (identically index). 
                For example above, would look like: [1, 1, 1, 1.920, 0.260, 0.215, 1]

    REQUIREMENTS: 
    -------------
        LIBRARIES: None

        CUSTOM FUNCTIONS REFERENCED:  None 

    USAGE: 
    ------
        Called Within:  parse_rxns_and_rates() within parse_kpp_main() within make_GEOSCHEM_AllRxns_file() 

        Output Usage:  Placed in output dict 'kpp_dict' of parse_rxns_and_rates() 
                       within parse_kpp_main() within make_GEOSCHEM_AllRxns_file() 

    AUTHOR: 
    -------
        Prof. Jessica D. Haskins (jessica.haskins@autah.edu) GitHub: @jhaskinsPhD

    CHANGE LOG: 
    ----------
    01/18/2022    JDH Created for pyF0AM
    09/28/2023    JDH Copy/pasted from pyF0AM to GC_Emulator/parse_KPP, updated documentation.
    """

    ls = []
    cmpds = []
    stoich = []
    for sep in seps:  # Split input string at the input seperators:
        # ls =['TOLU','OH','TRO2','1.920CH2O','0.260GLYX','0.215MGLY','OH']
        ls = ls+rxn_str.replace(' ', '').split(sep)

    for grp in ls:  # Loop over all groupings in the list split at seps
        if len(grp) > 0:
            # If str begins with a letter, then its a tracer name (e.g. 'TOLU').
            if grp[0].isalpha():
                cmpds.append(grp)
                stoich.append(1)  # Save the compound, save the stoichiometry.

            else:  # Otherwise, loop over the chars in the grp until you find a letter.
                yld = grp[0]
                ind = 1
                # Set stoichiometry= to first char (not a letter!)
                found_yld = False

                while found_yld == False:  # loop til you find a letter.
                    # Found a letter or decimal place.
                    if grp[ind].isnumeric() or grp[ind] == '.':
                        # Add to str containing all stoichiometry.
                        yld = yld+grp[ind]
                    # Found beginning of compound name.
                    elif grp[ind].isalpha() or ind == len(grp)-1:
                        cmpds.append(grp[ind:])  # Store it
                        # and the #s making up the yield as a #
                        stoich.append(float(yld))
                        found_yld = True
                        break
                    ind = ind+1

    # Don't let 'hv' be a "compound" or have "stoichimetry"...
    if 'hv' in cmpds:
        stoich.pop(cmpds.index('hv'))
        cmpds.remove('hv')

    return cmpds, stoich

# A.2.3 L2-function used in parse_kpp_main() called within make_GEOSCHEM_AllRxns_file():
def _is_blank(line):
    """Function to determine if a line is blank or not after all spaces have been removed."""
    tf = True if len(line.replace(' ', '')) == 0 else False
    return tf

# A.2.3 L2-function used in parse_kpp_main() called within make_GEOSCHEM_AllRxns_file():
def parse_rxns_and_rates(line: str, kpp_dict: dict, file_location: dict):
    """Function used to parse reaction lines in KPP file and extract their info into dictionary, kpp_dict. 

    INPUTS:
    -------

        (1) line - STR of line with reaction, rate, and (maybe) a comment. These lines look like this:     
                      Reaction               Rate                              Comment
                "O3 + NO = NO2 + O2 :  GCARR_ac(3.00d-12, -1500.0d0); {2014/02/03; Eastham2014; SDE}"

        (2)  kpp_dict - DICT with the following Key/Value Pairs.
                KEYS           :           VALUES
                ----------------------------------------------------------------------
                KPP_line       : List of original line from KPP file
                rxn            : List of strings with full reaction only... 
                reactants      : *Nested* list containing all reactants per reaction.  
                rct_stoich     : *Nested* list of reactants stoichiometry per reaction.  
                products       : *Nested* list containing all products per reaction. 
                prd_stoich     : *Nested* list of products stoichiometry per reaction. 
                is_gas         : Boolean List of whether rxn is gas phase or not. 
                is_het         : Boolean List of whether rxn is heterogeneous or not. 
                is_photo       : Boolean list of whether rxn is a photolysis rxn or not. 
                is_3body       : Boolean list of whether rxn includes 3rd body reactant, M, or not. 
                rate           : List of strings with full rate for each reaction 
                rate_function  : List of strings with rate function for each reaction
                comments       : List of comments for each reaction 

                NOTE: All lists in the dict, kpp_dict, should be the same length upon input and output to 
                    this function, as they all are all "indexed "by the reaction they refer to. 

        (3) file_location - Dictionary with the following Key/Value Pairs: 

                KEYS           :           VALUES
                ----------------------------------------------------------------------
                in_header      : Boolean. True if currently parsing the header, otherwise False.   
                in_defvar      : Boolean. True if currently parsing the fortran variable defintion section
                in_deffix      : Boolean. True if currently parsing the "#deffix" section  
                in_sulf_rxns   : Boolean. True if currently parsing sulfate mod het reactions  
                in_gas_rxns    : Boolean. True if currently parsing Gas Phase reactions
                in_het_rxns    : Boolean. True if currently parsing Heterogeneous reactions
                in_photo_rxns  : Boolean. True if currently parsing Photolysis reactions

     OUTPUTS: 
     -------   
         (1)  kpp_dict       - Updated DICT with the same Key/Value Pairs as above, but outputted 
                             lists in kpp_dict will be +1 in len since they will now include info 
                             about the reaction contained on the input line. This function checks 
                             for len consistency @ end. 

        (2) file_location - Updated DICT with the same Key/Value Pairs as above. 

    REQUIREMENTS: 
    -------------
        LIBRARIES:           None

        CUSTOM FUNCTIONS:    check_kpp_dict()
                            _str_multi_replace()
                            sep_stoich_vs_tracer()
                            fortran2matlab_nums()
                            check_list_lens()

    USAGE: 
    -----  
        Called Within: parse_kpp_main() within make_GEOSCHEM_AllRxns_file(). 
        Output Usage: Output dict 'kpp_dict' used in parse_kpp_main() to sep info for dif rxn types into own dicts! 

    AUTHOR: 
    -------
        Prof. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD

    """
    # Check that the input dictionary, kpp_dict, has the required input keys & all lists are consistent lengths
    _check_kpp_dict(kpp_dict, check_lens=True)

    # Clean up line to maintain in output: 
    clean_line= utils.str_multi_replace(line, ['\n', '\t'], rep_all='').strip()
    
    #####################################################################################
    # Split Line in to reaction, rate, rate function, and comment & Clean up for MATLAB!
    #####################################################################################
    rxn_parts = line.split(':')  # Split on colon to get rxn only
    rate = rxn_parts[1].split(';')[0].strip() # Split on ';' to get rates, remove trailing/tailing spaces
        
    # Define Booleans that tell us what type of rxn it is: 
    rxn_with_M=True if '{+M}' in rxn_parts[0] else False  
    rxn_is_gas= file_location['in_gas_rxns']
    rxn_is_photo= file_location['in_photo_rxns']
    rxn_is_het=True if any ([file_location['in_het_rxns'],file_location['in_sulf_rxns']]) else False 

    # If this part of the line after the reaction has a comment, then split that on'{' and clean it up.
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
    rct_i, rct_stoich_i = _sep_stoich_vs_tracer(rxn.split('=')[0], seps=['+'])
    prd_i, prd_stoich_i = _sep_stoich_vs_tracer(rxn.split('=')[1], seps=['+'])

    # Clean up rates so the #s no longer have FORTRAN style numbers (does NOT
    # change operators), only #s and case ofcommon math functions: 
    rate=convert_fortran_nums(rate)

    # Get function name only from the rate constant defs/ The lines look like this: 
    # 'GCARR_ac(2.60e-13, 1300)' but some maybe sums of calls to multiple functions 
    # like so: GCARR_ac(2.60e-13, 1300) + GCARR_ac(2.60e-13, 1300) and some may 
    # contain exp(), log(), or even log10()... So use cray regex: 
    pattern = r'\b((?!exp\b|log\b|log10\b)[A-Za-z_]+[a-zA-Z0-9_]*)\('

    # Use a list comprehension to apply the regex to each string in rates
    functions = [match_i[:] for match_i in re.findall(pattern, rate)]     
        
    # ------------------------------------------------------------------------- 
    # Now save everything into its appropriate output list: 
    # ------------------------------------------------------------------------- 
    # Original KPP line(s) - minus any new line chars. 
    kpp_dict['KPP_line'].append(clean_line)
    # List of strings with full reaction only.
    kpp_dict['rxn'].append(rxn)
    # List of bools indicating reaction type of each reaction
    kpp_dict['is_gas'].append(rxn_is_gas)
    kpp_dict['is_het'].append(rxn_is_het)
    kpp_dict['is_photo'].append(rxn_is_photo)
    kpp_dict['is_3body'].append(rxn_with_M)  
    # *Nested* list containing all reactants for each reaction
    kpp_dict['reactants'].append(rct_i)
    # *Nested* list of reactants stoichiometry for each reaction
    kpp_dict['rct_stoich'].append(rct_stoich_i)
    # *Nested* list containing all products for each reaction
    kpp_dict['products'].append(prd_i)
    # *Nested* list of products stoichiometry for each reaction
    kpp_dict['prd_stoich'].append(prd_stoich_i)
    # List of strings with full rate for each reaction
    kpp_dict['rate'].append(rate)
    # *Nested* List of rate function(s) used for each reaction
    kpp_dict['rate_function'].append(functions)
    # List of comments for each reaction
    kpp_dict['comments'].append(comment)

    # Check that all new lists in kpp_dict are still the same length!
    _check_list_lens(kpp_dict)

    return kpp_dict

# A.2.2 L2-function used in parse_kpp_main() called within make_GEOSCHEM_AllRxns_file():
def fix_multiline_rxns(line, line_stripped, count, kppfile):
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
        CALLED WITHIN:  parse_kpp_main() via a call to make_GEOSCHEM_AllRxns_file() 

        OUTPUT USAGE: Updated output "line" next passed to parse_rxns_and_rates() within 
                        parse_kpp_main() called within make_GEOSCHEM_AllRxns_file(). 

    AUTHOR: 
    -------
        Prof. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD

    CHANGE LOG: 
    ----------
        09/28/2023    JDH Created 
        11/15/2024: JDH added "':' not in line_stripped: " to while statement 
                    to work with Alfie's custom version of 14.5.0... 

    """

    # Read more lines til you get the entire reaction, rate & comments in a single line!
    while line_stripped.endswith('+') == True or ':' not in line_stripped:
        nline = kppfile.readline()   # Read next line from file
        count = count+1                 # Update line counter variable
        # append this new line to our current line str.
        line = line + nline.strip()
        # Update line_stripped to not contain spaces!
        line_stripped = line.strip().replace(' ', '')

    # Now that you got the entire reaction on this line, make sure to replace any new line characters with a space!
    line = line.replace('\n', ' ')
    line_stripped = line_stripped.replace('\n', ' ')

    # Return the count, new line, and line_stripped vars to continue on.
    return line, line_stripped, count

# A.2.1 L2-function used in parse_kpp_main() called within make_GEOSCHEM_AllRxns_file():
def parse_var_defs(line: str, tracers: list = [], long_name: list = [], split_on: str = 'ignore'):
    """Function to parse variable defintion lines in the KPP File header. These lines look like: 
        'A3O2       = IGNORE; {CH3CH2CH2OO; Primary RO2 from C3H8}'  

    INPUTS: 
    -------
        (1) line      - STR of the current KPP file line you're parsing for info. 

        (2) tracers   - LIST containing names of all tracers that have already been defined 
                         by parsing previous lines. Only contains tracer abbreviation! 

        (3) long_name - LIST containing extra comment info about the tracer being defined that 
                         has been compiled by parsing previous lines. 

        (4) split_on  - (OPTIONAL) STR we want to split the input line on... 
                        It is case insensitive & space insensitive. Default is set to 
                        'ignore' to work with most recent GC KPP files. 

    OUTPUTS:
    -------- 
        tracers   - Updated LIST with this line's new defined variable in it.

        long_name - Updated LIST with this line's tracer info in it. 

    REQUIREMENTS: 
     -------------
        LIBRARIES:           None

        CUSTOM FUNCTIONS:    _str_multi_replace()

    USAGE: 
    -----  
        CALLED WITHIN:  parse_kpp_main() via a call within make_GEOSCHEM_AllRxns_file(). 

        OUTPUT USAGE:   Zipped into dict "var_info" passed as output of make_GEOSCHEM_AllRxns_file()

    AUTHOR: 
    -------
        Prof. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD

    CHANGE LOG: 
    ----------
        09/28/2023    JDH Created 
        10/30/2023    JDH updated documentation. 

    """

    if split_on in line.lower():
        # Figure out index of where to split line between def and comments. Case insensitive.
        ind = line.lower().index(split_on)

        # Remove spaces and "=" from line, to get only the tracer name being defined!
        tracer = utils.str_multi_replace(line[:ind],  [' ', '='], rep_all='')

        # Keep tracer name of this compound defined in a list.
        tracers.append(tracer)

        # Pull out the rest of the comment about this tracer using the index.
        tracer_info = line[ind+len(split_on):]

        # Remove any remaining weird characters we don't care about (";", "{", "}").
        if tracer_info[0] == ';':
            tracer_info = tracer_info[1:]
        tracer_info = utils.str_multi_replace(
            tracer_info,  ['{', '}'], rep_all='').strip()

        # Now add the "long name" slash comment about this tracer to a list!
        long_name.append(tracer_info)

    return tracers, long_name

# A.2 L1-function called within make_GEOSCHEM_AllRxns_file()
def parse_kpp_main(kppfile, out2yaml:bool=True, out2excel:bool=False, debug:bool=False):
    """Function to parse a KPP file and extract info about tracers, reactions,
    & rates in a pandas dataframe. Output can also be provided in other ways using 
    keywords. 

    INPUTS:  
    -------  
       (1) kppfile - STR containing the full path and name of a KPP file.  

       (2) out2yaml - BOOL indicating if outputs should be saved to a yaml file. 
       
       (3) out2excel - BOOL indicating if outputs shoudl be saved to an excel .xlsx sheet. 
              
       (4) debug- BOOL for printing additional output for debugging. Default is FALSE. 
       
    OUTPUTS:  
    -------- 
       (2) var_info - 
       
       (2) kpp_df - 

    REQUIREMENTS: 
    -------------
        LIBRARIES:           None

        CUSTOM FUNCTIONS:   parse_var_defs()
                            _is_blank()
                            fix_multiline_rxns()
                            parse_rxns_and_rates()
                            check_kpp_dict()

    USAGE: 
    -----  
        CALLED WITHIN: 

        OUTPUT USAGE: 

    AUTHOR: 
    -------
       Prof. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD
    """
    ###########################################################################
    # Initialize all vars to hold outputs/ keep track of where we are!
    ###########################################################################

    # Get # of lines to expect so we know when we've reached the end of the file!
    num_lines = sum(1 for line in open(kppfile))
    tracers = []
    long_name = []  # Make empty lists to hold tracers & tracer info!
    count = 0  # Initialize line counter variable.

    # Define boolean dictionary with key/values that will tell us what part of the
    # KPP file we are currently parsing.
    file_location = dict({
        'in_header': True, # Boolean. True on initiation & if currently parsing header.
        'in_defvar': False,# Boolean. True if currently parsing the Fortran variable definition section
        'in_deffix': False,  # Boolean. True if currently parsing the "#deffix" section
        'in_sulf_rxns': False,  # Boolean. True if currently parsing sulfate mod het reactions
        'in_gas_rxns': False,  # Boolean. True if currently parsing gas phase reactions
        'in_het_rxns': False,  # Boolean. True if currently parsing heterogenous reactions
        'in_photo_rxns': False  # Boolean. True if currently parsing photolysis reactions
    })

    # Define dictionary to hold lists of info about dif rxns!
    kpp_dict = dict({
        # List of original KPP file lines.
        'KPP_line': list([]),      
        # List of strings with full reaction only.
        'rxn': list([]),           
        # Boolean list of whether rxn is gas phase or not. 
        'is_gas': list([]),
        # Boolean List of whether rxn is heterogeneous or not.
        'is_het': list([]),
        # Boolean list of whether rxn is a photolysis rxn or not.
        'is_photo': list([]),
        # Boolean list of whether rxn is a 3-bdy rxn w/ M or not.
        'is_3body': list([]),
        # NESTED list containing all reactants for each reaction
        'reactants': list([]),
        # NESTED list of reactants stoichiometry for each reaction
        'rct_stoich': list([]),
        # NESTED list containing all products for each reaction
        'products': list([]),
        # NESTED list of products stoichiometry for each reaction
        'prd_stoich': list([]),
        # List of strings with full rate for each reaction
        'rate': list([]),
        # List of strings with rate function for each reaction
        'rate_function': list([]),
        # List of comments for each reaction
        'comments': list([])       
    })

    ##############################################################################
    # Open the KPP file in "read" mode & Loop over all lines until end is reached!
    ##############################################################################
    file = open(kppfile, 'r')

    while True:  # loop over each line of the KPP file.
    
        line = file.readline()  # Read next line from file
        
        # Strip all white space from current line.
        line_stripped = line.strip().replace(' ', '')
        
        # Update line counter variable now that you read in this line.
        count = count+1

        # --------------------------------------------------------------------------------------------------------------------------
        # BEFORE parsing the line, set file location booleans to FALSE if the current line indicates you're entering a new section!
        # --------------------------------------------------------------------------------------------------------------------------
        # You are no longer in the header section once #defvar appears on the line!
        file_location['in_header'] = False if '#defvar' in line_stripped.lower() else file_location['in_header']
        # You are no longer in the variable defintion section if #deffix appears in the current line!
        file_location['in_defvar'] = False if '#deffix' in line_stripped.lower() else file_location['in_defvar']
        # You are no longer in section defining fixed vars if #equations appears in the current line!
        file_location['in_deffix'] = False if '#equations' in line_stripped.lower() else file_location['in_deffix']
        # You are no longer in the sulfate_mod reaction section on the next line if 
        # all of these strs: '//', 'gas', 'phase', 'chemistry', & 'reactions' appear in the current line.
        file_location['in_sulf_rxns'] = False if all([item_i in line_stripped.lower() \
                                        for item_i in ['//', 'gas', 'phase', 'chemistry', 'reactions']]) \
                                        else file_location['in_sulf_rxns']
        # You are no longer in the gas phase reaction section on the next line if 
        # all of these strs: '//', 'heterogeneous','chemistry', & 'reactions' appear in the current line!
        file_location['in_gas_rxns'] = False if all([item_i in line_stripped.lower() \
                                       for item_i in ['//', 'heterogeneous', 'chemistry', 'reactions']]) \
                                       else file_location['in_gas_rxns']
        # You are no longer in the photolysis reaction section n the next line if 
        # all of these strs: '//', 'photolysis', & 'reactions' appear in the current line!
        file_location['in_het_rxns'] = False if all([item_i in line_stripped.lower() \
                                       for item_i in ['//', 'photolysis', 'reactions']]) \
                                       else file_location['in_het_rxns']
        
        # Get a list of the keys that are currently set to "TRUE" (useful for debugging).. 
        true_keys = '.'.join([key for key, value in file_location.items() if value])
        if debug == True: 
            print(f'Line: "{line}"\n\tFile Location: {true_keys} \n')

        # Break if reached end of file, don't attempt to parse line.
        if count > num_lines:
            break
        # ------------------------------------------------------------------------------------------
        # If you're in the part of the file where variable definitions are made  //
        # tracer names are listed then pass to ==> parse_var_defs()
        # ------------------------------------------------------------------------------------------
        if file_location['in_defvar'] == True or file_location['in_deffix'] == True:

            # Pass to parse_var_defs() & return updated lists with tracer name and long name.
            tracers, long_name = parse_var_defs(line, tracers=tracers, long_name=long_name, split_on='ignore')

        # -----------------------------------------------------------------------------------------
        # If you're in part of the file with chemical equations & rates ==> parse_rxns_and_rates()
        # -----------------------------------------------------------------------------------------
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

                    line, line_stripped, count = fix_multiline_rxns(line, line_stripped, count, file)

                # (After fixing if necessary), parse this reaction and add all info about this rxn to lists:
                kpp_dict = parse_rxns_and_rates(line, kpp_dict, file_location)

        # --------------------------------------------------------------------------------------------------------------------------
        # AFTER parsing the line, set file location booleans to TRUE if the current line indicates you're entering a new section!
        # --------------------------------------------------------------------------------------------------------------------------
        # You will be in the variable definition section on the next line if '#defvar' appears in the current line.
        file_location['in_defvar'] = True if '#defvar' in line_stripped.lower() else file_location['in_defvar']
        # You will be in the deffix section on the next line if #deffix appears in the current line.
        file_location['in_deffix'] = True if '#deffix' in line_stripped.lower() else file_location['in_deffix']
        # You will be in the sulfate_mod reaction section on the next line if
        # all these strings appear in the current line.
        file_location['in_sulf_rxns'] = True if all([item_i in line_stripped.lower() \
                                        for item_i in ['//', 'sulfate_mod', 'reactions']]) \
                                        else file_location['in_sulf_rxns']
        # You will be in the gas-phase reaction section on the next line if 
        # all these strings appear in the current line.
        file_location['in_gas_rxns'] = True if all([item_i in line_stripped.lower() \
                                            for item_i in ['//', 'gas', 'phase', 'chemistry', 'reactions']]) \
                                            else file_location['in_gas_rxns']
        # You will be in the heterogenous reaction section on the next line if 
        #all these strings appear in the current line.
        file_location['in_het_rxns'] = True if all([item_i in line_stripped.lower() \
                                            for item_i in ['//', 'heterogeneous', 'chemistry', 'reactions']]) \
                                            else file_location['in_het_rxns']
        # You will be in the photolysis reaction section on the next line if 
        # all these strings appear in the current line.
        file_location['in_photo_rxns'] = True if all([item_i in line_stripped.lower(
        ) for item_i in ['//', 'photolysis', 'reactions']]) else file_location['in_photo_rxns']

    # While loop will break once we found the end of the KPP file and have parsed all lines.
    file.close()  # Close the KPP file.

    # Build var_info, a dictionary connecting tracer names & long names.
    var_info = dict(zip(tracers, long_name))
    
    # --------------------------------------------------------------------------
    # Seperate out reactions from kpp_dict into dif dicts for each type of rxn.
    # --------------------------------------------------------------------------
    # First, just make sure kpp_dict has right keys & consistent list lengths!
    _check_kpp_dict(kpp_dict, check_lens=True, ignore_extras=False)

    # Make dataframe with all this information... 
    kpp_df = pd.DataFrame(kpp_dict)
    
    if out2yaml==True: 
        yaml_file=utils.check_filename(filename='kpp_mech', default_name= 'kpp_mech', ext='.yml', 
                           savepath=master_outdir, overwrite=True, return_full=True, quiet=True)
        # Write the dictionary to a YAML file
        with open(yaml_file, 'w') as f:   
            yaml.dump(kpp_dict, f)
        print(f'Output of () saved at: \n\t{yaml_file}')
    
    if out2excel==True: 
        excel_file=utils.check_filename(filename='kpp_mech', default_name= 'kpp_mech', ext='.xlsx', 
                           savepath=master_outdir, overwrite=True, return_full=True, quiet=True)
        # Write the dataframe to an excel file: 
        kpp_df.to_excel(excel_file)
        print(f'Output of parse_kpp() saved at: \n\t{excel_file}')
            
    return var_info, kpp_df, kpp_dict 


# %%---------------------------------------------------------------------------
# -------------- (A.1) create_jmap_dict() & Sub-Functions ---------------------
# -----------------------------------------------------------------------------

def make_unique_Jname(jname, current_names):
    """Function to make a new , unique nickname for a J-value currently named 'jname' 
    that does NOT appear in the input list 'current_names'

    INPUTS:
    ------
        (1) jname         - STR of current "nickname" of a jvalue (e.g. 'jNO3')
        (2) current_names - LIST of strs with all nicknames of UNIQUE j-values that are already "taken". 

    OUTPUTS: 
    ------- 
        (1) jname - STR of unique "nickname" of a jvalue that does NOT apprear in the list "current_names"

    REQUIREMENTS: 
    ------------ 
        LIBRARIES:           import re
        CUSTOM FUNCTIONS:    None.

    USAGE: 
    ------ 
        CALLED WITHIN: create_jmap_dict() within make_GEOSCHEM_AllRxns_file() within make_GC_mechanism()
        OUTPUT USAGE: Used to assign "nicknames" of J-values in GEOSChem_J.m and GEOSChem_AllRxns.m 

    AUTHOR: 
    ------
        Prof. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD

    CHANGE LOG: 
    --------- 
        11/06/23 - JDH Created. 
    """

    ct = 1  # Initialize counter variable

    while True:
        if jname in current_names:
            # Pull the "name" of the thing out. For "jCH3COOH", nm_only='CH3COOH'
            nm_only = re.findall(r'j([a-zA-Z0-9]+)', jname)
            if len(nm_only) == 0:
                print(jname)
                sys.exit()
            # Decide on a new suffix to try using (iterates through alphabet).
            new_let = '_'+chr(ord('`')+ct)
            jname = 'j'+nm_only[0]+new_let

        if jname not in current_names:
            break
        else:
            ct = ct+1  # Update counter var to try a new '_'+LETTER(ct)

    return jname

# A.1.1 L2-function used in create_jmap_dict() called within make_GEOSCHEM_AllRxns_file():
def read_FJX_file(file):
    """ Function to read in a FJX fj2j.dat file as a pandas dataframe. 

    INPUTS: 
    -------

       (1) file - full path to the FJX fj2j.dat file you'd like to read in. 

    OUTPUTS: 
    --------
       (1)  fjx_df - pandas dataframe with info from the FJX file read in. Includes these columns: 

                'FJX_Index'     -   INTEGER internal index used by Fast-JX, & the number 
                                    used to index PHOTOL(), the J-value rate function, in the KPP file.  

                'Quantum_Yield' -   FLOAT of the quantum yield of the reaction... A flat multiplier 
                                    applied to the first-order rate returned by Fast-JX for this reaction only. 

                'FJX_CrossSec'   -  STR of the cross-section (from FJX_spec.dat) which will be used 
                                    to calculate the first-order reaction rate. NOTE: Sometimes we 
                                    use a "similar" cross-section to define a photolysis rate we don't 
                                    have exact info about.

                'Photolysis_Rxn' -  STR of the actual photolysis reaction this is for... 

                'J_NickName'     -  Nickname of J-value... 

    REQUIREMENTS: 
    -------------
        LIBRARIES:                     import pandas as pd 

        CUSTOM FUNCTIONS REFERENCED:   None

    USAGE: 
    -----  
        CALLED WITHIN: create_jmap_dict() which is called within make_GEOSCHEM_AllRxns_file()

        OUTPUT USAGE: Output of 'fjx' containing info from the raw FJX file is used to build 
                        the output of create_jmap)dict() 'jdict' that gest passed to 
                        KPP_Dict_to_F0AM_IO(..., is_photo=True )in make_GEOSCHEM_AllRxns_file(). 

    AUTHOR: 
    -------
        Prof. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD

    CHANGE LOG: 
    ----------
        01/14/2022    JDH Created   
        10/03/2023    JDH modified documentation, moved FJX handling functions around... 

    """
    # Use pandas to read in the file, skipping the 1st two rows and with dif widths of file columns.
    fjx_df = pd.read_fwf(file, skip_rows=[1, 2], widths=[5, 10, 10, 34, 6, 8])

    fjx_df = fjx_df.drop(index=len(fjx_df)-1)  # Drop nan column.

    names = list(fjx_df.columns)  # Get original column names.

    # Rename all the columns to more intuitive names.
    fjx_df = fjx_df.rename(columns={names[0]: 'FJX_Index',     # See above
                                    # Name of the GEOS-Chem species that will undergo photolysis.
                                    names[1]: 'Photolysis_Of',
                                    # Indicator that the reaction is a photolysis rxn.
                                    names[2]: 'PHOTON',
                                    # Products of this specific photolysis reaction ...
                                    names[3]: 'Into',
                                    names[4]: 'Quantum_Yield',  # See above
                                    names[5]: 'FJX_CrossSec'})  # See above

    # Join parts of the columns together to make a coherant photolysis reaction
    fjx_df['Photolysis_Rxn'] = [str(fjx_df.loc[ind, 'Photolysis_Of']) + '+hv -> '+'+'.join([prd.replace('...', '')
                                                                                            for prd in str(fjx_df.loc[ind, 'Into']).split(' ') if len(prd) > 0]) for ind in range(0, len(fjx_df))]

    # Clean up cross section used column from spaces and other characters we don't need!
    fjx_df['FJX_CrossSec'] = fjx_df['FJX_CrossSec'].apply(lambda x: x.replace(
        ' ', '').replace('/', '').replace('(', '').replace(')', '').replace('-', ''))

    # Stick a j in front of the cross section used to define a "nickname"...
    fjx_df['J-NickName'] = fjx_df['FJX_CrossSec'].apply(lambda x: 'j'+str(x))

    # Now drop columns that we don't actually need.
    fjx_df = fjx_df.drop(columns=['Photolysis_Of', 'PHOTON', 'Into'])

    return fjx_df

# A.1 L1-function called within make_GEOSCHEM_AllRxns_file()
def create_jmap_dict(version: str, jmap_type='CrossSec_Match'):
    """Function that takes a GEOS-Chem version number, and reads in the appropriate 
    FJX file. Then it uses the PHOTOL(IND) and cross-section name used to 
    match the GEOS-Chem J-values to an analogous J-value that is already defined in F0AM.

    NOTE: This means that the photolysis values used between GEOS-Chem Classic and 
    in this GEOS-Chem Emulator are not the same. We are simply mapping the j-values 
    based on how they are set in GEOS-Chem classic to their closest F0AM alternative. 

    Inputs: 
    -------
       (1) version - string formatted as '13.3.4' corresponding to what GEOS-Chem 
                 version you are seeking J-Value matches to.

       (2) jmap_type - string corresponding to what mapping data  in the file, 
                 'FJX_cross_sect_to_jvalues.xlsx' you would like to use. 
                 Default is 'CrossSec_Match', which maps j-values exactly as GEOS-Chem does. 

                 Other options include...
                     'MCM', which maps the J-Values to their BEST match in the MCM, 
                            regardless of whether FJX & GEOS-Chem use that cross-section 

                     'RCIM', which maps J-values as the RCIM-F0AM version does, 
                             which uses more realistic J-Values that are pre-defined 
                             in F0AM, and does differ from those defined in the MCM. 

                 NOTE: Both 'RCIM' and 'MCM' options will not replicate perfectly GEOS-Chem's 
                 chemical mechanism, as they contain more detail than it woul, which is why 
                 the DEFAULT is to do the 'CrossSec_Match' match based off the cross-section used! 

    Outputs: 
    -------

        (1) jdict - A dictionary that has GEOS-Chem's PHOTOL(IND) as keys, and their 
                matched, corresponding J value in the MCM as its values. 
                This is passed as input to GEOSChem_Emulator.make_GEOSCHEM_AllRxns_file() to turn all 
                references of PHOTOL(IND) into the right J-value. 

    REQUIREMENTS: 
    -------------
        LIBRARIES:           import pandas as pd 
                             import os
                             import sys 
                             import numpy as np 

        CUSTOM FUNCTIONS:    read_FJX_file()

    USAGE: 
    -----  
        CALLED WITHIN: make_GCAllRxns_file()

        OUTPUT USAGE: Passed as input to KPP_Dict_to_F0AM_IO(...,is_photo=True, jdict=THIS_OUTPUT)
                     within make_GCAllRxns_file()  

    AUTHOR: 
    -------
        Dr. Jessica D. Haskins (jessica.haskins@utah.edu)   GitHub: @jhaskinsPhD 

    CHANGE LOG: 
    ----------
    8/10/2022: JDH Created 
    10/3/2023: JDH Modified documentation, moved FJX stuff to own file. 
    """
    # Retrieve the path of where the photolysis data live:
    photolysis_paths = utils.get_my_relative_paths(photo_only=True)

    # Read in dictionary defining what version of FJX files are used as CHEM_INPUTS for each dif GC version w/ gc version as key
    vfjx = pd.read_excel(
        photolysis_paths['FJX_input_by_GC-version.xlsx'], index_col=0).to_dict()['FAST-JX-INPUT']

    # If you have that version, then set the path to its J2J file ...
    if version in list(vfjx.keys()):
        j2j_file = os.path.join(
            photolysis_paths['FJX_Files'], 'FJX_j2j_'+vfjx[version]+'.dat')
        exists = os.path.exists(j2j_file)
        if exists is True:
            print('Using: FJX_j2j_'+vfjx[version] +
                  '.dat for GEOS-Chem version:  ', version)
        else:
            raise FileNotFoundError('We could not locate the correct FJX file to open for GEOS-Chem Version: '+version +
                                    ' at the following expected path:  \n' +
                                    '    '+j2j_file).with_traceback(sys.exc_info()[2])
    else:
        raise ValueError('Version ', version, ' of GEOS-Chem Classic was not found in the dictionary \n' +
                         'defining its corresponding FJX chem input file... \n\n ' +
                         'This likely means you are trying to Emulate a new version of GEOS-Chem not yet supported by our Emulator. \n'
                         'To proceed, you will need to figure out what FJX version is used for this version of GEOS-Chem and then: \n' +
                         '    (1) Open and edit the file: ".../GEOSChem_Emulator/Input_Files/Photolysis_Files/FJX_input_by_GC-version.xlsx" \n' +
                         "           - Add the GEOS-Chem version you're emulating in the FIRST column \n" +
                         '           - Add the version of the FJX file that GEOS-Chem version uses to the the SECOND column.' +
                         '    (2) Add the Raw FJX.dat file to the directory: "GEOSChem_Emulator/Input_Files/Photolysis_Files/FJX_Files/" if it does not exist already. \n\n' +
                         'Why? #1 Will make sure we know what FJX file to open/read when emulating this version of GEOS-Chem in our dictionary that popped this error and \n' +
                         '#2 Will make sure we have the actual FJX file to open and read.').with_traceback(sys.exc_info()[2])

    # Read in info about what Photolysis Rxns are Defined in this FJX File into a pandas dataframe!
    fjx = read_FJX_file(j2j_file)
    # Create an empty list-like column to hold J-values...
    fjx['MCM_J2Use'] = list([''])*len(fjx)

    # Get FJX-Cross Section --> MCM J-Value Mapping done in do_jval_mapping.py
    cs = pd.read_excel(
        photolysis_paths['FJX_cross_sect_to_jvalues.xlsx']).fillna('')

    # Assign J-Values in mechanism based on the Cross-Section used in FAST-JX:
    for ind in fjx.index:  # All data is defined in the dataframe fjx for this matching.
        # Cross-Section we need a J-Value for...
        csect_need = fjx.loc[ind, 'FJX_CrossSec']
        # Index of the Reaction in PHOTOL()
        fjx_indx = fjx.loc[ind, 'FJX_Index']

        # If You can't figure out what cross section is needed...
        if csect_need == '':
            # Figure out what reaction this is for...
            rx_i = fjx.loc[ind, 'Photolysis_Rxn'].split(
                '->')[0].replace(' ', '')

            # Some versions don't have names for NITP cross sections... Fix that.
            if any([nit_str in rx_i for nit_str in ['NIT+hv', 'NITs+hv']]):
                csect_need = 'NITP'
                fjx.at[ind, 'FJX_CrossSec'] = 'NITP'
                fjx.at[ind, 'J-NickName'] = 'j'+fjx.loc[ind,
                                                        'Photolysis_Rxn'].split('+')[0].replace(' ', '')

        # Get Index in dataframe, cs, where that cross-section defined:
        has = cs.index[cs['FJX_Cross_Section'] == csect_need]

        if len(np.unique(cs.loc[has, jmap_type])) == 1:  # If we have a match
            # Then use that in the df.
            fjx.at[ind, 'MCM_J2Use'] = np.unique(cs.loc[has, jmap_type])[0]

        # Fill in info about everything used in the output fjx_info dataframe!
        if fjx.loc[ind, 'Quantum_Yield'] != 1:
            fjx.at[ind, 'F0AM_Assignment'] = str(
                fjx.loc[ind, 'Quantum_Yield'])+'.*('+fjx.loc[ind, 'MCM_J2Use']+')'
        else:
            fjx.at[ind, 'F0AM_Assignment'] = fjx.loc[ind, 'MCM_J2Use']

    # Get a list of all duplicate j-nicknames in FJX
    unique_Jnames = set()
    dupe_Jnames = []
    for ji, j in enumerate(fjx.loc[:, 'J-NickName']):
        if j in unique_Jnames and j not in dupe_Jnames:
            dupe_Jnames.append(j)
        elif j not in unique_Jnames:
            unique_Jnames.add(j)

    for j in dupe_Jnames:
        # ind in fjx of all matches to a dupe name
        idx_matches = fjx.index[j == fjx.loc[:, 'J-NickName']]
        assigns = [fjx.loc[ix, 'F0AM_Assignment'] for ix in idx_matches]
        is_diff = [True if assignment == assigns[0]
                   else False for assignment in assigns]

        # If some of the assignments vary for the same J-Nickname...
        if not all(is_diff):
            # Loop over all mataches to this nickname
            for ii, ix in enumerate(idx_matches):
                # If it doesn't match the assignment of the one already defined (1st occurance) or is 1st occurance...
                if (is_diff[ii] == False) or (ii == 0):
                    # Update the name to be unique!
                    new = make_unique_Jname(j, list(fjx.loc[:, 'J-NickName']))
                    fjx.at[ix, 'J-NickName'] = new

    for ix in fjx.index:
        fjx.at[ix, 'PHOTOL(IND)'] = 'PHOTOL('+str(fjx.loc[ix, 'FJX_Index'])+')'
        fjx.at[ix, 'Info'] = '% PHOTOL('+str(fjx.at[ix, 'FJX_Index'])+')'+' used for: '+fjx.at[ix, 'Photolysis_Rxn'] + \
            ' based on '+fjx.at[ix, 'FJX_CrossSec'] + \
            ' crossection with QY='+str(fjx.at[ix, 'Quantum_Yield'])
        fjx.at[ix, 'GEOSChem_Jline'] = fjx.at[ix, 'J-NickName']+'='+fjx.at[ix, 'F0AM_Assignment']+'; % PHOTOL('+str(fjx.at[ix, 'FJX_Index'])+') \n' +\
            '% Used for: '+fjx.at[ix, 'Photolysis_Rxn']+' based on '+fjx.at[ix,
                                                                            'FJX_CrossSec']+' crossection with QY='+str(fjx.at[ix, 'Quantum_Yield'])

    # And also create a dictionary to hold results mapping PHOTOL(IND) => J-NickName in F0AM!
    jdict = dict(zip(fjx['PHOTOL(IND)'], ["'" + j + "'" for j in fjx['J-NickName']]))

    return jdict, fjx

# %%###########################################################################
# -----------------(A) CREATE THE GEOSCHEM_AllRxns.m File ---------------------
###############################################################################

def make_GEOSCHEM_AllRxns_file(kppfile: str, GC_version: str, 
                               jmap_type: str = 'CrossSec_Match',
                               include_het_rxns: bool = False, 
                               output_fname: str = '', output_dir: str = '',
                               verbose: bool = True, overwrite: bool = False):
    """Function to convert a GEOS-Chem KPP File to a F0AM Compliant Mechanism. 

   First, this function figures out which FJX File corresponds to this GC version & 
   builds a dictionary that allows us to map J values from GEOS-Chem to their MCM/F0AM 
   corresponding value. Then, we parse the KPP file to get list of rxns and rates we 
   need to write a F0AM file. Next, we format that output so we can pass it to the F0AM 
   IO files to write a preliminary mechanism file. If the user does NOT want to include 
   het rxns, then we find stuff only created in the preliminary mech from a het rxn and
   remove that chemistry & rebuild the prelim mech! Next, we build a custom header for 
   the GEOSChem_AllRxns.m file and writes the final mechanism with its header to the 
   output file. This function OUTPUTS two things: a list of all unique reaction rate 
   definations that must be defined in the GEOSChem_K.m file and a dictionary with 
   keys about each unique rate call and what rxns use that rate so that we can build 
   the reaction rates file needed in F0AM. 


    INPUTS: 
    -------
        (1) kppfile          - STR with full path to the KPP input file that's being parsed. 

        (2) GC_version       - STR formatted as '13.3.4' corresponding to what GEOS-Chem 
                              version you are seeking to create a F0AM file for. 

        (3) jmap_type        - STR indicating HOW the j-value mapping from GEOS-Chem to 
                                MCM j-values should be done. Default is 'CrossSecMatch'.

        (4) include_het_rxns - BOOL indicating whether or not heterogenous reactions should 
                                be included in the output mechanisms or not                   

        (5) output_fname     - (OPTIONAL) STR with name of the output file to write. Default is 'GEOSCHEM_AllRxns.m'

        (6) output_dir     - (OPTIONAL) STR with full path where the outputfile should be written. 
                                Default is set to: "Your_Path_To/GC_Emulator/Output_Files/"

        (7) overwrite        - (OPTIONAL) BOOL of whether or not to overwite the output file or 
                                create a new one if it exists already. Default is False. 

        (8) verbose          - (OPTIONAL) BOOL of whether or not to print extras/progress.

    OUTPUTS: 
    --------

    (1) all_rfuncts    - LIST of all rate functions used in the mechanism that must be defined
                        in GEOSChem_K.m for the mechanism to work. 

    (2) rxn_rate_dict - DICT with keys corresponding to unique rate functions called in mechanism, 
                        that point to values of a nested DICT with key/value pairs of:  
                            'rxns'           => LIST of rxns that use that specific rate call 
                            'rate_functions' => LIST of the rate function being called as a STR, 
                                                no parenthesis or arguemnts to rate
                            'rate_args'      => LIST of a LIST with the rate function arguements 
                                                passed to that in 'rate_function' that is passed to 
                                                make_GCRates_file to make 'GEOSChem_K.m'
    REQUIREMENTS: 
    -------------
        LIBRARIES:           import numpy as np 

        CUSTOM FUNCTIONS:    create_jmap_dict()
                             parse_kpp_main()
                             KPP_Dict_to_F0AM_IO()
                             eliminate_dead_ends()
                             fio.build_all_from_rxns()
                             make_header_for_F0AM_mech()

    USAGE: 
    -----  
        CALLED WITHIN: 

        OUTPUT USAGE: 

    AUTHOR: 
    -------
        Dr. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD

    CHANGE LOG: 
    ----------
            10/30/2023    JDH Created 
            11/6/2023     JDH added need for fjx_df to get info on photo rxns... 

    """

    # Read in the right FJX file for this GC-Version and build the j-mapping
    # dictionary to use to map KPP J's to F0AM J's.
    jmap_dict, fjx_df = create_jmap_dict(version=GC_version, jmap_type=jmap_type)

    # Parse the KPP file. Get list of rxns and rates we need to write a F0AM file.
    print('Parsing KPP file...') 
    var_info, kpp_df, kpp_dict = parse_kpp_main(kppfile)

    # Prep reactions & rates for F0AM: 
    #rxn_df=prep_df4f0am(rxn_df)
    print('Preparing for F0AM...')
    rxn_df=prep4f0am(kpp_df, fjx_df, jmap_dict)

    # Pull out info about dif rxn types & get it in the format we need to pass to build_mech()
    nhet_df = rxn_df[(rxn_df['is_het'] == 0)]
    het_df = rxn_df[(rxn_df['is_het'] == 1)]
 
    # Build a F0AM compliant mechanism with ONLY the het rxns ...
    het_info = fio.build_all_from_rxns(list(het_df['Reactions']), k_list=list(het_df['Rates']),
                                       sort_rxn_i=False, sort_rxn_list=False, verbose=verbose)
    
    # Build a F0AM compliant mechanism without all the het rxns ...
    nhet_info = fio.build_all_from_rxns(list(nhet_df['Reactions']),k_list=list(nhet_df['Rates']),
                                   sort_rxn_i=False, sort_rxn_list=False, verbose=verbose)
       
    # ---------------------------------------------------------------------------------------------
    # Build your mechanism either with het rxns in their own file  or without het rxns...
    # ----------------------------------------------------------------------------------------------
    if include_het_rxns is False:
        
        # Parse the non-het mech to remove reactions of stuff formed ONLY via a het rxn...
        info, rxns2drop, inds2drop = eliminate_dead_ends(het_info, nhet_info)
        [sp_list, rxns, fs, gs, rates, rct_cmpds, prd_cmpds, rct_ylds, prd_ylds]=info
        #[sps, rxns, rates, gs, fs, rctt, prdd, r_yld,p_yld]= info
        
        if verbose is True and len(rxns2drop) > 0:
            print('Dropping these gas/ photolysis rxns from mech because they ' +
                  'involve species only formed in het reactions, which are not included:\n')
            [print('    '+r) for r in rxns2drop]
            print('\n')

        # Rebuild your mech without these rections/ species they form...
        out = fio.build_all_from_rxns(rxns, k_list=rates, sort_rxn_i=False, sort_rxn_list=False, verbose=verbose)
        [sp_list, rxns, fs, gs, ks, rct_cmpds, prd_cmpds, rct_ylds, prd_ylds] = out
        
        # Add rxns, fs, gs, ks, to outputdf :
        df=pd.DataFrame({'Rxn':rxns,'Ks':ks,'Gs':gs,'Fs':fs})
        
        # Drop those same rxns from rxn df and use to get required k's & js. 
        mech_df=nhet_df.drop(nhet_df.index[inds2drop]).reset_index() 
        
        # Get a list of unique reaction rate functions (K's, not J's) used in the resulting mechanism: 
        ks_required=utils.flatten_nested_list(mech_df['Rate_Functions'][mech_df['is_photo']==False].to_list(), drop_dupes=True)
        js_required=list(np.unique(mech_df['Rate_Functions'][mech_df['is_photo']==True]))
            
    # else:
    #     # Unpack stuff from each mechanism: 
    #     [het_sp_list, het_rxns, fhet_s, het_gs, het_ks, het_rct_cmpds, het_prd_cmpds, het_rct_ylds, het_prd_ylds] = het_info
    #     [sp_list, rxns, fs, gs, ks, rct_cmpds, prd_cmpds, rct_ylds, prd_ylds] = nhet_info
        
    #     # Pull references to Species declared in main rxn file out of the het rxn file (e.g. duplicate declarations)
    #     [het_sp_list.remove(sp) for sp in sp_list  if sp in het_sp_list]

    #     # Make sure rates are in MATLAB syntax...
    #     het_ks=[_str_multi_replace(k, {'.FALSE.':'false', '.TRUE.':'true'}) for k in het_ks]
        
    #     # Get list of required Het  rxns... 
    #     het_ks_required=np.unique(het_df['rate_function'][(het_df['is_photo']==False) & (het_df['rate_function']!='None')])

    #     # # Make a mechanism title for the het rxns file:
    #     het_title, het_set_globals, het_call, het_jzero = make_header_for_F0AM_mech(het_rfuncts,[], kppfile)

    #     # het_fname=output_fname.split('.m')[0]+'_HetRxns.m'
    #     # fio.write_mech_to_file(het_fname, hsp_list, [''],hrxns, hks,hgs,hfs,mech_title=het_title,
    #     #                lines_to_add_before_mech=het_set_globals+het_call+het_jzero)

    # Make a header & add info to the (main) gas/photolysis rxns file
    mech_title, mech_comments = make_header_for_AllRxns(kppfile, GC_version, ks_required, js_required, jmap_type)

    # Now, FINALLY, Write  gas/photolysis  mechanism to a F0AM file:
    fio.write_mech_to_file(output_fname,  species_list=sp_list,  ro2_list=[''],
                           rxn_list=rxns,       k_list=ks,       g_list=gs,  f_list=fs,
                           mech_name=output_fname.replace('.m', ''), mech_title=mech_title,
                           mech_comments=mech_comments,
                           out_filepath=output_dir, overwrite=overwrite)

    # # Unique list of all rate functions used in Het and regualr rxns...
    # if include_het_rxns is True:
    #     all_rfuncts = np.unique(ks_required+het_ks_required).tolist()
    # else:
    #     all_rfuncts = np.unique(ks_required).tolist()

    return mech_df, ks_required,  js_required, fjx_df

# %%###########################################################################
# -----------------(B) CREATE THE GEOSCHEM_K.m File ----------------------------
###############################################################################
def make_GEOSCHEM_K_file(kppfile: str, GC_version: str, template_file: str, 
                         ks_required: list, js_required: list, 
                         rxn_rate_dict: dict,
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

        template_file - STRING containing the full path to the template file for GC_Rates_Template.txt

        ks_required   - LIST of all rate functions used in the mechanism that must be defined
        
        js_required   - LIST of all named J-values used in the mechanism that must be defined
        
        output_fname  - (OPTIONAL) STR with name of the outputfile to write. Default is 'GEOSCHEM_K.m'

        output_dir  - (OPTIONAL) STR with full path where the outputfile should be written. Default is to "GC_Emulator/Output_Files/"

        overwrite     - (OPTIONAL) BOOL of whether or not to overwite the outputfile or create a new one if it exists already. Default is False. 


    OUTPUTS:
    --------
        A F0AM compliant rate defintion file saved at output_dir+output_fname. 

    REQUIREMENTS: 
    -------------
        LIBRARIES:           from datetime import date  

        CUSTOM FUNCTIONS:    pylb.join_list_for_MATLAB()
                             pyMCM_utils.check_filename()

    USAGE: 
    -----  
        CALLED WITHIN: 

        OUTPUT USAGE: 


    AUTHOR: 
    -------
            Dr. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD

    CHANGE LOG: 
    ----------
        10/30/2023    JDH Created 

    """
    
    # Create a formatted list of all the required rate functions (K's) that must 
    # be imported & named J's. This info goes IN THE HEADER (commented out). Format 
    # it so both are in in a commented out comma delimited list with line breaks 
    # where  necessary (for long MATLAB lists that continue on multiple lines). 
    hdr_k_list= utils.join_list_for_MATLAB(',', [k+'()' for k in ks_required], # Add '()' after func name... 
                                           insert='                    ', 
                                           insert_skip=[1],comment=True, 
                                           comment_skip=[1]) 
    
    hdr_j_list= utils.join_list_for_MATLAB(',', js_required,
                                           insert='                    ', 
                                           insert_skip=[1],comment=True, 
                                           comment_skip=[1]) 
   
    # Now, create the correctly formatted list w/ same info that is actually used 
    # to IMPORT THESE RATE FUNCTION HANDLES in the actual MATLAB code 
    # unlike above, this is NOT commented/ in header of the file:  
    rate_list= utils.join_list_for_MATLAB(',', ks_required, insert='   ', 
                                          insert_skip=[1])
    
    ###########################################################################
    # Read in GC_Rates Template and customize the header used ...
    ###########################################################################
    # Initialize list of lines we'll write to the  output file...
    function_lines = list([])

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
        line = _edit_line(line, '**GC_VERSION**', GC_version) 
        
        # Add to HEADER what date this fiel was generated on in the function header
        line = _edit_line(line, '**THE_DATE**', str(date.today())) 
        
        # Add the path to the KPP file that was used to generate this mechanism
        line = _edit_line(line, '**KPP_FILE_PATH**', kppfile)  
        
        # Add to HEADER the number of species declared in the mechanism: 
        line = _edit_line(line, '**N_SPECIES**', n_species) 
        
        # Add to HEADER the number of reactions in the mechanism: 
        line = _edit_line(line, '**N_RXNS**', n_rxns) 
        
        # Add to HEADER  how j-values were mapped: 
        line = _edit_line(line, '**JMAP_TYPE**', jmap_type) 
        
        # Add to HEADER which named J-values are required for mechanism: 
        line = _edit_line(line, '**HDR_LIST_OF_REQUIRED_J_NAMES**', hdr_j_list) 
        
        # Add to HEADER which functions are required for rate constants: 
        line = _edit_line(line, '**HDR_LIST_OF_REQUIRED_RATE_FUNCTIONS**', hdr_k_list) 
        
        # Update actual MATLAB CODE to import only the handles of the rate
        # functions used in the mechanism in the rate dictionary 
        line = _edit_line(line, '**LIST_OF_RATE_FUNCTS_TO_IMPORT**', rate_list) 
    
        # Add the line that unpacks the handles from the rate dictionary for use
        line = _edit_line(line, '**LIST_OF_IMPORTED_FUNCTION_HANDLES**', rate_list) 
        
        # After modifying lines (if needed), save the line to output list!
        function_lines.append(line)

        count += 1  # Update counter.

        if count > num_lines:
            break  # Exit while loop if you reach end of file.

    # Close the template file now that you got all the good stuff out.
    inF.close()

    ###########################################################################
    # Write the GEOSCHEM_vX_X_K.m rates file:
    ###########################################################################

    # Decide what to call the output file & make sure the path to it exists...
    rate_file = utils.check_filename(filename=output_fname, default_name='GEOSCHEM_K.m', ext='.m',
                                     savepath=output_dir, overwrite=overwrite, return_full=True)

    # Open the Output file to Write
    outF = open(rate_file, "w")

    # Write each line to the output file
    [outF.write(ln) for ln in function_lines]

    outF.close()  # Close the output file
    print('F0AM Rate Defintion file saved at: '+rate_file)

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
    for key in fdict:  # Loop over all functions
        glbs_tf = False
        # list of lines we need to write for each function.
        ln_list = fdict[key]
        for i, line in enumerate(ln_list):
            outF.write(line+'\n')
            if glbs_tf is False and ')' in line:
                outF.write(get_globals+'\n')
                glbs_tf = True

    outF.close()  # Close the output file
    print('MATLAB compliant version of gckpp_Rates.F90 called in "GEOSCHEM_K.m" to define all rate functions saved at: '+rate_file)

    return


# %%###########################################################################
# -----------------(D) Wrapper Function to Make ALL Files -----------------------
###############################################################################
# D.1 L1-Function used in make_GC_mechanism() to check user inputs to make_GC_mechanism()
def _check_user_inputs(kppfile, rate_files, GC_version: str, jmap_type: str, output_dir: str,
                       template_paths: dict, photolysis_paths: dict):

    # Check that the output directory the user entered exists:
    if os.path.exists(output_dir) is False:
        raise NotADirectoryError("INPUT ERROR in make_GC_mechanism(): \n" +
                                 "The following output directory you listed to store the output files does not exist: \n" +
                                 "    output_dir='"+output_dir + "'").with_traceback(sys.exc_info()[2])

    # Check that the KPP files passed as input exist
    if os.path.isfile(kppfile) is False:
        raise FileNotFoundError("INPUT ERROR in make_GC_mechanism(): \n" +
                                "The path to the KPP file you entered does not exist: \n" +
                                "    kppfile='"+kppfile + "'").with_traceback(sys.exc_info()[2])

    # Check that the gckpp rates files passed as input are the right type & exist:
    if type(rate_files) not in [type('str'), type(list([1]))]:
        raise TypeError("INPUT ERROR in make_GC_mechanism(): \n" +
                        "Input for argument 'rate_files' must either be TYPE = STR or LIST: \n" +
                        "    (1) A STR with a path to a gckpp_Rates.F90 file \n " +
                        "    (2) A LIST of strings with paths to rate function files like: \n" +
                        "             ['/path/to/gckpp_Rates.F90', /path/to/gckpp_HetRates.F90'] \n" +
                        "Instead, you entered input of 'rate_files' of type="+str(type(rate_files))).with_traceback(sys.exc_info()[2])
    # If its a single rate file ...
    elif type(rate_files) == type('str') and os.path.isfile(rate_files) is False:
        raise FileNotFoundError("INPUT ERROR in make_GC_mechanism(): \n" +
                                "The path to the GCKPP rate file you entered does not exist: \n" +
                                "    rate_files='"+rate_files+"'").with_traceback(sys.exc_info()[2])
    # If its a list of rate files..
    elif type(rate_files) == type(list([1])) and not all([os.path.isfile(file) for file in rate_files]):
        err_files = ["    ("+str(i+1)+") '"+file+"'" for i,
                     file in enumerate(rate_files) if os.file.exists(file) is False]
        raise FileNotFoundError("INPUT ERROR in make_GC_mechanism(): \n" +
                                "The path to the following gckPPP rate file(s) you entered does not exist: \n" +
                                '\n'.join(err_files)).with_traceback(sys.exc_info()[2])

    # Check that the template files exist at these paths.
    for key in list(template_paths.keys()):
        if os.path.isfile(template_paths[key]) == False:
            raise FileNotFoundError("INPUT ERROR in make_GC_mechanism(): \n" +
                                    "The path to the "+key+" template file was not found where expected at: \n" +
                                    "    '"+template_paths[key]+"'\n" +
                                    "It should be located at: '.../GEOSChem_Emulator/Input_Files/Template_Files/"+key+"' + \n" +
                                    "and if you get this error, there may be an install issue!").with_traceback(sys.exc_info()[2])

    # Check that the Photolysis files exist at these paths.
    for key in list(photolysis_paths.keys()):
        if os.path.exists(photolysis_paths[key]) == False:
            raise FileNotFoundError("INPUT ERROR in make_GC_mechanism(): \n" +
                                    "The path to the "+key+" photolysis input file was not found where expected at: \n" +
                                    "    '"+photolysis_paths[key]+"'\n" +
                                    "It should be located at: '.../GEOSChem_Emulator/Input_Files/Photolysis_Files/"+key+"' \n" +
                                    "and if you get this error, there may be an install issue!").with_traceback(sys.exc_info()[2])

    # Check that the jmap_type is allowed:
    if jmap_type not in ['CrossSec_Match', 'MCM', 'RCIM']:
        # Be case insensitive & check for miss-spelling
        if jmap_type.lower() in ['crosssec_match', 'crossec_match', 'mcm', 'rcim']:
            if jmap_type.lower() in ['crosssec_match', 'crossec_match']:
                jmap_type = 'CrossSec_Match'
            if jmap_type.lower() == 'mcm':
                jmap_type = 'MCM'
            if jmap_type.lower() == 'rcim':
                jmap_type = 'RCIM'
        else:
            raise ValueError("INPUT ERROR in make_GC_mechanism(): \n" +
                             "Input for argument 'jmap_type' must be one of the following: \n" +
                             "    (1) 'CrossSec_Match' \n" +
                             "    (2) 'MCM' \n" +
                             "    (3) 'RCIM' \n" +
                             " For more info on these options try: >>utils.jmap_type_help()").with_traceback(sys.exc_info()[2])

    return jmap_type


def make_GC_mechanism(kppfile, rate_files, GC_version: str, jmap_type: str, include_het_rxns: bool,
                      output_dir: str, overwrite: bool = False):
    """Main function to write a F0AM version of a GEOS-Chem mechanism.

    INPUTS:

        (1) kppfile - string containing the path to a GEOS-Chem KPP file (path_to/fullchem.eqn')
        (2) rate_files- list of strings OR string contaiing the path to a GEOS-Chem a 'gckpp_Rates.F90' file. 
        (3) GC_version - Version of GEOS-Chem being Emulated (so chooses right J-files!). 
        (4) jmap_type - STRING indicating HOW you want J-values to be mapped. 
        (5) include_het_rxns-  bool indicating whether or not you want heterogenous rxns to be included in the resulting mech. 
        (6) output_dir -str containing the path where the output F0AM files should be located. 
        (7) overwrite - bool indicating whether or not output files should overwrite any exisiting files there or not. 

    """
    global master_outdir
    master_outdir=output_dir
    
    # Get the relative paths to input/templates, etc...
    GCEM_path, input_path, template_paths, photolysis_paths, jmap_type_help = utils.get_my_relative_paths()

    # Check the user's inputs to make sure all files we need exist and are located where we expect them.
    jmap_type = _check_user_inputs(kppfile=kppfile, rate_files=rate_files, GC_version=GC_version,
                                   jmap_type=jmap_type, output_dir=output_dir, template_paths=template_paths,
                                   photolysis_paths=photolysis_paths)

    # Write the GEOSChem_AllRxns.m file & get output we need to write GEOSChem_K.m file:
    mech_df, ks_required,  js_required, fjx_df = make_GEOSCHEM_AllRxns_file(kppfile=kppfile,
                                                                    GC_version=GC_version,
                                                                    jmap_type=jmap_type,
                                                                    include_het_rxns=include_het_rxns,
                                                                    output_fname='GEOSCHEM_AllRxns.m',
                                                                    output_dir=output_dir,
                                                                    overwrite=overwrite)

    # # Write the GEOSChem_K.m file
    # make_GEOSCHEM_K_file(kppfile=kppfile,
    #                      GC_version=GC_version,
    #                      template_file=template_paths['GEOSChem_K_template.txt'],
    #                      all_rfuncts=all_rfuncts,
    #                      rxn_rate_dict=rxn_rate_dict,
    #                      output_fname='GEOSCHEM_K.m',
    #                      output_dir=output_dir,
    #                      overwrite=overwrite)

    # # Create import_GC_rates.m file where all rate functions are defined which is
    # # referenced and used as is without needing modifications in the GEOSChem_K.m file.
    # if type(rate_files) == str:
    #     rate_files = [rate_files]

    # fdict = make_import_GC_rates_file(rate_files=rate_files,
    #                                   all_rfuncts=all_rfuncts,
    #                                   template_file=template_paths['import_GC_rates_template.txt'],
    #                                   include_het_rxns=include_het_rxns,
    #                                   output_fname='import_GC_rates.m',
    #                                   output_dir=output_dir,
    #                                   overwrite=overwrite)

    return fjx_df

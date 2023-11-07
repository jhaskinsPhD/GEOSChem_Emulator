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
import utils 

sys.path.append('/uufs/chpc.utah.edu/common/home/u6044586/python_scripts/modules/pyF0AM/')
import pyF0AM as fio 


# %%---------------------------------------------------------------------------
#-------------- (A.4) make_header_for_F0AM_mech() & Sub-Functions -------------
#------------------------------------------------------------------------------
def make_header_for_GEOSCem_J(photo_functs): 
    """Function to make a descriptive header for the GEOSChem_J.m file. 
    
       (4) photo_functs  - LIST of unique J-values that must be defined in F0AM 
                            for the resulting mechanism to work outputted from KPP_Dict_to_F0AM_IO()  
    
    """
    cmt_space1= '%                    '; cmt_space2='%              '
    j_info =    '%    j_values={' +','.join([f if np.mod(i,15)!=0 or i==0 else cmt_space2+f for i,f in enumerate(photo_functs) ] )+'};' 
    jzero='Jzero=0.*J3; % Set photolysis=0 for stratospheric only rxns \n' 


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
    r_functs=np.unique(r_functs).tolist()
    
    # Join the list of J-values/ Rate Functions used in this mech to put into the mechanism's header. 
    cmt_space1= '%                    '; cmt_space2='%              '
    rate_info = '%    rate_functions={' +','.join([f+ '()' if np.mod(i,5)!=0 or i==0 else cmt_space1+f+ '()'for i,f in enumerate(r_functs) ])+'};'  
    j_info =    '%    j_values={' +','.join([f if np.mod(i,15)!=0 or i==0 else cmt_space2+f for i,f in enumerate(photo_functs) ] )+'};' 
    
    # Make Title Containing details + list of functions this rxn file requires to run.  
    mech_title=''.join(['% F0AM Compliant GEOS-Chem Mechanism for version'+GC_version +' generated on '+str(date.today())+' from KPP file: ',  
                        '%    '+kppfile+ ' using jmap_type =' +jmap_type])
                       
    mech_comments=''.join(['% List of reaction rates that must be defined in "GEOSChem_K.m" for this mech to work:\n',
                           #rate_info, 
                           '\n\n',
                           '% List of photolysis rates that must be defined in "GEOSChem_J.m" for this mech to work:\n'])
                           #j_info]) 
    
    return mech_title, mech_comments

# A.4 L1-function called within make_GEOSCHEM_AllRxns_file() 
def eliminate_dead_ends(het_info, info):
    """ID Reactions from stuff that's ONLY formed via a Heterogenous reaction
       Basically cut dead end reactions that we don't need...
       
       Custom Functions: 
           flatten_nested_list  --> in pyMCM/F0AM_Tools/utils.py
           enforce_list_inds    --> in pyMCM/F0AM_Tools/utils.py
           find_in_list         --> in pyMCM/F0AM_Tools/utils.py
           drop_dupes           --> in pyMCM/F0AM_Tools/utils.py
           """
    [hspecies, hrxn, hf, hg, hk, hrct_cmpds, hprd_cmpds, hrct_ylds, hprd_ylds]= het_info
    [species, rxn, f, g, k, rct_cmpds, prd_cmpds, rct_ylds, prd_ylds]= info
    
    # A list of all products formed in any het rxns ... 
    het_prds=utils.flatten_nested_list(hprd_cmpds, drop_dupes=True);
    
    growing=True; start=True 
    while growing: 
    
        if start==True: # Initialize everything on first pass through. 
            rct_a= rct_cmpds; prd_a=prd_cmpds; rxn2drop=[]; start=False; known=[]
        else: 
            # Pretend you dropped all those reactions... whats a dead end now? 
            keep=[v for v in range(0, len(rxn)) if v not in rxn2drop]
            [rct_a,prd_a]=utils.enforce_list_inds(keep, [rct_cmpds, prd_cmpds])
        
        # Get a unique list of reactants and products that are still left in the non-het rxns.
        prds=utils.flatten_nested_list(prd_a, drop_dupes=True);
                                          
        sz0=len(rxn2drop)             
                 
        # Find all products in het rxns that are not formed from other reactions. 
        dead=list()
        [dead.append(sp) for sp in het_prds if (sp not in prds) and (sp not in ['H2O', 'RO2', 'hv']+dead)] 

        # Loop over all species ID'd as only forming through a HET rxn. 
        for cmpd in dead:
            # Find all reactions in the gas+photo rxns that this het only product is a REACTANT in. 
            ind1,args= utils.find_in_list(cmpd, rct_cmpds, match_if_contains=True)
            
            # Add the index of this rxn (in gas+photolysis rxns) to drop. 
            [rxn2drop.append(i) for i in ind1 if i not in rxn2drop ]
            
            # Add it to the list of compounds you know you need to drop. 
            if cmpd not in known: known.append(cmpd)
            
        sz1=len(rxn2drop)
        if sz0==sz1: break # list stopped growing, so stop! 
        
    keep=[v for v in range(0, len(rxn)) if v not in rxn2drop]
    
    # Keep the reactions that aren't dead ends. 
    [rxns, ks,gs,fs,rctt,prdd,r_yld,p_yld]=utils.enforce_list_inds(keep,[rxn,k,g,f,rct_cmpds,prd_cmpds,rct_ylds,prd_ylds])
    
    # Make a new species & RO2 list using stuff in reactants and products. 
    sps=utils.drop_dupes(utils.flatten_nested_list(rctt, drop_dupes=True)+utils.flatten_nested_list(prdd, drop_dupes=True)) 
    [sps.remove(item) for item in ['hv','RO2'] if item in sps]
    sps.sort()
     
    # Return all the good stuff you'd need to write a mech. 
    to_drop=[rxn[l] for l in rxn2drop]

    return [sps, rxns, ks,gs,fs,rctt,prdd,r_yld,p_yld], to_drop
#------------------------------------------------------------------------------
#-------------- (A.3) KPP_Dict_to_F0AM_IO() & Sub-Functions -------------------
#------------------------------------------------------------------------------
# A.3.3 L2-function called within KPP_Dict_to_F0AM_IO within make_GEOSCHEM_AllRxns_file() 
def _append_to_ndict_list(ddict_in:dict, key1:str, key2:str, values, allow_dupes:bool=False): 
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
    ddict=ddict_in.copy() 
    
    # Make sure the input is a list so the fuction works. 
    if type(values) !=list: values=[values]
    
    # If key 2 doesn't exist, add the values as a list underneath it. 
    if key2 not in list(ddict[key1].keys()): 
        ddict[key1][key2]=values 
    else:
        # Don't worry about dupes, just append the values to the nested dict. 
        if allow_dupes is True: 
            ddict[key1][key2]=ddict[key1][key2]+values
        else: 
            # Loop over values and only add those that aren't in there already.
            for v in values: 
                if v not in ddict[key1][key2]: 
                    ddict[key1][key2]=ddict[key1][key2]+[v]
    
    return ddict

# A.3.2 L2-function called within KPP_Dict_to_F0AM_IO within make_GEOSCHEM_AllRxns_file() 
def get_all_MCM_Js_needed(J_in:str, fjx_df:pd.DataFrame, all_Js:list=[]): 
    """Function to parse and extract a unique list of MCM J-values referenced in input string, 
    'j_in' and append these to a "master list" of all unique Js referenced in the mechanism. 
    
    Because SOME MCM j-values are defined as a combination of scale_factors/ different, we 
    can't just use the list of rates we pass to build_mech() to keep track of all the 
    individual J-values that must be defined in the MCM for the resulting mechanism to 
    work... So this function exists to pull those out and add to a master list, "all_Js". 
    
    INPUTS:
    -------
        (1) J_in - A STRING containing a MCM j-value that (may) be a combo of 
                   scale factors/ multiple Js. 
               
        (2) fjx_df - A Pandas Dataframe outputted from create_jmap_dict() with info on all J-rxns 
       
        (3) all_Js - (OPTIONAL) A LIST containing all the unique MCM J-values used in 
                     rate definitions. Values in J_in are APPENDED to this list if they're new!
        
    OUTPUTS: 
    ------- 
        all_Js - An updated LIST of all the individual MCM J-values referenced in 
                what we'll pass to build_mech() including anything new defined in the input string, 
                "j_in" after its J-Values are seperated from scale factors & MATLAB math... 
                
    EXAMPLE: 
    --------
        Parses jMVK_a='0.5.*(J22+J15)', extracts ['J22', 'J15'], and adds them to the output list 'all_Js'
            OR
        Parses jMVK_b='J54.*4.6', extracts  ['J54'] , and adds it to the output list 'all_Js'
                             
    USAGE: 
    -----  
        CALLED WITHIN: convert_FJX_PhotoRxns_to_MCM() within KPP_Dict_to_F0AM_IO() within make_GEOSCHEM_AllRxns_file()
                
        OUTPUT USAGE: Output 'all_Js' is directly passed as output of convert_FJX_PhotoRxns_to_MCM
     
    REQUIREMENTS: 
    ------------- 
        LIBRARIES:           import re 
                             import numpy as np 
                             import sys 
                                             
        CUSTOM FUNCTIONS:    None
          
    AUTHOR: 
    -------
         Prof. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD
         
    CHANGE LOG: 
    ----------
         1/14/2022    JDH Created   
         10/2/2023    JDH Modified documentation & changed from string split to regex, added error pop! 
                                     
    """
    # Pull out the MCM assignment for this J-value name from fjx_df that we need to parse.. 
    fjx_ind=fjx_df.index[fjx_df['PHOTOL(IND)']==J_in]
    J_value=fjx_df.loc[fjx_ind[0],'F0AM_Assignment']
    
    # Define the Regex string to find J-value matches. 'Jregex' will hit str matches to the letter "J" 
    # followed by one or more numbers 0-9. Match will end when something other than 0-9 is encountered.
    # Won't match just "J" or return blank items in a list with spaces! 
    J_pattern1 = r'J[0-9]+'
    J_pattern2 = r'Jzero'# Catches "Jzero"... 
    J_pattern3 = r'Jn[0-9]+'# Catches "Jn3" etc. 

    # Use regex's "Findall()" to return a list of all matches to MCM type J-values in the string "J_in"
    # That follow the format defined in 'J_pattern'
    indv_Jlist=re.findall(J_pattern1,J_value)+re.findall(J_pattern2,J_value)+re.findall(J_pattern3,J_value)

    # Pop an error if no J-Values were found in "Jlist". 
    if len(indv_Jlist)==0: 
        raise ValueError("The function get_all_MCM_Js_needed() could not find any J-values in "+ 
                         "this string: '"+J_value+"' which is the assignment for '"+J_in+"'. \n\n"+ 
                         "This is either your fault (bad input) or means that the regex format string \n"+
                         "used in the function 'get_all_MCM_Js_needed()' needs to be modified to catch this instance!").with_traceback(sys.exc_info()[2])
    else: 
        # Remove any duplicates from list!  
        indv_Jlist= np.unique(indv_Jlist).tolist() 
    
    # Add all (new) unique J-values used in 'J_in' to the master list of all Js used in new MCM mech!  
    all_Js=np.unique(all_Js+indv_Jlist).tolist() 
        
    return all_Js

# A.3.1 L2-function called within KPP_Dict_to_F0AM_IO within make_GEOSCHEM_AllRxns_file() 
def convert_FJX_PhotoRxns_to_MCM(photo_dict, jmap_dict, fjx_df): 
    """Function used to format photolysis reactions and convert the rates for 
    photolysis reactions into lists we can use to build an MCM compliant mechanism. 
    Specifically, this is where the actual mapping of FJX/KPP photolysis rates into 
    their corresponding MCM photolysis rate is done. 
    
    INPUTS: 
    ------- 
       (1) photo_dict -         Dictionary of info about photolysis reactions generated by 
                                the function "parse_KPP" (a sub-selection of the dict, "rinfo"). 
                    
       (2) jmap_dict  -         Mapping Dictionary with keys corresponding to FJX phototlysis rates and 
                                values corresponding to F0AM corresponding rate values.
                                
       (3) fjx_df     -         Pandas df outputted from create_jmap_dict() 
       
    OUTPUTS: 
    --------
        (1) mapped_F0AM_JRates - List of Jrates to use in F0AM, with FJX values converted to 
                                 their F0AM counterpart and formatted with k(:,i) in front of them 
                                 and with a ';' behind! Should be in same order as values in 
                                 photo_dict! 
         
        (2) unq_MCM_Js_req -     A list of all individual J values that must be defined in F0AM 
                                 for our new mechanism to work since they're referenced in mech.
                 
    REQUIREMENTS: 
    -------------
        LIBRARIES:           import sys 
                             import numpy as np 
            
        CUSTOM FUNCTIONS:    get_all_MCM_Js_needed() 
        
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
    
    # Check that all rate functions in user input photo_dict['rate_function'] are
    # indeed "PHOTOL" (aka that all rates used in what we think are the 
    # photolysis reactions are actually photolysis rates. Otherwise pop an error! 
    unq_rates= np.unique(photo_dict['rate_function'])
    if len(unq_rates)!=1 or unq_rates[0]!='PHOTOL': 
        other_rates='\n'.join(['    '+rate for rate in unq_rates if rate !='PHOTOL'])
        raise ValueError("Not all rate functions in the 'photo_dict['rate_functions]' list passed to the function, \n" + 
                         "convert_FJX_PhotoRxns_to_MCM(), are 'PHOTOL' as expected.+ \n"+ 
                         "The following rate functions were also found: "+other_rates +
                         "\n\n  This likely means something weird is happening with: \n" + 
                         "    (1) How you detect a rxn is a photolysis rxn in the function, parse_KPP(). \n" + 
                         "                   OR \n   "+
                         "    (2) That something weird is happening in the function, _enforce_list_dict_inds(), \n"
                         "        which is used to sub-select photolysis info in the dict, 'rinfo', within the function, parse_KPP()").with_traceback(sys.exc_info()[2])
    
    mapped_F0AM_JRates=list([]) # Initialize empty list to hold mapped J-rates to use in F0AM! 
    unq_MCM_Js_req= list([]) # Initialize list of unique MCM Js required to be defined to be an empty list, we'll fill. 
    
    # Loop over all rates used in photolysis reactions within the input "photo_dict": 
    for j, rate in enumerate(photo_dict['rate']): 
        
        # Use the KPP rate parsed as the "key" to look up the j-value match in "jdict". 
        jkey=photo_dict['rate'][j]  # Both keys & rates should be formatted as 'PHOTOL(IND)'.
                
        # If this photolysis rate does NOT haev a correponding match to an MCM j-value in 'jdict, pop an error
        if jkey not in list(jmap_dict.keys()) or photo_dict['rate_args'][j][0] == '': 
            raise ValueError("In function, convert_FJX_rates_to_MCM(), we could not find a MCM J-Value \n" + 
                                 "match in 'jdict' corresponding to: '"+jkey + 
                                 "'\n\nThis likely means you need to update 'jdict' to account for new MCM/FJX mapping j-value changes!")
        
        else: #Otherwise, use the mapped value in J-Dict to define rate in MCM!  
                
            # Define the j-value rate that should be used for this reaction! 
            jind=fjx_df.index[fjx_df['PHOTOL(IND)']==jkey.replace(' ','')]
            j_info=fjx_df.loc[jind[0],'Info']
            mapped_F0AM_JRates.append('k(:,i) = '+jmap_dict[jkey].replace('EXP(','exp(')+';'+j_info)
                
            # Parse the MCM j-value you just added and add the individual 
            # Js used in that rate (since some defs include scale factors like '54.*J3') 
            # to a "MASTER" list of all J-s that must be defined in MCM for resulting mech to work! 
            unq_MCM_Js_req=get_all_MCM_Js_needed(rate,fjx_df,all_Js=unq_MCM_Js_req)
            
    return mapped_F0AM_JRates, unq_MCM_Js_req
   
# A.3 L1-function called within make_GEOSCHEM_AllRxns_file() 
def KPP_Dict_to_F0AM_IO(rxn_dict, is_photo:bool=False, jmap_dict:dict=dict({}),
                        fjx_df=pd.DataFrame(),make_rate_dict:bool=False): 
    """ Function to take a dictionary outputted from parse_KPP with info about 
    the rxns/rates and convert it into the format the build_mech() function 
    needs to create a F0AM compliant mechanism file.
    
    INPUTS: 
    -------
        (1) rxn_dict   - DICTIONARY/sub-dictionary outputted by parse_KPP with rxn info. 
        
        (2) is_photo   - (OPTIONAL) BOOLEAN telling you whether these rxns are photolysis or not 
                         (Because this triggers extra J-value mapping functionality if so!). 
                 
        (3) jmap_dict - (ONLY USED if is_photo=True) DICT with info on how to map j-values! 
                          
        (4) fjx_df    - (ONLY USED if is_photo=True) Pandas Df outputted from create_jmap_dict() 
    
        (5) make_rate_dict - (OPTIONAL) BOOLEAN if you want to make/output the rxn_rate dictionary! 
   
    OUTPUTS: 
    --------  
        (1) rxns       - LIST of reactions, formatted as F0AM expects ready for build_mech(). 
        
        (2) rxn_index  - LIST of indexes telling you where these rxns appeared in the original KPP File. 
        
        (3) rates      - LIST of rates for each reaction, formatted as F0AM Expects, ready for build_mech(). 
       
        NOTE: These 3 output lists should be equal length: rxns, rxn_index, and rates. 
        
        (4) unq_functs - LIST of unique rate function names that must be defined in F0AM's GEOSChem_J.m' 
                         or 'GEOSChem_K.m' file for the resulting mechanism to work! 
                             
     REQUIREMENTS: 
     -------------
         LIBRARIES:           
             
         CUSTOM FUNCTIONS:    _append_to_ndict_list()
                             convert_FJX_PhotoRxns_to_MCM()
     
     USAGE: 
     -----  
         CALLED WITHIN: make_GEOSCHEM_AllRxns_file()
             
         OUTPUT USAGE: 
                    
    AUTHOR 
    -------
         Prof. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD
         
    CHANGE LOG: 
    ----------
         01/14/2022    JDH Created   
         10/03/2023    JDH Modified documentation & changed to dict inputs, added FJX input parsing.   
    """ 
    
    # Add comments about each rxn back to rxn string & encase in "'" that we require to call "build mech" 
    rxns=["'"+rxn_dict['rxn'][i]+"'; % "+rxn_dict['comments'][i] if  
          rxn_dict['comments'][i] != '' else "'"+rxn_dict['rxn'][i]+"';" for i in range(0, len(rxn_dict['rxn']))] 
     
    # Get the index where this reaction should go (in totoal): 
    rxn_index= rxn_dict['origin_index']
    
    # If you're not parsing photolysis reactions... 
    if is_photo is False: 
        # Make reaction ates have "k(:,i) = " in front. As the Build-Mech function expects.  
        rates=['k(:,i) = '+rr+';' for rr in rxn_dict['rate']]; 
        
        # Get a list of all unique functions called by the rates that must be 
        # defined in F0AM for out new mech to work!  
        unq_functs=np.unique(rxn_dict['rate_function']).tolist() 
        
        if make_rate_dict is True: 
            # Define a dictionary with each unique rate function "Name" as a key with values corresponding
            # to whatever reactions that rate is used for & info about the rate/args! (only for Gas+M rxns)
            rxn_rate_dict=dict({})
            for i in range(0,len(rxn_dict['rate'])): 
                
                # Only add to rxn_rate dict things with defined functions (aka that aren't HET(IND), PHOTOL(IND) or a #.)
                if all([baddie not in rxn_dict['rate_function'][i] for baddie in ['None','HET','PHOTOL']]):
                    
                    # Make the "key of the dict a stripped down string of what rate funct it references. 
                    kname_i=str(rxn_dict['rate'][i].strip()); 
                    if kname_i not in list(rxn_rate_dict.keys()): 
                        rxn_rate_dict[kname_i]=dict({'rxns':list([]),'rate_funct':list([]), 'rate_args':list([])});
                    # Add info about rxns, rate function, and rate args used to dict: 
                    rxn_rate_dict=_append_to_ndict_list(rxn_rate_dict, kname_i,'rxns',[rxn_dict['rxn'][i]], allow_dupes=False)
                    rxn_rate_dict=_append_to_ndict_list(rxn_rate_dict, kname_i,'rate_function',[rxn_dict['rate_function'][i]], allow_dupes=False)
                    rxn_rate_dict=_append_to_ndict_list(rxn_rate_dict, kname_i,'rate_args',[rxn_dict['rate_args'][i]], allow_dupes=False)
                
            return rxns, rxn_index, rates, unq_functs,rxn_rate_dict
        else: 
            return rxns, rxn_index, rates, unq_functs
        
    else: # If you ARE parsing photolysis reactions ... 
    
        # Map all FJX photolysis rates to their corresponding values in F0AM & get a 
        # list of all unique J values that must be defined in F0AM for mech to work! 
        # Includes appending 'k(:,i) = ' and appending ";" as is done for non-photolysis reactions! 
        rates, unq_functs = convert_FJX_PhotoRxns_to_MCM(rxn_dict, jmap_dict, fjx_df)
        
        return rxns, rxn_index, rates, unq_functs

#------------------------------------------------------------------------------
#-------------- (A.2) parse_kpp_main() & Sub-Functions -------------------------
#------------------------------------------------------------------------------

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
        Used to ensure that all lists in dictionary "rinfo" that contain info 
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
    list_lens= [len(list_dict[list_i]) for list_i in list(list_dict.keys())] 
    
    # If you have more than 1 unique value of a list length, then throw an error
    if len(np.unique(list_lens)) > 1: 
        
        # Pull out the name of the function & line number where this error was found! 
        err_function=sys._getframe(1).f_code.co_name # Name of parent function where error check initiated. 
        err_line=str(sys._getframe( 1 ).f_lineno) # Line # of parent function error check initiated
        
        # Make a string that will list all the lengths of the lists using their dict keys as names... 
        err_display='\n'.join(['    Len('+key+') = '+str(list_lens[i]) for i, key in enumerate(list(list_dict.keys()))])
        
        raise ValueError('Lists are not the same length in '+err_function+'() on Line # ' 
                         +err_line+': \n'+err_display).with_traceback(sys.exc_info()[2])
        
    return 

# A.2.3.3 L3-function used within parse_rxns_and_rates() within parse_kpp_main() within make_GEOSCHEM_AllRxns_file()
def _check_rinfo(rinfo:dict, check_lens:bool=False, ignore_extras: bool=False): 
    """Function to check if rinfo has all the right key/value pairs functions expect!
    Option to check the list lengths in rinfo as well... Returns nothing if "rinfo" is good!
    
    INPUTS: 
    -------
        (1) rinfo         - DICT with list names as key and corresponding lists as value. 
        
        (2) check_lens    - BOOL indicating wehther or not to check that the lists in the 
                            DICT, rinfo, stored as values are indeed the same length. 
                          
        (3) ignore_extras - BOOL of whether you'd like to check rinfo but ignore any 
                            extra keys you weren't expecting... Useful if rinfo gets modified.'
        
    OUTPUT: 
        (1) NONE if rinfo is all good and lists are same length, otherwise pops an error. 
        
    REQUIREMENTS: 
    -------------
         LIBRARIES:           import sys
             
         CUSTOM FUNCTIONS:    check_list_lens() (defined above)
     
    USAGE: 
    -----  
    Used to ensure that "rinfo" has all the keys we expect everywhere its used 
    in this function. Optionally will check the lengths of lists in "rinfo". 
    Nice to ensure that I haven't messed up the keys in one function or another.
    
         CALLED WITHIN: parse_rxns_and_rates() within parse_kpp_main() within make_GEOSCHEM_AllRxns_file()
             
    AUTHOR: 
    -------
        Prof. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD
    
    CHANGE LOG: 
    ----------
        09/28/2023    JDH Created 
        10/31/2023    JDH Moved to make_GCAllRxns.py and updated documentation. 
    
    """ 
    
    # Define a list of the expected keys in rinfo... 
    expected_keys=['rxn','reactants','rct_stoich','products','prd_stoich','rxn_is_gas', #'rxn_in_strat',
                   'rxn_is_het','rxn_is_photo','rxn_with_M','rate','rate_function',
                   'rate_args','comments'] 
    
    # Figure out what keys in rinfo are extra or missing... 
    missing=[x for x in expected_keys if x not in list(rinfo.keys())]
    extra=[x for x in list(rinfo.keys()) if x not in expected_keys]
    
    # If you have Missing and/or extra keys in your rinfo... pop an error! 
    if any([True if len(listy)!=0 else False for listy in [missing, extra]]):
        
        # Pull out the name of the function & line number where this error was found! 
        err_function=sys._getframe(1).f_code.co_name # Name of parent function where error check initiated. 
        err_line=str(sys._getframe( 1 ).f_lineno) # Line # of parent function error check initiated
        
        if len(missing)!=0 and len(extra)==0: # Only missing some keys... 
            raise ValueError('In '+err_function+'() on Line # ' +err_line+
                         ', the dictionary, rinfo, is missing the following required keys: \n    '+ 
                         ','.join(missing)+'\n').with_traceback(sys.exc_info()[2])
        elif len(extra)!=0 and len(missing)==0 and ignore_extras is False: # Only have extra keys... 
            raise ValueError('In '+err_function+'() on Line # ' +err_line+
                         ', the dictionary, rinfo, contains the following extra keys: \n    '+ 
                         ','.join(extra)+'\n').with_traceback(sys.exc_info()[2])
        elif len(missing)!=0 and len(extra)!=0 and ignore_extras is False : # Have missing keys and extra keys! 
            raise ValueError('In '+err_function+'() on Line # ' +err_line+
                         ', the dictionary, rinfo, is missing the following required keys: \n    '+ 
                         ','.join(missing)+'\n     and has the following extra keys: \n    '+ 
                         ','.join(extra)+'\n').with_traceback(sys.exc_info()[2])
            
    # Optionally, check if all the lengths of lists in rinfo are consistent!       
    if check_lens is True: 
        _check_list_lens(rinfo)
        
    return 

# A.2.3.2 L3-function used within parse_rxns_and_rates() within parse_kpp_main() within make_GEOSCHEM_AllRxns_file()
def _fortran2matlab_nums(line): 
    """Function to convert scientific notation & # types in a line from Fortran 
    to MATLAB compliant #s & will remove spaces.
    
    INPUTS: 
    ------  
        (1) line - STR of line containing some numerical FORTRAN notation 
        
    OUTPUTS: 
    --------
        (2) line - updated STR of line with FORTRAN math converted to MATLAB math. 
                   NOTE: outputted line does not contain amy spaces! 
    
    USAGE: 
    ------ 
        Called Within: parse_rxns_and_rates() within parse_KPP_main() within make_GEOSCHEM_AllRxns_file()
        
        Output Usage: Used in output of prase_kpp_main() for rinfo['rates']
        
    REQUIREMENTS: 
    -------------
        LIBRARIES:           None
        
        CUSTOM FUNCTIONS:    _str_multi_replace()
        
    AUTHOR: 
    -------
        Prof. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD
    
    CHANGE LOG: 
    ----------
    09/28/2023    JDH Created 
    """
    
    # Dictionary with keys in Fortran code and values in MATLAB code to 
    # Convert numbers/math used in Fortran Rates to MATLAB compliant values. 
    fortran2MATLAB=dict({
        'E+00': '', # Numbers raised to a 0 should be just #s in MATLAB. 
               'e00': '',
               'E0': '', 
               'd0': '', 
               '_dp': '',
               'E+':'e',  # for numbers raise to a postivie power
               'd+': 'e', 
               'E-':'e-', # for numbers raised to a negative power.
               'd-': 'e-', 
               'EXP(':'exp('
               }) 
    
    # Replace all keys in dict that appear in line with their values & strip spaces! 
    line=utils.str_multi_replace(line, fortran2MATLAB).strip() 
    
    return line 

# A.2.3.1 L3-function used within parse_rxns_and_rates() within parse_kpp_main() within make_GEOSCHEM_AllRxns_file()
def _sep_stoich_vs_tracer(rxn_str:str, seps:list=['+','=', '->', '-->', '<->'] ):
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
        
        Output Usage:  Placed in output dict 'rinfo' of parse_rxns_and_rates() 
                       within parse_kpp_main() within make_GEOSCHEM_AllRxns_file() 
           
    AUTHOR: 
    -------
        Prof. Jessica D. Haskins (jessica.haskins@autah.edu) GitHub: @jhaskinsPhD
    
    CHANGE LOG: 
    ----------
    01/18/2022    JDH Created for pyF0AM
    09/28/2023    JDH Copy/pasted from pyF0AM to GC_Emulator/parse_KPP, updated documentation.
    """
    
    ls=[]; cmpds=[]; stoich=[] 
    for sep in seps:  # Split input string at the input seperators: 
        ls=ls+rxn_str.replace(' ' , '').split(sep)   # ls =['TOLU','OH','TRO2','1.920CH2O','0.260GLYX','0.215MGLY','OH'] 
    
    for grp in ls: # Loop over all groupings in the list split at seps
        if len(grp) > 0:
            if grp[0].isalpha(): # If str begins with a letter, then its a tracer name (e.g. 'TOLU'). 
                cmpds.append(grp); stoich.append(1) # Save the compound, save the stoichiometry. 
            
            else: # Otherwise, loop over the chars in the grp until you find a letter. 
                yld=grp[0]; ind=1; found_yld=False # Set stoichiometry= to first char (not a letter!)
                
                while found_yld==False:  # loop til you find a letter. 
                    if grp[ind].isnumeric() or grp[ind]=='.': # Found a letter or decimal place. 
                        yld=yld+grp[ind];  # Add to str containing all stoichiometry. 
                    elif grp[ind].isalpha() or  ind==len(grp)-1: # Found beginning of compound name. 
                        cmpds.append(grp[ind:]) # Store it 
                        stoich.append(np.float64(yld)) # and the #s making up the yield as a #
                        found_yld=True
                        break
                    ind=ind+1
                    
    # Don't let 'hv' be a "compound" or have "stoichimetry"... 
    if 'hv' in cmpds: 
        stoich.pop(cmpds.index('hv') ); cmpds.remove('hv')
    
    return cmpds, stoich

# A.2.4 L2-function used within parse_kpp_main() within make_GEOSCHEM_AllRxns_file()
def _enforce_list_dict_inds(list_dict:dict, inds:list,  add_index:bool=False):
    """Return a dictionary containing only a sub-selection from the lists stored
    as values within the input dictionary, list_dict, using a list of the indices
    that'd you like to sub-select. Useful for splitting info about all reactions
    into individual lists pertaining only to reactions of 1 type! 
    
    INPUT: 
    -------
        (1) list_dict - DICTIONARY with equal length lists stored as values you'd like to 
                        sub-select items from.  
                
        (2) inds      - LIST containing indices of the values you want to sub-select from in list_dict.  
        
        (3) add_index - (OPTIONAL) Boolean of whether to add index used to sub-select data 
                         from list_dict within the output dictionary! Default is FALSE
    OUTPUT: 
    --------
        (1) out_dict - DICT with only the values at indices inputted stored as values.
                        and may contain a new key with orgin indexes in it. 
    
    REQUIREMENTS: 
    -------------
        LIBRARIES:           None
            
        CUSTOM FUNCTIONS:    check_rinfo()  (defined above)
                
    AUTHOR: 
    -------
        Prof. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD
         
    CHANGE LOG: 
    ----------
        01/14/2022    JDH Created for pyMCM 
        09/28/2023    JDH Modified to take input of dict instead of a list of lists for GC_Emulator. 
        10/02/2023     JDH Added "add_index" option to function
    """ 
    # Create a new dictionary with the same keys as the input "list_dict" containing empty lists! 
    out_dict=dict({})
    
    # Loop over all keys in input "list_dict". 
    for key in list(list_dict.keys()): 
        # Pull out the full "master_list" contained under this key. 
        master_list=list_dict[key]
        
        new_list=[]; # Create a NEW List to hold only items with index in the input "inds".
        
        # Fill the new list from the master list where the index matches "inds" 
        new_list= [master_list[i] for i in inds]
        
        # Save this new list in the output dictionary! 
        out_dict[key]=new_list 
        
        # Clear local variables so you don't mess up the next key/list pair. 
        del new_list, master_list
        
    # Check that all the keys are still there and that all lists are the same length! 
    _check_rinfo(out_dict, check_lens=True, ignore_extras=False)
    
    # Now store the index used to sub-select these values from list_dict if asked! 
    if add_index is True: 
        out_dict['origin_index']=inds 

    return out_dict

# A.2.3 L2-function used in parse_kpp_main() called within make_GEOSCHEM_AllRxns_file():
def _is_blank(line): 
    """Function to determine if a line is blank or not after all spaces have been removed."""
    tf=True if len(line.replace(' ',''))==0 else False
    return  tf

# A.2.3 L2-function used in parse_kpp_main() called within make_GEOSCHEM_AllRxns_file():
def parse_rxns_and_rates(line:str, rinfo:dict, file_location:dict): 
    """Function used to parse reaction lines in KPP file and extract their info into dictionary, rinfo. 
    
    INPUTS:
    -------
    
        (1) line - STR of line with reaction, rate, and (maybe) a comment. These lines look like this:     
                      Reaction               Rate                              Comment
                "O3 + NO = NO2 + O2 :  GCARR_ac(3.00d-12, -1500.0d0); {2014/02/03; Eastham2014; SDE}"
        
        (2)  rinfo - DICT with the following Key/Value Pairs.
                KEYS           :           VALUES
                ----------------------------------------------------------------------
                rxn            : List of strings with full reaction only... 
                reactants      : *Nested* list containing all reactants per reaction.  
                rct_stoich     : *Nested* list of reactants stoichiometry per reaction.  
                products       : *Nested* list containing all products per reaction. 
                prd_stoich     : *Nested* list of products stoichiometry per reaction.  
                rxn_is_het     : Boolean List of whether rxn is heterogeneous or not. 
                rxn_is_photo   : Boolean list of whether rxn is a photolysis rxn or not. 
                rxn_in_strat   : Boolean list of whether rxn is a stratosphere only rxn or not. 
                rxn_with_M     : Boolean list of whether rxn has a (now hidden) "+M" reactant or not.
                rate           : List of strings with full rate for each reaction 
                rate_function  : List of strings with rate function for each reaction
                rate_args      : *Nested* list with arguments to rate function for each reaction 
                comments       : List of comments for each reaction 
                
                NOTE: All lists in the dict, rinfo, should be the same length upon input and output to 
                    this function, as they all are all "indexed "by the reaction they refer to. 
                 
        (3) file_location - Dictionary with the following Key/Value Pairs: 
        
                KEYS           :           VALUES
                ----------------------------------------------------------------------
                past_header    : Boolean. True if we have passed the header, otherwise False.   
                in_defvar      : Boolean. True if currently parsing the fortran variable defintion section
                in_deffix      : Boolean. True if currently parsing the "#deffix" section  
                in_gas_rxns    : Boolean. True if currently parsing Gas Phase reactions
                in_het_rxns    : Boolean. True if currently parsing Heterogeneous reactions
                in_photo_rxns  : Boolean. True if currently parsing Photolysis reactions
     
     OUTPUTS: 
     -------   
         (1)  rinfo       - Updated DICT with the same Key/Value Pairs as above, but outputted 
                             lists in rinfo will be +1 in len since they will now include info 
                             about the reaction contained on the input line. This function checks 
                             for len consistency @ end. 
                     
        (2) file_location - Updated DICT  with the same Key/Value Pairs as above. 
    
    REQUIREMENTS: 
    -------------
        LIBRARIES:           None
            
        CUSTOM FUNCTIONS:    check_rinfo()
                            _str_multi_replace()
                            sep_stoich_vs_tracer()
                            fortran2matlab_nums()
                            check_list_lens()
    
    USAGE: 
    -----  
        Called Within: parse_kpp_main() within make_GEOSCHEM_AllRxns_file(). 
        Output Usage: Output dict 'rinfo' used in parse_kpp_main() to sep info for dif rxn types into own dicts! 
        
    AUTHOR: 
    -------
        Prof. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD
            
    CHANGE LOG: 
    ----------
        10/03/2023    JDH Created
        
    """
    # Check that the input dictionary, rinfo, has the required input keys & all lists are consistent lengths 
    _check_rinfo(rinfo, check_lens=True)

    #####################################################################################
    # Split Line in to reaction, rate, rate function, and comment & Clean up for MATLAB! 
    #####################################################################################
    rxn_parts =line.split(':')  # Split on colon to get rxn only
    rate=rxn_parts[1].split(';')[0]  #Split on ';' to get rates 
    
    # If this part of the line after the reaction has a comment, then split that on'{' and clean it up.  
    if '{' in rxn_parts[1]: 
        replace_me=dict({'{':'', '}': '', ';':',', '\n':''}) 
        comment=utils.str_multi_replace(rxn_parts[1].split('{')[1], replace_me) 
    else: comment= '' # Otherwise make comment blank! 
    
    # Update rinfo to keep track of what kind of reaction this is (photolysis, het, strat, with M, etc.) 
    rinfo['rxn_is_gas'].append(np.int64(file_location['in_gas_rxns']))     # 1 if it is a gas phase rxn, 0 if not.   
    rinfo['rxn_is_photo'].append(np.int64(file_location['in_photo_rxns'])) # 1 if it is a photolysis rxn, 0 if not.   
    rinfo['rxn_is_het'].append(np.int64(file_location['in_het_rxns']))     # 1 if it is a het rxn, 0 if not.      
    #rinfo['rxn_in_strat'].append(np.int64(is_strat_rxn))                  # 1 if it is a stratosphere only rxn, 0 if not.          
    rinfo['rxn_with_M'].append(1)  if '{+M}' in rxn_parts[0] else rinfo['rxn_with_M'].append(0)

    # Now remove spaces, and '{+M'}, '{' , '}' characters from the reaction
    rxn= utils.str_multi_replace(rxn_parts[0], ['{+M}', '{', '}',' '], rep_all='').strip()  
    
    # Split entire reaction into Reactants & Products & their stoichiometry:  
    rct_i, rct_stoich_i= _sep_stoich_vs_tracer(rxn.split('=')[0], seps=['+']) 
    prd_i, prd_stoich_i= _sep_stoich_vs_tracer(rxn.split('=')[1], seps=['+']) 
    
    # Clean up rates so the #s appear how they should in MATLAB & strip spaces from rate line.
    rate= _fortran2matlab_nums(rate) # (e.g. convert 1300.0d --> 1300.0)
    
    # Split rate into function and args!  Lines look like this: 'GCARR_ac(2.60e-13, 1300)'
    if '(' in rate and ')' in rate:
        # Before you do a string split to extract the function, find all the matches to "exp()"
        exp_groups=re.findall(r'[0-9+-e./*]*exp\(.+\)',rate)
    
        # If there is a match in the string for "exp", remove it and check for "(" ")" again! 
        if len(exp_groups)>0: 
            # Remove all exp stuff from the rate! 
            rate_rev=rate
            for g in exp_groups: 
                if g in rate_rev:
                    rate_rev=rate_rev.replace(g,'')
                    
            # If after replacing these items, there are still parenthesis, then split it TO GET FUNCTION, ARGS
            if '(' in rate_rev and ')' in rate_rev and len(rate_rev)>0: 
                function=rate_rev.split('(')[0] # function= 'GCARR_ac'
                arg_list=rate_rev.split('(')[1].split(')')[0].split(',') # arg_list= list([2.60e-13, 1300])
            else: 
                function='None'; arg_list=list([])
        else: 
            function=rate.split('(')[0] # function= 'GCARR_ac'
            arg_list=rate.split('(')[1].split(')')[0].split(',') # arg_list= list([2.60e-13, 1300])
    else:
        function='None'; arg_list=list([])
     
    # Turn rate function calls into strings! 
    # If the first char is a letter (not a #!) 
    if rate[0].isalpha() is True and 'PHOTOL' not in rate and 'HET' not in rate: 
        # Then turn it into a STRING that will be recognized by MATLAB... 
        rate="'K_"+rate+"'"
             
             
    # Now save everything into its appropriate output list  
    rinfo['rxn'].append(rxn)                 # List of strings with full reaction only... 
    rinfo['reactants'].append(rct_i)         # *Nested* list containing all reactants per reaction.  
    rinfo['rct_stoich'].append(rct_stoich_i) # *Nested* list of reactants stoichiometry per reaction.  
    rinfo['products'].append(prd_i)          # *Nested* list containing all products per reaction. 
    rinfo['prd_stoich'].append(prd_stoich_i) # *Nested* list of products stoichiometry per reaction.  
    rinfo['rate'].append(rate)               # List of strings with full rate for each reaction 
    rinfo['rate_function'].append(function)  # List of strings with rate function for each reaction
    rinfo['rate_args'].append(arg_list)      # *Nested* list with arguments to ratre function for each reaction 
    rinfo['comments'].append(comment)        # List of comments for each reaction 

    # Check that all new lists in rinfo are still the same length! 
    _check_list_lens(rinfo)
    
    return rinfo

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
                
    """    
    
    # Read more lines til you get the entire reaction, rate & comments in a single line!  
    while line_stripped.endswith('+')==True:  
        nline = kppfile.readline()   # Read next line from file 
        count=count+1                 # Update line counter variable  
        line = line+ nline.strip()    # append this new line to our current line str.  
        line_stripped=line.strip().replace(' ', '') # Update line_stripped to not contain spaces! 
    
    # Now that you got the entire reaction on this line, make sure to replace any new line characters with a space! 
    line=line.replace('\n', ' '); line_stripped=line_stripped.replace('\n', ' ') 
    
    # Return the count, new line, and line_stripped vars to continue on. 
    return line, line_stripped, count

# A.2.1 L2-function used in parse_kpp_main() called within make_GEOSCHEM_AllRxns_file():
def parse_var_defs(line:str, tracers:list =[], long_name:list=[],split_on:str='ignore'): 
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
        ind=line.lower().index(split_on) 
        
        # Remove spaces and "=" from line, to get only the tracer name being defined! 
        tracer=utils.str_multi_replace(line[:ind],  [' ', '='], rep_all='') 
        
        tracers.append(tracer)  # Keep tracer name of this compound defined in a list.  
        
        # Pull out the rest of the comment about this tracer using the index. 
        tracer_info=line[ind+len(split_on):] 
        
        # Remove any remaining weird characters we don't care about (";", "{", "}"). 
        if tracer_info[0]==';': tracer_info=tracer_info[1:]  
        tracer_info= utils.str_multi_replace(tracer_info,  ['{', '}' ], rep_all='').strip() 
        
        #Now add the "long name" slash comment about this tracer to a list! 
        long_name.append(tracer_info) 
        
    return tracers, long_name

# A.2 L1-function called within make_GEOSCHEM_AllRxns_file() 
def parse_kpp_main(kppfile): 
    """Function to parse a KPP file and extract info about tracers, reactions & rates.  
     
    INPUTS:  
    -------  
       (1)  kppfile - STR containing the full path and name of a KPP file.  
     
 
    OUTPUTS:  
    -------- 
    
    
    REQUIREMENTS: 
    -------------
        LIBRARIES:           None
            
        CUSTOM FUNCTIONS:    parse_var_defs()
                            _is_blank()
                            fix_multiline_rxns()
                            check_rinfo()
                            parse_rxns_and_rates()
                            _enforce_list_dict_inds()
    
    USAGE: 
    -----  
        CALLED WITHIN: 
            
        OUTPUT USAGE: 
    
    AUTHOR: 
    -------
       Prof. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD
             
    CHANGE LOG:  
    ----------
       01/14/2022   JDH Created   
       10/03/2023   JDH modified documentation, moved FJX handling to its own function. 
        
    """ 
    ###########################################################################
    # Initialize all vars to hold outputs/ keep track of where we are!  
    ###########################################################################
    
    # Get # of lines to expect so we know when we've reached the end of the file! 
    num_lines = sum(1 for line in open(kppfile))  
    tracers=[]; long_name=[]; # Make empty lists to hold tracers & tracer info!  
    count = 0  # Initialize line counter variable. 
    
    # Define boolean dictionary with key/values that will tell us what part of the KPP file we are currently parsing.   
    file_location=dict({
        'in_header'    :True,  # Boolean. True on initiation & if currently parsing header. 
        'in_defvar'    :False, # Boolean. True if currently parsing the Fortran variable definition section   
        'in_deffix'    :False, # Boolean. True if currently parsing the "#deffix" section  
        'in_gas_rxns'  :False, # Boolean. True if currently parsing gas phase reactions  
        'in_het_rxns'  :False, # Boolean. True if currently parsing heterogenous reactions  
        'in_photo_rxns':False  # Boolean. True if currently parsing photolysis reactions
        })
    
    # Define  dictionary to hold lists of info about dif rxns! 
    rinfo=dict({
        'rxn':list([]),           # List of strings with full reaction only... 
        'rxn_is_gas':list([]),    # Boolean list of whether rxn is gas phase or not. .
        'rxn_is_het':list([]),    # Boolean List of whether rxn is heterogeneous or not. 
        'rxn_is_photo':list([]),  # Boolean list of whether rxn is a photolysis rxn or not. 
        #'rxn_in_strat':list([]),  # Boolean list of whether rxn is a stratosphere only rxn or not. 
        'rxn_with_M':list([]),    # Boolean list of whether rxn has a (now hidden) "+M" reactant or not.
        'reactants':list([]),     # NESTED list containing all reactants per reaction. 
        'rct_stoich':list([]),    # NESTED list of reactants stoichiometry per reaction.  
        'products':list([]),      # NESTED list containing all products per reaction. 
        'prd_stoich':list([]),    # NESTED list of products stoichiometry per reaction.  
        'rate':list([]),          # List of strings with full rate for each reaction 
        'rate_function':list([]), # List of strings with rate function for each reaction
        'rate_args':list([]),     # NESTED list with arguments to rate function for each reaction 
        'comments':list([])       # List of comments for each reaction 
        })
    
    ##############################################################################
    # Open the KPP file in "read" mode & Loop over all lines until end is reached!
    ##############################################################################
    file = open(kppfile, 'r')  
    
    while True:  # loop over each line of the KPP file. 
        line = file.readline() # Read next line from file 
        line_stripped=line.strip().replace(' ', '') # Strip all white space from current line. 
        count=count+1 # Update line counter variable now that you read in this line. 
        
        #--------------------------------------------------------------------------------------------------------------------------
        # BEFORE parsing the line, set file location booleans to FALSE if the current line indicates you're entering a new section!
        #--------------------------------------------------------------------------------------------------------------------------       
        # You are no longer in the header section once #defvar appears on the line! 
        file_location['in_header']   = False if '#defvar' in line_stripped.lower() else file_location['in_header']   
        # You are no longer in the variable defintion section if #deffix appears in the current line! 
        file_location['in_defvar']   = False if '#deffix' in line_stripped.lower() else file_location['in_defvar'] 
        # You are no longer in section defining fix vars if #equations appears in the current line!   
        file_location['in_deffix']   = False if '#equations' in line_stripped.lower() else file_location['in_deffix'] 
        # You are no longer in the gas phase reaction section once both '//' and 'heterogeneous' appear in the current line!
        file_location['in_gas_rxns'] = False if all([item_i in line_stripped.lower()  for item_i in ['//','heterogeneous']]) else file_location['in_gas_rxns']  
        # You are no longer in the photolysis reaction section once both '//' and 'photolysis' appear in the current line! 
        file_location['in_het_rxns'] = False if all([item_i in line_stripped.lower()  for item_i in ['//','photolysis']]) else file_location['in_het_rxns']  
           
        # Break if reached end of file, don't attempt to parse line. 
        if count>num_lines: break  
    
        #------------------------------------------------------------------------------------------
        # If you're in the part of the file where variable definitions are made ==> parse_var_defs()
        #------------------------------------------------------------------------------------------
        if file_location['in_defvar']==True or file_location['in_deffix']==True:  
            
            # Pass to parse_var_defs() & return updated lists with tracer name and long name.
            tracers, long_name = parse_var_defs(line, tracers=tracers, long_name=long_name, split_on='ignore'); 
         
        #-----------------------------------------------------------------------------------------    
        # If you're in part of the file with chemical equations & rates ==> parse_rxns_and_rates()
        #-----------------------------------------------------------------------------------------
        if any(cond ==True for cond in [file_location['in_gas_rxns'], file_location['in_photo_rxns'], file_location['in_het_rxns']]):  
            
            # Don't parse commented out or blank lines.  
            if (('//'!= line[0:2]) and (_is_blank(line)==False)): 
                
                # Decide if you've been given a line that contains a partial reaction & fix it before parsing as rxn.  
                if ':' not in line: 
                    # These lines won't contain a ':' and look like this: 
                    #   ln #     Content of Line      
                    #-------------------------------------------------------------
                    #   ln 42   'MVKOHOO + HO2 = 0.360MCO3 + 0.360GLYC +'   **(Would trigger when parsing this line)** 
                    #   ln 43   '  0.335MVKHP + 0.050MGLY + 0.050CH2O :  GCARR(2.12E-13, 0.0E+00, 1300.0); {2019/11/06; Bates2019; KHB}' 
                    line, line_stripped, count= fix_multiline_rxns(line, line_stripped, count, file)
            
                # (After fixing if necessary), parse this reaction and add all info about this rxn to lists: 
                rinfo=parse_rxns_and_rates(line, rinfo, file_location)

        #--------------------------------------------------------------------------------------------------------------------------
        # AFTER parsing the line, set file location booleans to TRUE if the current line indicates you're entering a new section!
        #--------------------------------------------------------------------------------------------------------------------------
        # You will be in the variable definition section on the next line if '#defvar' appears in the current line. 
        file_location['in_defvar']     = True if '#defvar' in line_stripped.lower() else file_location['in_defvar']
        # You will be in the deffix section on the next line if #deffix appears in the current line. 
        file_location['in_deffix']     = True if '#deffix' in line_stripped.lower() else file_location['in_deffix']
        # You will be in the gas-phase reaction section on the next line if #equations appears in the current line. 
        file_location['in_gas_rxns']   = True if '#equations' in line_stripped.lower() else file_location['in_gas_rxns'] 
        # You will be in the heterogenous reaction section on the next line if '//' and 'heterogenous' appear in the current line. 
        file_location['in_het_rxns']   = True if all([itemi in line_stripped.lower() for itemi in ['//','heterogeneous',]]) else file_location['in_het_rxns']
        # You will be in the photolysis reaction section on the next line if '//' and 'photolysis' appear in the current line. 
        file_location['in_photo_rxns'] = True if all([itemi in line_stripped.lower()  for itemi in ['//','photolysis']]) else file_location['in_photo_rxns'] 
             
    # While loop will break once we found the end of the KPP file and have parsed all lines. 
    file.close() # Close the KPP file. 
     
    # Build Var info, a dictionary connecting tracer names & long names. 
    var_info= dict(zip(tracers, long_name)) 
     
    #--------------------------------------------------------------------------
    # Seperate out reactions from rinfo into dif dicts for each type of rxn. 
    #--------------------------------------------------------------------------
    # First, just make sure rinfo has right keys & consistent list lengths! 
    _check_rinfo(rinfo, check_lens=True, ignore_extras=False)
    
    # Build list of Indexes in lists under rinfo that correspond to dif rxn types: 
    gas_inds    = [i for i,grx in enumerate(rinfo['rxn_is_gas'])   if grx==1] 
    het_inds    = [i for i,hrx in enumerate(rinfo['rxn_is_het'])   if hrx==1] 
    nhet_inds   = [i for i,hrx in enumerate(rinfo['rxn_is_het'])   if hrx==0] 
    photo_inds  = [i for i,jrx in enumerate(rinfo['rxn_is_photo']) if jrx==1] 
    M_inds      = [i for i,mrx in enumerate(rinfo['rxn_with_M'])   if mrx==1] 
    #strat_inds = [i for i,srx in enumerate(rinfo['rxn_is_strat']) if srx==1] 

    # Use indices to sep info about each dif reaction type into their own dict with the same format as rinfo. 
    gas_dict=  _enforce_list_dict_inds(rinfo, gas_inds,   add_index=True)
    het_dict=  _enforce_list_dict_inds(rinfo, het_inds,   add_index=True)
    nhet_dict= _enforce_list_dict_inds(rinfo, nhet_inds,  add_index=True)
    photo_dict=_enforce_list_dict_inds(rinfo, photo_inds, add_index=True)
    Mrxn_dict= _enforce_list_dict_inds(rinfo, M_inds,     add_index=True)
    gasM_dict=_enforce_list_dict_inds(rinfo, np.sort(gas_inds+M_inds),add_index=True)
    #start_dict=_enforce_list_dict_inds(rinfo, strat_inds,add_index=True)
    
    return var_info, rinfo, gas_dict, het_dict, nhet_dict, photo_dict, Mrxn_dict, gasM_dict, # strat_dict

#------------------------------------------------------------------------------
#-------------- (A.1) create_jmap_dict() & Sub-Functions -----------------------
#------------------------------------------------------------------------------
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
    
    ct =1; # Initialize counter variable
    
    while True: 
        if jname in current_names: 
            # Pull the "name" of the thing out. For "jCH3COOH", nm_only='CH3COOH'
            nm_only=re.findall(r'j([a-zA-Z0-9]+)',jname)
            if len(nm_only)==0: 
                print(jname);sys.exit()
            # Decide on a new suffix to try using (iterates through alphabet).
            new_let='_'+chr(ord('`')+ct)
            jname='j'+nm_only[0]+new_let
            
        if jname not in current_names: 
            break
        else: 
            ct=ct+1  # Update counter var to try a new '_'+LETTER(ct) 
   
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
    fjx_df= pd.read_fwf(file, skip_rows=[1,2], widths=[5,10,10,34,6,8]) 
    
    fjx_df= fjx_df.drop(index=len(fjx_df)-1); # Drop nan column. 
    
    names=list(fjx_df.columns);  # Get original column names. 
    
    # Rename all the columns to more intuitive names. 
    fjx_df=fjx_df.rename(columns={names[0]:'FJX_Index',     # See above
                            names[1]:'Photolysis_Of', # Name of the GEOS-Chem species that will undergo photolysis.
                            names[2]:'PHOTON',        # Indicator that the reaction is a photolysis rxn.
                            names[3]:'Into',          # Products of this specific photolysis reaction ...
                            names[4]:'Quantum_Yield', # See above
                            names[5]:'FJX_CrossSec'}) # See above
        
    # Join parts of the columns together to make a coherant photolysis reaction 
    fjx_df['Photolysis_Rxn']= [str(fjx_df.loc[ind, 'Photolysis_Of'])+ '+hv -> '+'+'.join( [prd.replace('...','') for prd in str(fjx_df.loc[ind, 'Into']).split(' ') if len(prd) >0]) for ind in range(0, len(fjx_df))]
    
    # Clean up cross section used column from spaces and other characters we don't need! 
    fjx_df['FJX_CrossSec']=fjx_df['FJX_CrossSec'].apply(lambda x: x.replace(' ', '').replace('/', '').replace('(', '').replace(')', '').replace('-', ''))
    
    # Stick a j in front of the cross section used to define a "nickname"... 
    fjx_df['J-NickName']= fjx_df['FJX_CrossSec'].apply(lambda x: 'j'+str(x))

    # Now drop columns that we don't actually need. 
    fjx_df=fjx_df.drop(columns=[ 'Photolysis_Of', 'PHOTON','Into'])
    
    return fjx_df

# A.1 L1-function called within make_GEOSCHEM_AllRxns_file() 
def create_jmap_dict(version:str, jmap_type='CrossSec_Match'): 
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
    photolysis_paths=utils.get_my_relative_paths(photo_only=True)
    
    # Read in dictionary defining what version of FJX files are used as CHEM_INPUTS for each dif GC version w/ gc version as key 
    vfjx=pd.read_excel(photolysis_paths['FJX_input_by_GC-version.xlsx'], index_col=0).to_dict()['FAST-JX-INPUT']
    
    if version in list(vfjx.keys()): # If you have that version, then set the path to its J2J file ... 
        j2j_file= os.path.join(photolysis_paths['FJX_Files'],'FJX_j2j_'+vfjx[version]+'.dat')
        exists=os.path.exists(j2j_file)
        if exists is True: 
            print('Using: FJX_j2j_'+vfjx[version]+'.dat for GEOS-Chem version:  ', version )
        else: 
            raise FileNotFoundError('We could not locate the correct FJX file to open for GEOS-Chem Version: '+version+
                                    ' at the following expected path:  \n' + 
                                    '    '+j2j_file).with_traceback(sys.exc_info()[2])
    else: 
        raise ValueError('Version ',version,' of GEOS-Chem Classic was not found in the dictionary \n'+
                         'defining its corresponding FJX chem input file... \n\n '+ 
                         'This likely means you are trying to Emulate a new version of GEOS-Chem not yet supported by our Emulator. \n'
                         'To proceed, you will need to figure out what FJX version is used for this version of GEOS-Chem and then: \n' +
                         '    (1) Open and edit the file: ".../GEOSChem_Emulator/Input_Files/Photolysis_Files/FJX_input_by_GC-version.xlsx" \n'+
                         "           - Add the GEOS-Chem version you're emulating in the FIRST column \n"+
                         '           - Add the version of the FJX file that GEOS-Chem version uses to the the SECOND column.' + 
                         '    (2) Add the Raw FJX.dat file to the directory: "GEOSChem_Emulator/Input_Files/Photolysis_Files/FJX_Files/" if it does not exist already. \n\n' +
                         'Why? #1 Will make sure we know what FJX file to open/read when emulating this version of GEOS-Chem in our dictionary that popped this error and \n'+ 
                         '#2 Will make sure we have the actual FJX file to open and read.').with_traceback(sys.exc_info()[2])
        
    # Read in info about what Photolysis Rxns are Defined in this FJX File into a pandas dataframe!        
    fjx= read_FJX_file(j2j_file);             
    fjx['MCM_J2Use']= list([''])*len(fjx)  # Create an empty list-like column to hold J-values... 
    
    # Get FJX-Cross Section --> MCM J-Value Mapping done in do_jval_mapping.py
    cs=pd.read_excel(photolysis_paths['FJX_cross_sect_to_jvalues.xlsx']).fillna('') 

    # Assign J-Values in mechanism based on the Cross-Section used in FAST-JX: 
    for ind in fjx.index: # All data is defined in the dataframe fjx for this matching.
        csect_need=fjx.loc[ind,'FJX_CrossSec'] # Cross-Section we need a J-Value for... 
        fjx_indx= fjx.loc[ind,'FJX_Index'] # Index of the Reaction in PHOTOL()
        
        if csect_need == '': # If You can't figure out what cross section is needed... 
            # Figure out what reaction this is for... 
            rx_i= fjx.loc[ind,'Photolysis_Rxn'].split('->')[0].replace(' ','')
            
            # Some versions don't have names for NITP cross sections... Fix that. 
            if any([nit_str in rx_i for nit_str in ['NIT+hv','NITs+hv']]) :
                csect_need='NITP'
                fjx.at[ind,'FJX_CrossSec']='NITP'
                fjx.at[ind,'J-NickName']= 'j'+fjx.loc[ind,'Photolysis_Rxn'].split('+')[0].replace(' ', '')
                
        # Get Index in dataframe, cs, where that cross-section defined: 
        has=cs.index[cs['FJX_Cross_Section']==csect_need]
        
        if len(np.unique(cs.loc[has,jmap_type]))==1: # If we have a match 
                fjx.at[ind,'MCM_J2Use']=np.unique(cs.loc[has,jmap_type])[0]# Then use that in the df. 
        
        # Fill in info about everything used in the output fjx_info dataframe! 
        if fjx.loc[ind,'Quantum_Yield'] != 1: 
            fjx.at[ind, 'F0AM_Assignment'] = str(fjx.loc[ind,'Quantum_Yield'])+'.*('+fjx.loc[ind,'MCM_J2Use']+')'
        else: 
            fjx.at[ind, 'F0AM_Assignment'] =fjx.loc[ind,'MCM_J2Use']
            
    # Get a list of all duplicate j-nicknames in FJX
    unique_Jnames=set(); dupe_Jnames=[] 
    for ji,j in enumerate(fjx.loc[:,'J-NickName']): 
        if j in unique_Jnames and j not in dupe_Jnames: 
            dupe_Jnames.append(j)
        elif j not in unique_Jnames: 
            unique_Jnames.add(j) 

    for j in dupe_Jnames: 
        idx_matches= fjx.index[j==fjx.loc[:,'J-NickName']] # ind in fjx of all matches to a dupe name 
        assigns=[fjx.loc[ix,'F0AM_Assignment'] for ix in idx_matches] 
        is_diff=[True if assignment==assigns[0] else False for assignment in assigns]
        
        if not all(is_diff):# If some of the assignments vary for the same J-Nickname...  
            for ii,ix in enumerate(idx_matches): # Loop over all mataches to this nickname
                if (is_diff[ii]==False) or (ii==0): # If it doesn't match the assignment of the one already defined (1st occurance) or is 1st occurance...
                    new=make_unique_Jname(j,list(fjx.loc[:,'J-NickName'])) # Update the name to be unique! 
                    fjx.at[ix,'J-NickName']=new
    
    for ix in fjx.index: 
        fjx.at[ix,'PHOTOL(IND)']= 'PHOTOL('+str(fjx.loc[ix,'FJX_Index'])+')'
        fjx.at[ix, 'Info'] = '% PHOTOL('+str(fjx.at[ix,'FJX_Index'])+')'+' used for: '+fjx.at[ix,'Photolysis_Rxn']+' based on '+fjx.at[ix,'FJX_CrossSec']+' crossection with QY='+str(fjx.at[ix,'Quantum_Yield'])
        fjx.at[ix,'GEOSChem_Jline']=   fjx.at[ix,'J-NickName']+'='+fjx.at[ix,'F0AM_Assignment']+'; % PHOTOL('+str(fjx.at[ix,'FJX_Index'])+') \n'+\
                                       '% Used for: '+fjx.at[ix,'Photolysis_Rxn']+' based on '+fjx.at[ix,'FJX_CrossSec']+' crossection with QY='+str(fjx.at[ix,'Quantum_Yield'])
    
    # And also create a dictionary to hold results mapping PHOTOL(IND) => J-NickName in F0AM! 
    jdict=dict(zip(fjx['PHOTOL(IND)'],["'" +j +"'" for j in fjx['J-NickName']]))
    
    return jdict, fjx

###############################################################################
#-----------------(A) CREATE THE GEOSCHEM_AllRxns.m File ----------------------
###############################################################################

def make_GEOSCHEM_AllRxns_file(kppfile:str, GC_version:str, jmap_type: str= 'CrossSec_Match', 
                include_het_rxns:bool=False, output_fname:str= '',output_dir:str= '',
                verbose:bool=True, overwrite:bool=False): 
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
    var_info, rinfo, gas_dict, het_dict, nhet_dict, photo_dict, Mrxn_dict,gasM_dict=parse_kpp_main(kppfile) 
            
    # Pull out info about dif rxn types & get it in the format we need to pass to build_mech()
    nhet_rxns, nhet_index, nhet_rates, nhet_functs= KPP_Dict_to_F0AM_IO(nhet_dict,is_photo=False)
    gas_rxns,   gas_index,  gas_rates,  gas_functs= KPP_Dict_to_F0AM_IO(gas_dict, is_photo=False)
    het_rxns,   het_index,  het_rates,  het_functs= KPP_Dict_to_F0AM_IO(het_dict, is_photo=False)
    M_rxns,      M_index,   M_rates,    M_functs= KPP_Dict_to_F0AM_IO(het_dict, is_photo=False)
    gasM_rxns, gasM_index, gasM_rates, gasM_functs,rxn_rate_dict= KPP_Dict_to_F0AM_IO(gasM_dict, is_photo=False, make_rate_dict=True)
    
    # Same as above, but also do j-value mapping to translate FJX/KPP j-values into their F0AM counterparts! 
    photo_rxns, photo_index, photo_rates, photo_functs= KPP_Dict_to_F0AM_IO(photo_dict,is_photo=True,jmap_dict=jmap_dict, fjx_df=fjx_df)
    
    rxns=gas_rxns+M_rxns+photo_rxns 
    rates=gas_rates+M_rates+photo_rates
    r_functs=gasM_functs
    
    #---------------------------------------------------------------------------------------------
    # Build your mechanism either with het rxns in their own file  or without het rxns... 
    #----------------------------------------------------------------------------------------------
    if include_het_rxns is False: 
        # Build a F0AM compliant mechanism without all the het rxns ... 
        het_info= fio.build_all_from_rxns(het_rxns, k_list=het_rates, sort_rxn_i=False, sort_rxn_list=False, verbose=verbose) 
        info = fio.build_all_from_rxns(rxns, k_list=rates, sort_rxn_i=False, sort_rxn_list=False, verbose=verbose) 
        
        # Now parse the non-het mech to remove reactions of stuff only formed via a het rxn... 
        [sps, rxns, rates,gs,fs,rctt,prdd,r_yld,p_yld], rxn2drop=eliminate_dead_ends(het_info,info)
        
        if verbose is True and len(rxn2drop) > 0: 
                print('Dropping these gas/ photolysis rxns from mech because they '+
                      'involve species only formed in het reactions, which are not included:\n')
                [print('    '+r) for r in rxn2drop]; print('\n')
            
        # Rebuild your mech without these rections/ species they form... 
        out= fio.build_all_from_rxns(rxns, k_list=rates, sort_rxn_i=False, sort_rxn_list=False, verbose=verbose)
        [sp_list, rxns, fs, gs, ks, rct_cmpds, prd_cmpds, rct_ylds, prd_ylds]= out
         
    #else: 
        # # Build a mechanism for the non= het rxns and rates ... 
        #out= fio.build_all_from_rxns(rxns, k_list=rates, sort_rxn_i=True, sort_rxn_list=False, verbose=verbose) 
        # [sp_list, rxns, fs, gs, ks, rct_cmpds, prd_cmpds, rct_ylds, prd_ylds]= out
        
        # # Build a mechanism from the het rxns and rates ... 
        # het_out= fio.build_all_from_rxns(het_rxns, k_list=het_rates, sort_rxn_i=True, sort_rxn_list=False, verbose=verbose) 
        # [hsp_list, hrxns, hfs, hgs, hks, hrct_cmpds, hprd_cmpds, hrct_ylds, hprd_ylds]= het_out
        
        # # Comment out all halogen Het Rxns... 
        # rmk=list()
        # hals=['Br','Cl','I','SALC','SALA','AERI','HI','ISALA','ISALC', 'HOCl',
        #       'HCl', 'HOI','SALACL','ICl','IONO','BrSALA','BrSALC','HOBr','HBr','Cl2',
        #       'BrCl','ClNO2','ClNO3']
        # for kk in range(0,len(hrxns)): 
        #     if any([hal in hrct_cmpds[kk]+hprd_cmpds[kk] for hal in hals]):
        #         hrxns[kk]='%'+hrxns[kk]
        #         hfs[kk]='%'+hfs[kk]
        #         hgs[kk]='%'+hgs[kk]
        #         hks[kk]='%'+hks[kk]
        #         rmk.append(hks[kk].split('=')[1].split('(')[0].strip())
                
        # het_rfuncts=[h for h in het_rfuncts if h not in rmk+['HET']]
        
        # # Pull references to Species declared in main rxn file out of the het rxn file (e.g. duplicate declarations)
        # [hsp_list.remove(sp) for sp in sp_list  if sp in hsp_list]
        
        # # Make sure rates are in MATLAB syntax... 
        # hks=[_str_multi_replace(k, {'.FALSE.':'false', '.TRUE.':'true'}) for k in hks]
        
        # # Make a mechanism title for the het rxns file: 
        # het_title, het_set_globals, het_call, het_jzero = make_header_for_F0AM_mech(het_rfuncts,[], kppfile)
        
        # het_fname=output_fname.split('.m')[0]+'_HetRxns.m'
        # fio.write_mech_to_file(het_fname, hsp_list, [''],hrxns, hks,hgs,hfs,mech_title=het_title, 
        #                lines_to_add_before_mech=het_set_globals+het_call+het_jzero) 
    
    # Make a header & add info to the gas/photolysis rxns file 
    mech_title, mech_comments = make_header_for_AllRxns(kppfile, GC_version, rates, photo_rates,jmap_type)

    # Now, FINALLY, Write  gas/photolysis  mechanism to a F0AM file:  
    fio.write_mech_to_file(output_fname,  species_list=sp_list,  ro2_list=[''],
                           rxn_list=rxns,       k_list=ks,       g_list=gs,  f_list=fs,
                           mech_name=output_fname.replace('.m',''), mech_title=mech_title, 
                           mech_comments=mech_comments, 
                           out_filepath=output_dir, overwrite=overwrite)
    
    # Unique list of all rate functions used in Het and regualr rxns... 
    if include_het_rxns is True: 
        all_rfuncts=np.unique(r_functs+het_functs).tolist()
    else: 
        all_rfuncts=np.unique(r_functs).tolist()
        
    return all_rfuncts, rxn_rate_dict, fjx_df

# %%###########################################################################
#-----------------(B) CREATE THE GEOSCHEM_K.m File ----------------------------
###############################################################################

def make_GEOSCHEM_K_file(kppfile:str, GC_version:str, template_file:str, all_rfuncts:list, rxn_rate_dict:dict,
                      output_fname:str='GEOSCHEM_K.m', output_dir:str='', overwrite:bool=False):
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
    
        kppfile       - STR with full path to the KPP input file that's being parsed. 
        
        GC_version    - STR formatted as '13.3.4' corresponding to what GEOS-Chem 
                        version you are seeking to create a F0AM file for. 
                        
        template_file - STRING containing the full path to the template file for GC_Rates_Template.txt
        
        all_rfuncts   - LIST of all rate functions used in the mechanism that mustl be defined! 
        
        rxn_rate_dict - DICT with keys corresponding to unique rate functions called in mechanism, 
                        that point to values of a nested DICT with key/value pairs of:  
                            'rxns' => LIST of rxns that use that specific rate call 
                            'rate_functions' => LIST of the rate function being called as a STR, no parenthesis or arguemnts to rate
                            'rate_args' => LIST of a LIST with the rate function arguements passed to that in 'rate_function'. 
                        that was generated as output from the function make_GEOSCHEM_AllRxns_file() and passed here. 
    
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
    ###########################################################################
    # Read in GC_Rates Template and customize the header used ... 
    ###########################################################################
    function_lines=list([]); # Initialize list of lines we'll write to the  output file...
    
    # Figure out how many lines there are in the "Template" file so you know when you reach the end. 
    num_lines = sum(1 for line in open(template_file))  
    
    # Using readlines(), open the template file and read in line by line. 
    inF= open(template_file, 'r')
    lines = inF.readlines()
    count = 0; # Initialize line counter variable. 
    
    for line in lines:
        # On most lines we don't need to modify anything, but there are several 
        # "Trigger" strings where we need to insert info about the mechanism...  
        
        # Update the template's mechanism name in the function header
        if '**MECHANISM_NAME**' in line.strip(): 
            mech_name='GEOSCHEM_v'+ GC_version.replace('.','_')
            line=line.replace('**MECHANISM_NAME**', mech_name)
            
        # Update the date the file was generated on in function header: 
        if '**THE_DATE**' in line.strip():
            line=line.replace('**THE_DATE**', str(date.today()))

        # Update the line to List the KPP file that was used to generate this function in MATLAB:
        if '**KPP_FILENAME**' in line.strip(): 
            line=line.replace('**KPP_FILENAME**', kppfile)
                              
        # Update the line that lists the total number of rate constants defined in this file: 
        if '**NUMBER_OF_RATE_CONSTANTS**' in line.strip():
            line=line.replace('**NUMBER_OF_RATE_CONSTANTS**',str(len(all_rfuncts)))
        
        # Add a list to the header of all the rate functions required to be imported! 
        if '**LIST_OF_RATE_FUNCTS**' in line.strip():
            # Create the list of rate_functions required: 
            all_rfuncts=[r+'()' for r in all_rfuncts if r!='None']
            req_rates='rate_functions={'+utils.join_list_for_MATLAB(',', all_rfuncts, 
                            insert='                    ', insert_skip=[1],
                            comment=True, comment_skip=[1])+'};'  
            
            # And update the line: 
            line=line.replace('**LIST_OF_RATE_FUNCTS**',req_rates)
            
        if '**IMPORT_RATES_CALL**' in line.strip(): 
            # Create the line that is used to import all the GEOS-Chem Rates Functions 
            get_rates='['+utils.join_list_for_MATLAB(',', all_rfuncts, insert= '   ', insert_skip=[1])+'] = out{:};'
            line=line.replace('**IMPORT_RATES_CALL**',get_rates)
        
        if '**INSERT_RATE_DEFINITIONS**' in line.strip(): 
            # Create the lines that define rates with names in this file like these: 
            # Normal F0AM Example |   What ours will look like: 
            # --------------------|-------------------------------------
            # i=i+1;              |    i=i+1;            
            # Knames{i} = 'KDEC'; |    Knames{i} = 'GCARR_ac(1300,5)';
            # krx(:,i) = 1e6;     |    krx(:,i) = GCARR_ac(Met,1300,5);
            
            iline="i=i+1;\n" # Define "iline" for the 1st thing that goes in a rate def! 

            rate_def_ls=list([]) 
            for k_i,kname_i in enumerate(list(rxn_rate_dict.keys())): 
                                
                # Create the line that defines the "KNames" 
                if len(rxn_rate_dict[kname_i]['rxns'])==1:# If Kname_i is used in only 1 reactions. 
                    kdec="Knames{i} = "+kname_i+"; % for: "+rxn_rate_dict[kname_i]['rxns'][0]+"\n"
                
                else: # If this rate is used in MULTIPLE reactions.. 
                
                    # Format the list of rxns that use this rate nicely to display in MATLAB 
                    rx_ls='\n'.join(['%    ' +r for r in rxn_rate_dict[kname_i]['rxns']])
                    
                    # Craft the declariation line & line stating info about what rxns use this rate: 
                    kdec0="Knames{i} = "+kname_i+";\n"
                    kdec1="% "+kname_i+" is used as the rate in these reactions: \n"+ rx_ls 
                    kdec=kdec0+kdec1+" \n"

                # Write actual call to function with args and add 'MET' as input!
                kdef= "krx(:,i) = "+rxn_rate_dict[kname_i]['rate_function'][0]+"(Met, "+','.join(rxn_rate_dict[kname_i]['rate_args'][0])+");\n"

                # Add this indv rate defintion to the list of rate_defs
                rate_def_ls.append(iline+kdec+kdef)
                
            # Join all rate_defs into a single string.... 
            rate_defs='\n'.join(rate_def_ls)
            
            # Replace the line with out insertion! 
            line=line.replace('**INSERT_RATE_DEFINITIONS**', rate_defs)
                    
        # After modifying lines (if needed), save the line to output list! 
        function_lines.append(line)
        
        count += 1 # Update counter. 
        
        if count>num_lines: break  # Exit while loop if you reach end of file. 
    
    inF.close() # Close the template file now that you got all the good stuff out. 
    
    ###########################################################################
    # Write the GEOSCHEM_vX_X_K.m rates file: 
    ###########################################################################     
    
    # Decide what to call the output file & make sure the path to it exists... 
    rate_file= utils.check_filename(filename=output_fname, default_name='GEOSCHEM_K.m', ext='.m',
                              savepath=output_dir, overwrite=overwrite, return_full=True)
    
    # Open the Output file to Write
    outF = open(rate_file, "w");     
    
    # Write each line to the output file 
    [outF.write(ln) for ln in function_lines] 

    outF.close() # Close the output file  
    print('F0AM Rate Defintion file saved at: '+rate_file) 
        
    return 


# %%###########################################################################
#-----------------(C) CREATE THE import_GC_Rates.m File -----------------------
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
    master_funct_lines=list([])
    
    # Figure out how many lines there are in the "Template" file so you know when you reach the end. 
    num_lines = sum(1 for line in open(template_file))  
    
    # Using readlines(), open the template file and read in line by line. 
    file = open(template_file, 'r')
    lines = file.readlines()
     
    count = 0; # Initialize line counter variable. 
    
    # Loop over each line... 
    for line in lines:
        # On most lines we don't need to modify anything, but there are several 
        # "Trigger" strings where we need to insert info about the specific 
        # rates that are being defined in this file where we need to make modifications. 
    
        # Update the template's Commented Out Example to contain an example for importing "all" the rates used in this mech.
        if 'CMT_INSERT_FUNCTION_NAMES' in line.strip(): 
            cmt_list_names=utils.join_list_for_MATLAB(',', all_rfuncts, insert='    %       ')
            cmt_list_names=cmt_list_names[12:] # remove first '    %       '
            line=line.replace('CMT_INSERT_FUNCTION_NAMES', cmt_list_names)

        # Update the actual function where rate_keys are defined to include each 
        # of the actual functions referenced in this file with their names. 
        if 'INSERT_INDV_FUNCTION_NAMES' in line.strip():
            rfuncts_as_strs=["'"+rate+"'" for rate in all_rfuncts]
            list_names=utils.join_list_for_MATLAB(',', rfuncts_as_strs, insert='                ')
            list_names=list_names[16:] # Remove first spaces... 
            line=line.replace('INSERT_INDV_FUNCTION_NAMES', list_names)
            
        # Update the master function where rate handles are defined as keys in dict! 
        if 'INSERT_INDV_FUNCTION_HANDLES' in line.strip(): 
            handles=['@'+funct for funct in all_rfuncts]
            list_handles=utils.join_list_for_MATLAB(',', handles,insert='                ')
            list_handles=list_handles[16:] # Remove first spaces... 
            line=line.replace('INSERT_INDV_FUNCTION_HANDLES', list_handles)            
            
        # After modifying lines (if needed) in master function, save each to output list! 
        master_funct_lines.append(line)
        
        count += 1 # Update counter. 
        
        if count>num_lines: break  # Exit while loop if you reach end of file. 
    
    # While loop will break once we found the end of the template file and have parsed all lines. 
    file.close() # Close the template file. 
    
    return master_funct_lines

# C.3 L1-function called inside make_import_GC_rates_file()     
def redo_function_header(line:str, spl_str:str='FUNCTION'): 
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
        l=line.split(spl_str)[1].strip() # e.g: ['', 'ARRPLUS_ade( a0, d0, e0 ) RESULT( k )']
        
        # Get the "result" output var ("k") isolated from the word "RESULT" and "(", and ")" 
        l=l.split('RESULT') #e.g: ['ARRPLUS_ade( a0, d0, e0 )' , '( k )']
        output_var= utils.str_multi_replace(l[1].strip(), ['(', ')'], rep_all='').strip() # e.g: 'k'
        
        # Now remove spaces from the function title & arguements: 
        funct_title_args=l[0].strip() # e.g:  ARRPLUS_ade(a0,d0,e0)'  
        
        # Add "Met" as an input arg so we can always retrieve the Met Vars in MATLAB as inputs to this function: 
        title_ls= funct_title_args.split('(') # e.g: ['ARRPLUS_ade', 'a0,d0,e0)']
        new_title_args=title_ls[0]+'(Met,'+title_ls[1] # e.g: 'ARRPLUS_ade(Met,a0,d0,e0)']
        
        # Add it all back together how MATLAB expects: 
        matlab_dec= 'function ['+output_var+'] = '+new_title_args.replace(' ', '') 
         
    else:  
        # Turn FORTRAN:  
        #   'REAL(kind=dp) FUNCTION FALL( A0, B0, C0, A1, B1, C1, CF )' 
        # Into MATLAB:  
        #   'function [FALL] = FALL(Met,A0,B0,C0,A1,B1,C1,CF)'
        
        # Split the input line at the word "FUNCTION" to just get the rest & remove spaces
        l=line.split(spl_str)[1].strip() # e.g: 'FALL(A0,B0,C0,A1,B1,C1,CF)'
        
        # Split the function on the "(" character: 
        title_ls=l.split('(') # e.g. ['FALL', 'A0,B0,C0,A1,B1,C1,CF)']
        
        # Isolate just the function name (to get the output var) & remove spaces: 
        output_var=title_ls[0].strip() # e.g: 'FALL'
        
        # Insert "Met" as an input arg so we can always retrieve the Met Vars in MATLAB as inputs to this function: 
        new_title_args=title_ls[0]+'( Met,'+title_ls[1]
        
        # Add it all back together how MATLAB expects: 
        matlab_dec= 'function ['+output_var+'] = '+new_title_args.replace(' ', '')  
 
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
    not_spaces=[i for i,ll in enumerate(line) if ll !=' '] 
    
    # If what you want to insert is not equal to the very first NON-Space character: 
    if line[not_spaces[0]] != insert:  
        # Then re-shuffle your line to be: (1) whatever came before the spaces, 
        # plus (2) what you want to insert, and then plus (3) whatever came after it...
        line= line[0:not_spaces[0]]+ insert + line[not_spaces[0]:] 
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
    not_spaces=[i for i,ll in enumerate(line) if ll !=' '] 
    if len(not_spaces) > 0:  
        tf=False if line[not_spaces[0]] not in ['%', '!'] else True 
    else:  
        tf=False 
     
    return tf 

# C.1 Function used in to create the import_GC_rates.m output file: 
def make_import_GC_rates_file(rate_files:list, all_rfuncts:list, template_file:str,
                               include_het_rxns: bool= False, output_fname:str='import_GC_rates', 
                               output_dir:str='', overwrite: bool=False):  
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
    mechnm=os.path.basename(rate_files[0]).split('.')[0] 

    # -------------------------------------------------------------------------------------
    # ------     Convert all the Fortran functions to matlab functions!  ------------ 
    # -------------------------------------------------------------------------------------
    key=''; all_lines= list(); fdict=dict({})
    for file in rate_files: 
        count = 0; pass_go=False  # Initialize counter 
        
        # Count # of lines in file the original file so we know when we reached the end.  
        num_lines = sum(1 for line in open(file, encoding="UTF-8", errors='ignore' )) 
        
        # Open the Input file to Read  
        inF = open(file, 'r',encoding="UTF-8", errors='ignore' ) 
        
        while True:  # Loop over each line of the KPP rate file. 
         
            line = inF.readline()  # Read next line from file 
            line=line.replace('\n', '') #Remove new line characters.  
            valid_starts=['RATE-LAW FUNCTIONS FOR GAS-PHASE REACTIONS','BEGIN RATE LAW FUNCTIONS','ARRHENIUS FUNCTIONS']
            valid_ends=['End INLINED Rate Law Functions','END MODULE '+mechnm+'_RateLawFuncs',
                        'RATE-LAW FUNCTIONS FOR HETEROGENEOUS REACTIONS']
    
            if any([item in line.upper() for item in valid_starts]): pass_go=True # Boolean of whether to start writing file or not.   
             
            if pass_go is True:  
                if len(line) > 0 :  
                    if '::' in line and 'PARAMETER' in line: # keep parameter definitions.  
                        line= insert_after_indent(line, line.split('::')[1].strip() +'&&baddie&&') 
                        line=line.split('&&baddie&&')[0] 
                    # Comment out all var definitions. Don't need them in MATLAB.  
                    if '::' in line and is_comment(line) is False: line= insert_after_indent(line, '%') 
                     
                    # Replace "END FUNCTION" with 'end' so it won't think the "end"line is a line with a function to fix...  
                    if 'END FUNCTION' in line and is_comment(line) is False: line= 'end'  
                     
                    # Replace "END SUBROUTINE" with 'end' so it won't think the "end"line is a line with a function to fix...  
                    if 'END SUBROUTINE' in line and is_comment(line) is False: line= 'end'  
                    
                    # Make sure we comment out "intrinsic variables.. b/c I have no idea what that means & neither does MATLAB 
                    if 'INTRINSIC' in line and is_comment(line) is False: line= insert_after_indent(line, '%') 
                     
                    # Move things around at function declarations to make MATLAB happy. 
                    if 'FUNCTION' in line and is_comment(line) is False : line= redo_function_header(line) 
                     
                    # Only after fixing functions can we make sure REAL declarations are commented out.  
                    if 'REAL' in line and is_comment(line) is False: line = '%'+line 
                    
                    # Now also comment out declariations of Character vars... 
                    if 'CHARACTER' in line and is_comment(line) is False: line = '%'+line 
                    
                    if 'WRITE' in line.upper() and is_comment(line) is False:
                        line="disp('"+ "'".join(line.split("'")[-2:])+")"
                     
                    # Check for use of state global vars in Fortran, which use a '%' but in middle of the line, not at beginning
                    if (is_comment(line) is False) and (' %' in line): line=line.replace(' %','')
                    if (is_comment(line) is False) and ('%' in line): line=line.replace('%','')
                    
                    # If its a one line "if" statement then need to add an "end" for MATLAB.  
                    if 'IF' in line and 'THEN' not in line and 'END IF' not in line and \
                        '&' not in line and 'ENDIF' not in line and 'DIFF' not in line: 
                        ln_ls=line.split(')')
                        test=ln_ls[0]
                        if len(ln_ls)==1: 
                            line =test+'); end '
                        elif len(ln_ls)>1:
                            contents=')'.join(ln_ls[1:]); sub=False
                            for cmt in ['!','%']:
                                if cmt in contents: 
                                    content=contents.split('!')
                                    line =test+');'+ content[0]+'; end %'+' '.join(content[1:])
                                    sub=True
                            if sub== False:  line =test+');'+contents+'; end '
                        else: 
                            print(line, len(ln_ls))
                    
                    # Move things around at function declarations to make MATLAB happy. 
                    if 'SUBROUTINE' in line and is_comment(line) is False : line= redo_function_header(line,spl_str='SUBROUTINE') 
                
                    # Otherwise, FORTRAN --> MATLAB is a series of operator substitutions...  
                    fortran_2_matlab= dict({ '!':'%',   
                                             '**':'.^',  
                                             '*' :'.*',  
                                             '/':'./', 
                                             'EXP':'exp',  
                                             'temp':'TEMP',  
                                             'MAX' :'max', 
                                             'LOG10':'log10',  
                                             '&':'...',  
                                             'ELSE IF':'elseif', 'ELSEIF':'elseif',
                                             'ELSE':'else', 
                                             'END MODULE':'% END MODULE', 
                                             'ENDIF':'end',  'ENDif':'end', 'END IF':'end',
                                             'IF':'if',  
                                             'DBLE':'double',   
                                             'e+':'e',  
                                             'E-':'e-', 
                                             'E+':'e',  
                                             'RETURN':'',  
                                             'SQRT':'sqrt',  
                                             'CALL': "", # Not a great solution just yet .. 
                                             'ENDDO':'end',  
                                             '.TRUE.': 'true', '.FALSE.': 'false',
                                             '.eq.':'==',  '.EQ.':'==',
                                             '.not.':'~',  '.NOT.':'~',
                                             '.lt.':'<',   '.LT.':'<', 
                                             '.gt.':'>',   '.GT.':'>',  
                                             '.ne.':'~=',  '.NE.':'~=', 
                                             '.ge.':'>=',  '.GE.':'>=', 
                                             '.le.':'<=',  '.LE.':'<=',  
                                             '.and.':'&&', '.AND.':'&&', 
                                             '.or.': '||' , '.OR.': '||'}) 
                     
                    line= utils.str_multi_replace(line,fortran_2_matlab) 
                     
                    # Remove things like double precision and 'Then' statements - unnecessary.  
                    line = utils.str_multi_replace(line,['_dp', 'THEN'], rep_all='') 
                     
                    # Fix anything my substitutitons messed up. 
                    line= utils.str_multi_replace(line,dict({'../':'./', '..*': '.*', '..^':'.^', 
                                                       ' .* ':'.*', ' ./ ':'./',
                                                       '); ) ;': '));', '); ;': ');'})) 
                    
                    # Make sure all the lines are suppressing their content in MATLAB by adding a ';' where appropriate.  
                    if 'end' not in line and 'function' not in line and is_comment(line) is False and line[-3:]!= '...':  
                        if '%' not in line:  
                            if len(line.replace(' ', ''))!=0 : line=line+';' # Make sure stuff is on line if place a ;
                        else:  
                            cmt_i=line.index('%') 
                            line=line[0:cmt_i]+';'+line[cmt_i:] 
                    
                    # Collect all lines from every function in a dictionary .... 
                    if not any([True if e in line else False for e in valid_ends ]): 
                        if line.strip() not in ['%', ';', '']: # Don't store empties! 
                            if 'function [' in line: 
                                key=line.split( '=')[1].split('(')[0].replace(' ','')
                                fdict[key]=[line]
                            elif line!='end' and key != '':
                                if line.strip()[0]=='%': line= '    % '+line.strip()[1:].strip()
                                all_lines.append(line)
                            elif line =='end' and key!= '': 
                                fdict[key]=fdict[key]+all_lines+[line]
                                key=''; all_lines= list()

            count=count+1 # Update counter var.  
             
            # Break if reached end of file or end of the functions for rate constants.  
            if count>num_lines or any([True if e in line else False for e in valid_ends ]): break 
        inF.close() # Close the input file  
    
    # Make sure the output file path exists... Set default and and ext incase users don't!  
    rate_file= utils.check_filename(filename=output_fname, default_name='import_GC_rates', ext='.m', 
                              savepath=output_dir, overwrite=overwrite, return_full=True)
    
    outF = open(rate_file, "w"); # Open the Output file to Write  
    
    # -------------------------------------------------------------------------------------
    # ------Write the Wrapper Function for All rates to the Output Rates File! ------------ 
    # -------------------------------------------------------------------------------------

    # Create the wrapper function that needs to be written to the "rates" file
    wrapper_funct_lines= make_GCrates_wrapper_funct(template_file, all_rfuncts)
    
    # Write this "Wrapper" function to the top of the output rates file. 
    for ln in wrapper_funct_lines: 
        if '\n' in ln: 
            outF.write(ln)
        else: 
            outF.write(ln+'\n')
    
    # Write a line to retrieve all global vars inside each function. 
    get_globals='    % Pass input struct, "Met" to function "extract_MetVars()" to get vars this function may depend on...\n'+ \
                '    [TEMP, PRESS, NUMDEN, H2O, TEMP_OVER_K300, K300_OVER_TEMP,INV_TEMP,SR_TEMP,RELHUM]= extract_MetVars(Met);\n' 
    for key in fdict: # Loop over all functions 
        glbs_tf=False
        ln_list=fdict[key] # list of lines we need to write for each function. 
        for i,line in enumerate(ln_list): 
            outF.write(line+'\n')
            if glbs_tf is False and ')' in line:
                outF.write(get_globals+'\n'); glbs_tf=True

    outF.close() # Close the output file  
    print('MATLAB compliant version of gckpp_Rates.F90 called in "GEOSCHEM_K.m" to define all rate functions saved at: '+rate_file ) 
    
    return


# %%###########################################################################
#-----------------(D) Wrapper Function to Make ALL Files -----------------------
###############################################################################
# D.1 L1-Function used in make_GC_mechanism() to check user inputs to make_GC_mechanism()
def _check_user_inputs(kppfile, rate_files, GC_version:str, jmap_type:str, output_dir:str,
                       template_paths:dict, photolysis_paths:dict): 
    
    # Check that the output directory the user entered exists: 
    if os.path.exists(output_dir) is False: 
        raise NotADirectoryError("INPUT ERROR in make_GC_mechanism(): \n"+
                         "The following output directory you listed to store the output files does not exist: \n"+ 
                         "    output_dir='"+output_dir+ "'").with_traceback(sys.exc_info()[2])
        
    # Check that the KPP files passed as input exist
    if os.path.isfile(kppfile) is False:    
        raise FileNotFoundError("INPUT ERROR in make_GC_mechanism(): \n"+
                         "The path to the KPP file you entered does not exist: \n"+ 
                         "    kppfile='"+kppfile+ "'").with_traceback(sys.exc_info()[2])
        
    # Check that the gckpp rates files passed as input are the right type & exist: 
    if type(rate_files) not in [type('str'), type(list([1]))]: 
        raise TypeError("INPUT ERROR in make_GC_mechanism(): \n"+
                        "Input for argument 'rate_files' must either be TYPE = STR or LIST: \n"+
                        "    (1) A STR with a path to a gckpp_Rates.F90 file \n " + 
                        "    (2) A LIST of strings with paths to rate function files like: \n"+
                        "             ['/path/to/gckpp_Rates.F90', /path/to/gckpp_HetRates.F90'] \n" + 
                        "Instead, you entered input of 'rate_files' of type="+str(type(rate_files))).with_traceback(sys.exc_info()[2])  
    # If its a single rate file ... 
    elif type(rate_files) == type('str') and os.path.isfile(rate_files) is False: 
        raise FileNotFoundError("INPUT ERROR in make_GC_mechanism(): \n"+
                                 "The path to the GCKPP rate file you entered does not exist: \n"+ 
                                 "    rate_files='"+rate_files+"'").with_traceback(sys.exc_info()[2])   
    # If its a list of rate files.. 
    elif type(rate_files) == type(list([1])) and not all([os.path.isfile(file) for file in rate_files]): 
        err_files=["    ("+str(i+1)+") '"+file+"'" for i,file in enumerate(rate_files) if os.file.exists(file) is False]
        raise FileNotFoundError("INPUT ERROR in make_GC_mechanism(): \n"+
                                 "The path to the following gckPPP rate file(s) you entered does not exist: \n"+ 
                                 '\n'.join(err_files)).with_traceback(sys.exc_info()[2])    

    # Check that the template files exist at these paths.
    for key in list(template_paths.keys()): 
        if os.path.isfile(template_paths[key]) ==False: 
            raise FileNotFoundError("INPUT ERROR in make_GC_mechanism(): \n"+
                                     "The path to the "+key+" template file was not found where expected at: \n"+ 
                                     "    '"+template_paths[key]+"'\n"+ 
                                     "It should be located at: '.../GEOSChem_Emulator/Input_Files/Template_Files/"+key+"' + \n"+
                                     "and if you get this error, there may be an install issue!").with_traceback(sys.exc_info()[2]) 
    
    # Check that the Photolysis files exist at these paths.
    for key in list(photolysis_paths.keys()): 
        if os.path.exists(photolysis_paths[key]) ==False: 
            raise FileNotFoundError("INPUT ERROR in make_GC_mechanism(): \n"+
                                     "The path to the "+key+" photolysis input file was not found where expected at: \n"+ 
                                     "    '"+photolysis_paths[key]+"'\n"+ 
                                     "It should be located at: '.../GEOSChem_Emulator/Input_Files/Photolysis_Files/"+key+"' \n"+
                                     "and if you get this error, there may be an install issue!").with_traceback(sys.exc_info()[2]) 
    
    # Check that the jmap_type is allowed: 
    if jmap_type not in ['CrossSec_Match','MCM','RCIM']: 
        # Be case insensitive & check for miss-spelling
        if jmap_type.lower() in ['crosssec_match','crossec_match','mcm','rcim']: 
            if jmap_type.lower() in ['crosssec_match','crossec_match']:jmap_type='CrossSec_Match'
            if jmap_type.lower()=='mcm': jmap_type='MCM'
            if jmap_type.lower()=='rcim': jmap_type='RCIM'
        else: 
            raise ValueError("INPUT ERROR in make_GC_mechanism(): \n"+
                            "Input for argument 'jmap_type' must be one of the following: \n" +
                            "    (1) 'CrossSec_Match' \n" + 
                            "    (2) 'MCM' \n" + 
                            "    (3) 'RCIM' \n" +             
                            " For more info on these options try: >>utils.jmap_type_help()").with_traceback(sys.exc_info()[2]) 
            
    return jmap_type

def make_GC_mechanism(kppfile, rate_files, GC_version:str, jmap_type:str, include_het_rxns:bool, 
                      output_dir:str, overwrite:bool=False): 
    
    # Get the relative paths to input/templates, etc... 
    GCEM_path, input_path, template_paths, photolysis_paths, jmap_type_help = utils.get_my_relative_paths()
    
   # Check the user's inputs to make sure all files we need exist and are located where we expect them. 
    jmap_type= _check_user_inputs(kppfile=kppfile, rate_files=rate_files, GC_version=GC_version, 
                      jmap_type=jmap_type, output_dir=output_dir, template_paths=template_paths, 
                      photolysis_paths=photolysis_paths)
    
    # Write the GEOSChem_AllRxns.m file & get output we need to write GEOSChem_K.m file:
    all_rfuncts,rxn_rate_dict, fjx_df=make_GEOSCHEM_AllRxns_file(kppfile = kppfile, 
                                                         GC_version = GC_version, 
                                                         jmap_type = jmap_type, 
                                                         include_het_rxns = include_het_rxns,   
                                                         output_fname = 'GEOSCHEM_AllRxns.m',
                                                         output_dir=output_dir, 
                                                         overwrite=overwrite)
    
    # Write the GEOSChem_K.m file 
    make_GEOSCHEM_K_file(kppfile=kppfile, 
                         GC_version = GC_version, 
                         template_file = template_paths['GEOSChem_K_template.txt'],
                         all_rfuncts = all_rfuncts, 
                         rxn_rate_dict = rxn_rate_dict,
                         output_fname='GEOSCHEM_K.m', 
                         output_dir=output_dir, 
                         overwrite=overwrite)

    # Create import_GC_rates.m file where all rate functions are defined which is 
    # referenced and used as is without needing modifications in the GEOSChem_K.m file. 
    if type(rate_files)==str: rate_files=[rate_files]
    
    fdict=make_import_GC_rates_file(rate_files=rate_files, 
                                    all_rfuncts=all_rfuncts, 
                                    template_file= template_paths['import_GC_rates_template.txt'],
                                    include_het_rxns=include_het_rxns, 
                                    output_fname='import_GC_rates.m',
                                    output_dir=output_dir, 
                                    overwrite=overwrite)
    
    return fjx_df




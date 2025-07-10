#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 23 18:29:55 2025

@author: u6044586
"""
import re 
import numpy as np 
from itertools import chain
import utils

# -------------- Prep/Write to F0AM  & Sub-Functions---------------------

# Helper function called within prep4f0am()
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

# Called within prep4f0am() to format the f-line of the F0AM input: 
def format_f_line(data, fixed_vars): 
    """Helper function to prepare/format the F-line for F0AM""" 
    # Initialize a dictionary with keys as the (unique) species involved in 
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
                                 
    rct_dict = {};prd_dict = {};fstrs=[]
    
    def accumulate_stoich(cmpd:str, stoich, net_dict:dict):
        # Local Helper function to accumulate stoichiometry for input compound 
        # in dictionary or assign it to 0 if no value exists. 
        if cmpd not in net_dict:
            net_dict[cmpd] = 0.0
        net_dict[cmpd] += stoich
        return net_dict
    
    # Accumulate reactant stoichiometry  
    for cmpd, stoich in zip(data['rct_cmpds'], data['rct_stoich']):
        #Rct stoich negative, as rxn consumes this compound
        rct_dict=accumulate_stoich(cmpd, -stoich, rct_dict)
        
    # Accumulate product stoich 
    for cmpd, stoich in  zip(data['prd_cmpds'], data['prd_stoich']):
        #Prd stoich positive, as rxn creates this compound
        prd_dict=accumulate_stoich(cmpd, stoich, prd_dict)
        
    for cmpd, stoich in chain(rct_dict.items(), prd_dict.items()):
        if stoich!=0: 
            operator='+' if stoich > 0 else '-'
            if cmpd not in fixed_vars: # Don't adjust concs of fixed vars! 
                fstrs.append(f"f{cmpd}(i)=f{cmpd}(i){operator}{np.abs(stoich)};")       
                
    data['F_Line'] =' '.join(fstrs)
    
    return data 


# Called within prep4f0am() to add "Met" as arg to rate functions assigned in GEOSChem_K.m. 
def add_met_to_rates(string): 
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


# Top level function called within make_GC_mechanism() to convert the raw output
# dict form parsing a input files to get it ready/formatted  for F0AM: 
def prep4f0am(kpp_dict, inds2drop, fixed_vars, jmap_dict): 
    
    
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
    maybe_add_to_het_js=[];
    
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
                # Create a name for the rate using the rxn: 
                rate_name='_'.join(data['rct_cmpds'])+'_to_'+'_'.join(data['prd_cmpds'])
                
                # Make sure its not so long it makes MATLAB angry... 
                if len(rate_name) > 60: 
                    # Split the rate name into its components
                    parts = rate_name.split('_')
                
                    # Remove elements from the end until the length constraint is satisfied
                    while len(rate_name) > 60 and len(parts) > 3:
                        parts.pop()  # Remove the last element
                        rate_name = '_'.join(parts)  # Recreate the string

                # Since its a func, we'll make it a named rate in the GEOSChem_AllRxns.m 
                #(e.g. written as a string) that is later defined in the 
                # GEOSChem_K.m & import_rates_GC.m files... 
                data['K_Line']=f"k(:,i) = {rate_name};" 
                
                if any([fvar in data['rct_cmpds'] for  fvar in fixed_vars]): 
                    if 'N2' in data['rct_cmpds']: 
                        rate=f'({rate}).*0.78.*Met.M'
                    if 'O2' in data['rct_cmpds']: 
                        rate=f'({rate}).*0.21.*Met.M'  
                    if 'H2' in data['rct_cmpds']: 
                        rate=f'({rate}).*5.5e-7.*Met.M'    
                    
                # If it IS a function, then add 'Met' as an argument for all rate 
                # constant functions so we can indeed use Temp/ M/RH as vars in
                # the rate constant functions & add it to key/value pair: 
                if data['is_het']==False: 
                    if "'"+rate+"'" not in list(foam_dict['nhet']['rate_dict'].keys()): 
                        foam_dict['nhet']['rate_dict']["'"+rate_name+"'"]=add_met_to_rates(rate)
                else: 
                    if "'"+rate+"'" not in list(foam_dict['het']['rate_dict'].keys()): 
                        foam_dict['het']['rate_dict']["'"+rate_name+"'"]=add_met_to_rates(rate)
                    

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
            if rate not in list(jmap_dict['Inds_to_JNames'].keys()):
                raise ValueError("In function, prep4f0am(), we could not find a J-Value defined in F0AM \n" +
                                 "within 'jmap_dict' corresponding to: '"+rate +
                                 "'\n\nThis likely means you need to update 'jmap_dict' to account for new F0AM/FJX mapping j-value changes!")
            else:  
                # Otherwise, use the mapped value in J-Dict to define rate in F0AM. 
               
                # Pull out comment on where the value used in PHOTOL(IND) comes from:
                # j-nickname of PHOTOL(IND) is: 
                j_name= jmap_dict['Inds_to_JNames'][rate].replace("'", '')
                j_info= jmap_dict['JNames_to_Vals'][j_name]['Info']
                
                # Add the info on the rxn / prepre the rate for F0AM as we did abmove. 
                data['K_Line']=f"k(:,i)={j_name}; {j_info}"
                
                # Likewise, keep up with the Named Js too. 
                data['Function']=jmap_dict['Inds_to_JNames'][rate]
        
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
        data=format_f_line(data,fixed_vars)
        
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
                
                # Add required js from this rxn we're commenting out to list to check @ end... 
                # to see if we need to add it to req'ks' or req_js' for het rxns... 
                if data['is_photo'] ==True: 
                    maybe_add_to_het_js.append(data['Function']) 
                
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
        
    # Move any non-het J's that are commented out b/c are photolysis of a het 
    # species to foam['het']['req_js'] if they aren't otherwise used in the 
    # nhet mech & aren't in the het mech... 
    for j in maybe_add_to_het_js: 
        if (j not in foam_dict['nhet']['req_js']) and (j not in foam_dict['het']['req_js']):
            foam_dict=_append_to_ndict_list(foam_dict, 'het', 'req_js', j, allow_dupes= False)
       
    return foam_dict


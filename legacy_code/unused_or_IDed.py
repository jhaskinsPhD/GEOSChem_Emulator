#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 13:05:24 2023

@author: u6044586
"""

if '**INSERT_RATE_DEFINITIONS**' in line.strip():
    # Create the lines that define rates with names in this file like these:
    # Normal F0AM Example |   What ours will look like:
    # --------------------|-------------------------------------
    # i=i+1;              |    i=i+1;
    # Knames{i} = 'KDEC'; |    Knames{i} = 'GCARR_ac(1300,5)';
    # krx(:,i) = 1e6;     |    krx(:,i) = GCARR_ac(Met,1300,5);

    iline = "i=i+1;\n"  # Define "iline" for the 1st thing that goes in a rate def!

    rate_def_ls = list([])
    for k_i, kname_i in enumerate(list(rxn_rate_dict.keys())):

        # Create the line that defines the "KNames"
        # If Kname_i is used in only 1 reactions.
        if len(rxn_rate_dict[kname_i]['rxns']) == 1:
            kdec = "Knames{i} = "+kname_i+"; % for: " + \
                rxn_rate_dict[kname_i]['rxns'][0]+"\n"

        else:  # If this rate is used in MULTIPLE reactions..

            # Format the list of rxns that use this rate nicely to display in MATLAB
            rx_ls = '\n'.join(
                ['%    ' + r for r in rxn_rate_dict[kname_i]['rxns']])

            # Craft the declariation line & line stating info about what rxns use this rate:
            kdec0 = "Knames{i} = "+kname_i+";\n"
            kdec1 = "% "+kname_i+" is used as the rate in these reactions: \n" + rx_ls
            kdec = kdec0+kdec1+" \n"

        # Write actual call to function with args and add 'MET' as input!
        kdef = "krx(:,i) = "+rxn_rate_dict[kname_i]['rate_function'][0] + \
            "(Met, " + \
            ','.join(rxn_rate_dict[kname_i]['rate_args'][0])+");\n"

        # Add this indv rate defintion to the list of rate_defs
        rate_def_ls.append(iline+kdec+kdef)

    # Join all rate_defs into a single string....
    rate_defs = '\n'.join(rate_def_ls)

    # Replace the line with out insertion!
    line = line.replace('**INSERT_RATE_DEFINITIONS**', rate_defs)

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

        Output Usage: Used in output of prase_kpp_main() for kpp_dict['rates']

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
    fortran2MATLAB = dict({
        'E+00': '',  # Numbers raised to a 0 should be just #s in MATLAB.
        'd+00':'',
        'e00':'',
        'E0': '',
        'd0': '',
        '_dp': '',
        'E+': 'e',  # for numbers raise to a postivie power
        'd+': 'e',
        'E-': 'e-',  # for numbers raised to a negative power.
        'd-': 'e-',    
        'EXP(': 'exp('
    })

    # Replace all keys in dict that appear in line with their values & strip spaces!
    line = utils.str_multi_replace(line, fortran2MATLAB).strip()

    return line    
            
# A.3.2 L2-function called within KPP_Dict_to_F0AM_IO within make_GEOSCHEM_AllRxns_file()
def get_all_MCM_Js_needed(J_in: str, fjx_df: pd.DataFrame, all_Js: list = []):
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
    fjx_ind = fjx_df.index[fjx_df['PHOTOL(IND)'] == J_in]
    J_value = fjx_df.loc[fjx_ind[0], 'F0AM_Assignment']

    # Define the Regex string to find J-value matches. 'Jregex' will hit str matches to the letter "J"
    # followed by one or more numbers 0-9. Match will end when something other than 0-9 is encountered.
    # Won't match just "J" or return blank items in a list with spaces!
    J_pattern1 = r'J[0-9]+'
    J_pattern2 = r'Jzero'  # Catches "Jzero"...
    J_pattern3 = r'Jn[0-9]+'  # Catches "Jn3" etc.
    J_pattern4 = r'Jv[0-9]+'  # Catches "Jv3" etc.

    # Use regex's "Findall()" to return a list of all matches to MCM type J-values in the string "J_in"
    # That follow the format defined in 'J_pattern'
    indv_Jlist = re.findall(J_pattern1, J_value)+re.findall(J_pattern2, J_value) + \
        re.findall(J_pattern3, J_value)+re.findall(J_pattern4, J_value)

    # Pop an error if no J-Values were found in "Jlist".
    if len(indv_Jlist) == 0:
        raise ValueError("The function get_all_MCM_Js_needed() could not find any J-values in " +
                         "this string: '"+J_value+"' which is the assignment for '"+J_in+"'. \n\n" +
                         "This is either your fault (bad input) or means that the regex format string \n" +
                         "used in the function 'get_all_MCM_Js_needed()' needs to be modified to catch this instance!").with_traceback(sys.exc_info()[2])
    else:
        # Remove any duplicates from list!
        indv_Jlist = np.unique(indv_Jlist).tolist()

    # Add all (new) unique J-values used in 'J_in' to the master list of all Js used in new MCM mech!
    all_Js = np.unique(all_Js+indv_Jlist).tolist()

    return all_Js

def sort_my_js(js_used): 
    """Function to parse a list of Js used in mechanism and "sort" them numerically. 
    Puts MCM default J-values first and then any "custom" Jns in them added from other mechs. """
    
    # Get the # of the Jns out first as numbers and sort them. 
    jns=sorted([np.int64(item.replace('Jn','')) for item in js_used if 'Jn' in item and item != 'Jzero'])
    # Now get the # of the J values and sort them as well.
    rjs=sorted([np.int64(item.replace('J',''))  for item in js_used if 'Jn' not in item and item != 'Jzero'])
    
    # Decide if you need to add "Jzero" back in or not. 
    jzero=['Jzero'] if 'Jzero' in js_used else []
    
    # Now combine them back into a nice sorted list. 
    my_js=['J'+str(rj) for rj in rjs]+jzero+['Jn'+str(jn) for jn in jns] 
    
    # Check that you return a list the same length as you started with. 
    if len(my_js) != len(js_used):
        print('In KPP_to_F0AM/sort_my_js()... Not writing full list of rates to file header.'); 
        sys.exit()
    
    return my_js

# A.2.4 L2-function used within parse_kpp_main() within make_GEOSCHEM_AllRxns_file()
def _enforce_list_dict_inds(list_dict: dict, inds: list,  add_index: bool = False):
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
    out_dict = dict({})

    # Loop over all keys in input "list_dict".
    for key in list(list_dict.keys()):
        # Pull out the full "master_list" contained under this key.
        master_list = list_dict[key]

        # Create a NEW List to hold only items with index in the input "inds".
        new_list = []

        # Fill the new list from the master list where the index matches "inds"
        new_list = [master_list[i] for i in inds]

        # Save this new list in the output dictionary!
        out_dict[key] = new_list

        # Clear local variables so you don't mess up the next key/list pair.
        del new_list, master_list

    # Check that all the keys are still there and that all lists are the same length!
    _check_rinfo(out_dict, check_lens=True, ignore_extras=False)

    # Now store the index used to sub-select these values from list_dict if asked!
    if add_index is True:
        out_dict['origin_index'] = inds

    return out_dict

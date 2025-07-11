#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 23 14:30:02 2025

@author: u6044586
"""
import re 
import os 
import sys 
import numpy as np
import pandas as pd 

import utils_class as utils


def jmap_type_help(): 
    """Helper function to return to a confused uer, information about options 
    for the GEOS-CHem EMulator input argument, 'jmap_type'. 
    
    AUTHOR: 
    -------
        Prof. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD
         
    """
    # Get the relative path to the photolysis help file on your computer. 
    my_paths  = utils.get_my_relative_paths()
    
    #Open the help/readme document & print each line to the screen. 
    file = open(my_paths['FJX']['jmap_type_help'], 'r')
    lines = file.readlines()
    for line in lines:
        print(line)
        
    file.close()
    return 

# Helper function called within create_jmap_dict() to find right FJX file to read in: 
def find_my_FJX_file(photolysis_paths:str, GC_version:str, verbose:bool=False): 
    """
    Function that takes a GEOS-Chem version number, and provides the path to the 
    apprpriate FJX file to read in: 
    """ 
    # Read in dictionary defining what version of FJX files are used for each dif 
    # GEOS-Chem version (In the output dict, vfjx, the GEOS-Chem verison is the key
    vfjx = pd.read_excel(photolysis_paths['FJX_input_by_GC-version.xlsx'], 
                         index_col=0).to_dict()['FAST-JX-INPUT']

    # If you have info on what version of FJX files to use for this version of 
    # GEOS_Chem, in the dict, vfjx then use it: 
    if GC_version in list(vfjx.keys()):
       
        #######################################################################################
        # Define path to where the j2J file that should be used for photolysis mapping lives: 
        j2j_file = os.path.join(photolysis_paths['FJX_Files'], 'FJX_j2j_'+vfjx[GC_version]+'.dat')
        #######################################################################################
        
        # Check that that file actually exists: 
        exists = os.path.exists(j2j_file)
        
        # If the file does exist , tell them & add it to the dict, user_inp for later: 
        if exists is True:
            if verbose is True: 
                print(f'Using: "FJX_j2j_{vfjx[GC_version]}.dat" for GEOS-Chem version {GC_version} photolysis mapping.')
            return j2j_file
        else:
            raise FileNotFoundError('We could not locate the correct FJX file to open for GEOS-Chem Version: '+GC_version +
                                    ' at the following expected path:  \n' +
                                    '    '+j2j_file).with_traceback(sys.exc_info()[2])
            return
    else:
        raise ValueError('Version ', GC_version, ' of GEOS-Chem Classic was not found in the dictionary \n' +
                         'defining its corresponding FJX chem input file... \n\n ' +
                         'This likely means you are trying to Emulate a new version of GEOS-Chem not yet supported by our Emulator. \n'
                         'To proceed, you will need to figure out what FJX version is used for this version of GEOS-Chem and then: \n' +
                         '    (1) Open and edit the file: ".../GEOSChem_Emulator/Input_Files/Photolysis_Files/FJX_input_by_GC-version.xlsx" \n' +
                         "           - Add the GEOS-Chem version you're emulating in the FIRST column \n" +
                         '           - Add the version of the FJX file that GEOS-Chem version uses to the the SECOND column.' +
                         '    (2) Add the Raw FJX.dat file to the directory: "GEOSChem_Emulator/Input_Files/Photolysis_Files/FJX_Files/" if it does not exist already. \n\n' +
                         'Why? #1 Will make sure we know what FJX file to open/read when emulating this version of GEOS-Chem in our dictionary that popped this error and \n' +
                         '#2 Will make sure we have the actual FJX file to open and read.').with_traceback(sys.exc_info()[2])
        return
       
# Helper function called within create_jmap_dict() to read in FJX file ID'd in  find_my_FJX_file(). 
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
    
    # -------------------------------------------------------------------------
    # Join parts of the columns together to make a coherant photolysis reaction
    # -------------------------------------------------------------------------
    # Create an empty list to store all coherant photolysis reactions built from pars of file: 
    photo_rxns = []
    
    # Iterate over each index in the DataFrame
    for ind in range(len(fjx_df)):
        
        # Get the name of the compound being photolyzed at the current index and convert to string
        photo_of = str(fjx_df.loc[ind, 'Photolysis_Of'])
        
        # Split 'Into' column at the current index into products, filter, and join them with '+'
        all_prods = fjx_df.loc[ind, 'Into'].split(' ')
        prods = '+'.join([prd_i.replace('...', '') for prd_i in all_prods if len(prd_i) > 0])
        
        # Construct the reaction and append to the list
        rxn_i = photo_of + '+hv -> ' + prods
        photo_rxns.append(rxn_i)

    # Assign the constructed list to the new DataFrame column
    fjx_df['Photolysis_Rxn'] = photo_rxns
    
    # Clean up cross section used column from spaces and other characters we don't need!
    fjx_df['FJX_CrossSec'] = fjx_df['FJX_CrossSec'].str.replace(r'[\s/()\-\']', '', regex=True)
   
    # Stick a j in front of the cross section used to define a "nickname"...
    fjx_df['J-NickName'] = 'j' + fjx_df['FJX_CrossSec']
   
    # Now drop columns that we don't actually need.
    fjx_df = fjx_df.drop(columns=['Photolysis_Of', 'PHOTON', 'Into'])
    
    return fjx_df    

# Helper function called within create_jmap_dict() to ensure J-Nicknames created for us in F0AM Are unique. 
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

    """

    ct = 1  # Initialize counter variable

    while True:
        # If you already have something with this name: 
        if jname in current_names:
            
            # Pull the "name" of the thing being photolyzed out. For "jCH3COOH", nm_only='CH3COOH'
            nm_only = re.findall(r'j([a-zA-Z0-9]+)', jname)
            
            # De-bugging check to make sure we can ID what's being photolyzed... 
            if len(nm_only) == 0:
                print("Error in make_unique_Jname(). Couldn't ID what's being photolyzed. Exiting...")
                print(f"\n\t{jname}")
                sys.exit()
                
            # Decide on a new suffix to try using (iterates through alphabet).
            new_letter = '_'+chr(ord('`')+ct)
            jname = 'j'+nm_only[0]+new_letter
        
        #Once you find a new unique jname, go ahead and quite! 
        if jname not in current_names:
            break
        else:
            ct = ct+1  # Update counter var to try a new '_'+LETTER(ct) combo... 

    return jname

# Top level function called to crate the dictionaries/data fames used to do J-Value mapping: 
def create_jmap_dict(user_inp, tracer_info):
    """Function that takes a GEOS-Chem version number, and reads in the appropriate 
    FJX file. Then it uses the PHOTOL(IND) and cross-section name used to 
    match the GEOS-Chem J-values to an analogous J-value that is already defined in F0AM.

    NOTE: This means that the photolysis values used between GEOS-Chem Classic and 
    in this GEOS-Chem Emulator are not the same. We are simply mapping the j-values 
    based on how they are set in GEOS-Chem classic to their closest F0AM alternative. 

    Inputs: 
    -------
    (1) user_inp -  DICT containing user intputs. At minimum, must have the following keys: 
        
       'Photolysis_Files' - DICT with paths to input photolysis files (should be 
                            auto-assigned correctly if GCEM installed right). 
        
       'GC_version' -  STR formatted as '13.3.4' corresponding to what GEOS-Chem 
                       version you are seeking J-Value matches to/ trying to Emulate

       'jmap_type' - STR corresponding to what mapping data in the file, 
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

        (1) jdict - A nested dictionary with 2 top level keys: 
            
            'Inds_to_JNames'- DICT that has GEOS-Chem's PHOTOL(IND) as STR keys, and their 
                              matched, corresponding J- Nickname used in F0AM as its value. 
                              Contents used to assign named J-alues in the GEOSChem_Rxns.m file. 
                              (e.g.  'PHOTOL(22)':'jNO2') 
            'JNames_to_Vals'-   DICT that has a STR J-Nicknames used in the GEOSChem_Rxns.m file as 
                               the key and a STR with the assignment that should be used in GEOSChem_J.m 
                               to define that in terms of some J value that F0AM has predefined. 
                               
        (2) fjx_df - Adataframe with a ton of info on the J-values/mapping performed. 

    REQUIREMENTS: 
    -------------

        CUSTOM FUNCTIONS:    
            hv.find_my_FJX_file()
            hv.read_FJX_file()
            hv.make_unique_Jname()

    AUTHOR: 
    -------
        Dr. Jessica D. Haskins (jessica.haskins@utah.edu)   GitHub: @jhaskinsPhD 

    """
    
    # Pull out dict from user input dict  defining where the photolysis files live: 
    photolysis_paths=user_inp['input_paths']['Photolysis_Files']
    
    # And what type of J-Mapping was requested: 
    jmap_type=user_inp['jmap_type']
    
    ###########################################################################
    # (1) Read in the FJX file for this GC version as a pandas dataframe 
    ###########################################################################
    # Figure out what FJX file to use with this GEOS-Chem version/
    j2j_file= find_my_FJX_file(photolysis_paths=photolysis_paths,
                               GC_version=user_inp['GC_version'], 
                               verbose= user_inp['verbose']) 
    
    # Also, keep name of FJX file used to generate stuff in user_inp dict:
    user_inp['FJX_version']=j2j_file
    
    # Read that file in as a df. Tells you what CROSS-SECTION  is used as the
    # basis for each reaction in the mechanism. fjx columns include: 
    #      'FJX_Index'      -  INTEGER index used in FJX & KPP for index in PHOTOL()
    #      'Quantum_Yield'  -  FLOAT of the quantum yield of the reaction... 
    #      'FJX_CrossSec'   -  STR of the cross-section (from FJX_spec.dat) which will be used 
    #                          to calculate the first-order reaction rate. 
    #      'Photolysis_Rxn' -  STR of the actual photolysis reaction this is for... 
    #      'J_NickName'     -  Nickname of J-value... 
    fjx_df = read_FJX_file(j2j_file)
    
    ###########################################################################
    # Read in DataFrame with Photolysis Mapping: (Cross-Sections to F0AM J-Values) 
    ###########################################################################
    # Tells you when this cross section is used, what pre-defined F0AM J-value should 
    # be used? ( decided in parrt by jmap_type)This file must be edited when new 
    # cross sections are added to FJX... 
    cs = pd.read_excel(photolysis_paths['FJX_cross_sect_to_jvalues.xlsx']).fillna('')

    ###########################################################################
    # Assign J-Values in mechanism based on the Cross-Section used in FAST-JX:
    ########################################################################### 
    # Create an empty list-like column to hold J-values that will be used. 
    fjx_df['J_for_CrossSec'] = list([''])*len(fjx_df) 
    
    # Loop over all photolysis rates used in the mechanism: 
    for ind in fjx_df.index:  
        # Figure out what cross-section is used here / we need to assign a J-Value for...
        csect_need = fjx_df.loc[ind, 'FJX_CrossSec']
        
        #----------------------------------------------------------------------
        # -- Quick Error Check we found what cross-sect is used for this rxn --
        #----------------------------------------------------------------------
        if csect_need == '':
            # Figure out what reaction this is for...
            rx_i = fjx_df.loc[ind, 'Photolysis_Rxn'].split('->')[0].replace(' ', '')

            # Some versions ofFJX files don't have names for NITP cross sections... 
            # Fix that by looking for strings that'd use the NITP cross section: 
            if any([nit_str in rx_i for nit_str in ['NIT+hv', 'NITs+hv']]):
                csect_need = 'NITP'
                fjx_df.at[ind, 'FJX_CrossSec'] = 'NITP'
                fjx_df.at[ind, 'J-NickName'] = 'j'+fjx_df.loc[ind,'Photolysis_Rxn'].split('+')[0].replace(' ', '')
        #----------------------------------------------------------------------
        
        # Get Index in dataframe, cs, where that cross-section is defined:
        has = cs.index[cs['FJX_Cross_Section'] == csect_need]
        
        # If we have a match in the cross-section dataframe: 
        if len(np.unique(cs.loc[has, jmap_type])) == 1:  
            # Then use the J- value assigned to that cross-section for this rxn
            fjx_df.at[ind, 'J_for_CrossSec'] = np.unique(cs.loc[has, jmap_type])[0]

        # Multiply the J used for this cross-section by the Quantum Yield 
        # that's noted in the FJX file for this rxn. Place that info in output df! 
        q_yield=fjx_df.loc[ind, 'Quantum_Yield']
        j_for_crossec=fjx_df.loc[ind, 'J_for_CrossSec']
        if q_yield != 1:
            fjx_df.at[ind, 'J_for_F0AM'] = f'{str(q_yield)}.*({j_for_crossec})'
        else:
            fjx_df.at[ind, 'J_for_F0AM'] = j_for_crossec
            
    ###########################################################################
    # Assign Unique Nicknames for Each dif. Cross-Section: 
    ########################################################################### 
    # Get a list of all duplicate j-nicknames in FJX
    unique_Jnames = set(); dupe_Jnames = []
    for ji, j in enumerate(fjx_df.loc[:, 'J-NickName']):
        if j in unique_Jnames and j not in dupe_Jnames:
            dupe_Jnames.append(j)
        elif j not in unique_Jnames:
            unique_Jnames.add(j)
            
    # Loop oveer all duplicate J-Names: 
    for j in dupe_Jnames:
        # Get ind in fjx of all matches to a dupe name
        idx_matches = fjx_df.index[j == fjx_df.loc[:, 'J-NickName']]
        # Figure out what J-value / quantum yield this specific match points to: 
        assigns = [fjx_df.loc[ix, 'J_for_F0AM'] for ix in idx_matches]
        # Figure out if this assignment is the same as the first occurance or not. 
        is_diff = [True if assignment == assigns[0]
                   else False for assignment in assigns]

        # If some of the assignments vary for the same J-Nickname (e.g. same name points to dif values)
        if not all(is_diff):
            # Loop over all mataches to this nickname
            for ii, ix in enumerate(idx_matches):
                # If it doesn't match the assignment of the one already defined (1st occurance) or is 1st occurance...
                if (is_diff[ii] == False) or (ii == 0):
                    # Update the name to be unique!
                    new = make_unique_Jname(j, list(fjx_df.loc[:, 'J-NickName']))
                    fjx_df.at[ix, 'J-NickName'] = new
                    
    ###########################################################################
    # Now, create output dictionary mapping what's in kpp file (PHOTOL(INDX)) o 
    # the J-Name used in the F0AM 'GESOChem_Rxns.m' file & the J value assigned 
    # to that name in the 'F0AM GEOSChem_J.m 'file: 
    ########################################################################### 
    for ix in fjx_df.index:
        # Pull out info on this row only to make crafting strings readable...  
        fjx_ind_i= fjx_df.loc[ix, 'FJX_Index']
        photo_rxn_i= fjx_df.loc[ix, 'Photolysis_Rxn'].replace(' ', '')
        cx_i =fjx_df.loc[ix, 'FJX_CrossSec']
        what_is= '(aka '+tracer_info[cx_i]+')' if cx_i in list(tracer_info.keys()) else '' 
        qy_i= fjx_df.loc[ix, 'Quantum_Yield']
        name_i= fjx_df.loc[ix, 'J-NickName']
        
        # Add coulmn that tells us what this reaction's index is in the KPP File: 
        fjx_df.at[ix, 'PHOTOL(IND)'] = f"PHOTOL({  str(fjx_ind_i)  })"
        
        # Add column telling F0AM User what photolysis rxn / PHOTOL_Index / 
        # cross-Sec/ QY was used here: 
        fjx_df.at[ix, 'Info']=f"% PHOTOL({ str(fjx_ind_i) }) = {name_i} & used for {photo_rxn_i}" 
                                   # based on cross-section of {cx_i} {what_is} with QY={str(qy_i)}."  

    # Create a dictionary to hold results of J-Mapping: 
    jdict={
        # 'Inds_to_JNames' maps PHOTOL(IND) in KPP to the J-NickName in F0AM used in GEOSChem_Rxns.m'
        'Inds_to_JNames': dict(zip(fjx_df['PHOTOL(IND)'], ["'" + j + "'" for j in fjx_df['J-NickName']])),
        
        # 'JNames_to_Vals' maps the J-Nickname to the F0AM assignment used in 'GEOSChem_J.m'
        'JNames_to_Vals':  {fjx_df.loc[ix, 'J-NickName']: 
                                 {'Assignment': fjx_df.loc[ix,'J_for_F0AM'], 
                                  'Info':fjx_df.loc[ix, 'Info'],
                                  'Photo_Indx': fjx_df.loc[ix, 'FJX_Index']}  for ix in fjx_df.index}
        }
    
    
    return jdict, fjx_df,user_inp




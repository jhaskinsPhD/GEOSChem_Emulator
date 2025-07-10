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
from prep_for_foam import prep4f0am 
import photolysis as hv 
import utils


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
def make_import_GC_rates_file(user_inp: dict, all_rfuncts: list,  
                              output_fname: str = 'import_GC_rates_v1452',
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
    # Unpack vars from user_inp dict referenced herein: 
    template_file=user_inp['input_paths']['Template_Files']['import_GC_rates']
    rate_files=user_inp['rate_files']
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
        '    [temp, press, numden, h2o,temp_over_k300, k300_over_temp,inv_temp,sr_temp,relhum]= extract_MetVars(Met);\n'
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
                
                
                if i!=0 and line.strip()[0]!="%": 
                    rfuncts_in_it=[funct for funct in all_rfuncts if funct+'(' in line.strip()]
                    line=line.lower()
                    if len(rfuncts_in_it) >0: 
                        for r in rfuncts_in_it: 
                            line= line.replace(r.lower()+'(', r+'(Met,')
                    outF.write(line+'\n')
                else: 
                    outF.write(line+'\n')
                if glbs_tf is False and ')' in line:
                    outF.write(get_globals+'\n')
                    glbs_tf = True

    outF.close()  # Close the output file
    print('MATLAB compliant version of gckpp_Rates.F90 called in "GEOSCHEM_K.m" to define all rate functions saved at: \n\t'+rate_file)
    
    # print('MISSING:') 
    missing=[print(f'\t{f}()') for f in all_rfuncts if f not in written] 
    return missing





# %%###########################################################################
# ---------- CREATE THE GEOSCHEM_K.m and GEOSChem_J.m Files -------------------
###############################################################################
def make_GEOSChem_J_file(user_inp:dict, 
                         foam_dict: dict, 
                         fjx_df, 
                         jmap_dict:dict, 
                         output_fname: str = 'GEOSCHEM_J_v1452.m',
                         output_dir: str = '', overwrite: bool = False):  
    
    # Unpack vars from user_inp dict referenced herein: 
    kppfile= user_inp['kppfile'] 
    GC_version= user_inp['GC_version'] 
    fjx_version= user_inp['FJX_version'] 
    template_file=user_inp['input_paths']['Template_Files']['J']
    include_het_rxns= user_inp['include_het_rxns']
    
    # Create a formatted list of all the required named J's. This info goes 
    #IN THE HEADER (commented out). Format it so its in a commented out comma 
    # delimited list with line breaks where  necessary (for long MATLAB lists 
    # that continue on multiple lines). 
    nhet_hdr_j_list= utils.join_list_for_MATLAB(',', foam_dict['nhet']['req_js'],
                                           insert='    ',  adj_ln1_width=33,
                                           insert_skip=[1],comment=True, 
                                           comment_skip=[1]) 
    het_hdr_j_list= utils.join_list_for_MATLAB(',', foam_dict['het']['req_js'],
                                       insert='    ',  adj_ln1_width=33,
                                       insert_skip=[1],comment=True, 
                                       comment_skip=[1]) 
   
    # Create lists to hold the lines assiging Js to het species, 
    # with direct MCM analogoues , & others (seperately) 
    mcm_js=[]; jvs=[]; jns=[]; het_js=[] 
    
    # Loop over all the Js required for (all) Rxns: 
    all_jnames= [j_name.replace("'",'') for j_name in foam_dict['nhet']['req_js']+foam_dict['het']['req_js']]
    all_jinds=[jmap_dict['JNames_to_Vals'][j_name]['Photo_Indx'] for j_name in all_jnames] 
    # Zip the lists together and sort the jnames based on the values in PHOTOL() 
    sorted_jnames = [item for ind, item in sorted(zip(all_jinds, all_jnames))  ]
    
    for j_name in sorted_jnames:
        # Pull out assigned value rom jmpa_dict: 
        assignment=jmap_dict['JNames_to_Vals'][j_name]['Assignment'] 
        info=jmap_dict['JNames_to_Vals'][j_name]['Info'] 
        
        # Fix name of assignment we're using so it is a sub-field of a struct named 'Jmcm'... 
        #( e.g. 'J1' turns into 'Jmcm.J1' and 'J13+J14' turns into 'Jmcm.J13+Jmcm.J14' etc.) 
        if 'Jn' in assignment: 
            assignment=assignment.replace('Jn','Jmcm.Jn') 
            mcm_analog=False; is_Jn=True
        elif 'Jv' in assignment:
            assignment=assignment.replace('Jv','Jmcm.Jv') 
            mcm_analog=False; is_Jn=False
        elif 'J' in assignment: 
            assignment=assignment.replace('J','Jmcm.J')
            mcm_analog=True
        else: 
            raise ValueError(f'IDK what to do with this J-assignment: {j_name}=assignment; {info}').with_traceback(sys.exc_info()[2])

        # Format this line for the F0AM GEOSChem_J. file:      
        line=f"J.{j_name} = {assignment}; {info}" 
        
        # Organize J's by whether het, have an MCM analogue or not... 
        if "'"+j_name+"'" in foam_dict['het']['req_js']: 
            # Comment out assignment for het J's, but still add to file: 
            if include_het_rxns is False: line='% '+line
            # Since its photolysis of a het species, keep it in its own list:
            het_js.append(line); 
        elif mcm_analog==True: 
            mcm_js.append(line)
        elif is_Jn==True: 
            jns.append(line)
        else:
            jvs.append(line)
            
    
    # Figure out how many lines there are in the "Template" file so you know when 
    # you reach the end.
    num_lines = sum(1 for line in open(template_file))
    
    # Using readlines(), open the template file and read in line by line.
    inF = open(template_file, 'r')
    lines = inF.readlines()
    count = 0  # Initialize line counter variable.
    
    jlines_to_write=[] 
    
    # Change a few parts of the file to read dif if het rxns are True/False: 
    if include_het_rxns==False:
        het_rx='FALSE'; not_or_blank1= ''; not_or_blank2= 'NOT'; extra=', but are retained here for reference'
    else: 
        het_rx='TRUE'; not_or_blank1= 'NOT'; not_or_blank2= 'indeed '; extra=' and/or GEOSChem_HetRxns.m file'
        
        
    for line in lines:
        
        # Add to HEADER what version of GEOS-Chem this mechanism is based on to header: 
        line = utils._edit_line(line, '**GC_VERSION**', GC_version) 
        
        line = utils._edit_line(line, '**KPP_FILE_PATH**', kppfile) 
        
        # Add to HEADER what date this fiel was generated on in the function header
        line = utils._edit_line(line, '**THE_DATE**', str(date.today())) 
        
        # Add to HEADER what version of GEOS-Chem this mechanism is based on to header: 
        line = utils._edit_line(line, '**FJX_FILEPATH**', fjx_version) 

        # Add lines with the MCM direct analogue J-Values here: 
        line = utils._edit_line(line, '**HDR_LIST_OF_REQUIRED_J_NAMES**', nhet_hdr_j_list) 

        # Add lines with the MCM direct analogue J-Values here: 
        line = utils._edit_line(line, '**INSERT_NHET_MCM_Js**', '\n'.join(mcm_js)) 
        
        # Add lines with the non-MCM direct analogue J-Values here: 
        line = utils._edit_line(line, '**INSERT_NHET_OTHER_Js**', '\n'.join(jns)+'\n'+'\n'.join(jvs)) 
        
        # Edit line(s) about comment on Het rxns being commented out or not... 
        line = utils._edit_line(line, '**NOT_OR_BLANK1**', not_or_blank1) 
        line = utils._edit_line(line, '**INCLUDE_HET_RXNS**', het_rx) 
        line = utils._edit_line(line, '**NOT_OR_BLANK2**', not_or_blank2) 
        line = utils._edit_line(line, '**EXTRA**', extra) 
        
        # Add lines with the non-MCM direct analogue J-Values here: 
        line = utils._edit_line(line, '**INSERT_HET_Js**', '\n'.join(het_js)) 
        
        # Append possibly modified line to list to write to output file: 
        jlines_to_write.append(line) 
        
        count += 1  # Update counter.

        if count > num_lines:
            break  # Exit while loop if you reach end of file.
            
    # Close the template file: 
    inF.close() 
    
    # Decide what to call the output file & make sure the path to it exists...
    j_file = utils.check_filename(filename=output_fname, default_name='GEOSCHEM_J_v1452.m', 
                                         ext='.m',savepath=output_dir, 
                                         overwrite=overwrite, return_full=True)

    # Open the Output file to Write
    out_JF = open(j_file, "w")

    # Write the header lines to the file: 
    [out_JF.write(line) for line in jlines_to_write]
    
    out_JF.close()  # Close the output file
    print('GEOSCHEM_J.m file for F0AM saved at: \n\t'+j_file)
    
    return user_inp

def edit_K_template(user_inp:dict, ks_required:list, rxn_rate_dict:dict,
                    GC_Rxn_File:str):
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
    (1) user_inp - DICT with a minumum of the following keys: 
        'kppfile'       - STR with full path to the KPP input file mech was built from. 

        'GC_version'    - STR formatted as '13.3.4' corresponding to what GEOS-Chem 
                        version you are seeking to create a F0AM file for. 
                        
        'Template_files' - DICT pointing to paths to template files. 

    (2) ks_required - Unique list of all K's required in file
    
    (3) rxn_rate_dict  - DICT with reaction rate names in F0AM as keys, assignment in 
                         GEOSChem_k.m file as values (including "Met" as arg!) 

    (4) GC_Rxn_File    - STR with name of GEOSChem_Rxns.m file this rates file 
                         corresponds to. 


    OUTPUTS:
    --------
        A F0AM compliant rate defintion file saved at output_dir+output_fname. 

    REQUIREMENTS: 
    -------------
        LIBRARIES:           from datetime import date  

        CUSTOM FUNCTIONS:    utils.join_list_for_MATLAB()
                             utils.check_filename()
                             utils._edit_line() 


    AUTHOR: 
    -------
            Dr. Jessica D. Haskins (jessica.haskins@utah.edu) GitHub: @jhaskinsPhD


    """
    # Create a formatted list of all the required rate functions (K's) that must 
    # be imported & named J's. This info goes IN THE HEADER (commented out). Format 
    # its in a commented out comma delimited list with line breaks 
    # where necessary (for long MATLAB lists that continue on multiple lines). 
    
    k_functs= [k+'()' for k in ks_required] # Add '()' after func name... 
    hdr_k_list= utils.join_list_for_MATLAB(',', k_functs,
                                           insert='        ', 
                                           adj_ln1_width=38,
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
    # Read in GC_Rates Template and sub in specific info for our version of the file.
    ###########################################################################
    # Initialize list of lines we'll write to the  output file...
    all_k_lines=list()
    
    # Get path to K file template: 
    template_file= user_inp['input_paths']['Template_Files']['K'] 
    
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
        line = utils._edit_line(line, '**GC_VERSION**', user_inp['GC_version'] ) 
        
        # Add to HEADER what date this fiel was generated on in the function header
        line = utils._edit_line(line, '**THE_DATE**', str(date.today())) 
        
        # Add to HEADER what date this fiel was generated on in the function header
        line = utils._edit_line(line, '**KPP_FILE_PATH**', user_inp['kppfile']) 
       
        # Add to HEADER what version of GEOS-Chem this mechanism is based on to header: 
        line = utils._edit_line(line, '**GC_RXN_FILE**', GC_Rxn_File) 
        
        # Add to HEADER which functions are required for rate constants: 
        line = utils._edit_line(line, '**HDR_LIST_OF_REQUIRED_RATE_FUNCTIONS**', hdr_k_list) 
    
        # Update actual MATLAB CODE to import only the handles of the rate
        # functions used in the mechanism in the rate dictionary 
        line = utils._edit_line(line, '**LIST_OF_RATE_FUNCTS_TO_IMPORT**', import_rate_list) 
    
        # Add the line that unpacks the handles from the rate dictionary for use
        line = utils._edit_line(line, '**LIST_OF_IMPORTED_FUNCTION_HANDLES**', rate_list) 
        
        # Add number of Ks required to create correct len arrays in F0AM. 
        line = utils._edit_line(line, '**NUMBER_OF_RATE_CONSTANTS**', str(int(len(all_k_names)))) 
        
        # Add all lines defining rates to the file in the rate def section: 
        line = utils._edit_line(line, '**INSERT_RATE_DEFINITIONS**',rate_defs ) 
       
        # Add (now possibly modified) line to list of lines to write: 
        all_k_lines.append(line)
        
        count += 1  # Update counter.

        if count > num_lines:
            break  # Exit while loop if you reach end of file.

    # Close the template file now that you got all the good stuff out.
    inF.close()
   
    return all_k_lines


def make_GEOSChem_K_file(user_inp:dict, foam_dict: dict, output_fname: str = 'GEOSCHEM_K_v1452.m',
                         output_dir: str = '', overwrite: bool = False):
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
    (1) user_inp - DICT with a minumum of the following keys: 
        'kppfile'       - STR with full path to the KPP input file mehch was built from. 

        'GC_version'    - STR formatted as '13.3.4' corresponding to what GEOS-Chem 
                        version you are seeking to create a F0AM file for. 

        'Template_Files'- DICT containing STRING w/ the full path to the 
                          GEOSChem_K_template.txt template files

    (2) foam_dict - DICT with F0AM mechanism info 
    
    (3) output_fname  - (OPTIONAL) STR with name of the outputfile to write. 
                        Default is 'GEOSCHEM_K.m'

    (4) output_dir    - (OPTIONAL) STR with full path where the outputfile should 
                         be written. Default is to "GC_Emulator/Output_Files/"

    (5) overwrite     - (OPTIONAL) BOOL of whether or not to overwite the output 
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
    ###########################################################################
    #   Write the  GEOSCHem_K.m File (for all Non-Heterogeneous rxns): 
    ###########################################################################
    # Pull out only the reaction rates names in F0AM/ assingments used in the non-het mech: 
    nhet_rxn_rates=foam_dict['nhet']['rate_dict']

    # Names of unique) rate functions required in non-het mechanism: 
    nhet_ks_req=foam_dict['nhet']['req_ks']
    
    # Read in the K template & edit with info on all the non-het rxn rates: 
    nhet_k_lines= edit_K_template(user_inp=user_inp, 
                                  ks_required=nhet_ks_req, 
                                  rxn_rate_dict=nhet_rxn_rates,
                                  GC_Rxn_File=user_inp['GEOSCHEM_GasRxns'])
    
    
    # Decide what to call the output file & make sure the path to it exists...
    nhet_k_file = utils.check_filename(filename=output_fname, default_name='GEOSChem_K_v1452.m', 
                                         ext='.m',savepath=output_dir, 
                                         overwrite=overwrite, return_full=True)
    
    # Open the Output file to Write
    out_gaskF = open(nhet_k_file, "w")
    
    # Write the (customized) lines to the file: 
    [out_gaskF.write(line) for line in nhet_k_lines]
    
    # Close the new KK file and print where it is: 
    out_gaskF.close()
    print('GEOSCHEM_K.m file defining reaction rates for non-het F0AM rxns saved at: \n\t'+nhet_k_file)
    
    # Add path to user_inp dict: 
    user_inp['GEOSCHEM_K']=nhet_k_file
   
    ###########################################################################
    #   (IF ASKED... ) Write the  GEOSCHem_K.m File (for al-Heterogeneous rxns): 
    ###########################################################################
    if user_inp['include_het_rxns'] is True: 
        # Pull out only the reaction rates names in F0AM/ assingments used in the non-het mech: 
        het_rxn_rates=foam_dict['het']['rate_dict']

        # Names of unique) rate functions required in non-het mechanism: 
        het_ks_req=foam_dict['het']['req_ks']
        
        # Read in the K template & edit with info on all the non-het rxn rates: 
        het_k_lines= edit_K_template(user_inp=user_inp, 
                                      ks_required=het_ks_req, 
                                      rxn_rate_dict=het_rxn_rates,
                                      GC_Rxn_File=user_inp['GEOSCHEM_HetRxns'])
        
        
        # Decide what to call the output file & make sure the path to it exists...
        het_k_file = utils.check_filename(default_name='GEOSChem_HET_K_v1452.m', 
                                             ext='.m',savepath=output_dir, 
                                             overwrite=overwrite, return_full=True)
        
        # Open the Output file to Write
        out_hetkF = open(het_k_file, "w")
        
        # Write the (customized) lines to the file: 
        [out_hetkF.write(line) for line in het_k_lines]
        
        # Close the new KK file and print where it is: 
        out_hetkF.close()
        print('GEOSChem_HET_K.m file defining rxn rates used in heterogeneous F0AM rxn rates saved at: \n\t'+het_k_file)
    
        # Add path to user_inp dict: 
        user_inp['GEOSCHEM_HET_K']=het_k_file
    
    return user_inp
 
                
# %%###########################################################################
# ----------------- CREATE THE GEOSCHEM_Rxns.m File ----------------------------
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


def make_GEOSChem_Rxns_file(user_inp:dict, foam_dict:dict, output_dir: str = '', 
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
                                      template_file=user_inp['input_paths']['Template_Files']['GasRxns'], 
                                      het_species_list= foam_dict['het']['sp_list'])
                                                                     
    # Decide what to call the output file & make sure the path to it exists...
    gas_rxns_file = utils.check_filename(default_name='GEOSCHEM_GasRxns_v1452.m', 
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
    
    # Add path to user_inp: 
    user_inp['GEOSCHEM_GasRxns']=gas_rxns_file
    

    if user_inp['include_het_rxns'] is True:
        ###########################################################################
        #   Write the  GEOSCHem_HetRxns.m File (for all Het rxns): 
        ###########################################################################
        het_header=edit_GCRxns_template(user_inp, foam_dict['het'], 
                                           has_excluded_rxns=has_excluded_rxns, 
                                           template_file=user_inp['input_paths']['Template_Files']['HetRxns'])
            
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
        
        # Add path to user_inp: 
        user_inp['GEOSCHEM_HetRxns']=het_rxns_file
        
    return user_inp
        

# %%###########################################################################
# ----------------- Top level functs referenced in make_GC_mchanism()----------
############################################################################### 
# Function used in make_GC_mechanism() to remove het only rxns in non-het mech: 
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

# Function used in make_GC_mechanism() to check user inputs to make_GC_mechanism()
def check_user_inputs(kppfile:str, rate_files, GC_version: str, jmap_type: str, 
                       output_dir: str, include_het_rxns:bool, verbose:bool) :
    """ Function to check all user inputs to the top level GEOS-Chem Emulator function""" 
   
    # Get the relative paths to input/templates, etc...
    my_paths = utils.get_my_relative_paths()
    
    # Check that the output directory the user entered exists:
    if os.path.exists(output_dir) is False:
        raise NotADirectoryError("The following output directory you listed to store the output files does not exist: \n" +
                                 f"\toutput_dir = '{output_dir}'").with_traceback(sys.exc_info()[2])

    # Check that the KPP files passed as input exist
    if os.path.isfile(kppfile) is False:
        raise FileNotFoundError("The path to the input KPP file you entered does not exist: \n" +
                                f"\tkppfile = '{kppfile}'").with_traceback(sys.exc_info()[2])

    # Check that the input reaction rates files passed as input are the right type & exist:
    if type(rate_files) not in [type('str'), type(list([1]))]:
        raise TypeError("Input for argument 'rate_files' must either be TYPE = STR or LIST: \n" +
                        "\t1) A STR with a path to a single .F90 file containing reaction rates \n" +
                        "\t(2) A LIST of strings with paths to rate function files like: \n" +
                        "\t\t['/path/to/gckpp_Rates.F90', '/path/to/fullchem_RateLawFuncs.F90'] \n" +
                        f"You entered input for 'rate_files' of type={str(type(rate_files))}").with_traceback(sys.exc_info()[2])
    
    # If its a single rate file convert it to a len 1 list from here out.  
    elif type(rate_files) == type('str'): 
        rate_files=[rate_files] 
   
    # Check that each file in the (now list), rate_files does indeed exist: 
    if not all([os.path.isfile(file) for file in rate_files]):
        err_files = [f"\t({str(i+1)}) '{file}'" for i,file in enumerate(rate_files) if os.path.exists(file) is False]
        raise FileNotFoundError("The path to the following reaction rate file(s) you entered do not exist: \n" +
                                '\n'.join(err_files)).with_traceback(sys.exc_info()[2])

    # Check that the Template & Phototlysis files exist at the expected paths: 
    for ftype in ['Template_Files', 'Photolysis_Files']:
        these_paths= my_paths[ftype]
        for key in list(these_paths.keys()):
            if os.path.exists(these_paths[key]) == False:
                raise FileNotFoundError(f"The path to the {key} {ftype.split('_')[0].lower()}"+
                                        "file was not found where expected at: \n" +
                                        f"\t'{these_paths[key]}'\n" +
                                        f"It should be located at: '.../GEOSChem_Emulator/Input_Files/{ftype}/{key}'" 
                                        ).with_traceback(sys.exc_info()[2])

    # Check that the user input for "jmap_type" is allowed/ recognized: :
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
            raise ValueError("Input for argument 'jmap_type' must be one of the following: \n" +
                             "    (1) 'CrossSec_Match' \n" +
                             "    (2) 'MCM' \n" +
                             "    (3) 'RCIM' \n" +
                             " For more info on these options do >>hv.jmap_type_help()").with_traceback(sys.exc_info()[2])

    # Store all the user input in a nice verified dictionary (passed to many functions)
    user_inp={'kppfile':kppfile, 
              'rate_files':rate_files,
              'GC_version': GC_version,
              'jmap_type': jmap_type,
              'include_het_rxns':include_het_rxns,
              'input_paths': my_paths, 
              'output_dir':output_dir,
              'verbose':verbose}
        
    return user_inp

# %%###########################################################################
# -----------------Main Wrapper Function to Make ALL Files --------------------
###############################################################################
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
    

    # Parse the KPP file. Get list of rxns and rates we need to write a F0AM file.
    print('Parsing KPP file...') 
    tracer_info, fixed_vars, kpp_dict = read_kpp(kppfile, output_dir=output_dir)

    
    # Read in the right FJX file for this GC-Version and build the j-mapping
    # dictionary & FJX dataframes to use in mapping KPP J's to F0AM J's.
    jmap_dict, fjx_df, user_inp = hv.create_jmap_dict(user_inp, tracer_info)
    
    
    # If NOT including heterogeneous reactions in output mech, we should be able to 
    # eliminate some reactions from the mechanism (commenting them out). (e.g. Why
    # photolyze ClNO2 if there is none ever made b/c formation is het only?).
    inds2drop=[] # Inialize list of inds to drop from kpp_dict... 
    if include_het_rxns is False:
        inds2drop=eliminate_dead_ends(kpp_dict, verbose=user_inp['verbose'])
    
    # Prep reactions & rates for F0AM (e.g. turn kpp_dict info into format F0AM 
    # The output "foam dict" has 2 top level keys: 'het' & 'nhet' w/ same info 
    # about stuff that will go in the het/nonhet mechanism files (if het is requested): 
    print('Preparing for F0AM...')
    foam_dict=prep4f0am(kpp_dict, inds2drop, fixed_vars, jmap_dict)

    # # Save Tracer info to an outputfile too: 
    # Check that the filename is available & return the full path to write the file: 
    tracer_file=utils.check_filename(default_name= 'TracerInfo', ext='.csv', 
                        savepath=output_dir, overwrite=overwrite, return_full=True, quiet=True)
    
    # Compile all info about tracers into a df to write to an output file: 
    tracers=fixed_vars+sorted(list(tracer_info.keys())); 
    tinfo =[tracer_info[key] for key in tracers]
    is_fixed =[True if key in fixed_vars else False for key in tracers]
    het_only_sp=[True if key in foam_dict['nhet']['sp_list'] else False for key in tracers] 
    
    tracer_df=pd.DataFrame({'Tracers': tracers, 'Information':tinfo, 'Is_Fixed' :is_fixed, 'Het_Only_Species':het_only_sp,})
    tracer_df.to_csv(tracer_file)


    # Write the GEOSChem_GasRxns.m file (&, if asked, the GEOSChem_HetRxns.m file). 
    user_inp=make_GEOSChem_Rxns_file(user_inp=user_inp,
                            foam_dict=foam_dict, 
                            output_dir=output_dir,
                            overwrite=overwrite)
    
    # Write the GEOSChem_K.m file  
    user_inp=make_GEOSChem_K_file(user_inp=user_inp,
                         foam_dict=foam_dict,
                         output_fname='GEOSCHEM_K_v1452.m',
                         output_dir=output_dir,
                         overwrite=overwrite)
        
    # Write the GEOSChem_J.m file  
    user_inp=make_GEOSChem_J_file(user_inp=user_inp,
                         foam_dict=foam_dict,
                         fjx_df=fjx_df,
                         jmap_dict=jmap_dict, 
                         output_fname='GEOSCHEM_J_v1452.m',
                         output_dir=output_dir,
                         overwrite=overwrite)

    # # Create import_GC_rates.m file where all rate functions are defined which is
    # # referenced and used as is without needing modifications in the GEOSChem_K.m file.
    all_rfuncts=foam_dict['nhet']['req_ks']
    fdict = make_import_GC_rates_file(user_inp=user_inp,
                                      all_rfuncts=all_rfuncts,
                                      output_fname='import_GC_rates_v1452.m',
                                      output_dir=output_dir,
                                      overwrite=overwrite)


    
    return foam_dict,kpp_dict,inds2drop,fixed_vars, fjx_df, jmap_dict ,




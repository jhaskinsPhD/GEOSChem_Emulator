#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 13:05:24 2023

@author: u6044586
"""

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



#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 12 12:01:38 2025

@author: u6044586
"""

import os
import shutil


def list_directories(path):
    ls= [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]
    return  ls

def create_directory(path, directory_name):
    # Construct the full path for the new directory
    new_directory_path = os.path.join(path, directory_name)
    
    # Check if the directory already exists
    if not os.path.exists(new_directory_path):
        # Create the directory
        os.makedirs(new_directory_path)

    return new_directory_path


# Get a list of all GEOS-Chem versions I have: 
fldr='/uufs/chpc.utah.edu/common/home/haskins-group1/users/jhask/GEOSChem/'
vers=list_directories(fldr)
vers.remove('GC_RunDirs') # Remove run dir folder... 

# Get version numbers only from folder names: 
vnums_only=[v.replace('GCClassic_','') for v in vers] 

# Create dirs to hold kpp files for each version if they don't already exist: 
inpt_pth='/uufs/chpc.utah.edu/common/home/u6044586/python_scripts/modules/GEOSChem_Emulator/Input_Files/GC_Mech_Files'
dest_pths={}
for v in vnums_only: 
    pth_i=create_directory(inpt_pth, v)
    dest_pths['GCClassic_'+v]=pth_i
    
for v in vers: 
    # Set path to kpp folder for this version: 
    kpp=os.path.join(fldr,v,'src/GEOS-Chem/KPP/fullchem')
    cloud_j= os.path.join(fldr,v,'src/Cloud-J/tables/FJX_j2j.dat')
    
    # Copy kpp .eqn file and gckpp rates over to folders: 
    shutil.copy(os.path.join(kpp,'fullchem.eqn'), os.path.join(dest_pths[v],'fullchem.eqn'))
    shutil.copy(os.path.join(kpp,'gckpp_Rates.F90'), os.path.join(dest_pths[v],'gckpp_Rates.F90'))
    shutil.copy(os.path.join(kpp,'fullchem_RateLawFuncs.F90'), os.path.join(dest_pths[v],'fullchem_RateLawFuncs.F90'))
    if v not in ['GCClassic_v14.2.3','GCClassic_v14.0.0']:
        shutil.copy(cloud_j, os.path.join(dest_pths[v],'FJX_j2j.dat'))

    

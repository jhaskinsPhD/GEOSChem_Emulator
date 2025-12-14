#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 10 21:26:24 2025

@author: u6044586
"""

from rdkit import Chem

class Reaction:
    def __init__(self, reaction_str, name_map, lump_map=None):
        """
        reaction_str: e.g., 'TOLU + OH = TRO2 + HO2'
        name_map: dict of tracer_name -> SMILES or InChI
        lump_map: dict of lumped_name -> list of component names
        """
        self.reaction_str = reaction_str
        self.name_map = name_map
        self.lump_map = lump_map or {}
        self.reactants, self.products = self.parse_reaction()
    
    def parse_reaction(self):
        # Parse the reaction string into reactants and products lists
        reactant_part, product_part = self.reaction_str.split('=')
        reactant_names = [r.strip() for r in reactant_part.split('+')]
        product_names = [p.strip() for p in product_part.split('+')]
        # Map to compounds (Mol objects)
        reactant_mols = [self.get_molecule(name) for name in reactant_names]
        product_mols = [self.get_molecule(name) for name in product_names]
        return reactant_mols, product_mols
    
    def get_molecule(self, name):
        # Handle lumped groups
        if name in self.lump_map:
            # optionally, handle lumped groups as dummy molecules
            return self.create_lumped_molecule(name)
        smiles_or_inchi = self.name_map.get(name)
        if smiles_or_inchi is None:
            raise ValueError(f"Unknown compound: {name}")
        mol = Chem.MolFromSmiles(smiles_or_inchi) or Chem.MolFromInchi(smiles_or_inchi)
        return mol
    
    def create_lumped_molecule(self, lump_name):
        # Create dummy molecule with label (or from constituent compounds)
        # For simplicity, create a molecule with a label
        from rdkit.Chem import rdmolops
        mol = Chem.MolFromSmiles("C")  # Placeholder
        mol.SetProp("Lumped", lump_name)
        return mol
    
    def get_participating_compounds(self):
        # Return set/list of compound names (or SMILES, etc.)
        names = []
        for mol in self.reactants + self.products:
            name = mol.GetProp("Lumped") if mol.HasProp("Lumped") else None
            if name:
                names.append(name)
            else:
                # In real code, store a property with the name
                pass
        return set(names)
    
    def is_equivalent(self, other):
        """Compare reactions based on molecule identities, ignoring order and naming."""
        # Compare sets of mol InChI or SMILES
        def mols_id(mols):
            return set(Chem.MolToInchi(mol) for mol in mols)
        return (mols_id(self.reactants) == mols_id(other.reactants) and
                mols_id(self.products) == mols_id(other.products))

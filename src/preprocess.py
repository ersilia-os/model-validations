import os
import pandas as pd
from rdkit import Chem

def is_processable_smiles(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            return True
        else:
            return False
    except Exception as e:
        print("Error:", e)
        return False    

def get_inchikey(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    else:
        # Generate InChIKey
        inchi_key = Chem.InchiToInchiKey(Chem.MolToInchi(mol))
        return inchi_key  

def canonical_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    else:
        # Generate canonical SMILES
        canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=False)
        return canonical_smiles    

def standardize_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)

    if mol is not None:
        try:
            mol = standardise.run(mol)  # standardise.run expects a list of molecules
        except:
            mol = None

    if mol is not None:
        smiles_standard = Chem.MolToSmiles(mol)
    else:
        smiles_standard = None

    return smiles_standard

def preprocess_eos_data(eos_df, smi_col="Smiles", act_col="activity10"):
    print("Inside preprocess_eos_data function...")

    # Drop rows with missing InChI Keys
    eos_df = eos_df[eos_df[smi_col].apply(is_processable_smiles)]
    initial_len = len(eos_df)

    # Generate InChI Keys
    eos_df['inchikey'] = eos_df[smi_col].apply(get_inchikey)

    # Drop rows with missing InChI Keys
    eos_df.dropna(subset=['inchikey'], inplace=True)
    eliminated = initial_len - len(eos_df)
    print("Number of SMILES eliminated:", eliminated)

    # Rename columns
    eos_df.rename(columns={smi_col: "smiles", act_col: "activity"}, inplace=True)

    # Select relevant columns
    eos_df = eos_df[["smiles", "inchikey", "activity"]]

    return eos_df


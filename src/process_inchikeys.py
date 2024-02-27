# process_inchikeys.py
import os
import pandas as pd
from rdkit import Chem

def preprocess_eos_data(eos_df, smi_col="Smiles", act_col="activity10"):
    print("Inside preprocess_eos_data function...")
    
    # Generate InChI Keys
    inchikeys = []
    for smi in eos_df[smi_col]:
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            inchikey = Chem.MolToInchiKey(mol)
            inchikeys.append(inchikey)
    
    # Drop rows with missing InChI Keys
    initial_len = len(eos_df)
    eos_df['inchikey'] = inchikeys
    eos_df.dropna(subset=['inchikey'], inplace=True)
    eliminated = initial_len - len(eos_df)
    print("Number of SMILES eliminated:", eliminated)
    
    # Rename columns
    eos_df.rename(columns={smi_col: "smiles", act_col: "activity"}, inplace=True)
    
    # Select relevant columns
    eos_df = eos_df[["smiles", "inchikey", "activity"]]
    
    return eos_df


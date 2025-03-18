import pandas as pd
import numpy as np

import Bio
from Bio.PDB import PDBParser, mmcifio
from Bio.PDB import alphafold_db
import os

def fetch_alphafold(uniprot_id, rename = None):
    '''Fetches the AlphaFold 3D structure of a protein from the AlphaFoldDB.'''
    # fetch alphafold file
    try:
        prediction = alphafold_db.get_predictions(uniprot_id)
        af_dict = next(prediction)
    except: 
        return False
    alphafold_file_path = alphafold_db.download_cif_for(af_dict, directory="alphafold_files")

    if not alphafold_file_path:
        return False
    if rename:
        os.rename(alphafold_file_path, f"alphafold_files/{rename}.cif")
        alphafold_file_path = f"alphafold_files/{rename}.cif"

    return True

def main():

    # open tsv file
    pdb_data = pd.read_csv("low_coverage_pdbs.tsv", sep="\t")
    print(f"Number of PDB IDs: {len(pdb_data)}")

    # fetch alphafold files
    failed = []
    for index, row in pdb_data.iterrows():
        print(f"{index + 1}. ({row["UniProt ID"]}) ", end = "")
        result = fetch_alphafold(row["UniProt ID"], rename = row["UniProt ID"])
        print(result)
        if not result:
            failed.append(row["UniProt ID"])

    print(f"Done. Number of Failed Fetches: {len(failed)}")
    # save failed as list
    with open("failed_alphafold_fetches.txt", "w") as f:
        for item in failed:
            f.write(f"{item}\n")

    #print(fetch_alphafold("Q9NYG8"))

if __name__ == "__main__":
    main()

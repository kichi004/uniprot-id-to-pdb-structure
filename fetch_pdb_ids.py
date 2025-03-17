import Bio.PDB
import pandas as pd
import time

def fetch_pdb_id(pdb_id):
    '''Fetches the 3D structure of a protein from the PDB.'''
    pdbl = Bio.PDB.PDBList()
    pdb_file_path = pdbl.retrieve_pdb_file(pdb_id, pdir="pdb_files", file_format="mmCif")
    time.sleep(1)
    return pdb_file_path

def main():
    threhold = 0.7

    # organize found pdb ids by coverage
    pdb_data = pd.read_csv("result_pdbs.tsv", sep="\t")
    # sort by coverage
    pdb_data = pdb_data.sort_values(by="Coverage", ascending=False)
    # filter into two dataframes
    high_coverage_pdbs = pdb_data[pdb_data["Coverage"] >= threhold]
    low_coverage_pdbs = pdb_data[pdb_data["Coverage"] < threhold]
    
    # get list of pdb ids
    pdb_ids_list = high_coverage_pdbs["PDB ID"].tolist()
    print(f"Number of High Coverage PDB IDs: {len(pdb_ids_list)}")
    print(f"Number of Low Coverage PDB IDs: {len(low_coverage_pdbs)}")

    # fetch pdb ids
    count = 0
    for pdb_id in pdb_ids_list[0:5]:
        fetch_pdb_id(pdb_id)
        print(f"{count + 1}. Downloaded {pdb_id}")
        count += 1

if __name__ == "__main__":
    main()
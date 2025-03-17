import Bio.PDB
import pandas as pd
import time
import os

def fetch_pdb_id(pdb_id, chain_id = None, rename = None, start = None, end = None):
    '''Fetches the 3D structure of a protein from the PDB and slicing it for a specific part.'''
    # fetch pdb file
    pdbl = Bio.PDB.PDBList()
    pdb_file_path = pdbl.retrieve_pdb_file(pdb_id, pdir="pdb_files", file_format="mmCif")

    if start and end and chain_id:
        sliced_file_path = slice_pdb_file(pdb_file_path, chain_id, start, end)
        os.remove(pdb_file_path)
        pdb_file_path = sliced_file_path
    if rename:
        os.rename(pdb_file_path, f"pdb_files/{rename}.cif")
        pdb_file_path = f"pdb_files/{rename}.cif"

    time.sleep(1)
    return pdb_file_path

def slice_pdb_file(pdb_filepath, chain_id, start, end):
    '''Slices, modifying a mmCif file to only represent a specific region.'''
    # parse pdb file
    parser = Bio.PDB.MMCIFParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_filepath)
    model = structure[0]
    chain = model[chain_id]

    # remove non-target residues
    for residue in chain:
        if residue.id[1] < start or residue.id[1] > end:
            chain.detach_child(residue.id)
    # remove non-target chains
    for chain in model:
        if chain.id != chain_id:
            model.detach_child(chain.id)

    io = Bio.PDB.MMCIFIO()
    io.set_structure(structure)
    sliced_pdb_filepath = pdb_filepath.replace(".cif", f"_{chain_id}_{start}_{end}.cif")
    io.save(sliced_pdb_filepath)
    return sliced_pdb_filepath

def main():
    threhold = 0.7

    # organize found pdb ids by coverage
    pdb_data = pd.read_csv("result_pdbs.tsv", sep="\t")
    # sort by coverage
    pdb_data = pdb_data.sort_values(by="Coverage", ascending=False)
    # filter into two dataframes
    high_coverage_pdbs = pdb_data[pdb_data["Coverage"] >= threhold]
    low_coverage_pdbs = pdb_data[pdb_data["Coverage"] < threhold]
    
    # print dfs
    print(f"Number of High Coverage PDB IDs: {len(high_coverage_pdbs)}")
    print(f"Number of Low Coverage PDB IDs: {len(low_coverage_pdbs)}")

    # output separated dfs
    high_coverage_pdbs.to_csv("high_coverage_pdbs.tsv", sep="\t", index=False)
    low_coverage_pdbs.to_csv("low_coverage_pdbs.tsv", sep="\t", index=False)

    # iterate through high coverage pdbs and fetch them
    count = 0
    for _, row in high_coverage_pdbs.iterrows():
        pdb_id = row["PDB ID"]
        uniprot_id = row["UniProt ID"]
        chain_id = row["Chain ID"]
        start = row["PDB Start"]
        end = row["PDB End"]
        print(f"{count + 1}. ({uniprot_id}) ", end = "")
        fetch_pdb_id(pdb_id, chain_id, rename = uniprot_id, start = start, end = end)
        count += 1

if __name__ == "__main__":
    main()
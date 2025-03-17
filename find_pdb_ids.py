import requests
import Bio.PDB
from Bio import Align, SeqIO
import time
from io import StringIO

def read_list(file_path):
    '''Reads a list of UniProt IDs from a file.'''
    with open(file_path, "r") as f:
        uniprot_ids = f.read().splitlines()
    return uniprot_ids

def fetch_uniprot_sequence(uniprot_id):
    uniprot_url = f"https://rest.uniprot.org/uniprotkb/accession/{uniprot_id}"
    headers = {"Accept": "application/fasta"}
    response = requests.get(uniprot_url, headers=headers)
    if response.status_code != 200:
        return None  # UniProt ID not found or error.

    uniprot_sequence = str(SeqIO.read(StringIO(response.text), "fasta").seq)
    return uniprot_sequence

def identify_best_structure(uniprot_id):
    '''Makes API call to EBI to fetch the PDB IDs for the best structures for a given UniProt ID.'''
    ebi_url = f"https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/{uniprot_id}"
    response = requests.get(ebi_url)
    if response.status_code != 200:
        # print(f"Best Structure not Found for {uniprot_id}")
        return {"uniprot_id": uniprot_id, "pdb_id": "____", "coverage": 0.0, "chain_id": "_", 
                "pdb_start": -1, "pdb_end": -1, "unp_start": -1, "unp_end": -1}
    
    # grab first PDB ID, coverage, start and end positions
    pdb_data = response.json()
    pdb_id = pdb_data[uniprot_id][0]["pdb_id"]
    coverage = pdb_data[uniprot_id][0]["coverage"]
    chain_id = pdb_data[uniprot_id][0]["chain_id"]
    start = pdb_data[uniprot_id][0]["start"]
    end = pdb_data[uniprot_id][0]["end"]
    unp_start = pdb_data[uniprot_id][0]["unp_start"]
    unp_end = pdb_data[uniprot_id][0]["unp_end"]

    time.sleep(1)
    return {"uniprot_id": uniprot_id, "pdb_id": pdb_id, "coverage": coverage, "chain_id": chain_id,
            "pdb_start": start, "pdb_end": end, "unp_start": unp_start, "unp_end": unp_end}

def main():
    uniprot_ids_list = read_list("lrr_benchmark_proteins_list.txt")
    print(f"Number of UniProt IDs: {len(uniprot_ids_list)}")

    # fetch pdb ids
    pdb_data_list = []
    for uniprot_id in uniprot_ids_list:
        print(f"{len(pdb_data_list) + 1}. ", end = "")
        pdb_data = identify_best_structure(uniprot_id)
        if pdb_data:
            pdb_data_list.append(pdb_data)
            print(pdb_data)
    print(f"Number of PDB IDs fetched: {len(pdb_data_list)}")

    # output list as tsv
    with open("result_pdbs.tsv", "w") as f:
        f.write("UniProt ID\tPDB ID\tCoverage\tChain ID\tPDB Start\tPDB End\tUniprot Start\tUniprot End\n")
        for pdb_data in pdb_data_list:
            f.write(f"{pdb_data['uniprot_id']}\t{pdb_data['pdb_id']}\t{pdb_data['coverage']}\t{pdb_data['chain_id']}\t{pdb_data['pdb_start']}\t{pdb_data['pdb_end']}\t{pdb_data['unp_start']}\t{pdb_data['unp_end']}\n")

if __name__ == "__main__":
    main()
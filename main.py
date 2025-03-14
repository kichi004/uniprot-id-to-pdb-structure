import requests
import Bio.PDB
from Bio import Align, SeqIO
import time
from io import StringIO

def interpret_tsv(tsv_data):
    pdb_codes = []
    tracked_uniprot_id = None
    tracked_pdb_ids = []
    for line in tsv_data.split("\n"):

        if line:
            line_uniprot_id = line.split("\t")[0]
            line_pdb_id = line.split("\t")[1]
        if not tracked_uniprot_id:
            tracked_uniprot_id = line_uniprot_id
            tracked_pdb_ids = [line_pdb_id]
        elif tracked_uniprot_id != line_uniprot_id:
            pdb_codes.append({"uniprot_id": tracked_uniprot_id, "pdb_ids": tracked_pdb_ids})
            tracked_uniprot_id = line_uniprot_id
            tracked_pdb_ids = [line_pdb_id]
        else:
            tracked_pdb_ids.append(line_pdb_id)
    return pdb_codes

def fetch_uniprot_sequence(uniprot_id):
    uniprot_url = f"https://rest.uniprot.org/uniprotkb/accession/{uniprot_id}"
    headers = {"Accept": "application/fasta"}
    response = requests.get(uniprot_url, headers=headers)
    if response.status_code != 200:
        return None  # UniProt ID not found or error.

    uniprot_sequence = str(SeqIO.read(StringIO(response.text), "fasta").seq)
    return uniprot_sequence

def fetch_best_structures(uniprot_id):
    '''Makes API call to EBI to fetch the PDB IDs for the best structures for a given UniProt ID.'''
    ebi_url = f"https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/{uniprot_id}"
    response = requests.get(ebi_url)
    if response.status_code != 200:
        print(f"Error fetching data for {uniprot_id}")
        return None  # UniProt ID not found or error.
    
    # grab first PDB ID, coverage, start and end positions
    pdb_data = response.json()
    pdb_id = pdb_data[uniprot_id][0]["pdb_id"]
    coverage = pdb_data[uniprot_id][0]["coverage"]
    start = pdb_data[uniprot_id][0]["start"]
    end = pdb_data[uniprot_id][0]["end"]

    time.sleep(1)
    return {"uniprot_id": uniprot_id, "pdb_id": pdb_id, "coverage": coverage, "start": start, "end": end}

def main():
    data = fetch_best_structures("O00206")
    print(data)

if __name__ == "__main__":
    main()
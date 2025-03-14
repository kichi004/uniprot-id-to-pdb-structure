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

def find_best_pdb(uniprot_id, pdb_ids):
    '''Determines the best PDB structure for a given UniProt ID. Prioritizes sequence coverage, deposition date, and resolution.'''
    
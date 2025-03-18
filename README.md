# uniprot-id-to-pdb-structure
Set of Python scripts to identify and download the best candidate PDB structure for a given UniProt ID based on sequence coverage. 

**Finding PDB IDs**

`find_pdb_ids.py` uses the EBI best structures API on a list of UniProt IDs to generate a table (organized as tsv) of the best PDB structure for a given UniProt ID based on sequence coverage. 

The outputted table contains the UniProt ID, PDB ID, Chain ID, Sequence Coverage, and the Start/End residues for both the PDB and UniProt sequence.

**Fetch PDB Structures**

`fetch_pdb_ids.py` uses the Biopython PDB module to fetch mmCif PDB models of high coverage (>70%) listed in the table generated from above. 

It then uses the `MMCIFParser` to isolate the chain and region of interest in the mmCIF file. The script separates the uncollected and collected structures into the `low_coverage_pdbs.tsv` and the `high_coverage_pdbs.tsv`. 

**Fetch AlphaFold Structures**

`fetch_alphafold.py` uses the Biopython alphafolddb module to fetch AlphaFold structures from the AlphaFold database. 

This is run on the table of proteins without high coverage PDB structures, `low_coverage_pdbs.tsv`. For proteins without an entry in the database, the script returns a list of failed AF fetches as well.

**Run FoldSeek Search**

1. `foldseek easy-search 2bnh.cif foldseek_db/ aln temp_dir/`


import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO, pairwise2
import requests
from sys import stderr
import pandas as pd
from typing import List, Tuple, Dict
import config

def _load_reference_sequence_records() -> Dict[str, SeqRecord]:
    reference_sequence_records = {}
    with open(config.sgd_fasta_path) as f:
        for record in SeqIO.parse(f, "fasta"):
            reference_sequence_records[record.id] = record
    return reference_sequence_records

def _load_uniprot_sequence_records() -> Dict[str, SeqRecord]:
    reference_sequence_records = {}
    with open(config.uniprot_fasta_path) as f:
        for record in SeqIO.parse(f, "fasta"):
            parts = record.id.split("|")
            assert len(parts) == 3
            reference_sequence_records[parts[1]] = record
    return reference_sequence_records


preload_reference_sequence_records = _load_reference_sequence_records()
preload_uniprot_sequence_records = _load_uniprot_sequence_records()


def _get_reference_sequence_and_description(
        systematic_gene_name: str, uniprot_id: str
) -> Tuple[Seq, str]:
    
    seq = None
    description = None
    if uniprot_id in preload_uniprot_sequence_records:
        seq_record = preload_uniprot_sequence_records[uniprot_id]
        seq = str(seq_record.seq)
        description = seq_record.description
        
    if systematic_gene_name in preload_reference_sequence_records:
        seq_record = preload_reference_sequence_records[systematic_gene_name]
        if seq is None:
            print(f"Using SGD sequence for protein {uniprot_id}", file=stderr)
            seq = str(seq_record.seq).rstrip('*')
        description = seq_record.description

    if seq is None or description is None:
        print(
            f'Did not find reference sequence for protein {systematic_gene_name}. '
            f'Protein will not be inserted!', file=stderr
        )
    return seq, description


def _get_pdb_filename(uniprot_id: str) -> str:
    return f'{config.pdb_file_prefix}{uniprot_id}.pdb'

def _download_pdb(uniprot_id: str) -> None:
    pdb_file_path = _get_pdb_filename(uniprot_id)
    if not os.path.exists(pdb_file_path):
        url = f'https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v2.pdb'
        r = requests.get(url)
        if not r.status_code == 200:
            print(
                f'Failed to get 3D structure of protein {uniprot_id}. Request status code {r.status_code}. '
                f'Inserting protein anyway.', file=stderr
            )
            return

        with open(pdb_file_path, 'wb') as f:
            f.write(r.content)


def _get_pdb_sequence(uniprot_id: str) -> Seq:
    pdb_file_path = _get_pdb_filename(uniprot_id)
    if not os.path.exists(pdb_file_path):
        return

    records = list(SeqIO.parse(pdb_file_path, 'pdb-seqres'))
    if len(records) != 1:
        print(f'More than one record for protein {uniprot_id}. Ids of records:', file=stderr)
        for record in records:
            print(f'  {record.id}', file=stderr)
        return

    record = records[0]
    return record.seq


def _get_mapping(reference_sequence: Seq, pdb_sequence: Seq) -> List[int]:
    #print("\nPAIR", str(reference_sequence), str(pdb_sequence))
    alignments = pairwise2.align.globalxx(reference_sequence, pdb_sequence)
    alignment = alignments[0]

    mapping = []
    i = 0
    for reference_element, pdb_element in zip(alignment[0], alignment[1]):
        if reference_element == '-':
            if pdb_element != '-':
                i += 1
        else:
            if pdb_element == '-':
                mapping.append(-1)
            else:
                mapping.append(i)
                i += 1
    return mapping


def get_protein_info(systematic_gene_name: str, uniprot_id: str) -> Tuple[str, List[int], str]:
    reference_sequence, description = _get_reference_sequence_and_description(systematic_gene_name, uniprot_id)
    _download_pdb(uniprot_id)
    pdb_sequence = _get_pdb_sequence(uniprot_id)

    mapping = None
    if reference_sequence is not None and pdb_sequence is not None:
        mapping = _get_mapping(reference_sequence, pdb_sequence)

    reference_sequence = None if reference_sequence is None else str(reference_sequence)
    return reference_sequence, mapping, description


def _load_map_systematic_gene_name_to_uniprot_id() -> dict:
    df = pd.read_csv(config.sgd_gene_table_path, sep='\t', skiprows=range(7), header=None, names=['col{}'.format(x) for x in range(17)],
                     usecols=[10, 16])

    map_sys_gene_name_to_uniprot_id = {}
    for i in range(len(df)):
        try:
            gene_name = df.iloc[i, 0].partition('|')[0]
            uniprot_id = df.iloc[i, 1][10:]
            map_sys_gene_name_to_uniprot_id[gene_name] = uniprot_id
        except TypeError:
            # happens if systematic_gene_name is 'nan' (not filled) and is interpreted as float (float can't be sliced)
            pass

    return map_sys_gene_name_to_uniprot_id


preloaded_map_systematic__name_to_uniprot_id = _load_map_systematic_gene_name_to_uniprot_id()


def has_valid_systematic_gene_name_and_uniprot_id(sys_gene_name: str, uniprot_id: str) -> bool:
    if sys_gene_name not in preloaded_map_systematic__name_to_uniprot_id:
        print(
            f'{sys_gene_name} is not in reference table. Can not check pair {sys_gene_name}, {uniprot_id}. '
            f'Inserting protein anyway.', file=stderr
        )
        return True
    if preloaded_map_systematic__name_to_uniprot_id[sys_gene_name] != uniprot_id:
        print(
            f'Our pair ({sys_gene_name},{uniprot_id}) does not match the expected pair:'
            f'({sys_gene_name},{preloaded_map_systematic__name_to_uniprot_id[sys_gene_name]}). '
            f'Protein will not be inserted!', file=stderr
        )
        return False
    return True

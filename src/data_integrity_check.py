"""Several functions for checking data from the Excel file"""

from sys import stderr

def compare_excel_sequence_length_and_reference_sequence_length(
    uniprot_id: str,
    excel_length: int,
    reference_sequence: str
) -> None:
    if len(reference_sequence) != excel_length:
        print(
            f'Protein {uniprot_id} has different sequence length in excel {excel_length} '
            f'than length of reference sequence {len(reference_sequence)}. '
            f'Inserting protein anyway.', file=stderr
        )


def _get_reference_sequence(uniprot_id: str, db_connection) -> str:
    cursor = db_connection.cursor()
    cursor.execute(f"SELECT protein_sequence FROM mtmod_proteins WHERE uniprot_id = '{uniprot_id}'")
    result = cursor.fetchone()

    assert cursor.fetchone() is None
    return None if result is None else result[0]


def has_modifications_on_correct_aminoacids(
    uniprot_id: str,
    allowed_amino_acids: str,
    position: int,
    modification_type : str,
    db_connection    
) -> bool:
    sequence = _get_reference_sequence(uniprot_id, db_connection)

    if sequence is None:
        print(
            f"'Position {position} of modification {modification_type} in {uniprot_id} can not be verified, protein not in db. "
            f'Might be skipped due to an error in previous phase. Modification will not be inserted!.', file=stderr
        )
        return False

    if position > len(sequence):
        print(
            f'Position of modification {position} is outside of protein {uniprot_id} with sequence length {len(sequence)}. '
            f'Modification will not be inserted!', file=stderr
        )
        return False

    amino_acid = sequence[position - 1]
    if allowed_amino_acids != "*" and amino_acid not in allowed_amino_acids:
        print(
            f'Protein {uniprot_id} has modification {modification_type} '
            f'on amino acid {amino_acid}, position {position}, '
            f'allowed {allowed_amino_acids}.', file=stderr
        )
    return True

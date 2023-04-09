"""The main script for parsing the database in Excel format and storing it in SQLite database which needs to be created beforehand."""

import sqlite3
import json
import pandas as pd
import re
import config
from sys import stderr, stdout

from tqdm import tqdm
from reference_db import get_protein_info, \
    has_valid_systematic_gene_name_and_uniprot_id
from data_integrity_check import has_modifications_on_correct_aminoacids,   \
    compare_excel_sequence_length_and_reference_sequence_length


def fill_proteins(db_connection) -> None:
    cursor = db_connection.cursor()

    proteome_sheet = pd.read_excel(
        config.excel_path,
        sheet_name='24_Reference mt proteome'
    )
    names_dict = {
        'Systematic gene name' : 'sys_gene_name',
        'Standard gene name' : 'std_gene_name',
        'Uniprot_ID' : 'uniprot_id',
        'Protein names': 'prot_name',
        'Gene names' : 'gene_names',
        'Length': 'protein_length'
    }
    proteome_sheet.rename(columns=names_dict, inplace=True)

    
    sql_query = """
        INSERT INTO mtmod_proteins
          (uniprot_id, systematic_gene_name, standard_gene_name, 
           protein_name, gene_names, description, 
           protein_sequence, mapping)
          VALUES (?, ?, ?, ?, ?, ?, ?, ?)
    """

    for row in tqdm(proteome_sheet.itertuples(),
                    total=len(proteome_sheet), desc='Populating proteins',
                    file=stdout):

        if not has_valid_systematic_gene_name_and_uniprot_id(row.sys_gene_name, row.uniprot_id):
            continue

        protein_sequence, mapping, description = get_protein_info(row.sys_gene_name, row.uniprot_id)

        if protein_sequence is None:
            continue

        compare_excel_sequence_length_and_reference_sequence_length(row.uniprot_id, row.protein_length, protein_sequence)

        # make sure that std gene name is the first in gene name list
        name_list = row.gene_names.split()
        if len(name_list) == 0 or name_list[0] != row.std_gene_name:
            old_name_list = name_list.copy()
            if row.std_gene_name in name_list:
                name_list.remove(row.std_gene_name)
            name_list.insert(0, row.std_gene_name)
            if old_name_list != name_list:
                print(f"Gene name list changed from {old_name_list} to {name_list} uniprot_id={row.uniprot_id}", file=stderr)
        gene_names = " ".join(name_list)
        
        if mapping is not None:
            mapping = json.dumps(mapping)
        try:
            cursor.execute(sql_query,
                           (row.uniprot_id, row.sys_gene_name,
                            row.std_gene_name, row.prot_name,
                            gene_names, description, protein_sequence, mapping))
        except Exception as e:
            print(e, f"uniprot_id={row.uniprot_id}, sys_gene_name={row.sys_gene_name}", file=stderr)

            
    db_connection.commit()


def fill_sources(db_connection) -> None:
    cursor = db_connection.cursor()

    source_sheet = pd.read_excel(
        config.excel_path,
        sheet_name='1_References'
    )

    names_dict = {
    'Short reference:' : 'source_id',
    'Full reference' : 'source_description', 
    'DOI' : 'source_url',
    'Annotation' : 'annotation'    
    }
    #print(",".join(source_sheet.columns) + "*")
    source_sheet.rename(columns=names_dict, inplace=True)
    #print(",".join(source_sheet.columns) + "*")
    
    sql_query = """
        INSERT INTO mtmod_source
          (source_id, source_description, source_url, annotation, position)
          VALUES (?, ?, ?, ?, ?)
    """

    position = 0;
    for row in tqdm(source_sheet.itertuples(), total=len(source_sheet), desc='Populating sources', file=stdout):

        # skip rows without ID (empty etc)
        if pd.isna(row.source_id):
            continue

        try:
            cursor.execute(sql_query, (str(row.source_id).strip(),
                                       str(row.source_description).strip(),
                                       str(row.source_url).strip(),
                                       str(row.annotation).strip(),
                                       position))
            position += 1
        except Exception as e:
            print(e, f"{row.source_id}", file=stderr)

    db_connection.commit()

def fill_modifications(db_connection) -> None:
    # read table of modifications
    modifications_df = pd.read_csv(config.modifications_csv_path)
    for row in modifications_df.itertuples():
        if row.Code == "multiple":
            continue
        print("Starting modification", row.Sheet, row.Code, row.Amino_acids)
        fill_one_modification(db_connection, row.Sheet, row.Code, row.Amino_acids)

def fill_one_modification(db_connection, sheet_name, modification_type, allowed_amino_acids) -> None:
    cursor = db_connection.cursor()

    proteome_sheet = pd.read_excel(
        config.excel_path,
        sheet_name=sheet_name
    )

    sql_query_modifications = """
        INSERT INTO mtmod_modifications 
          (uniprot_id, position, modification_type)
          VALUES (?, ?, ?)
    """

    sql_query_modifications_source = """
        INSERT INTO mtmod_modification_source
          (modification_id, source_id)
          VALUES (?, ?)
    """

    all_columns = list(proteome_sheet.columns)
    uniprot_column = 2
    site_column = 3
    last_fixed_column = 4
    if all_columns[uniprot_column] != "Uniprot_ID":
        raise ValueError(f"Bad column names {all_columns} (column {uniprot_column} should be the Uniprot ID")
    if "site" not in all_columns[site_column]:
        raise ValueError(f"Bad column names {all_columns} (column {site_column} should be the site number")
    if not all_columns[last_fixed_column].startswith("Annotation"):
        raise ValueError(f"Bad column names {all_columns} (column {last_fixed_column+1} should start with Annotation)")

    # mapping from column number to source id
    source_dict = {}
    for (index, column_name) in list(enumerate(all_columns))[last_fixed_column + 1:]:
        if not str(column_name).startswith("Unnamed"):
            source_dict[index] = str(column_name).strip()

    for row in tqdm(proteome_sheet.itertuples(name=None), total=len(proteome_sheet), desc='Populating modifications', file=stdout):

        
        # skip one value in the row, which is the index of the row
        uniprot_id = str(row[uniprot_column+1]).strip()

        # skip empty positions
        position = row[site_column+1]
        if pd.isna(position):
            print(f"Position missing in row {row[0]} {(uniprot_id, modification_type)}", file=stderr)
            continue
        modification = (uniprot_id, position, modification_type)

        if isinstance(position, str) and re.match(r'[a-zA-Z][.-]', position):
            position = int(position[2:])
        if not isinstance(position, int):
            print(f"Position {position} not int {modification}", file=stderr)
            position = int(position)

        if not has_modifications_on_correct_aminoacids(uniprot_id, allowed_amino_acids, position, modification_type, db_connection):
            continue

        try:
            cursor.execute(sql_query_modifications, (uniprot_id, position, modification_type))
            last_id = cursor.lastrowid
        except Exception as e:
            print(f' did not add modification {modification}.', file=stderr)
            print(e, file=stderr)
            continue

        for (column_number, column_name) in source_dict.items():
            (use_this_ref, ref_id) = process_ref(row[column_number + 1], column_name, modification, cursor)
        
            if use_this_ref:
                try:
                    cursor.execute(sql_query_modifications_source, (last_id, ref_id))
                except Exception as e:
                    print(f' did not add modification source pair (modification, ref_id)', file=stderr)
                    print(e, file=stderr)


    db_connection.commit()

def process_ref(value, column_name, modification, cursor):
    use_this_ref = False
    value = str(value).strip()
    ref_id = column_name
    if value.startswith('YES'):
        # YES means take name of column as reference for this modification
        use_this_ref = True
    elif value != "nan" and value != '':
        # other noneempty value is typically error, unless column name if Other,
        #   then take value in excel cell as reference for this modification
        if column_name == "Other":
            use_this_ref = True
            ref_id = value
        else: # print warning
            print(f"Bad string in reference list {value} for {modification}", file=stderr)

    if use_this_ref:
        # check that ref is valid
        sql_query = "SELECT source_id FROM mtmod_source WHERE source_id=?"
        cursor.execute(sql_query, (ref_id,))
        results = cursor.fetchall()
        if not results:
            print(f'did not find {ref_id} from {modification} among sources', file=stderr)
            
    return (use_this_ref, ref_id)
    

def main() -> None:
    # connect to a database
    db_connection = sqlite3.connect(config.database_path)
    
    fill_proteins(db_connection)
    fill_sources(db_connection)
    fill_modifications(db_connection)

    # close db connection
    db_connection.close()
    

if __name__ == "__main__":
    main()

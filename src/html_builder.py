"""Script generating website for individual proteins as well as the main page.
It uses the SQLite database of proteins and modifications."""

import subprocess
from tqdm import tqdm
import os
import sys
from typing import List, Tuple
import json
import sqlite3
import shutil
import argparse
import inspect

import jinja2
import markupsafe
import pandas as pd

import config


def get_all_uniprot_ids() -> List[str]:
        cnx = create_context()
        cursor = cnx.cursor()
        cursor.execute("SELECT DISTINCT uniprot_id FROM mtmod_proteins")
        uniprot_ids = [row[0] for row in cursor]
        return uniprot_ids

def get_protein_info():
        cnx = create_context()
        cnx.row_factory = sqlite3.Row
        cursor = cnx.cursor()
        cursor.execute(
                "SELECT uniprot_id, protein_name, systematic_gene_name, standard_gene_name, gene_names, description, protein_sequence, mapping "
                "FROM mtmod_proteins ORDER BY uniprot_id"
        )

        # get all rows to a dictionary with uniprot_id as key
        result = {row['uniprot_id']:dict(row) for row in cursor}
        
        # fix mapping, compute protein length
        for id, row in result.items():
                if row['mapping'] is not None:
                        row['mapping'] = json.loads(row['mapping'])
                row['length'] = len(row['protein_sequence'])
                
        return result

def get_modification_counts(protein_info, modification_list):
        """Add modification counts to protein_info"""
        
        cnx = create_context()
        cursor = cnx.cursor()
        cursor.execute(
                "SELECT COUNT(*),uniprot_id, modification_type FROM mtmod_modifications GROUP BY uniprot_id, modification_type"
                )

        count_dict = {(row[1], row[2]):row[0] for row in cursor}
        for uniprot_id, row in protein_info.items():
                for mod_record in modification_list:
                        mod_type =  mod_record["Code"]
                        key = (uniprot_id, mod_type)
                        if key in count_dict:
                                row[mod_type] = count_dict[key]
                        else:
                                row[mod_type] = 0
        

def get_modifications_for_protein(protein_id) -> List[Tuple[int, str, int]]:
        cnx = create_context()
        cursor = cnx.cursor()
        cursor.execute("SELECT position, modification_type, modification_id FROM mtmod_modifications "
                       "WHERE uniprot_id = ? ORDER BY position", (protein_id,)
        )
        uniprot_modifications = cursor.fetchall()
        return uniprot_modifications


def create_context():
        # connect to a database
        return sqlite3.connect(config.database_path)


def get_jinja_template(filename):
        environment = jinja2.Environment(
                loader=jinja2.FileSystemLoader("templates/"),
                autoescape=jinja2.select_autoescape()
        )
        return environment.get_template(filename)
        

def get_sources_for_modification_id(
        modification_id
) -> List[Tuple[str, str]]:
        cnx = create_context()
        cursor = cnx.cursor()
        cursor.execute(
                "SELECT source_id FROM mtmod_modification_source WHERE modification_id = ?", (modification_id,)
        )
        source_ids = cursor.fetchall()
        sources = []
        for source_id in source_ids:
                cursor.execute(
                        "SELECT source_description, source_url FROM mtmod_source WHERE source_id = ?", (source_id[0],)
                )
                sources.extend(cursor.fetchall())
        return sources



def create_protein_page(protein_info, verbose, modification_df):
        uniprot_id = protein_info["uniprot_id"]
        modifications = get_modifications_for_protein(uniprot_id)
        
        mapping = protein_info["mapping"]
                
        sources = []
        positions2d = []
        positions3d = []
        multiple = {}

        # 1-indexed to 0-indexed
        positions = [mod[0] - 1 for mod in modifications]
        for pos, pos_type, mod_id in modifications:
                pos -= 1
                mod_sources = get_sources_for_modification_id(mod_id)
                for source_name, source_url in mod_sources:
                        sources.append({'name':source_name, 'url':source_url,
                                        'pos':pos, 'type':pos_type})
                if positions.count(pos) > 1:
                        if pos in multiple:
                           multiple[pos].append(pos_type)
                        else:
                           multiple[pos] = [pos_type]
                        pos_type = "multiple"
                positions2d.append((pos, pos_type))
                if mapping is not None and mapping[pos] != -1:
                        positions3d.append((mapping[pos], pos_type))

        # sort and add to each source its order within the same position
        sources.sort(key=lambda x: x["pos"])
        last = -1
        count = 0
        for source in sources:
                if source["pos"] != last:
                        last = source["pos"]
                        count = 0;
                count += 1
                source['count'] = count
                        

        different_modifications = []
        single_chars = list(protein_info['protein_sequence'])
        for pos, pos_type in positions2d:
                modification_type = pos_type
                if pos in multiple:
                    modification_type = " ".join(multiple[pos])
                single_chars[pos] = (
                        # f'<div class="tooltip"><span style="color: {modification_df.loc[pos_type,"Text_color"]}; '
                        # f'background-color:{modification_df.loc[pos_type,"Color"]}">{single_chars[pos]}</span>'
                        # f'<span class="tooltiptext">{pos + 1}, {pos_type}</span></div>'

                         f'<a href="#source{pos + 1}_1" data-bs-toggle="tooltip" data-bs-placement="top" title="{pos + 1}, {modification_type}">'
                         f'<span style="color: {modification_df.loc[pos_type,"Text_color"]}; '
                         f'background-color:{modification_df.loc[pos_type,"Color"]}">'
                         f'{single_chars[pos]}</span></a>'

                )
                if pos_type not in different_modifications:
                        different_modifications.append(pos_type)
        
        hasStructure = True
        if not os.path.exists(f"{config.pdb_file_prefix}/{uniprot_id}.pdb"):
                print(
                        "Structure with uniprot_id",
                        uniprot_id,
                        "not in 3d_structures folder. Will generate page without 3d structure.", file=sys.stderr
                )
                hasStructure = False
                template = get_jinja_template("protein_page.html")
                with open(f"{config.web_output_dir}/{uniprot_id}.html", "w") as text_file:
                        print(template.render(protein_info = protein_info,
                                single_chars = single_chars,
                                sources = sources,
                                hasStructure = hasStructure,
                                modifications = different_modifications,
                                modification_df = modification_df ), file = text_file)
                return

        
        script = f"pymol -c pymol_script.pml -- ../{config.pdb_file_prefix}"
        script += " {}"
        for _ in positions3d:
                script += " {}"

        pos3d_str = []
        for pos, pos_type in positions3d:
                pos3d_str.append(f'{pos+1}."{modification_df.loc[pos_type,"Color"]}"')

        script = script.format(uniprot_id, *pos3d_str)

        if verbose:
                print("Printing the generated pymol script...")
                print(script)
        else:
                script += " > /dev/null"

        subprocess.call(script, shell=True, cwd="pdb_to_html")

        with open("./pdb_to_html/representation.txt", 'r') as file:
                representation = file.read()

        with open("./pdb_to_html/pdb_file.txt", 'r') as file:
                pdb_file = file.read()

        template = get_jinja_template("protein_page.html")
        with open(f"{config.web_output_dir}/{uniprot_id}.html", "w") as text_file:
                print(template.render(protein_info = protein_info,
                                single_chars = single_chars,
                                sources = sources,
                                hasStructure = hasStructure,
                                pdb_file = pdb_file,
                                representation = representation,
                                modifications = different_modifications,
                                modification_df = modification_df ), file = text_file)


def get_modification_list(modification_df):
        modification_list = []
        for (i, row) in enumerate(modification_df.itertuples()):
                if row.Index != "multiple":
                        modification_list.append({'Column':i+5, 'Code':row.Index, 'Full_name':row.Full_name, 'Tiny':row.Tiny})

        return modification_list

                

def main(verbose=False, debug=False):
        """build index.html, database.html and protein pages
        
        verbose: print scripts fr indovidual proteins
        debug: consider only proteins P31380 P00360 P18963 A5Z2X5
        """
        
        print('Building protein browser...')


        uniprot_ids = get_all_uniprot_ids()
        if len(uniprot_ids) == 0:
                print("Error: no uniprot ids loaded from database!", file=sys.stderr)
                exit(1)

        modification_df = pd.read_csv(config.modifications_csv_path)
        modification_df.set_index('Code', inplace=True)
                
        modification_list = get_modification_list(modification_df)
        protein_info = get_protein_info()
        get_modification_counts(protein_info, modification_list)
                
        if not os.path.exists(config.pdb_file_prefix):
                print(f"Error: could not find folder {config.pdb_file_prefix}!", file=sys.stderr)
                exit(1)

        if not os.path.exists("pdb_to_html"):
                print("Error: could not find folder pdb_to_html containing script!", file=sys.stderr)
                exit(1)
                
        if not os.path.exists(config.web_output_dir):
                os.makedirs(config.web_output_dir)


        protein_list = [protein_info[uniprot_id] for uniprot_id in uniprot_ids]

        template = get_jinja_template("index.html")
        with open(f"./{config.web_output_dir}/index.html", "w") as text_file:
                print(template.render(protein_list=protein_list), file=text_file)

        template = get_jinja_template("database.html")
        with open(f"./{config.web_output_dir}/database.html", "w") as text_file:
                print(template.render(protein_list=protein_list, modification_list = modification_list), file=text_file)

        shutil.copytree("web_include", f"{config.web_output_dir}/include", dirs_exist_ok=True)
                
        if debug: # only several proteins
                pages_todo = ["P31380", "P00360", "P18963", "A5Z2X5"]
        else:
                # not debug - do all proteins
                pages_todo = uniprot_ids
                
        for uniprot_id in tqdm(pages_todo,  desc='Uniprot IDs', file=sys.stdout):
                create_protein_page(protein_info[uniprot_id], verbose, modification_df)



if __name__ == "__main__":
        parser = argparse.ArgumentParser(description=inspect.getdoc(main))
        parser.add_argument(
                "-v", dest="verbose",  action='store_true')
        parser.add_argument(
                "-d", dest="debug",  action='store_true')
        args = parser.parse_args()
        main(** vars(args))

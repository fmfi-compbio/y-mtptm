# The source code for building the Yeast Mitochondrial PTM Database at http://compbio.fmph.uniba.sk/y-mtptm/

This repository contains scripts necessary to build the y-mtPTM database located at http://compbio.fmph.uniba.sk/y-mtptm/
Assuch, it is not really meant to be run by others, but it may help in setting up similar projects.
Once generated, the whole website consists only of HTML, CSS and JavaScript files and does not require server-side computation.

## How to run

Prerequisites:
* The software was tested only on linux systems
* Python version at least 3.8 is required

Also install 
* Pymol (e.g. via [conda](https://pymol.org/conda/) or as a ubuntu package)
* [SQLite](https://sqlite.org/index.html)
* Data files are downloaded via `wget` but alternatives such as `curl` can also be used

Downloading repository and setting up Python libraries:
```bash
git clone  git@github.com:fmfi-compbio/y-mtptm.git

# change to src dir
cd y-mtptm/src

# setup virtual environment
python3 -m venv env
# activate the environment
source env/bin/activate
# install required libraries to the environment
python3 -m pip install -r requirements.txt


# all subsequent commands should be run in src folder
# and with environment activated
```

Downloading necessary data files:
```bash
# S. cerevisiae proteome fasta from Uniprot
cd ../data/uniprot
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000002311/UP000002311_559292.fasta.gz
gunzip *.gz

# proteome fasta and gene names for SGD 
cd ../sgd
wget http://sgd-archive.yeastgenome.org/sequence/S288C_reference/orf_protein/archive/orf_trans_all_R64-3-1_20210421.fasta.gz
wget http://sgd-archive.yeastgenome.org/curation/literature/archive/gene_association.sgd.20210510.gaf.gz
gunzip *.gz

cd ../../src
```

Building the SQLite database from the Excel file:
* This will rewrite database provided in the repository but also download necessary pdb files
```bash

# remove the original database
rm ../data/excel/ymtptm.db
# create database schema
sqlite3 ../data/excel/ymtptm.db < create_db.sql
# fill database with data, download pdb files
# this takes some time
python3 excel_parser.py  2> excel.err > excel.log
# pdb files will be in ../data/pdb
# excel.err will contain some warnings
```

Building website from the SQLite database
```bash
python3 html_builder.py 
```

This final step creates html files in `../web`; these files can be then viewed in a browser locally or placed on a webserver.


## The main files

* `data/excel/*.xlsx` Human readable Excel file with PTM data
* `data/excel/ymtptm.db` data in SQLite format
* `data/excel/modifications.csv` configuration file with all considered PTM types
* `data/pdb` folder for pdb files downloaded by excel_parser.py
* `data/sgd` folder for files from the SGD database
* `data/uniprot` folder for files from the Uniprot database
* `web` folder for the resulting website
* `src/config.py` filenames of data files
* `src/create_db.sql` SQLite database schema
* `src/excel_parser.py` script for converting database from Excel to SQLite
* `src/html_builder.py` script for building website from SQLite database
* `src/templates` HTML templates for jinja library
* `src/web_include` images and CSS files used on the website directly
* `src/pdb_to_html` files needed to convert pdb files for GLmol via pymol


The following files contain code originating from the [GLmol package](https://github.com/biochem-fan/GLmol):

* `src/pdb_to_html/pymol2glmol.py`
* `src/templates/protein_page.html`
* `src/web_include/glmol.js`

## Authors and license

The database content, source code and design of the website were created by Bronislava Brejová, Veronika Vozáriková, Ivan Agarský, Hana Derková, Matej Fedor, Dominika Harmanová, Lukáš Kiss, Andrej Korman, Martin Pašen, Filip Brázdovič, Jozef Nosek, Tomáš Vinař, Ľubomír Tomáška.

All of the authors were affiliated at the Faculty of Natural Sciences or the Faculty of Mathematics, Physics and Informatics, Comenius University in Bratislava.

Several files listed above contain modified code from the [GLmol package](https://github.com/biochem-fan/GLmol) by biochem_fan, which is distributed under the GNU Lesser General Public License. Therefore the whole repository y-mtptm is also distributed under this license. 
CREATE TABLE `mtmod_proteins` (
  `uniprot_id` varchar(64) NOT NULL,
  `systematic_gene_name` text NOT NULL,
  `standard_gene_name` text NOT NULL,  
  `protein_name` text NOT NULL,
  `gene_names` text NOT NULL,
  `description` text NOT NULL,
  `protein_sequence` text NOT NULL,
  `mapping` text,
  PRIMARY KEY (`uniprot_id`), UNIQUE (`systematic_gene_name`) 
);


CREATE TABLE `mtmod_source` (
  `source_id` varchar(40) NOT NULL,
  `source_description` text NOT NULL,
  `source_url` text,
  PRIMARY KEY (`source_id`)
);

CREATE TABLE `mtmod_modifications` (
  `modification_id` INTEGER PRIMARY KEY AUTOINCREMENT,
  `uniprot_id` varchar(64) NOT NULL,
  `position` int NOT NULL,
  `modification_type` text NOT NULL,
  CONSTRAINT `fk_uniprot_id` FOREIGN KEY (`uniprot_id`)
    REFERENCES mtmod_proteins(`uniprot_id`)
    ON UPDATE CASCADE ON DELETE CASCADE
  /* PRIMARY KEY (`modification_id`)  */
);


CREATE TABLE `mtmod_modification_source` (
  `modification_id` int NOT NULL,
  `source_id` varchar(40) NOT NULL,
  CONSTRAINT `fk_modification_id` FOREIGN KEY (`modification_id`)
    REFERENCES mtmod_modifications(`modification_id`)
    ON UPDATE CASCADE ON DELETE CASCADE,
  CONSTRAINT `fk_source_id` FOREIGN KEY (`source_id`)
    REFERENCES mtmod_source(`source_id`)
    ON UPDATE CASCADE ON DELETE CASCADE
); 


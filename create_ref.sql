-- Create reference database for lipid network

-- Show databases and tables
.databases
.tables

-- open an existing database
.open HiVE.db

-- Drop database
drop database HiVE;

-- Add master list table
create table master(
  Source TEXT,
  Species TEXT,
  Organ TEXT,
  Matrix TEXT,
  Tx TEXT,
  standard_name TEXT,
  Name TEXT,
  Major_Class TEXT,
  Class TEXT,
  Abbreviation TEXT,
  Saturation TEXT,
  Highest_Saturation_Degree TEXT,
  Cluster TEXT,
  Subclass TEXT,
  Chain_Length TEXT,
  Chain_Type TEXT,
  Subclass_2 TEXT,
  synthesis_pathway TEXT,
  label TEXT,
  label_saturation TEXT
);

-- Add lipid network nodes
create table nodes(
  Node TEXT,
  Label TEXT,
  ID TEXT,
  ID_string TEXT,
  Abbreviation TEXT,
  Representative_name TEXT,
  Major_class TEXT,
  Synthesis_pathway TEXT,
  DOI TEXT,
  SMILES TEXT,
  LipidMAP TEXT
);

-- Add lipid network edges
create table edges(
  'from' TEXT,
  'to' TEXT,
  Level TEXT,
  Major_class TEXT,
  Synthesis_pathway TEXT,
  Leading_name TEXT,
  Enzyme_id TEXT,
  Ensembl TEXT,
  hgnc_symbol TEXT,
  PMCID TEXT,
  URL TEXT,
  DOI TEXT
);

-- Add lipid class edges
create table lipid_edges_pw(
  'x' TEXT,
  'y' TEXT,
  synthesis_pathway TEXT
);

-- Add lipid class ratio comparisons
create table lipid_class_ratios(
  'class1' TEXT,
  'class2' TEXT,
  comp_type TEXT,
  label TEXT
);

-- import .txt file containing annotations into respective tables
.mode csv
.separator ","
.import --csv --skip 1 ref/archive/0_masterlist.csv master

-- nodes
.import --csv --skip 1 nodes.csv nodes

-- edges
.import --csv --skip 1 edges.csv edges

-- pathway edges
.import --csv --skip 1 archive/0_lipidnet_pw_edges.csv lipid_edges_pw

-- class ratios
.import --csv --skip 1 lipid_comps.csv lipid_class_ratios

-- add rows to an existing table

-- Remove table
drop table lipid_nodes;

-- Export the results of a query (run line-by-line)
-- .show shows changes to the database parameters
.show
.headers on
.mode csv
.output test.csv
select * from lipid_nodes;
.output stdout

-- quit SQLite
.quit
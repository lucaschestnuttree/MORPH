# Overview of the Repository

```
├── README.md          <- The top-level README for developers using this project.
├── data
│   ├── 1-raw            <- The original, immutable data dump.
│   └── 2-processed      <- exploitable databases and predictions
│
├── requirements
│   ├── requirements-operator.txt (conda environment for general python): Python 3.10.18
│   ├── requirements-hail.txt (conda environment to run hail): Python 3.10.18
│   ├── requirements-pymut.txt (conda environment to run pymut): Python 3.10.19
│   ├── environment-operator.yml
│   ├── environment-hail.yml
│   └── environment-pymut.yml
│                       
├── src                <- Source code for use in this project.
    ├── 1-data_generation Data downloads, cleaning, harmonization, filtering operations
    └── 2-mapping   <- scripts to map & visualize data
```

# Script Description

## 1-data_generation
- 1-fetch_data.py (fetches data from AlphaMissense).
- 2-parse_AM_data.py (parses alphamissense .tsv(s) by Uniprot Id).
- 3-gnom_1_file.py (fetches gnomAD data using hail starting from alphamissense hg38 data for given Uniprot Id).
- 4-missense_gnomXmiss_1_file.py (merges alphamissense hg38 and gnomAD data & uses alphamissense aa_substitutions and ensemble-VEP to annotate entries).
- 5--tar_amissenseXgnome tar and gz the data in batches and make a lookup table.

## 2-mapping
- map_2_structure: scripts to map values to structures b-factors/additional chains.
    - all_aa_AM_2_structure.py: maps alphamissense pathogenicity.
    - pymut_gnomad_af_2_structure.py maps gnomad allele frequency.
    - pymut_gnomad_anc_2_structure maps gnomad ancestry specific allele frequency.
- pymol: scripts for visualization of map_2_structure produced pdbs in pymol.
    - color_am.py Tamas Hegedus @ hegelab.org function to visualize AM score in pymol.
    - color_gnom.py my custom functions to visualize gnomad allele frequency and ancestry specific data in pymol.
# Data description

## 1-raw
### missense_db
Available AlphaMissense data from Cheng et al https://www.science.org/doi/10.1126/science.adg7492

## 2-processed
### missense_jul
- aa_subsitutitons AlphaMissense .tsv for all aa substitutions parsed by Uniprot
- hg38 AlphaMissense .tsv for all hg38 substitutions parsed by Uniprot
- gnomAD_AF gnomAD data for corresponding ranges in hg38 matching to hg38 AlphaMissense

### missenseXgnomeAD
- merged missense_jul records keeping only missense variants

### missenseXgnomeAD_tars_gz
- merged records tar and gz in batches

### pdb
- a test pdb

### stuff
- a collection of very seriously thought out items

### z_empty 
- records where gnomAD_AF gnomAD data is empty with our workflow
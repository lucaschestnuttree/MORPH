# Advanced usage

[Back to the project overview](../README.md)

This documentation is intended to share workflows to rebuild datasets or map prepared data onto structures using your local server and/or machine.

Examples below use Bash from the repository root. Replace example paths with your own files.

Prepared data are available from the public [Zenodo record 18163874](https://zenodo.org/records/18163874).

## Environment setup

The repository contains paired Conda YAML definitions and Conda package exports for three workflow stages:

| Environment | Role                                         | Conda YAML                                                           | Conda package export                                                   |
| ----------- | -------------------------------------------- | -------------------------------------------------------------------- | ---------------------------------------------------------------------- |
| `operator`  | Download, parse, and annotate data           | [environment-operator.yml](../requirements/environment-operator.yml) | [requirements-operator.txt](../requirements/requirements-operator.txt) |
| `hail`      | Query gnomAD with Hail                       | [environment-hail.yml](../requirements/environment-hail.yml)         | [requirements-hail.txt](../requirements/requirements-hail.txt)         |
| `pymut-env` | Map scores and build alternative side chains | [environment-pymut.yml](../requirements/environment-pymut.yml)       | [requirements-pymut.txt](../requirements/requirements-pymut.txt)       |

The YAML files are portable environment definitions with package names and channels. The `.txt` files are Conda exports with pinned package versions, including Hail in `requirements-hail.txt`; they are most useful for rebuilding on a similar Linux system, not as portable pip input. Use the `pymut-env` environment for structure mapping, and clone/install [pyMUT](https://github.com/gerdos/pyMUT) before running the mapping scripts:

```bash
conda activate pymut-env
git clone https://github.com/gerdos/pyMUT.git ../pyMUT
export PYTHONPATH="$(cd ../pyMUT && pwd)"
```

For local structure viewing, install open-source PyMOL from conda-forge either in your active environment or in a dedicated environment:

```bash
# Install into the active environment
conda install conda-forge::pymol-open-source

# Or create a dedicated PyMOL environment
conda create -n pymol -c conda-forge pymol-open-source
conda activate pymol
```

## Wrapper configuration

The provided SLURM/Bash wrappers are cluster-specific. Before using them, obtain Hail and pyMUT and update the configuration files for your machine:

- `src/1-data_generation/server_wrappers/__config.sh`: `CONDA_SH`, `HT_PATH`, `LOG_DIR`, `VEP_SIF`, `CACHE_DIR`, and `CLINVAR_VCF`
- `src/2-map_2_structure/server_wrappers/__config.sh`: `CONDA_SH`

Submit wrappers with `sbatch` from the matching `server_wrappers` directory. The direct Python commands below do not use these configuration files.

## Input requirements

The gnomAD workflows expect:

- A per-protein merged TSV, normally named `<UniProt>_amissenseXgnomAD.tsv`. Output names are derived from the input TSV filename.
- The target protein in **chain A** of a PDB, with its sequence corresponding to the selected UniProt reference.
- **PDB residue numbers matching the one-based UniProt positions** used in `protein_variant`. A truncated protein can be used only when the remaining residues retain their canonical numbering. Automatic offset correction is not implemented.
- Continuous residues for pyMUT side-chain construction. Insertion codes are skipped by the mapping scripts' residue maps; normalize inputs explicitly before use.

The overall-frequency script requires `protein_variant`, `am_class`, `allele_count`, and `allele_frequency`. Category mapping also requires `uniprot_id`, and frequency columns named `<category>_allele_frequency`. It accepts `allele_frequency` as a category for the overall-frequency column. An absent `intron` column is created internally.

The notebook additionally validates that the UniProt entry is reviewed and human, that prepared data exists, and that every chain A CA residue has the matching amino acid at the same UniProt position. A truncated PDB must retain its canonical residue numbers; a fragment renumbered from 1 is rejected. A fragment whose sequence occurs more than once in UniProt is also rejected because its numbering cannot be verified from sequence alone.

The command-line gnomAD mappers check the reference amino acids at retained variant positions and fail on mismatches or no overlap. They cannot verify positions absent from the TSV, so validate the full chain against UniProt before command-line use.

## Map prepared data to a structure

### Overall allele frequency

```bash
python src/2-map_2_structure/pymut_gnomad_af_2_structure.py \
  data/2-processed/amissenseXgnomAD/P10635_amissenseXgnomAD.tsv \
  data/2-processed/pdb/P10635.pdb \
  data/2-processed/morph_outputs/af
```

Output: `data/2-processed/morph_outputs/af/P10635_AF_gnomad.pdb`.

The script ranks up to three retained variant rows by overall frequency at each position and builds variant chains `V`, `W`, and `X`. It writes the highest-ranked frequency bin to chains `A`, `V`, `W`, and `X`; the B-factor values on `W` and `X` therefore do not encode their individual variant frequencies.

### Category-specific frequency

```bash
python src/2-map_2_structure/pymut_gnomad_anc_2_structure.py \
  data/2-processed/amissenseXgnomAD/P10635_amissenseXgnomAD.tsv \
  data/2-processed/pdb/P10635.pdb \
  data/2-processed/morph_outputs/categories \
  nfe,sas,afr,fin,amr
```

Output: `data/2-processed/morph_outputs/categories/P10635_category_gnomad.pdb`.

For each category, the script selects the most frequent retained substitution at each position. The current script assigns categories in input order to chains `T`, `U`, `V`, `W`, `X`, `Y`, and `Z`, supporting up to seven categories. Read the generated PDB remarks to recover the assignments. Chain A receives the maximum frequency bin across the selected categories.

The extraction pipeline includes `nfe`, `fin`, `mid`, `remaining`, `ami`, `asj`, `afr`, `eas`, `sas`, `amr`, `XX`, and `XY` frequency columns. `XX` and `XY` are sex categories. Availability depends on the input table; an unavailable category currently produces a warning and is skipped during frequency calculation.

### Mean AlphaMissense score per residue

```bash
mkdir -p data/2-processed/morph_outputs/am
python src/2-map_2_structure/all_aa_AM_2_structure.py \
  data/2-processed/missense_jul/aa_substitutions/P10635_aa_substitutions.tsv \
  data/2-processed/pdb/P10635.pdb \
  data/2-processed/morph_outputs/am/P10635_AM.pdb
```

This mapper requires `protein_variant` and `am_pathogenicity` columns. It averages the supplied rows at each position, rounds to two decimal places, and writes that value to matching ATOM records.

It currently matches residue numbers across **all protein chains**, rather than restricting annotations to chain A. Use an input containing only the intended target protein, or explicitly prepare the chains you intend to annotate. Unscored residues retain their input B-factors, so a generic color scale can mix original values with mapped scores.

## View in PyMOL

From a PyMOL session whose working directory is the repository root:

```pymol
run src/3-pymol/color_gnom.py
load data/2-processed/morph_outputs/af/P10635_AF_gnomad.pdb, morph_af_model
morph_af sel=morph_af_model
```

For a category-mapped file:

```pymol
load data/2-processed/morph_outputs/categories/P10635_category_gnomad.pdb, morph_category_model
morph_anc sel=morph_category_model
```

To create selections for the reference residue and alternative side chains at one position:

```pymol
pos_var sel=morph_af_model, 123
```

For an AlphaMissense-mapped file:

```pymol
run src/3-pymol/color_am.py
load data/2-processed/morph_outputs/am/P10635_AM.pdb, morph_am_model
coloram sel=morph_am_model
```

The commands select and color different frequency ranges by default. `morph_af sel=morph_af_model, 1` broadens the highlighted range to include bin 4; `morph_anc sel=morph_category_model, 1` includes bin 5. See [color_gnom.py](../src/3-pymol/color_gnom.py) for the display settings.

## Understand the annotations

The gnomAD mapping scripts store the following integer bins in the PDB B-factor field:

| Bin | Allele-frequency interval                         |
| --- | ------------------------------------------------- |
| 0   | `[1e-9, 1e-8)`                                    |
| 1   | `[1e-8, 1e-7)`                                    |
| 2   | `[1e-7, 1e-6)`                                    |
| 3   | `[1e-6, 1e-5)`                                    |
| 4   | `[1e-5, 1e-4)`                                    |
| 5   | `[1e-4, 1e-3)`                                    |
| 6   | `[1e-3, 1e-2)`                                    |
| 7   | `[1e-2, 0.51)`                                    |
| 8   | `[0.51, 1.0)`                                     |
| 42  | No retained frequency annotation at this position |

The upper endpoint is exclusive in the current implementation. **Known boundary issue:** a frequency of exactly `1.0` falls through to bin 0. Values outside the defined intervals also fall through to bin 0. Verify such entries in the source TSV before interpreting them.

The variant chains are composite visualization constructs. Their side chains are generated using pyMUT rotamers, with a fallback for detected clashes within the same chain. They do not represent validated mutant structures, observed haplotypes, or a full assessment of contacts with other chains or ligands.

## Rebuild the prepared data

The data-generation scripts cover these steps:

| Step              | Script                                                                                    | Purpose                                                                                                       |
| ----------------- | ----------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------- |
| Fetch             | [\_0-fetch.sh](../src/1-data_generation/server_wrappers/_0-fetch.sh)                      | Download the two AlphaMissense tables                                                                         |
| Fetch preparation | [0-fetch_AM.py](../src/1-data_generation/0-fetch_AM.py)                                   | Decompress and prepare the downloaded tables                                                                  |
| Parse             | [1-parse_AM_data.py](../src/1-data_generation/1-parse_AM_data.py)                         | Split `aa_substitutions` and `hg38` tables into per-UniProt files                                             |
| gnomAD            | [2-gnomad_1_file.py](../src/1-data_generation/2-gnomad_1_file.py)                         | Extract overall and category frequencies from a compatible gnomAD Hail table using the hg38 ranges            |
| Merge             | [3-missense_gnomXmiss_1_file.py](../src/1-data_generation/3-missense_gnomXmiss_1_file.py) | Merge AlphaMissense and gnomAD records, run VEP with ClinVar annotations, and filter/format the merged tables |
| Package           | [tar_amissenseXgnome.sh](../src/1-data_generation/server_wrappers/tar_amissenseXgnome.sh) | Package per-protein TSVs into archives and generate the UniProt lookup table                                  |

The separate [PyMOL helpers](../src/3-pymol/) color mapped structures.

### 1. Download and prepare the raw tables

From the repository root, submit the fetch job on a cluster that can reach Zenodo:

```bash
cd src/1-data_generation/server_wrappers
mkdir -p ../../../log ../../../log/gnomad ../../../log/gnomXmiss
sbatch ./_0-fetch.sh
```

The fetch job downloads the two required AlphaMissense files, activates the `operator` environment, decompresses them, removes the first three introductory lines, and writes the prepared tables to `data/1-raw/missense_db/`. Existing prepared files are reused. Complete `.gz` files staged there are prepared without downloading them again. Slurm writes job output to `log/`, and the wrapper also saves a timestamped log there. Wait for this job to finish before submitting parsing.

### 2. Submit parsing to Slurm

After the two raw TSV files are present, run:

```bash
sbatch ./_1-parse_AM_data.sh
```

The parsing job checks that the raw tables exist and writes Slurm logs to the repository's `log/` directory. The job writes per-protein tables to `data/2-processed/missense_jul/aa_substitutions/` and `data/2-processed/missense_jul/hg38/`, then prints the number of files written for each table and the total.

After both `.gz` files have been downloaded, the preparation and parsing Python scripts can also be run directly in the `operator` environment from this wrapper directory:

```bash
python ../0-fetch_AM.py ../../..
python ../1-parse_AM_data.py ../../..
```

Set `ALPHAMISSENSE_RECORD_ID` before submitting the fetch wrapper to select a different Zenodo record. The preparation Python script accepts an optional raw-data directory as its second positional argument; the parser accepts `--missense-dir` and repeatable `--tsv` options. Relative paths resolve against the project root supplied to each Python script. The preparation script processes downloaded files; it does not download them.

The gnomAD extraction and merging stages require a compatible gnomAD Hail table, Hail with Java/Spark, an Ensembl VEP Singularity container and cache, and a ClinVar VCF. Set their locations in `src/1-data_generation/server_wrappers/__config.sh`; the current wrappers target the gnomAD joint v4.1 sites table.

### 3. Run the gnomAD and merge arrays

From the same wrapper directory, submit each stage after the previous stage and its array tasks finish:

```bash
sbatch ./_2-gnomad-batch.sh
# Wait for its gnomAD array tasks to finish.
sbatch ./_3-gnomXmiss_sbatch.sh
# Wait for its merge array tasks to finish.
bash ./tar_amissenseXgnome.sh
```

The gnomAD and merge launchers each create a job list and submit a Slurm array. The archive script packages the completed merge outputs.

Data-generation shell scripts and `__config.sh` are in `src/1-data_generation/server_wrappers/`; the Python scripts and array ID lists are in `src/1-data_generation/`. Set external locations in each group's `__config.sh`, then adapt Slurm partition, QoS, and resources for your cluster. Submit jobs with `sbatch` from the wrapper directory. The `amissense_check.sh`, `output_counts.sh`, and `output_empties.sh` scripts check local data; their reports and categorized TSV copies are written under `log/`. The `rosetta_*_old.py` files are legacy implementations; the pyMUT scripts are the current mapping workflow.

Record the source releases, genome build, VEP/cache version, ClinVar date, filtering choices, package versions, and checksums with any regenerated dataset.



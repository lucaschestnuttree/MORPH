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

The mapping workflow expects:

- A per-protein merged TSV, normally named `<UniProt>_amissenseXgnomAD.tsv`. Output names are derived from the input TSV filename.
- The target protein in **chain A** of a PDB, with its sequence corresponding to the selected UniProt reference.
- **PDB residue numbers matching the one-based UniProt positions** used in `protein_variant`. A truncated protein can be used only when the remaining residues retain their canonical numbering. Automatic offset correction is not implemented.
- Continuous residues for pyMUT side-chain construction. Insertion codes are skipped by the mapping scripts' residue maps; normalize inputs explicitly before use.

## Map prepared data to a structure

### Overall allele frequency

```bash
conda activate pymut-env

python src/2-map_2_structure/pymut_gnomad_af_2_structure.py \
  data/2-processed/amissenseXgnomAD/P10635_amissenseXgnomAD.tsv \
  data/2-processed/pdb/P10635.pdb \
  data/2-processed/morph_outputs/af
```

Output: `data/2-processed/morph_outputs/af/P10635_AF_gnomad.pdb`.

The script ranks up to three retained variant rows by overall frequency at each position and builds variant chains `V`, `W`, and `X`. It writes the highest-ranked frequency bin to chains `A`, `V`, `W`, and `X`; the B-factor values on `W` and `X` therefore do not encode their individual variant frequencies.

### Category-specific frequency

```bash
conda activate pymut-env

python src/2-map_2_structure/pymut_gnomad_anc_2_structure.py \
  data/2-processed/amissenseXgnomAD/P10635_amissenseXgnomAD.tsv \
  data/2-processed/pdb/P10635.pdb \
  data/2-processed/morph_outputs/categories \
  nfe,sas,afr,fin,amr
```

Output: `data/2-processed/morph_outputs/categories/P10635_category_gnomad.pdb`.

For each category, the script selects the most frequent retained substitution at each position. The current script assigns categories in input order to chains `T`, `U`, `V`, `W`, `X`, `Y`, and `Z`, supporting up to seven categories. Read the generated PDB remarks to recover the assignments. Chain A receives the maximum frequency bin across the selected categories.

The extraction pipeline includes `nfe`, `fin`, `mid`, `remaining`, `ami`, `asj`, `afr`, `eas`, `sas`, `amr`, `XX`, and `XY` frequency columns. `XX` and `XY` are sex categories. Availability depends on the input table; an unavailable category currently produces a warning and is skipped during frequency calculation.

## View in PyMOL

Install the local helper scripts before running these commands. In PyMOL, open **Plugin > Plugin Manager > Install Plugin > Install from local file**, then select `MORPH\src\3-pymol\color_gnom.py`. Restart PyMOL.

After installation, start a PyMOL session.

For an overall-frequency structure:

```pymol
load data/2-processed/morph_outputs/af/P10635_AF_gnomad.pdb, P10635_AF
morph_af sel=P10635_AF
```

or with a broader threshold:

```pymol
load data/2-processed/morph_outputs/af/P10635_AF_gnomad.pdb, P10635_AF
morph_af sel=P10635_AF, 1
```

For a category-mapped file:

```pymol
load data/2-processed/morph_outputs/categories/P10635_category_gnomad.pdb, P10635_cat
morph_anc sel=P10635_cat
```

or with a broader threshold:

```pymol
load data/2-processed/morph_outputs/categories/P10635_category_gnomad.pdb, P10635_cat
morph_anc sel=P10635_cat, 1
```

The threshold controls which variant and side-chain residues are selected by their integer B-factor bins.

For `morph_af`, the default threshold selects bins 5 to 8 `[1e-4, 1.0]`; threshold `1` expands this to bins 4 to 8 `[1e-5, 1.0]`. Here each bin has a distinct color as seen in table below.

For `morph_anc`, the default selects bins 6 to 8 (approximately `[1e-3, 1.0]`); threshold `1` expands this to bins 5 to 8 (approximately `[1e-4, 1.0]`). Chain A is still colored across all bins. See [color_gnom.py](../src/3-pymol/color_gnom.py) for the display settings.

## Understand the annotations

The mapping scripts store integer frequency bins in the PDB B-factor field. Both `morph_af` and `morph_anc` use the same bin definitions; their threshold settings select different subsets of those bins:

| Bin | Allele-frequency interval        | `morph_af` color                           | `morph_anc` use                               |
| --- | -------------------------------- | ------------------------------------------ | --------------------------------------------- |
| 0   | `[1e-9, 1e-8)`                   | light grey                                 | not bin-colored                               |
| 1   | `[1e-8, 1e-7)`                   | grey                                       | not bin-colored                               |
| 2   | `[1e-7, 1e-6)`                   | darker grey                                | not bin-colored                               |
| 3   | `[1e-6, 1e-5)`                   | dark grey                                  | not bin-colored                               |
| 4   | `[1e-5, 1e-4)`                   | light red at threshold `1`; otherwise grey | not bin-colored                               |
| 5   | `[1e-4, 1e-3)`                   | red                                        | lighter category-chain color at threshold `1` |
| 6   | `[1e-3, 1e-2)`                   | dark red                                   | category-chain color                          |
| 7   | `[1e-2, 0.51)`                   | blue                                       | category-chain color                          |
| 8   | `[0.51, 1.0)`                    | purple                                     | category-chain color                          |
| 42  | No retained frequency annotation | pale yellow                                | not bin-colored                               |

Variant chains are composite visualization constructs. Their side chains are generated using pyMUT rotamers, with a fallback for detected clashes within the same chain. They do not represent validated mutant structures, observed haplotypes, or a full assessment of contacts with other chains or ligands.

## Rebuild the prepared data

Before submitting jobs:

1. Create Conda environments from the YAML & requirements.
2. Install or obtain Hail with Java/Spark, a compatible Hail table, VEP Singularity image/cache, and a ClinVar VCF for the gnomAD workflow.
3. Edit src/1-data_generation/server_wrappers/\_\_config.sh with the Conda initialization, data paths, Hail table, VEP image/cache, ClinVar VCF, and Slurm settings.

| Step   | Slurm/Bash wrapper                                                                          | Python script                                                                               | Environment and purpose                          |
| ------ | ------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------- | ------------------------------------------------ |
| Fetch  | [`_0-fetch.sh`](../src/1-data_generation/server_wrappers/_0-fetch.sh)                       | [`0-fetch_AM.py`](../src/1-data_generation/0-fetch_AM.py)                                   | `operator`; download and prepare AlphaMissense   |
| Parse  | [`_1-parse_AM_data.sh`](../src/1-data_generation/server_wrappers/_1-parse_AM_data.sh)       | [`1-parse_AM_data.py`](../src/1-data_generation/1-parse_AM_data.py)                         | `operator`; create per-UniProt tables            |
| gnomAD | [`_2-gnomad-batch.sh`](../src/1-data_generation/server_wrappers/_2-gnomad-batch.sh)         | [`2-gnomad_1_file.py`](../src/1-data_generation/2-gnomad_1_file.py)                         | `hail`; extract overall and category frequencies |
| Merge  | [`_3-gnomXmiss_sbatch.sh`](../src/1-data_generation/server_wrappers/_3-gnomXmiss_sbatch.sh) | [`3-missense_gnomXmiss_1_file.py`](../src/1-data_generation/3-missense_gnomXmiss_1_file.py) | `operator` plus VEP/ClinVar; merge and annotate  |
| -      | package completed per-protein outputs                                                       |

From `src/1-data_generation/server_wrappers/`, submit each stage after the previous one completes:

```bash
sbatch ./_0-fetch.sh
sbatch ./_1-parse_AM_data.sh
sbatch ./_2-gnomad-batch.sh
sbatch ./_3-gnomXmiss_sbatch.sh
```

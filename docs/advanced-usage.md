# Advanced usage

[Back to the project overview](../README.md)

Use these workflows to map prepared data onto structures or to rebuild the datasets. The guided notebook is the intended entry point for users who prefer a step-by-step interface.

Examples below use Bash and run from the repository root. Replace the example paths with your own files. The public prepared-data release is pending, so the mapping commands currently require locally available per-protein TSVs.

## Environment setup

The repository contains three development environment definitions:

| Environment | Intended role | Definition |
| --- | --- | --- |
| `operator` | Downloading, parsing, and annotating data | [environment-operator.yml](../requirements/environment-operator.yml) |
| `hail` | Querying the gnomAD Hail table | [environment-hail.yml](../requirements/environment-hail.yml) |
| `pymut-env` | Mapping scores and alternative side chains | [environment-pymut.yml](../requirements/environment-pymut.yml) |

These files are development snapshots. Each YAML currently includes a machine-specific `prefix:`; remove or adapt it before creating an environment elsewhere. The Hail YAML does not yet declare Hail itself or all Python dependencies used by the extraction script. Several `requirements-*.txt` entries point to local Conda build paths (`file:///home/conda/...`), so those exports are not portable pip installation instructions.

For structure mapping, a minimal starting environment is:

```bash
conda create -n morph-mapping -c conda-forge python=3.10 "numpy<2" pandas biopython
conda activate morph-mapping
```

The gnomAD mapping scripts also import [pyMUT](https://github.com/gerdos/pyMUT). Clone it and make its directory importable, following its upstream library usage:

```bash
git clone https://github.com/gerdos/pyMUT.git ../pyMUT
export PYTHONPATH="$(cd ../pyMUT && pwd)${PYTHONPATH:+:$PYTHONPATH}"
python -c "import numpy, pandas, Bio, pymut; print('Mapping imports OK')"
```

This is a setup starting point, rather than a validated release environment. Record the pyMUT commit and package versions used for reproducible work. PyMOL helpers run inside PyMOL, which needs its own installation.

## Input requirements

The gnomAD workflows expect:

- A per-protein merged TSV, normally named `<UniProt>_amissenseXgnomAD.tsv`. Output names are derived from the input TSV filename.
- The target protein in **chain A** of a PDB, with its sequence corresponding to the selected UniProt reference.
- **PDB residue numbers matching the one-based UniProt positions** used in `protein_variant`. A truncated protein can be used only when the remaining residues retain their canonical numbering. Automatic offset correction is not implemented.
- Continuous residues for pyMUT side-chain construction. Insertion codes are skipped by the mapping scripts' residue maps; normalize inputs explicitly before use.

The overall-frequency script requires `protein_variant`, `am_class`, `allele_count`, and `allele_frequency`. Category mapping also requires `uniprot_id`, and frequency columns named `<category>_allele_frequency`. It accepts `allele_frequency` as a category for the overall-frequency column. An absent `intron` column is created internally.

The notebook additionally validates that the UniProt entry is reviewed and human, and that prepared data exists for that accession. The command-line scripts do not perform that complete validation; prepare and verify the inputs before running them.

## Map prepared data to a structure

### Overall allele frequency

```bash
python src/2-mapping/map_2_structure/pymut_gnomad_af_2_structure.py \
  data/2-processed/amissenseXgnomAD/P10635_amissenseXgnomAD.tsv \
  data/2-processed/pdb/P10635.pdb \
  data/2-processed/morph_outputs/af
```

Output: `data/2-processed/morph_outputs/af/P10635_AF_gnomad.pdb`.

The script ranks up to three retained variant rows by overall frequency at each position and builds variant chains `V`, `W`, and `X`. It writes the highest-ranked frequency bin to chains `A`, `V`, `W`, and `X`; the B-factor values on `W` and `X` therefore do not encode their individual variant frequencies.

### Category-specific frequency

```bash
python src/2-mapping/map_2_structure/pymut_gnomad_anc_2_structure.py \
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
python src/2-mapping/map_2_structure/all_aa_AM_2_structure.py \
  data/2-processed/missense_jul/aa_substitutions/P10635_aa_substitutions.tsv \
  data/2-processed/pdb/P10635.pdb \
  data/2-processed/morph_outputs/am/P10635_AM.pdb
```

This mapper requires `protein_variant` and `am_pathogenicity` columns. It averages the supplied rows at each position, rounds to two decimal places, and writes that value to matching ATOM records.

It currently matches residue numbers across **all protein chains**, rather than restricting annotations to chain A. Use an input containing only the intended target protein, or explicitly prepare the chains you intend to annotate. Unscored residues retain their input B-factors, so a generic color scale can mix original values with mapped scores.

## View in PyMOL

From a PyMOL session whose working directory is the repository root:

```pymol
run src/2-mapping/pymol/color_gnom.py
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
run src/2-mapping/pymol/color_am.py
load data/2-processed/morph_outputs/am/P10635_AM.pdb, morph_am_model
coloram sel=morph_am_model
```

The commands select and color different frequency ranges by default. `morph_af sel=morph_af_model, 1` broadens the highlighted range to include bin 4; `morph_anc sel=morph_category_model, 1` includes bin 5. See [color_gnom.py](../src/2-mapping/pymol/color_gnom.py) for the display settings.

## Understand the annotations

The gnomAD mapping scripts store the following integer bins in the PDB B-factor field:

| Bin | Allele-frequency interval |
| --- | --- |
| 0 | `[1e-9, 1e-8)` |
| 1 | `[1e-8, 1e-7)` |
| 2 | `[1e-7, 1e-6)` |
| 3 | `[1e-6, 1e-5)` |
| 4 | `[1e-5, 1e-4)` |
| 5 | `[1e-4, 1e-3)` |
| 6 | `[1e-3, 1e-2)` |
| 7 | `[1e-2, 0.51)` |
| 8 | `[0.51, 1.0)` |
| 42 | No retained frequency annotation at this position |

The upper endpoint is exclusive in the current implementation. **Known boundary issue:** a frequency of exactly `1.0` falls through to bin 0. Values outside the defined intervals also fall through to bin 0. Verify such entries in the source TSV before interpreting them.

The variant chains are composite visualization constructs. Their side chains are generated using pyMUT rotamers, with a fallback for detected clashes within the same chain. They do not represent validated mutant structures, observed haplotypes, or a full assessment of contacts with other chains or ligands.

## Rebuild the prepared data

The numbered scripts implement the following stages:

| Stage | Script | Purpose |
| --- | --- | --- |
| 1 | [1-fetch_data.py](../src/1-data_generation/1-fetch_data.py) | Download AlphaMissense files from its upstream Zenodo record, decompress them, and remove introductory lines |
| 2 | [2-parse_AM_data.py](../src/1-data_generation/2-parse_AM_data.py) | Split AlphaMissense tables into per-UniProt files |
| 3 | [3-gnom_1_file.py](../src/1-data_generation/3-gnom_1_file.py) | Extract overall and category frequencies from a compatible gnomAD Hail table using the hg38 ranges |
| 4 | [4-missense_gnomXmiss_1_file.py](../src/1-data_generation/4-missense_gnomXmiss_1_file.py) | Merge AlphaMissense and gnomAD records, run VEP with ClinVar annotations, and filter/format the merged tables |
| 5 | [5--tar_amissenseXgnome.sh](../src/1-data_generation/5--tar_amissenseXgnome.sh) | Package per-protein TSVs into archives and generate the UniProt lookup table |

Stages 3–4 require a compatible gnomAD Hail table, Hail and its Java/Spark environment, an Ensembl VEP Singularity container and cache, and a ClinVar VCF. The current wrappers are configured for the gnomAD joint v4.1 sites table.

The accompanying `.sh` files are examples from the development cluster. Adapt `BASE_PATH`, Conda activation paths, data/container/cache paths, log locations, and SLURM partition/QoS/resource settings before submitting jobs. The `z_*.sh` files are local data-check utilities. The `rosetta_*_old.py` files are legacy implementations; the pyMUT scripts are the current mapping workflow.

Record the source releases, genome build, VEP/cache version, ClinVar date, filtering choices, package versions, and checksums with any regenerated dataset.

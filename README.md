# MORPH (Missense on Reference Protein Human)

**Human missense variation, mapped onto protein structure.**

MORPH (_Missense on Reference Protein Human_) connects gnomAD population frequencies and human protein structure predictions.

The **guided notebook(s) are the intended starting point for researchers who want to reproduce visuals for their protein of interest**.

The Python and shell scripts support advanced users who want to rebuild the underlying datasets or do proteome level analyses.

> **Release status:** The three workflow notebooks are available at the repository root, and the prepared dataset is publicly available on [Zenodo](https://zenodo.org/records/18163874).

## What you can explore

- **An optional structure prediction:** the notebook includes a Boltz-2 workflow when you need to generate an input structure.
- **Population variation in structural context:** map gnomAD missense allele frequencies onto a reference protein and explore where observed missense variants occur in 3D.
- **Differences across ancestry groups:** build category-specific views using the frequency columns available in the prepared data.
- **Alternative side chains:** Inspect alternative amino-acid side chains alongside the reference structure at a glance in PyMOL using pyMUT-generated chains.

## Choose your workflow

| Your goal                                          | Where to start                                                                  | Experience needed                                          |
| -------------------------------------------------- | ------------------------------------------------------------------------------- | ---------------------------------------------------------- |
| Generate a structure                               | [1-Boltz2.ipynb](1-Boltz2.ipynb)                                                | Google Colab GPU runtime                                   |
| Map prepared variation onto a structure            | [2-MORPH.ipynb](2-MORPH.ipynb)                                                  | Google Colab CPU runtime                                   |
| Explore mapped structures                          | [3-pymol.ipynb](3-pymol.ipynb)                                                  | PyMOL                                                      |
| Map prepared variant data onto your own structures | [Advanced usage guide](docs/advanced-usage.md#map-prepared-data-to-a-structure) | Python, command line, and environment setup                |
| Rebuild or extend the prepared datasets            | [Data-generation pipeline](docs/advanced-usage.md#rebuild-the-prepared-data)    | Bioinformatics, Hail, Ensembl VEP, and Linux/HPC workflows |

## Start with the notebook

The workflow is split into three notebooks: [1-Boltz2.ipynb](1-Boltz2.ipynb) generates an optional input structure, [2-MORPH.ipynb](2-MORPH.ipynb) validates and maps variation, and [3-pymol.ipynb](3-pymol.ipynb) walks through visualization in PyMOL. The original combined notebook remains available as `MORPH_complete.ipynb`.

### The intended workflow

1. **Choose your protein.** Enter its UniProt accession and validate the reviewed human entry.
2. **Provide a structure.** Upload a suitable PDB, or generate and download one through `1-Boltz2.ipynb`.
3. **Retrieve prepared data.** The notebook uses a UniProt lookup table to select the relevant Zenodo archive.
4. **Choose a mapping mode.** Use overall allele frequency or select genetic ancestry categories.
5. **Download and explore.** Export the annotated PDB and follow `3-pymol.ipynb` or the [visualization helpers](docs/advanced-usage.md#view-in-pymol).

## Prepared data

**Public dataset: [Zenodo record 18163874](https://zenodo.org/records/18163874).**

| File                            | Purpose                                                                                                             |
| ------------------------------- | ------------------------------------------------------------------------------------------------------------------- |
| `lookup_amissenseXgnomAD.tsv`   | Maps each available UniProt accession to an archive and its per-protein TSV                                         |
| `amissenseXgnomAD_part*.tar.gz` | Batches of per-protein variant tables combining AlphaMissense data, gnomAD frequencies, and VEP/ClinVar annotations |

## What the outputs represent

| Mapping                    | Output                          | Meaning                                                                                         |
| -------------------------- | ------------------------------- | ----------------------------------------------------------------------------------------------- |
| Overall gnomAD frequency   | `<UniProt>_AF_gnomad.pdb`       | Frequency bins on the reference chain, with up to three ranked variant chains (`V`, `W`, `X`)   |
| Selected gnomAD categories | `<UniProt>_category_gnomad.pdb` | One variant chain per selected category; chain/category assignments are recorded in PDB remarks |
| AlphaMissense              | User-named PDB                  | Mean AlphaMissense score per residue, calculated from the supplied variant table                |

The mapped PDB **B-factor field stores visualization values**. For gnomAD outputs these are frequency bins; for AlphaMissense outputs they are mean scores. Keep the original structure if you need its original B-factors or model-confidence values.

Variant chains combine substitutions across positions to help display alternative side chains. They are visualization constructs, rather than experimentally determined mutant proteins or observed haplotypes. A position without retained frequency data should not be interpreted as proof that variation is absent.

See the [advanced guide](docs/advanced-usage.md#understand-the-annotations) for the bin scheme, chain behavior, and current implementation limitations.

## Repository layout

```text
README.md                         Project overview and notebook entry point
docs/advanced-usage.md            Command-line workflows and annotation details
requirements/                    Development environment definitions and exports
src/
1-Boltz2.ipynb                      Optional Boltz-2 structure prediction notebook
2-MORPH.ipynb                       Variation mapping notebook
3-pymol.ipynb                       PyMOL visualization walkthrough
MORPH_complete.ipynb                Original combined notebook
   1-data_generation/             Download, harmonize, annotate, and package data
   2-map_2_structure/             Map prepared values and substitutions to PDBs
   3-pymol/                        PyMOL coloring and selection helpers
data/                            Local raw and processed data (ignored by Git)
log/                             Local execution logs (ignored by Git)
```

## Sources and acknowledgements

MORPH builds on the following datasets and tools:

- **AlphaMissense:** Cheng et al., [_Accurate proteome-wide missense variant effect prediction with AlphaMissense_](https://www.science.org/doi/10.1126/science.adg7492), with [precomputed predictions on Zenodo](https://zenodo.org/records/10813168).
- **gnomAD:** [population variant frequencies](https://gnomad.broadinstitute.org/).
- **UniProt:** [reference protein sequences and metadata](https://www.uniprot.org/).
- **Boltz-2:** [optional biomolecular structure prediction](https://github.com/jwohlwend/boltz).
- **pyMUT:** [rotamer-based side-chain mutation](https://github.com/gerdos/pyMUT).
- **Ensembl VEP and ClinVar:** annotations used by the data-generation pipeline.
- **Tamas Hegedus, hegelab.org:** author credited in the supplied AlphaMissense PyMOL coloring helper.

When reporting results, cite the upstream datasets and tools used, and the public MORPH Zenodo record linked above.

**Project license:** MIT. Upstream datasets and dependencies retain their own licenses and attribution requirements.

## License

Copyright (c) 2026 Lucas Chataigner

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

# MORPH (Missense on Reference Protein Human)

**Human missense variation, mapped onto protein structure.**

MORPH (_Missense on Reference Protein Human_) connects gnomAD population frequencies and human protein structure predictions.

The **guided notebook(s) are the intended starting point for those who want to reproduce visuals for their protein of interest with little code knowledge**.

The Python and shell scripts support users who want to run jobs on their machines, rebuild the underlying datasets, or do proteome level analyses.

## Start with the notebook(s)

1. **Boltz2.ipynb** [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/lucaschestnuttree/MORPH/blob/main/1-Boltz2.ipynb) - `Generate an input structure with Boltz-2`

2. **MORPH.ipynb** [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/lucaschestnuttree/MORPH/blob/main/2-MORPH.ipynb) - `Validate and map variations`

3. **PyMOL.ipynb** [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/lucaschestnuttree/MORPH/blob/main/3-PyMOL.ipynb) - `Walk through of PyMOL visualization.`

## What you can explore

- **An optional structure prediction:** includes a Boltz-2 notebook when you need to generate an input structure.
- **Population variation in structural context:** map gnomAD missense allele frequencies onto a reference protein prediction and explore where observed missense variants occur in 3D.
- **Differences across ancestry groups:** build category-specific views using the frequency columns available in the prepared data.
- **Alternative side chains:** Inspect alternative amino-acid side chains alongside the reference structure at a glance in PyMOL using pyMUT-generated chains.

## Choose your workflow

| Your goal                                 | Where to start                                                               | Experience needed                                          |
| ----------------------------------------- | ---------------------------------------------------------------------------- | ---------------------------------------------------------- |
| Generate a structure                      | [1-Boltz2.ipynb](1-Boltz2.ipynb)                                             | Google Colab GPU runtime                                   |
| Map prepared variation onto a structure   | [2-MORPH.ipynb](2-MORPH.ipynb)                                               | Google Colab CPU runtime                                   |
| Explore mapped structures                 | [3-PyMOL.ipynb](3-PyMOL.ipynb)                                               | PyMOL                                                      |
| Map prepared variant data on your machine | [Advanced usage guide](docs/advanced-usage.md)                               | Python, command line, and environment setup                |
| Rebuild or extend the prepared datasets   | [Data-generation pipeline](docs/advanced-usage.md#rebuild-the-prepared-data) | Bioinformatics, Hail, Ensembl VEP, and Linux/HPC workflows |

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

The mapped PDB **B-factor field stores visualization values**. For gnomAD outputs these are frequency bins. Keep the original structure if you need its original B-factors or model-confidence values.

Variant chains combine substitutions across positions to help display alternative side chains. They are visualization constructs, rather than experimentally determined mutant proteins or observed haplotypes. A position without retained frequency data should not be interpreted as proof that variation is absent.

See the [advanced guide](docs/advanced-usage.md#understand-the-annotations) for the bin scheme, chain behavior, and current implementation limitations.

## Repository layout

```text
README.md                         Project overview and notebook entry point
docs/advanced-usage.md            Command-line workflows and annotation details
requirements/                    Development environment definitions and exports
src/
1-Boltz2.ipynb                      Boltz-2 structure prediction notebook
2-MORPH.ipynb                       Variation mapping notebook
3-PyMOL.ipynb                       PyMOL visualization walkthrough
   1-data_generation/             Download, harmonize, annotate, and package data
   2-map_2_structure/             Map prepared values and substitutions to PDBs
   3-PyMOL/                        PyMOL coloring and selection helpers
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

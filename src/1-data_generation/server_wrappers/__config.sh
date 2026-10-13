#!/usr/bin/env bash
# External locations used by this folder's scripts.
# Edit these defaults for your cluster, or export a variable before running or
# submitting a wrapper. Environment variables take precedence over these defaults.
# All project-owned paths are derived by the wrappers from their script location.

# Conda initialization script (environment names remain in the wrappers).
CONDA_SH="${CONDA_SH:-/pasteur/appa/homes/lchataig/miniconda3/etc/profile.d/conda.sh}"

# gnomAD Hail table and scratch directory used for Hail logs.
HT_PATH="${HT_PATH:-/pasteur/helix/projects/GGS_JASS/WKD_Lucas/gnomad.joint.v4.1.sites.ht}"
LOG_DIR="${LOG_DIR:-/pasteur/appa/scratch/lchataig}"

# Ensembl VEP container, cache and ClinVar annotation file.
VEP_SIF="${VEP_SIF:-/pasteur/helix/projects/GGS_JASS/WKD_Lucas/ensemble-VEP/vep.sif}"
CACHE_DIR="${CACHE_DIR:-/pasteur/helix/projects/GGS_JASS/WKD_Lucas/ensemble-VEP/vep_data}"
CLINVAR_VCF="${CLINVAR_VCF:-/pasteur/helix/projects/GGS_JASS/WKD_Lucas/ensemble-VEP/custom_annotations/clinvar.vcf.gz}"


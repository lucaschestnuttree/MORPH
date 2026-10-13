#!/usr/bin/env bash
# External locations used by the wrappers in this folder.
# Edit this default for your cluster, or export CONDA_SH before running or
# submitting a wrapper. An environment variable takes precedence over this default.
# All project-owned paths are derived by the wrappers from their script location.

# Conda initialization script (the pymut-env environment name stays in the wrappers).
CONDA_SH="${CONDA_SH:-/pasteur/appa/homes/lchataig/miniconda3/etc/profile.d/conda.sh}"

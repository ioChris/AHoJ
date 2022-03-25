#!/bin/bash

# Executes apoholo_J.py inside required conde env

CONDA_SH="/usr/share/miniconda/etc/profile.d/conda.sh"
if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    export CONDA_SH="$HOME/miniconda3/etc/profile.d/conda.sh"
fi
source "$CONDA_SH"

conda activate ahoj
PYTHON="$CONDA_PREFIX/bin/python"
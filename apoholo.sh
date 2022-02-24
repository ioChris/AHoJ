#!/bin/bash

# Executes apoholo_J.py inside required conde env

source ~/miniconda3/etc/profile.d/conda.sh

conda activate ahoj
exec python apoholo_J.py "$@"
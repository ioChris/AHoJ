#!/bin/bash

# TODO: consider syncing to tmp dir and then swapping

TARGET_DIR="$1"  # dir used by AHoJ with potentially extra data (may be by running process from webapp)

echo "Synchronizing directory '$TARGET_DIR' with PDB"

rsync -rlpt -v -z --delete rsync.ebi.ac.uk::pub/databases/pdb/data/structures/divided/mmCIF/ "$TARGET_DIR"

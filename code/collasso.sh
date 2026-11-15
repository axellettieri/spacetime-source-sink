#!/bin/bash
# Continues past the turnaround to the collapse (same parameters as futuro2). Needs parametri3d_originale.json in the folder.
set -e
python universo3d.py collasso --np ${N:-512} --ng ${N:-512} --L ${L:-1000} --passi ${P:-1100} --par parametri3d_originale.json --lam inf --afin 8 --collasso
python riassunto3d.py | grep -A4 '^collasso'

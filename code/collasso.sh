#!/bin/bash
# Prosegue oltre l'inversione fino al collasso (stessi parametri di futuro2). Serve parametri3d_originale.json nella cartella.
set -e
python universo3d.py collasso --np ${N:-512} --ng ${N:-512} --L ${L:-1000} --passi ${P:-1100} --par parametri3d_originale.json --lam inf --afin 8 --collasso
python riassunto3d.py | grep -A4 '^collasso'

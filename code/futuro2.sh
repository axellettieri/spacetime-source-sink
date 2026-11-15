#!/bin/bash
# Future with the absorption written as in the original model (does not feed on itself). Run in the folder of the previous runs.
set -e
python adatta3d.py diffuso.json originale
python universo3d.py futuro2 --np ${N:-512} --ng ${N:-512} --L ${L:-1000} --passi ${P:-1100} --par parametri3d_originale.json --lam inf --afin 8
python riassunto3d.py | grep -A2 futuro

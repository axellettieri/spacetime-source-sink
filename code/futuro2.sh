#!/bin/bash
# Futuro con l'assorbimento scritto come nel modello originale (non si autoalimenta). Lanciare nella cartella delle corse precedenti.
set -e
python adatta3d.py diffuso.json originale
python universo3d.py futuro2 --np ${N:-512} --ng ${N:-512} --L ${L:-1000} --passi ${P:-1100} --par parametri3d_originale.json --lam inf --afin 8
python riassunto3d.py | grep -A2 futuro

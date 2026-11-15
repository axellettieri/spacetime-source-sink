#!/bin/bash
# Full uncertainty chain with the emergent clock (full Planck + DESI + Pantheon+), then the summary.
# If interrupted, rerun the same command: it resumes from the last save (every 20 steps).
set -e
NP=${NP:-16}; PASSI=${PASSI:-3000}
python mcmc_orologio.py fine_L250_s1.json 10 $PASSI $NP
python analizza_orologio.py fine_L250_s1.json 10 | tee analisi_orologio.txt

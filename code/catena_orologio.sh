#!/bin/bash
# Catena completa delle incertezze con l'orologio emergente (Planck completo + DESI + Pantheon+), poi il riassunto.
# Se si interrompe, rilanciare lo stesso comando: riprende dall'ultimo salvataggio (ogni 20 passi).
set -e
NP=${NP:-16}; PASSI=${PASSI:-3000}
python mcmc_orologio.py fine_L250_s1.json 10 $PASSI $NP
python analizza_orologio.py fine_L250_s1.json 10 | tee analisi_orologio.txt

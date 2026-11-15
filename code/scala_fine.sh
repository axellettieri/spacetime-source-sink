#!/bin/bash
# Orologio degli eventi a risoluzione piu' fine, su 4 GPU in parallelo (una corsa per GPU).
# Scatole piu' piccole = celle piu' piccole = aloni di galassie risolti meglio. Poi adatta ogni orologio ai dati.
# Servono nella cartella: fit3.py, i dati scaricati (scarica_dati.py) e results3.json (python fit3.py).
set -e
N=${N:-512}; P=${P:-1500}; SOGLIE=3,10,30,100,200
[ -f results3.json ] || python fit3.py
i=0
for CFG in "500 1" "250 1" "125 1" "250 2"; do
  set -- $CFG; NOME=fine_L$1_s$2
  [ -f $NOME.json ] || CUDA_VISIBLE_DEVICES=$((i % ${NGPU:-4})) python universo3d.py $NOME --np $N --ng $N --L $1 --passi $P --xinit 0.3107 --dc 50 --soglie $SOGLIE --seme $2 > $NOME.log 2>&1 &
  i=$((i+1))
  [ $((i % ${NGPU:-4})) -eq 0 ] && wait
done
wait
for F in fine_L*.json; do echo "=== $F"; python adatta3d.py $F originale | tee ${F%.json}_adattamento.txt; mv parametri3d_originale.json ${F%.json}_parametri.json; done

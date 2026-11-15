#!/bin/bash
# Catena completa su GPU. Ogni corsa salva lo stato ogni 100 passi: se si interrompe, rilanciare questo script e riprende.
set -e
N=${N:-512}; L=${L:-1000}; P=${P:-600}
[ -f orologio.json ]   || python universo3d.py orologio --np $N --ng $N --L $L --passi $P --xinit 0.3107 --dc 20 --soglie 1.5,2,3,5,10
python adatta3d.py orologio.json
[ -f diffuso.json ]    || python universo3d.py diffuso   --np $N --ng $N --L $L --passi $P --par parametri3d.json --lam inf
for LAM in 300 100 30 10; do
[ -f locale$LAM.json ] || python universo3d.py locale$LAM --np $N --ng $N --L $L --passi $P --par parametri3d.json --lam $LAM
done
[ -f futuro.json ]     || python universo3d.py futuro    --np $N --ng $N --L $L --passi $((P*3/2)) --par parametri3d.json --lam inf --afin 6
python riassunto3d.py

# Source–sink model: analysis code

English guide to the scripts. The original Italian notes are in `LEGGIMI.txt`. Comments, docstrings and screen output are in English; file, variable and function names are still in Italian (see the glossary at the end), so that the archived code matches the published results exactly.

## Requirements

- Python 3.10 or later
- numpy and scipy (everything)
- emcee (`mcmc.py`, `mcmc_cmb.py`, `mcmc_orologio.py`)
- camb (`cmb_full.py`, `mcmc_cmb.py`, `cmb_orologio.py`, `mcmc_orologio.py`)
- astropy and pandas (`ambiente.py`)
- CuPy with an NVIDIA GPU (`universo3d.py`; without it the script falls back to NumPy, which is fine only for small test runs)

Run all scripts from this folder.

## Data fits

| Command | What it does | Time |
|---|---|---|
| `python scarica_dati.py` | Downloads the public data once (Pantheon+, DESI DR2, Union3, DES, planck-lite-py, 2M++ map) | about 210 MB |
| `python fit3.py` | Main fit: supernovae + DESI + compressed Planck priors; writes `results3.json` | 30 s |
| `python profile.py` | Quick error bars on turnaround and collapse times | 1 min |
| `python mcmc.py cum_n3 2500` | Full sampling with the compressed Planck priors | 5 min |
| `python growth.py` | Structure growth, fσ8 and S8 | 10 s |
| `python controllo_crescita.py` | Source of the gravitational potential in the three versions of the model | 10 s |
| `python check_growth.py` | Accuracy of the effective-fluid approximation used for the CMB | 10 s |
| `python orologi.py` | Alternative event clocks and energy budget | 10 min |
| `python shells.py` | Distance-shell test for local expansion | 10 s |
| `python ambiente.py` | Environment test with the 2M++ density and velocity map | 1 min |
| `python versione_locale.py` | Local version of the model | short |
| `python cmb_full.py LCDM U` | Full CMB with CAMB, ΛCDM | 10–25 min |
| `python cmb_full.py cum_n3 U` | Full CMB, sink proportional to accumulated mass | 20–30 min |
| `python cmb_full.py prop_n3 U` | Same, sink proportional to created space | 20–30 min |
| `python cmb_full.py w0wa U` | Same, empirical w0wa curve | 20–30 min |
| `python mcmc_cmb.py cum_n3 3000 8` | Long sampling with 8 processes; can be stopped and resumed | about 14 h |
| `python traccia_w0wa.py [cum\|prop\|const]` | The forced line in the (w0, wa) plane predicted by the model | 1 min |

`profile.py`, `mcmc.py`, `growth.py`, `controllo_crescita.py` and `orologi.py` read `results3.json`, so run `fit3.py` first.

`fit.py` and `fit2.py` are the first exploratory fits (supernovae + DESI only). They are superseded by `fit3.py` and kept for reference.

If `scarica_dati.py` cannot download planck-lite-py, clone it into this folder: `git clone https://github.com/heatherprince/planck-lite-py`.

### Other supernova samples

Set the environment variable `SNSET` to `union3`, `desy5` or `dovekie`. Without it, Pantheon+ is used.

```
# Linux
SNSET=union3 python fit3.py ; SNSET=union3 python cmb_full.py LCDM U ; SNSET=union3 python cmb_full.py cum_n3 U
```

```powershell
# PowerShell
$env:SNSET='union3'; python fit3.py; python cmb_full.py LCDM U; python cmb_full.py cum_n3 U; Remove-Item Env:SNSET
```

Results go to `results3_<sample>.json` and `cmb2_<model>_U_<sample>.json`.

## GPU N-body simulations (emergent clock)

| Script | What it does |
|---|---|
| `universo3d.py` | Particle-mesh N-body simulation with gravity. Events (matter collapsing) create space and masses absorb it. The event clock is not imposed: the simulation produces it. |
| `adatta3d.py` | Fits the model with the simulation's clock to supernovae, DESI and Planck; writes `parametri3d.json`. |
| `corsa3d.sh` | Full chain: clock, fit, diffuse run, local runs (300, 100, 30, 10 Mpc), future. |
| `riassunto3d.py` | Final summary table. |
| `futuro2.sh` | Future with the original sink form (turnaround in about 29 billion years). |
| `collasso.sh` | Like `futuro2.sh` but continues past the turnaround to the collapse (`--collasso` option). About 1–2 h on a V100. |
| `scala_fine.sh` | Event clock at fine resolution (boxes of 500, 250, 125 Mpc) and fit. Designed for 4 GPUs; with fewer, use `NGPU=1 ./scala_fine.sh`. About 1 h per run. |

They need the data already downloaded with `scarica_dati.py`.

Quick test (a few minutes): `N=128 L=400 P=200 ./corsa3d.sh`. Full run: `./corsa3d.sh` (512³ particles, about 9 GB of GPU memory).

## Full Planck with the emergent clock

```
python cmb_orologio.py fine_L250_s1.json 10 500
```

Full Planck + DESI + supernovae with the clock from a simulation file and a collapse threshold. Takes about 30–40 minutes on one core. Other samples work with `SNSET=union3 ...` as above.

### Uncertainty chain

`catena_orologio.sh` runs `mcmc_orologio.py` (3000 steps, 32 walkers) and then `analizza_orologio.py`.

- `NP` sets the number of parallel processes (16 on a Ryzen 1700). The chain resumes automatically if interrupted.
- The result goes to `analisi_orologio.txt`: parameter intervals, H0, matter fraction today, and the (w0, wa) point.
- It needs the same data, planck-lite-py and CAMB as `mcmc_cmb.py`. `fine_L250_s1.json` and the starting `.json` files are included.
- Other supernova samples (optional): `SNSET=union3 python mcmc_orologio.py fine_L250_s1.json 10 3000 16`

## Graph experiments (`grafo/`)

Exploratory toy models on graphs without coordinates (see section 6 of the article).

| Script | Experiment |
|---|---|
| `percolazione.py`, `esperimento1.py` | Growth of a causal set with source and sink |
| `esperimento2.py` | Gas on a ring: expansion, halt and collapse |
| `esperimento3.py` | Triangle network: does created space stay local? |
| `esperimento4.py` | Matter and events with conserved energy: singularity to singularity |
| `esperimento5.py` | Decay, unstable singularity and critical mass: cycles |
| `esperimento6.py` | Maximum signal speed and isotropy; horizon from space creation |
| `esperimento7.py` | Moving clock: emergent time dilation |
| `esperimento8.py` | Length contraction of an object held together by signals (radar rule) |

Example for experiment 8:

- `python esperimento8.py 300000 200 1 0.12 0.2,0.4,0.6,0.8` takes about 5 minutes.
- For a larger object, `python esperimento8.py 1200000 150 1 0.12 0.4,0.6,0.8` takes about 15–20 minutes.

## Included result files

| File | Content |
|---|---|
| `cmb2_LCDM_U*.json` | ΛCDM reference fits with full Planck, for each supernova sample |
| `fine_L250_s1.json` | Event clock from the fine-resolution simulation (250 Mpc box) |
| `cmb_orologio_fine_L250_s1_10.json` | Full-Planck fit with that clock |
| `nucleo_lente.npz` | Lensing kernel used by `cmb_full.py` |

## Glossary of Italian names

| Italian | English | Italian | English |
|---|---|---|---|
| orologio | clock | soglia | threshold |
| sorgente | source | pozzo | sink |
| materia | matter | spazio | space |
| evento, eventi | event, events | collasso | collapse |
| collassata | collapsed | inversione | turnaround |
| catena | chain | passi | steps |
| camminatori | walkers | processi | processes |
| corsa | run | riassunto | summary |
| adatta, adattamento | fit (verb), fit (noun) | analizza | analyse |
| scarica_dati | download data | crescita | growth |
| controllo | check | ambiente | environment |
| futuro | future | scala fine | fine scale |
| universo | universe | traccia | trace |
| versione locale | local version | cono | (light) cone |
| gusci | shells | orologi | clocks |
| grafo | graph | esperimento | experiment |
| percolazione | percolation | rete | network |
| vertice, vicini | vertex, neighbours | lato, legame | edge, link |
| seme | random seed | stelle, stellare | stars, stellar |
| originale | original (sink form) | diluito | diluted |
| ritmo | rate | massa | mass |
| oggi | today | fine | end |
| alfa | alpha | ripresa | restart |

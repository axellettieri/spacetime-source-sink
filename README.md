# Space created by events instead of dark energy

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24250555.svg)](https://doi.org/10.5281/zenodo.24250555)

**A source–sink cosmological model, tested against Planck, DESI DR2 and supernova data**

Pasquale Axel Lettieri · Independent researcher, Novedrate (CO), Italy

## The idea

Spacetime is treated as a causal graph of events. Every event that dissipates energy creates space, and masses reabsorb it. What we call dark energy is then not a substance. It is the balance of space created and not yet reabsorbed.

The model turns this idea into two terms inside general relativity:

- a **source**, proportional to the rate of events (star formation, or matter collapsing into halos)
- a **sink**, proportional to the absorbing mass

## Main results

- With full Planck 2018, DESI DR2 and Pantheon+, the model improves the ΛCDM χ² by 8.3–8.8 points. That is as much as the w0wa parametrization, with one parameter fewer.
- The nominal preference is 2.9–3.6 σ, depending on the supernova sample (Pantheon+, Union3, DES-Dovekie).
- Two independent event clocks give the same fit: the observed star-formation history, and a clock computed from N-body simulations and Press–Schechter theory.
- **Testable prediction.** The model is forced onto a line in the (w0, wa) plane: wa ≈ −0.41 − 1.61 (1 + w0).
- **Future.** If absorption continues, expansion stops in about 30 billion years and the universe recollapses.
- **Negative results.** The model does not explain the Hubble tension, and its local version is excluded.
- **Graph experiments.** Exploratory two-dimensional toys show expansion–collapse cycles with conserved energy, and emergent time dilation and length contraction.

## Contents

| Folder | Content |
|---|---|
| `paper/article_en.md` | The article (English) |
| `paper/articolo_it.md` | The article (Italian, original version) |
| `code/` | All scripts: data fits, full Planck with CAMB, MCMC, GPU N-body simulations, graph experiments (`code/grafo/`) |

How to run the code is described in [`code/README.md`](code/README.md) (English) and `code/LEGGIMI.txt` (Italian). Comments and screen output are in English; file and variable names are in Italian, with a glossary in `code/README.md`.

Quick start:

```
pip install -r requirements.txt
cd code
python scarica_dati.py    # downloads public data, about 210 MB
python fit3.py            # main fit: supernovae + DESI + compressed Planck, about 30 s
```

## Status

This is a draft. It has not been peer reviewed, and the calculations have not yet been verified independently.

## Use of AI

The analysis code and the text were developed with the assistance of an AI model (Claude, Anthropic). The hypothesis, the research direction and all choices are the author's. The author takes responsibility for the content. Several errors were found and corrected during the work, and others may remain.

## License

- Code: [MIT](LICENSE)
- Article text and figures (`paper/`): [CC BY 4.0](LICENSE-CC-BY-4.0.md)

## How to cite

Lettieri, P. A. (2026). *Space created by events instead of dark energy: a source–sink cosmological model confronted with the data*. Zenodo. https://doi.org/10.5281/zenodo.24250555 (all versions; v1.0: https://doi.org/10.5281/zenodo.24250556)

See also [CITATION.cff](CITATION.cff).

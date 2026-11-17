# Space created by events instead of dark energy: a source–sink model confronted with the data

**Pasquale Axel Lettieri** ([ORCID 0009-0001-0756-6467](https://orcid.org/0009-0001-0756-6467)) · Independent researcher, Novedrate (CO), Italy · 15 November 2026 · draft

## Abstract

A model that explains the accelerated expansion without dark energy as a separate substance: in its place there is space created by matter events and reabsorbed by masses. The model fits the cosmological data better than the standard ΛCDM model, with a nominal preference between 2.9 and 3.6 sigma depending on the supernova sample. It has one parameter more than ΛCDM and one fewer than the w0wa parametrization, which it matches.

The event rate was taken in two independent ways: from the observed star-formation history, and from the gravitational collapse of matter, computed with N-body simulations and with Press–Schechter theory. The two clocks give the same fit to Planck 2018, DESI DR2 and three supernova samples (Pantheon+, Union3, DES-Dovekie).

If masses keep absorbing space after the space created by recent events is used up, the model predicts that the expansion stops in about 30 billion years and that the universe recollapses in about 45–50 billion years. The model does not explain the Hubble tension, and its local version is excluded by the data.

A series of exploratory numerical experiments on graphs shows that a local creation rule can produce a diffuse effect; that a closed energy budget yields cycles of expansion and collapse; and that a moving object made of internal signals slows its own time and contracts, as in relativity, with deviations that shrink as the object grows relative to the graph's grain. These experiments are two-dimensional toys and are not a derivation of the model.

The model makes numerical predictions for the full DESI data and for Euclid, deposited before their release (section 4.10). The analysis has not been independently verified.

## 1. Introduction

The starting hypothesis is that spacetime is not a given container but the product of physical events: every event that dissipates energy creates space, and masses reabsorb it. In this reading dark energy is not a substance: it is the balance of space created and not yet reabsorbed.

The hypothesis arose qualitatively from the picture of the universe as a causal graph, whose nodes are events and whose links are cause–effect relations (Bombelli, Lee, Meyer and Sorkin, 1987). This work translates it into a minimal cosmological model, testable against data, within general relativity.

The observational context favours models of this kind. In 2025 the DESI collaboration reported a preference for time-varying dark energy, between 2.8 and 4.2 sigma depending on the supernova sample, described with the w0wa parametrization. The form preferred by the data is a dark energy that first grows and then declines. A model with a source that runs out and an absorption that continues produces exactly this form.

There were three initial goals: explain the accelerated expansion, explain the Hubble tension through locally non-uniform expansion, and obtain a recollapsing universe. The first and third are supported, under the conditions discussed in the text; the second is not.

## 2. The model

The model is an interacting vacuum that exchanges energy with matter, in a flat universe described by general relativity. The created space has energy density ρx and behaves as vacuum energy; matter has density ρm. The exchange Q conserves total energy:

```latex
\dot\rho_x = Q, \qquad \dot\rho_m + 3H\rho_m = -Q, \qquad H^2 = \frac{8\pi G}{3}\,(\rho_m+\rho_r+\rho_x)
```

The exchange is the difference between a source, proportional to the event rate, and a sink, proportional to the absorbing mass:

```latex
Q = \alpha\,\psi(z)\,(1+z)^3 \;-\; \gamma\,M(z), \qquad M(z)=\int \psi\,dt
```

Here ψ is the event rate per comoving volume, α the space created per event and γ the space absorbed per unit mass and time. Initially ρx is zero: all of it is created by events. There are two free parameters for this component (α and γ), against one for ΛCDM and three for w0wa.

**Two event clocks.** The rate ψ was taken in two independent ways.

| Clock | What counts as an event | Where it comes from |
|---|---|---|
| Stellar | Star formation | Observed history (Madau and Dickinson 2014) |
| Emergent | Matter entering a collapsed region for the first time | N-body simulations of this work and Press–Schechter theory |

With the emergent clock, the source is proportional to the matter collapsing per unit time and the sink to the fraction of matter already collapsed.

**Variants tested.** A source not diluted by the expansion is excluded by the data. Three forms of the sink were tested (proportional to the created space, constant, proportional to the accumulated mass): they fit today's data equally well but lead to different fates. The form proportional to the created space gives eternal expansion; the other two lead to collapse.

**Perturbations.** The exchange is treated as diffuse: matter loses or gains energy uniformly, with no momentum exchange in its own rest frame. Structure growth follows the corresponding equation; in the cosmic microwave background spectrum the model enters as an effective fluid, with a correction to the lensing amplitude computed from the growth equations (about 0.98).

## 3. Data and method

The model is compared with three families of data: supernovae, baryon acoustic oscillations and the cosmic microwave background.

| Data | Detail |
|---|---|
| Pantheon+ supernovae | 1590 supernova light curves with the official selection (redshift above 0.01), full covariance matrix |
| Union3 and DES-Dovekie supernovae | Alternative samples, used as a check |
| DESI DR2 | 13 baryon acoustic oscillation measurements |
| Planck 2018 | Full TT, TE, EE spectrum (plik-lite likelihood) with a prior on the optical depth; alternatively the compressed distance priors for quick comparisons |
| Structure growth | 7 fσ8 measurements (6dFGS, SDSS MGS, BOSS DR12, eBOSS DR16), from the compilation of Alam et al. 2021, used as a check |

**Two levels of comparison.** The quick comparison uses the compressed Planck priors and takes seconds; it is used to explore variants. The full comparison computes the CMB spectrum with CAMB and uses the sound horizon computed by CAMB; the main results refer to it.

**Statistics.** For each model the χ² minimum is found and compared with that of ΛCDM on the same data. The difference is converted into nominal sigma for the number of extra parameters. The models are not exactly nested, so the sigma values are indicative. For the stellar-clock model a full Markov chain (32 walkers, 3000 steps) was also run to derive the uncertainties on the future.

**Simulations.** The emergent clock comes from particle-mesh N-body simulations: 512³ particles, boxes from 125 to 1000 Mpc, on a GPU. Without events, large-scale structure growth reproduces linear theory to within 2%.

## 4. Results

### 4.1 Main comparison

With full Planck, DESI DR2 and Pantheon+, the source–sink model improves the ΛCDM χ² by 8.3–8.8 points, as much as w0wa, with one parameter fewer.

| Model | Extra parameters vs ΛCDM | χ² | Difference from ΛCDM |
|---|---|---|---|
| ΛCDM | 0 | 2008.70 | 0 |
| w0wa | 2 | 2000.01 | −8.7 |
| Source–sink, stellar clock, sink on accumulated mass | 1 | 2000.36 | −8.3 |
| Source–sink, stellar clock, sink on created space | 1 | 1999.94 | −8.8 |
| Source–sink, emergent clock (simulation) | 1, plus the clock threshold | 2000.03 | −8.7 |
| Source–sink, emergent clock (theory) | 1, plus the minimum mass | 2000.31 | −8.4 |

The model's expansion history matches that of the best-fit w0wa to within 0.2% up to redshift 5. The resulting Hubble constant is 67.7 km/s/Mpc.

### 4.2 Both ingredients are needed

In the quick comparison the source alone worsens the fit (+9.1) and the sink alone improves it only slightly (from −1.4 to −2.1); together they give −7.3 to −7.7. The exact shape of the clock is not decisive: a generic source that declines with expansion gives −6.7.

### 4.3 The emergent clock

The clock computed from the collapse of matter works as well as the observed one, provided events start early: between 7 and 11% of matter collapsed at redshift 3 and between 55 and 61% today. With Press–Schechter theory this corresponds to halos above 10^10–10^11 solar masses, i.e. those hosting galaxies; with minimum masses of 10^12 or more the fit is lost. In the simulations the result holds at all resolutions and seeds tested, but the density threshold that produces it changes with resolution: the particle-mesh method is not enough to fix the minimum mass.

### 4.4 Other supernova samples

| Supernovae | Stellar clock | Emergent clock | Nominal sigma |
|---|---|---|---|
| Pantheon+ | −8.3 | −8.7 | 2.9 |
| Union3 | −13.0 | −12.7 | 3.6 |
| DES-Dovekie | −10.9 | −11.0 | 3.3 |

Counting the emergent-clock threshold as a parameter too, the sigma values drop to 2.5, 3.1 and 2.9.

### 4.5 Negative results

- **Hubble tension.** Not explained: the model gives H0 = 67.7, like ΛCDM.
- **Faster local expansion.** The velocity data of nearby galaxies show no sign of faster local expansion, to within 1%; the 8% needed for the tension is excluded at more than 15 sigma. No environmental dependence (−0.6 ± 0.6%), measured with the 2M++ density and velocity map (Carrick et al. 2015).
- **Local version in simulation.** If created space stays near the events, local expansion increases, but by at most 3.6 percentage points within 75 Mpc and 0.2 within 300 Mpc: at least ten times less than needed.
- **Local exchange.** Local absorption of space by masses is excluded; local return of energy to all matter is excluded by the limits on the Sun's mass variation; return to dark matter alone is disfavoured at about 3 sigma.
- **Average of events instead of accumulation.** A component proportional to the average event density in the visible past worsens the fit by 104 points.

### 4.6 Energy budget

The source converts about 9% of matter into space, about 4.5 kg for every kg of stars formed. This is about 3000 times the energy that stars dissipate: the event acts as a trigger, not as the energy source, and the model does not explain this factor.

### 4.7 Matter, gravitational lensing and structure growth

The model predicts that matter today is 39–40% of the total, against 30% in ΛCDM, but this does not conflict with independent lensing and growth measurements. The matter returned by the sink is distributed uniformly: it raises the mean density and dilutes the contrasts in the same proportion, so their product, which is what lensing measures, stays that of ΛCDM.

| Redshift | Matter relative to ΛCDM | Contrast relative to ΛCDM | Lensing source relative to ΛCDM |
|---|---|---|---|
| 0 | 1.31 | 0.78 | 1.03 |
| 0.2 | 1.10 | 0.92 | 1.01 |
| 0.5 | 0.99 | 1.01 | 1.00 |
| 1 | 0.97 | 1.04 | 1.01 |
| 2 | 0.99 | 1.03 | 1.01 |

Values for the stellar clock; with the emergent clock they differ by less than 0.5% in the last column. Weighted over the typical weak-lensing distances, the signal amplitude is 1.004 times that of ΛCDM (w0wa: 1.017).

| Measurement | Observed value | Where the model stands |
|---|---|---|
| Weak lensing, KiDS-Legacy ([arXiv:2503.19441](https://arxiv.org/abs/2503.19441)) | S8 = 0.815 (+0.016, −0.021) | Like ΛCDM with Planck parameters (about 0.81–0.82): consistent |
| Weak lensing and clustering, DES Year 6 ([arXiv:2601.14559](https://arxiv.org/abs/2601.14559)) | S8 = 0.789 ± 0.012 | Like ΛCDM: about 2 sigma above, the same mild tension as the standard model |
| Growth rate, 7 fσ8 measurements |  | χ² = 5.8 against 6.3 for ΛCDM |

Two consequences. First, the check is passed: the higher matter fraction is not a problem. Second, the model is nearly indistinguishable from w0wa in these measurements too: lensing, galaxy clustering and growth rate all see the same combination of matter and contrast. Cluster counts were not computed; they are expected to change little, because the added matter is diffuse and recent, but this remains to be checked.

### 4.8 A testable prediction: the line in the (w0, wa) plane

In plain words: dark-energy measurements are published as two numbers, w0 (what it is like today) and wa (how different it was in the past). Here the model is not replaced by w0wa, but translated into those two numbers to compare it with published measurements. Translated, it says that the effect grew in the past and is declining today.

The model has a single free parameter more than ΛCDM, the sink strength. As it varies, the (w0, wa) point moves along a line, like a train on its track, while w0wa can occupy the whole plane. If future measurements fall off that line, the model is ruled out.

Since the model has one parameter fewer than w0wa, its points lie on a line. For each sink strength the (w0, wa) pair reproducing the model's expansion history was found (to within 0.5% up to redshift 2.5). With the stellar clock and the sink on accumulated mass, the line is:

```latex
w_a \simeq -0.41 \;-\; 1.61\,(1+w_0)
```

With the other two sink forms the slope is −1.46 (constant sink) and −1.23 (sink on created space); the starting point is the same. The line does not pass through ΛCDM (w0 = −1, wa = 0): without the sink the model gives w0 = −0.97 and wa = −0.46, and is excluded by the data.

| Supernovae | Model: w0, wa | Free w0wa on the same data: w0, wa |
|---|---|---|
| Pantheon+ | −0.82, −0.70 | −0.86, −0.51 |
| DES-Dovekie | −0.81, −0.72 | −0.83, −0.63 |
| DES Year 5 (earlier version) | −0.76, −0.79 | −0.77, −0.77 |
| Union3 | −0.72, −0.85 | −0.70, −0.98 |

Quick comparison (compressed Planck). In all four samples the point preferred by free w0wa lies within 0.2 in wa of the model's line, i.e. within the current uncertainty on wa (about ±0.2). With Pantheon+ the model fixes wa between −0.75 and −0.65, a range about four times narrower than that of free w0wa.

This is the prediction that upcoming data can falsify: when the uncertainty on wa drops below 0.1, the measured point must lie on the line. The direction along which current data are uncertain has a different slope (about −2.9), so the two lines diverge away from their crossing point, near w0 = −0.78.

### 4.9 Uncertainties with the emergent clock

A full Markov chain (emcee, 32 walkers, 7 parameters, 5040 steps) with the emergent clock from the fine-resolution simulation (250 Mpc box, threshold 10), the full Planck spectrum, DESI DR2 and Pantheon+ gives these intervals:

| Quantity | Central value | 68% interval | 95% interval |
|---|---|---|---|
| χ² improvement over ΛCDM | −8.48 | | |
| H0 (km/s/Mpc) | 67.9 | 67.3 – 68.6 | 66.7 – 69.4 |
| Matter today (fraction) | 0.38 | 0.34 – 0.41 | 0.30 – 0.44 |
| Equivalent w0 | −0.84 | −0.88 – −0.80 | −0.94 – −0.76 |
| Equivalent wa | −0.61 | −0.67 – −0.54 | −0.73 – −0.46 |
| ln γ (sink strength) | 0.32 | −0.10 – 0.58 | −1.33 – 0.78 |

The model fixes wa with an uncertainty of about ±0.065, against about ±0.2 for free w0wa on the same data. The central point lies close to the line of section 4.8, which for w0 = −0.84 gives wa ≈ −0.67. Extending the chain from 3000 to 5040 steps changed the values by less than 0.01: the results are stable. The kept steps are 24 times the autocorrelation time, slightly fewer than the recommended 30–50; the central values are solid, the ends of the 95% intervals should be taken as indicative.

### 4.10 Predictions for upcoming data: full DESI and Euclid

The following predictions were computed on 17 November 2026, before the release of the cosmological results of the full five-year DESI survey and of Euclid DR1 (expected between late 2026 and 2027), and deposited with a certified date on Zenodo. They come from 400 parameter combinations drawn from the chain of section 4.9; for comparison, the ΛCDM model that best describes the same current data is shown.

**The (w0, wa) point.** The model predicts w0 = −0.84 (−0.88 – −0.80) and wa = −0.61 (−0.67 – −0.54), on the line of section 4.8. A combined measurement that places wa outside this interval by more than twice its own uncertainty, or that brings the point back to ΛCDM (w0 = −1, wa = 0), rules the model out.

**Baryon acoustic oscillation distances (DESI).** Ratios between distance and sound horizon r_d at the redshifts of the DESI samples. 68% intervals.

| Sample | Quantity | Model | ΛCDM | Model − ΛCDM | DESI DR2 measured |
|---|---|---|---|---|---|
| BGS (z = 0.295) | DV/r_d | 7.915 (7.878 – 7.954) | 7.949 | −0.43% | 7.942 ± 0.076 |
| LRG1 (z = 0.510) | DM/r_d | 13.271 (13.214 – 13.327) | 13.331 | −0.45% | 13.588 ± 0.168 |
| LRG1 (z = 0.510) | DH/r_d | 22.400 (22.321 – 22.479) | 22.559 | −0.70% | 21.863 ± 0.429 |
| LRG2 (z = 0.706) | DM/r_d | 17.420 (17.361 – 17.479) | 17.504 | −0.48% | 17.351 ± 0.180 |
| LRG2 (z = 0.706) | DH/r_d | 19.979 (19.913 – 20.045) | 20.061 | −0.41% | 19.455 ± 0.334 |
| LRG3+ELG1 (z = 0.934) | DM/r_d | 21.688 (21.627 – 21.751) | 21.780 | −0.42% | 21.576 ± 0.162 |
| LRG3+ELG1 (z = 0.934) | DH/r_d | 17.514 (17.467 – 17.563) | 17.515 | −0.01% | 17.641 ± 0.201 |
| ELG2 (z = 1.321) | DM/r_d | 27.781 (27.717 – 27.844) | 27.856 | −0.27% | 27.601 ± 0.325 |
| ELG2 (z = 1.321) | DH/r_d | 14.130 (14.107 – 14.155) | 14.059 | +0.51% | 14.176 ± 0.225 |
| QSO (z = 1.484) | DM/r_d | 29.987 (29.923 – 30.053) | 30.050 | −0.21% | 30.512 ± 0.764 |
| QSO (z = 1.484) | DH/r_d | 12.967 (12.950 – 12.985) | 12.884 | +0.65% | 12.817 ± 0.518 |
| Lya (z = 2.330) | DM/r_d | 38.980 (38.916 – 39.042) | 38.971 | +0.02% | 38.989 ± 0.532 |
| Lya (z = 2.330) | DH/r_d | 8.709 (8.701 – 8.718) | 8.637 | +0.83% | 8.632 ± 0.101 |

The model predicts shorter distances than ΛCDM up to z ≈ 1 (between −0.3 and −0.7%) and a longer radial distance DH at high redshift (+0.5% at z = 1.3, +0.8% at z = 2.3, in the Lyman-α forest of quasars).

**Growth rate of structure, fσ8 (DESI and Euclid).**

| Redshift | Model | ΛCDM | Model − ΛCDM |
|---|---|---|---|
| 0.29 | 0.464 (0.457 – 0.471) | 0.464 | −0.2% |
| 0.51 | 0.469 (0.462 – 0.476) | 0.468 | +0.2% |
| 0.71 | 0.461 (0.454 – 0.468) | 0.458 | +0.7% |
| 0.93 | 0.442 (0.436 – 0.449) | 0.437 | +1.2% |
| 1.00 | 0.436 (0.430 – 0.442) | 0.430 | +1.4% |
| 1.20 | 0.415 (0.409 – 0.421) | 0.408 | +1.8% |
| 1.40 | 0.393 (0.388 – 0.399) | 0.385 | +2.0% |
| 1.65 | 0.367 (0.362 – 0.371) | 0.359 | +2.1% |

The model predicts the same growth as ΛCDM up to z ≈ 0.5 and increasingly higher growth beyond, up to +2% at z = 1.65, in the range covered by Euclid.

**Gravitational lensing.** The product σ8·Ωm, which sets the amplitude of the weak-lensing signal, is 0.251 (0.245 – 0.258) against 0.244 for ΛCDM (+3%). The amplitude σ8 of matter alone, instead, is much lower: 0.67 (0.62 – 0.73) against 0.81, because the matter returned by the sink is diffuse and dilutes the contrasts. The S8 value published by surveys is derived assuming ΛCDM, so it should not be compared directly with this σ8.

**What can decide.** The differences from ΛCDM are small: less than 1% in the distances and up to 2% in the growth, comparable to the error of a single measured point. The decision will come from the combination of all points, i.e. the (w0, wa) point, and above all from full DESI and Euclid together. The final DESI data also include those of DR2, already used to fit the model: the fully independent test is Euclid, together with the new part of the DESI data.

## 5. The future of the expansion

With the sink proportional to the accumulated mass, the model predicts that the expansion stops in about 30 billion years and that the universe recollapses. Two independent clocks give the same dates.

|  | Stellar clock (full chain) | Emergent clock (simulation) |
|---|---|---|
| Created space used up | in 17 billion years (13–26) | in 16 |
| Expansion stops | in 31 (26–41) | in 29 |
| Maximum size | 3.2 times today's (2.6–4.7) | 3.0 |
| Collapse | in 49 (42–62) | in 45, end of simulation at one tenth of today's size |

68% intervals in parentheses. These values replace the rougher ones from the quick comparison.

**The prediction is conditional.** The expansion stops because the energy of created space drops below zero: masses keep absorbing after the space created by the events counted in the model is used up. If absorption stopped at zero, a matter-only universe expanding forever would remain. Today's data do not distinguish the two cases. In the starting hypothesis all space, including that of the primordial expansion, is born from events, and absorption beyond zero is natural; but the model does not count primordial events, and this part remains an assumption.

**The simulated collapse.** In the N-body simulation carried past the turnaround, 99.9% of matter has collapsed by the end. Matter grows by about 80 times during the cycle, because the sink keeps returning energy to matter while the space deficit grows. No bounce appears: the model describes the collapse, not the restart.

## 6. Exploratory experiments on graphs

The model's equations are not derived from the causal graph. Eight numerical experiments explore whether the qualitative ideas of the hypothesis hold on graphs without coordinates. They are toys, almost all two-dimensional, with hand-chosen rules: they show that a behaviour is possible, not that it is the real one.

| Experiment | Question | Result |
|---|---|---|
| 1. Growth of a causal set with source and sink | What separates events? | The graph width follows the link probability; the graph collapses to a chain only if the sink is present |
| 2. Gas on a ring | Does the rule "events create space, masses absorb it" produce expansion? | Expansion, halt and collapse without imposed equations; non-uniform expansion that averages out over large regions |
| 3. Triangle network without coordinates | Does space created at a point stay there? | Without relaxation new vertices create no distance; with weak relaxation space stays local; with strong relaxation it spreads everywhere (93%) |
| 4. Matter and events on the network, energy conserved | Are the two singularities equal? | From 4 vertices to about 11,000 and back to 4; identical mass and energy at start and end, but 300 particles at the start and a single object at the end |
| 5. Decay | What restarts a cycle? | Unstable singularity: cycles with the same maximum within 2%. Critical mass: cycles without singularity, maxima varying by 15%. Decay of every hot object: irregular equilibrium, no collapse |
| 6. Signals on the network | Is the maximum speed the same in every direction? What does space creation do? | Regular network: 10% difference between directions; random network: 2.4%, consistent with noise. Space creation delays signals and produces a horizon |
| 7. Moving clock | Does a moving object perform fewer internal events? | Yes: the slowdown approaches the relativistic factor as the object grows (2.03 against 2.29 at nine tenths of the maximum speed) |
| 8. Object held together by signals | Does a moving object contract? | Yes, if its length along the motion is set by signals rather than imposed. At eight tenths of the maximum speed the factor is 0.631 for a 97-link object, against 0.600 in relativity; the deviation shrinks with size (12.5%, 9.0%, 5.1% for 25, 49, 97 links) |

Three useful indications for the hypothesis.

- **Local rule, diffuse effect.** The data exclude a component concentrated near events. On the graph, creation can be local and the effect diffuse, if the geometry relaxes much faster than space is created. The experiment does not say why this should be so.
- **Cycles.** With a closed energy budget the universe returns to the minimal state with the same mass and the same energy. To restart, an instability is needed that triggers only in the extreme state, or when a critical mass is exceeded.
- **Time and motion.** Even on a graph with a single global step, the time of an object made of internal signals slows down when the object moves. Length contraction also appears, if the object's length along the motion is set by signals rather than imposed (experiment 8). Both effects approach the relativistic ones the larger the object is compared with the graph's grain. The graph, however, has a preferred frame, its own: a theorem by Bombelli, Henson and Sorkin shows that a network with a finite number of neighbours per node cannot represent a discrete spacetime without picking one, whereas causal sets sprinkled at random into spacetime avoid this.

## 7. Discussion

### 7.1 Strengths

- **Few parameters with a physical meaning.** Space created per event and space absorbed per mass, against the two descriptive parameters of w0wa.
- **Two independent clocks.** Observed star formation and computed gravitational collapse give the same fit, on three supernova samples.
- **A distinguishing prediction.** A created-space energy that becomes negative and an expansion that stops, unlike ΛCDM.

### 7.2 Limitations

- **The signal is not new.** The statistical preference is largely the same that DESI reports for w0wa. The model offers a physical interpretation of it; it does not discover it, nor strengthen it beyond saving one parameter.
- **The significance is below the discovery threshold**, and depends on the supernova sample.
- **The mechanism is missing.** There is no theory saying why an event should create space, nor justifying a created energy about 3000 times larger than the dissipated one.
- **The emergent clock has a parameter chosen on the data** (the collapse threshold or minimum mass), although its value falls in the physically expected range.
- **Perturbations are treated approximately.** The effective fluid and the lensing correction hold for the diffuse exchange; the integrated Sachs–Wolfe effect, CMB lensing reconstruction and cluster counts have not been fully compared with data.
- **The link with the causal graph is qualitative.** The experiments of section 6 do not derive the equations of section 2.
- **Few distinguishing predictions.** Expansion, lensing and growth coincide with those of w0wa to within 1–2%. What distinguishes the model is the forced line in the (w0, wa) plane of section 4.8, the link with the history of structures, and the future, which is not observable.
- **No independent verification.** Analysis and code were produced in work sessions assisted by an artificial intelligence model; several errors were found and corrected along the way, and others may remain.

### 7.3 Relation to the literature

| Work | Relation to this model |
|---|---|
| DESI DR2 (arXiv:2503.14738) | Source of the preference for evolving dark energy; our w0wa values reproduce theirs |
| Interacting vacuum (broad literature) | The model is a special case, with the exchange tied to an astrophysical clock |
| Croker et al., cosmologically coupled black holes (arXiv:2405.12282) | Same idea of a source tied to star formation, without a sink; our source-only variant reproduces their H0 near 70 |
| Li and Zhang, cosmic sign reversal (arXiv:2506.18477) | They reconstruct from DESI DR2 data, without assuming a form, the energy exchange between dark matter and dark energy: the exchange reverses direction over time, from matter to dark energy in the past and the opposite today. This is the behaviour of our model, in which the source prevails first and the sink later |
| Luu, Qiu and Tye, lifespan of the universe (arXiv:2506.24011) | Collapse predicted by a different model (a light field plus a negative cosmological constant), with a total age of about 33 billion years |
| Sorkin, everpresent cosmological constant | Dark energy from causal sets, from random fluctuations instead of matter events |
| Bombelli, Henson and Sorkin, discreteness without symmetry breaking (arXiv:gr-qc/0605006) | Explains why the networks of section 6, with few neighbours per node, necessarily have a preferred frame, and why causal sets sprinkled at random into spacetime do not |

Proposals of a varying ratio between space and time in the early universe (varying speed of light, Hořava–Lifshitz gravity, causal dynamical triangulations) are close to the qualitative ideas discussed in section 6, but were not used in the calculations.

## 8. Conclusions and future work

Space created by matter events and reabsorbed by masses, instead of dark energy as a separate substance, describes current cosmological data as well as the best parametrization in use, with one parameter fewer and a physical reading. The result holds when the observed clock is replaced by a computed one. The preference over the standard model is between 2.9 and 3.6 nominal sigma: a hint, not a discovery.

The starting hypothesis comes out partly confirmed and partly refuted. The confirmed part is accumulation: created space persists, and the data accept this form while rejecting the alternatives tested. The refuted part is the locality of the effect: expansion is not faster where there are more events, and the Hubble tension remains unexplained.

Future work, in order of priority:

1. Independent verification of the main calculations by a cosmologist.
1. Comparison of the predictions of section 4.10 with the full DESI data and with Euclid, as soon as they are released.
1. Comparison with the integrated Sachs–Wolfe effect and with CMB lensing reconstruction.
1. Emergent clock with a real halo finder, to fix the minimum mass without looking at the data.
1. A causal set growing event by event in more dimensions, with matter made of chains of events: the step that would link the equations to the graph, and an open problem in research too.

## Glossary

| Term | Meaning |
|---|---|
| Dark energy | Whatever makes the expansion of the universe accelerate. Its nature is unknown; only its effect is measured. |
| ΛCDM | The standard model of cosmology: dark energy constant in time (Λ) plus cold dark matter (CDM). |
| w0wa | Not a theory but a description: it lets dark energy change in time and uses two numbers to say how. w0 says what it is like today, wa how much it changes going into the past. With w0 = −1 and wa = 0 one recovers ΛCDM. |
| Free parameter | A number the theory does not fix and that is derived from the data. For an equal fit, fewer parameters is better. |
| χ² (chi-squared) | A measure of how far a model is from the data: the lower, the better the fit. What matters here is the difference from ΛCDM. |
| Sigma | How unlikely it is that a difference is due to chance. At 3 sigma one speaks of a hint, at 5 of a discovery. |
| Redshift (z) | A measure of how far back in time what we observe is: z = 0 is today, z = 1 about 8 billion years ago, z = 3 about 11.5. |
| Source | In the model, the creation of space by events. |
| Sink | In the model, the reabsorption of space by masses. |
| Created space | The space produced by events and not yet reabsorbed. In the model it takes the place of dark energy: it produces the same effect without being a separate substance. |
| Event | A physical occurrence that dissipates energy. In the calculations: the birth of stars, or matter falling into an already dense region. |
| Event clock | The curve giving how many events occur at each epoch. Stellar: from observations of stars. Emergent: computed from the clustering of matter. |
| Diluted, undiluted source | Two ways of counting events as the universe expands. Diluted: counts events per physical volume, which grows with the expansion. Undiluted: counts them per volume that follows the expansion. The data accept only the first. |
| Interacting vacuum | The mathematical form of the model: created space takes energy from matter when it is born and returns it when reabsorbed. |
| Diffuse, local exchange | Diffuse: energy is taken from and returned to matter uniformly everywhere. Local: only near events or masses. The data exclude the local version. |
| Collapsed matter, halo | Matter that gravity has gathered into a stable clump. A halo is one such clump; large ones host galaxies. |
| Collapse threshold | In the simulations, the density above which a region is considered collapsed. |
| N-body simulation | A calculation following millions of particles attracting each other by gravity, to see how structures form. |
| Hubble constant (H0) | The expansion rate of the universe today. |
| Hubble tension | The disagreement between H0 measured nearby (about 73) and that inferred from the early universe (about 67.5). |
| Planck | Satellite that measured the oldest light in the universe (cosmic microwave background). |
| DESI | Instrument that measures the distance of millions of galaxies using an imprint left by sound in the early universe (baryon acoustic oscillations). |
| Supernovae | Stellar explosions of known luminosity, used to measure distances. Pantheon+, Union3 and DES are three different compilations. |
| Weak gravitational lensing, S8 | Matter slightly distorts the images of distant galaxies. S8 is the number summarizing how much matter there is and how clustered it is. |
| Contrast | How much denser than average a region is. |
| Mpc (megaparsec) | Unit of distance: about 3.3 million light years. |
| Turnaround, collapse | Turnaround: the moment the expansion stops. Collapse: the contraction that follows. |
| Causal graph, causal set | Representation of the universe as a network: nodes are events, links say which event can influence which. |
| Vertex, link | On the experiment graphs: a vertex is an atom of space, a link connects two neighbouring vertices. |
| Relaxation | In the graph experiments, the network settling into a regular shape after vertices have been added. |

## 9. Reproducibility and references

All results can be reproduced with the code in the archive modello_sorgente_pozzo_v18.zip (`code/` folder); graph experiments 4 to 7 were rerun on a second machine with identical numbers. The code is in the `code/` folder of this repository; the full work log, with the corrections made along the way, is available on request. Italian version: [articolo_it.md](articolo_it.md).

| Part | Scripts |
|---|---|
| Model and quick comparison | fit3.py |
| Full Planck, stellar clock | cmb_full.py, mcmc_cmb.py |
| Full Planck, emergent clock | cmb_orologio.py |
| N-body simulations | universo3d.py, adatta3d.py, corsa3d.sh, collasso.sh, scala_fine.sh |
| Graph experiments | grafo/percolazione.py, grafo/esperimento1.py … esperimento8.py |
| Uncertainty chain, emergent clock | mcmc_orologio.py, analizza_orologio.py, catena_orologio.sh |

**References**

- Ahmed, Dodelson, Greene, Sorkin (2004), Everpresent Λ, Physical Review D 69, 103523, [arXiv:astro-ph/0209274](https://arxiv.org/abs/astro-ph/0209274).
- Alam et al., eBOSS Collaboration (2021), Completed SDSS-IV extended Baryon Oscillation Spectroscopic Survey: cosmological implications from two decades of spectroscopic surveys at the Apache Point Observatory, Physical Review D 103, 083533, [arXiv:2007.08991](https://arxiv.org/abs/2007.08991). Source of the fσ8 compilation.
- Beutler et al. (2012), The 6dF Galaxy Survey: z ≈ 0 measurements of the growth rate and σ8, Monthly Notices of the Royal Astronomical Society 423, 3430, [arXiv:1204.4725](https://arxiv.org/abs/1204.4725).
- Bombelli, Henson, Sorkin (2009), Discreteness without symmetry breaking: a theorem, Modern Physics Letters A 24, 2579, [arXiv:gr-qc/0605006](https://arxiv.org/abs/gr-qc/0605006).
- Bombelli, Lee, Meyer, Sorkin (1987), Space-time as a causal set, Physical Review Letters 59, 521.
- Carrick, Turnbull, Lavaux, Hudson (2015), Cosmological parameters from the comparison of peculiar velocities with predictions from the 2M++ density field, Monthly Notices of the Royal Astronomical Society 450, 317, [arXiv:1504.04627](https://arxiv.org/abs/1504.04627). Map used in the version distributed by [pvhub](https://github.com/KSaid-1/pvhub).
- Croker et al. (2024), DESI dark energy time evolution is recovered by cosmologically coupled black holes, JCAP 10 (2024) 094, [arXiv:2405.12282](https://arxiv.org/abs/2405.12282).
- DES Collaboration (2026), Dark Energy Survey Year 6 results: cosmological constraints from galaxy clustering and weak lensing, accepted by Physical Review D, [arXiv:2601.14559](https://arxiv.org/abs/2601.14559).
- DESI Collaboration (2025), DESI DR2 results II: measurements of baryon acoustic oscillations and cosmological constraints, Physical Review D 112, [arXiv:2503.14738](https://arxiv.org/abs/2503.14738).
- Eisenstein, Hu (1998), Baryonic features in the matter transfer function, Astrophysical Journal 496, 605, [arXiv:astro-ph/9709112](https://arxiv.org/abs/astro-ph/9709112).
- Foreman-Mackey, Hogg, Lang, Goodman (2013), emcee: the MCMC hammer, PASP 125, 306, [arXiv:1202.3665](https://arxiv.org/abs/1202.3665).
- Lewis, Challinor, Lasenby (2000), Efficient computation of CMB anisotropies in closed FRW models (CAMB), Astrophysical Journal 538, 473, [arXiv:astro-ph/9911177](https://arxiv.org/abs/astro-ph/9911177).
- Li, Zhang (2025), Cosmic sign-reversal: non-parametric reconstruction of interacting dark energy with DESI DR2, [arXiv:2506.18477](https://arxiv.org/abs/2506.18477).
- Luu, Qiu, Tye (2025), The lifespan of our universe, JCAP (2025), [arXiv:2506.24011](https://arxiv.org/abs/2506.24011).
- Madau, Dickinson (2014), Cosmic star-formation history, Annual Review of Astronomy and Astrophysics 52, 415, [arXiv:1403.0007](https://arxiv.org/abs/1403.0007).
- Planck Collaboration (2020), Planck 2018 results. VI. Cosmological parameters, Astronomy & Astrophysics 641, A6, [arXiv:1807.06209](https://arxiv.org/abs/1807.06209).
- Popovic et al. (2025), The Dark Energy Survey supernova program: a reanalysis of cosmology results and evidence for evolving dark energy with an updated Type Ia supernova calibration (DES-Dovekie), [arXiv:2511.07517](https://arxiv.org/abs/2511.07517).
- Press, Schechter (1974), Formation of galaxies and clusters of galaxies by self-similar gravitational condensation, Astrophysical Journal 187, 425.
- Prince, Dunkley (2019), Data compression in cosmology: a compressed likelihood for Planck data (planck-lite-py), Physical Review D 100, 083502, [arXiv:1909.05869](https://arxiv.org/abs/1909.05869).
- Rideout, Sorkin (2000), A classical sequential growth dynamics for causal sets, Physical Review D 61, 024002, [arXiv:gr-qc/9904062](https://arxiv.org/abs/gr-qc/9904062).
- Rubin et al. (2025), Union through UNITY: cosmology with 2,000 SNe using a unified Bayesian framework (Union3), Astrophysical Journal, [arXiv:2311.12098](https://arxiv.org/abs/2311.12098).
- Scolnic et al. (2022), The Pantheon+ analysis: the full data set and light-curve release, Astrophysical Journal 938, 113, [arXiv:2112.03863](https://arxiv.org/abs/2112.03863).
- Wright et al. (2025), KiDS-Legacy: cosmological constraints from cosmic shear with the complete Kilo-Degree Survey, Astronomy & Astrophysics 703, A158, [arXiv:2503.19441](https://arxiv.org/abs/2503.19441).

References checked on 15 November 2026 against arXiv, NASA ADS or the journal pages. The 7 fσ8 measurements are those summarized by Alam et al. (2021): 6dFGS (Beutler et al. 2012), SDSS MGS, BOSS DR12 and the three eBOSS DR16 samples (luminous red galaxies, emission-line galaxies, quasars).

# Independent SiPM charge-to-energy transfer

This experiment connects frozen blind alpha fits to deposited energy. It uses
the current primary-descendant avalanche charge in equivalent pe, excluding
dark ancestry. It does not equate charge, source energy, and deposited energy.
The cumulative material light-yield table already includes quenching.

The physics and waveform fits are unchanged: fixed 95/351 ns, 2 ns samples,
2.048 us local waveforms, adequate analog dynamic range, and unknown K <= 2.
`model_relative_residual` and `joint_conserved_nll` are the main comparison.
No ADC clipping rejection is applied. No new ROOT files are produced.

## Data and frozen evaluation

* 280 isolated optical parents fit four prespecified energy mappings.
* 280 other isolated parents select the mapping by equal-group relative MSE,
  including a unit penalty for out-of-domain energy predictions.
* 320 disjoint parents supply the existing 1600-frame confirmation design:
  1200 pairs at 8/40/200 ns, 320 single events, and 80 empty controls.
* These confirmation waveforms were already published for charge/timing;
  this is a new energy endpoint, not a newly unseen transport campaign.
* Every source-energy group retains its events. Inverse conversion outside
  the observed calibration energy domain returns an explicit invalid value.
* The primary endpoint requires both events within 10% deposited energy and
  min(10 ns, separation/4) timing, with no extra output. Secondary energy
  tolerances are 5% and 20%. The truth target threshold remains 500 pe.
* Pair comparisons use 2000 paired, stratified mother-block bootstrap draws;
  reused optical parents are not independent samples across conditions.

An additional true-charge conversion diagnostic is explicitly truth-assisted.
It is neither blind reconstruction nor a rigorous performance upper bound.
The final selection and single-event transfer criterion are in `protocol.json`.
Calibration, source, protocol, and energy-model hashes are preserved separately.

## Reproduction

From the project root, using the project Python environment:

```bash
python pileup/experiments/sipm_energy_calibration_v2_20261028/run.py collect
python pileup/experiments/sipm_energy_calibration_v2_20261028/run.py calibrate
python pileup/experiments/sipm_energy_calibration_v2_20261028/run.py predict
python pileup/experiments/sipm_energy_calibration_v2_20261028/run.py evaluate
```

The first three stages refuse to overwrite existing completed outputs. A clean
checkout requires the local native input exports and prior frozen fits; these
large inputs are not distributed with GitHub. Each stage must finish before
the next starts. `run.py` uses the actual source snapshot in this directory.
ROOT identity is traced through original hashes and metadata; this round hashes
the derived truth files rather than rehashing every original ROOT.

Local results: `jupyter/artifacts/sipm_energy_calibration_v2_20261028/`.
Local provenance and validation: `jupyter/run_records/sipm_energy_calibration_v2_20261028/`.
Local research notebook: `jupyter/codex_0007_SiPM_Energy_Calibration.ipynb`.

The first directory without `v2` is preserved as an interrupted collection
attempt: float32 ROOT source energy failed an exact equality check. Version 2
uses a 2e-7 relative / 1e-9 MeV absolute source-identity tolerance and preserves
both ROOT and manifest energy values. No physical parameter was changed.

## Completed result

The prespecified rule selected `yield_scale`. The independent single-event
energy-transfer control failed: for source energies 1.625, 2.5, 5.486 and 8 MeV,
conversion of true primary charge has deposited-energy median biases of
+13.52%, +5.91%, -1.57% and -5.95%, respectively (80/80 finite in each group).
The required absolute median bias was at most 5% in every group.

Selection data expose an edge-domain problem. The saturating-yield model has
small conditional biases, but 12/80 high-energy selection charges exceed its
inverse domain. The unit rejection penalty makes its loss 0.038214, compared
with 0.005401 for the biased scaled-yield model. Do not remove these failures or
retune on the final evaluation; wider calibration guard bands and intermediate
energy points are needed before a new independent confirmation.

With this unqualified energy calibration, strict pair recovery within 10% E is
21.33% versus 24.83% for combinations containing weak alpha, Model GLS versus
joint NLL. The paired difference is +3.50 percentage points, with mother-block
95% bootstrap interval [+0.67,+6.50]. Other pairs give 47.33% versus 49.50%.
These results diagnose combined calibration and reconstruction error, not final
energy performance. The old 20%-charge endpoint is a different criterion.

Five numbered figures cover calibration, single-event controls, the three-gap
map, paired comparison, and absolute energy-count spectra. No default fitting
method was replaced. Eight numerical tests, independent-parent identities,
truth-isolated inference, and consistency with all 3200 old blind outputs pass.

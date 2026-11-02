"""Replay the frozen numeric pipeline into new directories; original files stay intact."""
from pathlib import Path
import argparse
import multiprocessing as mp
import os
import re
import shutil
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--label', required=True, help='New directory name, letters/digits/underscore/hyphen')
    parser.add_argument('--prepare-only', action='store_true')
    args = parser.parse_args()
    if not re.fullmatch(r'[A-Za-z0-9_-]+', args.label):
        parser.error('Use a single simple directory name')
    base = Path(__file__).resolve().parent
    root = base.parents[2]
    os.environ['GAGG_RESIDUAL_SOURCES_ROOT'] = str(root)
    os.environ['OPENBLAS_NUM_THREADS'] = '1'
    os.environ['OMP_NUM_THREADS'] = '1'
    sys.path.insert(0, str(base/'source'))
    from pileup.analysis import alpha_single_residual_sources_study as study
    # This project runner is verified on Linux/WSL; fork preserves the explicit
    # output roots in workers while every calculation uses the frozen sources.
    mp.set_start_method('fork', force=True)
    study.verify()
    original_art, original_rec = study.ART, study.REC
    destination = base/'reproductions'/args.label
    art = original_art/'reproductions'/args.label
    rec = original_rec/'reproductions'/args.label
    assert not any(p.exists() for p in (destination, art, rec)), 'Choose a new label'
    for p in (destination, art, rec):
        p.mkdir(parents=True)
    names = ['template.npz','energy-model.json','response-parameters.json','mean-correction.npz',
             'noise-calibration.json','sensor-config.json','charge-map-spatial12.npz','spatial-profile.npz',
             'mothers.json','protocol.json','design-freeze.json','input-freeze.json','bin/age.so']
    for name in names:
        p = destination/name
        p.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(base/name, p)
    shutil.copytree(original_art/'photons', art/'photons')
    for name in ('voltage','diagnostics','features','pilot','replay','figures'):
        (art/name).mkdir()
    copied = [destination/name for name in names] + list((art/'photons').glob('*.npz'))
    freeze = {str(p.relative_to(root)): study.sha(p) for p in copied}
    study.save(destination/'reproduction-inputs.json', dict(parent=str(base.relative_to(root)),
        original_source=str((base/'source').relative_to(root)), copied_sha256=freeze,
        output_label=args.label, original_outputs_modified=False, start_method='fork'))
    print('Prepared reproduction:', destination, flush=True)
    if args.prepare_only:
        return
    study.EXP, study.ART, study.REC = destination, art, rec
    for stage in ('pilot','production','replay','features','analyze'):
        for path, digest in freeze.items():
            assert study.sha(root/path)==digest, path
        if stage in ('pilot','production','replay'):
            study.generate(stage)
        elif stage=='features':
            study.features()
        else:
            study.analyze()
    result = study.read(destination/'result-summary.json')
    reference = study.read(base/'result-summary.json')
    assert result==reference, 'Numerical reproduction differs; retain both results for inspection'
    study.save(destination/'reproduction-equality.json', dict(all_summary_values_equal=True))


if __name__=='__main__':
    main()

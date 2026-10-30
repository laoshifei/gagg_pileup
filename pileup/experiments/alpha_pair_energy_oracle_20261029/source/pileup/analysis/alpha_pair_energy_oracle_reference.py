"""Matched-mother isolated reference, with no extra fitting or tuning.

Isolated/native and shared/offline digitizations differ in implementation,
windows and noise draws. This reference is descriptive, not a causal recovery
contrast or a rigorous performance bound.
"""
from collections import defaultdict
import csv

import numpy as np

from .alpha_pair_energy_oracle_study import EXPERIMENT, OUTPUT, RECORD, load, lines, save, verify


def main():
    verify('reference-freeze.json')
    cases = load(OUTPUT/'inputs/cases.json')
    predictions = {(r['method'], r['index']): r for r in lines(OUTPUT/'predictions.jsonl')}
    scores = {(r['method'], r['index']): r for r in lines(OUTPUT/'scores.jsonl')}
    with np.load(OUTPUT/'inputs/evaluation.labels.npz') as f: energy = f['energy_MeV']
    single = {}
    for i, c in enumerate(cases):
        if c['fold'] == 1 and c['kind'] == 'single':
            root = c['mother_roots'][0]
            assert root not in single
            single[root] = i
    assert len(single) == 280
    rows = []
    for i, c in enumerate(cases):
        if c['fold'] != 1 or c['kind'] != 'pair': continue
        ids = [single[root] for root in c['mother_roots']]
        for j, index in enumerate(ids): assert energy[i, j] == energy[index, 0]
        isolated = all(scores['fixed_true_time', index]['frame_success_10'] for index in ids)
        shared = bool(scores['fixed_true_time', i]['frame_success_10'])
        rows.append(dict(index=i, condition=c['condition'], gap_ns=c['gap_ns'],
            isolated_both_energy10=bool(isolated), shared_fixed_time_pair_energy10=shared,
            lost_after_pairing=bool(isolated and not shared), gained_after_pairing=bool(shared and not isolated)))
    grouped = defaultdict(list)
    for row in rows: grouped[row['condition'], row['gap_ns']].append(row)
    summary = []
    for (condition, gap), rr in grouped.items():
        summary.append(dict(condition=condition, gap_ns=gap, pairs=len(rr),
            isolated_both_energy10=sum(r['isolated_both_energy10'] for r in rr),
            shared_fixed_time_pair_energy10=sum(r['shared_fixed_time_pair_energy10'] for r in rr),
            lost_after_pairing=sum(r['lost_after_pairing'] for r in rr), gained_after_pairing=sum(r['gained_after_pairing'] for r in rr)))
    totals = []
    for weak, group in [(True, 'weak'), (False, 'strong')]:
        rr = [r for r in rows if ('weak' in r['condition']) == weak]
        totals.append(dict(group=group, pairs=len(rr), **{key: sum(r[key] for r in rr)
            for key in ['isolated_both_energy10', 'shared_fixed_time_pair_energy10', 'lost_after_pairing', 'gained_after_pairing']}))
    for name, data in [('isolated-reference-frames.csv', rows), ('isolated-reference-summary.csv', summary)]:
        with (OUTPUT/name).open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(data[0]), lineterminator='\n'); writer.writeheader(); writer.writerows(data)
    result = dict(pairs=900, single_mothers=280, method='fixed_true_time for both shared pair and matched isolated mothers',
        totals=totals, summary=summary, true_energy_identity_verified=True,
        limit='Native isolated versus offline shared: different digitizer/window/noise realizations; descriptive reference, not a recovery-only causal contrast or strict upper bound')
    save(OUTPUT/'isolated-reference.json', result)
    report = load(OUTPUT/'audited-report.json'); report['isolated_reference'] = result
    save(OUTPUT/'audited-report.json', report)
    print('ORACLE_ISOLATED_REFERENCE', totals, flush=True)


if __name__ == '__main__': main()

"""Independent accounting and descriptive mother-block resampling of controls."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import numpy as np

ROOT=Path(os.environ.get("GAGG_CAUSAL_REPORT_ROOT",Path(__file__).resolve().parents[2])).resolve()
SINGLE="alpha_pair_causal_controls_20261101"
SHARED="alpha_pair_causal_shared_20261101"


def read(path):return json.loads(Path(path).read_text())


def main():
    checked=0
    for name,manifests in ((SINGLE,("source-config-freeze.json","input-freeze.json")),
                           (SHARED,("freeze.json","input-freeze.json","report-source-freeze.json","plot-source-freeze.json"))):
        exp=ROOT/"pileup/experiments"/name
        for manifest in manifests:
            data=read(exp/manifest)
            for rel,digest in data.get("files",data).items():
                if isinstance(digest,str) and len(digest)==64:
                    assert hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()==digest,rel
                    checked+=1
    single=ROOT/"jupyter/artifacts"/SINGLE
    shared=ROOT/"jupyter/artifacts"/SHARED
    for path in (single/"single").glob("*.npz"):
        with np.load(path) as f:
            assert f["event_signal_mV"].shape==(8,1024)
            assert f["analog_mV"].shape==(8,1024)
            assert all(np.isfinite(f[k]).all() for k in f.files)
    assert len(list((single/"single").glob("*.npz")))==64
    assert len(list((shared/"cases").glob("*.npz")))==20
    replay=[]
    for directory,ids,main_folder in ((single,[0,48],"single"),(shared,[1],"cases")):
        for identifier in ids:
            with np.load(directory/"pilot"/f"{identifier:04d}.npz") as pilot, np.load(directory/main_folder/f"{identifier:04d}.npz") as full:
                for key in pilot.files:
                    np.testing.assert_array_equal(pilot[key],full[key][:2])
            replay.append(dict(study=directory.name,id=identifier,first_two_repeats_exact=True))
    rows=read(shared/"rows.json");rng=np.random.default_rng(20261101731)
    intervals=[]
    for condition,gap in (("high_weak",8.),("high_weak",40.),("high_weak",200.),("high_weak",1200.),("weak_high",40.)):
        group=[r for r in rows if (r["condition"],r["gap_ns"])==(condition,gap)]
        joint=[];isolated=[]
        for row in group:
            weak=row["parents"].index(row["block"])
            with np.load(shared/"cases"/f"{row['id']:04d}.npz") as f:
                assert all(np.isfinite(f[k]).all() for k in f.files)
                joint.append(f["on_shared"][:,weak].sum(-1))
                isolated.append(f["on_isolated"][:,weak].sum(-1))
        joint=np.array(joint);isolated=np.array(isolated)
        observed=100*(joint.mean(-1)/isolated.mean(-1)-1)
        np.testing.assert_allclose(observed,[r["actual_weak_area_change_percent"] for r in group],atol=1e-10)
        mothers=rng.integers(0,4,(10000,4));repeats=rng.integers(0,8,(10000,4,8))
        actual=100*(joint[mothers[:,:,None],repeats].mean(-1)/isolated[mothers[:,:,None],repeats].mean(-1)-1)
        prediction=np.array([r["predicted_weak_area_change_percent"] for r in group])
        error=(prediction[mothers]-actual).mean(-1)
        intervals.append(dict(condition=condition,gap_ns=gap,
            error_point_percentage_points=float(np.mean(prediction-observed)),
            error_percentile95_percentage_points=np.quantile(error,[.025,.975]).tolist(),
            resamples=10000,mother_pairs=4,digitizer_repeats=8))
    record=dict(frozen_hashes_checked=checked,finite_arrays=True,pilot_repeat_reproduction=replay,
        production_digitizations=64*8+20*8*6,cost_pilot_digitizations=2*2+2*6,
        single_optical_mothers=16,shared_optical_mothers=8,shared_mothers_subset_of_single=True,
        uncertainty=intervals,
        interpretation="Descriptive nested percentile intervals from four selected historical mother pairs and paired digital repeats; no independent transport/systematic uncertainty, not a universal performance interval. Gates remain the prespecified point-estimate rules.")
    destination=ROOT/"jupyter/run_records"/SHARED/"delivery-audit.json"
    destination.write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(record,ensure_ascii=False,indent=2))


if __name__=="__main__":main()

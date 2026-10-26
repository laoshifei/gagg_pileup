"""Resume frozen inference after restoring the optional CPU ML runtime."""
from pathlib import Path
import json
import os
import sys

FOLDER=Path(__file__).resolve().parent
PROJECT=FOLDER.parents[2]
for key in ('GAGG_CONDITIONAL_PROJECT','GAGG_RESOLUTION_PROJECT','GAGG_PAIR_PROJECT','GAGG_COMPONENT_PROJECT'):
    os.environ[key]=str(PROJECT)
os.environ.setdefault('MPLCONFIGDIR','/tmp/gagg-conditional-mpl')
local=Path(json.loads((PROJECT/'jupyter/run_records/alpha_pair_conditional_20261022/local-cache.json').read_text())['directory'])
sys.path.insert(0,str(local/'source' if (local/'source').exists() else FOLDER/'inference-source'))


def replay():
    import numpy as np
    import torch
    from pileup.analysis import alpha_pair_conditional_study as study
    from pileup.analysis.alpha_pair_ordered import OrderedNet,predict
    study.check_freeze();torch.set_num_threads(1)
    with np.load(study.OLD/'assets/shared.waveforms.npz') as f:y=f['voltage_mV'][:64]
    rows=[]
    for name in study.MODEL_NAMES:
        folder=FOLDER/'models'/name
        model=OrderedNet(study.load(FOLDER/'settings.json'),head='physical')
        model.load_state_dict(torch.load(folder/'model.pt',map_location='cpu',weights_only=True));model.eval()
        fits,_=predict(model,y,study.load(folder/'selection.json')['threshold'])
        expected=study.load(study.OLD/'predictions/shared'/(name+'.json'))['fits'][:64]
        difference=dict(times_ns=0.,charge_pe=0.,count_probabilities=0.)
        for a,b in zip(fits,expected):
            assert len(a['times_ns'])==len(b['times_ns'])
            for key in difference:
                if not len(a[key]):continue
                d=float(np.max(np.abs(np.asarray(a[key])-np.asarray(b[key]))))
                difference[key]=max(difference[key],d)
                np.testing.assert_allclose(a[key],b[key],rtol=2e-6,atol=1e-4)
        rows.append(dict(model=name,frames=64,max_absolute_difference=difference))
    study.save(study.RECORD/'restored-neural-replay.json',dict(passed=True,torch=torch.__version__,
        rows=rows,truth_inputs=[],weights_and_thresholds_changed=False,
        reference='Previously saved blind outputs, original batch of64; not new test labels'))
    print('RESTORED_MODEL_REPLAY_PASSED',rows,flush=True)


def main():
    from pileup.analysis import alpha_pair_conditional_study as study
    from pileup.analysis.alpha_pair_conditional_report import report,figures
    replay()
    if not (study.cache()/'test.seeds.json').exists():study.test_seeds()
    study.infer('test',workers=4)
    report('test');figures()
    print('INDEPENDENT_TEST_COMPLETE',flush=True)


if __name__=='__main__':main()

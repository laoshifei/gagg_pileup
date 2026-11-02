"""Check the rebuilt main executable against four current-environment controls.

Added after the user's explicit choice to keep the current environment. Uses
separate ROOTs/macros and preserves the original 16-control frozen experiment.
"""
from pathlib import Path
from datetime import datetime,timezone
import json
import os
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent
os.environ['GAGG_TRANSPORT_ENV_ROOT']=str(HERE.parents[2])
sys.path.insert(0,str(HERE/'analysis-source'))
import numpy as np
from pileup.analysis.alpha_transport_environment_study import ROOT,EXP,ART,REC,read,save,sha,load_root,verify
from pileup.analysis.alpha_transport_environment import relocate_macro,macro_semantics


def main():
    verify();assert not (EXP/'main-replay.json').exists()
    main_binary=ROOT/'g4/build/gagg';assert main_binary.is_file()
    source={}
    for old,digest in read(EXP/'source-freeze.json').items():
        p=ROOT/old;q=ROOT/'g4'/p.relative_to(EXP/'source/g4')
        assert sha(q)==digest,q;source[str(q.relative_to(ROOT))]=digest
    test=subprocess.run(['ctest','--test-dir',str(ROOT/'g4/build'),'--output-on-failure'],capture_output=True,text=True)
    (REC/'main-ctest.log').write_text(test.stdout+test.stderr);assert test.returncode==0
    linked=subprocess.run(['ldd',str(main_binary)],capture_output=True,text=True)
    (REC/'main-ldd.log').write_text(linked.stdout+linked.stderr)
    assert linked.returncode==0 and 'not found' not in linked.stdout+linked.stderr
    save(EXP/'main-build.json',dict(binary='g4/build/gagg',binary_sha256=sha(main_binary),ctest_passed=True,
        all_libraries_resolved=True,source_sha256=source,
        source_identical_to_frozen_experiment=True,
        configure=['cmake','-S','g4','-B','g4/build','-DCMAKE_BUILD_TYPE=Release',
            '-DGeant4_DIR=/opt/geant4/geant4-v1132/lib/cmake/Geant4','-DROOT_DIR=/opt/root/root63410/cmake'],
        build=['cmake','--build','g4/build','-j','4'],
        archive=read(REC/'main-build-migration.json')))
    plan=[r for r in read(EXP/'manifest.json')['runs'] if r['id'] in (0,4,8,12)]
    plan.sort(key=lambda r:r['id']);assert len(plan)==4
    save(EXP/'main-replay-plan.json',dict(created_utc=datetime.now(timezone.utc).isoformat(),
        user_direction='未必可用，你在现在的环境中让g4跑通',
        purpose='Verify current-environment determinism and main-build usability;not old/new population equivalence',
        references=plan,main_binary_sha256=sha(main_binary),new_independent_mothers=0,max_workers=1))
    save(EXP/'main-replay-source-freeze.json',{
        str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__).resolve(),EXP/'main-replay-plan.json',EXP/'main-build.json']})
    env=os.environ.copy();env.update(read(EXP/'environment.json')['current_environment']['g4_dataset_environment'])
    rows=[];last=0.;arrays_checked=0
    for task in plan:
        while True:
            delay=1.05-(time.time()-last)
            if delay>0:time.sleep(delay)
            output=ROOT/'g4/data'/(datetime.now().strftime('%Y%m%d_%Hh%Mm%Ss')+'.root')
            if not output.exists():break
            time.sleep(1.05)
        old=load_root(ROOT/task['root'])
        macro=EXP/'macros'/f"main_replay{task['id']:02d}_{output.stem}.mac"
        body=relocate_macro(old[0]['macro_contents'],output,ROOT/'g4/configs/sipm_hamamatsu_s13360_6050pe_v2.json')
        assert macro_semantics(body)==macro_semantics(old[0]['macro_contents']);macro.write_text(body)
        log=REC/'logs'/f'main_{output.stem}.log';started=time.monotonic();last=time.time()
        with log.open('w') as stream:
            result=subprocess.run([str(main_binary),str(macro)],cwd=ROOT/'g4/build',env=env,stdout=stream,stderr=subprocess.STDOUT)
        row=dict(id=task['id'],source_energy_MeV=task['source_energy_MeV'],seed=task['seed'],
            reference_root=task['root'],root=str(output.relative_to(ROOT)),macro=str(macro.relative_to(ROOT)),
            macro_sha256=sha(macro),log=str(log.relative_to(ROOT)),started_epoch=last,
            seconds=time.monotonic()-started,returncode=result.returncode,role='same_environment_main_build_repeat')
        save(REC/'attempts'/f'main_{output.stem}.json',row)
        assert result.returncode==0,row
        new=load_root(output)
        assert new[0]['random_engine_state']==old[0]['random_engine_state']
        assert new[0]['macro_contents']==body and new[0]['random_seed']==task['seed']
        assert json.loads(new[0]['sipm_config_json'])==json.loads(old[0]['sipm_config_json'])
        assert new[1]==old[1] and new[4]==old[4]
        for part in (2,3):
            for key in old[part]:np.testing.assert_array_equal(old[part][key],new[part][key]);arrays_checked+=1
        row.update(root_sha256=sha(output),root_bytes=output.stat().st_size,exact_selected_arrays=True,
                   edep_MeV=new[1]['edep_total_MeV'],status='completed')
        save(EXP/'runinfo'/f'main_{task["id"]:02d}.json',new[0])
        save(REC/'attempts'/f'main_{output.stem}.json',row);rows.append(row)
        print('Main build replay',task['id'],'exact',round(row['seconds'],2),'s',flush=True)
    save(EXP/'main-replay.json',dict(passed=True,rows=rows,exact_arrays=arrays_checked,
        native_replays=4,new_independent_mothers=0,main_binary='g4/build/gagg'))
    save(EXP/'current-environment-decision.json',dict(user_direction='Use the current environment;do not require the old environment',
        main_build_working=True,current_same_seed_reproducible=True,old_environment_bitwise_matching=False,
        old_16_controls_preserved=True,new_native_ROOT_count=20,new_independent_mothers=0,
        next='UseGeant4 11.3.2/ROOT6.34.10 for future work and record the migration;establish current-environment response controls before interpreting old-model transfer as independent confirmation',
        confirmation_design_status='336-mother design remains a draft for a separate current-environment campaign;no such fresh sample generated'))
    print('Current main build runs and repeats exactly;20new control ROOTs total',flush=True)


if __name__=='__main__':main()

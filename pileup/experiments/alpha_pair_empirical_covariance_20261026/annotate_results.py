"""Post-run convergence audit and balanced reading notes for the delivered book."""
from pathlib import Path
import argparse
import json
import os
import sys

folder=Path(__file__).resolve().parent
project=folder.parents[2]
sys.path.insert(0,str(folder/'report-source'))
for name in ('EMPIRICAL','WEAK','NULL','CONDITIONAL','RESIDUAL','RESOLUTION','PAIR','COMPONENT'):
    os.environ['GAGG_'+name+'_PROJECT']=str(project)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--notebook',action='store_true')
    args=parser.parse_args()
    output=project/'jupyter/artifacts/alpha_pair_weak_diagnostics_20261026'
    cases=json.loads((project/'jupyter/artifacts/alpha_pair_null_gate_20261023/test/labels.cases.json').read_text())
    rows=[json.loads(s) for s in (output/'diagnostics.jsonl').read_text().splitlines()]
    result=[]
    for gap in (8,40,200):
        for method in ('fixed_true_times','truth_start_joint','truth_start_joint_ols','oracle_subtract_other'):
            fits=[r['diagnostics'][method] for r in rows if cases[r['index']]['condition']=='high_weak'
                  and cases[r['index']]['gap_ns']==gap]
            result.append(dict(gap_ns=gap,method=method,frames=len(fits),
                converged=sum(bool(f.get('converged')) for f in fits),
                max_evaluations=max(f['evaluations'] for f in fits)))
    (output/'convergence-audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print('CONVERGENCE_AUDITED',len(result),'groups',flush=True)
    if args.notebook:
        import nbformat
        from pileup.analysis.alpha_pair_conditional_notebook import atomic_notebook
        from pileup.analysis.alpha_pair_resolution_data import digest
        path=project/'jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb'
        previous=path.read_bytes();nb=nbformat.reads(previous.decode(),as_version=4)
        identifier='weak-noise-convergence-20261026'
        tradeoff_identifier='weak-noise-tradeoff-20261026'
        if any(c.get('id')==identifier for c in nb.cells):
            assert any(c.get('id')==tradeoff_identifier for c in nb.cells)
            print('RESULT_NOTES_ALREADY_PRESENT',flush=True);return
        index=next(i for i,c in enumerate(nb.cells) if '<a id="empirical-pulse-covariance-20261026">' in c.source)
        original=json.loads(json.dumps(nb.cells))
        lines=['**优化器收敛核验**：这里的converged只表示满足局部数值停止条件，不证明全局最优或Q/t物理正确。',
            '', '| Δ / ns | 真值启动GLS正常停止 / 对数 | 最大函数评估次数 |', '| ---: | ---: | ---: |']
        for r in result:
            if r['method']=='truth_start_joint':
                lines.append(f"| {r['gap_ns']} | {r['converged']}/{r['frames']} | {r['max_evaluations']} |")
        lines+=['', '40与200 ns的全部真值启动GLS均正常停止，最大评估次数分别25和17，低于120的上限；不能把这些失败直接归为耗尽迭代预算。8 ns有1例达到上限，仍保留原结果，没有选择性重跑。完整[收敛核验](artifacts/alpha_pair_weak_diagnostics_20261026/convergence-audit.json)另列固定时刻、OLS及oracle对照。']
        cell=nbformat.v4.new_markdown_cell('\n'.join(lines));cell['id']=identifier
        nb.cells.insert(index,cell)
        report=json.loads((project/'jupyter/artifacts/alpha_pair_empirical_covariance_20261026/development-report.json').read_text())
        successes={r['method']:f"{r['success']}/{r['eligible_frames']}" for r in report['summary']
            if r['condition']=='weak_high' and r['gap_ns']==8.}
        text=('**局部改善与整体选择分开读**：弱→强、8 ns，旧AR同门槛样本量对照为'+successes['old_ar']+
            '；经验对角、16模式×0.5、16模式×0.9、32模式×0.9分别为'+
            '、'.join(successes[k] for k in ('diagonal','rank16_half','rank16_full','rank32_full'))+
            '。这是明确的开发局部收益；同时，强→弱及较强组合存在损失，四种方案均未通过整体选用条件。\n\n'+
            '表中的old_ar只作对照，不参与新候选选择，所以其eligible为false不表示控制失败。'+
            '额外输出按通常10 ns一对一匹配定义，包括时刻偏差超出容限的输出，不必然意味着真实没有物理事件。')
        tradeoff=nbformat.v4.new_markdown_cell(text);tradeoff['id']=tradeoff_identifier
        destination=next(i for i,c in enumerate(nb.cells) if '<a id="weak-noise-conclusion-20261026">' in c.source)
        nb.cells.insert(destination,tradeoff)
        assert [c for c in nb.cells if c.get('id') not in (identifier,tradeoff_identifier)]==original
        atomic_notebook(path,nb,previous)
        record=project/'jupyter/run_records/alpha_pair_empirical_covariance_20261026/notebook-annotations.json'
        record.write_text(json.dumps(dict(passed=True,existing_cells_preserved=True,inserted_indices=[index,destination],
            script_sha256=digest(Path(__file__)),notebook_sha256=digest(path)),indent=2)+'\n')
        print('RESULT_NOTES_DELIVERED',flush=True)


if __name__=='__main__':main()

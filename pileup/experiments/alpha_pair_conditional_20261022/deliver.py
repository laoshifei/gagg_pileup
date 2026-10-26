"""Deliver the local notebooks and execute only the newly appended cells."""
from pathlib import Path
import base64
from collections import Counter
import gzip
import hashlib
import json
import os
import re
import sys

FOLDER=Path(__file__).resolve().parent
PROJECT=FOLDER.parents[2]
for key in ('GAGG_CONDITIONAL_PROJECT','GAGG_RESOLUTION_PROJECT','GAGG_PAIR_PROJECT','GAGG_COMPONENT_PROJECT'):
    os.environ[key]=str(PROJECT)
sys.path.insert(0,str(FOLDER/'inference-source'))


def add_result_interpretation(nb, study, first):
    """Describe frozen outcomes; this postprocessing never changes inference."""
    import nbformat
    report=study.load(study.OUTPUT/'test-report.json')
    summary=report['summary']
    with gzip.open(study.OUTPUT/'test-scores.jsonl.gz','rt') as stream:
        scores=[json.loads(line) for line in stream]
    counts=[]
    for method in ('blind_fit','ml_1','ml_2','old_hybrid_1','selective_1','selective_2','conditional_fit_only','conditional_fit'):
        rows=[r for r in scores if r['condition']=='high_low' and r['gap_ns']==200. and r['method']==method]
        assert len(rows)==40
        histogram=Counter(r['count'] for r in rows)
        counts.append(dict(method=method,frames=40,output_K0=histogram[0],output_K1=histogram[1],
            output_K2=histogram[2],strict_pair_success=sum(r['success'] for r in rows)))
    errors=[{k:r[k] for k in ('method','condition','success','frames','q_error_p16_p50_p84','q_error_std')}
        for r in summary if r['kind']=='single' and r['condition']=='low' and r['method'] in
        ('ml_1','old_hybrid_1','mean_hybrid_1','correlated_hybrid_1','selective_1')]
    study.save(study.OUTPUT/'test-diagnostics.json',dict(
        report_sha256=study.digest(study.OUTPUT/'test-report.json'),
        purpose='Post-test explanation only; no parameter selection or inference change',
        high_low_200ns_counts=counts,low_single_conditional_errors=errors))
    count={r['method']:r for r in counts}
    def result(method,condition,gap=None):
        return next(r for r in summary if r['method']==method and r['condition']==condition and r['gap_ns']==gap)
    low=result('selective_1','low')
    interpretation=f'''**独立结果的具体解读**：中／中组实际Q约9400 pe/事件。Δ=8 ns时，原始网络两个种子分别19/40、31/40成功，固定模板重拟合均0/40，新均值＋GLS重拟合均40/40；Δ=12–200 ns各已测点也均40/40。高／高组在8–200 ns各点同样为40/40。40/40的逐点Wilson95区间约91.2%–100%，不是绝对正确保证。中／中8 ns仅修正均值而不作相关加权时，种子1为1/40；仅使用旧波形候选起点的新GLS为0/40，而加入网络起点为40/40。这支持继续研究条件均值、相关涨落和搜索起点的共同作用，尚不能给出普适8 ns阈值。

**两类退化必须保留**：低／低组30 ns时，旧混合方法为40/40，新GLS为34/40；低电荷原生单事件的新GLS电荷相对误差中位数为{100*low['q_error_p16_p50_p84'][1]:.2f}%，成功{low['success']}/{low['frames']}，只修正均值则80/80且中位误差约−0.68%。因此相关加权并非对所有Q都合适。先强后弱200 ns时，网络种子1有{count['ml_1']['output_K1']}/40帧只报一个事件，种子2为{count['ml_2']['output_K1']}/40；保持网络K的重拟合无法补出缺失槽，严格成功只有0/40和4/40。允许重新选择K的条件拟合为36/40，原盲拟合为39/40；这两者的双事件结果不能忽略单事件误拆控制。

新GLS及条件K选择的空、低／中／高单事件各组额外输出均0/80；旧盲拟合的高电荷单事件为34/80。0/80的Wilson95上限约4.58%。不将测试后的优势条件拼接成一个新方法，也不把两个网络种子的差异当作置信区间。

**下一轮优先项**：先在开发数据上修正低Q的相关权重偏差，并让残差触发新候选、重新比较K，解决强后弱的漏检；再做新的光学位置／源能量和α/γ测试。上面的误差诊断是在冻结测试后对已有输出的解释，没有修改本轮参数或重新选择门槛。仅波形起点对照继承旧模板选出的单个K=2起点，不能代表充分多起点纯拟合的上限。'''
    md=nbformat.v4.new_markdown_cell(interpretation)
    md['id']='conditional-20261022-results-interpretation'
    code=nbformat.v4.new_code_cell('''conditional_diagnostics = json.loads((conditional_output / "test-diagnostics.json").read_text())
print("先强后弱、200 ns：输出事件数与整对成功；各方法同一40对")
display(pd.DataFrame(conditional_diagnostics["high_low_200ns_counts"]))
print("低电荷原生单事件：匹配条件下的Q相对误差及标准差，须与整帧成功数一起读")
display(pd.DataFrame(conditional_diagnostics["low_single_conditional_errors"]))''')
    code['id']='conditional-20261022-results-diagnostics'
    ids={md.id,code.id}
    nb.cells[first:]=[c for c in nb.cells[first:] if c.get('id') not in ids]
    position=next(i+1 for i,c in enumerate(nb.cells) if c.cell_type=='markdown' and
        '## 44. 本轮结果、限制和下一步' in c.source)
    nb.cells[position:position]=[md,code]


def main():
    import nbformat
    from IPython.core.interactiveshell import InteractiveShell
    from IPython.utils.capture import capture_output
    from pileup.analysis import alpha_pair_conditional_study as study
    from pileup.analysis.alpha_pair_conditional_notebook import deliver,catalog,atomic_notebook
    assert study.load(study.RECORD/'science-validation.json')['passed']
    path=PROJECT/'jupyter/codex_0004_Alpha_Pileup_ML_Comparison_20261014.ipynb'
    nb=nbformat.read(path,as_version=4)
    if not any(c.get('id','').startswith('conditional-20261022-') for c in nb.cells):deliver()
    cat=nbformat.read(PROJECT/'jupyter/Run_Catalog.ipynb',as_version=4)
    if not any(c.get('id')=='conditional-catalog-20261022' for c in cat.cells):catalog()
    previous=path.read_bytes();nb=nbformat.reads(previous.decode(),as_version=4)
    validation=study.load(study.RECORD/'notebook-validation.json');first=validation['preserved_cells']
    original=json.loads(json.dumps(nb.cells[:first]))
    add_result_interpretation(nb,study,first)
    InteractiveShell.instance();namespace={};executed=0
    os.chdir(PROJECT/'jupyter')
    for cell in nb.cells[first:]:
        if cell.cell_type!='code':continue
        with capture_output(stdout=True,stderr=True,display=True) as captured:
            exec(compile(cell.source,'<conditional-notebook>','exec'),namespace)
        assert not captured.stderr,captured.stderr
        cell.outputs=[]
        if captured.stdout:cell.outputs.append(nbformat.v4.new_output('stream',name='stdout',text=captured.stdout))
        for output in captured.outputs:
            cell.outputs.append(nbformat.v4.new_output('display_data',data=output.data,metadata=output.metadata))
        executed+=1;cell.execution_count=executed
    assert json.loads(json.dumps(nb.cells[:first]))==original
    anchors=[];images=[];forbidden=[]
    for cell in nb.cells[first:]:
        if cell.cell_type=='markdown':
            anchors.extend(re.findall(r'id="codex-figure-([^"]+)"',cell.source))
            forbidden.extend(re.findall(r'\]\([^)]*\.(?:pdf|ipynb|md)(?:[?#][^)]*)?\)',cell.source,flags=re.I))
        else:
            for output in cell.outputs:
                if 'image/png' in output.get('data',{}):
                    value=output.data['image/png'];value=''.join(value) if isinstance(value,list) else value
                    images.append(hashlib.sha256(base64.b64decode(value)).hexdigest())
                assert output.output_type!='error'
    assert len(anchors)==len(set(anchors))==4 and len(images)==4 and not forbidden
    expected={study.digest(p) for p in (study.OUTPUT/'figures').glob('*.png')}
    assert set(images)==expected
    atomic_notebook(path,nb,previous)
    validation.update(actual_new_code_cells_executed=executed,new_outputs_regenerated=True,
        historical_cells_unchanged=True,new_figures_unique=True,prohibited_new_links=[],
        appended_cells=len(nb.cells)-first,post_test_interpretation_added=True,
        notebook_sha256=study.digest(path),delivery_script_sha256=study.digest(Path(__file__)))
    study.save(study.RECORD/'notebook-validation.json',validation)
    print(json.dumps(validation,ensure_ascii=False),flush=True)


if __name__=='__main__':main()

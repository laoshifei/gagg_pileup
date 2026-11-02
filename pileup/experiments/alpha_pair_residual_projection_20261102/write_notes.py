"""Append the completed diagnostic to existing Chinese project records."""
from pathlib import Path
import base64
import csv
import json

ROOT=Path(__file__).resolve().parents[3]
EXP=Path(__file__).resolve().parent
NAME=EXP.name
ART=ROOT/'jupyter/artifacts'/NAME
REC=ROOT/'jupyter/run_records'/NAME
NOTEBOOK='codex_0007_SiPM_Energy_Calibration.ipynb'


def read(path):return json.loads(path.read_text())
def rows(name):return list(csv.DictReader((EXP/'tables'/name).open()))
def label(condition,gap):return ('强→弱' if condition=='high_weak' else '弱→强')+f' {float(gap):g} ns'
def table(headers,body):
    return '| '+' | '.join(headers)+' |\n| '+' | '.join('---' for _ in headers)+' |\n'+''.join('| '+' | '.join(map(str,row))+' |\n' for row in body)
def md(text,id):return dict(cell_type='markdown',id=id,metadata={},source=text.splitlines(keepends=True))


def figure(section,name):
    source=("from pathlib import Path\nfrom IPython.display import display, Image\n"
        "root = Path.cwd()\nif root.name == 'jupyter':\n    root = root.parent\n"
        f"display(Image(filename=str(root / 'jupyter/artifacts/{NAME}/figures/{name}.png')))\n")
    encoded=base64.b64encode((ART/'figures'/(name+'.png')).read_bytes()).decode()
    return dict(cell_type='code',id=f'residual-projection-figure-{section}',metadata={},execution_count=None,
        source=source.splitlines(keepends=True),outputs=[dict(output_type='display_data',metadata={},
        data={'image/png':encoded,'text/plain':[f'图 0007-{section}-1（已生成PNG）']})])


def append_notebook(name,cells):
    path=ROOT/'jupyter'/name;before=read(REC/(name+'.before'));current=read(path)
    if current!=before:raise AssertionError(f'Notebook changed since backup: {name}')
    current['cells'].extend(cells)
    path.write_text(json.dumps(current,ensure_ascii=False,indent=1)+'\n')


def main():
    result=read(EXP/'result-summary.json');search=result['search']
    projection=[r for r in result['projections'] if r['weight']=='full']
    nonlinear=result['nonlinear']
    def fit(gap,variant):
        return next(r for r in nonlinear if r['condition']=='high_weak' and r['gap_ns']==gap
                    and r['method']=='spatial12' and r['variant']==variant)
    short=[]
    for gap in (8.,40.):
        raw,removed=fit(gap,'true_E_start'),fit(gap,'delete_strong_isolated')
        short.append(f"{gap:g} ns弱E偏差{raw['weak_bias_percent']:+.2f}%→{removed['weak_bias_percent']:+.2f}%，RMSE {raw['weak_rmse_percent']:.2f}%→{removed['weak_rmse_percent']:.2f}%")
    overview=(f"复用20条8次平均波形，完成120组线性投影与88条正式非线性对照。40条真实E起点对照中，"
        f"{search['lower_cost']}条找到超过预设数值容差的更低代价，最大能量变化{search['maximum_energy_change_percent']:.5f}%；"
        "这项对照未支持原搜索起点是本批主要偏差来源，也不证明全局最优。"
        "空间模型精确删除强脉冲孤立响应残差后，"+'；'.join(short)+"。删除使用真实E和祖先分量，只用于误差诊断。")
    limitation=("仅四个已暴露历史母对，五个条件，每条均值由8次数字化组成；重复不增加独立输运样本，均值不是无噪声理论期望。"
        "强/弱孤立残差同时包含固定光学母事件实现的差异和有限数字化涨落，不能全部归为普遍模板系统偏差。"
        "共享与孤立数字化的随机数消耗不同，四项互相关，删除后的波形是人工构造，不能把其收益当作可实现算法效率。"
        "真值投影、真值起点和残差删除均标为oracle；权重对照仍给定真实K/t。未测盲漏检、误拆分、误合并、能谱选择效率或高率性能。"
        "固定95/351 ns、快比例0.3、29 ns恢复代理、14400微胞和12×12开发先验；未传播先验、能标及器件系统误差。")
    next_step=("下一步优先利用原有开发母事件研究强脉冲孤立形状的逐事件变化及其与能量的关联，"
        "在独立开发划分上冻结低维形状/噪声处理方案，再用孤立与受控强弱双脉冲检查它是否吸收真实弱信号。"
        "固定权重对照须作为完整候选方案验证，不按本批每个母事件或间隔选择最优权重。"
        "本轮不调整恢复强度、不升级默认重建、不增加独立输运；后续改变模型或权重时另建实验快照。")
    projection_table=table(['条件','模型','强孤立','弱孤立','共享修正','暗/电子学','总投影','非线性弱E偏差/%'],[
        [label(r['condition'],r['gap_ns']),r['method'],*[f"{100*r[k+'_log']:+.2f}" for k in
         ('strong_isolated','weak_isolated','shared_correction','dark_electronics')],
         f"{100*r['mean_weak_delta_log']:+.2f}",f"{r['mean_weak_old_fit_bias_percent']:+.2f}"] for r in projection])
    variants=[('original_3starts','原full/三起点'),('diagonal_3starts','对角/三起点'),('unit_fir_3starts','单位FIR/三起点')]
    weight_table=table(['强→弱间隔','空间模型权重','弱E偏差/%','弱E RMSE/%','严格成功','边界','未收敛'],[
        [f'{gap:g} ns',title,f"{fit(gap,var)['weak_bias_percent']:+.2f}",f"{fit(gap,var)['weak_rmse_percent']:.2f}",
         f"{fit(gap,var)['strict']}/4",fit(gap,var)['boundary'],fit(gap,var)['nonconverged']]
        for gap in (8.,40.) for var,title in variants])
    deletion_variants=[('true_E_start','原平均波形'),('delete_strong_isolated','删强孤立残差'),
        ('delete_weak_isolated','删弱孤立残差'),('delete_shared_correction','删共享修正残差'),('delete_dark_electronics','删暗/电子学余项')]
    deletion_table=table(['强→弱间隔','真值辅助输入','弱E偏差/%','弱E RMSE/%','严格成功','边界','未收敛'],[
        [f'{gap:g} ns',title,f"{fit(gap,var)['weak_bias_percent']:+.2f}",f"{fit(gap,var)['weak_rmse_percent']:.2f}",
         f"{fit(gap,var)['strict']}/4",fit(gap,var)['boundary'],fit(gap,var)['nonconverged']]
        for gap in (8.,40.) for var,title in deletion_variants])
    split=rows('symmetric-model-change.csv');split_body=[]
    for condition,gap in [('high_weak',8),('high_weak',40),('high_weak',200),('high_weak',1200),('weak_high',40)]:
        group=[r for r in split if r['condition']==condition and float(r['gap_ns'])==gap and r['weight']=='full']
        split_body.append([label(condition,gap),*[f"{100*sum(float(r[key]) for r in group)/4:+.2f}" for key in
            ('residual_change_log','sensitivity_change_log','total_change_log')]])
    split_table=table(['条件','平均残差变化项','敏感方向变化项','空间−均匀总变化'],split_body)
    methods=("沿用上一阶段的原始电压、空间先验、孤立均值能标、各自Q→光量映射及噪声标定，不重新数字化。"
        "观测窗0–2048 ns、2 ns采样，保留完整1024点；真实时间为128 ns及128+间隔。输入路径与原种子映射见cases.json、input-freeze.json和父实验。\n\n"
        "令A为真E处孤立模型，M为真E处共享模型之和，I/S为物理孤立/共享祖先信号，Y为完整含噪电压、b为固定背景。"
        "四项为I强−A强、I弱−A弱、(ΣS−ΣI)−(M−ΣA)、Y−b−ΣS；逐重复严格相加为Y−b−M。末项含背景期望失配，不只是一份白噪声。\n\n"
        "θ=ln E，J=W·T·∂M/∂θ；用SVD计算δθ=J⁺W·T·r。T是继承的FIR变换，W固定在原完整平均波形的电压初值上。"
        "强/弱导数耦合被共同求解，也单独核对消去强导数后的弱方向。原始电压影响函数用于累计时间图，避免把白化后的索引误称为局部物理时间贡献。"
        "各项只在线性δlnE中可加；exp(δlnE)−1用于总变化，不能分别指数化再相加。另给出真实能量边界下的线性有界解。\n\n"
        "三种固定权重为原full、只保留边际方差的diagonal及unit FIR；后者不是原始电压OLS。"
        "40条full真E单起点、16条空间短间隔权重三起点、32条空间短间隔单项删除真E起点，共88条。"
        "继承logE边界[0.1183134541,8.6611928940] MeV、每起点80次评估及原求解容差。"
        "比较权重时仍使用原三个电压起点；删除项时权重固定为原波形full、起点固定为真E，与同起点原波形对照。"
        "不同权重的代价绝对值不横向排名；非线性删除效果不相加。")
    validation=("代数/变换、梯度张成与正交残差控制通过；120组分解的最大电压闭合误差1.14×10⁻¹³ mV，"
        "上一阶段代价复算最大相对差9.66×10⁻¹⁵，小幅实际模型扰动的线性logE误差最大1.83×10⁻⁵。"
        "正式88条＋先导4条＋独立输出目录重放4条，共96条非线性拟合；先导/正式及重放/正式的能量、代价、收敛、边界和起点逐项一致。"
        f"正式未收敛{result['nonconverged']}条，所有行均保留；边界与每组强能量偏差/RMSE详见tables/。"
        "线性贡献的10,000次母对/重复嵌套重采样在方法和四项间配对，梯度与权重固定，区间只作四母对描述。")
    # A small table includes the paired uncertainty, without treating repeat draws as mothers.
    intervals_table=table(['条件','空间−均匀弱δlnE×100','配对描述性95%范围'],[
        [label(r['condition'],r['gap_ns']),f"{100*r['full_spatial_minus_uniform_total_log']:+.2f}",
         f"[{100*r['paired_log_ci95'][0]:+.2f}, {100*r['paired_log_ci95'][1]:+.2f}]"] for r in result['paired_intervals']])
    sections=[
        md('## 93. 完整残差投影：研究问题与冻结对照（2026-11-02）\n\n'+overview+'\n\n'+methods+'\n','residual-projection-93'),
        md('## 94. 四项残差沿弱能量方向的贡献\n\n'+projection_table+'\n'
           '表中四项与总投影的单位均为100×δlnE，最后一列才是相对E百分比；分别对四母等权平均。大偏差时一阶反演与非线性结果不同，逐母差异及边界也影响均值，不能只对组平均log变化取指数。\n\n'
           +intervals_table+'\n图 0007-94-1：(a)均匀、(b)空间模型在full权重下的四项贡献。正负柱分开堆叠，黑横线为线性总量，黑点为原非线性拟合的平均log能量变化。\n','residual-projection-94'),
        figure(94,'fig_0007_94_1_residual_contributions'),
        md('## 95. 完整观测窗与模型变化的分解\n\n'
           '第一母对block 5、强→弱40 ns作固定示例。全波形幅度残差较小，仍可在消去强能量方向后投影成很大的弱能量误差。导数相关系数及条件数是本权重下的局部几何，不是校准的Fisher界或普适物理分辨极限。\n\n'
           '对两个原始电压影响算子P和残差r，使用对称恒等式：Δ空间−Δ均匀 = ½(P均匀+P空间)(r空间−r均匀) + ½(P空间−P均匀)(r空间+r均匀)。这是描述性代数分解，不是唯一因果归因；下表单位仍为100×δlnE。\n\n'
           +split_table+'\n图 0007-95-1：(a)完整均值电压及真E处模型，(b)空间四项残差，(c)原始电压到弱logE的影响函数，(d)四项沿真实时间累计的贡献。虚线标出给定时刻，所有面板保留0–2048 ns。\n','residual-projection-95'),
        figure(95,'fig_0007_95_1_full_window_projection'),
        md('## 96. 有限搜索与固定权重对照\n\n'
           f"40条真E起点对照有{search['lower_cost']}条低于原三起点代价、{search['higher_cost']}条高于原代价（阈值max(10⁻⁵,10⁻⁷×旧代价)），"
           f"最大E变化{search['maximum_energy_change_percent']:.5f}%。这只排查了一种额外起点，不证明全局最优。\n\n"
           +weight_table+'\n权重对照使用同一平均电压及原三起点，不使用真实E初始化；真实K/t仍给定。严格成功要求两个E_dep均在±10%且收敛，4条均值的比例不是单次盲效率。完整强分量误差与SD见实验表。\n\n'
           '图 0007-96-1：(a)强→弱8 ns、(b)40 ns，四条母事件曲线展示各固定权重的弱E误差；灰带为±10%。包括原均匀full作为参照，不逐母选最佳权重。\n','residual-projection-96'),
        figure(96,'fig_0007_96_1_weight_contrasts'),
        md('## 97. 真值辅助的单项残差删除\n\n'+deletion_table+'\n'
           '删除强孤立项，相当于在原共享电压中替换该母事件的强孤立响应残差，而非直接取出一个可观测的真实强分量；该操作需要真E与祖先信号。其余三项分别删除，保留原full权重及真实E起点。各项相关且拟合非线性，改善幅度不能求和，也不能称为新重建算法。\n\n'
           '图 0007-97-1：(a)8 ns、(b)40 ns，原波形及四种单项删除的逐母弱E误差；灰带为±10%，图中全部输入均属于真值辅助诊断。\n','residual-projection-97'),
        figure(97,'fig_0007_97_1_oracle_deletions'),
        md('## 98. 验证、限制与下一步\n\n'+validation+'\n\n'+limitation+'\n\n'+next_step+'\n\n'
           '沿用项目L89与冻结源码，不新增外部文献结论。配置的共享库/media/sf_works/papers仍不存在，未核原文或改变人工审核状态；共享库恢复后核对条件形状/协方差的既有证据。\n\n'
           '实际模型、协议、输入哈希、源码、运行器和表保存在pileup/experiments/'+NAME+'；本地数组和图在jupyter/artifacts/同名目录，日志在jupyter/run_records/同名目录。此次分析没有新ROOT或新数字化，也没有Poisson源率扫描。\n','residual-projection-98')]
    append_notebook(NOTEBOOK,sections)
    append_notebook('Research_Plan_20260921.ipynb',[md('## 31. 残差投影完成：强孤立形状与权重的下一步（2026-11-02）\n\n'
        +overview+'\n\n'+next_step+'\n\n'+limitation+'\n\n完整结果在codex_0007第93–98节、四幅固定编号图及采用记录L90。\n','residual-projection-plan-31')])
    permission=read(ROOT/'jupyter/run_records/permissions_20261102/permission-restoration.json')
    append_notebook('Run_Catalog.ipynb',[md('## 2026-11-02：权限恢复、同步与残差诊断\n\n'
        f"拷贝增加的执行位已按Git记录恢复，共{len(permission['tracked'])}个受跟踪文件；另恢复{len(permission['local_data'])}个本地笔记/数据文件。"
        '保留原有可执行入口、构建程序及虚拟环境；清单在permissions_20261102/permission-restoration.json。恢复后内容无差异，没有空提交；git push origin main返回Everything up-to-date，开始时HEAD为402f792。\n\n'
        +overview+'\n\n'+validation+'\n\n'
        '本轮复用父实验20条件的160条已有数字化及相应祖先数组，不生成ROOT或新数字化。输入原种子、母ID、时间与沉积能量见cases.json及父实验；新拟合确定性，bootstrap种子731100–731104。源码/配置/输入先冻结，最多8进程、BLAS/OMP单线程，g++按build.json构建。运行日志及阶段计时在同名run_records目录；4条成本先导通过预设成本门。\n\n'
        '新模块最初与历史residual_study/report名称碰撞，已将根目录历史模块逐字节恢复到402f792，并将新模块改为residual_projection_*。实际运行一直使用本实验冻结源码，保留其原内部模块名；namespace-revision.json记录对应关系，未改变计算或冻结输入。\n\n'
        '96条非线性拟合分别为正式88、先导4和新目录重放4。无新源率/宏；已有ROOT条件以原RunInfo与父实验记录为准。新阶段完成后只作本地提交，未自动上传后续修改。\n','residual-projection-catalog')])
    experiment_readme=('# 弱能量残差投影、权重与搜索诊断\n\n'+overview+'\n\n## 冻结方案\n\n'+methods+
        '\n\n## 残差投影\n\n'+projection_table+'\n四项及总量为100×δlnE，最后一列为相对E百分比；大偏差和边界下使用非线性对照解释。\n\n'+intervals_table+
        '\n## 固定权重的非线性对照\n\n'+weight_table+'\n不同权重的代价绝对值不直接比较；全部输入均给定真实K/t。\n\n## 真值辅助删除\n\n'+deletion_table+
        '\n所有删除输入都使用真实能量与祖先信号；不作为新算法效率。五种条件及两模型的强E误差、弱误差SD、边界和未收敛完整保存在tables/。\n\n## 验证与复现\n\n'+validation+
        '\n\n从项目根目录运行。source/保留实际Python/C++、JSON配置及本地NPZ模型数组，bin/age.so可按build.json中的相对命令构建；依赖版本见environment.json。模型数组、输入电压和祖先数组不随Git分发，须在input-freeze.json指定的相对位置提供。\n\n'
        '```bash\npython3 -B pileup/experiments/'+NAME+'/run.py verify\npython3 -B pileup/experiments/'+NAME+'/run_report.py\npython3 -B pileup/experiments/'+NAME+'/validate_delivery.py\n```\n\n'
        '完整计算入口为run.py execute，包括数值检查、投影、成本先导、88条正式拟合及4条重放。它要求一个具备相同输入、尚无本实验输出的工作副本；已有输出目录会被拒绝覆盖，无须切换Git提交。报告可从已有数组重建表和图。源码起点402f792，参数和输入在看本轮输出前冻结；报告源码与最终交付哈希分别记录，最终交付冻结不代表预注册。\n\n'
        '新可复用模块名为alpha_pair_residual_projection{,_study,_report}.py。冻结快照内部使用准备时的study/report模块名；namespace-revision.json记录根目录命名修复与严格内容映射，历史20261023模块已恢复，计算代码未改动。\n\n## 限制与下一步\n\n'
        +limitation+'\n\n'+next_step+'\n\n研究过程追加到本地Notebook0007第93–98节，图0007-94-1至97-1。共享文献库不可用，沿用项目L89，不新增文献性能主张；采用记录L90记录依据范围。\n')
    (EXP/'README.md').write_text(experiment_readme)
    for name in ('README.md','literature-use.md'):
        path=ROOT/name if name=='README.md' else ROOT/'docs'/name
        assert path.read_bytes()==(REC/(name+'.before')).read_bytes(),name
    root_append=('\n\n2026-11-02继续完成[弱能量残差投影与权重/搜索诊断](pileup/experiments/'+NAME+'/README.md)：'+overview+
        ' 在本批强→弱40 ns条件，原均匀共享响应残差与强孤立残差部分抵消，修正共享响应后需要单独处理强脉冲形状/权重。'
        '结果限于四个历史母对的真值辅助诊断，没有默认算法升级。完整结果在本地codex_0007第93–98节、四幅编号图及[采用记录L90](docs/literature-use.md#residual-projection-20261102)。'+next_step+'\n')
    with (ROOT/'README.md').open('a') as f:f.write(root_append)
    literature=('\n\n<a id="residual-projection-20261102"></a>\n## L90 — 2026-11-02：残差方向而非面积比解释弱能量偏差\n\n'
        '**项目结果**：'+overview+'\n\n**采用与区分**：将强/弱孤立响应、共享修正和暗/电子学余项分开，以真E处的联合logE导数检查残差投影；四项在线性log能量中可加，但非线性删除收益不相加。'
        '原均匀模型在强→弱40 ns有强孤立项与共享修正项部分抵消；共享面积控制改善后，旧抵消减小，能量结果仍需单独检验。\n\n'+projection_table+
        '\n四项/总量单位100×δlnE，最后一列为相对E百分比。均匀→空间变化还包括敏感方向变化，完整对称分解保存在实验表。\n\n**权重与输入边界**：'+weight_table+
        '\n权重对照固定原电压导出的三个起点，仍给定真实K/t；oracle初始化/删除使用真E。所有结果均保留失败与边界，不从每条件挑选最优方法。\n\n**核验**：'+validation+
        '\n\n**限制与下一步**：'+limitation+' '+next_step+
        '\n\n**文献核验层级**：依据L89、项目冻结模型和本轮代数/数值对照。配置的共享库/media/sf_works/papers仍不可用，未新增外部论文结论、核原文或改变审核状态。共享库恢复后，定点核对强脉冲条件形状、相关噪声与弱脉冲可辨识性的既有精读依据。'
        '\n\n位置：`pileup/experiments/'+NAME+'/`，根目录可复用模块`alpha_pair_residual_projection.py`、`alpha_pair_residual_projection_study.py`、`alpha_pair_residual_projection_report.py`；本地codex_0007第93–98节、图0007-94-1至97-1、研究计划第31节。\n')
    with (ROOT/'docs/literature-use.md').open('a') as f:f.write(literature)
    (REC/'delivery-summary.txt').write_text(overview+'\n\n'+next_step+'\n')
    print('Appended10Notebook0007cells,plan31,catalog,L90;oldcells preserved')


if __name__=='__main__':main()

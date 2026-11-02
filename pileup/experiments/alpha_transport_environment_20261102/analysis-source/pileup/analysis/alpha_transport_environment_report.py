"""Append the current-environment build and transport checks to project records."""
from __future__ import annotations

import csv
import json
import os

os.environ.setdefault('MPLCONFIGDIR', '/tmp/gagg_transport_mpl')
os.environ.setdefault('IPYTHONDIR', '/tmp/gagg_transport_ipython')
os.environ.setdefault('JUPYTER_RUNTIME_DIR', '/tmp/gagg_transport_jupyter')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import nbformat
from nbclient import NotebookClient
import numpy as np

from .alpha_transport_environment_study import ROOT, EXP, ART, REC, NAME, read, save, sha, verify


def table(headers, rows):
    return '\n'.join(['| ' + ' | '.join(headers) + ' |',
                      '| ' + ' | '.join(['---'] * len(headers)) + ' |'] +
                     ['| ' + ' | '.join(map(str, row)) + ' |' for row in rows])


def write_csv(name, rows):
    (EXP / 'tables').mkdir(exist_ok=True)
    with (EXP / 'tables' / f'{name}.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def figures(rows):
    plt.rcParams.update({'font.size': 10, 'figure.dpi': 140, 'axes.grid': True, 'grid.alpha': .2})
    paths = {}
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    for group, energy in enumerate((1.625, 2.5, 5.486, 8.)):
        selected = [r for r in rows if r['group'] == group]
        x = [r['id'] for r in selected]
        axes[0].scatter(x, [1000 * (r['current']['edep_MeV'] - r['old']['edep_MeV'])
                            for r in selected], label=f'{energy:g} MeV source')
        axes[1].scatter(x, [r['photons']['relative_count_change_percent'] for r in selected])
    for axis in axes:
        axis.axhline(0, color='black', lw=.8)
        axis.set(xlabel='Fixed same-seed control ID', xticks=[0, 4, 8, 12, 15])
    axes[0].set(ylabel='Current - historical deposited energy / keV', title='(a) Deposited energy')
    axes[0].legend(fontsize=8)
    axes[1].set(ylabel='Arriving photon count change / %', title='(b) All recorded photons')
    fig.suptitle('Figure 0007-120-1: Cross-environment differences (16 engineering controls)')
    path = ART / 'figures/fig_0007_120_1_environment_differences.png'
    fig.savefig(path); plt.close(fig); paths['120'] = path
    with np.load(ART / 'comparisons/events.npz') as data:
        time_us = data['native_time_ns'] / 1000
        fig, axes = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)
        for group, (axis, energy) in enumerate(zip(axes.flat, (1.625, 2.5, 5.486, 8.))):
            inset = axis.inset_axes([.43, .37, .53, .53])
            for label, caption, color, style in [('old', 'Historical 11.4.2', 'tab:blue', '-'),
                                                ('current', 'Current 11.3.2', 'tab:orange', '--')]:
                mean = np.mean([data[f'case{i}_{label}_signal_mV']
                                for i in range(group * 4, group * 4 + 4)], axis=0)
                axis.plot(time_us, mean, style, color=color, lw=1.1, label=caption)
                inset.plot(time_us, mean, style, color=color, lw=1.1)
            axis.set(xlim=(0, 20), xlabel='Time in native window / us', ylabel='Native signal / mV',
                     title=f'({chr(97 + group)}) Source {energy:g} MeV; mean of 4 controls')
            axis.legend(loc='upper left', fontsize=7)
            inset.set(xlim=(4.98, 6.6)); inset.tick_params(labelsize=7)
            inset.set_title('Pulse detail / us', fontsize=8)
        fig.suptitle('Figure 0007-121-1: Complete native SiPM waveforms and pulse detail')
        path = ART / 'figures/fig_0007_121_1_native_waveforms.png'
        fig.savefig(path); plt.close(fig); paths['121'] = path
    return paths


def main():
    verify()
    result = read(EXP / 'result-summary.json')
    replay = read(EXP / 'main-replay.json')
    smoke = read(EXP / 'main-smoke.json')
    assert replay['passed'] and smoke['passed']
    rows = result['rows']
    paths = figures(rows)
    flat = [dict(id=r['id'], source_energy_MeV=r['source_energy_MeV'], seed=r['seed'],
                 old_edep_MeV=r['old']['edep_MeV'], current_edep_MeV=r['current']['edep_MeV'],
                 delta_edep_keV=1000 * (r['current']['edep_MeV'] - r['old']['edep_MeV']),
                 photon_count_change_percent=r['photons']['relative_count_change_percent'],
                 native_signal_rmse_mV=r['native_signal_rmse_mV'],
                 exact_sorted_photons=r['photons']['sorted_exact_equal']) for r in rows]
    write_csv('cross_environment', flat)
    runs = read(EXP / 'manifest.json')['runs']
    inventory = []
    for kind, entries in [('historical_seed_control', runs), ('main_build_repeat', replay['rows']),
                          ('documented_single_alpha', [smoke])]:
        for r in entries:
            inventory.append(dict(role=kind, root=r['root'], macro=r['macro'], seed=r['seed'],
                                  source_energy_MeV=r['source_energy_MeV'], edep_MeV=r['edep_MeV'],
                                  sha256=r['root_sha256'], log=r['log']))
    write_csv('root_catalog', inventory)
    max_energy = max(abs(r['delta_edep_keV']) for r in flat)
    changes = [r['photon_count_change_percent'] for r in flat]
    key = ('当前 Ubuntu/WSL 的 Geant4 11.3.2、ROOT 6.34.10 已完成主工程和实验快照构建。'
           '16条跨环境同种子对照全部正常输出；主工程在当前环境重跑四档各一例，'
           '28组光子/电压数组逐项一致。文档中的单α命令也已通过，共新增21个有效原生ROOT。')
    caution = (f'旧 Fedora 环境为 Geant4 11.4.2、ROOT 6.40.04。16对初始随机状态全部相同，'
               f'但沉积能量、产生光子数、排序后的完整到达光子记录和原生信号均未逐事件精确复现。'
               f'最大沉积能量差为{max_energy:.5f} keV，到达光子数相对差为'
               f'[{min(changes):+.3f}%, {max(changes):+.3f}%]。'
               '这些是同种子工程对照，不能据此断定总体响应有偏差，或把差异单独归因于Geant4版本；'
               '数据集、编译器和平台也发生了变化。')
    next_text = ('用户明确要求使用当前环境，因此后续以本机11.3.2/6.34.10继续。'
                 '336个正式强母＋8个成本先导的草案保留，尚未运行；后续须另建实验快照，'
                 '将问题明确为旧冻结模型在当前环境的新母事件迁移检验。先检查四档孤立响应与读出一致性，'
                 '再根据冻结方案验证逐母投影MSE、均值不确定度和孤立E控制。'
                 '不把跨环境差异当作独立同环境的均值修正收益；不使用已暴露对照调修正。')
    check_table = table(['检查', '结果'], [
        ['冻结G4源码与主工程', '39文件内容完全相同；未改物理参数'],
        ['主工程／快照CTest', '各1/1通过'],
        ['宏路径及光子记录针对性测试', '2项通过'],
        ['跨环境初始随机状态', '16/16相同'],
        ['跨环境E_dep／产生光子数／完整到达记录', '分别0/16逐项相同'],
        ['当前环境主工程重跑', '4/4通过，28组数组完全一致'],
        ['文档单α命令', '1/1通过'],
        ['新ROOT／新增独立光学母样本', '21／0']])
    commands = '''source /opt/root/root63410/bin/thisroot.sh
source /opt/geant4/geant4-v1132/bin/geant4.sh
cmake -S g4 -B g4/build -DCMAKE_BUILD_TYPE=Release \\
  -DGeant4_DIR=/opt/geant4/geant4-v1132/lib/cmake/Geant4 \\
  -DROOT_DIR=/opt/root/root63410/cmake
cmake --build g4/build -j 4
ctest --test-dir g4/build --output-on-failure
cd g4/build
./gagg ../macros/run_single_alpha_sipm.mac'''
    sections = [
        ('119', f'''## 119. 使用当前环境重建Geant4（2026-11-02）

{key}

研究问题是运行环境迁移与可重复性，而非重新选择物理参数。输入为L93训练128母中按固定路径哈希选出的16个历史母（四档源能量1.625/2.5/5.486/8 MeV各4），并复用其实际宏、三组种子和RunInfo。源能量不等于沉积能量；两者逐例分开记录。

旧可执行文件依赖GLIBC_2.38、GLIBCXX_3.4.32、Qt6及旧Geant4库，本机无法加载。复制来的CMakeCache还指向/media/sf_works路径。旧g4/build已完整移入.local_archive/build_migration_20261102/g4_build_from_previous_machine，之后在g4/build重新配置构建。39个物理源码/配置文件与既有冻结实验逐个哈希相同；无物理模型修改。

当前g++为11.4.0；G4EMLOW为8.6.1、PhotonEvaporation为6.1、G4PARTICLEXS为4.1，旧记录分别为8.8、6.1.2、4.2。完整外部数据集及环境在environment.json；这些依赖变化没有与编译器/平台影响单独分离。
'''),
        ('120', f'''## 120. 跨环境同种子对照的适用范围

{caution}

{check_table}

预先冻结的门要求初始随机状态、E_dep、产生光子数及完整到达光子记录逐项相同，是保守的精确重放检查，并非统计等效性检验。该门未通过，不等于当前程序不能运行或已证明物理总体失真。比较还核对了原生时间网格、信号波形及电荷；16条旧ROOT均按哈希保留，实际宏只改变输出/配置文件路径。

图 0007-120-1：(a)逐对E_dep差（当前减旧值，keV）；(b)所有记录面到达光子数的相对差。横轴为冻结对照ID，颜色区分源能量；没有把16对当作新母样本或计算显著性结论。
'''),
        ('121', f'''## 121. 主工程重放与日常运行入口

主工程使用相同当前依赖，按预定ID 0/4/8/12各重跑一例。每例比较四列光子到达时间/位置/能量和三组采样时间、信号电压、含噪电压，合计28组数组逐项完全一致；初始随机状态、初级粒子记录和SiPM标量也相同。ROOT文件自身包含不同路径/时间等记录，不要求整个ROOT二进制哈希相同。

随后直接执行run_single_alpha_sipm.mac，源能量5.486 MeV，G4/雪崩/读出种子261010201/261010301/261010401，得到{smoke['root']}；实际E_dep={smoke['edep_MeV']:.9f} MeV，记录{smoke['photons']}个光子、{smoke['samples']}个原生电压采样。该默认宏也是既有固定种子，只作可用性检查。

从项目根目录执行（其他机器按安装位置调整两个/opt前缀）：

```bash
{commands}
```

每例仅一个20 μs WaveformEvent和SiPMEvent，中央正入射单α在5 μs发射；这是固定单事件，不对应一个泊松注入率。G4内部串行，对照最多4进程，任务启动间隔至少1.05 s；所有ROOT使用根目录g4/data下的秒级文件名且不覆盖旧数据。

图 0007-121-1：(a–d)四档源能量的完整20 μs原生SiPM信号，插图展示脉冲细节。各环境曲线仅为相同4个种子的平均，未画独立统计置信带；横轴是原始时间窗时刻，未对齐首光子。原生SiPM输出与L93的离线数字化分开，本轮没有新增离线数字化或能量拟合。
'''),
        ('122', f'''## 122. 环境迁移记录与下一步

{next_text}

当前结论限于程序可用、输出结构正确及四例同环境重复一致；没有验证跨环境总体统计等效性、盲漏检、误拆分、误合并、能量分辨率、能谱选择效率或高率性能。21个新文件包括16条旧种子跨环境对照、4条同环境主工程重放和1条默认宏检查，新增独立研究母样本计0。

源码、宏、种子、RunInfo和表格在pileup/experiments/{NAME}，root_catalog.csv逐条给出21个ROOT的实际宏与哈希。运行基于Git b2735ea及当时未提交的环境检查代码，protocol.json保留状态，design-freeze.json及后续独立冻结文件记录实际源码。旧结果和原始停步协议完整保留；current-environment-decision.json记载用户要求，environment-status.json汇总21条完成状态。

本次不新增文献性能主张或原PDF核查，沿用L91已有证据边界；环境核验为项目计算结果。L95、研究计划第36节与Run_Catalog同步记录。旧Notebook计算和图像保留，未更换默认重建模型。
''')]
    names = ('codex_0007_SiPM_Energy_Calibration.ipynb', 'Research_Plan_20260921.ipynb', 'Run_Catalog.ipynb')
    notebooks = {}
    for name in names:
        notebooks[name] = read(ROOT / 'jupyter' / name)
        assert notebooks[name] == read(REC / (name + '.before')), name
    assert [len(notebooks[name]['cells']) for name in names] == [193, 104, 65]
    cells = []
    for number, body in sections:
        cells.append(nbformat.v4.new_markdown_cell(body, id=f'transport-environment-{number}'))
        if number in paths:
            source = ('from pathlib import Path\nfrom IPython.display import display, Image\n'
                      "root = Path.cwd()\nif root.name == 'jupyter':\n    root = root.parent\n"
                      f"display(Image(filename=str(root / '{paths[number].relative_to(ROOT).as_posix()}')))\n")
            cells.append(nbformat.v4.new_code_cell(source, id=f'transport-environment-figure-{number}'))
    new_nb = nbformat.v4.new_notebook(cells=cells)
    NotebookClient(new_nb, timeout=120, kernel_name='python3',
                   resources={'metadata': {'path': str(ROOT)}}).execute()
    notebooks[names[0]]['cells'].extend(new_nb['cells'])
    body = '## 36. 当前Geant4环境已跑通，后续按环境迁移检验推进（2026-11-02）\n\n'
    body += key + '\n\n' + caution + '\n\n' + next_text
    body += '\n\n详见codex_0007第119–122节、图0007-120-1与121-1、L95；无新离线数字化或能量拟合。'
    notebooks[names[1]]['cells'].append(nbformat.v4.new_markdown_cell(body, id='transport-environment-plan-20261102'))
    mapping = table(['类型', 'ROOT', '实际宏', 'G4种子', '源能量/MeV', 'E_dep/MeV'], [
        [r['role'], r['root'], r['macro'], r['seed'], r['source_energy_MeV'], f"{r['edep_MeV']:.9f}"]
        for r in inventory])
    body = ('## 当前Ubuntu/WSL环境重建、同种子对照与主工程重放（2026-11-02）\n\n' + key + '\n\n' +
            caution + '\n\n' + mapping + '\n\n' +
            '每ROOT一个20 μs窗、固定单α在5 μs中央正入射、无背景和泊松源率。全部完成且保留，未替换失败种子。'
            '16个原宏只替换路径，雪崩/读出种子仍取原RunInfo；重放4例沿用相应三种种子。'
            '默认宏三种种子261010201/261010301/261010401。所有实际宏、配置JSON、RunInfo及SHA256在实验目录；'
            '39份源码一致、编译与数据集版本见environment.json/main-build.json。'
            'G4串行、批量最多4核，启动间隔≥1.05 s，原始日志在jupyter/run_records/' + NAME + '/logs。'
            'ROOT只存g4/data，图/数组在jupyter/artifacts/' + NAME + '。\n\n' + next_text)
    notebooks[names[2]]['cells'].append(nbformat.v4.new_markdown_cell(body, id='transport-environment-runs-20261102'))
    for name, nb in notebooks.items():
        nbformat.validate(nbformat.from_dict(nb))
        (ROOT / 'jupyter' / name).write_text(json.dumps(nb, ensure_ascii=False, indent=1) + '\n')
    overview = f'''# 当前Geant4环境重建与同种子检查

{key}

{caution}

{check_table}

从项目根目录运行主工程：

```bash
{commands}
```

旧build完整保留在.local_archive/build_migration_20261102/g4_build_from_previous_machine。
当前二进制在g4/build/gagg，冻结实验二进制在本目录bin/gagg，均需在相同外部依赖下使用；构建产物、ROOT及Notebook不随Git分发。

本实验source/g4保留39份实际源码/配置，analysis-source保留检查与报告源码；build.json记录独立快照构建命令。外部依赖本身不包含在快照中。原16条控制的协议、种子及哈希在protocol.json、plan.json、design-freeze.json，main-replay/main-smoke各有单独计划、结果和源码冻结。tables/root_catalog.csv给出21个ROOT与实际宏/日志/哈希，runinfo保存全部RunInfo；同环境比较为数组逐项相同，不要求ROOT文件哈希一致。

```bash
python3 -B pileup/experiments/{NAME}/run.py verify
python3 -B pileup/experiments/{NAME}/validate.py
```

validate.py只读核对已有数据，不重复生产。run.py的preflight/complete用于原16条任务的续跑，已完成时拒绝重复；replay_main.py与smoke_main.py也保护已完成结果。新生产应另建实验与输出，不能重写既有快照。

{next_text}

本阶段为运行环境检查，没有新增独立母样本、离线数字化或能量拟合。旧的跨环境精确复现门仍记为未通过；用户随后指定当前环境，这一选择和21条最终状态分别记录在current-environment-decision.json与environment-status.json。未检验总体分布等效或实际探测器性能。本地codex_0007第119–122节保留两幅固定编号图；L95记录研究影响。
'''
    (EXP / 'README.md').write_text(overview)
    literature = ROOT / 'docs/literature-use.md'
    old_literature = (REC / 'literature-use.md.before').read_text()
    assert literature.read_text() == old_literature and '## L95 ' not in old_literature
    literature.write_text(old_literature + f'''

<a id="transport-environment-20261102"></a>
## L95 — 2026-11-02：当前环境跑通，历史同种子精确复现与新确认分开

**项目结果。** {key}

{caution}

**方法与条件。** 原L93训练128母中按冻结路径哈希选四档各4例，源能量1.625/2.5/5.486/8 MeV，E_dep单独记录。39份物理源码/配置完全相同，实际宏只改输出/配置路径；同G4/雪崩/读出种子。每ROOT一个20 μs窗、一个5 μs中央正入射α；串行G4，批量最多4核，启动≥1.05 s。16对属于既有种子的环境对照，主工程4例和默认宏1例也不计新增独立母。外部G4EMLOW/PhotonEvaporation/PARTICLEXS等版本随环境变化，完整数据集版本及源码哈希保留。

{check_table}

**核验范围。** 主工程重放比较4列光子＋3组电压/时间数组，每例7组、共28组逐项一致，初级粒子/随机状态/SiPM标量也一致。CTest检查SiPM响应；2项Python测试验证宏只改路径和光子记录顺序/位模式判别。全部21条输出及原始16条输入的文件哈希、单窗结构、实际宏和种子核对；旧build完整归档保留。此次没有离线数字化或能量反演，没有评价盲漏检/误拆并、E分辨率、能谱选择效率或高率性能。

**采用与待验证。** {next_text} 跨环境精确门未通过与本机运行成功同时保留，后续需要总体响应和旧模型迁移的独立证据。336正式＋8先导仅为冻结草案，不计已提交运行，也不保证统计精度。

**证据层级。** 本条为项目环境与计算记录；没有新增论文结论、联网检索或原PDF核查，沿用L91已保存的文献证据及审核边界。文献本身未支持本机版本间统计等效性，不能借用其结论作此保证。

位置：`pileup/experiments/{NAME}/`；模块`pileup/analysis/alpha_transport_environment.py`、`alpha_transport_environment_study.py`与`alpha_transport_environment_report.py`。本地codex_0007第119–122节、图0007-120-1与121-1，研究计划第36节、Run_Catalog及ROOT逐条表同步。全部旧笔记计算/图像保留。
''')
    readme = ROOT / 'README.md'
    old_readme = (REC / 'README.md.before').read_text()
    assert readme.read_text() == old_readme
    paragraph = ('需要 CMake、C++17 编译器、带 UI/可视化组件的 Geant4、ROOT 的 Core/RIO/Tree，以及 nlohmann_json ≥ 3.2'
                 '（读取 SiPM 配置；优先使用独立 CMake 包，也支持 ROOT 自带头文件）。既有生产数据使用 Geant4 11.3.2、ROOT 6.34.10；'
                 '2026-10-09 本机使用 Geant4 11.4.2、ROOT 6.40.04，主工程、两个冻结实验和阻止本领工具已重新构建。'
                 '环境检查记录保存在本地 `jupyter/Run_Catalog.ipynb`，历史结果仍按各自的运行记录解释。')
    replacement = ('需要 CMake、C++17 编译器、带 UI/可视化组件的 Geant4、ROOT 的 Core/RIO/Tree，以及 nlohmann_json ≥ 3.2'
                   '（读取 SiPM 配置；优先使用独立 CMake 包，也支持 ROOT 自带头文件）。'
                   '**2026-11-02 当前 Ubuntu/WSL 已用 Geant4 11.3.2、ROOT 6.34.10 和 g++ 11.4.0 重建主工程并跑通单α。** '
                   '2026-10-09 起的部分历史数据来自 Fedora 的 Geant4 11.4.2 / ROOT 6.40.04；历史条件以对应 RunInfo 和记录为准。'
                   '旧构建目录已完整归档至 `.local_archive/build_migration_20261102/g4_build_from_previous_machine`。'
                   '当前环境四例同种子重放的28组光子/电压数组逐项一致，跨环境不能逐事件精确复现；'
                   '详见[迁移验证记录](pileup/experiments/' + NAME + '/README.md)。')
    assert paragraph in old_readme
    updated = old_readme.replace(paragraph, replacement, 1)
    old_command = '''从项目根目录运行，先按本机安装方式加载 Geant4/ROOT 环境：

```bash
cmake -S g4 -B g4/build
cmake --build g4/build -j 2
cd g4/build
./gagg ../macros/run_all_pileup.mac
```'''
    assert old_command in updated
    updated = updated.replace(old_command, '从项目根目录执行以下本机已验证的构建和单α命令；其他机器按安装位置调整 `/opt` 前缀：\n\n```bash\n' + commands + '\n```', 1)
    updated += f'\n\n2026-11-02完成[当前Geant4环境重建与验证](pileup/experiments/{NAME}/README.md)：{key} {caution} 后续按用户要求采用当前环境，独立母事件的冻结模型确认另立环境迁移实验；当前336＋8设计尚未运行。详见[采用记录L95](docs/literature-use.md#transport-environment-20261102)及本地codex_0007第119–122节。\n'
    readme.write_text(updated)
    save(EXP / 'environment-status.json', dict(main_build_working=True, current_same_seed_reproducible=True,
         documented_command_passed=True, native_root_count=21, independent_new_mothers=0,
         old_environment_exact_matching=False, fresh_confirmation_run=False,
         current_environment='Geant4 11.3.2 / ROOT 6.34.10 / g++ 11.4.0, Ubuntu/WSL',
         next='Current-environment response controls and a separate frozen-model transfer experiment',
         historical_decision_record='current-environment-decision.json'))
    outputs = [*paths.values(), EXP / 'README.md', readme, literature, EXP / 'environment-status.json',
               *[ROOT / 'jupyter' / name for name in names], *sorted((EXP / 'tables').glob('*.csv'))]
    save(EXP / 'report-output-freeze.json', {str(p.relative_to(ROOT)): sha(p) for p in outputs})
    print('Wrote 2 figures, sections 119–122, plan 36, 21-ROOT catalog, README and L95', flush=True)


if __name__ == '__main__':
    main()

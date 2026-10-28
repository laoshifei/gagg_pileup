# Physical count decision and ML candidate ablation

在相同缓存拟合候选上分离事件数判决：原生先强后弱200 ns的学习型门槛拒绝13个严格正确K2，改用物理似然判决后严格恢复7/40→20/40；40 ns仍0/40。追加局部精修／网络起点＋物理判决，原生均241/360，高斯分别500/600、501/600。门槛修复通过开发控制，但网络相对局部精修的增益未达两域各5个百分点要求；保持原默认流程。

| 流程 | 原生360对 | 已公开连续Q高斯600对 |
| --- | ---: | ---: |
| baseline | 241/360（66.94%） | 500/600（83.33%） |
| learned_gate | 226/360（62.78%） | 482/600（80.33%） |
| local_physical | 241/360（66.94%） | 500/600（83.33%） |
| mlp_physical | 241/360（66.94%） | 501/600（83.50%） |
| mlp_fixed_gate（门槛敏感性） | 241/360（66.94%） | 501/600（83.50%） |

## 范围与数据边界

纯α，固定双τ 95／351 ns，完整1024×2 ns模拟电压，充分动态范围，未知局部K≤2。原生Q为实现雪崩等效pe，高斯Q是均值／协方差条件参数；没有沉积能标。保持物理响应、噪声模型、K0／K1与源数据不变。

输入来自上一实验alpha_pair_risk_training_v2_20261027：9600训练／960验证所得relative网络及已生成的候选；本轮未重新训练、推断网络、优化或生成样本。280原生门槛帧、760原生评价、1000连续Q高斯评价在本轮均为开发复用。训练／验证／高斯生成种子基数270310001／270310002／270310003；权重选择与母事件分组、输入哈希完整继承于本地input-provenance.json。没有新的G4宏或ROOT。当前代码HEAD及依赖记录于environment.json，未提交的先前研究与本轮源码分别冻结保留。

## 对照与判决

baseline为原候选与原物理门槛；learned_gate为与上轮mlp_nll逐条相同的混合候选及学习型门槛；local_physical从缓存trial0取原K2起点的70步局部精修；mlp_physical额外包含网络起点的70步精修，均保留原初值。两种新物理流程以同280校准帧重新校准同一种NLL差规则；mlp_fixed_gate沿用原数值门槛，仅作敏感性。一个与两个起点不是等预算结构比较。

单事件存在需要Q≥500且NLL1<NLL0。K2需两个Q均≥500、max(NLL1−NLL2,0)>0；观测初拟合单Q按1500／5000／12500分档，加一上尾秩≤.025。并列保留，不读取网络logit、不加Q正则或重复源率先验；阈下最低代价候选不被悄悄换掉。详见protocol.json与gate-calibrations.json。

严格整对正确：每Q±20%、每t±min(10 ns,gap/4)、无额外。普通漏检／额外等使用10 ns时间匹配；时间匹配后的相对Q偏差／sigma68和时间误差、误拆／误并、未收敛均保留。强制K2只诊断候选。预定门槛修复和网络增益要求分开，后者未通过。不改默认算法。

原生先强后弱200 ns：三个库各20/40正确K2，学习型门槛只保留7，物理保留20；40 ns三个库均3/40正确候选且被拒，最终0/40。原生control额外0/0/0/1/0（每组80），高斯local为0/0/0/2/1、mlp为0/0/0/2/0；高斯mid比基线多一帧，未声称全面降低误拆。网络相对local只在高斯长间隔强→弱新增1个成功，不能证明短间隔改善。

配对bootstrap2000次、种子270327062；原生条件内按实际母事件块，高斯独立抽样；已公开开发结果不作物理确认。只有局部K≤2／固定τ／纯α，不提供最高计数率或普适最小间隔。

## 复现与输出

从项目根目录使用项目Python环境运行；设置OPENBLAS_NUM_THREADS=1、OMP_NUM_THREADS=1。run.py使用inference-source，report.py使用report-source，deliver-v2.py使用delivery-source-v2；每层校验冻结哈希，不需要切Git版本。初版交付单次约62 MB写入被共享盘返回长度异常阻断，原笔记未替换；v2只改为256 KiB分块、完整回读后替换，失败源码和临时文件保留，不改变任何计算。

```bash
bash jupyter/run_gagg_python.sh pileup/experiments/alpha_pair_physical_gate_20261027/report.py analyze
bash jupyter/run_gagg_python.sh pileup/experiments/alpha_pair_physical_gate_20261027/report.py figures
bash jupyter/run_gagg_python.sh pileup/experiments/alpha_pair_physical_gate_20261027/deliver-v2.py validate
```

原始prepare／pilot／infer为一次性运行器，拒绝覆盖原目录／缓存。已有缓存可重复analyze／figures；重建缓存时应复制实验及依赖输入到独立项目输出树，保持本记录。Notebook与文档交付亦拒绝覆盖用户后续修改。

本地结果jupyter/artifacts/alpha_pair_physical_gate_20261027/，记录jupyter/run_records/alpha_pair_physical_gate_20261027/；其report／summary／scores／forced-K2／stage-diagnostics／paired-comparison等文件保留全部结果。Notebook0004第93–97节、固定图0004-93-1／93-2／94-1／95-1／96-1，40／200 ns第0对完整0–2048 ns波形，缺失第二输出为零。

6项数值检查、32帧先导、2040帧代价重放（相对差0）、7040主评分／8800含敏感性评分、2880已知K2诊断；旧基线与学习型评分精确复现。五图目视检查、Notebook新单元实际执行，旧455单元原样保留；最终完整校验记录见final-validation.json。输入和各层源码哈希追溯，未收敛候选不删除，推断阻止读取标签。缓存重放11.48 s不是完整拟合吞吐。

## 文献采用与后续

已有Liddick2012 DOI 10.1016/j.nima.2011.12.034（硅DSSD、40 MSPS）仅支持将检测与参数精度分开、控制假双；AI精读待人工审核，门选记号有未决项。Morad2024 DOI 10.1109/TNS.2024.3371573（CLYC、1 GS/s）仅借鉴联合指标，人工审核553为限定确认。当前实验不是其复现，未重查PDF或修改文献库；具体位置、条件、疑问与采用边界见项目docs/literature-use.md的L62。

下一步分析弱Q歧义、证据与区间覆盖，公开误拆—恢复取舍；尚需新母事件、噪声／位置、α/γ、free-τ及高率多事件验证。当前网络不是已证实更好的默认解法。

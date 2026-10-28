# Weak-pulse mean-score diagnostics

Known-count development on 44 previously exposed frozen-model observations. See protocol.json, freeze.json, input provenance, separate primary/oracle predictions and advancement-decision.json. No default change.

40ns的严格整对恢复，原NLL／扩大搜索NLL／固定一次权重GLS／迭代GLS分别为2/16、2/16、2/16、2/16。所有失败与未收敛样本均保留在分母中。

同一组在固定真实时刻／固定真实电荷／从全部真实参数启动后自由拟合时，严格成功分别为4/16、7/16、6/16；这三项使用额外真值信息，只是误差来源诊断。

200ns的严格整对恢复，原NLL／扩大搜索NLL／固定一次权重GLS／迭代GLS分别为10/16、10/16、11/16、11/16。所有失败与未收敛样本均保留在分母中。

同一组在固定真实时刻／固定真实电荷／从全部真实参数启动后自由拟合时，严格成功分别为13/16、14/16、13/16；这三项使用额外真值信息，只是误差来源诊断。

扩大搜索使4/32条的J降低超过0.05；从真值出发仍有4/32条找到比扩大搜索更低超过0.05的J。有限起点和有限预算不能证明全局最优；J更低也不等于Q和t更接近真值。0.05只用于描述数值差异，不是统计显著性阈值。

共有7/32条出现：真值启动后的自由拟合严格成功，扩大搜索候选失败，但成功解的J反而高出0.05以上。其中40ns为4条、200ns为3条。这些反例说明同一代价的排名与本次逐事件误差并不一致；增加搜索起点无法保证选中更准确的那个解，也不能据此判定物理上一定不可分辨。

没有方法同时通过预先冻结的四项推进条件，因此本轮不增加112条确认样本、不改变默认重建流程。这是44条已暴露模型样本上的开发筛查结论，不是物理不可分辨极限，也没有证明其他优化器无效。

Reproduction: run.py tests; run.py pilot --workers 2; run.py infer --workers 2; run.py oracle --workers 2. Existing outputs are protected; use an isolated copy with the recorded relative input tree for a clean rerun. report.py regenerates scoring/figures from saved predictions; delivery.py appends to the preserved 570-cell baseline and checks it before writing. Models/source are archived; switching Git revisions is unnecessary.

The interrupted pre-inference Notebook backup remains under the sibling *_prepare_interrupted_005415 run record. The original Notebook was unchanged. No new ROOT or random draws were generated.

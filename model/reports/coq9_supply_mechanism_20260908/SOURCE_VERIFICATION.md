# Source verification — CoQ9 static supply mechanism

核验范围仅限指定的本地一级技术来源：固定 model.xml、runner/helper/context 源码、两条旧 t=0 trajectory 和本轮输出。未使用网络或新增求解。

## 覆盖

- Claims checked: 26/26
- Supported: 21
- Partially supported with explicit limits: 5
- Unsupported / contradicted / unchecked: 0

## 核心来源结论

1. model.xml 的 Q9/Q9H2 两行相加后，氧化还原循环项互相抵消；R385 是原模型唯一总 Q 净生成项。runtime helper 增加 source `+1` 与 dilution `-1`，并约束 dilution=`alpha*growth`。所以 `R385 + source = dilution = alpha*growth` 是模型结构事实，不是原生生物学验证。
2. Condition 05 的 source flux 明显高于 FeasibilityTol，且与 imposed demand 的残差为 `3.59e-15`。它仅证明人工 source 能在这个静态 t=0 优化问题中替代 R385 净供给，不证明真实 Q9 pool、耗尽、周转或长期 rescue。
3. 当前 GPR 中，YALI1A14736g KO 只关闭 R305；YALI1A21711g KO 只关闭 R2062。对应 reaction KOs 重现旧 gene-KO t=0 growth。
4. R1889 关闭后的 R1977/R740 通量属于一条 pFBA 解；不能升级为唯一性或 FVA 必需性结论。

## 限定项

- Condition 02 / R558 的 bound exception 使严格全-bound gate 失败；依赖该条件的 R1889 边际量化应标为 provisional，但 Q9 balance 与定性机制比较仍可报告。
- 两个基因在本地模型中可能带 accession/EC 等数据库注释；这里的“native function uncharacterized/unverified”指没有原生实验功能验证，不应解释为没有任何数据库注释。
- restoration fingerprint 覆盖 reaction IDs、stoichiometry、bounds 和 GPR，但不覆盖 objective 及所有额外 solver constraints；“恢复”声明限于已记录字段。
- manifest 的 `started_utc` 在求解结束后写入，不是精确开始时间。

状态：内部 runtime-only 机制判断 `CONDITIONAL PASS`；严格数值全通过 `FAIL`；生物学发表与正式模型/GPR curation 均 `NOT CLEARED`。

# F38820 候选功能元数据：独立来源核查

核查时间：2026-09-11T00:51:03.395165+00:00。审阅者：`/root/audit_gpr`。

**结论：本轮候选文字与数字符合已审原始来源，可用于用户授权的候选功能元数据标记。** 本核查不采纳 GPR、原生定位或 essentiality 变更；此前 33 项研究审计分母不重算。本轮没有新 BLAST、结构搜索、GEM 求解或模型数据修改。

目标 **YALI1F38820g** 的原生正式基因名尚未核实；可标记为“偏 Vph1-like 的 V0 a 亚基候选，仍需实验确认”，功能为质子转运／复合体装配候选。其功能判断基于序列同源及 **AlphaFold 预测**；不是目标蛋白的实验功能确认。

## 实际核查

1. **候选状态与边界。** 直接读取新增 `gene_function_annotations.json`：名称明示 candidate／experimental confirmation required；notes 保留 provisional、正式名未建立、实验确认待完成及 W29 液泡与 Golgi／内体定位未确认。另明确一对一同源关系、功能互换和生长必需性未确认。未将 VPH1 作为已确立的 Yarrowia 正式基因名。
2. **序列身份。** 直接读原始 GenPept AOW07942.1 的 VERSION、source 和 ORIGIN：CLIB89(W29)，804 aa，完整序列 SHA 为 `bbe56eb89a0ab7edad9e0e65450a1574c7dee90acd6695134b4008b6026a0083`，与本轮声明及此前固定 A0A1H6PMT1 v1 序列一致。GenPept 标记 conceptual translation，不能据此称原生功能实验。
3. **BLAST 数字。** 直接读原 XML 的目标查询及 P32563／P37296 HSP。**YOR270C／VPH1／P32563**（酿酒酵母实验确认的液泡 V0 a 亚基）为 446/832＝53.61% identity，查询范围9–798，覆盖790/804＝98.26%。**YMR054W／STV1／P37296**（实验确认的 Golgi／内体 V0 a 亚基）有两个 HSP：297/681 与87/116，查询范围9–606及686–798，两者不重叠；合并已报告比对片段为384/797＝48.18%，查询覆盖711/804＝88.43%。43.61% 只对应 Stv1 第一段，75% 只对应第二段；48.18% 也不是全长逐位 identity。新数据保存分段原始分子／分母，没有混用这些量。
4. **结构数字。** 直接读已有 US-align 原始输出：目标归一化 TM-score 对 Vph1／Stv1 为0.84079／0.83102；参考坐标长度归一化为0.91765／0.93202；对齐715／710残基，RMSD2.51／2.54 Å。新数据逐项一致。归一化后排序会改变，两个结构都相似，支持保留“区室仍未确定”的限制，不能把小幅结构分数差变成功能定位结论。这里 RMSD 指原始程序打印的最小二乘 RMSD，未混用导出 TM 优化变换的 RMSD。
5. **应用步骤的静态边界。** 只读检查 `apply_curated_gene_function_annotations`：先核版本、目标、旧显示名、预期 UniProt 与必要候选限制；通过后仅赋值 `gene.name` 和 `gene.notes`。所查函数没有修改 gene ID、身份 annotation、反应、GPR、区室或 essentiality。此为静态检查；实际构建、导出、幂等性及全模型差异验证由根任务完成，不在本来源审计中冒充已执行。

“序列比较偏向 Vph1；预测结构同时接近两类；原生定位及具体生理功能尚需实验确认”是当前允许的解释。未发现需修订的源级文字或数字错误。

<details>
<summary>所读文件身份与定位</summary>

原 XML 定位：查询 `YALI1F38820g`，P32563 单 HSP、P37296 两个 HSP。结构定位：两份 `target_a` 原输出的长度、Aligned length 与 TM-score 行。GenPept 定位：`VERSION AOW07942.1` 记录。下面 SHA 绑定本次实际读取版本；后续代码或数据变化需由执行记录另行标明。

- `data/reference_build/curation/gene_function_annotations.json`
  - SHA256: `a85a9acebc2c316e6e3133995968ce87d6e7362ccac4f13faade771ce79afd48`
- `scripts/gem_annotate/genes.py`
  - SHA256: `01981083d18b23d1a6e8c428d793ca77348eb323b7b759fa2e30763211c9fce0`
- `artifacts/r794_r795_gpr_review_20260909/sources/ncbi_target_proteins.gp`
  - SHA256: `e105e1b3ab2be40f1aff21760f7eb363388e6d902152956358c605d85c2a9a47`
- `artifacts/r794_r795_gpr_review_20260909/continuation_20260910/sources/old_rid_response.xml`
  - SHA256: `ec0b8f34d2f68f51b327a1481227e46b4c803a8f67c0b54a3b1288bba77977d7`
- `artifacts/r794_r795_gpr_review_20260909/continuation_20260910/structure/results/target_a__6o7t_a.txt`
  - SHA256: `afec9318d57137964e48311458644a049282eba785b11dc9dac8339ff83f2bd9`
- `artifacts/r794_r795_gpr_review_20260909/continuation_20260910/structure/results/target_a__6o7u_a.txt`
  - SHA256: `0ead9d9b97508b06feec3fef60e9746e59df3c66eb942849eccec33d961ff8bc`

</details>

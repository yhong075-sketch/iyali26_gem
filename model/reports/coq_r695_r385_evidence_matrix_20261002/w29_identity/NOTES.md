# W29 身份核验补充

本轮只解析保存的序列、原始 BLAST 表、AlphaFold 元数据/PDB/PAE 和当前候选 XML，并重新计算哈希、覆盖率与置信度摘要。没有重新运行 BLAST、结构比对、AlphaFold、构建或求解。具体核验时间与完整身份见 [identity.json](identity.json)。23/23 个有界检查通过；这不是 23 个生物学命题被实验证明。

| 系统 ID — 名称/功能与证据 | 已核实 W29 序列 | 历史 BLAST（酿酒酵母参考；一致性/目标覆盖/参考覆盖） | 当前候选中的作用 |
|---|---|---|---|
| YALI1E18269g — 原生正式名称未核实；COQ7 家族去甲氧基泛醌羟化酶候选（项目整理注释；底层 W29 记录为同源推断） | XP_503973.1，192 aa，CLIB89(W29) | P41735 / CAT5：58.421% / 97.396% / 81.545%；E=2.44e-75 | R695 的催化候选，与 COQ9 候选组成条件性功能 AND |
| YALI1F34675g — 原生正式名称未核实；COQ9 家族脂质结合/底物呈递辅助候选（项目整理注释；底层 W29 记录为同源推断） | XP_505941.3，258 aa，CLIB89(W29) | Q05779 / COQ9：23.985% / 96.899% / 98.846%；E=9.46e-21 | R695 辅助依赖；不表示其催化氧插入 |
| YALI1B20835g — 原生正式名称未核实；COQ3 家族 CoQ O-甲基转移酶候选（项目整理注释；底层 W29 记录为同源推断） | AOW01767.1，367 aa，CLIB89(W29) | P27680 / COQ3：47.143% / 73.025% / 89.103%；E=1.46e-82 | R385 保留单基因；同一候选也关联 R715 |

三条序列均重新确认 GenBank ORIGIN、原生 FASTA、缓存 UniProt、AlphaFold 元数据和 PDB CA 残基序列完全一致。BLAST 是 2026-09-18 的 BLASTP 2.17.0+、9 候选对 11 预选 reviewed 参考面板；本轮只核实保存对齐片段、匹配数和边界，未执行新检索。覆盖率按单 HSP 的未加 gap 序列边界除以各自全长计算，不能拿对齐长度直接相除。该面板不排除未检查的旁系蛋白。

## 必须保留的身份缺口

**COQ3 当前模型/数据中的 XP_500950.3 未被这套历史序列证据验证。** 历史 BLAST 与 AlphaFold 均绑定 AOW01767.1；当前 XML 的 geneProduct 和 `data/coq9_gene_evidence.tsv` 列 XP_500950.3，同时含 Q6CEG2。本轮检查范围内没有二者与已验证 W29 序列完全一致的桥接材料，不能替换 accession 或把它们写成已验证同一序列。AOW01767.1 文件及历史 manifest 的哈希一致，但 identity_retrieval.json 未含其精确获取时刻，保持未知。该缺口限制当前 RefSeq 注释的身份声明，不阻塞以 AOW01767.1 为明确对象报告 COQ3 家族候选。

三个 W29 UniProt 缓存条目原本均为 unreviewed、protein existence=“Inferred from homology”；2026-09-18 保存的当前查询显示因“不属于 reference proteome”而停用。这不等于功能被反证，也不能称为 reviewed 的原生实验注释。W29 一致性不自动认证 PO1f 序列、定位或原生酶活。

## AlphaFold 与 GPR 的外推上限

复用的是 2022-06-01 创建、2026-09-18 获取的既有 **AlphaFold 预测**：AlphaFold Monomer v2.0 pipeline，AFDB 文件 release 6。三者完整序列匹配；COQ7/COQ9/COQ3 平均 pLDDT 分别 93.44/77.31/67.83，平均 PAE 3.85/11.77/16.04 Å，pLDDT<50 占 1.04%/16.67%/31.34%。COQ3 的 AFDB 原始 globalMetricValue 为 67.81，与 PDB CA 均值 67.827 小有差异，均原样保留，未混写。

历史 COQ7/COQ9 单链与人源 7SSS 同家族链的 USalign TM-score（按目标全长/实验有坐标链长度）分别 0.82825/0.91435 与 0.64390/0.90159；实验链为 E=COQ7、A=COQ9。数字只是读取历史结果，本轮未重新拟合。高参考归一化得分不代表完整目标蛋白均被验证，单链叠合也不验证结合界面。COQ3 没有本子任务所查材料中的对应实验结构比较。

当前候选 XML 确认 `R695 = YALI1E18269g and YALI1F34675g`、`R385 = YALI1B20835g`。R695 AND 的前提是原生膜内长侧链底物的可及性需要辅助呈递；不是两个酶都催化羟化，也不是原生缺失时通量为零的实验证明。既有 Nicoll 2024 审核记录的祖先 COQ7 在非异戊二烯化底物/0.05% DDM 中单独活跃、加入 COQ9 后效率约增加 1.5 倍，反驳普遍的本征催化必需性，但不证明 W29 长链膜底物下 COQ9 可省略。此反证从既有原始来源审核继承，本身份子任务未重审论文。R385 原生底物氧化态及两端连接不能由 COQ3 同源或 AlphaFold 自动补齐。

## 已直接读取的证据

- `artifacts/coq_system_gpr_review_20260918/sources/XP_503973.1.gb`、`XP_505941.3.gb`、`AOW01767.1.gb`：VERSION、strain、locus_tag、ORIGIN。
- 同目录 `native_sequences.fasta`、`reference_panel.fasta`、`native_vs_reviewed.blast.tsv`、`native_cached_records.json`、`sources/Sc_COQ_references.json`：原始序列与对齐。
- 同目录 `sources/A0A1H6Q2E8_af.json`、`A0A1D8NQ60_af.json`、`A0A1D8N802_af.json` 及对应 v6 PDB/PAE：既有预测身份与置信度。
- 同目录 `coq79_structure_comparison.json`、四份 `YALI1*_vs_7SSS_*.txt`、`structure_audit.json`：历史结构分数与独立审核范围。
- `artifacts/coq_gpr_completion_20261001/coq89_audit/REPORT_zh.md`、`source_claim_audit.json`、`native_identity_reuse.json` 及 `gene_decisions.tsv`：功能 AND、反证与历史身份限定。
- `artifacts/coq_literature_pipeline_20261002/coq_literature_pipeline.xml`：实际候选 GPR、notes 与 geneProduct 交叉引用；`data/coq9_gene_evidence.tsv`：当前整理注释范围。

完整输入文件 SHA、版本、原始 BLAST 参数与本次检查项保存在 identity.json；未改上述输入文件。

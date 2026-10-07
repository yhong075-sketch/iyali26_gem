# Escher 基因显示来源调查（只读）

核查时间：2026-09-09T21:59:25.132840+00:00。核查范围：静态参考 1074 个基因、现有注释来源、当前 core 图 GPR 可读性。没有联网、文献重审、序列比对、结构计算、求解或模型/GPR修改。

## 结论

- 静态参考 1074 个 geneProduct 中，1054 个 name 为 `COBRAProtein…`，15 个为 `G_YALI…`；只有 5 个非占位名称。`COBRAProtein` 不可作为已核实名显示，模型 EC/UniProt 交叉引用本身也不等于本次核实的蛋白功能。
- `scenario_0.json.gz` 有 1075 个基因，额外一个来自运行时质粒互补覆盖；本图 1074 基因身份应从 reference.xml 取得，不能把运行时覆盖混入静态基因总数。
- 现成最容易安全复用的是 `data/coq9_gene_evidence.tsv`：48 行，13 行带限定的整理注释，35 行明确为 model/GPR assignment only。保留 `evidence_scope`，不要把候选名称、跨版本桥接和原生位点验证合并。
- 本地 UniProt 快照有广覆盖的候选功能，但大部分未审校，不能为了填空升级为已核实功能。无足够功能证据时显示“无已核实名称；功能未核实；模型角色：所连反应；model/GPR assignment only”。

## 可复用字段与可靠程度

| 来源 | 可读字段 | 安全复用方式 |
|---|---|---|
| reference.xml | gene ID、原始 name、annotation URI、GPR | ID 与 GPR 直接复用；占位 name 屏蔽；5 个非占位名字仍保留“模型所载名称”身份。 |
| data/coq9_gene_evidence.tsv | yali1_gene、established_name、protein_function、evidence_status、evidence_scope、sources、accession | 以系统 ID 精确连接；优先作为带证据范围的显示层。assigned_role/established_name 有时是功能描述，不必硬叫基因符号。 |
| coq9_pipeline…/gene_evidence_matrix.tsv | gene_id、name_or_symbol、function、evidence_status、literature_review_outcome | 更明确保留“no established Yarrowia gene name; COQn candidate”。fitness_summary 是已整理的分实验结果，不能转为本次独立 essential 真值。 |
| uniprot_UP000182444.json | primaryAccession、entryType、genes.orfNames、proteinDescription、comments、evidenceCode、entryAudit | 7894 条均为 UniProt TrEMBL 未审校记录；仅去除 YALI1 与下划线的格式差异，可唯一精确匹配 1058/1074 位点。可作“未审校数据库功能候选”，不可自动写成 curated annotation 或 experimentally verified。 |
| uniprot_UP000001300.json | 同上 | 6454 条，667 条 reviewed。对参考 ID 不跨版本直接匹配 5 条，其中 reviewed 2 条。已有 safe crosswalk 可找 123 个 reviewed 记录，但只是已有跨版本身份桥接，不能一键认证 123 个 W29 位点功能。 |
| gene_locus_tag_map.json | _meta.crosswalk_fingerprint、_meta.model_gene_fingerprint、lookup | 属于已有映射缓存；格式归一化与跨版本映射分开记录。需遵守 identity exclusions；不能仅凭相同数字后缀或模型现有 accession 绕过排除。 |
| curated_gene_annotation_overrides.csv | gene_id、UniProt/KEGG/NCBI/RefSeq/EC、case_id、evidence_path | 当前仅 1 条身份纠正，不是通用功能表；只可按对应 case 复用。 |
| data/kegg/yli_genes.tsv | KEGG gene ID、特征类别、坐标、YALI2 标签与描述 | 含 YALI2 新版本标识；不要通过 YALI1 数字后缀直接连接。需要现成已核实跨版本桥接，否则本轮留空。 |
| batch_annotation_results.csv / annotation_audit.csv | 反应自动/LLM 注释 | 属于反应注释，不能作为基因名称/原生蛋白功能证据。 |

## 可直接保留的限定注释条目

下表仅复用现有证据整理，不表示本次重新验证蛋白身份或生化活性。全部默认 `curated annotation` 并显示原 evidence_scope；其中功能候选仍写“候选”。

| 系统 ID | 名称或名称状态 | 简要蛋白功能 | 证据范围／重要限定 |
|---|---|---|---|
| YALI1A08781g | COQ6 候选；原生标准名未在本轮重新核实 | candidate flavin-dependent CoQ ring hydroxylase; native electron donor unresolved | Existing curated candidate/annotation and local crosswalk; pathway role checked where discussed, not new native-locus enzymology or sequence identity certification |
| YALI1A14736g | 无已核实标准名（Complex III core/MPP-like 成员注释） | Complex III core/MPP-like member | Native Yarrowia CIII structural membership via local YALI1-to-YALI0 crosswalk and PDB8ABF entity5; no new sequence alignment or single-KO essentiality validation; 本地跨版本映射支持的原生CIII成员；不等同于原生单基因绝对必需 |
| YALI1A21711g | NUPM 候选；accession 桥接未闭合 | NUPM candidate; exact accession bridge not newly closed | Local crosswalk and archived CI structural candidate; PDB6YJ4 directly confirms Q6CGB4 NUPM, not the entire YALI1-to-accession bridge |
| YALI1B20527g | COQ8 候选；原生标准名未在本轮重新核实 | candidate accessory ATPase/kinase-like protein and CoQ complex stability | Existing curated candidate/annotation and local crosswalk; pathway role checked where discussed, not new native-locus enzymology or sequence identity certification |
| YALI1B20835g | COQ3 候选；原生标准名未在本轮重新核实 | candidate CoQ O-methyltransferase | Existing curated candidate/annotation and local crosswalk; pathway role checked where discussed, not new native-locus enzymology or sequence identity certification |
| YALI1C25352g | COQ5 候选；原生标准名未在本轮重新核实 | candidate CoQ ring C-methyltransferase | Existing curated candidate/annotation and local crosswalk; pathway role checked where discussed, not new native-locus enzymology or sequence identity certification |
| YALI1C26017g | COQ1 候选；原生标准名未在本轮重新核实 | candidate polyprenyl diphosphate synthase | Existing curated candidate/annotation and local crosswalk; pathway role checked where discussed, not new native-locus enzymology or sequence identity certification |
| YALI1D11769g | 无已核实标准名（细胞色素 c 携带体注释） | Cytochrome c carrier, not cytochrome c1 | Model-to-accession plus official cytochrome-c family annotation; enzyme-subunit vs carrier distinction |
| YALI1E18269g | COQ7 候选；原生标准名未在本轮重新核实 | candidate demethoxyubiquinone hydroxylase | Existing curated candidate/annotation and local crosswalk; pathway role checked where discussed, not new native-locus enzymology or sequence identity certification |
| YALI1F08349g | COQ2 候选；原生标准名未在本轮重新核实 | candidate 4-hydroxybenzoate polyprenyltransferase | Existing curated candidate/annotation and local crosswalk; pathway role checked where discussed, not new native-locus enzymology or sequence identity certification |
| YALI1F32476g | NDH2（整理注释） | NDH2: external alternative NADH dehydrogenase | Model/local locus mapping; native species NDH2 external orientation and retargeting experiments |
| YALI1F34625g | COQ4 候选；原生标准名未在本轮重新核实 | candidate CoQ assembly/decarboxylation role; native mechanism unresolved | Existing curated candidate/annotation and local crosswalk; pathway role checked where discussed, not new native-locus enzymology or sequence identity certification |
| YALI1F34675g | COQ9 候选；原生标准名未在本轮重新核实 | candidate lipid binding/substrate presentation and COQ7 accessory function | Existing curated candidate/annotation and local crosswalk; pathway role checked where discussed, not new native-locus enzymology or sequence identity certification |

当前 core 图与这 13 条注释的交集仅有上表的 Complex III core/MPP-like 成员、细胞色素 c 携带体、NDH2 3 个位点。其余基因不能仅因连到熟悉酶反应而赋予已验证蛋白功能。

## 已有证据冲突必须保留

| 系统 ID | 名称 | 功能 | 显示结论 |
|---|---|---|---|
| YALI1C24124g | ICL1 候选 | 异柠檬酸裂解酶候选 | 当前 W29 记录 A0A1D8NBI0 为 458 aa；已整理的 reviewed P41555 为 YALI0C16885g、540 aa。现有独立审阅明确把当前 458 aa 功能维持为 model/GPR assignment only。不可因跨版本表能连到 P41555，就继承为实验验证的 ICL1。 |
| YALI1E18171g | 无已核实名称 | 功能未核实 | 已有 identity exclusion 指出当前 W29 ORF 与所重叠 CLIB122 ORF 不是同一蛋白；禁止经该 crosswalk 传播名称或功能。 |

## 当前 core GPR 与压缩显示建议

- 当前 4 张 core 图去重 47 条反应，121 个不同基因；分类：单基因 20、纯 OR 13、无 GPR 2、纯 AND 11、混合 AND/OR 1。此统计直接来自参考 XML。
- 最长是 R851：23 个 OR 候选；R171：14 基因 AND 复合项 OR 两个单基因项，合计 16 基因；R304：12 基因 AND；R305：11 基因 AND。图上强行展开会遮挡通路。
- R1311 和 R1889 原始 GPR 为空，显示“无 GPR”；不能给 R1889 移植另一个复合体反应的 GPR。
- R2206 有 3 个 AND 基因，且含非标准 YALI 前缀 ID。解析应依照 geneProductRef 或全量 model.genes 精确词元；只匹配 YALI 前缀会漏基因。
- 画布反应标记可为“R851 · 23 基因 OR”“R305 · 11 基因 AND”“R171 · 16 基因 / 混合规则”。这只是索引摘要，不是等价缩写公式。点击后呈现完整原始带括号 GPR，并按原始 AND/OR 分组逐项列出所有系统 ID。
- 单基因、双基因可直接显示系统 ID；常用名作为补充，不能替换唯一系统 ID。各基因详细表按“系统 ID—名称/未核实—功能/未核实—证据级别—模型所连反应”呈现。
- 图上 AND 表示模型要求共同满足的布尔条件；OR 表示模型视为替代项，不能由此宣称实验上可互换。混合规则绝不能压成全部 AND、全部 OR 或无括号列表。
- 用反应颜色显示条件化状态时，基因数量摘要不能改变颜色判定逻辑。inactive、essential 要分别保留其原始实验/模型判定来源、培养条件、阈值与未知值；注释文件本身不提供当前模型全基因标签。
- 特别避免由“AND 中某基因 essential”推出全反应 essential，或由“OR 项里有 essential 基因”推出另一 OR 项应删除；这些都是未授权科学改写。

## 后台来源身份

以下为本次只读源文件 SHA256，便于后续记录显示来源版本。

| 路径 | SHA256 |
|---|---|
| `artifacts/r608_engineering_20260907/inputs/reference.xml` | `bc2aac8fecd8f2f5f20de7bb3c988bf46b3a5831e525f556498ed51159bc1bee` |
| `data/coq9_gene_evidence.tsv` | `4318e3212f548582e053f81f02cd46dfa6e3892a2587251c252993d21cd9e8c0` |
| `artifacts/coq9_pipeline_20260909/input/evidence/gene_evidence_matrix.tsv` | `10e4b68498849c679240a2930908ab5ab5274827c58420917d98b26f212aab59` |
| `artifacts/r608_engineering_20260907/research/cache/data/uniprot_UP000182444.json` | `e5b0a04874079b4057ffe25dadcb6b812cba8c96b227f41187187b27715753a3` |
| `artifacts/r608_engineering_20260907/research/cache/data/uniprot_UP000001300.json` | `89f2dc9c77ff1aa5f1b4c245567285e5e98bebc0916ab94ed25b43e6210a3a4b` |
| `artifacts/r608_engineering_20260907/research/cache/data/gene_locus_tag_map.json` | `fee341c644fac7a50562d846d0c34fe6608e883d5178fb55422ce42f6181a8c7` |
| `artifacts/r608_engineering_20260907/research/state/essentiality/repository/curated_locus_identity_exclusions.csv` | `883f1aac417838e923f12a00285545a093b985e90fd87ecfd9e6f29adbeee408` |
| `artifacts/r608_engineering_20260907/research/state/essentiality/repository/curated_gene_annotation_overrides.csv` | `ca375a30afb96200da86bdb7a3de68c5169b14adc786920e0cc709374d1583b6` |
| `artifacts/next_evidence_independent_audit_20260908/REPORT.md` | `a4b2fa3cbcbb60eaf8f812246e9c40269d6f54445cf1f04d5cb971105301cbec` |
| `data/kegg/yli_genes.tsv` | `33cb8ae3689f9774f66e2c22e4c3caeab9fa55db89941fd4aa791be5c0d5bf4f` |

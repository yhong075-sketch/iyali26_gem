# R2251 / R2087 独立来源审计

审计日期：2026-09-11。审计者独立读取原模型、封存反应快照及 WT 通量、旧筛查目标行、原始 UniProt 缓存、UniSave 历史记录、当前 UniProt 返回、AlphaFold API/PDB/PAE 和 PMC7003347 正文；没有依据主报告的结论进行投票。新增优化求解、BLAST、结构预测或集群作业均为 0。此文件只记录本次范围，不接受科学模型变更。

## 蛋白身份及审计口径

| W29 系统 ID | 名称与蛋白功能 | 证据边界 |
|---|---|---|
| YALI1E14619g | 原生正式名未核实；amidase signature 家族酰胺水解酶候选；A0A1D8NI23，540 aa | 未审校数据库/同源预测；具体苯乙酰胺活性仅模型赋值 |
| YALI1E41276g | W29 原生正式名未直接核实；对应文献 YlAMD1 的乙酰胺利用酰胺酶候选；A0A1D8NL97，539 aa | 数据库序列与 CLIB122 YALI0E34771g 对应蛋白完全一致；论文的 YlAMD1 遗传实验在 YB-392，不能直接称 W29 实验确认 |
| YALI1F21875g | 原生正式名未核实；amidase 候选，泛 EC 3.5.1.4；A0A1D8NNQ1，552 aa | 未审校数据库自动注释/同源预测；具体苯乙酰胺活性仅模型赋值 |

文献涉及的 YALI0E34771g — YlAMD1 — 主要乙酰胺利用功能有 YB-392 删除/回补遗传学支持；YALI0E11847g — 正式名未核实 — amidase 家族候选，在该论文条件下删除未出现同类明显缺陷，不能推断其完全无酶活。下表“支持”限于各条明确的声明层级。

## 审计覆盖

`total claims 19 | audited 19 | supported 15 | unresolved 4 | contradicted 0 | unchecked 0`

4 项未决是 B7–B10 的科学证据缺口；“已审计”表示检查了来源及适用范围，不表示 19 项科学功能全部已确认。已完成主报告最终措辞核对，详见本文末；报告补充说明按同一 19 项口径审计，不重复扩大分母。

| 编号 | 被审计声明 | 直接打开的来源与定位 | 判定与限制 |
|---|---|---|---|
| M1 | 两条是同一胞质计量列，均为正向 [0,1000] | `model_metadata_trna.xml` 的 R_R2087/R_R2251；旧 `reaction_snapshot.json` 同名键 | supported；两条均消耗 m32/m1931、生成 m38/m1932，各系数绝对值 1。仅是模型的物种身份/计量，非完整化学式已验证 |
| M2 | R2251 为 E14619 单基因；R2087 为 E41276 OR F21875 | XML `geneProductAssociation` 及旧快照 `gpr` | supported；独立 XML 递归布尔检查与主代理 8 状态真值表一致 |
| M3 | 完整底物/产物行强制两条均零 | 独立从 XML 全 2315 反应重建；`compartment_results.json` 的 m1931/m1932、round 1 certificates | supported；底物行只有 −v2087−v2251=0，产物行只有 v2087+v2251=0；两条非负，所以分别为零。R144/R2252/R2088 也分别由仅有 −1 消耗项的 m298/m1933/m2052 完整行强制零。无需 WT 数值或新求解，不推广到改连接/边界后的模型 |
| M4 | 旧 WT 完整见证仍可行且生长为 1.8718823069403 | 旧 `results.json#/runs/WT`；旧快照全计量/边界；XML biomass_C objective | supported；独立重算 1877 行最大残差 8.058630802595945e−15，最大边界违反 0。原 WT optimal 是历史记录，本次未新增最优性证明 |
| M5 | 三项单敲、三基因联合关闭以及两反应关闭均保留旧 WT 见证 | 独立 XML GPR 求值；旧 WT 的全部受影响反应通量 | supported；E14619 单敲关闭 R2251/R2252，另两单敲不新关闭反应；三者联合关闭 R2087/R2088/R2251/R2252，这四项旧通量均为零。关闭两主反应同理，不需要修改通量 |
| M6 | 历史筛查三目标均预测非必需且标签为 unlabelled | `screen_test_metadata_trna_20260910/screen_predictions.tsv` 和 `essentiality_per_gene.tsv` 的三目标行 | supported；三者 KO/WT≈1，四个保存阈值均非必需；不能计为实验真阴性或已证实的 FN |
| M7 | 当前保存的两个原始输入和旧目录模型已有同样规则/计量 | 直接打开 `data/iyali26.xml`、`data/iyli21.xml`、`model.xml` 的两条 reaction | supported；只证明现存输入已有，不证明最早引入日期、原作者依据或修复动机 |
| C1 | 当前化学元数据不足以通过元素/电荷验收 | 直接打开 XML 的 m1931/m1932/m32/m38 species 属性 | supported；m1931 chemicalFormula 为字符串 Phenylacetamide，m1932 缺公式，m38 为 H3N/charge+1。按保存电荷计算反应净+1；这只是注释状态，不以解析出的伪元素 Ph 断言真实元素不平衡或推定修复方式 |
| B1 | 三条 W29 序列可由确切 locus 及版本定位 | 原始 `uniprot_UP000182444.json`；本轮三份 `_cached.json` 与 `_v23/v28/v29.txt` | supported；缓存提取与原缓存对象完全相同，UniSave 历史序列与缓存逐字相同；不以 UniProt 其他菌株交叉引用代替 W29 locus |
| B2 | 三条当前 UniProt 记录已 inactive | 三份 `_current.json` 的 `inactiveReason` | supported；均显示 DELETED / Not part of a reference proteome。这是数据库状态，不是蛋白不存在或序列作废的证据 |
| B3 | 三者有泛 amidase 家族证据，NNQ1 另有泛 EC3.5.1.4 注释 | 三份历史/缓存记录的 comments、features、ECO、GO | supported；TrEMBL，蛋白存在等级 Inferred from homology。ARBA/PIRSR/Pfam 及自动 GO 不等于目标酶学实验；未见原生细胞定位条目 |
| B4 | W29–CLIB122 桥接需分别表述 | 三份目标缓存与 Q6C676/Q6C3H6/Q6C1G7 当前 JSON 的完整序列 | supported；NI23–Q6C676 为 539/540、W29 第342位 Y 对 CLIB122 N；NL97–Q6C3H6 为539/539；NNQ1–Q6C1G7 为552/552。前者不能写完全一致 |
| B5 | 论文支持 YlAMD1 在 YB-392 的主要乙酰胺利用功能 | `PMC7003347.xml`，DOI 10.1186/s12934-020-1292-9，Par9–13、Par18，Fig1/4/6 图注 | supported；2.3 g/L乙酰胺平板等条件，删除与质粒回补支持；Par12 仍有少量慢生长且保留杂质/次要路径解释，因此不能写“完全不能生长”或绝对必需 |
| B6 | 论文中 E11847 删除没有同类明显乙酰胺生长缺陷 | 同文 Par11–12；NS996 对照 | supported；只限试验菌株/培养条件。它不是对该蛋白所有底物活性为零的证明；该文未检验本次第三目标的独立苯乙酰胺活性 |
| B7 | 论文 YB-392 克隆与 W29 NL97 是同一完整蛋白序列 | 同文 Par18；W29 与 CLIB122 记录 | unverified；论文说克隆来自 YB-392 DNA，当前审计未取得并核对该论文克隆的完整序列。W29=CLIB122 不能补上 YB-392 这一跳 |
| B8 | 三者能分别水解苯乙酰胺 | 同文实验底物；历史 UniProt；KEGG R02540 | unsupported by audited sources；乙酰胺不是苯乙酰胺，泛 amidase 注释与反应数据库条目不能证明三目标底物特异性；不等于无此活性 |
| B9 | 三者为 W29 胞质酶 | 历史 UniProt comments/features/GO、论文方法、AlphaFold 记录 | unverified；模型使用 C_cy，不构成实验定位。当前来源没有足够定位证据 |
| B10 | 三者具有可独立补偿的该步 OR 关系，或应改为共同 AND | 同文删除/回补；数据库；模型布尔规则；结构记录 | unverified；来源未建立三者在同底物/同区室的独立替代，亦未建立共同必需复合体。不能以提升必需性匹配率作为 AND 依据 |
| S1 | 已有与三目标完全同序列的 AlphaFold 预测可复用，置信度有明确来源 | 三份 AlphaFold API JSON、model_v6.pdb、predicted_aligned_error_v6.json | supported；独立核对每条 PDB 全部 Cα 残基序列、API序列及目标一致。预测模型置信度不确认催化底物、定位或复合体关系 |

## 独立结构核对

均为 **既有 AlphaFold 预测**，API 标注 `AlphaFold Monomer v2.0 pipeline`、模型创建日期 2022-06-01、分发文件 version 6；本轮为 2026-09-11 获取，不是本次运行预测。未运行结构叠合、结合实验或新 BLAST，不能称已凭结构确认家族/底物功能。

| 目标 accession | 长度 | API globalMetricValue | PDB Cα pLDDT 均值 | pLDDT 最小值 | pLDDT<50 残基数 | PAE均值/最大值 Å |
|---|---:|---:|---:|---:|---:|---:|
| A0A1D8NI23 | 540 | 94.06 | 94.0427222222 | 35.62 | 2 | 4.0871193416 / 29 |
| A0A1D8NL97 | 539 | 96.00 | 95.9743413729 | 60.38 | 0 | 3.4816002974 / 28 |
| A0A1D8NNQ1 | 552 | 96.31 | 96.3370833333 | 70.81 | 0 | 3.3905823356 / 24 |

API 与 PDB 的均值略不同，保留口径，不混用。PAE 矩阵尺寸分别为 540²、539²、552²；均值是本次对下载矩阵的算术重算，不代表结合或催化正确性。结构均为单体预测，未提供本次靶向底物/复合体判断。

## 身份、重算与局限

- 本次独立重建 XML 计量列，2315 个反应和封存快照逐项相同；目标边界/GPR 从 XML 与快照分别核对，全部 WT 可行性使用封存 SD-Leu/PO1f 有效快照的完整边界。
- 主代理 `static_verification.json` 的 12 项输入 SHA 和 3 项历史 loader 源码 SHA 均重新校验相符。仅能说明这些文件身份匹配，不是恢复全部历史 dirty 环境。
- 当前模型 SHA256：`d274bad3050e3c9220a8b6287eae847f3bf1334892284d565a6c4d96b38135a0`。
- 旧反应快照 SHA256：`d9eb0a95af7cd6944b37f0667623c565f2a0999dc82b303a2d4d54774ca4bde3`；旧 results SHA256：`3e38ddf833e0f84c7bdeb21cccd3166c58f81b37ae09adab0f65c6291f163abb`。
- 原始 UniProt W29 缓存 SHA256：`e5b0a04874079b4057ffe25dadcb6b812cba8c96b227f41187187b27715753a3`。目标序列 SHA256 分别为 `a3c09c01edafe295fa2d8dd61b120b48fa05fed1fcae9b01f97ce69e1a325cd0`、`6108dfd127107f203f872532796a523f4fb6ca96c7b44f43a1a005476af5e092`、`bcb18dd33db365d088aa7d45b90e1cd273732b64da8ec63140cc134805fab1c3`。
- 本轮来源文件与 `sources/retrieval_curl.json`、`sources/retrieval_details.json` 的成功记录 SHA 均相符。初次 SSL 失败仍有记录，没有将失败当成无证据或无模型。
- `protein_evidence.json` 的全部来源 SHA 与目标序列、PDB 均值、低置信区和预测位点字段重新核对相符。三组标注位点均为 K/S/S，来源是 PIRSR 自动预测，不是 AlphaFold 新发现的实验活性位点；pLDDT 高不改变该证据等级。
- 论文审计直接打开 XML 正文、方法和图注，并独立目视本轮取得的 Fig.1/4 平板图片：中央 NS995 在乙酰胺上相对两侧明显受损，较高浓度及较久培养仍有可见生长；不估计未经量化的效应大小。未目检 PCR/回补原图、补充克隆文件或底层实验数据。对遗传学结论的支持限于作者报告、上述两图及其描述的对照；并非重现实验。Fig.1 SHA256 为 `7e24fb42a761861f126518171f8ff91925523f76a7f61a5eb7d1e4bc56456a83`，Fig.4 为 `a6a3535feee916fa83eb3c76e93aa4f07f860c39f808af189a4dce4c522a8f82`。
- 生长最优值继承旧 WT 的 optimal 报告；本次仅重算可行性。旧 WT 若为该固定问题的最优解，则收紧至这些 KO/反应关闭后同一见证仍可行，最优生长上下界夹定为原值。未核对新对偶证书，未声称本次重新求解。

## 主报告核对

已直接通读 `REPORT.md` 并对照本轮 `static_verification.json` 与 `protein_evidence.json`。YlAMD1 的 YB-392 遗传证据与 W29/CLIB122 参考序列桥接分开，乙酰胺与苯乙酰胺分开，残余生长未抹去；本模型完全行阻断证明与真实生物学分开；不把 AlphaFold 置信度当作催化、定位或 OR/AND 证据；化学计算只陈述保存元数据局限。未发现需要阻塞当前只读审查的来源冲突。

最终核对版本 SHA256：`REPORT.md` 为 `eb10a35f9228bbf57ea548a1666ad3b796344a877b793146b333156fd9290b9f`；`static_verification.json` 为 `b223fcb3f397eaa00746756bd740fa1a7683592e7ac9ea8512c8e3eb9258c5e8`；`protein_evidence.json` 为 `ca975db195a577d77d4218f9b3bbe017cb9a2e71fb48d817fc90f803e7904658`。最终蛋白证据中的所有来源 SHA 再次相符，新增的 `search_scope.json` 说明无关 locus Q6C5W4 的初始抓取不用于结论，与本审计实际使用范围一致。

对合并的数学说明亦核对成立：若三基因状态为 a、b、c，在暂不考虑当前供给断点时，两列所设边界的名义总容量是 `1000·a + 1000·(b OR c)`；共同单列、统一固定上界无法自动保持所有基因状态的名义容量。当前完整模型中两列实际可行通量仍为 0，名义 2000 不能称作本轮实现或计算出的通量。此属于 M2/M3 的条件性数学说明，不是独立生物学发现。

核对后保留四个科学未决项 B7–B10；不得将“报告审计完成”理解为接受候选 GPR。审计代码只执行读取、静态断言、数值重算，未重跑优化或预测。

# R794/R795 GPR 来源审查

2026-09-09；按 TASK.md 进行定点只读审查。输入完整 SHA、摘录 XML、历史表原始行和来源定位见 `provenance_sources.json`。本次未加载求解器、未运行 GEM 计算、未修改模型、GPR、病例账本或 Git。此文件是来源子任务的待独立审核交付，不能把本作者核对计为独立审计。

## 结论

该 12-AND OR 2-AND 规则已存在于出版商公开的 iYli21 `mmc6.xml`，本轮取回附件与本地 `data/iyli21.xml` 完整字节一致；它不是近期 essentiality 修模才新增。它与本地更早命名的 iYali v4.1.2 对应反应存在具体布尔结构差异：旧文件是13成员的共同 AND；按现有跨版本ID表转换后，恰好为现12成员加上 YALI1F13017g。iYli21 新增了 YALI1F38820g，并把原共同成员 YALI1F13017g 移到与新增成员组成的独立替代分支。这个差异是可核对的软件/来源事实；尚没有找到当时为何这样改的逐反应原始证据。

本来源链不能单独证明当前 OR 生物学错误，也不能证明旧13成员规则正确，更不能支持为了让基因必需而直接删分支。其价值是把待查问题缩小到“原共同成员为何退出共同依赖、新成员为何被赋予独立替代能力”。

| 关键位点 | 名称、功能及证据状态 | 来源角色 |
|---|---|---|
| YALI1F13017g（本地表映射 YALI0F09405g） | 本来源子任务未独立确认通用名称；蛋白功能以主任务序列/文献审查为准。此处证据等级仅 model/GPR assignment only，不从反应名反推蛋白功能 | 在本地 iYali v4.1.2 为共同13成员之一；在 iYli21 为后2成员分支之一 |
| YALI1F38820g（本地表映射 YALI0F31119g） | 本来源子任务未独立确认通用名称或蛋白功能；证据等级仅 model/GPR assignment only | 不在本地 iYali v4.1.2 两个泵反应内；在 iYli21 与上述成员共同组成后分支 |

此前工作消息中提出“末尾两个可能是 a 亚基同工型”仅为未检验假设，已撤回；主任务序列证据提示二者属于不同亚基类别。此审查不保留该猜想作为可用结论，不给精确新 GPR。

## 可核对的沿革

1. **图源暂定参考** P1 (`reference.xml`) 的 R794/R795 分别在100436、100470行。均为12-AND OR 2-AND；R794 `[0,1000]`、R795 `[0,0]`。这里的身份是图源暂定参考，文件复制或“canonical”称谓不代表正式模型批准。
2. **构建输入** P2 (`data/iyli21.xml`) 16713、16747行具有同一 GPR 和上述边界。P14 `scripts/update_model.py:28` 指向该输入。最早能从本仓 `git log --reverse --all -- data/iyli21.xml` 定位的提交是 `05ae6e16955a180cdb6416143312ad806ada0440`，日期2026-03-19T14:12:05-07:00；`git show` 取出的完整文件 SHA 与当前 P2 相同。故规则至少在该时点已存在于本地构建输入。该事实不等于完整历史执行环境已恢复。
3. **本地 iYali v4.1.2 对照** P3 31733、31782行的 y001085/y001086 均为13成员 AND，二者均允许正向通量。两反应自带 `Confidence Level: 2`、`NOTES: Check`，引用 PMID11278748 和 PMID11836511；这些都是模型内注释，不是对精确GPR的实验认证。按 P5 ID表转换，旧成员集与现14成员集只差新增的 YALI1F38820g；原13成员中的 YALI1F13017g 在现模型被移出共同组。YALI0E16192g 没有在该表猜补映射，只保留字面一致的旧ID。
4. **独立模型对照** P4 是本地 Yeast-GEM v9.0.2，物种 S. cerevisiae S288C。其 r_1085/r_1086 各有两个大量共享成员的 AND 分支；不存在本题这种“多数亚基一组 OR 完全不相交的两成员一组”结构。该文件可作为布尔结构对照，不能移植为 Y. lipolytica 的正确规则，且不代表本任务已核对其所有亚基或最新在线版本。
5. **近期 curation 来源排查的范围** P11/P12 两份 `gpr_isozyme_additions.csv` 字节相同；与 P13 `curated_model_patches.csv` 中，均没有本题反应或末尾两个位点条目。这仅排除本次点查的三份表，不能宣称搜索过所有历史修改来源。

2022年 iYli21 原始论文为 Guo 等，DOI [10.1016/j.csbj.2022.05.018](https://doi.org/10.1016/j.csbj.2022.05.018)。初次PMC网页、出版商网页和附件访问失败的记录保留；后续主任务通过系统curl成功取回[出版商 mmc6.xml](https://ars.els-cdn.com/content/image/1-s2.0-S200103702200174X-mmc6.xml)，本子任务重新读取、解析并计算哈希：2,614,400字节，SHA256 `6974b7588f2a6c60ba2cde2f26e20d3aba1334d0d501572bf01cee47eda86631`，与本地 P2 完整字节一致。R794/R795的GPR、各自边界也完全一致。

随后本子任务成功打开[Europe PMC提供的论文全文XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9136261/fullTextXML)：其DOI、题名、PII `S2001-0370(22)00174-X` 与论文一致，`supplementary-material id=m000511` 明确列出 Supplementary data 6，`href=mmc6.xml`。出版商下载URL含对应PII和mmc6.xml，因此附件对应性有原文元数据支持；该PMC XML给的是相对文件名，并非完整CDN URL。摘录、全文XML哈希及读取时间保存于 `paper_attachment_linkage`。已经不再把“原补充XML未取回”作为缺口；但论文内外尚未核实本题精确GPR的原始修订理由、依据及具体经手人。


## 历史计算已经提示的限制

P6 旧33组 isozyme ledger 中 EGC-4b3a988fa17f 绑定 `39f4cae11c3f270400c8a227c78b6af3ed412e85b1ade6cb604b0f85c3d8b1d9`，记录11个FN，单KO比值和全部关联反应关闭比值均为1，标为 `noncausal_redundancy_signal` / `reclassified_noncausal`。人工决定为 pending，pipeline patch 为 none，回归为 not_run。

P7/P8 后来的2026-08-06诊断绑定另一模型 `0f3a6c2b151e945b3461d3fa85f04575f8e8570ba817ed2879013aec91f62415`，将涉及泵的11个FN归为 inactive，保存的关闭关联反应生长比值约1。R794 标为 `growth_objective_exclusion`，R795 标为 `hard_bound_closed`。该案已不在 P10 的24组 isozyme表中。这里核实的是保存的交付物，不是本次重新计算。

P9 记录旧计算使用 SD-Leu（完整SHA在JSON）、Gurobi13.0.2 / COBRA0.30.0，Git提交35c959b3032b14661653a5bdd8eb2f10c11d5495、dirty=true。manifest没有独立 runtime strain 字段，不能回填PO1f或默认配置。因此这些旧数值不能直接声称在图源 bc2aac8f 和其目标培养条件已复现。它们足以限制“只修 GPR 就能让这些基因 essential”的承诺。

## 供独立审阅的原子声明

| Claim ID | 声明 | 直接来源及定位 | 证据类型 | 本作者核对 | 独立审计 |
|---|---|---|---|---|---|
| PROV1 | 图源两反应与本地原始 iYli21 输入的 GPR、各自边界相同 | P1/P2 两反应XML | 代码/数据直接静态核对 | supported | unchecked |
| PROV2 | 2026-03-19本仓最早可查 iYli21 文件与当前P2字节相同 | `earliest_local_history` | Git对象直接核对 | supported | unchecked |
| PROV3 | 本地iYali v4.1.2是共同13AND；按现有映射，F13017从共同组移出且新增F38820 | P3两反应、P5、`identifier_delta` | 数据库映射与结构比较；非序列身份认证 | supported（仅限给定映射） | unchecked |
| PROV4 | 本地Yeast-GEM两分支共享大量必需成员，不同于12 OR 2不相交结构 | P4两反应XML | 模型对照，非生物学证明 | supported | unchecked |
| PROV5 | 历史相关病例记录全反应关闭后生长比值约1，且未获得patch/回归认证 | P6案行、P7目标行、P8、P9 | 历史交付物核实；非本次复现 | supported（仅限所载历史输入） | unchecked |
| PROV6 | 取回的出版商mmc6与本地iYli21完整字节相同，原文元数据将mmc6列为该论文补充材料6 | P15、paper_attachment_linkage | 公开补充文件与原文元数据直接核对 | supported | unchecked |
| PROV7 | 已找到该GPR原始修订理由和生物证据 | 尚未定位逐反应理由 | 未核实 | unverified | unchecked |

本子任务审计覆盖：`7 total | 0 independently audited | 0 independently supported | 0 independently unresolved | 0 contradicted | 7 unchecked`。最终总报告须由另一个审阅者重新打开来源，再统计正式覆盖；不要将本作者的6条 supported 自动算作独立通过。

建议：来源结论可纳入审查，精确 GPR 修改保持提案待定。缺口是原始逐反应整理说明、目标蛋白的序列与亚基/区室证据，以及适用暂定参考下反应对生长依赖的验证。前两者限制归属结论，最后一者限制 essential 效果承诺；它们不阻塞本次来源审查收口。

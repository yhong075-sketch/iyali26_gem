# 检索边界、来源访问与执行记录

## 已做及未做

输入为既有coverage_1290.tsv；只将其中1282个 `identified_absent_from_model` 作为筛选总体。8个 `mapping_conflict` 按原始字段单列。读取全体既有蛋白注释，宽关键词筛选命中131条，再结合代谢、运输、辅因子成熟与复合体的可定位线索人工选20项。关键词不是分类器；未命中、未入选或旧 `non_metabolic` 均不作为无代谢功能证据。

官方UniProt共检查21条记录，其中20条进入最终列表；最初YALI1D30420g — 无已核实基因名 — 肽运输候选（自动／旧注释）未入最终20，因为底物范围及具体模型步骤较不明确，后由证据更直接的PMR1线索替入。该未入选记录只保留为原始注释线索，没有进行单案精查或重分类。20项之外的1262项未完成逐基因文献检索。前五之外的15项仅浅查，故1282项中共1277项未精查。

按最终排序精查5项，而非把最初选出的前五固定不变。官方PDB9个entity只做原生成员／身份桥接浅查，最终四个复合体I案进一步对照原文表格和方法。这里不存在1282个基因均已完成原生功能验证的主张。

## 检索来源与式样

执行日期2026-09-07。使用公开网页检索、PubMed元数据／摘要、PMC与期刊原文、官方UniProt JSON及RCSB实体JSON。只用原生Y. lipolytica实验支持原生功能；异物种同源或数据库只交叉核查身份／发现候选。未以搜索结果数量作置信度。未查询付费私人论文库或把无法读取的来源记为已读全文。

前段查询围绕已有系统ID／旧ID、NADH复合体成员、硫辛酸附着／合成、铁硫簇、细胞色素成熟与离子运输展开。此处是主题记录，不伪装成完整逐条检索日志。末段如下检索式实际运行，用于支持／反向证据及条件核查：

- `"Disruption of PMR1" "Yarrowia" "Materials and Methods"`
- `"Molecular cloning of YlPMR1" Ca2 EGTA strain`
- `"Yarrowia" "NUZM" deletion knockout`
- `"Yarrowia" "NIDM" deletion knockout`
- `"Yarrowia" "NIMM" deletion knockout`
- `"Yarrowia" "NUJM" deletion knockout`
- `"Disruption of PMR1" "Ca2+" "10 mM" SMS397A CS3`
- `"Yarrowia" "PMR1" transport assay ATPase activity purified`

UniProt查询原文及取回时间／SHA记录于[取回记录](/Users/david/.codex/worktrees/55a8/iyali26_gem/artifacts/missing_metabolic_function_followup_20260907/uniprot_retrieval.json)；另外的PMR1条目单独保存。PDB查询实体15、25、26、27、28、29、33、39、42。序列来源、entry/sequence版本、长度、SHA、源行及作者链ID均保存于[排序数据](/Users/david/.codex/worktrees/55a8/iyali26_gem/artifacts/missing_metabolic_function_followup_20260907/ranked20.json)。追加单条下载未保存独立的精确请求时间，登记日期和保存的内容SHA；不杜撰秒级时间。数据获取使用公开只读请求，没有提交私人输入到外部服务。

检索到的Candida同名亚基删除研究因物种不同，不作为本任务原生必要性证据。其他原生复合体I亚基的删除论文只能提供检索背景，不能代替最终四个目标自身的实验。商业重组蛋白目录和数据库中一般酶名不作为催化实验证明。多数据库对同一结构的映射不算独立重复。

## 实际读取范围与访问限制

| 来源 | 实际访问层级 | 本轮未完成 |
|---|---|---|
| Grba & Hirst 2020，DOI 10.1038/s41594-020-0473-x | 完整HTML可打开，读取相关Results、Methods、Extended Data Table 1及正文限制 | 未重解原始密度，未逐一检查全部补充文件／引用；结构样本与目标筛选菌株等同性未验证 |
| Angerer等2011，DOI 10.1042/BJ20110359 | 完整HTML可打开，读取方法、Table I、相关结果与Figure 5说明 | 未重分析原始MS；定位表的推测部分保留为推测 |
| Park等1998，PMID9461422 | PubMed元数据及完整摘要 | Gene原文、构建和重复／剂量细节未取得；期刊页面访问失败 |
| Sohn等1998，PMID9852022 | PubMed摘要；搜索索引返回的PMC Methods／Results／Discussion片段 | PMC直接打开反复出现浏览器验证；ASM访问受限，Europe PMC全文接口404；未完整阅读PDF |
| UniProt、RCSB | 官方原始JSON、明确字段与序列比较 | 没有BLAST、AlphaFold、靶菌株全序列重建，也没有把登记序列覆盖等同于密度覆盖 |

上述访问失败是公开网页／接口限制，不是自动审批拒绝。没有由此请求用户再批准；使用可取得材料并下调结论范围。关键来源的题名、DOI、定位与独立性登记于[source_registry.json](/Users/david/.codex/worktrees/55a8/iyali26_gem/artifacts/missing_metabolic_function_followup_20260907/source_registry.json)。本报告不是系统综述，不声称检索穷尽或已完成偏倚评估。

## 实际采用的技能和规则

- [govern-agentic-research](/Users/david/.codex/skills/govern-agentic-research/SKILL.md)：已读SKILL，采用中性问题、只读边界、条件匹配、主张溯源和人类闸门；总协调另行审计，本任务不自称完成独立审计。
- [gene-identity-function](/Users/david/.codex/skills/gene-identity-function/SKILL.md)：已读SKILL；每个候选分开记录系统ID、名称、蛋白功能及证据状态；分清identity、function和model assignment。
- [academic-research-suite](/Users/david/.codex/skills/academic-research-suite/SKILL.md)：已读SKILL；进一步阅读deep-research WORKFLOW的路由／事实核查相关部分及[来源验证角色](/Users/david/.codex/skills/academic-research-suite/ars/deep-research/agents/source_verification_agent.md)；采用原文定位、元数据与命题核对。没有运行全套ARS流水线、S2 API验证或模拟独立审稿。阅读ARS中的deep-research流程不表示启用了另一个Deep research插件。
- [ponytail](/Users/david/.codex/plugins/cache/ponytail/ponytail/4.9.0/skills/ponytail/SKILL.md)：已读SKILL，复用既有静态XML解析器，以标准库及已有openpyxl做表格读取／真实输入断言，没有新装依赖或引入框架。
- [Spreadsheets](/Users/david/.codex/plugins/cache/openai-primary-runtime/spreadsheets/26.905.11957/skills/spreadsheets/SKILL.md)：已读SKILL和已配置依赖位置；只读TSV／S2原始表、保留标签和行定位。交付Markdown／JSON，未制作新的XLSX，故不声称完成新工作簿公式重算／视觉QA。

PDF、gene-annotation-auditor、essentiality专用代理或patch builder均未调用。本轮没有新模型修改，四种中文essentiality命令的登记／接受流程不由此次筛选触发。已检查本任务祖先目录、本目录和输出路径规则，也读了源项目规则、PROJECT_STATE相关内容和BENCHMARK_CONTRACT／参考清单；历史任务的许可没有用来扩大本包范围。

## 完成检查及下一步边界

已执行静态输入与真实数据断言：1290=1282+8；20条唯一且属于既有缺失组；源行与旧ID一致；最终五个ID顺序一致；9个登记蛋白序列一致；R1889 GPR空；R2062/R573系数相等而R573为0…0。最终验证另外核对原始字段未变、五档案／24主张／11来源引用完整、输入SHA未改变，以及Git修改仅位于本次新增目录。检查记录见[validation.json](/Users/david/.codex/worktrees/55a8/iyali26_gem/artifacts/missing_metabolic_function_followup_20260907/validation.json)。

独立来源审计0/24，尚未进行；化学平衡／微物种审计、SD-Leu求解与删除回归、目标菌株序列验证以及直接原生运输／单亚基活性核查均未完成，也未据此提议激活补丁。静态输出可以供总协调审核，不替代后续验证。本轮求解器调用0、模型/GPR/旧标签/共享状态写入0，提交和推送0。

开始时间2026-09-07 22:56:10 UTC，最长主动研究窗口60分钟，停止线23:56:10 UTC。实际结束和总耗时写入[validation.json](/Users/david/.codex/worktrees/55a8/iyali26_gem/artifacts/missing_metabolic_function_followup_20260907/validation.json)；到交付即停止，不用剩余时间继续扩张候选集。

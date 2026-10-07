# R334 / R2081 模型层限定审查

核验日期：2026-09-11。对象为状态入口登记的当前已发布 `model_metadata_trna.xml`，不是根目录旧 `model.xml` 或未接受的 GPR 假设版。只读解析 SBML、构建规则和已有数值记录；本次优化求解 0、模型/数据变更 0。本文仅覆盖模型机制，原生蛋白身份与生物学 OR/AND 依据交给独立来源检索，不因模型结果代替生物学证据。

已核对空的全局 `/Users/david/.codex/AGENTS.md`、各实际祖先、根目录 AGENTS 及目标目录；目标 artifacts 和 scripts/data 路径无其他 override。历史 loader 源码目录的 `artifacts/r608_engineering_20260907/code/AGENTS.md` 已读，其特定 FN 批量命令/模块工作流不由本次限定审查触发。外部冻结源所在 `/private/tmp/worland-fixed-audit-20260910` 及祖先无附加 AGENTS。本次不执行构建、完整静态扫描、screen 或集群作业。

## 结论

当前模型把同一胞质谷氨酰胺合成步骤表示为两条精确重复反应，使用不一致的 GPR。目标基因单敲关闭了 R334，但 R2081 的另一个 OR 成员使同一步骤保持可用。已有 WT、基因单敲、R2081 单关及联合关闭四份记录，足以支持“当前 SD-Leu / PO1f 模型中的重复反应保留解释该次基因单敲非必需预测”。本次重新核对旧见证的守恒和边界，不称重新求解或验证原生基因非必需性。

## 基因与反应身份

| 系统 ID | 已核实名称/符号 | 简要蛋白功能与证据等级 | 全部当前模型关联 |
|---|---|---|---|
| YALI1F00821g | 本模型层未核实名；XML 的 COBRAProtein367 是占位名 | 谷氨酰胺合成酶候选（model/GPR assignment only；旧实验来源表的同源功能注释不是本次原生酶活验证） | R334、R2081 |
| YALI1D16151g | 本模型层未核实名；XML 的 COBRAProtein1001 是占位名 | 谷氨酰胺合成步骤候选成员（model/GPR assignment only；其独立催化能力尚未由本层证据证明） | R2081 |

XML 的身份指针分别为 Q6C3E0 / CAG77624 / EC 6.3.1.2，以及 A0A1D8NED6 / NCBI 2911272 / XP_502772.1。这里只确认这些字段存在，不把跨菌株映射、数据库身份与蛋白功能视为已实验验证。定位：`model_metadata_trna.xml:168807`、`:177323`。

| 反应 | 模型名称 | GPR | 边界 |
|---|---|---|---|
| R334 | glutamine synthetase | YALI1F00821g | [0, 1000] |
| R2081 | L-Glutamate:ammonia ligase (ADP-forming) | YALI1D16151g OR YALI1F00821g | [0, 1000] |

两列都严格等于：

`m38[C_cy] + m50[C_cy] + m141[C_cy] → m35[C_cy] + m130[C_cy] + m143[C_cy]`

即 XML 标签中的 ammonium_H3N + L-glutamate + ATP → phosphate + L-glutamine + ADP，系数全部 1，全部在胞质 C_cy，无其他水/质子物种。反应与 GPR 原始定位：`model_metadata_trna.xml:78409`、`:78492`、`:155059`、`:155140`。

本次另对上述六个物种作定点公式/电荷求和：m38 标签为 ammonium_H3N，公式 H3N、charge +1；其余五物种 charge 均 0。两反应的 C/H/N/O/P 原子残差全部 0，但按“产物减反应物”计的电荷残差均为 −1。这是两列共有的化学元数据不一致，不能称当前反应化学已完全验证；它与 OR/AND 所表达的蛋白依赖是不同问题。未凭名称替换氮物种、增加 H+ 或改动计量，后续若处理须审查关联物种的统一化学表示，不能作为本轮 GPR 判断的捷径。

令 F、D 分别表示上表两个基因可用，当前规则为 `R334=F`、`R2081=F∨D`：

| F | D | R334 | R2081 | 两列上界之和 |
|---:|---:|---|---|---:|
| 1 | 1 | 保留 | 保留 | 2000 |
| 0 | 1 | 关闭 | 保留 | 1000 |
| 1 | 0 | 保留 | 保留 | 2000 |
| 0 | 0 | 关闭 | 关闭 | 0 |

这是完整布尔逻辑的直接推导，不是新的 KO 生长计算。现有静态扫描 `gpr_results.json` 的 `gene_checks[467]`、`gene_checks[812]` 与 `duplicate_groups[3]` 给出相同结果；D 单敲没有改变任何关联反应的边界，F 单敲仅关闭 R334。旧实验比较把 F 标为 essential/FN；D 是 unlabelled，不能把它补成实验非必需或新负例。

## 已有数值证据核对

实际历史输入为同一发布模型、冻结 SD-Leu 培养和 PO1f 运行时配置，最大化 biomass_C。历史使用 COBRA 0.30.0、optlang 1.8.3、Gurobi 13.0.1；Threads=1、TimeLimit=30 s、FeasibilityTol=OptimalityTol=1e-7、Method=0、Presolve=0。历史诊断运行时间为 2026-09-11 17:05:11–17:05:12 UTC，定向对照为 17:08:05 UTC。原文件保留历史 dirty 状态和原始返回；本次未重建完整历史环境。

| 已保存条件 | 生长原始值 | R334 通量 | R2081 通量 | 原始状态 |
|---|---:|---:|---:|---|
| WT | 1.8718823069403 | 0 | 3.30531432815666 | optimal |
| F 单敲 | 1.8718823069402977 | 0 | 3.3053143281566704 | optimal |
| 仅关闭 R2081 | 1.8718823069402988 | 2.9562120216769134 | 0 | optimal |
| F 单敲并关闭 R2081 | 0 | 0 | 0 | optimal |

数据定位：`nonessential_diagnosis_20260911/results.json` 的 `/runs/WT`、`/runs/YALI1F00821g`，以及 `controls.json` 的 `/runs/YALI1F00821g__R2081__reaction_only`、`/runs/YALI1F00821g__R2081__double`。后二者是反应关闭实验，不是 D 单基因 KO 的记录。

本次直接读取每份完整通量向量及 `reaction_snapshot.json`，按各条件实际关闭反应设置 [0,0] 后，重新以 `math.fsum` 计算全部物种行残差与界限违反。四份见证最大质量残差依次为 8.058630802595945e-15、1.3542143209761826e-13、4.355828306552845e-14、0；边界违反均 0。`controls.json` 绑定的前序 results SHA 与当前文件一致，历史诊断所记录的三个实际导入源码文件 SHA 全部仍匹配。

`model_static_screen_20260911/input_verification.json` 还保存了历史 loader 重建全部 2315 反应快照匹配以及 22 份旧见证检查；本次读取其中涉及本问题的条目，不重复运行全模型审计。当前 XML 两目标反应与旧快照计量、边界和 GPR 相同。

解释限制：剩余 R2081 通量约 3.305，远低于 1000，上界减少没有限制这份见证；不能据此称 WT 与 KO 整个可行域相同。WT 已使用 R2081，不能叫作 KO 诱导上调。旧联合关闭记录的 optimum=0 支持其当时求解结果；本次守恒核对本身只验证可行性，不能单凭全零向量重新证明最优值为零。保存的最优状态、定向对照和既有输入核验共同支持当前模型机制，不能外推全部培养条件或生物学真实不可替代性。

## 构建来源

1. 原始 `data/iyali26.xml:8336` 与 `:32196` 已有这两条反应、相同规则和 [0,1000]。原始两列也精确重复，但都比当前多一个产物 `m10[C_cy]`（H+）；重复关系和 OR 不是 metadata 阶段新引入。
2. 当前主构建 `scripts/gem_annotate/main.py:122` 读取 starting model；`model_metadata_trna.build.json` 明确本次原始输入就是该文件。`main.py:459` 的 isozyme 添加表在封存 `data/reference_build/curation/gpr_isozyme_additions.csv` 中没有这两反应或这两基因的记录，本次没有发现另一个目标 GPR 覆写条目。
3. `data/metadata_reaction_selection.json:4943` 与 `:7265` 的 fields 都只有 stoichiometry，before/after GPR 完全相同；`scripts/gem_annotate/reaction_selection.py:62` 仅写入 fields 指定项目。构建记录对两者也只记录 stoichiometry applied。
4. 因此可核实的直接来源是原始模型继承。本文没有确认最早作者何时、为何建立该 OR；不能把当前构建代码中的一般 isozyme 说明当作这个历史 OR 的特异实验证据。

## 外部 AND 的可用范围与候选修改约束

本次直接打开三份已冻结 iYLI647 源文件，其 GLNS 均为胞质正向 [0,1000]，规则 `YALI0D13024g AND YALI0F00506g`（外部旧版本系统 ID；两者名称为空，功能在这些文件里只达到该谷氨酰胺合成反应的 model/GPR assignment only）。corr 第一版带括号、后二版不带，逻辑不变。GLNS 两成员共同必需未由三个相同模型规则独立验证；这也不是 corr→corr_3 中新增的 AND 修复。跨版本 ID 的完整生物学确认不在本模型层完成。

不能直接采纳 AND，理由是：

- 外部模型是另一项模型赋值，三版本沿袭不等于三个独立证据；未提供两蛋白共同催化/复合体依赖的直接实验依据。
- 已存外部命中 essential 的结果只说明该外部网络的预测，比较培养只映射了 28/35 项供给，保留不同生物量和化学表示，不能由命中倒推本地 AND 正确。
- 只把 R2081 改成 AND 会令 F KO 关闭两条，但 D KO 仍由 F-only R334 保留合成；这不能称为已实现完整两成员共同依赖。
- 只删 R334、保留 R2081 OR，F KO 仍可用 D 分支；去重复本身不解决该 GPR 判定。反过来，直接将两列合并为一列 [0,1000] 会将 WT 合计上界从 2000 降为 1000，不能声称全可行域等价。

后续候选应由生物学证据分流：若 D 能独立催化同区室同一步骤，OR 可成立，应审查 F 的实验上下文/特异需求；若只有 F 有该催化功能，应在两列中一致移除不支持的 D 成员；若存在两成员共同必需证据，应同时处理两列的依赖关系或经容量核对后去重复。以上均为待审批候选，没有实施。验收至少包括目标布尔规则、重复列处理后容量影响、固定条件 WT/两个目标单敲以及相关既有正确结果；不能只验收一个 FN 是否被翻转。

## 后台身份记录

当前 Git HEAD `ff36d87eb6c8f933dfc4413f43a2fcefcf0eea27`。旧构建记录 HEAD 为 `b9723026d08e28de42fc72c28104456e23638ed7` 且 dirty；当前部分构建源码已继续编辑，不能把当前 HEAD/源码当作历史完整执行环境。本次未切换、提交、推送或清理。

| 文件 | 本次读取 SHA-256 |
|---|---|
| model_metadata_trna.xml | d274bad3050e3c9220a8b6287eae847f3bf1334892284d565a6c4d96b38135a0 |
| data/iyali26.xml | 5c8c199e2c5b622e97daf2b3500f763f83519fb598702a11dd153052c6a99f9d |
| data/metadata_reaction_selection.json | d994792f374732897c839172917558ace854a347a0a47f4b760b17572530c134 |
| scripts/gem_annotate/main.py | f78e1e56ef088f861a9ee032c5db4c779d0b4e84a7a9465472923426c98fe67e |
| scripts/gem_annotate/reaction_selection.py | aeec4d6855179466b825f3d9f0d0e3552d46401798d66d1af4b5cd5ad1ee372c |
| model_metadata_trna.build.json | 802e82ae6c02e64536c75c66707d0bd1def762edb98b1967de3f38d9a3ace6e5 |
| nonessential_diagnosis_20260911/results.json | 3e38ddf833e0f84c7bdeb21cccd3166c58f81b37ae09adab0f65c6291f163abb |
| nonessential_diagnosis_20260911/controls.json | d97daac8d41a2dbed3ba11c4b49d6b8d3632bc717ac3a837796a7393b99873df |
| nonessential_diagnosis_20260911/reaction_snapshot.json | d9eb0a95af7cd6944b37f0667623c565f2a0999dc82b303a2d4d54774ca4bde3 |
| model_static_screen_20260911/input_verification.json | 142a807b6918c689e26195c082bd05feab57b316948c746a968ffcae84fc2b59 |
| model_static_screen_20260911/gpr_results.json | bfdb497c9a5ec257f675e47b69338411077a697ac08b69810661e1e91a232d17 |
| /private/tmp/worland-fixed-audit-20260910/iYLI647_corr.json | 96d1ee6bc5bbb38a447490b8c343aebd41846448218f96315ef80e6b00360bfe |
| /private/tmp/worland-fixed-audit-20260910/iYLI647_corr_2.json | 0bc227319eab7d18d96e035dc1298d3d6e60dfa7c89bd8acbba68a7398ce646c |
| /private/tmp/worland-fixed-audit-20260910/iYLI647_corr_3.json | 329be540c099409c2c7b76ee581a23f86eefdbbffe97e9b03b39b5e9c014b5d2 |

实际培养文件 SHA 为 `ed176d26a373f98cc413ed2e32a71f5f060a06e343f90f7db25cd32eff268e85`；菌株 profile SHA 为 `35307853a477d0b8540919acc6cd18d922e1e010ce98fb355316172a15048383`；运行时 simulation fingerprint 为 `77f5ce6f4f1f57a0125e015fda6c060493dc5fc3dc5da734029a6f936ea0517d`。它们是旧实际运行记录中的值，不用默认参数回填；本次没有重新实例化运行条件。完整路径、培养供给、菌株操作、源码清单与历史 dirty 文本仍保存于 results.json。

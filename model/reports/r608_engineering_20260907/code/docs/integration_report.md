# Final pipeline integration report

本次工作从本地 `main` 的 `36bb6f0735e4c6458bd53c0ceb01952b116b8be7`
创建 `integration/final-pipeline`。审计 inventory 与实施方案在任何代码整合之前
已交付。未更新远程 refs、未合并回 main、未删除 worktree/branch、未覆盖其他
工作区未提交内容。

## 1. Integration inventory / Git decisions

| Branch / worktree | 功能、主要代码及测试 | 依赖与重叠 | 决策 |
|---|---|---|---|
| main `36bb6f0` | canonical annotation/chemistry/build、B-group biomass、CoQ9 runtime、calibration、`src/iyali26_flow` sampling/analysis；原有 265 tests | 外部 research workspace；canonical model 与运行时副本 | 作为基底，保留完整功能 |
| `codex/canonical-model-refresh` / `codex/workspace-cleanup` `52160ef` | 确定性 rebuild、路径解析、legacy forwarding；project paths/SBML/compat tests | 两分支同一个 commit，均为 main 祖先 | 已包含，不重复 merge |
| `codex/canonical-b-group-model` `96026c5` | B-group translation、patch runner、run registry、research paths、flow | 56 个 patch-equivalent commits；剩余历史 hygiene 无新增功能 | 不合并旧历史 |
| `codex/coq9-dilution-runtime` `8f839b8` | runtime CoQ9 dilution/dFBA；CoQ9 tests | main 已有随机 alpha、区间对齐、解释 flux 等后续修复 | 已包含，保留 main 最新实现 |
| `lipid-unlump` / `codex/lipid-unlump-a014` `2aa6a0c` | strict-sn candidate、moiety ledger、CoA/R1521/fatty-acyl-CoA audit、VLCFA correction；六组模块 tests | 与 main patches、biomass、化学约定重叠；候选与审计 activation blocked | 手动语义整合，从最新后继 `994b09b` 选择成熟代码 |
| `codex/dfba-new-fn` `2446e2a` | baseline/candidate dFBA、HPCC sharding、WT/rescue diagnostics；new-FN tests | 共用实验标签 loader、strict-sn fingerprint；独立动态培养基合同 | 手动整合最终版本；复用 main loader |
| `codex/r989-gpr-main-worktree` / `origin/lipid-unlump` `994b09b` | 候选 R989、R1521、R39/R40 修订及证据限制 | 全模型 fingerprint；暂定 mapping 不能成为生产确认 | 仅候选内保留；不替换 canonical model |
| `codex/module-agent-workflow` `35c959b` 已提交历史 | 旧 canonical/history 与 CI indentation | 57 个 patch-equivalent commits，CI 已在 main | 不 merge 旧历史 |
| 同上 worktree 的未提交功能 | `module_workflow.py`、JSON/schema、四个 agent TOML、CI、8 tests | 研究工作依赖/角色/交接验证器，不执行 GEM simulation；jsonschema 已存在于 lock | 审核后复制独立功能，追加依赖与 AGENTS 段落 |
| detached `bae7` `694e6ac` | 早期 annotation/C16:1/chemistry 历史 | 45 个 patch-equivalent commits，无独有必需功能 | 不纳入旧基线 |
| `origin/patch-1` / `upstream/main` | 历史 CI 修复 | 已进入 main；剩余空行差异 | 不重复整合 |
| `upstream/gh-pages` | 发布的生成结果 | 非 pipeline 源码 | 不纳入 |

三个旧 `/private/tmp` worktree 目录已不存在，但登记和 branch 都保留。
当前原工作区有未跟踪研究数据、连接配置、图片和演示文件；a014 有监控脚本；
lipid worktree 有未提交 dFBA 草稿和其他修改。这些内容未被重置或顺带提交。
未提交 module 功能的逐文件 SHA 和原 Git 状态见
[`integration_source_snapshot.json`](integration_source_snapshot.json)。

选择性整合排除了旧 remix/unlump 实现、教学 prototype、巨型参考模型、
生成候选 XML、debug 和大量历史报告。只保留运行合同和测试实际需要的数据，
两份 handoff 历史报告明确作为回归 fixture。`data/iyli21.xml` 仅用于锁定的
VLCFA 证据验证，canonical 原始输入仍为 `data/iyali26.xml`。

## 2. Architecture and usage

统一入口：`python -m scripts.gem_annotate.pipeline`。

```text
data/iyali26.xml + explicitly configured research inputs
  → existing gem_annotate canonical build (optional --rebuild, offline)
      annotation / chemical identity and balance / duplicate cleanup
      → medium / biomass / gene annotation / CoQ9 chemistry
      → FVA / curated gap fill / final curated patches and chemical checks
      → canonical B-group tRNA biomass / deterministic SBML
  → canonical baseline load and structural/solver validation
  → optional strict-sn candidate on a copy, gated by exact source SHA
  → objective, chemistry and boundary regression validation
  → deterministic candidate export
  → optional baseline/candidate dFBA comparison
  → hash-bound comparison summary and pipeline manifest
```

使用已有 canonical，不需要为每次模拟重新注释：

```bash
python -m scripts.gem_annotate.pipeline \
  --lipid strict-sn --output-dir outputs/lipid-001
```

从原始输入运行完整构建：

```bash
python -m scripts.gem_annotate.pipeline \
  --rebuild --research-root "$IYALI26_RESEARCH_ROOT" \
  --lipid strict-sn --output-dir outputs/rebuilt-lipid-001
```

连接现有 dFBA 比较算法：

```bash
python -m scripts.gem_annotate.pipeline \
  --lipid strict-sn --output-dir outputs/dfba-001 \
  --experimental "$IYALI26_RESEARCH_ROOT/state/essentiality/repository/consensus_essential_genes.csv" \
  --dynamic-medium data/media/po1f_csm_leu_dfba.csv \
  --hours 24 --step-hours 0.1 --initial-biomass 0.05 --growth-cutoff 0.01
```

所有输出目录必须是新路径。`manifest.json` 记录源和输出 SHA、软件版本、
Git commit、校验结果和失败原因。dFBA 子目录额外保留源码/输入 hashes、
参数、模型指纹、结果表及科学比较门。非最优求解或比较门失败返回非零。
短时入口 smoke 只验证软件衔接，不代表 24 小时或全基因验证。

CoQ9 runtime 与 new-FN dFBA 的储备、培养基、端点和求解策略不同，不能互换：
前者使用 runtime CoQ9 reserve 和 dynamic doublings，后者比较两个显式 SBML
的 biomass gain。两者复用基础代码并保留独立科学合同。原有
`iyali26-flow phase1/analyze/train-fmpe` 与 calibration API 原样保留；
随机采样继续使用显式 seed，capacity profile 继续校验模型指纹。
模块工作流验证器是可选管理工具，不是生物学运行步骤。

Public API 兼容：旧 `scripts.gem_annotate`、`scripts/update_model.py`、
candidate builder/writer、dFBA `compare_models` 入口保留。新增
`canonical_copy=True` / `--canonical-copy` 用于在新路径重建同一 canonical
biomass 表示；拒绝已有文件及实验 overlay 混入。新增 `--offline` 完整控制
注释联网。new-FN `main(argv=None)` 支持统一入口直接调用。

## 3. Dependencies and invariants

| 阶段 | 输入/输出和顺序要求 | 不变量与验证 |
|---|---|---|
| canonical build | raw model + research contracts → baseline SBML；化学身份先于平衡，medium 先于 FVA，最终 tRNA 层在 curated structure 后 | 原有 chemistry、GPR、tRNA、SBML tests；raw SHA 不变；离线缓存缺失不得联网；缺源抛错 |
| strict-sn candidate | 只能接受冻结 canonical SHA；在副本上扩展脂质 | source-copy tests；1134 新反应/状态；全部新增反应守恒；模板 bounds/GPR 与 biomass 权重保留；重复 candidate 输入拒绝 |
| validation | baseline/candidate → audit records；先校验再导出/模拟 | IDs 唯一，有限且有序 bounds，objective 及方向保持，原内部反应不得变成物质来源；不把未知化学式算作 balanced；已平衡反应不得退化 |
| runtime simulation | 模型副本/context + 显式培养基/strain/settings | 求解状态先于 flux；context 退出恢复；WT 必须有可比增长；无 silent infeasible success |
| CoQ9 analysis | completed WT/KO → cutoff calls | 未完成、非有限、WT 无增长 → undetermined；汇总拒绝错误计为 FN；schema 1.9 记录 call policy |
| sampling/calibration | SHA-pinned profile + data + seed | 原有 flow/calibration tests 验证 seed、resume、输入 SHA 和门控；不把旧 profile 套到新模型 |
| export/report | 新输出目录 → atomic deterministic SBML / JSON/TSV | writer 不重复实现；源码/输入 hashes；CLI/HPCC 不覆盖已存在 run；硬链接/符号链接 alias 检查 |

canonical builder 在其独立新加载模型上逐步修改。strict-sn 调整脂质与 biomass
并应用候选 CoQ/GPR 修订，因此只能在 canonical 构建完成之后运行，不能把
canonical 的后续 generic annotation 再应用到 candidate。CoA/VLCFA/R1521
审计独立运行，不偷偷修改 canonical 或解除证据门。

## 4. Conflicts resolved

- **旧核心覆盖风险**：旧 lipid `main.py` 会丢失 main 的 microspecies、CoQ9、
  tRNA、路径和导出语义。保留 main orchestration，仅接入候选步骤。
- **CoA 与 shared patches**：仅加入完整 inactive curation unit，未引入旧
  `apply_all_patches`、旧 CoA charge patch 或 remix。默认 blocked 语义保留。
- **共享实验 loader**：保留 main 更完整的正例/证据/strain API；补空 ID 与
  矛盾标签诊断，不降低原重复 ID 拒绝规则。
- **确定性导出重复实现**：candidate public wrapper 转调既有 atomic SBML writer。
- **隐藏状态与副本**：ledger 原先改动全局 solver，现仅设置新加载的模型；
  strict-sn 显式验证 `max` 目标方向，并为候选分配独立的区室名称字典，避免
  COBRA 部分深拷贝造成来源元数据变更。
- **路径与 frozen evidence**：去掉执行默认值对个人 worktree/HPCC home 的依赖；
  历史 evidence 中的原路径仍是 provenance 文本，内容 hash 仍验证。
- **VLCFA stale tests**：旧测试把特定 builder 字符串当作必需，和 main 已升级
  writer/顺序矛盾。改成可观察的 optional API、source-copy 和独立性检查；
  未修改化学、生长或 GPR expected values。
- **失败分类**：修复已存在的 CoQ9 截断结果分类错误，以及 new-FN diagnostics
  输出目录缺失/覆盖问题。完整模型/模拟测试验证真实行为。
- **重建接口**：普通新输出过去可能省略 canonical tRNA 表示；新增显式
  canonical-copy 模式，保留旧默认 API，并用 raw rebuild SHA 检查验证。
- **离线与资源**：补 KEGG 和 isozyme 遗漏开关；FVA 使用单进程，避免 CLI
  rebuild 派生整个构建进程并复制大型注释数据。删除无用 DEBUG 计数输出。

## 5. HPCC invocation

在干净 checkout 中先生成候选模型，再提交明确绑定输入的任务：

```bash
export REPO="$PWD"
export EXPECTED_GIT_COMMIT="$(git rev-parse HEAD)"
export CANDIDATE_MODEL=/absolute/path/to/candidate.xml
export EXPECTED_CANDIDATE_SHA256="$(sha256sum "$CANDIDATE_MODEL" | awk '{print $1}')"
export RUNS_ROOT=/absolute/path/to/new-fn-runs
export IYALI26_RESEARCH_ROOT=/absolute/path/to/research
sbatch scripts/hpcc_dfba_new_fn.slurm
```

源码、默认 baseline 和培养基的 SHA 从明确的 Git commit 解析，candidate 有
独立 hash。分片/聚合及 WT/rescue 模式保留原固定参数门；账户、partition 和
资源使用 sbatch 参数覆盖。脚本不会自行提交/重提任务。本次不提交真实集群作业。

## 6. Validation evidence

审计基线：main **265 passed，0 failed，319 warnings**；Python 3.13.5、
COBRApy 0.30.0、pytest 9.0.2，显式配置现有 research workspace。

逐功能 targeted：

| 范围 | 结果 |
|---|---|
| Module workflow | 8 passed；CLI 与 offline lock check 通过 |
| CoA/VLCFA/handoff + strict-sn/ledger + 两个输出/路径检查 | 120 passed，1 intentional skip，314 subtests；无失败 |
| new-FN + main experimental loader | 29 passed，1 expected infeasible warning |
| CoQ9 endpoint/status/summary | 22 passed，1 expected infeasible warning |
| offline + annotation/locus | 15 passed，无 warnings |
| unified entry + offline + compatibility/capacity | 17 passed；覆盖真实 candidate 两次构建、toy dFBA、源保护、chemistry/objective/infeasible |
| HPCC 输出保护 | 14 passed，2 subtests，1 expected infeasible warning；shell syntax 通过 |
| solver 隔离及原始 no-op regression | 2 passed，原测试 expected 保持不变 |
| 最后 candidate 边界及统一入口复核 | 4 passed，3 subtests；90 个已解释 warnings |

警告不是 blanket ignore：历史 GEM 的 polymer/wildcard formula 无法按简单分子式
解析；另外有明确触发 infeasible 的负例、COBRA Group warning 和 memote
SQLAlchemy deprecation。一个已有 skip 对应 frozen source 尚非 anionic CoA
约定，相关 audit 仍报告 blocked。

真实 raw rebuild → strict-sn → validation → export 已完成：

| 项目 | baseline | candidate |
|---|---|---|
| Reactions | 2313 | 3413 |
| Metabolites | 1877 | 3006 |
| Objective | biomass_C, max | biomass_C, max |
| Solver status | optimal | optimal |
| Objective value | 1.3875914754870777 | 1.386113657777711 |
| 可检测的内部不平衡反应 | 331 | 328 |
| 化学不可完整判定的反应 | 449 | 425 |

baseline SHA：`bc2aac8fecd8f2f5f20de7bb3c988bf46b3a5831e525f556498ed51159bc1bee`。
candidate SHA：`d5bc83aadf07bfbd3e56b81188b30c96ee9c0858b4434020174b769aefff29e1`。
生长值用于软件回归，不作实验准确性或生产接受声明。

两次独立 raw rebuild 的 baseline SBML 完全相同，两次 candidate SBML 也完全
相同。完整入口还完成了真实模型、1 个实验参考基因、0.1 小时的 dFBA pilot：
0 新增 FN，比较门通过，但 `production_gate_passed=false`。
紧凑验证记录见 [`integration_smoke_validation.json`](integration_smoke_validation.json)。

第一轮全套测试：415 passed、1 failed、1 skipped、793 warnings、314 subtests。
唯一失败为 `test_empty_curated_patch_table_does_not_mutate_model`。
分类为 **A：整合引入的全局 solver 泄漏**；表现为 **E：既有 GLPK 数值舍入**。
独立对照证明 main 和 integration 的补丁函数源码、模型内容相同；空补丁前
GLPK 连续求解已有约 `4.663e-15` 差异，Gurobi 没有该差异，确定性模型字节
保持相同。修复 ledger 全局状态后，原测试与新增隔离测试通过，未改 expected。
证据见 [`integration_solver_isolation.json`](integration_solver_isolation.json)。

修复后的完整 `python -m pytest -q`：**421 passed，0 failed，1 skipped，
837 warnings，316 subtests passed，464.02 秒**。最后的 objective-direction /
compartment-copy 边界修复另经 4 项 targeted tests、3 个 subtests 验证。

| 失败类别 | 本次处理 |
|---|---|
| A integration regression | 第一轮全局 solver 泄漏已修复；最终无失败 |
| B existing failure | 配置环境下 main 基线无失败；main 已有的 CoQ9 截断分类与诊断输出缺陷一并修复 |
| C external data/dependency | 本次所需输入和 solver 可用，无缺失失败；research workspace 仍是显式运行依赖 |
| D stale test/benchmark | 两项 VLCFA 旧实现字符串检查已根据新架构改为行为检查；没有调整生物学/数值 expected |
| E environment/version | 独立复现 GLPK 最后几位的既有数值差异；隔离 solver 后原 no-op 测试通过 |

额外检查：新入口与既有 flow CLI import/help 通过；无 circular import；
HPCC shell syntax 通过；新增运行代码没有个人 worktree 或临时目录默认值；
无 tracked cache/debug/generated candidate；`git diff --check` 通过。
原始工作区 Git 状态保持一致，module 工作区 14 个来源文件 SHA 保持一致，
main 仍为 `36bb6f0`。`src/iyali26_flow`、legacy forwarding、`model.xml` 和
canonical raw `data/iyali26.xml` 均没有内容差异。

## 7. Remaining scientific boundaries

- strict-sn、CoA biochemical-pH migration、connected-component chemistry、
  cardiolipin expansion，以及候选暂定 GPR 保留其现有 blockers。
  软件整合不等于批准生物学模型发布。
- 当前模型仍有历史不平衡与未定义/聚合分子式；报告明确区分已知不平衡和
  无法判定。没有为了让测试变绿重新定义这些物种。
- 全 24 小时全基因 dFBA 和真实 HPCC 集群资源验证不属于短时入口 smoke；
  需使用明确参数、完整实验输入和独立输出运行，结果仍经过科学门。
- CoQ9 output schema 从 1.8 升为 1.9；严格 schema 消费者需要重新运行分析，
  不会把历史截断结果自动认证为完整比较。

## 8. Final recommendation

**SAFE TO MERGE INTO MAIN** — 结论针对本报告所述的软件整合和保留门控的
研究工具。没有未修复的 integration regression。没有自动将任何 blocked 或
provisional 生物学候选激活到 canonical 模型，也没有执行 main 合并或远程发布。

Integration branch：`integration/final-pipeline`。
持久 worktree：项目旁的 `iyali26_gem_integration`。
保留所有原 feature branches/worktrees 和未提交研究工作。

功能提交：

```text
1fafcdb Integrate validated metabolic module workflow controls
69befe3 Prevent essentiality calls from incomplete CoQ9 simulations
a2a998a Integrate guarded lipid candidates and chemistry audit workflows
8ea0cb1 Integrate portable dFBA comparison and failure diagnostics
dd7a4ec Connect canonical rebuild, lipid validation and dFBA pipeline
574b73b Keep lipid ledger solver selection local to its model
c240edf Preserve prior HPCC run status and serialize paths safely
770118c Protect lipid candidate objective and source compartment metadata
```

本报告和紧凑证据随最后一个 documentation commit 保存。

# 候选构建与覆盖修正

本轮已实际修改候选构建机制，原始基线 XML 未改。新增 `scripts/gem_annotate/energy_candidates.py`、`data/energy_candidate_repairs.json`、独立入口 `scripts/build_energy_candidates.py` 和构建测试；在原 `main.py`/`cli.py` 增加显式默认 E0 的候选开关，并在 `reaction_selection.py` 添加持久化优先级保护。未运行完整历史构建；执行的是固定完成态参考的最终候选阶段和当前可执行函数的窄范围阶段重放。

## 已实际重放的覆盖点

输入使用当前基线的已注释物种与基因目录，删除四个目标列后，执行真实 `add_gap_fill_reactions`（仅原表四行及对应方向整理）、`fix_proton_water_balance`、`apply_metadata_reaction_selection`、E3 候选补丁、再次 metadata、SBML 导出/重载。35 份反应阶段快照包含完整计量、边界、GPR、字段来源/顺序及生效 notes；见 `build_stage_snapshots.json` 与 `build_stage_reaction_diff.tsv`。

| 反应 | charge-aware balance 后 | metadata selection 后首次错误覆盖 |
|---|---|---|
| R_PGAM1_PhosHydro | PEP 水解、[0,1000] | 计量反向成 PEP 合成，边界恢复 [-1000,1000] |
| R_NTP3pp | 无额外 H⁺ 的 GTP 水解、[0,1000] | 反向成 GTP 合成并加入反应物 H⁺，边界恢复双向 |
| R_NTP7 | 无额外 H⁺ 的 UTP 水解、[0,1000] | 加入产物 H⁺，边界恢复双向 |
| r0242 | DHA+Pi→DHAP+H₂O+2H⁺、双向 | 删除产物 2H⁺ |

R2010 同步保存，阶段间未变。覆盖来自 `data/metadata_reaction_selection.json` 明确列出的 `bounds`/`stoichiometry` 旧 metadata 定义；不是普通名称字段本身造成。为保持默认基线，本轮不全局重写旧选择数据，而使**显式候选 > metadata 化学/边界**：补丁先逐条核对完整反应与物种签名，再统一写入；任何冲突在第一次 mutation 前失败。候选定义及其物种约定存为 SBML notes，后续 metadata 在任何编辑前验证所有 guard，并跳过这些已保护的化学/边界字段。普通未受保护反应仍沿用原选择逻辑。

独立入口和主构建最后阶段共用同一候选函数，非诊断求解末尾 override。导出必须新建路径；重载比较全模型反应、物种、目标及完整实际求解器定义，拒绝非稳态自定义约束、被改动的稳态行界和额外矩阵变化。重复应用用绝对目标定义，不累加 H⁺；E3→E4→E5 与直接 E5 等价，并禁止退回较弱候选。

## 实际版本与来源界限

E0/E1/E2/E3/E4/E5 分别保存在 `candidates/`，每份有独立构建记录。E0 保持**加载后的定义**，不声称原 XML 字节相同。SBML writer 将 20 个原无显式 fbc:charge 的 tRNA 私有残基写为 charge=0；这是 COBRA 读入缺失属性时的默认表示，原 XML 的 notes 也没有提供该电荷证据。加载后的共享物种定义相同，但原始未记载状态仍须保留，不能解释为新化学验收。原始基线完整 SHA 保持不变。

E1 限制三个已知合成方向；E2 去掉 NTP3pp/NTP7 的额外 H⁺并给 r0242 增加产物 2H⁺；E3 合并两组。E4 累计修订 NDP1（禁合成并去掉反应物 H⁺）和 R72（保留 IP7 水解方向）。E5 累计修订 CAT2p（保留**过氧化物酶式乙醇氧化**方向）与 OAADCm（保留草酰乙酸脱羧方向，并去掉中性物种约定下多余的产物 H⁺）。均未修改共享物种、GPR、培养或正常 ATP 维护。

CAT2p 实際为 ethanol[cy]+H₂O₂[pe]→acetaldehyde[cy]+2H₂O[pe]，没有 O₂。这与 catalase dismutation 是不同反应；官方 EC1.11.1.6 的 peroxidatic activity 支持活性类别，不能确认模型基因或混合区室抽象。最终 E5 notes 已据此写明。

历史执行配置分别封存于 `build_specs/E0_E3_executed.json`、`E4_executed.json` 和 `E5_final_executed.json`，均实际核对对应构建记录 SHA。E4 曾多写一个未经核实的 R72 KEGG R05202 URL，它不是方向判断依据；当前整理数据及最终 E5 已删除，使用独立审核打开的 EC3.6.1.52/EC2.7.4.24 来源。已执行 E4 XML 不回写，保留原始记录并以此说明订正。正式生物学接纳、原生功能/定位与整体能量验收由主报告分别裁定。

## 保留的软件失败与验证

1. 初次阶段重放的窄模型缺少 proton 约定要求的其他区室物种，balance 函数按设计拒绝；保留 `trace_build_stages_attempt1.py`、`stage_replay.log` 及第一批局部缓存。第二次保留全部当前物种，仅处理五条目标反应，阶段验证完成，无优化/网络调用。
2. 初次求解器导出核验错误使用变量的显示字符串（包含 `0`/`0.0` 等界格式）；改为稳定变量 ID。
3. 随后的严格矩阵核验发现 Gurobi `model.copy()` 经 LP 中间格式将 R1372 小系数舍入，36 个稳态行存在真实微小系数差。严格 guard 保留，不放宽 1e−7 求解容差。实际候选生成始终从新读 XML 开始，未用该复制路径；round-trip 测试改为与真实入口相同的新读取，并新增篡改实际矩阵的拒绝检查。前两次测试失败日志保留。
4. 初稿 E5 的 CAT2p rationale 错写为 dismutation，独立审核在优化前指出；旧文件及配置保存在 `E5_preflight_rationale_error.*`，未作为执行输入。修正说明后新建最终 E5，原已执行 E0–E4 未改。

最终测试命令和状态见 `build_commands.txt`、`build_tests_final2.log`。本分工未调用任何优化，不占主任务的 150 次优化预算；运行时 guard 同时禁止网络。科学结果仍以主任务对最终 XML 的实际重新求解为准。

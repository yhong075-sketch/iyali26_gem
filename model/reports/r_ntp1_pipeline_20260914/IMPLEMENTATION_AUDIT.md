# R_NTP1 管线实施独立审计

核查日期：2026-09-14。审计者：`/root/ntp1_integration_audit`。范围：直接比较本轮 `before/` 副本与实际整理、选择函数和测试；核对已有审查的本地/外部证据边界及真实构建入口顺序。本代理没有修改实现、整理或模型，也没有运行优化。

**实施审查通过：8/8 项支持，无阻塞项。最终完整构建输出尚待独立核验，本节不把软件测试通过视为生长、必需性或原生酶功能验证。**

| ID | 核查声明 | 独立直接核查与判定 |
|---|---|---|
| I1 | 修改限于授权的 R_NTP1 计量/方向保留例外 | 逐项比较旧、新 JSON，反应集合相同，唯一变化的反应条目是 R_NTP1；其 before、species、GPR 和旧水解来源 notes 均相同。supported。 |
| I2 | 复用实际构建路径，不会被最终 metadata 恢复为旧可逆式 | `main.py` 先执行已有 gap-fill 水解整理，后调用 `apply_metadata_reaction_selection`；R_NTP1 现在 `fields=[]`、`after=before`，通过精确前提检查后保留四物种水解式及 `[0,1000]`。本轮 main、patches、上游方向整理文件 SHA 与修改前一致。supported（静态核实）。 |
| I3 | 旧历史选择身份和计数没有悄悄混用 | 旧 JSON 完整快照 SHA 与 before.json 一致；amendment 指向该快照，显式记录旧计数与当前例外。独立重算当前 fields 得到 stoichiometry=201、bounds=6、gpr=9，与声明一致；其余排除条目保留。supported。 |
| I4 | 化学检查与已知反态拒绝覆盖目标行为 | 新测试明确断言 ATP+H₂O→ADP+Pi 的四个系数、`[0,1000]`、质量/电荷平衡及原 GPR；构造旧五物种可逆状态，要求 conflict 并保持原反应/notes。不是把任意旧 XML 当成已授权迁移。supported。 |
| I5 | 通用选择逻辑只增加所需证据写入 | Python 实现 diff 唯一新增行是无选择字段分支的可选 `preserve_notes` 更新，位于全部已有化学/GPR/物种前提检查之后；其余反应选择行为不变。新增键使用 `ntp1_*`，不覆盖旧 `gap_fill_*` 来源。supported。 |
| I6 | 回归测试覆盖受影响通用行为 | 直接读取测试源码、tests.log、tests_execution.json：7 项通过，包含幂等、SBML 往返、局部冲突拒绝、布尔等价、tRNA biomass 保留和新增 NTP1 行为。本代理核验已运行记录，未重复运行整个测试组。记录注明第 2 次执行修复了 notes XHTML 箭头表示，没有静默改变科学字段。supported。 |
| I7 | 新证据标签没有升级原生 GPR 结论 | preserve_notes 明示“native GPR activity unverified; growth and essentiality effects not tested”，并明确取代旧无耦联 ATP 合成/H⁺形式。原生酶归属与 ATP 底物活性没有由本次化学修正获得确认。supported。 |
| I8 | 当前输入与旧比较模型身份固定 | 独立核对旧比较模型 SHA 为 `b83bbc65a00b6dc049501f8c22e82414d6133fb2c0b1845d24f4106e19a4fab0`，与 before.json 相同。使用具体已固定构建作比较，而非依文件名“最新”推定科研基线。supported。 |

覆盖：`total claims 8 | audited 8 | supported 8 | unresolved 0 | contradicted 0 | unchecked 0`。分母仅为上述实施声明；最终完整导出验收与未求解的生物学结果不在此分母内。

本次直接核对版本：

- `data/metadata_reaction_selection.json`：`abfc4a94a03b4ec8da9e3e2fc3ba2e4f04806ef802b6e93c442274065ca0cb24`
- `scripts/gem_annotate/reaction_selection.py`：`1e51bf3a8fbf5afa8d4af3c399b728e0850a3cfbede33231faf51ecba1bf77ec`
- `tests/test_reaction_selection.py`：`8fd70f497fb02858f61bdb864ed57e608b7b6d1232f43c9275d2dcc2d9d20aab`

科学来源复用前轮 [LOCAL_AUDIT.md](../r_ntp1_direction_review_20260914/LOCAL_AUDIT.md)、[SOURCE_AUDIT.md](../r_ntp1_direction_review_20260914/SOURCE_AUDIT.md) 与报告；本代理本轮直接打开这些文件，未重新开展文献检索或把外部蛋白证据限制算作新的实验支持。

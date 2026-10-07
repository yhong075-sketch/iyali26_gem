# R_NTP1 本地来源独立审计

核查日期：2026-09-14。审计人：独立本地来源代理 `/root/ntp1_local_audit`。直接读取固定 SBML、整理记录、metadata 选择、旧保存向量及当前相关构建代码；未运行 LP、未重跑 screen、未更改模型或科学整理数据。根代理完成报告，本代理独立核对其第 1、3、4、5 节的本地事实。

**覆盖结果：10/10 项有直接证据或明确限定；无阻塞本地事实。** 其中第 9 项为未实施候选的代数推论，不是模拟观察。蛋白原生功能及方向热力学依据交由外部来源审计；不计入本地通过数。

| 编号 | 审核主张 | 直接核查与判定 |
|---|---|---|
| L1 | 截图对应模型身份 | `DIRECTION.md` 记载 `model_metadata_trna.xml`；实测 SHA-256 为 `d274bad3050e3c9220a8b6287eae847f3bf1334892284d565a6c4d96b38135a0`，与旧记录一致。通过。 |
| L2 | 化学计量、区室、边界 | 直接读 SBML 反应 `R_R_NTP1`：ADP + Pi + H⁺ ⇌ ATP + H₂O，五种物质均 `C_cy`，边界 `[-1000,1000]`；无另一侧膜质子或独立高能供体计量。通过。 |
| L3 | 元素和电荷平衡 | ATP、ADP、Pi、水在固定模型均是中性式；H⁺ 为 H/+1。产物减反应物为 H = −1、电荷 = −1；C/N/O/P 为零。当前式不是平衡式，不能套用另一套 ATP⁴⁻/ADP³⁻ 物种默认值。通过。 |
| L4 | 当前 GPR | 原始 FBC 为 `YALI1A09766g or YALI1E41893g`。通过。 |
| L5 | 两基因名称及功能等级 | `YALI1A09766g` 正式名本轮未核实、EC 7.1.2.2 质子转运 ATPase 相关候选；`YALI1E41893g` 正式名本轮未核实、EC 3.1.3.2 酸性磷酸酶候选。固定 SBML 只有 `COBRAProtein663`/`COBRAProtein94` 占位名。EC 来源为 gap-fill 表，对应的 native 催化能力、定位、单基因独立催化、同工酶关系仍未决；仅属 model/GPR assignment。报告的限定正确。 |
| L6 | 22 个旧保存观察 | 从底层 `results.json` 的 12 个 run 与 `controls.json` 的 10 个 run 独立提取 R_NTP1；22 值全有限、均大于旧阈值 1e−8，与方向报告全部值逐一相同；WT = 157.55576349404305。仅复读这一反应的旧通量，不重新证明完整向量可行或最优。通过。 |
| L7 | 维护相加不是精确零循环 | 直接 SBML 列加法得到 `S_R_NTP1 + S_xMAINTENANCE = −H⁺[C_cy]`，未删去质子。xMAINTENANCE 自身按中性式平衡；边界 `[7.8625,1000]`。不能由两列声称完整封闭能量循环。通过。 |
| L8 | 旧整理与当前版本来源 | active direction curation 要求 reverse 后水解、边界 `[0,1000]`；metadata selection 明确选择当前五物种可逆式并保留 superseded notes，证据状态原文是用户选择版本、不是新的化学或原生功能验证。当前 `main.py` 先调用 gap-fill（约 310 行），后调用 metadata selection（约 613 行）；selection 检查已知 before/after 后应用目标字段。不能直接称无来源软件回归。通过。 |
| L9 | 候选化学表示与方向需共同处理 | 保留当前中性物种，ATP + H₂O → ADP + Pi 为平衡的水解式。仅改边界不修正 H/charge 不平衡；只去掉额外 H⁺但仍可逆，合成列恰好是 xMAINTENANCE 的负列。此为代数候选推论，未实施、未求解，不能当作 GPR 活性已验证。通过。 |
| L10 | 报告声明边界与输入保护 | 报告第 1、3、4、5 节区分本地注释不平衡、EC 映射、历史观察、候选与未决生物学；没有承诺改后 essential。所有 12 个检查输入前后完整 SHA 一致。通过。 |

旧 22 个向量源位于 `artifacts/iyli647_screen_20260910/nonessential_diagnosis_20260911/results.json` 与 `controls.json`；二者关联 SHA 本次核对。旧培养为保存记录的 SD-Leu，启用 `po1f_sd_leu_accrispr_v1`；本次没有运行时配置或求解参数。历史软件、solver 参数、培养记录和 dirty 状态按旧记录保存在 `LOCAL_EVIDENCE.json`，未从当前默认值回填。

复核入口：从仓库运行 `python3 artifacts/r_ntp1_direction_review_20260914/local_audit.py`。本次已执行通过，解析后直接断言物种、残差、固定输入、22 个旧值及输入不变，并生成 [LOCAL_EVIDENCE.json](LOCAL_EVIDENCE.json)。该脚本无求解器导入、无模型写入，仅写自己的证据 JSON。完整输入与脚本 SHA 见该 JSON；本审计没有重建完整历史软件环境。

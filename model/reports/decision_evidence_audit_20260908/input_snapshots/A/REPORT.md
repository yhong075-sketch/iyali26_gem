# CoQ9 finite-reserve depletion — REPORT

## 首屏：六个决策问题

1. **要区分的机制或决策是什么？** 关闭 R385、切断 Q9/Q9H2 总池的原模型净新合成后，检验现有人工有限 reserve 是否真的随生长扣减并最终限制增长；不是再次检验仍能合成 Q9 的 KO。
2. **新做了什么，复用了什么？** 新增一条 `po1f_nonlimiting / no additional gene KO / R385 closed / alpha=1e-4 / pool=1 / dt=0.03125 h / 0–1 h` 轨迹，共 32 个推进区间；有效背景仍含 `po1f_sd_leu_accrispr_v1` overlay。复用旧五条件中“R385+source 同关时 exact-zero growth”和“固定 source 可恢复 t=0 growth”的结果；旧八轨迹、三步长和 R558 条件均未重跑。
3. **最决定性的原始证据是什么？** 保存的人工 Q9 reserve 从 `1.0000000000000002e-6 mmol/L` 降到 0；累计 raw source withdrawal 为 `1.0000000000012337e-6 mmol/L`，实际库存扣减为 `1.0000000000000002e-6 mmol/L`，首次保存为 0 在 `0.5 h`。该事件的 raw pre-cap inventory 为 `-1.2333196759e-18 mmol/L`，由库存下界得到 0；本轨迹没有正的 `<=1e-12 mmol/L` snap loss。biomass 从 `0.01` 增至 `0.0200000000000397 gDW/L`，较参数导出的编码上限 `0.020000000000000004 gDW/L` 超出 `3.97e-14 gDW/L`，在所报数值容差内；最大累计预算残差为 `3.97e-18 mmol/L`。
4. **哪个解释被支持、削弱或仍无法区分？** 从 `0.5 h` 起，solver status 为 optimal、所选 raw biomass/source flux 均精确为 0，保存的 reserve 同时为 0；没有记录 glucose exact-zero 或 infeasible。该结果表明 Q9 reserve 约束在这个运行时编码中足以强制零增长，但单条轨迹不能排除所有其他模型约束或未保存的有限营养库存，也不能外推真实细胞中的 Q9 pool、周转或长期表型。
5. **当前允许及禁止的声明是什么？** 可声明该运行时实现以 `v_reserve=alpha*mu`、`Q(t)+alpha[B(t)-B0]=Q0` 与 `B<=B0+Q0/alpha` 为编码目标，而且本条所选数值轨迹在下表残差/超量内与这些关系一致。不能把这些等式当实验事实、细胞死亡、真实 CoQ9 定量校准、旧 STOP/时间步收敛的修复，或正式模型/GPR 依据。
6. **唯一建议的下一动作是什么？** 将本轨迹作为“有限 reserve 分支按预算工作”的运行时实现验证证据送回总控，进入提案审查；本轮无需再计算或科学变更。

## 结果摘要

| 指标 | 结果 |
|---|---:|
| run termination | `completed` |
| backend solver statuses | 128 `optimal` |
| advanced intervals | 32 / 32 |
| backend optimize calls | 128 / 128 |
| wall seconds | 3.54264 / 600 |
| saved Q9 zero time | 0.5 h |
| first optimal raw-exact-zero growth | 0.5 h |
| final biomass | 0.0200000000000397 gDW/L |
| dynamic doublings | 1.000000000002864 |
| integrated raw source withdrawal | 1.0000000000012337e-6 mmol/L |
| actual inventory deduction | 1.0000000000000002e-6 mmol/L |
| max absolute `source-alpha*mu` | 3.59e-15 mmol gDW^-1 h^-1 |
| max absolute reconstructed Q9-row residual | 3.04e-15 mmol gDW^-1 h^-1 |
| max absolute reconstructed Q9H2-row residual | 3.00e-15 mmol gDW^-1 h^-1 |
| max absolute total-Q residual | 3.59e-15 mmol gDW^-1 h^-1 |
| max absolute cumulative budget residual | 3.97e-18 mmol/L |
| max biomass-cap overage | 3.97e-14 gDW/L |
| max selected-solution bound violation | 1.78e-10 mmol gDW^-1 h^-1 |
| max selected-solution mass-balance residual | 4.48e-12 mmol gDW^-1 h^-1 |

`raw_q9_end_before_min_and_zero_clamp_mmol_L` 与保存值分开保留；runner 先限制扣减不超过库存并实施非负下界，再把正的 `<=1e-12 mmol/L` 剩余库存置零。本次耗尽事件来自前者，没有正的 snap loss。exact-zero 是编码事件，不是死亡判定。

## 资源、异常和证据边界

- Gurobi/pFBA readback、输入 SHA、实际调用层级、fallback 次数和模型恢复记录在 `run_manifest.json`。
- 旧 R558 数值例外原样保留：YALI1D05440g — no established gene name found in the evidence inspected for this report — native protein function uncharacterized（uncharacterized；model role only: GPR-assigned to R558, whose reaction annotation is myo-inositol 1-phosphatase）。旧 condition 02 的 R558 raw flux=`-1.8631115666087675e-9`、bounds=`[0,1000]`；本轮不改阈值、不重算、不据此暂停本机制实验。本轮另观察到 R558 的 bound violation 为 `1.7813790807036378e-10`，低于现有 `FeasibilityTol=1e-9`，如实保留而不调整阈值。
- 本轮未施加额外单基因 KO；“WT”仅是 runner 的 no-KO 标签，不表示未修改的野生型，因为 effective context 启用了 `po1f_sd_leu_accrispr_v1` overlay。R385 是 reaction-level runtime closure，不构成蛋白功能断言。
- `po1f_nonlimiting` 的 uracil state 为 NaN、uptake flux 可读。动态轨迹还保存 glucose，但没有保存其他 12 个有限营养池的库存状态，因此“未见 glucose 耗尽”不能扩大为“已排除所有营养限制”。
- 结果保持 `runtime_only / sensitivity_only_not_calibrated`。pFBA 只给一条 parsimonious optimum，不是 FVA。

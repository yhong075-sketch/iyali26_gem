# CoQ9 supply mechanism：静态 t=0 结果与审计

## 范围

本实验只比较五个相互独立的初始状态优化问题；没有推进动态时间、没有重跑旧轨迹，也没有持久修改 canonical model.xml、GPR、bounds、标签或共享状态。所有反应关闭、t=0 medium cap 和 Q9 source/dilution 都只发生在可丢弃的内存模型副本中，并逐项记录。结论均为 `runtime_only / sensitivity_only_not_calibrated`。

## 结果

| 条件 | 状态 | growth (h⁻¹) | 分类 | R385 | reserve source | dilution |
|---|---|---:|---|---:|---:|---:|
| 01_finite_R305_reaction_KO | optimal | 0.481691968964 | positive | 4.81691968964e-05 | 0 | 4.81691968964e-05 |
| 02_nonlimiting_R2062_reaction_KO | optimal | 1.46480464495 | positive | 0.000146480464495 | 0 | 0.000146480464495 |
| 03_nonlimiting_YALI1A21711g_KO_plus_R1889 | optimal | 1.28732167922 | positive | 0.000128732167922 | 0 | 0.000128732167922 |
| 04_nonlimiting_WT_R385_closed_source_closed | optimal | 0 | exact_zero | 0 | 0 | 0 |
| 05_nonlimiting_WT_R385_closed_fixed_source | optimal | 1.46507605681 | positive | 0 | 0.000146507605677 | 0.000146507605681 |

## 预注册问题的回答

- R305 reaction KO 对旧 finite YALI1A14736g gene KO 的 t=0 growth：`isclose=True`。
- R2062 reaction KO 对旧 nonlimiting YALI1A21711g gene KO 的 t=0 growth：`isclose=True`。
- 在 YALI1A21711g KO 已关闭 R2062 后再关闭 R1889，静态最优 growth 变化为 `-0.177482965727 h⁻¹`；该 pFBA 解中其余非零 Q9→Q9H2 入口：`R740,R1977`。这只描述一条 parsimonious optimum，不是 FVA 或路径唯一性结论。
- R385 关闭且 source 关闭时 growth 分类为 `exact_zero`；打开固定 source 上界后为 `positive`。`static_initial_rescue=True`。这只检验人工 source 在 t=0 能否满足 imposed demand，不验证储备耗尽、周转或长期 rescue。
- 每个 optimal 条件均直接审计 Q9 与 Q9H2 两行，并检查 `v_R385 + v_source = v_dilution = alpha*mu`；详见 `q9_balance.tsv`。

## 已保存通量的机制解释

### R305 关闭后 ATP 从哪里来

Condition 01 中 R305、R304、R2206 和 R171 ATP synthase 的通量均为 0，因此该解没有使用所建模的经典 Q9→复合体 III→复合体 IV→ATP synthase 路径。主要 ATP 生成来自底物水平磷酸化：R642 phosphoglycerate kinase `17.7569176346`、R694 pyruvate kinase `17.4322272855`、R741 succinate-CoA ligase (ADP-forming) `14.6420807896 mmol/gDW/h`。其中约 `14.4920636538` 的线粒体 ATP 经 R815 ADP/ATP transporter 转到胞质。它们支持了 biomass ATP 消耗 `11.1222675634`、NGAM `7.8625` 以及糖代谢和核苷酸等耗能反应。

这解释的是保存的一条 pFBA 最优解，不证明这些旁路在细胞内具有同样容量。该条件仍有 oxygen uptake `11.7463144068`，但主要来自 R1883 long-chain alcohol oxidase `9.4051713989` 和 R688 pyridoxine oxidase `1.9971811701`；氧摄取本身不等于氧化磷酸化已开启。

### 关闭 R1889 后如何改道

Condition 02 中 `R1889 22.7234413278 + R740 3.08708670243 = R305 25.8105280302`。Condition 03 再关闭 R1889 后，所保存的解改为 `R1977 32.5732197723 + R740 5.23231154252 = R305 37.8055313148`：R1977 使用 FADH2 还原 Q9，R740 使用 succinate 还原 Q9。与此同时 R171 ATP synthase 从 `55.3772684642` 降到 `44.5337218700`，growth 从 `1.46480464495` 降到 `1.28732167922 h⁻¹`，下降 `12.1165%`。这些数值支持“该最优解发生电子入口重排并损失部分模型能力”，不支持“R1977/R740 是唯一旁路”。

## 基因证据边界

- YALI1A14736g — no established gene name — native protein function uncharacterized（model/GPR assignment only；model role R305）。
- YALI1A21711g — no established gene name — native protein function uncharacterized（model/GPR assignment only；model role R2062）。
- Reaction KO 不构成蛋白功能断言；零 pFBA flux 不代表 FVA 证明反应不可用或不重要。

## 数值例外

Condition 02 的求解器状态是 `optimal`，但导出的 R558（myo-inositol 1-phosphatase）净通量为 `-1.8631115666087675e-9 mmol/gDW/h`，而模型 bounds 为 `[0,1000]`。因此独立 bound 检查得到 `1.8631115666087675e-9` 的下界超量，略高于冻结的 `FeasibilityTol=1e-9`；这是 11,575 条通量中唯一超过该容差的 bound violation，不能静默截为 0。

`optimal` 是求解器对内部优化问题的终止状态；这里的 bound violation 是对 COBRA 导出的净反应通量与名义反应边界所做的独立复核。未保存 forward/reverse 原始 primal 与 Gurobi quality attributes，不能事后指定误差来源。该项与 Q9/Q9H2 无关：同一条件的两行残差为 `±8.88e-16`、总 Q 和 dilution coupling 残差为 0，旧/新 growth 也精确重现。因此它是必须披露的非阻塞数值例外；严格的“全部反应 bounds 通过”判定仍为 False，依赖 condition 02 的边际量化保持 provisional。

## 执行审计

- 条件完成：5/5；native backend calls：10/12；总 wall：10.955 s。
- 所有条件 fresh context 与恢复检查通过：True。
- 所有 optimal 条件 Q-row 与 growth coupling residual 通过 FeasibilityTol：True。
- 所有反应 bounds 严格通过 FeasibilityTol：False（仅 condition 02 / R558；见上）。
- 状态判定：内部静态机制结论 `CONDITIONAL PASS`；严格全-bound 数值验收 `FAIL`；不授权生物学发表结论或正式模型/GPR curation。

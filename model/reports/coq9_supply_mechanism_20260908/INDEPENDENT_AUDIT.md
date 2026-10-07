# Independent audit — CoQ9 static supply mechanism

本审计仅重读并复算已有 JSON/TSV 与固定输入；没有调用优化器，也没有修改模型或科学输入。

## 覆盖与结论

| 项目 | 覆盖 | 结果 |
|---|---:|---|
| 条件身份 | 5/5 | 与 manifest 一致 |
| reaction flux | 11,575/11,575 | 唯一且为有限数值 |
| solver readback | 5/5 | 一致 |
| backend calls | 10/10 | 每条件 2 次，全部 optimal，无异常 |
| input/driver SHA | 10/10 | 当前文件匹配 |
| bound-change records | 18/18 | 最终目标 bounds 与 flux 表一致 |
| gene-KO impacts | 2/2 | 各只关闭声明反应 |
| old t=0 comparisons | 2/2 | 通过 abs `1e-8` / rel `1e-6` |
| Q9/Q9H2 terms | 95/95 | 十个 metabolite rows 完整 |
| source/dilution coupling | 5/5 | 通过 `FeasibilityTol=1e-9` |

总体判定为 `CONDITIONAL PASS`：五条件的 Q9 机制问题均可解释，唯一严格数值例外是 condition 02 的 R558 bound overage；因此不能声称整包所有反应 bounds 均严格通过。

## 独立复算

- R305 reaction KO 与旧 YALI1A14736g gene-KO t=0 growth 一致。
- R2062 reaction KO 与旧 YALI1A21711g gene-KO t=0 growth 一致。
- 所有 Q9/Q9H2 单行残差最大 `2.22e-15`，总 Q/direct identity 最大残差 `3.60e-15`；5/5 满足 `v_R385 + v_source = v_dilution = alpha*mu` 至少达到 `1e-9`。
- Condition 01–04 的 source flux 严格为 0。Condition 05 source=`1.46507605677471e-4`、ub=`0.0032`，与 dilution 相差 `3.59e-15`。
- Condition 02→03 只额外关闭 R1889；growth 下降 `0.177482965727011 h^-1`（`12.1165%`），仍为正；该 pFBA 解使用 R1977=`32.5732` 与 R740=`5.23231`。
- R385/source 同关得到 `optimal + exact_zero`，不是 infeasible；开放固定 source 后得到 positive growth。这只支持静态 t=0 人工 source rescue。

## 唯一严格数值例外

Condition 02 的 R558 flux=`-1.86311156660877e-9`，lb=`0`，故 lower-bound overage=`1.86311156660877e-9 > FeasibilityTol 1e-9`。它是 11,575 条 flux 中唯一超过该容差的 bound violation。Q9 rows、coupling、旧 growth 重现均不受影响，所以本轮不追加求解、不改阈值、不把整组科学比较判为失败；原值和 `all_bounds_within_FeasibilityTol=false` 保留。

## 证据边界

- YALI1A14736g — no established gene name — native protein function uncharacterized/unverified at the requested evidence level（model/GPR assignment only；model role R305）。
- YALI1A21711g — no established gene name — native protein function uncharacterized/unverified at the requested evidence level（model/GPR assignment only；model role R2062）。
- Reaction KO 不创建蛋白功能断言；pFBA 是一条 parsimonious optimum，不是 FVA。
- 结论保持 `runtime_only / sensitivity_only_not_calibrated`。

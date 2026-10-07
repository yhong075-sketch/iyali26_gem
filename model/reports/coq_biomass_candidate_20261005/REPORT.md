# 已选定 α＝10⁻⁴ 的 CoQ9 生物量候选

2026-10-05，用户明确回复“using the 10^-4”，选定 **α＝0.0001 mmol/gDW，未校准的暂定假设**。

有效交付为 [E5_coq9_biomass_alpha_1e-4_validated.xml](E5_coq9_biomass_alpha_1e-4_validated.xml)。输入固定为原历史E5，不合并其他CoQ化学或R1889 GPR候选。共享构建入口在 `biomass_C` 增加 `m468[C_mi]: −0.0001`，原有生物量系数及R385边界[0,1000]保持：

\[
v_{R385}=10^{-4}v_{biomass_C}.
\]

新增项表示随新生物量保留的CoQ9，并非每次电子传递消耗CoQ。数值选择不构成原生含量测量或其他CoQ机制的验收。

本轮通过实际XML导出／重新读取验证：完整模型和求解器定义一致；非目标计量、边界、GPR及原有注释保持；新增CoQ来源说明经XML实体解码后语义完全一致。构建记录见[build manifest](E5_coq9_biomass_alpha_1e-4_validated.build.json)，含输入／输出完整SHA、代码身份、选择来源和既有dirty状态。保护器记录本轮优化0次、网络0次。此前同α的[静态结果](../coq_biomass_growth_20261005/REPORT.md)和[dFBA曲线](../coq_biomass_dfba_extended_20261005/REPORT.md)为既有结果，本轮未重跑。

有限独立核对3/3项通过：实际输入／输出／实现SHA；唯一模型系数变更及原注释保持；相同静态培养／菌株运行时条件下，与此前α＝10⁻⁴保存优化问题的完整签名（含稀疏约束矩阵）完全一致。审核未新增求解；不是本轮重新验证生长或必需性。

首次文件 `E5_coq9_biomass_alpha_1e-4.xml` 因新增JSON注释的引号被XML实体编码、旧检查按字符串比较而未通过；模型方程没有失败。已保留该文件及[失败标记](E5_coq9_biomass_alpha_1e-4.failed.json)。仅对新增JSON说明采用实体解码后数据比较，未放宽其他注释、模型或求解器检查；重新构建得到上述有效交付。未覆盖原E5、改变默认构建、提交或推送。

实际成功命令：

```sh
.venv/bin/python -B scripts/build_coq_biomass_candidate.py --output artifacts/coq_biomass_candidate_20261005/E5_coq9_biomass_alpha_1e-4_validated.xml --alpha 0.0001 --alpha-source 'User selected 1e-4 mmol/gDW on 2026-10-05; uncalibrated provisional assumption, not a measured physiological value.'
```

后续交付说明：用户已授权提交／推送此候选。新增8次优化的结果见[生长依赖、R305化学和封闭ATP检查](../coq_candidate_checks_20261005/REPORT.md)：CoQ生长耦合及所测封闭ATP检查通过，R305的H和电荷残差仍均为−2。本次发布有效XML、构建入口及该次检查记录；上文静态扫描／dFBA链接为本地历史研究记录，失败XML仅在本地保留。重新构建时须指定新的输出路径，不能覆盖已发布XML。

发布前使用已提交的共享代码及本次暂存的构建入口／数据／测试，在临时目录完成合成软件检查和无求解重建，生成XML与交付文件字节完全一致；不依赖本地其他未提交修改，无新增GEM优化。记录见[publish_check.json](publish_check.json)。

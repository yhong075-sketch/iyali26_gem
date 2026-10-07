# V-ATPase 三共同 AND 假设版：实现独立审计

审计时间：2026-09-11T01:04:40.014294+00:00。审核者：独立子代理 `/root/vatpase_gpr_entry`。范围为用户明确选择的候选实现，源码只读；本代理仅写本报告。未执行构建、优化或新生物学检索，不改写既往证据审计。

**结论：当前实现静态通过。** 代码将三共同依赖显式作为待验证假设接入真实完整构建，默认不启用；没有把假设录入 accepted essentiality 病例或提高为实验功能证据。生长效应仍须根任务的实际运行结果判定。

## 被审计的身份和规则

| 系统 ID | 功能与证据边界 | 本次模型角色 |
|---|---|---|
| YALI1D00581g；原生正式名称未核实 | V1 D 中央转轴亚基候选，已有序列/AlphaFold 预测支持 | 共同 AND |
| YALI0E16192g；原生正式名称未核实 | CLIB122 V1 F 中央转轴亚基候选；W29 活性位点/伪基因对应冲突未决 | 保留这个旧 ID 作为共同 AND 假设，不自动映射新 ID |
| YALI1F38820g；原生正式名称未核实 | 偏 Vph1-like 的 V0 a 亚基候选，仍需实验确认；区室和替代性未决 | 共同 AND 假设，保留候选功能标记 |

设 G0 为原 R794/R795 两支路规则，数据中的新规则精确为 `YALI1D00581g AND YALI0E16192g AND YALI1F38820g AND (G0)`。原内部 OR 完整保留；不是全部基因直接 AND，也没有将后支路另一成员改成独立 OR。

## 静态核查清单

本表分母只含实现主张：**10 项 | 已审 10 | 静态支持 10 | 未决 0 | 反证 0 | 未审 0**。这不代表原生功能或 essentiality 得到验证。

| 项 | 源码/数据核查 | 判定 |
|---|---|---|
| 1 | JSON 只列 R794/R795 和三个准确系统 ID；status 为 provisional_hypothesis；来源文件均存在 | 支持 |
| 2 | helper 用既有 boolean_key 验证 after 等于 required AND before，并禁止改变基因集合或使用模型不存在基因 | 支持 |
| 3 | 所有靶标的 GPR、界限、计量及物种 formula/charge/compartment 先完整核对，再开始 mutation | 支持 |
| 4 | 末版增加全旧/全新判定；两条反应混合前后态拒绝，避免自动补齐部分迁移 | 支持 |
| 5 | 自有假设备注若已存在但不一致则拒绝；其他 notes 保留；已达目标规则和 notes 时幂等 | 支持 |
| 6 | CLI store_true 默认关闭，只有显式 flag 调用；要求 offline/no-solve/metadata 且禁止混合其他实验 overlay | 支持 |
| 7 | 强制 no-solve 后的已有输出/符号链接/原始输入/model.xml 路径检查，要求新输出 | 支持 |
| 8 | ordinary 排除候选，但 canonical_build 显式纳入候选，因此仍通过全部 20 个 tRNA 耦联步骤及同样 biomass notes；这是完整构建语义，不是发布批准 | 支持 |
| 9 | 先完成 metadata reaction selection 与 CoQ9，二者全部完成才应用假设，之后直接导出；不会被旧选择规则覆盖 | 支持 |
| 10 | 构建清单保存 enabled/evidence_status，JSON 自动纳入现有 data/reference_build 完整 SHA；未修改 accepted essentiality 门控 | 支持 |

实际前版 `model_metadata_trna_vph1like.xml` 的 R794/R795 仅有四项 metadata 选择 notes，无新 helper 自有 notes 键冲突。其 R795 仍为 [0,0]，候选数据保留该界限。

## 检查证据与限度

本代理读取了测试源码和根任务已生成日志：`tests.log` 记录初版 16 项相关测试通过；末版 guard 后 `tests_after_guard.log` 记录 3 项候选测试通过。候选测试覆盖两条规则的 2^14 布尔状态、三单敲关停规则、SBML 往返/幂等、第二目标冲突不部分修改、混合状态拒绝、默认关闭和新输出检查。本代理未重跑这些测试，不能称独立执行复现。

初次源码检查发现逐反应允许旧/新会容忍混合状态、自有 notes 可覆盖同名值；根任务已增加相应拒绝 guard，本报告核对的是修正后的快照。实际完整构建后仍需比较整模型差异，确认除两条 GPR 和假设备注外所有字段一致，尤其 tRNA、生物量、R795 界限和既有功能候选注释。screen 运行及生长结论由根任务实际产物另记，本报告不提前给出通过判断。

## 审计快照

当前代码 HEAD：`ff36d87eb6c8f933dfc4413f43a2fcefcf0eea27`；工作区存在未提交修改，以下 SHA 标识实际审计内容，不声称完整历史环境等同于 HEAD。

- `data/reference_build/curation/vatpase_gpr_hypothesis.json`: `b41e5cc151c575dd9f31c4f454b55f60c1c053b9e0f69915d9844f3ece200e6b`
- `scripts/gem_annotate/patches.py`: `cd19ecc729809c9661eec423b741a61313c3470cdae4017383750383c5566b5d`
- `scripts/gem_annotate/main.py`: `f78e1e56ef088f861a9ee032c5db4c779d0b4e84a7a9465472923426c98fe67e`
- `scripts/gem_annotate/cli.py`: `79fff4611d5cb2f10bd74ff450e0180eabee6bc30823d9625374611d694088b2`
- `scripts/gem_annotate/reaction_selection.py`: `aeec4d6855179466b825f3d9f0d0e3552d46401798d66d1af4b5cd5ad1ee372c`
- `scripts/gem_annotate/coq9.py`: `b587b3867a488fd992a9a9411a93cac8eec954bb1e59dbc35f9f0cea7e78cab5`
- `tests/test_vatpase_gpr_hypothesis.py`: `795b53dddc31264dfec34b5c38d5be3400b644e4eac921a9db8ef4c89f122ce5`
- `model_metadata_trna_vph1like.xml`: `f77f07f96c8f19d0507112813110dc2afe613ddb31c912c0700927f036b23736`
- `artifacts/vatpase_gpr_hypothesis_20260911/tests.log`: `617b7dbccd99ca61f439f75357bf54b92108a25a7cfc1f6d2d7fbadca734b4bf`
- `artifacts/vatpase_gpr_hypothesis_20260911/tests_after_guard.log`: `bc7095f32bfd3020d28d1a97d88391d6cdbe5683a66bb7251be10883d121e950`

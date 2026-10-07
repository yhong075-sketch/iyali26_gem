# R_NTP1 最终构建独立核验

核验日期：2026-09-14。审计者：`/root/ntp1_integration_audit`。直接用 Python 标准库 ElementTree 解析旧比较模型与新导出 XML，独立计算计量、元素/电荷残差、通量边界、反应及全 XML 差异，并核对构建清单与输入哈希；没有导入求解器、重新构建或执行 LP/screen。

**最终输出核验通过：8/8 项支持，无阻塞项。新模型只有 R_NTP1 的化学、方向及其 notes 变化；本次没有验证生长或必需性变化。**

| ID | 核查声明 | 直接核查与判定 |
|---|---|---|
| O1 | 导出真正使用四物种 ATP 水解 | 直接读取 `R_R_NTP1`：ATP 与 H₂O 系数均 −1，ADP 与 Pi 均 +1；无 H⁺物种参与。四物种均 C_cy。supported。 |
| O2 | 方向以真实化学过程为准 | 从实际 SBML bound 参数解析下界 0、上界 1000，reversible=false；对应 ATP+H₂O→ADP+Pi，而非仅按通量符号称水解。supported。 |
| O3 | 当前实际物种下元素及电荷平衡 | 对实际 SBML 中性公式逐元素重算，C/H/N/O/P 残差全部 0，电荷残差 0；没有改用另一套离子公式。supported。 |
| O4 | 没有附带更改全模型结构 | 逐反应标准化 XML 比较，唯一差异 ID 为 R_R_NTP1；移除该反应节点后，整个 XML 的标签、属性、非空文本和子节点顺序完全相同。因而其他反应、物种、基因、目标、参数、区室及模型元数据保持。此比较忽略无语义的空白缩进，不声称两文件字节相同。supported。 |
| O5 | 本反应 GPR 和来源 annotation 未改 | 独立比较 R_NTP1 的 FBC geneProductAssociation 与 annotation 子树，二者与旧模型完全相同。当前 OR 赋值没有由此次修正升级为独立原生催化的实验证据。supported。 |
| O6 | 修正依据和限制实际写入模型 | 直接读取 notes，含新 ntp1_curation_id、用户授权、取代旧五物种可逆式的说明、KEGG/IUBMB 来源及 native GPR activity unverified / growth and essentiality effects not tested；上游 active/reverse 水解整理来源保留。supported。 |
| O7 | 实际输出身份和旧输入保护成立 | 实测新输出 SHA 与 .build.json 一致；旧比较模型与原始 data/iyali26.xml SHA 与 before.json 一致。独立逐个检查 before.json 中排除 4 个授权改动项后的 65 个文件，全部 SHA 不变，包括旧模型、旧证据和他人既有工作范围内的受保护记录。supported。 |
| O8 | 本轮完整构建符合记录范围 | 直接读取 build_execution.json：第一次完整构建成功，62.4394 秒；命令 offline/no-solve、CoQ9 metadata。实际 .build.json 为 requested_build_complete=true，201/6/9 字段计数，runtime strain overlay=null，独立 V-ATPase 假设关闭；本代理核对实际 metadata 整理、上游方向整理、reaction_selection.py、main.py SHA 与运行清单一致。没有把 solver 接口名称当作已经求解。supported。 |

覆盖：`total claims 8 | audited 8 | supported 8 | unresolved 0 | contradicted 0 | unchecked 0`。该分母只覆盖本轮输出验收；真实蛋白功能、胞内 ΔG、生长和必需性未纳入通过项。软件实施审查另见 [IMPLEMENTATION_AUDIT.md](IMPLEMENTATION_AUDIT.md)。

固定文件身份：

- 对照 `model_metadata_trna_r1025_gpr.xml`：`b83bbc65a00b6dc049501f8c22e82414d6133fb2c0b1845d24f4106e19a4fab0`。
- 新输出 `model_metadata_trna_ntp1_hydrolysis.xml`：`9adfb0f6187360770f1ca38ebbdb7a9e68d99705c6868a91a8e7f583cd7e7e36`。
- 原始构建输入 `data/iyali26.xml`：`5c8c199e2c5b622e97daf2b3500f763f83519fb598702a11dd153052c6a99f9d`。

当前输出属于本次授权的普通完整构建，不将其称为正式发布，也不将另一份既有 V-ATPase 假设模型的条件默认为本次输入。本审计独立解析原 XML 后得到上述结论，没有从根代理生成的比较结论倒填核验结果。

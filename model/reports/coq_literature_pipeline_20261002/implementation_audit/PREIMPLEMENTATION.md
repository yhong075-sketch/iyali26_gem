# 管线接入的只读预审

2026-10-02。复用前轮生物学与候选验收，不新增检索、比对、预测、构建或求解。本文件是实现建议，不是本轮构建已完成的声明。

最小接口为一个默认关闭的 `--coq-literature-revision`。在构建入口先解析其有效依赖为既有 CoQ9 functional GPR 与 C5 GPR，再统一检查 offline、no-solve、metadata、新输出及禁止混合实验 overlay。尾部按既有 R695、R39、文献 R19/R18 顺序应用。构建记录同时保留用户请求、有效依赖、实际原始输入及 `partial_literature_candidate_redox_gap` 状态。直接 API 与 CLI 必须经历同一完整参考阶段，不能因 API 未指定 canonical_copy 而漏掉必要 tRNA 表示。

原始输入应在重计算前锁定；前轮成品仅用于比较与参考身份，不作为原始 starting-model。现有 curation 的 source 字段是前轮审核参考；实际执行 input 由本轮 build manifest 明确区分。

导出根因：COBRApy reader 对有 FBC 插件的物种直接调用 getCharge，未先检查字段是否存在，使成品中的缺失值在重读时变成零。tRNA split 在完整链内新建的20个残基未设置 charge，内存为 None；writer 仅写入非 None 电荷，因此直接首次导出应保持缺省。原始 SBML 根无 SBO，首次 writer 加入 SBO；下一次读取再输出时 document annotation 触发 meta_。这些差异发生于成品读写；新完整链应优先使用原 writer，不能预先添加无必要的全局修复。以上为源码核实后的条件推断，待实际输出确认。

独立验收范围：CLI/API默认关闭与依赖顺序；守卫在读取/写出前拒绝非法组合；原始输入、代码和 curation 身份；完整 SBML 对前轮已审核 v2 逐树比较（只忽略排版、前缀、属性顺序，科学属性不得归一化）；原始缺失 charge 不补零；目标 GPR、公式/电荷与供需保持；候选 partial 状态与未闭合结构明确；无求解/网络泄漏。前轮序列和 AlphaFold 预测证据按既有身份继承。

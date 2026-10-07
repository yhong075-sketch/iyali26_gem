# V-ATPase 12共同AND + 两a候选OR：一次全基因screen

用户2026-09-21明确授权“用新的候选 做一次screen test”。以上轮给出的14位点候选为科学测试假设：共同12 AND (YALI1E12482g OR YALI1F38820g)，同时用于R794/R795，仅为检验假设后果，不构成原生区室或W29 F身份被接受。基因功能/序列/局限沿用 ../vatpase_branch_alphafold_20260918/web_server/GPR_REVIEW.md。

固定底稿为model_metadata_trna_r539_alphafold_labeled.xml，SHA cff96f40368b3de3d502956cf33dda246e049becfa43c1bba12fbd8ded6acd9a；不以历史三基因AND诊断模型和人工二肽供给条件代替。通过现有apply_metadata_reaction_selection读取任务内候选整理JSON生成独立候选SBML，正式整理表/默认管线/已有模型不改。检查两个目标GPR和记录notes以外的模型字段、计量、基因集合与全部边界保持；静态单KO核验12个共同成员关闭两条GPR、两个a单缺失仍真、二者双缺失为假，无双敲LP。

复用已记录SD-Leu、PO1f运行时overlay及冻结context loader，biomass_C最大化，Gurobi原Presolve=0/容差1e-7、1线程、60秒/LP。候选一次WT+1073单KO=1074主LP；未修改底稿同条件对照另1074，合计上限2148，每模型最多600秒。保留R794[0,1000]、R795[0,0]及其余原边界，不加供给、不开放反应。主分类用未舍入KO/本模型WT<0.15，附1/5/10%。失败/非optimal/缺失/非有限/负growth/守恒或边界违反>=1e-7立即保留原值并停止，不重试、不改变基准。

保存两个实际输入SHA、整理规范、构建检查、每条原始growth/status/GPR传导/残差、目标14基因及WT完整通量、全1073基因比对、软件代码及运行条件指纹。实验正例只复用已有冻结标签/ID作开发参考，未标注保持未知，不宣称独立验证或分类准确率。独立代理审核原始输出、输入变化及结果比较。无FVA、结构预测、集群、提交、推送或正式模型接纳。

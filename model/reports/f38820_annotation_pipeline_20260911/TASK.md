# F38820候选功能标记接入

用户明确授权：将YALI1F38820g标记为偏Vph1-like的V0 a亚基候选，注明仍需实验确认，并入管线。该授权接受候选功能元数据，不把液泡定位、同源一对一关系、GPR或essentiality变成已确认事实。

顺序：复用已审序列/AlphaFold证据，在正常整理数据目录保存候选标记；现有genes模块增加最小名称/notes应用步骤，接入真实构建；执行针对性元数据幂等/导出/边界测试；一次offline/no-solve完整构建输出到新文件；比较全部数学语义及其他元数据，独立审阅实际数据和改动后交付。

输入：发布model_metadata_trna.xml仅作前后比较，构建仍读取data/iyali26.xml和已封存reference_pipeline_restore_20260909/research。实际完整SHA记录before.json和新build.json。保留旧发布模型与全部旧证据及同事dirty工作。

成功：新模型gene.name/notes同时出现candidate及experimental confirmation required；不伪造VPH1正式名称，来源和证据级别可读、导出可保留；只改变目标功能元数据，模型数学和其他对象注释保持。

资源：最多30分钟；一次完整构建，限定验证和必要失败修复；GEM优化0、联网0、集群0、新序列/结构计算0、Git提交/推送0。缺输入/不一致/身份guard失败时只停止对应构建，不自动换基线或放宽标准。此前研究审计分母不重算。

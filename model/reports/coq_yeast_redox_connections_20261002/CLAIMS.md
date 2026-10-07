# CoQ 末端反应待审查声明

范围：固定当前 Yarrowia lipolytica 管线局部候选；酵母参照用于反应识别，不将数据库和跨物种证据升级为原生长链底物验证。请独立判断，不预设支持结论。

|ID|待核声明|来源与位置|主任务当前证据判断|
|---|---|---|---|
|C01|YeastCyc 将 RXN3O-75 表示为 DMQ6H2 + 通用还原供体 + O2 → DMeQ6H2 + 通用氧化供体 + H2O，关联 YOR125C—CAT5/COQ7（去甲氧基CoQ羟化酶，数据库注释）。这是数据库反应，不是自由醌醇底物纯酶实验。|primary_sources/YeastCyc_COQ7.html 与 YeastCyc_RXN3O_75.html；yeast_model/REPORT.md|数据库内容可核查；机制外推保留限制|
|C02|Poon1999 在酵母线粒体、短链demethyl-Q3类似物体系检测末步O甲基化；YOL096C—COQ3（CoQ O甲基转移酶，酵母实验支持）缺失失活/回补恢复，NADH必需。作者把还原力解释为生成hydroquinone；分析前用Ce(IV)氧化产物。纯化酶对预还原底物的直接要求来自细菌UbiG实验，不冒称纯酵母酶。|Poon1999 DOI10.1074/jbc.274.31.21665，p21667 Methods、p21668 Fig3旁文字、p21670 Discussion|支持醌醇末步表示及酵母COQ3角色；纯酶与氧化态局限明确|
|C03|RXN3O-102及当前EC2.1.1.64定义末步 DMeQnH2 + SAM → QnH2 + SAH（具体H+随微物种约定）。因此R385的醌醇版本有明确参考依据；YALI1B20835g仅保留为既有COQ3家族/O甲基转移酶候选，原生功能未实验确证。|YeastCyc COQ3；https://iubmb.qmul.ac.uk/enzyme/EC2/1/1/64.html；既有coq_system_gpr_review_20260918/coq_gene_gpr_review.tsv|反应家族支持，Q9原生外推有限|
|C04|IUBMB2024将真核羟化酶分为EC1.14.13.253醌/NADH净式，EC1.14.99.60当前注释为原核醌醇步骤。其酶学基础主要是脊椎/人蛋白，不能把它当作酵母自由醌醇底物的直接否定。Nicoll2024的对应比较终点是NADH消耗，不能扩写为所有条件下无产物。|IUBMB两页；Nicoll2024 DOI10.1038/s41929-023-01087-z Fig5/Results；Lu2013|分类事实支持；物种和终点限制必须保留|
|C05|目前尚未查实可赋予确定酵母基因和电子受体的独立自由DMQH2→DMQ氧化反应。Lu2013证据是人酶结合态底物介导双铁还原，不足以独立证明自由DMQH2氧化酶；DMQH2+O2→DMeQ+H2O仅为机制推演候选。|Lu2013 DOI10.1021/bi301674p，Results substrate-mediated reduction；本轮限定检索记录|未解决，不声明酶不存在|
|C06|不能把成熟QH2氧化直接复制给DMQH2。Padilla2004见DMQ6-only酵母呼吸不足，但bc1装配缺陷使作者明确不能区分DMQ6是否是功能性bc1底物；Bradley2020对YLR290C—COQ11（推定氧化还原酶/CoQ相关调控蛋白，功能注释与假说）提出的是成熟Q6H2氧化假说，未证明DMQH2底物或电子受体。|Padilla2004 DOI10.1074/jbc.M400001200 Discussion；Bradley2020 PMC7196636 Discussion|不足以赋值，也不能当作纯酶底物阴性证明|
|C07|固定Yeast-GEM9.1.1中r_0963使用醌→醌；r_0532使用DMeQ6醌+SAM+H+→Q6H2+SAH，元素残差0、电荷残差−2。完整物种扫描未找到三种所需中间醌醇；没有可照搬的DMQH2氧化连接。|yeast_model/model_extract.json、verification.json及固定官方XML|限定版本文件事实，不外推所有酵母模型|
|C08|可提出待审连续路线：R18→DMQ9H2→R695醌醇版本→DMeQ9H2→R385醌醇版本→Q9H2→现有成熟Q氧化。R695对自由DMQ9H2的适用性及实际供体仍未确立，因此本轮不能称这条Q9路线已被生化验证或模型已闭合。|C01–C07；当前候选仅R305已连接成熟Q池|整体Q9反应候选未解决；本轮不修改模型|

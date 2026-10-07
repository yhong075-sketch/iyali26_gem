# 酿酒酵母反应参照：氧化态连接审查

范围：只读核对 Yeast-GEM 9.1.1 固定提交、仓库 9.0.2、YeastCyc 22.5 和原始酵母实验；没有改模型、GPR 或运行求解。日期 2026-10-02。目标是寻找 DMQ9H₂ 接往末端 CoQ 合成的依据，不预设必须存在独立氧化酶。计划先逐计量比较，再核对数据库基因与原始实验，完成限定取证即停止。

## 结论

酵母参考给出两条不同层次的信息：**YeastCyc 已把下游合成画成醌醇连续路线；Yeast-GEM 沿用醌类中间体，不能提供缺失氧化反应。** 数据库路线足以提出重审 R695 化学形式的候选，但 generic electron donor 及实验限制仍不能直接生成已验证的 NADH-GPR 净反应。

|目标|Yeast-GEM 9.1.1 / 9.0.2 实际表示（均 m 线粒体，界限 0–1000）|YeastCyc 22.5 表示|
|---|---|---|
|R18 C-甲基化|r_0021：DDMQ6（醌）+ SAM → DMQ6（醌）+ SAH + H⁺|RXN3O-54：DDMQ6H₂ + SAM → DMQ6H₂ + SAH + H⁺|
|R695 C-羟化|r_0963：DMQ6 + NADH + H⁺ + O₂ → DMeQ6 + NAD + H₂O|RXN3O-75：DMQ6H₂ + 还原型通用电子载体 + O₂ → DMeQ6H₂ + 氧化型通用电子载体 + H₂O|
|R385 末端O-甲基化|r_0532：DMeQ6（醌）+ SAM + H⁺ → Q6H₂ + SAH|RXN3O-102 / EC2.1.1.64：DMeQ6H₂ + SAM → Q6H₂ + SAH + H⁺|

这里 DMeQ 指去甲基辅酶Q，DMQ 指去甲氧基辅酶Q，两者不是同一物质。S. cerevisiae 天然侧链主要 Q6，不能把其物种ID原样放进 Q9 模型。

**文件级反证：** Yeast-GEM r_0532 混用醌底物与醌醇产物。由文件实际分子式、电荷求和，元素守恒但电荷净残差为 −2；完整物种扫描没有 C37H56O3、C38H58O3、C38H58O4 三种所需中间醌醇。因此不能把该模型作为氧化态与电子平衡已解决的证据。r_0021、r_0022、r_0963 元素/电荷残差均为0，但守恒不等于生物化学已验证。

成熟辅酶Q的氧化已经存在：酵母 r_0439（ubiquinol:ferricytochrome c reductase）消耗成熟 Q6H₂，生成 Q6；目标模型 R305 消耗成熟 Q9H₂。该事实不支持它们具有 DMQH₂ 底物特异性。r_0439 在文件中另带非整数质子系数，不作为本次可复制化学模板。

## 基因与证据层次

|系统ID|已核实名称与功能|证据和模型角色|
|---|---|---|
|YML110C|COQ5；辅酶Q环C-甲基转移酶|YeastCyc curated annotation；对应 RXN3O-54。模型将其包含于下述六成员AND，不能据此称每个成员直接催化此步。|
|YOR125C|CAT5，别名COQ7；辅酶Q去甲氧基中间体羟化酶|YeastCyc curated annotation；RXN3O-75 使用醌醇和通用供体。该数据库方程不单独证明纯酵母酶对自由DMQ6H₂的体外活性。|
|YOL096C|COQ3；辅酶Q生物合成O-甲基转移酶|酵母实验支持O-甲基化作用和线粒体内膜基质侧定位；对应 RXN3O-102。|
|YDR204W|COQ4；辅酶Q生物合成复合体相关蛋白|本次仅记录模型/GPR归属，未新增催化作用判断。|
|YGR255C|COQ6；辅酶Q生物合成黄素羟化酶|本次仅记录模型/GPR归属，未独立重审其位点功能。|
|YLR201C|COQ9；辅酶Q生物合成复合体相关脂质结合蛋白|本次仅记录模型/GPR归属，未新增原生依赖判断。|

Yeast-GEM r_0021、r_0022、r_0963、r_0532 共享 `YDR204W and YGR255C and YLR201C and YML110C and YOL096C and YOR125C`。r_0021/r_0022 引用 PMID15792955（复合体论文），r_0963/r_0532 未带直接PubMed引用；整个复合体AND不是逐反应底物特异性的实验确认，不直接移植。

## 原始证据与可用来源

Poon等1999用短链demethyl-Q3类似物和离体酵母线粒体测末步O-甲基化：删COQ3丢失活性，回补恢复；加入NADH才检测到活性。讨论将NADH解释为产生hydroquinone的还原力，且纯化细菌UbiG需还原型底物。故它支持末步以醌醇甲基化表示；没有鉴定连接DMQH₂→DMQ的独立氧化酶。物种、短链类似物和粗线粒体条件不能忽略。实际已读PDF第21668、21670页文字；PDF截图工具失败，未宣称视觉核验图3/5。

- [Yeast-GEM 固定提交文件](https://github.com/SysBioChalmers/yeast-GEM/blob/2d594ae1c4a2d550ccef120d96a58c7bbf586255/model/yeast-GEM.xml)，本地文件完整SHA在 model_extract.json，核对与既有官方获取记录一致。
- [YeastCyc COQ5](https://pathway.yeastgenome.org/gene?id=YML110C&orgid=YEAST)：RXN3O-54链接及方程已读。
- [YeastCyc CAT5/COQ7](https://pathway.yeastgenome.org/gene?id=YOR125C&orgid=YEAST)：RXN3O-75链接及方程已读；独立反应详情页无法取回，原始支持文献仍未确定。
- [YeastCyc COQ3末步反应](https://pathway.yeastgenome.org/YEAST/NEW-IMAGE?object=RXN3O-102&type=REACTION-IN-PATHWAY)：ID与EC及参考文献已读；方程另由COQ3基因页检索结果核对。
- [Poon等1999，JBC 274:21665–21672](https://www.biochemistry.ucla.edu/Faculty/CClarke/pdf/21665.pdf)，PMID10419476，DOI10.1074/jbc.274.31.21665，第21668页 Fig.3旁文字，第21670页 Discussion。

审查限制：本包为原始取证供总协调独立审计，不把本代理自检称独立来源审计。网络网页通过web工具读取；本地HTTP抓取域名解析失败保存在 retrieval.json，未保存完整网页快照。对新酶没有作身份或序列推断，不需要扩展AlphaFold。未解决的具体问题是COQ7自由醌醇底物/供体净计量，不能为闭合通路编造氧化酶。

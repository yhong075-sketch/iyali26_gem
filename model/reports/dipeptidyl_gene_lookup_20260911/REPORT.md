# YALI1 二肽基肽酶候选查询

核验日期：2026-09-11。范围为已有注释、编号映射及原始文献的查阅；没有进行新的酶身份鉴定、BLAST、结构推断、求解或模型修改。

| YALI1 ID | 名称与功能 | YALI0 对应 ID | 证据等级及局限 |
|---|---|---|---|
| YALI1B04274g（原表 YALI1_B04274g） | 正式基因名未核实；二肽基肽酶 IV 候选 | YALI0B02838g | 已发表 W29 注释及映射；UniProt A0A1D8N680 未审阅，名称来自 ProtNLM、S9B 家族/液泡膜为自动推断。机构报告将旧 ID 称为 putative yylDAP。没有据此确认原生催化活性、底物或定位。 |
| YALI1B25603g（原表 YALI1_B25603g） | 正式基因名未核实；二肽基肽酶 III 候选，注释功能为从至少四残基肽链释放 N 端二肽 | YALI0B19580g | 已发表 W29 注释及映射；UniProt A0A1D8N8I3 未审阅，EC3.4.14.4、M49 家族及细胞质位置为自动推断。没有原生功能或具体底物验证。 |

## 来源与解释边界

1. Magnan et al. (2016), DOI 10.1371/journal.pone.0162363，S2 Table / DOI https://doi.org/10.1371/journal.pone.0162363.s004。本地原表 sheet1 第1276和2144行逐格核对；对映射编号去下划线仅为显示归一化。
2. UniProt 既有缓存 A0A1D8N680（entry38，sequence1，887aa）、A0A1D8N8I3（entry31，sequence1，707aa）。完整条目、来源和序列SHA见 source_extracts.json。新 API 请求失败，不称最新实时复核。两条原表编码坐标长度与后续 UniProt 蛋白长度存在版本差异；本轮不据此推断完整序列跨版本相同。
3. IPN 项目20060695机构报告：https://www.sappi.ipn.mx/cgpi/archivos_anexo/20060695_3686.pdf。PDF第21页描述同源检索，第22页PCR图注包括MI-IPN-1和W29，第24页图2caption称 putative yylDAP (YALI0B02838g)。读取的是Web提取文本，PDF截图失败；不宣称已目视核实图中数据。该材料不等同同行评审的基因功能验证。第二报告20071132无法取得，排除。
4. Hernández-Montañez et al. (2007), DOI 10.1111/j.1574-6968.2006.00578.x，PMID17227470。已读摘要支持解脂耶氏酵母三种细胞内蛋白酶活性及氨肽酶纯化；摘要本身不提供具体YALI1编号。全文未读取，不能断言全文完全没有编号或基因证据。

## 本轮结论

能提供上述两个已注释的YALI1候选。其中B04274与早期yylDAP候选的联系更直接；B25603是另一类可能释放二肽的酶候选。候选注释不是实验确认，不确定其是否生成Gly-Asp、Gly-Glu、Ala-Gly或Gly-Pro；也不确定其与R795所表示液泡质子泵的生理依赖。不得直接把候选添加为这四种产物的生成GPR。独立核查范围和未决项见AUDIT.md。

# 四条二肽水解底物来源：独立来源审计

2026-09-23。审核者 `/root/dipeptide_source_audit`。只读原 XML、历史运行 manifest、原始来源与本轮提取；本轮新增 LP/模型变更为 0。以下为独立来源核查，不把审查者与主代理的同意作为证据。

## 固定输入与复核方法

- 9 月 21 日候选：`../screen_vatpase_candidate_20260921/candidate.xml`，SHA256 `0fcc2f6ff05124b91977c6f96f78cbe5906a3ccbf3f936f95990afe02b540353`；原 XML 2314 反应、1877 物种。
- 实际 SD-Leu/PO1f 运行清单：`../screen_vatpase_candidate_20260921/candidate/run_manifest.json`，SHA256 `348ae99577e76c851340a871f20cb200807c7f38f29b06cef42825294efbeb14`。使用其 `base_bounds` 和 `medium.active_medium`，不把 XML 默认边界充作运行时培养基。
- 审核者从原 XML 独立重建 12 个目标二肽物种的全部反应邻接、计量和三池合并守恒，与 `model_evidence.json` 四条链逐项断言匹配；该提取文件记录的 5 个来源 SHA 均重新核算匹配。
- `../screen_vatpase_candidate_20260921/mechanism_details.json` 的四链机制与独立重建一致。模型代谢物名称 `Gly_Glu` 对应这里的 Gly-Glu；未因连字符/下划线差别遗漏其 3 个区室物种。
- 原生功能候选仅复用 `../dipeptidyl_function_prediction_20260911/REFERENCE_AUDIT.md` 的既有限定；本轮没有重做序列/结构鉴定，也不据此接受新 GPR。

## 核心声明与判定

| ID | 有限范围内的声明 | 直接核查依据 | 判定与限制 |
|---|---|---|---|
| C1 | R2021、R2029、R2034、R2039 依次消耗 Gly-Asp、Gly-Glu、Ala-Gly、Gly-Pro 和水，产生对应游离氨基酸；实际计量全部在液泡 | 原 XML 四反应及相应 species 的 `compartment=C_va` | supported。反应名称写有 cytosol，但不能用名称覆盖实际区室。模型反应赋值不等于已验证的酶/区室生物学。 |
| C2 | 当前固定模型没有这四种二肽的内源净生成反应 | 穷举每种二肽的胞外、胞质、液泡三物种的全部反应；每链只有 exchange、两段 transport、单向 hydrolysis | supported，限定该模型。物种节点存在不等于有初始库存；没有建模来源也不证明真实细胞不能产生。 |
| C3 | 当前运行培养基不提供这四种二肽摄入 | R2018/R2026/R2031/R2036 均为 `[-0.0,1000]`；计量为消耗胞外二肽，负向才是摄入；四者不在 active_medium | supported。不能误写为 `[0,0]`；分泌方向仍允许。 |
| C4 | 在这组结构与边界下，四条水解通量必为 0 | 三个区室的二肽守恒相加，依次得到 `-v_EX-v_hydrolysis=0`；两项均非负 | supported，静态守恒结论，无新求解。运输可逆不构成二肽来源；开放其他不改变该守恒的反应也不会凭空供料。 |
| C5 | 9 月 14 日 pool 是独立人工连续来源，不是证明原生已有二肽 | `../dipeptide_pool_diagnostic_20260914/run.py:100–108`：新增 `Reaction(lower_bound=0,upper_bound=1)`，只生成对应胞质二肽，注释明确 artificial source；旧报告说明非有限初始浓度 | supported。当前 9 月 21 日候选没有这 4 条二肽 POOL；模型存在其他脂质 xPOOL，不能笼统声称完全没有 POOL。 |
| C6 | 所引用的标准 CSM-Leu 补充配方未列这四种二肽 | Sunrise Science CSM-Leu，货号 1005-010，“Recipe in mg/L”，详见来源 S1 | supported，仅配方级证据。页面为单 Leu 缺失配方；不是 Leu-Lys 配方，也不是实际培养批次的化学测定。 |
| C7 | 解脂耶氏酵母具有细胞内二肽基氨肽酶活性的实验报道，支持内源肽加工可作为候选来源 | Hernández-Montañez 等 2007，出版商原始摘要，S2 | supported：物种级活性事实；“可作为候选来源”为受限推断。最高活性见含 peptone 培养基、稳定期，非当前 PO1f/SD-Leu 四二肽测量。 |
| C8 | 这四种精确的游离二肽在 PO1f/SD-Leu 的液泡中已有并持续供应，足以支撑模型水解通量 | 上述来源及既有功能审查不提供该条件、四产物和区室的直接测量 | unverified。未证实不等于不存在；不能据此接纳无质量来源的永久 pool 或把浓度/通量填为默认正值。 |
| C9 | 当前 R1363 关闭是独立的液泡水供给断点 | 独立穷举原 XML 物种 `m1384[C_va]` 全部邻接，仅有 R1363 系数 +1 和四条水解系数 −1；运行清单 R1363 `[0,0]`，四条水解 `[0,1000]` | supported。由守恒得四非负水解之和为 0；即使有二肽输入，在保持该水连接关闭且不增加其他水来源时仍阻塞。不能把该建模缺口当作细胞不含水。 |

**审计覆盖：9 total | 9 audited | 8 supported | 1 unresolved | 0 contradicted | 0 unchecked。** C8 是明确保留的问题，不能把 9 项均已检查写成 9 项科学结论均成立。

## 实际打开的外部来源

**S1：供应商原始配方，非研究论文。** [Sunrise Science CSM-Leu Powder, 10 grams](https://sunrisescience.com/shop/growth-media/amino-acid-supplement-mixtures/csm-formulations/csm-leu-powder-10-grams/)，2026-09-23 读取。货号 1005-010；定位：Recipe in mg/L / Amino Acids & Supplements。列有单体氨基酸与腺嘌呤、尿嘧啶，不列这四种二肽；说明与 YNB、氮源、葡萄糖配成培养基。未检验实际培养瓶、杂质或非标准补料。

**S2：同行评议原始研究，当前仅核摘要。** Hernández-Montañez Z 等，*The intracellular proteolytic system of Yarrowia lipolytica and characterization of an aminopeptidase*，FEMS Microbiology Letters 268(2):178–186 (2007)，DOI [10.1111/j.1574-6968.2006.00578.x](https://doi.org/10.1111/j.1574-6968.2006.00578.x)，PMID 17227470。2026-09-23 实际打开 [出版商摘要](https://academic.oup.com/femsle/article-abstract/268/2/178/656541)，跳转至 `oup.silverchair-cdn.com/article-minimal/656541`；定位 Abstract 第 2–4 段。论文在细胞可溶提取物中测三种蛋白酶活性，并以底物显色凝胶得到各自活性条带；DAP 的存在不能替代四精确产物、胞内浓度或液泡定位证据。97 kDa、Lys-pNA 动力学等具体数值属于进一步纯化的 aminopeptidase，不能套用到 DAP。全文未获得，本轮不推断样本量/重复数、具体菌株或未展示的底物谱。

PubMed 直接打开出现验证页，DOI 直接打开失败；随后成功取得 S2 出版商摘要。没有绕过访问限制，也没有把只检索到的标题当作内容验证。

1989 PMID 2649495 和后续 1997 PMID 9353927 的原始摘要另经检索得到：前者支持分泌蛋白前体发生逐二肽加工，后者提示成熟蛋白分泌不必依赖该加工。其全文及精确产品未在本轮核验，不属于核心 9 项证据的新增直接产品/通量证明，不能外推为这四种二肽必定存在或生长必需。

## 对本轮报告的复核

已读取本目录 `REPORT.md` 全文。模型、配方与历史人工源部分均在上述审计范围；2007 来源只用于物种级细胞内蛋白酶活性和候选供给机制，未把酶活性写成四精确游离产物或液泡供应已证实。报告将内源供给标为候选并保留 PO1f/SD-Leu/区室限制，适当。全文未获取，关于原论文没有证明精确产物的句子应明确限定为“本轮读取的摘要未提供这些证据”，不能声称查尽未读全文。审核者已向主代理提出这一措辞修正。

## 接受的解释范围

可以说：“生物上二肽不一定需要外加；蛋白或肽的加工是可能来源。当前模型却仅编码了这四种二肽的外源利用链，所用 media 又没有允许其摄入，因此在模型中没有可持续供给。”

不能说：“酵母本来就有足量这四种液泡二肽”“所有培养基都必须额外提供它们”“模型缺失来源证明生物没有内源生成”“自由二肽 pool 已有实验依据”。这一区分不授权新的模型/培养修改。

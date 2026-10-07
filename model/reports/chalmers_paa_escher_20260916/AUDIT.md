# Chalmers PAM/PAA Escher 图：独立审计

核验日期：2026-09-16。审计者 `/root/source_audit`。范围仅为两份已保存官方模型和既有已审来源的静态复核，以及本轮生成地图的对应关系；没有新搜文献、运行比对/结构预测/代谢求解或修改模型。本文件是新图的审计，不覆盖或重写2026-09-14来源审计。

## 固定输入

| 输入 | 模型内部版本 | 完整SHA256 |
|---|---|---|
| `../phe_paa_sc_w29_20260914/sources/Chalmers_Yl.xml` | iYali v4.1.2 | `c0b54165301bfba7edc7083efee1106309f049494f3a3f35fb436306d4fc86ae` |
| `../phe_paa_sc_w29_20260914/sources/Chalmers_Sc.xml` | Yeast-GEM v9.1.1 | `30842b15eefb0ef7e36cbdea86a9efddfacf69a871c8b054165faa9af9f6c8eb` |

两者与2026-09-14已审输入完全相同；来源和原研究的证据层级复用该目录SOURCE_AUDIT.md。XML中的`R_`/`M_`前缀与图中的较短ID如有归一化，应保留可逆对应，不借图改动原模型。

## 模型事实复核

| 图应保留的事实 | 本次独立核查 |
|---|---|
| iYali的胞质PAA没有模型去路 | `M_s_1321`完整行仅`+v_R_y000185`，该反应边界[0,1000]。物种为非边界物种，标准稳态直接强制其通量为0。 |
| Yeast-GEM的PAM和PAA分支仍断开 | `s_3867`完整行仅`−v_r_4227`；`s_1321`完整行为`v_r_0185+v_r_4227`。`r_4227`边界[−1000,1000]仍不能解除PAM单行守恒，两通量均被强制为0。 |
| 上游选定转氨反应可逆 | `R_y000851`/`r_0851`均为Phe[c]+2-oxoglutarate[c] ⇌ phenylpyruvate[c]+glutamate[c]，边界[−1000,1000]，四项系数均1。 |
| 苯丙酮酸脱羧到苯乙醛是正向 | `R_y000854`/`r_0854`：phenylpyruvate[c]+H+[c] → phenylacetaldehyde[c]+CO2[c]，边界[0,1000]。 |
| PAA形成反应须保留辅因子和质子计量 | `R_y000185`/`r_0185`：phenylacetaldehyde[c]+NAD[c]+H2O[c] → PAA[c]+NADH[c]+2H+[c]，边界[0,1000]。 |
| 苯乙醛胞质/胞外/线粒体支路不是PAA运输 | `R_y002002`/`r_2002`为[c]⇌[e]；`R_y002003`/`r_2003`为[c]⇌[m]，均[−1000,1000]。`R_y002001`/`r_2001`仅消耗[e]苯乙醛，边界[0,1000]，即排到模型外。 |
| 存在其他苯乙醛消耗分支 | 两模型`0169`为胞质NADH还原生成2-phenylethanol；`0170`为线粒体同类步骤；Sc还含胞质NADPH步骤`0171`，均正向边界。若图有意省略，不能让读者以为这些全不存在。 |
| 图不是完整苯丙氨酸/苯丙酮酸网络 | 两模型还含`2117`（Phe/pyruvate转氨）、`2118`（Trp/phenylpyruvate转氨）及`0938`（prephenate生成phenylpyruvate）。仅画`0851`允许，但须明确选定支路范围。 |

红色或“阻断”标记只能表示上述完整物种行推导的模型稳态约束，不能称本轮FBA通量、实验没有活性或真实生物学不存在去路。黑/蓝的其他反应同样不是已求得非零可行通量。

基因名称和功能层级沿用已审记录：Sc YGL202W/ARO8为芳香族氨基酸转氨酶，YDR380W/ARO10为2-氧代酸脱羧酶，YMR169C/ALD3及YMR170C/ALD2为醛脱氢酶；YDR242W/AMD2仅为推定酰胺酶，其PAM GPR是模型赋值而不是底物活性实证。Yarrowia旧位点按模型原ID呈现；未经核实的正式名称保持未知，GPR集合及OR仅表示模型赋值。运输/交换反应没有GPR，不添加PDR12或W29候选来替代空规则。

## 地图与显示核对

已独立读取`w29-map.json`（8反应）和`sc-map.json`（10反应），用源XML重建各反应的计量、FBC边界和完整布尔GPR后逐项核对。18个反应的67项计量、可逆性、GPR、103段连接端点及`evidence.json`的区室/名称/边界/系数均通过；`M_`/`R_`只作可追溯ID归一化。每个反应的代谢物端点与其计量集合一一对应。两地图涵盖PAA/PAM与苯乙醛[c]/[e]/[m]全部关联反应，完整行与源模型完全相同。

详情中未核实的蛋白正式名称保持未知并标`model/GPR assignment only`，AMD2明确为推定酰胺水解酶且PAM特异活性未确认。地图没有新增运输GPR，也没有把无GPR变成PDR12候选赋值。

首次目视两份PNG后，向根代理提出以下显示修正：将苯乙醛胞外运输线绕开PAA节点；修正Sc辅因子标签重叠；将Sc的PAM注释收窄为“没有独立生成或摄入”，PAA注释收窄为“没有独立排出或继续降解”，以免否认可逆r_4227在局部方向上的生成/消耗能力。这些修改均已纳入最终地图。

最终两份PNG已再次目视：PAA在单独红色分支，胞外/线粒体运输线不再穿过PAA节点；Sc两个胞质2-苯乙醇分支以“同一胞质池”注明重复绘图节点，NADH/NADPH标签已分开；PAM可逆箭头、上下游限制措辞和完整守恒式均可见；单底物交换旁标“→ 模型外”。红色是源模型守恒证明，青色是已有运输/交换，灰色是其余步骤；没有把这些颜色冒充求解所得通量。

最终布局修改后重新执行上述源XML→地图静态断言，18个反应及5类目标池的完整行再次通过；HTML嵌入的两套数据与最终`evidence.json`完全一致，后者所含地图与两个独立地图JSON完全一致。显示样式通过反应类别着色，未用虚构flux overlay。

根代理`render-check.json`记录Escher 1.8.1/headless Chrome的两模型渲染、目标反应详情、布局拖动、切换后编辑保留、下载JSON计量/GPR保留，以及空错误列表。本审计读取并核对该记录的科学文字与来源；**交互操作由根代理执行，本独立审计不声称亲自重复了浏览器点击测试**。

## 最终覆盖与文件锁定

| 审计项 | 范围 | 结果 |
|---|---|---|
| 1 | 两个输入SHA与已审版本身份 | 通过 |
| 2 | iYali PAA完整行及稳态断点 | 通过 |
| 3 | Yeast-GEM PAM/PAA完整行及可逆边界下断点 | 通过 |
| 4 | 0851转氨与0854脱羧的计量、方向 | 通过 |
| 5 | PAA形成步骤辅因子及2H+系数 | 通过 |
| 6 | 苯乙醛c/e/m运输、胞外交换的底物与方向 | 通过 |
| 7 | 苯乙醛还原分支及完整目标行覆盖 | 通过 |
| 8 | 上游仅选一条转氨支路、2-苯乙醇下游未展开的范围提示 | 通过 |
| 9 | 18反应、67计量项、103段连接端点的源模型一致性 | 通过 |
| 10 | 保存边界/可逆性与全部布尔GPR一致性 | 通过 |
| 11 | 已核实名与未知名称分开；AMD2/PAM及OR证据边界 | 通过 |
| 12 | 最终两PNG的方向、标记、节点连线和文字修正 | 通过，独立目视 |
| 13 | HTML嵌入数据、evidence与独立地图JSON对应 | 通过，独立静态核对 |
| 14 | 交互测试记录及科学详情表述 | 通过，已读根侧执行记录，未独立重跑UI |

**14项已审/14项；14项在所列验证层级通过，0项未决、0项矛盾、0项未检查。** 该覆盖仅指图与既定源模型事实的一致性；生物学功能、GPR成员独立催化、化学整理、全模型可行通量及模型修改接受均不在本图审计中得到新增验证。

最终文件SHA256：

| 文件 | SHA256 |
|---|---|
| w29-map.json | `82288435547df26d5a6d95093d3943bd907f5a5ae2931ead40415e8a63b008a6` |
| sc-map.json | `4a5f4812a0f32be69c6029a25b97a695589654f26fe5c7a0b6b0313cc713d430` |
| evidence.json | `0d97e6614c7a6dcbb7a01d9b8a7c3e21f7e1e6c71dc491e5a16780c979b339de` |
| chalmers-paa-escher.html | `2c61b1d23313c375718666077e0202f0c19a122e13e0d251233befdfe5f37025` |
| w29-escher.svg | `b1ddcd635fd7f78422aa0f5951647b3dd2fc484318c0845ac334004c6d04df8e` |
| sc-escher.svg | `478219784b540a533f4be6ec15111c2430aa4791d086dbb95ba42f04dc524843` |
| w29-escher.png | `7db9b4a4aba524b1488f8ae0ee8153a72ddc79a771e11a79b111f79e03e2b9e3` |
| sc-escher.png | `a87655640c50c94179f0b757180059e215be8fb4902cc0c1a8112d8d72f0fff9` |
| render-check.json | `f0111d6a584d612979d6747cb9d589d02216f9d5a6a6f1e18cd7fa4261d4ca2c` |

## English edition addendum — 2026-09-16

The Chinese section above is preserved unchanged. Its pre-addendum SHA256 was `e1f3c405062d0f2e327c3f29bbbc68537ec68f8e6b1ac6f38ea3dfa62af0e4b3`; all nine locked Chinese artifact hashes were independently rechecked and remain unchanged. This addendum audits translation and presentation only. No literature search, sequence comparison, structure work, flux calculation or model change was performed.

The English maps retain the exact 18 reaction identifiers, all metabolite coefficients, reaction reversibility, Boolean GPRs, gene identifiers, node identities and segment connection endpoints. Layout coordinates and displayed text may differ. The English evidence retains the source model identities, complete pool rows, blocked/transport classifications, stored bounds, original model reaction names and equation terms. The English HTML's embedded datasets exactly match `evidence-en.json`, whose two maps exactly match the English map JSON files.

The final wording preserves the relevant limits: the W29 PAM pathway is **not included in this model**, rather than claimed absent biologically; Sc lacks **independent** PAM production/input and an independent PAA outlet; reversible hydrolysis remains constrained to zero by the full steady-state balances. Arrows indicate permitted directions, not calculated flux. Gene details retain model-assignment-only status where appropriate. AMD2 remains a probable amidase with phenylacetamide-specific activity unconfirmed, and OR does not become evidence of independently validated catalysis.

Both final English PNGs were independently inspected. Phenylacetaldehyde export remains visually distinct from the PAA branch, compartment labels and cofactor labels are readable, the repeated cytosolic pools are identified by compartment and the shared-ID legend, and both Sc balance equations remain visible. No source-conflicting display issue remains.

`render-check-en.json` reports 8/10 rendered reactions, matching scientific reaction details, drag support, preservation of layout edits on switching, preservation of stoichiometry/GPR on download and no browser errors. This independent auditor read that execution record; **the browser actions themselves were performed by the root agent and were not independently repeated here**.

| Additional audit item | Verification level | Result |
|---|---|---|
| EN1 — Scientific fields and connections unchanged by translation | Independent static comparison of both map/evidence pairs and locked Chinese hashes | Pass |
| EN2 — Scoped pathway/flux/GPR/gene-evidence wording | Independent reading of English map labels and evidence fields against the locked source conclusions | Pass |
| EN3 — Final English presentation | Independent visual inspection of both final PNGs | Pass |
| EN4 — Embedded data and browser-check record | Independent embedded-data equality check; review of root-performed browser results | Pass at the stated levels |

**Final coverage: 18/18 declared items audited; 18 supported at their stated verification levels, 0 unresolved, 0 contradicted, 0 unchecked.** The English additions do not expand biological or whole-model validation beyond the Chinese audit.

| English artifact | SHA256 |
|---|---|
| w29-map-en.json | `c34e518ad7128a62217843b5b02676cfdcdbd1764343e1bccc96c5789c418f5b` |
| sc-map-en.json | `c110011bcbf0573e1dceaa62411567bbf0dca821b577bbe12286976eec473cf3` |
| evidence-en.json | `a7bda1426ce14763cd2906f2ed91d267461cd289d244572fcfe552c6df8c4b17` |
| chalmers-paa-escher-en.html | `296821072f763ffbcd4bc247af46902d33ae28e8f5d145deb9630485d1dc786b` |
| w29-escher-en.svg | `b8659ec57865ef396872470390df6a38c868bfab51dea0066e4ef8683efb4098` |
| sc-escher-en.svg | `74d43bf3a6255f026baccf6f4967a39b35ea26e461882aede5a36e707483ea5e` |
| w29-escher-en.png | `f625601a6447795d6426071b751a086ef64a0e59f1b370c6b7eb953dfef07275` |
| sc-escher-en.png | `334a893e2eee8936b5d5924e96a85872c895d6ebf1011a390d74b9d2710ee6eb` |
| render-check-en.json | `a6bdd82fc6b1db39b8af970244f83adb77ec161332f16dc13efc7b48fb5a05c2` |

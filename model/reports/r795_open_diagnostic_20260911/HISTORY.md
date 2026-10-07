# R795 关闭边界的来源追溯

核验时间：2026-09-11 01:20 UTC（当地 2026-09-10）。只读调查；未求解、未改动模型／代码／整理数据、未联网、未写 Git。当前 HEAD：`ff36d87eb6c8f933dfc4413f43a2fcefcf0eea27`。

**结论：当前关闭设置继承自导入本仓库时已经关闭的原始模型。本仓库可读取的历史没有给出“原作者为什么关闭”的直接解释。** 最早可核实时间为 2026-03-19；不能把本轮开放实验的结果倒推成当年的设计动机。

## 最早可核实边界

1. `0004d768d739f7519742ebfc08ecc290aed48df3`，2026-03-19 13:10:19 -07:00，提交信息 `feat: add initial structure for the model repository`：`model.xml` 是 95 反应模板，不含 R795。
2. `05ae6e16955a180cdb6416143312ad806ada0440`，2026-03-19 14:12:05 -07:00，提交信息 `chore: Updating structure, adding new update script`：新增 `data/iyli21.xml`；此文件与当时 `model.xml` 字节相同。R795 已名为 `V-ATPase, vacuole`，上下界均引用 `FB3N0`，其数值为 0；该反应无 notes 或 annotation。
3. 当时新增的 `scripts/update_model.py` 第 11、20–29 行只从 `data/iyli21.xml` 读取并写回 `model.xml`；更新函数均是注释占位。未见关闭 R795 的设置语句或理由。
4. 当前 `data/iyali26.xml` 来自后续输入文件整理；在 2026-07-23 的 `d180d902781964ed4f5aa4925fe9ad5052bc14cf` 及平行历史 `e123e04b302f4b0719f543608c924b76c1b71c32` 中即与当前文件完整 SHA 相同，R795 仍为 [0, 0]。

补充范围检查：针对 `model.xml`、`data/iyli21.xml`、`data/iyali26.xml` 执行 `git log --all --format=%H -- <上述三个路径>`，取得 45 个路径相关 revision；读取各 revision 可访问的这三个路径并按 Git blob 去重，得到 23 个不同 XML，其中 22 个含 R795。逐一解析反应边界及对应参数，22 个全部为 [0, 0]；余下一个是上述 95 反应模板。**这覆盖的是本仓库目前可读取的这些历史对象，不包括未提交文件、上游未取得的版本或不可达 Git 对象。**

## 当前构建如何保留它

- 原输入 [`data/iyali26.xml`](../../data/iyali26.xml:16747) 的 R795 两个边界均是 `FB3N0`；参数定义在第 1903 行，值为 0。
- 旧目录模型 [`model.xml`](../../model.xml:98580) 与已发布 [`model_metadata_trna.xml`](../../model_metadata_trna.xml:101055) 也均为 [0, 0]。
- 当前管线 [`main.py`](../../scripts/gem_annotate/main.py:122) 从指定输入读入这些边界。
- [`data/metadata_reaction_selection.json`](../../data/metadata_reaction_selection.json:11971) 对 R795 的 `fields` 只有 `stoichiometry`；`before.bounds` 与 `after.bounds` 均为 `[0.0, 0.0]`。这一轮把胞质 H+ 系数从 −1 改为 −2，未选择修改边界。7 项真正选择 bounds 的反应是 R2041、R_CAT2p、R_NDP1、R_NTP1、R_NTP3pp、R_NTP7、R_PGAM1_PhosHydro，不含 R795。
- 共享应用函数 [`reaction_selection.py`](../../scripts/gem_annotate/reaction_selection.py:61) 只应用 `fields` 指定的字段；第 71–72 行虽有通用 bounds 赋值，但 R795 不进入该分支。两个实际构建记录的 `reaction_selection.items[item=R795].fields` 均只记录 stoichiometry applied。
- 最新 GPR 假设整理数据 [`vatpase_gpr_hypothesis.json`](../../data/reference_build/curation/vatpase_gpr_hypothesis.json:71) 把 R795 的 [0, 0] 列为前置条件。应用函数 [`patches.py`](../../scripts/gem_annotate/patches.py:79) 校验边界，最终只赋值 GPR 和 notes。因此本次共同 AND 也未把原本开放的 R795 关闭。

## “为什么关闭”目前能说到哪里

在上述原始反应记录、最初导入脚本和提交说明中，没有找到与 R795 关闭关联的理由。`HISTORY.rst` 仅有第一版发布的模板条目；最初 README 未给出这条反应的处理理由。对现存 scripts/data 的精确 R795 名称检索以及 scripts、README、HISTORY、docs 和 metadata 选择数据的 Git `-G 'R795|V.ATPase, vacuole'` 检索，也未发现更早的关闭说明。

因此：

- **可确定：** 关闭在最初导入含 R795 的模型时已存在，后续构建沿用。
- **尚未确定：** 上游为何采用零边界，以及这是人工策划、条件设定、历史修补还是导出遗留。
- **不能当作已知原因：** 为了避免 ATP／质子循环、认为液泡泵不存在、根据本次 essentiality 结果关闭、或者本次 GPR 修订意外关闭。前两种未找到直接历史证据，后两种与已核验时间链不符。

要进一步确定原作者意图，需要导入前的模型来源、上游 curation 记录或作者说明。当前文件曾名 `iyli21.xml` 只是溯源线索，不能仅据文件名确认具体上游版本或论文。本轮未获取这些外部资料。

## 交给独立审计的原子声明

| ID | 声明 | 直接来源 | 本记录的证据状态 |
|---|---|---|---|
| RH-1 | 首个含 R795 的可核实导入文件已经 [0, 0] | `05ae6e1:data/iyli21.xml` 及同提交参数；父历史模板无 R795 | 本次源文件核实，待独立审计 |
| RH-2 | 22 个含 R795 的不同历史 XML blob 全为 [0, 0] | 上述 45 revision／23 blob 的定点路径扫描 | 本次静态历史检查，待独立审计 |
| RH-3 | 当前 metadata 选择对 R795 仅改计量，不改 bounds | `metadata_reaction_selection.json` R795；共享应用函数；两个 build.json | 本次源文件核实，待独立审计 |
| RH-4 | 最新 GPR 假设仍保留 R795 的 [0, 0] | 假设整理数据、前置条件检查及最终赋值 | 本次源文件核实，待独立审计 |
| RH-5 | 已检查范围未提供原始关闭动机 | 原始反应无备注；导入提交及脚本；限定范围文本／历史检索 | 有范围限制的未找到，不是证明所有来源都无说明 |

本记录不把自身检查计作独立审计；审计覆盖率由根任务的独立审计记录报告。

## 输入身份

| 来源 | 完整 SHA-256 |
|---|---|
| `0004d768:model.xml`（95 反应模板） | `671c33157ebf07874ea3c5b1928cfd011d127272f62bc6c830da8d58c9975ea9` |
| `05ae6e1:data/iyli21.xml`／`05ae6e1:model.xml` | `6974b7588f2a6c60ba2cde2f26e20d3aba1334d0d501572bf01cee47eda86631` |
| 当前 `data/iyali26.xml` | `5c8c199e2c5b622e97daf2b3500f763f83519fb598702a11dd153052c6a99f9d` |
| 当前 `model.xml` | `576a284ee86f0b96c802ea2e4445a862da49463022962471bfc841c469fcb5f2` |
| 当前 `model_metadata_trna.xml` | `d274bad3050e3c9220a8b6287eae847f3bf1334892284d565a6c4d96b38135a0` |
| `data/metadata_reaction_selection.json` | `d994792f374732897c839172917558ace854a347a0a47f4b760b17572530c134` |
| `scripts/gem_annotate/reaction_selection.py` | `aeec4d6855179466b825f3d9f0d0e3552d46401798d66d1af4b5cd5ad1ee372c` |
| `scripts/gem_annotate/main.py` | `f78e1e56ef088f861a9ee032c5db4c779d0b4e84a7a9465472923426c98fe67e` |
| `data/reference_build/curation/vatpase_gpr_hypothesis.json` | `b41e5cc151c575dd9f31c4f454b55f60c1c053b9e0f69915d9844f3ece200e6b` |
| `scripts/gem_annotate/patches.py` | `cd19ecc729809c9661eec423b741a61313c3470cdae4017383750383c5566b5d` |
| `model_metadata_trna_vatpase_and_hypothesis.build.json` | `c180546fbd84820f4b0499afa47d46f7be252fb4b8b066fc0e74de40c16120d5` |

现有工作区含其他任务的未提交改动。本调查读取的是上述哈希对应工作文件与指定 Git 对象；不声称 HEAD 本身重建了完整当前工作环境。

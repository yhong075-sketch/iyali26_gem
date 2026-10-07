# iYali26 GEM 静态基线评价合同

状态：`provisional`；2026-09-05。机器可读身份见 [baseline_manifest_20260905.json](baseline_manifest_20260905.json)，项目入口见 [STATE.md](../STATE.md)。本次授权仅创建文档；下文计算命令 **未执行**，计算状态为 `awaiting_compute_authorization`。文件存在、哈希核实、报告声称完成和本次计算复现是不同证据层级。

## 1. 参考对象、覆盖率与召回率

本合同固定的历史静态参照是 `ce6d752991bb-20260826T072416Z` 的 B-group／PO1f screen。实际输入 XML 为 `model.trna_biomass_group_b.xml`（`07dbdc7b…`；manifest 的 `artifacts.historical_executed_model`）；不能把 canonical（`bc2aac8f…`；`artifacts.canonical_model`）写成历史已执行模型，也不能用当前目录模型自动替换它们。完整SHA与路径只在manifest维护。日期、文件名和分支名称均不决定基线。

| 口径 | 历史交付物中的值 | 含义 |
|---|---:|---|
| 实验正例 | 1,612 | consensus-essential 正例参考，不是完整正负标签 |
| 模型内实验正例 | 322 | 正例列表与模型基因的交集 |
| 覆盖率 | 322 / 1,612 = 19.98% | 正例评价交集覆盖率 |
| 未进入当前模型评价交集的正例 | 1,290 | **1290个未进入当前模型评价交集，具体原因待分类**；原程序标记`outside_gem_scope`，不能当非必需或直接归因于某个缺失机制 |
| 10% 主阈值 | TP 67、FN 255 | 分母是322，不是1,612 |
| 召回率 | 67 / 322 = 20.81% | 已覆盖实验正例中被模型判为必需的比例 |
| WT | 1.4650760568106092 h⁻¹ | 既有静态FBA结果，本轮没有重算 |
| 模型／实际筛选基因 | 1,075 / 1,074 | 模型数含运行时质粒LEU2伪基因；该伪基因不进入KO筛选 |

历史1%、5%、10%、15%阈值分别为 TP/FN：57/265、63/259、67/255、79/243。上述结果已核对交付物，不是本次复现；6个或更少对照的检查不能更新全局召回率。

## 2. 数据来源、条件及依赖限制

- 原始正例工作簿 `42003_2023_4996_MOESM10_ESM.xlsx` 的 `Sheet1` 只有两列：A1 为 ID、B1 为 Function，A2:B1613 为1,612条记录。标准化CSV有5列；新增的来源／置信度等字段不意味着原表有逐基因三来源投票数据。
- 本地文件的SHA可核实，出版商原件的字节身份未核实；已有README对Supplementary Data 7/8的定位不一致，保持未决，不自动选定一个来源。guide QC及sample到gene结果的完整血缘同样未建立，见manifest的`assay_source_audit`与`gaps.claim_limitation`中的CL01。
- 项目保存的文献与来源审计称该consensus取Cas9、Cas12a、transposon三类screen中至少2/3来源支持；这是文献报告的定义，不是本次重建的逐基因投票。对应manifest条目为`consensus_dependency_sources`、`consensus_dependency_source_audit`；工作簿和映射依赖分别见`consensus_source_workbook`、`experimental_positive`、`locus_mapping_source`。
- 将这1,612个ID与已存在的Cas9/Cas12a结果交叉核对，模式为：**898 essential/essential；536 essential/nonessential；176 nonessential/essential；2 essential/该基因行缺失**。这是两份assay结果的对应关系，不是独立三来源一致性验证。transposon逐基因calls尚未恢复；不得由“consensus”名称补出投票、独立性或具体缺失原因。
- Cas9/Cas12a原始双sheet工作簿保留7,854／7,795行，标准化长表15,649行，来源工作簿身份见manifest的 `artifacts.assay_source_workbook`（`ba1eca8f…`）。连续fitness、assay calls和consensus正例有来源依赖，不能把同源数据的内部匹配当独立外部验证。
- 本静态历史run的 `assay_fitness_enabled=false`。该结果验证正例召回，不等于用原始guide-level acCRISPR数据完成校准；guide质量、切割效率／QC和独立外部检验的证据缺口保持明确。
- 此基准已经用于项目诊断与历史修模：manifest中的`curation_ledger`、`curated_patches`、`legacy_patch_evidence`分别指向持久案件账本、已整理模型补丁表和legacy补丁证据。它目前应作为开发／回归参考；与模型开发的独立性尚未建立，不能称未参与开发的独立测试集。
- 使用冻结的SD-Leu medium和PO1f strain profile：先加载SBML，再应用medium，再应用既有strain overlay。overlay只作用于一次性内存模型：关闭R612、对R45采用质粒互补的运行时GPR、uracil uptake设为1000；不写回任何XML。此1000是静态非限制供给条件，不是实测摄取上限；leucine uptake保持关闭，有限铁供给保持开启。
- Gene ID、mapping、培养基、实验标签及profile按manifest中的路径和SHA使用；不能删减实验CSV以伪装成子集运行，不能补写缺失标签，不能自动换表、换模型或换条件。

## 3. 历史代码实际采用的评价规则

已恢复工作树：`/private/tmp/iyali26_coq9_wp12_20260904`。历史dirty checkout的5个有记录源码SHA（`patches.py`、`sbml.py`、`strain_overlay.py`、`trna_biomass_pipeline.py`、`validate_essential_genes.py`）逐项匹配；这仅恢复已记录文件的内容身份，不能恢复整个dirty工作区、所有导入依赖或完整历史环境。manifest的 `code_identity.source_files` 区分历史已记录与本次补充固定的文件；补充SHA不能冒充历史记录。

依据该工作树 `scripts/gem_annotate/validate_essential_genes.py`：

| 行号／函数 | 原有行为，必须保留 |
|---|---|
| 520–564 `run_single_gene_deletions` | 指定solver；WT必须optimal、objective非None且0.1 ≤ WT ≤ 2.0 h⁻¹。NaN/±Inf不通过区间检查。基因排序，`processes=1`。 |
| 549–560 | KO缺status默认optimal；非optimal、growth=None或NaN均归为0；其余负growth截为0。KO/WT用该growth除以WT。没有显式拒绝正Inf。 |
| 590、618、645 | 四个阈值全部使用未舍入比值的严格 `<`；等于阈值不判essential。没有epsilon。 |
| 572–621 `build_per_gene_table` | 预测表之外的ID可被标为outside_gem_scope；不能把只跑6个KO的预测表当作全量结果。实验正例按阈值分TP/FN；未标注模型基因不是TN。 |
| 636–653 `make_summary` | 分母是实验正例与实际模型基因的交集。若KO整行缺失，NaN比较可落入FN；原函数没有完整返回行数门。因此本任务不调用该汇总函数。 |
| 153–197、698–700 | positive-only允许缺essential列时把输入行作为正例；不把未标注基因补成负例。若输入本身含负标签，逐基因函数仍有FP/TN分支，标志不会自动删除负标签。正例总报告不报accuracy、MCC、precision、specificity。 |
| 418–423 | 连续assay诊断才把ratio舍入9位；essentiality判定保留原值。本历史run没有启用连续assay。 |

原KO失败转零可能使实验正例被计为TP，不能静默改写此历史规则。已导出的1,074条历史KO记录全部optimal，ratio没有缺失或非有限值；这只能说明归一化交付物中未见显式异常行，不能排除原始None／NaN／负growth已被归零，因而不能据此证明历史67/255未受影响。未来小批重放额外检查原始求解返回值，发现异常便停止，并保留原始值和能够取得的原函数legacy结果；这不修改旧标签或旧计数。

历史软件记录为 Python 3.13.5、COBRApy 0.30.0、optlang 1.8.3、pandas 2.3.3、gurobipy 13.0.1。此运行链没有显式记录或设定完整solver容差、Threads、Seed、Method参数；**历史容差未核实**。代码中的FVA用`FLUX_EPS=1e-9`不是本静态essentiality判定的阈值容差。不得把CoQ9 WP2另一个运行的参数移植成这次历史设置。

### XML差异与等价性边界

对历史执行XML与canonical XML进行忽略格式空白的逐节点静态比较：144,217个节点匹配，只有21处属性变化：SBML根新增 `metaid="meta_"`；20个 `M_trna_biomass_residue_R*` species 的缺省 `fbc:charge` 变成显式 `0`。完整差异在manifest的 `xml_comparison`。

所检查的notes／annotation子树未见差异；根metaid是元数据变化；后20项charge是化学表示变化，不能概括为“只有无意义注释”。XML所存反应计量、bounds、GPR、objective未发现其他差异；实际加载后的LP尚未核查，未求解，不能宣布完全优化等价或科学等价。通过6个对照也只能支持这6个对照在本合同条件下的有限数值一致性。

## 4. 缺口分类与唯一下一任务

| 分类 | 当前内容 | 停止范围 |
|---|---|---|
| `blocking_current_task` | 当前文档整理无此类未解决项。计算尚待授权是task状态，单独记录。未来计算前若必需输入／源码SHA、软件版本、固定ID或历史对照记录不匹配，即形成该次计算阻塞。 | 只停止受影响任务，不阻塞已核实资料的文档整理。 |
| `claim_limitation` | 来源依赖、guide QC／逐基因transposon calls、独立性；**1290个未进入当前模型评价交集，具体原因待分类**；完整历史dirty环境、化学与优化等价性未建立。 | 限制结论，不能据此声称已独立校准、全面复现或生物学缺失。 |
| `other_workstream_dependency` | 旧inactive／isozyme案件版本绑定；CoQ9 exact-zero及14个非有限参数记为null问题；同事lipid候选activation门。 | 分别管理，不将其混入本静态小批任务或修改模型来绕过。 |

**唯一下一任务：固定6 KO的两阶段静态小批重放。** 先在历史实际执行XML运行WT+6 KO；所有检查通过，才对canonical XML运行同样WT+6 KO。原CLI没有基因子集参数，`--batch-size`只控制证据案件批次；使用原函数的`excluded_gene_ids = 全部模型ID − 固定6 ID`即可，不改源码。

固定控制集合如下；顺序和ID不得根据新结果替换。选择依据是历史表和冻结模型中的氨基酸／中心碳代谢，覆盖零TP、部分TP、接近10%阈值FN、高比值单酶和OR对照；不做CoQ9或lipid反应干预。这些geneProduct的name均为COBRAProtein编号，**本轮未核实正式基因symbol**；下表功能是模型GPR／反应注释，不宣称原生蛋白功能已实验验证。

| 固定Gene ID／模型功能 | 反应 | 历史KO/WT；10%分类 | 历史per_gene行／XML反应起始行 |
|---|---|---|---|
| YALI1B19079g；methionine adenosyltransferase，合成S-adenosylmethionine | R544 | 0；TP | 417／87272 |
| YALI1B20289g；saccharopine dehydrogenase，生成lysine | R718 | 0.061205762303123405；TP | 440／96404 |
| YALI1E15659g；methionine synthase | R65、R545 | 0.10231467743348883；FN | 1430／64284、87332 |
| YALI1A15879g；histidine合成中的histidinol DH、phosphoribosyl-AMP cyclohydrolase、phosphoribosyl-ATP pyrophosphatase | R389、R651、R652 | 0.10881714957959472；FN | 146／80839、93019、93061 |
| YALI1F12842g；pyruvate kinase | R694 | 0.958744421521669；FN | 1980／95373 |
| YALI1B09025g；ribose-5-phosphate isomerase，模型含OR规则 | R712 | 0.999999999999993；FN | 306／96074，OR位于96126 |

输入只取manifest锁定的历史XML、canonical XML、medium、strain profile、正例表、历史per_gene／summary／run manifest和恢复源码。动作只是在新加载的内存模型上调用原context loader、原WT／单KO函数和既有运行时overlay。`unittest.mock`包装器将固定6 KO按预定ID排序逐个交给原`single_gene_deletion`引擎，观察原始返回值；此预声明的逐个派发用于首个异常后立即停止，不改变原引擎、单KO求解参数或历史归零公式。可解析的部分返回表仍交回原函数归一化并保留；缺列／空表等无法归一化时明确标记legacy不可得，不伪造结果。不得调用full CLI、`make_summary`、FVA、dFBA、全量KO或重建管线。

预算为总600秒、单计算进程、两阶段合计最多14次primary FBA调用、每个阶段最多1 WT+6 KO、只尝试一次；该上限不是库内部数值迭代次数。`Threads=1`是**本次预先声明的资源设置，不是已知历史设置**；不覆写FeasibilityTol、OptimalityTol等容差。每次求解前记录Gurobi Presolve和optlang presolve；若optlang为`auto`则停止，因为该配置在非optimal时可能触发库层隐式重试，不自动改值。监督进程仅负责600秒硬停止，不执行模型计算。

验收：固定ID完整唯一；原始返回值含ids/growth/status；无缺行、重复、非optimal、None、NaN、Inf或负growth；WT通过历史有效性要求；原函数legacy输出保持可见。WT、KO growth和KO/WT分别以`abs_tol=1e-8, rel_tol=1e-6`比较；四阈值调用仍按原值严格`<0.01/0.05/0.10/0.15`且必须一致。第二阶段同时对比历史参考和刚通过的第一阶段。差异门不是新essentiality阈值。任何失败立即停止后续阶段，不调模型、培养基、GPR、阈值或对照集合。

产出只为stdout中的一个JSON对象，包含原始返回、legacy结果、软件／声明参数、比较结果和失败状态。非有限值保留为`"NaN"`、`"+Infinity"`、`"-Infinity"`字符串，None保持null，禁止JSON NaN常量。不得创建输出文件、写回状态／模型／ledger。未来明确获批运行时，可另行授权捕获stdout；本轮不运行或捕获。

## 5. 完整命令：未执行，须另行取得计算授权

解释器已核实为项目现有`.venv/bin/python`；`-B`关闭字节码写入。命令不依赖当前目录，不创建脚本。导入和solver的原始日志被送入已有`/dev/null`设备以避免污染JSON及暴露许可证信息；只记录白名单环境字段。监督进程捕获内存中的阶段检查点，在失败或超时后只输出最后已获得的事实，不伪造未完成结果。

```bash
'/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/.venv/bin/python' -B - <<'PY'
import json
import subprocess
import sys

worker = r'''
import csv
import importlib.metadata as metadata
import json
import math
import os
from pathlib import Path
import platform
import sys
from hashlib import sha256
from unittest.mock import patch

channel = os.fdopen(os.dup(1), "w", buffering=1)
sink = os.open(os.devnull, os.O_WRONLY)
os.dup2(sink, 1)
os.dup2(sink, 2)
os.close(sink)
report = {"status": "running", "phase": "preflight", "stages": [],
          "primary_solves_requested_upper_bound": 0,
          "purpose": "bounded_static_reference_check",
          "historical_full_environment_reproduction": "not_established",
          "logs_suppressed": True}

def safe(value):
    if isinstance(value, dict):
        return {str(k): safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [safe(v) for v in value]
    if isinstance(value, (set, frozenset)):
        return [safe(v) for v in sorted(value, key=str)]
    if hasattr(value, "item"):
        return safe(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return "NaN" if math.isnan(value) else ("+Infinity" if value > 0 else "-Infinity")
    return value

def emit():
    print(json.dumps(safe(report), ensure_ascii=False, allow_nan=False), file=channel, flush=True)

def require(condition, message):
    if not condition:
        report["failure_reason"] = message
        raise RuntimeError("declared_gate_failed")

def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()

def finite_nonnegative(value):
    try:
        return value is not None and math.isfinite(float(value)) and float(value) >= 0
    except (TypeError, ValueError):
        return False

emit()
try:
    manifest_path = Path("/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/docs/baseline_manifest.json")
    manifest = json.loads(manifest_path.read_text())
    report["manifest_sha256"] = digest(manifest_path)
    require(manifest["schema_version"] == "1.0", "Unexpected manifest schema")
    require(manifest["baseline_id"] == "iyali26_static_reference_provisional_20260905", "Unexpected baseline")
    artifacts = manifest["artifacts"]
    required = ("historical_executed_model", "canonical_model", "experimental_positive",
                "medium", "strain_profile", "historical_per_gene", "historical_run_manifest",
                "historical_summary")
    report["input_hashes"] = {}
    for key in required:
        item = artifacts[key]
        observed = digest(item["path"])
        report["input_hashes"][key] = observed
        require(observed == item["sha256"], "Input SHA mismatch: " + key)
    code = manifest["code_identity"]
    worktree = Path(code["recovered_worktree"]).resolve()
    require(str(worktree) == "/private/tmp/iyali26_coq9_wp12_20260904", "Unexpected recovered worktree")
    require(bool(code["source_files"]), "Missing source-file inventory")
    report["source_hashes"] = {}
    for item in code["source_files"]:
        path = Path(item["path"]).resolve()
        require(path == (worktree / item["relative_path"]).resolve(), "Source path mismatch")
        observed = digest(path)
        report["source_hashes"][item["relative_path"]] = {
            "sha256": observed, "historically_hash_recorded": item["historically_hash_recorded"]}
        require(observed == item["sha256"], "Source SHA mismatch: " + item["relative_path"])
    actual_software = {"python": platform.python_version()}
    for package in ("cobra", "optlang", "pandas", "gurobipy"):
        actual_software[package] = metadata.version(package)
    report["software"] = actual_software
    report["python_executable"] = sys.executable
    require(actual_software == code["expected_software"], "Software version mismatch")
    selected = ["YALI1B19079g", "YALI1B20289g", "YALI1E15659g",
                "YALI1A15879g", "YALI1F12842g", "YALI1B09025g"]
    task = manifest["next_task"]
    require(task["selected_gene_ids"] == selected, "Predeclared control list changed")
    require(task["max_primary_solves"] == 14 and task["wall_time_seconds"] == 600, "Budget changed")
    require(task["comparison_abs_tol"] == 1e-8 and task["comparison_rel_tol"] == 1e-6, "Comparison gate changed")
    report["selected_gene_ids"] = selected
    report["new_declared_resource_settings"] = {"Threads": 1, "processes": 1, "wall_time_seconds": 600}
    cuts = (0.01, 0.05, 0.10, 0.15)
    old_manifest = json.loads(Path(artifacts["historical_run_manifest"]["path"]).read_text())
    for old_key, key in (("model", "historical_executed_model"), ("experimental", "experimental_positive"),
                         ("medium", "medium"), ("strain_profile", "strain_profile")):
        require(old_manifest["inputs"][old_key]["sha256"] == artifacts[key]["sha256"], "Historical input binding mismatch: " + key)
    old_summary = json.loads(Path(artifacts["historical_summary"]["path"]).read_text())
    with Path(artifacts["historical_per_gene"]["path"]).open(newline="") as handle:
        old_rows = [r for r in csv.DictReader(handle, delimiter="\t") if r["gene_id"] in selected]
    require(len(old_rows) == 6 and len({r["gene_id"] for r in old_rows}) == 6, "Historical controls missing or duplicated")
    reference = {r["gene_id"]: r for r in old_rows}
    require(all(r["ko_status"] == "optimal" and r["in_model"] == "True" and
                r["experimental_essential"] == "True" and finite_nonnegative(r["ko_growth"]) and
                finite_nonnegative(r["ko_growth_ratio"]) for r in old_rows), "Invalid historical control record")
    call_columns = ("essential_at_1pct", "essential_at_5pct", "essential_at_10pct", "essential_at_15pct")
    historical_calls = {}
    for row in old_rows:
        require(all(row.get(column) in ("True", "False") for column in call_columns), "Missing or invalid historical cutoff calls")
        recorded = [row[column] == "True" for column in call_columns]
        require(recorded == [float(row["ko_growth_ratio"]) < cutoff for cutoff in cuts], "Historical recorded calls disagree with raw-ratio rule")
        historical_calls[row["gene_id"]] = recorded
    require(finite_nonnegative(old_summary["wt_growth"]), "Invalid historical WT")
    os.environ["IYALI26_RESEARCH_ROOT"] = manifest["workspace"]["research_root"]
    report["research_root"] = os.environ["IYALI26_RESEARCH_ROOT"]
    sys.path.insert(0, str(worktree))
    from scripts.gem_annotate import validate_essential_genes as validator
    from scripts.gem_annotate.essentiality_simulation_context import load_effective_simulation_context
    require(Path(validator.__file__).resolve() == worktree / "scripts/gem_annotate/validate_essential_genes.py", "Wrong module imported")
    experimental = validator.load_experimental(Path(artifacts["experimental_positive"]["path"]), positive_only=True)
    positive_ids = set(experimental.loc[experimental["essential"].eq(True), "gene_id"])
    require(set(selected) <= positive_ids, "Control absent from positive reference")
    report["phase"] = "preflight_passed"
    emit()

    def compare(stage, wt, rows, ref_wt, ref_rows, ref_calls, label):
        checks = {"wt": math.isclose(wt, float(ref_wt), abs_tol=1e-8, rel_tol=1e-6), "genes": {}}
        for row in rows:
            target = ref_rows[row["gene_id"]]
            ratio = float(row["ko_growth_ratio"])
            ref_ratio = float(target["ko_growth_ratio"])
            checks["genes"][row["gene_id"]] = {
                "growth": math.isclose(float(row["ko_growth"]), float(target["ko_growth"]), abs_tol=1e-8, rel_tol=1e-6),
                "ratio": math.isclose(ratio, ref_ratio, abs_tol=1e-8, rel_tol=1e-6),
                "calls": [ratio < c for c in cuts], "reference_calls": ref_calls[row["gene_id"]]}
        stage.setdefault("comparisons", {})[label] = checks
        emit()
        require(checks["wt"] and all(x["growth"] and x["ratio"] and x["calls"] == x["reference_calls"]
                                    for x in checks["genes"].values()), "Numerical or cutoff mismatch: " + label)

    for key in ("historical_executed_model", "canonical_model"):
        stage = {"model_artifact": key, "status": "running"}
        report["stages"].append(stage)
        report["phase"] = key
        emit()
        context = load_effective_simulation_context(
            model_path=artifacts[key]["path"], media_path=artifacts["medium"]["path"],
            strain_profile_path=artifacts["strain_profile"]["path"])
        stage["loaded_input_hashes"] = {key: context.canonical_model_sha256,
                                       "medium": context.medium_sha256,
                                       "strain_profile": context.strain_profile_sha256}
        emit()
        require(all(value == artifacts[name]["sha256"] for name, value in stage["loaded_input_hashes"].items()),
                "Loaded input SHA mismatch; no solve dispatched")
        model = context.model
        model.solver = "gurobi"
        model.solver.problem.Params.Threads = 1
        stage["context"] = context.provenance()
        stage["solver_parameters"] = {name: getattr(model.solver.problem.Params, name)
                                      for name in ("Threads", "FeasibilityTol", "OptimalityTol", "Method", "Seed", "Presolve")}
        stage["optlang_presolve"] = model.solver.configuration.presolve
        require(stage["optlang_presolve"] != "auto", "optlang auto presolve permits implicit retry; stop without changing it")
        stage["cobra_tolerance"] = model.tolerance
        all_ids = {gene.id for gene in model.genes}
        require(set(selected) <= all_ids, "Predeclared control missing from model")
        require(context.strain_overlay_enabled and context.active_medium.get("R1354") == 1000.0 and
                context.active_medium.get("R1189", 0) > 0 and context.active_medium.get("R1219", 0) == 0,
                "Historical runtime medium/strain invariant failed")
        raw_tables = []
        original_deletion = validator.single_gene_deletion
        original_optimize = model.optimize

        def observe_wt(*args, **kwargs):
            require(model.solver.configuration.presolve != "auto", "Unexpected auto presolve before WT")
            report["primary_solves_requested_upper_bound"] += 1
            require(report["primary_solves_requested_upper_bound"] <= 14, "Primary solve budget exceeded")
            emit()
            solution = original_optimize(*args, **kwargs)
            stage["raw_wt"] = {"status": solution.status, "objective_value": solution.objective_value}
            emit()
            if solution.status == "optimal" and finite_nonnegative(solution.objective_value) and 0.1 <= float(solution.objective_value) <= 2.0:
                checks = {"historical_record": math.isclose(float(solution.objective_value), float(old_summary["wt_growth"]), abs_tol=1e-8, rel_tol=1e-6)}
                if key == "canonical_model":
                    checks["first_stage_replay"] = math.isclose(float(solution.objective_value), report["stages"][0]["legacy"]["wt_growth"], abs_tol=1e-8, rel_tol=1e-6)
                stage["early_wt_comparisons"] = checks
                emit()
                if not all(checks.values()):
                    stage["legacy_unavailable_reason"] = "Stopped after WT comparison, before original function requested KO"
                    require(False, "WT comparison failed; no KO dispatched")
            return solution

        def observe_deletion(*args, **kwargs):
            require(len(args) == 1 and set(kwargs) == {"processes", "gene_list"} and
                    kwargs["processes"] == 1 and list(kwargs["gene_list"]) == sorted(selected), "KO request exceeded declared scope")
            frames = []
            stage["raw_deletion"] = {"dispatch": "predeclared single-KO calls to unchanged engine", "calls": []}
            for gene_id in sorted(selected):
                require(model.solver.configuration.presolve != "auto", "Unexpected auto presolve before KO")
                report["primary_solves_requested_upper_bound"] += 1
                require(report["primary_solves_requested_upper_bound"] <= 14, "Primary solve budget exceeded")
                stage["requested_ko"] = gene_id
                emit()
                raw = original_deletion(args[0], gene_list=[gene_id], processes=1)
                stage["raw_deletion"]["calls"].append({"requested_gene_id": gene_id,
                    "columns": list(raw.columns), "rows": [
                        {"index": index, "values": row.to_dict()} for index, row in raw.iterrows()]})
                emit()
                if not {"ids", "growth", "status"} <= set(raw.columns) or raw.empty:
                    stage["legacy_unavailable_reason"] = "Missing raw columns or empty return; do not fabricate normalized rows"
                    require(False, "Raw KO malformed or missing; no further KO dispatched")
                frames.append(raw.copy(deep=True))
                valid_id = (len(raw) == 1 and isinstance(raw.iloc[0]["ids"], (set, frozenset)) and raw.iloc[0]["ids"] == {gene_id})
                if not valid_id or not all(raw["status"].eq("optimal")) or not all(finite_nonnegative(v) for v in raw["growth"]):
                    stage["early_raw_gate_failure"] = "Raw KO wrong/duplicate ID, nonoptimal, None/NaN/Inf or negative growth"
                    break
                growth = float(raw.iloc[0]["growth"])
                ratio = growth / float(stage["raw_wt"]["objective_value"])
                observed_calls = [ratio < cutoff for cutoff in cuts]
                references = [("historical_record", reference[gene_id], historical_calls[gene_id])]
                if key == "canonical_model":
                    first = report["stages"][0]
                    first_row = next(r for r in first["legacy"]["predictions"] if r["gene_id"] == gene_id)
                    first_calls = first["comparisons"]["historical_record"]["genes"][gene_id]["calls"]
                    references.append(("first_stage_replay", first_row, first_calls))
                early_checks = {}
                for label, target, target_calls in references:
                    early_checks[label] = {"growth": math.isclose(growth, float(target["ko_growth"]), abs_tol=1e-8, rel_tol=1e-6),
                                           "ratio": math.isclose(ratio, float(target["ko_growth_ratio"]), abs_tol=1e-8, rel_tol=1e-6),
                                           "calls": observed_calls, "reference_calls": target_calls}
                stage.setdefault("early_ko_comparisons", {})[gene_id] = early_checks
                emit()
                if not all(c["growth"] and c["ratio"] and c["calls"] == c["reference_calls"] for c in early_checks.values()):
                    stage["early_comparison_gate_failure"] = "KO numerical or cutoff mismatch; do not dispatch later KO"
                    break
            combined = validator.pd.concat(frames, ignore_index=True)
            raw_tables.append(combined)
            return combined

        with patch.object(model, "optimize", observe_wt), patch.object(validator, "single_gene_deletion", observe_deletion):
            predictions, wt = validator.run_single_gene_deletions(model, "gurobi", excluded_gene_ids=all_ids - set(selected))
        rows = predictions.to_dict(orient="records")
        stage["legacy"] = {"wt_growth": wt, "predictions": rows}
        emit()
        require("early_raw_gate_failure" not in stage, "Raw KO anomaly: partial legacy retained, later KO not dispatched")
        require("early_comparison_gate_failure" not in stage, "KO comparison failed: partial legacy retained, later KO not dispatched")
        require(len(raw_tables) == 1, "Unexpected deletion call count")
        raw = raw_tables[0]
        require({"ids", "growth", "status"} <= set(raw.columns), "Raw KO columns missing")
        raw_ids = [validator._gene_id_from_deletion_row(index, row) for index, row in raw.iterrows()]
        require(len(raw) == 6 and len(set(raw_ids)) == 6 and set(raw_ids) == set(selected), "Raw KO rows missing, duplicated or unexpected")
        require(all(isinstance(ids, (set, frozenset)) and len(ids) == 1 for ids in raw["ids"]), "Raw KO IDs are not singleton deletions")
        require(all(status == "optimal" for status in raw["status"]), "Raw KO nonoptimal: retain legacy result and stop")
        require(all(finite_nonnegative(value) for value in raw["growth"]), "Raw KO None/NaN/Inf/negative: retain legacy result and stop")
        require(len(rows) == 6 and {r["gene_id"] for r in rows} == set(selected) and
                all(finite_nonnegative(r["ko_growth"]) and finite_nonnegative(r["ko_growth_ratio"]) for r in rows), "Invalid legacy result shape or numbers")
        compare(stage, wt, rows, old_summary["wt_growth"], reference, historical_calls, "historical_record")
        if key == "canonical_model":
            first = report["stages"][0]["legacy"]
            first_calls = {gene_id: report["stages"][0]["comparisons"]["historical_record"]["genes"][gene_id]["calls"] for gene_id in selected}
            compare(stage, wt, rows, first["wt_growth"], {r["gene_id"]: r for r in first["predictions"]}, first_calls, "first_stage_replay")
        stage["status"] = "passed"
        emit()
    report["status"] = "passed"
    report["phase"] = "complete"
    report["claim_limit"] = "Only these six controls and WT; no global recall, calibration or complete model equivalence claim"
except Exception as exc:
    report["status"] = "stopped"
    report.setdefault("failure_reason", "Execution error; no automatic retry or substitution")
    report["exception_type"] = type(exc).__name__
    if report["stages"] and report["stages"][-1]["status"] == "running":
        report["stages"][-1]["status"] = "stopped"
emit()
'''

process = subprocess.Popen([sys.executable, "-B", "-c", worker], stdout=subprocess.PIPE,
                           stderr=subprocess.DEVNULL, text=True)
timed_out = False
try:
    output, _ = process.communicate(timeout=600)
except subprocess.TimeoutExpired:
    timed_out = True
    process.kill()
    output, _ = process.communicate()
checkpoints = []
for line in output.splitlines():
    try:
        checkpoints.append(json.loads(line))
    except json.JSONDecodeError:
        pass
result = checkpoints[-1] if checkpoints else {"status": "stopped", "stages": [], "failure_reason": "No complete worker checkpoint"}
if timed_out or process.returncode != 0 or result.get("status") == "running":
    result["status"] = "stopped"
    result["failure_reason"] = "600-second wall-clock limit" if timed_out else "Worker did not complete normally"
    for stage in result.get("stages", []):
        if stage.get("status") == "running":
            stage["status"] = "stopped"
result["worker_returncode"] = process.returncode
result["timed_out"] = timed_out
print(json.dumps(result, ensure_ascii=False, allow_nan=False))
sys.exit(0 if result["status"] == "passed" else 1)
PY
```

本命令文本的语法检查不构成模型检查。只有未来取得计算授权并实际执行、保留stdout结果后，才可为该次限定检查记录“已执行”；任何当前文档不得预先写成已通过。

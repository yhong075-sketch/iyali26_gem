# 本轮实际执行的主要命令

以下命令均已在项目根目录执行成功；不是待实现的参数示例。构建/验证故意拒绝覆盖已有输出。记录中的路径现已有产物；再次执行同一路径会拒绝覆盖，不应先删除旧结果。配置的培养/菌株身份固定，模型输入必须由命令显式提供。两个 fresh_validation 共用本轮独立预算；audit_saved_run 只读历史预算及保存记录，不消耗优化次数。

## 重建开放候选并验证实际新文件

```bash
TMPDIR="$PWD/artifacts/dipeptide_chemistry_routes_20260924/tmp" .venv/bin/python -B scripts/build_vacuole_candidate.py --source artifacts/atp_candidate_repair_20260924/candidates/E5.xml --source-sha256 43042d16f874a91d61f12f9f2b65838cbc264821d98bcf6e108300bf787007be --output artifacts/dipeptide_chemistry_routes_20260924/E5_vacuole_open_rebuilt.xml --enabled

TMPDIR="$PWD/artifacts/dipeptide_chemistry_routes_20260924/tmp" .venv/bin/python -B -m scripts.validate_vacuole_supply --mode fresh_validation --config artifacts/dipeptide_chemistry_routes_20260924/config.json --baseline-model artifacts/atp_candidate_repair_20260924/candidates/E5.xml --baseline-sha256 43042d16f874a91d61f12f9f2b65838cbc264821d98bcf6e108300bf787007be --candidate-model artifacts/dipeptide_chemistry_routes_20260924/E5_vacuole_open_rebuilt.xml --build-manifest artifacts/dipeptide_chemistry_routes_20260924/E5_vacuole_open_rebuilt.build.json --expected-diff four_connections --output artifacts/dipeptide_chemistry_routes_20260924/fresh_rebuilt --budget artifacts/dipeptide_chemistry_routes_20260924/solver_budget.json
```

## 化学候选构建和新优化

```bash
TMPDIR="$PWD/artifacts/dipeptide_chemistry_routes_20260924/tmp" .venv/bin/python -B scripts/build_dipeptide_chemistry.py --source artifacts/dipeptide_chemistry_routes_20260924/E5_vacuole_open_rebuilt.xml --source-sha256 e4806a1a7fcd69af1b5ec7c2e014890828548dd1e77630db2b3a2a96dbfc9a39 --patch data/dipeptide_chemistry_patch.json --output artifacts/dipeptide_chemistry_routes_20260924/E5_vacuole_open_chemistry.xml --enabled

TMPDIR="$PWD/artifacts/dipeptide_chemistry_routes_20260924/tmp" .venv/bin/python -B -m scripts.validate_vacuole_supply --mode fresh_validation --config artifacts/dipeptide_chemistry_routes_20260924/config.json --baseline-model artifacts/dipeptide_chemistry_routes_20260924/E5_vacuole_open_rebuilt.xml --baseline-sha256 e4806a1a7fcd69af1b5ec7c2e014890828548dd1e77630db2b3a2a96dbfc9a39 --candidate-model artifacts/dipeptide_chemistry_routes_20260924/E5_vacuole_open_chemistry.xml --build-manifest artifacts/dipeptide_chemistry_routes_20260924/E5_vacuole_open_chemistry.build.json --expected-diff chemistry_only --output artifacts/dipeptide_chemistry_routes_20260924/fresh_chemistry --budget artifacts/dipeptide_chemistry_routes_20260924/solver_budget.json > artifacts/dipeptide_chemistry_routes_20260924/fresh_chemistry.log 2>&1
```

## 保存记录核查和无优化检查

```bash
TMPDIR="$PWD/artifacts/dipeptide_chemistry_routes_20260924/tmp" .venv/bin/python -B -m scripts.validate_vacuole_supply --mode audit_saved_run --config artifacts/vacuole_open_supply_20260924/config.json --baseline-model artifacts/atp_candidate_repair_20260924/candidates/E5.xml --baseline-sha256 43042d16f874a91d61f12f9f2b65838cbc264821d98bcf6e108300bf787007be --candidate-model artifacts/vacuole_open_supply_20260924/E5_vacuole_open.xml --build-manifest artifacts/vacuole_open_supply_20260924/E5_vacuole_open.build.json --saved-run artifacts/vacuole_open_supply_20260924/run --source-archive artifacts/dipeptide_chemistry_routes_20260924/previous_source_snapshot --output artifacts/dipeptide_chemistry_routes_20260924/audit_previous_saved_run --budget artifacts/vacuole_open_supply_20260924/budget.json

TMPDIR="$PWD/artifacts/dipeptide_chemistry_routes_20260924/tmp" .venv/bin/python -B -m unittest tests.test_dipeptide_chemistry tests.test_vacuole_explicit_inputs > artifacts/dipeptide_chemistry_routes_20260924/engineering_tests.log 2>&1

TMPDIR="$PWD/artifacts/dipeptide_chemistry_routes_20260924/tmp" .venv/bin/python -B artifacts/dipeptide_chemistry_routes_20260924/audit_chemical_candidate.py
```

另已在 `execution_limits(no_solve=True, allow_network=False)` 内执行 `tests.test_reaction_selection` 的7项既有测试，结果见 `metadata_regression_tests.log`。结构解析、序列比较、路线参数检查和独立审计的执行记录分别随对应子目录保存。没有重新执行历史114次矩阵，没有新增FVA。

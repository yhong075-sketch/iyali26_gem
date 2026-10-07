# 本轮已执行的分析

工作目录为项目根目录。以下分析命令均已实际执行；不调用代谢优化、不修改模型。原始公共来源的检索由同目录记录脚本及retrieval logs保存；远程BLAST只提交过一个RID，重新执行提交脚本会拒绝重复提交。

```bash
.venv/bin/python -B artifacts/glypro_function_localization_20260924/sequence/analyze_target.py
.venv/bin/python artifacts/glypro_function_localization_20260924/localization/sequence_localization_heuristic.py
```

前者运行3个全局仿射比对及既有AlphaFold模型与实验参考5M4G的序列引导Kabsch叠合，包含最小得分与重建自检；后者仅输出描述性KD19及有限PTS模式，不是专门定位预测器。软件版本、评分矩阵和原始输入SHA在analysis.json和来源清单中。

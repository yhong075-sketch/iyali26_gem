# Chalmers 苯乙酸分支 · Escher

打开 `chalmers-paa-escher.html`，在顶部切换 W29 和酿酒酵母。无需联网；请将 HTML 与同目录的 escher.min.js 一起保留。

可缩放、拖动画布、点击反应查看方程/保存边界/GPR与证据等级，或开启“编辑布局”拖动节点与标签。切换视图会在当前页面保留布局修改；下载地图 JSON 或导出 SVG 保存修改。

- `w29-map.json`、`sc-map.json`：可导入 Escher 的地图。
- `w29-escher.svg/png`、`sc-escher.svg/png`：初始布局的矢量图和图片。
- 红色：完整守恒行证明的稳态阻断；青色：原模型已有运输/交换；灰色：其他反应。没有叠加模拟或实测通量。

源模型为2026-09-14已审计的Chalmers iYali v4.1.2与Yeast-GEM v9.1.1固定快照；本次没有查询后续版本。上游只显示0851转氨支路，2-苯乙醇后续去路未展开；PAM/PAA与苯乙醛三个区室的完整反应行已核对。W29图没有虚构PAM节点或缺失反应。基因赋值是源模型GPR，不能当作底物活性验证。

来源与完整SHA见 provenance.json / evidence.json；独立核对见 AUDIT.md，浏览器交互检查见 render-check.json。未运行FBA、未修改任何原模型、培养或GPR。Escher 1.8.1许可证随附。

模型来源：[Chalmers W29 iYali](https://github.com/SysBioChalmers/Yarrowia_lipolytica_W29-GEM)；[Chalmers Yeast-GEM](https://github.com/SysBioChalmers/yeast-GEM)。固定版本身份保存在 provenance.json 与 evidence.json，图件整理日期2026-09-16。

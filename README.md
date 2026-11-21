# TCM Database —— 草药-成分-靶点-疾病 四层数据库网站

独立于网络药理学分析平台的数据检索网站，数据取自 **TCMSP** 官方数据库全库：

| 层级 | 规模 | 说明 |
| --- | --- | --- |
| 草药 | 502 味 | TCMSP 全库中药（含药类分类 / 拼音 / 拉丁名） |
| 成分 | 2585 个 | OB ≥ 30% 且 DL ≥ 0.18 筛选，含 ADME 参数 |
| 靶点 | 572 个 | UniProt 标准化基因符号，附 SVM / RF 预测分 |
| 疾病 | 445 种 | TCMSP 靶点-疾病关联（DrugBank 验证子集） |

## 功能

- **中药 / 成分 / 靶点 / 疾病** 四个检索页：搜索 + 筛选 + 详情卡；
- **靶点预测页**：TCMSP SVM / RF 双模型预测分检索（约 1.7 万条成分-靶点预测对），
  支持预测分阈值过滤、实验验证筛选与综合分排行图；
- **蛋白三维结构**：全部靶点在线查看 AlphaFold 预测模型（3Dmol.js 可旋转缩放），
  结构文件按需下载；UniProt 未覆盖的靶点自动按蛋白名在线解析登录号；
- **四层联动**：疾病反查靶点、成分、来源中药；成分看靶点与相关疾病；
- **四层网络图**：按中药或按疾病构建可交互网络（含疾病层开关）；
- 成分详情含 PubChem 结构图与结构信息。

## 本地运行

```bash
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
streamlit run app.py        # http://localhost:8502
```

## 部署

推送到 GitHub 后，在 <https://share.streamlit.io> 新建应用：
仓库选本仓库、分支 `main`、主文件 `app.py` 即可（数据文件已入库，无需额外配置）。

## 数据说明

- `data/network.json`：四层网络（由原项目 `data_prep.py` 从 `tcmsp_full.json` 构建）；
- `data/target_info.json`：UniProt 蛋白注释；
- `data/compound_structures.json`：PubChem 结构信息。

> 数据库预测结果仅供科研检索，不构成临床建议。

"""TCM Database —— 草药-成分-靶点-疾病 四层数据库网站（独立项目）。

数据来源：TCMSP 官方数据库全库抓取（502 味中药），
  成分经 OB≥30% / DL≥0.18 筛选，靶点经 UniProt 标准化为基因符号，
  疾病为 TCMSP 靶点-疾病关联（DrugBank 验证子集）。
本网站为纯数据浏览与检索平台，运行 network.json 单一数据文件。

启动：streamlit run app.py
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Set

import pandas as pd
import requests
import streamlit as st
import streamlit.components.v1 as components
from st_link_analysis import st_link_analysis
from st_link_analysis.component.styles import EdgeStyle, NodeStyle

HERE = Path(__file__).parent
NETWORK_FILE = HERE / "data" / "network.json"
TARGET_INFO_FILE = HERE / "data" / "target_info.json"
STRUCTURES_FILE = HERE / "data" / "compound_structures.json"
PDB_CACHE = HERE / "data" / "pdb_cache"   # AlphaFold PDB 本地缓存（不入 git）

st.set_page_config(page_title="TCM Database · 草药-成分-靶点-疾病",
                   page_icon="🌿", layout="wide")

# ------------------------------------------------------------------ 数据层
@st.cache_resource(show_spinner="正在加载 TCMSP 全库数据…")
def load_db() -> dict:
    """加载四层网络并建立双向索引。"""
    net = json.loads(NETWORK_FILE.read_text(encoding="utf-8"))
    node_by_id: Dict[str, dict] = {}
    herbs: Dict[str, dict] = {}          # 药名 -> node data
    compounds: Dict[str, dict] = {}      # MOL_ID -> node data
    targets: Dict[str, dict] = {}        # 基因符号 -> node data
    diseases: Dict[str, dict] = {}       # 疾病名 -> node data
    for n in net["nodes"]:
        d = n["data"]
        node_by_id[d["id"]] = d
        name = d.get("name", "")
        if d["label"] == "HERB":
            herbs[name] = d
        elif d["label"] == "COMPOUND":
            compounds[d["id"]] = d
        elif d["label"] == "TARGET":
            targets[name] = d
        elif d["label"] == "DISEASE":
            diseases[name] = d

    herb_compounds: Dict[str, List[str]] = {}    # 药名 -> [MOL_ID]
    compound_targets: Dict[str, List[str]] = {}  # MOL_ID -> [基因符号]
    target_compounds: Dict[str, List[str]] = {}  # 基因符号 -> [MOL_ID]
    target_diseases: Dict[str, List[str]] = {}   # 基因符号 -> [疾病名]
    disease_targets: Dict[str, List[str]] = {}   # 疾病名 -> [基因符号]
    edge_scores: Dict[tuple, dict] = {}          # (MOL_ID, 基因) -> svm/rf/validated
    for e in net["edges"]:
        d = e["data"]
        label = d.get("label")
        if label == "CONTAINS":
            herb_compounds.setdefault(
                str(d["source"]).replace("HERB::", ""), []).append(d["target"])
        elif label == "TARGETS":
            comp = str(d["source"])
            tgt = str(d["target"]).replace("TGT::", "")
            compound_targets.setdefault(comp, []).append(tgt)
            target_compounds.setdefault(tgt, []).append(comp)
            if d.get("svm") is not None:
                edge_scores[(comp, tgt)] = {
                    "svm": d.get("svm"), "rf": d.get("rf"),
                    "validated": d.get("validated", False)}
        elif label == "INDICATES":
            tgt = str(d["source"]).replace("TGT::", "")
            dis = str(d["target"]).replace("DIS::", "")
            target_diseases.setdefault(tgt, []).append(dis)
            disease_targets.setdefault(dis, []).append(tgt)

    target_info = {}
    if TARGET_INFO_FILE.exists():
        target_info = json.loads(TARGET_INFO_FILE.read_text(encoding="utf-8"))
    structures = {}
    if STRUCTURES_FILE.exists():
        structures = json.loads(STRUCTURES_FILE.read_text(encoding="utf-8"))
    return {
        "node_by_id": node_by_id, "herbs": herbs, "compounds": compounds,
        "targets": targets, "diseases": diseases,
        "herb_compounds": {k: list(dict.fromkeys(v)) for k, v in herb_compounds.items()},
        "compound_targets": {k: list(dict.fromkeys(v)) for k, v in compound_targets.items()},
        "target_compounds": {k: list(dict.fromkeys(v)) for k, v in target_compounds.items()},
        "target_diseases": {k: sorted(set(v)) for k, v in target_diseases.items()},
        "disease_targets": {k: list(dict.fromkeys(v)) for k, v in disease_targets.items()},
        "edge_scores": edge_scores,
        "target_info": target_info, "structures": structures,
        "meta": net.get("meta", {}),
    }


db = load_db()


# ------------------------------------------------- 蛋白三维结构（AlphaFold）
@st.cache_data(show_spinner="正在从 AlphaFold 获取蛋白结构…")
def load_alphafold_pdb(accession: str) -> str:
    """按 UniProt 登录号取 AlphaFold 预测的 PDB 文本（带本地缓存）。"""
    PDB_CACHE.mkdir(exist_ok=True)
    cache_file = PDB_CACHE / f"AF-{accession}-F1.pdb"
    if cache_file.exists():
        text = cache_file.read_text(encoding="utf-8", errors="ignore")
        if "ATOM" in text:
            return text
    try:
        resp = requests.get(
            f"https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v6.pdb",
            timeout=45)
    except Exception:
        return ""
    if resp.status_code == 200 and "ATOM" in resp.text:
        try:
            cache_file.write_text(resp.text, encoding="utf-8")
        except Exception:
            pass  # 云端只读文件系统时仅失去磁盘缓存
        return resp.text
    return ""


@st.cache_data(show_spinner=False)
def resolve_accession(target_name: str, protein_hint: str = "") -> str:
    """靶点名 → UniProt 登录号（target_info 未覆盖时按蛋白名全文检索）。

    TCMSP 蛋白命名与 UniProt 官方名常有差异（如
    “5-hydroxytryptamine 1A receptor” vs “5-hydroxytryptamine receptor 1A”），
    用全文检索比精确字段匹配命中率高。
    """
    if protein_hint:
        return protein_hint
    try:
        resp = requests.get(
            "https://rest.uniprot.org/uniprotkb/search",
            params={"query": f"{target_name} AND organism_id:9606 AND reviewed:true",
                    "size": 1, "fields": "accession"},
            timeout=20)
        if resp.status_code == 200:
            results = resp.json().get("results", [])
            if results:
                return results[0]["primaryAccession"]
    except Exception:
        pass
    return ""


def render_protein_3d(pdb_text: str, height: int = 430) -> None:
    """用 3Dmol.js 渲染蛋白三维结构（PDB 内容内联，浏览器无需额外联网）。"""
    html = f"""
<div id="protein-viewer-{height}" style="height:{height}px; width:100%; background:#ffffff;"></div>
<script src="https://3dmol.org/build/3Dmol-min.js"></script>
<script>
  const viewer = $3Dmol.createViewer(document.getElementById("protein-viewer-{height}"),
                                     {{backgroundColor: "white"}});
  viewer.addModel({json.dumps(pdb_text)}, "pdb");
  viewer.setStyle({{}}, {{cartoon: {{color: "spectrum"}}}});
  viewer.zoomTo();
  viewer.render();
</script>
"""
    components.html(html, height=height + 25)


def target_3d_block(picked: str, height: int = 430) -> None:
    """靶点三维结构展示块：解析登录号 → 取 PDB → 渲染。"""
    info = db["target_info"].get(picked.upper(), {})
    acc = info.get("accession") or resolve_accession(picked, info.get("accession", ""))
    if not acc:
        st.info("未能确定该靶点的 UniProt 登录号（多为非人源蛋白或自定义名称），"
                "暂无可靠三维结构。")
        return
    with st.spinner(f"正在获取 {picked}（{acc}）的 AlphaFold 结构…"):
        pdb_text = load_alphafold_pdb(acc)
    if pdb_text:
        st.markdown(f"**AlphaFold 预测模型**（UniProt：`{acc}`），可旋转缩放：")
        render_protein_3d(pdb_text, height=height)
        st.download_button("下载 PDB 结构文件", data=pdb_text,
                           file_name=f"AF-{acc}-F1-model_v6.pdb",
                           mime="chemical/x-pdb", key=f"dl_{acc}_{height}")
    else:
        st.info("AlphaFold DB 暂无该蛋白的预测结构。")

# 侧边栏
with st.sidebar:
    st.title("🌿 TCM Database")
    st.caption("**草药 → 成分 → 靶点 → 疾病** 四层数据库")
    m = db["meta"]
    st.metric("收录中药", len(db["herbs"]))
    st.metric("活性成分", len(db["compounds"]))
    st.metric("预测靶点", len(db["targets"]))
    st.metric("关联疾病", len(db["diseases"]))
    st.caption(f"OB ≥ {m.get('ob_min', 30):g}% · DL ≥ {m.get('dl_min', 0.18):g} 筛选\n\n"
               "数据来源：TCMSP（成分 / 靶点 / SVM·RF 预测分）+ UniProt（基因符号与蛋白注释）"
               " + PubChem（化学结构）。结果为数据库预测，不构成临床建议。")

tab_home, tab_herb, tab_comp, tab_tgt, tab_pred, tab_dis, tab_net = st.tabs(
    ["📊 总览", "🌿 中药", "🧪 成分", "🎯 靶点", "🔬 靶点预测（含 3D 结构）",
     "🩺 疾病", "🕸️ 四层网络图"])


def _compound_row(mid: str) -> dict:
    d = db["compounds"].get(mid, {})
    return {
        "成分": d.get("name", mid), "编号": mid,
        "OB(%)": d.get("ob"), "DL": d.get("dl"),
        "分子量": d.get("mw", ""), "靶点数": d.get("n_targets", len(
            db["compound_targets"].get(mid, []))),
        "来源中药数": d.get("n_herbs", ""),
    }


def _disease_chips(names: List[str], limit: int = 60) -> str:
    shown = names[:limit]
    more = f"　…等 {len(names)} 种" if len(names) > limit else ""
    return "、".join(shown) + more


# ------------------------------------------------------------------ 总览
with tab_home:
    st.header("TCMSP 全库 · 草药-成分-靶点-疾病四层数据库")
    st.markdown(
        "本网站整合 **TCMSP** 数据库全部 **502 味中药**的四层数据：\n\n"
        "- **草药层**：502 味中药（含药类分类、拼音与拉丁名）；\n"
        "- **成分层**：口服生物利用度 OB ≥ 30% 且类药性 DL ≥ 0.18 的活性成分"
        f"（{len(db['compounds'])} 个，含分子量、AlogP、Caco-2、血脑屏障等 ADME 参数）；\n"
        "- **靶点层**：成分预测作用靶点，已经 UniProt 标准化为基因符号"
        f"（{len(db['targets'])} 个，附 SVM / RF 预测分）；\n"
        "- **疾病层**：靶点-疾病关联（DrugBank 验证子集，"
        f"{len(db['diseases'])} 种疾病）。")
    st.divider()
    c1, c2 = st.columns(2)
    with c1:
        st.subheader("关联疾病最多的靶点 Top 20")
        top_t = sorted(db["target_diseases"], key=lambda t: -len(db["target_diseases"][t]))[:20]
        st.bar_chart(pd.DataFrame(
            {"关联疾病数": {t: len(db["target_diseases"][t]) for t in top_t}}))
    with c2:
        st.subheader("活性成分最多的中药 Top 20")
        top_h = sorted(db["herb_compounds"], key=lambda h: -len(db["herb_compounds"][h]))[:20]
        st.bar_chart(pd.DataFrame(
            {"活性成分数": {h: len(db["herb_compounds"][h]) for h in top_h}}))
    st.divider()
    st.subheader("疾病-靶点关联数排行 Top 20")
    top_d = sorted(db["disease_targets"], key=lambda d: -len(db["disease_targets"][d]))[:20]
    st.bar_chart(pd.DataFrame(
        {"关联靶点数": {d: len(db["disease_targets"][d]) for d in top_d}}))
    st.caption("以上为数据库预测结果统计，仅供科研检索参考，不构成临床结论。")

# ------------------------------------------------------------------ 中药
with tab_herb:
    st.header("中药浏览")
    c1, c2 = st.columns([3, 1])
    kw = c1.text_input("搜索中药（名称 / 拼音 / 拉丁名）",
                       placeholder="例如：甘草 / gancao / Glycyrrhiza")
    cats = sorted({v.get("category", "") for v in db["herbs"].values()
                   if v.get("category")})
    cat = c2.selectbox("按药类筛选", ["全部"] + cats)
    view = {h: d for h, d in db["herbs"].items()
            if cat == "全部" or d.get("category") == cat}
    if kw.strip():
        needle = kw.strip().lower()
        view = {h: d for h, d in view.items()
                if needle in h.lower()
                or needle in str(d.get("pinyin", "")).lower()
                or needle in str(d.get("latin", "")).lower()}
    st.caption(f"共 **{len(view)}** 味中药")
    if view:
        rows = [{"中药": h, "拼音": d.get("pinyin", ""), "拉丁名": d.get("latin", ""),
                 "药类": d.get("category", ""), "活性成分数": len(
                     db["herb_compounds"].get(h, []))}
                for h, d in sorted(view.items())]
        st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch",
                     height=min(28 * (len(rows) + 1), 460))
        picked = st.selectbox("查看中药详情", sorted(view))
        mids = db["herb_compounds"].get(picked, [])
        tgt_set = {t for m in mids for t in db["compound_targets"].get(m, [])}
        dis_set = {d2 for t in tgt_set for d2 in db["target_diseases"].get(t, [])}
        h1, h2, h3 = st.columns(3)
        h1.metric("活性成分", len(mids))
        h2.metric("作用靶点", len(tgt_set))
        h3.metric("关联疾病", len(dis_set))
        d = view[picked]
        st.markdown(
            f"**拼音**：{d.get('pinyin') or '—'}　**拉丁名**：{d.get('latin') or '—'}　"
            f"**药类**：{d.get('category') or '—'}")
        with st.expander(f"活性成分（{len(mids)} 个，按靶点数排序）", expanded=True):
            if mids:
                rows = [_compound_row(m) for m in mids]
                st.dataframe(pd.DataFrame(rows).sort_values(
                    "靶点数", ascending=False), hide_index=True, width="stretch",
                    height=min(28 * (len(rows) + 1), 460))
        with st.expander(f"关联疾病（{len(dis_set)} 种）"):
            st.write(_disease_chips(sorted(dis_set)) if dis_set else "—")

# ------------------------------------------------------------------ 成分
with tab_comp:
    st.header("活性成分检索")
    c1, c2, c3 = st.columns(3)
    kw = c1.text_input("搜索成分（名称 / 编号）",
                       placeholder="例如：quercetin / MOL000098 / 黄芩苷")
    ob_min = c2.slider("OB(%) 下限", 0.0, 100.0, 0.0, 5.0)
    dl_min = c3.slider("DL 下限", 0.0, 1.0, 0.0, 0.05)
    view = {mid: d for mid, d in db["compounds"].items()
            if (d.get("ob", 0) or 0) >= ob_min and (d.get("dl", 0) or 0) >= dl_min}
    if kw.strip():
        needle = kw.strip().lower()
        view = {mid: d for mid, d in view.items()
                if needle in str(d.get("name", "")).lower()
                or needle in str(d.get("en", "")).lower()
                or needle in mid.lower()}
    st.caption(f"共 **{len(view)}** 个成分")
    if view:
        rows = [_compound_row(mid) for mid in sorted(view)]
        st.dataframe(pd.DataFrame(rows).sort_values(
            ["靶点数", "OB(%)"], ascending=False), hide_index=True,
            width="stretch", height=min(28 * (len(rows) + 1), 460))
        picked = st.selectbox("查看成分详情", sorted(view))
        d = view[picked]
        struct = db["structures"].get(picked, {})
        c1, c2 = st.columns([1, 2])
        with c1:
            if struct.get("img_url"):
                st.image(struct["img_url"], caption=d.get("name", picked))
            else:
                st.info("PubChem 暂无该成分的结构图。")
        with c2:
            st.markdown(f"### {d.get('name', picked)}")
            st.markdown(f"**编号**：`{picked}`　**英文名**：{d.get('en') or '—'}")
            m1, m2, m3, m4 = st.columns(4)
            m1.metric("OB(%)", d.get("ob", "—"))
            m2.metric("DL", d.get("dl", "—"))
            m3.metric("分子量", d.get("mw") or struct.get("mw") or "—")
            m4.metric("AlogP", d.get("alogp", "—"))
            st.markdown(
                f"**Caco-2**：{d.get('caco2') or '—'}　**BBB**：{d.get('bbb') or '—'}　"
                f"**H键供体**：{d.get('hdon') or '—'}")
            if struct:
                st.markdown(f"**分子式**：{struct.get('formula') or '—'}　"
                            f"**SMILES**：`{struct.get('smiles') or '—'}`")
                if struct.get("cid"):
                    st.markdown(f"[PubChem 页面（CID {struct['cid']}）]"
                                f"(https://pubchem.ncbi.nlm.nih.gov/compound/{struct['cid']})")
        herbs_of = db["node_by_id"].get(picked, {}).get("herbs") or [
            h for h, ms in db["herb_compounds"].items() if picked in ms]
        st.markdown(f"**来源中药（{len(herbs_of)} 味）**：{'、'.join(herbs_of) or '—'}")
        tids = db["compound_targets"].get(picked, [])
        scores = db["edge_scores"]
        with st.expander(f"预测靶点（{len(tids)} 个，含 SVM / RF 预测分）", expanded=True):
            if tids:
                rows = []
                for t in tids:
                    s = scores.get((picked, t), {})
                    info = db["target_info"].get(t.upper(), {})
                    rows.append({
                        "靶点": t, "蛋白名": info.get("protein", ""),
                        "SVM": s.get("svm", ""), "RF": s.get("rf", ""),
                        "实验验证": "✔" if s.get("validated") else "",
                        "关联疾病数": len(db["target_diseases"].get(t, [])),
                    })
                st.dataframe(pd.DataFrame(rows), hide_index=True,
                             width="stretch", height=min(28 * (len(rows) + 1), 460))
        dis_of = sorted({d2 for t in tids for d2 in db["target_diseases"].get(t, [])})
        with st.expander(f"经靶点关联的疾病（{len(dis_of)} 种）"):
            st.write(_disease_chips(dis_of) if dis_of else "—")

# ------------------------------------------------------------------ 靶点
with tab_tgt:
    st.header("靶点检索")
    kw = st.text_input("搜索靶点（基因符号 / 蛋白名）",
                       placeholder="例如：TNF / AKT1 / interleukin-6")
    view = db["targets"]
    if kw.strip():
        needle = kw.strip().lower()
        view = {t: d for t, d in view.items()
                if needle in t.lower()
                or needle in str(db["target_info"].get(t.upper(), {})
                                 .get("protein", "")).lower()}
    st.caption(f"共 **{len(view)}** 个靶点")
    if view:
        rows = [{"靶点": t, "蛋白名": db["target_info"].get(t.upper(), {}).get("protein", ""),
                 "关联成分数": len(db["target_compounds"].get(t, [])),
                 "关联疾病数": len(db["target_diseases"].get(t, []))}
                for t in sorted(view)]
        st.dataframe(pd.DataFrame(rows).sort_values(
            "关联成分数", ascending=False), hide_index=True, width="stretch",
            height=min(28 * (len(rows) + 1), 460))
        picked = st.selectbox("查看靶点详情", sorted(view))
        info = db["target_info"].get(picked.upper())
        if info:
            st.markdown(f"### {picked}　{info.get('protein', '')}")
            m1, m2, m3 = st.columns(3)
            m1.metric("UniProt 登录号", info.get("accession") or "—")
            m2.metric("氨基酸长度", info.get("length") or "—")
            m3.metric("关联疾病", len(db["target_diseases"].get(picked, [])))
            if info.get("accession"):
                st.markdown(f"[在 UniProt 打开](https://www.uniprot.org/uniprotkb/"
                            f"{info['accession']})　"
                            f"[AlphaFold 三维结构](https://alphafold.ebi.ac.uk/entry/"
                            f"{info['accession']})")
            if info.get("locations"):
                st.markdown("**亚细胞定位**：" +
                            "　".join(f"`{x}`" for x in info["locations"]))
            if info.get("function"):
                st.caption(info["function"])
        else:
            st.markdown(f"### {picked}")
            st.caption("UniProt 中暂无该符号的人源已审阅记录。")
        with st.expander("蛋白三维结构（AlphaFold 预测模型，可旋转缩放）"):
            target_3d_block(picked)
        cids = db["target_compounds"].get(picked, [])
        with st.expander(f"作用该靶点的成分（{len(cids)} 个）", expanded=True):
            rows = [_compound_row(c) for c in cids]
            if rows:
                st.dataframe(pd.DataFrame(rows).sort_values(
                    "OB(%)", ascending=False), hide_index=True, width="stretch",
                    height=min(28 * (len(rows) + 1), 460))
        dis_of = db["target_diseases"].get(picked, [])
        with st.expander(f"关联疾病（{len(dis_of)} 种）", expanded=True):
            st.write("、".join(dis_of) if dis_of else "—")

# ------------------------------------------------------------------ 靶点预测
with tab_pred:
    st.header("靶点预测（TCMSP SVM / RF 模型）")
    st.caption("TCMSP 对每个成分-靶点对给出 **SVM**（支持向量机）与 **RF**（随机森林）"
               "两个独立模型的预测概率分（0-1），并标注是否经**实验验证**（DrugBank）。"
               "选择成分后按预测分排序展示作用靶点，并可直接查看靶点蛋白的三维结构。")

    scores = db["edge_scores"]
    comp_pool = sorted({c for (c, _t) in scores})
    n_pairs = len(scores)
    n_valid = sum(1 for v in scores.values() if v.get("validated"))
    m1, m2, m3 = st.columns(3)
    m1.metric("预测成分-靶点对", f"{n_pairs:,}")
    m2.metric("覆盖成分", len(comp_pool))
    m3.metric("实验验证对", f"{n_valid:,}（{n_valid / n_pairs:.0%}）")

    options = {f"{db['compounds'][c].get('name', c)}（{c}）": c
               for c in comp_pool}
    pick_label = st.selectbox("选择成分（按名称搜索）", sorted(options))
    picked_comp = options[pick_label]

    c1, c2, c3 = st.columns(3)
    svm_min = c1.slider("SVM 预测分下限", 0.0, 1.0, 0.0, 0.05)
    rf_min = c2.slider("RF 预测分下限", 0.0, 1.0, 0.0, 0.05)
    only_validated = c3.checkbox("仅显示实验验证靶点")

    tids = db["compound_targets"].get(picked_comp, [])
    rows = []
    for t in tids:
        s = scores.get((picked_comp, t))
        if s is None:
            s = {"svm": None, "rf": None, "validated": False}
        if s["svm"] is not None and (s["svm"] < svm_min or s["rf"] < rf_min):
            continue
        if only_validated and not s.get("validated"):
            continue
        info = db["target_info"].get(t.upper(), {})
        rows.append({
            "靶点": t,
            "蛋白名": info.get("protein", ""),
            "SVM 预测分": s["svm"] if s["svm"] is not None else "—",
            "RF 预测分": s["rf"] if s["rf"] is not None else "—",
            "综合分": (round((s["svm"] + s["rf"]) / 2, 3)
                     if s["svm"] is not None else "—"),
            "实验验证": "✔" if s.get("validated") else "",
            "关联疾病数": len(db["target_diseases"].get(t, [])),
        })
    st.markdown(f"**{db['compounds'][picked_comp].get('name', picked_comp)}** 的预测靶点："
                f"{len(tids)} 个（当前筛选后 {len(rows)} 个）")
    if rows:
        st.dataframe(pd.DataFrame(rows).sort_values(
            "综合分", ascending=False, na_position="last"),
            hide_index=True, width="stretch",
            height=min(28 * (len(rows) + 1), 460))
        chart_rows = [r for r in rows if r["综合分"] != "—"]
        if chart_rows:
            top = sorted(chart_rows, key=lambda r: -r["综合分"])[:15]
            st.bar_chart(pd.DataFrame(
                {"综合预测分": {r["靶点"]: r["综合分"] for r in top}}))
            st.caption("综合预测分 = (SVM + RF) / 2，柱状图展示前 15 个靶点。")

        st.divider()
        st.subheader("预测靶点的蛋白三维结构")
        pred_t = st.selectbox("选择靶点查看 3D 结构",
                              [r["靶点"] for r in rows])
        target_3d_block(pred_t, height=480)
        if tids:
            dis_of = sorted({d2 for t in tids
                             for d2 in db["target_diseases"].get(t, [])})
            with st.expander(f"该成分经靶点关联的疾病（{len(dis_of)} 种）"):
                st.write(_disease_chips(dis_of) if dis_of else "—")
    else:
        st.info("当前筛选条件下没有靶点，请调低预测分下限或取消验证过滤。")

# ------------------------------------------------------------------ 疾病
with tab_dis:
    st.header("疾病检索（四层反查）")
    kw = st.text_input("搜索疾病（英文名，不区分大小写）",
                       placeholder="例如：diabetes / cardiac / leukemia")
    view = sorted(db["diseases"],
                  key=lambda d: -len(db["disease_targets"].get(d, [])))
    if kw.strip():
        needle = kw.strip().lower()
        view = [d for d in view if needle in d.lower()]
    st.caption(f"共 **{len(view)}** 种疾病（未输入关键词时按关联靶点数排序）")
    if view:
        options = [f"{d}　关联靶点 {len(db['disease_targets'].get(d, []))} 个"
                   for d in view[:300]]
        pick = st.selectbox("选择疾病", options)
        picked = pick.split("　关联靶点")[0]
        tids = db["disease_targets"].get(picked, [])
        cids = {c for t in tids for c in db["target_compounds"].get(t, [])}
        herb_count: Dict[str, int] = {}
        for c in cids:
            for h in db["node_by_id"].get(c, {}).get("herbs") or []:
                herb_count[h] = herb_count.get(h, 0) + 1
        st.markdown(f"### {picked}")
        m1, m2, m3 = st.columns(3)
        m1.metric("关联靶点", len(tids))
        m2.metric("相关成分", len(cids))
        m3.metric("来源中药", len(herb_count))
        with st.expander(f"关联靶点（{len(tids)} 个，按关联成分数排序）", expanded=True):
            if tids:
                rows = [{"靶点": t,
                         "蛋白名": db["target_info"].get(t.upper(), {}).get("protein", ""),
                         "关联成分数": len(db["target_compounds"].get(t, [])),
                         "关联疾病数": len(db["target_diseases"].get(t, []))}
                        for t in tids]
                st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch",
                             height=min(28 * (len(rows) + 1), 460))
        with st.expander(f"来源中药（{len(herb_count)} 味，按覆盖成分数排序）"):
            if herb_count:
                rows = [{"中药": h, "覆盖成分数": n}
                        for h, n in sorted(herb_count.items(), key=lambda kv: -kv[1])]
                st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch",
                             height=min(28 * (len(rows) + 1), 460))
                top = dict(sorted(herb_count.items(), key=lambda kv: -kv[1])[:15])
                st.bar_chart(pd.DataFrame({"覆盖成分数": top}).T)
        with st.expander(f"相关成分（{len(cids)} 个）"):
            rows = [_compound_row(c) for c in sorted(cids)]
            if rows:
                st.dataframe(pd.DataFrame(rows).sort_values(
                    "来源中药数", ascending=False), hide_index=True,
                    width="stretch", height=min(28 * (len(rows) + 1), 460))

# ------------------------------------------------------------------ 网络图
with tab_net:
    st.header("四层网络可视化")
    st.caption("中药（橙）→ 成分（蓝）→ 靶点（绿）→ 疾病（红）。"
               "点击节点高亮邻居，可拖拽调整布局。")
    mode = st.radio("构建方式", ["按中药", "按疾病"], horizontal=True)
    show_disease = st.checkbox("包含疾病层", value=True)
    max_dis = st.slider("疾病节点上限（按关联数取前 N）", 5, 60, 20)

    nodes: List[dict] = []
    edges: List[dict] = []

    def add_compound_targets(mid: str) -> None:
        d = db["compounds"].get(mid)
        if d and not any(n["data"]["id"] == mid for n in nodes):
            nodes.append({"data": {"id": mid, "label": "COMPOUND",
                                   "name": d.get("name", mid)}})
        for t in db["compound_targets"].get(mid, []):
            tid = f"TGT::{t}"
            if not any(n["data"]["id"] == tid for n in nodes):
                nodes.append({"data": {"id": tid, "label": "TARGET", "name": t}})
            edges.append({"data": {"id": f"E::{mid}->{tid}", "source": mid,
                                   "target": tid, "label": "TARGETS"}})

    if mode == "按中药":
        picks = st.multiselect("选择中药（可多选）", sorted(db["herb_compounds"]),
                               default=sorted(db["herb_compounds"])[:2])
        for h in picks:
            hid = f"HERB::{h}"
            nodes.append({"data": {"id": hid, "label": "HERB", "name": h}})
            for mid in db["herb_compounds"].get(h, []):
                edges.append({"data": {"id": f"E::{hid}->{mid}", "source": hid,
                                       "target": mid, "label": "CONTAINS"}})
                add_compound_targets(mid)
    else:
        dis_pool = sorted(db["disease_targets"],
                          key=lambda d: -len(db["disease_targets"][d]))
        pick = st.selectbox("选择疾病", [f"{d}　关联靶点 {len(db['disease_targets'][d])} 个"
                                    for d in dis_pool])
        picked_dis = pick.split("　关联靶点")[0]
        did = f"DIS::{picked_dis}"
        nodes.append({"data": {"id": did, "label": "DISEASE", "name": picked_dis}})
        for t in db["disease_targets"].get(picked_dis, []):
            tid = f"TGT::{t}"
            nodes.append({"data": {"id": tid, "label": "TARGET", "name": t}})
            edges.append({"data": {"id": f"E::{tid}->{did}", "source": tid,
                                   "target": did, "label": "INDICATES"}})
            for c in db["target_compounds"].get(t, [])[:15]:
                add_compound_targets(c)

    if show_disease and mode == "按中药":
        tgt_ids = {n["data"]["id"] for n in nodes if n["data"]["label"] == "TARGET"}
        dis_count: Dict[str, int] = {}
        for t in tgt_ids:
            for d2 in db["target_diseases"].get(str(t).replace("TGT::", ""), []):
                dis_count[d2] = dis_count.get(d2, 0) + 1
        for d2, _ in sorted(dis_count.items(), key=lambda kv: -kv[1])[:max_dis]:
            did = f"DIS::{d2}"
            nodes.append({"data": {"id": did, "label": "DISEASE", "name": d2}})
        tgt_set = {str(n["data"]["id"]).replace("TGT::", "") for n in nodes
                   if n["data"]["label"] == "TARGET"}
        for t in tgt_set:
            for d2 in db["target_diseases"].get(t, []):
                if f"DIS::{d2}" in {n["data"]["id"] for n in nodes}:
                    edges.append({"data": {"id": f"E::TGT::{t}->DIS::{d2}",
                                           "source": f"TGT::{t}",
                                           "target": f"DIS::{d2}",
                                           "label": "INDICATES"}})

    if not nodes:
        st.warning("请先选择中药或疾病。")
    else:
        node_styles = [
            NodeStyle("HERB", "#E67E22", "name", "store"),
            NodeStyle("COMPOUND", "#2A629A", "name", "science"),
            NodeStyle("TARGET", "#27AE60", "name", "medical_services"),
            NodeStyle("DISEASE", "#C0392B", "name", "sick"),
        ]
        edge_styles = [
            EdgeStyle("CONTAINS", caption="label", directed=True),
            EdgeStyle("TARGETS", caption="label", directed=True),
            EdgeStyle("INDICATES", caption="label", directed=True),
        ]
        st_link_analysis({"nodes": nodes, "edges": edges},
                         layout="cose", node_styles=node_styles,
                         edge_styles=edge_styles, height=640)
        st.caption(f"当前网络：{len(nodes)} 个节点、{len(edges)} 条边。")

import streamlit as st
import pandas as pd
import numpy as np
import streamlit_ketcher as sk
from streamlit.components.v1 import html

from drugcomp import Drug
from diffdock import get_diffdock, get_binding_prob
from prob_success import compute_prob_clin_success
from mol_viewer import run_wrapper

df_filename = "data/belka_smiles_test.csv"

# ── Page config ───────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="MOracle",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="collapsed",
)

# ── Global CSS ────────────────────────────────────────────────────────────────
st.markdown("""
<style>
#MainMenu, footer, header { visibility: hidden; }

.stApp { background: #f1f5f9; }

.block-container {
    padding-top: 1.8rem !important;
    padding-bottom: 3rem !important;
    max-width: 1440px !important;
}

/* Selectbox — indigo tint inside the purple options expander */
[data-testid="stSelectbox"] > div > div {
    border-radius: 8px !important;
    border-color: #a5b4fc !important;
    background: #e0e7ff !important;
    color: #3730a3 !important;
}
[data-testid="stSelectbox"] [data-testid="stWidgetLabel"] p {
    color: #4338ca !important;
    font-weight: 600 !important;
    font-size: 12px !important;
}

/* Mode toggle radio — full-width 50/50 pills, no bullet */
[data-testid="stRadio"] { margin-bottom: 0 !important; }
[data-testid="stRadio"] > div {
    display: flex !important;
    width: 100% !important;
    gap: 8px !important;
}
[data-testid="stRadio"] > div label {
    flex: 1 !important;
    display: flex !important;
    align-items: center !important;
    justify-content: center !important;
    background: #e0e7ff !important;
    border: 1.5px solid #a5b4fc !important;
    border-radius: 8px !important;
    padding: 9px 0 !important;
    font-size: 13px !important;
    font-weight: 500 !important;
    color: #4338ca !important;
    cursor: pointer !important;
}
/* Hide the radio circle indicator */
[data-testid="stRadio"] > div label > div:first-child {
    display: none !important;
}
/* Remove the left gap the hidden bullet used to occupy */
[data-testid="stRadio"] > div label p {
    margin: 0 !important;
    padding: 0 !important;
}
[data-testid="stRadio"] > div label:has(input:checked) {
    background: #fef9c3 !important;
    border-color: #fde047 !important;
    color: #854d0e !important;
    font-weight: 600 !important;
}

/* Expander — default (About MOracle) */
[data-testid="stExpander"] {
    border: 1px solid #e2e8f0 !important;
    border-radius: 10px !important;
    background: white !important;
    box-shadow: 0 1px 3px rgba(0,0,0,0.04) !important;
}

/* Options expander — purple, only when it contains selectboxes */
[data-testid="stExpander"]:has([data-testid="stSelectbox"]) {
    border: 1.5px solid #a5b4fc !important;
    background: #eef2ff !important;
    box-shadow: 0 1px 6px rgba(79,70,229,0.10) !important;
}
[data-testid="stExpander"]:has([data-testid="stSelectbox"]) summary p {
    font-size: 12px !important;
    color: #4338ca !important;
    font-weight: 600 !important;
}
[data-testid="stExpander"]:has([data-testid="stSelectbox"]) summary svg {
    fill: #4338ca !important;
}
/* Equalise all label text inside the options expander */
[data-testid="stExpander"]:has([data-testid="stSelectbox"]) [data-testid="stMarkdownContainer"] p,
[data-testid="stExpander"]:has([data-testid="stSelectbox"]) [data-testid="stMarkdownContainer"] div {
    color: #4338ca !important;
    font-weight: 600 !important;
    font-size: 12px !important;
}

/* Wipe all backgrounds and borders inside the options expander */
[data-testid="stExpander"]:has([data-testid="stSelectbox"]) details,
[data-testid="stExpander"]:has([data-testid="stSelectbox"]) details > div,
[data-testid="stExpander"]:has([data-testid="stSelectbox"]) details div {
    background: transparent !important;
    box-shadow: none !important;
    border: none !important;
}
/* Restore just the dropdown control */
[data-testid="stExpander"]:has([data-testid="stSelectbox"]) [data-testid="stSelectbox"] > div > div {
    background: #e0e7ff !important;
    border: 1px solid #a5b4fc !important;
    border-radius: 8px !important;
    box-shadow:0 2px 8px rgba(79,70,229,0.18)
}

            /* Bordered containers */
[data-testid="stVerticalBlockBorderWrapper"] {
    background: white !important;
    border-radius: 12px !important;
    border-color: #e2e8f0 !important;
    box-shadow: 0 1px 4px rgba(0,0,0,0.05) !important;
    padding: 10px !important;
    box-sizing: border-box !important;
}
[data-testid="stVerticalBlockBorderWrapper"] > div {
    background: transparent !important;
}


/* Search input inside library — compact */
[data-testid="stVerticalBlockBorderWrapper"] [data-testid="stTextInput"] input {
    font-size: 12px !important;
    padding: 6px 10px !important;
    border-radius: 6px !important;
    border-color: #e2e8f0 !important;
}

[data-testid="stDataFrameResizable"] {
    border-radius: 8px !important;
    overflow: hidden !important;
}
</style>
""", unsafe_allow_html=True)

# ── Session state ─────────────────────────────────────────────────────────────
if "failed_retrieval"       not in st.session_state:
    st.session_state.failed_retrieval       = False
if "df"                     not in st.session_state:
    st.session_state.df                     = pd.read_csv(df_filename, sep=",", header=0)
if "selected_chemical"      not in st.session_state:
    st.session_state.selected_chemical      = st.session_state.df.iloc[0]
if "selected_chemical1"     not in st.session_state:
    st.session_state.selected_chemical1     = st.session_state.df.iloc[0]
if "count"                  not in st.session_state:
    st.session_state.count                  = 0
if "selected_index1"        not in st.session_state:
    st.session_state.selected_index1        = 0
if "selected_source_option" not in st.session_state:
    st.session_state.selected_source_option = "Belka"
if "main_chemical"          not in st.session_state:
    st.session_state.main_chemical          = st.session_state.df.iloc[0]


# ── UI helpers ────────────────────────────────────────────────────────────────

def render_hero():
    st.markdown("""
    <div style="background:#eef2ff; padding:22px 28px 18px; border-radius:16px;
                margin-bottom:20px; border:1.5px solid #a5b4fc;
                box-shadow:0 2px 8px rgba(79,70,229,0.18);">
      <div style="display:flex; align-items:center; justify-content:space-between;
                  flex-wrap:wrap; gap:16px;">
        <div style="display:flex; align-items:center; gap:14px;">
          <div style="width:46px; height:46px; border-radius:12px;
                      background:linear-gradient(135deg,#4f46e5,#7c3aed);
                      display:flex; align-items:center; justify-content:center;
                      font-size:24px; box-shadow:0 4px 12px rgba(79,70,229,0.35);">🧪</div>
          <div>
            <div style="font-size:26px; font-weight:800; color:#1e1b4b;
                        letter-spacing:-0.03em; line-height:1.1;">MOracle</div>
            <div style="font-size:12px; color:#6d6f6c; margin-top:2px;">
              Explainable clinical viability scoring in drug discovery
            </div>
          </div>
        </div>
      </div>
    </div>
    """, unsafe_allow_html=True)


def section_label(icon, title, right=""):
    """Card-level label: small caps inside a bordered container."""
    right_html = (
        f'<span style="font-size:10px;font-weight:700;color:#64748b;text-transform:uppercase;letter-spacing:0.09em;">{right}</span>'
        if right else ""
    )
    st.markdown(
        f'<div style="display:flex;align-items:center;justify-content:space-between;'
        f'margin:14px 0 6px;padding-bottom:7px;border-bottom:1px solid #e2e8f0;">'
        f'<div style="display:flex;align-items:center;gap:7px;">'
        f'<span style="font-size:12px;line-height:1;">{icon}</span>'
        f'<span style="font-size:10px;font-weight:700;color:#64748b;'
        f'text-transform:uppercase;letter-spacing:0.09em;">{title}</span>'
        f'</div>'
        f'{right_html}'
        f'</div>',
        unsafe_allow_html=True,
    )


def column_header(title, subtitle=""):
    """Column-level header: one tier above section_label."""
    sub = (f'<div style="font-size:11px;color:#94a3b8;margin-top:3px;">{subtitle}</div>'
           if subtitle else "")
    st.markdown(
        f'<div style="margin-bottom:14px;padding-bottom:10px;border-bottom:1.5px solid #e2e8f0;">'
        f'<div style="font-size:14px;font-weight:700;color:#1e293b;letter-spacing:-0.01em;">{title}</div>'
        f'{sub}</div>',
        unsafe_allow_html=True,
    )


def render_library(source_df, event_key, search_key, current_smiles, height=260,
                   show_sim_score=False):
    """
    Polished molecule library with search/filter.

    Returns (selected_row_series | None, filtered_df).
    Caller should check if selected_row_series is not None and update session state.
    """
    # ── Header row: label + count badge ──────────────────────────────────────
    query = st.text_input(
        "filter", placeholder="🔍  Filter by name…",
        label_visibility="collapsed", key=search_key
    )

    if query:
        mask = source_df["Name"].astype(str).str.contains(query, case=False, na=False)
        filtered = source_df[mask].reset_index(drop=True)
    else:
        filtered = source_df.reset_index(drop=True)

    n = len(filtered)
    count_label = f"{n} molecule{'s' if n != 1 else ''}"
    if query:
        count_label += f" matching <em style='color:#4338ca'>{query}</em>"

    # Currently selected indicator
    sel_name = ""
    matching = source_df[source_df["Smiles"].astype(str) == str(current_smiles)]
    if not matching.empty:
        sel_name = str(matching.iloc[0]["Name"])

    st.markdown(f"""
    <div style="display:flex; align-items:center; justify-content:space-between;
                margin:4px 0 6px; font-size:10px; color:#94a3b8; flex-wrap:wrap; gap:4px;">
      <span>{count_label}</span>
    </div>
    """, unsafe_allow_html=True)

    # ── Build display DataFrame ───────────────────────────────────────────────
    cols_to_show: list[str]
    if show_sim_score and "Sim. Score" in filtered.columns:
        display_df = pd.DataFrame({
            "Molecule": filtered["Name"].astype(str),
            "Similarity": filtered["Sim. Score"],
            "Phase": filtered["Trial Phase"].astype(str),
        })
        col_cfg = {
            "Molecule":   st.column_config.TextColumn("Molecule",   width="medium"),
            "Similarity": st.column_config.NumberColumn("Sim.", format="%.2f", width="small"),
            "Phase":      st.column_config.TextColumn("Phase",      width="small"),
        }
        cols_to_show = ["Molecule", "Similarity", "Phase"]
    else:
        smiles_series = filtered["Smiles"].astype(str)
        display_df = pd.DataFrame({
            "ID":    filtered["Name"].astype(str),
            "SMILES": smiles_series.apply(
                lambda s: s[:44] + "…" if len(s) > 44 else s
            ),
        })
        col_cfg = {
            "ID":     st.column_config.TextColumn("ID",    width="small"),
            "SMILES": st.column_config.TextColumn("SMILES", width="large"),
        }
        cols_to_show = ["ID", "SMILES"]

    event = st.dataframe(
        display_df[cols_to_show],
        height=height,
        key=event_key,
        on_select="rerun",
        selection_mode="single-row",
        hide_index=True,
        use_container_width=True,
        column_config=col_cfg,
    )

    selected_row = None
    if len(event.selection["rows"]) > 0:
        row = filtered.iloc[event.selection["rows"][0]]
        if str(row.get("Smiles", row.get("smiles", ""))) != str(current_smiles):
            selected_row = row

    return selected_row, filtered


def get_similar_drugs_data(smiles):
    try:
        curr_drug = Drug(smiles)
        similar_drugs = curr_drug.get_similar_drugs()
        st.session_state.failed_retrieval = False
        df = pd.DataFrame.from_dict([{
            "Name": drug.name,
            "Smiles": drug.smiles,
            "Sim. Score": np.round(float(sim), 2),
            "ID": drug.chembl_id,
            "Trial Phase": drug.max_phase,
            "Ind.": drug.indication
        } for drug, sim in similar_drugs])
        df["Trial Phase"].fillna("Preclinical", inplace=True)
        return df.sort_values(by="Sim. Score", ascending=False)
    except Exception:
        st.session_state.failed_retrieval = True
        return pd.DataFrame(columns=["Name", "Smiles", "Sim. Score", "ID", "Trial Phase", "Ind."])


def _hex_to_rgba(hex_color, alpha):
    h = hex_color.lstrip("#")
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    return f"rgba({r},{g},{b},{alpha})"


def _band_color(value, low_thresh, high_thresh):
    try:
        v = float(value)
        if v >= high_thresh:
            return "#15803d"
        if v >= low_thresh:
            return "#b45309"
        return "#b91c1c"
    except (TypeError, ValueError):
        return "#64748b"


def _verdict(value, low_thresh, high_thresh):
    try:
        v = float(value)
        if v >= high_thresh:
            return "HIGH"
        if v >= low_thresh:
            return "MODERATE"
        return "LOW"
    except (TypeError, ValueError):
        return "N/A"


def _metric_card(title, value, verdict, color):
    """Render a single metric card using st.markdown with no leading whitespace."""
    bg  = _hex_to_rgba(color, 0.07)
    bdr = _hex_to_rgba(color, 0.22)
    vbdr = _hex_to_rgba(color, 0.30)
    card = (
        f'<div style="background:{bg};border:1px solid {bdr};border-radius:10px;'
        f'padding:16px 10px;text-align:center;">'
        f'<div style="font-size:9px;font-weight:700;color:#64748b;text-transform:uppercase;'
        f'letter-spacing:0.1em;margin-bottom:8px;">{title}</div>'
        f'<div style="font-size:34px;font-weight:800;color:{color};line-height:1;margin-bottom:8px;">{value}</div>'
        f'<span style="background:white;color:{color};border:1px solid {vbdr};'
        f'font-size:9px;font-weight:700;padding:2px 9px;border-radius:20px;'
        f'text-transform:uppercase;letter-spacing:0.06em;">{verdict}</span>'
        f'</div>'
    )
    st.markdown(card, unsafe_allow_html=True)


def display_molecule_data(label, chemical,
                          prob_success=None, confidence=None, prob=None):
    name   = chemical["Name"]
    smiles = chemical["Smiles"]

    ps_str   = str(prob_success) if prob_success is not None else "—"
    p_str    = str(prob)         if prob         is not None else "—"
    conf_str = str(confidence)   if confidence   is not None else "—"

    ps_color   = _band_color(prob_success, 0.08, 0.15)
    ps_verdict = _verdict(prob_success, 0.08, 0.15)
    p_color    = _band_color(prob, 0.05, 0.10)
    p_verdict  = _verdict(prob, 0.05, 0.10)

    # Row label (omitted when label is None)
    if label:
        st.markdown(
            f'<div style="display:flex;align-items:center;gap:7px;margin-bottom:10px;">'
            f'<div style="width:6px;height:6px;border-radius:50%;background:#6366f1;flex-shrink:0;"></div>'
            f'<span style="font-size:10px;font-weight:700;color:#64748b;text-transform:uppercase;'
            f'letter-spacing:0.09em;">{label}</span>'
            f'<span style="font-size:10px;color:#94a3b8;">&nbsp;·&nbsp;{name}</span>'
            f'</div>',
            unsafe_allow_html=True,
        )

    c1, c2, c3 = st.columns(3)
    with c1:
        _metric_card("Clinical Success",    ps_str,   ps_verdict, ps_color)
    with c2:
        _metric_card("Binding Probability", p_str,    p_verdict,  p_color)
    with c3:
        _metric_card("DiffDock Score",      conf_str, "DOCKING",  "#4338ca")

    st.markdown(
        f'<div style="font-size:10px;color:#94a3b8;margin-top:8px;word-break:break-all;">'
        f'<span style="font-weight:600;color:#64748b;">SMILES</span>&nbsp;'
        f'<span style="font-family:\'SF Mono\',\'Fira Mono\',monospace;color:#475569;">{smiles}</span>'
        f'</div>',
        unsafe_allow_html=True,
    )


def _fmt(value, decimals=6):
    """Format a float to string, returning 'N/A' for None."""
    return str(round(value, decimals)) if value is not None else "N/A"


def display_comparison(name_p, name_r,
                       ps_p, ps_r,
                       conf_p, conf_r,
                       prob_p, prob_r,
                       protein):
    """Render full-width head-to-head comparison panel."""

    def _try_float(v):
        try:
            return float(v)
        except (TypeError, ValueError):
            return None

    def _winner(vp, vr):
        fp, fr = _try_float(vp), _try_float(vr)
        if fp is None and fr is None:
            return "tie"
        if fp is None:
            return "ref"
        if fr is None:
            return "primary"
        if fp > fr:
            return "primary"
        if fr > fp:
            return "ref"
        return "tie"

    w_ps   = _winner(ps_p,   ps_r)
    w_prob = _winner(prob_p, prob_r)
    w_conf = _winner(conf_p, conf_r)

    p_wins = sum(1 for w in [w_ps, w_prob, w_conf] if w == "primary")
    r_wins = sum(1 for w in [w_ps, w_prob, w_conf] if w == "ref")

    if p_wins > r_wins:
        vtext = f"Primary leads {p_wins}–{r_wins}"
        vbg, vcol, vbdr = "#f0fdf4", "#15803d", "#bbf7d0"
    elif r_wins > p_wins:
        vtext = f"Reference leads {r_wins}–{p_wins}"
        vbg, vcol, vbdr = "#fffbeb", "#b45309", "#fde68a"
    else:
        vtext = "All tied"
        vbg, vcol, vbdr = "#f8fafc", "#64748b", "#e2e8f0"

    BEST = (
        '<span style="background:#15803d;color:white;font-size:8px;font-weight:700;'
        'padding:2px 6px;border-radius:20px;text-transform:uppercase;'
        'letter-spacing:0.05em;margin-left:6px;">BEST</span>'
    )

    def _cell(val_str, verdict_str, color, is_winner):
        bg  = "#f0fdf4" if is_winner else "#f8fafc"
        bdr = "#bbf7d0" if is_winner else "#e2e8f0"
        bdg = BEST if is_winner else ""
        return (
            f'<div style="flex:1;background:{bg};border:1px solid {bdr};'
            f'border-radius:8px;padding:14px 16px;">'
            f'<div style="font-size:22px;font-weight:800;color:{color};line-height:1;">{val_str}</div>'
            f'<div style="margin-top:5px;display:flex;align-items:center;">'
            f'<span style="font-size:9px;font-weight:700;color:{color};'
            f'text-transform:uppercase;letter-spacing:0.08em;">{verdict_str}</span>'
            f'{bdg}</div>'
            f'</div>'
        )

    def _conf_cell(val_str, is_winner):
        bg  = "#f0fdf4" if is_winner else "#f8fafc"
        bdr = "#bbf7d0" if is_winner else "#e2e8f0"
        bdg = BEST if is_winner else ""
        return (
            f'<div style="flex:1;background:{bg};border:1px solid {bdr};'
            f'border-radius:8px;padding:14px 16px;">'
            f'<div style="font-size:22px;font-weight:800;color:#4338ca;line-height:1;">{val_str}</div>'
            f'<div style="margin-top:5px;display:flex;align-items:center;">'
            f'<span style="font-size:9px;font-weight:700;color:#6366f1;'
            f'text-transform:uppercase;letter-spacing:0.08em;">DOCKING</span>'
            f'{bdg}</div>'
            f'</div>'
        )

    def _row(label, sublabel, cell_p, cell_r):
        return (
            f'<div style="display:flex;gap:10px;align-items:stretch;margin-bottom:8px;">'
            f'<div style="width:150px;flex-shrink:0;display:flex;flex-direction:column;'
            f'justify-content:center;padding:10px 0;">'
            f'<div style="font-size:11px;font-weight:700;color:#1e293b;">{label}</div>'
            f'<div style="font-size:10px;color:#94a3b8;margin-top:2px;">{sublabel}</div>'
            f'</div>'
            f'{cell_p}{cell_r}'
            f'</div>'
        )

    def _fmt_band(v, lo, hi):
        f = _try_float(v)
        if f is None:
            return (str(v) if v is not None else "N/A"), "#64748b", "N/A"
        return str(round(f, 4)), _band_color(f, lo, hi), _verdict(f, lo, hi)

    ps_p_v, ps_p_c, ps_p_l = _fmt_band(ps_p,   0.08, 0.15)
    ps_r_v, ps_r_c, ps_r_l = _fmt_band(ps_r,   0.08, 0.15)
    pb_p_v, pb_p_c, pb_p_l = _fmt_band(prob_p, 0.05, 0.10)
    pb_r_v, pb_r_c, pb_r_l = _fmt_band(prob_r, 0.05, 0.10)
    cf_p_v = _fmt(conf_p) if conf_p is not None else "N/A"
    cf_r_v = _fmt(conf_r) if conf_r is not None else "N/A"

    np_s = str(name_p)[:18] + ("…" if len(str(name_p)) > 18 else "")
    nr_s = str(name_r)[:18] + ("…" if len(str(name_r)) > 18 else "")

    hdr = (
        f'<div style="display:flex;align-items:center;justify-content:space-between;'
        f'margin-bottom:16px;padding-bottom:10px;border-bottom:1.5px solid #e2e8f0;">'
        f'<div style="display:flex;align-items:center;gap:8px;">'
        f'<span style="font-size:12px;">⚖️</span>'
        f'<span style="font-size:10px;font-weight:700;color:#64748b;'
        f'text-transform:uppercase;letter-spacing:0.09em;">Head-to-Head Comparison</span>'
        f'</div>'
        f'<div style="display:flex;align-items:center;gap:12px;">'
        f'<span style="font-size:10px;font-weight:700;color:#64748b;text-transform:uppercase;letter-spacing:0.09em;">Protein: {protein}</span>'
        f'<span style="background:{vbg};color:{vcol};border:1px solid {vbdr};'
        f'font-size:11px;font-weight:700;padding:3px 12px;border-radius:20px;">{vtext}</span>'
        f'</div>'
        f'</div>'
    )
    col_hdrs = (
        f'<div style="display:flex;gap:10px;margin-bottom:8px;">'
        f'<div style="width:150px;flex-shrink:0;"></div>'
        f'<div style="flex:1;font-size:9px;font-weight:700;color:#475569;'
        f'text-transform:uppercase;letter-spacing:0.07em;padding:0 16px;">'
        f'PRIMARY MOLECULE</div>'
        f'<div style="flex:1;font-size:9px;font-weight:700;color:#475569;'
        f'text-transform:uppercase;letter-spacing:0.07em;padding:0 16px;">'
        f'REFERENCE MOLECULE</div>'
        f'</div>'
    )
    row_ps   = _row("Probability of Clinical Success",    "Similarity to approved molecules based on ChEMBL data",
                    _cell(ps_p_v, ps_p_l, ps_p_c, w_ps   == "primary"),
                    _cell(ps_r_v, ps_r_l, ps_r_c, w_ps   == "ref"))
    row_prob = _row("Binding Probability", "Estimate of successful protein binding probability given by DiffDock",
                    _cell(pb_p_v, pb_p_l, pb_p_c, w_prob == "primary"),
                    _cell(pb_r_v, pb_r_l, pb_r_c, w_prob == "ref"))
    row_conf = _row("DiffDock Score",      "Confidence of 3D docking prediction given by DiffDock",
                    _conf_cell(cf_p_v, w_conf == "primary"),
                    _conf_cell(cf_r_v, w_conf == "ref"))

    st.markdown(
        f'<div style="padding:4px 2px;">{hdr}{col_hdrs}{row_ps}{row_prob}{row_conf}</div>',
        unsafe_allow_html=True,
    )


def display_single_analysis(name, protein, ps, conf, prob):
    """Single-molecule analysis panel matching the comparative row/cell layout."""

    def _try_float(v):
        try:
            return float(v)
        except (TypeError, ValueError):
            return None

    def _fmt_band(v, lo, hi):
        f = _try_float(v)
        if f is None:
            return (str(v) if v is not None else "N/A"), "#64748b", "N/A"
        return str(round(f, 4)), _band_color(f, lo, hi), _verdict(f, lo, hi)

    def _cell(val_str, verdict_str, color):
        return (
            f'<div style="flex:1;background:#f8fafc;border:1px solid #e2e8f0;'
            f'border-radius:8px;padding:14px 16px;">'
            f'<div style="font-size:22px;font-weight:800;color:{color};line-height:1;">{val_str}</div>'
            f'<div style="margin-top:5px;">'
            f'<span style="font-size:9px;font-weight:700;color:{color};'
            f'text-transform:uppercase;letter-spacing:0.08em;">{verdict_str}</span>'
            f'</div></div>'
        )

    def _conf_cell(val_str):
        return (
            f'<div style="flex:1;background:#f8fafc;border:1px solid #e2e8f0;'
            f'border-radius:8px;padding:14px 16px;">'
            f'<div style="font-size:22px;font-weight:800;color:#4338ca;line-height:1;">{val_str}</div>'
            f'<div style="margin-top:5px;">'
            f'<span style="font-size:9px;font-weight:700;color:#6366f1;'
            f'text-transform:uppercase;letter-spacing:0.08em;">DOCKING</span>'
            f'</div></div>'
        )

    def _row(label, sublabel, cell):
        return (
            f'<div style="display:flex;gap:10px;align-items:stretch;margin-bottom:8px;">'
            f'<div style="width:150px;flex-shrink:0;display:flex;flex-direction:column;'
            f'justify-content:center;padding:10px 0;">'
            f'<div style="font-size:11px;font-weight:700;color:#1e293b;">{label}</div>'
            f'<div style="font-size:10px;color:#94a3b8;margin-top:2px;">{sublabel}</div>'
            f'</div>'
            f'{cell}'
            f'</div>'
        )

    ps_v,  ps_c,  ps_l  = _fmt_band(ps,   0.08, 0.15)
    pb_v,  pb_c,  pb_l  = _fmt_band(prob, 0.05, 0.10)
    cf_v = _fmt(conf) if conf is not None else "N/A"
    name_s = str(name)[:18] + ("…" if len(str(name)) > 18 else "")

    hdr = (
        f'<div style="display:flex;align-items:center;justify-content:space-between;'
        f'margin-bottom:16px;padding-bottom:10px;border-bottom:1.5px solid #e2e8f0;">'
        f'<div style="display:flex;align-items:center;gap:8px;">'
        f'<span style="font-size:12px;">📊</span>'
        f'<span style="font-size:10px;font-weight:700;color:#64748b;'
        f'text-transform:uppercase;letter-spacing:0.09em;">Numerical Prediction Analysis</span>'
        f'</div>'
        f'<span style="font-size:10px;font-weight:700;color:#64748b;text-transform:uppercase;letter-spacing:0.09em;">Protein: {protein}</span>'
        f'</div>'
    )
    col_hdr = (
        f'<div style="display:flex;gap:10px;margin-bottom:8px;">'
        f'<div style="width:150px;flex-shrink:0;"></div>'
        f'<div style="flex:1;font-size:9px;font-weight:700;color:#475569;'
        f'text-transform:uppercase;letter-spacing:0.07em;padding:0 16px;">'
        f'PRIMARY MOLECULE</div>'
        f'</div>'
    )
    row_ps   = _row("Probability of Clinical Success", "Similarity to approved molecules based on ChEMBL data",
                    _cell(ps_v, ps_l, ps_c))
    row_prob = _row("Binding Probability", "Estimate of successful protein binding probability given by DiffDock",
                    _cell(pb_v, pb_l, pb_c))
    row_conf = _row("DiffDock Score", "Confidence of 3D docking prediction given by DiffDock",
                    _conf_cell(cf_v))

    st.markdown(
        f'<div style="padding:4px 2px;">{hdr}{col_hdr}{row_ps}{row_prob}{row_conf}</div>',
        unsafe_allow_html=True,
    )


def _maybe_add_molecule(new_smiles):
    if new_smiles and new_smiles not in st.session_state.df["Smiles"].tolist():
        st.session_state.count += 1
        new_row = pd.DataFrame({
            "Name": f"new_mol{st.session_state.count}",
            "Smiles": new_smiles,
            "ID": "", "Trial Phase": "", "Ind.": ""
        }, index=[0])
        st.session_state.df = pd.concat([new_row, st.session_state.df], ignore_index=True)
        st.rerun()


# ── Hero + mode selector ──────────────────────────────────────────────────────

render_hero()

# with st.expander("ℹ  About MOracle"):
#     st.write("""
#         MOracle scores drug molecules for clinical viability by combining DiffDock binding
#         predictions with similarity-weighted historical trial outcomes from ChEMBL.
#         - **Single Molecule** — deep-dive on one compound.
#         - **Comparative** — side-by-side analysis of two molecules.
#         - **Proteins supported**: BRD4, sEH, HSA (Belka-v1 dataset).
#     """)

with st.expander("⚙  Options", expanded=False):
    st.markdown(
        '<div style="font-size:12px;font-weight:600;color:#a5b4fc;margin-bottom:6px;">Analysis Mode</div>',
        unsafe_allow_html=True,
    )
    layout_mode = st.radio(
        "Mode", ["Single Molecule", "Comparative"],
        horizontal=True, label_visibility="collapsed"
    )
    if layout_mode == "Single Molecule":
        protein_mode = st.selectbox("Protein Target", options=["BRD4", "sEH", "HSA"])
    else:
        _ctrl1, _ctrl2 = st.columns(2)
        with _ctrl1:
            protein_comp = st.selectbox(
                "Protein Target", options=["BRD4", "sEH", "HSA"], key="protein_comp"
            )
        with _ctrl2:
            st.session_state.selected_source_option = st.selectbox(
                "Reference Molecule Database",
                options=["Belka", "ChEMBL"],
                index=["Belka", "ChEMBL"].index(st.session_state.selected_source_option),
            )


# ── Single Molecule Mode ──────────────────────────────────────────────────────

if layout_mode == "Single Molecule":

    # ── Library ───────────────────────────────────────────────────────────────
    with st.container(border=True):
        section_label("🔍", "Molecule Library")
        selected_row, _ = render_library(
            source_df=st.session_state.df,
            event_key="data1",
            search_key="search1",
            current_smiles=st.session_state.selected_chemical["Smiles"],
            height=260,
        )
    if selected_row is not None:
        st.session_state.selected_chemical = selected_row
        st.rerun()

    chemical = st.session_state.selected_chemical.copy()
    st.session_state.main_chemical = chemical

    # ── Editor ────────────────────────────────────────────────────────────────
    with st.container(border=True):
        section_label("✏️", "Molecule Editor")
        new_smiles = sk.st_ketcher(chemical["Smiles"], key=str(chemical["Name"]) + "_")
    _maybe_add_molecule(new_smiles)

    similar_drugs_df = get_similar_drugs_data(chemical["Smiles"])
    prob_success = np.round(compute_prob_clin_success(similar_drugs_df), 2)
    filename, confidence = get_diffdock(protein_mode, chemical["Smiles"])
    prob = get_binding_prob(protein_mode, chemical["Smiles"])

    # ── 3D Docking View ───────────────────────────────────────────────────────
    with st.container(border=True):
        section_label("🔬", "DiffDock - 3D Docking View", right=f"Protein: {protein_mode}")
        html(run_wrapper(filename), height=620)

    # ── Analysis Summary ──────────────────────────────────────────────────────
    with st.container(border=True):
        display_single_analysis(
            chemical["Name"], protein_mode,
            prob_success, confidence, prob
        )


# ── Comparative Mode ──────────────────────────────────────────────────────────

else:

    col1, spacer, col2 = st.columns([1, 0.03, 1])

    # ── Left: primary molecule ────────────────────────────────────────────────
    with col1:
        column_header("Primary Molecule", "Select and edit the molecule under investigation")

        with st.container(border=True):
            section_label("🔍", "Molecule Library")
            selected_row, _ = render_library(
                source_df=st.session_state.df,
                event_key="data1",
                search_key="search1",
                current_smiles=st.session_state.selected_chemical["Smiles"],
                height=240,
            )
        if selected_row is not None:
            st.session_state.selected_chemical = selected_row
            st.rerun()

        chemical = st.session_state.selected_chemical.copy()
        st.session_state.main_chemical = chemical

        with st.container(border=True):
            section_label("✏️", "Molecule Editor", right=chemical["Name"])
            new_smiles = sk.st_ketcher(chemical["Smiles"], key=str(chemical["Name"]) + "_")
        _maybe_add_molecule(new_smiles)

        similar_drugs_df = get_similar_drugs_data(chemical["Smiles"])
        prob_success = np.round(compute_prob_clin_success(similar_drugs_df), 2)
        filename, confidence = get_diffdock(protein_comp, chemical["Smiles"])
        prob = get_binding_prob(protein_comp, chemical["Smiles"])

        with st.container(border=True):
            section_label("🔬", "3D Docking View", right=f"Protein {protein_comp}")
            html(run_wrapper(filename), height=620)

    # ── Spacer ────────────────────────────────────────────────────────────────
    with spacer:
        st.markdown(
            '<div style="width:1px;background:#e2e8f0;min-height:600px;margin:0 auto;"></div>',
            unsafe_allow_html=True,
        )

    # ── Right: reference molecule ─────────────────────────────────────────────
    with col2:
        column_header("Reference Molecule", "Compare against Belka dataset or ChEMBL similar drugs")

        if st.session_state.selected_source_option == "Belka":
            with st.container(border=True):
                section_label("🔍", "Molecule Library")
                selected_row_ref, _ = render_library(
                    source_df=st.session_state.df,
                    event_key=str(st.session_state.main_chemical["Name"]) + "_ref2",
                    search_key="search_ref2",
                    current_smiles=st.session_state.selected_chemical1["Smiles"],
                    height=240,
                )
            if selected_row_ref is not None:
                st.session_state.selected_chemical1 = selected_row_ref
                st.rerun()

            chemical_ref = st.session_state.selected_chemical1.copy()
            similar_drugs_df_ref = get_similar_drugs_data(chemical_ref["Smiles"])
            prob_success_ref = np.round(compute_prob_clin_success(similar_drugs_df_ref), 2)

        else:  # ChEMBL
            st.session_state.similar_drugs_df = get_similar_drugs_data(
                st.session_state.main_chemical["Smiles"]
            )
            with st.container(border=True):
                section_label("🔍", "Similar Drugs · ChEMBL")
                selected_row_ref, _ = render_library(
                    source_df=st.session_state.similar_drugs_df,
                    event_key=str(st.session_state.main_chemical["Name"]) + "_chembl2",
                    search_key="search_chembl2",
                    current_smiles=st.session_state.selected_chemical1.get(
                        "Smiles", st.session_state.selected_chemical1.get("smiles", "")
                    ),
                    height=240,
                    show_sim_score=True,
                )
            if selected_row_ref is not None:
                st.session_state.selected_chemical1 = selected_row_ref
                st.rerun()

            chemical_ref = st.session_state.selected_chemical1.copy()
            prob_success_ref = "N/A"

        with st.container(border=True):
            section_label("✏️", "Molecule Editor", right=chemical_ref["Name"])
            new_smiles_ref = sk.st_ketcher(
                chemical_ref["Smiles"], key=str(chemical_ref["Name"]) + "_1"
            )
        _maybe_add_molecule(new_smiles_ref)

        filename_ref, confidence_ref = get_diffdock(protein_comp, chemical_ref["Smiles"])
        prob_ref = get_binding_prob(protein_comp, chemical_ref["Smiles"])

        with st.container(border=True):
            section_label("🔬", "3D Docking View", right=f"Protein {protein_comp}")
            html(run_wrapper(filename_ref), height=620)

    # ── Full-width comparison panel ────────────────────────────────────────────
    with st.container(border=True):
        display_comparison(
            chemical["Name"],  chemical_ref["Name"],
            prob_success,      prob_success_ref,
            confidence,        confidence_ref,
            prob,              prob_ref,
            protein_comp,
        )

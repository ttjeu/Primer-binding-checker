import re
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional

import pandas as pd
import streamlit as st
from io import BytesIO
import streamlit.components.v1 as components
from openpyxl.utils import get_column_letter


# -----------------------------
# 1) Primer definitions (IUPAC + explicit variants both supported)
# -----------------------------
@dataclass(frozen=True)
class Primer:
    name: str
    seq: str  # may include IUPAC
    variants: Optional[List[str]] = None  # if provided, use these


# From your image (fixed set)
PRIMERS: List[Primer] = [
    Primer(
        name="337F",
        seq="GACTCCTACGGGAGGCWGCAG",
        variants=[
            "GACTCCTACGGGAGGCAGCAG",
            "GACTCCTACGGGAGGCTGCAG",
        ]
    ),
    Primer(name="518F", seq="CCAGCAGCCGCGGTAATACG"),
    Primer(name="785F", seq="GGATTAGATACCCTGGTA"),
    Primer(name="800R", seq="TACCAGGGTATCTAATCC"),
    Primer(
        name="907R",
        seq="CCGTCAATTCMTTTRAGTTT",
        variants=[
            "CCGTCAATTCATTTAAGTTT",
            "CCGTCAATTCATTTGAGTTT",
            "CCGTCAATTCCTTTAAGTTT",
            "CCGTCAATTCCTTTGAGTTT",
        ]
    ),
    Primer(name="1100R", seq="GGGTTGCGCTCGTTG"),
]

# IUPAC map (extend if needed)
IUPAC: Dict[str, str] = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "W": "AT",  # weak
    "M": "AC",  # amino
    "R": "AG",  # purine
    "Y": "CT",  # pyrimidine
    "N": "ACGT",
}

COMP: Dict[str, str] = {
    "A": "T", "C": "G", "G": "C", "T": "A",
    "W": "W",  # complement of W(AT) is W
    "M": "K",  # complement of M(AC) is K(GT) - not used unless primer contains it after rc
    "R": "Y",
    "Y": "R",
    "N": "N",
    # if K appears (from complement), handle it:
    "K": "M",  # K(GT) <-> M(AC)
}


def sanitize_seq(s: str) -> str:
    s = s.upper()
    s = re.sub(r"[^ACGTNWMRYK]", "", s)  # keep IUPAC we use
    return s


def revcomp(seq: str) -> str:
    seq = sanitize_seq(seq)
    return "".join(COMP.get(b, "N") for b in reversed(seq))


def expand_iupac(seq: str, limit: int = 128) -> List[str]:
    """
    Expand IUPAC into explicit variants. limit is a safety cap.
    """
    seq = sanitize_seq(seq)
    variants = [""]

    for b in seq:
        choices = IUPAC.get(b)
        if choices is None:
            choices = "N"
        new_vars = []
        for v in variants:
            for c in choices:
                new_vars.append(v + c)
                if len(new_vars) > limit:
                    # Stop expanding further; fallback to partial set
                    return new_vars[:limit]
        variants = new_vars

    return variants


# -----------------------------
# 2) Matching / scoring
# -----------------------------
@dataclass
class MatchResult:
    primer: str
    call: str
    orientation: str
    start_1b: int
    end_1b: int
    mismatches: int
    matched_seq: str
    primer_used: str


def best_match_in_seq(
    seq: str,
    primer_variants: List[str],
) -> Tuple[int, int, str, str]:
    """
    Returns:
    (best_start0, best_mismatches, matched_seq, primer_used)
    """
    best = (-1, 10**9, "", "")
    n = len(seq)

    for pv in primer_variants:
        L = len(pv)
        if L == 0 or n < L:
            continue

        for start in range(0, n - L + 1):
            window = seq[start:start + L]
            mm = 0

            for i in range(L):
                if pv[i] != window[i]:
                    mm += 1

            cand = (start, mm, window, pv)

            if best[0] == -1:
                best = cand
            else:
                if cand[1] < best[1]:
                    best = cand

    return best




def find_best_for_primer(
    seq: str,
    primer: Primer,
) -> MatchResult:
    seq = sanitize_seq(seq)

    if primer.variants and len(primer.variants) > 0:
        f_vars = [sanitize_seq(v) for v in primer.variants]
    else:
        f_vars = expand_iupac(primer.seq)

    r_vars = [revcomp(v) for v in f_vars]

    f_best = best_match_in_seq(seq, f_vars)
    r_best = best_match_in_seq(seq, r_vars)

    pick = f_best
    orientation = "F"

    if r_best[0] != -1:
        if f_best[0] == -1 or r_best[1] < f_best[1]:
            pick = r_best
            orientation = "R"

    start0, mm, window, primer_used = pick

    if start0 < 0:
        start0, mm, window, primer_used = 0, 0, "", ""

    call = "Found" if mm == 0 else "Mismatch"

    L = len(window) if window else (len(primer_used) if primer_used else len(sanitize_seq(primer.seq)))

    return MatchResult(
        primer=primer.name,
        call=call,
        orientation=orientation,
        start_1b=start0 + 1,
        end_1b=start0 + L,
        mismatches=mm,
        matched_seq=window,
        primer_used=primer_used if primer_used else sanitize_seq(primer.seq),
    )




# -----------------------------
# 3) Visualization (track)
# -----------------------------
def color_mismatch_html(primer_used: str, matched_seq: str) -> str:
    """
    primer_used vs matched_seq를 비교해서 mismatch 문자만 빨간색으로 표시한 HTML 문자열 반환
    """
    if not primer_used or not matched_seq or len(primer_used) != len(matched_seq):
        return matched_seq

    out = []
    for p, s in zip(primer_used, matched_seq):
        if p == s:
            out.append(s)
        else:
            out.append(f"<span style='color:#ff4b4b;font-weight:700'>{s}</span>")
    return "".join(out)

def build_excel_bytes(summary_df: pd.DataFrame, detail_df: pd.DataFrame) -> bytes:
    output = BytesIO()

    with pd.ExcelWriter(output, engine="openpyxl") as writer:
        summary_df.to_excel(writer, sheet_name="Summary", index=False)
        detail_df.to_excel(writer, sheet_name="Detail", index=False)

        for sheet_name in ["Summary", "Detail"]:
            ws = writer.book[sheet_name]

            for col_idx, column_cells in enumerate(ws.columns, start=1):
                header = str(column_cells[0].value) if column_cells[0].value is not None else ""
                column_letter = get_column_letter(col_idx)

                max_length = 0
                for cell in column_cells:
                    cell_value = "" if cell.value is None else str(cell.value)
                    if len(cell_value) > max_length:
                        max_length = len(cell_value)

                if header in ["primer_used", "matched_seq", "sites"]:
                    adjusted_width = min(max(max_length + 2, 20), 80)
                else:
                    adjusted_width = min(max(max_length + 2, 12), 30)

                ws.column_dimensions[column_letter].width = adjusted_width

    output.seek(0)
    return output.getvalue()

def parse_sequence_and_header_info(text: str):
    """
    FASTA header 예:
    >H260310-R01_G01_S260207_1604_1_907R.ab1 1056
    """
    text = str(text).strip()
    lines = [line.strip() for line in text.splitlines() if line.strip()]

    header_name = ""
    seq_lines = lines

    if lines and lines[0].startswith(">"):
        header_name = lines[0][1:].split()[0]   # 공백 전까지만
        seq_lines = lines[1:]

    seq = sanitize_seq("".join(seq_lines))

    sample_name = ""
    read_label = ""

    if header_name:
        base = header_name.replace(".ab1", "")
        parts = base.split("_")

        # 예: H260310-R01_G01_S260207_1604_1_907R
        if len(parts) >= 6:
            sample_name = "_".join(parts[1:4])   # G01_S260207_1604
            read_label = parts[-1]              # 907R

    return {
        "sequence": seq,
        "sample_name": sample_name,
        "read_label": read_label,
        "header_name": header_name,
    }

def make_result_filename(order_no: str, sample_name: str, read_label: str, idx: int = 1):
    parts = [order_no]
    if sample_name:
        parts.append(sample_name)
    if read_label:
        parts.append(read_label)
    else:
        parts.append(f"result{idx}")
    return "_".join(parts) + ".xlsx"

def color_call_badge_html(call_value: str) -> str:
    if call_value == "Found":
        bg = "#0f9d58"
    else:
        bg = "#d97706"

    return (
        f"<span style='background:{bg};color:white;padding:3px 8px;"
        f"border-radius:999px;font-weight:700'>{call_value}</span>"
    )

def color_yesno_badge_html(value: str) -> str:
    bg = "#b91c1c" if value == "Yes" else "#374151"
    return (
        f"<span style='background:{bg};color:white;padding:3px 8px;"
        f"border-radius:999px;font-weight:700'>{value}</span>"
    )


def find_all_matches_in_seq(
    seq: str,
    primer_variants: List[str],
    orientation_label: str,
    max_mismatches: int = 2,
):
    hits = []
    n = len(seq)

    for pv in primer_variants:
        L = len(pv)
        if L == 0 or n < L:
            continue

        for start in range(0, n - L + 1):
            window = seq[start:start + L]
            mm = 0

            for i in range(L):
                if pv[i] != window[i]:
                    mm += 1

            if mm <= max_mismatches:
                hits.append({
                    "orientation": orientation_label,
                    "start_1b": start + 1,
                    "end_1b": start + L,
                    "mismatches": mm,
                    "primer_used": pv,
                    "matched_seq": window,
                })

    unique_hits = []
    seen = set()
    for h in hits:
        key = (
            h["orientation"],
            h["start_1b"],
            h["end_1b"],
            h["matched_seq"]
        )
        if key not in seen:
            seen.add(key)
            unique_hits.append(h)

    unique_hits = sorted(
        unique_hits,
        key=lambda x: (x["mismatches"], x["start_1b"])
    )

    return unique_hits


def find_multi_sites_for_primer(
    seq: str,
    primer: Primer,
    max_mismatches: int = 2,
):
    seq = sanitize_seq(seq)

    if primer.variants and len(primer.variants) > 0:
        f_vars = [sanitize_seq(v) for v in primer.variants]
    else:
        f_vars = expand_iupac(primer.seq)

    r_vars = [revcomp(v) for v in f_vars]

    f_hits = find_all_matches_in_seq(
        seq, f_vars, "F",
        max_mismatches=max_mismatches
    )
    r_hits = find_all_matches_in_seq(
        seq, r_vars, "R",
        max_mismatches=max_mismatches
    )

    all_hits = sorted(
        f_hits + r_hits,
        key=lambda x: (x["mismatches"], x["start_1b"])
    )

    return all_hits

# -----------------------------
# 4) Streamlit UI
# -----------------------------
st.set_page_config(page_title="Primer Binding Checker (Demo)", layout="wide")
st.title("Primer Binding Checker (Demo)")

with st.sidebar:
    st.header("Primer Binding check")

st.subheader("Input sequences")

if "input_rows" not in st.session_state:
    st.session_state.input_rows = pd.DataFrame(
        [
            {"order_no": "", "sequence": ""},
        ]
    )

input_df_editor = st.data_editor(
    st.session_state.input_rows,
    num_rows="dynamic",
    use_container_width=True,
    key="input_editor",
    column_config={
        "order_no": st.column_config.TextColumn("Order No"),
        "sequence": st.column_config.TextColumn("Sequence"),
    }
)

run = st.button("Run")

if run:
    input_df = input_df_editor.copy()

    input_df["order_no"] = input_df["order_no"].fillna("").astype(str).str.strip()
    input_df["sequence"] = input_df["sequence"].fillna("").astype(str).str.strip()

    input_df = input_df[
        (input_df["order_no"] != "") &
        (input_df["sequence"] != "")
    ].copy()

    if len(input_df) == 0:
        st.error("Please enter at least one valid row with Order No / Sequence.")
        st.stop()

    all_summary_rows = []
    all_detail_rows = []
    grouped_output = {}

    for _, row in input_df.iterrows():
        order_no = row["order_no"]

        parsed = parse_sequence_and_header_info(row["sequence"])
        seq = parsed["sequence"]
        sample_name = parsed["sample_name"]
        read_label = parsed["read_label"]

        if len(seq) == 0:
            continue

        result_summary_rows = []
        result_detail_rows = []

        best_results = [find_best_for_primer(seq, p) for p in PRIMERS]
        best_df = pd.DataFrame([r.__dict__ for r in best_results])
        best_df = best_df[
            ["primer", "call", "orientation", "start_1b", "end_1b", "mismatches", "primer_used", "matched_seq"]
        ].sort_values(
            by=["mismatches", "primer"],
            ascending=[True, True]
        ).reset_index(drop=True)

        for p in PRIMERS:
            hits = find_multi_sites_for_primer(
                seq,
                p,
                max_mismatches=2,
            )

            hit_count = len(hits)
            multi_possible = "Yes" if hit_count >= 2 else "No"

            pos_text = ", ".join(
                [f"{h['orientation']}:{h['start_1b']}-{h['end_1b']}" for h in hits[:10]]
            )

            summary_row = {
                "order_no": order_no,
                "sample_name": sample_name,
                "read_label": read_label,
                "primer": p.name,
                "hit_count": hit_count,
                "multi_possible": multi_possible,
                "sites": pos_text if pos_text else "-"
            }

            result_summary_rows.append(summary_row)
            all_summary_rows.append(summary_row)

            for idx, h in enumerate(hits, start=1):
                detail_row = {
                    "order_no": order_no,
                    "sample_name": sample_name,
                    "read_label": read_label,
                    "primer": p.name,
                    "hit_no": idx,
                    "orientation": h["orientation"],
                    "start_1b": h["start_1b"],
                    "end_1b": h["end_1b"],
                    "mismatches": h["mismatches"],
                    "primer_used": h["primer_used"],
                    "matched_seq": h["matched_seq"],
                }
                result_detail_rows.append(detail_row)
                all_detail_rows.append(detail_row)

        if order_no not in grouped_output:
            grouped_output[order_no] = []

        grouped_output[order_no].append({
            "sample_name": sample_name,
            "read_label": read_label,
            "best_df": best_df,
            "summary_df": pd.DataFrame(result_summary_rows),
            "detail_df": pd.DataFrame(result_detail_rows),
        })

    for order_no, blocks in grouped_output.items():
        st.markdown(f"## Order No: {order_no}")

        order_summary_df = pd.concat([b["summary_df"] for b in blocks], ignore_index=True)
        order_detail_df = pd.concat([b["detail_df"] for b in blocks], ignore_index=True)

        order_summary_df = order_summary_df.sort_values(
            by=["sample_name", "read_label", "primer"]
        ).reset_index(drop=True)

        order_detail_df = order_detail_df.sort_values(
            by=["sample_name", "read_label", "primer", "hit_no"]
        ).reset_index(drop=True)

        order_excel_bytes = build_excel_bytes(order_summary_df, order_detail_df)

        if len(blocks) == 1 and blocks[0]["sample_name"]:
            sample_name = blocks[0]["sample_name"]
            read_label = blocks[0]["read_label"]
            file_name = make_result_filename(order_no, sample_name, read_label, 1)
        else:
            file_name = f"{order_no}.xlsx"

        st.download_button(
            label=f"Download Excel - {order_no}",
            data=order_excel_bytes,
            file_name=file_name,
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            key=f"excel_{order_no}"
        )

        for idx, block in enumerate(blocks, start=1):
            sample_name = block["sample_name"]
            read_label = block["read_label"]
            best_df = block["best_df"]
            summary_df = block["summary_df"]
            detail_df = block["detail_df"]

            if sample_name or read_label:
                st.markdown(f"### Result {idx} - {sample_name} {read_label}".strip())
            else:
                st.markdown(f"### Result {idx}")

            st.markdown("#### Best binding result")
            best_df_show = best_df.copy()
            best_df_show["call"] = best_df_show["call"].apply(color_call_badge_html)
            best_df_show["matched_seq"] = best_df_show.apply(
                lambda row: color_mismatch_html(row["primer_used"], row["matched_seq"]),
                axis=1
            )

            best_html = best_df_show.to_html(escape=False, index=False)

            components.html(
                f"""
                <style>
                    table {{
                        width: 100%;
                        border-collapse: collapse;
                        font-family: Arial, sans-serif;
                        font-size: 14px;
                        color: white;
                        background-color: #0e1117;
                    }}
                    th, td {{
                        border: 1px solid #31333f;
                        padding: 8px 10px;
                        vertical-align: middle;
                        text-align: left;
                        white-space: nowrap;
                    }}
                    th {{
                        background-color: #1c1f2b;
                        font-weight: 600;
                    }}
                </style>
                {best_html}
                """,
                height=320,
                scrolling=True
            )

            st.markdown("#### Multi binding summary")

            summary_df_show = summary_df.sort_values(
                by=["hit_count", "primer"],
                ascending=[False, True]
            ).reset_index(drop=True)

            summary_df_show["multi_possible"] = summary_df_show["multi_possible"].apply(color_yesno_badge_html)

            summary_html = summary_df_show.to_html(escape=False, index=False)

            components.html(
                f"""
                <style>
                    table {{
                        width: 100%;
                        border-collapse: collapse;
                        font-family: Arial, sans-serif;
                        font-size: 14px;
                        color: white;
                        background-color: #0e1117;
                    }}
                    th, td {{
                        border: 1px solid #31333f;
                        padding: 8px 10px;
                        vertical-align: middle;
                        text-align: left;
                        white-space: nowrap;
                    }}
                    th {{
                        background-color: #1c1f2b;
                        font-weight: 600;
                    }}
                </style>
                {summary_html}
                """,
                height=260,
                scrolling=True
            )

            st.markdown("#### Multi binding detail")
            if len(detail_df) > 0:
                detail_df_show = detail_df.copy()
                detail_df_show["matched_seq"] = detail_df_show.apply(
                    lambda row: color_mismatch_html(row["primer_used"], row["matched_seq"]),
                    axis=1
                )

                detail_html = detail_df_show.to_html(escape=False, index=False)

                components.html(
                    f"""
                    <style>
                        table {{
                            width: 100%;
                            border-collapse: collapse;
                            font-family: Arial, sans-serif;
                            font-size: 14px;
                            color: white;
                            background-color: #0e1117;
                        }}
                        th, td {{
                            border: 1px solid #31333f;
                            padding: 8px 10px;
                            vertical-align: middle;
                            text-align: left;
                            white-space: nowrap;
                        }}
                        th {{
                            background-color: #1c1f2b;
                            font-weight: 600;
                        }}
                    </style>
                    {detail_html}
                    """,
                    height=320,
                    scrolling=True
                )
            else:
                st.info("No detailed hits for this result.")
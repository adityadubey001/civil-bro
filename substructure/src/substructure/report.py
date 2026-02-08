"""
PDF report generator for bridge substructure design.

Produces a comprehensive calculation report using reportlab.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, List

from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_RIGHT
from reportlab.platypus import (
    SimpleDocTemplate,
    Paragraph,
    Spacer,
    Table,
    TableStyle,
    PageBreak,
    KeepTogether,
)


# ── colours ──────────────────────────────────────────────────────────────
_DARK_BLUE = colors.HexColor("#003366")
_MID_BLUE = colors.HexColor("#336699")
_HEADER_BG = colors.HexColor("#003366")
_HEADER_FG = colors.white
_ROW_EVEN = colors.HexColor("#F2F2F2")
_ROW_ODD = colors.white
_GREEN = colors.HexColor("#C6EFCE")
_YELLOW = colors.HexColor("#FFEB9C")
_RED = colors.HexColor("#FFC7CE")

PAGE_W, PAGE_H = A4
_USABLE_W = PAGE_W - 40 * mm  # 20 mm each side


# ── styles ───────────────────────────────────────────────────────────────

def _build_styles() -> dict:
    ss = getSampleStyleSheet()
    custom: dict = {"Normal": ss["Normal"]}
    custom["Title"] = ParagraphStyle(
        "Title", parent=ss["Normal"], fontSize=18, leading=22,
        textColor=_DARK_BLUE, alignment=TA_CENTER, spaceAfter=6 * mm,
        fontName="Helvetica-Bold",
    )
    custom["H1"] = ParagraphStyle(
        "H1", parent=ss["Normal"], fontSize=14, leading=18,
        textColor=_DARK_BLUE, spaceBefore=6 * mm, spaceAfter=3 * mm,
        fontName="Helvetica-Bold",
    )
    custom["H2"] = ParagraphStyle(
        "H2", parent=ss["Normal"], fontSize=11, leading=14,
        textColor=_MID_BLUE, spaceBefore=4 * mm, spaceAfter=2 * mm,
        fontName="Helvetica-Bold",
    )
    custom["CoverTitle"] = ParagraphStyle(
        "CoverTitle", parent=ss["Normal"], fontSize=24, leading=30,
        textColor=_DARK_BLUE, alignment=TA_CENTER,
        fontName="Helvetica-Bold",
    )
    custom["CoverSub"] = ParagraphStyle(
        "CoverSub", parent=ss["Normal"], fontSize=14, leading=18,
        textColor=_MID_BLUE, alignment=TA_CENTER,
    )
    custom["Small"] = ParagraphStyle(
        "Small", parent=ss["Normal"], fontSize=8, leading=10,
    )
    return custom


# ── table helper ─────────────────────────────────────────────────────────

def _tbl(data: list[list], col_widths=None, font_size: int = 9) -> Table:
    """Build a styled Table with alternating row colours."""
    if col_widths is None:
        n_cols = max(len(r) for r in data)
        col_widths = [_USABLE_W / n_cols] * n_cols

    t = Table(data, colWidths=col_widths, repeatRows=1)
    style_cmds = [
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, -1), font_size),
        ("BACKGROUND", (0, 0), (-1, 0), _HEADER_BG),
        ("TEXTCOLOR", (0, 0), (-1, 0), _HEADER_FG),
        ("ALIGN", (1, 0), (-1, -1), "RIGHT"),
        ("ALIGN", (0, 0), (0, -1), "LEFT"),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("GRID", (0, 0), (-1, -1), 0.4, colors.grey),
        ("TOPPADDING", (0, 0), (-1, -1), 2),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
    ]
    for i in range(1, len(data)):
        bg = _ROW_EVEN if i % 2 == 0 else _ROW_ODD
        style_cmds.append(("BACKGROUND", (0, i), (-1, i), bg))
    t.setStyle(TableStyle(style_cmds))
    return t


def _util_color(val: float) -> colors.Color:
    if val > 1.0:
        return _RED
    if val > 0.9:
        return _YELLOW
    return _GREEN


def _fnum(val, dp: int = 1, unit: str = "") -> str:
    """Format a number safely."""
    if val is None:
        return "N/A"
    try:
        s = f"{float(val):.{dp}f}"
    except (TypeError, ValueError):
        return str(val)
    return f"{s} {unit}".strip()


def _fpct(val) -> str:
    """Format utilisation as percentage."""
    if val is None:
        return "N/A"
    try:
        return f"{float(val):.1%}"
    except (TypeError, ValueError):
        return str(val)


def _ok(val) -> str:
    try:
        return "OK" if float(val) <= 1.0 else "NOT OK"
    except Exception:
        return "N/A"


def _safe(obj, attr, default=None):
    """Safely get attribute."""
    return getattr(obj, attr, default)


def _avail(val) -> bool:
    """Check if a result value is usable."""
    return val is not None and val != "not yet computed"


# ── page footer ──────────────────────────────────────────────────────────

def _footer(canvas, doc):
    canvas.saveState()
    canvas.setFont("Helvetica", 8)
    canvas.drawString(20 * mm, 12 * mm, "Bridge Substructure Design - IRC 6 & IRC 112")
    canvas.drawRightString(PAGE_W - 20 * mm, 12 * mm, f"Page {doc.page}")
    canvas.restoreState()


# ── public API ───────────────────────────────────────────────────────────

def generate_report(output_path: str, config: dict, results: dict) -> None:
    """
    Generate PDF design report.

    Args:
        output_path: Path to output PDF file.
        config: Parsed input configuration dict.
        results: Dict mapping module names to result dataclass objects.
    """
    styles = _build_styles()
    story: List = []

    _cover_page(story, styles, config)
    story.append(PageBreak())

    _dashboard(story, styles, results)
    story.append(PageBreak())

    _input_summary(story, styles, config)
    story.append(PageBreak())

    _geometry_section(story, styles, config, results)
    story.append(PageBreak())

    _loads_section(story, styles, results)
    story.append(PageBreak())

    _combinations_section(story, styles, results)
    story.append(PageBreak())

    _pier_cap_section(story, styles, results)
    story.append(PageBreak())

    _pier_section(story, styles, results)
    story.append(PageBreak())

    _pile_capacity_section(story, styles, results)
    story.append(PageBreak())

    _pile_design_section(story, styles, results)
    story.append(PageBreak())

    _pilecap_section(story, styles, results)
    story.append(PageBreak())

    _boq_section(story, styles, results)

    doc = SimpleDocTemplate(
        str(output_path),
        pagesize=A4,
        leftMargin=20 * mm,
        rightMargin=20 * mm,
        topMargin=25 * mm,
        bottomMargin=25 * mm,
    )
    doc.build(story, onFirstPage=_footer, onLaterPages=_footer)


# ── 1. Cover page ────────────────────────────────────────────────────────

def _cover_page(story: list, styles: dict, config: dict):
    proj = config.get("project", {})
    story.append(Spacer(1, 60 * mm))
    story.append(Paragraph("BRIDGE SUBSTRUCTURE DESIGN", styles["CoverTitle"]))
    story.append(Spacer(1, 10 * mm))
    story.append(Paragraph(proj.get("name", ""), styles["CoverSub"]))
    story.append(Paragraph(f"Pier: {proj.get('pier_id', '')}", styles["CoverSub"]))
    story.append(Spacer(1, 15 * mm))
    story.append(Paragraph("Per IRC:6-2017 &amp; IRC:112-2020", styles["CoverSub"]))
    story.append(Spacer(1, 20 * mm))
    info_data = [
        ["Designer", proj.get("designer", "")],
        ["Checker", proj.get("checker", "")],
        ["Date", proj.get("date", "")],
    ]
    t = Table(info_data, colWidths=[40 * mm, 80 * mm])
    t.setStyle(TableStyle([
        ("FONTSIZE", (0, 0), (-1, -1), 11),
        ("FONTNAME", (0, 0), (0, -1), "Helvetica-Bold"),
        ("ALIGN", (0, 0), (-1, -1), "CENTER"),
    ]))
    story.append(t)
    story.append(Spacer(1, 30 * mm))
    story.append(Paragraph(
        "Generated by <b>substructure-design v0.1.0</b>",
        ParagraphStyle("footer", parent=styles["Normal"], alignment=TA_CENTER, fontSize=9),
    ))


# ── 2. Dashboard ─────────────────────────────────────────────────────────

def _dashboard(story: list, styles: dict, results: dict):
    story.append(Paragraph("DESIGN SUMMARY DASHBOARD", styles["Title"]))

    checks: list[list] = [["Check", "Demand", "Capacity", "Util.", "Status"]]

    pc = results.get("pier_cap_design")
    pier = results.get("pier_design")
    pcap = results.get("pile_capacity")
    pd = results.get("pile_design")
    pcd = results.get("pilecap_design")

    util_cells: list[tuple[int, float]] = []  # (row, util_value) for colouring

    def _add(label, demand, capacity, util_val):
        row_idx = len(checks)
        checks.append([
            label,
            demand,
            capacity,
            _fpct(util_val),
            _ok(util_val),
        ])
        if util_val is not None:
            try:
                util_cells.append((row_idx, float(util_val)))
            except Exception:
                pass

    if _avail(pc):
        _add("Pier Cap Flexure (pier face)",
             _fnum(pc.Mu_pier, 1, "kN.m"), _fnum(pc.Mu_capacity_pier, 1, "kN.m"),
             pc.flexure_util_pier)
        _add("Pier Cap Flexure (curtailment)",
             _fnum(pc.Mu_curt, 1, "kN.m"), _fnum(pc.Mu_capacity_curt, 1, "kN.m"),
             pc.flexure_util_curt)
        _add("Pier Cap Shear",
             _fnum(pc.VEd_pier, 1, "kN"),
             _fnum(min(pc.VRds, pc.VRd_max) if pc.shear_reinf_reqd else pc.VRdc_pier, 1, "kN"),
             pc.shear_util)
        _add("Pier Cap Crack Width",
             _fnum(pc.crack_width, 3, "mm"), _fnum(pc.crack_width_limit, 2, "mm"),
             pc.crack_width / pc.crack_width_limit if pc.crack_width_limit else None)

    if _avail(pier):
        _add("Pier Biaxial Bending", "-", "-", pier.util_biaxial)
        _add("Pier Slenderness",
             _fnum(pier.lambda_, 1), _fnum(pier.lambda_lim, 1),
             pier.lambda_ / pier.lambda_lim if pier.lambda_lim else None)
        _add("Pier Crack Width",
             _fnum(pier.crack_width, 3, "mm"), "0.30 mm",
             pier.crack_width / 0.3 if pier.crack_width else None)

    if _avail(pcap):
        _add("Pile Compression",
             _fnum(pcap.max_compression_pile.P, 1, "kN"),
             _fnum(pcap.pile_capacity_compression, 1, "kN"),
             pcap.compression_util)
        _add("Pile Uplift",
             _fnum(pcap.max_uplift_pile.P, 1, "kN"),
             _fnum(pcap.pile_capacity_uplift, 1, "kN"),
             pcap.uplift_util)

    if _avail(pd):
        _add("Pile P-M Interaction", "-", "-", pd.util_PM)
        _add("Pile Shear",
             _fnum(pd.VEd, 1, "kN"), _fnum(pd.VRdc, 1, "kN"),
             pd.shear_util)
        _add("Pile Crack Width",
             _fnum(pd.crack_width, 3, "mm"), _fnum(pd.crack_width_limit, 2, "mm"),
             pd.crack_width / pd.crack_width_limit if pd.crack_width_limit else None)

    if _avail(pcd):
        _add("Pilecap Flexure (Long)",
             _fnum(pcd.Mu_long, 1, "kN.m"), _fnum(pcd.Mu_cap_long, 1, "kN.m"),
             pcd.flexure_util_long)
        _add("Pilecap Flexure (Trans)",
             _fnum(pcd.Mu_trans, 1, "kN.m"), _fnum(pcd.Mu_cap_trans, 1, "kN.m"),
             pcd.flexure_util_trans)
        _add("Pilecap Punching (Pier)",
             _fnum(pcd.VEd_punch_pier, 1, "kN"), _fnum(pcd.VRdc_punch_pier, 1, "kN"),
             pcd.punch_util_pier)
        _add("Pilecap Punching (Pile)",
             _fnum(pcd.VEd_punch_pile, 1, "kN"), _fnum(pcd.VRdc_punch_pile, 1, "kN"),
             pcd.punch_util_pile)
        _add("Pilecap Crack Width",
             _fnum(pcd.crack_width, 3, "mm"), _fnum(pcd.crack_width_limit, 2, "mm"),
             pcd.crack_width / pcd.crack_width_limit if pcd.crack_width_limit else None)

    cw = [55 * mm, 30 * mm, 30 * mm, 25 * mm, 25 * mm]
    t = _tbl(checks, cw, font_size=8)

    # colour the util / status columns
    extra = []
    for row_idx, uval in util_cells:
        bg = _util_color(uval)
        extra.append(("BACKGROUND", (3, row_idx), (4, row_idx), bg))
    if extra:
        t.setStyle(TableStyle(extra))

    story.append(t)


# ── 3. Input summary ─────────────────────────────────────────────────────

def _input_summary(story: list, styles: dict, config: dict):
    story.append(Paragraph("INPUT SUMMARY", styles["Title"]))

    proj = config.get("project", {})
    story.append(Paragraph("Project", styles["H2"]))
    story.append(_tbl([
        ["Parameter", "Value"],
        ["Name", proj.get("name", "")],
        ["Pier ID", proj.get("pier_id", "")],
        ["Designer", proj.get("designer", "")],
    ], [60 * mm, _USABLE_W - 60 * mm]))

    sup = config.get("superstructure", {})
    story.append(Paragraph("Superstructure", styles["H2"]))
    story.append(_tbl([
        ["Parameter", "Value"],
        ["Deck width", _fnum(sup.get("deck_width"), 1, "m")],
        ["Type", str(sup.get("type", ""))],
        ["Left span (c/c brg)", _fnum(sup.get("left_span_cc_bearing"), 2, "m")],
        ["Right span (c/c brg)", _fnum(sup.get("right_span_cc_bearing"), 2, "m")],
        ["Num girders", str(sup.get("num_girders", ""))],
        ["Girder spacing", _fnum(sup.get("girder_spacing"), 2, "m")],
        ["Depth (incl slab)", _fnum(sup.get("depth_incl_slab"), 3, "m")],
        ["Wearing coat", _fnum(sup.get("wearing_coat_thickness"), 0, "mm")],
    ], [60 * mm, _USABLE_W - 60 * mm]))

    pier_cfg = config.get("pier", {})
    story.append(Paragraph("Pier", styles["H2"]))
    story.append(_tbl([
        ["Parameter", "Value"],
        ["Type", pier_cfg.get("type", "")],
        ["Diameter / width (bottom)", _fnum(pier_cfg.get("diameter_bottom"), 2, "m")],
        ["Diameter / width (top)", _fnum(pier_cfg.get("diameter_top"), 2, "m")],
    ], [60 * mm, _USABLE_W - 60 * mm]))

    pc_cfg = config.get("pier_cap", {})
    story.append(Paragraph("Pier Cap", styles["H2"]))
    story.append(_tbl([
        ["Parameter", "Value"],
        ["Width (long)", _fnum(pc_cfg.get("width_long"), 2, "m")],
        ["Length (trans)", _fnum(pc_cfg.get("length_trans"), 2, "m")],
        ["Depth max", _fnum(pc_cfg.get("depth_max"), 2, "m")],
        ["Depth min", _fnum(pc_cfg.get("depth_min"), 2, "m")],
    ], [60 * mm, _USABLE_W - 60 * mm]))

    fdn = config.get("foundation", {})
    story.append(Paragraph("Foundation", styles["H2"]))
    story.append(_tbl([
        ["Parameter", "Value"],
        ["Pile diameter", _fnum(fdn.get("pile_diameter"), 2, "m")],
        ["Piles (long x trans)", f"{fdn.get('piles_long', 0)} x {fdn.get('piles_trans', 0)}"],
        ["Spacing factor", str(fdn.get("pile_spacing_factor", ""))],
        ["Pilecap thickness", _fnum(fdn.get("pilecap_thickness"), 2, "m")],
        ["Pile capacity (comp)", _fnum(fdn.get("pile_capacity_top"), 0, "tonnes")],
        ["Pile capacity (uplift)", _fnum(fdn.get("pile_uplift_capacity"), 0, "tonnes")],
    ], [60 * mm, _USABLE_W - 60 * mm]))

    mat = config.get("materials", {})
    story.append(Paragraph("Materials", styles["H2"]))
    story.append(_tbl([
        ["Parameter", "Value"],
        ["Pier fck", _fnum(mat.get("pier", {}).get("fck"), 0, "MPa")],
        ["Pier cap fck", _fnum(mat.get("pier_cap", {}).get("fck"), 0, "MPa")],
        ["Pile / pilecap fck", _fnum(mat.get("pile_pilecap", {}).get("fck"), 0, "MPa")],
        ["Steel fyk", _fnum(mat.get("steel", {}).get("fyk"), 0, "MPa")],
        ["Exposure", mat.get("exposure", "")],
    ], [60 * mm, _USABLE_W - 60 * mm]))

    seis = config.get("seismic", {})
    story.append(Paragraph("Seismic", styles["H2"]))
    story.append(_tbl([
        ["Parameter", "Value"],
        ["Zone", str(seis.get("zone", ""))],
        ["Importance factor", str(seis.get("importance_factor", ""))],
        ["Response reduction R", str(seis.get("response_reduction", ""))],
        ["Soil type", str(seis.get("soil_type", ""))],
    ], [60 * mm, _USABLE_W - 60 * mm]))


# ── 4. Geometry ──────────────────────────────────────────────────────────

def _geometry_section(story: list, styles: dict, config: dict, results: dict):
    story.append(Paragraph("GEOMETRY RESULTS", styles["Title"]))

    geom = results.get("geometry")
    if not _avail(geom):
        story.append(Paragraph("Geometry results not available.", styles["Normal"]))
        return

    story.append(Paragraph("Levels", styles["H2"]))
    lvl = config.get("levels", {})
    story.append(_tbl([
        ["Level", "Elevation (m)"],
        ["FRL", _fnum(lvl.get("frl"), 3)],
        ["GL", _fnum(lvl.get("gl"), 3)],
        ["Pilecap top", _fnum(geom.pilecap_top_level, 3)],
        ["Pilecap bottom", _fnum(geom.pilecap_bottom_level, 3)],
        ["Pier cap top", _fnum(geom.piercap_top_level, 3)],
        ["Pier cap bottom", _fnum(geom.piercap_bottom_level, 3)],
    ], [60 * mm, _USABLE_W - 60 * mm]))

    story.append(Paragraph("Pier Section", styles["H2"]))
    ps = geom.pier_section
    story.append(_tbl([
        ["Property", "Value"],
        ["Type", ps.pier_type],
        ["Width long", _fnum(ps.width_long, 3, "m")],
        ["Width trans", _fnum(ps.width_trans, 3, "m")],
        ["Area", _fnum(ps.area, 4, "m2")],
        ["Ixx", _fnum(ps.inertia_xx, 4, "m4")],
        ["Iyy", _fnum(ps.inertia_yy, 4, "m4")],
        ["Pier height", _fnum(geom.pier_height, 3, "m")],
    ], [60 * mm, _USABLE_W - 60 * mm]))

    story.append(Paragraph("Pile Arrangement", styles["H2"]))
    ph = [["Pile", "Row", "Col", "X (m)", "Y (m)"]]
    for i, pc in enumerate(geom.pile_coordinates):
        ph.append([str(i), str(pc.row), str(pc.col),
                    _fnum(pc.x, 3), _fnum(pc.y, 3)])
    story.append(_tbl(ph, [20 * mm, 20 * mm, 20 * mm, 40 * mm, 40 * mm]))

    story.append(Paragraph("Pile Cap Dimensions", styles["H2"]))
    story.append(_tbl([
        ["Parameter", "Value"],
        ["Length (long)", _fnum(geom.pilecap_length_long, 3, "m")],
        ["Width (trans)", _fnum(geom.pilecap_width_trans, 3, "m")],
        ["Thickness", _fnum(geom.pilecap_thickness, 3, "m")],
    ], [60 * mm, _USABLE_W - 60 * mm]))


# ── 5. Loads ─────────────────────────────────────────────────────────────

def _fv_row(label: str, fv) -> list:
    return [label, _fnum(fv.P, 1), _fnum(fv.FL, 1), _fnum(fv.FT, 1),
            _fnum(fv.ML, 1), _fnum(fv.MT, 1)]

_FV_HDR = ["Load", "P (kN)", "FL (kN)", "FT (kN)", "ML (kN.m)", "MT (kN.m)"]
_FV_CW = [40 * mm, 24 * mm, 24 * mm, 24 * mm, 28 * mm, 28 * mm]


def _loads_section(story: list, styles: dict, results: dict):
    story.append(Paragraph("LOADS SUMMARY", styles["Title"]))

    loads = results.get("loads")
    if not _avail(loads):
        story.append(Paragraph("Load results not available.", styles["Normal"]))
        return

    # Dead loads
    story.append(Paragraph("Dead Loads", styles["H2"]))
    dl_rows = [_FV_HDR]
    dl_rows.append(_fv_row("DL super LHS", loads.dl_super_lhs))
    dl_rows.append(_fv_row("DL super RHS", loads.dl_super_rhs))
    dl_rows.append(_fv_row("DL pier cap", loads.dl_piercap))
    dl_rows.append(_fv_row("DL pier", loads.dl_pier))
    dl_rows.append(_fv_row("Total DL", loads.total_dl_at_pier_base))
    story.append(_tbl(dl_rows, _FV_CW))

    # SIDL
    story.append(Paragraph("SIDL", styles["H2"]))
    sidl_rows = [_FV_HDR]
    sidl_rows.append(_fv_row("SIDL LHS", loads.sidl_lhs))
    sidl_rows.append(_fv_row("SIDL RHS", loads.sidl_rhs))
    sidl_rows.append(_fv_row("WC LHS", loads.wc_lhs))
    sidl_rows.append(_fv_row("WC RHS", loads.wc_rhs))
    story.append(_tbl(sidl_rows, _FV_CW))

    # Live load cases
    story.append(Paragraph("Live Load Cases (at pier base)", styles["H2"]))
    ll_rows = [_FV_HDR]
    for llc in loads.ll_cases:
        ll_rows.append(_fv_row(llc.name, llc.force_at_pier_base))
    story.append(_tbl(ll_rows, _FV_CW))

    # Wind
    story.append(Paragraph("Wind Loads", styles["H2"]))
    w_rows = [_FV_HDR]
    w_rows.append(_fv_row("Wind super LHS", loads.wind_on_super_lhs))
    w_rows.append(_fv_row("Wind super RHS", loads.wind_on_super_rhs))
    w_rows.append(_fv_row("Wind on pier", loads.wind_on_pier))
    w_rows.append(_fv_row("Wind on pier cap", loads.wind_on_piercap))
    story.append(_tbl(w_rows, _FV_CW))

    # Seismic
    story.append(Paragraph("Seismic Loads", styles["H2"]))
    story.append(Paragraph(
        f"Ah_long = {_fnum(loads.Ah_long, 4)} &nbsp;&nbsp; Ah_trans = {_fnum(loads.Ah_trans, 4)}",
        styles["Normal"],
    ))
    s_rows = [_FV_HDR]
    s_rows.append(_fv_row("Seismic super (L)", loads.seismic_super_long))
    s_rows.append(_fv_row("Seismic super (T)", loads.seismic_super_trans))
    s_rows.append(_fv_row("Seismic pier (L)", loads.seismic_pier_long))
    s_rows.append(_fv_row("Seismic pier (T)", loads.seismic_pier_trans))
    story.append(_tbl(s_rows, _FV_CW))

    # Temperature
    story.append(Paragraph("Temperature", styles["H2"]))
    t_rows = [_FV_HDR, _fv_row("Temp (long)", loads.temp_long)]
    story.append(_tbl(t_rows, _FV_CW))


# ── 6. Load combinations ─────────────────────────────────────────────────

def _combinations_section(story: list, styles: dict, results: dict):
    story.append(Paragraph("LOAD COMBINATIONS", styles["Title"]))

    combos = results.get("combinations")
    if not _avail(combos):
        story.append(Paragraph("Combination results not available.", styles["Normal"]))
        return

    story.append(Paragraph(
        f"Total combinations: {len(combos.all_combinations)}", styles["Normal"]))

    story.append(Paragraph("Governing Combinations (at pier base)", styles["H2"]))
    hdr = ["Case", "Name", "P (kN)", "ML (kN.m)", "MT (kN.m)"]
    cw = [25 * mm, 55 * mm, 30 * mm, 30 * mm, 30 * mm]
    rows = [hdr]
    for label, combo in [
        ("Max P", combos.max_P),
        ("Min P", combos.min_P),
        ("Max ML", combos.max_ML),
        ("Max MT", combos.max_MT),
    ]:
        fv = combo.forces_pier_base
        rows.append([label, combo.name[:30],
                     _fnum(fv.P, 1), _fnum(fv.ML, 1), _fnum(fv.MT, 1)])
    story.append(_tbl(rows, cw))

    # ULS Basic governing
    story.append(Paragraph("ULS Basic Governing", styles["H2"]))
    gub = combos.governing_uls_basic
    fv = gub.forces_pier_base
    story.append(_tbl([
        ["Parameter", "Value"],
        ["Name", gub.name],
        ["P", _fnum(fv.P, 1, "kN")],
        ["FL", _fnum(fv.FL, 1, "kN")],
        ["FT", _fnum(fv.FT, 1, "kN")],
        ["ML", _fnum(fv.ML, 1, "kN.m")],
        ["MT", _fnum(fv.MT, 1, "kN.m")],
    ], [50 * mm, _USABLE_W - 50 * mm]))

    # ULS Seismic governing
    story.append(Paragraph("ULS Seismic Governing", styles["H2"]))
    gus = combos.governing_uls_seismic
    fv = gus.forces_pier_base
    story.append(_tbl([
        ["Parameter", "Value"],
        ["Name", gus.name],
        ["P", _fnum(fv.P, 1, "kN")],
        ["FL", _fnum(fv.FL, 1, "kN")],
        ["FT", _fnum(fv.FT, 1, "kN")],
        ["ML", _fnum(fv.ML, 1, "kN.m")],
        ["MT", _fnum(fv.MT, 1, "kN.m")],
    ], [50 * mm, _USABLE_W - 50 * mm]))


# ── 7. Pier cap design ───────────────────────────────────────────────────

def _pier_cap_section(story: list, styles: dict, results: dict):
    story.append(Paragraph("PIER CAP DESIGN", styles["Title"]))

    cap = results.get("pier_cap_design")
    if not _avail(cap):
        story.append(Paragraph("Results not available.", styles["Normal"]))
        return

    story.append(Paragraph("Section", styles["H2"]))
    story.append(_tbl([
        ["Parameter", "Value"],
        ["Width", _fnum(cap.width, 0, "mm")],
        ["Depth at pier", _fnum(cap.depth_at_pier, 0, "mm")],
        ["Depth at curtailment", _fnum(cap.depth_at_curtailment, 0, "mm")],
        ["Cover", _fnum(cap.cover, 0, "mm")],
        ["d_eff (pier)", _fnum(cap.d_eff_pier, 0, "mm")],
        ["d_eff (curt)", _fnum(cap.d_eff_curt, 0, "mm")],
    ], [60 * mm, _USABLE_W - 60 * mm]))

    story.append(Paragraph("Flexure - Pier Face", styles["H2"]))
    story.append(_tbl([
        ["Parameter", "Value"],
        ["Mu (design)", _fnum(cap.Mu_pier, 1, "kN.m")],
        ["Mu (capacity)", _fnum(cap.Mu_capacity_pier, 1, "kN.m")],
        ["xu/d", _fnum(cap.xu_d_pier, 4)],
        ["Ast required", _fnum(cap.Ast_reqd_pier, 0, "mm2")],
        ["Ast provided", _fnum(cap.Ast_provided_pier, 0, "mm2")],
        ["Utilisation", _fpct(cap.flexure_util_pier)],
    ], [60 * mm, _USABLE_W - 60 * mm]))

    story.append(Paragraph("Flexure - Curtailment", styles["H2"]))
    story.append(_tbl([
        ["Parameter", "Value"],
        ["Mu (design)", _fnum(cap.Mu_curt, 1, "kN.m")],
        ["Mu (capacity)", _fnum(cap.Mu_capacity_curt, 1, "kN.m")],
        ["xu/d", _fnum(cap.xu_d_curt, 4)],
        ["Ast required", _fnum(cap.Ast_reqd_curt, 0, "mm2")],
        ["Ast provided", _fnum(cap.Ast_provided_curt, 0, "mm2")],
        ["Utilisation", _fpct(cap.flexure_util_curt)],
    ], [60 * mm, _USABLE_W - 60 * mm]))

    story.append(Paragraph("Shear", styles["H2"]))
    story.append(_tbl([
        ["Parameter", "Value"],
        ["VEd (pier)", _fnum(cap.VEd_pier, 1, "kN")],
        ["VRd,c (pier)", _fnum(cap.VRdc_pier, 1, "kN")],
        ["Shear reinf required?", "Yes" if cap.shear_reinf_reqd else "No"],
        ["VRd,s", _fnum(cap.VRds, 1, "kN")],
        ["VRd,max", _fnum(cap.VRd_max, 1, "kN")],
        ["Shear util", _fpct(cap.shear_util)],
    ], [60 * mm, _USABLE_W - 60 * mm]))

    story.append(Paragraph("Crack Width", styles["H2"]))
    story.append(_tbl([
        ["Parameter", "Value"],
        ["wk", _fnum(cap.crack_width, 4, "mm")],
        ["Limit", _fnum(cap.crack_width_limit, 2, "mm")],
        ["Status", "OK" if cap.crack_ok else "NOT OK"],
    ], [60 * mm, _USABLE_W - 60 * mm]))

    story.append(Spacer(1, 4 * mm))
    story.append(Paragraph(f"Overall Status: {cap.status}", styles["H1"]))


# ── 8. Pier design ───────────────────────────────────────────────────────

def _pier_section(story: list, styles: dict, results: dict):
    story.append(Paragraph("PIER DESIGN", styles["Title"]))

    pier = results.get("pier_design")
    if not _avail(pier):
        story.append(Paragraph("Results not available.", styles["Normal"]))
        return

    story.append(Paragraph("Section Properties", styles["H2"]))
    story.append(_tbl([
        ["Parameter", "Value"],
        ["Type", pier.pier_type],
        ["Diameter", _fnum(pier.diameter, 0, "mm")],
        ["fck", _fnum(pier.fck, 0, "MPa")],
        ["Cover", _fnum(pier.cover, 0, "mm")],
        ["n_bars", str(pier.n_bars)],
        ["Bar dia", _fnum(pier.bar_dia, 0, "mm")],
        ["Ast provided", _fnum(pier.Ast_provided, 0, "mm2")],
        ["rho_l", _fnum(pier.rho_l * 100, 3, "%")],
    ], [60 * mm, _USABLE_W - 60 * mm]))

    story.append(Paragraph("Slenderness", styles["H2"]))
    story.append(_tbl([
        ["Parameter", "Value"],
        ["lambda", _fnum(pier.lambda_, 1)],
        ["lambda_lim", _fnum(pier.lambda_lim, 1)],
        ["Second order required?", "Yes" if pier.second_order else "No"],
    ], [60 * mm, _USABLE_W - 60 * mm]))

    story.append(Paragraph("ULS Check", styles["H2"]))
    story.append(_tbl([
        ["Parameter", "Value"],
        ["Governing combo", pier.governing_combo_name[:50]],
        ["P applied", _fnum(pier.P_applied, 1, "kN")],
        ["ML applied", _fnum(pier.ML_applied, 1, "kN.m")],
        ["MT applied", _fnum(pier.MT_applied, 1, "kN.m")],
        ["ML (incl 2nd order)", _fnum(pier.ML_with_2nd_order, 1, "kN.m")],
        ["MT (incl 2nd order)", _fnum(pier.MT_with_2nd_order, 1, "kN.m")],
        ["M capacity (XX)", _fnum(pier.M_capacity_xx, 1, "kN.m")],
        ["M capacity (YY)", _fnum(pier.M_capacity_yy, 1, "kN.m")],
        ["Uniaxial util XX", _fpct(pier.util_uniaxial_xx)],
        ["Uniaxial util YY", _fpct(pier.util_uniaxial_yy)],
        ["Biaxial util", _fpct(pier.util_biaxial)],
    ], [60 * mm, _USABLE_W - 60 * mm]))

    story.append(Paragraph("SLS Check", styles["H2"]))
    story.append(_tbl([
        ["Parameter", "Value"],
        ["sigma_c (rare)", _fnum(pier.sigma_c_rare, 1, "MPa")],
        ["sigma_c (QP)", _fnum(pier.sigma_c_qp, 1, "MPa")],
        ["sigma_s (rare)", _fnum(pier.sigma_s_rare, 1, "MPa")],
        ["Crack width", _fnum(pier.crack_width, 4, "mm")],
        ["SLS OK?", "Yes" if pier.sls_ok else "No"],
    ], [60 * mm, _USABLE_W - 60 * mm]))

    story.append(Paragraph("Ductile Detailing", styles["H2"]))
    story.append(_tbl([
        ["Parameter", "Value"],
        ["Plastic hinge Lp", _fnum(pier.Lp, 0, "mm")],
        ["Spiral dia", _fnum(pier.spiral_dia, 0, "mm")],
        ["Spiral spacing", _fnum(pier.spiral_spacing, 0, "mm")],
        ["omega_wd", _fnum(pier.omega_wd, 4)],
        ["omega_wd required", _fnum(pier.omega_wd_required, 4)],
        ["Ductile OK?", "Yes" if pier.ductile_ok else "No"],
    ], [60 * mm, _USABLE_W - 60 * mm]))

    story.append(Spacer(1, 4 * mm))
    story.append(Paragraph(f"Overall Status: {pier.status}", styles["H1"]))


# ── 9. Pile capacity ─────────────────────────────────────────────────────

def _pile_capacity_section(story: list, styles: dict, results: dict):
    story.append(Paragraph("PILE CAPACITY", styles["Title"]))

    pcap = results.get("pile_capacity")
    if not _avail(pcap):
        story.append(Paragraph("Results not available.", styles["Normal"]))
        return

    story.append(Paragraph("Pile Configuration", styles["H2"]))
    story.append(_tbl([
        ["Parameter", "Value"],
        ["Number of piles", str(pcap.num_piles)],
        ["Pile diameter", _fnum(pcap.pile_diameter, 3, "m")],
        ["Sum x^2", _fnum(pcap.sum_x2, 3, "m2")],
        ["Sum y^2", _fnum(pcap.sum_y2, 3, "m2")],
        ["Compression capacity", _fnum(pcap.pile_capacity_compression, 1, "kN")],
        ["Uplift capacity", _fnum(pcap.pile_capacity_uplift, 1, "kN")],
    ], [60 * mm, _USABLE_W - 60 * mm]))

    story.append(Paragraph("Max Compression Pile", styles["H2"]))
    mp = pcap.max_compression_pile
    story.append(_tbl([
        ["Parameter", "Value"],
        ["Combination", pcap.max_compression_combo[:50]],
        ["Pile index", str(mp.pile_index)],
        ["Location", f"x={_fnum(mp.pile_coord.x, 3)} m, y={_fnum(mp.pile_coord.y, 3)} m"],
        ["P", _fnum(mp.P, 1, "kN")],
        ["Utilisation", _fpct(pcap.compression_util)],
    ], [60 * mm, _USABLE_W - 60 * mm]))

    story.append(Paragraph("Max Uplift Pile", styles["H2"]))
    up = pcap.max_uplift_pile
    story.append(_tbl([
        ["Parameter", "Value"],
        ["Combination", pcap.max_uplift_combo[:50]],
        ["Pile index", str(up.pile_index)],
        ["Location", f"x={_fnum(up.pile_coord.x, 3)} m, y={_fnum(up.pile_coord.y, 3)} m"],
        ["P", _fnum(up.P, 1, "kN")],
        ["Utilisation", _fpct(pcap.uplift_util)],
    ], [60 * mm, _USABLE_W - 60 * mm]))

    # Summary of all combos — show first 20
    story.append(Paragraph("Combination Summary (first 20)", styles["H2"]))
    hdr = ["Combination", "Category", "Max P (kN)", "Min P (kN)"]
    cw = [55 * mm, 30 * mm, 40 * mm, 40 * mm]
    rows = [hdr]
    for pr in pcap.all_pile_results[:20]:
        rows.append([pr.combo_name[:30], pr.category,
                     _fnum(pr.max_pile_P, 1), _fnum(pr.min_pile_P, 1)])
    story.append(_tbl(rows, cw, font_size=7))

    story.append(Spacer(1, 4 * mm))
    story.append(Paragraph(f"Overall Status: {pcap.status}", styles["H1"]))


# ── 10. Pile design ──────────────────────────────────────────────────────

def _pile_design_section(story: list, styles: dict, results: dict):
    story.append(Paragraph("PILE STRUCTURAL DESIGN", styles["Title"]))

    pd = results.get("pile_design")
    if not _avail(pd):
        story.append(Paragraph("Results not available.", styles["Normal"]))
        return

    story.append(Paragraph("Section", styles["H2"]))
    story.append(_tbl([
        ["Parameter", "Value"],
        ["Diameter", _fnum(pd.diameter, 0, "mm")],
        ["Cover", _fnum(pd.cover, 0, "mm")],
        ["d_eff", _fnum(pd.d_eff, 1, "mm")],
        ["fck", _fnum(pd.fck, 0, "MPa")],
        ["n_bars", str(pd.n_bars)],
        ["Bar dia", _fnum(pd.bar_dia, 0, "mm")],
        ["Ast provided", _fnum(pd.Ast_provided, 0, "mm2")],
        ["rho_l", _fnum(pd.rho_l * 100, 3, "%")],
    ], [60 * mm, _USABLE_W - 60 * mm]))

    story.append(Paragraph("Applied Loads (governing combo)", styles["H2"]))
    story.append(_tbl([
        ["Parameter", "Value"],
        ["Governing combo", pd.governing_combo_name[:50]],
        ["P applied", _fnum(pd.P_applied, 1, "kN")],
        ["H applied", _fnum(pd.H_applied, 1, "kN")],
        ["M applied", _fnum(pd.M_applied, 1, "kN.m")],
    ], [60 * mm, _USABLE_W - 60 * mm]))

    story.append(Paragraph("Design Checks", styles["H2"]))
    story.append(_tbl([
        ["Check", "Value"],
        ["P-M interaction util", _fpct(pd.util_PM)],
        ["VEd", _fnum(pd.VEd, 1, "kN")],
        ["VRd,c", _fnum(pd.VRdc, 1, "kN")],
        ["Shear util", _fpct(pd.shear_util)],
        ["Crack width", _fnum(pd.crack_width, 4, "mm")],
        ["Crack limit", _fnum(pd.crack_width_limit, 2, "mm")],
        ["Crack OK?", "Yes" if pd.crack_ok else "No"],
        ["Curtailment depth", _fnum(pd.curtailment_depth, 2, "m")],
    ], [60 * mm, _USABLE_W - 60 * mm]))

    story.append(Spacer(1, 4 * mm))
    story.append(Paragraph(f"Overall Status: {pd.status}", styles["H1"]))


# ── 11. Pilecap design ───────────────────────────────────────────────────

def _pilecap_section(story: list, styles: dict, results: dict):
    story.append(Paragraph("PILE CAP DESIGN", styles["Title"]))

    pcd = results.get("pilecap_design")
    if not _avail(pcd):
        story.append(Paragraph("Results not available.", styles["Normal"]))
        return

    story.append(Paragraph("Dimensions", styles["H2"]))
    story.append(_tbl([
        ["Parameter", "Value"],
        ["Length (long)", _fnum(pcd.length_long, 0, "mm")],
        ["Width (trans)", _fnum(pcd.width_trans, 0, "mm")],
        ["Thickness", _fnum(pcd.thickness, 0, "mm")],
        ["d_eff", _fnum(pcd.d_eff, 1, "mm")],
        ["fck", _fnum(pcd.fck, 0, "MPa")],
        ["Cover", _fnum(pcd.cover, 0, "mm")],
    ], [60 * mm, _USABLE_W - 60 * mm]))

    story.append(Paragraph("Flexure - Longitudinal", styles["H2"]))
    story.append(_tbl([
        ["Parameter", "Value"],
        ["Mu", _fnum(pcd.Mu_long, 1, "kN.m")],
        ["Mu capacity", _fnum(pcd.Mu_cap_long, 1, "kN.m")],
        ["Ast required", _fnum(pcd.Ast_reqd_long, 0, "mm2")],
        ["Ast provided", _fnum(pcd.Ast_prov_long, 0, "mm2")],
        ["Utilisation", _fpct(pcd.flexure_util_long)],
    ], [60 * mm, _USABLE_W - 60 * mm]))

    story.append(Paragraph("Flexure - Transverse", styles["H2"]))
    story.append(_tbl([
        ["Parameter", "Value"],
        ["Mu", _fnum(pcd.Mu_trans, 1, "kN.m")],
        ["Mu capacity", _fnum(pcd.Mu_cap_trans, 1, "kN.m")],
        ["Ast required", _fnum(pcd.Ast_reqd_trans, 0, "mm2")],
        ["Ast provided", _fnum(pcd.Ast_prov_trans, 0, "mm2")],
        ["Utilisation", _fpct(pcd.flexure_util_trans)],
    ], [60 * mm, _USABLE_W - 60 * mm]))

    story.append(Paragraph("One-Way Shear", styles["H2"]))
    story.append(_tbl([
        ["Direction", "VEd (kN)", "VRd,c (kN)", "Util"],
        ["Longitudinal", _fnum(pcd.VEd_long, 1), _fnum(pcd.VRdc_long, 1),
         _fpct(pcd.shear_util_long)],
        ["Transverse", _fnum(pcd.VEd_trans, 1), _fnum(pcd.VRdc_trans, 1),
         _fpct(pcd.shear_util_trans)],
    ], [40 * mm, 35 * mm, 35 * mm, 35 * mm]))

    story.append(Paragraph("Punching Shear", styles["H2"]))
    story.append(_tbl([
        ["Location", "VEd (kN)", "VRd,c (kN)", "u (mm)", "Util"],
        ["Pier", _fnum(pcd.VEd_punch_pier, 1), _fnum(pcd.VRdc_punch_pier, 1),
         _fnum(pcd.u_pier, 0), _fpct(pcd.punch_util_pier)],
        ["Pile", _fnum(pcd.VEd_punch_pile, 1), _fnum(pcd.VRdc_punch_pile, 1),
         _fnum(pcd.u_pile, 0), _fpct(pcd.punch_util_pile)],
    ], [30 * mm, 30 * mm, 30 * mm, 30 * mm, 30 * mm]))

    story.append(Paragraph("Crack Width", styles["H2"]))
    story.append(_tbl([
        ["Parameter", "Value"],
        ["wk", _fnum(pcd.crack_width, 4, "mm")],
        ["Limit", _fnum(pcd.crack_width_limit, 2, "mm")],
        ["OK?", "Yes" if pcd.crack_ok else "No"],
    ], [60 * mm, _USABLE_W - 60 * mm]))

    story.append(Spacer(1, 4 * mm))
    story.append(Paragraph(f"Overall Status: {pcd.status}", styles["H1"]))


# ── 12. BOQ ──────────────────────────────────────────────────────────────

def _boq_section(story: list, styles: dict, results: dict):
    story.append(Paragraph("BILL OF QUANTITIES", styles["Title"]))

    boq = results.get("boq")
    if not _avail(boq):
        story.append(Paragraph("BOQ not available.", styles["Normal"]))
        return

    hdr = ["Element", "Concrete (m3)", "Steel (kg)", "Conc. Cost (INR)",
           "Steel Cost (INR)", "Total (INR)"]
    cw = [28 * mm, 26 * mm, 24 * mm, 28 * mm, 28 * mm, 28 * mm]
    rows = [hdr]
    for item in boq.items:
        rows.append([
            item.element,
            _fnum(item.concrete_volume, 2),
            _fnum(item.steel_weight, 1),
            f"{item.concrete_cost:,.0f}",
            f"{item.steel_cost:,.0f}",
            f"{item.total_cost:,.0f}",
        ])
    rows.append([
        "TOTAL",
        _fnum(boq.total_concrete, 2),
        _fnum(boq.total_steel, 1),
        "", "",
        f"{boq.total_cost:,.0f}",
    ])
    story.append(_tbl(rows, cw))

    story.append(Spacer(1, 5 * mm))
    story.append(Paragraph(
        f"Total Concrete: {boq.total_concrete:.2f} m3 &nbsp;&nbsp;|&nbsp;&nbsp; "
        f"Total Steel: {boq.total_steel:.1f} kg ({boq.total_steel / 1000:.2f} tonnes) &nbsp;&nbsp;|&nbsp;&nbsp; "
        f"Total Cost: INR {boq.total_cost:,.0f}",
        styles["Normal"],
    ))

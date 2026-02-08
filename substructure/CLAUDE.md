# CLAUDE.md — Substructure Design Streamlit App

## What This Project Is

A Streamlit web app for bridge substructure design per IRC 6:2017 (loads) and IRC 112:2020 (concrete design). It wraps a CLI-based design engine (~12,000 lines) with an interactive UI.

The engine was originally built as a CLI tool (source: `/Users/adityadubey/Documents/substructure-design-handoff/`). The engine modules were copied into `src/substructure/` and a new Streamlit UI layer was built on top at `src/ui/`.

There is a **separate** beam/girder design app at `/Users/adityadubey/Documents/beam-design/` — same author, same design code (IRC 112), but for superstructure (girders) not substructure (piers/piles). The substructure app's UI patterns were modeled after that beam app.

## How to Run

```bash
cd /Users/adityadubey/Documents/substructure-app
streamlit run app.py
```

Dependencies: `pip install streamlit pyyaml numpy click reportlab matplotlib`

## Architecture

```
app.py                    # Streamlit entry point
src/
├── ui/                   # NEW — UI layer (1,981 lines total)
│   ├── state.py          # Session state defaults + build_config_from_state()
│   ├── runner.py         # Engine wrapper — calls all 11 design steps
│   ├── yaml_mode.py      # YAML upload/paste input mode
│   ├── wizard.py         # 8-step guided input wizard
│   └── results.py        # Dashboard (6 metric cards) + 6 detail tabs
├── substructure/         # COPIED from handoff — the design engine (~12K lines)
│   ├── loads.py          # 1908 lines — DL, SIDL, LL, wind, seismic, braking, temp
│   ├── load_combinations.py  # IRC 6 Tables B.1-B.4, 105 combinations
│   ├── pier_cap_design.py    # Flexure, shear, crack width (tapered cantilever)
│   ├── pier_design.py        # Biaxial P-M interaction, slenderness, ductile detailing
│   ├── interaction.py        # P-M diagram generation (circular/rectangular)
│   ├── pile_capacity.py      # Pile group distribution, IRC 78 capacity
│   ├── pile_design.py        # Pile structural design (IS:2911 moment method)
│   ├── pilecap_design.py     # Flexure, punching shear (4 perimeters)
│   ├── geometry.py           # Derived dimensions, pile grid, section props
│   ├── materials.py          # IRC 112 concrete/steel properties
│   ├── creep.py              # IRC 112 Annex A-2 creep coefficient
│   ├── input_parser.py       # YAML validation
│   ├── boq.py                # Bill of quantities + cost
│   ├── report.py             # PDF generation (not used by Streamlit app yet)
│   └── utils.py              # Section properties, unit converters, IRC table cache
config/
├── sample_input.yaml     # Complete reference input file
└── irc_tables.yaml       # IRC code tables (load factors, material props, etc.)
```

## Design Pipeline (runner.py)

```
YAML config dict
    → geometry.calculate_geometry()
    → materials.get_concrete_properties() per component
    → creep.calculate_creep()
    → loads.calculate_loads()
    → load_combinations.calculate_combinations()   # 105 combos
    → pier_cap_design.design_pier_cap()
    → pier_design.design_pier()
    → pile_capacity.calculate_pile_capacity()
    → pile_design.design_pile()
    → pilecap_design.design_pilecap()
    → boq.calculate_boq()
    → returns dict of result dataclasses
```

## Two Input Modes

1. **YAML Upload** (default) — upload file, paste text, or load the bundled sample. Has validate + run buttons.
2. **Wizard** — 8 steps: Project/Levels → Superstructure → Pier → Pier Cap/Bearings → Foundation → Materials/Creep → Seismic/Wind → Review/Run

## Results Display

- **Dashboard**: 6 metric cards (Pier Cap, Pier, Pile Capacity, Pile Design, Pile Cap, Cost) with OK/FAIL delta
- **6 Detail Tabs**: Pier Cap, Pier, Piles, Pile Cap, Loads, BOQ — each with sub-metrics and tables

## Current Status — All 5 Checks Pass

| Check | Utilization | Status |
|---|---|---|
| Pier Cap | flex 0.803, shear 0.928, crack 0.265mm | OK |
| Pier | biaxial 0.884 | OK |
| Pile Capacity | compression 0.991 | OK |
| Pile Design | P-M 0.602 | OK |
| Pile Cap | pier punch 0.758, pile punch 0.995 | OK |
| Cost | INR 15.8L | — |

**Tests**: Run `pytest tests/` — 15 validation tests covering all design checks.

## Critical Conventions (READ THESE)

### Engine Conventions
- **Units**: Geometry in SI (metres). Design modules convert to mm internally.
- **Status strings**: `"OK"` / `"NOT OK"` — always these exact strings.
- **Attribute naming inconsistency**: `CombinationResult.name` (not `.combo_name`), but `PileCombinationResult.combo_name` (not `.name`).
- **Steel properties**: `materials.get_steel_properties(fyk, Es)` — `Es` is in **GPa** (not MPa). Default 200.0.
- **Load combination categories**: `'uls_basic'`, `'uls_seismic'`, `'sls_rare'`, `'sls_frequent'`, `'sls_quasi'`
- **Use `python3`** not `python`. Needs `PYTHONPATH=src` for direct module execution.
- **irc_tables.yaml path**: `utils.load_irc_tables()` searches relative to `utils.py` — resolves to `config/irc_tables.yaml`. Do not move `config/`.

### UI Conventions
- Session state keys are flat strings (not nested dicts). See `_DEFAULTS` in `state.py` for the full list.
- `build_config_from_state()` assembles the nested YAML-compatible config dict from flat state keys.
- `runner.run_design(config)` returns a dict with string keys. Each value is the raw engine dataclass (or None if that step failed).
- `safe_attr(obj, attr, default)` is used throughout results.py for null-safe attribute access.
- Results page uses `st.metric()` with `delta="OK"/"FAIL"` and `delta_color="normal"/"inverse"` for color coding.

## Known Gaps / Future Work

1. **PDF report not wired up** — `src/substructure/report.py` exists and works (CLI uses it), but the Streamlit app doesn't expose it yet. Next step: add a download button on the results page that calls `report.generate_report()`.
2. ~~**No tests**~~ — **FIXED**: `tests/` directory with 15 pytest validation tests.
3. ~~**Pier cap VEd 25% lower than Excel**~~ — **FIXED**: Added optional `reactions_uls` and `reactions_sls` in bearings section. Users can now provide per-bearing reactions from FE analysis.
4. **Pier biaxial util 0.884 vs Excel 0.658** — hand-calc wind forces are larger than Excel's SOFiSTiK FE-computed forces. Not fixable without FE analysis.
5. **Cost 12% high** — conservative steel quantities from iterative crack width control, punching reinforcement, auto-sized spiral.
6. ~~**fywd for shear design**~~ — **FIXED**: Now uses `fywd = 0.8 * min(fyk, 500)` = 400 MPa per IRC 112. Shear util: 0.776 → 0.928.
7. **Modules not implemented** — bearing_design, STM, fatigue, construction stage, ancillary (pedestal/AD block).
8. **Wizard bearing coordinates** — currently fixed at 4 bearings per side. Could add +/- buttons to dynamically add/remove bearings.

## Deployment

**GitHub**: https://github.com/adityadubey001/civil-bro
**Streamlit Cloud**: Deployed at `substructure/app.py` (separate from beam design app at root `app.py`)

To deploy updates:
1. Push changes to the `civil-bro` repo
2. Streamlit Cloud auto-redeploys from `main` branch

## Key Decisions Made

1. **Hosted in civil-bro repo** — substructure app lives in `substructure/` subfolder alongside the beam design app.
2. **Engine modified for IRC compliance** — fywd shear fix applied to `pier_cap_design.py`.
3. **No Pydantic models** — the engine uses plain dataclasses and dicts (not Pydantic like the beam app). The config is a plain dict matching the YAML structure.
4. **PDF report deferred** — focus was on interactive UI first. The CLI's report.py can be wired up later.
5. **Two input modes** — YAML upload for power users, wizard for guided input. Both produce the same config dict and use the same runner.

## Related Projects

| Project | Path | What |
|---|---|---|
| Beam design app | `/Users/adityadubey/Documents/beam-design/` | IRC 112 girder design (Streamlit) |
| Substructure handoff | `/Users/adityadubey/Documents/substructure-design-handoff/` | Original CLI tool + HANDOFF.md with full technical details |
| Reference Excel | `beam-design/reference-files/23B002-SDN-SUP-01-R2.xlsm` | Excel workbook the engine replicates |

## Full Technical Reference

For detailed formulas, five major fixes, exact IRC code references, and complete Excel comparison, read:
`/Users/adityadubey/Documents/substructure-design-handoff/HANDOFF.md`

---

## Session Summary (2026-02-08)

### Completed This Session
1. **Deployed to Streamlit Cloud** via civil-bro repo (`substructure/app.py`)
2. **Fixed fywd shear design** — now uses 400 MPa per IRC 112 (was 478 MPa). Shear util: 0.776 → 0.928
3. **Added validation tests** — 15 pytest tests in `tests/test_design_validation.py`
4. **Added per-bearing reactions** — optional `reactions_uls`/`reactions_sls` in YAML bearings section

### Remaining Work
| # | Item | Notes |
|---|------|-------|
| 1 | PDF report | Wire up download button to `report.generate_report()` |
| 5 | Cost 12% high | Conservative steel quantities |
| 7 | Missing modules | bearing_design, STM, fatigue, etc. |
| 8 | Wizard dynamic bearings | Add +/- buttons for bearing count |

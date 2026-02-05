# RC Beam Design Application

A Q&A-style web application for designing simply supported reinforced concrete beams per IS 456:2000.

## Features

- **Q&A Flow**: Step-by-step guided input collection
- **IS 456:2000 Compliance**: Design calculations per Indian Standard
- **Auto-sizing**: Preliminary dimensions from span/depth ratios
- **Singly/Doubly Reinforced**: User choice for reinforcement type
- **Visual Output**: Cross-section diagrams with Matplotlib
- **PDF Reports**: Downloadable design reports with ReportLab

## Quick Start

### 1. Install Dependencies

```bash
cd /Users/adityadubey/Documents/beam-design
pip install -r requirements.txt
```

### 2. Run the Application

```bash
streamlit run app.py
```

The app will open in your browser at `http://localhost:8501`

## Design Workflow

1. **Geometry**: Enter span length (optional: custom dimensions)
2. **Materials**: Select concrete and steel grades
3. **Loads**: Enter dead load and live load
4. **Design Options**: Choose singly/doubly reinforced preference
5. **Exposure**: Select exposure condition for cover requirements
6. **Review**: Verify all inputs
7. **Results**: View design output, diagrams, and download PDF

## Project Structure

```
beam-design/
├── app.py                 # Main Streamlit application
├── requirements.txt       # Python dependencies
├── src/
│   ├── core/              # Calculation engine
│   │   ├── beam_design.py # Main orchestrator
│   │   ├── flexure.py     # Flexural design
│   │   ├── shear.py       # Shear design
│   │   └── serviceability.py
│   ├── codes/             # Code provisions
│   │   └── is456.py       # IS 456:2000
│   ├── models/            # Pydantic data models
│   │   ├── inputs.py
│   │   └── outputs.py
│   ├── ui/                # Streamlit UI
│   │   └── session_state.py
│   └── reports/           # Report generation
│       ├── pdf_generator.py
│       └── diagrams/
│           └── cross_section.py
└── reference-files/       # Reference documents
```

## Design Approach

### Flexural Design (IS 456, Clause 38)
- Limit state method with rectangular stress block
- xu_max/d = 0.46 for Fe500 steel
- Minimum reinforcement per Clause 26.5.1.1

### Shear Design (IS 456, Clause 40)
- Design shear strength τc from Table 19
- Maximum shear stress τc,max from Table 20
- Stirrup spacing from Clause 40.4

### Deflection Check (IS 456, Clause 23.2.1)
- Span/depth ratio method
- Modification factors for tension and compression steel

## Example

For a 6m span beam with 10 kN/m dead load and 5 kN/m live load:

- **Section**: 200 × 350 mm
- **Main Steel**: 3-25φ (Ast = 1473 mm²)
- **Stirrups**: 2L-8φ @ 225mm c/c
- **All checks pass**

## Code References

- IS 456:2000 - Plain and Reinforced Concrete
- IS 875 (Part 1-5) - Loading Standards

## Notes

- Self-weight is calculated automatically from section dimensions
- Cover is determined from exposure condition (Table 16)
- Design iterates automatically if initial section is inadequate

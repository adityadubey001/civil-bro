"""
PDF report generator for IRC 112 Bridge Girder Design using ReportLab.
"""

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import mm, cm
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
    PageBreak, Image, KeepTogether
)
from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_RIGHT
from datetime import datetime
from typing import Optional
import io

from src.models.inputs import BridgeGirderInput
from src.models.outputs import BridgeGirderOutput, DesignStatus
from src.reports.diagrams.cross_section import generate_cross_section


class PDFReportGenerator:
    """
    Generate PDF design reports for IRC 112 bridge girder design.
    """

    def __init__(self):
        self.styles = getSampleStyleSheet()
        self._setup_custom_styles()

    def _setup_custom_styles(self):
        """Set up custom paragraph styles."""
        self.styles.add(ParagraphStyle(
            name='TitleStyle',
            parent=self.styles['Title'],
            fontSize=18,
            spaceAfter=30,
            textColor=colors.HexColor('#2c3e50')
        ))

        self.styles.add(ParagraphStyle(
            name='SectionHeader',
            parent=self.styles['Heading2'],
            fontSize=14,
            spaceBefore=20,
            spaceAfter=10,
            textColor=colors.HexColor('#34495e'),
            borderWidth=1,
            borderColor=colors.HexColor('#3498db'),
            borderPadding=5
        ))

        self.styles.add(ParagraphStyle(
            name='SubSection',
            parent=self.styles['Heading3'],
            fontSize=11,
            spaceBefore=10,
            spaceAfter=5,
            textColor=colors.HexColor('#7f8c8d')
        ))

        self.styles.add(ParagraphStyle(
            name='BodyTextSmall',
            parent=self.styles['BodyText'],
            fontSize=9
        ))

        self.styles.add(ParagraphStyle(
            name='CodeRef',
            parent=self.styles['BodyText'],
            fontSize=8,
            textColor=colors.HexColor('#7f8c8d'),
            fontName='Helvetica-Oblique'
        ))

    def generate_report(
        self,
        inputs: BridgeGirderInput,
        outputs: BridgeGirderOutput,
        project_name: str = "IRC 112 Bridge Girder Design",
        engineer_name: Optional[str] = None,
        include_diagrams: bool = True,
    ) -> bytes:
        """
        Generate complete PDF report.

        Args:
            inputs: Design input parameters
            outputs: Design calculation results
            project_name: Project name for report header
            engineer_name: Engineer name (optional)
            include_diagrams: Whether to include cross-section diagram

        Returns:
            PDF file content as bytes
        """
        buffer = io.BytesIO()
        doc = SimpleDocTemplate(
            buffer,
            pagesize=A4,
            rightMargin=2*cm,
            leftMargin=2*cm,
            topMargin=2*cm,
            bottomMargin=2*cm
        )

        # Build story (content)
        story = []

        # Title page
        story.extend(self._build_title_page(project_name, outputs, engineer_name))

        # Input summary
        story.extend(self._build_input_summary(inputs, outputs))

        # Load analysis
        story.extend(self._build_load_analysis_section(outputs))

        # Flexural design details
        story.extend(self._build_flexure_section(outputs))

        # Shear design details
        story.extend(self._build_shear_section(outputs))

        # SLS Stress checks
        story.extend(self._build_sls_stress_section(outputs))

        # Crack width check
        if outputs.crack_width:
            story.extend(self._build_crack_width_section(outputs))

        # Deflection check
        story.extend(self._build_deflection_section(outputs))

        # Reinforcement summary
        story.extend(self._build_reinforcement_summary(outputs))

        # Cross-section diagram
        if include_diagrams:
            story.extend(self._build_diagram_section(outputs))

        # Build PDF
        doc.build(story)

        buffer.seek(0)
        return buffer.getvalue()

    def _build_title_page(self, project_name, outputs, engineer_name):
        """Build title page elements."""
        elements = []

        # Title
        elements.append(Paragraph(
            "IRC 112 BRIDGE GIRDER DESIGN REPORT",
            self.styles['TitleStyle']
        ))

        # Project name
        elements.append(Paragraph(
            f"<b>{project_name}</b>",
            self.styles['Heading2']
        ))

        elements.append(Spacer(1, 20))

        # Status badge
        status = outputs.overall_status.value.upper()
        status_color = '#27ae60' if status == 'PASS' else '#f39c12' if status == 'WARNING' else '#e74c3c'
        elements.append(Paragraph(
            f'<font color="{status_color}"><b>DESIGN STATUS: {status}</b></font>',
            self.styles['Heading3']
        ))

        elements.append(Spacer(1, 30))

        # Quick summary table
        summary_data = [
            ['Section', f'{outputs.girder_width:.0f} × {outputs.girder_depth:.0f} mm'],
            ['Main Steel', outputs.flexure.tension_bar_arrangement],
            ['Stirrups', f'{outputs.shear.stirrup_legs}L-{outputs.shear.stirrup_diameter}φ @ {outputs.shear.spacing_provided:.0f}mm c/c'],
            ['Design Code', outputs.design_code],
        ]

        summary_table = Table(summary_data, colWidths=[5*cm, 8*cm])
        summary_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (0, -1), colors.HexColor('#ecf0f1')),
            ('FONTNAME', (0, 0), (0, -1), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, -1), 10),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.gray),
            ('TOPPADDING', (0, 0), (-1, -1), 8),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 8),
        ]))
        elements.append(summary_table)

        elements.append(Spacer(1, 30))

        # Report metadata
        meta_data = [
            ['Report Date', datetime.now().strftime('%Y-%m-%d %H:%M')],
            ['Engineer', engineer_name or 'Not specified'],
        ]
        meta_table = Table(meta_data, colWidths=[3*cm, 6*cm])
        meta_table.setStyle(TableStyle([
            ('FONTSIZE', (0, 0), (-1, -1), 9),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
        ]))
        elements.append(meta_table)

        elements.append(PageBreak())

        return elements

    def _build_input_summary(self, inputs, outputs):
        """Build input summary section."""
        elements = []

        elements.append(Paragraph("1. INPUT DATA", self.styles['SectionHeader']))

        # Geometry
        elements.append(Paragraph("1.1 Geometry", self.styles['SubSection']))
        geo_data = [
            ['Parameter', 'Value', 'Unit'],
            ['Effective Span', f'{outputs.span_length:.1f}', 'm'],
            ['Girder Type', outputs.girder_type, '-'],
            ['Web Width', f'{outputs.girder_width:.0f}', 'mm'],
            ['Overall Depth', f'{outputs.girder_depth:.0f}', 'mm'],
            ['Effective Depth', f'{outputs.effective_depth:.0f}', 'mm'],
            ['Clear Cover', f'{outputs.clear_cover:.0f}', 'mm'],
        ]
        if outputs.flange_width:
            geo_data.append(['Flange Width', f'{outputs.flange_width:.0f}', 'mm'])
            geo_data.append(['Flange Depth', f'{outputs.flange_depth:.0f}', 'mm'])

        elements.append(self._create_data_table(geo_data))
        elements.append(Spacer(1, 10))

        # Deck parameters
        elements.append(Paragraph("1.2 Deck Parameters", self.styles['SubSection']))
        deck_data = [
            ['Parameter', 'Value', 'Unit'],
            ['Deck Width per Girder', f'{inputs.deck.width:.2f}', 'm'],
            ['Deck Thickness', f'{inputs.deck.thickness:.0f}', 'mm'],
            ['Wearing Coat', f'{inputs.deck.wearing_coat:.0f}', 'mm'],
        ]
        elements.append(self._create_data_table(deck_data))
        elements.append(Spacer(1, 10))

        # Materials
        elements.append(Paragraph("1.3 Materials", self.styles['SubSection']))
        mat_data = [
            ['Parameter', 'Value', 'Unit'],
            ['Concrete Grade', outputs.concrete_grade, '-'],
            ['Steel Grade', outputs.steel_grade, '-'],
            ['fck', f'{outputs.fck:.0f}', 'MPa'],
            ['fy', f'{outputs.fy:.0f}', 'MPa'],
            ['fcd (0.446×fck)', f'{outputs.fcd:.2f}', 'MPa'],
            ['fyd (fy/1.15)', f'{outputs.fyd:.2f}', 'MPa'],
            ['Ecm', f'{outputs.ecm:.0f}', 'MPa'],
        ]
        elements.append(self._create_data_table(mat_data))
        elements.append(Spacer(1, 10))

        # Exposure
        elements.append(Paragraph("1.4 Exposure Condition", self.styles['SubSection']))
        exp_val = inputs.exposure_condition
        if hasattr(exp_val, 'value'):
            exp_val = exp_val.value
        elements.append(Paragraph(
            f"Exposure: {exp_val.replace('_', ' ').title()}",
            self.styles['BodyText']
        ))

        return elements

    def _build_load_analysis_section(self, outputs):
        """Build load analysis section with IRC 6 results."""
        elements = []

        elements.append(Paragraph("2. LOAD ANALYSIS", self.styles['SectionHeader']))
        elements.append(Paragraph("Reference: IRC 6:2017", self.styles['CodeRef']))

        loads = outputs.loads

        # Dead loads
        elements.append(Paragraph("2.1 Dead Loads", self.styles['SubSection']))
        dl_data = [
            ['Component', 'Value', 'Unit'],
            ['Girder Self-Weight', f'{loads.sw_girder:.2f}', 'kN/m'],
            ['Deck Slab', f'{loads.sw_deck:.2f}', 'kN/m'],
            ['Wearing Coat', f'{loads.sw_wearing_coat:.2f}', 'kN/m'],
            ['Total Dead Load', f'{loads.total_dead_load:.2f}', 'kN/m'],
        ]
        elements.append(self._create_data_table(dl_data))
        elements.append(Spacer(1, 10))

        # Dead load effects
        elements.append(Paragraph("2.2 Dead Load Effects", self.styles['SubSection']))
        dle_data = [
            ['Effect', 'Value', 'Unit'],
            ['Bending Moment (Mdl)', f'{loads.dead_load_moment:.2f}', 'kNm'],
            ['Shear Force (Vdl)', f'{loads.dead_load_shear:.2f}', 'kN'],
        ]
        elements.append(self._create_data_table(dle_data))
        elements.append(Spacer(1, 10))

        # IRC 6 Live load analysis
        elements.append(Paragraph("2.3 Live Load Analysis (IRC 6)", self.styles['SubSection']))
        ll_header = ['Vehicle Type', 'BM (kNm)', 'SF (kN)', 'Impact', 'BM+IF (kNm)', 'SF+IF (kN)', 'Governs']
        ll_data = [ll_header]

        for ll_result in loads.live_load_results:
            governs = ""
            if ll_result.is_governing_bm:
                governs += "BM "
            if ll_result.is_governing_sf:
                governs += "SF"
            ll_data.append([
                ll_result.load_type,
                f'{ll_result.max_bm_without_impact:.0f}',
                f'{ll_result.max_sf_without_impact:.0f}',
                f'{ll_result.impact_factor:.2f}',
                f'{ll_result.max_bm_with_impact:.0f}',
                f'{ll_result.max_sf_with_impact:.0f}',
                governs.strip() or '-'
            ])

        ll_table = Table(ll_data, colWidths=[3*cm, 2*cm, 2*cm, 1.5*cm, 2.2*cm, 2.2*cm, 2*cm])
        ll_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#3498db')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, -1), 8),
            ('ALIGN', (1, 0), (-1, -1), 'CENTER'),
            ('ALIGN', (0, 0), (0, -1), 'LEFT'),
            ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.gray),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#f8f9fa')]),
            ('TOPPADDING', (0, 0), (-1, -1), 4),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 4),
        ]))
        elements.append(ll_table)
        elements.append(Spacer(1, 10))

        # Governing cases
        elements.append(Paragraph("2.4 Governing Cases & Load Combinations", self.styles['SubSection']))
        gov_data = [
            ['Item', 'Value', 'Unit'],
            ['Governing BM Case', loads.governing_bm_case, '-'],
            ['Governing SF Case', loads.governing_sf_case, '-'],
            ['Max Live Load Moment', f'{loads.governing_live_moment:.2f}', 'kNm'],
            ['Max Live Load Shear', f'{loads.governing_live_shear:.2f}', 'kN'],
        ]
        elements.append(self._create_data_table(gov_data))
        elements.append(Spacer(1, 10))

        # Factored loads (ULS)
        elements.append(Paragraph("2.5 ULS Load Combination (IRC 112 Table B.1)", self.styles['SubSection']))
        elements.append(Paragraph(
            "Mu = 1.35×Mdl + 1.50×Mll | Vu = 1.35×Vdl + 1.50×Vll",
            self.styles['BodyText']
        ))
        uls_data = [
            ['Design Force', 'Calculation', 'Value', 'Unit'],
            ['Factored Moment (Mu)',
             f'1.35×{loads.dead_load_moment:.0f} + 1.50×{loads.governing_live_moment:.0f}',
             f'{loads.factored_moment:.2f}', 'kNm'],
            ['Factored Shear (Vu)',
             f'1.35×{loads.dead_load_shear:.0f} + 1.50×{loads.governing_live_shear:.0f}',
             f'{loads.factored_shear:.2f}', 'kN'],
            ['Service Moment (Ms)',
             f'{loads.dead_load_moment:.0f} + {loads.governing_live_moment:.0f}',
             f'{loads.service_moment:.2f}', 'kNm'],
        ]
        uls_table = Table(uls_data, colWidths=[4*cm, 5.5*cm, 2.5*cm, 1.5*cm])
        uls_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#3498db')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, -1), 9),
            ('ALIGN', (1, 0), (-1, -1), 'CENTER'),
            ('ALIGN', (0, 0), (0, -1), 'LEFT'),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.gray),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#f8f9fa')]),
        ]))
        elements.append(uls_table)

        return elements

    def _build_flexure_section(self, outputs):
        """Build flexural design section per IRC 112 Cl. 10.3."""
        elements = []

        elements.append(Paragraph("3. FLEXURAL DESIGN", self.styles['SectionHeader']))
        elements.append(Paragraph("Reference: IRC 112:2020 Cl. 10.3", self.styles['CodeRef']))

        flexure = outputs.flexure

        # Design forces
        elements.append(Paragraph("3.1 Design Moment", self.styles['SubSection']))
        moment_data = [
            ['Parameter', 'Value', 'Unit'],
            ['Design Moment (Mu)', f'{flexure.design_moment:.2f}', 'kNm'],
            ['Moment Capacity Limit (Mu,lim)', f'{flexure.moment_capacity_limit:.2f}', 'kNm'],
            ['Moment Capacity (provided)', f'{flexure.moment_capacity:.2f}', 'kNm'],
            ['Utilization Ratio', f'{flexure.utilization_ratio:.2%}', '-'],
        ]
        elements.append(self._create_data_table(moment_data))
        elements.append(Spacer(1, 10))

        # Section analysis
        elements.append(Paragraph("3.2 Section Analysis", self.styles['SubSection']))
        section_data = [
            ['Parameter', 'Value', 'Unit'],
            ['Neutral Axis Depth (xu)', f'{flexure.neutral_axis_depth:.1f}', 'mm'],
            ['xu/d ratio', f'{flexure.xu_d_ratio:.3f}', '-'],
            ['xu,max/d limit', '0.456', '-'],
            ['Lever Arm (z)', f'{flexure.lever_arm:.1f}', 'mm'],
        ]
        elements.append(self._create_data_table(section_data))
        elements.append(Spacer(1, 10))

        # Tension reinforcement
        elements.append(Paragraph("3.3 Tension Reinforcement", self.styles['SubSection']))
        reinf_data = [
            ['Parameter', 'Value', 'Unit'],
            ['Required Ast', f'{flexure.required_ast:.0f}', 'mm²'],
            ['Provided Ast', f'{flexure.provided_ast:.0f}', 'mm²'],
            ['Bar Arrangement', flexure.tension_bar_arrangement, '-'],
            ['Steel %', f'{flexure.pt_provided:.2f}', '%'],
        ]
        elements.append(self._create_data_table(reinf_data))
        elements.append(Spacer(1, 10))

        # Compression reinforcement if doubly reinforced
        if flexure.is_doubly_reinforced:
            elements.append(Paragraph("3.4 Compression Reinforcement", self.styles['SubSection']))
            comp_data = [
                ['Parameter', 'Value', 'Unit'],
                ['Required Asc', f'{flexure.required_asc:.0f}', 'mm²'],
                ['Provided Asc', f'{flexure.provided_asc:.0f}', 'mm²'],
                ['Bar Arrangement', flexure.compression_bar_arrangement or '-', '-'],
            ]
            elements.append(self._create_data_table(comp_data))
            elements.append(Spacer(1, 10))

        # Code checks
        elements.append(Paragraph("3.5 Code Checks", self.styles['SubSection']))
        check_data = [
            ['Check', 'Status'],
            ['Under-reinforced (xu/d ≤ 0.456)', 'PASS' if flexure.is_under_reinforced else 'FAIL'],
            ['Minimum reinforcement (Cl. 16.5.1)', 'PASS' if flexure.min_ast_satisfied else 'FAIL'],
            ['Maximum reinforcement (Cl. 16.5.1)', 'PASS' if flexure.max_ast_satisfied else 'FAIL'],
        ]
        elements.append(self._create_check_table(check_data))

        return elements

    def _build_shear_section(self, outputs):
        """Build shear design section per IRC 112 Cl. 10.4."""
        elements = []

        elements.append(Paragraph("4. SHEAR DESIGN", self.styles['SectionHeader']))
        elements.append(Paragraph("Reference: IRC 112:2020 Cl. 10.3.2 & 10.3.3", self.styles['CodeRef']))

        shear = outputs.shear

        # Design forces
        elements.append(Paragraph("4.1 Design Shear", self.styles['SubSection']))
        shear_data = [
            ['Parameter', 'Value', 'Unit'],
            ['Design Shear (Vu)', f'{shear.design_shear:.2f}', 'kN'],
            ['Shear at d from support', f'{shear.shear_at_d_from_support:.2f}', 'kN'],
        ]
        elements.append(self._create_data_table(shear_data))
        elements.append(Spacer(1, 10))

        # Shear capacity
        elements.append(Paragraph("4.2 Shear Capacity (IRC 112 Cl. 10.3.2)", self.styles['SubSection']))
        cap_data = [
            ['Parameter', 'Value', 'Unit'],
            ['Size effect factor (k)', f'{shear.size_effect_factor_k:.3f}', '-'],
            ['Reinforcement ratio (ρl)', f'{shear.reinforcement_ratio_rho:.4f}', '-'],
            ['Concrete shear capacity (VRd,c)', f'{shear.concrete_shear_capacity:.2f}', 'kN'],
            ['Min concrete capacity (VRd,c,min)', f'{shear.min_concrete_shear_capacity:.2f}', 'kN'],
            ['Max shear capacity (VRd,max)', f'{shear.max_shear_capacity:.2f}', 'kN'],
        ]
        elements.append(self._create_data_table(cap_data))
        elements.append(Spacer(1, 10))

        # Stirrups
        elements.append(Paragraph("4.3 Shear Reinforcement", self.styles['SubSection']))
        stirrup_data = [
            ['Parameter', 'Value', 'Unit'],
            ['Stirrup diameter', f'{shear.stirrup_diameter}', 'mm'],
            ['Number of legs', f'{shear.stirrup_legs}', '-'],
            ['Stirrup area (Asv)', f'{shear.stirrup_area:.0f}', 'mm²'],
            ['Required spacing', f'{shear.spacing_required:.0f}', 'mm'],
            ['Provided spacing', f'{shear.spacing_provided:.0f}', 'mm'],
            ['Max spacing (0.75d)', f'{shear.spacing_max:.0f}', 'mm'],
        ]
        elements.append(self._create_data_table(stirrup_data))
        elements.append(Spacer(1, 10))

        # Code checks
        elements.append(Paragraph("4.4 Code Checks", self.styles['SubSection']))
        check_data = [
            ['Check', 'Status'],
            ['Shear adequacy (Vu ≤ VRd,max)', 'PASS' if shear.shear_adequacy_check else 'FAIL'],
            ['Minimum shear reinforcement', 'PASS' if shear.min_shear_reinf_check else 'FAIL'],
        ]
        elements.append(self._create_check_table(check_data))

        return elements

    def _build_sls_stress_section(self, outputs):
        """Build SLS stress check section per IRC 112 Cl. 12.2.1."""
        elements = []

        elements.append(Paragraph("5. SLS STRESS CHECKS", self.styles['SectionHeader']))
        elements.append(Paragraph("Reference: IRC 112:2020 Cl. 12.2.1", self.styles['CodeRef']))

        sls = outputs.sls_stress

        # Elastic analysis
        elements.append(Paragraph("5.1 Cracked Section Analysis", self.styles['SubSection']))
        elastic_data = [
            ['Parameter', 'Value', 'Unit'],
            ['Service Moment (Ms)', f'{sls.rare_combination_moment:.2f}', 'kNm'],
            ['Modular Ratio (m = Es/Ecm)', f'{sls.modular_ratio:.2f}', '-'],
            ['Elastic NA Depth (x)', f'{sls.elastic_na_depth:.1f}', 'mm'],
            ['Cracked MI (Icr)', f'{sls.cracked_moment_of_inertia:.2e}', 'mm⁴'],
        ]
        elements.append(self._create_data_table(elastic_data))
        elements.append(Spacer(1, 10))

        # Rare combination stress check
        elements.append(Paragraph("5.2 Rare Combination Stresses", self.styles['SubSection']))
        elements.append(Paragraph("Limit: σc ≤ 0.48×fck, σs ≤ 0.8×fy", self.styles['CodeRef']))
        rare_data = [
            ['Stress', 'Calculated', 'Limit', 'Utilization', 'Status'],
            ['Concrete (σc)', f'{sls.concrete_stress:.2f} MPa', f'{sls.concrete_stress_limit:.2f} MPa',
             f'{sls.concrete_utilization:.1%}', 'PASS' if sls.concrete_stress_check else 'FAIL'],
            ['Steel (σs)', f'{sls.steel_stress:.2f} MPa', f'{sls.steel_stress_limit:.2f} MPa',
             f'{sls.steel_utilization:.1%}', 'PASS' if sls.steel_stress_check else 'FAIL'],
        ]
        rare_table = Table(rare_data, colWidths=[3*cm, 3*cm, 3*cm, 2.5*cm, 2*cm])
        rare_table.setStyle(self._get_stress_table_style(rare_data))
        elements.append(rare_table)
        elements.append(Spacer(1, 10))

        # Quasi-permanent combination (to avoid non-linear creep)
        if sls.qp_moment is not None:
            elements.append(Paragraph("5.3 Quasi-Permanent Combination", self.styles['SubSection']))
            elements.append(Paragraph(
                "Limit: σc ≤ 0.36×fck (to avoid non-linear creep per Cl. 12.2.1(3))",
                self.styles['CodeRef']
            ))
            elements.append(Paragraph(
                "Note: ψ2 = 0 for traffic loads per IRC 112 Table B.3, so Mqp = Mdl only",
                self.styles['BodyTextSmall']
            ))
            qp_data = [
                ['Parameter', 'Value', 'Unit'],
                ['Quasi-Permanent Moment (Mqp)', f'{sls.qp_moment:.2f}', 'kNm'],
                ['Concrete Stress (σc,qp)', f'{sls.qp_concrete_stress:.2f}', 'MPa'],
                ['Stress Limit (0.36×fck)', f'{sls.qp_concrete_stress_limit:.2f}', 'MPa'],
                ['Check', 'PASS' if sls.qp_concrete_stress_check else 'FAIL', '-'],
            ]
            elements.append(self._create_data_table(qp_data))

        return elements

    def _build_crack_width_section(self, outputs):
        """Build crack width check section per IRC 112 Cl. 12.3.4."""
        elements = []

        elements.append(Paragraph("6. CRACK WIDTH CHECK", self.styles['SectionHeader']))
        elements.append(Paragraph("Reference: IRC 112:2020 Cl. 12.3.4", self.styles['CodeRef']))

        cw = outputs.crack_width

        # Parameters
        elements.append(Paragraph("6.1 Crack Spacing Parameters", self.styles['SubSection']))
        param_data = [
            ['Parameter', 'Value', 'Unit'],
            ['Steel stress under QP (σs)', f'{cw.steel_stress_qp:.2f}', 'MPa'],
            ['Clear cover (c)', f'{cw.clear_cover:.0f}', 'mm'],
            ['Bar diameter (φ)', f'{cw.bar_diameter:.0f}', 'mm'],
            ['Bar spacing (s)', f'{cw.bar_spacing:.0f}', 'mm'],
            ['Effective reinf. ratio (ρp,eff)', f'{cw.effective_reinforcement_ratio:.4f}', '-'],
        ]
        elements.append(self._create_data_table(param_data))
        elements.append(Spacer(1, 10))

        # Crack width calculation
        elements.append(Paragraph("6.2 Crack Width Calculation", self.styles['SubSection']))
        elements.append(Paragraph(
            "wk = sr,max × (εsm - εcm) per Eq. 12.5",
            self.styles['CodeRef']
        ))
        calc_data = [
            ['Parameter', 'Value', 'Unit'],
            ['Max crack spacing (sr,max)', f'{cw.max_crack_spacing:.1f}', 'mm'],
            ['Strain difference (εsm - εcm)', f'{cw.strain_difference:.6f}', '-'],
            ['Calculated crack width (wk)', f'{cw.calculated_crack_width:.3f}', 'mm'],
            ['Allowable crack width', f'{cw.allowable_crack_width:.2f}', 'mm'],
        ]
        elements.append(self._create_data_table(calc_data))
        elements.append(Spacer(1, 10))

        # Check result
        status = 'PASS' if cw.crack_width_check else 'FAIL'
        status_color = '#27ae60' if status == 'PASS' else '#e74c3c'
        elements.append(Paragraph(
            f'<font color="{status_color}"><b>CRACK WIDTH CHECK: {status}</b></font>',
            self.styles['Heading3']
        ))

        return elements

    def _build_deflection_section(self, outputs):
        """Build deflection check section per IRC 112 Cl. 12.4."""
        elements = []

        elements.append(Paragraph("7. DEFLECTION CHECK", self.styles['SectionHeader']))
        elements.append(Paragraph("Reference: IRC 112:2020 Cl. 12.4", self.styles['CodeRef']))

        defl = outputs.deflection

        # Section properties
        elements.append(Paragraph("7.1 Section Properties", self.styles['SubSection']))
        prop_data = [
            ['Parameter', 'Value', 'Unit'],
            ['Gross MI (Ig)', f'{defl.gross_moment_of_inertia:.2e}', 'mm⁴'],
            ['Centroid from bottom', f'{defl.centroid_from_bottom:.1f}', 'mm'],
            ['Cracking Moment (Mcr)', f'{defl.cracking_moment:.2f}', 'kNm'],
            ['Effective MI (Ie)', f'{defl.effective_moment_of_inertia:.2e}', 'mm⁴'],
        ]
        elements.append(self._create_data_table(prop_data))
        elements.append(Spacer(1, 10))

        # Deflection results
        elements.append(Paragraph("7.2 Deflection Results", self.styles['SubSection']))
        defl_data = [
            ['Parameter', 'Value', 'Unit'],
            ['Calculated Deflection', f'{defl.calculated_deflection:.2f}', 'mm'],
            ['Allowable (Total, L/250)', f'{defl.allowable_deflection_total:.2f}', 'mm'],
            ['Allowable (Live, L/800)', f'{defl.allowable_deflection_live:.2f}', 'mm'],
        ]
        elements.append(self._create_data_table(defl_data))
        elements.append(Spacer(1, 10))

        # Check result
        status = 'PASS' if defl.deflection_check_passed else 'FAIL'
        status_color = '#27ae60' if status == 'PASS' else '#e74c3c'
        elements.append(Paragraph(
            f'<font color="{status_color}"><b>DEFLECTION CHECK: {status}</b></font>',
            self.styles['Heading3']
        ))

        return elements

    def _build_reinforcement_summary(self, outputs):
        """Build final reinforcement summary section."""
        elements = []

        elements.append(Paragraph("8. REINFORCEMENT SUMMARY", self.styles['SectionHeader']))

        # Main reinforcement table
        summary_data = [
            ['Reinforcement Type', 'Details'],
            ['Bottom Steel (Tension)', outputs.flexure.tension_bar_arrangement],
        ]

        if outputs.flexure.is_doubly_reinforced:
            summary_data.append(['Top Steel (Compression)', outputs.flexure.compression_bar_arrangement or '-'])

        summary_data.append([
            'Stirrups',
            f'{outputs.shear.stirrup_legs}L-{outputs.shear.stirrup_diameter}φ @ {outputs.shear.spacing_provided:.0f}mm c/c'
        ])

        if outputs.side_face_reinforcement.required:
            summary_data.append(['Side Face', outputs.side_face_reinforcement.arrangement or '-'])

        summary_table = Table(summary_data, colWidths=[5*cm, 9*cm])
        summary_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#34495e')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
            ('FONTNAME', (0, 0), (-1, -1), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, -1), 10),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.gray),
            ('TOPPADDING', (0, 0), (-1, -1), 8),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 8),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.HexColor('#ecf0f1'), colors.white]),
        ]))
        elements.append(summary_table)

        return elements

    def _build_diagram_section(self, outputs):
        """Build cross-section diagram."""
        elements = []

        elements.append(PageBreak())
        elements.append(Paragraph("APPENDIX: REINFORCEMENT DRAWING", self.styles['SectionHeader']))

        try:
            # Generate cross-section diagram
            img_bytes = generate_cross_section(
                width=outputs.girder_width,
                depth=outputs.girder_depth,
                cover=outputs.clear_cover,
                tension_bars=outputs.flexure.tension_bar_count,
                tension_bar_dia=outputs.flexure.tension_bar_diameter,
                stirrup_dia=outputs.shear.stirrup_diameter,
                compression_bars=outputs.flexure.compression_bar_count if outputs.flexure.is_doubly_reinforced else None,
                compression_bar_dia=outputs.flexure.compression_bar_diameter if outputs.flexure.is_doubly_reinforced else None,
                return_figure=False  # Get PNG bytes
            )

            # Convert to ReportLab Image
            img = Image(io.BytesIO(img_bytes), width=12*cm, height=16*cm)
            elements.append(img)

        except Exception as e:
            elements.append(Paragraph(
                f"<i>Diagram generation failed: {str(e)}</i>",
                self.styles['BodyText']
            ))

        return elements

    def _create_data_table(self, data):
        """Create a formatted data table."""
        # Determine column widths based on number of columns
        if len(data[0]) == 3:
            col_widths = [5*cm, 4*cm, 2*cm]
        elif len(data[0]) == 4:
            col_widths = [4*cm, 4*cm, 3*cm, 2*cm]
        else:
            col_widths = None

        table = Table(data, colWidths=col_widths)
        table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#3498db')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, -1), 9),
            ('ALIGN', (1, 0), (-1, -1), 'CENTER'),
            ('ALIGN', (0, 0), (0, -1), 'LEFT'),
            ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.gray),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#f8f9fa')]),
            ('TOPPADDING', (0, 0), (-1, -1), 6),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 6),
        ]))
        return table

    def _create_check_table(self, data):
        """Create a check status table with color coding."""
        table = Table(data, colWidths=[8*cm, 3*cm])
        style = [
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#34495e')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, -1), 9),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('ALIGN', (-1, 0), (-1, -1), 'CENTER'),
            ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.gray),
            ('TOPPADDING', (0, 0), (-1, -1), 6),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 6),
        ]

        # Color code PASS/FAIL
        for i, row in enumerate(data[1:], 1):
            if row[1] == 'PASS':
                style.append(('TEXTCOLOR', (1, i), (1, i), colors.HexColor('#27ae60')))
                style.append(('FONTNAME', (1, i), (1, i), 'Helvetica-Bold'))
            elif row[1] == 'FAIL':
                style.append(('TEXTCOLOR', (1, i), (1, i), colors.HexColor('#e74c3c')))
                style.append(('FONTNAME', (1, i), (1, i), 'Helvetica-Bold'))

        table.setStyle(TableStyle(style))
        return table

    def _get_stress_table_style(self, data):
        """Get table style for stress check tables with PASS/FAIL coloring."""
        style = [
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#3498db')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, -1), 9),
            ('ALIGN', (1, 0), (-1, -1), 'CENTER'),
            ('ALIGN', (0, 0), (0, -1), 'LEFT'),
            ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.gray),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#f8f9fa')]),
            ('TOPPADDING', (0, 0), (-1, -1), 5),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 5),
        ]

        # Color code Status column
        for i, row in enumerate(data[1:], 1):
            if row[-1] == 'PASS':
                style.append(('TEXTCOLOR', (-1, i), (-1, i), colors.HexColor('#27ae60')))
                style.append(('FONTNAME', (-1, i), (-1, i), 'Helvetica-Bold'))
            elif row[-1] == 'FAIL':
                style.append(('TEXTCOLOR', (-1, i), (-1, i), colors.HexColor('#e74c3c')))
                style.append(('FONTNAME', (-1, i), (-1, i), 'Helvetica-Bold'))

        return TableStyle(style)

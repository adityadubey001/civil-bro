"""
PDF report generator using ReportLab.
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

from src.models.inputs import BeamDesignInput
from src.models.outputs import BeamDesignOutput, DesignStatus
from src.reports.diagrams.cross_section import generate_cross_section


class PDFReportGenerator:
    """
    Generate PDF design reports using ReportLab.
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

    def generate_report(
        self,
        inputs: BeamDesignInput,
        outputs: BeamDesignOutput,
        project_name: str = "RC Beam Design",
        engineer_name: Optional[str] = None,
    ) -> bytes:
        """
        Generate complete PDF report.

        Args:
            inputs: Design input parameters
            outputs: Design calculation results
            project_name: Project name for report header
            engineer_name: Engineer name (optional)

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
        story.extend(self._build_input_summary(inputs))

        # Results summary
        story.extend(self._build_results_summary(outputs))

        # Flexural design details
        story.extend(self._build_flexure_section(outputs))

        # Shear design details
        story.extend(self._build_shear_section(outputs))

        # Serviceability check
        story.extend(self._build_serviceability_section(outputs))

        # Build PDF
        doc.build(story)

        buffer.seek(0)
        return buffer.getvalue()

    def _build_title_page(self, project_name, outputs, engineer_name):
        """Build title page elements."""
        elements = []

        # Title
        elements.append(Paragraph(
            "RC BEAM DESIGN REPORT",
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
            ['Section', f'{outputs.beam_width:.0f} × {outputs.beam_depth:.0f} mm'],
            ['Main Steel', outputs.flexure.tension_bar_arrangement],
            ['Stirrups', f'{outputs.shear.stirrup_legs}L-{outputs.shear.stirrup_diameter}φ @ {outputs.shear.spacing_at_support:.0f}mm c/c'],
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

    def _build_input_summary(self, inputs):
        """Build input summary section."""
        elements = []

        elements.append(Paragraph("1. INPUT DATA", self.styles['SectionHeader']))

        # Geometry
        elements.append(Paragraph("1.1 Geometry", self.styles['SubSection']))
        geo_data = [
            ['Parameter', 'Value', 'Unit'],
            ['Clear Span', f'{inputs.geometry.span_length:.1f}', 'm'],
            ['Support Type', 'Simply Supported', '-'],
        ]
        if inputs.geometry.width:
            geo_data.append(['Custom Width', f'{inputs.geometry.width:.0f}', 'mm'])
        if inputs.geometry.depth:
            geo_data.append(['Custom Depth', f'{inputs.geometry.depth:.0f}', 'mm'])

        elements.append(self._create_data_table(geo_data))
        elements.append(Spacer(1, 10))

        # Materials
        elements.append(Paragraph("1.2 Materials", self.styles['SubSection']))
        mat_data = [
            ['Parameter', 'Value', 'Unit'],
            ['Concrete Grade', inputs.materials.concrete_grade.value, '-'],
            ['Steel Grade', inputs.materials.steel_grade.value, '-'],
            ['fck', f'{inputs.materials.fck:.0f}', 'MPa'],
            ['fy', f'{inputs.materials.fy:.0f}', 'MPa'],
        ]
        elements.append(self._create_data_table(mat_data))
        elements.append(Spacer(1, 10))

        # Loads
        elements.append(Paragraph("1.3 Loading", self.styles['SubSection']))
        load_data = [
            ['Load Type', 'Value', 'Unit'],
            ['Dead Load (SDL)', f'{inputs.dead_load:.1f}', 'kN/m'],
            ['Live Load', f'{inputs.live_load:.1f}', 'kN/m'],
            ['Total (unfactored)', f'{inputs.dead_load + inputs.live_load:.1f}', 'kN/m'],
        ]
        elements.append(self._create_data_table(load_data))
        elements.append(Spacer(1, 10))

        # Exposure
        elements.append(Paragraph("1.4 Exposure Condition", self.styles['SubSection']))
        elements.append(Paragraph(
            f"Exposure: {inputs.exposure_condition.replace('_', ' ').title()}",
            self.styles['BodyText']
        ))

        return elements

    def _build_results_summary(self, outputs):
        """Build results summary section."""
        elements = []

        elements.append(Paragraph("2. FINAL SECTION", self.styles['SectionHeader']))

        # Section dimensions
        dim_data = [
            ['Parameter', 'Symbol', 'Value', 'Unit'],
            ['Beam Width', 'b', f'{outputs.beam_width:.0f}', 'mm'],
            ['Overall Depth', 'D', f'{outputs.beam_depth:.0f}', 'mm'],
            ['Effective Depth', 'd', f'{outputs.effective_depth:.1f}', 'mm'],
            ['Clear Cover', 'c', f'{outputs.clear_cover:.0f}', 'mm'],
            ['Self Weight', 'sw', f'{outputs.self_weight:.2f}', 'kN/m'],
        ]
        elements.append(self._create_data_table(dim_data))
        elements.append(Spacer(1, 10))

        # Design loads
        elements.append(Paragraph("2.1 Design Loads", self.styles['SubSection']))
        load_data = [
            ['Load', 'Value', 'Unit'],
            ['Factored Dead Load', f'{outputs.factored_dead_load:.2f}', 'kN/m'],
            ['Factored Live Load', f'{outputs.factored_live_load:.2f}', 'kN/m'],
            ['Total Factored Load (wu)', f'{outputs.total_factored_load:.2f}', 'kN/m'],
        ]
        elements.append(self._create_data_table(load_data))

        return elements

    def _build_flexure_section(self, outputs):
        """Build flexural design section."""
        elements = []

        elements.append(Paragraph("3. FLEXURAL DESIGN", self.styles['SectionHeader']))

        flexure = outputs.flexure

        # Design forces
        elements.append(Paragraph("3.1 Design Moment", self.styles['SubSection']))
        moment_data = [
            ['Parameter', 'Value', 'Unit'],
            ['Design Moment (Mu)', f'{flexure.design_moment:.2f}', 'kNm'],
            ['Moment Capacity (Mu,lim)', f'{flexure.moment_capacity_limit:.2f}', 'kNm'],
            ['Moment Capacity (provided)', f'{flexure.moment_capacity:.2f}', 'kNm'],
            ['Utilization', f'{flexure.utilization_ratio:.2f}', '-'],
        ]
        elements.append(self._create_data_table(moment_data))
        elements.append(Spacer(1, 10))

        # Reinforcement
        elements.append(Paragraph("3.2 Tension Reinforcement", self.styles['SubSection']))
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
            elements.append(Paragraph("3.3 Compression Reinforcement", self.styles['SubSection']))
            comp_data = [
                ['Parameter', 'Value', 'Unit'],
                ['Required Asc', f'{flexure.required_asc:.0f}', 'mm²'],
                ['Provided Asc', f'{flexure.provided_asc:.0f}', 'mm²'],
                ['Bar Arrangement', flexure.compression_bar_arrangement, '-'],
            ]
            elements.append(self._create_data_table(comp_data))
            elements.append(Spacer(1, 10))

        # Code checks
        elements.append(Paragraph("3.4 Code Checks", self.styles['SubSection']))
        check_data = [
            ['Check', 'Status'],
            ['Under-reinforced section', 'PASS' if flexure.is_under_reinforced else 'FAIL'],
            ['Minimum reinforcement', 'PASS' if flexure.min_ast_satisfied else 'FAIL'],
            ['Maximum reinforcement', 'PASS' if flexure.max_ast_satisfied else 'FAIL'],
        ]
        elements.append(self._create_check_table(check_data))

        return elements

    def _build_shear_section(self, outputs):
        """Build shear design section."""
        elements = []

        elements.append(Paragraph("4. SHEAR DESIGN", self.styles['SectionHeader']))

        shear = outputs.shear

        # Design forces
        elements.append(Paragraph("4.1 Design Shear", self.styles['SubSection']))
        shear_data = [
            ['Parameter', 'Value', 'Unit'],
            ['Design Shear (Vu)', f'{shear.design_shear:.2f}', 'kN'],
            ['Shear at d (Vu,d)', f'{shear.shear_at_critical_section:.2f}', 'kN'],
            ['Nominal stress (τv)', f'{shear.nominal_shear_stress:.3f}', 'MPa'],
            ['Concrete strength (τc)', f'{shear.concrete_shear_strength:.3f}', 'MPa'],
            ['Max stress (τc,max)', f'{shear.max_shear_stress:.2f}', 'MPa'],
        ]
        elements.append(self._create_data_table(shear_data))
        elements.append(Spacer(1, 10))

        # Stirrups
        elements.append(Paragraph("4.2 Shear Reinforcement", self.styles['SubSection']))
        stirrup_data = [
            ['Parameter', 'Value', 'Unit'],
            ['Stirrup diameter', f'{shear.stirrup_diameter}', 'mm'],
            ['Number of legs', f'{shear.stirrup_legs}', '-'],
            ['Spacing at support', f'{shear.spacing_at_support:.0f}', 'mm'],
            ['Spacing at center', f'{shear.spacing_at_center:.0f}', 'mm'],
        ]
        elements.append(self._create_data_table(stirrup_data))
        elements.append(Spacer(1, 10))

        # Code checks
        elements.append(Paragraph("4.3 Code Checks", self.styles['SubSection']))
        check_data = [
            ['Check', 'Status'],
            ['Maximum shear stress', 'PASS' if shear.max_shear_check else 'FAIL'],
            ['Minimum reinforcement', 'PASS' if shear.min_shear_reinf_check else 'FAIL'],
        ]
        elements.append(self._create_check_table(check_data))

        return elements

    def _build_serviceability_section(self, outputs):
        """Build serviceability section."""
        elements = []

        elements.append(Paragraph("5. DEFLECTION CHECK", self.styles['SectionHeader']))

        serv = outputs.serviceability

        # L/d ratio check
        elements.append(Paragraph("5.1 Span/Depth Ratio Method", self.styles['SubSection']))
        ld_data = [
            ['Parameter', 'Value', 'Unit'],
            ['Basic L/d ratio', f'{serv.basic_ld_ratio:.0f}', '-'],
            ['Modification factor (kt)', f'{serv.modification_factor_tension:.2f}', '-'],
            ['Modification factor (kc)', f'{serv.modification_factor_compression:.2f}', '-'],
            ['Allowable L/d', f'{serv.allowable_ld_ratio:.2f}', '-'],
            ['Actual L/d', f'{serv.actual_ld_ratio:.2f}', '-'],
        ]
        elements.append(self._create_data_table(ld_data))
        elements.append(Spacer(1, 10))

        # Check result
        status = 'PASS' if serv.deflection_check_passed else 'FAIL'
        status_color = '#27ae60' if status == 'PASS' else '#e74c3c'
        elements.append(Paragraph(
            f'<font color="{status_color}"><b>DEFLECTION CHECK: {status}</b></font>',
            self.styles['Heading3']
        ))

        return elements

    def _create_data_table(self, data):
        """Create a formatted data table."""
        table = Table(data, colWidths=[5*cm, 4*cm, 2*cm])
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

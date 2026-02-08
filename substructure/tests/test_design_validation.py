"""Validation tests for substructure design against reference values.

These tests verify the design engine produces correct results for the
sample_input.yaml configuration. Reference values are from Excel validation.
"""
import sys
from pathlib import Path

import pytest
import yaml

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from ui.runner import run_design


@pytest.fixture(scope="module")
def design_results():
    """Run the design once and share results across all tests."""
    config_path = Path(__file__).parent.parent / "config" / "sample_input.yaml"
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)
    return run_design(config)


class TestPierCapDesign:
    """Pier cap design validation tests."""

    def test_flexure_utilization(self, design_results):
        """Flexure utilization at pier face should be ~0.80."""
        pc = design_results["pier_cap"]
        assert pc.flexure_util_pier == pytest.approx(0.803, abs=0.01)

    def test_shear_utilization(self, design_results):
        """Shear utilization with fywd=400 MPa should be ~0.93."""
        pc = design_results["pier_cap"]
        # Updated after fywd fix: was 0.776 with fyd=478, now ~0.928 with fywd=400
        assert pc.shear_util == pytest.approx(0.928, abs=0.02)

    def test_crack_width(self, design_results):
        """Crack width should be ~0.265mm (< 0.3mm limit)."""
        pc = design_results["pier_cap"]
        assert pc.crack_width == pytest.approx(0.265, abs=0.01)
        assert pc.crack_width < 0.3  # IRC 112 limit

    def test_status_ok(self, design_results):
        """Pier cap design should pass all checks."""
        pc = design_results["pier_cap"]
        assert pc.status == "OK"


class TestPierDesign:
    """Pier design validation tests."""

    def test_biaxial_utilization(self, design_results):
        """Biaxial P-M utilization should be ~0.88."""
        pier = design_results["pier"]
        assert pier.util_biaxial == pytest.approx(0.884, abs=0.02)

    def test_status_ok(self, design_results):
        """Pier design should pass all checks."""
        pier = design_results["pier"]
        assert pier.status == "OK"


class TestPileCapacity:
    """Pile capacity validation tests."""

    def test_compression_utilization(self, design_results):
        """Pile compression utilization should be ~0.99."""
        pile_cap = design_results["pile_capacity"]
        assert pile_cap.compression_util == pytest.approx(0.991, abs=0.02)

    def test_status_ok(self, design_results):
        """Pile capacity should pass all checks."""
        pile_cap = design_results["pile_capacity"]
        assert pile_cap.status == "OK"


class TestPileDesign:
    """Pile structural design validation tests."""

    def test_pm_utilization(self, design_results):
        """P-M utilization should be ~0.60."""
        pile = design_results["pile_design"]
        assert pile.util_PM == pytest.approx(0.602, abs=0.02)

    def test_status_ok(self, design_results):
        """Pile design should pass all checks."""
        pile = design_results["pile_design"]
        assert pile.status == "OK"


class TestPilecapDesign:
    """Pile cap design validation tests."""

    def test_pier_punching_utilization(self, design_results):
        """Pier punching shear utilization should be ~0.76."""
        pilecap = design_results["pilecap"]
        assert pilecap.punch_util_pier == pytest.approx(0.758, abs=0.02)

    def test_pile_punching_utilization(self, design_results):
        """Pile punching shear utilization should be ~0.995."""
        pilecap = design_results["pilecap"]
        assert pilecap.punch_util_pile == pytest.approx(0.995, abs=0.02)

    def test_status_ok(self, design_results):
        """Pile cap design should pass all checks."""
        pilecap = design_results["pilecap"]
        assert pilecap.status == "OK"


class TestBOQ:
    """Bill of quantities validation tests."""

    def test_total_cost_reasonable(self, design_results):
        """Total cost should be ~15.8 lakhs INR."""
        boq = design_results["boq"]
        # Cost in lakhs (1 lakh = 100,000 INR)
        cost_lakhs = boq.total_cost / 100000
        assert cost_lakhs == pytest.approx(15.8, abs=2.0)


class TestAllChecksPass:
    """Verify all design checks pass."""

    def test_all_status_ok(self, design_results):
        """All five design checks should return OK status."""
        components = ["pier_cap", "pier", "pile_capacity", "pile_design", "pilecap"]
        for comp in components:
            result = design_results[comp]
            assert result.status == "OK", f"{comp} failed with status: {result.status}"

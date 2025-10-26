"""Tests for geometry module."""

import pytest
import numpy as np
from lh2sim.geometry import (
    cyl_v_to_h,
    cylinder_cross_section_area,
    cylinder_lateral_surface_area,
    horizontal_cylinder_liquid_surface_area,
)


class TestCylVToH:
    """Tests for cylinder volume to height conversion."""
    
    def test_empty_tank(self):
        """Test with zero volume."""
        R, L = 1.0, 5.0
        V = 0.0
        H = cyl_v_to_h(V, R, L)
        assert H == 0.0
    
    def test_full_tank(self):
        """Test with full volume."""
        R, L = 1.0, 5.0
        V = np.pi * R**2 * L
        H = cyl_v_to_h(V, R, L)
        assert abs(H - 2 * R) < 1e-3
    
    def test_half_full_tank(self):
        """Test with half-full tank."""
        R, L = 1.0, 5.0
        V = 0.5 * np.pi * R**2 * L
        H = cyl_v_to_h(V, R, L)
        # For horizontal cylinder, half volume should give height = R
        assert abs(H - R) < 0.1  # Allow some tolerance due to iteration
    
    def test_quarter_full_tank(self):
        """Test with quarter-full tank."""
        R, L = 1.0, 5.0
        V = 0.25 * np.pi * R**2 * L
        H = cyl_v_to_h(V, R, L)
        assert 0 < H < R
    
    def test_realistic_dimensions(self):
        """Test with realistic LH2 tank dimensions."""
        R, L = 1.5, 10.0  # Trailer tank size
        V = 0.3 * np.pi * R**2 * L  # 30% full
        H = cyl_v_to_h(V, R, L)
        assert 0 < H < 2 * R
        assert not np.isnan(H)


class TestCylinderGeometry:
    """Tests for basic cylinder geometry functions."""
    
    def test_cross_section_area(self):
        """Test cross-sectional area calculation."""
        R = 1.0
        A = cylinder_cross_section_area(R)
        assert abs(A - np.pi) < 1e-10
    
    def test_lateral_surface_area(self):
        """Test lateral surface area calculation."""
        R, L = 1.0, 5.0
        A = cylinder_lateral_surface_area(R, L)
        assert abs(A - 2 * np.pi * R * L) < 1e-10
    
    def test_liquid_surface_area_empty(self):
        """Test liquid surface area for empty tank."""
        R, L = 1.0, 5.0
        V = 0.0
        A = horizontal_cylinder_liquid_surface_area(V, R, L)
        assert A == 0.0
    
    def test_liquid_surface_area_full(self):
        """Test liquid surface area for full tank."""
        R, L = 1.0, 5.0
        V = np.pi * R**2 * L
        A = horizontal_cylinder_liquid_surface_area(V, R, L)
        assert A == 0.0  # No interface when completely full
    
    def test_liquid_surface_area_half_full(self):
        """Test liquid surface area for half-full tank."""
        R, L = 1.0, 5.0
        V = 0.5 * np.pi * R**2 * L
        A = horizontal_cylinder_liquid_surface_area(V, R, L)
        # At half full, chord length = 2R, so area = 2R * L
        expected = 2 * R * L
        assert abs(A - expected) < 0.5  # Allow tolerance due to iteration

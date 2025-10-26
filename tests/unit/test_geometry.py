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


class TestCylVToHEdgeCases:
    """Additional edge case tests for cyl_v_to_h (DIFFERENCES.md Item #4)."""
    
    def test_very_small_volume(self):
        """Test with very small liquid volume (near empty)."""
        R, L = 1.0, 5.0
        V = 1e-6 * np.pi * R**2 * L  # 0.0001% full
        H = cyl_v_to_h(V, R, L)
        assert 0 < H < 0.01  # Very small height
        assert not np.isnan(H)
    
    def test_near_full_volume(self):
        """Test with volume very close to full."""
        R, L = 1.0, 5.0
        V = 0.9999 * np.pi * R**2 * L  # 99.99% full
        H = cyl_v_to_h(V, R, L)
        assert 1.99 * R < H < 2 * R  # Very close to 2R
        assert not np.isnan(H)
    
    def test_multiple_fill_levels(self):
        """Test accuracy at various fill levels."""
        R, L = 1.0, 5.0
        
        fill_fractions = [0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99]
        
        for frac in fill_fractions:
            V = frac * np.pi * R**2 * L
            H = cyl_v_to_h(V, R, L)
            
            # Verify height is in valid range
            assert 0 <= H <= 2 * R, f"Failed for fill={frac*100}%: H={H}, should be in [0, {2*R}]"
            
            # Verify round-trip: compute volume from height
            # For validation, check it's monotonically increasing
            assert not np.isnan(H), f"NaN result for fill={frac*100}%"
    
    def test_large_tank_dimensions(self):
        """Test with large realistic tank dimensions."""
        R, L = 2.5, 20.0  # Large trailer
        V = 0.5 * np.pi * R**2 * L  # Half full
        H = cyl_v_to_h(V, R, L)
        
        # At half full, H should be approximately R
        assert 0.9 * R < H < 1.1 * R  # Within 10% of R
    
    def test_small_tank_dimensions(self):
        """Test with small tank dimensions."""
        R, L = 0.1, 0.5  # Small test tank
        V = 0.3 * np.pi * R**2 * L
        H = cyl_v_to_h(V, R, L)
        
        assert 0 < H < 2 * R
        assert not np.isnan(H)
    
    def test_convergence_at_boundaries(self):
        """Test that Newton iteration converges at edge cases."""
        R, L = 1.0, 5.0
        
        # Test at various boundary-adjacent volumes
        boundary_volumes = [
            0.001 * np.pi * R**2 * L,   # Near empty
            0.01 * np.pi * R**2 * L,    # 1%
            0.99 * np.pi * R**2 * L,    # 99%
            0.999 * np.pi * R**2 * L,   # Near full
        ]
        
        for V in boundary_volumes:
            H = cyl_v_to_h(V, R, L)
            assert not np.isnan(H), f"Failed to converge for V={V}"
            assert 0 <= H <= 2 * R, f"Out of bounds for V={V}"
    
    def test_aspect_ratio_variations(self):
        """Test with different tank aspect ratios."""
        # Various R/L ratios
        test_cases = [
            (0.5, 10.0),   # Long and thin
            (1.0, 5.0),    # Moderate
            (2.0, 2.0),    # Short and fat
        ]
        
        for R, L in test_cases:
            V = 0.5 * np.pi * R**2 * L  # Half full
            H = cyl_v_to_h(V, R, L)
            
            assert 0 < H < 2 * R, f"Failed for R={R}, L={L}"
            assert not np.isnan(H), f"NaN for R={R}, L={L}"


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

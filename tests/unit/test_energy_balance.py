"""
Tests for multi-node energy balance module.

Tests the helper functions for multi-node energy balance calculations
following the MATLAB LLNL implementation.
"""

import pytest
import numpy as np
from lh2sim.simulation.energy_balance import (
    generate_boundary_layer_grid,
    compute_surface_temperature,
    compute_latent_heat,
    compute_natural_convection_heat_transfer_coefficient,
    compute_wall_convection_nusselt,
    compute_condensation_rate,
    compute_temperature_from_internal_energy,
)


class TestBoundaryLayerGrid:
    """Test exponential grid generation for boundary layers."""

    def test_single_node_grid(self):
        """Test grid with single node."""
        n = 1
        kappa = 0.1  # W/m/K
        tmin = 10.0  # s
        c = 9000.0  # J/kg/K
        rho = 70.0  # kg/m³
        
        l, l12 = generate_boundary_layer_grid(n, kappa, tmin, c, rho)
        
        assert len(l) == 1
        assert len(l12) == 1
        assert l[0] > 0
        assert l12[0] > 0

    def test_multi_node_grid(self):
        """Test grid with multiple nodes."""
        n = 5
        kappa = 0.1
        tmin = 10.0
        c = 9000.0
        rho = 70.0
        
        l, l12 = generate_boundary_layer_grid(n, kappa, tmin, c, rho)
        
        assert len(l) == n
        assert len(l12) == n
        
        # Grid should be increasing (moving away from interface)
        assert np.all(np.diff(l) > 0)
        assert np.all(np.diff(l12) > 0)
        
        # Check exponential spacing
        # l12[i+1] / l12[i] should be approximately constant
        ratios = l12[1:] / l12[:-1]
        assert np.allclose(ratios, ratios[0], rtol=1e-10)

    def test_first_node_at_interface(self):
        """Test that first node is closest to interface."""
        n = 3
        kappa = 0.1
        tmin = 10.0
        c = 9000.0
        rho = 70.0
        
        l, l12 = generate_boundary_layer_grid(n, kappa, tmin, c, rho)
        
        # First node should be smallest distance
        assert l[0] < l[1]
        assert l[0] < l[2]


class TestSurfaceTemperature:
    """Test surface temperature correlation."""

    def test_surface_temperature_at_critical_point(self):
        """Test that surface temperature equals critical temperature at critical pressure."""
        T_c = 32.938  # K
        p_c = 1.2884e6  # Pa
        lambda_ = 1.5
        
        Ts = compute_surface_temperature(p_c, T_c, p_c, lambda_)
        
        assert np.isclose(Ts, T_c)

    def test_surface_temperature_below_critical(self):
        """Test surface temperature below critical point."""
        T_c = 32.938  # K
        p_c = 1.2884e6  # Pa
        lambda_ = 1.5
        p_v = 0.5 * p_c  # Half critical pressure
        
        Ts = compute_surface_temperature(p_v, T_c, p_c, lambda_)
        
        assert Ts < T_c
        assert Ts > 13.804  # Above triple point

    def test_surface_temperature_monotonic(self):
        """Test that surface temperature increases with pressure."""
        T_c = 32.938
        p_c = 1.2884e6
        lambda_ = 1.5
        
        # Use pressure range above minimum to avoid clamping effects
        pressures = np.linspace(2e5, p_c, 10)  # Start at 2 bar instead of 1 bar
        temperatures = [compute_surface_temperature(p, T_c, p_c, lambda_) for p in pressures]
        
        # Should be monotonically increasing (allowing for numerical precision)
        diffs = np.diff(temperatures)
        assert np.all(diffs >= -1e-10), "Temperature should be non-decreasing with pressure"


class TestLatentHeat:
    """Test latent heat of vaporization."""

    def test_latent_heat_positive(self):
        """Test that latent heat is positive for typical LH2 temperatures."""
        # Note: MATLAB polynomial is valid for T > ~17K based on empirical testing
        T_values = np.linspace(18.0, 30.0, 10)
        
        for T in T_values:
            qh = compute_latent_heat(T)
            assert qh > 0, f"Latent heat should be positive at T={T}K, got {qh}"

    def test_latent_heat_decreases_with_temperature(self):
        """Test that latent heat decreases with increasing temperature."""
        T_low = 20.0  # K
        T_high = 30.0
        
        qh_low = compute_latent_heat(T_low)
        qh_high = compute_latent_heat(T_high)
        
        assert qh_low > qh_high, "Latent heat should decrease with temperature"

    def test_latent_heat_realistic_values(self):
        """Test that latent heat is in realistic range for LH2."""
        T = 20.0  # K
        qh = compute_latent_heat(T)
        
        # LH2 latent heat is typically 400-450 kJ/kg
        assert 300e3 < qh < 500e3, f"Latent heat {qh/1e3:.1f} kJ/kg is outside expected range"


class TestNaturalConvection:
    """Test natural convection correlations."""

    def test_zero_delta_T_gives_zero_h(self):
        """Test that zero temperature difference gives zero heat transfer coefficient."""
        kappa = 0.1
        g = 9.81
        beta = 0.01
        cp = 9000.0
        rho = 70.0
        mu = 1e-5
        delta_T = 0.0
        length_scale = 1.0
        
        h = compute_natural_convection_heat_transfer_coefficient(
            kappa, g, beta, cp, rho, mu, delta_T, length_scale
        )
        
        assert h == 0.0

    def test_positive_h_for_positive_delta_T(self):
        """Test that positive temperature difference gives positive heat transfer coefficient."""
        kappa = 0.1
        g = 9.81
        beta = 0.01
        cp = 9000.0
        rho = 70.0
        mu = 1e-5
        delta_T = 5.0
        length_scale = 1.0
        
        h = compute_natural_convection_heat_transfer_coefficient(
            kappa, g, beta, cp, rho, mu, delta_T, length_scale
        )
        
        assert h > 0


class TestWallConvectionNusselt:
    """Test wall convection Nusselt number correlations."""

    def test_vertical_wall_nusselt(self):
        """Test Nusselt number for vertical wall."""
        g = 9.81
        beta = 0.01
        delta_T = 5.0
        length_scale = 1.0
        nu = 1e-6
        Pr = 0.7
        
        Nu = compute_wall_convection_nusselt(g, beta, delta_T, length_scale, nu, Pr, "vertical_wall")
        
        assert Nu >= 0.68, "Nusselt number should be at least the minimum value"

    def test_horizontal_plate_nusselt(self):
        """Test Nusselt number for horizontal plate."""
        g = 9.81
        beta = 0.01
        delta_T = 5.0
        length_scale = 1.0
        nu = 1e-6
        Pr = 0.7
        
        Nu = compute_wall_convection_nusselt(
            g, beta, delta_T, length_scale, nu, Pr, "horizontal_plate"
        )
        
        assert Nu > 0

    def test_zero_delta_T_gives_minimum_nusselt(self):
        """Test that zero temperature difference gives minimum Nusselt number."""
        g = 9.81
        beta = 0.01
        delta_T = 0.0
        length_scale = 1.0
        nu = 1e-6
        Pr = 0.7
        
        Nu = compute_wall_convection_nusselt(g, beta, delta_T, length_scale, nu, Pr, "vertical_wall")
        
        assert Nu == 1.0

    def test_invalid_geometry_raises_error(self):
        """Test that invalid geometry raises ValueError."""
        with pytest.raises(ValueError, match="Unknown geometry"):
            compute_wall_convection_nusselt(9.81, 0.01, 5.0, 1.0, 1e-6, 0.7, "invalid")


class TestCondensationRate:
    """Test condensation rate calculation."""

    def test_zero_heat_flow_gives_zero_condensation(self):
        """Test that zero heat flows give zero condensation."""
        Q_L = 0.0
        Q_V = 0.0
        qh = 400e3
        
        J_cd = compute_condensation_rate(Q_L, Q_V, qh)
        
        assert J_cd == 0.0

    def test_positive_heat_flow_gives_condensation(self):
        """Test that heat flow to interface gives condensation (positive rate)."""
        Q_L = 1000.0  # W from liquid
        Q_V = 1000.0  # W from vapor
        qh = 400e3  # J/kg
        
        J_cd = compute_condensation_rate(Q_L, Q_V, qh)
        
        # Both heat flows to interface → condensation
        assert J_cd < 0  # Negative because -(Q_L + Q_V) / qh

    def test_zero_latent_heat_gives_zero_condensation(self):
        """Test that zero latent heat gives zero condensation."""
        Q_L = 1000.0
        Q_V = 1000.0
        qh = 0.0
        
        J_cd = compute_condensation_rate(Q_L, Q_V, qh)
        
        assert J_cd == 0.0


class TestTemperatureFromInternalEnergy:
    """Test temperature computation from internal energy."""

    def test_liquid_temperature_in_valid_range(self):
        """Test that liquid temperature is in valid range for typical LH2."""
        # Typical LH2 internal energy
        u = 200e3  # J/kg
        
        T = compute_temperature_from_internal_energy(u, "liquid")
        
        assert 13.804 <= T <= 32.93, f"Temperature {T}K outside valid range"

    def test_vapor_temperature_from_internal_energy(self):
        """Test vapor temperature calculation."""
        # For ideal gas: u = c_v * T
        c_v = 6490.0  # J/kg/K
        T_expected = 25.0  # K
        u = c_v * T_expected
        
        T = compute_temperature_from_internal_energy(u, "vapor")
        
        assert np.isclose(T, T_expected)

    def test_invalid_fluid_type_raises_error(self):
        """Test that invalid fluid type raises ValueError."""
        with pytest.raises(ValueError, match="Unknown fluid type"):
            compute_temperature_from_internal_energy(100e3, "invalid")

    def test_liquid_temperature_clamped_to_valid_range(self):
        """Test that liquid temperature is clamped to physical limits."""
        # Very low internal energy
        u_low = 0.0
        T_low = compute_temperature_from_internal_energy(u_low, "liquid")
        assert T_low >= 13.804
        
        # Very high internal energy
        u_high = 1e6
        T_high = compute_temperature_from_internal_energy(u_high, "liquid")
        assert T_high <= 32.93


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

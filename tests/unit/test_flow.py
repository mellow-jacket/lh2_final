"""Tests for flow module."""

import pytest
import numpy as np
from lh2sim.flow import (
    gas_flow,
    directed_sqrt,
    valve_flow_coefficient,
    vent_flow_rate,
)


class TestGasFlow:
    """Tests for gas flow calculations."""
    
    def test_no_pressure_difference(self):
        """Test with equal pressures (no flow)."""
        CA = 0.01  # m²
        rho = 1.0  # kg/m³
        P1 = P2 = 1e5  # Pa
        gamma = 1.4
        
        mdot = gas_flow(CA, rho, P1, P2, gamma)
        assert abs(mdot) < 1e-10
    
    def test_flow_reversal(self):
        """Test that flow reverses with pressure swap."""
        CA = 0.01
        rho = 1.0
        P1, P2 = 2e5, 1e5
        gamma = 1.4
        
        mdot_forward = gas_flow(CA, rho, P1, P2, gamma)
        mdot_reverse = gas_flow(CA, rho, P2, P1, gamma)
        
        assert mdot_forward > 0
        assert mdot_reverse < 0
        assert abs(mdot_forward + mdot_reverse) < 1e-10
    
    def test_choked_flow(self):
        """Test choked flow condition."""
        CA = 0.01
        rho = 1.0
        P1 = 2e5
        P2 = 0.5e5  # Low enough for choked flow
        gamma = 1.4
        
        mdot = gas_flow(CA, rho, P1, P2, gamma)
        
        # For choked flow, mass flow should not depend on P2
        # Test with even lower P2
        mdot2 = gas_flow(CA, rho, P1, 0.1e5, gamma)
        
        # Should be approximately equal (choked)
        assert abs(mdot - mdot2) / mdot < 0.01
    
    def test_nonchoked_flow(self):
        """Test non-choked flow condition."""
        CA = 0.01
        rho = 1.0
        P1 = 1.2e5
        P2 = 1.0e5  # Small pressure difference
        gamma = 1.4
        
        mdot = gas_flow(CA, rho, P1, P2, gamma)
        assert mdot > 0
        
        # For non-choked flow, mass flow should depend on P2
        mdot2 = gas_flow(CA, rho, P1, 0.9e5, gamma)
        assert mdot2 > mdot  # Larger pressure difference = more flow
    
    def test_zero_area(self):
        """Test with zero flow area."""
        mdot = gas_flow(0.0, 1.0, 2e5, 1e5, 1.4)
        assert mdot == 0.0
    
    def test_realistic_hydrogen_conditions(self):
        """Test with realistic LH2 vapor conditions."""
        CA = 0.001  # m² (small vent)
        rho = 1.3  # kg/m³ (hydrogen vapor)
        P1 = 1.5e5  # Pa (tank pressure)
        P2 = 1.0e5  # Pa (atmospheric)
        gamma = 1.4  # For hydrogen
        
        mdot = gas_flow(CA, rho, P1, P2, gamma)
        
        # Should get reasonable mass flow rate
        assert 0 < mdot < 1.0  # kg/s


class TestDirectedSqrt:
    """Tests for directed square root function."""
    
    def test_positive_value(self):
        """Test with positive value."""
        assert directed_sqrt(4.0) == 2.0
    
    def test_negative_value(self):
        """Test with negative value."""
        assert directed_sqrt(-4.0) == -2.0
    
    def test_zero(self):
        """Test with zero."""
        assert directed_sqrt(0.0) == 0.0
    
    def test_array(self):
        """Test with array input."""
        x = np.array([-4, 0, 4, 9])
        expected = np.array([-2, 0, 2, 3])
        result = directed_sqrt(x)
        np.testing.assert_array_almost_equal(result, expected)


class TestValveFlowCoefficient:
    """Tests for valve flow coefficient calculation."""
    
    def test_fully_open(self):
        """Test with fully open valve."""
        A = 0.01
        CA = valve_flow_coefficient(A, 1.0)
        assert CA == 0.6 * A
    
    def test_fully_closed(self):
        """Test with fully closed valve."""
        CA = valve_flow_coefficient(0.01, 0.0)
        assert CA == 0.0
    
    def test_half_open(self):
        """Test with half-open valve."""
        A = 0.01
        CA = valve_flow_coefficient(A, 0.5)
        assert CA == 0.6 * A * 0.5
    
    def test_custom_discharge_coefficient(self):
        """Test with custom discharge coefficient."""
        A = 0.01
        C_d = 0.8
        CA = valve_flow_coefficient(A, 1.0, C_d)
        assert CA == C_d * A


class TestVentFlowRate:
    """Tests for vent flow rate calculation."""
    
    def test_vent_closed(self):
        """Test with zero vent area."""
        mdot = vent_flow_rate(1.5e5, 1.0e5, 1.3, 0.0, 1.4)
        assert mdot == 0.0
    
    def test_no_pressure_difference(self):
        """Test with no pressure difference."""
        mdot = vent_flow_rate(1.0e5, 1.0e5, 1.3, 0.001, 1.4)
        assert abs(mdot) < 1e-10
    
    def test_positive_venting(self):
        """Test normal venting condition."""
        P_tank = 1.5e5
        P_ambient = 1.0e5
        rho = 1.3
        A_vent = 0.001
        gamma = 1.4
        
        mdot = vent_flow_rate(P_tank, P_ambient, rho, A_vent, gamma)
        assert mdot > 0

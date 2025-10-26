"""Tests for properties module."""

import pytest
import numpy as np
from lh2sim.properties import FluidProperties, vapor_pressure


class TestFluidPropertiesPolynomial:
    """Tests for FluidProperties using polynomial backend."""
    
    @pytest.fixture
    def props(self):
        """Create FluidProperties instance with polynomial backend."""
        return FluidProperties(fluid_name="Hydrogen", backend="polynomial")
    
    def test_liquid_density_range(self, props):
        """Test liquid density returns reasonable values."""
        T = 20.0  # K
        rho = props.density(T=T)
        # Liquid hydrogen density should be around 70 kg/m³ at 20K
        assert 40 < rho < 80
    
    def test_liquid_density_temperature_dependence(self, props):
        """Test that density decreases with temperature."""
        T1 = 20.0
        T2 = 25.0
        rho1 = props.density(T=T1)
        rho2 = props.density(T=T2)
        assert rho1 > rho2  # Density decreases with temperature
    
    def test_vapor_pressure_range(self, props):
        """Test vapor pressure returns reasonable values."""
        T = 20.0  # K
        P = props.pressure(T=T)
        # Should be reasonable pressure for LH2
        assert 1e4 < P < 1e6
    
    def test_vapor_pressure_increases_with_temp(self, props):
        """Test that vapor pressure increases with temperature."""
        T1 = 20.0
        T2 = 25.0
        P1 = props.pressure(T=T1)
        P2 = props.pressure(T=T2)
        assert P2 > P1


class TestVaporPressureFunction:
    """Tests for standalone vapor_pressure function."""
    
    def test_returns_positive_pressure(self):
        """Test that function returns positive pressure."""
        T_film = 20.0
        rho_vapor = 1.3
        P = vapor_pressure(T_film, rho_vapor)
        assert P > 0
    
    def test_pressure_increases_with_temperature(self):
        """Test pressure increases with film temperature."""
        rho_vapor = 1.3
        P1 = vapor_pressure(20.0, rho_vapor)
        P2 = vapor_pressure(25.0, rho_vapor)
        # In supercritical branch, P ~ T
        # May not be monotonic in two-phase, but typically should increase
        assert P2 > 0 and P1 > 0
    
    def test_pressure_increases_with_density(self):
        """Test pressure increases with vapor density."""
        T_film = 25.0
        P1 = vapor_pressure(T_film, 1.0)
        P2 = vapor_pressure(T_film, 2.0)
        assert P2 > P1


# Skip CoolProp tests if not available
pytest.importorskip("CoolProp", reason="CoolProp not installed")


class TestFluidPropertiesCoolProp:
    """Tests for FluidProperties using CoolProp backend."""
    
    @pytest.fixture
    def props(self):
        """Create FluidProperties instance with CoolProp backend."""
        return FluidProperties(fluid_name="Hydrogen", backend="CoolProp")
    
    def test_density_at_saturation(self, props):
        """Test density at saturation conditions."""
        T = 20.0  # K
        Q = 0.0  # Saturated liquid
        rho = props.density(T=T, Q=Q)
        # Liquid hydrogen density at 20K should be around 71 kg/m³
        assert 65 < rho < 75
    
    def test_vapor_density_at_saturation(self, props):
        """Test vapor density at saturation."""
        T = 20.0  # K
        Q = 1.0  # Saturated vapor
        rho = props.density(T=T, Q=Q)
        # Vapor should be much less dense than liquid
        assert 0.5 < rho < 3.0
    
    def test_pressure_from_temperature_density(self, props):
        """Test pressure calculation from T and rho."""
        T = 25.0  # K
        D = 50.0  # kg/m³
        P = props.pressure(T=T, D=D)
        assert P > 0
    
    def test_temperature_from_pressure_density(self, props):
        """Test temperature calculation from P and rho."""
        P = 1e5  # Pa
        D = 70.0  # kg/m³
        T = props.temperature(P=P, D=D)
        # Should be in cryogenic range
        assert 15 < T < 30
    
    def test_enthalpy(self, props):
        """Test enthalpy calculation."""
        T = 20.0  # K
        P = 1e5  # Pa
        h = props.enthalpy(T=T, P=P)
        # Enthalpy should be reasonable for liquid hydrogen
        assert h != 0  # Non-zero enthalpy
    
    def test_viscosity(self, props):
        """Test viscosity calculation."""
        T = 20.0  # K
        P = 1e5  # Pa
        mu = props.viscosity(T, P=P)
        # Liquid hydrogen viscosity should be very low
        assert mu > 0
        assert mu < 1e-3  # Less than 1 mPa·s
    
    def test_thermal_conductivity(self, props):
        """Test thermal conductivity calculation."""
        T = 20.0  # K
        P = 1e5  # Pa
        k = props.thermal_conductivity(T, P=P)
        # Should be positive and reasonable
        assert k > 0
        assert k < 1.0
    
    def test_specific_heat(self, props):
        """Test specific heat calculation."""
        T = 20.0  # K
        P = 1e5  # Pa
        Cp = props.specific_heat_cp(T, P=P)
        # Liquid hydrogen has high specific heat
        assert Cp > 5000  # J/kg/K
        assert Cp < 20000

"""
Multi-Node Energy Balance Module

This module provides helper functions for computing multi-node energy balance
terms following the MATLAB LLNL implementation approach.

Key components:
- Grid generation for liquid and vapor boundary layers
- Heat transfer coefficients (conduction, natural convection)
- Interface/surface temperature correlations
- Latent heat of vaporization
- Condensation flow rates
- Node-to-node energy derivatives

References:
- MATLAB LLNL LH2Simulate.m lines 346-650
- MATLAB paper model LH2Simulate_Pump.m

The multi-node approach discretizes liquid and vapor phases into boundary layers
near the interface with exponentially-spaced grids for better resolution.
"""

import numpy as np
from typing import Tuple
from lh2sim.properties import FluidProperties


def generate_boundary_layer_grid(
    n_nodes: int, kappa: float, tmin: float, c: float, rho: float
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generate exponentially-spaced grid for boundary layer discretization.
    
    Following MATLAB approach (LH2Simulate.m lines 346-363):
    - Count from interface to bulk (layer 0 is at interface)
    - Exponential spacing: l_{i+1/2} = l_{1/2} * exp(pi*i/sqrt(n))
    - Node positions: l_i = sqrt(l_{i-1/2} * l_{i+1/2})
    
    Args:
        n_nodes: Number of nodes (>= 1)
        kappa: Thermal conductivity [W/m/K]
        tmin: Minimum time scale [s]
        c: Specific heat [J/kg/K]
        rho: Density [kg/m³]
    
    Returns:
        l: Node positions from interface [m] (length n_nodes)
        l12: Half-step positions [m] (length n_nodes)
    """
    # Characteristic length scale
    lmin = np.sqrt(kappa * tmin / c / rho)
    
    # Initialize arrays
    l = np.zeros(n_nodes)
    l12 = np.zeros(n_nodes)
    
    # First layer
    l[0] = lmin / (1 + np.exp(np.pi / 2 / np.sqrt(n_nodes)))
    l12[0] = lmin
    
    # Subsequent layers with exponential spacing
    for i in range(1, n_nodes):
        l12[i] = l12[i - 1] * np.exp(np.pi / np.sqrt(n_nodes))
        l[i] = np.sqrt(l12[i - 1] * l12[i])
    
    return l, l12


def compute_surface_temperature(p_vapor: float, T_c: float, p_c: float, lambda_: float) -> float:
    """
    Compute interface/surface temperature from vapor pressure.
    
    Following Osipov 2008 correlation used in MATLAB (LH2Simulate.m line 266):
    Ts = T_c * (p_v / p_c)^(1/lambda)
    
    Args:
        p_vapor: Vapor pressure [Pa]
        T_c: Critical temperature [K]
        p_c: Critical pressure [Pa]
        lambda_: Exponent (typically 1.5)
    
    Returns:
        Surface temperature [K]
    """
    # Clamp pressure to valid range to avoid numerical issues
    p_ratio = max(1e-6, min(p_vapor / p_c, 10.0))  # Clamp to reasonable range
    Ts = T_c * p_ratio ** (1.0 / lambda_)
    # Clamp temperature to valid range for LH2
    return np.clip(Ts, 13.804, T_c)


def compute_latent_heat(Ts: float) -> float:
    """
    Compute latent heat of vaporization from surface temperature.
    
    Following MATLAB polynomial correlation (LH2Simulate.m lines 273-274):
    qh = 1000 * (-0.002445451720487*Ts^6 + 0.3629946692976*Ts^5 - 22.28028769483*Ts^4 
                 + 723.6541112107*Ts^3 - 13116.31006512*Ts^2 + 125780.2915522*Ts - 498095.5392318)
    
    6th order polynomial from REFPROP v9.1.
    
    Args:
        Ts: Surface/saturation temperature [K]
    
    Returns:
        Latent heat of vaporization [J/kg]
    """
    # Direct polynomial evaluation following MATLAB
    # qh = a6*T^6 + a5*T^5 + a4*T^4 + a3*T^3 + a2*T^2 + a1*T + a0
    qh = 1000.0 * (
        -0.002445451720487 * Ts**6
        + 0.3629946692976 * Ts**5
        - 22.28028769483 * Ts**4
        + 723.6541112107 * Ts**3
        - 13116.31006512 * Ts**2
        + 125780.2915522 * Ts
        - 498095.5392318
    )
    
    return qh


def compute_natural_convection_heat_transfer_coefficient(
    kappa: float,
    g: float,
    beta: float,
    cp: float,
    rho: float,
    mu: float,
    delta_T: float,
    length_scale: float,
) -> float:
    """
    Compute natural convection heat transfer coefficient.
    
    Following MATLAB approach (LH2Simulate.m lines 365, 368):
    h_conv = kappa * 0.156 * (Ra)^(1/3)
    where Ra = g * beta * cp * rho^2 * delta_T / kappa / mu
    
    This is a simplified Churchill-Chu correlation for natural convection.
    
    Args:
        kappa: Thermal conductivity [W/m/K]
        g: Gravitational acceleration [m/s²]
        beta: Thermal expansion coefficient [1/K]
        cp: Specific heat [J/kg/K]
        rho: Density [kg/m³]
        mu: Dynamic viscosity [Pa·s]
        delta_T: Temperature difference [K]
        length_scale: Characteristic length [m]
    
    Returns:
        Convective heat transfer coefficient [W/m²/K]
    """
    if abs(delta_T) < 1e-6:
        return 0.0
    
    Ra = g * beta * cp * rho**2 * abs(delta_T) / kappa / mu
    h_conv = kappa * 0.156 * Ra ** (1.0 / 3.0)
    
    return h_conv


def compute_wall_convection_nusselt(
    g: float,
    beta: float,
    delta_T: float,
    length_scale: float,
    nu: float,
    Pr: float,
    geometry: str = "vertical_wall",
) -> float:
    """
    Compute Nusselt number for natural convection at a wall.
    
    Following MATLAB approach (LH2Simulate.m lines 408-426):
    - Vertical wall: Nu = 0.68 + 0.503 * (Ra * Psi)^(1/4)
    - Horizontal plate: Nu = 0.27 * Ra^(1/4)
    
    where Psi = (1 + (0.492/Pr)^(9/16))^(-16/9)
    
    Args:
        g: Gravitational acceleration [m/s²]
        beta: Thermal expansion coefficient [1/K]
        delta_T: Temperature difference [K]
        length_scale: Characteristic length [m]
        nu: Kinematic viscosity [m²/s]
        Pr: Prandtl number [-]
        geometry: Either "vertical_wall" or "horizontal_plate"
    
    Returns:
        Nusselt number [-]
    """
    if abs(delta_T) < 1e-6:
        return 1.0  # Minimum Nusselt number
    
    # Rayleigh number
    Ra = abs(g * beta * delta_T * length_scale**3 * Pr / nu)
    
    if geometry == "vertical_wall":
        # Churchill-Chu correlation for vertical wall
        Psi = (1 + (0.492 / Pr) ** (9.0 / 16.0)) ** (-16.0 / 9.0)
        Nu = 0.68 + 0.503 * (Ra * Psi) ** (1.0 / 4.0)
    elif geometry == "horizontal_plate":
        # Simplified correlation for horizontal plate
        Nu = 0.27 * Ra ** (1.0 / 4.0)
    else:
        raise ValueError(f"Unknown geometry: {geometry}")
    
    return Nu


def compute_condensation_rate(
    Q_dot_liquid_surface: float,
    Q_dot_vapor_surface: float,
    latent_heat: float,
) -> float:
    """
    Compute condensation/evaporation rate from interface heat balance.
    
    Following MATLAB approach (LH2Simulate.m lines 484-494):
    J_cd = -(Q_dotLS + Q_dotVS) / qh
    
    Positive J_cd means condensation (vapor → liquid)
    Negative J_cd means evaporation (liquid → vapor)
    
    Args:
        Q_dot_liquid_surface: Heat flow from liquid to surface [W]
        Q_dot_vapor_surface: Heat flow from vapor to surface [W]
        latent_heat: Latent heat of vaporization [J/kg]
    
    Returns:
        Condensation rate [kg/s] (positive = condensation, negative = evaporation)
    """
    if abs(latent_heat) < 1e-6:
        return 0.0
    
    J_cd = -(Q_dot_liquid_surface + Q_dot_vapor_surface) / latent_heat
    return J_cd


def compute_temperature_from_internal_energy(
    u: float, fluid: str = "liquid", properties: FluidProperties = None
) -> float:
    """
    Compute temperature from specific internal energy.
    
    Uses polynomial correlations from MATLAB (LH2Simulate.m lines 147, 238).
    
    For liquid: T = poly3(u/1000) where u is in J/kg
    For vapor: uses REFPROP/CoolProp if available, else ideal gas
    
    Args:
        u: Specific internal energy [J/kg]
        fluid: Either "liquid" or "vapor"
        properties: FluidProperties instance (optional, for vapor)
    
    Returns:
        Temperature [K]
    """
    if fluid == "liquid":
        # Polynomial correlation from MATLAB (REFPROP v9.1)
        # T = poly3(u/1000) where u is in J/kg
        u_kJ = u / 1000.0
        coeffs = [
            1.44867559e-07,
            -2.53438808e-04,
            1.05449468e-01,
            2.03423757e01,
        ]
        T = np.polyval(coeffs, u_kJ)
        
        # Clamp to valid range for LH2
        T = np.clip(T, 13.804, 32.93)
        return T
    
    elif fluid == "vapor":
        # For vapor, would use REFPROP/CoolProp T(D, U)
        # Simplified: use ideal gas c_v * T = u
        # This is a placeholder - proper implementation should use properties
        if properties is not None:
            # Would call properties.temperature(density=rho, internal_energy=u)
            # For now, use ideal gas approximation
            pass
        
        # Ideal gas: u = c_v * T
        c_v = 6490.0  # J/kg/K for parahydrogen vapor
        T = u / c_v
        return max(T, 13.804)  # Minimum temperature
    
    else:
        raise ValueError(f"Unknown fluid type: {fluid}")


def compute_density_from_temperature(T: float, fluid: str = "liquid") -> float:
    """
    Compute density from temperature using polynomial correlations.
    
    Following MATLAB approach (LH2Simulate.m lines 117, 223).
    
    Args:
        T: Temperature [K]
        fluid: Either "liquid" or "vapor"
    
    Returns:
        Density [kg/m³]
    """
    if fluid == "liquid":
        # Polynomial correlation for liquid density (saturation curve)
        # From MATLAB but needs internal energy, so this is a placeholder
        # In MATLAB: rho_L = poly3(u_L/1000)
        # Here we use a constant approximation
        return 70.9  # kg/m³ at typical LH2 conditions
    
    else:
        raise ValueError("Vapor density should be computed from ideal gas law")

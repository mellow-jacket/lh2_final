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


def compute_interface_heat_flows(
    T_liquid_bulk: float,
    T_vapor_bulk: float,
    Ts: float,
    dTs_dt: float,
    rho_liquid: float,
    rho_vapor: float,
    kappa_liquid: float,
    kappa_vapor: float,
    cp_liquid: float,
    cp_vapor: float,
    cv_vapor: float,
    mu_liquid: float,
    mu_vapor: float,
    beta_liquid: float,
    beta_vapor: float,
    g: float,
    interface_area: float,
    n_liquid_nodes: int,
    n_vapor_nodes: int,
    tmin_liquid: float,
    tmin_vapor: float,
) -> Tuple[float, float, float, float]:
    """
    Compute heat transfer at the liquid-vapor interface.
    
    Following MATLAB approach (LH2Simulate.m lines 323-387):
    - Generates boundary layer grids for liquid and vapor
    - Computes conduction and convection heat transfer coefficients
    - Computes heat flows from liquid and vapor to interface
    - Returns both conduction and convection components
    
    Args:
        T_liquid_bulk: Bulk liquid temperature [K]
        T_vapor_bulk: Bulk vapor temperature [K]
        Ts: Interface/surface temperature [K]
        dTs_dt: Time derivative of surface temperature [K/s]
        rho_liquid: Liquid density [kg/m³]
        rho_vapor: Vapor density [kg/m³]
        kappa_liquid: Liquid thermal conductivity [W/m/K]
        kappa_vapor: Vapor thermal conductivity [W/m/K]
        cp_liquid: Liquid specific heat (constant pressure) [J/kg/K]
        cp_vapor: Vapor specific heat (constant pressure) [J/kg/K]
        cv_vapor: Vapor specific heat (constant volume) [J/kg/K]
        mu_liquid: Liquid dynamic viscosity [Pa·s]
        mu_vapor: Vapor dynamic viscosity [Pa·s]
        beta_liquid: Liquid thermal expansion coefficient [1/K]
        beta_vapor: Vapor thermal expansion coefficient [1/K]
        g: Gravitational acceleration [m/s²]
        interface_area: Interface area [m²]
        n_liquid_nodes: Number of liquid boundary layer nodes
        n_vapor_nodes: Number of vapor boundary layer nodes
        tmin_liquid: Minimum time scale for liquid [s]
        tmin_vapor: Minimum time scale for vapor [s]
    
    Returns:
        Tuple of (Q_dot_LS, Q_dot_VS, h_LS, h_VS):
            Q_dot_LS: Heat flow from liquid to surface [W]
            Q_dot_VS: Heat flow from vapor to surface [W]
            h_LS: Liquid-surface heat transfer coefficient [W/m²/K]
            h_VS: Vapor-surface heat transfer coefficient [W/m²/K]
    """
    # Generate boundary layer grids
    l_liquid, l12_liquid = generate_boundary_layer_grid(
        n_liquid_nodes, kappa_liquid, tmin_liquid, cp_liquid, rho_liquid
    )
    l_vapor, l12_vapor = generate_boundary_layer_grid(
        n_vapor_nodes, kappa_vapor, tmin_vapor, cp_vapor, rho_vapor
    )
    
    # Liquid-to-surface heat transfer
    h_LS_cond = kappa_liquid / l12_liquid[0]
    h_LS_conv = compute_natural_convection_heat_transfer_coefficient(
        kappa_liquid, g, beta_liquid, cp_liquid, rho_liquid, mu_liquid,
        abs(T_liquid_bulk - Ts), l_liquid[0]
    )
    
    # Heat flows (conduction and convection)
    Q_dot_LS_cond = h_LS_cond * interface_area * (T_liquid_bulk - Ts) - l_liquid[0] * cp_liquid * rho_liquid * dTs_dt
    Q_dot_LS_conv = h_LS_conv * interface_area * (T_liquid_bulk - Ts) * (T_liquid_bulk > Ts)
    
    # Choose max for heating, conduction only for cooling
    if Q_dot_LS_conv > 0:
        Q_dot_LS = max(Q_dot_LS_conv, Q_dot_LS_cond)
        h_LS = Q_dot_LS / (interface_area * (T_liquid_bulk - Ts) + 1e-10)
    else:
        Q_dot_LS = Q_dot_LS_cond
        h_LS = h_LS_cond
    
    # Vapor-to-surface heat transfer
    h_VS_cond = kappa_vapor / l12_vapor[0]
    h_VS_conv = compute_natural_convection_heat_transfer_coefficient(
        kappa_vapor, g, beta_vapor, cp_vapor, rho_vapor, mu_vapor,
        abs(Ts - T_vapor_bulk), l_vapor[0]
    )
    
    # Heat flows (conduction and convection)
    Q_dot_VS_cond = h_VS_cond * interface_area * (T_vapor_bulk - Ts) - l_vapor[0] * cv_vapor * rho_vapor * dTs_dt
    Q_dot_VS_conv = h_VS_conv * interface_area * (T_vapor_bulk - Ts) * (Ts > T_vapor_bulk)
    
    # Choose min for cooling, conduction only for heating
    if Q_dot_VS_conv < 0:
        Q_dot_VS = min(Q_dot_VS_conv, Q_dot_VS_cond)
        h_VS = Q_dot_VS / (interface_area * (T_vapor_bulk - Ts) + 1e-10)
    else:
        Q_dot_VS = Q_dot_VS_cond
        h_VS = h_VS_cond
    
    return Q_dot_LS, Q_dot_VS, h_LS, h_VS


def compute_wall_heat_transfer(
    T_wall: float,
    T_liquid_bulk: float,
    T_vapor_bulk: float,
    rho_liquid: float,
    rho_vapor: float,
    kappa_liquid: float,
    kappa_vapor: float,
    mu_liquid: float,
    mu_vapor: float,
    beta_liquid: float,
    beta_vapor: float,
    Pr_liquid: float,
    Pr_vapor: float,
    g: float,
    tank_height: float,
    liquid_height: float,
    tank_radius: float,
    tank_cross_section_area: float,
) -> Tuple[float, float]:
    """
    Compute heat transfer between tank wall and fluid phases.
    
    Following MATLAB approach (LH2Simulate.m lines 388-434):
    - Churchill-Chu correlation for natural convection at vertical walls
    - Simplified correlation for horizontal plates (bottom)
    - Separate treatment for liquid and vapor regions
    
    Args:
        T_wall: Wall temperature [K]
        T_liquid_bulk: Bulk liquid temperature [K]
        T_vapor_bulk: Bulk vapor temperature [K]
        rho_liquid: Liquid density [kg/m³]
        rho_vapor: Vapor density [kg/m³]
        kappa_liquid: Liquid thermal conductivity [W/m/K]
        kappa_vapor: Vapor thermal conductivity [W/m/K]
        mu_liquid: Liquid dynamic viscosity [Pa·s]
        mu_vapor: Vapor dynamic viscosity [Pa·s]
        beta_liquid: Liquid thermal expansion coefficient [1/K]
        beta_vapor: Vapor thermal expansion coefficient [1/K]
        Pr_liquid: Liquid Prandtl number [-]
        Pr_vapor: Vapor Prandtl number [-]
        g: Gravitational acceleration [m/s²]
        tank_height: Total tank height [m]
        liquid_height: Current liquid height [m]
        tank_radius: Tank radius [m]
        tank_cross_section_area: Tank cross-sectional area [m²]
    
    Returns:
        Tuple of (Q_dot_WL, Q_dot_WV):
            Q_dot_WL: Heat flow from wall to liquid [W]
            Q_dot_WV: Heat flow from wall to vapor [W]
    """
    vapor_height = tank_height - liquid_height
    
    # Wall-to-vapor heat transfer
    nu_vapor = mu_vapor / rho_vapor
    Nu_vapor = compute_wall_convection_nusselt(
        g, beta_vapor, T_wall - T_vapor_bulk, vapor_height, nu_vapor, Pr_vapor, "vertical_wall"
    )
    h_WV = Nu_vapor * kappa_vapor / max(vapor_height, 0.01)  # Avoid division by zero
    
    # Vapor region areas: top + side
    A_vapor_top = tank_cross_section_area
    A_vapor_side = 2 * np.pi * tank_radius * vapor_height
    Q_dot_WV = h_WV * (T_wall - T_vapor_bulk) * (A_vapor_top + A_vapor_side)
    
    # Wall-to-liquid heat transfer
    if liquid_height > 1e-6:
        nu_liquid = mu_liquid / rho_liquid
        
        # Side wall (vertical)
        Nu_liquid_side = compute_wall_convection_nusselt(
            g, beta_liquid, T_wall - T_liquid_bulk, liquid_height, nu_liquid, Pr_liquid, "vertical_wall"
        )
        h_WL_side = Nu_liquid_side * kappa_liquid / liquid_height
        A_liquid_side = 2 * np.pi * tank_radius * liquid_height
        
        # Bottom (horizontal plate)
        Nu_liquid_bottom = compute_wall_convection_nusselt(
            g, beta_liquid, T_wall - T_liquid_bulk, tank_radius / 2, nu_liquid, Pr_liquid, "horizontal_plate"
        )
        h_WL_bottom = Nu_liquid_bottom * kappa_liquid / (tank_radius / 2)
        A_liquid_bottom = tank_cross_section_area
        
        Q_dot_WL = (T_wall - T_liquid_bulk) * (h_WL_bottom * A_liquid_bottom + h_WL_side * A_liquid_side)
    else:
        Q_dot_WL = 0.0
    
    return Q_dot_WL, Q_dot_WV


def compute_environmental_heat_leak(
    liquid_volume: float, correlation_coeffs: Tuple[float, float, float] = None
) -> float:
    """
    Compute environmental heat leak to tank from surroundings.
    
    Following MATLAB approach (LH2Simulate.m line 645):
    Q_dot_EW = a*V_L^2 + b*V_L + c
    
    This is a correlation for a specific tank geometry and insulation.
    Default coefficients are for LLNL 3,300 gallon vertical Dewar.
    
    Args:
        liquid_volume: Liquid volume [m³]
        correlation_coeffs: Tuple of (a, b, c) for quadratic correlation
                           Default is for LLNL 3,300 gallon Dewar
    
    Returns:
        Heat leak from environment [W]
    """
    if correlation_coeffs is None:
        # Default LLNL Dewar correlation (LH2Simulate.m line 645)
        a = -7.462776654302e-02  # W/m⁶
        b = 4.445867251697e00  # W/m³
        c = 3.108170556297e01  # W
    else:
        a, b, c = correlation_coeffs
    
    Q_dot_EW = a * liquid_volume**2 + b * liquid_volume + c
    return max(Q_dot_EW, 0.0)  # Heat leak is always positive (into tank)


def compute_wall_temperature_derivative(
    T_wall: float,
    Q_dot_environment: float,
    Q_dot_to_liquid: float,
    Q_dot_to_vapor: float,
    wall_mass: float,
    wall_specific_heat: float,
) -> float:
    """
    Compute time derivative of wall temperature.
    
    Energy balance on wall thermal mass:
    M_w * c_w * dT_w/dt = Q_dot_env - Q_dot_WL - Q_dot_WV
    
    Args:
        T_wall: Current wall temperature [K]
        Q_dot_environment: Heat flow from environment to wall [W]
        Q_dot_to_liquid: Heat flow from wall to liquid [W]
        Q_dot_to_vapor: Heat flow from wall to vapor [W]
        wall_mass: Wall thermal mass [kg]
        wall_specific_heat: Wall specific heat [J/kg/K]
    
    Returns:
        dT_wall/dt [K/s]
    """
    if wall_mass < 1e-6 or wall_specific_heat < 1e-6:
        return 0.0
    
    dT_wall_dt = (Q_dot_environment - Q_dot_to_liquid - Q_dot_to_vapor) / (wall_mass * wall_specific_heat)
    return dT_wall_dt

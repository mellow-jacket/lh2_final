"""
Flow Module

Provides fluid flow calculations including:
- Choked and non-choked gas flow through orifices
- Valve flow models
- Mass flow rate calculations
"""

import numpy as np


def gas_flow(CA, rho, P1, P2, gamma):
    """
    Calculate mass flow rate through an orifice for compressible gas flow.

    Handles both choked and non-choked flow conditions based on pressure ratio
    and gas properties. Based on isentropic flow relations.

    Parameters
    ----------
    CA : float
        Effective orifice area (C_d * A) [m²]
    rho : float
        Gas density at upstream conditions [kg/m³]
    P1 : float
        Upstream pressure [Pa]
    P2 : float
        Downstream pressure [Pa]
    gamma : float
        Ratio of specific heats (Cp/Cv) [-]

    Returns
    -------
    mdot : float
        Mass flow rate [kg/s]
        Positive for flow from P1 to P2

    Notes
    -----
    The function automatically handles flow reversal when P1 < P2.

    Choked flow occurs when the pressure ratio reaches the critical value:
    P2/P1 < ((gamma+1)/2)^(gamma/(gamma-1))

    References
    ----------
    Original MATLAB implementation: LLNL_model/gasFlow.m
    """
    # Handle flow reversal
    if P1 < P2:
        return -gas_flow(CA, rho, P2, P1, gamma)

    # Avoid division by zero
    if P1 <= 0 or CA <= 0:
        return 0.0

    # Critical pressure ratio for choked flow
    pressure_ratio = P2 / P1
    critical_ratio = ((gamma + 1) / 2) ** (gamma / (gamma - 1))

    if pressure_ratio < (1.0 / critical_ratio):
        # Choked flow
        # mdot = CA * sqrt(gamma * rho * P1 * (2/(gamma+1))^((gamma+1)/(gamma-1)))
        exponent = (gamma + 1) / (gamma - 1)
        term = (2.0 / (gamma + 1)) ** exponent
        mdot = CA * np.sqrt(gamma * rho * P1 * term)
    else:
        # Non-choked (subsonic) flow
        # mdot = CA * sqrt(2 * rho * P1 * (gamma/(gamma-1)) *
        #                  [(P2/P1)^(2/gamma) - (P2/P1)^((gamma+1)/gamma)])
        term1 = pressure_ratio ** (2.0 / gamma)
        term2 = pressure_ratio ** ((gamma + 1) / gamma)
        factor = (gamma / (gamma - 1)) * (term1 - term2)

        if factor < 0:
            # Numerical issue - should not happen for valid pressure ratios
            mdot = 0.0
        else:
            mdot = CA * np.sqrt(2.0 * rho * P1 * factor)

    return mdot


def directed_sqrt(x):
    """
    Compute signed square root: sign(x) * sqrt(abs(x)).

    Helper function used in flow calculations to preserve flow direction.

    Parameters
    ----------
    x : float or array
        Input value(s)

    Returns
    -------
    result : float or array
        Signed square root
    """
    return np.sign(x) * np.sqrt(np.abs(x))


def valve_flow_coefficient(A_valve, lambda_valve, C_d=0.6):
    """
    Calculate effective valve flow area.

    Parameters
    ----------
    A_valve : float
        Full valve area [m²]
    lambda_valve : float
        Valve opening fraction [0-1]
    C_d : float, optional
        Discharge coefficient, default 0.6

    Returns
    -------
    CA : float
        Effective flow area [m²]
    """
    return C_d * A_valve * lambda_valve


def vent_flow_rate(P_tank, P_ambient, rho_vapor, A_vent, gamma, C_d=0.6):
    """
    Calculate vent flow rate from tank to ambient.

    Parameters
    ----------
    P_tank : float
        Tank pressure [Pa]
    P_ambient : float
        Ambient pressure [Pa]
    rho_vapor : float
        Vapor density in tank [kg/m³]
    A_vent : float
        Vent area [m²]
    gamma : float
        Ratio of specific heats [-]
    C_d : float, optional
        Discharge coefficient, default 0.6

    Returns
    -------
    mdot_vent : float
        Vent mass flow rate [kg/s]
    """
    CA = C_d * A_vent
    return gas_flow(CA, rho_vapor, P_tank, P_ambient, gamma)

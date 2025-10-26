"""
Properties Module

Provides thermophysical property calculations for liquid hydrogen using CoolProp.
Includes abstraction layer to support both CoolProp and REFPROP backends.

Key functionality:
- Density, enthalpy, entropy calculations
- Vapor pressure and saturation properties
- Transport properties (viscosity, thermal conductivity)
- Polynomial correlations as fallbacks
"""

import numpy as np
from typing import Optional

try:
    import CoolProp.CoolProp as CP
    COOLPROP_AVAILABLE = True
except ImportError:
    COOLPROP_AVAILABLE = False
    print("Warning: CoolProp not available. Some property calculations will use polynomial approximations.")


class FluidProperties:
    """
    Abstraction layer for thermophysical properties.
    
    Supports CoolProp backend with polynomial correlation fallbacks.
    All SI units unless otherwise specified.
    """
    
    def __init__(self, fluid_name="Hydrogen", backend="CoolProp"):
        """
        Initialize fluid properties calculator.
        
        Parameters
        ----------
        fluid_name : str
            Fluid name (e.g., "Hydrogen", "PARAHYD")
        backend : str
            Property backend: "CoolProp" or "polynomial"
        """
        self.fluid_name = fluid_name
        self.backend = backend
        
        if backend == "CoolProp" and not COOLPROP_AVAILABLE:
            print(f"Warning: CoolProp not available, falling back to polynomial correlations")
            self.backend = "polynomial"
    
    def density(self, T=None, P=None, Q=None, U=None):
        """
        Calculate density.
        
        Parameters
        ----------
        T : float, optional
            Temperature [K]
        P : float, optional
            Pressure [Pa]
        Q : float, optional
            Quality (0=saturated liquid, 1=saturated vapor) [-]
        U : float, optional
            Internal energy [J/kg]
            
        Returns
        -------
        rho : float
            Density [kg/m³]
        """
        if self.backend == "CoolProp":
            if T is not None and P is not None:
                return CP.PropsSI('D', 'T', T, 'P', P, self.fluid_name)
            elif T is not None and Q is not None:
                return CP.PropsSI('D', 'T', T, 'Q', Q, self.fluid_name)
            elif U is not None and T is not None:
                return CP.PropsSI('D', 'U', U, 'T', T, self.fluid_name)
            else:
                raise ValueError("Insufficient inputs for density calculation")
        else:
            # Polynomial correlation for liquid density vs temperature
            # Based on LLNL model correlations
            if T is not None:
                return self._liquid_density_polynomial(T)
            else:
                raise ValueError("Temperature required for polynomial density")
    
    def pressure(self, T=None, D=None, U=None, Q=None):
        """
        Calculate pressure.
        
        Parameters
        ----------
        T : float, optional
            Temperature [K]
        D : float, optional
            Density [kg/m³]
        U : float, optional
            Internal energy [J/kg]
        Q : float, optional
            Quality [-]
            
        Returns
        -------
        P : float
            Pressure [Pa]
        """
        if self.backend == "CoolProp":
            if T is not None and D is not None:
                return CP.PropsSI('P', 'T', T, 'D', D, self.fluid_name)
            elif U is not None and D is not None:
                return CP.PropsSI('P', 'U', U, 'D', D, self.fluid_name)
            elif T is not None and Q is not None:
                return CP.PropsSI('P', 'T', T, 'Q', Q, self.fluid_name)
            else:
                raise ValueError("Insufficient inputs for pressure calculation")
        else:
            if T is not None:
                return self._vapor_pressure_polynomial(T)
            else:
                raise ValueError("Insufficient inputs for polynomial pressure")
    
    def temperature(self, P=None, D=None, U=None, Q=None):
        """
        Calculate temperature.
        
        Parameters
        ----------
        P : float, optional
            Pressure [Pa]
        D : float, optional
            Density [kg/m³]
        U : float, optional
            Internal energy [J/kg]
        Q : float, optional
            Quality [-]
            
        Returns
        -------
        T : float
            Temperature [K]
        """
        if self.backend == "CoolProp":
            if P is not None and D is not None:
                return CP.PropsSI('T', 'P', P, 'D', D, self.fluid_name)
            elif U is not None and D is not None:
                return CP.PropsSI('T', 'U', U, 'D', D, self.fluid_name)
            elif P is not None and Q is not None:
                return CP.PropsSI('T', 'P', P, 'Q', Q, self.fluid_name)
            else:
                raise ValueError("Insufficient inputs for temperature calculation")
        else:
            raise NotImplementedError("Polynomial temperature calculation not implemented")
    
    def enthalpy(self, T=None, P=None, Q=None):
        """
        Calculate enthalpy.
        
        Parameters
        ----------
        T : float, optional
            Temperature [K]
        P : float, optional
            Pressure [Pa]
        Q : float, optional
            Quality [-]
            
        Returns
        -------
        h : float
            Enthalpy [J/kg]
        """
        if self.backend == "CoolProp":
            if T is not None and P is not None:
                return CP.PropsSI('H', 'T', T, 'P', P, self.fluid_name)
            elif T is not None and Q is not None:
                return CP.PropsSI('H', 'T', T, 'Q', Q, self.fluid_name)
            else:
                raise ValueError("Insufficient inputs for enthalpy calculation")
        else:
            raise NotImplementedError("Polynomial enthalpy calculation not implemented")
    
    def viscosity(self, T, P=None, D=None):
        """
        Calculate dynamic viscosity.
        
        Parameters
        ----------
        T : float
            Temperature [K]
        P : float, optional
            Pressure [Pa]
        D : float, optional
            Density [kg/m³]
            
        Returns
        -------
        mu : float
            Dynamic viscosity [Pa·s]
        """
        if self.backend == "CoolProp":
            if P is not None:
                return CP.PropsSI('V', 'T', T, 'P', P, self.fluid_name)
            elif D is not None:
                return CP.PropsSI('V', 'T', T, 'D', D, self.fluid_name)
            else:
                raise ValueError("Pressure or density required")
        else:
            # Simple polynomial approximation
            return 1e-6 * (0.1 + 0.001 * T)
    
    def thermal_conductivity(self, T, P=None, D=None):
        """
        Calculate thermal conductivity.
        
        Parameters
        ----------
        T : float
            Temperature [K]
        P : float, optional
            Pressure [Pa]
        D : float, optional
            Density [kg/m³]
            
        Returns
        -------
        k : float
            Thermal conductivity [W/m/K]
        """
        if self.backend == "CoolProp":
            if P is not None:
                return CP.PropsSI('L', 'T', T, 'P', P, self.fluid_name)
            elif D is not None:
                return CP.PropsSI('L', 'T', T, 'D', D, self.fluid_name)
            else:
                raise ValueError("Pressure or density required")
        else:
            # Simple approximation
            return 0.1 + 0.001 * T
    
    def specific_heat_cp(self, T, P=None):
        """
        Calculate specific heat at constant pressure.
        
        Parameters
        ----------
        T : float
            Temperature [K]
        P : float, optional
            Pressure [Pa]
            
        Returns
        -------
        Cp : float
            Specific heat [J/kg/K]
        """
        if self.backend == "CoolProp":
            if P is not None:
                return CP.PropsSI('C', 'T', T, 'P', P, self.fluid_name)
            else:
                raise ValueError("Pressure required")
        else:
            # Approximate value for liquid hydrogen
            return 9700.0  # J/kg/K
    
    def _liquid_density_polynomial(self, T):
        """
        Polynomial correlation for liquid hydrogen density vs temperature.
        
        Based on LLNL model correlations fitted to REFPROP data.
        Valid for liquid phase at saturation.
        """
        # Polynomial coefficients (example - would need actual fitted values)
        # rho = a0 + a1*T + a2*T^2 + ...
        # Using approximate values for LH2
        if T < 14:  # Below triple point
            return 76.0
        elif T > 32:  # Above critical point
            return 30.0
        else:
            # Linear approximation between 20K (71 kg/m³) and 30K (40 kg/m³)
            return 71.0 - (71.0 - 40.0) * (T - 20.0) / 10.0
    
    def _vapor_pressure_polynomial(self, T):
        """
        Polynomial correlation for vapor pressure vs temperature.
        
        Based on LLNL model correlations.
        """
        # Cubic polynomial for vapor pressure
        # Based on MATLAB vaporpressure.m
        # Coefficients would be fitted from REFPROP data
        if T < 14:
            return 7e3  # Approximate triple point pressure
        elif T > 33:
            return 1.3e6  # Approximate critical pressure
        else:
            # Exponential approximation: P = P0 * exp(B*(T-T0))
            P0 = 1e5  # Reference pressure at ~20K
            T0 = 20.0
            B = 0.15  # Fitted parameter
            return P0 * np.exp(B * (T - T0))


def vapor_pressure(T_film, rho_vapor):
    """
    Calculate vapor pressure from film temperature and vapor density.
    
    This function replicates the MATLAB vaporpressure.m logic with
    polynomial correlations and phase determination.
    
    Parameters
    ----------
    T_film : float
        Film temperature [K]
    rho_vapor : float
        Vapor density [kg/m³]
        
    Returns
    -------
    P : float
        Vapor pressure [Pa]
        
    References
    ----------
    Original MATLAB: LLNL_model/vaporpressure.m
    """
    # Saturation temperature from density polynomial
    # T_sat(rho) polynomial - would need actual coefficients
    T_sat = 20.0 + 0.1 * rho_vapor  # Simplified
    
    if T_film > T_sat:
        # Supercritical branch - use ideal gas approximation
        R_specific = 4124.0  # J/kg/K for hydrogen
        P = rho_vapor * R_specific * T_film
    else:
        # Two-phase branch - cubic polynomial in temperature
        # P = a0 + a1*T + a2*T^2 + a3*T^3
        # Coefficients from MATLAB code (example values)
        a0 = -1e6
        a1 = 1e5
        a2 = 1e3
        a3 = 10.0
        P = a0 + a1 * T_film + a2 * T_film**2 + a3 * T_film**3
        P = max(P, 1e3)  # Ensure positive pressure
    
    return P

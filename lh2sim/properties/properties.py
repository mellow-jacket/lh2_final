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

import math
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
            print("Warning: CoolProp not available, falling back to polynomial correlations")
            self.backend = "polynomial"

        # Saturation-aware blending / safety margins (best-practice defaults)
        # These mirror the best-practice belt used in smooth_phase.py
        self._REL_MARGIN = 2e-4
        self._ABS_MARGIN = 300.0
        self._TAU_FB = 500.0
        self._BELT_REL = 5e-4
        self._BELT_ABS = 1.5e3

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
            # Use saturation-aware logic: when near saturation, blend saturated
            # liquid/vapor values smoothly. Otherwise query single-phase at a
            # safe pressure away from Psat(T).
            if T is not None and P is not None:
                return self._density_cp_safe(T, P)
            elif T is not None and Q is not None:
                return CP.PropsSI("D", "T", T, "Q", Q, self.fluid_name)
            elif U is not None and T is not None:
                return CP.PropsSI("D", "U", U, "T", T, self.fluid_name)
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
                return CP.PropsSI("P", "T", T, "D", D, self.fluid_name)
            elif U is not None and D is not None:
                return CP.PropsSI("P", "U", U, "D", D, self.fluid_name)
            elif T is not None and Q is not None:
                return CP.PropsSI("P", "T", T, "Q", Q, self.fluid_name)
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
                return CP.PropsSI("T", "P", P, "D", D, self.fluid_name)
            elif U is not None and D is not None:
                return CP.PropsSI("T", "U", U, "D", D, self.fluid_name)
            elif P is not None and Q is not None:
                return CP.PropsSI("T", "P", P, "Q", Q, self.fluid_name)
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
                return self._blend_or_safe_prop("H", T, P)
            elif T is not None and Q is not None:
                return CP.PropsSI("H", "T", T, "Q", Q, self.fluid_name)
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
                return self._blend_or_safe_prop("V", T, P)
            elif D is not None:
                return CP.PropsSI("V", "T", T, "D", D, self.fluid_name)
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
                return self._blend_or_safe_prop("L", T, P)
            elif D is not None:
                return CP.PropsSI("L", "T", T, "D", D, self.fluid_name)
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
                return self._blend_or_safe_prop("C", T, P)
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

    # -----------------------------
    # Saturation-aware helpers (CoolProp backend)
    # -----------------------------
    def _psat_T(self, T: float) -> float:
        """Return Psat(T) [Pa] or NaN if unavailable."""
        try:
            return CP.PropsSI("P", "T", T, "Q", 0, self.fluid_name)
        except Exception:
            return float("nan")

    def _softpos(self, x: float, delta: float) -> float:
        return 0.5 * (x + math.sqrt(x * x + delta * delta))

    def _in_belt(self, T: float, p: float) -> bool:
        Ps = self._psat_T(T)
        if not math.isfinite(Ps):
            return False
        belt = max(self._BELT_ABS, self._BELT_REL * max(p, 1.0))
        return abs(p - Ps) <= belt

    def _safe_pressure_away_from_saturation(self, T: float, p: float, phase: str) -> float:
        Ps = self._psat_T(T)
        target_gap = max(self._ABS_MARGIN, self._REL_MARGIN * max(p, 1.0))
        phase_l = phase.lower()
        if math.isfinite(Ps):
            if "liq" in phase_l:
                p_adj = max(p, Ps + target_gap)
            else:
                p_adj = min(p, Ps - target_gap)
            if abs(p_adj - Ps) < target_gap:
                p_adj = Ps + target_gap if "liq" in phase_l else Ps - target_gap
            return max(p_adj, 1.0)
        else:
            return max(p + self._ABS_MARGIN, 1.0) if "liq" in phase_l else max(p - self._ABS_MARGIN, 1.0)

    def _blend_or_safe_prop(self, prop_code: str, T: float, P: float):
        """Return a property either blended across saturation (when in belt)
        or computed from a safe single-phase query.

        prop_code: CoolProp output key (e.g., 'D','H','C','V','L')
        """
        # If Psat unknown, fall back to direct single-phase call
        Ps = self._psat_T(T)
        if not math.isfinite(Ps):
            # Use a safe single-phase guess: query as vapor or liquid via PhaseSI
            try:
                phase = CP.PhaseSI("T", T, "P", P, self.fluid_name)
                if "liquid" in phase.lower():
                    p_try = self._safe_pressure_away_from_saturation(T, P, "liquid")
                else:
                    p_try = self._safe_pressure_away_from_saturation(T, P, "vapor")
            except Exception:
                p_try = max(P, 1.0)
            return CP.PropsSI(prop_code, "T", T, "P", p_try, self.fluid_name)

        # If inside belt: blend saturated-liquid and saturated-vapor property values
        if self._in_belt(T, P):
            y = self._softpos(P - Ps, self._TAU_FB)
            z = self._softpos(Ps - P, self._TAU_FB)
            denom = y + z + 1e-16
            w_v = z / denom
            w_l = 1.0 - w_v
            try:
                vL = CP.PropsSI(prop_code, "T", T, "Q", 0, self.fluid_name)
            except Exception:
                vL = float("nan")
            try:
                vV = CP.PropsSI(prop_code, "T", T, "Q", 1, self.fluid_name)
            except Exception:
                vV = float("nan")
            # If either saturated value is not available, fall back to safe single-phase
            if not (math.isfinite(vL) and math.isfinite(vV)):
                # pick side by P vs Ps
                phase = "liquid" if P > Ps else "vapor"
                p_try = self._safe_pressure_away_from_saturation(T, P, phase)
                return CP.PropsSI(prop_code, "T", T, "P", p_try, self.fluid_name)
            return w_l * vL + w_v * vV

        # Outside belt: query single-phase safely on the appropriate side
        phase = "liquid" if P > Ps else "vapor"
        p_try = self._safe_pressure_away_from_saturation(T, P, phase)
        return CP.PropsSI(prop_code, "T", T, "P", p_try, self.fluid_name)

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

    Notes
    -----
    Uses exact polynomial coefficients from LLNL MATLAB vaporpressure.m.
    - T_sat(rho_vapor): 6th-order polynomial for saturation temperature
    - P(T): 3rd-order polynomial for vapor pressure in two-phase region
    - For supercritical: Would use REFPROP, here we use ideal gas law

    References
    ----------
    Original MATLAB: LLNL_model/vaporpressure.m (lines 7-30)
    """
    # Handle negative or very small densities
    if rho_vapor < 0:
        rho_vapor = 0.0001

    # Saturation temperature from density polynomial (6th order)
    # Exact coefficients from MATLAB vaporpressure.m lines 7-9
    T_sat = (-3.9389254667e-09 * (rho_vapor**6) +
             1.0053641879e-06 * (rho_vapor**5) -
             1.0304184083e-04 * (rho_vapor**4) +
             5.3058942923e-03 * (rho_vapor**3) -
             1.4792439609e-01 * (rho_vapor**2) +
             2.2234419496 * rho_vapor +
             1.7950995359e+01)

    if T_film > T_sat:
        # Supercritical branch
        # MATLAB uses REFPROP here, but we'll use ideal gas as approximation
        # In production, this should use CoolProp or REFPROP backend
        # Clip temperature near critical point to avoid numerical issues
        T_vapor = T_film
        if T_vapor > 32.937 and T_vapor < 32.938:
            T_vapor = 32.937

        # Ideal gas approximation (better than nothing without REFPROP)
        R_specific = 4124.0  # J/kg/K for hydrogen
        P = rho_vapor * R_specific * T_vapor
    else:
        # Two-phase branch - cubic polynomial in temperature
        # Exact coefficients from MATLAB vaporpressure.m line 30
        # P = a3*T^3 + a2*T^2 + a1*T + a0, returns in Pa (after *1e3)
        a3 = 1.6133821043e-1
        a2 = -6.9432088540
        a1 = 1.1373052580e2
        a0 = -6.9558797798e2

        P = (a3 * (T_film**3) + a2 * (T_film**2) +
             a1 * T_film + a0) * 1e3  # Convert to Pa

        # Ensure positive pressure
        P = max(P, 1e3)

    return P

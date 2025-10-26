"""
Control Module

Provides control logic for liquid hydrogen transfer operations:
- Pressure-driven control strategies
- Pump-driven control strategies
- Vent control with hysteresis
- Fill regime management
"""

import numpy as np
from typing import Dict, Any
from dataclasses import dataclass


@dataclass
class ControlOutputs:
    """Container for control outputs."""

    lambda_E: float  # Transfer valve opening fraction [0-1]
    lambda_V: float  # Vaporizer valve opening fraction [0-1]
    ST_vent_state: int  # Source tank vent state (0=closed, 1=open)
    ET_vent_state: int  # End tank vent state (0=closed, 1=open)


class PressureDrivenControl:
    """
    Pressure-driven control for LH2 transfer.

    Based on paper_model/LH2Control.m with improvements.
    Controls transfer valve, vaporizer valve, and vent valves based on
    liquid level and tank pressures.

    Parameters
    ----------
    params : dict
        Dictionary of control parameters including:
        - H: End tank height [m]
        - p_ST_slow, p_ST_fast, p_ST_final: ST pressure setpoints [Pa]
        - p_ET_low, p_ET_high, p_ET_final: ET pressure limits [Pa]
    """

    def __init__(self, params: Dict[str, Any]):
        """Initialize pressure-driven controller."""
        self.params = params

        # Internal state for hysteresis
        self.ST_vent_state = 0
        self.ET_vent_state = 0
        self.vap_valve_state = 0.0

    def compute_control(
        self,
        h_L2: float,
        p_ST: float,
        p_ET: float,
        ET_fill_complete: bool = False,
        ST_vent_complete: bool = False,
    ) -> ControlOutputs:
        """
        Compute control outputs based on current state.

        Parameters
        ----------
        h_L2 : float
            Liquid height in end tank [m]
        p_ST : float
            Source tank (trailer) pressure [Pa]
        p_ET : float
            End tank (dewar) pressure [Pa]
        ET_fill_complete : bool
            Flag indicating end tank is full
        ST_vent_complete : bool
            Flag indicating source tank venting complete

        Returns
        -------
        outputs : ControlOutputs
            Control valve positions and vent states
        """
        P = self.params
        H = P["H"]  # End tank height

        # Check if explicit ST vent thresholds are provided (for single-tank venting)
        has_explicit_ST_vent = "p_ST_vent_open" in P and "p_ST_vent_close" in P

        # Determine fill regime based on liquid height
        h_fraction = h_L2 / H

        if h_fraction < 0.05:
            # Slow fill regime
            lambda_E = 0.5 * (1 - float(ET_fill_complete))
            lambda_V = (1 - float(ET_fill_complete)) * self._get_vaporizer_valve_state(p_ST, P["p_ST_slow"])
            if has_explicit_ST_vent:
                ST_vent_state = self._get_ST_vent_state_explicit(p_ST)
            else:
                ST_vent_state = self._get_ST_vent_state(p_ST, P["p_ST_slow"])

        elif h_fraction < 0.70:
            # Fast fill regime
            lambda_E = 1.0 * (1 - float(ET_fill_complete))
            lambda_V = (1 - float(ET_fill_complete)) * self._get_vaporizer_valve_state(p_ST, P["p_ST_fast"])
            if has_explicit_ST_vent:
                ST_vent_state = self._get_ST_vent_state_explicit(p_ST)
            else:
                ST_vent_state = self._get_ST_vent_state(p_ST, P["p_ST_fast"])

        elif h_fraction < 0.85:
            # Reduced fast fill regime
            lambda_E = 0.8 * (1 - float(ET_fill_complete))
            lambda_V = (1 - float(ET_fill_complete)) * self._get_vaporizer_valve_state(p_ST, P["p_ST_slow"])
            if ET_fill_complete:
                ST_vent_state = 1 if p_ST > P["p_ST_final"] else 0
            else:
                if has_explicit_ST_vent:
                    ST_vent_state = self._get_ST_vent_state_explicit(p_ST)
                else:
                    ST_vent_state = self._get_ST_vent_state(p_ST, P["p_ST_slow"])

        else:
            # Topping regime
            lambda_E = 0.7 * (1 - float(ET_fill_complete))
            lambda_V = (1 - float(ET_fill_complete)) * self._get_vaporizer_valve_state(p_ST, P["p_ST_slow"])
            if ET_fill_complete:
                ST_vent_state = int(not ST_vent_complete)
            else:
                ST_vent_state = 0

        # ET vent control
        if ET_fill_complete:
            ET_vent_state = 1 if p_ET > P["p_ET_final"] else 0
        else:
            ET_vent_state = self._get_ET_vent_state(p_ET)

        # Override ST vent when fill is complete
        if ET_fill_complete:
            ST_vent_state = 1 if p_ST > P["p_ST_final"] else 0

        # Update internal state
        self.ST_vent_state = ST_vent_state
        self.ET_vent_state = ET_vent_state

        return ControlOutputs(
            lambda_E=lambda_E,
            lambda_V=lambda_V,
            ST_vent_state=ST_vent_state,
            ET_vent_state=ET_vent_state,
        )

    def _get_ST_vent_state(self, p_ST: float, threshold: float) -> int:
        """
        Determine source tank vent state with hysteresis.

        Parameters
        ----------
        p_ST : float
            Source tank pressure [Pa]
        threshold : float
            Pressure threshold [Pa]

        Returns
        -------
        state : int
            Vent state (0=closed, 1=open)
        """
        # Hysteresis: ±10% deadband around threshold
        # Open vent when pressure exceeds upper threshold
        # Close vent when pressure drops below lower threshold
        upper_threshold = threshold + 0.1 * threshold  # threshold * 1.1
        lower_threshold = threshold - 0.1 * threshold  # threshold * 0.9
        
        if p_ST > upper_threshold:
            # Pressure above upper limit - open vent
            return 1
        elif p_ST < lower_threshold:
            # Pressure below lower limit - close vent
            return 0
        else:
            # Within deadband - maintain current state (hysteresis)
            return self.ST_vent_state

    def _get_ET_vent_state(self, p_ET: float) -> int:
        """
        Determine end tank vent state with hysteresis.

        Parameters
        ----------
        p_ET : float
            End tank pressure [Pa]

        Returns
        -------
        state : int
            Vent state (0=closed, 1=open)
        """
        P = self.params

        if p_ET < P["p_ET_low"]:
            # Below low limit - close vent
            return 0
        elif p_ET > P["p_ET_high"]:
            # Above high limit - open vent
            return 1
        else:
            # Within range - maintain current state
            return self.ET_vent_state

    def _get_ST_vent_state_explicit(self, p_ST: float) -> int:
        """
        Determine supply tank vent state with explicit open/close thresholds.
        
        Used for single-tank venting scenarios with explicit pressure thresholds.
        
        Parameters
        ----------
        p_ST : float
            Supply tank pressure [Pa]
            
        Returns
        -------
        state : int
            Vent state (0=closed, 1=open)
        """
        P = self.params
        
        if p_ST > P["p_ST_vent_open"]:
            # Pressure above open threshold - open vent
            return 1
        elif p_ST < P["p_ST_vent_close"]:
            # Pressure below close threshold - close vent
            return 0
        else:
            # Within hysteresis band - maintain current state
            return self.ST_vent_state

    def _get_vaporizer_valve_state(self, p_ST: float, p_set: float) -> float:
        """
        Determine vaporizer valve opening to maintain ST pressure.

        Proportional control to maintain p_ST at p_set.

        Parameters
        ----------
        p_ST : float
            Source tank pressure [Pa]
        p_set : float
            Desired pressure setpoint [Pa]

        Returns
        -------
        state : float
            Valve opening fraction [0-1]
        """
        # ±2% deadband
        if p_ST < p_set - 0.02 * p_set:
            # Pressure too low - open valve proportionally
            state = 10.0 * (p_set - p_ST) / p_set
            state = np.clip(state, 0.0, 1.0)
        elif p_ST > p_set + 0.02 * p_set:
            # Pressure too high - close valve
            state = 0.0
        else:
            # Within deadband - maintain current state
            state = self.vap_valve_state

        self.vap_valve_state = state
        return state


class PumpDrivenControl:
    """
    Pump-driven control for LH2 transfer.

    Based on paper_model/LH2Control_pump.m.
    Controls pump speed, transfer valve, and vent valves.

    Parameters
    ----------
    params : dict
        Dictionary of control parameters including:
        - H: End tank height [m]
        - PumpMassTransferSlow, PumpMassTransferFast: Pump flow rates [kg/s]
        - p_ET_low, p_ET_high, p_ET_final: ET pressure limits [Pa]
    """

    def __init__(self, params: Dict[str, Any]):
        """Initialize pump-driven controller."""
        self.params = params

        # Internal state
        self.ET_vent_state = 0
        self.ST_vent_state = 0

    def compute_control(
        self,
        h_L2: float,
        p_ET: float,
        ET_fill_complete: bool = False,
        p_ST: float = None,  # Optional for ST vent control
    ) -> ControlOutputs:
        """
        Compute control outputs for pump-driven mode.

        Parameters
        ----------
        h_L2 : float
            Liquid height in end tank [m]
        p_ET : float
            End tank pressure [Pa]
        ET_fill_complete : bool
            Flag indicating end tank is full
        p_ST : float, optional
            Supply tank pressure [Pa] (for ST vent control if enabled)

        Returns
        -------
        outputs : ControlOutputs
            Control valve positions and vent states
        """
        P = self.params
        H = P["H"]

        # Determine fill regime
        h_fraction = h_L2 / H

        if h_fraction < 0.85:
            # Fast fill with pump
            # Scale transfer valve by ratio of pump speeds
            lambda_E = (P["PumpMassTransferSlow"] / P["PumpMassTransferFast"]) * (1 - float(ET_fill_complete))
        else:
            # Topping regime
            lambda_E = 0.7 * (1 - float(ET_fill_complete))

        # No vaporizer in pump mode
        lambda_V = 0.0

        # ST vent control (if p_ST provided and thresholds set)
        if p_ST is not None and "p_ST_vent_open" in P and "p_ST_vent_close" in P:
            ST_vent_state = self._get_ST_vent_state(p_ST)
        else:
            # Default: no ST vent in pump mode
            ST_vent_state = 0

        # ET vent control
        if ET_fill_complete:
            ET_vent_state = 1 if p_ET > P["p_ET_final"] else 0
        else:
            ET_vent_state = self._get_ET_vent_state(p_ET)

        self.ET_vent_state = ET_vent_state
        self.ST_vent_state = ST_vent_state

        return ControlOutputs(
            lambda_E=lambda_E,
            lambda_V=lambda_V,
            ST_vent_state=ST_vent_state,
            ET_vent_state=ET_vent_state,
        )
    
    def _get_ST_vent_state(self, p_ST: float) -> int:
        """
        Determine supply tank vent state with hysteresis (for single-tank scenarios).
        
        Parameters
        ----------
        p_ST : float
            Supply tank pressure [Pa]
            
        Returns
        -------
        state : int
            Vent state (0=closed, 1=open)
        """
        P = self.params
        
        if p_ST > P["p_ST_vent_open"]:
            # Pressure above open threshold - open vent
            return 1
        elif p_ST < P["p_ST_vent_close"]:
            # Pressure below close threshold - close vent
            return 0
        else:
            # Within hysteresis band - maintain current state
            return self.ST_vent_state

    def _get_ET_vent_state(self, p_ET: float) -> int:
        """Determine end tank vent state with hysteresis."""
        P = self.params

        if p_ET < P["p_ET_low"]:
            return 0
        elif p_ET > P["p_ET_high"]:
            return 1
        else:
            return self.ET_vent_state

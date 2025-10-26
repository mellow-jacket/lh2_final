"""
LH2 Simulation Module

This module provides the core simulation engine for LH2 transfer operations.
It implements mass and energy balance equations, event detection, and ODE integration.

Key components:
- SimulationState: State vector management
- Simulator: Main simulation orchestrator
- mass_balance: Mass balance equations
- energy_balance: Energy balance equations (stub for future)

This is a simplified initial implementation focusing on mass balance.
Energy balance and detailed heat transfer will be added in future iterations.
"""

from dataclasses import dataclass, field
from typing import Optional, Tuple, Callable
import numpy as np
from scipy.integrate import solve_ivp

from lh2sim.parameters import ScenarioConfig
from lh2sim.properties import FluidProperties
from lh2sim.geometry import cyl_v_to_h, cylinder_cross_section_area
from lh2sim.flow import gas_flow
from lh2sim.control import PressureDrivenControl, PumpDrivenControl, ControlOutputs


@dataclass
class SimulationState:
    """
    State vector for LH2 transfer simulation.
    
    For a simplified two-tank system:
    - Each tank has liquid mass and vapor mass
    - Each tank has liquid internal energy and vapor internal energy
    
    This is a simplified state for initial implementation.
    Full implementation would include:
    - Wall temperatures
    - Multi-zone discretization
    - Interface temperatures
    
    Attributes:
        m_L_supply: Liquid mass in supply tank [kg]
        m_v_supply: Vapor mass in supply tank [kg]
        m_L_receiver: Liquid mass in receiver tank [kg]
        m_v_receiver: Vapor mass in receiver tank [kg]
        U_L_supply: Liquid internal energy in supply tank [J]
        U_v_supply: Vapor internal energy in supply tank [J]
        U_L_receiver: Liquid internal energy in receiver tank [J]
        U_v_receiver: Vapor internal energy in receiver tank [J]
    """
    m_L_supply: float
    m_v_supply: float
    m_L_receiver: float
    m_v_receiver: float
    U_L_supply: float
    U_v_supply: float
    U_L_receiver: float
    U_v_receiver: float
    
    def to_array(self) -> np.ndarray:
        """Convert state to numpy array for ODE solver."""
        return np.array([
            self.m_L_supply,
            self.m_v_supply,
            self.m_L_receiver,
            self.m_v_receiver,
            self.U_L_supply,
            self.U_v_supply,
            self.U_L_receiver,
            self.U_v_receiver,
        ])
    
    @classmethod
    def from_array(cls, arr: np.ndarray) -> 'SimulationState':
        """Create state from numpy array."""
        return cls(
            m_L_supply=arr[0],
            m_v_supply=arr[1],
            m_L_receiver=arr[2],
            m_v_receiver=arr[3],
            U_L_supply=arr[4],
            U_v_supply=arr[5],
            U_L_receiver=arr[6],
            U_v_receiver=arr[7],
        )


@dataclass
class SimulationResult:
    """
    Results from a simulation run.
    
    Attributes:
        t: Time array [s]
        states: Array of state vectors over time
        success: Whether simulation completed successfully
        message: Status message from solver
        nfev: Number of function evaluations
        njev: Number of Jacobian evaluations
    """
    t: np.ndarray
    states: np.ndarray
    success: bool
    message: str
    nfev: int
    njev: int
    
    def get_state_at(self, index: int) -> SimulationState:
        """Get simulation state at a specific time index."""
        return SimulationState.from_array(self.states[:, index])


class Simulator:
    """
    Main simulator for LH2 transfer operations.
    
    This orchestrates the simulation by:
    1. Setting up initial conditions from scenario
    2. Defining the ODE system (mass and energy balances)
    3. Running the integration with event detection
    4. Returning results
    
    Attributes:
        config: Scenario configuration
        properties: Fluid properties calculator
        controller: Control logic (pressure-driven or pump-driven)
    """
    
    def __init__(self, config: ScenarioConfig):
        """
        Initialize simulator with scenario configuration.
        
        Args:
            config: Complete scenario configuration
        """
        self.config = config
        
        # Initialize property calculator
        backend = config.property_backend
        if backend == "CoolProp":
            self.properties = FluidProperties("ParaHydrogen", backend="coolprop")
        else:
            self.properties = FluidProperties("ParaHydrogen", backend="polynomial")
        
        # Initialize controller based on transfer mode
        if config.transfer.mode == "pressure_driven":
            # Create parameter dict for PressureDrivenControl
            control_params = {
                'H': config.receiver_tank.length_or_height,
                'p_ST_slow': config.supply_tank.max_working_pressure * 0.5,
                'p_ST_fast': config.supply_tank.max_working_pressure * 0.7,
                'p_ST_final': config.supply_tank.max_working_pressure * 0.9,
                'p_ET_low': config.receiver_tank.max_working_pressure * 0.9,
                'p_ET_high': config.receiver_tank.max_working_pressure * 1.05,
                'p_ET_final': config.receiver_tank.max_working_pressure * 1.0,
            }
            self.controller = PressureDrivenControl(control_params)
        else:  # pump_driven
            # Create parameter dict for PumpDrivenControl
            control_params = {
                'H': config.receiver_tank.length_or_height,
                'PumpMassTransferSlow': config.transfer.pump_flow_rate * config.physics.rho_liquid * 0.5,
                'PumpMassTransferFast': config.transfer.pump_flow_rate * config.physics.rho_liquid * 1.0,
                'p_ET_low': config.receiver_tank.max_working_pressure * 0.9,
                'p_ET_high': config.receiver_tank.max_working_pressure * 1.05,
                'p_ET_final': config.receiver_tank.max_working_pressure * 1.0,
            }
            self.controller = PumpDrivenControl(control_params)
        
        # Store cross-sectional areas for geometry calculations
        self.A_supply = cylinder_cross_section_area(config.supply_tank.radius)
        self.A_receiver = cylinder_cross_section_area(config.receiver_tank.radius)
    
    def _compute_initial_state(self) -> SimulationState:
        """
        Compute initial state from scenario configuration.
        
        Returns:
            Initial simulation state
        """
        # Supply tank initial conditions
        V_supply = self.config.supply_tank.volume
        fill_supply = self.config.supply_tank.initial_fill_fraction
        T_L_supply = self.config.supply_tank.initial_liquid_temp
        T_v_supply = self.config.supply_tank.initial_vapor_temp
        p_supply = self.config.supply_tank.initial_pressure
        
        # Get densities at specified pressure
        rho_L_supply = self.properties.density(T_L_supply, p_supply)
        # For vapor, use ideal gas law to be consistent with derivative calculation
        rho_v_supply = p_supply / (self.config.physics.R_vapor * T_v_supply)
        
        # Compute masses based on volumes
        V_L_supply = V_supply * fill_supply
        V_v_supply = V_supply * (1 - fill_supply)
        m_L_supply = rho_L_supply * V_L_supply
        m_v_supply = rho_v_supply * V_v_supply
        
        # Compute internal energies (simplified - use specific energy approximation)
        # For now, use c_liquid * T as a simple approximation for liquid
        # and c_v * T for vapor (ignoring pressure work terms)
        u_L_supply = self.config.physics.c_liquid * T_L_supply
        u_v_supply = self.config.physics.c_v_vapor * T_v_supply
        U_L_supply = m_L_supply * u_L_supply
        U_v_supply = m_v_supply * u_v_supply
        
        # Receiver tank initial conditions
        V_receiver = self.config.receiver_tank.volume
        fill_receiver = self.config.receiver_tank.initial_fill_fraction
        T_L_receiver = self.config.receiver_tank.initial_liquid_temp
        T_v_receiver = self.config.receiver_tank.initial_vapor_temp
        p_receiver = self.config.receiver_tank.initial_pressure
        
        # Get densities at specified pressure
        rho_L_receiver = self.properties.density(T_L_receiver, p_receiver)
        # For vapor, use ideal gas law to be consistent
        rho_v_receiver = p_receiver / (self.config.physics.R_vapor * T_v_receiver)
        
        # Compute masses based on volumes
        V_L_receiver = V_receiver * fill_receiver
        V_v_receiver = V_receiver * (1 - fill_receiver)
        m_L_receiver = rho_L_receiver * V_L_receiver
        m_v_receiver = rho_v_receiver * V_v_receiver
        
        # Compute internal energies (simplified)
        u_L_receiver = self.config.physics.c_liquid * T_L_receiver
        u_v_receiver = self.config.physics.c_v_vapor * T_v_receiver
        U_L_receiver = m_L_receiver * u_L_receiver
        U_v_receiver = m_v_receiver * u_v_receiver
        
        return SimulationState(
            m_L_supply=m_L_supply,
            m_v_supply=m_v_supply,
            m_L_receiver=m_L_receiver,
            m_v_receiver=m_v_receiver,
            U_L_supply=U_L_supply,
            U_v_supply=U_v_supply,
            U_L_receiver=U_L_receiver,
            U_v_receiver=U_v_receiver,
        )
    
    def _derivatives(self, t: float, y: np.ndarray) -> np.ndarray:
        """
        Compute time derivatives for ODE system.
        
        This is a simplified implementation focusing on mass transfer.
        Full implementation would include:
        - Detailed heat transfer
        - Multi-zone discretization
        - Wall thermal dynamics
        - Interface dynamics
        
        Args:
            t: Current time [s]
            y: Current state vector
            
        Returns:
            Time derivatives dy/dt
        """
        # Parse state
        state = SimulationState.from_array(y)
        
        # Compute current tank properties
        # Supply tank
        V_supply = self.config.supply_tank.volume
        V_L_supply = state.m_L_supply / self.config.physics.rho_liquid if state.m_L_supply > 0 else 0
        V_v_supply = V_supply - V_L_supply
        
        # Approximate temperatures from internal energies (simplified)
        T_L_supply = self.config.supply_tank.initial_liquid_temp
        T_v_supply = self.config.supply_tank.initial_vapor_temp
        
        # Approximate pressures (simplified - use ideal gas for vapor)
        if V_v_supply > 1e-6:
            p_supply = state.m_v_supply * self.config.physics.R_vapor * T_v_supply / V_v_supply
        else:
            p_supply = self.config.supply_tank.max_working_pressure
        
        # Receiver tank
        V_receiver = self.config.receiver_tank.volume
        V_L_receiver = state.m_L_receiver / self.config.physics.rho_liquid if state.m_L_receiver > 0 else 0
        V_v_receiver = V_receiver - V_L_receiver
        
        T_L_receiver = self.config.receiver_tank.initial_liquid_temp
        T_v_receiver = self.config.receiver_tank.initial_vapor_temp
        
        if V_v_receiver > 1e-6:
            p_receiver = state.m_v_receiver * self.config.physics.R_vapor * T_v_receiver / V_v_receiver
        else:
            p_receiver = self.config.receiver_tank.max_working_pressure
        
        # Compute receiver liquid height for control
        if self.config.receiver_tank.geometry == "vertical_cylinder":
            h_L_receiver = V_L_receiver / self.A_receiver
        else:  # horizontal_cylinder
            h_L_receiver = cyl_v_to_h(
                V=V_L_receiver,
                R=self.config.receiver_tank.radius,
                L=self.config.receiver_tank.length_or_height
            )
        
        # Get control outputs (different signatures for pressure vs pump driven)
        ET_fill_complete = h_L_receiver >= (0.99 * self.config.receiver_tank.length_or_height)
        
        if self.config.transfer.mode == "pressure_driven":
            controls = self.controller.compute_control(
                h_L2=h_L_receiver,
                p_ST=p_supply,
                p_ET=p_receiver,
                ET_fill_complete=ET_fill_complete,
            )
        else:  # pump_driven
            controls = self.controller.compute_control(
                h_L2=h_L_receiver,
                p_ET=p_receiver,
                ET_fill_complete=ET_fill_complete,
            )
        
        # Compute transfer flow rate
        if self.config.transfer.mode == "pressure_driven":
            # Pressure-driven flow
            CA = self.config.transfer.transfer_valve_area * controls.lambda_E
            rho_L = self.config.physics.rho_liquid
            mdot_transfer = gas_flow(
                CA=CA,
                rho=rho_L,
                P1=p_supply,
                P2=p_receiver,
                gamma=self.config.physics.gamma_vapor
            )
        else:  # pump_driven
            # Pump-driven flow (volumetric flow rate * density)
            Q_pump = self.config.transfer.pump_flow_rate
            # For pump-driven, controls.lambda_E represents pump fraction
            mdot_transfer = Q_pump * self.config.physics.rho_liquid * controls.lambda_E
        
        # Mass balance derivatives (simplified - no evaporation/condensation yet)
        dm_L_supply_dt = -mdot_transfer  # Loses liquid
        dm_v_supply_dt = 0.0  # Simplified - no phase change
        dm_L_receiver_dt = mdot_transfer  # Gains liquid
        dm_v_receiver_dt = 0.0  # Simplified - no phase change
        
        # Energy balance derivatives (simplified - using specific energy)
        u_transfer = self.config.physics.c_liquid * T_L_supply
        dU_L_supply_dt = -mdot_transfer * u_transfer
        dU_v_supply_dt = 0.0
        dU_L_receiver_dt = mdot_transfer * u_transfer
        dU_v_receiver_dt = 0.0
        
        # Return derivatives
        return np.array([
            dm_L_supply_dt,
            dm_v_supply_dt,
            dm_L_receiver_dt,
            dm_v_receiver_dt,
            dU_L_supply_dt,
            dU_v_supply_dt,
            dU_L_receiver_dt,
            dU_v_receiver_dt,
        ])
    
    def run(
        self,
        rtol: float = 1e-6,
        atol: float = 1e-8,
        max_step: Optional[float] = None,
    ) -> SimulationResult:
        """
        Run the simulation.
        
        Args:
            rtol: Relative tolerance for ODE solver
            atol: Absolute tolerance for ODE solver
            max_step: Maximum step size [s] (None for automatic)
            
        Returns:
            SimulationResult with time series and final state
        """
        # Get initial conditions
        state0 = self._compute_initial_state()
        y0 = state0.to_array()
        
        # Time span
        t_span = (0.0, self.config.t_final)
        
        # Run ODE integration
        kwargs = {
            'fun': self._derivatives,
            't_span': t_span,
            'y0': y0,
            'method': 'BDF',  # Good for stiff problems
            'rtol': rtol,
            'atol': atol,
            'dense_output': False,
        }
        
        # Only add max_step if it's not None
        if max_step is not None:
            kwargs['max_step'] = max_step
        
        sol = solve_ivp(**kwargs)
        
        return SimulationResult(
            t=sol.t,
            states=sol.y,
            success=sol.success,
            message=sol.message,
            nfev=sol.nfev,
            njev=sol.njev if sol.njev is not None else 0,
        )


# Export key classes
__all__ = [
    "SimulationState",
    "SimulationResult",
    "Simulator",
]

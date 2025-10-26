"""
LH2 Simulation Parameters Module

This module provides parameter configuration classes for LH2 transfer simulation scenarios.
Follows the design from project_instructions.md with validation and clear structure.

Key components:
- TankParameters: Tank geometry and initial conditions
- PhysicsParameters: Thermophysical constants
- TransferParameters: Transfer control parameters
- ScenarioConfig: Complete scenario configuration

References MATLAB implementations for parameter structure but uses clean Python dataclasses.
"""

from dataclasses import dataclass, field
from typing import Literal, Optional
import numpy as np


@dataclass
class TankParameters:
    """
    Parameters for a single tank (supply or receiving).

    Attributes:
        name: Tank identifier (e.g., "ST", "ET", "Trailer", "Dewar")
        geometry: Tank geometry type
        volume: Total tank volume [m³]
        radius: Tank radius [m]
        length_or_height: Cylinder length (horizontal) or height (vertical) [m]
        initial_pressure: Initial pressure [Pa]
        initial_liquid_temp: Initial liquid temperature [K]
        initial_vapor_temp: Initial vapor temperature [K]
        initial_fill_fraction: Initial fill fraction [0-1]
        max_working_pressure: Maximum working pressure [Pa]
        vent_area: Vent valve effective area [m²]
        wall_mass: Tank wall mass for thermal mass modeling [kg] (optional)
        heat_leak_liquid: Heat leak to liquid [W]
        heat_leak_vapor: Heat leak to vapor [W]
        heat_leak_wall: Heat leak to wall [W] (if wall_mass > 0)
    """

    name: str
    geometry: Literal["horizontal_cylinder", "vertical_cylinder", "sphere"]
    volume: float
    radius: float
    length_or_height: float

    # Initial conditions
    initial_pressure: float
    initial_liquid_temp: float
    initial_vapor_temp: float
    initial_fill_fraction: float

    # Operating limits
    max_working_pressure: float

    # Vent parameters
    vent_area: float

    # Thermal parameters
    wall_mass: float = 0.0  # kg, 0 means no wall thermal mass
    heat_leak_liquid: float = 0.0  # W
    heat_leak_vapor: float = 0.0  # W
    heat_leak_wall: float = 0.0  # W

    def __post_init__(self):
        """Validate parameters after initialization."""
        if self.volume <= 0:
            raise ValueError(f"Tank volume must be positive, got {self.volume}")
        if self.radius <= 0:
            raise ValueError(f"Tank radius must be positive, got {self.radius}")
        if self.length_or_height <= 0:
            raise ValueError(f"Tank length/height must be positive, got {self.length_or_height}")
        if not 0 <= self.initial_fill_fraction <= 1:
            raise ValueError(f"Initial fill fraction must be in [0,1], got {self.initial_fill_fraction}")
        if self.initial_pressure <= 0:
            raise ValueError(f"Initial pressure must be positive, got {self.initial_pressure}")
        if self.max_working_pressure <= 0:
            raise ValueError(f"Max working pressure must be positive, got {self.max_working_pressure}")
        if self.vent_area < 0:
            raise ValueError(f"Vent area must be non-negative, got {self.vent_area}")


@dataclass
class PhysicsParameters:
    """
    Thermophysical constants and correlations for hydrogen.

    These are typically constant for a given simulation and follow REFPROP/CoolProp values.
    For ParaHydrogen at typical LH2 conditions.

    Attributes:
        T_critical: Critical temperature [K]
        p_critical: Critical pressure [Pa]
        rho_liquid: Liquid density at reference conditions [kg/m³]
        c_liquid: Liquid specific heat [J/kg/K]
        kappa_liquid: Liquid thermal conductivity [W/m/K]
        mu_liquid: Liquid dynamic viscosity [Pa·s]
        R_vapor: Vapor gas constant [J/kg/K]
        c_v_vapor: Vapor specific heat at constant volume [J/kg/K]
        gamma_vapor: Ratio of specific heats for vapor [-]
        mu_vapor: Vapor dynamic viscosity [Pa·s]
        kappa_vapor: Vapor thermal conductivity [W/m/K]
        g: Gravitational acceleration [m/s²]
        lambda_: Exponent for surface temperature correlation [-]
        n_liquid_nodes: Number of liquid boundary layer nodes (for multi-node energy balance)
        n_vapor_nodes: Number of vapor boundary layer nodes (for multi-node energy balance)
        tmin_liquid: Minimum time scale for liquid boundary layer [s]
        tmin_vapor: Minimum time scale for vapor boundary layer [s]
        wall_specific_heat: Wall material specific heat [J/kg/K]
        wall_initial_temperature: Initial wall temperature [K]
        beta_liquid: Liquid thermal expansion coefficient [1/K]
        beta_vapor: Vapor thermal expansion coefficient [1/K]
    """

    # Critical properties
    T_critical: float = 32.938  # K
    p_critical: float = 1.2884e6  # Pa (186.49 psi)

    # Liquid properties (reference at ~1 bar, 20K)
    rho_liquid: float = 70.9  # kg/m³
    c_liquid: float = 9702.5  # J/kg/K
    kappa_liquid: float = 0.10061  # W/m/K
    mu_liquid: float = 13.54e-6  # Pa·s

    # Vapor properties
    R_vapor: float = 4124.0  # J/kg/K
    c_v_vapor: float = 6490.0  # J/kg/K
    gamma_vapor: float = 5.0 / 3.0  # dimensionless
    mu_vapor: float = 0.98e-6  # Pa·s
    kappa_vapor: float = 0.0166  # W/m/K

    # Other physics constants
    g: float = 9.81  # m/s²
    lambda_: float = 1.5  # Exponent for surface temperature correlation (Ts = Tc * (P/Pc)^(1/lambda))

    # Multi-node energy balance parameters
    n_liquid_nodes: int = 3  # Number of liquid boundary layer nodes (start with 3 for speed)
    n_vapor_nodes: int = 3  # Number of vapor boundary layer nodes (start with 3 for speed)
    tmin_liquid: float = 10.0  # s, minimum time scale for liquid boundary layer
    tmin_vapor: float = 10.0  # s, minimum time scale for vapor boundary layer

    # Wall thermal properties (for receiver tank wall thermal mass)
    wall_specific_heat: float = 500.0  # J/kg/K, typical for stainless steel
    wall_initial_temperature: float = 300.0  # K, ambient temperature
    
    # Thermal expansion coefficients (for natural convection calculations)
    beta_liquid: float = 0.025  # 1/K, typical for LH2
    beta_vapor: float = 1.0 / 20.0  # 1/K, ~1/T for ideal gas at 20K

    @property
    def c_p_vapor(self) -> float:
        """Vapor specific heat at constant pressure [J/kg/K]."""
        return self.c_v_vapor + self.R_vapor
    
    @property
    def Pr_liquid(self) -> float:
        """Liquid Prandtl number [-]."""
        return self.mu_liquid * self.c_liquid / self.kappa_liquid
    
    @property
    def Pr_vapor(self) -> float:
        """Vapor Prandtl number [-]."""
        return self.mu_vapor * self.c_p_vapor / self.kappa_vapor

    def __post_init__(self):
        """Validate physical parameters."""
        if self.T_critical <= 0:
            raise ValueError("Critical temperature must be positive")
        if self.p_critical <= 0:
            raise ValueError("Critical pressure must be positive")
        if self.gamma_vapor <= 1:
            raise ValueError("Gamma must be > 1 for ideal gas")
        if self.n_liquid_nodes < 1:
            raise ValueError("Number of liquid nodes must be at least 1")
        if self.n_vapor_nodes < 1:
            raise ValueError("Number of vapor nodes must be at least 1")


@dataclass
class TransferParameters:
    """
    Parameters for transfer control and valves/pipes.

    Attributes:
        mode: Transfer mode (pressure-driven or pump-driven)
        transfer_valve_area: Effective area of main transfer valve [m²]
        pipe_length: Transfer line length [m]
        pipe_diameter: Transfer line diameter [m]
        pipe_roughness: Pipe wall roughness [m]
        vaporizer_area: Vaporizer valve effective area [m²] (for pressure-driven)
        vaporizer_coefficient: Vaporizer flow coefficient [m²·√(Pa·kg/m³)]
        vaporizer_time_constant: Time constant for vaporizer boil-off [s]
        pump_flow_rate: Pump volumetric flow rate [m³/s] (for pump-driven, base rate)
        pump_flow_slow: Pump mass flow rate for slow fill [kg/s] (for pump-driven)
        pump_flow_fast: Pump mass flow rate for fast fill [kg/s] (for pump-driven)
        pump_flow_topping: Pump mass flow rate for topping [kg/s] (for pump-driven)
        pump_efficiency: Pump efficiency [-] (for pump-driven)
        transfer_time_constant: Time constant for transfer flow lag [s]

        # Control thresholds
        slow_fill_threshold: Height threshold for slow fill [m]
        fast_fill_threshold: Height threshold for fast fill [m]
        topping_fill_threshold: Height threshold for topping [m]
        ET_vent_open_threshold: Pressure ratio to open ET vent [-]
        ET_vent_close_threshold: Pressure ratio to close ET vent [-]
        ST_vent_open_threshold: Pressure ratio to open ST vent [-]
        ST_vent_close_threshold: Pressure ratio to close ST vent [-]
    """

    mode: Literal["pressure_driven", "pump_driven"]

    # Transfer system geometry
    transfer_valve_area: float
    pipe_length: float = 10.0  # m
    pipe_diameter: float = 0.05  # m
    pipe_roughness: float = 1.5e-5  # m (smooth steel)

    # Mode-specific parameters
    vaporizer_area: float = 0.0  # m² (only for pressure-driven)
    vaporizer_coefficient: float = 1.0  # Flow coefficient for vaporizer
    vaporizer_time_constant: float = 10.0  # s, time constant for vaporizer boil-off
    transfer_time_constant: float = 1.0  # s, time constant for transfer flow lag
    pump_flow_rate: float = 0.0  # m³/s (only for pump-driven, base rate)
    pump_flow_slow: float = 0.0  # kg/s (pump-driven slow fill rate)
    pump_flow_fast: float = 0.0  # kg/s (pump-driven fast fill rate)
    pump_flow_topping: float = 0.0  # kg/s (pump-driven topping rate)
    pump_efficiency: float = 0.8  # [-] (only for pump-driven)

    # Control thresholds (as fractions of receiver tank height)
    slow_fill_threshold: float = 0.2  # 20% full
    fast_fill_threshold: float = 0.5  # 50% full
    topping_fill_threshold: float = 0.9  # 90% full

    # Vent control (pressure ratios relative to MWP)
    ET_vent_open_threshold: float = 1.05  # Open at 105% of target
    ET_vent_close_threshold: float = 0.95  # Close at 95% of target
    ST_vent_open_threshold: float = 1.05
    ST_vent_close_threshold: float = 0.95

    def __post_init__(self):
        """Validate transfer parameters."""
        if self.transfer_valve_area <= 0:
            raise ValueError("Transfer valve area must be positive")
        if self.pipe_length <= 0:
            raise ValueError("Pipe length must be positive")
        if self.pipe_diameter <= 0:
            raise ValueError("Pipe diameter must be positive")

        if self.mode == "pressure_driven" and self.vaporizer_area <= 0:
            raise ValueError("Pressure-driven mode requires positive vaporizer area")
        if self.mode == "pump_driven":
            if self.pump_flow_rate <= 0 and (self.pump_flow_slow <= 0 or self.pump_flow_fast <= 0):
                raise ValueError("Pump-driven mode requires positive pump flow rate(s)")

        # Validate thresholds
        if not (0 < self.slow_fill_threshold < self.fast_fill_threshold < self.topping_fill_threshold < 1):
            raise ValueError("Fill thresholds must be in order: 0 < slow < fast < topping < 1")


@dataclass
class ScenarioConfig:
    """
    Complete configuration for an LH2 transfer simulation scenario.

    This is the top-level configuration object that contains all parameters
    needed to run a simulation.

    Attributes:
        name: Scenario name/identifier
        description: Human-readable description
        supply_tank: Supply tank (source) parameters
        receiver_tank: Receiver tank (destination) parameters
        physics: Thermophysical constants
        transfer: Transfer control parameters
        t_final: Simulation end time [s]
        property_backend: Property calculation backend to use
        use_multinode_energy_balance: Enable multi-node energy balance (experimental)
    """

    name: str
    description: str
    supply_tank: TankParameters
    receiver_tank: TankParameters
    physics: PhysicsParameters
    transfer: TransferParameters
    t_final: float
    property_backend: Literal["CoolProp", "polynomial"] = "CoolProp"
    use_multinode_energy_balance: bool = False  # Enable experimental multi-node mode

    def __post_init__(self):
        """Validate complete scenario configuration."""
        if self.t_final <= 0:
            raise ValueError("Simulation time must be positive")

        # Ensure tank names are distinct
        if self.supply_tank.name == self.receiver_tank.name:
            raise ValueError("Supply and receiver tanks must have different names")

        # Validate initial pressures are reasonable
        if self.supply_tank.initial_pressure >= self.supply_tank.max_working_pressure:
            raise ValueError("Supply tank initial pressure exceeds max working pressure")
        if self.receiver_tank.initial_pressure >= self.receiver_tank.max_working_pressure:
            raise ValueError("Receiver tank initial pressure exceeds max working pressure")


def create_trailer_to_dewar_scenario() -> ScenarioConfig:
    """
    Create default Trailer-to-Dewar scenario (pressure-driven).

    Based on LLNL reference model: horizontal trailer to vertical dewar.

    Returns:
        ScenarioConfig for trailer-to-dewar transfer
    """
    # Constants
    bar_to_pa = 1e5

    # Supply tank (Trailer - horizontal cylinder)
    supply_tank = TankParameters(
        name="Trailer",
        geometry="horizontal_cylinder",
        volume=18.0,  # m³
        radius=1.0,  # m
        length_or_height=18.0 / (np.pi * 1.0**2),  # L = V/A
        initial_pressure=1.5 * bar_to_pa,  # Closer to saturation pressure
        initial_liquid_temp=21.0,  # K
        initial_vapor_temp=21.5,  # K (slightly warmer)
        initial_fill_fraction=0.9,
        max_working_pressure=10.0 * bar_to_pa,
        vent_area=0.001,  # m²
        heat_leak_liquid=200.0,  # W
        heat_leak_vapor=40.0,  # W
    )

    # Receiver tank (Dewar - vertical cylinder)
    receiver_tank = TankParameters(
        name="Dewar",
        geometry="vertical_cylinder",
        volume=18.0,  # m³
        radius=1.0,  # m
        length_or_height=18.0 / (np.pi * 1.0**2),  # H = V/A
        initial_pressure=1.2 * bar_to_pa,  # Closer to supply pressure
        initial_liquid_temp=20.0,  # K (slightly colder)
        initial_vapor_temp=20.5,  # K
        initial_fill_fraction=0.1,  # Start nearly empty
        max_working_pressure=10.0 * bar_to_pa,
        vent_area=0.001,  # m²
        wall_mass=5000.0,  # kg (with wall thermal mass)
        heat_leak_liquid=0.0,  # W
        heat_leak_vapor=0.0,  # W
        heat_leak_wall=100.0,  # W
    )

    # Physics parameters (default para-hydrogen)
    physics = PhysicsParameters()

    # Transfer parameters (pressure-driven)
    transfer = TransferParameters(
        mode="pressure_driven",
        transfer_valve_area=0.0001,  # m² (reduced for realistic flow)
        vaporizer_area=0.00005,  # m² (reduced proportionally)
        pipe_length=10.0,  # m
        pipe_diameter=0.05,  # m
    )

    return ScenarioConfig(
        name="Trailer-to-Dewar (Pressure-Driven)",
        description="LLNL reference scenario: horizontal trailer to vertical dewar, pressure-driven transfer",
        supply_tank=supply_tank,
        receiver_tank=receiver_tank,
        physics=physics,
        transfer=transfer,
        t_final=3600.0,  # 1 hour
        property_backend="CoolProp",
    )


def create_pump_driven_scenario() -> ScenarioConfig:
    """
    Create pump-driven transfer scenario.

    Similar geometry to pressure-driven but uses pump instead of pressure difference.
    Based on paper model pump-driven variant.

    Returns:
        ScenarioConfig for pump-driven transfer
    """
    # Start with pressure-driven scenario
    config = create_trailer_to_dewar_scenario()

    # Modify for pump-driven
    config.name = "Trailer-to-Dewar (Pump-Driven)"
    config.description = "Pump-driven transfer from horizontal trailer to vertical dewar"

    # Update transfer parameters with regime-specific flow rates
    # Typical LH2 transfer: 100-1000 kg/hr depending on regime
    rho_liquid = config.physics.rho_liquid  # ~70.9 kg/m³

    config.transfer = TransferParameters(
        mode="pump_driven",
        transfer_valve_area=0.0001,  # m² (not used in pump mode but required)
        pump_flow_rate=0.001,  # m³/s base rate
        pump_flow_slow=rho_liquid * 0.0005,  # kg/s (~35 kg/s = ~127 kg/hr)
        pump_flow_fast=rho_liquid * 0.002,  # kg/s (~142 kg/s = ~511 kg/hr)
        pump_flow_topping=rho_liquid * 0.0003,  # kg/s (~21 kg/s = ~76 kg/hr)
        pump_efficiency=0.8,
        pipe_length=10.0,  # m
        pipe_diameter=0.05,  # m
    )

    return config


def create_single_tank_venting_scenario() -> ScenarioConfig:
    """
    Create a simplified single-tank venting scenario with heat input.
    
    A single tank starts 50% full at ~10 bar. External heat input (simulating
    vaporizer) causes pressure to rise. Vent valve opens at 10 bar setpoint
    to regulate pressure, then closes when pressure falls.
    
    This demonstrates:
    - Pressure control via vent cycling at 10 bar setpoint
    - Vent control with hysteresis around setpoint
    - Liquid level decrease over time as vapor is vented
    - Energy balance during heating and venting
    
    Returns:
        ScenarioConfig for single-tank venting (supply tank only)
    """
    # Constants
    bar_to_pa = 1e5
    
    # Main tank (vertical cylinder at 50% fill, near setpoint pressure)
    # Start at 10 bar so we can see vent cycling immediately
    main_tank = TankParameters(
        name="Main Tank",
        geometry="vertical_cylinder",
        volume=10.0,  # m³ (medium-sized tank)
        radius=1.0,  # m
        length_or_height=10.0 / (np.pi * 1.0**2),  # H = V/A ≈ 3.18 m
        initial_pressure=10.5 * bar_to_pa,  # Above setpoint - will vent immediately
        initial_liquid_temp=30.0,  # K (saturated at 10.5 bar)
        initial_vapor_temp=29.5,  # K (saturated equilibrium)
        initial_fill_fraction=0.5,  # 50% full as specified
        max_working_pressure=15.0 * bar_to_pa,  # Safety limit
        vent_area=0.002,  # m² (moderate vent size for good regulation)
        heat_leak_liquid=0.0,  # W (no heat input for this demo)
        heat_leak_vapor=0.0,  # W
    )
    
    # Dummy receiver tank (not used, but required by ScenarioConfig structure)
    # Set it to be very large and nearly empty
    dummy_tank = TankParameters(
        name="Dummy",
        geometry="vertical_cylinder",
        volume=1000.0,  # Very large
        radius=5.0,
        length_or_height=1000.0 / (np.pi * 5.0**2),
        initial_pressure=1.0 * bar_to_pa,  # Atmospheric
        initial_liquid_temp=20.0,
        initial_vapor_temp=20.0,
        initial_fill_fraction=0.01,  # Nearly empty
        max_working_pressure=20.0 * bar_to_pa,
        vent_area=0.0,  # No venting
        heat_leak_liquid=0.0,
        heat_leak_vapor=0.0,
    )
    
    # Physics parameters (default para-hydrogen)
    physics = PhysicsParameters()
    
    # Transfer parameters - use pump mode to disable vaporizer
    # The heat leak will drive pressure rise, and vent will regulate
    rho_liquid = physics.rho_liquid
    transfer = TransferParameters(
        mode="pump_driven",  # No vaporizer - use heat leak instead
        transfer_valve_area=1e-10,  # Essentially closed (no transfer)
        pump_flow_rate=1e-10,  # Minimal (no actual pumping)
        pump_flow_slow=rho_liquid * 1e-10,
        pump_flow_fast=rho_liquid * 1e-10,
        pump_flow_topping=rho_liquid * 1e-10,
        pipe_length=1.0,
        pipe_diameter=0.01,
        # Set vent thresholds for 10 bar setpoint
        ST_vent_open_threshold=10.0 * bar_to_pa,  # Open at 10 bar
        ST_vent_close_threshold=9.9 * bar_to_pa,  # Close at 9.9 bar
        # Dummy tank thresholds (very high - won't vent)
        ET_vent_open_threshold=20.0 * bar_to_pa,
        ET_vent_close_threshold=19.0 * bar_to_pa,
    )
    
    return ScenarioConfig(
        name="Single Tank Venting at 10 Bar",
        description="Heat input causes pressure rise, vent regulates at 10 bar setpoint, demonstrates liquid level decrease",
        supply_tank=main_tank,
        receiver_tank=dummy_tank,
        physics=physics,
        transfer=transfer,
        t_final=1800.0,  # 30 minutes
        property_backend="CoolProp",
    )


# Export key classes
__all__ = [
    "TankParameters",
    "PhysicsParameters",
    "TransferParameters",
    "ScenarioConfig",
    "create_trailer_to_dewar_scenario",
    "create_pump_driven_scenario",
    "create_single_tank_venting_scenario",
]

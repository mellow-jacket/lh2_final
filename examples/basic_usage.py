"""
Basic Usage Examples for LH2 Simulation Package

Demonstrates the core functionality of geometry, flow, properties, and control modules.
"""

import numpy as np
import sys
sys.path.insert(0, '/home/runner/work/lh2_final/lh2_final')

from lh2sim.geometry import cyl_v_to_h, cylinder_cross_section_area
from lh2sim.flow import gas_flow, vent_flow_rate
from lh2sim.properties import FluidProperties, vapor_pressure
from lh2sim.control import PressureDrivenControl, ControlOutputs


def example_geometry():
    """Demonstrate geometry calculations."""
    print("=" * 60)
    print("GEOMETRY MODULE EXAMPLES")
    print("=" * 60)
    
    # Horizontal cylinder (trailer tank)
    R = 1.5  # m, radius
    L = 10.0  # m, length
    V = 20.0  # m³, liquid volume
    
    print(f"\nHorizontal Cylinder Tank:")
    print(f"  Radius: {R} m")
    print(f"  Length: {L} m")
    print(f"  Liquid Volume: {V} m³")
    
    # Calculate liquid height
    height = cyl_v_to_h(V, R, L)
    print(f"  Liquid Height: {height:.3f} m")
    
    # Calculate cross-sectional area
    A = cylinder_cross_section_area(R)
    print(f"  Cross-sectional Area: {A:.3f} m²")
    
    # Fill percentage
    total_volume = A * L
    fill_percent = (V / total_volume) * 100
    print(f"  Fill Percentage: {fill_percent:.1f}%")


def example_flow():
    """Demonstrate flow calculations."""
    print("\n" + "=" * 60)
    print("FLOW MODULE EXAMPLES")
    print("=" * 60)
    
    # Gas flow through orifice
    print("\nChoked/Non-choked Flow:")
    CA = 0.001  # m², effective flow area
    rho = 1.3  # kg/m³, hydrogen vapor density
    P1 = 1.5e5  # Pa, upstream pressure
    P2 = 1.0e5  # Pa, downstream pressure
    gamma = 1.4  # ratio of specific heats for hydrogen
    
    mdot = gas_flow(CA, rho, P1, P2, gamma)
    print(f"  Orifice Area: {CA} m²")
    print(f"  Upstream Pressure: {P1/1e5:.2f} bar")
    print(f"  Downstream Pressure: {P2/1e5:.2f} bar")
    print(f"  Mass Flow Rate: {mdot*3600:.2f} kg/hr")
    
    # Vent flow
    print("\nVent Flow:")
    P_tank = 1.4e5  # Pa
    P_ambient = 1.0e5  # Pa
    A_vent = 0.0005  # m²
    
    mdot_vent = vent_flow_rate(P_tank, P_ambient, rho, A_vent, gamma)
    print(f"  Tank Pressure: {P_tank/1e5:.2f} bar")
    print(f"  Ambient Pressure: {P_ambient/1e5:.2f} bar")
    print(f"  Vent Area: {A_vent*1e4:.2f} cm²")
    print(f"  Vent Flow Rate: {mdot_vent*3600:.2f} kg/hr")


def example_properties():
    """Demonstrate property calculations."""
    print("\n" + "=" * 60)
    print("PROPERTIES MODULE EXAMPLES")
    print("=" * 60)
    
    # Using polynomial backend (no CoolProp required)
    props = FluidProperties("Hydrogen", backend="polynomial")
    
    print("\nLiquid Hydrogen Properties (Polynomial Backend):")
    T = 20.0  # K
    print(f"  Temperature: {T} K")
    
    rho = props.density(T=T)
    print(f"  Liquid Density: {rho:.2f} kg/m³")
    
    P = props.pressure(T=T)
    print(f"  Vapor Pressure: {P/1e5:.3f} bar")
    
    # Vapor pressure calculation
    print("\nVapor Pressure Calculation:")
    T_film = 22.0  # K
    rho_vapor = 1.5  # kg/m³
    
    P_vapor = vapor_pressure(T_film, rho_vapor)
    print(f"  Film Temperature: {T_film} K")
    print(f"  Vapor Density: {rho_vapor} kg/m³")
    print(f"  Vapor Pressure: {P_vapor/1e5:.3f} bar")


def example_control():
    """Demonstrate control logic."""
    print("\n" + "=" * 60)
    print("CONTROL MODULE EXAMPLES")
    print("=" * 60)
    
    # Control parameters
    params = {
        'H': 10.0,  # m, end tank height
        'p_ST_slow': 1.2e5,  # Pa, slow fill pressure
        'p_ST_fast': 1.5e5,  # Pa, fast fill pressure
        'p_ST_final': 1.1e5,  # Pa, final ST pressure
        'p_ET_low': 1.0e5,  # Pa, ET low pressure limit
        'p_ET_high': 1.3e5,  # Pa, ET high pressure limit
        'p_ET_final': 1.2e5,  # Pa, final ET pressure
    }
    
    controller = PressureDrivenControl(params)
    
    print("\nPressure-Driven Control:")
    print(f"  End Tank Height: {params['H']} m")
    
    # Different fill scenarios
    scenarios = [
        ("Slow Fill (3%)", 0.3, 1.2e5, 1.0e5),
        ("Fast Fill (50%)", 5.0, 1.5e5, 1.1e5),
        ("Topping (90%)", 9.0, 1.3e5, 1.15e5),
    ]
    
    for name, h_L2, p_ST, p_ET in scenarios:
        outputs = controller.compute_control(h_L2, p_ST, p_ET)
        print(f"\n  Scenario: {name}")
        print(f"    Liquid Height: {h_L2:.1f} m ({h_L2/params['H']*100:.0f}%)")
        print(f"    ST Pressure: {p_ST/1e5:.2f} bar")
        print(f"    ET Pressure: {p_ET/1e5:.2f} bar")
        print(f"    Transfer Valve: {outputs.lambda_E*100:.0f}% open")
        print(f"    Vaporizer Valve: {outputs.lambda_V*100:.0f}% open")
        print(f"    ST Vent: {'OPEN' if outputs.ST_vent_state else 'CLOSED'}")
        print(f"    ET Vent: {'OPEN' if outputs.ET_vent_state else 'CLOSED'}")


def example_integrated():
    """Demonstrate integrated calculation."""
    print("\n" + "=" * 60)
    print("INTEGRATED EXAMPLE: Tank Venting Scenario")
    print("=" * 60)
    
    # Tank parameters
    R = 1.5  # m
    L = 10.0  # m
    V_liquid = 15.0  # m³
    
    # Current state
    T = 21.0  # K
    P_tank = 1.4e5  # Pa
    
    print(f"\nInitial Conditions:")
    print(f"  Tank: R={R}m, L={L}m")
    print(f"  Liquid Volume: {V_liquid} m³")
    print(f"  Temperature: {T} K")
    print(f"  Pressure: {P_tank/1e5:.2f} bar")
    
    # Calculate liquid height
    h = cyl_v_to_h(V_liquid, R, L)
    A = cylinder_cross_section_area(R)
    fill_percent = (V_liquid / (A * L)) * 100
    
    print(f"\nGeometry:")
    print(f"  Liquid Height: {h:.3f} m")
    print(f"  Fill Level: {fill_percent:.1f}%")
    
    # Get properties
    props = FluidProperties("Hydrogen", backend="polynomial")
    rho_liquid = props.density(T=T)
    
    # Estimate vapor space
    V_vapor = A * L - V_liquid
    m_liquid = rho_liquid * V_liquid
    
    print(f"\nMass Distribution:")
    print(f"  Liquid Density: {rho_liquid:.1f} kg/m³")
    print(f"  Liquid Mass: {m_liquid:.1f} kg")
    print(f"  Vapor Space: {V_vapor:.2f} m³")
    
    # Calculate vent flow if needed
    P_ambient = 1.0e5  # Pa
    A_vent = 0.0005  # m²
    rho_vapor = 1.3  # kg/m³
    gamma = 1.4
    
    if P_tank > P_ambient:
        mdot_vent = vent_flow_rate(P_tank, P_ambient, rho_vapor, A_vent, gamma)
        print(f"\nVenting Required:")
        print(f"  Vent Flow Rate: {mdot_vent*3600:.2f} kg/hr")
        print(f"  Pressure Relief: {(P_tank - P_ambient)/1e3:.1f} kPa")
    else:
        print(f"\nNo Venting Required")


def main():
    """Run all examples."""
    print("\n")
    print("╔" + "=" * 58 + "╗")
    print("║" + " " * 10 + "LH2 SIMULATION PACKAGE EXAMPLES" + " " * 16 + "║")
    print("╚" + "=" * 58 + "╝")
    
    example_geometry()
    example_flow()
    example_properties()
    example_control()
    example_integrated()
    
    print("\n" + "=" * 60)
    print("Examples completed successfully!")
    print("=" * 60 + "\n")


if __name__ == "__main__":
    main()

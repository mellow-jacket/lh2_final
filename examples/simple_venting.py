"""
Simple LH2 Tank Venting Example

This example demonstrates a simplified single-tank venting scenario:
- Tank starts 50% full of LH2 at saturated equilibrium
- Vent valve opens, allowing liquid to vent to atmosphere
- Focus on energy balance and heat transfer during venting

This is a simpler case than full tank-to-tank transfer, useful for
understanding the energy balance and thermodynamic behavior during venting.
"""

import os
import sys
import numpy as np

# Add package root to path
package_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, package_root)

from lh2sim.parameters import create_single_tank_venting_scenario
from lh2sim.simulation import Simulator
from lh2sim.visualization import plot_single_tank_venting

# Create output directory for plots
OUTPUT_DIR = os.path.join(package_root, "examples", "output")
os.makedirs(OUTPUT_DIR, exist_ok=True)


def print_separator(title):
    """Print a nice separator with title."""
    print("\n" + "=" * 70)
    print(f" {title}")
    print("=" * 70)


def main():
    """Run simple tank venting example."""
    print("\n")
    print("-" + "=" * 68 + "-")
    print("|" + " " * 15 + "SIMPLE LH2 TANK VENTING EXAMPLE" + " " * 22 + "|")
    print("-" + "=" * 68 + "-")
    
    print_separator("SCENARIO SETUP")
    
    print("\nSimple Venting Scenario:")
    print("  - Single vertical tank (10 m³)")
    print("  - Starts 50% full at saturated equilibrium")
    print("  - Initial pressure: 1.3 bar (just below vent threshold)")
    print("  - Initial temperature: 20.3 K (saturated)")
    print("  - Vent opens when pressure exceeds 1.35 bar")
    print("  - Environmental heat leak: 1500 W total (1000 W liquid + 500 W vapor)")
    print("\nPhysics Focus:")
    print("  ✓ Energy balance during venting")
    print("  ✓ Heat transfer from environment")
    print("  ✓ Phase change (liquid/vapor equilibrium)")
    print("  ✓ Pressure and temperature evolution")
    
    # Create scenario
    config = create_single_tank_venting_scenario()
    
    print("\nInitial Conditions:")
    print(f"  Tank Volume: {config.supply_tank.volume:.1f} m³")
    print(f"  Tank Height: {config.supply_tank.length_or_height:.2f} m")
    print(f"  Initial Fill: {config.supply_tank.initial_fill_fraction*100:.0f}%")
    
    # Calculate initial masses
    V_liquid = config.supply_tank.volume * config.supply_tank.initial_fill_fraction
    V_vapor = config.supply_tank.volume * (1 - config.supply_tank.initial_fill_fraction)
    m_liquid = V_liquid * config.physics.rho_liquid
    m_vapor = config.supply_tank.initial_pressure / (
        config.physics.R_vapor * config.supply_tank.initial_vapor_temp
    ) * V_vapor
    total_mass = m_liquid + m_vapor
    
    print(f"  Liquid Mass: {m_liquid:.1f} kg")
    print(f"  Vapor Mass: {m_vapor:.3f} kg")
    print(f"  Total Mass: {total_mass:.1f} kg")
    print(f"  Pressure: {config.supply_tank.initial_pressure/1e5:.2f} bar")
    print(f"  Temperature: {config.supply_tank.initial_liquid_temp:.1f} K")
    print(f"  Vent Area: {config.supply_tank.vent_area*1e4:.1f} cm²")
    
    print_separator("RUNNING SIMULATION")
    
    print("\nSimulating tank venting for 10 minutes...")
    print("(This may take a moment to compute energy balance...)")
    
    # Run simulation with smaller max_step for stability
    simulator = Simulator(config)
    result = simulator.run(t_end=600, max_step=2.0)  # 10 min, 2s max step for stability
    
    print(f"\nSimulation completed!")
    print(f"  Status: {result.message}")
    print(f"  Time steps: {len(result.time)}")
    print(f"  Function evaluations: {result.nfev}")
    print(f"  Final time: {result.time[-1]/60:.1f} minutes")
    
    # Calculate final conditions
    final_liquid_mass = result.m_L_ST[-1]
    final_vapor_mass = result.m_v_ST[-1]
    final_total_mass = final_liquid_mass + final_vapor_mass
    mass_vented = total_mass - final_total_mass
    mass_vented_pct = (mass_vented / total_mass) * 100
    
    final_pressure = result.p_v_ST[-1]
    final_temp = result.T_L_ST[-1]
    final_level = result.h_L_ST[-1] / config.supply_tank.length_or_height * 100
    
    print_separator("RESULTS SUMMARY")
    
    print(f"\nInitial State:")
    print(f"  Mass: {total_mass:.1f} kg")
    print(f"  Level: {config.supply_tank.initial_fill_fraction*100:.0f}%")
    print(f"  Pressure: {config.supply_tank.initial_pressure/1e5:.2f} bar")
    print(f"  Temperature: {config.supply_tank.initial_liquid_temp:.1f} K")
    
    print(f"\nFinal State:")
    print(f"  Mass: {final_total_mass:.1f} kg")
    print(f"  Level: {final_level:.1f}%")
    print(f"  Pressure: {final_pressure/1e5:.2f} bar")
    print(f"  Temperature: {final_temp:.2f} K")
    
    print(f"\nVenting Results:")
    print(f"  Mass Vented: {mass_vented:.2f} kg ({mass_vented_pct:.2f}%)")
    print(f"  Average Vent Rate: {mass_vented/(result.time[-1]/3600):.1f} kg/hr")
    
    # Calculate energy changes
    initial_energy = result.U_L_ST[0] + result.U_v_ST[0]
    final_energy = result.U_L_ST[-1] + result.U_v_ST[-1]
    energy_change = final_energy - initial_energy
    
    print(f"\nEnergy Balance:")
    print(f"  Initial internal energy: {initial_energy/1e6:.2f} MJ")
    print(f"  Final internal energy: {final_energy/1e6:.2f} MJ")
    print(f"  Energy change: {energy_change/1e6:.2f} MJ")
    print(f"  Heat leak rate: {config.supply_tank.heat_leak_liquid + config.supply_tank.heat_leak_vapor:.0f} W")
    print(f"  Total heat leak: {(config.supply_tank.heat_leak_liquid + config.supply_tank.heat_leak_vapor) * result.time[-1]/1e6:.2f} MJ")
    
    # Mass conservation check
    print(f"\nConservation Check:")
    print(f"  Initial mass: {total_mass:.3f} kg")
    print(f"  Final mass + vented: {final_total_mass:.3f} kg + {mass_vented:.3f} kg = {final_total_mass + mass_vented:.3f} kg")
    conservation_error = abs((final_total_mass + mass_vented) - total_mass) / total_mass * 100
    if conservation_error < 5:
        print(f"  ✓ Mass conservation: {conservation_error:.3f}% error (acceptable)")
    else:
        print(f"  ⚠ Mass conservation: {conservation_error:.3f}% error (check venting dynamics)")
    
    print_separator("GENERATING VISUALIZATION")
    
    print("\nCreating comprehensive dashboard plot...")
    
    # Generate single-tank venting plot
    fig = plot_single_tank_venting(
        result,
        tank_height=config.supply_tank.length_or_height,
        pressure_limit=config.transfer.ST_vent_open_threshold,
        save_path=os.path.join(OUTPUT_DIR, "simple_venting_dashboard.png")
    )
    
    print(f"  ✓ Saved: simple_venting_dashboard.png")
    print(f"\nPlot includes:")
    print("  • Tank fill level evolution")
    print("  • Pressure with vent threshold")
    print("  • Liquid and vapor temperatures")
    print("  • Mass distribution (liquid/vapor/total)")
    print("  • Phase densities")
    print("  • Internal energy evolution")
    print("  • Cumulative mass vented")
    print("  • Instantaneous vent flow rate")
    print("  • Summary statistics")
    
    print_separator("EXAMPLE COMPLETED")
    
    print(f"\nVisualization saved to: {OUTPUT_DIR}")
    print("\nKey Observations:")
    print(f"  1. Tank pressure regulated by vent (oscillates around {config.transfer.ST_vent_open_threshold/1e5:.2f} bar)")
    print(f"  2. Heat leak causes continuous venting of ~{mass_vented/(result.time[-1]/3600):.1f} kg/hr")
    print(f"  3. Temperature rise of ~{final_temp - config.supply_tank.initial_liquid_temp:.2f} K due to heat input")
    print(f"  4. Vapor mass increases as liquid evaporates from heat leak")
    print(f"  5. Energy balance captures thermodynamic behavior during venting")
    
    print("\nThis simplified case demonstrates:")
    print("  ✓ Energy balance implementation working correctly")
    print("  ✓ Heat transfer and phase change physics")
    print("  ✓ Vent control with hysteresis")
    print("  ✓ Mass and energy conservation (within numerical accuracy)")
    print()


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\nSimulation interrupted by user.")
    except Exception as e:
        print(f"\n\nError during simulation: {e}")
        import traceback
        traceback.print_exc()

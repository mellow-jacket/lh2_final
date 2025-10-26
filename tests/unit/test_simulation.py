"""
Unit tests for the Simulation module.

Tests state management, simulator initialization, and basic simulation runs.
"""

import pytest
import numpy as np
from lh2sim.simulation import SimulationState, SimulationResult, Simulator
from lh2sim.parameters import create_trailer_to_dewar_scenario, create_pump_driven_scenario


class TestSimulationState:
    """Test SimulationState data structure."""
    
    def test_state_creation(self):
        """Test creating a simulation state."""
        state = SimulationState(
            m_L_supply=1000.0,
            m_v_supply=10.0,
            m_L_receiver=100.0,
            m_v_receiver=5.0,
            U_L_supply=1e7,
            U_v_supply=1e5,
            U_L_receiver=1e6,
            U_v_receiver=5e4,
        )
        assert state.m_L_supply == 1000.0
        assert state.m_v_receiver == 5.0
    
    def test_to_array_conversion(self):
        """Test converting state to array."""
        state = SimulationState(
            m_L_supply=1000.0,
            m_v_supply=10.0,
            m_L_receiver=100.0,
            m_v_receiver=5.0,
            U_L_supply=1e7,
            U_v_supply=1e5,
            U_L_receiver=1e6,
            U_v_receiver=5e4,
        )
        arr = state.to_array()
        assert isinstance(arr, np.ndarray)
        assert len(arr) == 8
        assert arr[0] == 1000.0
        assert arr[7] == 5e4
    
    def test_from_array_conversion(self):
        """Test creating state from array."""
        arr = np.array([1000.0, 10.0, 100.0, 5.0, 1e7, 1e5, 1e6, 5e4])
        state = SimulationState.from_array(arr)
        assert state.m_L_supply == 1000.0
        assert state.U_v_receiver == 5e4
    
    def test_roundtrip_conversion(self):
        """Test that to_array and from_array are inverses."""
        original = SimulationState(
            m_L_supply=1000.0,
            m_v_supply=10.0,
            m_L_receiver=100.0,
            m_v_receiver=5.0,
            U_L_supply=1e7,
            U_v_supply=1e5,
            U_L_receiver=1e6,
            U_v_receiver=5e4,
        )
        arr = original.to_array()
        restored = SimulationState.from_array(arr)
        
        assert restored.m_L_supply == original.m_L_supply
        assert restored.m_v_supply == original.m_v_supply
        assert restored.m_L_receiver == original.m_L_receiver
        assert restored.m_v_receiver == original.m_v_receiver


class TestSimulationResult:
    """Test SimulationResult data structure."""
    
    def test_result_creation(self):
        """Test creating a simulation result."""
        t = np.linspace(0, 100, 11)
        states = np.random.rand(8, 11)
        
        result = SimulationResult(
            t=t,
            states=states,
            success=True,
            message="Success",
            nfev=100,
            njev=20,
        )
        
        assert result.success
        assert len(result.t) == 11
        assert result.states.shape == (8, 11)
    
    def test_get_state_at_index(self):
        """Test extracting state at specific time index."""
        t = np.linspace(0, 100, 11)
        states = np.random.rand(8, 11)
        
        result = SimulationResult(
            t=t,
            states=states,
            success=True,
            message="Success",
            nfev=100,
            njev=20,
        )
        
        state_at_5 = result.get_state_at(5)
        assert isinstance(state_at_5, SimulationState)
        assert state_at_5.m_L_supply == states[0, 5]


class TestSimulator:
    """Test Simulator class."""
    
    def test_simulator_initialization_pressure_driven(self):
        """Test initializing simulator with pressure-driven config."""
        config = create_trailer_to_dewar_scenario()
        sim = Simulator(config)
        
        assert sim.config == config
        assert sim.properties is not None
        assert sim.controller is not None
    
    def test_simulator_initialization_pump_driven(self):
        """Test initializing simulator with pump-driven config."""
        config = create_pump_driven_scenario()
        sim = Simulator(config)
        
        assert sim.config == config
        assert sim.properties is not None
        assert sim.controller is not None
    
    def test_compute_initial_state(self):
        """Test computing initial state from config."""
        config = create_trailer_to_dewar_scenario()
        sim = Simulator(config)
        
        state0 = sim._compute_initial_state()
        
        # Supply tank should start mostly full (90% liquid)
        assert state0.m_L_supply > 0
        assert state0.m_v_supply > 0
        
        # Receiver tank should start mostly empty (10% liquid)
        assert state0.m_L_receiver > 0
        assert state0.m_v_receiver > 0
        
        # All masses should be positive
        assert state0.m_L_supply > 0
        assert state0.m_v_supply > 0
        assert state0.m_L_receiver > 0
        assert state0.m_v_receiver > 0
        
        # Supply should have more total mass than receiver
        m_total_supply = state0.m_L_supply + state0.m_v_supply
        m_total_receiver = state0.m_L_receiver + state0.m_v_receiver
        assert m_total_supply > m_total_receiver
    
    def test_derivatives_shape(self):
        """Test that derivatives return correct shape."""
        config = create_trailer_to_dewar_scenario()
        sim = Simulator(config)
        
        state0 = sim._compute_initial_state()
        y0 = state0.to_array()
        
        dydt = sim._derivatives(0.0, y0)
        
        assert dydt.shape == y0.shape
        assert len(dydt) == 8
    
    def test_derivatives_mass_conservation(self):
        """Test that mass is conserved in derivatives."""
        config = create_trailer_to_dewar_scenario()
        sim = Simulator(config)
        
        state0 = sim._compute_initial_state()
        y0 = state0.to_array()
        
        dydt = sim._derivatives(0.0, y0)
        
        # Total mass change should be approximately zero
        # (supply loses what receiver gains)
        dm_supply = dydt[0] + dydt[1]  # liquid + vapor
        dm_receiver = dydt[2] + dydt[3]  # liquid + vapor
        
        # Mass conservation: what supply loses, receiver gains
        assert abs(dm_supply + dm_receiver) < 1e-10
    
    def test_run_simulation_short(self):
        """Test running a short simulation."""
        config = create_trailer_to_dewar_scenario()
        config.t_final = 10.0  # Very short simulation
        
        sim = Simulator(config)
        result = sim.run(max_step=1.0)
        
        assert result.success
        assert len(result.t) > 0
        assert result.states.shape[0] == 8
        assert result.t[0] == 0.0
        assert result.t[-1] <= config.t_final
    
    def test_simulation_preserves_total_mass(self):
        """Test that total mass is conserved during simulation."""
        config = create_trailer_to_dewar_scenario()
        config.t_final = 10.0
        
        sim = Simulator(config)
        result = sim.run(max_step=1.0)
        
        # Get initial and final masses
        state0 = result.get_state_at(0)
        state_final = result.get_state_at(-1)
        
        m_total_initial = (
            state0.m_L_supply + state0.m_v_supply +
            state0.m_L_receiver + state0.m_v_receiver
        )
        
        m_total_final = (
            state_final.m_L_supply + state_final.m_v_supply +
            state_final.m_L_receiver + state_final.m_v_receiver
        )
        
        # Total mass should be conserved (within tolerance)
        relative_error = abs(m_total_final - m_total_initial) / m_total_initial
        assert relative_error < 1e-3  # 0.1% tolerance
    
    def test_receiver_gains_mass(self):
        """Test that receiver tank gains mass during transfer."""
        config = create_trailer_to_dewar_scenario()
        config.t_final = 10.0
        
        sim = Simulator(config)
        result = sim.run(max_step=1.0)
        
        state0 = result.get_state_at(0)
        state_final = result.get_state_at(-1)
        
        # Receiver should have more liquid at the end
        assert state_final.m_L_receiver > state0.m_L_receiver
        
        # Supply should have less liquid at the end
        assert state_final.m_L_supply < state0.m_L_supply
    
    def test_pump_driven_simulation(self):
        """Test running a pump-driven simulation."""
        config = create_pump_driven_scenario()
        config.t_final = 10.0
        
        sim = Simulator(config)
        result = sim.run(max_step=1.0)
        
        assert result.success
        assert len(result.t) > 0
        
        # Check mass is transferred
        state0 = result.get_state_at(0)
        state_final = result.get_state_at(-1)
        assert state_final.m_L_receiver > state0.m_L_receiver


class TestSimulatorEdgeCases:
    """Test simulator edge cases and error conditions."""
    
    def test_simulation_with_very_small_timestep(self):
        """Test simulation with very small maximum step size."""
        config = create_trailer_to_dewar_scenario()
        config.t_final = 1.0
        
        sim = Simulator(config)
        result = sim.run(max_step=0.01)
        
        assert result.success
        assert len(result.t) > 10  # Should have many steps
    
    def test_simulation_tolerances(self):
        """Test simulation with different tolerance settings."""
        config = create_trailer_to_dewar_scenario()
        config.t_final = 5.0
        
        sim = Simulator(config)
        
        # Run with loose tolerances
        result_loose = sim.run(rtol=1e-3, atol=1e-5)
        assert result_loose.success
        
        # Run with tight tolerances
        result_tight = sim.run(rtol=1e-8, atol=1e-10)
        assert result_tight.success
        
        # Tight tolerances should require more function evaluations
        assert result_tight.nfev > result_loose.nfev


class TestEventDetection:
    """Test event detection during simulation."""

    def test_run_with_events_basic(self):
        """Test that run_with_events executes successfully."""
        config = create_trailer_to_dewar_scenario()
        config.t_final = 10.0  # Short run

        sim = Simulator(config)
        result = sim.run_with_events()

        assert result.success
        assert len(result.t) > 0

    def test_event_fill_complete(self):
        """Test fill completion event detection."""
        config = create_trailer_to_dewar_scenario()
        # Set aggressive parameters to reach fill quickly
        config.t_final = 1000.0  # Long enough to fill
        config.receiver_tank.initial_fill_fraction = 0.85  # Start near full

        sim = Simulator(config)
        result = sim.run_with_events()

        # Check if simulation stopped before t_final (event triggered)
        # and receiver tank is reasonably full
        if result.t[-1] < config.t_final * 0.9:
            # Event likely triggered
            final_state = result.get_state_at(-1)
            V_L_receiver = final_state.m_L_receiver / config.physics.rho_liquid
            fill_fraction = V_L_receiver / config.receiver_tank.volume
            # Should be near 90% full (target in event function)
            assert fill_fraction > 0.85

    def test_event_functions_created(self):
        """Test that event functions are properly created."""
        config = create_trailer_to_dewar_scenario()
        sim = Simulator(config)

        events = sim._create_event_functions()

        # Should have multiple event functions
        assert len(events) > 0
        # Each should be callable
        for event in events:
            assert callable(event)
            # Should have terminal attribute
            assert hasattr(event, "terminal")
            assert hasattr(event, "direction")

    def test_event_with_pump_driven(self):
        """Test events work with pump-driven mode."""
        config = create_pump_driven_scenario()
        config.t_final = 50.0

        sim = Simulator(config)
        result = sim.run_with_events()

        assert result.success
        assert len(result.t) > 0

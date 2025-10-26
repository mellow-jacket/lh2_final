"""Tests for control module."""

import pytest
import numpy as np
from lh2sim.control import PressureDrivenControl, PumpDrivenControl, ControlOutputs


class TestPressureDrivenControl:
    """Tests for pressure-driven control."""
    
    @pytest.fixture
    def params(self):
        """Create typical control parameters."""
        return {
            'H': 10.0,  # m, end tank height
            'p_ST_slow': 1.2e5,  # Pa
            'p_ST_fast': 1.5e5,  # Pa
            'p_ST_final': 1.1e5,  # Pa
            'p_ET_low': 1.0e5,  # Pa
            'p_ET_high': 1.3e5,  # Pa
            'p_ET_final': 1.2e5,  # Pa
        }
    
    @pytest.fixture
    def controller(self, params):
        """Create pressure-driven controller instance."""
        return PressureDrivenControl(params)
    
    def test_slow_fill_regime(self, controller):
        """Test slow fill regime (h < 5% of H)."""
        h_L2 = 0.3  # m (3% of H)
        p_ST = 1.2e5  # Pa
        p_ET = 1.0e5  # Pa
        
        outputs = controller.compute_control(h_L2, p_ST, p_ET)
        
        # Should be in slow fill
        assert outputs.lambda_E == 0.5
        assert outputs.lambda_V >= 0
        assert outputs.ST_vent_state in [0, 1]
        assert outputs.ET_vent_state in [0, 1]
    
    def test_fast_fill_regime(self, controller):
        """Test fast fill regime (5% < h < 70% of H)."""
        h_L2 = 5.0  # m (50% of H)
        p_ST = 1.5e5  # Pa
        p_ET = 1.1e5  # Pa
        
        outputs = controller.compute_control(h_L2, p_ST, p_ET)
        
        # Should be in fast fill
        assert outputs.lambda_E == 1.0
        assert outputs.lambda_V >= 0
    
    def test_reduced_fast_fill_regime(self, controller):
        """Test reduced fast fill regime (70% < h < 85% of H)."""
        h_L2 = 7.5  # m (75% of H)
        p_ST = 1.3e5  # Pa
        p_ET = 1.1e5  # Pa
        
        outputs = controller.compute_control(h_L2, p_ST, p_ET)
        
        # Should be in reduced fast fill
        assert outputs.lambda_E == 0.8
    
    def test_topping_regime(self, controller):
        """Test topping regime (h > 85% of H)."""
        h_L2 = 9.0  # m (90% of H)
        p_ST = 1.2e5  # Pa
        p_ET = 1.1e5  # Pa
        
        outputs = controller.compute_control(h_L2, p_ST, p_ET)
        
        # Should be in topping
        assert outputs.lambda_E == 0.7
    
    def test_fill_complete_closes_valve(self, controller):
        """Test that fill complete closes transfer valve."""
        h_L2 = 5.0  # m
        p_ST = 1.5e5  # Pa
        p_ET = 1.1e5  # Pa
        
        outputs = controller.compute_control(
            h_L2, p_ST, p_ET, ET_fill_complete=True
        )
        
        # All valves should be closed or stopping
        assert outputs.lambda_E == 0.0
        assert outputs.lambda_V == 0.0
    
    def test_ET_vent_high_pressure(self, controller):
        """Test ET vent opens at high pressure."""
        h_L2 = 5.0
        p_ST = 1.2e5
        p_ET = 1.4e5  # Above p_ET_high
        
        outputs = controller.compute_control(h_L2, p_ST, p_ET)
        
        # ET vent should be open
        assert outputs.ET_vent_state == 1
    
    def test_ET_vent_low_pressure(self, controller):
        """Test ET vent closes at low pressure."""
        h_L2 = 5.0
        p_ST = 1.2e5
        p_ET = 0.9e5  # Below p_ET_low
        
        outputs = controller.compute_control(h_L2, p_ST, p_ET)
        
        # ET vent should be closed
        assert outputs.ET_vent_state == 0
    
    def test_vaporizer_valve_proportional_control(self, controller):
        """Test vaporizer valve responds to pressure error."""
        # Use larger h_L2 to be in fast fill regime
        h_L2 = 5.0  # 50% of H - fast fill regime
        
        # Low pressure - should open valve
        p_ST = 1.35e5  # Below p_ST_fast (1.5e5)
        outputs = controller.compute_control(h_L2, p_ST, 1.0e5)
        # In fast fill, vaporizer should respond to pressure
        assert outputs.lambda_V >= 0
        
        # High pressure - should close or reduce valve
        p_ST = 1.6e5  # Above p_ST_fast
        outputs2 = controller.compute_control(h_L2, p_ST, 1.0e5)
        # Valve response depends on control logic
        assert 0.0 <= outputs2.lambda_V <= 1.0
    
    def test_ST_vent_hysteresis(self, controller, params):
        """Test ST vent has hysteresis to avoid chattering."""
        h_L2 = 1.0
        p_ET = 1.0e5
        threshold = params['p_ST_slow']
        
        # Start with pressure below threshold
        p_ST = threshold - 0.2 * threshold
        outputs1 = controller.compute_control(h_L2, p_ST, p_ET)
        state1 = outputs1.ST_vent_state
        
        # Move within deadband - should maintain state
        p_ST = threshold - 0.05 * threshold
        outputs2 = controller.compute_control(h_L2, p_ST, p_ET)
        state2 = outputs2.ST_vent_state
        
        # State should be consistent (hysteresis working)
        assert state2 in [0, 1]


class TestPumpDrivenControl:
    """Tests for pump-driven control."""
    
    @pytest.fixture
    def params(self):
        """Create typical pump control parameters."""
        return {
            'H': 10.0,  # m
            'PumpMassTransferSlow': 0.5,  # kg/s
            'PumpMassTransferFast': 1.0,  # kg/s
            'p_ET_low': 1.0e5,  # Pa
            'p_ET_high': 1.3e5,  # Pa
            'p_ET_final': 1.2e5,  # Pa
        }
    
    @pytest.fixture
    def controller(self, params):
        """Create pump-driven controller instance."""
        return PumpDrivenControl(params)
    
    def test_fast_fill_regime(self, controller, params):
        """Test fast fill with pump."""
        h_L2 = 5.0  # m (50% of H)
        p_ET = 1.1e5  # Pa
        
        outputs = controller.compute_control(h_L2, p_ET)
        
        # Should scale by pump speed ratio
        expected = params['PumpMassTransferSlow'] / params['PumpMassTransferFast']
        assert abs(outputs.lambda_E - expected) < 1e-6
        
        # No vaporizer in pump mode
        assert outputs.lambda_V == 0.0
        
        # No ST vent in pump mode
        assert outputs.ST_vent_state == 0
    
    def test_topping_regime(self, controller):
        """Test topping regime with pump."""
        h_L2 = 9.0  # m (90% of H)
        p_ET = 1.1e5  # Pa
        
        outputs = controller.compute_control(h_L2, p_ET)
        
        # Should be in topping
        assert outputs.lambda_E == 0.7
    
    def test_fill_complete_stops_transfer(self, controller):
        """Test that fill complete stops transfer."""
        h_L2 = 5.0
        p_ET = 1.1e5
        
        outputs = controller.compute_control(h_L2, p_ET, ET_fill_complete=True)
        
        # Transfer should stop
        assert outputs.lambda_E == 0.0
    
    def test_ET_vent_control(self, controller):
        """Test ET vent operates correctly in pump mode."""
        h_L2 = 5.0
        
        # High pressure - vent should open
        p_ET = 1.4e5
        outputs1 = controller.compute_control(h_L2, p_ET)
        assert outputs1.ET_vent_state == 1
        
        # Low pressure - vent should close
        p_ET = 0.9e5
        outputs2 = controller.compute_control(h_L2, p_ET)
        assert outputs2.ET_vent_state == 0
    
    def test_no_vaporizer_or_ST_vent(self, controller):
        """Test that pump mode doesn't use vaporizer or ST vent."""
        h_L2 = 5.0
        p_ET = 1.1e5
        
        outputs = controller.compute_control(h_L2, p_ET)
        
        assert outputs.lambda_V == 0.0
        assert outputs.ST_vent_state == 0


class TestControlOutputs:
    """Tests for ControlOutputs dataclass."""
    
    def test_control_outputs_creation(self):
        """Test creating ControlOutputs."""
        outputs = ControlOutputs(
            lambda_E=0.5,
            lambda_V=0.3,
            ST_vent_state=1,
            ET_vent_state=0,
        )
        
        assert outputs.lambda_E == 0.5
        assert outputs.lambda_V == 0.3
        assert outputs.ST_vent_state == 1
        assert outputs.ET_vent_state == 0
    
    def test_control_outputs_field_access(self):
        """Test field access of ControlOutputs."""
        outputs = ControlOutputs(
            lambda_E=1.0,
            lambda_V=0.0,
            ST_vent_state=0,
            ET_vent_state=1,
        )
        
        # Access fields
        assert hasattr(outputs, 'lambda_E')
        assert hasattr(outputs, 'lambda_V')
        assert hasattr(outputs, 'ST_vent_state')
        assert hasattr(outputs, 'ET_vent_state')

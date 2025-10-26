"""
Unit tests for the Parameters module.

Tests configuration classes, validation, and scenario creation.
"""

import pytest
import numpy as np
from lh2sim.parameters import (
    TankParameters,
    PhysicsParameters,
    TransferParameters,
    ScenarioConfig,
    create_trailer_to_dewar_scenario,
    create_pump_driven_scenario,
    create_single_tank_venting_scenario,
)


class TestTankParameters:
    """Test TankParameters validation and creation."""
    
    def test_valid_tank_creation(self):
        """Test creating a valid tank configuration."""
        tank = TankParameters(
            name="TestTank",
            geometry="vertical_cylinder",
            volume=10.0,
            radius=1.0,
            length_or_height=10.0,
            initial_pressure=3e5,
            initial_liquid_temp=20.0,
            initial_vapor_temp=21.0,
            initial_fill_fraction=0.5,
            max_working_pressure=10e5,
            vent_area=0.001,
        )
        assert tank.name == "TestTank"
        assert tank.volume == 10.0
        assert tank.initial_fill_fraction == 0.5
    
    def test_negative_volume_raises_error(self):
        """Test that negative volume raises ValueError."""
        with pytest.raises(ValueError, match="volume must be positive"):
            TankParameters(
                name="BadTank",
                geometry="vertical_cylinder",
                volume=-10.0,
                radius=1.0,
                length_or_height=10.0,
                initial_pressure=3e5,
                initial_liquid_temp=20.0,
                initial_vapor_temp=21.0,
                initial_fill_fraction=0.5,
                max_working_pressure=10e5,
                vent_area=0.001,
            )
    
    def test_invalid_fill_fraction_raises_error(self):
        """Test that fill fraction outside [0,1] raises ValueError."""
        with pytest.raises(ValueError, match="fill fraction must be in"):
            TankParameters(
                name="BadTank",
                geometry="vertical_cylinder",
                volume=10.0,
                radius=1.0,
                length_or_height=10.0,
                initial_pressure=3e5,
                initial_liquid_temp=20.0,
                initial_vapor_temp=21.0,
                initial_fill_fraction=1.5,  # Invalid
                max_working_pressure=10e5,
                vent_area=0.001,
            )
    
    def test_wall_mass_optional(self):
        """Test that wall mass defaults to zero."""
        tank = TankParameters(
            name="TestTank",
            geometry="vertical_cylinder",
            volume=10.0,
            radius=1.0,
            length_or_height=10.0,
            initial_pressure=3e5,
            initial_liquid_temp=20.0,
            initial_vapor_temp=21.0,
            initial_fill_fraction=0.5,
            max_working_pressure=10e5,
            vent_area=0.001,
        )
        assert tank.wall_mass == 0.0


class TestPhysicsParameters:
    """Test PhysicsParameters properties and validation."""
    
    def test_default_physics_parameters(self):
        """Test default physics parameters are reasonable."""
        params = PhysicsParameters()
        assert params.T_critical > 0
        assert params.p_critical > 0
        assert params.gamma_vapor > 1
    
    def test_c_p_vapor_property(self):
        """Test that c_p = c_v + R."""
        params = PhysicsParameters()
        expected_cp = params.c_v_vapor + params.R_vapor
        assert abs(params.c_p_vapor - expected_cp) < 1e-10
    
    def test_invalid_gamma_raises_error(self):
        """Test that gamma <= 1 raises ValueError."""
        with pytest.raises(ValueError, match="Gamma must be > 1"):
            PhysicsParameters(gamma_vapor=0.9)


class TestTransferParameters:
    """Test TransferParameters validation."""
    
    def test_valid_pressure_driven_config(self):
        """Test creating valid pressure-driven transfer config."""
        params = TransferParameters(
            mode="pressure_driven",
            transfer_valve_area=0.002,
            vaporizer_area=0.0005,
        )
        assert params.mode == "pressure_driven"
        assert params.vaporizer_area > 0
    
    def test_valid_pump_driven_config(self):
        """Test creating valid pump-driven transfer config."""
        params = TransferParameters(
            mode="pump_driven",
            transfer_valve_area=0.002,
            pump_flow_rate=0.005,
        )
        assert params.mode == "pump_driven"
        assert params.pump_flow_rate > 0
    
    def test_pressure_driven_requires_vaporizer_area(self):
        """Test that pressure-driven mode requires vaporizer area."""
        with pytest.raises(ValueError, match="vaporizer area"):
            TransferParameters(
                mode="pressure_driven",
                transfer_valve_area=0.002,
                vaporizer_area=0.0,  # Invalid for pressure-driven
            )
    
    def test_pump_driven_requires_flow_rate(self):
        """Test that pump-driven mode requires pump flow rate."""
        with pytest.raises(ValueError, match="pump flow rate"):
            TransferParameters(
                mode="pump_driven",
                transfer_valve_area=0.002,
                pump_flow_rate=0.0,  # Invalid for pump-driven
            )
    
    def test_fill_thresholds_must_be_ordered(self):
        """Test that fill thresholds must be in increasing order."""
        with pytest.raises(ValueError, match="Fill thresholds must be in order"):
            TransferParameters(
                mode="pressure_driven",
                transfer_valve_area=0.002,
                vaporizer_area=0.0005,
                slow_fill_threshold=0.5,  # Wrong order
                fast_fill_threshold=0.2,
                topping_fill_threshold=0.9,
            )


class TestScenarioConfig:
    """Test complete scenario configuration."""
    
    def test_valid_scenario_creation(self):
        """Test creating a valid complete scenario."""
        supply = TankParameters(
            name="Supply",
            geometry="horizontal_cylinder",
            volume=18.0,
            radius=1.0,
            length_or_height=5.73,
            initial_pressure=3e5,
            initial_liquid_temp=20.0,
            initial_vapor_temp=21.0,
            initial_fill_fraction=0.9,
            max_working_pressure=10e5,
            vent_area=0.001,
        )
        
        receiver = TankParameters(
            name="Receiver",
            geometry="vertical_cylinder",
            volume=18.0,
            radius=1.0,
            length_or_height=5.73,
            initial_pressure=3e5,
            initial_liquid_temp=20.0,
            initial_vapor_temp=21.0,
            initial_fill_fraction=0.1,
            max_working_pressure=10e5,
            vent_area=0.001,
        )
        
        physics = PhysicsParameters()
        
        transfer = TransferParameters(
            mode="pressure_driven",
            transfer_valve_area=0.002,
            vaporizer_area=0.0005,
        )
        
        config = ScenarioConfig(
            name="TestScenario",
            description="Test configuration",
            supply_tank=supply,
            receiver_tank=receiver,
            physics=physics,
            transfer=transfer,
            t_final=1800.0,
        )
        
        assert config.name == "TestScenario"
        assert config.t_final == 1800.0
    
    def test_duplicate_tank_names_raise_error(self):
        """Test that duplicate tank names raise ValueError."""
        supply = TankParameters(
            name="Tank",  # Same name
            geometry="horizontal_cylinder",
            volume=18.0,
            radius=1.0,
            length_or_height=5.73,
            initial_pressure=3e5,
            initial_liquid_temp=20.0,
            initial_vapor_temp=21.0,
            initial_fill_fraction=0.9,
            max_working_pressure=10e5,
            vent_area=0.001,
        )
        
        receiver = TankParameters(
            name="Tank",  # Same name
            geometry="vertical_cylinder",
            volume=18.0,
            radius=1.0,
            length_or_height=5.73,
            initial_pressure=3e5,
            initial_liquid_temp=20.0,
            initial_vapor_temp=21.0,
            initial_fill_fraction=0.1,
            max_working_pressure=10e5,
            vent_area=0.001,
        )
        
        with pytest.raises(ValueError, match="must have different names"):
            ScenarioConfig(
                name="TestScenario",
                description="Test configuration",
                supply_tank=supply,
                receiver_tank=receiver,
                physics=PhysicsParameters(),
                transfer=TransferParameters(
                    mode="pressure_driven",
                    transfer_valve_area=0.002,
                    vaporizer_area=0.0005,
                ),
                t_final=1800.0,
            )


class TestScenarioFactories:
    """Test scenario factory functions."""
    
    def test_create_trailer_to_dewar_scenario(self):
        """Test creating default trailer-to-dewar scenario."""
        config = create_trailer_to_dewar_scenario()
        
        assert config.supply_tank.name == "Trailer"
        assert config.receiver_tank.name == "Dewar"
        assert config.supply_tank.geometry == "horizontal_cylinder"
        assert config.receiver_tank.geometry == "vertical_cylinder"
        assert config.transfer.mode == "pressure_driven"
        assert config.t_final > 0
        
        # Trailer should start mostly full
        assert config.supply_tank.initial_fill_fraction > 0.8
        
        # Dewar should start mostly empty
        assert config.receiver_tank.initial_fill_fraction < 0.2
    
    def test_create_pump_driven_scenario(self):
        """Test creating pump-driven scenario."""
        config = create_pump_driven_scenario()
        
        assert config.transfer.mode == "pump_driven"
        assert config.transfer.pump_flow_rate > 0
        assert config.supply_tank.name == "Trailer"
        assert config.receiver_tank.name == "Dewar"
    
    def test_scenarios_are_valid(self):
        """Test that factory scenarios pass validation."""
        # Should not raise any errors
        config1 = create_trailer_to_dewar_scenario()
        config2 = create_pump_driven_scenario()
        config3 = create_single_tank_venting_scenario()
        
        # Verify they have reasonable values
        assert config1.t_final > 0
        assert config2.t_final > 0
        assert config3.t_final > 0
        assert config1.supply_tank.volume > 0
        assert config2.supply_tank.volume > 0
        assert config3.supply_tank.volume > 0
    
    def test_create_single_tank_venting_scenario(self):
        """Test creating single-tank venting scenario."""
        config = create_single_tank_venting_scenario()
        
        # Should be pressure-driven mode (natural venting behavior)
        assert config.transfer.mode == "pressure_driven"
        
        # Main tank should be 50% full
        assert config.supply_tank.initial_fill_fraction == 0.5
        
        # Main tank should be vertical cylinder
        assert config.supply_tank.geometry == "vertical_cylinder"
        
        # Dummy tank should be very large
        assert config.receiver_tank.volume > config.supply_tank.volume
        
        # Vent threshold should be reasonable for LH2 (around 10 bar for this demo)
        assert 8e5 < config.transfer.ST_vent_open_threshold < 15e5  # 8-15 bar
        
        # Transfer valve should be essentially closed (no transfer)
        assert config.transfer.transfer_valve_area < 0.001
        
        # Should have ST vent parameters set
        assert config.transfer.ST_vent_open_threshold > 0
        assert config.transfer.ST_vent_close_threshold > 0

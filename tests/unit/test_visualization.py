"""
Unit tests for visualization module.
"""

import pytest
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for testing
import matplotlib.pyplot as plt

from lh2sim.visualization import (
    plot_tank_levels,
    plot_pressures,
    plot_temperatures,
    plot_masses,
    plot_densities,
    plot_summary_dashboard
)


@pytest.fixture
def simple_result():
    """Create a simple simulation result for testing."""
    n_points = 10
    time = np.linspace(0, 600, n_points)  # 10 minutes
    
    # Create a mock object with the required attributes
    class MockResult:
        pass
    
    result = MockResult()
    result.time = time
    result.h_L_ET = np.linspace(0.5, 8.0, n_points)
    result.m_L_ST = np.linspace(10000, 8000, n_points)
    result.m_L_ET = np.linspace(100, 2100, n_points)
    result.m_v_ST = np.linspace(50, 45, n_points)
    result.m_v_ET = np.linspace(5, 15, n_points)
    result.p_v_ST = np.linspace(1.2e5, 1.1e5, n_points)
    result.p_v_ET = np.linspace(1.0e5, 1.3e5, n_points)
    result.T_L_ST = np.linspace(20.0, 20.5, n_points)
    result.T_L_ET = np.linspace(20.0, 22.0, n_points)
    result.T_v_ST = np.linspace(20.5, 21.0, n_points)
    result.T_v_ET = np.linspace(21.0, 23.0, n_points)
    result.rho_L_ST = np.linspace(70.0, 69.5, n_points)
    result.rho_L_ET = np.linspace(70.0, 68.0, n_points)
    result.rho_v_ST = np.linspace(1.2, 1.1, n_points)
    result.rho_v_ET = np.linspace(1.0, 1.5, n_points)
    
    return result


class TestPlotTankLevels:
    """Tests for plot_tank_levels function."""
    
    def test_basic_plot(self, simple_result, tmp_path):
        """Test basic tank levels plot."""
        tank_heights = (3.0, 10.0)
        fig = plot_tank_levels(simple_result, tank_heights)
        
        assert fig is not None
        assert len(fig.axes) == 3  # 2x2 subplot with bottom spanning both columns
        
        plt.close(fig)
    
    def test_save_plot(self, simple_result, tmp_path):
        """Test saving plot to file."""
        tank_heights = (3.0, 10.0)
        save_path = tmp_path / "tank_levels.png"
        
        fig = plot_tank_levels(simple_result, tank_heights, save_path=str(save_path))
        
        assert save_path.exists()
        plt.close(fig)
    
    def test_custom_figsize(self, simple_result):
        """Test custom figure size."""
        tank_heights = (3.0, 10.0)
        figsize = (8, 6)
        
        fig = plot_tank_levels(simple_result, tank_heights, figsize=figsize)
        
        assert fig is not None
        plt.close(fig)


class TestPlotPressures:
    """Tests for plot_pressures function."""
    
    def test_basic_plot(self, simple_result):
        """Test basic pressure plot."""
        fig = plot_pressures(simple_result)
        
        assert fig is not None
        assert len(fig.axes) == 1
        
        plt.close(fig)
    
    def test_with_pressure_limits(self, simple_result):
        """Test pressure plot with limits."""
        pressure_limits = {
            'p_ET_high': 1.4e5,
            'p_ET_low': 1.0e5
        }
        
        fig = plot_pressures(simple_result, pressure_limits=pressure_limits)
        
        assert fig is not None
        plt.close(fig)
    
    def test_save_plot(self, simple_result, tmp_path):
        """Test saving pressure plot."""
        save_path = tmp_path / "pressures.png"
        
        fig = plot_pressures(simple_result, save_path=str(save_path))
        
        assert save_path.exists()
        plt.close(fig)


class TestPlotTemperatures:
    """Tests for plot_temperatures function."""
    
    def test_basic_plot(self, simple_result):
        """Test basic temperature plot."""
        fig = plot_temperatures(simple_result)
        
        assert fig is not None
        assert len(fig.axes) == 2
        
        plt.close(fig)
    
    def test_save_plot(self, simple_result, tmp_path):
        """Test saving temperature plot."""
        save_path = tmp_path / "temperatures.png"
        
        fig = plot_temperatures(simple_result, save_path=str(save_path))
        
        assert save_path.exists()
        plt.close(fig)


class TestPlotMasses:
    """Tests for plot_masses function."""
    
    def test_basic_plot(self, simple_result):
        """Test basic mass plot."""
        fig = plot_masses(simple_result)
        
        assert fig is not None
        assert len(fig.axes) == 4  # 2x2 subplot
        
        plt.close(fig)
    
    def test_save_plot(self, simple_result, tmp_path):
        """Test saving mass plot."""
        save_path = tmp_path / "masses.png"
        
        fig = plot_masses(simple_result, save_path=str(save_path))
        
        assert save_path.exists()
        plt.close(fig)


class TestPlotDensities:
    """Tests for plot_densities function."""
    
    def test_basic_plot(self, simple_result):
        """Test basic density plot."""
        fig = plot_densities(simple_result)
        
        assert fig is not None
        assert len(fig.axes) == 2
        
        plt.close(fig)
    
    def test_save_plot(self, simple_result, tmp_path):
        """Test saving density plot."""
        save_path = tmp_path / "densities.png"
        
        fig = plot_densities(simple_result, save_path=str(save_path))
        
        assert save_path.exists()
        plt.close(fig)


class TestPlotSummaryDashboard:
    """Tests for plot_summary_dashboard function."""
    
    def test_basic_dashboard(self, simple_result):
        """Test basic dashboard creation."""
        tank_heights = (3.0, 10.0)
        
        fig = plot_summary_dashboard(simple_result, tank_heights)
        
        assert fig is not None
        assert len(fig.axes) == 9  # 3x3 grid
        
        plt.close(fig)
    
    def test_dashboard_with_limits(self, simple_result):
        """Test dashboard with pressure limits."""
        tank_heights = (3.0, 10.0)
        pressure_limits = {'p_ET_high': 1.4e5}
        
        fig = plot_summary_dashboard(
            simple_result, 
            tank_heights,
            pressure_limits=pressure_limits
        )
        
        assert fig is not None
        plt.close(fig)
    
    def test_save_dashboard(self, simple_result, tmp_path):
        """Test saving dashboard."""
        tank_heights = (3.0, 10.0)
        save_path = tmp_path / "dashboard.png"
        
        fig = plot_summary_dashboard(
            simple_result, 
            tank_heights,
            save_path=str(save_path)
        )
        
        assert save_path.exists()
        plt.close(fig)


class TestPlotEdgeCases:
    """Tests for edge cases in plotting."""
    
    def test_single_point_result(self):
        """Test plotting with single data point."""
        class MockResult:
            pass
        
        result = MockResult()
        result.time = np.array([0.0])
        result.h_L_ET = np.array([5.0])
        result.m_L_ST = np.array([10000.0])
        result.m_L_ET = np.array([100.0])
        result.m_v_ST = np.array([50.0])
        result.m_v_ET = np.array([5.0])
        result.p_v_ST = np.array([1.2e5])
        result.p_v_ET = np.array([1.0e5])
        result.T_L_ST = np.array([20.0])
        result.T_L_ET = np.array([20.0])
        result.T_v_ST = np.array([20.5])
        result.T_v_ET = np.array([21.0])
        result.rho_L_ST = np.array([70.0])
        result.rho_L_ET = np.array([70.0])
        result.rho_v_ST = np.array([1.2])
        result.rho_v_ET = np.array([1.0])
        
        tank_heights = (3.0, 10.0)
        
        # Should not raise errors
        fig1 = plot_tank_levels(result, tank_heights)
        plt.close(fig1)
        
        fig2 = plot_pressures(result)
        plt.close(fig2)
        
        fig3 = plot_summary_dashboard(result, tank_heights)
        plt.close(fig3)

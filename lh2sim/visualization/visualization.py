"""
Visualization Module for LH2 Simulation Results

This module provides plotting utilities for visualizing simulation results,
including time series of pressures, temperatures, mass flows, and tank levels.
"""

import warnings
import numpy as np
import matplotlib.pyplot as plt
from typing import Optional, Tuple
from ..simulation import SimulationResult


def _safe_xlim(ax, time_min):
    """
    Safely set xlim without triggering warnings for single-point or zero-duration data.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to set xlim on
    time_min : np.ndarray
        Time array in minutes
    """
    max_time = time_min[-1] if len(time_min) > 0 else 0
    if max_time > 1e-6:  # Only set xlim if we have meaningful time range
        ax.set_xlim([0, max_time])
    # Otherwise, let matplotlib auto-scale


def plot_tank_levels(
    result: SimulationResult,
    tank_heights: Tuple[float, float],
    figsize: Tuple[float, float] = (12, 8),
    save_path: Optional[str] = None,
) -> plt.Figure:
    """
    Plot liquid levels and transfer flow rate.

    Parameters
    ----------
    result : SimulationResult
        Simulation results to plot
    tank_heights : tuple of float
        (H_ST, H_ET) - heights of supply and end tanks (m)
    figsize : tuple of float, optional
        Figure size in inches (width, height)
    save_path : str, optional
        Path to save the figure. If None, figure is not saved.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The created figure
    """
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)
    fig.suptitle("Tank Levels and Transfer Flow", fontsize=14, fontweight="bold")

    time_min = result.time / 60.0  # Convert to minutes
    H_ST, H_ET = tank_heights

    # Supply Tank level
    ax = fig.add_subplot(gs[0, 0])
    # For horizontal cylinder, use diameter for percentage
    if hasattr(result, "h_L_ST"):
        level_pct_ST = (result.h_L_ST / (2 * H_ST)) * 100
    else:
        # Fallback if height not available
        level_pct_ST = np.zeros_like(time_min)

    ax.plot(time_min, level_pct_ST, "b-", linewidth=2, label="ST level")
    ax.set_ylabel("Supply Tank Level [%]", fontsize=11)
    ax.set_xlabel("Time [min]", fontsize=11)
    _safe_xlim(ax, time_min)
    ax.set_ylim([0, 100])
    ax.legend(loc="best")
    ax.grid(True, alpha=0.3)

    # End Tank level with regime markers
    ax = fig.add_subplot(gs[0, 1])
    level_pct_ET = (result.h_L_ET / H_ET) * 100

    ax.plot(time_min, level_pct_ET, "r-", linewidth=2, label="ET level")
    # Add regime markers
    ax.axhline(y=70, color="gold", linestyle="--", linewidth=1, label="70% (Regime change)")
    ax.axhline(y=80, color="orange", linestyle="--", linewidth=1, label="80% (Regime change)")
    ax.axhline(y=90, color="purple", linestyle="--", linewidth=1, label="90% (Regime change)")

    ax.set_ylabel("End Tank Level [%]", fontsize=11)
    ax.set_xlabel("Time [min]", fontsize=11)
    _safe_xlim(ax, time_min)
    ax.set_ylim([0, 100])
    ax.legend(loc="best", fontsize=8)
    ax.grid(True, alpha=0.3)

    # Transfer flow rate (bottom row, spanning both columns)
    ax = fig.add_subplot(gs[1, :])
    # Calculate transfer flow rate from mass change
    if len(result.m_L_ET) > 1:
        dt = np.diff(result.time)
        dm_ET = np.diff(result.m_L_ET)
        flow_rate = np.zeros_like(time_min)
        flow_rate[1:] = dm_ET / dt * 60  # kg/min
    else:
        flow_rate = np.zeros_like(time_min)

    ax.plot(time_min, flow_rate, "g-", linewidth=2, label="Transfer flow")
    ax.set_ylabel("Transfer Flow Rate [kg/min]", fontsize=11)
    ax.set_xlabel("Time [min]", fontsize=11)
    _safe_xlim(ax, time_min)
    if np.max(flow_rate) > 0:
        ax.set_ylim([0, np.max(flow_rate) * 1.1])
    ax.legend(loc="best")
    ax.grid(True, alpha=0.3)

    # Adjust layout to avoid suptitle overlap
    # Suppress tight_layout warning for GridSpec layouts
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="This figure includes Axes that are not compatible with tight_layout")
        plt.tight_layout(rect=[0, 0, 1, 0.96])

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")

    return fig


def plot_pressures(
    result: SimulationResult,
    pressure_limits: Optional[dict] = None,
    figsize: Tuple[float, float] = (10, 6),
    save_path: Optional[str] = None,
) -> plt.Figure:
    """
    Plot pressure evolution in both tanks.

    Parameters
    ----------
    result : SimulationResult
        Simulation results to plot
    pressure_limits : dict, optional
        Dictionary with 'p_ET_high' and 'p_ET_low' keys for vent limits (Pa)
    figsize : tuple of float, optional
        Figure size in inches (width, height)
    save_path : str, optional
        Path to save the figure

    Returns
    -------
    fig : matplotlib.figure.Figure
        The created figure
    """
    fig, ax = plt.subplots(figsize=figsize)

    time_min = result.time / 60.0

    # Convert to bar
    p_ST_bar = result.p_v_ST / 1e5
    p_ET_bar = result.p_v_ET / 1e5

    ax.plot(time_min, p_ST_bar, "b-", linewidth=2, label="Supply Tank")
    ax.plot(time_min, p_ET_bar, "r-", linewidth=2, label="End Tank")

    # Add pressure limit lines if provided
    if pressure_limits:
        if "p_ET_high" in pressure_limits:
            p_high = pressure_limits["p_ET_high"] / 1e5
            ax.axhline(y=p_high, color="red", linestyle="--", linewidth=1, label=f"ET High Limit ({p_high:.2f} bar)")
        if "p_ET_low" in pressure_limits:
            p_low = pressure_limits["p_ET_low"] / 1e5
            ax.axhline(y=p_low, color="orange", linestyle="--", linewidth=1, label=f"ET Low Limit ({p_low:.2f} bar)")

    ax.set_xlabel("Time [min]", fontsize=11)
    ax.set_ylabel("Pressure [bar]", fontsize=11)
    ax.set_title("Tank Pressures", fontsize=14, fontweight="bold")
    _safe_xlim(ax, time_min)
    ax.legend(loc="best")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")

    return fig


def plot_temperatures(
    result: SimulationResult, figsize: Tuple[float, float] = (10, 8), save_path: Optional[str] = None
) -> plt.Figure:
    """
    Plot temperature evolution in both tanks.

    Parameters
    ----------
    result : SimulationResult
        Simulation results to plot
    figsize : tuple of float, optional
        Figure size in inches (width, height)
    save_path : str, optional
        Path to save the figure

    Returns
    -------
    fig : matplotlib.figure.Figure
        The created figure
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize)
    fig.suptitle("Temperature Evolution", fontsize=14, fontweight="bold")

    time_min = result.time / 60.0

    # Supply Tank temperatures
    ax1.plot(time_min, result.T_L_ST, "b-", linewidth=2, label="Liquid")
    ax1.plot(time_min, result.T_v_ST, "r--", linewidth=2, label="Vapor")

    ax1.set_ylabel("Supply Tank Temperature [K]", fontsize=11)
    ax1.set_xlabel("Time [min]", fontsize=11)
    _safe_xlim(ax1, time_min)
    ax1.legend(loc="best")
    ax1.grid(True, alpha=0.3)

    # End Tank temperatures
    ax2.plot(time_min, result.T_L_ET, "b-", linewidth=2, label="Liquid")
    ax2.plot(time_min, result.T_v_ET, "r--", linewidth=2, label="Vapor")

    ax2.set_ylabel("End Tank Temperature [K]", fontsize=11)
    ax2.set_xlabel("Time [min]", fontsize=11)
    _safe_xlim(ax2, time_min)
    ax2.legend(loc="best")
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")

    return fig


def plot_masses(
    result: SimulationResult, figsize: Tuple[float, float] = (12, 8), save_path: Optional[str] = None
) -> plt.Figure:
    """
    Plot mass evolution in both tanks (liquid and vapor phases).

    Parameters
    ----------
    result : SimulationResult
        Simulation results to plot
    figsize : tuple of float, optional
        Figure size in inches (width, height)
    save_path : str, optional
        Path to save the figure

    Returns
    -------
    fig : matplotlib.figure.Figure
        The created figure
    """
    fig, axes = plt.subplots(2, 2, figsize=figsize)
    fig.suptitle("Mass Distribution", fontsize=14, fontweight="bold")

    time_min = result.time / 60.0

    # Supply Tank liquid mass
    ax = axes[0, 0]
    ax.plot(time_min, result.m_L_ST, "b-", linewidth=2)
    ax.set_ylabel("ST Liquid Mass [kg]", fontsize=11)
    ax.set_xlabel("Time [min]", fontsize=11)
    _safe_xlim(ax, time_min)
    ax.grid(True, alpha=0.3)

    # End Tank liquid mass
    ax = axes[0, 1]
    ax.plot(time_min, result.m_L_ET, "r-", linewidth=2)
    ax.set_ylabel("ET Liquid Mass [kg]", fontsize=11)
    ax.set_xlabel("Time [min]", fontsize=11)
    _safe_xlim(ax, time_min)
    ax.grid(True, alpha=0.3)

    # Supply Tank vapor mass
    ax = axes[1, 0]
    ax.plot(time_min, result.m_v_ST, "b--", linewidth=2)
    ax.set_ylabel("ST Vapor Mass [kg]", fontsize=11)
    ax.set_xlabel("Time [min]", fontsize=11)
    _safe_xlim(ax, time_min)
    ax.grid(True, alpha=0.3)

    # End Tank vapor mass
    ax = axes[1, 1]
    ax.plot(time_min, result.m_v_ET, "r--", linewidth=2)
    ax.set_ylabel("ET Vapor Mass [kg]", fontsize=11)
    ax.set_xlabel("Time [min]", fontsize=11)
    _safe_xlim(ax, time_min)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")

    return fig


def plot_densities(
    result: SimulationResult, figsize: Tuple[float, float] = (10, 8), save_path: Optional[str] = None
) -> plt.Figure:
    """
    Plot density evolution in both tanks.

    Parameters
    ----------
    result : SimulationResult
        Simulation results to plot
    figsize : tuple of float, optional
        Figure size in inches (width, height)
    save_path : str, optional
        Path to save the figure

    Returns
    -------
    fig : matplotlib.figure.Figure
        The created figure
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize)
    fig.suptitle("Density Evolution", fontsize=14, fontweight="bold")

    time_min = result.time / 60.0

    # Vapor densities
    ax1.plot(time_min, result.rho_v_ST, "b-", linewidth=2, label="Supply Tank")
    ax1.plot(time_min, result.rho_v_ET, "r-", linewidth=2, label="End Tank")

    ax1.set_ylabel("Vapor Density [kg/m³]", fontsize=11)
    ax1.set_xlabel("Time [min]", fontsize=11)
    _safe_xlim(ax1, time_min)
    ax1.legend(loc="best")
    ax1.grid(True, alpha=0.3)

    # Liquid densities
    ax2.plot(time_min, result.rho_L_ST, "b-", linewidth=2, label="Supply Tank")
    ax2.plot(time_min, result.rho_L_ET, "r-", linewidth=2, label="End Tank")

    ax2.set_ylabel("Liquid Density [kg/m³]", fontsize=11)
    ax2.set_xlabel("Time [min]", fontsize=11)
    _safe_xlim(ax2, time_min)
    ax2.legend(loc="best")
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")

    return fig


def plot_summary_dashboard(
    result: SimulationResult,
    tank_heights: Tuple[float, float],
    pressure_limits: Optional[dict] = None,
    figsize: Tuple[float, float] = (16, 10),
    save_path: Optional[str] = None,
) -> plt.Figure:
    """
    Create a comprehensive dashboard with all key metrics.

    Parameters
    ----------
    result : SimulationResult
        Simulation results to plot
    tank_heights : tuple of float
        (H_ST, H_ET) - heights of supply and end tanks (m)
    pressure_limits : dict, optional
        Dictionary with pressure limit keys
    figsize : tuple of float, optional
        Figure size in inches (width, height)
    save_path : str, optional
        Path to save the figure

    Returns
    -------
    fig : matplotlib.figure.Figure
        The created figure
    """
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)

    time_min = result.time / 60.0
    H_ST, H_ET = tank_heights

    # Top row: Tank levels
    ax1 = fig.add_subplot(gs[0, 0])
    level_pct_ET = (result.h_L_ET / H_ET) * 100
    ax1.plot(time_min, level_pct_ET, "r-", linewidth=2)
    ax1.axhline(y=70, color="gold", linestyle="--", linewidth=0.8, alpha=0.7)
    ax1.axhline(y=80, color="orange", linestyle="--", linewidth=0.8, alpha=0.7)
    ax1.axhline(y=90, color="purple", linestyle="--", linewidth=0.8, alpha=0.7)
    ax1.set_ylabel("ET Level [%]", fontsize=10)
    ax1.set_xlabel("Time [min]", fontsize=10)
    _safe_xlim(ax1, time_min)
    ax1.set_ylim([0, 100])
    ax1.grid(True, alpha=0.3)
    ax1.set_title("End Tank Level", fontsize=11)

    # Pressures
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.plot(time_min, result.p_v_ST / 1e5, "b-", linewidth=2, label="ST")
    ax2.plot(time_min, result.p_v_ET / 1e5, "r-", linewidth=2, label="ET")
    if pressure_limits:
        if "p_ET_high" in pressure_limits:
            ax2.axhline(y=pressure_limits["p_ET_high"] / 1e5, color="red", linestyle="--", linewidth=0.8, alpha=0.7)
    ax2.set_ylabel("Pressure [bar]", fontsize=10)
    ax2.set_xlabel("Time [min]", fontsize=10)
    _safe_xlim(ax2, time_min)
    ax2.legend(loc="best", fontsize=8)
    ax2.grid(True, alpha=0.3)
    ax2.set_title("Tank Pressures", fontsize=11)

    # Transfer flow rate
    ax3 = fig.add_subplot(gs[0, 2])
    if len(result.m_L_ET) > 1:
        dt = np.diff(result.time)
        dm_ET = np.diff(result.m_L_ET)
        flow_rate = np.zeros_like(time_min)
        flow_rate[1:] = dm_ET / dt * 60  # kg/min
    else:
        flow_rate = np.zeros_like(time_min)
    ax3.plot(time_min, flow_rate, "g-", linewidth=2)
    ax3.set_ylabel("Flow [kg/min]", fontsize=10)
    ax3.set_xlabel("Time [min]", fontsize=10)
    _safe_xlim(ax3, time_min)
    ax3.grid(True, alpha=0.3)
    ax3.set_title("Transfer Flow Rate", fontsize=11)

    # Middle row: Temperatures
    ax4 = fig.add_subplot(gs[1, 0])
    ax4.plot(time_min, result.T_L_ST, "b-", linewidth=2, label="Liquid")
    ax4.plot(time_min, result.T_v_ST, "b--", linewidth=2, label="Vapor")
    ax4.set_ylabel("Temperature [K]", fontsize=10)
    ax4.set_xlabel("Time [min]", fontsize=10)
    _safe_xlim(ax4, time_min)
    ax4.legend(loc="best", fontsize=8)
    ax4.grid(True, alpha=0.3)
    ax4.set_title("Supply Tank Temperatures", fontsize=11)

    ax5 = fig.add_subplot(gs[1, 1])
    ax5.plot(time_min, result.T_L_ET, "r-", linewidth=2, label="Liquid")
    ax5.plot(time_min, result.T_v_ET, "r--", linewidth=2, label="Vapor")
    ax5.set_ylabel("Temperature [K]", fontsize=10)
    ax5.set_xlabel("Time [min]", fontsize=10)
    _safe_xlim(ax5, time_min)
    ax5.legend(loc="best", fontsize=8)
    ax5.grid(True, alpha=0.3)
    ax5.set_title("End Tank Temperatures", fontsize=11)

    # Densities
    ax6 = fig.add_subplot(gs[1, 2])
    ax6.plot(time_min, result.rho_v_ST, "b-", linewidth=1.5, label="ST vapor")
    ax6.plot(time_min, result.rho_v_ET, "r-", linewidth=1.5, label="ET vapor")
    ax6.set_ylabel("Vapor Density [kg/m³]", fontsize=10)
    ax6.set_xlabel("Time [min]", fontsize=10)
    _safe_xlim(ax6, time_min)
    ax6.legend(loc="best", fontsize=8)
    ax6.grid(True, alpha=0.3)
    ax6.set_title("Vapor Densities", fontsize=11)

    # Bottom row: Masses
    ax7 = fig.add_subplot(gs[2, 0])
    ax7.plot(time_min, result.m_L_ST, "b-", linewidth=2, label="Liquid")
    ax7.plot(time_min, result.m_v_ST, "b--", linewidth=2, label="Vapor")
    ax7.set_ylabel("Mass [kg]", fontsize=10)
    ax7.set_xlabel("Time [min]", fontsize=10)
    _safe_xlim(ax7, time_min)
    ax7.legend(loc="best", fontsize=8)
    ax7.grid(True, alpha=0.3)
    ax7.set_title("Supply Tank Masses", fontsize=11)

    ax8 = fig.add_subplot(gs[2, 1])
    ax8.plot(time_min, result.m_L_ET, "r-", linewidth=2, label="Liquid")
    ax8.plot(time_min, result.m_v_ET, "r--", linewidth=2, label="Vapor")
    ax8.set_ylabel("Mass [kg]", fontsize=10)
    ax8.set_xlabel("Time [min]", fontsize=10)
    _safe_xlim(ax8, time_min)
    ax8.legend(loc="best", fontsize=8)
    ax8.grid(True, alpha=0.3)
    ax8.set_title("End Tank Masses", fontsize=11)

    # Total system mass (conservation check)
    ax9 = fig.add_subplot(gs[2, 2])
    total_mass = result.m_L_ST + result.m_v_ST + result.m_L_ET + result.m_v_ET
    ax9.plot(time_min, total_mass, "k-", linewidth=2)
    ax9.set_ylabel("Total Mass [kg]", fontsize=10)
    ax9.set_xlabel("Time [min]", fontsize=10)
    _safe_xlim(ax9, time_min)
    ax9.grid(True, alpha=0.3)
    ax9.set_title("Total System Mass", fontsize=11)

    fig.suptitle("LH2 Transfer Simulation Dashboard", fontsize=16, fontweight="bold", y=0.995)

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")

    return fig


__all__ = [
    "plot_tank_levels",
    "plot_pressures",
    "plot_temperatures",
    "plot_masses",
    "plot_densities",
    "plot_summary_dashboard",
]

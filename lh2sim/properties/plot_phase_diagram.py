"""Plot a hydrogen phase diagram (saturation curve) using the project's property package.

This script attempts to use CoolProp via `lh2sim.properties` when available. If CoolProp
is not installed, it falls back to a Clausius--Clapeyron approximation calibrated at
normal boiling point (20.271 K at 1 atm) so the plot is still useful for visualization.

Output:
  - phase_diagram_h2.png saved into repository root

Usage:
  python plot_phase_diagram.py

Dependencies:
  - numpy, matplotlib
  - optional: CoolProp (recommended for accurate saturation curve). If not installed,
    the script produces an approximate curve and prints a note.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
# Use GridSpec to reserve columns for colorbars (more robust than append_axes)

ROOT = Path(__file__).resolve().parent
OUTFILE = ROOT / "phase_diagram_h2.png"

try:
    # import the project's property abstraction
    from lh2sim.properties import COOLPROP_AVAILABLE

    coolprop_ok = bool(COOLPROP_AVAILABLE)
except Exception:
    coolprop_ok = False

if coolprop_ok:
    try:
        from CoolProp.CoolProp import PropsSI

        def saturation_pressure_coolprop(T: float) -> float:
            """Return saturation pressure [Pa] at temperature T [K] using CoolProp."""
            # Q=0 (saturated liquid) at given T returns saturation pressure
            return PropsSI("P", "T", T, "Q", 0, "Hydrogen")

        # Try to obtain critical point from CoolProp
        try:
            Tcrit = PropsSI("Tcrit", "Hydrogen")
            Pcrit = PropsSI("pcrit", "Hydrogen")
        except Exception:
            # Fallback values (literature): Tc ~ 33.145 K, Pc ~ 1.296e6 Pa
            Tcrit = 33.145
            Pcrit = 1.296e6
    except Exception:
        coolprop_ok = False

if not coolprop_ok:
    print("CoolProp not available through lh2sim.properties — using Clausius-Clapeyron fallback.")

    # Physical constants and calibration point
    M_molar = 2.01588e-3  # kg/mol for H2
    R_universal = 8.314462618  # J/mol/K
    R_specific = R_universal / M_molar  # J/kg/K

    # Approximate latent heat of vaporization for H2 near boiling (~20 K)
    # Source: order-of-magnitude reference; for accuracy use CoolProp/REFPROP
    L_v = 4.49e5  # J/kg (approx)

    # Calibration at 1 atm boiling point for H2
    T_boil = 20.271  # K
    P_boil = 101325.0  # Pa

    # Derived constant for Clausius-Clapeyron integration (integrated form)
    # ln P = -L/(R_specific * T) + C  -> C = ln(P) + L/(R*T)
    C_cc = math.log(P_boil) + L_v / (R_specific * T_boil)

    # Critical point approximate (literature)
    Tcrit = 33.145
    Pcrit = 1.296e6

    def saturation_pressure_cc(T: float) -> float:
        """Clausius-Clapeyron exponential approximation for Psat(T)."""
        return math.exp(-L_v / (R_specific * T) + C_cc)


def build_curve(n_points: int = 400):
    # Build temperature array from 10 K up to just below critical
    T_min = 10.0
    T_max = float(Tcrit) * 0.999
    T = np.linspace(T_min, T_max, n_points)

    P = np.empty_like(T)
    for i, Ti in enumerate(T):
        if coolprop_ok:
            try:
                P[i] = saturation_pressure_coolprop(Ti)
            except Exception:
                P[i] = np.nan
        else:
            P[i] = saturation_pressure_cc(Ti)

    return T, P


def build_density_field(T_min=10.0, T_max=None, P_min=1e2, P_max=None, nT=140, nP=140):
    """Build a T-P grid and compute density [kg/m^3] at each point.

    Uses CoolProp if available (recommended). If not available, falls back to
    an ideal-gas estimate rho = P / (R_specific * T) which is useful only for
    vapor-phase regions.
    Returns (T_grid, P_grid, rho_grid) with P in Pa.
    """
    if T_max is None:
        T_max = float(Tcrit) * 1.02
    if P_max is None:
        P_max = float(Pcrit) * 1.5

    T_vals = np.linspace(T_min, T_max, nT)
    # Use log-spaced pressures from P_min to P_max
    P_vals = np.logspace(np.log10(P_min), np.log10(P_max), nP)

    T_mesh, P_mesh = np.meshgrid(T_vals, P_vals, indexing="xy")

    rho_mesh = np.full_like(T_mesh, np.nan, dtype=float)

    # Try to compute with CoolProp if available. If available, compute a
    # phase-aware density field by comparing local pressure to the saturated
    # pressure at that temperature and asking for the saturated vapor/liquid
    # densities as appropriate. This avoids mixing single-phase and two-phase
    # behavior which produced spurious low-density "gaps" when masked.
    cp_local_ok = False
    try:
        from CoolProp.CoolProp import PropsSI

        cp_local_ok = True
    except Exception:
        cp_local_ok = False

    # Compute a 1D saturated pressure along the mesh temperatures (if possible)
    T_vals = T_vals  # (already defined)
    try:
        if cp_local_ok:
            # Use CoolProp to obtain Psat(T) robustly
            P_sat_1d = np.array([PropsSI("P", "T", float(Ti), "Q", 0, "Hydrogen") for Ti in T_vals])
        else:
            # Use the module-level Clausius-Clapeyron fallback
            P_sat_1d = np.array([saturation_pressure_cc(float(Ti)) for Ti in T_vals])
    except Exception:
        # If anything fails, fill with NaNs so we fall back gracefully later
        P_sat_1d = np.full_like(T_vals, np.nan, dtype=float)

    # Tile to mesh shape for comparison with P_mesh
    P_sat_mesh = np.tile(P_sat_1d[np.newaxis, :], (P_mesh.shape[0], 1))

    if cp_local_ok:
        # Fill each cell with either saturated-vapor or saturated-liquid density
        for i in range(P_mesh.shape[0]):
            for j in range(P_mesh.shape[1]):
                Ti = float(T_mesh[i, j])
                Pi = float(P_mesh[i, j])
                P_sat_local = P_sat_mesh[i, j]
                try:
                    if np.isfinite(P_sat_local) and Pi < P_sat_local:
                        # Vapor-like region: use saturated vapor density at (T)
                        rho_mesh[i, j] = PropsSI("D", "T", Ti, "Q", 1, "Hydrogen")
                    else:
                        # Liquid-like region (including Pi >= Psat): saturated liquid density
                        rho_mesh[i, j] = PropsSI("D", "T", Ti, "Q", 0, "Hydrogen")
                except Exception:
                    rho_mesh[i, j] = np.nan
    else:
        # Fallback: ideal-gas for vapor-like region; leave liquid region as NaN
        M_molar = 2.01588e-3  # kg/mol
        R_universal = 8.314462618
        R_specific = R_universal / M_molar
        # Identify vapor-like cells where P < P_sat (if P_sat available), else assume vapor
        if np.all(np.isnan(P_sat_1d)):
            vapor_mask = np.ones_like(P_mesh, dtype=bool)
        else:
            vapor_mask = P_mesh < P_sat_mesh

        # Compute ideal-gas for vapor cells
        rho_mesh[vapor_mask] = P_mesh[vapor_mask] / (R_specific * T_mesh[vapor_mask])
        # Liquid-like entries remain NaN (so they'll be masked in plotting)

    return T_mesh, P_mesh, rho_mesh


def plot_phase_diagram(save_path: Optional[Path] = None):
    T, P = build_curve()

    # local flag to avoid assigning to the module-level `coolprop_ok`
    cp_ok = coolprop_ok
    # Convert pressure to GPa for x-axis
    P_GPa = P / 1e9

    # Compute a property to map to color: use saturated vapor density (rho_v)
    rho = np.empty_like(P)
    if cp_ok:
        try:
            from CoolProp.CoolProp import PropsSI

            for i, Ti in enumerate(T):
                try:
                    # saturated vapor density (Q=1)
                    rho[i] = PropsSI("D", "T", float(Ti), "Q", 1, "Hydrogen")
                except Exception:
                    rho[i] = np.nan
        except Exception:
            # fallback to ideal-gas estimate
            cp_ok = False

    if not cp_ok:
        # Use ideal-gas rho = P / (R_specific * T) as an approximate vapor density
        M_molar = 2.01588e-3  # kg/mol for H2
        R_universal = 8.314462618  # J/mol/K
        R_specific = R_universal / M_molar
        rho = P / (R_specific * T)

    # Compute density field for background visualization
    T_mesh, P_mesh, rho_mesh = build_density_field(T_min=T.min(), T_max=T.max(), P_min=max(1.0, P.min()), P_max=max(P.max() * 1.2, Pcrit * 1.2), nT=160, nP=160)

    # Convert P mesh to GPa for plotting
    P_mesh_GPa = P_mesh / 1e9

    # Create a figure with a 1x3 GridSpec: main plot + main colorbar + liquid colorbar
    # Reserving dedicated columns for the two colorbars is robust across figure
    # sizes and backends (no fragile manual padding required).
    from matplotlib.colors import LogNorm

    fig = plt.figure(figsize=(10, 6), constrained_layout=True)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 0.04, 0.04])
    ax = fig.add_subplot(gs[0, 0])
    cax_main = fig.add_subplot(gs[0, 1])
    cax_liq = fig.add_subplot(gs[0, 2])

    # Use LogNorm for coloring density across many decades
    # Mask invalid / NaN cells so they render as transparent instead of as the
    # lowest colormap color (which looked like "gaps" previously).
    rho_plot = np.ma.masked_invalid(rho_mesh)
    # Avoid non-positive values for LogNorm by computing vmin/vmax over finite data
    finite_vals = rho_plot.compressed() if rho_plot.count() > 0 else np.array([1e-12])
    vmin = finite_vals.min() if finite_vals.size > 0 else 1e-12
    vmax = finite_vals.max() if finite_vals.size > 0 else vmin * 1e3

    cmap = plt.get_cmap("plasma")
    # Ensure masked values are transparent
    cmap.set_bad((0, 0, 0, 0))

    # Unified normalization for all density-based artists
    norm = LogNorm(vmin=max(vmin, 1e-12), vmax=vmax)

    pcm = ax.pcolormesh(P_mesh_GPa, T_mesh, rho_plot, cmap=cmap, norm=norm, shading="auto")
    # Main colorbar attached to reserved axis
    cbar = fig.colorbar(pcm, cax=cax_main)
    cbar.set_label("Density [kg/m^3] (log scale)")

    # --- Overlay a separate 'hot' scale for liquid-region density ---
    try:
        # Reconstruct P_sat along mesh T values and create a liquid mask
        T_vals_mesh = T_mesh[0, :]
        P_sat_on_T = np.interp(T_vals_mesh, T, P)
        P_sat_mesh = np.tile(P_sat_on_T[np.newaxis, :], (P_mesh.shape[0], 1))
        liquid_mask = (P_mesh >= P_sat_mesh) & np.isfinite(rho_plot)

        # Create a masked array that shows only liquid cells
        rho_liq = np.ma.masked_where(~liquid_mask, rho_plot)

        if rho_liq.count() > 0:
            # Compute a normalization specifically for liquid densities
            liq_vals = rho_liq.compressed()
            vmin_liq = liq_vals.min()
            vmax_liq = liq_vals.max()

            cmap_liq = plt.get_cmap("hot")
            cmap_liq.set_bad((0, 0, 0, 0))
            norm_liq = LogNorm(vmin=max(vmin_liq, 1e-12), vmax=vmax_liq)

            # Draw liquid overlay on top of the background (slightly transparent)
            liq_pcm = ax.pcolormesh(P_mesh_GPa, T_mesh, rho_liq, cmap=cmap_liq, norm=norm_liq, shading="auto", alpha=0.95, zorder=34)

            # Place the liquid colorbar in the reserved right-side axis so it does
            # not overlap the main colorbar. We reserved `cax_liq` in the GridSpec
            # above; attach the colorbar there. If this axis is somehow missing or
            # invalid, fall back to attaching to the main axes.
            try:
                fig.colorbar(liq_pcm, cax=cax_liq)
                cax_liq.yaxis.set_label_text("Liquid density [kg/m^3]")
            except Exception:
                # Fallback: attach to main axes (may overlap)
                fig.colorbar(liq_pcm, ax=ax, pad=0.02, fraction=0.03).set_label("Liquid density [kg/m^3]")
    except Exception:
        # If anything fails here, continue without the liquid overlay
        pass

    # Overlay the saturation curve: draw a thick white line then a thinner black line
    # so it remains visible over the density colormap. We'll also compute an explicit
    # saturation contour from the saturated pressure curve interpolated onto the
    # mesh, and use it to contour-fill the liquid region for clearer phase separation.
    sat_x = P / 1e9
    ax.plot(sat_x, T, color="white", lw=4, zorder=30)
    ax.plot(sat_x, T, color="black", lw=1.25, zorder=31)

    # Build a mesh-aligned saturated-pressure field by interpolating the 1D
    # saturation curve (P vs T) onto the T values used in the mesh. Then compute
    # phi = log10(P_mesh) - log10(P_sat_mesh); the phi==0 contour corresponds to
    # the saturation boundary (P == P_sat(T)). This is robust on a log-pressure
    # mesh and draws a crisp contour line even when the colormap is busy.
    try:
        # T values used for mesh (columns)
        T_vals_mesh = T_mesh[0, :]
        # Interpolate P_sat(T) onto mesh T values. Use np.interp which handles
        # out-of-bounds by clamping to ends; values outside the curve will not
        # produce a valid contour and are masked later by phi.
        P_sat_on_T = np.interp(T_vals_mesh, T, P)
        # Tile the 1D P_sat onto full mesh shape (rows = pressure direction)
        P_sat_mesh = np.tile(P_sat_on_T[np.newaxis, :], (P_mesh.shape[0], 1))

        # Compute phi in log10 space to handle wide pressure ranges
        with np.errstate(divide="ignore", invalid="ignore"):
            phi = np.log10(P_mesh) - np.log10(P_sat_mesh)

        # Contour-fill the liquid region (phi >= 0) with a subtle shading
        # placed above the density colormap but below saturation line markers.
        # Use two levels: below 0 (gas) and >=0 (liquid); select a soft color.
        cf = ax.contourf(P_mesh_GPa, T_mesh, phi, levels=[0.0, phi.max() if np.isfinite(phi).any() else 1.0], colors=["#000000"], alpha=0.07, zorder=29)

        # Draw a strong saturation contour (phi == 0) as white then black edge for visibility
        cs_w = ax.contour(P_mesh_GPa, T_mesh, phi, levels=[0.0], colors=["white"], linewidths=3.2, zorder=36)
        cs_b = ax.contour(P_mesh_GPa, T_mesh, phi, levels=[0.0], colors=["black"], linewidths=1.0, zorder=37)
    except Exception:
        # Fallback to the simple plotted line if contouring fails for any reason
        pass

    # Overlay sample saturation points colored by saturated vapor density (and saturated
    # liquid if available) so the saturation behavior is highlighted. Do not create a
    # separate colorbar for these markers — we use the main density colorbar.
    try:
        rho_sat_v = None
        rho_sat_l = None
        if coolprop_ok:
            from CoolProp.CoolProp import PropsSI

            # saturated vapor (Q=1) and saturated liquid (Q=0)
            rho_sat_v = np.array([PropsSI("D", "T", float(Ti), "Q", 1, "Hydrogen") for Ti in T])
            rho_sat_l = np.array([PropsSI("D", "T", float(Ti), "Q", 0, "Hydrogen") for Ti in T])
        else:
            # fallback: saturated vapor via ideal-gas estimate; saturated liquid unknown
            M_molar = 2.01588e-3
            R_universal = 8.314462618
            R_specific = R_universal / M_molar
            rho_sat_v = P / (R_specific * T)
            rho_sat_l = None

        # Plot saturated vapor points without a white outline so markers don't
        # create a halo when drawn densely. Use the same cmap/norm as the
        # background so marker colors match the colorbar.
        sc = ax.scatter(sat_x, T, c=rho_sat_v, cmap=cmap, norm=norm, s=50, edgecolors="none", linewidths=0.0, zorder=38)

        # If saturated liquid density is available, plot it as squares without an outline
        if rho_sat_l is not None:
            ax.scatter(sat_x, T, c=rho_sat_l, cmap=cmap, norm=norm, s=50, marker="s", edgecolors="none", linewidths=0.0, alpha=0.95, zorder=33)
    except Exception:
        pass

    ax.set_xscale("log")
    ax.set_xlabel("Pressure [GPa]")
    ax.set_ylabel("Temperature [K]")
    ax.set_title("Hydrogen Phase Diagram — Density field with saturation boundary")

    # Mark critical point (no legend). Annotate using axes-fraction text so the
    # label is always readable and won't be obscured by the colorbar.
    ax.scatter(Pcrit / 1e9, Tcrit, color="red", zorder=45)
    ax.annotate(
        f"Tc={Tcrit:.3f} K\nPc={Pcrit/1e9:.3e} GPa",
        xy=(Pcrit / 1e9, Tcrit),
        xycoords="data",
        xytext=(0.92, 0.92),
        textcoords="axes fraction",
        ha="right",
        va="top",
        arrowprops=dict(arrowstyle="->", color="red"),
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.85),
        zorder=50,
    )
    ax.grid(which="both", linestyle="--", alpha=0.35)

    # Add 'Gas' / 'Liquid' labels in axis fraction coordinates so they are always visible
    ax.text(0.12, 0.88, "Gas", fontsize=14, color="white", transform=ax.transAxes, bbox=dict(facecolor="black", alpha=0.45), zorder=40)
    # Place 'Liquid' label in the bottom-right axis fraction so it's visually
    # located over the liquid region (high P, low T).
    ax.text(0.78, 0.12, "Liquid", fontsize=14, color="white", transform=ax.transAxes, bbox=dict(facecolor="black", alpha=0.45), zorder=40)

    # constrained_layout=True was used when creating the figure so avoid
    # calling tight_layout (they can interfere). If a save_path was provided,
    # write the image to disk; otherwise return the Figure for embedding in
    # other code (importing this module and calling the function).
    if save_path is not None:
        fig.savefig(save_path, dpi=300)
        print(f"Saved phase diagram to {save_path} (x-axis Pressure in GPa; color = density field)")

    return fig


if __name__ == "__main__":
    # When executed as a script, save to the module's OUTFILE for backwards
    # compatibility.
    plot_phase_diagram(save_path=OUTFILE)

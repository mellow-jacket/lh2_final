"""
Geometry Module

Provides geometric calculations for cylindrical tanks including:
- Volume to height conversions for horizontal cylinders
- Surface area calculations
- Cross-sectional area computations
"""

import numpy as np


def cyl_v_to_h(V, R, L):
    """
    Convert liquid volume to liquid height in a horizontal cylinder.
    
    Based on MATLAB cylVToH.m function from LLNL model.
    Uses Newton iteration to solve for height given volume.
    
    Parameters
    ----------
    V : float
        Liquid volume [m³]
    R : float
        Cylinder radius [m]
    L : float
        Cylinder length [m]
        
    Returns
    -------
    H : float
        Liquid height from bottom of cylinder [m]
        
    Notes
    -----
    For a horizontal cylinder, the relationship between volume and height
    is nonlinear. This function iteratively solves for the height parameter.
    
    References
    ----------
    Original MATLAB implementation: LLNL_model/cylVToH.m
    """
    # Cross-sectional area of cylinder
    A = np.pi * R**2
    
    # Cross-sectional area of liquid (s = V/L in MATLAB code)
    s = V / L
    
    # Handle edge cases
    if s <= 0:
        return 0.0
    if s >= A:
        return 2 * R
    
    # Special case: exactly half full
    if abs(s - A / 2) < 1e-10:
        return R
    
    # Determine which half of cylinder we're in
    if s > A / 2:
        sup = abs(s - A)
        b = 1
    else:
        sup = s
        b = 1
    
    # Initial guess for x (distance from center to liquid surface)
    # Start with a better initial guess based on approximate geometry
    x = 0.1 * R
    max_iter = 100
    tolerance = 1e-4
    
    # Newton iteration
    # Function: f(x) = sup - b*R²*atan(sqrt(R²-x²)/x) + x*sqrt(R²-x²)
    for _ in range(max_iter):
        x_old = x
        
        # Ensure x is in valid range
        if abs(x) >= R:
            x = 0.99 * R * np.sign(x)
        
        # Compute intermediate values
        R2 = R * R
        x2 = x * x
        
        # Check if x is valid for sqrt
        if R2 - x2 < 0:
            x = 0.99 * R
            x2 = x * x
        
        y = np.sqrt(R2 - x2)  # sqrt(R²-x²)
        
        # Function value
        f = sup - b * R2 * np.arctan(y / x) + x * y
        
        # Derivative
        # alpha = -1/sqrt(R²-x²) - sqrt(R²-x²)/x²
        alpha = -1.0 / y - y / x2
        
        # supd = -b*alpha*R²/(1+(y/x)²) + sqrt(R²-x²) - x²/sqrt(R²-x²)
        supd = -b * alpha * R2 / (1 + (y / x)**2) + y - x2 / y
        
        # Newton update
        if abs(supd) < 1e-10:
            break
        x = x - f / supd
        
        # Clamp x to valid range
        x = np.clip(x, -0.999 * R, 0.999 * R)
        
        # Check convergence
        error = abs((x_old - x) / x_old) if abs(x_old) > 1e-10 else 0
        if error < tolerance:
            break
    
    # Convert x parameter to height H
    if s >= A / 2:
        H = R + x
    else:
        H = R - x
    
    return H


def cylinder_cross_section_area(R):
    """
    Calculate cross-sectional area of a cylinder.
    
    Parameters
    ----------
    R : float
        Cylinder radius [m]
        
    Returns
    -------
    A : float
        Cross-sectional area [m²]
    """
    return np.pi * R**2


def cylinder_lateral_surface_area(R, L):
    """
    Calculate lateral surface area of a cylinder.
    
    Parameters
    ----------
    R : float
        Cylinder radius [m]
    L : float
        Cylinder length [m]
        
    Returns
    -------
    A : float
        Lateral surface area [m²]
    """
    return 2.0 * np.pi * R * L


def horizontal_cylinder_liquid_surface_area(V, R, L):
    """
    Calculate liquid-vapor interface area in a horizontal cylinder.
    
    Parameters
    ----------
    V : float
        Liquid volume [m³]
    R : float
        Cylinder radius [m]
    L : float
        Cylinder length [m]
        
    Returns
    -------
    A : float
        Liquid-vapor interface area [m²]
    """
    # For horizontal cylinder, interface area is approximately rectangular
    # with width = cylinder length and height = chord length at liquid level
    h = cyl_v_to_h(V, R, L)
    
    # Chord length at height h
    if h <= 0 or h >= 2 * R:
        return 0.0
    
    # Distance from center
    y = h - R
    chord_length = 2.0 * np.sqrt(R**2 - y**2)
    
    return chord_length * L

"""
================================================================================
1D SHALLOW WATER EQUATIONS (SWE) SOLVER
================================================================================
Developer: Gabriel Thomas Scarlett (2015 / Modernized 2026)
Domain: Coastal Engineering / Hydrodynamics

MATHEMATICAL MODEL:
-------------------
This solver implements the depth-integrated Navier-Stokes equations, known 
as the Shallow Water Equations. 

1. CONTINUITY (Mass Conservation):
   ∂ζ/∂t + ∂q/∂x = 0

2. MOMENTUM (Newton's Second Law):
   ∂q/∂t + ∂/∂x(q²/h) + gh(∂ζ/∂x) + τ_b/ρ = 0

COORDINATE SYSTEM & VARIABLES:
------------------------------
Reference: Still Water Level (SWL) is the vertical datum (Z = 0).
Vertical:  Z-axis is positive UPWARDS.

   Z = ζ   ~~~~~~~~~~~~~~~~  Free Surface (Dynamic Topography)
                ^          ^
                | ζ        |
   Z = 0      ----------------  SWL / Local Geoid (Reference Datum)
                |          | 
                | hs       | h (Total Water Depth)
                v          |
   Z = -hs    ________________  Bed (Bathymetry)

Definitions:
   ζ (zeta) : Free surface elevation relative to SWL [m]
   hs       : Distance from SWL down to the bed (Still water depth) [m]
   h        : Total water depth (hs + ζ) [m]
   q        : Discharge / Flux per unit width (u * h) [m²/s]
   u        : Depth-averaged velocity [m/s]
   τ_b      : Bed shear stress (calculated via Chezy coefficient) [N/m²]

ASSUMPTIONS & LIMITATIONS:
--------------------------
1. INVISCID FLOW: The model is inviscid (Euler-based), meaning internal fluid 
   viscosity and turbulent eddy viscosity are neglected. Energy dissipation 
   is accounted for solely through the Bed Shear Stress (τ_b) term.
   
2. HYDROSTATIC PRESSURE: Vertical accelerations are assumed negligible compared 
   to gravity. Pressure distribution is assumed to be linear with depth.

3. UNIFORM VELOCITY: The velocity profile is assumed uniform over the depth. 
   Three-dimensional effects like bottom boundary layer logs are averaged.

4. ABSENCE OF CORIOLIS: The Coriolis force (Earth's rotation) is neglected. 
   In 1D channel flow, Coriolis affects the lateral tilt of the water surface 
   but does not drive the longitudinal discharge. 

5. TIDAL RESOURCE ASSESSMENT NOTE: 
   While this 1D model is sufficient for localized channel characterization, 
   larger "Macro" scale tidal assessments must incorporate:
   - 2D/3D Depth-averaged equations.
   - Coriolis Force (essential for Kelvin Waves and large-scale circulation).
   - Wind stress and Atmospheric Pressure gradients.
   - Interaction with the Geoid (for global tidal constituents).

NUMERICAL SCHEME:
-----------------
- Spatial Discretization: 2nd Order Central Differences.
- Temporal Integration: 4th Order Runge-Kutta (RK4).
- Implementation: Python / NumPy Vectorized.
================================================================================
"""

import numpy as np
import matplotlib.pyplot as plt
import time

class Config:
    """Simulation Parameters"""
    imax = 101
    lx = 1000.0
    dx = lx / (imax - 1)
    dt = 2.8
    tend = 5000.0
    
    g = 9.81
    rho = 1000.0
    ch = 40.0
    cf = g / (ch**2)
    
    h0 = 5.0      # Reference depth
    S0 = 1.0 / lx # Bed Slope

# =============================================================================
# PHYSICS ENGINE
# =============================================================================

def get_derivatives(zeta, q, hs, config):
    """Calculates time derivatives dz/dt and dq/dt using Central Differences."""
    # Local depth and velocity
    h = hs + zeta
    u = q / h
    
    dzdt = np.zeros_like(zeta)
    dqdt = np.zeros_like(q)
    
    # 1. Continuity: d(zeta)/dt = -d(q)/dx
    dzdt[1:-1] = -(q[2:] - q[:-2]) / (2.0 * config.dx)
    
    # 2. Momentum: d(q)/dt = -d(q^2/h)/dx - g*h*d(zeta)/dx - Friction
    advection = -( (q[2:]**2 / h[2:]) - (q[2-2:-2]**2 / h[:-2]) ) / (2.0 * config.dx)
    pressure  = -config.g * h[1:-1] * (zeta[2:] - zeta[:-2]) / (2.0 * config.dx)
    friction  = -(config.cf * config.rho * u[1:-1] * np.abs(u[1:-1])) / config.rho
    
    dqdt[1:-1] = advection + pressure + friction
    
    # Simple extrapolation for boundary derivatives
    dzdt[0], dzdt[-1] = dzdt[1], dzdt[-2]
    dqdt[0], dqdt[-1] = dqdt[1], dqdt[-2]
    
    return dzdt, dqdt

def rk4_step(zeta, q, hs, config):
    """Performs one 4th-order Runge-Kutta time step."""
    dt = config.dt
    
    # k1
    dz1, dq1 = get_derivatives(zeta, q, hs, config)
    
    # k2
    z2, q2 = zeta + (dt/2)*dz1, q + (dt/2)*dq1
    dz2, dq2 = get_derivatives(z2, q2, hs, config)
    
    # k3
    z3, q3 = zeta + (dt/2)*dz2, q + (dt/2)*dq2
    dz3, dq3 = get_derivatives(z3, q3, hs, config)
    
    # k4
    z4, q4 = zeta + dt*dz3, q + dt*dq3
    dz4, dq4 = get_derivatives(z4, q4, hs, config)
    
    # Weighted average update
    zeta_new = zeta + (dt/6.0) * (dz1 + 2*dz2 + 2*dz3 + dz4)
    q_new    = q    + (dt/6.0) * (dq1 + 2*dq2 + 2*dq3 + dq4)
    
    # Apply Boundary Conditions
    zeta_new[0] = -config.S0 * 0.0
    zeta_new[-1] = -config.S0 * config.lx
    q_new[0] = q_new[1]
    q_new[-1] = q_new[-2]
    
    return zeta_new, q_new

# =============================================================================
# MAIN SOLVER
# =============================================================================

""" def run_simulation():
    c = Config()
    x = np.linspace(0, c.lx, c.imax)
    
    # Initialization
    # hs: Depth from SWL down to Bed. (positive magnitude)
    hs = c.h0 + c.S0 * x
    # zeta: Elevation relative to SWL (starts as a slope)
    zeta = -c.S0 * x
    q = np.zeros(c.imax)
    
    t = 0.0
    print(f"Running SWE Solver... End Time: {c.tend}s")
    
    start_time = time.time()
    while t < c.tend:
        zeta, q = rk4_step(zeta, q, hs, c)
        t += c.dt
        
        if round(t) % 500 == 0:
            mid = c.imax // 2
            h_mid = hs[mid] + zeta[mid]
            u_mid = q[mid] / h_mid
            print(f"Time: {t:6.1f}s | Depth: {h_mid:8.5f}m | Vel: {u_mid:8.5f}m/s")
            
    print(f"Done. Simulation took {time.time() - start_time:.2f}s")
    return x, zeta, hs, q """

def run_simulation():
    c = Config()
    x = np.linspace(0, c.lx, c.imax)
    
    # --- GAUSSIAN HUMP PARAMETERS ---
    hump_amplitude = 1.5     # Height of the hump (meters)
    hump_center = 500.0      # Position of the crest (meters)
    hump_width = 60.0        # Spread of the hump
    
    # Define the hump shape
    hump = hump_amplitude * np.exp(-0.5 * ((x - hump_center) / hump_width)**2)
    
    # --- INITIALIZATION ---
    # hs (Still Water Depth): The distance from SWL down to the bed.
    # We subtract the hump because the bed is rising TOWARD the SWL.
    hs = (c.h0 + c.S0 * x) - hump
    
    # zeta: Elevation relative to SWL
    zeta = -c.S0 * x
    
    # q: Discharge (Initial guess)
    q = np.zeros(c.imax)
    
    t = 0.0
    print(f"Running Task 1: Flow over Hump | Amplitude: {hump_amplitude}m")
    
    start_time = time.time()
    while t < c.tend:
        zeta, q = rk4_step(zeta, q, hs, c)
        t += c.dt
        
        # Check for stability (Stop if the model blows up)
        if np.any(np.isnan(zeta)):
            print("Instability detected! Decrease dt.")
            break

        if round(t) % 500 == 0:
            mid = c.imax // 2
            h_mid = hs[mid] + zeta[mid]
            u_mid = q[mid] / h_mid
            print(f"Time: {t:6.1f}s | Mid-Depth: {h_mid:8.5f}m | Mid-Vel: {u_mid:8.5f}m/s")
            
    print(f"Done. Simulation took {time.time() - start_time:.2f}s")
    return x, zeta, hs, q

def plot_results(x, zeta, hs, q):
    h = hs + zeta
    u = q / h
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    
    # Top Plot: Elevations
    ax1.plot(x, zeta, 'b', label='Surface (zeta)')
    ax1.plot(x, -hs, 'r', label='Bed (-hs)')
    ax1.axhline(0, color='k', linestyle='--', label='SWL (Z=0)')
    ax1.fill_between(x, -hs, zeta, color='skyblue', alpha=0.3)
    ax1.set_ylabel('Elevation (m)')
    ax1.set_title('Shallow Water Flow: Z-coordinate up from SWL')
    ax1.legend()
    
    # Bottom Plot: Velocity
    ax2.plot(x, u, 'g', label='Velocity (u)')
    ax2.set_ylabel('Velocity (m/s)')
    ax2.set_xlabel('Distance (m)')
    ax2.set_ylim(0, 5)
    ax2.legend()
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    x, zeta, hs, q = run_simulation()
    plot_results(x, zeta, hs, q)
"""
Numerical Solution to the 1D Shallow Water Equations (SWE)
Original FORTRAN Code: Gabriel Thomas Scarlett, University of Edinburgh (2015)
Python Port & Modernization: 2026

Governing Equations (Depth-Integrated Navier-Stokes):
1. Continuity: d(zeta)/dt + d(q)/dx = 0
2. Momentum:   d(q)/dt + d(q^2/h)/dx + g*h*d(zeta)/dx + (tau_b/rho) = 0

Where:
    zeta = Free surface elevation (m)
    h    = Total water depth (hs + zeta) (m)
    q    = Discharge / Flux (u * h) (m^2/s)
    hs   = Still water level (distance from datum to bed) (m)
    tau_b= Bed shear stress (N/m^2)
"""

import numpy as np
import matplotlib.pyplot as plt
import time

# =============================================================================
# 1. PHYSICAL & NUMERICAL PARAMETERS
# =============================================================================
imax = 101          # Number of spatial nodes
lx   = 1000.0       # Domain length (meters)
g    = 9.81         # Acceleration due to gravity (m/s^2)
ch   = 40.0         # Chezy friction coefficient (m^1/2 / s)
h0   = 5.0          # Reference initial water depth (meters)
rho  = 1000.0       # Density of water (kg/m^3)
tend = 5000.0       # Total simulation time (seconds)
dt   = 2.8          # Time step (seconds) - Must satisfy Courant (CFL) condition

# Derived Parameters
S0 = 1.0 / lx       # Constant bed gradient (Slope)
cf = g / (ch**2)    # Dimensionless bed friction coefficient
dx = lx / (imax-1)  # Spatial step size (m)

# =============================================================================
# 2. INITIALIZATION (NumPy Arrays)
# =============================================================================
# x is our spatial coordinate array from 0 to 1000m
x = np.linspace(0, lx, imax)

# zeta: Free surface elevation relative to the datum. 
# In FORTRAN: zeta(i) = -S0 * x. Initialized as a slope.
zeta = -S0 * x

# hs: Still water depth (distance from bed to datum).
# In FORTRAN: hs(i) = h0 + S0 * x
hs = h0 + S0 * x

#hs = h0 - S0*x # (Bed rises from 0 to 1m)

# h: Total depth (hs + zeta). At start, it is uniform (h0).
#h = hs + zeta
#h = (h0 - S0*x) + (-S0*x) #= 5.0 - 2.0 = 3.0m

# q: Local discharge (flux). Initially zero (static water).
q = np.zeros(imax)

# =============================================================================
# 3. THE PHYSICS ENGINE (Derivative Calculation)
# =============================================================================
def compute_derivatives(z_loc, q_loc, hs_loc):
    """
    This function replaces the 'calc' subroutine logic.
    It calculates the time-derivatives (dz/dt and dq/dt) using 
    Vectorized Central Differences.
    
    Vectorization Explained:
    Instead of 'do i=2, imax-1', we use slices:
    [1:-1] = The interior nodes (index 1 to 99)
    [2:]   = The 'East' nodes (i+1 in FORTRAN)
    [:-2]  = The 'West' nodes (i-1 in FORTRAN)
    """
    
    # Calculate local total depth and velocity
    h_loc = hs_loc + z_loc
    u = q_loc / h_loc
    
    # Initialize derivative arrays with zeros
    dzdt = np.zeros_like(z_loc)
    dqdt = np.zeros_like(q_loc)
    
    # --- A. CONTINUITY EQUATION (Mass Conservation) ---
    # dz/dt = -dq/dx
    # Central difference: (q_east - q_west) / (2 * dx)
    dzdt[1:-1] = -(q_loc[2:] - q_loc[:-2]) / (2.0 * dx)
    
    # --- B. MOMENTUM EQUATION (Conservation of Momentum) ---
    # 1. Advection Term: d(q^2/h)/dx
    advection = -( (q_loc[2:]**2 / h_loc[2:]) - (q_loc[:-2]**2 / h_loc[:-2]) ) / (2.0 * dx)
    
    # 2. Pressure/Surface Gradient Term: g * h * d(zeta)/dx
    pressure = -g * h_loc[1:-1] * (z_loc[2:] - z_loc[:-2]) / (2.0 * dx)
    
    # 3. Friction Term (Bed Shear Stress): taubx / rho
    # Based on Chezy: taubx = cf * rho * u * |u|
    friction = -(cf * rho * u[1:-1] * np.abs(u[1:-1])) / rho
    
    # Combine terms to get dq/dt
    dqdt[1:-1] = advection + pressure + friction
    
    # --- C. BOUNDARY CONDITIONS (Internal extrapolation) ---
    # Your FORTRAN code set edges to the value of the nearest neighbor
    dzdt[0] = dzdt[1]; dzdt[-1] = dzdt[-2]
    dqdt[0] = dqdt[1]; dqdt[-1] = dqdt[-2]
    
    return dzdt, dqdt

# =============================================================================
# 4. RUNGE-KUTTA 4th ORDER TIME INTEGRATION
# =============================================================================
print(f"Starting Simulation: {tend}s, dt={dt}s")
print(f"{'Time (s)':>10} | {'Mid-Depth (m)':>15} | {'Mid-Velocity (m/s)':>18}")
print("-" * 50)

start_wall_time = time.time()
t = 0.0

# Store history for potential steady-state check
while t < tend:
    # --- RK4 FOUR-STEP PROCESS ---
    # k1 = f(t, y)
    dz1, dq1 = compute_derivatives(zeta, q, hs)
    
    # k2 = f(t + dt/2, y + k1*dt/2)
    z2, q2 = zeta + (dt/2.0)*dz1, q + (dt/2.0)*dq1
    dz2, dq2 = compute_derivatives(z2, q2, hs)
    
    # k3 = f(t + dt/2, y + k2*dt/2)
    z3, q3 = zeta + (dt/2.0)*dz2, q + (dt/2.0)*dq2
    dz3, dq3 = compute_derivatives(z3, q3, hs)
    
    # k4 = f(t + dt, y + k3*dt)
    z4, q4 = zeta + dt*dz3, q + dt*dq3
    dz4, dq4 = compute_derivatives(z4, q4, hs)
    
    # --- FINAL WEIGHTED UPDATE ---
    # zeta_new = zeta + dt/6 * (k1 + 2k2 + 2k3 + k4)
    zeta_change = (dt/6.0) * (dz1 + 2.0*dz2 + 2.0*dz3 + dz4)
    q_change    = (dt/6.0) * (dq1 + 2.0*dq2 + 2.0*dq3 + dq4)
    
    zeta += zeta_change
    q    += q_change
    
    # --- ENFORCE RIGID BOUNDARY CONDITIONS ---
    # Inlet/Outlet elevations (zeta) remain fixed as per original 'calc'
    zeta[0] = -S0 * x[0]
    zeta[-1] = -S0 * x[-1]
    
    # Boundary condition for Discharge q (Zero gradient / flow continuity)
    q[0] = q[1]
    q[-1] = q[-1] # Fixed as per FORTRAN update loop logic
    
    t += dt
    
    # --- DIAGNOSTICS & PRINTING ---
    if round(t) % 500 == 0:
        h_current = hs + zeta
        imid = imax // 2
        velocity = q[imid] / h_current[imid]
        print(f"{t:10.1f} | {h_current[imid]:15.6f} | {velocity:18.6f}")

end_wall_time = time.time()
print("-" * 50)
print(f"Simulation Complete in {end_wall_time - start_wall_time:.4f} seconds.")

# --- 5. VISUALIZATION (Corrected Coordinate System) ---
h_final = hs + zeta
velocity_final = q / h_final

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

# Plot 1: Physical Geometry (Surface and Bed relative to Datum Z=0)
# We plot zeta as the surface and -hs as the bed to show the downward slope
ax1.plot(x, zeta, 'b', linewidth=2, label='Water Surface Elevation (zeta)')
ax1.plot(x, -hs, 'saddlebrown', linewidth=2, label='Seabed Elevation (-hs)')

# Fill the area between them to represent the water body
ax1.fill_between(x, -hs, zeta, color='skyblue', alpha=0.4, label='Water Column (h=5m)')

ax1.set_ylabel('Elevation relative to Datum (m)')
ax1.set_title('1D Shallow Water Flow: Uniform Parallel Downward Slope')
ax1.legend(loc='upper right')
ax1.grid(True, linestyle='--', alpha=0.6)

# Plot 2: Velocity Profile
ax2.plot(x, velocity_final, 'r', linewidth=2, label='Velocity (u)')
ax2.set_xlabel('Distance along Channel (m)')
ax2.set_ylabel('Velocity (m/s)')

# Fix the Y-axis so we don't see the 1e-12 numerical noise
# We set the limits around the theoretical value 2.828
ax2.set_ylim(0, 5) 

ax2.legend(loc='upper right')
ax2.grid(True, linestyle='--', alpha=0.6)

plt.tight_layout()
plt.show()
"""
================================================================================
1D SHALLOW WATER EQUATIONS (SWE) SOLVER
================================================================================
Developer: Gabriel Thomas Scarlett (2015 / Modernised 2026)
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
import warnings
import time

class SWESolver1D:
    def __init__(self, 
                 hump_amplitude=1.5,
                 node_count=101, 
                 domain_length=1000.0,
                 target_cfl=0.5):
        """
        Initializes the physical domain and numerical parameters.
        
        Args:
            hump_amplitude (float): Height of the Gaussian hump [m].
            node_count (int): Number of spatial nodes (imax).
            domain_length (float): Length of the channel [m].
        """
        # --- Physical Constants ---
        self.g = 9.81
        self.rho = 1000.0
        self.ch = 40.0              # Chezy Coefficient
        self.cf = self.g / (self.ch**2) # Friction Coefficient
        self.h0 = 5.0               # Nominal Depth [m]
        
        # --- Numerical Parameters ---
        self.imax = node_count
        self.lx = domain_length
        self.dx = self.lx / (self.imax - 1)
        self.dt = 2.8               # Time Step [s] (CFL sensitive)
        
        # --- Spatial Grid & Bathymetry ---
        self.x = np.linspace(0, self.lx, self.imax)
        self.S0 = 1.0 / self.lx     # Bed Slope
        
        # Gaussian Hump Definition
        hump_center = 500.0
        hump_width = 60.0
        hump = hump_amplitude * np.exp(-0.5 * ((self.x - hump_center) / hump_width)**2)
        
        # hs: Distance from SWL down to the bed (positive magnitude)
        # Bed falls away with S0*x, but rises with the hump.
        self.hs = (self.h0 + self.S0 * self.x) - hump
        
        # --- Initial State Variables ---
        # zeta: Surface elevation relative to SWL (Z=0)
        self.zeta = -self.S0 * self.x

        # Starting with a uniform flux based on inflow_velocity
        self.q = np.full(self.imax, 0.0)
        
        self.t = 0.0

        self.target_cfl = target_cfl
        self.dt = 0.1 # Initial guess, will be updated immediately

    def _update_dt(self):
        """
        Calculates the maximum allowable time step based on the 
        CFL condition: dt = (CFL * dx) / (u + sqrt(gh))
        """
        h = self.hs + self.zeta
        u = np.abs(self.q / h) 
        c = np.sqrt(self.g * h) # the celerity (wave speed)
        
        # Find the maximum 'information speed' across the whole domain
        max_speed = np.max(u + c)
        
        # Calculate new dt
        self.dt = (self.target_cfl * self.dx) / max_speed

    def _compute_derivatives(self, z_loc, q_loc):
        """
        Calculates dz/dt and dq/dt with Artificial Viscosity 
        to prevent 'jagged' profiles in near-critical flow.
        """
        h_loc = self.hs + z_loc
        u_loc = q_loc / h_loc
        
        dzdt = np.zeros_like(z_loc)
        dqdt = np.zeros_like(q_loc)
        
        # --- 1. Standard Physics Terms ---
        dzdt[1:-1] = -(q_loc[2:] - q_loc[:-2]) / (2.0 * self.dx)
        
        advection = -( (q_loc[2:]**2 / h_loc[2:]) - (q_loc[:-2]**2 / h_loc[:-2]) ) / (2.0 * self.dx)
        pressure  = -self.g * h_loc[1:-1] * (z_loc[2:] - z_loc[:-2]) / (2.0 * self.dx)
        friction  = -(self.cf * self.rho * u_loc[1:-1] * np.abs(u_loc[1:-1])) / self.rho
        
        # --- 2. NEW: Artificial Viscosity (Numerical Smoothing) ---
        # nu_art is the viscosity coefficient. 
        # Increase this slightly if the profile is still jagged.
        nu_art = 0.05
        
        # 2nd-order diffusion (Laplacian) of discharge and elevation
        # This damps high-frequency oscillations (the wiggles)
        diffusion_q = nu_art * (q_loc[2:] - 2*q_loc[1:-1] + q_loc[:-2]) / self.dx
        diffusion_z = nu_art * (z_loc[2:] - 2*z_loc[1:-1] + z_loc[:-2]) / self.dx

        # --- 3. Combine ---
        dqdt[1:-1] = advection + pressure + friction + diffusion_q
        dzdt[1:-1] += diffusion_z
        
        # BCs
        dzdt[0], dzdt[-1] = dzdt[1], dzdt[-2]
        dqdt[0], dqdt[-1] = dqdt[1], dqdt[-2]
        
        return dzdt, dqdt

    def _rk4_step(self):
        """Executes a single 4th-order Runge-Kutta time integration step."""
        # k1
        dz1, dq1 = self._compute_derivatives(self.zeta, self.q)
        
        # k2
        z2, q2 = self.zeta + (self.dt/2.0)*dz1, self.q + (self.dt/2.0)*dq1
        dz2, dq2 = self._compute_derivatives(z2, q2)
        
        # k3
        z3, q3 = self.zeta + (self.dt/2.0)*dz2, self.q + (self.dt/2.0)*dq2
        dz3, dq3 = self._compute_derivatives(z3, q3)
        
        # k4
        z4, q4 = self.zeta + self.dt*dz3, self.q + self.dt*dq3
        dz4, dq4 = self._compute_derivatives(z4, q4)
        
        # Final Weighted Update
        self.zeta += (self.dt/6.0) * (dz1 + 2.0*dz2 + 2.0*dz3 + dz4)
        self.q    += (self.dt/6.0) * (dq1 + 2.0*dq2 + 2.0*dq3 + dq4)
        
        # Enforce Hard Boundary Conditions
        self.zeta[0] = -self.S0 * self.x[0]   # Fixed Inlet Elevation
        self.zeta[-1] = -self.S0 * self.x[-1] # Fixed Outlet Elevation
        self.q[0] = self.q[1]                # Zero-gradient inflow
        self.q[-1] = self.q[-2]              # Zero-gradient outflow

    def run(self, tend=5000.0, target_cfl=None, verbose=True):
        """
        Runs the simulation until the target end time is reached.
        Uses the instance's target_cfl unless a new one is provided here.
        """
        # If the user provides a cfl in the run call, update the instance variable
        if target_cfl is not None:
            self.target_cfl = target_cfl
            
        if verbose:
            print(f"Adaptive Simulation started (Target CFL: {self.target_cfl})")
        
        start_wall = time.time()
        max_fr_detected = 0.0  # Track the peak Froude number
        
        while self.t < tend:
            # --- 1. Physics Calculations ---
            h_current = self.hs + self.zeta
            u_current = np.abs(self.q / h_current)
            celerity = np.sqrt(self.g * h_current)
            
            # --- 2. Froude Number Check ---
            # Fr = u / c
            fr_array = u_current / celerity
            current_max_fr = np.max(fr_array)
            if current_max_fr > max_fr_detected:
                max_fr_detected = current_max_fr
            
            # --- 3. Adaptive dt (CFL) ---
            max_signal_speed = np.max(u_current + celerity)
            new_dt = (self.target_cfl * self.dx) / max_signal_speed
            
            if self.t + new_dt > tend:
                new_dt = tend - self.t
            self.dt = new_dt
            
            # --- 4. Step & Stability ---
            self._rk4_step()
            self.t += self.dt
            
            if np.any(np.isnan(self.zeta)):
                raise ValueError(f"Instability at t={self.t:.2f}s. Max Fr reached: {max_fr_detected:.2f}")

        # --- 5. Post-Simulation Warnings ---
        if max_fr_detected >= 1.0:
            warnings.warn(f"\n[WARNING] Supercritical flow detected! (Max Fr = {max_fr_detected:.2f})\n"
                          "Central difference schemes are dispersive in this regime.\n"
                          "Results may contain numerical oscillations (jags).", UserWarning)
        elif max_fr_detected > 0.8:
            print(f"Note: Near-critical flow detected (Max Fr = {max_fr_detected:.2f}).")

        end_wall = time.time()
        if verbose:
            print(f"Convergence reached at T={self.t:.1f}s. Max Fr: {max_fr_detected:.2f}")

    def get_results(self):
        """
        Returns the final steady-state profiles.
        Returns: (x, zeta, h, u)
        """
        h_final = self.hs + self.zeta
        u_final = self.q / h_final
        return self.x, self.zeta, h_final, u_final

# ================================================================================
# TEST SCRIPT (Only runs if this file is executed directly)
# ================================================================================
if __name__ == "__main__":
    # Example usage: Flow over a 1.0m hump
    solver = SWESolver1D(hump_amplitude=1.7)
    solver.run(tend=5000, target_cfl=0.5)
    x, zeta, h, u = solver.get_results()
    
    import matplotlib.pyplot as plt
    plt.plot(x, u)
    plt.title("Steady State Velocity Profile")
    plt.show()
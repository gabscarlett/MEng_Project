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

FLOW CONDITIONS & REGIMES:
--------------------------
1. SUBCRITICAL FLOW: The model is designed for subcritical regimes (Froude < 1.0).
   The solver monitors the Froude Number (Fr = u / √(gh)) locally. 
   
2. BERNOULLI EFFECT: The solver captures non-linear surface deformation over 
   bathymetric obstructions (e.g., Gaussian humps), where localised 
   acceleration results in a free-surface "dip."

3. ADAPTIVE STABILITY: Time-integration is governed by the CFL condition. 
   The time step (dt) is updated dynamically based on the combined advective 
   velocity and wave celerity (u + √(gh)).

BOUNDARY CONDITIONS (FLUX-FORCED / TIDAL):
------------------------------------------
The model utilizes a "Flume-Type" boundary configuration to simulate a forced 
tidal strait or channel:

1. INLET (x = 0): Dirichlet for Discharge (q). A constant mass flux is 
   prescribed based on a target 'Reference Velocity' (u_ref). The surface 
   elevation (ζ) is allowed to float (Neumann), allowing the developement of 
   a 'Backwater Curve' required to overcome friction and blockage.

2. OUTLET (x = L): Dirichlet for Surface Elevation (ζ). The sea level is 
   pinned to the datum (Z = 0) to represent a downstream oceanic control. 
   The discharge (q) is transmissive (Neumann).

3. PARAMETRIC NOTE: Due to the backwater effect, u_ref physically manifests 
   at the outlet where h is at the reference depth. Inlet velocity will 
   naturally be lower than u_ref to maintain mass continuity.

ASSUMPTIONS & LIMITATIONS:
--------------------------
1. INVISCID FLOW: Internal fluid and turbulent eddy viscosity are neglected. 
   Energy dissipation is handled solely through Bed Shear Stress (τ_b).
   
2. HYDROSTATIC PRESSURE: Vertical accelerations are assumed negligible.

3. UNIFORM VELOCITY: The vertical velocity profile is assumed uniform.

4. ABSENCE OF CORIOLIS: Neglected for 1D longitudinal channel flow.


NUMERICAL SCHEME & STABILIZATION:
---------------------------------
1. DISCRETIZATION: 2nd Order Central Differences for spatial gradients and 
   4th-order Runge-Kutta (RK4) for temporal integration.
   
2. ADAPTIVE STEPPING: Time-step (dt) is dynamically updated via the CFL 
   condition (Target CFL = 0.5) to ensure stability as information celerity 
   increases over bathymetric obstructions.

3. NUMERICAL STABILIZATION (ARTIFICIAL VISCOSITY): 
   To mitigate dispersive oscillations (Gibbs phenomenon) and "jagged" profiles 
   inherent in non-dissipative central difference schemes, especially as the 
   Froude number approaches unity, a 2nd-order Artificial Viscosity term 
   is implemented.
   
   - Coefficient: nu_art = 0.05
   - Conceptual Basis: This provides explicit numerical diffusion, functionally 
     equivalent to the implicit dissipation found in 1st-order Upwind schemes, 
     but with the surgical control required to preserve the higher-order 
     accuracy of the central difference discretization.

4. IMPLEMENTATION: Fully vectorized Python / NumPy for performance parity 
   with legacy compiled FORTRAN architectures.
================================================================================
"""

import numpy as np
import warnings
import time

class SWESolver1D:
    def __init__(self, 
                 hump_amplitude=1.5,
                 reference_velocity=2.0, 
                 node_count=501, 
                 domain_length=5000.0,
                 target_cfl=0.5):
        """
        Initializes the physical domain and numerical parameters.
        
        Args:
            hump_amplitude (float): Height of the Gaussian hump [m].
            inflow_velocity (float): Initial velocity across the domain [m/s].
            node_count (int): Number of spatial nodes (imax).
            domain_length (float): Length of the channel [m].
        """
        # --- Physical Constants ---
        self.g = 9.81
        self.rho = 1000.0
        self.ch = 40.0              # Chezy Coefficient
        self.cf = self.g / (self.ch**2) # Friction Coefficient
        self.h0 = 30.0               # Nominal Depth [m]
        
        # --- Numerical Parameters ---
        self.imax = node_count
        self.lx = domain_length
        self.dx = self.lx / (self.imax - 1)
        self.dt = 2.8               # Time Step [s] (CFL sensitive)
        
        # --- Spatial Grid & Bathymetry ---
        self.x = np.linspace(0, self.lx, self.imax)
        self.S0 = 0.0 # flat channel
        
        # Gaussian Hump Definition
        hump_center = domain_length / 2
        hump_width = 100.0
        hump = hump_amplitude * np.exp(-0.5 * ((self.x - hump_center) / hump_width)**2)
        
        # hs: Distance from SWL down to the bed (positive magnitude)
        # Bed rises up from h0 towards Z=0
        self.hs = self.h0 - hump
        
        # --- Initial State Variables ---
        
        self.zeta = np.zeros(self.imax) # Start with a perfectly flat surface

        # Forcing: inlet flux based on velocity
        self.q_in = reference_velocity * self.h0
        
        # q: Discharge per unit width (u * h)
        # Start with a uniform flux based on inflow_velocity
        self.q = np.full(self.imax, self.q_in)
        
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
        
        # --- 2. Artificial Viscosity (Numerical Smoothing) ---
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
        
        # 1. INLET (Upstream Control): 
        # We force the mass flux (The Tide).
        # We let the elevation float so the "piling up" can happen naturally.
        self.q[0] = self.q_in
        self.zeta[0] = self.zeta[1] 

        # 2. OUTLET (The Downstream Anchor):
        # We fix the Sea Level (Z=0). This prevents the channel from draining.
        self.zeta[-1] = 0.0 
        self.q[-1] = self.q[-2]

    def run(self, tend=5000.0, target_cfl=None, verbose=True):
        """
        Runs the simulation until the target end time is reached.
        Uses the instance's target_cfl unless a new one is provided here.
        Monitors the Froude Number and provides warnings near and above critical.
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
    
    def plot_results(self, x, zeta, h, u):
        import matplotlib.pyplot as plt
        hs = h - zeta
        #h = hs + zeta
        #u = q / h
        
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

# ================================================================================
# TEST SCRIPT (Only runs if this file is executed directly)
# ================================================================================
if __name__ == "__main__":
    # Example usage: Flow over a 1.0m hump
    solver = SWESolver1D(hump_amplitude=10.7, reference_velocity=3.0)
    solver.run(tend=1000, target_cfl=0.5)
    x, zeta, h, u = solver.get_results()
    # Print inlet and outlet velocities
    print(f"inlet velocity={u[0]:.2f} m/s")
    print(f"outlet velocity={u[-1]:.2f} m/s")
    # Plot the result
    solver.plot_results(x, zeta, h, u)
    

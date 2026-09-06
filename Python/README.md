# 1D Shallow Water Equations (SWE) Solver 🌊

**Developer:** Gabriel Thomas Scarlett (2015 / Modernized 2026)  
**Domain:** Coastal Engineering / Hydrodynamics  

**Modern Python refactor** of the original Fortran solver `SWE_CHANNEL.f95` (see `../Fortran/SWE/`).  

A robust, fully vectorized Python solver for the 1D depth-integrated Navier-Stokes equations, commonly known as the Shallow Water Equations. This model simulates free-surface flow over varying bathymetry, specifically designed to model subcritical tidal channel dynamics and backwater curves.

---

## 🚀 Features

* **High-Order Numerical Scheme:** Utilizes 2nd-order Central Differences for spatial gradients and 4th-order Runge-Kutta (RK4) for temporal integration.
* **Adaptive Time Stepping:** Dynamically updates the time step ($\Delta t$) using the Courant-Friedrichs-Lewy (CFL) condition to ensure stability across varying wave celerities.
* **Artificial Viscosity:** Implements a 2nd-order Laplacian smoothing term to mitigate Gibbs phenomenon and jagged profiles near critical flow regimes ($Fr \approx 1$).
* **Froude Monitoring:** Actively tracks the local Froude number, warning the user if the flow transitions to a supercritical regime.
* **High Performance:** Written in fully vectorized `NumPy` to achieve performance parity with legacy compiled FORTRAN architectures.

---

## 📦 Dependencies

To run this solver, you will need Python 3.x and the following libraries:
* `numpy` (for matrix/vector operations)
* `matplotlib` (for visualizing the results)

```bash
pip install numpy matplotlib
```

---

## 💻 Quick Start

You can run the solver with just a few lines of code. By default, it simulates flow over a Gaussian bathymetric hump.

```python
from SWESolver1D import SWESolver1D

# Initialize the solver with physical/numerical parameters
solver = SWESolver1D(
    hump_amplitude=10.7,   # Height of the bathymetric hump [m]
    reference_velocity=3.0, # Target reference velocity [m/s]
    node_count=501,        # Number of spatial nodes
    domain_length=5000.0   # Length of the channel [m]
)

# Run the simulation 
solver.run(tend=1000, target_cfl=0.5)

# Retrieve data arrays
x, zeta, h, u = solver.get_results()

# Print metrics
print(f"Inlet Velocity: {u[0]:.2f} m/s")
print(f"Outlet Velocity: {u[-1]:.2f} m/s")

# Visualize the free surface, bed, and velocity profile
solver.plot_results(x, zeta, h, u)
```

---

## 🧮 Mathematical Model

This solver models the Conservation of Mass (Continuity) and Conservation of Momentum:

**1. Continuity (Mass Conservation):**
$$ \frac{\partial \zeta}{\partial t} + \frac{\partial q}{\partial x} = 0 $$

**2. Momentum (Newton's Second Law):**
$$ \frac{\partial q}{\partial t} + \frac{\partial}{\partial x}\left(\frac{q^2}{h}\right) + gh\frac{\partial \zeta}{\partial x} + \frac{\tau_b}{\rho} = 0 $$

Where:
* $\zeta$: Free surface elevation relative to SWL [m]
* $q$: Discharge / Flux per unit width ($u \times h$) [m²/s]
* $h$: Total water depth [m]
* $\tau_b$: Bed shear stress (calculated via Chezy coefficient) [N/m²]
* $\rho$: Fluid density [kg/m³]
* $g$: Acceleration due to gravity [m/s²]

---

## 📐 Coordinate System & Variables

The model uses the **Still Water Level (SWL)** as the vertical datum ($Z = 0$), with the Z-axis positive **upwards**.

```text
   Z = ζ   ~~~~~~~~~~~~~~~~  Free Surface (Dynamic Topography)
                ^          ^
                | ζ        |
   Z = 0      ----------------  SWL / Local Geoid (Reference Datum)
                |          | 
                | hs       | h (Total Water Depth)
                v          |
   Z = -hs    ________________  Bed (Bathymetry)
```

---

## 🌊 Flow Conditions & Boundary Forcing

### Regimes
* **Subcritical Flow:** Designed for $Fr < 1.0$. The Bernoulli effect is naturally captured, demonstrating non-linear surface deformation ("dips") over obstructions due to localized acceleration.

### Boundary Conditions ("Flume-Type")
The domain acts as a forced tidal strait or channel:
1. **Inlet ($x = 0$):** Fixed discharge ($q$) based on the target reference velocity. Surface elevation ($\zeta$) floats via a Neumann condition, allowing a natural backwater curve to develop.
2. **Outlet ($x = L$):** Pinned sea-level ($\zeta = 0$) acting as a downstream oceanic control. Discharge is transmissive (Neumann).

---

## ⚠️ Assumptions & Limitations

1. **Inviscid Flow:** Internal fluid and turbulent eddy viscosity are neglected. Energy dissipation is handled strictly through Bed Shear Stress.
2. **Hydrostatic Pressure:** Vertical accelerations are assumed negligible.
3. **Uniform Velocity:** The vertical velocity profile is assumed uniform (depth-averaged).
4. **No Coriolis Effect:** Rotational forces are neglected for 1D longitudinal channel flow.

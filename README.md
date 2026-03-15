```markdown
# MEng_Project: Ocean Wave Propagation in Shallow Waters

This repository contains the **Fortran** numerical code developed for my MEng Thesis, **"Ocean Wave Propagation in Shallow Waters"** (2015), at the University of Edinburgh.

## Description
The project involved creating a numerical model from scratch to solve **Sobey’s phase-resolving integral wave evolution equations**, as well as a solver for the standard **Shallow Water Equations (SWE)**. The model was designed to predict the behaviour of ocean waves as they travel through coastal waters, providing solutions for both shallow and intermediate water depths with increased accuracy compared to conventional models.

## Key Features
*   **Numerical Method:** Solves equations using a **central differencing scheme**.
*   **Mixed Derivatives:** Includes a matrix inversion method to separate mixed temporal and spatial derivative terms.
*   **Time Integration:** Uses **Runge-Kutta 4th order (RK4)** for marching the solution forward in time.
*   **Validation:** Tested against sloshing waves in a tank and the dispersion of an initial mound (the Cauchy-Poisson problem).
*   **Potential Applications:** Tsunami propagation modelling and predicting tidal current velocities for site resource assessments.

## Repository Structure
*   **`/SWE`**: Contains the finite difference solver for the standard 1HD Shallow Water Equations.
*   **`/Sobey`**: Contains the finite difference solver for Sobey's phase-resolving integral wave evolution equations.
*   **`MEng Thesis - Gabriel T Scarlett.pdf`**: The full thesis document providing mathematical derivations, numerical implementation details, and result discussions.

## Academic Context
*   **Author:** Gabriel T. Scarlett
*   **Degree:** Master of Mechanical Engineering with Renewable Energy
*   **Submission Date:** April 2015
*   **Supervisor:** Professor Alistair Borthwick
```

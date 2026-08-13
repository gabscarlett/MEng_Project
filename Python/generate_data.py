import numpy as np
import pandas as pd
import time
from SWESolver1D import SWESolver1D

def run_parametric_study():
    # --- Configuration ---
    amplitudes = np.linspace(0.0, 1.75, 6)   # Hump heights
    velocities = np.linspace(0.5, 3.0, 6)   # Tidal forcing
    node_count = 201                        # Option A: High Resolution
    t_end = 2500.0                          # Sufficient for steady state
    
    all_results = []
    
    print(f"Starting Data Factory: {len(amplitudes) * len(velocities)} scenarios.")
    start_total = time.time()

    for amp in amplitudes:
        for u_in in velocities:
            print(f" Simulating: Amp={amp:.2f}m, U_in={u_in:.2f}m/s...", end="")
            
            # Initialize the solver
            # target_cfl=0.5 and nu_art=0.05 inside the class
            solver = SWESolver1D(hump_amplitude=amp, 
                                 reference_velocity=u_in, 
                                 node_count=node_count)
            
            try:
                # Run to steady state
                solver.run(tend=t_end, verbose=False)
                
                # Extract profiles
                x, zeta, h, u = solver.get_results()
                
                # Calculate Froude number for the record
                celerity = np.sqrt(9.81 * h)
                fr = u / celerity
                
                # Append each node as a row
                for i in range(len(x)):
                    all_results.append({
                        'x': x[i],
                        'hump_amp': amp,
                        'u_inflow': u_in,
                        'target_zeta': zeta[i],
                        'target_u': u[i],
                        'target_h': h[i],
                        'froude': fr[i]
                    })
                print(f" Done. Max Fr: {np.max(fr):.2f}")
                
            except Exception as e:
                print(f" FAILED: {e}")

    # --- Save to Disk ---
    df = pd.DataFrame(all_results)
    df.to_csv("tidal_surrogate_data.csv", index=False)
    
    end_total = time.time()
    print("-" * 30)
    print(f"Dataset complete: {len(df)} samples saved to 'tidal_surrogate_data.csv'.")
    print(f"Total processing time: {end_total - start_total:.2f}s")

if __name__ == "__main__":
    run_parametric_study()
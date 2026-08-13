import numpy as np
import pandas as pd
import time
from SWESolver1D import SWESolver1D

def generate_hydro_dataset():
    # --- 1. Define Parameter Space ---
    # 10 values for amplitude and 10 for velocity = 100 simulations total
    amplitudes = np.linspace(0.0, 1.7, 10)
    u_refs = np.linspace(0.1, 3.0, 10) # Start at 0.1 to avoid div by zero in empty channel
    
    dataset = []
    sim_count = 1
    total_sims = len(amplitudes) * len(u_refs)
    
    print(f"Starting data generation: {total_sims} simulations planned.")
    start_wall = time.time()

    # --- 2. Production Loop ---
    for amp in amplitudes:
        for u_in in u_refs:
            # Initialize solver for this specific scenario
            solver = SWESolver1D(hump_amplitude=amp, reference_velocity=u_in)
            
            try:
                # Run to steady state (5000s is usually plenty for convergence)
                # We use a lower verbose setting to keep the console clean
                solver.run(tend=5000, target_cfl=0.5, verbose=False)
                
                # Extract final converged profiles
                x, zeta, h, u = solver.get_results()
                
                # Check for physical sanity (Froude check)
                # We know Fr > 1.0 causes wiggles; we'll tag the data
                max_fr = np.max(u / np.sqrt(9.81 * h))
                
                # Log only if the simulation was stable (no NaNs)
                if not np.any(np.isnan(zeta)):
                    for i in range(len(x)):
                        dataset.append({
                            'x': x[i],
                            'hump_amp': amp,
                            'u_ref': u_in,
                            'zeta': zeta[i],
                            'u': u[i],
                            'h': h[i],
                            'max_fr': max_fr
                        })
                    
                    if sim_count % 10 == 0:
                        print(f"Progress: {sim_count}/{total_sims} simulations complete.")
                
            except ValueError as e:
                print(f"Skipping scenario [Amp:{amp}, U:{u_in}]: {e}")
            
            sim_count += 1

    # --- 3. Save to Disk ---
    df = pd.DataFrame(dataset)
    
    # Professional Portfolio Tip: Save as Parquet for speed, or CSV for readability
    filename = "tidal_surrogate_data.csv"
    df.to_csv(filename, index=False)
    
    end_wall = time.time()
    print("-" * 50)
    print(f"DATA GENERATION COMPLETE")
    print(f"Total points collected: {len(df)}")
    print(f"File saved: {filename}")
    print(f"Total execution time: {end_wall - start_wall:.2f} seconds")

if __name__ == "__main__":
    generate_hydro_dataset()
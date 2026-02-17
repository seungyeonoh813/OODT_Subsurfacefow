import matlab.engine
import numpy as np
import matplotlib.pyplot as plt
import os
import time

T = 10.0            
dt = 0.1            
max_steps = int(T / dt)
test_k = 100.0      
test_init_p = 100.0 
test_u = 150.0      

time_log = []
max_p_log = []
prod_log = []

print("MATLAB engine starting...")
start_time = time.time()
eng = matlab.engine.start_matlab()
eng.addpath(os.getcwd(), nargout=0)
print(f"Engine started in {time.time() - start_time:.2f} seconds.")

try:
    print("Initialize Environment (is_reset=True)")
    
    res = eng.random_initial_k(
        float(0.0), float(dt), True, float(test_init_p), float(test_k), nargout=1 
    )
    p_array = np.array(res['pressure']).flatten() / 1e5
    prod = float(res['prod_rate'])
    
    current_time = 0.0

    print(f"Running steps (Total {max_steps} steps)")
    for step in range(max_steps):
        current_time += dt
        
        res = eng.random_initial_k(
            float(test_u), float(dt), False, 0.0, 0.0, nargout=1
        )
        
        p_array = np.array(res['pressure']).flatten() / 1e5
        max_p = np.max(p_array)
        prod = float(res['prod_rate'])
        
        time_log.append(current_time)
        max_p_log.append(max_p)
        prod_log.append(abs(prod))
        
        if (step + 1) % 10 == 0:
            print(f"Step {step+1}/{max_steps} | Max P: {max_p:.2f} Bar | Prod Rate: {abs(prod):.2f} m3/day")

finally:
    eng.quit()
    print("MATLAB engine closed.")

plt.figure(figsize=(10, 8))
plt.subplot(2, 1, 1)
plt.plot(time_log, max_p_log, 'b-', linewidth=2)
plt.axhline(test_init_p, color='k', linestyle='--', label='Initial Pressure')
plt.title(f"Pressure (Injection = {test_u} $m^3/day$)")
plt.xlabel("Time (Days)")
plt.ylabel("Max Pressure (Bar)")
plt.grid(True)
plt.legend()

plt.subplot(2, 1, 2)
plt.plot(time_log, prod_log, 'r-', linewidth=2)
plt.title("Production")
plt.xlabel("Time (Days)")
plt.ylabel("Production Rate ($m^3/day$)")
plt.grid(True)

plt.tight_layout()
plt.show()
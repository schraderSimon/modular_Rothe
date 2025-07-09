#!/usr/bin/env python3
# --------------------------------------------------------------------
import numpy as np
from HOscillator_1D import HOscillator_1D
from RotheSolver     import RotheSolver     # includes the I/O helpers

# ---------- initial wave-function parameters ------------------------
dt     = 0.01
omega  = 1.0
a=1/np.sqrt(3)
b=0.1
mu=1
q= 0
p_init = np.array([[a,b,mu,q], [a,b,-mu,q],[a,b,2*mu,q],[a,b,-2*mu,q]                  
                   ])          # p
c_init = np.array([1/np.sqrt(1),1,1,1])  # c

# ---------- build WF + solver ---------------------------------------
wf     = HOscillator_1D(p_init, c_init,
                        omega=omega, dt=dt, name="gauss_x2")
solver = RotheSolver()
solver.initialize_WF_object(wf)

# ---------- propagate t = 0 … 10 ------------------------------------
N_steps = int(10.0 / dt)            # 1 000 steps
solver.propagate_Nsteps_and_store(
        t0=0.0,
        N=N_steps,
        dt=dt,
        filename="HO_1D_trajectory.h5",   # will be created
        group_prefix="step",
        starting_guess_method="difference")        # groups: /step_000000 … /step_001000

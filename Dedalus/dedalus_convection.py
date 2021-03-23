

import time
import numpy as np
from dedalus import public as de

import logging
logger = logging.getLogger(__name__)


# Parameters
Lx = Ly = 256
Nx = Ny = 64
Lz = 128
Nz = 64
N2 = 1e-5
Q0 = 1e-8
ν = 1e-3
κ = 1e-3
f = 1e-4
dt = 60
stop_sim_time = 48 * 60 * 60

# Domain
x_basis = de.Fourier('x', Nx, interval=(0, Lx), dealias=3/2)
y_basis = de.Fourier('y', Ny, interval=(0, Ly), dealias=3/2)
z_basis = de.Chebyshev('z', Nz, interval=(-Lz, 0), dealias=3/2)
domain = de.Domain([x_basis, y_basis, z_basis], grid_dtype=np.float64)

# Problem
problem = de.IVP(domain, variables=['p','b','u','v','w','bz','uz','vz','wz'])
problem.parameters['N2'] = N2
problem.parameters['Q0'] = Q0
problem.parameters['ν'] = ν
problem.parameters['κ'] = κ
problem.parameters['f'] = f
problem.substitutions['Lap(A,Az)'] = "dx(dx(A)) + dy(dy(A)) + dz(Az)"
problem.substitutions['Adv(A,Az)'] = "u*dx(A) + v*dy(A) + w*Az"
problem.add_equation("dx(u) + dy(v) + wz = 0")
problem.add_equation("dt(b) - κ*Lap(b,bz)               = - Adv(b,bz)")
problem.add_equation("dt(u) - ν*Lap(u,uz) + dx(p) - f*v = - Adv(u,uz)")
problem.add_equation("dt(v) - ν*Lap(v,vz) + dy(p) + f*u = - Adv(v,vz)")
problem.add_equation("dt(w) - ν*Lap(w,wz) + dz(p) - b   = - Adv(w,wz)")
problem.add_equation("bz - dz(b) = 0")
problem.add_equation("uz - dz(u) = 0")
problem.add_equation("vz - dz(v) = 0")
problem.add_equation("wz - dz(w) = 0")
problem.add_bc("left(bz) = N2")
problem.add_bc("left(uz) = 0")
problem.add_bc("left(vz) = 0")
problem.add_bc("left(w) = 0")
problem.add_bc("-κ*right(bz) = Q0")
problem.add_bc("right(uz) = 0")
problem.add_bc("right(vz) = 0")
problem.add_bc("right(w) = 0", condition="(nx != 0) or (ny != 0)")
problem.add_bc("right(p) = 0", condition="(nx == 0) and (ny == 0)")

# Solver
solver = problem.build_solver("RK222")
solver.stop_sim_time = stop_sim_time

# Initial conditions
z = domain.grid(2)
b = solver.state['b']
b['g'] = N2 * z + 1e-3 * N2 * np.random.randn(*b['g'].shape)

# Analysis
out = solver.evaluator.add_file_handler('snapshots', iter=10, max_writes=10)
out.add_task("dy(w) - vz", name='x vorticity')
out.add_task("(u*u + v*v + w*w)/2", name='KE density')

# Main loop
try:
    logger.info('Starting loop')
    start_run_time = time.time()
    while solver.ok:
        solver.step(dt)
        if (solver.iteration-1) % 100 == 0:
            logger.info('Iteration: %i, Time: %e, dt: %e' %(solver.iteration, solver.sim_time, dt))
except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    end_run_time = time.time()
    logger.info('Iterations: %i' %solver.iteration)
    logger.info('Sim end time: %f' %solver.sim_time)
    logger.info('Run time: %f cpu-hr' %((end_run_time-start_run_time)/60/60*domain.dist.comm_cart.size))


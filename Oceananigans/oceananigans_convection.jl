# Free convection in Oceananigans.jl

# Grid setup

using Oceananigans

grid = RegularRectilinearGrid(size = (64, 64, 64),
                              x = (0, 256), y = (0, 256), z = (-128, 0),
                              topology = (Periodic, Periodic, Bounded))

# Boundary conditions and model setup

Q₀ = 1e-8 # Surface buoyancy flux (m² s⁻³)
N² = 1e-5 # Initial and bottom buoyancy gradient (s⁻²)

buoyancy_bcs = TracerBoundaryConditions(grid,
                                        top = BoundaryCondition(Flux, Q₀),
                                        bottom = BoundaryCondition(Gradient, N²))

model = IncompressibleModel(architecture = CPU(),
                            grid = grid,
                            advection = UpwindBiasedFifthOrder(),
                            tracers = :b,
                            buoyancy = BuoyancyTracer(),
                            coriolis = FPlane(f=1e-4),
                            closure = IsotropicDiffusivity(ν=1e-4, κ=1e-4),
                            boundary_conditions = (b=buoyancy_bcs,))

# Set initial condition
noise(z) = 1e-3 * (Q₀ * grid.Δz)^(1/3) * rand()

initial_buoyancy(x, y, z) = N² * z + noise(z)

set!(model, b=initial_buoyancy)

# Simulation setup

using Oceananigans.Units: minute, minutes, hours
using Oceananigans.Utils: prettytime

simulation_time(s) = @info "Time: $(prettytime(s.model.clock.time))"

simulation = Simulation(model,
                        Δt = 1minute,
                        stop_time = 48hours,
                        iteration_interval = 100,
                        progress = simulation_time)

# Output / diagnostics setup

using Oceananigans.AbstractOperations

u, v, w = model.velocities

k = (u^2 + v^2 + w^2) / 2
ξ = ∂y(w) - ∂z(v)

outputs = (kinetic_energy = ComputedField(k),
           x_vorticity = ComputedField(ξ))

simulation.output_writers[:fields] =
    NetCDFOutputWriter(model, outputs,
                       schedule = TimeInterval(30minutes),
                       mode = "c",
                       filepath = "oceananigans_convection.nc")

# Time-step!

run!(simulation)

# ImprovingScientificSoftware2021-HumanCentricAbstractions

Code for the talk "Human-centric Abstractions in Julia and Python for Geophysical Simulation" at the 2021 Improving Scientific Software conference.

```julia
using Oceananigans

grid = RegularRectilinearGrid(size = (256, 256, 256),
                              x = (0, 256), y = (0, 256), z = (-128, 0))

model = IncompressibleModel(grid = grid, 
closure = IsotropicDiffusivity(ν=1e-4, κ=1e-4))

simulation = Simulation(model, Δt=60, stop_time=3600)

run!(simulation)
```

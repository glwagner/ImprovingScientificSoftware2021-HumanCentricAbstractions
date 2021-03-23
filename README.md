# ImprovingScientificSoftware2021-HumanCentricAbstractions

Code for the talk "Human-centric Abstractions in Julia and Python for Geophysical Simulation" at the 2021 Improving Scientific Software conference.

This talk is about user interface design. Here's an example of an [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl) script:

```julia
using Oceananigans

grid = RegularRectilinearGrid(size = (256, 256, 256),
                              x = (-128, 128), y = (-128, 128), z = (-128, 128))

model = IncompressibleModel(grid = grid, 
                            closure = IsotropicDiffusivity(ν=1e-4, κ=1e-4))
                            
thermal_bubble(x, y, z) = 20 + 0.01 * exp(-(x^2 + y^2 + z^2) / 100)
set!(model, T=thermal_bubble)

simulation = Simulation(model, Δt=60, stop_time=3600)

run!(simulation)
```

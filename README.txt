
# Cart-Pole Simulation GUI

This project provides a GUI for simulating and controlling an inverted pendulum on a moving cart using Julia. The simulation supports PID, LQR, and Open Loop control strategies.

## Requirements

- **Julia**: Version 1.10.7
- **Required Packages**:
  - GLMakie
  - Observables
  - DifferentialEquations
  - ControlSystems
  - LinearAlgebra
  - GeometryBasics
- **Makie.jl**: Version v0.15. If you encounter errors with other versions, install version 0.15 specifically:

```julia
# In the Julia REPL:
] rm Makie
] add Makie@0.15.0
```

- Ensure `Project.toml` and `Manifest.toml` are set up for the project environment.

## Setup and Execution

1. Open a terminal or Julia REPL inside the project folder containing `Project.toml` and `Manifest.toml`.
2. Activate and instantiate the project environment:
   ```julia
   julia> ]
   pkg> activate .
   pkg> instantiate
   ```
3. Load and run the main script:
   ```julia
   julia> include("cartpole.jl")
   ```
4. The GUI will appear in a separate window or your current plotting environment.


## Notes

- Ensure you select a control mode and controller type when using Closed Loop.
- Modify parameters as needed before starting a new simulation.
- Use the Reset button to restore default settings.
- Compatibility: The GUI was developed and tested with Makie v0.15. Adjust figure or font sizes in the script for better visualization if needed.

## Conclusion

Experiment with different control parameters to understand the dynamics of an inverted pendulum system. For advanced usage, modify the project environment or explore alternative ODE solvers in the code.


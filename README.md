# Cart-Pole Simulation GUI

This project provides a GUI for simulating and controlling an inverted pendulum on a moving cart using Julia. The simulation supports PID, LQR, and Open Loop control strategies.

---

## Requirements

- **Julia**: Version 1.10.7 ([Download here](https://julialang.org/downloads/))
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

---

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

---

## Features

- **Control Modes**: Select from PID, LQR, or Open Loop strategies.
- **Parameter Customization**: Adjust control and simulation parameters directly from the GUI.
- **Reset Functionality**: Easily restore default settings between runs.
- **Visualization**: Observe real-time dynamics of the cart-pole system.

---

## Detailed Results

For a comprehensive analysis, the simulation generates a detailed report in PDF format. This report includes plots, key metrics, and insights into the system's behavior. You can find the PDF report [here](https://github.com/Melaku-Y/Control-of-Inverted-Pendulum-on-a-Cart/blob/main/Melaku_MESIHU_Final_Report.pdf).

If you want a detailed user manual, please take a look at this PDF [here](https://github.com/Melaku-Y/Control-of-Inverted-Pendulum-on-a-Cart/blob/main/GUI_User_Manual.pdf).



---

## Sample GUI

[here](https://github.com/Melaku-Y/Control-of-Inverted-Pendulum-on-a-Cart/blob/main/GUI.png)

---

## Notes

- Ensure you select a control mode and controller type when using Closed Loop.
- Modify parameters as needed before starting a new simulation.
- Compatibility: The GUI was developed and tested with Makie v0.15. Adjust figure or font sizes in the script for better visualization if needed.

---

## Conclusion

Experiment with different control parameters to understand the dynamics of an inverted pendulum system. For advanced usage, modify the project environment or explore alternative ODE solvers in the code.

---

## Thank You

Thank you for using the Cart-Pole Simulation GUI! We hope you find it helpful and informative.


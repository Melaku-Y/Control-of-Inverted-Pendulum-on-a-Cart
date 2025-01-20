module InvertedPendulumCartPole

using GLMakie
using Observables
using DifferentialEquations
using LinearAlgebra
using ControlSystems 
using GeometryBasics: Point2f0  # Add this import for Point2f0

# Define a two-dimensional point type alias (Float32) for convenience
const Point2f = Point2f0

#######################################################
# 1) DATA STRUCTURES
#######################################################
"""
    Cart

A structure to hold parameters of the Cart-Pole system:
- Physical parameters (masses, lengths, damping, gravity)
"""
struct Cart
    # Physical parameters
    M::Float64         # mass of the cart
    m::Float64         # mass of the pendulum bob
    b::Float64         # damping coefficient
    l::Float64         # pendulum length
    g::Float64         # gravitational acceleration

    # PID gains
    Kp::Float64
    Ki::Float64
    Kd::Float64

    # LQR gain (1x4 vector)
    K_lqr::Vector{Float64}

    # Geometry for visualization
    cart_width::Float64
    cart_height::Float64
    wheel_radius::Float64
    pendulum_length::Float64
    pendulum_mass_radius::Float64
end

"""
    SimParams

Holds the simulation parameters:
- `cart` is a `Cart` instance with system and controller parameters.
- `control_mode` is either "Open Loop" or "Closed Loop".
- `controller_type` is either "PID" or "LQR" (when in Closed Loop).
- `desired_pos` is the target cart position (for closed-loop control).
"""
struct SimParams
    cart::Cart
    control_mode::String
    controller_type::String
    desired_pos::Float64
end

#######################################################
# 2) INIT CART
#######################################################
"""
    init_cart()

Creates and returns a default `Cart` object with some preset physical
parameters, PID gains, and geometry values. The LQR gain is initially
set to a zero vector and will be computed/updated later.
"""
function init_cart()
    empty_Klqr = zeros(Float64, 4)  # placeholder, will be updated by compute_lqr_gain

    return Cart(
        # Physical params
        5.0,    # M (mass of cart)
        1.0,    # m (mass of pendulum bob)
        1.0,    # b (damping coefficient)
        2.5,    # l (pendulum length)
        9.81,   # g (gravitational acceleration)

        # PID gains
        50.0,   # Kp
        5.0,    # Ki
        10.0,   # Kd

        # LQR gain (placeholder)
        empty_Klqr,

        # Geometry
        3.0,    # cart_width
        2.0,    # cart_height
        0.4,    # wheel_radius
        2.5,    # pendulum_length
        0.3     # pendulum_mass_radius
    )
end

#######################################################
# LQR GAIN COMPUTATION
#######################################################
"""
    compute_lqr_gain(cart::Cart)

Computes the LQR gain for the given cart-pole system by constructing
the linearized state-space (A, B) matrices around the inverted position
(θ = π). Returns a 4-element vector (the row vector from LQR in 1x4 form).
"""
function compute_lqr_gain(cart::Cart)
    M = cart.M
    m = cart.m
    d = cart.b
    g = cart.g
    l = cart.l

    # State-space matrices (linearized about θ = π)
    A = [
        0   1           0           0
        0  -d/M       (m*g)/M       0
        0   0           0           1
        0  -d/(M*l)  ((M+m)*g)/(M*l) 0
    ]
    B = [0, 1/M, 0, 1/(M*l)]

    # LQR weighting matrices
    Q = Diagonal([1.0, 1.0, 10.0, 100.0])
    R = 0.0001

    # Compute LQR gain
    K_row = lqr(A, B, Q, R)  # result is 1×4 matrix
    return vec(K_row)        # convert to 4-element Vector
end

"""
    compute_lqr_gain_custom(cart::Cart, qvals::NTuple{4,Float64}, r_val::Float64)

Same as `compute_lqr_gain`, but uses user-specified Q and R values for
the LQR design. The Q diagonal is populated with `qvals` and R is `r_val`.
"""
function compute_lqr_gain_custom(cart::Cart, qvals::NTuple{4,Float64}, r_val::Float64)
    M = cart.M
    m = cart.m
    d = cart.b
    g = cart.g
    l = cart.l

    # State-space matrices
    A = [
        0   1           0           0
        0  -d/M       (m*g)/M       0
        0   0           0           1
        0  -d/(M*l)  ((M+m)*g)/(M*l) 0
    ]
    B = [0, 1/M, 0, 1/(M*l)]

    Q = Diagonal([qvals...])
    R = r_val

    K_row = lqr(A, B, Q, R)
    return vec(K_row)
end

#######################################################
# 3) DYNAMICS
#######################################################
"""
    cartpole_dynamics!(du, u, p::SimParams, t)

Defines the cart-pole dynamics:

States:
  u[1] = x   (cart position)
  u[2] = dx  (cart velocity)
  u[3] = θ   (pendulum angle, measured from the horizontal axis)
  u[4] = dθ  (angular velocity)
  u[5] = e_int (integral of position error, for PID)

Depending on the control mode and controller type in `p`:
- "Open Loop": F = 0
- "Closed Loop":
    - "PID": uses cart.Kp, Ki, Kd
    - "LQR": uses cart.K_lqr
"""
function cartpole_dynamics!(du, u, p::SimParams, t)
    # Unpack state
    x      = u[1]
    dx     = u[2]
    θ      = u[3]
    dθ     = u[4]
    e_int  = u[5]

    # Extract cart and physical parameters
    cart = p.cart
    M, m, b, l, g = cart.M, cart.m, cart.b, cart.l, cart.g

    # -----------------------------
    # 1) DETERMINE CONTROL FORCE F
    # -----------------------------
    if p.control_mode == "Open Loop"
        # Open-loop: no force, no integral error accumulation
        F = 0.0
        du[5] = 0.0
    else
        # Closed-loop control
        if p.controller_type == "PID"
            # Position error
            e = p.desired_pos - x

            # Basic PID (derivative term uses -dx as error derivative)
            F = cart.Kp*e + cart.Ki*e_int + cart.Kd*(-dx)

            # Clamp force to avoid excessive magnitudes
            F = clamp(F, -1000.0, 1000.0)

            # Integral of error for next step
            du[5] = e

        elseif p.controller_type == "LQR"
            # LQR uses full state feedback:
            #   states = [ x - x_ref, dx, θ - π, dθ ]
            x_ref = p.desired_pos
            state_error = [x - x_ref, dx, θ - π, dθ]

            # Compute control action: F = -K_lqr * state_error
            F = -dot(cart.K_lqr, state_error)

            # Clamp force
            F = clamp(F, -1000.0, 1000.0)

            # No integral error to track for LQR
            du[5] = 0.0
        else
            # If somehow neither PID nor LQR is selected, set force to zero
            F = 0.0
            du[5] = 0.0
        end
    end

    # -----------------------------
    # 2) EQUATIONS OF MOTION
    # -----------------------------
    if p.controller_type == "PID"
        # For the PID block, an alternate (slightly more explicit) EoM is used
        # I is a small moment of inertia 
        I = 0.006
        A = M + m
        B = b
        C = m*l*cos(θ)
        D = m*l*sin(θ)*(dθ^2)
        E = I + m*l^2

        denomθ = E - (C^2 / A)
        numθ   = m*g*l*sin(θ) - (C/A)*(F - B*dx + D)
        ddθ    = numθ / denomθ

        ddx    = (F - B*dx + D - C*ddθ) / A

        # Update derivatives
        du[1] = dx
        du[2] = ddx
        du[3] = dθ
        du[4] = ddθ

    else
        # LQR & 'else' case share the same standard cart-pole EoM
        sθ, cθ = sincos(θ)
        denom = M + m*sθ^2

        ddx = (F + m*sθ*(l*(dθ^2) + g*cθ) - b*dx) / denom
        ddθ = (
            -F*cθ
            - m*l*(dθ^2)*cθ*sθ
            - (M + m)*g*sθ
            + b*dx*cθ
        ) / (l * denom)

        # Update derivatives
        du[1] = dx
        du[2] = ddx
        du[3] = dθ
        du[4] = ddθ
    end
end

#######################################################
# 4) SIMULATE CART-POLE
#######################################################
"""
Sets up the initial conditions, ODE problem, and solves for the
cart-pole dynamics. Returns the solution (`sol`).
"""
function simulate_cartpole(x0, θ0, sim_time, control_mode, controller_type, cart::Cart, desired_pos)
    # State vector with an integral error set to 0 initially
    u0 = [x0, 0.0, θ0, 0.0, 0.0]
    tspan = (0.0, sim_time)

    # Store all relevant parameters in a SimParams
    params = SimParams(cart, control_mode, controller_type, desired_pos)

    println("Creating ODE problem...")
    prob = ODEProblem(cartpole_dynamics!, u0, tspan, params)
    println("Solving ODE problem...")
    sol = solve(prob, Tsit5())
    println("Solution obtained!")
    return sol
end

#######################################################
# 5) UPDATE VISUALIZATION
#######################################################
"""
    update_visualization!(ax, x, θ, cart::Cart)

Clears and re-draws the entire cart-pole system within the given
Makie Axis `ax`. The cart is drawn at position x, and the pendulum at angle θ.
"""
function update_visualization!(ax, x, θ, cart::Cart)
    # Clear existing drawings from the axis
    empty!(ax)

    # Unpack geometry
    cw = cart.cart_width
    ch = cart.cart_height
    wr = cart.wheel_radius
    al = cart.pendulum_length
    mr = cart.pendulum_mass_radius

    # Cart center vertical coordinate (so it sits above the wheels)
    y = wr + ch/2

    # Set axis limits and style
    limits!(ax, -20, 20, 0, 6.5)
    ax.aspect = DataAspect()
    hideydecorations!(ax)
    ax.xticks = -12:2:12
    ax.xgridvisible = false
    ax.ygridvisible = false
    hidespines!(ax)

    # Ground line
    lines!(ax, [-12, 12], [0, 0], color=:black, linewidth=3)

    # Cart body polygon
    cart_x = [x-cw/2, x+cw/2, x+cw/2, x-cw/2, x-cw/2]
    cart_y = [y-ch/2, y-ch/2, y+ch/2, y+ch/2, y-ch/2]
    poly!(ax, Point2f.(zip(cart_x, cart_y)), color=:blue)

    # Wheels (simple circles)
    t = range(0, 2π, length=50)

    wheel1_x = x - cw/2 .+ wr*cos.(t) .+ wr
    wheel1_y = wr*sin.(t) .+ wr
    wheel2_x = x + cw/2 - wr .+ wr*cos.(t)
    wheel2_y = wr*sin.(t) .+ wr
    poly!(ax, Point2f.(zip(wheel1_x, wheel1_y)), color=:black)
    poly!(ax, Point2f.(zip(wheel2_x, wheel2_y)), color=:black)

    # Pendulum rod (line from cart center to pendulum bob)
    lines!(ax, [x, x + al*sin(θ)], [y, y - al*cos(θ)], color=:brown, linewidth=3)

    # Pendulum mass (small circle at the rod tip)
    mass_x = x + al*sin(θ) .+ mr/2*cos.(t)
    mass_y = y - al*cos(θ) .+ mr/2*sin.(t)
    poly!(ax, Point2f.(zip(mass_x, mass_y)), color=:gray)
end

#######################################################
# 6) CREATE PARAM CONTROL (SLIDERS/TEXTBOX)
#######################################################
"""
Utility function to create a label, slider, textbox, and numeric label
in one row of a GridLayout. The slider and textbox both update the
same `observable`, so changes in one reflect in the other.

Returns `(slider, textbox)` so they can be further managed if needed.
"""
function create_parameter_control(layout, row, label_str, range, default, observable)
    # Label in first column
    Label(layout[row, 1], label_str)

    # Slider in second column
    slider = Slider(layout[row, 2], width=150, range=range, startvalue=default)

    # Textbox in third column
    textbox = Textbox(
        layout[row, 3],
        bordercolor_hover="Black",
        width=60,
        boxcolor="white",
        placeholder=string(default)
    )

    # Value label in fourth column (displays the slider value)
    value_label = Label(layout[row, 4], @lift(string(round($(slider.value), digits=3))))

    # Connect the slider’s value to the observable
    connect!(observable, slider.value)

    # When the textbox changes, parse as Float64 and try to set slider
    on(textbox.stored_string) do s
        try
            val = parse(Float64, s)
            # Only update slider if value is within slider range
            if val >= minimum(range) && val <= maximum(range)
                set_close_to!(slider, val)
            end
        catch
            # If parse fails or value is out of range, do nothing
        end
    end

    return slider, textbox
end

#######################################################
# 7) BUILD THE MAKIE GUI
#######################################################
"""
Builds and returns a GUI (Figure) containing:
Returns `(fig, sim_state)`.
"""
function create_makie_gui()
    # Main figure with a 2x1 grid layout
    fig = Figure(
        size=(1400, 800),
        backgroundcolor=:lightgray,
        font="Comic Sans",
        fontsize=14,
        padding=(10, 10, 10, 10)   # top, right, bottom, left
    )

    # Left control panel
    control_panel = fig[1, 1] = GridLayout(padding=(10, 10, 10, 10))

    # Right visualization panel
    viz_panel = fig[1, 2] = GridLayout(padding=(10, 10, 10, 10))

    # Error label for warnings or missing selections
    error_label = Label(control_panel[2, 3], "", color=:red, halign=:left)

    # Menus for Control Mode & Controller Type
    Label(control_panel[2, 1], "Control Mode:")
    control_mode_menu = Menu(
        control_panel[2, 2],
        width=110, height=15,
        cell_color_hover="Gray",
        options=["Select", "Open Loop", "Closed Loop"]
    )

    Label(control_panel[3, 1], "Controller Type:")
    controller_type_menu = Menu(
        control_panel[3, 2],
        width=110, cell_color_hover="Gray",
        options=["Select", "PID", "LQR"]
    )

    # Observables that hold parameter values
    current_pos_obs   = Observable(0.0)
    desired_pos_obs   = Observable(0.0)
    sim_time_obs      = Observable(30.0)
    initial_angle_obs = Observable(3.14)

    # Sliders & textboxes for starting position, desired position, sim time, etc.
    (curpos_slider, curpos_text) =
        create_parameter_control(control_panel, 4,
                                 "Starting Position:",
                                 -11.0:0.1:11.0, 0.0,
                                 current_pos_obs)

    (despos_slider, despos_text) =
        create_parameter_control(control_panel, 5,
                                 "Desired Position:",
                                 -11.0:0.1:11.0, 0.0,
                                 desired_pos_obs)

    (time_slider, time_text) =
        create_parameter_control(control_panel, 6,
                                 "Simulation Time:",
                                 1.0:1.0:100.0, 30.0,
                                 sim_time_obs)

    (initang_slider, initang_text) =
        create_parameter_control(control_panel, 7,
                                 "Initial Angle:",
                                 0.0:0.1:6.28, 3.14,
                                 initial_angle_obs)

    # PID Gains
    Kp_obs = Observable(50.0)
    Ki_obs = Observable(5.0)
    Kd_obs = Observable(10.0)

    (Kp_slider, Kp_text) = create_parameter_control(
        control_panel, 8, "Kp:", 0.0:1.0:300.0, 50.0, Kp_obs)
    (Ki_slider, Ki_text) = create_parameter_control(
        control_panel, 9, "Ki:", 0.0:0.5:150.0, 5.0, Ki_obs)
    (Kd_slider, Kd_text) = create_parameter_control(
        control_panel, 10, "Kd:", 0.0:1.0:300.0, 10.0, Kd_obs)

    # LQR Weights
    Q1_obs = Observable(1.0)
    Q2_obs = Observable(1.0)
    Q3_obs = Observable(10.0)
    Q4_obs = Observable(100.0)
    R_obs  = Observable(0.0001)

    (Q1_slider, Q1_text) = create_parameter_control(
        control_panel, 11, "Q1:", 0.0:1.0:200.0, 1.0, Q1_obs)
    (Q2_slider, Q2_text) = create_parameter_control(
        control_panel, 12, "Q2:", 0.0:1.0:200.0, 1.0, Q2_obs)
    (Q3_slider, Q3_text) = create_parameter_control(
        control_panel, 13, "Q3:", 0.0:1.0:200.0, 10.0, Q3_obs)
    (Q4_slider, Q4_text) = create_parameter_control(
        control_panel, 14, "Q4:", 0.0:1.0:500.0, 100.0, Q4_obs)
    (R_slider, R_text) = create_parameter_control(
        control_panel, 15, "R:", 1e-6:1e-5:1e-2, 0.0001, R_obs)

    # Simulation control buttons
    start_button = Button(
        control_panel[16, 1],
        cornerradius=15, buttoncolor_hover="Gray",
        label="Start Simulation"
    )
    reset_button = Button(
        control_panel[17, 1],
        cornerradius=15, buttoncolor_hover="Gray",
        label="Reset"
    )

    # --------------------------------------------------
    # Visualization Panel (Title + Axes + Plots)
    # --------------------------------------------------
    desc = """
      Inverted Pendulum Cart-Pole Control in Julia(PID & LQR)
    """
    Label(
        viz_panel[1, 1:2], desc, textsize=28, fontweight=:Bold, halign=:left,
        width=500, padding=(10, 10, 10, 10)
    )

    # Main animation axis
    animation_ax = Axis(
        viz_panel[2, 1:2],
        aspect=DataAspect(),
        backgroundcolor=:lightgray,
        titleweight=:bold,
        title="Cart-Pole Animation",
        titlesize=26,
        xzoomlock=true, yzoomlock=true,
        xpanlock=true, ypanlock=true,
        margin=10px
    )

    # Axes for plotting position/velocity and angle/angular velocity
    pos_vel_ax = Axis(
        viz_panel[3, 1],
        title="Position & Velocity",
        xlabel="Time (s)",
        ylabel="m, m/s",
        margin=10px
    )
    ang_vel_ax = Axis(
        viz_panel[3, 2],
        title="Angle & Angular Velocity",
        xlabel="Time (s)",
        ylabel="rad, rad/s",
        margin=10px
    )

    # Initialize a base Cart for the initial drawing
    base_cart = init_cart()
    update_visualization!(animation_ax, 0.0, π, base_cart)

    # Observables for plotting solution data
    t_data  = Observable(Float64[])
    x_data  = Observable(Float64[])
    dx_data = Observable(Float64[])
    θ_data  = Observable(Float64[])
    dθ_data = Observable(Float64[])

    # Plot lines for the timeseries
    lines!(pos_vel_ax, t_data, x_data,  color=:red,   label="Position (m)")
    lines!(pos_vel_ax, t_data, dx_data, color=:black, label="Velocity (m/s)")
    axislegend(pos_vel_ax)

    lines!(ang_vel_ax, t_data, θ_data,  color=:blue,   label="Angle (rad)")
    lines!(ang_vel_ax, t_data, dθ_data, color=:orange, label="AngVel (rad/s)")
    axislegend(ang_vel_ax)

    # Dictionary to keep track of simulation state
    sim_state = Dict(
        "running"         => false,
        "control_mode"    => "Select",
        "controller_type" => "Select"
    )

    # Update control mode and controller type when menus change
    on(control_mode_menu.selection) do sel
        sim_state["control_mode"] = sel
    end
    on(controller_type_menu.selection) do sel
        sim_state["controller_type"] = sel
    end

    #####################################
    # START SIMULATION BUTTON
    #####################################
    on(start_button.clicks) do _
        # Check user selections
        if sim_state["control_mode"] == "Select"
            error_label.text = "Please select a Control Mode!"
            return
        end
        if sim_state["control_mode"] == "Closed Loop"
            if sim_state["controller_type"] == "Select"
                error_label.text = "Please select a Controller Type (PID or LQR)!"
                return
            end
        end

        # Clear any previous error
        error_label.text = ""

        # Start simulation only if not already running
        if !sim_state["running"]
            try
                x0     = current_pos_obs[]
                desired= desired_pos_obs[]
                time   = sim_time_obs[]
                angle0 = initial_angle_obs[]

                # User-specified gains
                user_Kp = Kp_obs[]
                user_Ki = Ki_obs[]
                user_Kd = Kd_obs[]

                user_Q1 = Q1_obs[]
                user_Q2 = Q2_obs[]
                user_Q3 = Q3_obs[]
                user_Q4 = Q4_obs[]
                user_R  = R_obs[]

                # Compute new LQR gain using user-specified Q, R
                new_cart_lqr = compute_lqr_gain_custom(
                    base_cart,
                    (user_Q1, user_Q2, user_Q3, user_Q4),
                    user_R
                )

                # Create a new cart object with updated PID gains and LQR gain
                cart_updated = Cart(
                    base_cart.M, base_cart.m, base_cart.b, base_cart.l, base_cart.g,
                    user_Kp, user_Ki, user_Kd,
                    new_cart_lqr,
                    base_cart.cart_width, base_cart.cart_height,
                    base_cart.wheel_radius, base_cart.pendulum_length,
                    base_cart.pendulum_mass_radius
                )

                # Flag simulation as running
                sim_state["running"] = true
                start_button.label   = "Running..."

                # Solve the ODE
                sol = simulate_cartpole(
                    x0, angle0, time,
                    sim_state["control_mode"],
                    sim_state["controller_type"],
                    cart_updated,
                    desired
                )

                # Clear old data from plots
                empty!(t_data[])
                empty!(x_data[])
                empty!(dx_data[])
                empty!(θ_data[])
                empty!(dθ_data[])
                notify.([t_data, x_data, dx_data, θ_data, dθ_data])

                # Animation loop: sample solution 300 times over [0, time]
                @async begin
                    for tt in range(0, time, length=300)
                        if !sim_state["running"]
                            break
                        end

                        st = sol(tt)
                        xx, dxx, th, dth = st[1], st[2], st[3], st[4]

                        # Update the animation
                        update_visualization!(animation_ax, xx, th, cart_updated)

                        # Append data for plotting
                        push!(t_data[], tt)
                        push!(x_data[], xx)
                        push!(dx_data[], dxx)
                        push!(θ_data[], th)
                        push!(dθ_data[], dth)
                        notify.([t_data, x_data, dx_data, θ_data, dθ_data])

                        sleep(0.03)  # small pause for animation
                    end

                    # After the loop, mark simulation as done
                    sim_state["running"] = false
                    start_button.label   = "Start Simulation"
                end

            catch e
                # In case of any runtime error
                @warn "Simulation error: $e"
                sim_state["running"] = false
                start_button.label   = "Start Simulation"
            end
        end
    end

    #####################################
    # RESET BUTTON
    #####################################
    on(reset_button.clicks) do _
        # Stop any running simulation
        sim_state["running"] = false
        start_button.label   = "Start Simulation"

        # Reset Observables
        current_pos_obs[]   = 0.0
        desired_pos_obs[]   = 0.0
        sim_time_obs[]      = 30.0
        initial_angle_obs[] = 3.14

        Kp_obs[] = 50.0
        Ki_obs[] = 5.0
        Kd_obs[] = 10.0

        Q1_obs[] = 1.0
        Q2_obs[] = 1.0
        Q3_obs[] = 10.0
        Q4_obs[] = 100.0
        R_obs[]  = 0.0001

        control_mode_menu.selection    = "Select"
        sim_state["control_mode"]      = "Select"
        controller_type_menu.selection = "Select"
        sim_state["controller_type"]   = "Select"

        # Clear error label
        error_label.text = ""

        # Clear plot data
        empty!(t_data[])
        empty!(x_data[])
        empty!(dx_data[])
        empty!(θ_data[])
        empty!(dθ_data[])
        notify.([t_data, x_data, dx_data, θ_data, dθ_data])

        # Reset the main animation to default
        update_visualization!(animation_ax, 0.0, π, base_cart)

        # Reset all slider & textbox widgets
        set_close_to!(curpos_slider, 0.0)
        curpos_text.text[] = "0.0"

        set_close_to!(despos_slider, 0.0)
        despos_text.text[] = "0.0"

        set_close_to!(time_slider, 30.0)
        time_text.text[] = "30.0"

        set_close_to!(initang_slider, 3.14)
        initang_text.text[] = "3.14"

        set_close_to!(Kp_slider, 50.0)
        Kp_text.text[] = "50.0"

        set_close_to!(Ki_slider, 5.0)
        Ki_text.text[] = "5.0"

        set_close_to!(Kd_slider, 10.0)
        Kd_text.text[] = "10.0"

        set_close_to!(Q1_slider, 1.0)
        Q1_text.text[] = "1.0"

        set_close_to!(Q2_slider, 1.0)
        Q2_text.text[] = "1.0"

        set_close_to!(Q3_slider, 10.0)
        Q3_text.text[] = "10.0"

        set_close_to!(Q4_slider, 100.0)
        Q4_text.text[] = "100.0"

        set_close_to!(R_slider, 0.0001)
        R_text.text[] = "0.0001"
    end

    # Finally, display the figure
    display(fig)
    return fig, sim_state
end

#######################################################
# 8) MAIN
#######################################################
"""
    main()

Creates the Makie-based GUI by calling `create_makie_gui()` and returns
the figure and simulation state dictionary.
"""
function main()
    fig, sim_state = create_makie_gui()
    return fig, sim_state
end

# Run the main function to launch the GUI
fig, sim_state = main()
end  # module

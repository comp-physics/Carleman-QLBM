QCFD_SRC = ENV["QCFD_SRC"]  
include(QCFD_SRC * "LBM/cal_feq.jl")

function f_ini_feq(u0, rho0)
    _, _, w_value, e_value, _, _, _, _, a_v, b_v, c_v, d_v = lbm_const_sym()
    vx = u0; vy = 0.; wm = w_value; vcx = e_value; vcy = [0. , 0. , 0.]
    f_ini = cal_feq(rho0, vx, vy, wm, vcx, vcy, a_v, b_v, c_v, d_v)
    println("sum f_m =", sum(f_ini))
    return f_ini
end

function f_ini_test(uu)
    _, _, w_value, _, _, _, _, _, a_v, b_v, c_v, d_v = lbm_const_sym()
    f_ini = w_value .+ [-uu/2, 0, uu/2]
    return f_ini 
end

"""
    d1q3_state_from_velocity(u)

Build the single-site D1Q3 state used by the current repository from a
velocity-like parameter `u`. This is the per-site building block used for both
the original/legacy multigrid initial condition and the sinusoidal option.
"""
function d1q3_state_from_velocity(u)
    return f_ini_test(u)
end

raw"""
    d1q3_initial_condition_from_velocity_profile(velocity_profile)

Concatenate sitewise D1Q3 states into the interleaved multigrid state

```math
\phi = (f(x_1), f(x_2), \ldots, f(x_{N_x}))^\top.
```
"""
function d1q3_initial_condition_from_velocity_profile(velocity_profile)
    return reduce(vcat, (d1q3_state_from_velocity(u) for u in velocity_profile))
end

raw"""
    d1q3_sinusoidal_velocity_profile(nx; u_ini=0.1)

Return the paper-style sinusoidal D1Q3 velocity profile

```math
u_i = u_{\mathrm{ini}} \sin\left(\frac{2\pi i}{N_x}\right), \qquad i = 1,\dots,N_x.
```

In the notation of `LBE_theory.tex`, this corresponds to the spatial part of

```math
\mathbf{g}(r_x^\star,0) = \mathbf{g}^{\mathrm{eq}}\!\left(
u^\star = u^\star_{\mathrm{ini}} \sin\frac{2\pi r_x^\star}{N_x}
\right).
```
"""
function d1q3_sinusoidal_velocity_profile(nx; u_ini=0.1)
    return [u_ini * sin(2 * pi * rx / nx) for rx in 1:nx]
end

"""
    d1q3_legacy_velocity_profile(nx)

Return the original repository D1Q3 multigrid velocity profile.

- For `nx == 3`, this reproduces the historical regression-test profile
  `[0.12, 0.00, -0.08]`.
- Otherwise, it uses the original linear ramp from `0.12` to `-0.08` across the
  grid.
"""
function d1q3_legacy_velocity_profile(nx)
    if nx == 1
        return [0.12]
    end
    if nx == 3
        return [0.12, 0.00, -0.08]
    end
    return collect(range(0.12, -0.08, length=nx))
end

raw"""
    d1q3_multigrid_initial_condition(nx; initial_condition=:legacy, u_ini=0.1)

Select the D1Q3 multigrid initial condition.

Supported options:

- `initial_condition = :legacy` : original repository behavior
- `initial_condition = :sinusoidal` : sinusoidal profile

For the sinusoidal option, `u_ini` sets the amplitude in

```math
u_i = u_{\mathrm{ini}} \sin\left(\frac{2\pi i}{N_x}\right).
```
"""
function d1q3_multigrid_initial_condition(nx; initial_condition=:legacy, u_ini=0.1)
    velocity_profile = if initial_condition == :legacy
        d1q3_legacy_velocity_profile(nx)
    elseif initial_condition == :sinusoidal
        d1q3_sinusoidal_velocity_profile(nx; u_ini=u_ini)
    else
        error("Unsupported D1Q3 initial_condition=$(initial_condition). Supported values are :legacy and :sinusoidal.")
    end

    return d1q3_initial_condition_from_velocity_profile(velocity_profile)
end

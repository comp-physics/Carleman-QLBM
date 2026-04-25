# carleman_burgers.jl

using LinearAlgebra
using SparseArrays
using DifferentialEquations
using Interpolations
using Plots
using Printf
using LaTeXStrings

# -------------------------
# kronp.m from https://github.com/JuliaReach/CarlemanLinearization.jl -> Julia 
# -------------------------
function kronp(A, k::Int)
B = 1
for _ = 1:k
    B = kron(B, A)
end
return B
end

function main()
    # %% Simulation parameters [1]
    nx = 16
    nt = 4000
    nx_pde = 100
    nt_pde = 40000

    Re0 = 20.0
    L0  = 1.0
    U0  = 1 / sqrt(nx - 1)
    beta = 0.0
    f = 1.0
    T = 3.0

    F0_fun(t, x) = U0 .* exp.(-((x .- L0/4).^2) ./ (2*(L0/32)^2)) .* cos.(2*pi*t)

    N_max = 4
    ode_deg = 2

    Ns = collect(1:N_max)
    nu = U0 * L0 / Re0
    Tnl = L0 / U0
    t_plot = Tnl / 3

    x0, x1 = -L0/2, L0/2
    t0, t1 = 0.0, T

    dx = (x1 - x0) / (nx - 1)
    dt = (t1 - t0) / (nt - 1)
    xs = collect(range(x0, x1; length=nx))
    ts = collect(range(t0, t1; length=nt))

    nt_ode = nt * 10
    dt_ode = (t1 - t0) / (nt_ode - 1)
    ts_ode = collect(range(t0, t1; length=nt_ode))

    dx_pde = (x1 - x0) / (nx_pde - 1)
    dt_pde = (t1 - t0) / (nt_pde - 1)
    xs_pde = collect(range(x0, x1; length=nx_pde))
    ts_pde = collect(range(t0, t1; length=nt_pde))

    F0 = zeros(nt, nx)
    for it in 1:nt
        F0[it, :] .= F0_fun(ts[it], xs)
    end

    F1 = zeros(nx, nx)
    for i in 1:nx
        F1[i,i] = -2*nu/dx^2
    end
    for i in 1:nx-1
        F1[i, i+1] = nu/dx^2
        F1[i+1, i] = nu/dx^2
    end
    F1 .-= beta * I(nx)

    F2 = zeros(nx, nx^2)
    step = nx^2 + nx + 1
    vF2 = vec(F2)
    for idx in step:step:length(vF2)
        vF2[idx] = -1/(4*dx)
    end
    for idx in 2:step:length(vF2)
        vF2[idx] = +1/(4*dx)
    end
    F2 = reshape(vF2, nx, nx^2)

    F1[1, :] .= 0
    F1[end, :] .= 0
    F2[1, :] .= 0
    F2[end, :] .= 0

    u0(x) = -U0*sin(2*pi*f*x/L0)
    u0s = u0.(xs)

    F0_itps = [linear_interpolation(ts, F0[:,j]; extrapolation_bc=Line()) for j in 1:nx]
    F0_interp(t) = [F0_itps[j](t) for j in 1:nx]

    function burgers_odefun!(du, u, p, t)
        du .= F0_interp(t) .+ F1*u .+ F2*(kron(u,u))
        return nothing
    end

    function pde_rhs!(du, u, p, t)
        # enforce Dirichlet
        uL = 0.0
        uR = 0.0

        dxl = xs_pde[2] - xs_pde[1]

        dudx = zeros(length(u))
        for i in 2:length(u)-1
            dudx[i] = (u[i+1] - u[i-1])/(2*dxl)
        end

        flux = nu .* dudx .- 0.5 .* (u.^2)
        flux_x = zeros(length(u))
        for i in 2:length(u)-1
            flux_x[i] = (flux[i+1] - flux[i-1])/(2*dxl)
        end

        src = -beta .* u .+ F0_fun(t, xs_pde)
        du .= flux_x .+ src

        du[1] = 0.0
        du[end] = 0.0
        return nothing
    end

    C1_e = U0*dt/dx
    C2_e = 2*nu*dt/dx^2
    C1_ode = U0*dt_ode/dx
    C2_ode = 2*nu*dt_ode/dx^2
    C1_pde = U0*dt_pde/dx_pde
    C2_pde = 2*nu*dt_pde/dx_pde^2
    if C1_e > 1; error(@sprintf("C1_e = %.2f\n", C1_e)); end
    if C2_e > 1; error(@sprintf("C2_e = %.2f\n", C2_e)); end
    if C1_ode > 1; error(@sprintf("C1_ode = %.2f\n", C1_ode)); end
    if C2_ode > 1; error(@sprintf("C2_ode = %.2f\n", C2_ode)); end
    if C1_pde > 1; error(@sprintf("C1_pde = %.2f\n", C1_pde)); end
    if C2_pde > 1; error(@sprintf("C2_pde = %.2f\n", C2_pde)); end

    lambdas = eigvals(F1)
    lambdas = filter(!=(0.0), lambdas)
    lam = maximum(lambdas)

    f2n = opnorm(F2)
    f1n = opnorm(F1)
    f0n = 0.0
    for it in 1:nt
        f0n = max(norm(F0[it, :]), f0n)
    end
    R = (norm(u0s)*f2n + f0n/norm(u0s)) / abs(lam)

    if dt > 1/(N_max*f1n)
        error("Time step too large")
    end
    if f0n + f2n > abs(lam)
        @info "Perturbation too large"
    end

    println("Preparing Carleman matrix")
    dNs = zeros(Int, N_max)
    for N in Ns
        dNs[N] = Int((nx^(N+1) - nx) ÷ (nx - 1))
    end

    A = spzeros(Float64, dNs[end], dNs[end])

    Fs = hcat(reshape(F0_fun(1.0, xs), nx, 1), F1, F2)

    Inx = sparse(I(nx))

    for i in Ns
        for j in 0:min(ode_deg, N_max - i + 1)
            if i == 1 && j == 0
                continue
            end

            a0 = 1 + (nx^i - nx) ÷ (nx - 1)
            a1 = a0 + nx^i - 1
            b0 = 1 + (nx^(j+i-1) - nx) ÷ (nx - 1)
            b1 = b0 + nx^(j+i-1) - 1

            Aij = spzeros(Float64, nx^i, nx^(i+j-1))

            f0 = 1 + (nx^j - nx) ÷ (nx - 1) + 1
            f1 = f0 + nx^j - 1
            Fj = Fs[:, f0:f1]  # [1]

            for p in 1:i
                Ia = kronp(Inx, p-1)
                Ib = kronp(Inx, i-p)
                Aij .= Aij .+ kron(kron(Ia, sparse(Fj)), Ib)  # [1]
            end

            A[a0:a1, b0:b1] = Aij  # [1]
        end
    end

    ys_c_N = zeros(Float64, N_max, nt, dNs[end])
    for N in Ns
        dimN = dNs[N]
        A_N = A[1:dimN, 1:dimN]
        b_N = zeros(Float64, dimN)
        b_N[1:nx] .= F0_fun(1.0, xs)

        y0s = Float64[]
        for i in 1:N
            y0s = vcat(y0s, vec(kronp(u0s, i)))
        end

        @printf("Solving Carleman N=%d\n", N)
        ys = zeros(Float64, nt, dimN)
        ys[1, :] .= y0s

        for k in 1:(nt-1)
            # Rebuild the inhomogeneous part per time step [1]
            A_N_mod = copy(A_N)

            for i in 2:N
                a0 = 1 + (nx^i - nx) ÷ (nx - 1)
                a1 = a0 + nx^i - 1
                b0 = 1 + (nx^(i-1) - nx) ÷ (nx - 1)
                b1 = b0 + nx^(i-1) - 1

                Aij = spzeros(Float64, nx^i, nx^(i-1))
                Fj = reshape(F0_fun(ts[k], xs), nx, 1)  # [1]

                for p in 1:i
                    Ia = kronp(Inx, p-1)
                    Ib = kronp(Inx, i-p)
                    Aij .= Aij .+ kron(kron(Ia, sparse(Fj)), Ib)  # [1]
                end
                A_N_mod[a0:a1, b0:b1] = Aij
            end

            b_N[1:nx] .= F0_fun(ts[k], xs)
            ys[k+1, :] .= ys[k, :] .+ dt .* (A_N_mod*ys[k, :] .+ b_N)  # [1]
        end

        println("Done")
        ys_c_N[N, :, 1:dimN] .= real.(ys)
    end
    us_c_N = ys_c_N[:, :, 1:nx]  # [1]

    println("Solving direct Euler")
    us_e = zeros(Float64, nt, nx)
    us_e[1, :] .= u0s
    tmp = zeros(nx)
    for k in 1:(nt-1)
        burgers_odefun!(tmp, view(us_e, k, :), nothing, ts[k])
        us_e[k+1, :] .= us_e[k, :] .+ dt .* tmp
    end

    println("Solving \"exact\" ODE")
    prob_ode = ODEProblem(burgers_odefun!, u0s, (t0, t1))
    sol_ode = solve(prob_ode, Tsit5(); reltol=1e-10, abstol=1e-10, saveat=ts_ode)
    us_ode = reduce(hcat, sol_ode.u)'  # (nt_ode, nx)

    us_d = zeros(nt, nx)
    for j in 1:nx
        itp = linear_interpolation(ts_ode, us_ode[:, j]; extrapolation_bc=Line())
        us_d[:, j] .= itp.(ts)
    end

    
    println("Solving \"exact\" PDE")
    u0_pde = u0.(xs_pde)
    u0_pde[1] = 0.0
    u0_pde[end] = 0.0
    prob_pde = ODEProblem(pde_rhs!, u0_pde, (t0, t1))
    sol_pde = solve(prob_pde, Tsit5(); reltol=1e-6, abstol=1e-8, saveat=ts_pde)
    us_pde = reduce(hcat, sol_pde.u)'  # (nt_pde, nx_pde)

    us_pde_interp_temp = zeros(nt, nx_pde)
    for i in 1:nx_pde
        itp_t = linear_interpolation(ts_pde, us_pde[:, i]; extrapolation_bc=Line())
        us_pde_interp_temp[:, i] .= itp_t.(ts)
    end
    us_pde_interp = zeros(nt, nx)
    for k in 1:nt
        itp_x = linear_interpolation(xs_pde, us_pde_interp_temp[k, :]; extrapolation_bc=Line())
        us_pde_interp[k, :] .= itp_x.(xs)
    end

    # %% Calculate errors [1]
    dus_c_d_N = zeros(Float64, N_max, nt, nx)
    dus_rel_c_d_N = zeros(Float64, N_max, nt, nx)
    eps_c_d_N = zeros(Float64, N_max, nt)
    eps_rel_c_d_N = zeros(Float64, N_max, nt)

    dus_c_pde_N = zeros(Float64, N_max, nt, nx)
    dus_rel_c_pde_N = zeros(Float64, N_max, nt, nx)
    eps_c_pde_N = zeros(Float64, N_max, nt)
    eps_rel_c_pde_N = zeros(Float64, N_max, nt)

    dus_d_pde = zeros(Float64, nt, nx)
    dus_rel_d_pde = zeros(Float64, nt, nx)
    eps_d_pde = zeros(Float64, nt)
    eps_rel_d_pde = zeros(Float64, nt)

    dus_d_e = zeros(Float64, nt, nx)
    dus_rel_d_e = zeros(Float64, nt, nx)
    eps_d_e = zeros(Float64, nt)
    eps_rel_d_e = zeros(Float64, nt)

    for N in 1:N_max
        dus_c_d_N[N, :, :] .= us_c_N[N, :, :] .- us_d
        dus_rel_c_d = dus_c_d_N[N, :, :] ./ us_d
        dus_rel_c_d[isnan.(dus_rel_c_d)] .= 0.0
        dus_rel_c_d_N[N, :, :] .= dus_rel_c_d

        dus_c_pde_N[N, :, :] .= us_c_N[N, :, :] .- us_pde_interp
        dus_rel_c_pde = dus_c_pde_N[N, :, :] ./ us_pde_interp
        dus_rel_c_pde[isnan.(dus_rel_c_pde)] .= 0.0
        dus_rel_c_pde[isinf.(dus_rel_c_pde)] .= 0.0
        dus_rel_c_pde_N[N, :, :] .= dus_rel_c_pde

        dus_d_pde[:, :] .= us_d .- us_pde_interp
        dus_rel_d_pde[:, :] .= dus_d_pde ./ us_pde_interp
        dus_rel_d_pde[isnan.(dus_rel_d_pde)] .= 0.0
        dus_rel_d_pde[isinf.(dus_rel_d_pde)] .= 0.0

        dus_d_e[:, :] .= us_d .- us_e
        dus_rel_d_e[:, :] .= dus_d_e ./ us_e
        dus_rel_d_e[isnan.(dus_rel_d_e)] .= 0.0
        dus_rel_d_e[isinf.(dus_rel_d_e)] .= 0.0

        for k in 1:nt
            eps_c_d_N[N, k] = norm(vec(dus_c_d_N[N, k, :]))
            eps_rel_c_d_N[N, k] = norm(vec(dus_rel_c_d_N[N, k, :]), Inf)

            eps_c_pde_N[N, k] = norm(vec(dus_c_pde_N[N, k, :]))
            eps_rel_c_pde_N[N, k] = norm(vec(dus_rel_c_pde_N[N, k, :]), Inf)

            eps_d_pde[k] = norm(vec(dus_d_pde[k, :]))
            eps_rel_d_pde[k] = norm(vec(dus_rel_d_pde[k, :]), Inf)

            eps_d_e[k] = norm(vec(dus_d_e[k, :]))
            eps_rel_d_e[k] = norm(vec(dus_rel_d_e[k, :]), Inf)
        end
    end

    i_plot = findfirst(>=(t_plot), ts)
    i_start = Int(ceil(i_plot*3/4))

    p1 = plot(xs_pde, us_pde[1, :]; linestyle=:dash, color=:black, label="Initial condition")
    plot!(p1, xs_pde, F0_fun(1.0, xs_pde); linestyle=:dashdot, color=:black, label="Source shape")
    plot!(p1, xs, us_d[i_plot, :]; marker=:circle, color=:black, label=L"Direct Euler solution at $T_{nl}/3$")

    for N in (Ns[1], Ns[end])
        plot!(p1, xs, vec(ys_c_N[N, i_plot, 1:nx]); marker=:star5, label=latexstring("Carleman solution at \$T_{nl}/3\$, \$N=$N\$"))
    end
    ylims!(p1, (-maximum(abs.(us_pde[1, :])), maximum(abs.(us_pde[1, :]))))
    xlabel!(p1, L"x")
    ylabel!(p1, L"u")
    xlims!(p1, (x0, x1))

    p2 = plot(; yscale=:log10, title="Absolute error", xlabel=L"t", ylabel=L"\|\varepsilon_{\mathrm{abs}}\|_2")
    for N in Ns
        plot!(p2, ts, eps_c_d_N[N, :], label=latexstring("Carleman, \$N=$N\$"))
    end
    vline!(p2, [t_plot]; linestyle=:dot, color=:black, label="")
    ymin = min(minimum(eps_c_d_N[:, i_start:end]), minimum(eps_d_pde[i_start:end])) * 0.1
    ymax = maximum(eps_c_d_N[1, :]) * 10
    ylims!(p2, (ymin, ymax))

    p3 = plot(Ns, vec(maximum(eps_c_d_N; dims=2)); yscale=:log10, marker=:circle,
                title="Error convergence", xlabel=L"N", ylabel=latexstring("\\max_t \\|\\varepsilon_{\\mathrm{abs}}\\|_2"),
                label="Time-maximum error")

    Re_act = maximum(us_pde) * L0 / nu
    bigtitle = @sprintf("Forced VBE solution with Re=%.2f, n_x=%d, n_t=%d, R=%.2f", Re_act, nx, nt, R)

    plt = plot(p1, p2, p3; layout=@layout([a{0.55h}; b c]), size=(1000, 700), plot_title=bigtitle)
    display(plt)

    
    outname = @sprintf("vbe_re0_%.2f_N_%d_nx_%d_nt_%d_rev2.png", Re0, N_max, nx, nt)
    savefig(plt, outname)

        return nothing
    end

main()

# Running the production D2Q9 CLBE

```bash
cd "$HOME/Carleman-QLBM"
julia --project=. -e '
using Pkg; Pkg.instantiate()
ENV["QCFD_HOME"] = pwd()
ENV["QCFD_SRC"] = pwd() * "/src/"
include("src/CLBE/clbe_tg2d_run.jl")
Main.main(
    nx=3,
    ny=3,
    amplitude=0.02,
    local_n_time=10,
    l_plot=false,
    boundary_setup=false,
    coeff_method=:numerical,
    local_truncation_order=3,
    reference_model=:direct_lbe,
    integrator=:matrix_exponential,
    direct_lbe_integrator=:exponential_euler,
)
'
```

---

## Source file map

### Layer 1 — Driver and configuration

| File | Role | MPI / GPU notes |
|---|---|---|
| `src/CLBE/clbe_tg2d_run.jl` | Top-level D2Q9 CLBE driver. Builds the Taylor-Green initial condition, calls `prepare_d2q9_carleman_runtime`, calls `timeMarching_state_CLBM_sparse`, computes diagnostics and optionally saves HDF5 output. | **MPI**: domain decomposition of the `nx×ny` physical grid lives here — split the grid across ranks before calling the runtime builder. |
| `src/CLBE/clbe_config_2D.jl` | Global simulation parameters: `tau_value`, `dt`, `LX`/`LY`, `Q`, `D`, `truncation_order`, `poly_order`, `coeff_generation_method`. | No compute; edit to set problem size before scaling up. |

### Layer 2 — Carleman operator assembly

| File | Role | MPI / GPU notes |
|---|---|---|
| `src/CLBE/carleman_transferA_ngrid.jl` | Builds per-site collision coefficient matrices `F1`, `F2`, `F3` (numerical path via `numerical_carleman_coefficients`) and lifts them to the full n-point sparse operators `F1_ngrid`, `F2_ngrid`, `F3_ngrid` via `sparse_lift_local_collision`. Also provides `get_coeff_LBM_Fi_ngrid` and `build_sparse_ngrid_coefficients`. | **GPU**: the inner loops in `sparse_lift_local_collision` (orders 1, 2, 3) are the primary GPU kernel targets — each loop body is an independent index computation with no data dependency. |
| `src/CLBE/carleman_transferA.jl` | Assembles the full block-structured Carleman collision matrix `C` and streaming lift `S` via iterated Kronecker products (`sum_Kron_kth_identity`, `Kron_kth_identity`, `carleman_C`, `carleman_S`). Defines `carleman_V` (state lifting φ → V) and `carleman_C_dim` (lifted dimension). | **GPU**: the Kronecker product chains in `Kron_kth_identity` / `sum_Kron_kth_identity` are the second GPU target — replace with batched tensor-product kernels (e.g. `CUDA.jl` or `KernelAbstractions.jl`). |

### Layer 3 — Time marching

| File | Role | MPI / GPU notes |
|---|---|---|
| `src/CLBE/timeMarching.jl` | Sparse CLBE time loop (`timeMarching_state_CLBM_sparse`). Assembles the full generator `C_full = C − S` via `build_full_clbe_generator_sparse`, then advances with explicit Euler or Krylov matrix-exponential (`krylov_expv_affine`). Also contains `carleman_C_sparse` (block-by-block sparse assembly) and `carleman_S_sparse`. | **GPU**: the sparse matvec `C_sparse * V` inside the time loop is the primary GPU/MPI kernel — port `C_sparse` to `CUDA.CUSPARSE` or distribute across MPI ranks with `PartitionedArrays.jl`. **MPI**: the distributed sparse matvec across ranks maps directly to the loop at lines 110–122. |
| `src/LBE/direct_LBE.jl` | Direct semi-discrete n-point LBE reference (`timeMarching_direct_LBE_ngrid`). Euler or exponential-Euler (ETD/Krylov) step on `∂_t φ = −Sφ + F¹φ + F²φ⊗² + F³φ⊗³`. The apples-to-apples validation target for CLBE. | **GPU**: `kron(phi, phi)` and `kron(phi, kron(phi, phi))` in `direct_lbe_rhs_ngrid` are GPU-friendly tensor contractions — replace with batched outer-product kernels. |

### Layer 4 — LBM physics and utilities

| File | Role | MPI / GPU notes |
|---|---|---|
| `src/LBM/tg_d2q9_lbm_run.jl` | Pure numerical D2Q9 LBM Taylor-Green runner (`run_tg_d2q9_lbm`). Classical LBM baseline and validation reference. | Classical LBM streaming+collision loop — straightforward GPU port with `KernelAbstractions.jl`. |
| `src/LBM/lbm_const_sym.jl` | Symbolic D2Q9 weights, velocities, and equilibrium constants (`lbm_const_sym`). Used by the symbolic coefficient path. | Setup only; no runtime compute. |
| `src/LBM/cal_feq.jl` | Computes the equilibrium distribution `f_eq(ρ, u)`. Used to set the Taylor-Green initial condition. | Called once at initialization; no parallelism needed. |
| `src/CLBE/collision_sym.jl` | Symbolic BGK collision operator (`collision`). Entry point for the symbolic Carleman coefficient path (`coeff_method=:symbolic`). | Symbolic preprocessing only; runs once offline. |
| `src/CLBE/LBM_const_subs.jl` | Substitutes numerical values of `τ` into symbolic collision expressions. | Symbolic preprocessing only. |
| `julia_lib/matrix_kit.jl` | Utility: `Kron_kth` (iterated Kronecker product, dense). Used throughout the Carleman assembly. | Replace with sparse/GPU Kronecker (`Kron_kth_sparse` in `timeMarching.jl` already exists as the sparse counterpart). |

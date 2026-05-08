include("src/CLBE/cfull_norm_sweep_utils.jl")

const CFSU = CFullNormSweepUtils

default_h5_candidates() = [
    CFSU.data_output_path("d1q3_cfull_norm_vs_ngrid.h5"),
    CFSU.data_output_path("d2q9_cfull_norm_vs_grid.h5"),
]

function resolve_default_h5_path()
    existing = filter(isfile, default_h5_candidates())
    if length(existing) == 1
        return existing[1]
    elseif isempty(existing)
        error("No default HDF5 file found. Pass a path explicitly, e.g. `main(\"data/d1q3_cfull_norm_vs_ngrid.h5\")`.")
    else
        error("Multiple default HDF5 files found. Pass the desired file explicitly. Candidates: $(join(existing, ", "))")
    end
end

function main(h5_path::AbstractString=resolve_default_h5_path(); output_pdf::Union{Nothing,String}=nothing)
    return CFSU.plot_spectral_from_h5(h5_path; output_pdf=output_pdf)
end

if abspath(PROGRAM_FILE) == @__FILE__
    h5_path = length(ARGS) >= 1 ? ARGS[1] : resolve_default_h5_path()
    output_pdf = length(ARGS) >= 2 ? ARGS[2] : nothing
    main(h5_path; output_pdf=output_pdf)
end
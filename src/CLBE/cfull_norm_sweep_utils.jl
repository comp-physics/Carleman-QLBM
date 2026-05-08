module CFullNormSweepUtils

using Printf
using Random
using Dates
using HDF5
using PyPlot
using LinearAlgebra
using SparseArrays

export data_output_path,
       figure_output_path,
       power_iteration_spectral_norm,
       write_results_h5,
       load_results_h5,
       load_h5_metadata,
       summary_table_lines,
       write_summary_table,
       plot_norm_results,
       plot_spectral_from_h5,
       infer_model_name,
       run_main_if_script

function data_output_dir()
    outdir = joinpath(pwd(), "data")
    mkpath(outdir)
    return outdir
end

function tex_output_dir()
    if haskey(ENV, "TEXPATH")
        outdir = joinpath(ENV["TEXPATH"], "QC", "QCFD-CarlemanLBE", "figs")
    else
        outdir = data_output_dir()
    end
    mkpath(outdir)
    return outdir
end

data_output_path(filename::AbstractString) = joinpath(data_output_dir(), filename)
figure_output_path(filename::AbstractString) = joinpath(tex_output_dir(), filename)

function power_iteration_spectral_norm(A::SparseMatrixCSC; maxiter::Int=80, tol::Float64=1e-8, seed::Int=1234)
    n = size(A, 2)
    rng = MersenneTwister(seed)
    x = randn(rng, n)
    x ./= norm(x)

    sigma_prev = 0.0
    converged = false
    iterations = 0

    for iter = 1:maxiter
        y = A * x
        sigma = norm(y)
        iterations = iter
        if iszero(sigma)
            return 0.0, true, iter
        end

        z = transpose(A) * y
        z_norm = norm(z)
        if iszero(z_norm)
            return sigma, true, iter
        end

        x .= z ./ z_norm
        rel_change = abs(sigma - sigma_prev) / max(sigma, eps(Float64))
        if iter > 1 && rel_change < tol
            converged = true
            sigma_prev = sigma
            break
        end
        sigma_prev = sigma
    end

    return sigma_prev, converged, iterations
end

function infer_model_name(Q)
    q_int = try
        Int(Q)
    catch
        return string(Q)
    end
    return q_int == 3 ? "D1Q3" : q_int == 9 ? "D2Q9" : "Q=$q_int"
end

function _encode_h5_column(values)
    if isempty(values)
        return String[], "String"
    elseif all(v -> v isa Bool, values)
        return Bool[values...], "Bool"
    elseif all(v -> v isa Integer, values)
        return Int[values...], "Int"
    elseif all(v -> v isa AbstractFloat, values)
        return Float64[values...], "Float64"
    elseif all(v -> v isa Symbol, values)
        return String[string(v) for v in values], "Symbol"
    else
        return String[string(v) for v in values], "String"
    end
end

function _decode_h5_column(values, type_tag::AbstractString)
    if type_tag == "Bool"
        return Bool.(values)
    elseif type_tag == "Int"
        return Int.(values)
    elseif type_tag == "Float64"
        return Float64.(values)
    elseif type_tag == "Symbol"
        return Symbol.(String.(values))
    else
        return String.(values)
    end
end

_read_attr(attrs, key::AbstractString) = read(attrs[key])

function write_results_h5(path::AbstractString, results; metadata::Dict{String,<:Any}=Dict{String,Any}())
    h5open(path, "w") do h5
        attrs = attributes(h5)
        attrs["generated"] = Dates.format(now(), dateformat"yyyy-mm-dd HH:MM:SS")
        for (k, v) in metadata
            attrs[k] = v
        end

        g = create_group(h5, "results")
        if isempty(results)
            attrs["field_names_csv"] = ""
            return
        end

        field_names = String.(collect(keys(results[1])))
        attrs["field_names_csv"] = join(field_names, ",")
        gattrs = attributes(g)

        for name in field_names
            values = [getfield(r, Symbol(name)) for r in results]
            encoded, type_tag = _encode_h5_column(values)
            g[name] = encoded
            gattrs[name * "_type"] = type_tag
        end
    end
    println("Saved HDF5 results to: $path")
end

function load_results_h5(path::AbstractString; with_metadata::Bool=false)
    h5open(path, "r") do h5
        attrs = attributes(h5)
        field_names_csv = String(_read_attr(attrs, "field_names_csv"))
        field_names = isempty(field_names_csv) ? String[] : split(field_names_csv, ",")
        g = h5["results"]
        gattrs = attributes(g)

        columns = Dict{String,Vector}()
        for name in field_names
            raw = read(g[name])
            type_tag = String(_read_attr(gattrs, name * "_type"))
            columns[name] = _decode_h5_column(raw, type_tag)
        end

        results = NamedTuple[]
        if !isempty(field_names)
            nt_names = Tuple(Symbol.(field_names))
            nt_type = NamedTuple{nt_names}
            nrows = length(columns[first(field_names)])
            for i in 1:nrows
                row_values = Tuple(columns[name][i] for name in field_names)
                push!(results, nt_type(row_values))
            end
        end

        if with_metadata
            metadata = Dict{String,Any}()
            for key in keys(attrs)
                metadata[String(key)] = _read_attr(attrs, String(key))
            end
            return results, metadata
        end
        return results
    end
end

function load_h5_metadata(path::AbstractString)
    _, metadata = load_results_h5(path; with_metadata=true)
    return metadata
end

function _summary_cell(value)
    if value isa AbstractFloat
        return isfinite(value) ? @sprintf("%.8e", value) : "NaN"
    end
    return string(value)
end

function summary_table_lines(headers::Vector{String}, rows::Vector{<:AbstractVector{<:Any}}; left_align::Vector{Bool}=falses(length(headers)))
    string_rows = [[_summary_cell(v) for v in row] for row in rows]
    widths = [length(headers[j]) for j in eachindex(headers)]
    for row in string_rows
        for j in eachindex(headers)
            widths[j] = max(widths[j], length(row[j]))
        end
    end

    format_cell(j, s) = left_align[j] ? rpad(s, widths[j]) : lpad(s, widths[j])
    lines = String[]
    push!(lines, join([format_cell(j, headers[j]) for j in eachindex(headers)], "  "))
    for row in string_rows
        push!(lines, join([format_cell(j, row[j]) for j in eachindex(headers)], "  "))
    end
    return lines
end

function write_summary_table(path::AbstractString; title::AbstractString, intro_lines::Vector{String}=String[], headers::Vector{String}, rows::Vector{<:AbstractVector{<:Any}}, left_align::Vector{Bool}=falses(length(headers)))
    table_lines = summary_table_lines(headers, rows; left_align=left_align)
    open(path, "w") do io
        println(io, title)
        println(io, "Generated: ", Dates.format(now(), dateformat"yyyy-mm-dd HH:MM:SS"))
        for line in intro_lines
            println(io, line)
        end
        println(io)
        for line in table_lines
            println(io, line)
        end
    end
end

function plot_norm_results(results, output_pdf::AbstractString;
    x_values::Vector,
    x_labels=nothing,
    xlabel_text::AbstractString,
    title_text::AbstractString,
    include_inf::Bool=true,
    include_spectral::Bool=true,
    xscale_base=nothing,
    yscale_base::Real=10)

    close("all")
    figure(figsize=(7.6, 4.8))

    if include_inf
        inf_pairs = [(x_values[i], results[i].inf_norm) for i in eachindex(results) if hasproperty(results[i], :inf_norm) && isfinite(results[i].inf_norm)]
        if !isempty(inf_pairs)
            plot([p[1] for p in inf_pairs], [p[2] for p in inf_pairs], "--s",
                color="tab:red", linewidth=1.8, markersize=6,
                label="||C_full||∞")
        end
    end

    if include_spectral
        spectral_pairs = [(x_values[i], results[i].spectral_est) for i in eachindex(results) if hasproperty(results[i], :spectral_est) && isfinite(results[i].spectral_est)]
        if !isempty(spectral_pairs)
            plot([p[1] for p in spectral_pairs], [p[2] for p in spectral_pairs], "-o",
                color="tab:blue", linewidth=2.0, markersize=6,
                label="||C_full||₂ (power estimate)")
        end
    end

    if x_labels !== nothing
        xticks(x_values, x_labels)
    end
    if xscale_base !== nothing
        xscale("log", base=xscale_base)
    end
    yscale("log", base=yscale_base)
    xlabel(xlabel_text)
    ylabel("matrix norm")
    title(title_text)
    grid(true, which="both", alpha=0.3)
    legend(loc="best", fontsize=10)
    tight_layout()
    savefig(output_pdf, bbox_inches="tight")
    println("Saved figure to: $output_pdf")
end

function plot_spectral_from_h5(h5_path::AbstractString; output_pdf::Union{Nothing,AbstractString}=nothing)
    results, metadata = load_results_h5(h5_path; with_metadata=true)
    if isempty(results)
        error("No results found in HDF5 file: $h5_path")
    end

    q_val = get(metadata, "Q", missing)
    model_name = q_val === missing ? "Carleman" : infer_model_name(q_val)

    first_result = results[1]
    if hasproperty(first_result, :ngrid)
        x_values = [r.ngrid for r in results]
        x_labels = nothing
        xlabel_str = "ngrid"
        title_str = "$model_name spectral norm from cached HDF5"
        xscale_base = 2
    elseif hasproperty(first_result, :nx) && hasproperty(first_result, :ny)
        x_values = [r.nx * r.ny for r in results]
        x_labels = ["$(r.nx)×$(r.ny)" for r in results]
        xlabel_str = "grid size"
        title_str = "$model_name spectral norm from cached HDF5"
        xscale_base = nothing
    else
        error("Unsupported HDF5 result layout for plotting: expected `ngrid` or `nx`/`ny` fields.")
    end

    if output_pdf === nothing
        base = splitext(basename(h5_path))[1]
        output_pdf = joinpath(dirname(h5_path), base * "_spectral_only.pdf")
    end

    plot_norm_results(results, output_pdf;
        x_values=x_values,
        x_labels=x_labels,
        xlabel_text=xlabel_str,
        title_text=title_str,
        include_inf=false,
        include_spectral=true,
        xscale_base=xscale_base)
    return results
end

function run_main_if_script(caller_file::AbstractString, main_fn::Function)
    if abspath(PROGRAM_FILE) == abspath(caller_file)
        main_fn()
    end
    return nothing
end

end
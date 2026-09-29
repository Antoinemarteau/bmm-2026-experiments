# Experiment 2 — sparsity and cost of the rotation maps.
#
# Run from the repository root:  julia --project=. exp2_sparsity.jl
#
# For every (space, D, r) with D ∈ {2,3}, r ∈ {1,2,3,4}, over ALL π ∈ S_{D+1}:
#   • sparsity of the signed row map rotation_map(rc, π): number of "hit" rows
#     (rows with > 1 nonzero); for the trimmed space also single-term rows with
#     sign −1 (counted separately); average and max nnz per row; total nnz per
#     matrix vs n²;
#   • the per-π hit fraction grouped by the number of "min-breaking" faces of
#     π. Definition used: a face is a vertex subset F ⊆ {1,…,D+1} with |F| ≥ 2
#     (the faces that carry 1-form bubbles); F is min-breaking for π iff
#     π(min F) ≠ min(π(F)), i.e. π does not map the least vertex of F to the
#     least vertex of π(F).
#   • microbenchmark (@elapsed medians over 1000 reps after warmup, fixed
#     cyclic permutation π = (2,3,…,D+1,1)):
#       - rotation_map cold: first call on a fresh RotationCache (cache
#         construction excluded);
#       - rotation_map hot: memoised lookup (inner loop of 100 calls / rep);
#       - dense build: rotation_change_of_basis(b, π), one call / rep;
#       - dense n×n mat-vec (mul!, inner loop of 100 calls / rep), for scale.
#
# Output: results/exp2_sparsity.json and exp2_table.md
# Deterministic: mat-vec input vector from MersenneTwister(42).

using Gridap.Polynomials: BarycentricPΛBasis, BarycentricPmΛBasis,
                          RotationCache, rotation_map, rotation_change_of_basis
using Combinatorics: permutations, combinations
using LinearAlgebra
using Random
using Statistics: mean, median
using Printf

# The timings below are µs- and ns-scale, so the table is only interpretable
# next to the machine that produced it. Each probe is guarded so a missing
# introspection API degrades one field instead of losing the measurements.
_probe(f) = try
  string(f())
catch
  "unknown"
end

cpu_model()   = _probe(() -> strip(first(Sys.cpu_info()).model))
blas_threads() = _probe(LinearAlgebra.BLAS.get_num_threads)

machine_line() = string(cpu_model(), ", ", Sys.CPU_THREADS, " logical cores, Julia ", VERSION,
                        ", ", Threads.nthreads(), " Julia thread(s), ",
                        blas_threads(), " BLAS thread(s)")

# ── Minimal JSON writer (no external deps) ────────────────────────────────────
_json(io::IO, ::Nothing)          = print(io, "null")
_json(io::IO, x::Bool)            = print(io, x ? "true" : "false")
_json(io::IO, x::Integer)         = print(io, x)
_json(io::IO, x::AbstractFloat)   = print(io, isfinite(x) ? repr(x) : "\"$(x)\"")
_json(io::IO, s::AbstractString)  = print(io, '"', replace(s, "\\"=>"\\\\", "\""=>"\\\""), '"')
function _json(io::IO, v::AbstractVector{<:Pair})   # ordered JSON object
  print(io, '{')
  for (i, (k, x)) in enumerate(v)
    i > 1 && print(io, ',')
    _json(io, String(k)); print(io, ':'); _json(io, x)
  end
  print(io, '}')
end
function _json(io::IO, v::AbstractVector)
  print(io, '[')
  for (i, x) in enumerate(v)
    i > 1 && print(io, ',')
    _json(io, x)
  end
  print(io, ']')
end

# Vector-proxied 1-form bases (k = 1). flavor = :BMM selects the bases of the
# paper: support-indicator direction forms for the full space, bare barycentric
# monomials for the trimmed one; the rotation entries are then exactly ±1.
make_basis(name, D, r) =
  name == :rotating ? BarycentricPΛBasis( Val(D), Float64, r, 1; flavor=:BMM) :
                      BarycentricPmΛBasis(Val(D), Float64, r, 1; flavor=:BMM)

# Number of min-breaking faces of π (faces = vertex subsets of size ≥ 2).
function n_min_breaking_faces(π::Vector{Int}, D::Int)
  cnt = 0
  for s in 2:D+1, F in combinations(1:D+1, s)
    (π[minimum(F)] != minimum(π[F])) && (cnt += 1)
  end
  cnt
end
n_faces_ge2(D) = 2^(D + 1) - (D + 1) - 1   # subsets of size ≥ 2

# ── Microbenchmark helpers ────────────────────────────────────────────────────

const SINK = Ref(0)

function bench_cold(b, π; reps=1000)
  ts = Vector{Float64}(undef, reps)
  rotation_map(RotationCache(b), π)                 # warmup (compilation)
  for i in 1:reps
    rc = RotationCache(b)                           # fresh cache, not timed
    ts[i] = @elapsed (SINK[] += length(rotation_map(rc, π)))
  end
  median(ts)
end

function bench_hot(b, π; reps=1000, inner=100)
  rc = RotationCache(b)
  rotation_map(rc, π)                               # populate + warmup
  ts = Vector{Float64}(undef, reps)
  for i in 1:reps
    ts[i] = (@elapsed for _ in 1:inner
      SINK[] += length(rotation_map(rc, π))
    end) / inner
  end
  median(ts)
end

function bench_dense_build(b, π; reps=1000)
  rotation_change_of_basis(b, π)                    # warmup
  ts = Vector{Float64}(undef, reps)
  for i in 1:reps
    ts[i] = @elapsed (SINK[] += size(rotation_change_of_basis(b, π), 1))
  end
  median(ts)
end

function bench_matvec(C; reps=1000, inner=100)
  n = size(C, 1)
  rng = MersenneTwister(42)
  x = randn(rng, n); y = zeros(n)
  mul!(y, C, x)                                     # warmup
  ts = Vector{Float64}(undef, reps)
  for i in 1:reps
    ts[i] = (@elapsed for _ in 1:inner
      mul!(y, C, x)
    end) / inner
  end
  median(ts)
end

# ── The experiment ────────────────────────────────────────────────────────────

records = Vector{Any}()

for name in (:rotating, :trimmed), D in (2, 3), r in (1, 2, 3, 4)
  t0 = time()
  b = make_basis(name, D, r)
  n = length(b)
  rc = RotationCache(b)
  perms = [collect(p) for p in permutations(1:D+1)]

  hit_rows_tot = 0                  # rows with > 1 nonzero, summed over π
  neg_single_tot = 0                # trimmed: single-term rows with sign −1
  nnz_tot = 0
  max_nnz_row = 0
  max_hit_rows = 0
  # per-π hit fraction grouped by number of min-breaking faces
  by_mb = Dict{Int,Vector{Float64}}()

  for π in perms
    m = rotation_map(rc, π)
    hits = count(row -> length(row) > 1, m)
    negs = count(row -> length(row) == 1 && row[1][1] == -1.0, m)
    nnzπ = sum(length, m)
    hit_rows_tot += hits
    neg_single_tot += negs
    nnz_tot += nnzπ
    max_nnz_row = max(max_nnz_row, maximum(length, m))
    max_hit_rows = max(max_hit_rows, hits)
    push!(get!(by_mb, n_min_breaking_faces(π, D), Float64[]), hits / n)
  end

  nπ = length(perms)
  mb_table = [Pair{String,Any}[
      "min_breaking_faces" => k,
      "n_perms" => length(v),
      "mean_hit_fraction" => mean(v),
      "min_hit_fraction" => minimum(v),
      "max_hit_fraction" => maximum(v),
    ] for (k, v) in sort(collect(by_mb); by=first)]

  # microbenchmarks with the fixed cyclic permutation
  πc = [collect(2:D+1); 1]
  t_cold  = bench_cold(b, πc)
  t_hot   = bench_hot(b, πc)
  t_dense = bench_dense_build(b, πc)
  t_mv    = bench_matvec(rotation_change_of_basis(b, πc))

  rec = Pair{String,Any}[
    "space" => String(name), "D" => D, "r" => r, "n" => n,
    "n_permutations" => nπ,
    "n_faces_dim_ge1" => n_faces_ge2(D),
    "hit_rows_total" => hit_rows_tot,
    "hit_rows_mean_per_perm" => hit_rows_tot / nπ,
    "hit_rows_max_per_perm" => max_hit_rows,
    "hit_row_fraction_mean" => hit_rows_tot / (nπ * n),
    "neg_single_rows_total" => (name == :trimmed ? neg_single_tot : nothing),
    "neg_single_rows_mean_per_perm" => (name == :trimmed ? neg_single_tot / nπ : nothing),
    "nnz_mean_per_row" => nnz_tot / (nπ * n),
    "nnz_max_per_row" => max_nnz_row,
    "nnz_mean_per_perm" => nnz_tot / nπ,
    "n_squared" => n^2,
    "density_mean" => nnz_tot / (nπ * n^2),
    "hit_fraction_by_min_breaking_faces" => mb_table,
    "bench_perm" => πc,
    "t_rotation_map_cold_s" => t_cold,
    "t_rotation_map_hot_s" => t_hot,
    "t_dense_build_s" => t_dense,
    "t_dense_matvec_s" => t_mv,
    "walltime_s" => round(time() - t0, digits=2),
  ]
  push!(records, rec)
  @printf("done: %s D=%d r=%d n=%3d  hit%%=%5.1f  nnz/row=%.3f  cold=%.2fus hot=%.3fus dense=%.2fus mv=%.3fus (%.1fs)\n",
    name, D, r, n, 100hit_rows_tot / (nπ * n), nnz_tot / (nπ * n),
    1e6t_cold, 1e6t_hot, 1e6t_dense, 1e6t_mv, time() - t0)
end

# ── Write JSON + markdown table ───────────────────────────────────────────────

out = Pair{String,Any}[
  "experiment" => "exp2_sparsity",
  "description" => "Sparsity of rotation_map rows over all pi in S_{D+1}; hit fraction " *
                   "grouped by number of min-breaking faces (F with pi(min F) != min pi(F), |F|>=2); " *
                   "@elapsed-median microbenchmarks (1000 reps, warmup, cyclic pi).",
  "julia" => string(VERSION),
  "machine" => machine_line(),
  "records" => records,
]

resdir = joinpath(@__DIR__, "results")
mkpath(resdir)
open(joinpath(resdir, "exp2_sparsity.json"), "w") do io
  _json(io, out); println(io)
end

fmt_t(t) = t < 1e-6 ? @sprintf("%.0f ns", 1e9t) : @sprintf("%.2f µs", 1e6t)

open(joinpath(resdir, "exp2_table.md"), "w") do io
  println(io, "### Experiment 2 — sparsity and cost of the rotation maps\n")
  println(io, "Sparsity statistics aggregate over all (D+1)! permutations π. A \"hit\" row has ",
    "more than one nonzero; for the trimmed space, single-term rows with sign −1 are counted ",
    "separately (\"neg 1-rows\"). nnz counts signed entries of rotation_map(rc, π) (all entries ±1). ",
    "Timings are medians of 1000 @elapsed repetitions after warmup, fixed cyclic π = (2,…,D+1,1): ",
    "cold = first rotation_map call on a fresh RotationCache; hot = memoised call; dense = ",
    "rotation_change_of_basis build; mat-vec = dense n×n mul! (for scale).\n")
  println(io, "Timings measured on: ", machine_line(), ". They are strongly ",
    "machine-dependent; the sparsity columns are not.\n")
  println(io, "| space | D | r | n | hit rows: mean/π (frac) | max/π | neg 1-rows/π | nnz/row mean (max) | nnz/π vs n² | cold | hot | dense build | dense mat-vec |")
  println(io, "|---|---|---|---|---|---|---|---|---|---|---|---|---|")
  for rec in records
    d = Dict(rec)
    neg = d["neg_single_rows_mean_per_perm"]
    println(io, "| ", d["space"], " | ", d["D"], " | ", d["r"], " | ", d["n"],
      " | ", @sprintf("%.2f (%.1f%%)", d["hit_rows_mean_per_perm"], 100d["hit_row_fraction_mean"]),
      " | ", d["hit_rows_max_per_perm"],
      " | ", neg === nothing ? "—" : @sprintf("%.2f", neg),
      " | ", @sprintf("%.3f (%d)", d["nnz_mean_per_row"], d["nnz_max_per_row"]),
      " | ", @sprintf("%.1f / %d (%.1f%%)", d["nnz_mean_per_perm"], d["n_squared"], 100d["density_mean"]),
      " | ", fmt_t(d["t_rotation_map_cold_s"]),
      " | ", fmt_t(d["t_rotation_map_hot_s"]),
      " | ", fmt_t(d["t_dense_build_s"]),
      " | ", fmt_t(d["t_dense_matvec_s"]), " |")
  end
  println(io, "\n#### Hit fraction vs number of min-breaking faces of π\n")
  println(io, "Face = vertex subset F, |F| ≥ 2 (the faces carrying 1-form bubbles; ",
    "2D: 4 faces, 3D: 11). F is min-breaking for π iff π(min F) ≠ min π(F). ",
    "Entries: mean hit-row fraction over the permutations with that count (min–max in brackets); ",
    "— where no permutation has that count.\n")
  for (name, D) in ((:rotating, 2), (:trimmed, 2), (:rotating, 3), (:trimmed, 3))
    sel = [Dict(rec) for rec in records if Dict(rec)["space"] == String(name) && Dict(rec)["D"] == D]
    isempty(sel) && continue
    mbvals = sort(unique(vcat([[Dict(e)["min_breaking_faces"] for e in d["hit_fraction_by_min_breaking_faces"]]
                               for d in sel]...)))
    println(io, "\n**", name, ", D=", D, "** (columns: #min-breaking faces; #π in brackets)\n")
    hdr = ["r"; [string(k, " (", length([1 for π in permutations(1:D+1)
                   if n_min_breaking_faces(collect(π), D) == k]), "π)") for k in mbvals]]
    println(io, "| ", join(hdr, " | "), " |")
    println(io, "|", repeat("---|", length(hdr)))
    for d in sel
      tbl = Dict(Dict(e)["min_breaking_faces"] => Dict(e) for e in d["hit_fraction_by_min_breaking_faces"])
      cells = String[string(d["r"])]
      for k in mbvals
        if haskey(tbl, k)
          e = tbl[k]
          push!(cells, @sprintf("%.3f [%.3f–%.3f]", e["mean_hit_fraction"],
                                e["min_hit_fraction"], e["max_hit_fraction"]))
        else
          push!(cells, "—")
        end
      end
      println(io, "| ", join(cells, " | "), " |")
    end
  end
end

println("wrote ", joinpath(resdir, "exp2_sparsity.json"))

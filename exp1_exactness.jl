# Experiment 1 — exactness and representation property of the rotation
# change of basis C(π), for both the rotating (untrimmed P_rΛ¹) and the
# trimmed (P_r⁻Λ¹) bases.
#
# Run from the repository root:  julia --project=. exp1_exactness.jl
#
# For every (space, D, r) with D ∈ {2,3}, r ∈ {1,2,3,4}, over ALL π, τ ∈ S_{D+1}:
#   (a) every entry of C(π) lies in {−1, 0, +1} EXACTLY (Float64 ==, no tol);
#   (b) C(π)·C(invperm(π)) == I in exact integer arithmetic (entries rounded to
#       Int with an exact-roundtrip guard), plus the Float64 residual;
#   (c) the representation property C(π∘τ) vs products of C's — the correct
#       composition order is determined empirically, then verified exactly in
#       integer arithmetic over all 576 (2D: 36) pairs;
#   (d) pointwise pullback residual at 20 deterministic interior points,
#       against the numeric-pullback oracle of
#       test/PolynomialsTests/BarycentricPΛBases.jl.
#
# Output: results/exp1_exactness.json and exp1_table.md
#
# Deterministic: the only randomness is MersenneTwister(42) for the 20
# evaluation points (fixed seed).

using Gridap.Polynomials: BarycentricPΛBasis, BarycentricPmΛBasis,
                          bubble_entries, rotation_change_of_basis
using Gridap.Fields: Point
using Combinatorics: permutations, multinomial
using LinearAlgebra
using Random
using Printf

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
json_string(x) = (io = IOBuffer(); _json(io, x); String(take!(io)))

# ── Pullback oracles (identical formulas to the basis test suite) ─────────────
#
# A basis 1-form is written in the ambient barycentric frame as Σᵢ cᵢ dλⁱ, with
# N = D+1 coefficients. Its vector proxy on the reference simplex follows from
# dλ¹ = −Σⱼ dxʲ and dλʲ⁺¹ = dxʲ:  proxy[j] = c[j+1] − c[1].

# Gridap convention: λ = (1−Σx, x…), so V_1 = origin and V_{j+1} = e_j.
to_barycentric(x) = (1 - sum(x.data), x.data...)

# Vector proxy (as a plain Vector) of the ambient 1-form with coefficients c.
proxy(c::AbstractVector) = [c[j+1] - c[1] for j in 1:length(c)-1]

function rotating_psi(f::Vector{Int}, k::Int, α::Vector{Int}, N::Int)
  c = zeros(N)
  c[k] += 1.0
  supp = [i for i in 1:N if α[i] > 0]
  if !isempty(supp)
    w = (k in supp ? 1.0 : 0.0) / length(supp)
    for i in f; c[i] -= w; end
  end
  c
end

function oracle_eval(f::Vector{Int}, k::Int, α::Vector{Int}, x)   # untrimmed
  λ = to_barycentric(x); N = length(λ)
  val = multinomial(α...) * prod(λ[i]^α[i] for i in 1:N)   # Bernstein B_α
  val * proxy(rotating_psi(f, k, α, N))
end

function oracle_eval(f::Vector{Int}, e::Tuple{Int,Int}, α::Vector{Int}, x)  # trimmed
  λ = to_barycentric(x); N = length(λ)
  e1, e2 = e
  c = zeros(eltype(λ), N); c[e2] += λ[e1]; c[e1] -= λ[e2]   # Whitney form ϕ^e
  val = prod(λ[i]^α[i] for i in 1:N)      # BARE monomial λ^α (flavor :BMM)
  val * proxy(c)
end

# Vertex V_1 = origin, V_{j+1} = e_j (Gridap barycentric convention).
ref_vertices(D) = [[i == j + 1 ? 1.0 : 0.0 for j in 1:D] for i in 1:D+1]

function affine_map(p::Vector{Int}, D::Int)
  V = ref_vertices(D)
  A = x -> begin
    λ = to_barycentric(x)
    y = zeros(D)
    for i in 1:D+1, j in 1:D
      y[j] += λ[i] * V[p[i]][j]
    end
    Point(y...)
  end
  x0 = A(Point(zeros(D)...))
  J = zeros(D, D)
  for j in 1:D
    e = zeros(D); e[j] = 1.0
    J[:, j] = collect(A(Point(e...)).data) .- collect(x0.data)
  end
  A, J
end

# max over basis functions μ and points x of |J_a'·w_μ(A_a x) − Σ_ν C[μ,ν] w_ν(x)|
function pullback_error(b, π, a, pts, D)
  entries = bubble_entries(b)
  C = rotation_change_of_basis(b, π)
  A, J = affine_map(a, D)
  err = 0.0
  for (μ, t) in enumerate(entries)
    for x in pts
      lhs = transpose(J) * oracle_eval(t..., A(x))   # covariant Piola pullback
      rhs = zeros(D)
      for (ν, s) in enumerate(entries)
        c = C[μ, ν]
        c == 0.0 && continue
        rhs .+= c .* oracle_eval(s..., x)
      end
      err = max(err, maximum(abs.(lhs .- rhs)))
    end
  end
  err
end

# 20 deterministic interior points of the reference simplex (fixed seed).
function simplex_points(D, npts, seed)
  rng = MersenneTwister(seed)
  pts = Vector{Any}(undef, npts)
  for j in 1:npts
    e = -log.(rand(rng, D + 1))
    λ = e ./ sum(e)
    pts[j] = Point(λ[2:end]...)
  end
  pts
end

# ── The experiment ────────────────────────────────────────────────────────────

# Vector-proxied 1-form bases (k = 1). flavor = :BMM selects the bases of the
# paper (and of the oracles above): support-indicator direction forms for the
# full space, bare barycentric monomials for the trimmed one.
make_basis(name, D, r) =
  name == :rotating ? BarycentricPΛBasis( Val(D), Float64, r, 1; flavor=:BMM) :
                      BarycentricPmΛBasis(Val(D), Float64, r, 1; flavor=:BMM)

const NPTS = 20
const SEED = 42

records = Vector{Any}()

for name in (:rotating, :trimmed), D in (2, 3), r in (1, 2, 3, 4)
  t0 = time()
  b = make_basis(name, D, r)
  n = length(b)
  perms = [collect(p) for p in permutations(1:D+1)]

  # Dense C(π) for every π, plus exact Int copies.
  Cf = Dict{Vector{Int},Matrix{Float64}}()
  Ci = Dict{Vector{Int},Matrix{Int}}()
  entries_exact = true
  roundtrip_exact = true
  for π in perms
    C = rotation_change_of_basis(b, π)
    entries_exact &= all(c -> c == -1.0 || c == 0.0 || c == 1.0, C)
    K = round.(Int, C)
    roundtrip_exact &= all(Float64.(K) .== C)
    Cf[π] = C; Ci[π] = K
  end

  # (b) inverse: C(π)·C(invperm(π)) == I, exact Int arithmetic + Float residual.
  Id = Matrix{Int}(I, n, n)
  inverse_exact = true
  inverse_resid = 0.0
  for π in perms
    ιπ = invperm(π)
    inverse_exact &= (Ci[π] * Ci[ιπ] == Id)
    inverse_resid = max(inverse_resid, maximum(abs.(Cf[π] * Cf[ιπ] - I)))
  end

  # (c) representation: σ = π∘τ (σ[i] = π[τ[i]]). Test both matrix orders.
  n_pairs = 0; ok_CτCπ = 0; ok_CπCτ = 0
  comp_resid_CτCπ = 0.0; comp_resid_CπCτ = 0.0
  for π in perms, τ in perms
    n_pairs += 1
    σ = π[τ]
    (Ci[σ] == Ci[τ] * Ci[π]) && (ok_CτCπ += 1)
    (Ci[σ] == Ci[π] * Ci[τ]) && (ok_CπCτ += 1)
    comp_resid_CτCπ = max(comp_resid_CτCπ, maximum(abs.(Cf[σ] - Cf[τ] * Cf[π])))
    comp_resid_CπCτ = max(comp_resid_CπCτ, maximum(abs.(Cf[σ] - Cf[π] * Cf[τ])))
  end
  if ok_CτCπ == n_pairs
    comp_order, comp_exact, comp_resid = "C(pi∘tau) = C(tau)·C(pi)", true, comp_resid_CτCπ
  elseif ok_CπCτ == n_pairs
    comp_order, comp_exact, comp_resid = "C(pi∘tau) = C(pi)·C(tau)", true, comp_resid_CπCτ
  else
    comp_order = "NEITHER order holds for all pairs " *
                 "(C(tau)C(pi): $ok_CτCπ/$n_pairs, C(pi)C(tau): $ok_CπCτ/$n_pairs)"
    comp_exact = false
    comp_resid = min(comp_resid_CτCπ, comp_resid_CπCτ)
  end

  # (d) pullback residual at 20 points, all π (oracle: pullback along A_{π⁻¹}).
  pts = simplex_points(D, NPTS, SEED)
  pb_resid = 0.0
  for π in perms
    pb_resid = max(pb_resid, pullback_error(b, π, invperm(π), pts, D))
  end

  rec = Pair{String,Any}[
    "space" => String(name), "D" => D, "r" => r, "n" => n,
    "n_permutations" => length(perms), "n_pairs" => n_pairs,
    "entries_exact" => entries_exact,
    "int_roundtrip_exact" => roundtrip_exact,
    "inverse_exact" => inverse_exact,
    "inverse_float_residual" => inverse_resid,
    "composition_order" => comp_order,
    "composition_exact" => comp_exact,
    "composition_float_residual" => comp_resid,
    "pullback_n_points" => NPTS,
    "pullback_point_seed" => SEED,
    "pullback_residual" => pb_resid,
    "walltime_s" => round(time() - t0, digits=2),
  ]
  push!(records, rec)
  println("done: $name D=$D r=$r  n=$n  entries=$(entries_exact) inv=$(inverse_exact) " *
          "comp=$(comp_exact) pb=$(pb_resid)  ($(round(time()-t0, digits=1))s)")
end

# ── Write JSON + markdown table ───────────────────────────────────────────────

out = Pair{String,Any}[
  "experiment" => "exp1_exactness",
  "description" => "Exactness/inverse/representation of C(pi) over all pi,tau in S_{D+1}; " *
                   "pullback residual vs numeric oracle at 20 fixed points.",
  "julia" => string(VERSION),
  "records" => records,
]

resdir = joinpath(@__DIR__, "results")
mkpath(resdir)
open(joinpath(resdir, "exp1_exactness.json"), "w") do io
  _json(io, out); println(io)
end

fmt(x) = x == 0.0 ? "0" : @sprintf("%.2e", x)
open(joinpath(resdir, "exp1_table.md"), "w") do io
  println(io, "### Experiment 1 — exactness and representation property of C(π)\n")
  println(io, "All checks quantify over every π (and every pair (π, τ)) in S_{D+1}. ",
              "\"exact\" columns are verified in integer arithmetic after an exact ",
              "Float64→Int roundtrip; residuals are Float64. Pullback residual: max over ",
              "all π and all basis functions of the pointwise defect against the numeric ",
              "pullback oracle at 20 fixed interior points.\n")
  println(io, "| space | D | r | n | entries ∈ {−1,0,1} exact | inverse exact (resid) | composition exact (resid) | pullback residual |")
  println(io, "|---|---|---|---|---|---|---|---|")
  for rec in records
    d = Dict(rec)
    println(io, "| ", d["space"], " | ", d["D"], " | ", d["r"], " | ", d["n"],
      " | ", d["entries_exact"] ? "yes" : "NO",
      " | ", d["inverse_exact"] ? "yes" : "NO", " (", fmt(d["inverse_float_residual"]), ")",
      " | ", d["composition_exact"] ? "yes" : "NO", " (", fmt(d["composition_float_residual"]), ")",
      " | ", fmt(d["pullback_residual"]), " |")
  end
  orders = unique([Dict(rec)["composition_order"] for rec in records])
  println(io, "\nEmpirically determined composition order (uniform across all rows): ",
          join(orders, "; "), ".")
end

println("wrote ", joinpath(resdir, "exp1_exactness.json"))

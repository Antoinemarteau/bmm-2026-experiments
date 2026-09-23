# Experiment 4 — convergence under h-refinement on sorted vs scrambled meshes.
#
# L² projection of a smooth manufactured 1-form onto the rotating (full) P_rΛ¹
# and trimmed P_r⁻Λ¹ H(curl)-conforming spaces on simplexified Cartesian meshes
# of the unit square (D=2, n = 4,8,16,32 per side, r = 1,2,3) and the unit cube
# (D=3, n = 2,4[,8], r = 1,2).  Every mesh is built twice from the SAME
# connectivity: once with each cell's vertex list sorted (the fast path, no
# change of basis) and once with each row permuted pseudo-randomly by a
# deterministic per-cell generator (seed 42).  The two spaces are identical as
# sets, so errors and the projections themselves must agree to machine
# precision; we report |err_sorted − err_scrambled| and the pointwise
# ‖u_h^sorted − u_h^scrambled‖_L² (per-cell quadrature through both cells'
# affine parametrizations).
#
# Manufactured field, given by its vector proxy (non-polynomial, NONZERO curl —
# note the sign: the textbook pair (sin cos, cos sin) is a gradient, curl u = 0,
# which would make the curl column vacuous):
#   D=2: u = (sin(πx)cos(πy), −cos(πx)sin(πy)),
#        curl u = 2π sin(πx)sin(πy)
#   D=3: u = (sin(πx)cos(πy), −cos(πx)sin(πy), sin(πx)cos(πz)),
#        curl u = (0, −π cos(πx)cos(πz), 2π sin(πx)sin(πy))
# In form language u is the 1-form with these components and curl u is the
# proxy of du: the scalar coefficient of dx¹∧dx² in 2D, (du₂₃, −du₁₃, du₁₂) in 3D.
#
# Run from the repository root:  julia --project=. exp4_convergence.jl
#   EXP4_QUICK=1 runs a tiny case list (smoke test), written with a _quick suffix.
# Outputs: results/exp4_convergence.json, results/summary_4.md

using Gridap
using Gridap.Geometry
using Gridap.ReferenceFEs
using Gridap.FESpaces
using Gridap.CellData
using Gridap.Fields
using Gridap.Arrays: Table, array_cache, getindex!
using LinearAlgebra

const SEED = 42

# ─── Mesh builders ───────────────────────────────────────────────────────────

"Node coordinates + cell connectivity of a simplexified Cartesian unit-D-cube mesh."
function base_mesh(D::Int, n::Int)
  dom   = D == 2 ? (0.0, 1.0, 0.0, 1.0) : (0.0, 1.0, 0.0, 1.0, 0.0, 1.0)
  cells = D == 2 ? (n, n) : (n, n, n)
  model = simplexify(CartesianDiscreteModel(dom, cells))
  grid  = get_grid(model)
  coords = collect(get_node_coordinates(grid))
  conn   = [collect(Int, row) for row in get_cell_node_ids(grid)]
  coords, conn
end

"""
Deterministic pseudo-random permutation of 1:m for cell c (Fisher–Yates driven
by an xorshift64 stream whose state mixes the global seed and the cell id).
Dependency-free so the script runs with only the package's own deps.
"""
function cell_perm(c::Int, m::Int; seed::Int=SEED)
  s = UInt64(seed) * 0x9e3779b97f4a7c15 + UInt64(c) * 0xbf58476d1ce4e5b9
  p = collect(1:m)
  for i in m:-1:2
    s ⊻= s << 13; s ⊻= s >> 7; s ⊻= s << 17
    j = Int(mod(s, UInt64(i))) + 1
    p[i], p[j] = p[j], p[i]
  end
  p
end

"UnstructuredDiscreteModel from explicit connectivity rows (arbitrary orderings)."
function build_model(D::Int, coords, conn)
  reffe = LagrangianRefFE(Float64, D == 2 ? TRI : TET, 1)
  grid  = UnstructuredGrid(coords, Table(conn), [reffe], fill(Int8(1), length(conn)))
  UnstructuredDiscreteModel(grid)
end

# ─── Manufactured fields ─────────────────────────────────────────────────────

u_ex_2d(x)  = VectorValue(sin(π*x[1])*cos(π*x[2]),
                          -cos(π*x[1])*sin(π*x[2]))
du_ex_2d(x) = 2π*sin(π*x[1])*sin(π*x[2])

u_ex_3d(x)  = VectorValue(sin(π*x[1])*cos(π*x[2]),
                          -cos(π*x[1])*sin(π*x[2]),
                          sin(π*x[1])*cos(π*x[3]))
du_ex_3d(x) = VectorValue(0.0,
                          -π*cos(π*x[1])*cos(π*x[3]),
                          2π*sin(π*x[1])*sin(π*x[2]))

# Guard the hand-written curls against Gridap's curl (AD) at one point.
let x2 = Point(0.13, 0.27), x3 = Point(0.13, 0.27, 0.41)
  @assert abs(du_ex_2d(x2) - curl(u_ex_2d)(x2)) < 1e-12
  @assert maximum(abs, Tuple(du_ex_3d(x3) - curl(u_ex_3d)(x3))) < 1e-12
end

# ─── L² projection + errors ──────────────────────────────────────────────────

"""
L² projection of u_fun onto FESpace(model, rf, :HCurl): assemble the mass
matrix a(u,v) = ∫ u⋅v (= ∫ vol_coeff(u ∧ ⋆v) for the 1-forms) and rhs
l(v) = ∫ u_ex⋅v, solve, and return (uh, L² error, curl-seminorm error, ndofs).
"""
function l2_project(model, rf, u_fun, du_fun; deg)
  Ω  = Triangulation(model)
  dΩ = Measure(Ω, deg)
  V  = FESpace(model, rf, conformity=:HCurl)
  U  = TrialFESpace(V)
  u_cf  = CellField(u_fun, Ω, PhysicalDomain())
  du_cf = CellField(du_fun, Ω, PhysicalDomain())
  a(u, v) = ∫( u ⋅ v ) * dΩ
  l(v)    = ∫( u_cf ⋅ v ) * dΩ
  uh = solve(AffineFEOperator(a, l, U, V))
  e  = uh - u_cf
  errL2 = sqrt(sum(∫( e ⋅ e ) * dΩ))
  ec = curl(uh) - du_cf
  errCurl = sqrt(sum(∫( ec ⋅ ec ) * dΩ))
  uh, errL2, errCurl, num_free_dofs(V)
end

# ─── Pointwise sorted-vs-scrambled distance ──────────────────────────────────

"Affine geometry (o, J) of every cell: F(x̂) = o + J x̂."
cell_affines(coords, conn, D) = map(conn) do g
  o = collect(coords[g[1]].data)
  J = hcat([collect(coords[g[j+1]].data) .- o for j in 1:D]...)
  (o, J)
end

"""
‖uh_a − uh_b‖_L² where cell c of both models is the same physical simplex with
different reference parametrizations: quadrature points of cell c in model a
are mapped to physical space and pulled back through model b's affine map.
"""
function l2_diff(uh_a, uh_b, coords, conn_a, conn_b, D; deg)
  quad = Quadrature(D == 2 ? TRI : TET, deg)
  q̂ = get_coordinates(quad); w = get_weights(quad)
  da = Gridap.CellData.get_data(uh_a); db = Gridap.CellData.get_data(uh_b)
  ca = array_cache(da); cb = array_cache(db)
  ga = cell_affines(coords, conn_a, D); gb = cell_affines(coords, conn_b, D)
  acc = 0.0
  for c in 1:length(conn_a)
    oa, Ja = ga[c]; ob, Jb = gb[c]
    detJa = abs(det(Ja))
    xb = [Point((Jb \ (oa .+ Ja * collect(q.data) .- ob))...) for q in q̂]
    fa = getindex!(ca, da, c); fb = getindex!(cb, db, c)
    va = evaluate(fa, q̂); vb = evaluate(fb, xb)
    for (i, ω) in enumerate(w)
      dv = va[i] - vb[i]
      acc += ω * detJa * (dv ⋅ dv)
    end
  end
  sqrt(abs(acc))
end

# ─── Driver ──────────────────────────────────────────────────────────────────

# Vector-proxied P_rΛ¹ (full) and P_r⁻Λ¹ (trimmed) H(curl) elements; the
# default flavor :BMM gives the bases of the paper.
simplex(D) = D == 2 ? TRI : TET
make_reffe(name, D, r) = name == :full ? RotatingPΛRefFE(Float64, simplex(D), r) :
                                         TrimmedPΛRefFE(Float64, simplex(D), r)

"One (D, space, r) refinement study; returns a vector of per-level records."
function run_case(D, name, r, ns)
  u_fun  = D == 2 ? u_ex_2d  : u_ex_3d
  du_fun = D == 2 ? du_ex_2d : du_ex_3d
  deg = max(2r + 2, 6)
  rf  = make_reffe(name, D, r)
  levels = Vector{Vector{Pair{String,Any}}}()
  for n in ns
    coords, conn = base_mesh(D, n)
    conn_s = [sort(row) for row in conn]
    conn_x = [row[cell_perm(c, length(row))] for (c, row) in enumerate(conn_s)]
    m_s = build_model(D, coords, conn_s)
    m_x = build_model(D, coords, conn_x)
    t = @elapsed begin
      uh_s, eL2_s, eC_s, nd = l2_project(m_s, rf, u_fun, du_fun; deg=deg)
      uh_x, eL2_x, eC_x, _  = l2_project(m_x, rf, u_fun, du_fun; deg=deg)
      duh = l2_diff(uh_s, uh_x, coords, conn_s, conn_x, D; deg=deg)
    end
    push!(levels, Pair{String,Any}[
      "n" => n, "h" => 1.0 / n, "ncells" => length(conn_s), "ndofs" => nd,
      "errL2_sorted" => eL2_s, "errL2_scrambled" => eL2_x,
      "errCurl_sorted" => eC_s, "errCurl_scrambled" => eC_x,
      "abs_diff_errL2" => abs(eL2_s - eL2_x),
      "abs_diff_errCurl" => abs(eC_s - eC_x),
      "l2_norm_uh_diff" => duh,
      "quad_degree" => deg, "walltime_s" => t,
    ])
    println("  D=$D $name r=$r n=$n: ndofs=$nd  L2=$(eL2_s)  curl=$(eC_s)  " *
            "|Δerr|=$(abs(eL2_s - eL2_x))  ‖Δuh‖=$duh  ($(round(t, digits=1))s)")
    flush(stdout)
  end
  levels
end

getlv(lv, k) = last(lv[findfirst(p -> first(p) == k, lv)])
rates(levels, key) = [log2(getlv(levels[i], key) / getlv(levels[i+1], key))
                      for i in 1:length(levels)-1]

# ─── Minimal JSON writer (no external deps) ──────────────────────────────────

json(x::AbstractString) = "\"" * replace(x, "\\" => "\\\\", "\"" => "\\\"") * "\""
json(x::Bool)    = x ? "true" : "false"
json(x::Integer) = string(x)
json(x::Real)    = isfinite(x) ? string(Float64(x)) : "null"
json(v::AbstractVector{<:Pair}) = "{" * join([json(String(first(p))) * ":" * json(last(p)) for p in v], ",") * "}"
json(v::AbstractVector) = "[" * join(json.(v), ",") * "]"

# ─── Formatting helpers ──────────────────────────────────────────────────────

fnum(x) = x == 0 ? "0" : string(round(x, sigdigits=4))
frate(x) = string(round(x, digits=2))

# ─── Main ────────────────────────────────────────────────────────────────────

function main()
  # EXP4_QUICK=1 runs a tiny case list (smoke test of the full pipeline).
  quick = get(ENV, "EXP4_QUICK", "0") == "1"
  cases = quick ? [(2, :full, 1, [4, 8]), (2, :trimmed, 1, [4, 8]), (3, :full, 1, [2])] : [
    (2, :full,    1, [4, 8, 16, 32]),
    (2, :full,    2, [4, 8, 16, 32]),
    (2, :full,    3, [4, 8, 16, 32]),
    (2, :trimmed, 1, [4, 8, 16, 32]),
    (2, :trimmed, 2, [4, 8, 16, 32]),
    (2, :trimmed, 3, [4, 8, 16, 32]),
    (3, :full,    1, [2, 4, 8]),
    (3, :full,    2, [2, 4]),
    (3, :trimmed, 1, [2, 4, 8]),
    (3, :trimmed, 2, [2, 4]),
  ]

  results = Vector{Vector{Pair{String,Any}}}()
  for (D, name, r, ns) in cases
    println("Case D=$D space=$name r=$r  n=$ns"); flush(stdout)
    levels = run_case(D, name, r, ns)
    push!(results, Pair{String,Any}[
      "D" => D, "space" => String(name), "r" => r,
      "expected_L2_rate" => name == :full ? r + 1 : r,
      "levels" => levels,
      "ratesL2_sorted" => rates(levels, "errL2_sorted"),
      "ratesCurl_sorted" => rates(levels, "errCurl_sorted"),
    ])
  end

  outdir = joinpath(@__DIR__, "results")
  mkpath(outdir)
  suffix = quick ? "_quick" : ""

  # ── JSON ──
  top = Pair{String,Any}[
    "experiment" => "exp4_convergence",
    "description" => "L2 projection of a smooth 1-form onto full (rotating) P_rΛ¹ and trimmed P_r⁻Λ¹ under h-refinement; sorted vs per-cell scrambled vertex orderings (seed $SEED).",
    "field_2d" => "u = (sin(pi x)cos(pi y), -cos(pi x)sin(pi y)) (curl u = 2pi sin sin)",
    "field_3d" => "u = (sin(pi x)cos(pi y), -cos(pi x)sin(pi y), sin(pi x)cos(pi z)) (curl u = (0, -pi cos cos, 2pi sin sin))",
    "scramble" => "per-cell Fisher-Yates via xorshift64, state mixes seed $SEED and cell id",
    "norms" => "errL2 = ||u - u_h||_L2; errCurl = ||curl(u - u_h)||_L2 (= ||d(u - u_h)||_L2); l2_norm_uh_diff = ||u_h^sorted - u_h^scrambled||_L2 (per-cell quadrature)",
    "cases" => results,
  ]
  open(joinpath(outdir, "exp4_convergence$suffix.json"), "w") do io
    write(io, json(top), "\n")
  end

  # ── Markdown summary ──
  io = IOBuffer()
  print(io, """
  # Experiment 4 — convergence under h-refinement, sorted vs scrambled meshes

  L² projection of the smooth 1-form (written through its vector proxy)
  **u = (sin(πx)cos(πy), −cos(πx)sin(πy))** (D=2; curl u = 2π sin(πx)sin(πy))
  and **u = (sin(πx)cos(πy), −cos(πx)sin(πy), sin(πx)cos(πz))** (D=3)
  onto the full (rotating) **P_rΛ¹** and trimmed **P_r⁻Λ¹** H(curl) spaces on
  simplexified Cartesian meshes of the unit square/cube.

  Each mesh is assembled twice from the same connectivity: **sorted** (every
  cell's vertex list ascending — the no-change-of-basis fast path) and
  **scrambled** (every row permuted by a deterministic per-cell pseudo-random
  permutation, seed $SEED). Same space, different local orderings: the "|Δ err|"
  column is |err_sorted − err_scrambled| and "‖Δu_h‖" is the pointwise
  L² norm of u_h^sorted − u_h^scrambled (per-cell quadrature through both
  parametrizations). Both should be machine zero.

  Expected L² rates: **r+1** for full P_rΛ¹, **r** for trimmed P_r⁻Λ¹.
  The curl column is the curl-seminorm error ‖curl(u − u_h)‖ = ‖d(u − u_h)‖ of the same L²
  projection (no curl control is imposed by the projection; rates are reported
  as measured).

  """)
  for res in results
    D = getlv(res, "D"); sp = getlv(res, "space"); r = getlv(res, "r")
    levels = getlv(res, "levels")
    spname = sp == "full" ? "full P_rΛ¹ (rotating)" : "trimmed P_r⁻Λ¹"
    exp_rate = getlv(res, "expected_L2_rate")
    print(io, "## D=$D, $spname, r=$r  (expected L² rate $exp_rate)\n\n")
    print(io, "| n | ndofs | ‖u−u_h‖_L² | rate | ‖d(u−u_h)‖ | rate | \\|Δ err_L²\\| | ‖Δu_h‖_L² |\n")
    print(io, "|--:|------:|-----------:|:----:|-----------:|:----:|-------------:|----------:|\n")
    rL2 = getlv(res, "ratesL2_sorted"); rC = getlv(res, "ratesCurl_sorted")
    for (i, lv) in enumerate(levels)
      print(io, "| ", getlv(lv, "n"), " | ", getlv(lv, "ndofs"),
            " | ", fnum(getlv(lv, "errL2_sorted")),
            " | ", i == 1 ? "—" : frate(rL2[i-1]),
            " | ", fnum(getlv(lv, "errCurl_sorted")),
            " | ", i == 1 ? "—" : frate(rC[i-1]),
            " | ", fnum(getlv(lv, "abs_diff_errL2")),
            " | ", fnum(getlv(lv, "l2_norm_uh_diff")), " |\n")
    end
    print(io, "\n")
  end
  print(io, """
  Notes:
  - Quadrature degree max(2r+2, 6) for assembly, error integration, and the
    pointwise sorted-vs-scrambled comparison.
  - The manufactured field is deliberately NOT the gradient pair
    (sin cos, cos sin), whose curl vanishes identically; the sign
    flip on the second component gives a nonzero smooth curl so the curl column
    is informative.
  - D=3 uses n per side ∈ {2,4,8} for r=1 and {2,4} for r=2 to bound runtime;
    coarse 3D levels are pre-asymptotic, so 3D rates are indicative only.
  """)
  open(joinpath(outdir, "summary_4$suffix.md"), "w") do f
    write(f, String(take!(io)))
  end
  println("Wrote ", joinpath(outdir, "exp4_convergence$suffix.json"))
  println("Wrote ", joinpath(outdir, "summary_4$suffix.md"))
end

main()

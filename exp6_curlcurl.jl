# Experiment 6 — a curl–curl problem on sorted vs scrambled meshes.
#
# Definite Maxwell / H(curl) problem on the unit square (D=2) and unit cube
# (D=3), discretised with the rotating (full) P_rΛ¹ and trimmed P_r⁻Λ¹
# H(curl)-conforming spaces:
#
#   curl curl u + u = f   in Ω,      n × curl u = 0   on ∂Ω,
#
# weak form  ∫ curl u ⋅ curl v + ∫ u ⋅ v = ∫ f ⋅ v  for all v in the full
# conforming space (the boundary condition is natural, no essential condition
# is imposed). In form language: ∫⟨du,dv⟩ + ∫⟨u,v⟩ = ∫⟨f,v⟩.
#
# Manufactured solutions (vector proxies), chosen so that n × curl u = 0 holds
# exactly on ∂Ω and curl curl u is a multiple of u:
#   D=2: u = (sin(πx)cos(πy), −cos(πx)sin(πy))      (the field of exp4),
#        curl u = 2π sin(πx)sin(πy) vanishes on ∂Ω,  curl curl u = 2π² u,
#        f = (1 + 2π²) u.
#   D=3: u = (−sin(πx)cos(πy)cos(πz), cos(πx)sin(πy)cos(πz), 0),
#        curl u = (π cos sin sin, π sin cos sin, −2π sin sin cos),
#        every tangential component of curl u vanishes on ∂Ω,  curl curl u = 3π² u,
#        f = (1 + 3π²) u.
#
# Every mesh is assembled twice from the SAME connectivity, sorted and scrambled
# (per-cell deterministic permutation, seed 42), exactly as in exp4. Reported:
# ‖u − u_h‖_L², ‖curl(u − u_h)‖_L² with rates, and ‖u_h^sorted − u_h^scrambled‖_L².
# Expected rates: L²  r+1 (full) / r (trimmed);  curl  r (both).
#
# Run from the repository root:  julia --project=. exp6_curlcurl.jl
#   EXP6_QUICK=1 runs a tiny case list (smoke test), written with a _quick suffix.
# Outputs: results/exp6_curlcurl.json, results/summary_6.md

using Gridap
using Gridap.Geometry
using Gridap.ReferenceFEs
using Gridap.FESpaces
using Gridap.CellData
using Gridap.Fields
using Gridap.Arrays: Table, array_cache, getindex!
using LinearAlgebra

const SEED = 42

# ─── Mesh builders (as in exp4) ──────────────────────────────────────────────

function base_mesh(D::Int, n::Int)
  dom   = D == 2 ? (0.0, 1.0, 0.0, 1.0) : (0.0, 1.0, 0.0, 1.0, 0.0, 1.0)
  cells = D == 2 ? (n, n) : (n, n, n)
  model = simplexify(CartesianDiscreteModel(dom, cells))
  grid  = get_grid(model)
  coords = collect(get_node_coordinates(grid))
  conn   = [collect(Int, row) for row in get_cell_node_ids(grid)]
  coords, conn
end

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

function build_model(D::Int, coords, conn)
  reffe = LagrangianRefFE(Float64, D == 2 ? TRI : TET, 1)
  grid  = UnstructuredGrid(coords, Table(conn), [reffe], fill(Int8(1), length(conn)))
  UnstructuredDiscreteModel(grid)
end

# ─── Manufactured solutions ──────────────────────────────────────────────────

u_ex_2d(x)  = VectorValue(sin(π*x[1])*cos(π*x[2]),
                          -cos(π*x[1])*sin(π*x[2]))
du_ex_2d(x) = 2π*sin(π*x[1])*sin(π*x[2])
f_ex_2d(x)  = (1 + 2π^2) * u_ex_2d(x)

u_ex_3d(x)  = VectorValue(-sin(π*x[1])*cos(π*x[2])*cos(π*x[3]),
                           cos(π*x[1])*sin(π*x[2])*cos(π*x[3]),
                           0.0)
du_ex_3d(x) = VectorValue( π*cos(π*x[1])*sin(π*x[2])*sin(π*x[3]),
                           π*sin(π*x[1])*cos(π*x[2])*sin(π*x[3]),
                          -2π*sin(π*x[1])*sin(π*x[2])*cos(π*x[3]))
f_ex_3d(x)  = (1 + 3π^2) * u_ex_3d(x)

# Guard the hand-written curls against Gridap's curl (AD) at one point: a
# wrong du_ex would leave the solve untouched and only corrupt the curl column.
let x2 = Point(0.13, 0.27), x3 = Point(0.13, 0.27, 0.41)
  @assert abs(du_ex_2d(x2) - curl(u_ex_2d)(x2)) < 1e-12
  @assert maximum(abs, Tuple(du_ex_3d(x3) - curl(u_ex_3d)(x3))) < 1e-12
end

# ─── Curl–curl solve + errors ────────────────────────────────────────────────

"""
Solve ∫ curl u ⋅ curl v + ∫ u ⋅ v = ∫ f ⋅ v on FESpace(model, rf, :HCurl) and
return (uh, L² error, curl-seminorm error, ndofs).
"""
function solve_curlcurl(model, rf, u_fun, du_fun, f_fun; deg)
  Ω  = Triangulation(model)
  dΩ = Measure(Ω, deg)
  V  = FESpace(model, rf, conformity=:HCurl)
  U  = TrialFESpace(V)
  u_cf  = CellField(u_fun,  Ω, PhysicalDomain())
  du_cf = CellField(du_fun, Ω, PhysicalDomain())
  f_cf  = CellField(f_fun,  Ω, PhysicalDomain())
  a(u, v) = ∫( curl(u) ⋅ curl(v) + u ⋅ v ) * dΩ
  l(v)    = ∫( f_cf ⋅ v ) * dΩ
  uh = solve(AffineFEOperator(a, l, U, V))
  e  = uh - u_cf
  errL2 = sqrt(sum(∫( e ⋅ e ) * dΩ))
  ec = curl(uh) - du_cf
  errCurl = sqrt(sum(∫( ec ⋅ ec ) * dΩ))
  uh, errL2, errCurl, num_free_dofs(V)
end

# ─── Pointwise sorted-vs-scrambled distance (as in exp4) ─────────────────────

cell_affines(coords, conn, D) = map(conn) do g
  o = collect(coords[g[1]].data)
  J = hcat([collect(coords[g[j+1]].data) .- o for j in 1:D]...)
  (o, J)
end

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

function run_case(D, name, r, ns)
  u_fun  = D == 2 ? u_ex_2d  : u_ex_3d
  du_fun = D == 2 ? du_ex_2d : du_ex_3d
  f_fun  = D == 2 ? f_ex_2d  : f_ex_3d
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
      uh_s, eL2_s, eC_s, nd = solve_curlcurl(m_s, rf, u_fun, du_fun, f_fun; deg=deg)
      uh_x, eL2_x, eC_x, _  = solve_curlcurl(m_x, rf, u_fun, du_fun, f_fun; deg=deg)
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

# ─── Minimal JSON writer ─────────────────────────────────────────────────────

json(x::AbstractString) = "\"" * replace(x, "\\" => "\\\\", "\"" => "\\\"") * "\""
json(x::Bool)    = x ? "true" : "false"
json(x::Integer) = string(x)
json(x::Real)    = isfinite(x) ? string(Float64(x)) : "null"
json(v::AbstractVector{<:Pair}) = "{" * join([json(String(first(p))) * ":" * json(last(p)) for p in v], ",") * "}"
json(v::AbstractVector) = "[" * join(json.(v), ",") * "]"

fnum(x) = x == 0 ? "0" : string(round(x, sigdigits=4))
frate(x) = string(round(x, digits=2))

# ─── Main ────────────────────────────────────────────────────────────────────

function main()
  quick = get(ENV, "EXP6_QUICK", "0") == "1"
  cases = quick ? [(2, :full, 1, [4, 8]), (2, :trimmed, 1, [4, 8]), (3, :full, 1, [2])] : [
    (2, :full,    1, [4, 8, 16, 32]),
    (2, :full,    2, [4, 8, 16, 32]),
    (2, :full,    3, [4, 8, 16, 32]),
    (2, :trimmed, 1, [4, 8, 16, 32]),
    (2, :trimmed, 2, [4, 8, 16, 32]),
    (2, :trimmed, 3, [4, 8, 16, 32]),
    (3, :full,    1, [2, 4, 8]),
    (3, :full,    2, [2, 4, 8]),
    (3, :trimmed, 1, [2, 4, 8]),
    (3, :trimmed, 2, [2, 4, 8]),
  ]

  results = Vector{Vector{Pair{String,Any}}}()
  for (D, name, r, ns) in cases
    println("Case D=$D space=$name r=$r  n=$ns"); flush(stdout)
    levels = run_case(D, name, r, ns)
    push!(results, Pair{String,Any}[
      "D" => D, "space" => String(name), "r" => r,
      "expected_L2_rate" => name == :full ? r + 1 : r,
      "expected_curl_rate" => r,
      "levels" => levels,
      "ratesL2_sorted" => rates(levels, "errL2_sorted"),
      "ratesCurl_sorted" => rates(levels, "errCurl_sorted"),
    ])
  end

  outdir = joinpath(@__DIR__, "results")
  mkpath(outdir)
  suffix = quick ? "_quick" : ""

  top = Pair{String,Any}[
    "experiment" => "exp6_curlcurl",
    "description" => "curl curl u + u = f with natural boundary condition n x curl u = 0, full (rotating) P_rΛ¹ and trimmed P_r⁻Λ¹, h-refinement, sorted vs per-cell scrambled vertex orderings (seed $SEED).",
    "field_2d" => "u = (sin(pi x)cos(pi y), -cos(pi x)sin(pi y)), f = (1 + 2 pi^2) u",
    "field_3d" => "u = (-sin(pi x)cos(pi y)cos(pi z), cos(pi x)sin(pi y)cos(pi z), 0), f = (1 + 3 pi^2) u",
    "scramble" => "per-cell Fisher-Yates via xorshift64, state mixes seed $SEED and cell id",
    "norms" => "errL2 = ||u - u_h||_L2; errCurl = ||curl(u - u_h)||_L2 (= ||d(u - u_h)||_L2); l2_norm_uh_diff = ||u_h^sorted - u_h^scrambled||_L2 (per-cell quadrature)",
    "julia" => string(VERSION),
    "cases" => results,
  ]
  open(joinpath(outdir, "exp6_curlcurl$suffix.json"), "w") do io
    write(io, json(top), "\n")
  end

  io = IOBuffer()
  print(io, """
  # Experiment 6 — curl–curl problem, sorted vs scrambled meshes

  curl curl u + u = f on the unit square/cube with the natural boundary
  condition n × curl u = 0, discretised with the full (rotating) **P_rΛ¹** and
  trimmed **P_r⁻Λ¹** H(curl) spaces on simplexified Cartesian meshes.
  Manufactured solutions (vector proxies): **u = (sin(πx)cos(πy), −cos(πx)sin(πy))**,
  f = (1+2π²)u (D=2); **u = (−sin(πx)cos(πy)cos(πz), cos(πx)sin(πy)cos(πz), 0)**,
  f = (1+3π²)u (D=3). In both cases n × curl u = 0 on ∂Ω exactly.

  Each mesh is assembled twice from the same connectivity: **sorted** and
  **scrambled** (seed $SEED, as in exp4). "‖Δu_h‖" is the L² norm of
  u_h^sorted − u_h^scrambled (per-cell quadrature through both parametrizations)
  and must be machine zero.

  Expected rates: L² **r+1** (full) / **r** (trimmed); curl-seminorm **r** (both).

  """)
  for res in results
    D = getlv(res, "D"); sp = getlv(res, "space"); r = getlv(res, "r")
    levels = getlv(res, "levels")
    spname = sp == "full" ? "full P_rΛ¹ (rotating)" : "trimmed P_r⁻Λ¹"
    print(io, "## D=$D, $spname, r=$r  (expected rates: L² $(getlv(res, "expected_L2_rate")), curl $r)\n\n")
    print(io, "| n | ndofs | ‖u−u_h‖_L² | rate | ‖d(u−u_h)‖ | rate | ‖Δu_h‖_L² | time (s) |\n")
    print(io, "|--:|------:|-----------:|:----:|-----------:|:----:|----------:|--------:|\n")
    rL2 = getlv(res, "ratesL2_sorted"); rC = getlv(res, "ratesCurl_sorted")
    for (i, lv) in enumerate(levels)
      print(io, "| ", getlv(lv, "n"), " | ", getlv(lv, "ndofs"),
            " | ", fnum(getlv(lv, "errL2_sorted")),
            " | ", i == 1 ? "—" : frate(rL2[i-1]),
            " | ", fnum(getlv(lv, "errCurl_sorted")),
            " | ", i == 1 ? "—" : frate(rC[i-1]),
            " | ", fnum(getlv(lv, "l2_norm_uh_diff")),
            " | ", round(getlv(lv, "walltime_s"), digits=1), " |\n")
    end
    print(io, "\n")
  end
  print(io, "Environment: Julia $(VERSION); machine recorded in `environment.md`. ",
            "Quadrature degree max(2r+2, 6).\n")
  open(joinpath(outdir, "summary_6$suffix.md"), "w") do f
    write(f, String(take!(io)))
  end
  println("Wrote ", joinpath(outdir, "exp6_curlcurl$suffix.json"))
  println("Wrote ", joinpath(outdir, "summary_6$suffix.md"))
end

main()

# Experiment 3 — H(curl) conformity on scrambled simplicial meshes.
#
# Run from the repository root:  julia --project=. exp3_conformity.jl
#
# (a) Two-triangle / two-tet meshes exactly as in
#     test/FESpacesTests/RotatingPLambdaTests/PΛFESpaceTests.jl: for every
#     relative cell vertex ordering, the max tangential jump across the shared
#     facet over all global basis functions, with the rotation change of basis
#     ON (the FESpace) and OFF (negative control: raw covariant Piola
#     pushforward via get_cell_shapefuns_and_dof_basis(..., nothing, ...)).
#     Reported: max over orderings of both numbers, per (space, D, r).
# (b) Larger meshes: simplexified CartesianDiscreteModel (8×8 in 2D, 4×4×4 in
#     3D), rebuilt as an UnstructuredDiscreteModel in which every cell's
#     vertex ordering is permuted by a pseudo-random permutation drawn from a
#     single MersenneTwister(42) stream (one randperm per cell, in cell
#     order). Measured: max tangential jump over ALL interior facets and ALL
#     global basis functions, rotation ON vs OFF; assembled mass-matrix SPD
#     check (smallest eigenvalue); n_cells, n_dofs.
#     Grid: 2D r ∈ {1,2,3}, 3D r ∈ {1,2}, both spaces.
#
# Output: results/exp3_conformity.json and exp3_table.md
# Deterministic: the only randomness is the fixed-seed cell scrambling.

using Gridap
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap.Fields
using Gridap.Arrays: Table, Fill
using Combinatorics: permutations
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

# ── Two-cell mesh machinery (verbatim from test/…/PΛFESpaceTests.jl) ─────────

const X2D = [[0.0,0.0], [1.0,0.0], [0.0,1.0], [1.0,1.0]]
const X3D = [[0.0,0.0,0.0], [1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0], [1.0,1.0,1.0]]

function two_cell_model(D::Int, c1::Vector{Int}, c2::Vector{Int})
  X = D == 2 ? X2D : X3D
  coords = [Point(x...) for x in X]
  reffe  = LagrangianRefFE(Float64, D == 2 ? TRI : TET, 1)
  grid   = UnstructuredGrid(coords, Table([c1, c2]), [reffe], fill(Int8(1), 2))
  UnstructuredDiscreteModel(grid)
end

function cell_affine(model, c)
  X    = get_node_coordinates(get_grid(model))
  conn = get_cell_node_ids(get_grid(model))
  g = conn[c]
  D = length(first(X))
  o = collect(X[g[1]].data)
  J = hcat([collect(X[g[j+1]].data) .- o for j in 1:D]...)
  o, J
end

function interface_points_tangents(D::Int)
  X = D == 2 ? X2D : X3D
  ts = [X[2+j] .- X[2] for j in 1:D-1]
  wts = D == 2 ? [(0.15,), (0.5,), (0.83,)] :
                 [(0.2, 0.3), (0.5, 0.25), (0.1, 0.7), (1/3, 1/3)]
  pts = [X[2] .+ sum(w[j] .* ts[j] for j in 1:D-1) for w in wts]
  pts, ts
end

function interface_jump(model, cell_fields, cell_ids, D)
  pts, ts = interface_points_tangents(D)
  d = [collect(cell_ids[1]), collect(cell_ids[2])]
  shared = intersect(Set(d[1]), Set(d[2]))
  geo = [cell_affine(model, c) for c in 1:2]
  maxjump = 0.0
  for p in pts
    tr = Dict{Int,Vector{Vector{Float64}}}()
    for c in 1:2
      o, J = geo[c]
      x̂ = Point((J \ (p .- o))...)
      vals = evaluate(cell_fields[c], [x̂])
      for (l, gdof) in enumerate(d[c])
        v = collect(vals[1, l].data)
        tv = [dot(v, t) for t in ts]
        if gdof in shared
          push!(get!(tr, gdof, Vector{Float64}[]), tv)
        else
          maxjump = max(maxjump, maximum(abs, tv))
        end
      end
    end
    for (_, vs) in tr
      maxjump = max(maxjump, maximum(abs, vs[1] .- vs[2]))
    end
  end
  maxjump
end

# Vector-proxied P_rΛ¹ (rotating) and P_r⁻Λ¹ (trimmed) H(curl) elements; the
# default flavor :BMM gives the bases of the paper (rotation entries exactly ±1).
simplex(D) = D == 2 ? TRI : TET
make_reffe(name, D, r) = name == :rotating ? RotatingPΛRefFE(Float64, simplex(D), r) :
                                             TrimmedPΛRefFE(Float64, simplex(D), r)

# Rotation-ON cell basis fields (the FESpace's own) and rotation-OFF fields
# (raw covariant Piola pushforward, cell_changes = nothing — the negative
# control bypass pinned by the test file).
on_fields(V) = Gridap.CellData.get_data(get_fe_basis(V))
function off_fields(model, rf, ncells)
  cell_reffe = Fill(rf, ncells)
  cell_Jt = lazy_map(Broadcasting(∇), get_cell_map(get_grid(model)))
  cell_fields_raw, _ = FESpaces.get_cell_shapefuns_and_dof_basis(
    CoVariantPiolaMap(), model, cell_reffe, nothing, cell_Jt)
  cell_fields_raw
end

# ── Part (a): all two-cell relative orderings ─────────────────────────────────

const TRI1_ORDERINGS = ([1,2,3], [3,1,2], [2,1,3])
const TET_ORDERINGS = (
  ([1,2,3,4], [2,3,4,5]),   # both sorted (fast path, no change of basis)
  ([1,2,3,4], [3,2,4,5]),   # non-cyclic permutation of the shared face
  ([1,2,3,4], [4,3,2,5]),
  ([1,2,3,4], [5,4,3,2]),
  ([2,4,1,3], [5,4,3,2]),
  ([4,1,3,2], [3,2,5,4]),
)

two_cell_orderings(D) = D == 2 ?
  [(c1, collect(c2)) for c1 in TRI1_ORDERINGS for c2 in permutations([2,3,4])] :
  collect(TET_ORDERINGS)

function part_a()
  records = Vector{Any}()
  for name in (:rotating, :trimmed), (D, rs) in ((2, (1,2,3)), (3, (1,2))), r in rs
    t0 = time()
    rf = make_reffe(name, D, r)
    ords = two_cell_orderings(D)
    jon = 0.0; joff = 0.0
    for (c1, c2) in ords
      model = two_cell_model(D, c1, c2)
      V = FESpace(model, rf, conformity=:HCurl)
      ids = get_cell_dof_ids(V)
      jon  = max(jon,  interface_jump(model, on_fields(V), ids, D))
      joff = max(joff, interface_jump(model, off_fields(model, rf, 2), ids, D))
    end
    push!(records, Pair{String,Any}[
      "space" => String(name), "D" => D, "r" => r,
      "n_orderings" => length(ords),
      "max_jump_rotation_on" => jon,
      "max_jump_rotation_off" => joff,
      "walltime_s" => round(time() - t0, digits=2),
    ])
    @printf("(a) %s D=%d r=%d: %d orderings, jump ON=%.3e OFF=%.3e (%.1fs)\n",
      name, D, r, length(ords), jon, joff, time() - t0)
  end
  records
end

# ── Part (b): larger scrambled meshes ─────────────────────────────────────────

# Scramble every cell's vertex ordering with one MersenneTwister(42) stream
# (one randperm per cell, drawn in cell order — fixed and reproducible).
function scrambled_simplex_model(D::Int; seed=42)
  model0 = D == 2 ?
    simplexify(CartesianDiscreteModel((0.0,1.0,0.0,1.0), (8,8))) :
    simplexify(CartesianDiscreteModel((0.0,1.0,0.0,1.0,0.0,1.0), (4,4,4)))
  grid0  = get_grid(model0)
  coords = collect(get_node_coordinates(grid0))
  conn   = get_cell_node_ids(grid0)
  rng    = MersenneTwister(seed)
  rows   = [collect(conn[c])[randperm(rng, D + 1)] for c in 1:length(conn)]
  reffe  = LagrangianRefFE(Float64, D == 2 ? TRI : TET, 1)
  grid   = UnstructuredGrid(coords, Table(rows), [reffe], fill(Int8(1), length(rows)))
  UnstructuredDiscreteModel(grid)
end

# Max tangential jump over ALL interior facets and all global basis functions.
# Same trace logic as interface_jump, generalized to every 2-cell facet:
# shared dofs must have matching tangential traces from both cells; dofs of
# one cell only must have vanishing tangential trace on the facet.
function all_facet_jump(model, cell_fields, cell_ids, D)
  grid = get_grid(model)
  X    = get_node_coordinates(grid)
  topo = get_grid_topology(model)
  f2c  = get_faces(topo, D - 1, D)
  f2v  = get_faces(topo, D - 1, 0)
  ncells = num_cells(grid)
  geo  = [cell_affine(model, c) for c in 1:ncells]
  wts  = D == 2 ? [(0.15,), (0.5,), (0.83,)] :
                  [(0.2, 0.3), (0.5, 0.25), (0.1, 0.7), (1/3, 1/3)]
  maxjump = 0.0
  for f in 1:length(f2c)
    cells = f2c[f]
    length(cells) == 2 || continue
    vs = f2v[f]
    v1 = collect(X[vs[1]].data)
    ts = [collect(X[vs[1+j]].data) .- v1 for j in 1:D-1]
    pts = [v1 .+ sum(w[j] .* ts[j] for j in 1:D-1) for w in wts]
    d = [collect(cell_ids[c]) for c in cells]
    shared = intersect(Set(d[1]), Set(d[2]))
    vals = Vector{Any}(undef, 2)
    for ci in 1:2
      o, J = geo[cells[ci]]
      x̂s = [Point((J \ (p .- o))...) for p in pts]
      vals[ci] = evaluate(cell_fields[cells[ci]], x̂s)
    end
    for (ip, _) in enumerate(pts)
      tr = Dict{Int,Vector{Vector{Float64}}}()
      for ci in 1:2
        for (l, gdof) in enumerate(d[ci])
          v = collect(vals[ci][ip, l].data)
          tv = [dot(v, t) for t in ts]
          if gdof in shared
            push!(get!(tr, gdof, Vector{Float64}[]), tv)
          else
            maxjump = max(maxjump, maximum(abs, tv))
          end
        end
      end
      for (_, vv) in tr
        maxjump = max(maxjump, maximum(abs, vv[1] .- vv[2]))
      end
    end
  end
  maxjump
end

function mass_matrix_min_eig(model, V, r)
  trian = Triangulation(model)
  dΩ    = Measure(trian, 2r + 2)
  U     = TrialFESpace(V)
  a_mass(u, v) = ∫(u ⋅ v) * dΩ           # = ∫ vol_coeff(u ∧ ⋆v) for the 1-forms
  M  = Matrix(assemble_matrix(a_mass, U, V))
  asym = maximum(abs, M - M')
  λmin = eigvals(Symmetric((M + M') / 2), 1:1)[1]
  λmin, asym
end

function part_b()
  records = Vector{Any}()
  models = Dict(D => scrambled_simplex_model(D) for D in (2, 3))
  for name in (:rotating, :trimmed), (D, rs) in ((2, (1,2,3)), (3, (1,2))), r in rs
    t0 = time()
    model = models[D]
    ncells = num_cells(get_grid(model))
    rf  = make_reffe(name, D, r)
    V   = FESpace(model, rf, conformity=:HCurl)
    ids = get_cell_dof_ids(V)
    n_int = count(f -> length(f) == 2, get_faces(get_grid_topology(model), D - 1, D))
    jon  = all_facet_jump(model, on_fields(V), ids, D)
    joff = all_facet_jump(model, off_fields(model, rf, ncells), ids, D)
    λmin, asym = mass_matrix_min_eig(model, V, r)
    push!(records, Pair{String,Any}[
      "space" => String(name), "D" => D, "r" => r,
      "mesh" => D == 2 ? "simplexify(8x8), scrambled (seed 42)" :
                         "simplexify(4x4x4), scrambled (seed 42)",
      "n_cells" => ncells,
      "n_interior_facets" => n_int,
      "n_dofs" => num_free_dofs(V),
      "max_jump_rotation_on" => jon,
      "max_jump_rotation_off" => joff,
      "mass_min_eigenvalue" => λmin,
      "mass_spd" => λmin > 0,
      "mass_max_asymmetry" => asym,
      "walltime_s" => round(time() - t0, digits=2),
    ])
    @printf("(b) %s D=%d r=%d: %d cells, %d dofs, jump ON=%.3e OFF=%.3e, λmin=%.3e (%.1fs)\n",
      name, D, r, ncells, num_free_dofs(V), jon, joff, λmin, time() - t0)
  end
  records
end

# ── Run and write ─────────────────────────────────────────────────────────────

rec_a = part_a()
rec_b = part_b()

out = Pair{String,Any}[
  "experiment" => "exp3_conformity",
  "description" => "Tangential-jump conformity on scrambled simplicial meshes, rotation ON " *
                   "vs OFF (raw covariant Piola, no change of basis); mass-matrix SPD check.",
  "two_cell" => rec_a,
  "scrambled_mesh" => rec_b,
]

resdir = joinpath(@__DIR__, "results")
mkpath(resdir)
open(joinpath(resdir, "exp3_conformity.json"), "w") do io
  _json(io, out); println(io)
end

fmt(x) = x == 0.0 ? "0" : @sprintf("%.2e", x)

open(joinpath(resdir, "exp3_table.md"), "w") do io
  println(io, "### Experiment 3 — conformity on scrambled meshes\n")
  println(io, "#### (a) Two-cell meshes, all relative vertex orderings\n")
  println(io, "Max tangential jump across the shared facet over all global basis functions ",
    "and all facet sample points; \"max over orderings\" of that number, rotation ON ",
    "(the FESpace) vs OFF (raw covariant Piola pushforward with no change of basis — ",
    "negative control). 2D: 18 orderings (3 for cell 1 × 6 for cell 2); 3D: the 6 ",
    "scrambles of the test suite, incl. a non-cyclic shared-face permutation.\n")
  println(io, "| space | D | r | orderings | max jump, rotation ON | max jump, rotation OFF |")
  println(io, "|---|---|---|---|---|---|")
  for rec in rec_a
    d = Dict(rec)
    println(io, "| ", d["space"], " | ", d["D"], " | ", d["r"], " | ", d["n_orderings"],
      " | ", fmt(d["max_jump_rotation_on"]), " | ", fmt(d["max_jump_rotation_off"]), " |")
  end
  println(io, "\n#### (b) Larger scrambled meshes\n")
  println(io, "Simplexified Cartesian meshes (2D: 8×8, 3D: 4×4×4) with every cell's vertex ",
    "ordering permuted by a fixed pseudo-random permutation (MersenneTwister seed 42, one ",
    "randperm per cell). Jump = max over ALL interior facets, all global basis functions, ",
    "all facet sample points. λ_min = smallest eigenvalue of the assembled mass matrix ",
    "(rotation ON, quadrature degree 2r+2).\n")
  println(io, "| space | D | r | cells | interior facets | dofs | max jump ON | max jump OFF | mass λ_min | SPD |")
  println(io, "|---|---|---|---|---|---|---|---|---|---|")
  for rec in rec_b
    d = Dict(rec)
    println(io, "| ", d["space"], " | ", d["D"], " | ", d["r"], " | ", d["n_cells"],
      " | ", d["n_interior_facets"], " | ", d["n_dofs"],
      " | ", fmt(d["max_jump_rotation_on"]), " | ", fmt(d["max_jump_rotation_off"]),
      " | ", fmt(d["mass_min_eigenvalue"]), " | ", d["mass_spd"] ? "yes" : "NO", " |")
  end
end

println("wrote ", joinpath(resdir, "exp3_conformity.json"))

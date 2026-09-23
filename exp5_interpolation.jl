# Experiment 5 — convergence of the coefficient-functional (collocation)
# interpolant under h-refinement, on sorted and scrambled meshes.
#
# The DOFs are the duals of the rotating bases, represented by pointwise
# predofs at BB lattice nodes (paper §4.5).  They act only in interpolation,
# so this experiment is the numerical check that the interpolation operator
# converges: I_h u for the smooth manufactured 1-forms of exp4, errors
# ‖u − I_h u‖_L² and ‖curl(u − I_h u)‖_L², rates vs the Bramble–Hilbert
# prediction (L²: r+1 full, r trimmed).
#
# Run from the repository root:  julia --project=. exp5_interpolation.jl
# Outputs: results/exp5_interpolation.json, results/summary_5.md

using Gridap
using Gridap.Geometry
using Gridap.ReferenceFEs
using Gridap.FESpaces
using Gridap.CellData
using Gridap.Fields
using Gridap.Arrays: Table
using LinearAlgebra

const SEED = 42

# Minimal JSON writer (no external deps), as in exp4
json(x::AbstractString) = "\"" * replace(x, "\\" => "\\\\", "\"" => "\\\"") * "\""
json(x::Bool)    = x ? "true" : "false"
json(x::Integer) = string(x)
json(x::Real)    = isfinite(x) ? string(Float64(x)) : "null"
json(d::AbstractDict) = "{" * join([json(String(k)) * ":" * json(v) for (k,v) in d], ",") * "}"
json(v::AbstractVector) = "[" * join(json.(v), ",") * "]"
json(x, n) = json(x)


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

# Vector-proxied P_rΛ¹ (full) and P_r⁻Λ¹ (trimmed) H(curl) elements; the
# default flavor :BMM gives the bases of the paper.
simplex(D) = D == 2 ? TRI : TET
make_reffe(name, D, r) = name == :full ? RotatingPΛRefFE(Float64, simplex(D), r) :
                                         TrimmedPΛRefFE(Float64, simplex(D), r)

"Interpolate u_fun onto FESpace(model, rf, :HCurl); return (L², curl, ndofs)."
function interp_errors(model, rf, u_fun, du_fun; deg)
  Ω  = Triangulation(model)
  dΩ = Measure(Ω, deg)
  V  = FESpace(model, rf, conformity=:HCurl)
  uh = interpolate(u_fun, V)
  u_cf  = CellField(u_fun, Ω, PhysicalDomain())
  du_cf = CellField(du_fun, Ω, PhysicalDomain())
  e  = uh - u_cf
  errL2 = sqrt(abs(sum(∫( e ⋅ e ) * dΩ)))
  ec = curl(uh) - du_cf
  errCurl = sqrt(abs(sum(∫( ec ⋅ ec ) * dΩ)))
  errL2, errCurl, num_free_dofs(V)
end

function run_case(D, name, r, ns)
  u_fun  = D == 2 ? u_ex_2d  : u_ex_3d
  du_fun = D == 2 ? du_ex_2d : du_ex_3d
  deg = max(2r + 2, 6)
  rf  = make_reffe(name, D, r)
  levels = Vector{Dict{String,Any}}()
  for n in ns
    coords, conn = base_mesh(D, n)
    conn_s = [sort(row) for row in conn]
    conn_x = [row[cell_perm(c, length(row))] for (c, row) in enumerate(conn_s)]
    eL2_s, eC_s, nd = interp_errors(build_model(D, coords, conn_s), rf, u_fun, du_fun; deg=deg)
    eL2_x, eC_x, _  = interp_errors(build_model(D, coords, conn_x), rf, u_fun, du_fun; deg=deg)
    push!(levels, Dict("n"=>n, "h"=>1.0/n, "ndofs"=>nd,
      "errL2_sorted"=>eL2_s, "errL2_scrambled"=>eL2_x,
      "errCurl_sorted"=>eC_s, "errCurl_scrambled"=>eC_x))
    @info "done" D name r n eL2_s eL2_x
  end
  levels
end

rate(errs, hs) = [log(errs[i-1]/errs[i]) / log(hs[i-1]/hs[i]) for i in 2:length(errs)]

function main()
  results = Dict{String,Any}("seed"=>SEED, "cases"=>Dict{String,Any}())
  cases = [(2,:full,1,[4,8,16,32]), (2,:full,2,[4,8,16,32]), (2,:full,3,[4,8,16,32]),
           (2,:trimmed,1,[4,8,16,32]), (2,:trimmed,2,[4,8,16,32]), (2,:trimmed,3,[4,8,16,32]),
           (3,:full,1,[2,4,8]), (3,:full,2,[2,4,8]),
           (3,:trimmed,1,[2,4,8]), (3,:trimmed,2,[2,4,8])]
  for (D,name,r,ns) in cases
    key = "D$(D)_$(name)_r$(r)"
    levels = run_case(D, name, r, ns)
    hs   = [l["h"] for l in levels]
    eS   = [l["errL2_scrambled"] for l in levels]
    eC   = [l["errCurl_scrambled"] for l in levels]
    results["cases"][key] = Dict("levels"=>levels,
      "ratesL2"=>rate(eS,hs), "ratesCurl"=>rate(eC,hs),
      "expectedL2"=> name==:full ? r+1 : r)
    println("== $key  L2 rates: ", round.(rate(eS,hs), digits=2),
            "  curl rates: ", round.(rate(eC,hs), digits=2),
            "  expected L2: ", name==:full ? r+1 : r)
  end
  outdir = joinpath(@__DIR__, "results")
  mkpath(outdir)
  open(joinpath(outdir, "exp5_interpolation.json"),"w") do io; print(io, json(results, 2)); end
  # summary table
  open(joinpath(outdir, "summary_5.md"),"w") do io
    println(io, "| space | D | r | errL2 (finest, scrambled) | L2 rate | expected | curl rate | sorted-vs-scr rel. diff |")
    for (D,name,r,ns) in cases
      key = "D$(D)_$(name)_r$(r)"; c = results["cases"][key]; lv = c["levels"][end]
      rel = abs(lv["errL2_sorted"]-lv["errL2_scrambled"])/lv["errL2_scrambled"]
      println(io, "| $name | $D | $r | $(lv["errL2_scrambled"]) | $(round(c["ratesL2"][end],digits=2)) | $(c["expectedL2"]) | $(round(c["ratesCurl"][end],digits=2)) | $(round(rel,sigdigits=2)) |")
    end
  end
  println("wrote ", joinpath(outdir, "exp5_interpolation.json"), " and summary_5.md")
end

main()

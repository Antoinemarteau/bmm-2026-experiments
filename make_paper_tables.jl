# Post-processing: the data of the article's tables, as CSV.
#
# Reads the committed JSON in results/ and writes one CSV per table of the
# article that displays experiment data:
#
#   results/paper/table2_sparsity.csv     Table 2, Sect. 8.2  <- exp2_sparsity.json
#   results/paper/table3_conformity.csv   Table 3, Sect. 8.3  <- exp3_conformity.json
#   results/paper/table4_convergence.csv  Table 4, Sect. 8.4  <- exp4_convergence.json
#   results/paper/table5_curlcurl.csv     Table 5, Sect. 8.5  <- exp6_curlcurl.json
#
# Each file holds the rows, the columns and the rounding of the printed table,
# so that a cell of the article can be checked against the data without
# reading the JSON. Full precision stays in the JSON; the rounding rules are
# listed per column in the README.
#
# No experiment is re-run: the committed results are the only input. Table 1
# of the article compares strategies from the literature and carries no
# experiment data, and the article has no data-generated figure, so these four
# files cover every experimental number of the article except those quoted in
# the prose of Sect. 8.
#
# Usage: julia --project=. make_paper_tables.jl

using Printf

const RESDIR = joinpath(@__DIR__, "results")
const OUTDIR = joinpath(RESDIR, "paper")

# ── Minimal JSON reader (no external deps) ────────────────────────────────────
mutable struct _JsonIn
  s::String
  i::Int
end

_curr(r::_JsonIn) = r.s[r.i]

function _skipws!(r::_JsonIn)
  while r.i <= lastindex(r.s) && isspace(r.s[r.i])
    r.i = nextind(r.s, r.i)
  end
end

function _value(r::_JsonIn)
  _skipws!(r)
  c = _curr(r)
  c == '{' && return _object(r)
  c == '[' && return _array(r)
  c == '"' && return _string(r)
  c == 't' && (r.i += 4; return true)
  c == 'f' && (r.i += 5; return false)
  c == 'n' && (r.i += 4; return nothing)
  return _number(r)
end

function _object(r::_JsonIn)
  d = Dict{String,Any}()
  r.i += 1                                  # '{'
  _skipws!(r)
  _curr(r) == '}' && (r.i += 1; return d)
  while true
    _skipws!(r)
    k = _string(r)
    _skipws!(r); r.i += 1                   # ':'
    d[k] = _value(r)
    _skipws!(r)
    c = _curr(r); r.i += 1                  # ',' or '}'
    c == '}' && return d
  end
end

function _array(r::_JsonIn)
  v = Any[]
  r.i += 1                                  # '['
  _skipws!(r)
  _curr(r) == ']' && (r.i += 1; return v)
  while true
    push!(v, _value(r))
    _skipws!(r)
    c = _curr(r); r.i += 1                  # ',' or ']'
    c == ']' && return v
  end
end

function _string(r::_JsonIn)
  io = IOBuffer()
  r.i += 1                                  # opening quote
  while true
    c = r.s[r.i]; r.i = nextind(r.s, r.i)
    if c == '"'
      return String(take!(io))
    elseif c == '\\'
      e = r.s[r.i]; r.i = nextind(r.s, r.i)
      print(io, e == 'n' ? '\n' : e == 't' ? '\t' : e == 'r' ? '\r' : e)
    else
      print(io, c)
    end
  end
end

function _number(r::_JsonIn)
  j = r.i
  while r.i <= lastindex(r.s) && (isdigit(r.s[r.i]) || r.s[r.i] in "+-.eE")
    r.i = nextind(r.s, r.i)
  end
  t = r.s[j:prevind(r.s, r.i)]
  return any(c -> c in ".eE", t) ? parse(Float64, t) : parse(Int, t)
end

readjson(path) = _value(_JsonIn(read(path, String), 1))

# ── The column formats of the article ────────────────────────────────────────
dec1(x) = @sprintf("%.1f", x)
dec2(x) = @sprintf("%.2f", x)
dec3(x) = @sprintf("%.3f", x)
sci2(x) = @sprintf("%.1e", x)               # 2 significant digits: 2.8e-17
sci3(x) = @sprintf("%.2e", x)               # 3 significant digits: 8.97e-03
time_us(x) = x < 1 ? dec2(x) : dec1(x)      # microseconds, 2 significant digits
# The O(1) rotation-off jumps carry one digit more than the article prints:
# one of them is 0.745 up to float noise, which two digits would round to 0.74
# on one mesh and to 0.75 on the other.
jump_off(x) = dec3(x)

# Experiments 2 and 3 name the full space "rotating", 4 and 6 name it "full",
# which is the name the article uses.
paper_space(s) = s == "rotating" ? "full" : s

# The article's expected rates: r + 1 in L2 for the full space, r for the
# trimmed one, r in the curl seminorm for both.
expected_l2(space, r) = space == "full" ? r + 1 : r

# ── Table 2: sparsity and cost of the rotation maps ──────────────────────────
function table2(path)
  recs = readjson(joinpath(RESDIR, "exp2_sparsity.json"))["records"]
  open(path, "w") do io
    println(io, "space,D,r,n,multi_term_rows_mean_per_perm,nnz_per_row_mean,",
                "nnz_per_row_max,nnz_vs_n2_percent,map_cold_us,map_hot_ns")
    for d in recs
      # The article omits r = 1, whose maps are index permutations, and r = 2
      # in three dimensions.
      d["r"] >= 2 || continue
      (d["D"] == 3 && d["r"] == 2) && continue
      println(io, join([paper_space(d["space"]), d["D"], d["r"], d["n"],
                        dec2(d["hit_rows_mean_per_perm"]),
                        dec2(d["nnz_mean_per_row"]),
                        d["nnz_max_per_row"],
                        dec1(100 * d["density_mean"]),
                        time_us(1e6 * d["t_rotation_map_cold_s"]),
                        round(Int, 1e9 * d["t_rotation_map_hot_s"])], ","))
    end
  end
end

# ── Table 3: conformity on scrambled meshes ──────────────────────────────────
# Each row of the article pairs a two-cell case with the scrambled-mesh case of
# the same space and degree, printed as "two-cell / mesh"; the DOF count is the
# one of the mesh.
function table3(path)
  j  = readjson(joinpath(RESDIR, "exp3_conformity.json"))
  tc, sm = j["two_cell"], j["scrambled_mesh"]
  open(path, "w") do io
    println(io, "space,D,r,cells_two_cell,cells_mesh,dofs_mesh,",
                "jump_on_two_cell,jump_on_mesh,jump_off_two_cell,jump_off_mesh")
    for (a, b) in zip(tc, sm)
      key(x) = (x["space"], x["D"], x["r"])
      key(a) == key(b) || error("exp3 records out of step: $(key(a)) vs $(key(b))")
      println(io, join([paper_space(a["space"]), a["D"], a["r"],
                        2, b["n_cells"], b["n_dofs"],
                        sci2(a["max_jump_rotation_on"]),
                        sci2(b["max_jump_rotation_on"]),
                        jump_off(a["max_jump_rotation_off"]),
                        jump_off(b["max_jump_rotation_off"])], ","))
    end
  end
end

# ── Table 4: L2 projection on scrambled meshes ───────────────────────────────
# The article tabulates the two-dimensional cases at 8, 16 and 32 cells per
# side, the rate between the two finest of those, and the largest
# sorted-versus-scrambled difference over all levels, the 4-cell one included.
# The three-dimensional cases appear in the prose of Sect. 8.4 instead.
function table4(path)
  cases = readjson(joinpath(RESDIR, "exp4_convergence.json"))["cases"]
  open(path, "w") do io
    println(io, "space,r,err_L2_n8,err_L2_n16,err_L2_n32,rate,expected_rate,",
                "sorted_vs_scrambled")
    for c in cases
      c["D"] == 2 || continue
      lv = Dict(l["n"] => l for l in c["levels"])
      e8, e16, e32 = (lv[n]["errL2_scrambled"] for n in (8, 16, 32))
      space = paper_space(c["space"])
      println(io, join([space, c["r"], sci3(e8), sci3(e16), sci3(e32),
                        dec1(log2(e16 / e32)),
                        expected_l2(space, c["r"]),
                        sci2(maximum(l["l2_norm_uh_diff"] for l in c["levels"]))], ","))
    end
  end
end

# ── Table 5: curl-curl problem on scrambled meshes ───────────────────────────
# One row per case, at its finest level, with both rates taken between the two
# finest levels and the largest sorted-versus-scrambled difference over all of
# them. The article prints the two-dimensional rates to one decimal and the
# pre-asymptotic three-dimensional ones to two.
function table5(path)
  cases = readjson(joinpath(RESDIR, "exp6_curlcurl.json"))["cases"]
  open(path, "w") do io
    println(io, "space,D,r,cells_per_side,dofs,err_L2,rate_L2,err_curl,rate_curl,",
                "expected_rate_L2,expected_rate_curl,sorted_vs_scrambled")
    for c in cases
      lv, space = c["levels"], paper_space(c["space"])
      fine, prev = lv[end], lv[end-1]
      rate = c["D"] == 2 ? dec1 : dec2
      println(io, join([space, c["D"], c["r"], fine["n"], fine["ndofs"],
                        sci3(fine["errL2_scrambled"]),
                        rate(log2(prev["errL2_scrambled"] / fine["errL2_scrambled"])),
                        sci3(fine["errCurl_scrambled"]),
                        rate(log2(prev["errCurl_scrambled"] / fine["errCurl_scrambled"])),
                        expected_l2(space, c["r"]), c["r"],
                        sci2(maximum(l["l2_norm_uh_diff"] for l in lv))], ","))
    end
  end
end

mkpath(OUTDIR)
for (name, build) in ("table2_sparsity.csv"    => table2,
                      "table3_conformity.csv"  => table3,
                      "table4_convergence.csv" => table4,
                      "table5_curlcurl.csv"    => table5)
  path = joinpath(OUTDIR, name)
  build(path)
  println("wrote ", path)
end

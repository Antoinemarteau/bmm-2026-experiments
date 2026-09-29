# Assemble results/summary_1_3.md from the per-experiment markdown
# fragments (exp1_table.md, exp2_table.md, exp3_table.md), which are written by
# the exp{1,2,3}_*.jl scripts. Run after those scripts, from the repo root:
#
#   julia --project=. make_summary.jl

using Pkg

resdir = joinpath(@__DIR__, "results")

# Read the environment back from the active project instead of naming versions
# here: a literal version in this file outlives the run that produced the
# tables and then misreports it.
function dep_version(name::AbstractString)
  try
    for (_, p) in Pkg.dependencies()
      p.name == name && return string(p.version)
    end
  catch
  end
  "unknown"
end

const GRIDAP_VERSION = dep_version("Gridap")

header = """
# Numerical experiments 1–3 — quantitative results

Package: Gridap.jl $(GRIDAP_VERSION), with the vector-proxied rotating P_rΛ¹ and
trimmed P_r⁻Λ¹ bases and their rotation change-of-basis calculus. Environment:
Julia $(VERSION), Gridap $(GRIDAP_VERSION); the machine is recorded in
`results/environment.md`. All scripts are deterministic; the only
pseudo-randomness is `MersenneTwister(42)` (evaluation points in Exp. 1,
mat-vec input in Exp. 2, cell scrambling in Exp. 3b).

Reproduce from the repo root:

```
julia --project=. exp1_exactness.jl
julia --project=. exp2_sparsity.jl
julia --project=. exp3_conformity.jl
```

Raw data: `results/exp{1_exactness,2_sparsity,3_conformity}.json`.

Notation: `n` = dimension of the local basis; C(π) = `rotation_change_of_basis`
matrix of the vertex relabelling π ∈ S_{D+1}; "rotating" = untrimmed P_rΛ¹,
"trimmed" = P_r⁻Λ¹.

---

"""

open(joinpath(resdir, "summary_1_3.md"), "w") do io
  print(io, header)
  for (i, frag) in enumerate(("exp1_table.md", "exp2_table.md", "exp3_table.md"))
    i > 1 && print(io, "\n---\n\n")
    print(io, read(joinpath(resdir, frag), String))
  end
end

println("wrote ", joinpath(resdir, "summary_1_3.md"))

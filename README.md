[![doi](https://zenodo.org/badge/doi/10.5281/zenodo.23049283.svg)](https://doi.org/10.5281/zenodo.23049283)

# bmm-2026-experiments

Experiments for the article "Rotating Bases for Finite Element Exterior
Calculus: Closed-Form Change of Basis Under Vertex Permutations", by Santiago
Badia, Jordi Manyer and Antoine Marteau.

Each script writes its raw data as JSON and a Markdown table into `results/`;
re-running overwrites them. `results/environment.md` records the machine and
software environment the committed results were produced on.

## Setup

The repository contains the `Project.toml` and `Manifest.toml` specifying the
code version the results were generated on, using the open source Gridap.jl
v0.20.10 available from the official Julia packages registries. Instantiating
reproduces that environment:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

To run against a local Gridap checkout instead, `Pkg.develop(path="/path/to/Gridap.jl")`.

## Run

Everything, in order, with the environment recorded and a log per experiment:

```sh
./run_all.sh              # --quick for a smoke test, --only 2,6 for a subset, --help for the rest
```

Or one at a time from the repository root, in any order except
`make_summary.jl`, which concatenates the tables of Experiments 1–3 and so runs
last:

```sh
julia --project=. env_stamp.jl            # writes results/environment.{md,json}
julia --project=. exp1_exactness.jl
julia --project=. exp2_sparsity.jl
julia --project=. exp3_conformity.jl
julia --project=. exp4_convergence.jl     # the 3D cases take a while; EXP4_QUICK=1 for a smoke test
julia --project=. exp5_interpolation.jl
julia --project=. exp6_curlcurl.jl        # idem; EXP6_QUICK=1 for a smoke test
julia --project=. make_summary.jl
```

## The experiments

| script | what it measures |
|---|---|
| `exp1_exactness.jl` | entries of $T(\pi)$ exactly in $\{-1,0,1\}$, $T(\pi)T(\pi^{-1})=I$, composition law, pointwise pullback residual; all $\pi,\tau$, $D\in\{2,3\}$, $r\le 4$ |
| `exp2_sparsity.jl` | hit rows, nnz per row, timings of the rotation map (cold/hot) vs dense build |
| `exp3_conformity.jl` | tangential jumps on two-cell and scrambled meshes, rotation on/off, mass SPD |
| `exp4_convergence.jl` | $L^2$ projection on sorted vs scrambled meshes, rates, sorted-vs-scrambled difference |
| `exp5_interpolation.jl` | convergence of the coefficient-functional interpolant, which the paper does not use |
| `exp6_curlcurl.jl` | curl–curl problem $\mathrm{curl}\,\mathrm{curl}\,u+u=f$, natural boundary condition, sorted vs scrambled meshes, rates and sorted-vs-scrambled difference |
| `make_summary.jl` | concatenates the Exp. 1–3 tables into `results/summary_1_3.md` |
| `make_paper_tables.jl` | exports the data behind the paper, see [Post-processing](#post-processing-the-tables-of-the-paper) |
| `env_stamp.jl` | records CPU, Julia, Gridap and BLAS state in `results/environment.md` |
| `run_all.sh` | runs the above in order, serially, with logs and a pass/fail report |

The scripts run on Gridap's vector-proxied implementation of the bases: the
1-forms are `VectorValue`s, the local bases are `BarycentricPΛBasis` (full
$P_r\Lambda^1$) and `BarycentricPmΛBasis` (trimmed $P_r^-\Lambda^1$) with
`flavor = :BMM`, the elements are `RotatingPΛRefFE` / `TrimmedPΛRefFE`, and the
rotation calculus is `Gridap.Polynomials.rotation_change_of_basis` and friends.
Form-language quantities are computed through their vector proxies:
$\langle u, v\rangle = u\cdot v$ and $du \leftrightarrow \mathrm{curl}\,u$
(a scalar in 2D, $(du_{23}, -du_{13}, du_{12})$ in 3D).

The timings of Experiment 2 are very machine-dependent: `exp2_table.md` names
the machine it was measured on in its own header.

## The results files

`results/` holds the committed output of the runs recorded in
`results/environment.md`; re-running overwrites it.

| file | written by |
|---|---|
| `environment.{md,json}` | `env_stamp.jl` |
| `exp1_exactness.json`, `exp1_table.md` | `exp1_exactness.jl` |
| `exp2_sparsity.json`, `exp2_table.md` | `exp2_sparsity.jl` |
| `exp3_conformity.json`, `exp3_table.md` | `exp3_conformity.jl` |
| `exp4_convergence.json`, `summary_4.md` | `exp4_convergence.jl` |
| `exp5_interpolation.json`, `summary_5.md` | `exp5_interpolation.jl` |
| `exp6_curlcurl.json`, `summary_6.md` | `exp6_curlcurl.jl` |
| `summary_1_3.md` | `make_summary.jl` |
| `paper/table{2,3,4,5}_*.csv` | `make_paper_tables.jl`, see below |

## Post-processing: the tables of the paper

The paper displays experiment data in four tables, Tables 2 to 5 of Sects. 8.2
to 8.5. `make_paper_tables.jl` reads the committed JSON and writes the data of
each of those tables as a CSV, one line per line of the printed table:

```sh
julia --project=. make_paper_tables.jl    # results/*.json -> results/paper/*.csv
```

| paper | CSV | from | rows |
|---|---|---|---|
| Table 2, Sect. 8.2 | `results/paper/table2_sparsity.csv` | `exp2_sparsity.json` | $r \ge 2$, without $D = 3, r = 2$ |
| Table 3, Sect. 8.3 | `results/paper/table3_conformity.csv` | `exp3_conformity.json` | each two-cell record with the scrambled-mesh record of the same space and degree, printed `two-cell / mesh` |
| Table 4, Sect. 8.4 | `results/paper/table4_convergence.csv` | `exp4_convergence.json` | the $D = 2$ cases, at 8, 16 and 32 cells per side |
| Table 5, Sect. 8.5 | `results/paper/table5_curlcurl.csv` | `exp6_curlcurl.json` | every case, at its finest level |

The columns the paper derives rather than tabulates as such:

- `nnz vs n²` is `density_mean` as a percentage, and `map, cold` / `map, hot`
  are `t_rotation_map_cold_s` / `t_rotation_map_hot_s` in microseconds and in
  nanoseconds.
- A rate is the base-2 logarithm of the ratio of the errors of the two finest
  levels.
- `sorted vs scrambled` is the largest `l2_norm_uh_diff` over *all* levels of
  the case, the coarsest included, not the one of the tabulated level.
- `expected` is not measured: it is $r + 1$ in $L^2$ for the full space, $r$
  for the trimmed one, and $r$ in the curl seminorm for both.
- Experiments 2 and 3 name the full space `rotating`, as Gridap does; the
  paper names it `full`, and so do the CSVs.

Other numbers quoted in the text from Section 8 are also found in the same
json files.

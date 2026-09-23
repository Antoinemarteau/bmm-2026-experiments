# bmm-2026-experiments

Experiments for the article on rotating bases for exterior calculus.

Each script writes its raw data as JSON and a Markdown table into `results/`;
re-running overwrites them.

## Setup

The repository contains the `Project.toml` and `Manifest.toml` specifying the
code version the results were generated on, using Gridap.jl v0.20.10. They can
be reproduced by instantiating reproduces the exact environment:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

To run against a local Gridap checkout instead, `Pkg.develop(path="/path/to/Gridap.jl")`.

## Run

From the repository root, in any order except `make_summary.jl`, which
concatenates the tables of Experiments 1–3 and so runs last:

```sh
julia --project=. exp1_exactness.jl
julia --project=. exp2_sparsity.jl
julia --project=. exp3_conformity.jl
julia --project=. exp4_convergence.jl     # the 3D cases take a while; EXP4_QUICK=1 for a smoke test
julia --project=. exp5_interpolation.jl
julia --project=. exp6_curlcurl.jl        # idem; EXP6_QUICK=1 for a smoke test
julia --project=. make_summary.jl
```

## The experiments

| script | paper | what it measures |
|---|---|---|
| `exp1_exactness.jl` | Sect. 7.1 | entries of $T(\pi)$ exactly in $\{-1,0,1\}$, $T(\pi)T(\pi^{-1})=I$, composition law, pointwise pullback residual; all $\pi,\tau$, $D\in\{2,3\}$, $r\le 4$ |
| `exp2_sparsity.jl` | Sect. 7.2, Table 2 | hit rows, nnz per row, timings of the rotation map (cold/hot) vs dense build |
| `exp3_conformity.jl` | Sect. 7.3, Table 3 | tangential jumps on two-cell and scrambled meshes, rotation on/off, mass SPD |
| `exp4_convergence.jl` | Sect. 7.4, Table 4 | $L^2$ projection on sorted vs scrambled meshes, rates, sorted-vs-scrambled difference |
| `exp5_interpolation.jl` | (not in the paper) | convergence of the coefficient-functional interpolant |
| `exp6_curlcurl.jl` | Sect. 7.5 | curl–curl problem $\mathrm{curl}\,\mathrm{curl}\,u+u=f$, natural boundary condition, sorted vs scrambled meshes, rates and sorted-vs-scrambled difference |
| `make_summary.jl` | | concatenates the Exp. 1–3 tables into `results/summary_1_3.md` |

The scripts run on Gridap's vector-proxied implementation of the bases: the
1-forms are `VectorValue`s, the local bases are `BarycentricPΛBasis` (full
$P_r\Lambda^1$) and `BarycentricPmΛBasis` (trimmed $P_r^-\Lambda^1$) with
`flavor = :BMM`, the elements are `RotatingPΛRefFE` / `TrimmedPΛRefFE`, and the
rotation calculus is `Gridap.Polynomials.rotation_change_of_basis` and friends.
Form-language quantities are computed through their vector proxies:
$\langle u, v\rangle = u\cdot v$ and $du \leftrightarrow \mathrm{curl}\,u$
(a scalar in 2D, $(du_{23}, -du_{13}, du_{12})$ in 3D).

The timings in Table 2 are very machine-dependent.

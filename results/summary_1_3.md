# Numerical experiments 1–3 — quantitative results

Package: Gridap.jl, branch with the vector-proxied rotating P_rΛ¹ and trimmed
P_r⁻Λ¹ bases and their rotation change-of-basis calculus. Environment:
Julia 1.12.7, Gridap 0.20.9. All scripts are deterministic; the only pseudo-randomness is
`MersenneTwister(42)` (evaluation points in Exp. 1, mat-vec input in Exp. 2,
cell scrambling in Exp. 3b).

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

### Experiment 1 — exactness and representation property of C(π)

All checks quantify over every π (and every pair (π, τ)) in S_{D+1}. "exact" columns are verified in integer arithmetic after an exact Float64→Int roundtrip; residuals are Float64. Pullback residual: max over all π and all basis functions of the pointwise defect against the numeric pullback oracle at 20 fixed interior points.

| space | D | r | n | entries ∈ {−1,0,1} exact | inverse exact (resid) | composition exact (resid) | pullback residual |
|---|---|---|---|---|---|---|---|
| rotating | 2 | 1 | 6 | yes | yes (0) | yes (0) | 1.11e-16 |
| rotating | 2 | 2 | 12 | yes | yes (0) | yes (0) | 2.12e-16 |
| rotating | 2 | 3 | 20 | yes | yes (0) | yes (0) | 2.98e-16 |
| rotating | 2 | 4 | 30 | yes | yes (0) | yes (0) | 3.89e-16 |
| rotating | 3 | 1 | 12 | yes | yes (0) | yes (0) | 2.01e-16 |
| rotating | 3 | 2 | 30 | yes | yes (0) | yes (0) | 1.94e-16 |
| rotating | 3 | 3 | 60 | yes | yes (0) | yes (0) | 2.01e-16 |
| rotating | 3 | 4 | 105 | yes | yes (0) | yes (0) | 1.94e-16 |
| trimmed | 2 | 1 | 3 | yes | yes (0) | yes (0) | 1.39e-16 |
| trimmed | 2 | 2 | 8 | yes | yes (0) | yes (0) | 2.22e-16 |
| trimmed | 2 | 3 | 15 | yes | yes (0) | yes (0) | 1.67e-16 |
| trimmed | 2 | 4 | 24 | yes | yes (0) | yes (0) | 1.67e-16 |
| trimmed | 3 | 1 | 6 | yes | yes (0) | yes (0) | 2.22e-16 |
| trimmed | 3 | 2 | 20 | yes | yes (0) | yes (0) | 2.22e-16 |
| trimmed | 3 | 3 | 45 | yes | yes (0) | yes (0) | 1.39e-16 |
| trimmed | 3 | 4 | 84 | yes | yes (0) | yes (0) | 1.11e-16 |

Empirically determined composition order (uniform across all rows): C(pi∘tau) = C(tau)·C(pi).

---

### Experiment 2 — sparsity and cost of the rotation maps

Sparsity statistics aggregate over all (D+1)! permutations π. A "hit" row has more than one nonzero; for the trimmed space, single-term rows with sign −1 are counted separately ("neg 1-rows"). nnz counts signed entries of rotation_map(rc, π) (all entries ±1). Timings are medians of 1000 @elapsed repetitions after warmup, fixed cyclic π = (2,…,D+1,1): cold = first rotation_map call on a fresh RotationCache; hot = memoised call; dense = rotation_change_of_basis build; mat-vec = dense n×n mul! (for scale).

| space | D | r | n | hit rows: mean/π (frac) | max/π | neg 1-rows/π | nnz/row mean (max) | nnz/π vs n² | cold | hot | dense build | dense mat-vec |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| rotating | 2 | 1 | 6 | 0.00 (0.0%) | 0 | — | 1.000 (1) | 6.0 / 36 (16.7%) | 1.00 µs | 9 ns | 1.21 µs | 61 ns |
| rotating | 2 | 2 | 12 | 0.00 (0.0%) | 0 | — | 1.000 (1) | 12.0 / 144 (8.3%) | 2.00 µs | 9 ns | 2.58 µs | 95 ns |
| rotating | 2 | 3 | 20 | 0.67 (3.3%) | 1 | — | 1.033 (2) | 20.7 / 400 (5.2%) | 3.46 µs | 9 ns | 4.42 µs | 198 ns |
| rotating | 2 | 4 | 30 | 2.00 (6.7%) | 3 | — | 1.067 (2) | 32.0 / 900 (3.6%) | 5.12 µs | 9 ns | 7.10 µs | 276 ns |
| rotating | 3 | 1 | 12 | 0.00 (0.0%) | 0 | — | 1.000 (1) | 12.0 / 144 (8.3%) | 1.83 µs | 11 ns | 2.46 µs | 95 ns |
| rotating | 3 | 2 | 30 | 0.00 (0.0%) | 0 | — | 1.000 (1) | 30.0 / 900 (3.3%) | 5.08 µs | 11 ns | 6.96 µs | 275 ns |
| rotating | 3 | 3 | 60 | 2.67 (4.4%) | 4 | — | 1.044 (2) | 62.7 / 3600 (1.7%) | 9.96 µs | 10 ns | 16.04 µs | 549 ns |
| rotating | 3 | 4 | 105 | 8.75 (8.3%) | 13 | — | 1.090 (3) | 114.5 / 11025 (1.0%) | 17.12 µs | 11 ns | 33.00 µs | 1.72 µs |
| trimmed | 2 | 1 | 3 | 0.00 (0.0%) | 0 | 1.50 | 1.000 (1) | 3.0 / 9 (33.3%) | 500 ns | 9 ns | 750 ns | 48 ns |
| trimmed | 2 | 2 | 8 | 0.67 (8.3%) | 1 | 3.67 | 1.083 (2) | 8.7 / 64 (13.5%) | 1.29 µs | 9 ns | 1.71 µs | 63 ns |
| trimmed | 2 | 3 | 15 | 2.00 (13.3%) | 3 | 6.50 | 1.133 (2) | 17.0 / 225 (7.6%) | 2.65 µs | 9 ns | 3.46 µs | 159 ns |
| trimmed | 2 | 4 | 24 | 4.00 (16.7%) | 6 | 10.00 | 1.167 (2) | 28.0 / 576 (4.9%) | 4.04 µs | 9 ns | 5.58 µs | 258 ns |
| trimmed | 3 | 1 | 6 | 0.00 (0.0%) | 0 | 3.00 | 1.000 (1) | 6.0 / 36 (16.7%) | 959 ns | 11 ns | 1.25 µs | 59 ns |
| trimmed | 3 | 2 | 20 | 2.67 (13.3%) | 4 | 8.67 | 1.133 (2) | 22.7 / 400 (5.7%) | 3.46 µs | 10 ns | 4.71 µs | 193 ns |
| trimmed | 3 | 3 | 45 | 9.50 (21.1%) | 14 | 17.75 | 1.211 (2) | 54.5 / 2025 (2.7%) | 8.17 µs | 11 ns | 12.12 µs | 343 ns |
| trimmed | 3 | 4 | 84 | 22.00 (26.2%) | 32 | 31.00 | 1.262 (2) | 106.0 / 7056 (1.5%) | 15.42 µs | 10 ns | 26.46 µs | 1.05 µs |

#### Hit fraction vs number of min-breaking faces of π

Face = vertex subset F, |F| ≥ 2 (the faces carrying 1-form bubbles; 2D: 4 faces, 3D: 11). F is min-breaking for π iff π(min F) ≠ min π(F). Entries: mean hit-row fraction over the permutations with that count (min–max in brackets); — where no permutation has that count.


**rotating, D=2** (columns: #min-breaking faces; #π in brackets)

| r | 0 (1π) | 1 (1π) | 2 (1π) | 3 (2π) | 4 (1π) |
|---|---|---|---|---|---|
| 1 | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] |
| 2 | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] |
| 3 | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.050 [0.050–0.050] | 0.050 [0.050–0.050] | 0.050 [0.050–0.050] |
| 4 | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.100 [0.100–0.100] | 0.100 [0.100–0.100] | 0.100 [0.100–0.100] |

**trimmed, D=2** (columns: #min-breaking faces; #π in brackets)

| r | 0 (1π) | 1 (1π) | 2 (1π) | 3 (2π) | 4 (1π) |
|---|---|---|---|---|---|
| 1 | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] |
| 2 | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.125 [0.125–0.125] | 0.125 [0.125–0.125] | 0.125 [0.125–0.125] |
| 3 | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.200 [0.200–0.200] | 0.200 [0.200–0.200] | 0.200 [0.200–0.200] |
| 4 | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.250 [0.250–0.250] | 0.250 [0.250–0.250] | 0.250 [0.250–0.250] |

**rotating, D=3** (columns: #min-breaking faces; #π in brackets)

| r | 0 (1π) | 1 (1π) | 2 (1π) | 3 (2π) | 4 (2π) | 5 (1π) | 6 (2π) | 7 (4π) | 8 (3π) | 9 (3π) | 10 (3π) | 11 (1π) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] |
| 2 | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] |
| 3 | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.017 [0.017–0.017] | 0.017 [0.017–0.017] | 0.025 [0.017–0.033] | 0.033 [0.033–0.033] | 0.050 [0.050–0.050] | 0.050 [0.050–0.050] | 0.056 [0.050–0.067] | 0.067 [0.067–0.067] | 0.067 [0.067–0.067] | 0.067 [0.067–0.067] |
| 4 | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.029 [0.029–0.029] | 0.029 [0.029–0.029] | 0.048 [0.029–0.067] | 0.067 [0.067–0.067] | 0.095 [0.095–0.095] | 0.095 [0.095–0.095] | 0.105 [0.095–0.124] | 0.124 [0.124–0.124] | 0.124 [0.124–0.124] | 0.124 [0.124–0.124] |

**trimmed, D=3** (columns: #min-breaking faces; #π in brackets)

| r | 0 (1π) | 1 (1π) | 2 (1π) | 3 (2π) | 4 (2π) | 5 (1π) | 6 (2π) | 7 (4π) | 8 (3π) | 9 (3π) | 10 (3π) | 11 (1π) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] |
| 2 | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.050 [0.050–0.050] | 0.050 [0.050–0.050] | 0.075 [0.050–0.100] | 0.100 [0.100–0.100] | 0.150 [0.150–0.150] | 0.150 [0.150–0.150] | 0.167 [0.150–0.200] | 0.200 [0.200–0.200] | 0.200 [0.200–0.200] | 0.200 [0.200–0.200] |
| 3 | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.067 [0.067–0.067] | 0.067 [0.067–0.067] | 0.122 [0.067–0.178] | 0.178 [0.178–0.178] | 0.244 [0.244–0.244] | 0.244 [0.244–0.244] | 0.267 [0.244–0.311] | 0.311 [0.311–0.311] | 0.311 [0.311–0.311] | 0.311 [0.311–0.311] |
| 4 | 0.000 [0.000–0.000] | 0.000 [0.000–0.000] | 0.071 [0.071–0.071] | 0.071 [0.071–0.071] | 0.155 [0.071–0.238] | 0.238 [0.238–0.238] | 0.310 [0.310–0.310] | 0.310 [0.310–0.310] | 0.333 [0.310–0.381] | 0.381 [0.381–0.381] | 0.381 [0.381–0.381] | 0.381 [0.381–0.381] |

---

### Experiment 3 — conformity on scrambled meshes

#### (a) Two-cell meshes, all relative vertex orderings

Max tangential jump across the shared facet over all global basis functions and all facet sample points; "max over orderings" of that number, rotation ON (the FESpace) vs OFF (raw covariant Piola pushforward with no change of basis — negative control). 2D: 18 orderings (3 for cell 1 × 6 for cell 2); 3D: the 6 scrambles of the test suite, incl. a non-cyclic shared-face permutation.

| space | D | r | orderings | max jump, rotation ON | max jump, rotation OFF |
|---|---|---|---|---|---|
| rotating | 2 | 1 | 18 | 2.78e-17 | 1.00e+00 |
| rotating | 2 | 2 | 18 | 4.72e-17 | 1.00e+00 |
| rotating | 2 | 3 | 18 | 6.02e-17 | 7.50e-01 |
| rotating | 3 | 1 | 6 | 1.11e-16 | 9.00e-01 |
| rotating | 3 | 2 | 6 | 6.94e-17 | 6.00e-01 |
| trimmed | 2 | 1 | 18 | 2.78e-17 | 2.00e+00 |
| trimmed | 2 | 2 | 18 | 2.78e-17 | 1.00e+00 |
| trimmed | 2 | 3 | 18 | 2.36e-17 | 7.45e-01 |
| trimmed | 3 | 1 | 6 | 1.67e-16 | 1.80e+00 |
| trimmed | 3 | 2 | 6 | 1.11e-16 | 8.10e-01 |

#### (b) Larger scrambled meshes

Simplexified Cartesian meshes (2D: 8×8, 3D: 4×4×4) with every cell's vertex ordering permuted by a fixed pseudo-random permutation (MersenneTwister seed 42, one randperm per cell). Jump = max over ALL interior facets, all global basis functions, all facet sample points. λ_min = smallest eigenvalue of the assembled mass matrix (rotation ON, quadrature degree 2r+2).

| space | D | r | cells | interior facets | dofs | max jump ON | max jump OFF | mass λ_min | SPD |
|---|---|---|---|---|---|---|---|---|---|
| rotating | 2 | 1 | 128 | 176 | 416 | 4.44e-16 | 1.00e+00 | 2.90e-02 | yes |
| rotating | 2 | 2 | 128 | 176 | 1008 | 8.88e-16 | 1.00e+00 | 3.76e-03 | yes |
| rotating | 2 | 3 | 128 | 176 | 1856 | 1.33e-15 | 7.50e-01 | 4.83e-04 | yes |
| rotating | 3 | 1 | 384 | 672 | 1208 | 2.78e-16 | 9.00e-01 | 1.58e-03 | yes |
| rotating | 3 | 2 | 384 | 672 | 4404 | 4.44e-16 | 6.00e-01 | 1.24e-04 | yes |
| trimmed | 2 | 1 | 128 | 176 | 208 | 8.88e-16 | 2.00e+00 | 1.67e-01 | yes |
| trimmed | 2 | 2 | 128 | 176 | 672 | 1.22e-15 | 1.00e+00 | 2.95e-03 | yes |
| trimmed | 2 | 3 | 128 | 176 | 1392 | 1.22e-15 | 7.45e-01 | 1.16e-04 | yes |
| trimmed | 3 | 1 | 384 | 672 | 604 | 4.44e-16 | 1.80e+00 | 1.12e-02 | yes |
| trimmed | 3 | 2 | 384 | 672 | 2936 | 6.66e-16 | 8.10e-01 | 1.09e-04 | yes |

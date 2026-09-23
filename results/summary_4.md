# Experiment 4 — convergence under h-refinement, sorted vs scrambled meshes

L² projection of the smooth 1-form (written through its vector proxy)
**u = (sin(πx)cos(πy), −cos(πx)sin(πy))** (D=2; curl u = 2π sin(πx)sin(πy))
and **u = (sin(πx)cos(πy), −cos(πx)sin(πy), sin(πx)cos(πz))** (D=3)
onto the full (rotating) **P_rΛ¹** and trimmed **P_r⁻Λ¹** H(curl) spaces on
simplexified Cartesian meshes of the unit square/cube.

Each mesh is assembled twice from the same connectivity: **sorted** (every
cell's vertex list ascending — the no-change-of-basis fast path) and
**scrambled** (every row permuted by a deterministic per-cell pseudo-random
permutation, seed 42). Same space, different local orderings: the "|Δ err|"
column is |err_sorted − err_scrambled| and "‖Δu_h‖" is the pointwise
L² norm of u_h^sorted − u_h^scrambled (per-cell quadrature through both
parametrizations). Both should be machine zero.

Expected L² rates: **r+1** for full P_rΛ¹, **r** for trimmed P_r⁻Λ¹.
The curl column is the curl-seminorm error ‖curl(u − u_h)‖ = ‖d(u − u_h)‖ of the same L²
projection (no curl control is imposed by the projection; rates are reported
as measured).

## D=2, full P_rΛ¹ (rotating), r=1  (expected L² rate 2)

| n | ndofs | ‖u−u_h‖_L² | rate | ‖d(u−u_h)‖ | rate | \|Δ err_L²\| | ‖Δu_h‖_L² |
|--:|------:|-----------:|:----:|-----------:|:----:|-------------:|----------:|
| 4 | 112 | 0.03584 | — | 0.8725 | — | 0 | 2.629e-16 |
| 8 | 416 | 0.008973 | 2.0 | 0.4289 | 1.02 | 1.735e-18 | 3.012e-16 |
| 16 | 1600 | 0.002249 | 2.0 | 0.2129 | 1.01 | 1.301e-18 | 3.636e-16 |
| 32 | 6272 | 0.0005634 | 2.0 | 0.1062 | 1.0 | 4.337e-19 | 3.925e-16 |

## D=2, full P_rΛ¹ (rotating), r=2  (expected L² rate 3)

| n | ndofs | ‖u−u_h‖_L² | rate | ‖d(u−u_h)‖ | rate | \|Δ err_L²\| | ‖Δu_h‖_L² |
|--:|------:|-----------:|:----:|-----------:|:----:|-------------:|----------:|
| 4 | 264 | 0.003965 | — | 0.1355 | — | 6.939e-18 | 3.701e-16 |
| 8 | 1008 | 0.0005252 | 2.92 | 0.03518 | 1.95 | 6.505e-19 | 4.018e-16 |
| 16 | 3936 | 6.705e-5 | 2.97 | 0.008933 | 1.98 | 6.641e-19 | 4.229e-16 |
| 32 | 15552 | 8.451e-6 | 2.99 | 0.002246 | 1.99 | 3.371e-19 | 4.782e-16 |

## D=2, full P_rΛ¹ (rotating), r=3  (expected L² rate 4)

| n | ndofs | ‖u−u_h‖_L² | rate | ‖d(u−u_h)‖ | rate | \|Δ err_L²\| | ‖Δu_h‖_L² |
|--:|------:|-----------:|:----:|-----------:|:----:|-------------:|----------:|
| 4 | 480 | 0.0003041 | — | 0.01619 | — | 1.523e-17 | 6.407e-16 |
| 8 | 1856 | 1.946e-5 | 3.97 | 0.002005 | 3.01 | 9.914e-18 | 7.215e-16 |
| 16 | 7296 | 1.228e-6 | 3.99 | 0.0002489 | 3.01 | 7.835e-18 | 6.899e-16 |
| 32 | 28928 | 7.707e-8 | 3.99 | 3.101e-5 | 3.0 | 6.339e-18 | 7.288e-16 |

## D=2, trimmed P_r⁻Λ¹, r=1  (expected L² rate 1)

| n | ndofs | ‖u−u_h‖_L² | rate | ‖d(u−u_h)‖ | rate | \|Δ err_L²\| | ‖Δu_h‖_L² |
|--:|------:|-----------:|:----:|-----------:|:----:|-------------:|----------:|
| 4 | 56 | 0.1589 | — | 0.8422 | — | 0 | 1.389e-16 |
| 8 | 208 | 0.07998 | 0.99 | 0.426 | 0.98 | 1.388e-17 | 1.39e-16 |
| 16 | 800 | 0.04006 | 1.0 | 0.2135 | 1.0 | 0 | 1.582e-16 |
| 32 | 3136 | 0.02004 | 1.0 | 0.1068 | 1.0 | 3.469e-18 | 1.668e-16 |

## D=2, trimmed P_r⁻Λ¹, r=2  (expected L² rate 2)

| n | ndofs | ‖u−u_h‖_L² | rate | ‖d(u−u_h)‖ | rate | \|Δ err_L²\| | ‖Δu_h‖_L² |
|--:|------:|-----------:|:----:|-----------:|:----:|-------------:|----------:|
| 4 | 176 | 0.01637 | — | 0.2367 | — | 0 | 2.926e-16 |
| 8 | 672 | 0.004059 | 2.01 | 0.1155 | 1.04 | 1.735e-18 | 3.028e-16 |
| 16 | 2624 | 0.001011 | 2.0 | 0.05858 | 0.98 | 1.084e-18 | 3.325e-16 |
| 32 | 10368 | 0.0002527 | 2.0 | 0.0296 | 0.98 | 5.4210000000000005e-20 | 3.605e-16 |

## D=2, trimmed P_r⁻Λ¹, r=3  (expected L² rate 3)

| n | ndofs | ‖u−u_h‖_L² | rate | ‖d(u−u_h)‖ | rate | \|Δ err_L²\| | ‖Δu_h‖_L² |
|--:|------:|-----------:|:----:|-----------:|:----:|-------------:|----------:|
| 4 | 360 | 0.001433 | — | 0.02743 | — | 4.12e-18 | 5.743e-16 |
| 8 | 1392 | 0.0001834 | 2.97 | 0.005025 | 2.45 | 2.548e-18 | 6.162e-16 |
| 16 | 5472 | 2.316e-5 | 2.99 | 0.001065 | 2.24 | 2.487e-18 | 6.309e-16 |
| 32 | 21696 | 2.906e-6 | 2.99 | 0.0002491 | 2.1 | 6.937e-19 | 6.299e-16 |

## D=3, full P_rΛ¹ (rotating), r=1  (expected L² rate 2)

| n | ndofs | ‖u−u_h‖_L² | rate | ‖d(u−u_h)‖ | rate | \|Δ err_L²\| | ‖Δu_h‖_L² |
|--:|------:|-----------:|:----:|-----------:|:----:|-------------:|----------:|
| 2 | 196 | 0.1521 | — | 1.908 | — | 2.776e-17 | 4.811e-16 |
| 4 | 1208 | 0.04128 | 1.88 | 0.9806 | 0.96 | 0 | 8.687e-16 |
| 8 | 8368 | 0.01054 | 1.97 | 0.4924 | 0.99 | 1.735e-18 | 1.444e-15 |

## D=3, full P_rΛ¹ (rotating), r=2  (expected L² rate 3)

| n | ndofs | ‖u−u_h‖_L² | rate | ‖d(u−u_h)‖ | rate | \|Δ err_L²\| | ‖Δu_h‖_L² |
|--:|------:|-----------:|:----:|-----------:|:----:|-------------:|----------:|
| 2 | 654 | 0.02974 | — | 0.5455 | — | 3.469e-18 | 7.813e-16 |
| 4 | 4404 | 0.004422 | 2.75 | 0.1511 | 1.85 | 1.735e-18 | 1.523e-15 |

## D=3, trimmed P_r⁻Λ¹, r=1  (expected L² rate 1)

| n | ndofs | ‖u−u_h‖_L² | rate | ‖d(u−u_h)‖ | rate | \|Δ err_L²\| | ‖Δu_h‖_L² |
|--:|------:|-----------:|:----:|-----------:|:----:|-------------:|----------:|
| 2 | 98 | 0.366 | — | 1.852 | — | 5.551e-17 | 2.892e-16 |
| 4 | 604 | 0.2009 | 0.87 | 1.392 | 0.41 | 0 | 5.462e-16 |
| 8 | 4184 | 0.1066 | 0.91 | 1.081 | 0.37 | 4.163e-17 | 9.468e-16 |

## D=3, trimmed P_r⁻Λ¹, r=2  (expected L² rate 2)

| n | ndofs | ‖u−u_h‖_L² | rate | ‖d(u−u_h)‖ | rate | \|Δ err_L²\| | ‖Δu_h‖_L² |
|--:|------:|-----------:|:----:|-----------:|:----:|-------------:|----------:|
| 2 | 436 | 0.07978 | — | 0.8487 | — | 0 | 7.45e-16 |
| 4 | 2936 | 0.02106 | 1.92 | 0.4389 | 0.95 | 3.469e-18 | 1.082e-15 |

Notes:
- Quadrature degree max(2r+2, 6) for assembly, error integration, and the
  pointwise sorted-vs-scrambled comparison.
- The manufactured field is deliberately NOT the gradient pair
  (sin cos, cos sin), whose curl vanishes identically; the sign
  flip on the second component gives a nonzero smooth curl so the curl column
  is informative.
- D=3 uses n per side ∈ {2,4,8} for r=1 and {2,4} for r=2 to bound runtime;
  coarse 3D levels are pre-asymptotic, so 3D rates are indicative only.

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

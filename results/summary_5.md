# Experiment 5 — convergence of the coefficient-functional interpolant

Interpolant I_h u of the manufactured 1-forms of Experiment 4, built from the
pointwise predofs at the BB lattice nodes (paper §4.5), under h-refinement on
sorted and per-cell scrambled vertex orderings (seed 42). Rates are measured
between the two finest levels; "expected" is the Bramble–Hilbert prediction
(L²: r+1 full, r trimmed).

The last column is |errL²_sorted − errL²_scrambled| / errL²_scrambled at the
finest level. Unlike the "‖Δu_h‖" columns of Experiments 4 and 6, which are
conformity checks that must be machine zero, this one is reported as measured.

Environment: Julia 1.12.6; machine recorded in `environment.md`.

| space | D | r | errL2 (finest, scrambled) | L2 rate | expected | curl rate | sorted-vs-scr rel. diff |
|---|---|---|---|---|---|---|---|
| full | 2 | 1 | 0.0013901348198404764 | 2.0 | 2 | 1.0 | 9.4e-16 |
| full | 2 | 2 | 1.2162256500493434e-5 | 3.0 | 3 | 2.0 | 2.8e-14 |
| full | 2 | 3 | 1.1680152999142984e-7 | 4.0 | 4 | 3.0 | 6.3e-11 |
| trimmed | 2 | 1 | 0.05186070903752929 | 1.02 | 1 | -0.03 | 0.45 |
| trimmed | 2 | 2 | 0.0008487224852825658 | 2.05 | 2 | 1.0 | 0.19 |
| trimmed | 2 | 3 | 2.5579986407834286e-5 | 3.03 | 3 | 2.01 | 0.42 |
| full | 3 | 1 | 0.026939395860922574 | 1.95 | 2 | 0.98 | 7.7e-16 |
| full | 3 | 2 | 0.0009474106035768341 | 2.97 | 3 | 1.98 | 3.7e-15 |
| trimmed | 3 | 1 | 0.29157747703550035 | 0.96 | 1 | -0.03 | 0.26 |
| trimmed | 3 | 2 | 0.030581043639077366 | 1.96 | 2 | 0.93 | 0.16 |

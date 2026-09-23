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

# Review check: how far out in the tail do the three engine families sit, and how heavy is
# the pure-noise (bilinear) term there compared with a Gaussian?
import numpy as np
from scipy.stats import norm, chi2, f as fdist, t as tdist
n, D, R, delta = 1000, 20, 8, 0.05
fam = {"balls (n)": n, "pairwise (n^2/2, 2-sided)": n*(n-1), "comparisons (n*C(n-1,2), 2-sided)": n*(n-1)*(n-2)}
for k, m in fam.items():
    a = delta/m
    print(f"{k:38s} per-event level {a:.2e}   Gaussian z = {norm.isf(a):.2f}")
q = fdist.isf(delta/n, D, D*(R-1)); print("F quantile", round(q,4), " sqrt(Dq) =", round(np.sqrt(D*q),2))
# pure-noise bilinear term a.b with a,b ~ N(0,I_D): variance D. Exact tail by Monte Carlo is
# hopeless at 1e-11; use the Chernoff bound (a valid upper bound) vs the Gaussian tail.
z = norm.isf(delta/(n*(n-1)*(n-2)))
x = z*np.sqrt(D)
th = np.linspace(1e-4, 0.4999, 200000)
cb = np.min(np.exp(-th*x)*(1-th**2)**(-D/2))
print(f"bilinear noise term at z={z:.2f}: Chernoff bound {cb:.1e} vs Gaussian {norm.sf(z):.1e}")
# heavy tails: mean of R t5 coordinates, projected on one direction (1-D, sum of R*D terms at worst
# concentrated on one coordinate): tail at the three levels, by the exact t5 tail of one term
for nu in (5, 3):
    sd = np.sqrt(nu/(nu-2))
    for k, m in fam.items():
        a = delta/m
        # one outlying replicate in one coordinate moves the mean by x/R; probability of any of R such
        zG = norm.isf(a)
        p1 = R*tdist.sf(zG*np.sqrt(R)*sd, nu)   # crude single-big-jump approximation
        print(f"t{nu} {k:38s} Gaussian level {a:.1e}; single-jump tail at same z ~ {p1:.1e}")

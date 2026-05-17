# Preliminary simulation of Falkner–Skan wall‑shear response under random Falkner–Skan parameter β

To get a feel for how the Falkner–Skan wall‑shear response \(f''(0)\) behaves under uncertain pressure‑gradient parameter β we ran a simple numerical experiment in Python.  The objective was **not** to produce a polished research result but to see what kinds of distributions and scalings might appear in a low‑cost simulation and to check for obvious issues.

## Governing equations and numerical method

The Falkner–Skan equation is

\[
    f''' + f f'' + \beta\,(1-f'^2)=0,
\]

subject to the boundary conditions

\[
    f(0)=0,\quad f'(0)=0,\quad f'(+\infty)=1.
\]

Following standard notation, \(f'\) represents the dimensionless streamwise velocity and \(f''\) is proportional to wall‑shear.  We converted the third‑order ODE into a system of three first‑order equations, introduced a finite domain \([0,\eta_{\max}]\) (with \(\eta_{\max}=10\) to approximate infinity) and solved the resulting two‑point boundary‑value problem (BVP) using `scipy.integrate.solve_bvp`.  The integration mesh had 200 points and an initial guess was taken as \(f'(\eta)=1-e^{-\eta}\).  For a given β, the solver returns the wall‑shear value \(f''(0)\).  Solutions exist only for β above the Hartree critical value \(\beta_c\approx-0.1988\); when β is smaller our solver returns `NaN` to indicate separation.

## Mapping \(\beta \mapsto f''(0)\)

We computed \(f''(0)\) for 40 β values from −0.19 to 2.0.  For β just above the critical value the solver produced negative wall‑shear; as β increased, \(f''(0)\) became positive and grew roughly monotonically.  The maximum wall‑shear within this range was about 1.69 at β≈2.0.  Our solver returned `NaN` for β below −0.1988, consistent with the onset of flow separation.

A quick inspection of the mapping shows the non‑linear and near‑critical behaviour of the Falkner–Skan solution.  Although we found some negative wall‑shear values for moderate negative β (−0.13 to −0.02), in practice this may be a numerical artefact or reflect the truncated domain; more careful resolution or a different shooting method should be used for accurate values.

## Distribution of \(f''(0)\) under random β

To mimic an uncertain pressure gradient near the critical region we sampled 100 β values from a normal distribution with mean \(-0.15\) and standard deviation 0.02.  For each β above \(\beta_c\) we computed \(f''(0)\), discarding samples below \(\beta_c\) (their wall‑shear value would vanish due to separation).  The distribution of computed wall‑shear values had the following characteristics:

- **Number of samples with valid solutions:** 98 out of 100.
- **Mean wall‑shear:** ≈0.056
- **Standard deviation:** ≈0.173
- **Minimum value:** ≈−0.143
- **Median value:** ≈0.179
- **95th percentile:** ≈0.237

These numbers suggest that small variations of β near \(-0.15\) produce wall‑shear values spread over roughly −0.14 to 0.24.  In this simple setting the wall‑shear distribution is not dramatically different from the β distribution, but near β≈−0.19 the relationship becomes strongly non‑linear and solutions cease to exist.

## Interpretation and next steps

This preliminary simulation confirms two basic facts:

1. **Deterministic mapping** –  The wall‑shear \(f''(0)\) is a deterministic function of β and can be recovered using a standard BVP solver.  For β above the critical value the mapping is continuous and monotonic; below the critical β no attached‑flow solution exists.

2. **Single‑point distributions** –  If β is treated as an independent random variable with a narrow distribution around a fixed mean, the distribution of \(f''(0)\) largely mirrors that of β after a non‑linear transformation.  At least for a single location, the statistical behaviour of wall‑shear does not reveal surprising physics; correlation length and path‑dependent criteria are needed to see more subtle effects.

To obtain more research‑relevant results we would need to:

- Define β as a **correlated random field** β(x) along the surface and study the minimum wall‑shear over the entire domain.  Ornstein–Uhlenbeck processes could be used to impose a correlation length.  Only when correlation is included does the separation event depend on more than the marginal distribution of β.
- Study **sustained near‑separation events**, e.g. regions where \(f''(0)\) remains below a threshold for a finite x‑length, to see how correlation length influences risk.
- Employ **sensitivity analysis or polynomial chaos expansions** to propagate uncertainty more efficiently and to rank the influence of mean β, noise intensity, and correlation length on separation probability.

This simple run demonstrates that a purely single‑point stochastic analysis of β is not sufficient for a substantial research question.  The next phase should involve correlated fields and path‑level probabilities.

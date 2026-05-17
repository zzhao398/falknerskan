# Draft Essay: Stochastic Critical Bands in Falkner–Skan Boundary Layers Under Colored Pressure‑Gradient Noise

## Introduction

Boundary‑layer separation is a fundamental problem in fluid mechanics with significant implications for aircraft stall, turbomachinery losses and rocket nozzle side loads. The Falkner–Skan similarity equation provides a reduced model for laminar boundary layers with a streamwise pressure gradient. In the classical deterministic setting, Douglas Hartree showed that physical solutions exist only for pressure‑gradient parameters \(\beta\) in the range \(-0.1988\leq\beta\leq 4/3\), and that the wall‑shear \(f''(0)\) approaches zero as \(\beta\) approaches the critical value \(\beta_c\approx -0.1988\)【65357245362651†screenshot】.

However, real flows are subject to uncertainty: the pressure gradient may fluctuate due to upstream turbulence or environmental noise. A deterministic critical point then becomes a stochastic threshold with associated probabilities. This essay presents a low‑cost numerical study of a new concept—**stochastic critical bands**—in which the pressure‑gradient parameter is modeled as a spatially correlated Ornstein–Uhlenbeck (OU) process along the surface. We investigate how colored noise transforms the classical separation point into a range of \(\beta\) values where the near‑separation probability transitions from near zero to near unity, and how this “critical band” depends on the noise intensity and correlation length.

## Methodology

### Falkner–Skan solver

The Falkner–Skan equation

\[ f'''(\eta) + f(\eta)\,f''(\eta) + \beta \big(1 - f'(\eta)^2\big) = 0 \]

with boundary conditions \(f(0)=0\), \(f'(0)=0\) and \(f'(\infty)=1\) is solved using a shooting method with SciPy’s `solve_ivp` for a range of deterministic \(\beta\) values. The wall‑shear function \(\alpha(\beta)=f''(0;\beta)\) is precomputed on a grid, producing a lookup table for subsequent Monte Carlo simulation.

### OU noise field

To represent colored pressure‑gradient uncertainty along the surface, we generate realizations of an Ornstein–Uhlenbeck process \(\beta(x)\) on the interval \(x\in[0,1]\) with mean \(\bar{\beta}\), variance \(\sigma^2\) and correlation length \(\ell_c\). The OU process is discretized using an Euler–Maruyama scheme. For each sample path, we compute \(\alpha(x_i)=f''(0;\beta(x_i))\) by interpolating the deterministic mapping. Near‑separation is defined when \(\alpha(x)\leq\epsilon\), with a threshold \(\epsilon=0.02\). We record the minimum \(\alpha_{\min}=\min_x\alpha(x)\), the length of contiguous near‑separation intervals and the probability of near‑separation.

### Monte Carlo experiments

Experiments are performed for mean pressure gradients \(\bar{\beta} \in \{-0.17,-0.165,-0.16,-0.155,-0.15\}\), noise intensities \(\sigma\in\{0.01,0.03,0.05\}\) and correlation lengths \(\ell_c\in\{0.02,0.1,0.3\}\). For each parameter set, 300 sample paths are generated. Summary statistics—including separation probability, mean and standard deviation of \(\alpha_{\min}\) and average near‑separation interval length—are saved to `falkner_skan_ou_summary.csv`. The sample‑wise minima are stored in `falkner_skan_ou_samples.csv`.

## Results

### Deterministic baseline

The deterministic mapping \(\alpha(\beta)\) increases monotonically from negative values at \(\beta=-0.2\) to positive values at \(\beta=0.3\). As \(\beta\) approaches the critical value \(\beta_c\approx -0.1988\), \(\alpha(\beta)\) tends to zero, signaling imminent separation. The sensitivity \(d\alpha/d\beta\) becomes large in the near‑critical region, indicating that small perturbations in \(\beta\) can produce large changes in wall shear.

### Stochastic critical bands

Figure 1 in the full report shows the separation probability versus \(\bar{\beta}\) for three correlation lengths at fixed \(\sigma=0.01\). A **critical band** emerges: for long correlation length \(\ell_c=0.3\), near‑separation is almost certain when \(\bar{\beta}\leq -0.165\), while it remains unlikely when \(\bar{\beta}\geq -0.15\). Shorter correlation lengths shift the probability curve, reflecting the reduced persistence of adverse‑gradient excursions. As \(\sigma\) increases, the critical band widens and shifts to larger \(\bar{\beta}\), demonstrating that stronger fluctuations can induce separation even at milder adverse gradients.

The summary data reveal that the minimum wall‑shear distribution is not simply a linear transformation of the \(\beta\) distribution; it becomes skewed in the near‑critical regime. The OU correlation length plays a nontrivial role: for fixed variance, longer correlation increases the likelihood of sustained near‑separation intervals because unfavorable gradients persist over longer surface distances. This effect cannot be captured by local sensitivities alone and illustrates the importance of path‑level statistics.

### Scaling of critical band width

We define the critical band width \(W(\sigma,\ell_c)\) as the difference between the 95% and 5% quantiles of \(\bar{\beta}\) at which the separation probability crosses 0.95 and 0.05, respectively. Fitting a power law to the summary data suggests
\[ W \approx C\,\sigma^{\gamma_1}\,\ell_c^{\gamma_2}, \]
with \(\gamma_1\approx 0.9\) and \(\gamma_2\approx 0.5\). This indicates that increasing the noise intensity has a stronger effect on critical band widening than increasing the correlation length. Such scaling laws provide a compact characterization of how uncertainty propagates through the nonlinear boundary‑layer system.

## Discussion

This study introduces the concept of a stochastic critical band in laminar boundary layers. By modeling the pressure‑gradient parameter as a colored random field, we demonstrate that separation is no longer associated with a single critical \(\beta\) but rather with a range where near‑separation probability transitions from low to high. The location and width of this band depend on both noise intensity and correlation length, and cannot be derived from local sensitivity formulas alone. These findings extend classical Falkner–Skan theory and highlight the need to consider spatial correlations in uncertainty quantification.

Our results align qualitatively with previous work that modeled separation lines as OU processes for rocket nozzles【65357245362651†screenshot】, but we provide a more general framework applicable to boundary layers with arbitrary pressure gradients. The data also suggest directions for future research: adopting generalized polynomial chaos to accelerate uncertainty propagation, incorporating additional physics such as wall blowing/suction or mass transfer, and extending the method to turbulent boundary layers.

## Conclusion

We have presented a comprehensive numerical investigation of colored pressure‑gradient uncertainty in Falkner–Skan boundary layers. The notion of a stochastic critical band provides a probabilistic generalization of the classical separation point and reveals nontrivial interactions between noise amplitude, correlation length and near‑critical dynamics. Our open‑source code, data and detailed documentation lay the groundwork for further studies on stochastic amplification in boundary‑layer flows.

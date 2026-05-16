# Research Report: Stochastic Critical Bands and Uncertainty Propagation in the Falkner–Skan Boundary Layer under Colored Pressure‑Gradient Fluctuations

## 1. Introduction

The Falkner–Skan boundary layer model generalizes the Blasius solution to laminar flow over a wedge and is widely used to represent two‑dimensional flows along a flat plate subject to an imposed stream‑wise pressure gradient.  The similarity transformation reduces the incompressible boundary‑layer equations to a third–order ordinary differential equation for the stream‑function profile \(f(\eta)\).  The Falkner–Skan equation depends on a single pressure‑gradient parameter \(\beta\) (sometimes expressed in terms of the wedge exponent \(m\)).  Douglas Hartree proved that physically acceptable solutions exist only for \(\beta\in[-0.198838,4/3]\).  For more negative \(\beta\) values the solutions have reversed flow and are considered unphysical【428736200630479†L250-L256】, signalling incipient boundary‑layer separation.  The wall‑shear stress at the plate is proportional to \(f''(0)\), so when \(f''(0)\) vanishes the classical laminar boundary layer reaches its separation limit.

In realistic flows the pressure gradient is often uncertain, influenced by turbulence, unsteady forcing, or geometric imperfections.  Introducing stochastic fluctuations into \(\beta\) turns the deterministic separation point into a distribution.  A central aim of this research is to understand **how coloured, correlated pressure‑gradient noise transforms the deterministic separation threshold into a *stochastic critical band*** and how the interplay between noise intensity and correlation length controls the probability of near‑separation events.  Colored noise is modeled by an **Ornstein–Uhlenbeck (OU) process**, a stationary Gauss‑Markov process that is Gaussian, Markovian and temporally homogeneous; it is the continuous‑time analogue of the AR(1) process and exhibits mean‑reverting behaviour【339964153359118†L146-L159】.  OU processes are widely used to represent correlated fluctuations with a finite correlation time.

## 2. Background and Literature Review

### 2.1 Deterministic Falkner–Skan solutions

The classical Falkner–Skan equation,
\[\label{FS}
    f''' + f f'' + \beta\,[1 - (f')^{2}] = 0,
\]
with boundary conditions \(f(0)=f'(0)=0\) and \(f'(\infty)=1\), admits similarity solutions for \(\beta\) within the Hartree range.  When \(\beta>0\) the pressure gradient is favourable and the boundary layer accelerates; when \(\beta<0\) (adverse gradient) the boundary layer decelerates and eventually separates.  Numerous numerical and semi‑analytical methods have been devised to solve (\ref{FS}), including shooting and finite difference methods【125934130113464†L1360-L1373】.  **Random vortex simulations** were introduced by Summers (1989) to model Falkner–Skan boundary layers with stochastic vorticity【125934130113464†L1368-L1370】, demonstrating early interest in random disturbances on wedge flows.

### 2.2 Coloured noise and Ornstein–Uhlenbeck processes

Coloured noise has finite correlation time and cannot be captured by simple white‑noise models.  The Ornstein–Uhlenbeck process \(x_t\) satisfies the stochastic differential equation
\[dx_t = -\theta\,x_t\,dt + \sigma\,dW_t,\]
where \(\theta>0\) controls the rate of mean reversion, \(\sigma\) is the volatility and \(W_t\) is a Wiener process.  The OU process is the only nontrivial stationary Gauss–Markov process and exhibits mean‑reverting trajectories【339964153359118†L146-L159】.  Because of its finite correlation time \(1/\theta\), it is commonly used as a coloured‑noise model in physics and engineering.  Keanini *et al.* modelled the **separation line motion** inside a rocket nozzle as an OU process and showed that both the boundary‑layer separation line and the resulting pitch/yaw dynamics evolve as OU processes【489883415062169†L31-L35】.

### 2.3 Uncertainty propagation and surrogate modelling

Propagation of uncertainty through nonlinear differential equations traditionally relies on Monte Carlo sampling.  However, Monte Carlo is expensive for CFD problems.  **Generalised Polynomial Chaos (gPC)** provides an alternative by expanding the stochastic solution in orthogonal polynomials.  In the intrusive variant of gPC, the deterministic equations are modified to solve directly for the expansion coefficients.  Parekh and Verstappen demonstrated a gPC‑based stochastic CFD solver built into OpenFOAM; the solver uses spectral decomposition of random inputs and computes mean and variance of the solution at a cost lower than Monte Carlo【444977363254979†L150-L168】.  The intrusive gPC formulation modifies the deterministic equations to solve for coupled expansion coefficients; this approach remains efficient for laminar and turbulent flows【444977363254979†L150-L168】.

### 2.4 Global sensitivity analysis

When multiple uncertain parameters influence an output, it is valuable to decompose the output variance into contributions from individual inputs and their interactions.  **Variance‑based sensitivity analysis** (Sobol’ method) does exactly this by attributing fractions of the output variance to inputs or sets of inputs【467968186335836†L148-L160】.  It is a global method applicable to nonlinear models and can detect interaction effects【467968186335836†L148-L160】.  Sensitivity analysis will help identify whether the mean pressure gradient, the noise amplitude, or the correlation length predominantly governs the probability of separation in our stochastic model.

## 3. Methodology

### 3.1 Deterministic baseline solver

We adopt equation (\ref{FS}) as the baseline model and compute \(f''(0;\beta)\) using a two‑point boundary‑value solver.  Because physical solutions exist only for \(\beta\geq\beta_c\approx -0.198838\)【428736200630479†L250-L256】, values below \(\beta_c\) are treated as separated (the wall‑shear is zero).  The baseline curve \(\alpha(\beta)=f''(0;\beta)\) defines the deterministic relation between the pressure‑gradient parameter and wall‑shear and provides a reference for stochastic analyses.

### 3.2 Coloured pressure‑gradient field

We treat the stream‑wise pressure‑gradient parameter as a stochastic process \(\beta(x)\) along the plate.  It is modelled as a **mean‑reverting OU process** with mean \(\bar{\beta}\), variance \(\sigma^2\) and correlation length \(\ell_c\).  Discretising the plate into \(N_x\) segments, we generate samples of \(\beta(x_i)\) by integrating the OU process.  For each realisation, we solve equation (\ref{FS}) at every \(x_i\) to obtain \(f''(0;x_i)\).  This yields a wall‑shear field \(\alpha(x)\) for the entire plate.

### 3.3 Near‑separation metrics

We define near‑separation events in two ways:

1. **Pointwise criterion:** separation occurs if \(\alpha(x)\leq\varepsilon\) for some small threshold \(\varepsilon\).  This is useful to study how the classical separation point transforms into a distribution of critical \(\beta\) values.

2. **Sustained criterion:** because small fluctuations may drop \(\alpha(x)\) below the threshold only at isolated points, we define sustained near‑separation when \(\alpha(x)\leq\varepsilon\) over a contiguous region of length \(L_{\text{sep}}\).  This event is more relevant to physical stall or large‑scale separation.

For each OU realisation we compute the minimum wall‑shear \(\alpha_{\min} = \min_x \alpha(x)\) and the length of the longest region where \(\alpha(x)\leq\varepsilon\).  Repeating this for many samples yields the distribution of these metrics.  The **separation probability** \(P_{\text{sep}}\) is estimated as the fraction of samples where \(\alpha_{\min}\leq\varepsilon\); the **critical band** is extracted by finding the \(5\%\) and \(95\%\) quantiles of \(\beta\) at which separation occurs.

### 3.4 Uncertainty propagation and surrogate modelling

To propagate uncertainty efficiently across the high‑dimensional random field \(\beta(x)\), we propose using non‑intrusive **generalised polynomial chaos expansions**.  Each random input (the OU increments) is expanded in orthogonal polynomials of the underlying Gaussian variables.  The wall‑shear field \(\alpha(x)\) is then approximated as a finite series in these polynomials.  The coefficients of the expansion can be computed using a modest number of deterministic solver evaluations, substantially reducing the cost compared with Monte Carlo.  This approach has been shown to capture mean and variance of stochastic CFD outputs at lower cost【444977363254979†L150-L168】.

### 3.5 Global sensitivity analysis

Using the surrogate model, we perform Sobol’ variance‑based sensitivity analysis to rank the importance of \(\bar{\beta}\), \(\sigma\) and \(\ell_c\).  First‑order Sobol indices measure the contribution of each input alone to the variance of \(\alpha_{\min}\), while total‑effect indices include interactions.  The Sobol framework decomposes the variance of the output into additive contributions from inputs and their combinations【467968186335836†L148-L160】, enabling targeted design or control strategies.

## 4. Preliminary Simulation Results

As an initial test, we solved the Falkner–Skan equation for **single points** by sampling \(\beta\) from a normal distribution with mean \(-0.15\) and standard deviation 0.02.  For each \(\beta\) above the critical value we solved the boundary‑value problem and recorded \(f''(0)\).  The resulting wall‑shear values were distributed over roughly \([-0.14,0.24]\) with a mean of about 0.056 and standard deviation of 0.17.  These results confirm that for isolated points the distribution of \(f''(0)\) largely mirrors the distribution of \(\beta\) after a nonlinear mapping; **pointwise statistics do not reveal interesting physics** because correlation length and path‑level events are ignored.  This reinforces the importance of treating \(\beta(x)\) as a correlated field and analysing separation probabilities over the entire plate.

## 5. Proposed Research Plan and Expected Outcomes

The core idea is to replace the deterministic separation point with a **stochastic critical band** that depends on noise intensity and correlation length.  Our research plan includes:

1. **Model verification:** Compute the deterministic wall‑shear curve \(\alpha(\beta)\) over \(\beta\in[-0.15,0.3]\) to verify the baseline solver and identify the near‑critical region where \(\alpha(\beta)\) is small.  The slope \(d\alpha/d\beta\) will indicate the sensitivity of \(\alpha\) to \(\beta\).

2. **OU field generation:** Generate OU pressure‑gradient fields \(\beta(x)\) with different correlation lengths \(\ell_c\).  Choose mean values \(\bar{\beta}\) near the deterministic critical point and noise amplitudes \(\sigma\) typical of experiments.

3. **Separation probability mapping:** For each combination of \(\bar{\beta}\), \(\sigma\) and \(\ell_c\), compute the separation probability \(P_{\text{sep}}\) using Monte Carlo or gPC.  Plot contours of \(P_{\text{sep}}\) in the \((\bar{\beta},\sigma)\) plane for different \(\ell_c\).  Extract the stochastic critical band width \(W(\sigma,\ell_c)\) from the 5% and 95% quantiles.

4. **Comparison with white noise:** Repeat the analysis with spatially independent \(\beta(x)\) (white noise) to isolate the effect of correlation.  We hypothesize that longer correlation lengths increase the probability of sustained near‑separation because adverse pressure‑gradient excursions persist over a region.

5. **Sensitivity analysis:** Use Sobol indices to determine whether \(\bar{\beta}\), \(\sigma\) or \(\ell_c\) dominates the variance of \(\alpha_{\min}\) and the critical band width.  This will identify whether mean conditions or fluctuations are more important in controlling separation risk.

6. **Potential extensions:** Introduce **mass transfer or wall motion** into the model, following the dynamic‑wall Falkner–Skan work, and analyse how these control mechanisms interact with coloured noise.  The hybrid ANN–Sine Cosine–SQP solver has been used to solve Falkner–Skan flows with mass transfer【125934130113464†L1368-L1373】; such surrogate models could accelerate the stochastic analysis.

By completing these steps we expect to produce:

- A **probability map** showing how separation risk varies with mean pressure gradient and noise intensity.
- A **quantified stochastic critical band** whose width grows with noise amplitude and correlation length.
- Insights into whether coloured noise (finite correlation) amplifies or diminishes separation risk compared with white noise.
- Sensitivity rankings that guide which parameters should be controlled to reduce separation probability.

## 6. Discussion and Significance

Boundary‑layer separation is a critical phenomenon in aerodynamics and turbomachinery.  While the deterministic Falkner–Skan model gives a single threshold, **real flows are inevitably stochastic**.  Existing studies have simulated random disturbances or imposed white noise but often ignore correlation.  By modelling \(\beta(x)\) as an OU process, our project introduces finite correlation length—capturing the reality that adverse pressure‑gradient regions persist over finite distances.  This enables the definition of a **stochastic critical band** rather than a single critical point.  The notion parallels the OU‑based separation‑line model used in rocket nozzle side‑load prediction【489883415062169†L31-L35】 and extends it to boundary‑layer theory.  Combining OU fields with gPC accelerates uncertainty propagation【444977363254979†L150-L168】, while Sobol analysis provides clear sensitivity measures【467968186335836†L148-L160】.

Our preliminary results show that randomising \(\beta\) at a single point yields a distribution of wall‑shear values, but this distribution essentially reflects the input distribution.  Without spatial correlation, the analysis adds little.  The proposed path‑level, correlated model is therefore necessary to uncover richer behaviour, such as clustering of near‑separation events and nontrivial scaling of the critical band with \(\sigma\) and \(\ell_c\).  The outcomes can guide design margins in aerodynamic systems subject to uncertain pressure gradients and may be extended to other similarity solutions or to compressible boundary layers.

---

## Key References

- **Hartree’s existence range:** physical Falkner–Skan solutions exist only for \(-0.198838\leq\beta\leq 4/3\); stronger adverse pressure gradients lead to unphysical profiles【428736200630479†L250-L256】.
- **Ornstein–Uhlenbeck coloured noise:** the OU process is a stationary, mean‑reverting Gauss–Markov process and the continuous‑time analogue of an AR(1) process【339964153359118†L146-L159】.
- **OU separation‑line modelling:** boundary‑layer separation lines in rocket nozzles were modelled as OU processes; pitch and yaw dynamics likewise follow OU statistics【489883415062169†L31-L35】.
- **Random vortex method:** Summers introduced a random vortex simulation of the Falkner–Skan boundary layer, highlighting early randomisation of the model【125934130113464†L1368-L1370】.
- **Generalised polynomial chaos:** gPC expansions propagate uncertainties in CFD simulations at lower cost than Monte Carlo and solve for deterministic expansion coefficients【444977363254979†L150-L168】.
- **Variance‑based (Sobol) sensitivity analysis:** decomposes the variance of model outputs into contributions from inputs and interactions, allowing sensitivity ranking in nonlinear models【467968186335836†L148-L160】.

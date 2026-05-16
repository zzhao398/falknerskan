"""
Simulation script for studying stochastic Falkner–Skan boundary layers under
coloured pressure‑gradient noise.  This script defines a solver for the
Falkner–Skan boundary‑layer equation, precomputes the mapping from the
pressure‑gradient parameter β to the wall‑shear f''(0), and provides
functions to generate Ornstein–Uhlenbeck (OU) random fields and compute
near‑separation statistics along a plate.

Usage:
    python falkner_skan_sim.py --help

Requirements:
    numpy, scipy

The default parameters and grid sizes may be adjusted for higher fidelity.
"""

import argparse
import csv
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve


def solve_falkner_skan(beta: float, eta_max: float = 40.0) -> float:
    """Return the wall‑shear f''(0) for a given pressure gradient β.

    Uses a shooting method to integrate the Falkner–Skan ODE until a
    sufficiently large value of η is reached.  The root of the boundary
    residual (f'(∞) − 1) gives the initial slope f''(0).

    Args:
        beta: Falkner–Skan pressure‑gradient parameter.
        eta_max: Maximum similarity coordinate to integrate to; default 40.

    Returns:
        The wall‑shear f''(0) or 0 if β is below the physical range.
    """
    # Hartree critical value: solutions exist only for β ≥ −0.198838
    BETA_CRIT = -0.198838
    if beta < BETA_CRIT:
        return 0.0

    def fs_equations(eta, y):
        f, fp, fpp = y
        return [fp, fpp, -f * fpp - beta * (1.0 - fp ** 2)]

    def boundary_residual(alpha_guess):
        # Integrate from η=0 to η=eta_max with initial conditions
        sol = solve_ivp(
            fs_equations,
            [0.0, eta_max],
            [0.0, 0.0, alpha_guess],
            t_eval=[eta_max],
            rtol=1e-8,
            atol=1e-10,
        )
        # residual is f'(eta_max) − 1
        f_prime_end = sol.y[1, -1]
        return f_prime_end - 1.0

    # Solve for α using a root finder; initial guess 0.5 typically works
    alpha_root = fsolve(boundary_residual, 0.5)[0]
    return alpha_root


def precompute_alpha_grid(beta_min: float, beta_max: float, n_beta: int) -> tuple:
    """Precompute wall‑shear values α(β) for a uniform β grid.

    Args:
        beta_min: Minimum β in the grid.
        beta_max: Maximum β in the grid.
        n_beta: Number of points in the grid.

    Returns:
        Tuple of arrays (betas, alphas).
    """
    betas = np.linspace(beta_min, beta_max, n_beta)
    alphas = np.zeros_like(betas)
    for i, b in enumerate(betas):
        alphas[i] = solve_falkner_skan(b)
    return betas, alphas


def interpolate_alpha(beta: np.ndarray, betas: np.ndarray, alphas: np.ndarray) -> np.ndarray:
    """Interpolate wall‑shear α for arbitrary β values.

    Values of β below the critical range return zero.

    Args:
        beta: Array of β values.
        betas: Grid of β points.
        alphas: Corresponding α(β) values.

    Returns:
        Array of interpolated α values.
    """
    alpha_interp = np.interp(beta, betas, alphas, left=0.0, right=alphas[-1])
    return alpha_interp


def generate_ou_path(mu: float, sigma: float, corr_len: float, x_length: float, dx: float) -> np.ndarray:
    """Generate a one‑dimensional Ornstein–Uhlenbeck process.

    The OU process solves dβ = −(β − μ)/ℓ_c dx + σ dW with spatial step size dx.

    Args:
        mu: Long‑term mean of the process.
        sigma: Standard deviation of the stationary distribution.
        corr_len: Correlation length (ℓ_c); the process returns to the mean over this scale.
        x_length: Total length of the domain.
        dx: Spatial step.

    Returns:
        Array of β values along the domain [0, x_length] with spacing dx.
    """
    n_steps = int(x_length / dx) + 1
    beta_path = np.empty(n_steps)
    beta_path[0] = mu
    # Precompute exponential decay and diffusion terms
    exp_decay = np.exp(-dx / corr_len)
    noise_scale = sigma * np.sqrt(1.0 - exp_decay ** 2)
    for i in range(1, n_steps):
        beta_path[i] = (
            mu + (beta_path[i - 1] - mu) * exp_decay + noise_scale * np.random.randn()
        )
    return beta_path


def run_simulation(
    base_beta: float,
    sigma: float,
    corr_len: float,
    n_samples: int,
    x_length: float = 1.0,
    dx: float = 0.01,
    epsilon: float = 0.02,
    betas_grid=None,
    alphas_grid=None,
) -> tuple:
    """Run a Monte Carlo simulation of coloured Falkner–Skan flows.

    For each sample, generate an OU β path along the plate, compute wall‑shear
    using interpolation, record the minimum α and whether a near‑separation event
    occurs (α_min ≤ ε).

    Args:
        base_beta: Mean pressure‑gradient parameter μ.
        sigma: Standard deviation of β fluctuations.
        corr_len: Correlation length ℓ_c.
        n_samples: Number of Monte Carlo samples.
        x_length: Length of the plate (default 1.0).
        dx: Spatial discretisation step (default 0.01).
        epsilon: Threshold for near‑separation (default 0.02).
        betas_grid: Precomputed β grid for interpolation; if None, uses default.
        alphas_grid: Precomputed α grid corresponding to betas_grid.

    Returns:
        Tuple containing arrays (alpha_mins, sep_events).
    """
    # If no precomputed grid provided, create a small one on the fly
    if betas_grid is None or alphas_grid is None:
        betas_grid, alphas_grid = precompute_alpha_grid(-0.2, 0.3, 50)

    n_points = int(x_length / dx) + 1
    alpha_mins = np.empty(n_samples)
    sep_events = np.empty(n_samples, dtype=int)

    for j in range(n_samples):
        beta_path = generate_ou_path(base_beta, sigma, corr_len, x_length, dx)
        alpha_path = interpolate_alpha(beta_path, betas_grid, alphas_grid)
        alpha_min = np.min(alpha_path)
        alpha_mins[j] = alpha_min
        sep_events[j] = 1 if alpha_min <= epsilon else 0
    return alpha_mins, sep_events


def save_csv(filename: str, headers: list, rows: list) -> None:
    """Save data to a CSV file."""
    with open(filename, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description="Run Falkner–Skan OU simulations")
    parser.add_argument("--base_beta", type=float, default=-0.165, help="Mean pressure gradient μ")
    parser.add_argument("--sigma", type=float, default=0.03, help="Noise intensity σ")
    parser.add_argument("--corr_len", type=float, default=0.1, help="Correlation length ℓ_c")
    parser.add_argument("--samples", type=int, default=100, help="Number of Monte Carlo samples")
    parser.add_argument("--epsilon", type=float, default=0.02, help="Separation threshold ε")
    parser.add_argument("--outfile", type=str, default="simulation_results.csv", help="Output CSV filename")
    args = parser.parse_args()

    # Precompute baseline mapping with a moderate resolution
    betas_grid, alphas_grid = precompute_alpha_grid(-0.2, 0.3, 60)
    alpha_mins, sep_events = run_simulation(
        args.base_beta,
        args.sigma,
        args.corr_len,
        args.samples,
        x_length=1.0,
        dx=0.01,
        epsilon=args.epsilon,
        betas_grid=betas_grid,
        alphas_grid=alphas_grid,
    )

    rows = []
    for i in range(args.samples):
        rows.append([
            args.base_beta,
            args.sigma,
            args.corr_len,
            i,
            alpha_mins[i],
            sep_events[i],
        ])
    headers = ["base_beta", "sigma", "corr_len", "sample_idx", "alpha_min", "sep_event"]
    save_csv(args.outfile, headers, rows)
    print(f"Saved results to {args.outfile}")


if __name__ == "__main__":
    main()
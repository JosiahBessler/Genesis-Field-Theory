from __future__ import annotations

import argparse
import itertools
import math
from pathlib import Path

import numpy as np


def mu_star_from_lambda(
    lam_hat: np.ndarray,
    a1: float = 0.0,
    a2: float = 0.0,
    b: float = 1.0,
    c: float = 1.0,
) -> np.ndarray:
    """
    Kernel-vacuum style per-mode minimizer.

    Solve elementwise for mu>0:
        a1 + c*lam_hat + 2*a2*mu + b*(1 + ln(mu)) = 0

    - lam_hat must be finite and >= 0 (dimensionless eigenvalue)
    - b > 0 for entropy-like convexity
    - a2 >= 0 for stability

    Returns strictly positive float64 array.
    """
    lam_hat = np.asarray(lam_hat, dtype=np.float64)
    if b <= 0:
        raise ValueError("b must be > 0.")
    if a2 < 0:
        raise ValueError("a2 must be >= 0.")

    # Initial guess ignoring a2:  mu = exp(-1 - (a1 + c*lam_hat)/b)
    mu = np.exp(-1.0 - (a1 + c * lam_hat) / b)
    mu = np.clip(mu, 1e-300, np.inf)

    if a2 == 0.0:
        return mu

    # Newton iterations (vectorized)
    for _ in range(60):
        f = a1 + c * lam_hat + 2.0 * a2 * mu + b * (1.0 + np.log(mu))
        fp = 2.0 * a2 + b / mu
        step = f / fp
        mu_new = mu - step
        mu_new = np.maximum(mu_new, 1e-300)
        if np.max(np.abs(step)) < 1e-12:
            mu = mu_new
            break
        mu = mu_new

    return mu


def lambda_s3(ell: int, R: float) -> float:
    return (ell * (ell + 2)) / (R * R)


def lambda_t3(nx: int, ny: int, nz: int, Lx: float, Ly: float, Lz: float) -> float:
    twopi2 = (2.0 * math.pi) ** 2
    return twopi2 * ((nx * nx) / (Lx * Lx) + (ny * ny) / (Ly * Ly) + (nz * nz) / (Lz * Lz))


def generate_model_a_lambdas(
    R: float,
    Lx: float,
    Ly: float,
    Lz: float,
    ell_max: int,
    n_max: int,
    lam_max: float | None = None,
) -> np.ndarray:
    """
    Enumerate eigenvalues for Xi = S^3(R) x T^3(Lx,Ly,Lz) under O = -Δ.
    Returns lambdas as float64, sorted.
    """
    lambdas: list[float] = []

    n_range = range(-n_max, n_max + 1)
    n_triples = list(itertools.product(n_range, n_range, n_range))

    for ell in range(0, ell_max + 1):
        lam_s = lambda_s3(ell, R)
        for nx, ny, nz in n_triples:
            # Skip exact zero mode to avoid lambda=0
            if ell == 0 and nx == 0 and ny == 0 and nz == 0:
                continue

            lam_t = lambda_t3(nx, ny, nz, Lx, Ly, Lz)
            lam = lam_s + lam_t

            if lam_max is not None and lam > lam_max:
                continue

            lambdas.append(lam)

    lambdas_arr = np.asarray(lambdas, dtype=np.float64)
    lambdas_arr.sort()
    return lambdas_arr


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="data/model_a.npz")

    # Geometry
    ap.add_argument("--R", type=float, default=1.0, help="S^3 radius")
    ap.add_argument("--Lx", type=float, default=1.0, help="T^3 length in x")
    ap.add_argument("--Ly", type=float, default=1.0, help="T^3 length in y")
    ap.add_argument("--Lz", type=float, default=1.0, help="T^3 length in z")
    ap.add_argument("--ell-max", type=int, default=40)
    ap.add_argument("--n-max", type=int, default=6)
    ap.add_argument("--lam-max", type=float, default=None, help="Optional cutoff in lambda")

    # Kernel-vacuum params (dimensionless)
    ap.add_argument("--a1", type=float, default=0.0)
    ap.add_argument("--a2", type=float, default=0.0)
    ap.add_argument("--b", type=float, default=1.0)
    ap.add_argument("--c", type=float, default=1.0)
    ap.add_argument("--lam-ref", choices=["median", "max"], default="median",
                    help="Reference scale for lambda to form lam_hat=lambda/lam_ref")

    # If you want to keep old placeholder mode, flip this flag:
    ap.add_argument("--placeholder-mu", action="store_true",
                    help="Use the old placeholder degeneracy weight mu = (ell+1) instead of kernel vacuum.")

    args = ap.parse_args()

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    # 1) Enumerate lambdas
    lambdas = generate_model_a_lambdas(
        R=args.R, Lx=args.Lx, Ly=args.Ly, Lz=args.Lz,
        ell_max=args.ell_max, n_max=args.n_max, lam_max=args.lam_max
    )

    if lambdas.size == 0:
        raise RuntimeError("No modes generated; check ell-max/n-max/lam-max.")

    # 2) Build mus
    if args.placeholder_mu:
        # Placeholder: use a constant mu=1 (purely for debugging), or replace with your earlier degeneracy logic.
        mus = np.ones_like(lambdas, dtype=np.float64)
    else:
        lam_ref = float(np.median(lambdas)) if args.lam_ref == "median" else float(np.max(lambdas))
        lam_hat = lambdas / lam_ref
        mus = mu_star_from_lambda(lam_hat, a1=args.a1, a2=args.a2, b=args.b, c=args.c)

    # 3) Hard guarantee: paired arrays
    assert lambdas.shape == mus.shape, f"shape mismatch: lambdas {lambdas.shape} vs mus {mus.shape}"
    assert np.all(np.isfinite(mus)) and np.all(mus > 0), "mu array must be finite and strictly positive."

    np.savez(out, lambdas=lambdas, mus=mus)

    print(f"Wrote {out} with {lambdas.size} modes")
    print(f"lambda range: {lambdas.min():.6g} .. {lambdas.max():.6g}")
    print(f"mu range: {mus.min():.6g} .. {mus.max():.6g}")
    if not args.placeholder_mu:
        print(f"lambda_ref ({args.lam_ref}) = {lam_ref:.6g}, using a1={args.a1}, a2={args.a2}, b={args.b}, c={args.c}")


if __name__ == "__main__":
    main()

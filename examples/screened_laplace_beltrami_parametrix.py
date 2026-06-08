"""Screened Laplace-Beltrami Green kernel by parametrix splitting.

This prototype is the next step after the sphere-only spectral Green operator:
it separates the kernel into a universal local singular model and a smooth
remainder,

    G_alpha(x, y) = Phi_alpha(d(x, y)) + R_alpha(x, y),

where Phi_alpha(r) = K0(sqrt(alpha) r)/(2*pi) is the tangent-plane screened
fundamental solution.  On a curved surface the distance in the singular model
is corrected with the local Gaussian curvature.  The remaining smooth part is
then approximated spectrally.

For validation we use the unit sphere, where the exact screened Green kernel is
known for 0 < alpha < 1/4.  This keeps the example measurable while the
parametrix idea itself is not sphere-specific.
"""

import os
from pathlib import Path

cache_root = Path("/private/tmp/surfgeopy_cache")
cache_root.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("XDG_CACHE_HOME", str(cache_root))
os.environ.setdefault("MPLCONFIGDIR", str(cache_root / "matplotlib"))

import matplotlib
import numpy as np
from numpy.polynomial.chebyshev import chebfit, chebval
from scipy import special

matplotlib.use("Agg")

import matplotlib.pyplot as plt


def exact_sphere_green(cosine: np.ndarray, alpha: float) -> np.ndarray:
    """Exact screened Green kernel on the unit sphere for 0 < alpha < 1/4."""
    if alpha <= 0.0 or alpha >= 0.25:
        raise ValueError("this closed form requires 0 < alpha < 1/4")

    degree = (-1.0 + np.sqrt(1.0 - 4.0 * alpha)) / 2.0
    return -special.lpmv(0, degree, -cosine) / (4.0 * np.sin(np.pi * degree))


def curvature_corrected_distance(
    cosine: np.ndarray,
    gaussian_curvature: float = 1.0,
) -> np.ndarray:
    """Approximate geodesic distance from the chord distance."""
    chord = np.sqrt(np.maximum(2.0 * (1.0 - cosine), np.finfo(float).tiny))
    return chord * (
        1.0
        + gaussian_curvature * chord**2 / 24.0
        + 3.0 * gaussian_curvature**2 * chord**4 / 640.0
    )


def screened_parametrix(
    cosine: np.ndarray,
    alpha: float,
    gaussian_curvature: float = 1.0,
) -> np.ndarray:
    """Curvature-corrected tangent-plane screened singular model."""
    distance = curvature_corrected_distance(cosine, gaussian_curvature)
    return special.k0(np.sqrt(alpha) * distance) / (2.0 * np.pi)


def planar_parametrix(cosine: np.ndarray, alpha: float) -> np.ndarray:
    """Uncorrected tangent-plane screened singular model."""
    chord = np.sqrt(np.maximum(2.0 * (1.0 - cosine), np.finfo(float).tiny))
    return special.k0(np.sqrt(alpha) * chord) / (2.0 * np.pi)


def fit_smooth_remainder(
    alpha: float,
    degree: int,
    n_samples: int = 2000,
    gaussian_curvature: float = 1.0,
) -> np.ndarray:
    """Fit the smooth kernel remainder in the variable x.y."""
    theta = np.linspace(np.pi, 1.0e-4, n_samples)
    cosine = np.cos(theta)
    exact = exact_sphere_green(cosine, alpha)
    singular = screened_parametrix(cosine, alpha, gaussian_curvature)
    return chebfit(cosine, exact - singular, degree)


def reconstructed_green(
    cosine: np.ndarray,
    alpha: float,
    remainder_coefficients: np.ndarray,
    gaussian_curvature: float = 1.0,
) -> np.ndarray:
    singular = screened_parametrix(cosine, alpha, gaussian_curvature)
    smooth = chebval(cosine, remainder_coefficients)
    return singular + smooth


def main() -> None:
    alpha = 0.1
    test_theta = np.geomspace(1.0e-3, np.pi, 1200)
    test_cosine = np.cos(test_theta)
    exact = exact_sphere_green(test_cosine, alpha)

    anchor = np.cos(0.8)
    anchored_exact = exact_sphere_green(anchor, alpha)
    anchored_planar = planar_parametrix(anchor, alpha)
    anchored_curved = screened_parametrix(anchor, alpha)

    planar = (
        planar_parametrix(test_cosine, alpha)
        + anchored_exact
        - anchored_planar
    )
    curved = (
        screened_parametrix(test_cosine, alpha)
        + anchored_exact
        - anchored_curved
    )

    near = test_theta <= 0.25
    planar_near_error = np.max(np.abs(planar[near] - exact[near]))
    curved_near_error = np.max(np.abs(curved[near] - exact[near]))

    print("Screened Laplace-Beltrami parametrix split")
    print("Kernel:                G_alpha = singular parametrix + smooth remainder")
    print(f"alpha:                 {alpha:.3f}")
    print(f"Near-field interval:   theta <= {0.25:.2f}")
    print(f"Planar near max error: {planar_near_error:.3e}")
    print(f"Curved near max error: {curved_near_error:.3e}")
    print("\nremainder degree    max kernel error    max relative error")

    degrees = [0, 2, 4, 6, 8, 10, 12, 16, 20, 30, 40]
    max_errors = []
    relative_errors = []

    for degree in degrees:
        coefficients = fit_smooth_remainder(alpha, degree)
        reconstructed = reconstructed_green(test_cosine, alpha, coefficients)
        error = np.abs(reconstructed - exact)
        max_error = np.max(error)
        max_relative_error = np.max(error / np.maximum(np.abs(exact), 1.0e-30))
        max_errors.append(max_error)
        relative_errors.append(max_relative_error)
        print(f"{degree:16d}    {max_error:16.3e}    {max_relative_error:.3e}")

    root = Path(__file__).resolve().parents[1]
    image_dir = root / "images"
    image_dir.mkdir(exist_ok=True)
    figure_path = image_dir / "screened_laplace_parametrix_remainder.png"

    plt.figure(figsize=(6.5, 4.5))
    plt.semilogy(degrees, max_errors, "-or", label="absolute")
    plt.semilogy(degrees, relative_errors, "-ob", label="relative")
    plt.xlabel("Chebyshev degree of smooth remainder", fontsize=12)
    plt.ylabel("Max kernel error", fontsize=12)
    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(figure_path, dpi=200)

    print(f"\nSaved figure:          {figure_path}")


if __name__ == "__main__":
    main()

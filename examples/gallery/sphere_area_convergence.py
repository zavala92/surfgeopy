"""Generate a sphere-area convergence table and figure."""

from __future__ import annotations

import csv
import os
import sys
from pathlib import Path
from time import perf_counter
from typing import Dict, Union

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from surfgeopy import IntegrationConfig, LevelSetSurface, SurfaceMesh, integrate

OUTPUT_DIR = ROOT / "docs" / "gallery"
MPL_CACHE_DIR = OUTPUT_DIR / ".matplotlib-cache"
MPL_CACHE_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_CACHE_DIR))
os.environ.setdefault("XDG_CACHE_HOME", str(MPL_CACHE_DIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

MESH_PATH = ROOT / "tests" / "mesh_test" / "sphere_N=104.mat"
REFERENCE_AREA = 4.0 * np.pi
Record = Dict[str, Union[float, int]]


def phi(x: np.ndarray) -> float:
    return x[0] ** 2 + x[1] ** 2 + x[2] ** 2 - 1.0


def grad_phi(x: np.ndarray) -> np.ndarray:
    return np.array([2.0 * x[0], 2.0 * x[1], 2.0 * x[2]])


def run_case(interpolation_degree: int) -> Record:
    mesh = SurfaceMesh.from_mat(str(MESH_PATH))
    surface = LevelSetSurface(mesh, phi, grad_phi)
    config = IntegrationConfig(
        interpolation_degree=interpolation_degree,
        refinement_level=1,
        integration_degree=14,
        quadrature_rule="Gauss_Legendre",
    )

    start = perf_counter()
    result = integrate(surface, lambda _: 1.0, config)
    seconds = perf_counter() - start
    absolute_error = abs(result.total - REFERENCE_AREA)
    return {
        "interpolation_degree": interpolation_degree,
        "value": result.total,
        "absolute_error": absolute_error,
        "relative_error": absolute_error / REFERENCE_AREA,
        "seconds": seconds,
        "n_quadrature_points": result.n_quadrature_points,
    }


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    records = [run_case(degree) for degree in (2, 4, 6)]

    csv_path = OUTPUT_DIR / "sphere_area_convergence.csv"
    with open(csv_path, "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=records[0].keys())
        writer.writeheader()
        writer.writerows(records)

    fig, ax = plt.subplots(figsize=(5.5, 3.6))
    ax.semilogy(
        [record["interpolation_degree"] for record in records],
        [record["relative_error"] for record in records],
        marker="o",
        color="#1f6f8b",
        linewidth=2,
    )
    ax.set_xlabel("Interpolation degree")
    ax.set_ylabel("Relative area error")
    ax.set_title("Unit sphere area convergence")
    ax.grid(True, which="both", linewidth=0.5, alpha=0.35)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "sphere_area_convergence.png", dpi=180)

    print(f"Wrote {csv_path}")
    print(f"Wrote {OUTPUT_DIR / 'sphere_area_convergence.png'}")


if __name__ == "__main__":
    main()

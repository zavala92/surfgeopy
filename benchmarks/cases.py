"""Reusable benchmark cases for surfgeopy."""

from dataclasses import dataclass, replace
from pathlib import Path
from time import perf_counter
from typing import Callable, Dict, List, Optional, Sequence, Union

import numpy as np

from surfgeopy import IntegrationConfig, LevelSetSurface, SurfaceMesh, integrate


ROOT = Path(__file__).resolve().parents[1]
BenchmarkValue = Union[float, int, str, None]


@dataclass(frozen=True)
class BenchmarkCase:
    """One reproducible integration benchmark configuration."""

    name: str
    surface: str
    mesh_path: Path
    level_set: Callable[[np.ndarray], float]
    gradient: Callable[[np.ndarray], np.ndarray]
    reference_value: Optional[float]
    interpolation_degree: int
    refinement_level: int
    integration_degree: int
    quadrature_rule: str
    lp_degree: float = float("inf")
    integrand: Callable[[np.ndarray], float] = lambda _: 1.0
    quantity: str = "area"

    def config(self) -> IntegrationConfig:
        return IntegrationConfig(
            interpolation_degree=self.interpolation_degree,
            lp_degree=self.lp_degree,
            refinement_level=self.refinement_level,
            integration_degree=self.integration_degree,
            quadrature_rule=self.quadrature_rule,
        )


def sphere_level_set(x: np.ndarray) -> float:
    return x[0] ** 2 + x[1] ** 2 + x[2] ** 2 - 1.0


def sphere_gradient(x: np.ndarray) -> np.ndarray:
    return np.array([2.0 * x[0], 2.0 * x[1], 2.0 * x[2]])


def torus_level_set(x: np.ndarray, major_radius: float = 2.0, minor_radius: float = 1.0) -> float:
    radius_xy = np.sqrt(x[0] * x[0] + x[1] * x[1])
    return (radius_xy - major_radius) ** 2 + x[2] ** 2 - minor_radius ** 2


def torus_gradient(x: np.ndarray, major_radius: float = 2.0) -> np.ndarray:
    radius_xy = np.sqrt(x[0] * x[0] + x[1] * x[1])
    return np.array([
        -2.0 * major_radius * x[0] / radius_xy + 2.0 * x[0],
        -2.0 * major_radius * x[1] / radius_xy + 2.0 * x[1],
        2.0 * x[2],
    ])


def sphere_area(radius: float = 1.0) -> float:
    return 4.0 * np.pi * radius ** 2


def torus_area(major_radius: float = 2.0, minor_radius: float = 1.0) -> float:
    return (2.0 * np.pi * major_radius) * (2.0 * np.pi * minor_radius)


CASES = (
    BenchmarkCase(
        name="sphere_n104_gl",
        surface="sphere",
        mesh_path=ROOT / "tests" / "mesh_test" / "sphere_N=104.mat",
        level_set=sphere_level_set,
        gradient=sphere_gradient,
        reference_value=sphere_area(),
        interpolation_degree=6,
        refinement_level=1,
        integration_degree=14,
        quadrature_rule="Gauss_Legendre",
    ),
    BenchmarkCase(
        name="sphere_n124_pullback",
        surface="sphere",
        mesh_path=ROOT / "meshes" / "SphereMesh_N=124_r=1.mat",
        level_set=sphere_level_set,
        gradient=sphere_gradient,
        reference_value=sphere_area(),
        interpolation_degree=6,
        refinement_level=0,
        integration_degree=14,
        quadrature_rule="Pull_back_Gauss",
    ),
    BenchmarkCase(
        name="torus_n260_gl",
        surface="torus",
        mesh_path=ROOT / "meshes" / "torus_260.mat",
        level_set=torus_level_set,
        gradient=torus_gradient,
        reference_value=torus_area(),
        interpolation_degree=4,
        refinement_level=0,
        integration_degree=10,
        quadrature_rule="Gauss_Legendre",
    ),
)

SUITES = {
    "quick": ("sphere_n104_gl",),
    "baseline": ("sphere_n104_gl", "torus_n260_gl"),
    "sphere": ("sphere_n104_gl", "sphere_n124_pullback"),
    "torus": ("torus_n260_gl",),
    "all": tuple(case.name for case in CASES),
}


def list_cases() -> List[BenchmarkCase]:
    return list(CASES)


def get_case(name: str) -> BenchmarkCase:
    for case in CASES:
        if case.name == name:
            return case
    raise KeyError(f"Unknown benchmark case: {name}")


def select_cases(suite: str = "quick", names: Sequence[str] = ()) -> List[BenchmarkCase]:
    selected_names = list(names) if names else list(SUITES[suite])
    return [get_case(name) for name in selected_names]


def run_case(case: BenchmarkCase) -> Dict[str, BenchmarkValue]:
    mesh = SurfaceMesh.from_mat(str(case.mesh_path))
    surface = LevelSetSurface(mesh, case.level_set, case.gradient)

    start = perf_counter()
    result = integrate(surface, case.integrand, case.config())
    seconds = perf_counter() - start

    absolute_error = None
    relative_error = None
    if case.reference_value is not None:
        absolute_error = abs(result.total - case.reference_value)
        relative_error = absolute_error / abs(case.reference_value)

    return {
        "name": case.name,
        "surface": case.surface,
        "quantity": case.quantity,
        "mesh": str(case.mesh_path.relative_to(ROOT)),
        "n_vertices": mesh.n_vertices,
        "n_faces": mesh.n_faces,
        "interpolation_degree": case.interpolation_degree,
        "lp_degree": "inf" if case.lp_degree == float("inf") else case.lp_degree,
        "refinement_level": case.refinement_level,
        "integration_degree": case.integration_degree,
        "quadrature_rule": case.quadrature_rule,
        "value": result.total,
        "reference_value": case.reference_value,
        "absolute_error": absolute_error,
        "relative_error": relative_error,
        "seconds": seconds,
        "n_quadrature_points": result.points.shape[0],
    }


def clone_case(case: BenchmarkCase, **overrides: object) -> BenchmarkCase:
    """Return a copy of a case with selected fields changed."""
    return replace(case, **overrides)


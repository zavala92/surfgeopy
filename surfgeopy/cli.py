"""Command-line helpers for installed ``surfgeopy`` environments."""

from __future__ import annotations

import argparse
import importlib.metadata
import math
import platform
import sys
from typing import Sequence

from .version import version as __version__


_RUNTIME_DEPENDENCIES = (
    "numpy",
    "scipy",
    "numba",
    "minterpy",
    "modepy",
    "matplotlib",
)


def _dependency_status(name: str) -> str:
    try:
        version = importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "missing"
    return f"ok {version}"


def _run_doctor() -> int:
    print(f"surfgeopy {__version__}")
    print(f"python {platform.python_version()}")
    print(f"executable {sys.executable}")
    print("dependencies")
    for dependency in _RUNTIME_DEPENDENCIES:
        print(f"  {dependency}: {_dependency_status(dependency)}")
    return 0


def _run_demo(args: argparse.Namespace) -> int:
    from .api import IntegrationConfig, LevelSetSurface, integrate

    surface = LevelSetSurface.sphere(
        mesh_refinement_level=args.mesh_refinement_level,
        radius=args.radius,
    )
    config = IntegrationConfig(
        interpolation_degree=args.interpolation_degree,
        integration_degree=args.integration_degree,
        quadrature_rule=args.quadrature_rule,
    )
    result = integrate(surface, lambda _: 1.0, config)
    exact = 4.0 * math.pi * args.radius**2

    print("surfgeopy sphere demo")
    print(f"area        {result.total:.12f}")
    print(f"exact       {exact:.12f}")
    print(f"abs error   {abs(result.total - exact):.3e}")
    print(f"quad points {result.n_quadrature_points}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="surfgeopy",
        description="Smoke-test and demo commands for surfgeopy installations.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"surfgeopy {__version__}",
    )

    subparsers = parser.add_subparsers(dest="command")

    doctor = subparsers.add_parser(
        "doctor",
        help="Show Python and dependency versions without importing numerical backends.",
    )
    doctor.set_defaults(func=lambda _: _run_doctor())

    demo = subparsers.add_parser(
        "demo",
        help="Integrate 1 over a built-in sphere and print the area error.",
    )
    demo.add_argument(
        "--mesh-refinement-level",
        type=int,
        default=1,
        help="Icosphere refinement level for the reference mesh.",
    )
    demo.add_argument(
        "--radius",
        type=float,
        default=1.0,
        help="Sphere radius.",
    )
    demo.add_argument(
        "--interpolation-degree",
        type=int,
        default=4,
        help="Polynomial interpolation degree.",
    )
    demo.add_argument(
        "--integration-degree",
        type=int,
        default=8,
        help="Quadrature degree.",
    )
    demo.add_argument(
        "--quadrature-rule",
        default="Gauss_Legendre",
        help="Reference quadrature rule name.",
    )
    demo.set_defaults(func=_run_demo)

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if not hasattr(args, "func"):
        parser.print_help()
        return 0
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())

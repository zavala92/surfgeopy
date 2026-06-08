"""Public package API for surfgeopy.

The package keeps top-level imports lazy so ``import surfgeopy`` does not
eagerly import the numerical backends used by specific workflows.
"""

from importlib import import_module

from .version import version as __version__

_SUBMODULES = {
    "api",
    "cli",
    "mesh_io",
    "quadrature_points",
    "quadrature_points_gl",
    "reference_quadrature",
    "remesh",
    "spectral",
    "surf_integration",
    "surface",
    "utils",
}

_SYMBOL_TO_MODULE = {
    # High-level API
    "SurfaceMesh": "api",
    "LevelSetSurface": "api",
    "IntegrationConfig": "api",
    "IntegrationResult": "api",
    "CompiledIntegrationScheme": "api",
    "DiagnosticIntegrationResult": "api",
    "IndicatorRefinementIteration": "api",
    "IndicatorRefinementResult": "api",
    "PatchGeometry": "api",
    "SurfaceGeometryResult": "api",
    "compile_integration": "api",
    "integrate": "api",
    "integrate_with_diagnostics": "api",
    "refine_by_indicator": "api",
    "surface_geometry": "api",
    # Reference quadrature
    "PULL_BACK_GAUSS": "reference_quadrature",
    "GAUSS_LEGENDRE": "reference_quadrature",
    "MODEPY_XIAO_GIMBUTAS": "reference_quadrature",
    "MODEPY_GRUNDMANN_MOELLER": "reference_quadrature",
    "MODEPY_VIOREANU_ROKHLIN": "reference_quadrature",
    "MODEPY_SIMPLEX_RULES": "reference_quadrature",
    "DEFAULT_QUADRATURE_RULE": "reference_quadrature",
    "ReferenceQuadrature": "reference_quadrature",
    "make_reference_quadrature": "reference_quadrature",
    # Surface projection helpers
    "ImplicitSurface": "surface",
    "ProjectionResult": "surface",
    "simplex_barycentric_coordinates": "surface",
    "affine_triangle_points": "surface",
    "project_triangle_nodes": "surface",
    # Legacy integration and quadrature helpers
    "integration": "surf_integration",
    "accumulate_surface_integrals": "surf_integration",
    "compute_surf_quadrature": "surf_integration",
    "compute_surf_geometry": "surf_integration",
    "quadrature_surf_tri": "surf_integration",
    "quadrature_split_surf_tri": "surf_integration",
    "quadrule_on_flat": "quadrature_points",
    "quadrule_on_simplex": "quadrature_points",
    "gauss_legendre_square": "quadrature_points_gl",
    "q_gauss_legendre": "quadrature_points_gl",
    # Remeshing helpers
    "subdivide": "remesh",
    "subdivide_conforming": "remesh",
    "faces_to_edges": "remesh",
    "unique_rows": "remesh",
    "hashable_rows": "remesh",
    "unique_ordered": "remesh",
    # Utility helpers
    "compute_norm": "utils",
    "_cross": "utils",
    "max_edge_length": "utils",
    "decimal_to_digits": "utils",
    "float_to_int": "utils",
    "pushforward": "utils",
    "pullback": "utils",
    "SimpleImplicitSurfaceProjection": "utils",
    "read_mesh_data": "utils",
    # Mesh I/O helpers
    "read_gmsh_mesh": "mesh_io",
    "write_gmsh_mesh": "mesh_io",
    # FFT-backed tensor Chebyshev validation helpers
    "chebyshev_lobatto_nodes": "spectral",
    "chebyshev_coefficients_2d": "spectral",
    "chebyshev_values_2d": "spectral",
    "chebyshev_derivative_coefficients_2d": "spectral",
    "dense_chebyshev_coefficients_2d": "spectral",
}

__all__ = ["__version__", *_SYMBOL_TO_MODULE]


def __getattr__(name):
    """Load public objects and submodules on first access."""
    if name in _SUBMODULES:
        module = import_module(f".{name}", __name__)
        globals()[name] = module
        return module

    module_name = _SYMBOL_TO_MODULE.get(name)
    if module_name is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

    module = import_module(f".{module_name}", __name__)
    value = getattr(module, name)
    globals()[name] = value
    return value


def __dir__():
    """Return public names without importing every submodule."""
    return sorted(set(globals()) | _SUBMODULES | set(__all__))

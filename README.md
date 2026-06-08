# Surfgeopy

[![License](https://img.shields.io/github/license/zavala92/surfgeopy?style=flat-square)](https://choosealicense.com/licenses/mit/)
[![Documentation Status](https://readthedocs.org/projects/surfgeopy/badge/?version=latest)](https://surfgeopy.readthedocs.io/en/latest/?badge=latest)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg?style=flat-square)](https://www.python.org/downloads/release/python-3100/)

`surfgeopy` helps you integrate scalar functions over smooth implicit surfaces.
You provide a level-set function, its gradient, and a triangular reference
mesh; `surfgeopy` builds high-order curved surface patches and returns
quadrature values, points, weights, and diagnostics.

The main workflow is high-order integration of smooth scalar functions on
implicit surfaces.

## Install

```bash
pip install surfgeopy
```

Check the installation and run a tiny built-in demo:

```bash
surfgeopy --version
surfgeopy doctor
surfgeopy demo
```

For local development, install the test extras:

```bash
pip install -e ".[test]"
```

## First Example

This example works after a normal package install; it uses the built-in
icosphere reference mesh and integrates `1` over the unit sphere.

```python
import numpy as np

from surfgeopy import IntegrationConfig, LevelSetSurface, integrate

surface = LevelSetSurface.unit_sphere(mesh_refinement_level=1)
config = IntegrationConfig(
    interpolation_degree=4,
    integration_degree=8,
    quadrature_rule="Gauss_Legendre",
)

result = integrate(surface, lambda _: 1.0, config)

print(f"area          = {result.total:.12f}")
print(f"exact 4*pi    = {4.0 * np.pi:.12f}")
print(f"abs error     = {abs(result.total - 4.0 * np.pi):.3e}")
print(f"quad points   = {result.n_quadrature_points}")
```

For a custom surface, replace `LevelSetSurface.unit_sphere(...)` with
`LevelSetSurface(mesh, phi, grad_phi)`.

## Overview

`surfgeopy` starts from a linear triangulated reference mesh and an implicit
surface representation

```text
Gamma = { x in R^3 : phi(x) = 0 }.
```

For each reference triangle, interpolation nodes are mapped to the surface by
closest-point projection. The resulting curved patch is represented by a
Minterpy Newton interpolant on a square parameter domain. Surface integrals are
then evaluated by high-order quadrature on the interpolated patch.

The numerical trajectory is:

```text
smooth high-order surface quadrature
-> reusable geometry and quadrature diagnostics
-> robust mesh and high-order interoperability
-> future surface-analysis workflows on implicit manifolds
```

## Core Idea

The central method is a cubical reparametrization of each surface triangle. For
a reference triangle `T_i`, let

```text
tau_i : Delta_2 -> T_i
pi_i  : T_i -> Gamma_i
sigma : square_2 -> Delta_2
```

where `pi_i` is the closest-point projection to the implicit surface and
`sigma` is the square-squeezing map. The curved surface patch is represented by

```text
varphi_i = pi_i o tau_i o sigma : square_2 -> Gamma_i.
```

This pulls the interpolation task from a triangle to the square. On the square,
`surfgeopy` uses Chebyshev-Lobatto nodes through Minterpy, which gives a stable
high-order tensor-product interpolation setting and avoids the instability of
poorly chosen interpolation nodes.

After the geometry map has been interpolated, the surface integral is evaluated
as

```text
int_Gamma f dS
  ~= sum_i int_square f(varphi_i(u))
       sqrt(det(D varphi_i(u)^T D varphi_i(u))) du.
```

## What Surfgeopy Currently Supports

Implemented and tested core functionality:

- implicit surfaces given by `phi` and `grad_phi`,
- triangulated reference meshes loaded from MATLAB `.mat` files,
- closest-point projection from reference mesh nodes to the implicit surface,
- square-squeezing and pullback/pushforward transformations,
- Minterpy/Chebyshev-Lobatto interpolation of curved surface patches,
- high-order quadrature for smooth surface integrals,
- reusable compiled quadrature schemes for faster repeated integrations,
- simplex quadrature rules from ModePy and tensor-product Gauss-Legendre rules,
- per-face quadrature values, points, weights, and offsets,
- diagnostic integration by comparing base and enriched configurations,
- indicator-based conforming reference-mesh refinement,
- differential geometry samples from the interpolated surface map,
- FFT-backed tensor Chebyshev validation utilities for square-patch
  experiments.

## Geometry And Quadrature Pipeline

The regular quadrature pipeline is:

1. Load a `SurfaceMesh`.
2. Define a `LevelSetSurface` with `phi` and `grad_phi`.
3. Select an `IntegrationConfig`.
4. Build a high-order curved patch on each reference triangle.
5. Evaluate quadrature points, weights, and per-face contributions.
6. Return an `IntegrationResult`.

```python
import numpy as np

from surfgeopy import IntegrationConfig, LevelSetSurface, SurfaceMesh, integrate


def phi(x: np.ndarray) -> float:
    return x[0] ** 2 + x[1] ** 2 + x[2] ** 2 - 1.0


def grad_phi(x: np.ndarray) -> np.ndarray:
    return np.array([2.0 * x[0], 2.0 * x[1], 2.0 * x[2]])


mesh = SurfaceMesh.icosphere(refinement_level=1)
surface = LevelSetSurface(mesh, phi, grad_phi)
config = IntegrationConfig(
    interpolation_degree=4,
    integration_degree=8,
    quadrature_rule="Gauss_Legendre",
)

result = integrate(surface, lambda _: 1.0, config)
print(result.total)  # approximately 4*pi
```

For repeated integrations on the same surface and configuration, precompute the
quadrature once and reuse it:

```python
from surfgeopy import compile_integration

scheme = compile_integration(surface, config)

area = scheme.integrate(lambda _: 1.0)
z_second_moment = scheme.integrate(lambda x: x[2] ** 2)

fast_z_second_moment = scheme.integrate(
    lambda points: points[:, 2] ** 2,
    vectorized=True,
)
```

The same interpolated surface map can be differentiated to inspect geometry:

```python
from surfgeopy import surface_geometry

geometry = surface_geometry(surface, config)
print(geometry.normal)
print(geometry.metric_tensor)
print(geometry.gaussian_curvature)
```

Available quadrature rules include:

- `"ModePy_VioreanuRokhlin"`: default Vioreanu-Rokhlin simplex rule from
  ModePy, mapped to the reference triangle and pulled back through
  square-squeezing.
- `"ModePy_XiaoGimbutas"`: Xiao-Gimbutas simplex rule from ModePy, useful for
  higher requested degrees when Vioreanu-Rokhlin rules are unavailable.
- `"ModePy_GrundmannMoeller"`: Grundmann-Moeller simplex rule from ModePy.
- `"Pull_back_Gauss"`: built-in simplex rule pulled back through
  square-squeezing.
- `"Gauss_Legendre"`: tensor-product Gauss-Legendre quadrature on the square.

## Examples

The examples are organized conceptually in `examples/README.md`.

Core geometry and quadrature:

- `examples/quickstart_sphere_area.py`
- `examples/Test_integration_on_the_whole_sphere.ipynb`
- `examples/Test_integration_on_torus.ipynb`
- `examples/Gauss_Bonnet_theorem_bench.ipynb`
- `examples/adapted_mesh_Gauss_Bonnet_benchmark.ipynb`

## Current Limitations

- The high-order geometry pipeline assumes an implicit representation and a
  usable gradient for closest-point projection.
- Adaptive refinement is available for the reference mesh, but fully automatic
  accuracy control is still under development.
- The package currently focuses on smooth scalar surface integration rather
  than general surface PDE solvers.

## Research Direction

The intended research direction is:

```text
implicit surface geometry
+ high-order cubical reparametrization
+ adaptive diagnostics
= reliable smooth-surface quadrature workflows
```

Near-term development should focus on:

- stronger convergence studies on analytic surfaces,
- better diagnostics and automatic configuration guidance,
- broader mesh interoperability,
- careful performance profiling before large-scale solvers are introduced.

## Install

Install from PyPI:

```bash
pip install surfgeopy
```

For development from a local checkout:

```bash
git clone https://github.com/zavala92/surfgeopy.git
cd surfgeopy
pip install -e ".[test]"
```

Avoid `python setup.py install`; modern `pip` builds from `pyproject.toml`.

Documentation: <https://surfgeopy.readthedocs.io>

For MATLAB-related notes, see [README_MATLAB.md](./README_MATLAB.md).

## Testing

Run the unit tests with:

```bash
pytest
```

Run benchmark entry points with:

```bash
python -m benchmarks.run_benchmarks --list
python -m benchmarks.run_benchmarks
```

## Credits And Contributors

This work was partly funded by the Center for Advanced Systems Understanding
(CASUS), financed by Germany's Federal Ministry of Education and Research
(BMBF), and by the Saxony Ministry for Science, Culture and Tourism (SMWK)
with tax funds on the basis of the budget approved by the Saxony State
Parliament.

Main code development:

- Gentian Zavalani (HZDR/CASUS) <gentian.zavalani@tu-dresden.de>

Mathematical foundation:

- Gentian Zavalani (HZDR/CASUS) <gentian.zavalani@tu-dresden.de>
- Oliver Sander (TU Dresden) <oliver.sander@tu-dresden.de>
- Michael Hecht (HZDR/CASUS) <m.hecht@hzdr.de>

Acknowledgement:

- Minterpy development team

## Reference

If you use `surfgeopy` in a program or publication, please cite:

```bibtex
@article{zavalani2025high,
  title={High-Order Integration on Regular Triangulated Manifolds Reaches Superalgebraic Approximation Rates Through Cubical Reparametrizations},
  author={Zavalani, Gentian and Sander, Oliver and Hecht, Michael},
  journal={SIAM Journal on Numerical Analysis},
  volume={63},
  number={6},
  pages={2454--2482},
  year={2025},
  publisher={SIAM}
}
```

## License

[MIT](LICENSE)

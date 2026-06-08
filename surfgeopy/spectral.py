"""FFT-backed tensor Chebyshev utilities for square patch experiments."""

from __future__ import annotations

import numpy as np
from numpy.polynomial.chebyshev import chebder, chebvander
from scipy.fft import dctn

__all__ = [
    "chebyshev_lobatto_nodes",
    "chebyshev_coefficients_2d",
    "chebyshev_values_2d",
    "chebyshev_derivative_coefficients_2d",
    "dense_chebyshev_coefficients_2d",
]


def chebyshev_lobatto_nodes(degree: int) -> np.ndarray:
    """Return Chebyshev-Lobatto nodes on ``[-1, 1]``."""
    if degree < 1:
        raise ValueError("degree must be at least 1")
    return np.cos(np.pi * np.arange(degree + 1) / degree)


def chebyshev_coefficients_2d(values: np.ndarray) -> np.ndarray:
    """Compute tensor Chebyshev coefficients using FFT-backed DCT-I.

    Parameters
    ----------
    values
        Samples on a square Chebyshev-Lobatto tensor grid. The first two axes
        must have equal length ``degree + 1``. Additional trailing axes are
        treated as vector components and transformed independently.
    """
    values = np.asarray(values, dtype=float)
    if values.ndim < 2:
        raise ValueError("values must have at least two dimensions")
    if values.shape[0] != values.shape[1]:
        raise ValueError("first two axes must form a square tensor grid")

    degree = values.shape[0] - 1
    if degree < 1:
        raise ValueError("grid degree must be at least 1")

    coefficients = dctn(values, type=1, axes=(0, 1)) / (degree * degree)
    coefficients[0, ...] *= 0.5
    coefficients[-1, ...] *= 0.5
    coefficients[:, 0, ...] *= 0.5
    coefficients[:, -1, ...] *= 0.5
    return coefficients


def chebyshev_values_2d(
    coefficients: np.ndarray,
    u: np.ndarray,
    v: np.ndarray,
) -> np.ndarray:
    """Evaluate a tensor Chebyshev expansion at ``(u, v)`` points."""
    coefficients = np.asarray(coefficients, dtype=float)
    if coefficients.ndim < 2:
        raise ValueError("coefficients must have at least two dimensions")
    if coefficients.shape[0] < 1 or coefficients.shape[1] < 1:
        raise ValueError("coefficients must have non-empty polynomial axes")

    u = np.asarray(u, dtype=float)
    v = np.asarray(v, dtype=float)
    if u.shape != v.shape:
        raise ValueError("u and v must have the same shape")

    flat_u = u.reshape(-1)
    flat_v = v.reshape(-1)
    basis_u = chebvander(flat_u, coefficients.shape[0] - 1)
    basis_v = chebvander(flat_v, coefficients.shape[1] - 1)
    values = np.einsum("pi,pj,ij...->p...", basis_u, basis_v, coefficients)
    return values.reshape(u.shape + coefficients.shape[2:])


def chebyshev_derivative_coefficients_2d(
    coefficients: np.ndarray,
    axis: int,
) -> np.ndarray:
    """Return tensor Chebyshev coefficients for a first derivative."""
    if axis not in (0, 1):
        raise ValueError("axis must be 0 or 1")
    return chebder(np.asarray(coefficients, dtype=float), axis=axis)


def dense_chebyshev_coefficients_2d(values: np.ndarray) -> np.ndarray:
    """Compute tensor Chebyshev coefficients by dense linear solve.

    This is intended for validation and benchmarking against
    :func:`chebyshev_coefficients_2d`.
    """
    values = np.asarray(values, dtype=float)
    if values.ndim < 2:
        raise ValueError("values must have at least two dimensions")
    if values.shape[0] != values.shape[1]:
        raise ValueError("first two axes must form a square tensor grid")

    degree = values.shape[0] - 1
    if degree < 1:
        raise ValueError("grid degree must be at least 1")

    nodes = chebyshev_lobatto_nodes(degree)
    uu, vv = np.meshgrid(nodes, nodes, indexing="ij")
    basis_u = chebvander(uu.reshape(-1), degree)
    basis_v = chebvander(vv.reshape(-1), degree)
    vandermonde = np.einsum("pi,pj->pij", basis_u, basis_v).reshape(
        (degree + 1) ** 2,
        (degree + 1) ** 2,
    )
    rhs = values.reshape((degree + 1) ** 2, -1)
    coefficients = np.linalg.solve(vandermonde, rhs)
    return coefficients.reshape((degree + 1, degree + 1) + values.shape[2:])

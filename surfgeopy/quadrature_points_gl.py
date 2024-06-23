# gauss_legendre.py
"""
gauss_legendre.py
-----------------
Gauss-Legendre quadrature and cubature rules implementation.
"""
import numpy as np
from numpy import linalg as LA
from .utils import *

__all__ = ['gauss_legendre_square', 'q_gauss_legendre']

def gauss_legendre_square(deg):
    """
    Gauss-Legendre cubature rule on the square [-1, 1]^2.

    Parameters
    ----------
    deg : int
        Algebraic degree of precision of the rule.

    Returns
    -------
    weights_ps : ndarray
        n-by-1 array of single or double, weights of quadrature points.
    quad_ps : ndarray
        n-by-2 array of single or double, natural coordinates of quadrature points.
    """
    
    # Generate 1D Gauss-Legendre quadrature points and weights
    q_points, q_weights = q_gauss_legendre(deg, np.array([-1, 1]))
    
    # Create a 2D grid of quadrature points
    q_points_meshgrid = np.meshgrid(q_points, q_points, indexing='xy')
    q2x_points, q2y_points = map(lambda x: x.reshape(-1), q_points_meshgrid)

    # Compute the weights for the 2D quadrature points
    weights_ps = np.kron(q_weights, q_weights)
    
    # Combine the x and y coordinates into a single array
    quad_ps = np.column_stack((q2x_points, q2y_points))

    return weights_ps, quad_ps


def q_gauss_legendre(n, Domain=np.array([-1, 1])):
    """
    Generate quadrature weights and nodes for Gauss-Legendre quadrature rule.

    An n-point Gauss-Legendre quadrature rule is of degree 2n-1; 
    that is, it integrates all polynomials up to degree 2n-1 exactly.

    Parameters
    ----------
    n : int
        Number of sample points and weights. It must be >= 1.

    Domain : ndarray, optional
        Refers to the domain of integration [a, b].
        The default domain is the bi-unit interval [-1, 1].

    Returns
    --------
    x : ndarray
        1-D ndarray containing the sample points.

    w : ndarray
        1-D ndarray containing the weights.

    Notes
    ------
    The function constructs the symmetric tridiagonal matrix and computes
    the eigenvalues and eigenvectors to obtain the quadrature points and weights.

    The next steps ensure that zero is zero, and the points and weights are symmetric.
    
    If the integration domain [a, b] is different from the default [-1, 1], 
    the points and weights are transformed accordingly.
    """
    # Construct the symmetric tridiagonal matrix
    A = np.divide(0*(np.arange(0, n, 1))/(np.arange(0, n, 1)+1), (2*np.arange(0, n, 1)+1)/(np.arange(0, n, 1)+1))
    B = np.sqrt(np.divide(np.arange(1, n, 1)/(np.arange(1, n, 1)+1), np.multiply((2*np.arange(0, n-1, 1)+1)/(np.arange(0, n-1, 1)+1), (2*np.arange(1, n, 1)+1)/(np.arange(1, n, 1)+1))))
    J = np.diagflat(B, 1) + np.diagflat(A) + np.diagflat(B, -1)

    # Compute the eigenvalues and eigenvectors
    q_points, v = LA.eigh(J)
    q_weights = np.zeros(int(len(q_points)))

    for i in range(int(len(q_points))):
        q_weights[i] = np.transpose(2*np.multiply(v[0, i], v[0, i]))

    # Ensure zero is zero and the points and weights are perfectly symmetric
    q_points[np.abs(q_points) < 10*np.finfo(float).eps] = 0
    q_points[int(np.ceil(q_points[-1]/2)+1):int(q_points[-1])] = -np.flipud(q_points[1:int(np.floor(q_points[-1]/2))])
    q_weights[int(np.ceil(q_weights[-1]/2)+1):int(q_weights[-1])] = np.flipud(q_weights[1:int(np.floor(q_weights[-1]/2))])

    # Transformation of points and weights if [a, b]~=[-1, 1]
    if not np.array_equal(Domain, np.array([-1, 1])):
        q_points = (Domain[1] - Domain[0])/2 * q_points + (Domain[0] + Domain[1])/2
        q_weights = (Domain[1] - Domain[0])/2 * q_weights

    return q_points, q_weights
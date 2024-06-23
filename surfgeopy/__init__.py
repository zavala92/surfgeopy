"""
This is the curved_integration package initialization file.

It imports all necessary submodules and defines the package's public API.
"""

# Import the version of the package
from .version import version as __version__

# Import all public symbols from submodules
from .surf_integration import *
from .remesh import *
from .quadrature_points import *
from .quadrature_points_gl import *
from .utils import *

__all__ = [
    "__version__"
] + (
    surf_integration.__all__ +
    remesh.__all__ +
    quadrature_points.__all__ +
    quadrature_points_gl.__all__ +
    utils.__all__
)

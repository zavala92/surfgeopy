Installation
============

``surfgeopy`` is designed to install from PyPI with ``pip``.

Requirements
------------

* Python 3.10 or newer
* ``numpy``
* ``scipy``
* ``numba``
* ``matplotlib``
* ``minterpy``
* ``modepy``
* ``pytest`` for tests when installing the ``test`` extra

Install From PyPI
-----------------

Create and activate a virtual environment before installing the package:

.. code-block:: bash

    python -m venv .venv
    source .venv/bin/activate

Install the package:

.. code-block:: bash

    pip install surfgeopy

Check the installed command-line entry point:

.. code-block:: bash

    surfgeopy --version
    surfgeopy doctor
    surfgeopy demo

Install From Source
-------------------

Use ``git`` to obtain the source code:

.. code-block:: bash

    git clone https://github.com/zavala92/surfgeopy.git
    cd surfgeopy

Install in editable mode with test tools:

.. code-block:: bash

    pip install -e ".[test]"

Run Tests
---------

After installation, run the test suite:

.. code-block:: bash

    pytest

The package depends on ``minterpy`` for polynomial interpolation and on
``modepy`` for the optional simplex quadrature families exposed through
``IntegrationConfig``.

The ``-e`` argument creates an editable install. Changes made in the source
tree are reflected when importing ``surfgeopy`` from the same environment.

.. warning::

    Avoid using the command ``python setup.py install`` to install
    ``surfgeopy``. Use ``pip install .`` or ``pip install -e ".[test]"`` so
    the build goes through ``pyproject.toml``.

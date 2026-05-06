Installation
============

``surfgeopy`` is currently installed from source.

Requirements
------------

* Python 3.10 or newer
* ``numpy``
* ``scipy``
* ``numba``
* ``matplotlib``
* ``minterpy``
* ``modepy``
* ``pytest`` for tests

Install From Source
-------------------

Use ``git`` to obtain the source code:

.. code-block:: bash

    git clone https://github.com/zavala92/surfgeopy.git

Create and activate a virtual environment before installing the package:

.. code-block:: bash

    python -m venv .venv
    source .venv/bin/activate

Install in editable mode:

.. code-block:: bash

    pip install -e .

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

    Avoid using the command ``python setup.py install`` to install ``surfgeopy``. This method is discouraged, as the presence of the ``setup.py`` file cannot be guaranteed in the ongoing development of the ``surfgeopy`` library.

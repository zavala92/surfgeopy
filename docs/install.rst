🛠️ Installation
===============


For this prototype implementation, we recommend installing `surfgeopy` by self-building from the source. Use the `git` version control system to obtain the source code:

.. code-block:: bash

    git clone https://codebase.helmholtz.cloud/interpol/surfgeopy.git

.. caution::

    **Switch to Your Virtual Environment:**

    Prior to installation, activate the desired virtual environment using either `conda` or `venv`.

Once inside the activated environment, proceed with the installation using [pip]:

.. code-block:: bash

    pip install -e .

The use of the `-e` argument ensures the installation creates symbolic links. This allows any modifications made by the user within the source folders to be reflected in the installed version when importing modules.

.. warning::

    Avoid using the command ``python setup.py install`` to install ``surfgeopy``. This method is discouraged, as the presence of the ``setup.py`` file cannot be guaranteed in the ongoing development of the ``surfgeopy`` library.

.. _install:

Installing orphans
====================

At some point you will be able to install `orphans` with `pip` or from source.
`orphans` depends upon `orphans`, call it even just an extension if you wish.


Installing with pip
-------------------

This is will work one day, believe me.

.. code-block:: bash

    pip install orphans



Installing from source
-----------------------

This is the recommended installation procedure in a conda environment.
It shall take care of the `afterglowpy` dependency cleanly.

.. code-block:: bash

    git clone git@gitlab.in2p3.fr:johan-bregeon/orphans.git
    cd orphans
    conda env create -f environment.yml
    conda activate orphans
    pip install -r requirements.txt
    pip install -e .


simphony-kratos
===============

The Kratos-CFD and Kratos-DEMPack engine-wrappers for the SimPhoNy framework (www.simphony-project.com)

.. image:: https://travis-ci.org/simphony/simphony-kratos.svg?branch=master
	:target: https://travis-ci.org/simphony/simphony-kratos
		: Build Status

Repository
----------

simphony-kratos is hosted on github: https://github.com/simphony/simphony-kratos

Requirements
------------

- `simphony-common`_ == 0.1.3

.. _simphony-common: https://github.com/simphony/simphony-common

Installation
------------

Installation is based on setuptools and requeires Python 2.7::

    # build and install
    python setup.py install

or::

    # build for in-place development
    python setup.py develop

Kratos Installation
~~~~~~~~~~~~~~~~~~~

Both Kratos-CFD and DEMPack wrappers use KratosMultiphysics. A suitable version to run with these wrappers can be downloaded here:

- `KratosMultiphysics`_

.. _KratosMultiphysics: https://web.cimne.upc.edu/users/croig/data/kratos-simphony.tgz

Once downloaded and extracted in $EXTRACT_DIR a couple of system variables need to be defined::

    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$EXTRACT_DIR/kratos/libs
    export PYTHONPATH=$PYTHONPATH:$EXTRACT_DIR/kratos
    
Usage
-----

To use any of the wrappers after the installation they need to be imported as plugins.

For Kratos-CFD::

  from simphony.engine import kratosCFD
    engine = kratos.CFDWrapper()
    
or for DEMPack::

  from simphony.engine import kratosDEM
    engine = kratos.DEMWrapper()

Testing
-------

To run the full test-suite run::

    python -m unittest discover


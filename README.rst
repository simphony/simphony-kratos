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

- `simphony-common`_ ~= 0.5.0

.. _simphony-common: https://github.com/simphony/simphony-common

Installation
------------

Installation is based on setuptools and requeire Python 2.7::

    # build and install
    python setup.py install

or::

    # build for in-place development
    python setup.py develop

Kratos Installation
~~~~~~~~~~~~~~~~~~~

Both Kratos-CFD and DEMPack wrappers use KratosMultiphysics. A suitable version to run with these wrappers can be downloaded here:

- `KratosMultiphysics`_

.. _KratosMultiphysics: https://github.com/KratosMultiphysics/Kratos/releases/download/v5.0-Simphony/kratos-simphony.tgz

This version is prepared for Ubuntu-12-04-LTS 64-bit. See http://kratos-wiki.cimne.upc.edu/index.php/Download for additional information on getting Kratos for a specific system.

Once downloaded and extracted in $EXTRACT_DIR a couple of system variables need to be defined::

    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$EXTRACT_DIR/kratos/libs
    export PYTHONPATH=$PYTHONPATH:$EXTRACT_DIR/kratos
    
Usage
-----

To use any of the wrappers after the installation they need to be imported as plugins.

For Kratos-CFD::

  from simphony.api import CUDS, Simulation
  sim = Simulation(cuds, "KRATOS_CFD", engine_interface=True)

or for DEMPack::

  from simphony.api import CUDS, Simulation
  sim = Simulation(cuds, "KRATOS_DEM", engine_interface=True)

Testing
-------

To run the full test-suite run::

    python -m unittest discover


""" This files provide a test example on how to use
Simphony along with the Kratos Multiphisics CFD wrapper

"""

import os

from simphony.core.cuba import CUBA

from simphony.api import CUDS, Simulation
from simphony.cuds.meta import api

# TODO: Utils now belong to probably another package
from simphony.engine import EngineInterface
from simphony.engine import kratos

# Add the problem path to the script
path = os.path.join(os.path.dirname(__file__), "CFD_exampleFluid")

cuds = CUDS(name='example_kratos_cfd_simulatiob')

# Integration time:
itime = api.IntegrationTime(name="cfd_integration_time")
itime.time = 0.0001
itime.step = 0.0025
itime.final = 0.0125  # 5 Kratos Timesteps
cuds.add(itime)

# Utils are used to read an existing Kratos model as raw data so we can
# initialize the correct simphony datasets
utils = kratos.CFD_Utils()

# Reads Kratos mpda as a simphony data.
kratos_model = utils.read_modelpart(path)

# Add the datasets readed from the conversor.
for dataset in kratos_model['datasets']:
    cuds.add(dataset)

# Add the boundary contitions from the conversor
for condition in kratos_model['conditions']:
    cuds.add(condition)

# Add the materials contitions from the conversor
for material in kratos_model['materials']:
    cuds.add(material)

# Create the simulation and run the problem
sim = Simulation(cuds, "KRATOS_CFD", engine_interface=True)
sim.run()

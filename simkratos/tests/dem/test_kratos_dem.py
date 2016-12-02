""" This files provide a test example on how to use
Simphony along with the Kratos Multiphisics DEM wrapper

"""

import os

from simphony.core.cuba import CUBA

from simphony.api import CUDS, Simulation
from simphony.cuds.meta import api

# TODO: Utils now belong to probably another package
from simphony.engine import EngineInterface
from simphony.engine import kratos

# Add the problem path to the script
pathParticles = os.path.join(os.path.dirname(__file__), "3balls")
pathSolid = os.path.join(os.path.dirname(__file__), "3ballsDEM_FEM_boundary")

cuds = CUDS(name='example_kratos_dem_somulation')

# Integration time:
itime = api.IntegrationTime(name="dem_integration_time")
itime.time = 0.0001
itime.step = 0.001
itime.final = 600 * itime.step
cuds.add(itime)

# Utils are used to read an existing Kratos model as raw data so we can
# initialize the correct simphony datasets
utils = kratos.DEM_Utils()

# Reads Kratos mpda as a simphony data.
kratos_model_particles = utils.read_modelpart_as_particles(pathParticles)
kratos_model_solid = utils.read_modelpart_as_mesh(pathSolid)

# Add all models
for kratos_model in [kratos_model_particles, kratos_model_solid]:
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
sim = Simulation(cuds, "KRATOS_DEM", engine_interface=True)
sim.run()

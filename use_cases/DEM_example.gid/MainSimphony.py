""" This files provide a test example on how to use
Simphony along with the Kratos Multiphisics CFD wrapper

"""

import os

from simphony.api import CUDS, Simulation
from simphony.cuds.meta import api

from simkratos.DEM.kratos_DEM_utils import DEM_Utils


def abs_path(relPath):
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), relPath)


# Path for the Kratos' MDPA
pathParticles = abs_path("3ballsDEM")
pathSolid = abs_path("3ballsDEM_FEM_boundary")

cuds = CUDS(name='example_kratos_dem_somulation')

# Integration time:
itime = api.IntegrationTime(name="dem_integration_time")
itime.time = 0.0001
itime.step = 0.001
itime.final = 60 * itime.step
cuds.add([itime])

# Utils are used to read an existing Kratos model as raw data so we can
# initialize the correct simphony datasets. We can also use manualy written
# datasets.
utils = DEM_Utils()

# Reads Kratos mpda as a simphony data.
model_particles = utils.read_modelpart_as_particles(pathParticles)
model_solid = utils.read_modelpart_as_mesh(pathSolid)

# Add all datasets from the particles to the DEM wrapper
for model in [model_particles, model_solid]:
    cuds.add(list(model['datasets']))
    cuds.add(list(model['conditions']))
    cuds.add(list(model['materials']))
    cuds.add([model['pe']])

# Create the simulation and run the problem
sim = Simulation(cuds, "KRATOS_DEM", engine_interface=True)
sim.run()

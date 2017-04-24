""" This files provide a test example on how to use
Simphony along with the Kratos Multiphisics CFD wrapper

"""

import os

from simphony.core.cuba import CUBA
from simphony.cuds.mesh import Mesh

from simphony.api import CUDS, Simulation
from simphony.cuds.meta import api

# TODO: Utils now belong to probably another package
from simphony.engine import EngineInterface
from simphony.engine import kratos

from simkratos.CFD.kratos_CFD_utils import CFD_Utils

# Add the problem path to the script
path = os.path.join(os.path.dirname(__file__), "CFD_exampleFluid")

cuds = CUDS(name='example_fluid_name')

# Integration time:
itime = api.IntegrationTime(name="md_nve_integration_time")
itime.time = 0.0001
itime.step = 0.0025
itime.final = 0.0125  # 5 Kratos Timesteps
cuds.add([itime])

# Utils are used to read an existing Kratos model as raw data so we can
# initialize the correct simphony datasets
cfd_utils = CFD_Utils()

# Reads kratos data so its interpretable by simphony
model_fluid = cfd_utils.read_modelpart(path)

# Add the datasets readed from the conversor.
# NOTE: That the 'mesh' variable is called that way because CFD is going to
# use meshes, but it represents a dataset object.
# Add all datasets from the fluid to the CFD wrapper
for model in [model_fluid]:
    cuds.add(list(model['datasets']))
    cuds.add(list(model['conditions']))
    cuds.add(list(model['materials']))
    cuds.add([model['pe']])

# Create the simulation and run the problem
sim = Simulation(cuds, "KRATOS_CFD", engine_interface=True)
sim.run()

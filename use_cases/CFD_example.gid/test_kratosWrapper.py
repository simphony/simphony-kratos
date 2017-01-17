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

# Add the problem path to the script
path = os.path.join(os.path.dirname(__file__), "CFD_exampleFluid")

cuds = CUDS(name='example_fluid_name')

# Integration time:
itime = api.IntegrationTime(name="md_nve_integration_time")
itime.time = 0.0001
itime.step = 0.0025
itime.final = 0.0125  # 5 Kratos Timesteps
cuds.add(itime)

# Utils are used to read an existing Kratos model as raw data so we can
# initialize the correct simphony datasets
utils = kratos.CFD_Utils()

# Reads kratos data so its interpretable by simphony
kratos_model = utils.read_modelpart(path)

# Add the datasets readed from the conversor.
# NOTE: That the 'mesh' variable is called that way because CFD is going to
# use meshes, but it represents a dataset object.
for mesh in kratos_model['meshes']:
    cuds.add(mesh)

# Add some boundary conditions readed from the conversor.
for bc in kratos_model['bcs']:
    conditionName = 'condition_'+bc['name']

    # Generate a regular bondary condition
    condition = api.Condition(name=conditionName)

    # Apply the data
    conditionData = condition.data
    conditionData[CUBA.VELOCITY] = bc['velocity']
    conditionData[CUBA.PRESSURE] = bc['pressure']
    condition.data = conditionData

    # Add it to the CUDS
    cuds.add(condition)

# Create the simulation and run the problem
sim = Simulation(cuds, "KRATOS_CFD", engine_interface=True)
sim.run()

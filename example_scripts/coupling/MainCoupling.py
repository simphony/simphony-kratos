""" This files provide a test example on how to use
Simphony for coupling CFD with DEM

"""
import os
import math

from simphony.api import CUDS, Simulation
from simphony.core.cuba import CUBA
from simphony.cuds.meta import api

# This line is optional and only adds some tools to read kratos models
# as simphony datasets rather than crating them in the script.
from simkratos.CFD.kratos_CFD_utils import CFD_Utils
from simkratos.DEM.kratos_DEM_utils import DEM_Utils


def simple_coupling(cuds, particles_dataset_name):
    viscosity = 1.0e-3
    something_shaddy = 6.0
    particle_dataset = cuds.get_by_name(name=particles_dataset_name)

    for particle in particle_dataset.iter(item_type=CUBA.PARTICLES):
        p_data = particle.data

        vel = p_data[CUBA.VELOCITY]
        rel_vel = p_data[CUBA.RELATIVE_VELOCITY]
        radius = p_data[CUBA.RADIUS]

        factor = something_shaddy * math.pi * viscosity * radius

        external_force = [factor * (x - y) for x, y in zip(rel_vel, vel)]

        particle.data[CUBA.EXTERNAL_APPLIED_FORCE] = external_force
        particle_dataset.update([particle])


def abs_path(relPath):
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), relPath)


# Path for the Kratos' MDPA
pathFluid = abs_path("fluid_prism")
pathParticles = abs_path("fibers_and_ballsDEM")

cuds = CUDS(name='Coupling')

# Integration time for the fluid solver:
CFDitime = api.IntegrationTime(name="CFD_Integration_time")
CFDitime.time = 0.0001
CFDitime.step = 0.005
CFDitime.final = 0.0125  # 1 Kratos Timesteps

# Integration time for the particle solver
DEMitime = api.IntegrationTime(name="DEM_Integration_time")
DEMitime.time = 0.0001
DEMitime.step = 0.005
DEMitime.final = 0.00125  # 5 Kratos Timesteps

COUiTime = api.IntegrationTime(name="COU_Inetegration_time")

cuds.add([COUiTime])

# Utils are used to read an existing Kratos model as raw data so we can
# initialize the correct simphony datasets. We can also use manualy written
# datasets.
cfd_utils = CFD_Utils()
dem_utils = DEM_Utils()

# Reads kratos data as simphony datasets
model_fluid = cfd_utils.read_modelpart(pathFluid)
model_particles = dem_utils.read_modelpart_as_particles(pathParticles)

# Add all datasets from the fluid to the CFD wrapper
for model in [model_fluid]:
    cuds.add(list(model['datasets']))
    cuds.add(list(model['conditions']))
    cuds.add(list(model['materials']))
    cuds.add([model['pe']])

# Add all datasets from the particles to the DEM wrapper
for model in [model_particles]:
    cuds.add(list(model['datasets']))
    cuds.add(list(model['conditions']))
    cuds.add(list(model['materials']))
    cuds.add([model['pe']])

# Get a reference to the datasets that we will use in the coupling
fluid_dataset_name = model_fluid['datasets'][0].name
particles_dataset_name = model_particles['datasets'][0].name

# Create the simulation and run the problem
simCFD = Simulation(cuds, "KRATOS_CFD", engine_interface=True)
simDEM = Simulation(cuds, "KRATOS_DEM", engine_interface=True)
simPRO = Simulation(cuds, "KRATOS_PRO", engine_interface=True)
simGID = Simulation(cuds, "KRATOS_GID", engine_interface=True)

for step in xrange(0, 100):

    # Set interval times
    CFDitime.final = 0.0125 * (step + 1)  # 5 Kratos Timesteps
    DEMitime.final = 0.0125 * (step + 1)  # 5 Kratos Timesteps

    # Make sure the timestep between wrappers are consistent
    COUiTime.time = 0.0125 * (step + 0)
    COUiTime.step = CFDitime.step
    COUiTime.final = CFDitime.final

    # Run the CFD
    simCFD.run()

    # Projects the velocity from the fluid to the particles
    simPRO.run()

    # Coupling
    simple_coupling(cuds, particles_dataset_name)

    # Make sure the timestep between wrappers are consistent
    COUiTime.time = 0.0125 * (step + 0)
    COUiTime.step = DEMitime.step
    COUiTime.final = DEMitime.final

    # Run the DEM
    simDEM.run()

    # Print using the GiD wrapper
    simGID.run()

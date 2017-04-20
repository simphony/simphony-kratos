""" This files provide a test example on how to use
Simphony for coupling CFD with DEM

"""
import os

from simphony.api import CUDS, Simulation
from simphony.core.cuba import CUBA
from simphony.cuds.meta import api

# This line is optional and only adds some tools to read kratos models
# as simphony datasets rather than crating them in the script.
from simkratos.CFD.kratos_CFD_utils import CFD_Utils
from simkratos.DEM.kratos_DEM_utils import DEM_Utils


def abs_path(relPath):
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), relPath)


# Path for the Kratos' MDPA
pathFluid = abs_path("fluid_prism")
pathParticles = abs_path("fibers_and_ballsDEM")

cuds = CUDS(name='Coupling')

# Integration time for the fluid solver:
CFDitime = api.IntegrationTime(name="CFD_Integration_time")
CFDitime.time = 0.0001
CFDitime.step = 0.01
CFDitime.final = 2  # 1 Kratos Timesteps

# Integration time for the particle solver
DEMitime = api.IntegrationTime(name="DEM_Integration_time")
DEMitime.time = 0.0001
DEMitime.step = 0.0025
DEMitime.final = 0.0125  # 5 Kratos Timesteps

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
fluid_dataset = model_fluid['datasets'][0].name
particles_dataset = model_particles['datasets'][0].name

# Create the simulation and run the problem
simCFD = Simulation(cuds, "KRATOS_CFD", engine_interface=True)
simDEM = Simulation(cuds, "KRATOS_DEM", engine_interface=True)
simPRO = Simulation(cuds, "KRATOS_PRO", engine_interface=True)

for step in xrange(0, 5):

    print("\t====== ITERATION BEG ======")
    # Set interval times
    CFDitime.final = 0.0125 * (step + 1)  # 5 Kratos Timesteps
    DEMitime.final = 0.0125 * (step + 1)  # 5 Kratos Timesteps

    print("\t==== CFD ITERATION BEG ====")
    # Make sure the timestep between wrappers are consistent
    COUiTime.time = 0.0125 * (step + 0)
    COUiTime.step = CFDitime.step
    COUiTime.final = CFDitime.final

    # simCFD.run()
    print("\t==== CFD ITERATION END ====")

    print("\t==== COUPLING ITER BEG ====")
    # Projects the velocity from the fluid to the particles
    simPRO.run()

    # print("Getting ", fluid_dataset)
    # fluid_dataset = cuds.get_by_name(fluid_dataset)
    # for node in fluid_dataset.iter(item_type=CUBA.NODE):
    #     print(node.data)
    #
    # # Simple coupling (Drag force)
    # print("Getting ", particles_dataset)
    # particles_dataset = cuds.get_by_name(particles_dataset)
    # for particle in particles_dataset.iter(item_type=CUBA.PARTICLES):
    #     print(particle.data)

    # for node in spheres_model_part.Nodes:
    #     velocity_factor = 0.02
    #     vx = -1.0 * velocity_factor * node.Y
    #     vy = velocity_factor * node.X
    #     vz = velocity_factor * 0.0
    #
    #     radius = node.GetSolutionStepValue(RADIUS)
    #     viscosity = 1.0e-3
    #     factor = 6.0 * math.pi * viscosity * radius  # falta densitat
    #     fx = factor * vx
    #     fy = factor * vy
    #     fz = factor * vz
    #
    #     node.SetSolutionStepValue(EXTERNAL_APPLIED_FORCE_X, fx)
    #     node.SetSolutionStepValue(EXTERNAL_APPLIED_FORCE_Y, fy)
    #     node.SetSolutionStepValue(EXTERNAL_APPLIED_FORCE_Z, fz)

    print("\t==== COUPLING ITER END ====")

    print("\t==== DEM ITERATION BEG ====")
    # Make sure the timestep between wrappers are consistent
    COUiTime.time = 0.0125 * (step + 0)
    COUiTime.step = DEMitime.step
    COUiTime.final = DEMitime.final

    simDEM.run()
    print("\t==== DEM ITERATION END ====")

    print("All steps of the iteration executed correctly.")

    print("\t====== ITERATION END ======")

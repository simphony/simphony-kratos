""" This files provide a test example on how to use
Simphony for coupling CFD with DEM

"""
import os
import math

from simphony.api import CUDS, Simulation
from simphony.core.cuba import CUBA
from simphony.cuds.meta import api

from simphony.core.data_container import DataContainer
from simphony.cuds.particles import Particles as SParticles
from simphony.cuds.particles import Particle as SParticle

# This line is optional and only adds some tools to read kratos models
# as simphony datasets rather than crating them in the script.
from simkratos.CFD.kratos_CFD_utils import CFD_Utils


def simple_coupling(cuds, particles_dataset_name):
    viscosity = 1.0e-3
    stokes_constant = 6.0
    particle_dataset = cuds.get_by_name(name=particles_dataset_name)

    for particle in particle_dataset.iter(item_type=CUBA.PARTICLES):
        p_data = particle.data

        vel = p_data[CUBA.VELOCITY]
        rel_vel = p_data[CUBA.RELATIVE_VELOCITY]
        radius = p_data[CUBA.RADIUS]

        factor = stokes_constant * math.pi * viscosity * radius

        external_force = [factor * (x - y) for x, y in zip(rel_vel, vel)]

        particle.data[CUBA.EXTERNAL_APPLIED_FORCE] = external_force
        particle_dataset.update([particle])


def generate_particles(smp_particles, smp_conditions, smp_materials, smp_pe):
    # Define the particle containers.
    particles = SParticles(name='particles')

    # Fill the data ( particle )
    data = DataContainer()
    data[CUBA.RADIUS] = 5e-5

    particles.add([
        SParticle(
            coordinates=(0.002, 0.004, 0.004),
            data=DataContainer(data),
        )
    ])

    material = api.Material(name='material_' + particles.name)
    materialData = material.data

    materialData[CUBA.DENSITY] = 950.0
    materialData[CUBA.YOUNG_MODULUS] = 35e9
    materialData[CUBA.POISSON_RATIO] = 0.20
    materialData[CUBA.FRICTION_COEFFICIENT] = 0.5773502691896257
    # materialData[CUBA.PARTICLE_COHESION] = 0.0
    materialData[CUBA.RESTITUTION_COEFFICIENT] = 0.02
    materialData[CUBA.ROLLING_FRICTION] = 0.01
    # materialData[CUBA.DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME] = \
    # "DEM_KDEMFabric"
    # materialData[CUBA.DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME] = \
    # "DEM_D_Hertz_viscous_Coulomb"
    # materialData[CUBA.CONTACT_TAU_ZERO] = 25
    # materialData[CUBA.CONTACT_SIGMA_MIN] = 5
    # materialData[CUBA.CONTACT_INTERNAL_FRICC] = 1

    material.data = materialData

    particlesData = particles.data
    particlesData[CUBA.MATERIAL] = material.name
    particles.data = particlesData

    # Pack the return objects
    smp_particles.append(particles)
    smp_materials.append(material)

    # Add the datasets that will be used by the wrapper
    smp_pe.data[CUBA.DATA_SET].append(particles.name)

    return {
        'datasets': smp_particles,
        'conditions': smp_conditions,
        'materials': smp_materials,
        'pe': smp_pe,
    }


def generate_fibers(smp_particles, smp_conditions, smp_materials, smp_pe):
    # Define the particle containers.
    fibers = SParticles(name='fibers')

    # Fill the data ( fiber )
    data = DataContainer()
    data[CUBA.RADIUS] = 1e-5

    fiberCoords = [
        (0.1166666666666666685, 0.5083333333333333037, 0.5166666666666666075),
        (0.1499999999999999944, 0.5250000000000000222, 0.5500000000000000444),
        (0.1833333333333333481, 0.5416666666666667407, 0.5833333333333332593),
        (0.2166666666666666741, 0.5583333333333333481, 0.6166666666666666963),
        (0.2500000000000000000, 0.5749999999999999556, 0.6499999999999999112),
        (0.2833333333333333259, 0.5916666666666665630, 0.6833333333333333481)
    ]

    fibers.add([
        SParticle(
            coordinates=fiberCoords[0],
            data=DataContainer(data),
        ),
        SParticle(
            coordinates=fiberCoords[1],
            data=DataContainer(data),
        ),
        SParticle(
            coordinates=fiberCoords[2],
            data=DataContainer(data),
        ),
        SParticle(
            coordinates=fiberCoords[3],
            data=DataContainer(data),
        ),
        SParticle(
            coordinates=fiberCoords[4],
            data=DataContainer(data),
        ),
        SParticle(
            coordinates=fiberCoords[5],
            data=DataContainer(data),
        )
    ])

    material = api.Material(name='material_' + fibers.name)
    materialData = material.data

    materialData[CUBA.DENSITY] = 1050.0
    materialData[CUBA.YOUNG_MODULUS] = 1.0e9
    materialData[CUBA.POISSON_RATIO] = 0.20
    materialData[CUBA.FRICTION_COEFFICIENT] = 0.9999999999999999
    # materialData[CUBA.PARTICLE_COHESION] = 0.0
    materialData[CUBA.RESTITUTION_COEFFICIENT] = 0.02
    materialData[CUBA.ROLLING_FRICTION] = 0.01
    # materialData[CUBA.FABRIC_COEFFICIENT] = 0.1
    # materialData[CUBA.DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME] = \
    # "DEM_KDEMFabric"
    # materialData[CUBA.DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME] = \
    # "DEM_D_Hertz_viscous_Coulomb"
    # materialData[CUBA.CONTACT_TAU_ZERO] = 25
    # materialData[CUBA.CONTACT_SIGMA_MIN] = 5
    # materialData[CUBA.CONTACT_INTERNAL_FRICC] = 1

    material.data = materialData

    fibersData = fibers.data
    fibersData[CUBA.MATERIAL] = material.name
    fibers.data = fibersData

    # Pack the return objects
    smp_particles.append(fibers)
    smp_materials.append(material)

    # Add the datasets that will be used by the wrapper
    smp_pe.data[CUBA.DATA_SET].append(fibers.name)

    return {
        'datasets': smp_particles,
        'conditions': smp_conditions,
        'materials': smp_materials,
        'pe': smp_pe,
    }


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
DEMitime.step = 5e-5
DEMitime.final = 0.0125  # 5 Kratos Timesteps

COUiTime = api.IntegrationTime(name="COU_Inetegration_time")

cuds.add([COUiTime])

# Utils are used to read an existing Kratos model as raw data so we can
# initialize the correct simphony datasets. We can also use manualy written
# datasets.
cfd_utils = CFD_Utils()

# Reads kratos data as simphony datasets
model_fluid = cfd_utils.read_modelpart(pathFluid)

# Dem is generated on the script
smp_particles = []
smp_conditions = []
smp_materials = []

dem_pe = api.GranularDynamics()
dem_pe.data[CUBA.DATA_SET] = []

model_particles = generate_particles(
    smp_particles, smp_conditions, smp_materials, dem_pe
)
model_particles = generate_fibers(
    smp_particles, smp_conditions, smp_materials, dem_pe
)

# Add all datasets from the fluid to the CFD wrapper
for model in [model_fluid]:
    cuds.add(list(model['datasets']))
    cuds.add(list(model['conditions']))
    cuds.add(list(model['materials']))
    cuds.add([model['pe']])

# Add all datasets from the particles / fibers to the DEM wrapper
for model in [model_particles]:
    cuds.add(list(model['datasets']))
    cuds.add(list(model['conditions']))
    cuds.add(list(model['materials']))
    cuds.add([model['pe']])

# Get a reference to the datasets that we will use in the coupling
fluid_dataset_name = model_fluid['datasets'][0].name
particles_dataset_name = model_particles['datasets'][0].name
fibers_dataset_name = model_particles['datasets'][1].name

# Create the simulation and run the problem
simCFD = Simulation(cuds, "KRATOS_CFD", engine_interface=True)
simDEM = Simulation(cuds, "KRATOS_DEM", engine_interface=True)
simPRO = Simulation(cuds, "KRATOS_PRO", engine_interface=True)
simGID = Simulation(cuds, "KRATOS_GID", engine_interface=True)

for step in xrange(0, 10):

    # Set interval times
    CFDitime.final = 0.0125 * (step + 1)  # 5 Kratos Timesteps
    DEMitime.final = 0.0125 * (step + 1)  # 5 Kratos Timesteps

    # Make sure the timestep between wrappers are consistent
    COUiTime.time = 0.0125 * (step + 0)
    COUiTime.step = CFDitime.step
    COUiTime.final = CFDitime.final

    # Run the CFD
    simCFD.run()

for step in xrange(0, 100):
    # Projects the velocity from the fluid to the particles
    simPRO.run()

    # Coupling ( For both the particle and the fiber)
    simple_coupling(cuds, particles_dataset_name)
    simple_coupling(cuds, fibers_dataset_name)

    # Make sure the timestep between wrappers are consistent
    COUiTime.time = 0.0125 * (step + 0)
    COUiTime.step = DEMitime.step
    COUiTime.final = 0.0125 * (step + 0) + COUiTime.step * 10

    # Run the DEM
    simDEM.run()

    # Print using the GiD wrapper
    simGID.run()

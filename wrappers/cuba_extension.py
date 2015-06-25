""" Provisional CUBA keywords

"""

from enum import Enum, unique


@unique
class CUBAExt(Enum):
    """ Provisional CUBA keywords specific for Kratos

    These are additional CUBA-Keywords for Kratos-CFD and
    Kratos-DEMPack that are not included in simphony-common yet.
    The proposed description for these new CUBA keywords is:

    - description: Mass of the node
    domain: [MD]
    key: NODAL_MASS
    name: NodalMass
    number: 100
    shape: [1]
    type: double

    - description: Displacement vector
    domain: [MD]
    key: DISPLACEMENT
    name: Displacement
    number: 101
    shape: [3]
    type: double

    - description: Vector with the element forces
    domain: [MD]
    key: TOTAL_FORCES
    name: TotalForces
    number: 102
    shape: [3]
    type: double

    - description: Friction coeficient with walls
    domain: [MD]
    key: WALL_FRICTION
    name: WallFriction
    number: 103
    shape: [1]
    type: double

    - description: Density of the particles
    domain: [MD]
    key: PARTICLE_DENSITY
    name: ParticleDensity
    number: 104
    shape: [1]
    type: double

    - description: Friction of the particles
    domain: [MD]
    key: PARTICLE_FRICTION
    name: ParticleFriction
    number: 105
    shape: [1]
    type: double

    - description: Cohesion of the particles
    domain: [MD]
    key: PARTICLE_COHESION
    name: ParticleCohesion
    number: 106
    shape: [1]
    type: double

    - description: Logarithm of the restitution coeficient
    domain: [MD]
    key: LN_OF_RESTITUTION_COEFF
    name: LdOfRestitutionCoeff
    number: 107
    shape: [1]
    type: double

    - description: Material of the particle
    domain: [MD]
    key: PARTICLE_MATERIAL
    name: ParticleMaterial
    number: 108
    shape: [1]
    type: double

    - description: Continuum law name
    domain: [MD]
    key: DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME
    name: ContinuumLawName
    number: 109
    shape: [1]
    type: string

    - description: Discontinuum law name
    domain: [MD]
    key: DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME
    name: DiscontinuumLawName
    number: 110
    shape: [1]
    type: string

    - description: Name of the fluid mesh
    domain: [DEM,CFD]
    key: FLUID_MESH_NAME
    name: FluidMeshName
    number: 109
    shape: [1]
    type: string

    - description: Name of the structure mesh
    domain: [DEM]
    key: STRUCTURE_MESH_NAME
    name: StructureMeshName
    number: 110
    shape: [1]
    type: string

"""

    NODAL_MASS = 0
    DISPLACEMENT = 1
    TOTAL_FORCES = 2
    WALL_FRICTION = 3
    PARTICLE_DENSITY = 4
    PARTICLE_FRICTION = 5
    PARTICLE_COHESION = 6
    LN_OF_RESTITUTION_COEFF = 7
    PARTICLE_MATERIAL = 8
    DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME = 9
    DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME = 10

    FLUID_MESH_NAME = 11
    STRUCTURE_MESH_NAME = 12 

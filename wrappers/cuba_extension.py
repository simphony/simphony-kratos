""" Provisional CUBA keywords

"""

from enum import Enum, unique


@unique
class CUBAExtension(Enum):
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

    - description: Reaction TOTAL_FORCES
    domain: [MD]
    key: REACTION
    name: Reaction
    number: 111
    shape: [1]
    type: double

    - description: Distance
    domain: [MD]
    key: REACTION
    name: Reaction
    number: 112
    shape: [1]
    type: double

    - description: Viscosity
    domain: [MD]
    key: VISCOSITY
    name: Viscosity
    number: 113
    shape: [1]
    type: double

    - description: Force of the body
    domain: [MD]
    key: BODY_FORCE
    name: BodyForce
    number: 114
    shape: [1]
    type: double

    - description: Flags of the node
    domain: [MD]
    key: FLAG_VARIABLE
    name: FlagVariable
    number: 115
    shape: [1]
    type: integer

    - description: Boolean for structure check
    domain: [MD]
    key: IS_STRUCTURE
    name: FlagVariable
    number: 116
    shape: [1]
    type: boolean

    - description: Boolean for slip check
    domain: [MD]
    key: IS_SLIP
    name: FlagVariable
    number: 117
    shape: [1]
    type: boolean

    - description: Boolean for imposed pressure check
    domain: [MD]
    key: IMPOSED_PRES
    name: FlagVariable
    number: 118
    shape: [1]
    type: boolean

    - description: Boolean for imposed velocity check
    domain: [MD]
    key: IMPOSED_VEL
    name: FlagVariable
    number: 119
    shape: [1]
    type: boolean

"""

    NODAL_MASS = 100
    DISPLACEMENT = 101
    TOTAL_FORCES = 102
    WALL_FRICTION = 103
    PARTICLE_DENSITY = 104
    PARTICLE_FRICTION = 105
    PARTICLE_COHESION = 106
    LN_OF_RESTITUTION_COEFF = 107
    PARTICLE_MATERIAL = 108
    DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME = 109
    DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME = 110

    REACTION = 111
    DISTANCE = 112
    VISCOSITY = 113
    BODY_FORCE = 114
    FLAG_VARIABLE = 115
    IS_STRUCTURE = 116
    IS_SLIP = 117

    IMPOSED_PRES = 118
    IMPOSED_VEL = 119

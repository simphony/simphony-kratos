""" TODO: Description was out of date

"""

# Kratos works with Python3 by default.
from __future__ import print_function, absolute_import, division

# Simphony Imports
from simphony.core.cuba import CUBA
from simphony.core.data_container import DataContainer

from simphony.cuds.mesh import Point as SPoint
from simphony.cuds.mesh import Face as SFace
from simphony.cuds.mesh import Cell as SCell
from simphony.cuds.particles import Particle as SParticle

# Wrapper Imports
from simkratos.kratosWrapper import KratosWrapper
from simkratos.DEM import DEM_explicit_solver_var as DEM_parameters

# Kratos Imports
import KratosMultiphysics as KRTS
import KratosMultiphysics.DEMApplication as KRTSDEM

import sphere_strategy as SolverStrategy
import DEM_procedures


class DEMWrapper(KratosWrapper):

    def __init__(self, use_internal_interface=True, **kwargs):
        super(DEMWrapper, self).__init__(use_internal_interface, **kwargs)

        self.time = 0
        self.step = 0
        self.substeps = 0
        self.particle_type = "SphericParticle3D"
        self.condition_type = "RigidFace3D3N"

        # The dictionary defines the relation between CUBA and
        # kratos variables

        self.variables_dictionary = {
            "RADIUS": [
                CUBA.RADIUS,
                KRTS.RADIUS
            ],
            "NODAL_MASS": [
                None,
                KRTS.NODAL_MASS
            ],
            "VELOCITY": [
                CUBA.VELOCITY,
                KRTS.VELOCITY,
                KRTS.VELOCITY_X,
                KRTS.VELOCITY_Y,
                KRTS.VELOCITY_Z
            ],
            "DISPLACEMENT": [
                None,
                KRTS.DISPLACEMENT,
                KRTS.DISPLACEMENT_X,
                KRTS.DISPLACEMENT_Y,
                KRTS.DISPLACEMENT_Z
            ]
        }

        self.initialize()

    def _load_cuds(self):
        """Load CUDS data into lammps engine."""
        cuds = self.get_cuds()
        if not cuds:
            return

        for component in cuds.iter(item_type=CUBA.MESH):
            self.add_dataset(component)

    def getNodalData(self, data, node, model):
        """ Extracts the node data

        Extracts the node data and puts in ina format readeable
        by the Simphony DataContainer

        """

        if model == "Particles":
            self.getSolutionStepVariable1D(data, node, "RADIUS")
            self.getSolutionStepVariable1D(data, node, "NODAL_MASS")
            self.getSolutionStepVariable3D(data, node, "VELOCITY")
            self.getSolutionStepVariable3D(data, node, "DISPLACEMENT")

    def setNodalData(self, data, node, model):
        """ Assembles the point data

        Assembles the point data and puts in ina format readeable
        by the Kratos ModelPart

        """

        if model == "Particles":
            self.setSolutionStepVariable1D(data, node, "RADIUS")
            self.setSolutionStepVariable1D(data, node, "NODAL_MASS")
            self.setSolutionStepVariable3D(data, node, "VELOCITY")
            self.setSolutionStepVariable3D(data, node, "DISPLACEMENT")

    def _setMeshData(self):
        " This probably needs to be done throug configuration"

        cLawString = "DEMContinuumConstitutiveLaw"
        dLawString = "DEMDiscontinuumConstitutiveLaw"

        self.SP[CUBA.DENSITY] = 2500.0
        self.SP[CUBA.YOUNG_MODULUS] = 1.0e5
        self.SP[CUBA.POISSON_RATIO] = 0.20
        self.SP[CUBA.ROLLING_FRICTION] = 0.01

        self.PARTICLE_FRICTION = 0.99
        self.PARTICLE_COHESION = 0.0
        self.LN_OF_RESTITUTION_COEFF = -1.6094379124341003
        self.PARTICLE_MATERIAL = 1
        self.WALL_FRICTION = 0.3
        self.DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME = cLawString
        self.DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME = dLawString

    def setElementData(self):
        self.element_properties.SetValue(
            KRTSDEM.PARTICLE_DENSITY,
            self.SP[CUBA.DENSITY]
        )
        self.element_properties.SetValue(
            KRTS.YOUNG_MODULUS,
            self.SP[CUBA.YOUNG_MODULUS]
        )
        self.element_properties.SetValue(
            KRTS.POISSON_RATIO,
            self.SP[CUBA.POISSON_RATIO]
        )
        self.element_properties.SetValue(
            KRTSDEM.PARTICLE_FRICTION,
            self.PARTICLE_FRICTION
        )
        self.element_properties.SetValue(
            KRTSDEM.PARTICLE_COHESION,
            self.PARTICLE_COHESION
        )
        self.element_properties.SetValue(
            KRTSDEM.LN_OF_RESTITUTION_COEFF,
            self.LN_OF_RESTITUTION_COEFF
        )
        self.element_properties.SetValue(
            KRTS.PARTICLE_MATERIAL,
            self.PARTICLE_MATERIAL
        )
        self.element_properties.SetValue(
            KRTSDEM.ROLLING_FRICTION,
            self.SP[CUBA.ROLLING_FRICTION]
        )
        self.element_properties.SetValue(
            KRTSDEM.DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME,
            self.DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME
        )
        self.element_properties.SetValue(
            KRTSDEM.DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME,
            self.DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME
        )

    def setConditionData(self):
        cmesh = self.get_dataset("Conditions")

        self.condition_properties.SetValue(
            KRTS.WALL_FRICTION,
            cmesh.data[CUBA.WALL_FRICTION]
        )

    def initialize(self):
        """ Initalizes the necessary kratos commponents

            Initalizes the necessary kratos commponents to
            execute Kratos' DEMPack solver
        """

        # DEMPack SubClasses
        self.procedures = DEM_procedures.Procedures(DEM_parameters)
        self.parallelutils = DEM_procedures.ParallelUtils()
        self.materialTest = DEM_procedures.MaterialTest()
        self.creator_destructor = KRTSDEM.ParticleCreatorDestructor()

        # Prepare ModelParts
        self.spheres_model_part = KRTS.ModelPart("Particles")
        self.rigid_face_model_part = KRTS.ModelPart("Conditions")
        self.mixed_model_part = KRTS.ModelPart("")
        self.cluster_model_part = KRTS.ModelPart("")
        self.DEM_inlet_model_part = KRTS.ModelPart("")
        self.mapping_model_part = KRTS.ModelPart("")
        self.contact_model_part = KRTS.ModelPart("")

        # Create solver
        self.solver = SolverStrategy.ExplicitStrategy(
            self.spheres_model_part,
            self.rigid_face_model_part,
            self.cluster_model_part,
            self.DEM_inlet_model_part,
            self.creator_destructor,
            DEM_parameters
        )

        # Prepare properties
        self.element_properties = KRTS.Properties(0)
        self.condition_properties = KRTS.Properties(1)

        self._setMeshData()
        self.setElementData()

        # Prepare variables
        self.solver.AddAdditionalVariables(
            self.spheres_model_part,
            DEM_parameters
        )
        self.procedures.AddCommonVariables(
            self.spheres_model_part,
            DEM_parameters
        )
        self.procedures.AddSpheresVariables(
            self.spheres_model_part,
            DEM_parameters
        )
        self.procedures.AddCommonVariables(
            self.rigid_face_model_part,
            DEM_parameters
        )
        self.procedures.AddRigidFaceVariables(
            self.rigid_face_model_part,
            DEM_parameters
        )

        # Set a search strategy
        self.solver.search_strategy = self.parallelutils.GetSearchStrategy(
            self.solver,
            self.spheres_model_part
        )

        self.solver.Initialize()

    def run(self):
        """ Run a step of the wrapper """

        fluid_particles = self.pcs
        solid_meshes = self.meshes

        cuds = self.get_cuds()

        self.spheres_model_part.GetMesh(len(fluid_particles))
        self.rigid_face_model_part.GetMesh(len(solid_meshes))

        fluid_properties = KRTS.PropertiesArray()
        meshNumber = 1
        meshDict = {}

        for particles in cuds.iter(item_type=CUBA.PARTICLE):

            group = meshNumber

            self.importKratosParticles(
                particles, self.spheres_model_part,
                group, self.particle_type
            )

            meshDict[particles.name] = meshNumber
            meshNumber += 1

        self.updateBackwardDicc()

        fluid_properties.append(self.element_properties)
        # solid_properties.append(self.condition_properties)

        self.spheres_model_part.SetProperties(fluid_properties)
        # self.rigid_face_model_part.SetProperties(solid_properties)

        SolverStrategy.AddDofs(self.spheres_model_part)

        self.solver.Initialize()

        self.dt = cuds.get_by_name('dem_integration_time').step

        # Start the simulation itself
        self.time = cuds.get_by_name('dem_integration_time').time
        self.final = cuds.get_by_name('dem_integration_time').final

        # Solve
        while self.time < self.final:

            self.dt = self.spheres_model_part.ProcessInfo.GetValue(
                KRTS.DELTA_TIME
            )

            cuds.get_by_name('dem_integration_time').step = self.dt

            self.spheres_model_part.ProcessInfo[KRTS.TIME] = self.time
            self.spheres_model_part.ProcessInfo[KRTS.DELTA_TIME] = self.dt
            self.spheres_model_part.ProcessInfo[KRTS.TIME_STEPS] = self.step

            self.rigid_face_model_part.ProcessInfo[KRTS.TIME] = self.time
            self.rigid_face_model_part.ProcessInfo[KRTS.DELTA_TIME] = self.dt
            self.rigid_face_model_part.ProcessInfo[KRTS.TIME_STEPS] = self.step

            self.solver.Solve()

            self.step += 1
            self.time = self.time + self.dt

        cuds.get_by_name('dem_integration_time').time = self.time
        cuds.get_by_name('dem_integration_time').final = self.final

        for particles in cuds.iter(item_type=CUBA.PARTICLE):

            group = meshDict[particles.name]

            self.exportKratosParticles(
                self.spheres_model_part,
                particles,
                group
            )

        self.updateForwardDicc()

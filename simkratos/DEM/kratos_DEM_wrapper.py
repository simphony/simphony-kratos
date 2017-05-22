""" TODO: Description was out of date

"""

# Kratos works with Python3 by default.
from __future__ import print_function, absolute_import, division

# Simphony Imports
from simphony.core.cuba import CUBA

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
            "DENSITY": [
                CUBA.DENSITY,
                KRTSDEM.PARTICLE_DENSITY
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
                CUBA.DELTA_DISPLACEMENT,
                KRTS.DISPLACEMENT,
                KRTS.DISPLACEMENT_X,
                KRTS.DISPLACEMENT_Y,
                KRTS.DISPLACEMENT_Z
            ],
            "EXTERNAL_APPLIED_FORCE": [
                CUBA.EXTERNAL_APPLIED_FORCE,
                KRTSDEM.EXTERNAL_APPLIED_FORCE,
                KRTSDEM.EXTERNAL_APPLIED_FORCE_X,
                KRTSDEM.EXTERNAL_APPLIED_FORCE_Y,
                KRTSDEM.EXTERNAL_APPLIED_FORCE_Z
            ]
        }

        self.properties_dictionary = {
            "PARTICLE_DENSITY": [
                CUBA.DENSITY,
                KRTSDEM.PARTICLE_DENSITY
            ],
            "YOUNG_MODULUS": [
                CUBA.YOUNG_MODULUS,
                KRTS.YOUNG_MODULUS
            ],
            "POISSON_RATIO": [
                CUBA.POISSON_RATIO,
                KRTS.POISSON_RATIO
            ],
            "PARTICLE_FRICTION": [
                CUBA.FRICTION_COEFFICIENT,
                KRTSDEM.PARTICLE_FRICTION
            ],
            "PARTICLE_COHESION": [
                None,
                KRTSDEM.PARTICLE_COHESION
            ],
            "ROLLING_FRICTION": [
                CUBA.ROLLING_FRICTION,
                KRTSDEM.PARTICLE_COHESION
            ],
            "DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME": [
                None,
                KRTSDEM.PARTICLE_COHESION
            ],
            "DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME": [
                None,
                KRTSDEM.PARTICLE_COHESION
            ]
        }

        # Set DEM parameters
        self.solver_strategy = SolverStrategy
        self.creator_destructor = KRTSDEM.ParticleCreatorDestructor()
        self.dem_fem_search = KRTSDEM.DEM_FEM_Search()
        self.procedures = DEM_procedures.Procedures(DEM_parameters)

        # Define some paths as they are nedded by some modules of them.
        # These paths will NOT be used in this wrapper.
        self.graphs_path = '.'

        # This should not be nedded for Simphony
        self.parallelutils = DEM_procedures.ParallelUtils()
        self.scheme = self.procedures.SetScheme()

        self.initialize()

    def _load_cuds(self):
        """Load CUDS data into kratos engine."""
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

        if model == "SpheresPart":
            self.getSolutionStepVariable1D(data, node, "RADIUS")
            self.getSolutionStepVariable1D(data, node, "NODAL_MASS")
            self.getSolutionStepVariable3D(data, node, "VELOCITY")
            self.getSolutionStepVariable3D(data, node, "DISPLACEMENT")
            self.getSolutionStepVariable3D(
                data, node, "EXTERNAL_APPLIED_FORCE"
            )

    def setNodalData(self, data, node, model):
        """ Assembles the point data

        Assembles the point data and puts in ina format readeable
        by the Kratos ModelPart

        """

        if model == "SpheresPart":
            self.setSolutionStepVariable1D(data, node, "RADIUS")
            self.setSolutionStepVariable1D(data, node, "NODAL_MASS")
            self.setSolutionStepVariable3D(data, node, "VELOCITY")
            self.setSolutionStepVariable3D(data, node, "DISPLACEMENT")
            self.setSolutionStepVariable3D(
                data, node, "EXTERNAL_APPLIED_FORCE"
            )

    def setDefaultElementData(self, propertyContainer):
        propertyContainer.SetValue(KRTSDEM.PARTICLE_DENSITY, 2500.0)
        propertyContainer.SetValue(KRTS.YOUNG_MODULUS, 1.0e5)
        propertyContainer.SetValue(KRTS.POISSON_RATIO, 0.20)
        propertyContainer.SetValue(KRTSDEM.PARTICLE_FRICTION, 0.99)
        propertyContainer.SetValue(KRTSDEM.PARTICLE_COHESION, 0.0)
        propertyContainer.SetValue(KRTS.PARTICLE_MATERIAL, 1)
        propertyContainer.SetValue(KRTSDEM.ROLLING_FRICTION, 0.99)

        propertyContainer.SetValue(
            KRTSDEM.DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME,
            "DEM_KDEM"
        )
        propertyContainer.SetValue(
            KRTSDEM.DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME,
            "DEM_D_Hertz_viscous_Coulomb"
        )

    def setProperties(self, material, propertyContainer):
        for propertyName in self.properties_dictionary:
            self.setProperty(material.data, propertyContainer, propertyName)

    def setConditionData(self):
        cmesh = self.get_dataset("RigidFacePart")

        self.condition_properties.SetValue(
            KRTS.WALL_FRICTION,
            cmesh.data[CUBA.WALL_FRICTION]
        )

    def initialize(self):
        """ Initalizes the necessary kratos commponents

            Initalizes the necessary kratos commponents to
            execute Kratos' DEMPack solver
        """

        # Prepare ModelParts
        self.spheres_model_part = KRTS.ModelPart("SpheresPart")
        self.rigid_face_model_part = KRTS.ModelPart("RigidFacePart")
        self.cluster_model_part = KRTS.ModelPart("ClusterPart")
        self.DEM_inlet_model_part = KRTS.ModelPart("DEMInletPart")
        self.mapping_model_part = KRTS.ModelPart("MappingPart")
        self.contact_model_part = KRTS.ModelPart("ContactPart")

        self.all_model_parts = DEM_procedures.SetOfModelParts([
            self.spheres_model_part,
            self.rigid_face_model_part,
            self.cluster_model_part,
            self.DEM_inlet_model_part,
            self.mapping_model_part,
            self.contact_model_part
        ])

        # Prepare the solver
        self.solver = SolverStrategy.ExplicitStrategy(
            self.all_model_parts,
            self.creator_destructor,
            self.dem_fem_search,
            self.scheme,
            DEM_parameters,
            self.procedures
        )

        # Honestly I don't know what this does
        self.procedures.AddAllVariablesInAllModelParts(
            self.solver,
            self.scheme,
            self.all_model_parts,
            DEM_parameters
        )

        # Read the modelparts.
        # This part is skipped and its done once the run is called.
        # Instead we:

        # ###################### #
        # ##     Simphony     ## #
        # ###################### #

        # Prepare properties
        self.kratos_props = {0: KRTS.Properties(0)}
        properties_array = KRTS.PropertiesArray()

        properties_array.append(self.kratos_props[0])

        self.spheres_model_part.SetProperties(properties_array)
        self.setDefaultElementData(properties_array[0])

        # Prepare variables
        self.solver.AddAdditionalVariables(
            self.spheres_model_part,
            DEM_parameters
        )

        self.spheres_model_part.AddNodalSolutionStepVariable(
            KRTSDEM.EXTERNAL_APPLIED_FORCE
        )

        self.procedures.solver = self.solver

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

        # ###################### #
        # ##     Simphony     ## #
        # ###################### #

        # Set Buffers for historical data
        self.procedures.SetUpBufferSizeInAllModelParts(
            self.spheres_model_part, 1,
            self.cluster_model_part, 1,
            self.DEM_inlet_model_part, 1,
            self.rigid_face_model_part, 1
        )

        # Adding dofs
        self.solver.AddDofs(self.spheres_model_part)
        self.solver.AddDofs(self.cluster_model_part)
        self.solver.AddDofs(self.DEM_inlet_model_part)

        # Set a search strategy
        self.solver.search_strategy = self.parallelutils.GetSearchStrategy(
            self.solver,
            self.spheres_model_part
        )

        # Prepare things before initialization
        self.solver.BeforeInitialize()
        self.parallelutils.Repart(self.spheres_model_part)
        self.bounding_box_time_limits = self.procedures.SetBoundingBoxLimits(
            self.all_model_parts,
            self.creator_destructor
        )

        # Finding the max id of the nodes
        max_Id = self.procedures.FindMaxNodeIdAccrossModelParts(
            self.creator_destructor,
            self.all_model_parts
        )

        self.creator_destructor.SetMaxNodeId(max_Id)

        # Initialize the solver
        self.solver.Initialize()

    def run(self):
        """ Run a step of the wrapper """

        fluid_particles = self.pcs
        solid_meshes = self.meshes

        cuds = self.get_cuds()

        self.spheres_model_part.GetMesh(len(fluid_particles))
        self.rigid_face_model_part.GetMesh(len(solid_meshes))

        meshNbr = 1
        meshDict = {}

        # Get the CFD pe
        if cuds.count_of(item_type=CUBA.GRANULAR_DYNAMICS) != 1:
            raise Exception("KratosDEM only allows one GRANULAR_DYNAMICS pe.")

        for gd_pe in cuds.iter(item_type=CUBA.GRANULAR_DYNAMICS):
            if len(gd_pe.data[CUBA.DATA_SET]) < 1:
                raise Exception("GD PE does not have any associated dataset")

            for name in gd_pe.data[CUBA.DATA_SET]:

                particles = cuds.get_by_name(name)
                group = meshNbr

                self.importKratosParticles(
                    particles, self.spheres_model_part,
                    group, self.particle_type
                )

                if(CUBA.MATERIAL in particles.data):
                    material = cuds.get_by_name(particles.data[CUBA.MATERIAL])
                    model_props = self.spheres_model_part.GetProperties()

                    if meshNbr not in self.kratos_props.keys():
                        self.kratos_props[meshNbr] = KRTS.Properties(meshNbr)
                        model_props[meshNbr] = self.kratos_props[meshNbr]
                        self.setDefaultElementData(model_props[meshNbr])

                    self.setProperties(material, model_props[meshNbr])

                meshDict[particles.name] = meshNbr
                meshNbr += 1

        self.updateBackwardDicc()

        self.solver.AddDofs(self.spheres_model_part)

        self.solver.Initialize()

        print("DEM Solver correctly initialized ...")

        if cuds.count_of(item_type=CUBA.INTEGRATION_TIME) < 0:
            raise Exception("No integran time")

        if cuds.count_of(item_type=CUBA.INTEGRATION_TIME) > 1:
            raise Exception("More than one integration time")

        iTime = [it for it in cuds.iter(item_type=CUBA.INTEGRATION_TIME)][0]

        self.dt = iTime.step

        # Start the simulation itself
        self.time = iTime.time
        self.final = iTime.final

        # Solve
        while self.time < self.final:

            self.dt = self.spheres_model_part.ProcessInfo.GetValue(
                KRTS.DELTA_TIME
            )

            iTime.step = self.dt

            self.spheres_model_part.ProcessInfo[KRTS.TIME] = self.time
            self.spheres_model_part.ProcessInfo[KRTS.DELTA_TIME] = self.dt
            self.spheres_model_part.ProcessInfo[KRTS.TIME_STEPS] = self.step

            self.rigid_face_model_part.ProcessInfo[KRTS.TIME] = self.time
            self.rigid_face_model_part.ProcessInfo[KRTS.DELTA_TIME] = self.dt
            self.rigid_face_model_part.ProcessInfo[KRTS.TIME_STEPS] = self.step

            self.solver.Solve()

            self.step += 1
            self.time = self.time + self.dt

        iTime.time = self.time
        iTime.final = self.final

        for gd_pe in cuds.iter(item_type=CUBA.GRANULAR_DYNAMICS):
            if len(gd_pe.data[CUBA.DATA_SET]) < 1:
                raise Exception("GD PE does not have any associated dataset")

            for name in gd_pe.data[CUBA.DATA_SET]:

                particles = cuds.get_by_name(name)
                group = meshDict[particles.name]

                self.exportKratosParticles(
                    self.spheres_model_part,
                    particles,
                    group
                )

        self.updateForwardDicc()

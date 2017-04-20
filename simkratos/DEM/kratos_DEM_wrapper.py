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

    def _setMeshData(self):
        " This probably needs to be done throug configuration"

        cLawString = "DEM_KDEM"
        dLawString = "DEM_D_Hertz_viscous_Coulomb"

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
        # self.element_properties.SetValue(
        #     KRTSDEM.LN_OF_RESTITUTION_COEFF,
        #     self.LN_OF_RESTITUTION_COEFF
        # )
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
        self.element_properties = KRTS.Properties(0)
        self.condition_properties = KRTS.Properties(1)

        element_properties_array = KRTS.PropertiesArray()

        self._setMeshData()
        self.setElementData()

        element_properties_array.append(self.element_properties)
        # solid_properties.append(self.condition_properties)

        self.spheres_model_part.SetProperties(element_properties_array)
        # self.rigid_face_model_part.SetProperties(solid_properties)

        # Prepare variables
        self.solver.AddAdditionalVariables(
            self.spheres_model_part,
            DEM_parameters
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

        # Constructing the inlet and initializing it
        # if (DEM_parameters.dem_inlet_option):
        #     self.DEM_inlet = DEM_Inlet(self.DEM_inlet_model_part)
        #     self.DEM_inlet.InitializeDEM_Inlet(
        #         self.spheres_model_part,
        #         self.creator_destructor,
        #         self.solver.continuum_type
        #     )

        # Enable this
        # self.DEMFEMProcedures = DEM_procedures.DEMFEMProcedures(
        #     DEM_parameters,
        #     self.graphs_path,
        #     self.spheres_model_part,
        #     self.rigid_face_model_part
        # )

    def run(self):
        """ Run a step of the wrapper """

        fluid_particles = self.pcs
        solid_meshes = self.meshes

        cuds = self.get_cuds()

        self.spheres_model_part.GetMesh(len(fluid_particles))
        self.rigid_face_model_part.GetMesh(len(solid_meshes))

        meshNumber = 1
        meshDict = {}

        # Get the CFD pe
        if cuds.count_of(item_type=CUBA.GRANULAR_DYNAMICS) != 1:
            raise Exception("KratosDEM only allows one GRANULAR_DYNAMICS pe.")

        for gd_pe in cuds.iter(item_type=CUBA.GRANULAR_DYNAMICS):
            if len(gd_pe.data[CUBA.DATA_SET]) < 1:
                raise Exception("GD PE does not have any associated dataset")

            for name in gd_pe.data[CUBA.DATA_SET]:

                particles = cuds.get_by_name(name)
                group = meshNumber

                self.importKratosParticles(
                    particles, self.spheres_model_part,
                    group, self.particle_type
                )

                meshDict[particles.name] = meshNumber
                meshNumber += 1

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

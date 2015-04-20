""" TODO: Description was out of date

"""

# Kratos works with Python3 by default.
from __future__ import print_function, absolute_import, division

# Simphony Imports
from simphony.core.cuba import CUBA
from simphony.core.data_container import DataContainer

from simphony.cuds.mesh import Point as SPoint
from simphony.cuds.mesh import Mesh as SMesh
from simphony.cuds.mesh import Face as SFace
from simphony.cuds.mesh import Cell as SCell

# Wrapper Imports
from wrappers.kratosWrapper import KratosWrapper

# Kratos Imports
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

import DEM_explicit_solver_var as DEM_parameters
import sphere_strategy as SolverStrategy
import DEM_procedures


class DEMPackWrapper(KratosWrapper):

    def __init__(self):
        KratosWrapper.__init__(self)
        self.time = 0
        self.step = 0

        # The dictionary defines the relation between CUBA and
        # kratos variables

        self.variables_dictionary = {
            "RADIUS": [
                CUBA.RADIUS,
                RADIUS
            ],
            "NODAL_MASS": [
                CUBA.NODAL_MASS,
                NODAL_MASS
            ],
            "VELOCITY": [
                CUBA.VELOCITY,
                VELOCITY,
                VELOCITY_X,
                VELOCITY_Y,
                VELOCITY_Z
            ],
            "DISPLACEMENT": [
                CUBA.DISPLACEMENT,
                DISPLACEMENT,
                DISPLACEMENT_X,
                DISPLACEMENT_Y,
                DISPLACEMENT_Z
            ],
            "TOTAL_FORCES": [
                CUBA.TOTAL_FORCES,
                TOTAL_FORCES,
                TOTAL_FORCES_X,
                TOTAL_FORCES_Y,
                TOTAL_FORCES_Z
            ]
        }

    def __addNodalVariablesToModelpart(self):
        """ Adds the DEMPack nodal variables

        Adds the DEMPack nodal variables to the particle and
        solid Kratos modelparts in order to be usable later
        while importing the mesh.

        """

        self.procedures.AddCommonVariables(
            self.spheres_model_part,
            DEM_parameters
        )
        self.procedures.AddSpheresVariables(
            self.spheres_model_part,
            DEM_parameters
        )

        SolverStrategy.AddAdditionalVariables(
            self.spheres_model_part,
            DEM_parameters
        )
        SolverStrategy.AddAdditionalVariables(
            self.rigid_face_model_part,
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

    def __getNodalData(self, data, node):
        """ Extracts the node data

        Extracts the node data and puts in ina format readeable
        by the Simphony DataContainer

        """

        self.__getSolutionStepVariable1D(data, node, "RADIUS")
        self.__getSolutionStepVariable1D(data, node, "NODAL_MASS")
        self.__getSolutionStepVariable3D(data, node, "VELOCITY")
        self.__getSolutionStepVariable3D(data, node, "DISPLACEMENT")
        self.__getSolutionStepVariable3D(data, node, "TOTAL_FORCES")

    def __setNodalData(self, data, node):
        """ Assembles the point data

        Assembles the point data and puts in ina format readeable
        by the Kratos ModelPart

        """

        self.__setSolutionStepVariable1D(data, node, "RADIUS")
        self.__setSolutionStepVariable1D(data, node, "NODAL_MASS")
        self.__setSolutionStepVariable3D(data, node, "VELOCITY")
        self.__setSolutionStepVariable3D(data, node, "DISPLACEMENT")
        self.__setSolutionStepVariable3D(data, node, "TOTAL_FORCES")

    def __exportKratosElements(self, src, dst, entitylist=None):
        """ Parses all kratos elements to simphony cells

        Iterates over all nodes in the kratos mesh (src) and
        converts them to simphony cells (dst). While doing this operation
        any point that has not currently mapped will have his uuid
        added in the 'id_map' of the weapper

        """

        for element in src.GetElements():

            point_list = [
                self.id_to_uuid_node_map[point.Id]
                for point in element.GetNodes()
            ]

            for node in element.GetNodes():

                data = {}

                self.__getNodalData(data, node)

                point = SPoint(
                    coordinates=(node.X, node.Y, node.Z),
                    data=DataContainer(data),
                    uid=self.id_to_uuid_node_map[node.Id]
                )

                dst.add_point(point)

            cell = SCell(
                points=point_list,
                data=DataContainer(data),
                uid=self.id_to_uuid_element_map[element.Id]
            )

            dst.add_cell(cell)

    def __exportKratosConditions(self, src, dst, entitylist=None):
        """ Parses all kratos conditions to simphony faces

        Iterates over all nodes in the kratos mesh ( src ) and
        converts them to simphony faces (dst). While doing this operation
        any point that has not currently mapped will have his uuid
        added in the 'id_map' of the weapper

        """

        for condition in src.GetConditions():

            point_list = [
                self.id_to_uuid_node_map[point.Id]
                for point in condition.GetNodes()
            ]

            for node in src.GetNodes():

                data = {}

                point = SPoint(
                    coordinates=(node.X, node.Y, node.Z),
                    data=DataContainer(data),
                    uid=self.id_to_uuid_node_map[node.Id]
                )

                try:
                    dst.add_point(point)
                except:
                    pass

            face = SFace(
                points=point_list,
                data=DataContainer(data),
                uid=self.id_to_uuid_condition_map[condition.Id]
            )

            dst.add_face(face)

    def __importKratosElements(self, src, dst, entitylist=None):
        """ Parses all simphony cells to kratos elements

        Iterates over all cells in the simphony mesh (src) and
        converts them to kratos SphericContinuumParticle3D elements (dst).
        While doing this operation any point/node pair that has not
        currently mapped will have  his uuid added in the 'id_map'
        of the wrapper

        """

        for element in src.iter_cells(entitylist):

            if element.uid not in self.uuid_to_id_element_map.keys():
                self.uuid_to_id_element_map.update(
                    {element.uid: self.free_id}
                )

                self.free_id += 1

            element_id = self.uuid_to_id_element_map[element.uid]

            for point in src.iter_points(element.points):

                if point.uid not in self.uuid_to_id_node_map.keys():
                    self.uuid_to_id_node_map.update(
                        {point.uid: self.free_id}
                    )

                    self.free_id += 1

                node_id = self.uuid_to_id_node_map[point.uid]

                data = point.data

                node = dst.CreateNewNode(
                    node_id,
                    point.coordinates[0],
                    point.coordinates[1],
                    point.coordinates[2])

                self.__setNodalData(data, node)

            dst.CreateNewElement(
                "SphericParticle3D",
                element_id,
                [self.uuid_to_id_node_map[p] for p in element.points],
                self.element_properties)

    def __importKratosConditions(self, src, dst, entitylist=None):
        """ Parses all simphony faces to kratos conditions

        Iterates over all faces in the simphony mesh (src) and
        converts them to kratos RigidFace3D3N conditions (dst).
        While doing this operation any point/node pair that has not
        currently mapped will have  his uuid added in the 'id_map'
        of the wrapper

        """

        for condition in src.iter_faces(entitylist):

            if condition.uid not in self.uuid_to_id_condition_map.keys():
                self.uuid_to_id_condition_map.update(
                    {condition.uid: self.free_id}
                )

                self.free_id += 1

            condition_id = self.uuid_to_id_condition_map[condition.uid]

            for point in src.iter_points(condition.points):

                if point.uid not in self.uuid_to_id_node_map.keys():
                    self.uuid_to_id_node_map.update(
                        {point.uid: self.free_id}
                    )

                    self.free_id += 1

                node_id = self.uuid_to_id_node_map[point.uid]

                if node_id not in [node.Id for node in dst.Nodes]:
                    node = dst.CreateNewNode(
                        node_id,
                        point.coordinates[0],
                        point.coordinates[1],
                        point.coordinates[2])

            dst.CreateNewCondition(
                "RigidFace3D3N",
                condition_id,
                [self.uuid_to_id_node_map[p] for p in condition.points],
                self.condition_properties)

    def read_modelpart(self, filename):
        """ Reads a Kratos formated modelpart NYI

        This adds partial support for the future FileIO
        """

        # f = open(filename, 'r')
        pass

    def write_modelpart(self, filename):
        """ Writes a Kratos formated modelpart

        This adds partial support for the future FileIO
        """

        pmesh = self.get_mesh("Particles")

        f = open(filename, 'w')

        # Variable information
        f.write('Begin ModelPartData\n')
        f.write('//  VARIABLE_NAME value\n')
        f.write('End ModelPartData\n\n')

        # Properties ( out shared values )
        f.write('Begin Properties 1\n')
        f.write('PARTICLE_DENSITY {}\n'.format(pmesh.data[PARTICLE_DENSITY]))
        f.write('YOUNG_MODULUS {}\n'.format(pmesh.data[YOUNG_MODULUS]))
        f.write('POISSON_RATIO {}\n'.format(pmesh.data[POISSON_RATIO]))
        f.write('PARTICLE_FRICTION {}\n'.format(pmesh.data[PARTICLE_FRICTION]))
        f.write('PARTICLE_COHESION {}\n'.format(pmesh.data[PARTICLE_COHESION]))
        f.write('LN_OF_RESTITUTION_COEFF {}\n'.format(
            pmesh.data[LN_OF_RESTITUTION_COEFF]))
        f.write('PARTICLE_MATERIAL {}\n'.format(pmesh.data[PARTICLE_MATERIAL]))
        f.write('ROLLING_FRICTION {}\n'.format(pmesh.data[ROLLING_FRICTION]))
        f.write('DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME {}\n'.format(
            pmesh.data[DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME]))
        f.write('DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME {}\n'.format(
            pmesh.data[DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME]))
        f.write('End Properties\n\n')

        # Nodes
        f.write('Begin Nodes\n')
        for point in src.iter_points(entitylist=None):
            f.write('{} {} {} {}\n').format(
                self.uuid_to_id_point_map[point.id],
                point.coordinates[0],
                point.coordinates[1],
                point.coordinates[2]
            )
        f.write('End Nodes\n\n')

        f.write('Begin Elements SphericParticle3D')
        for element in src.iter_cells(entitylist=None):
            f.write('{} {} {}\n').format(
                self.uuid_to_id_point_map[elemet.id],
                self.uuid_to_id_point_map[element.points[0]],
                self.uuid_to_id_point_map[element.points[1]]
            )
        f.write('End Elements\n\n')

        f.write('Begin NodalData RADIUS')
        for point in src.iter_points(entitylist=None):
            f.write('{} {} {} {}\n').format(
                self.uuid_to_id_point_map[point.id],
                0,
                point.data[CUBA.RADIUS],
            )
        f.write('End NodalData\n\n')

    def setMeshData(self, mesh):
        " This probably needs to be done throug configuration"

        cLawString = "DEMContinuumConstitutiveLaw"
        dLawString = "DEMDiscontinuumConstitutiveLaw"

        data = mesh.data

        data[CUBA.PARTICLE_DENSITY] = 2500.0
        data[CUBA.YOUNG_MODULUS] = 1.0e5
        data[CUBA.POISSON_RATIO] = 0.20
        data[CUBA.PARTICLE_FRICTION] = 0.99
        data[CUBA.PARTICLE_COHESION] = 0.0
        data[CUBA.LN_OF_RESTITUTION_COEFF] = -1.6094379124341003
        data[CUBA.PARTICLE_MATERIAL] = 1
        data[CUBA.ROLLING_FRICTION] = 0.01
        data[CUBA.WALL_FRICTION] = 0.3
        data[CUBA.DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME] = cLawString
        data[CUBA.DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME] = dLawString

        mesh.data = data

    def setElementData(self):
        pmesh = self.get_mesh("Particles")

        self.element_properties.SetValue(
            PARTICLE_DENSITY,
            pmesh.data[CUBA.PARTICLE_DENSITY]
        )
        self.element_properties.SetValue(
            YOUNG_MODULUS,
            pmesh.data[CUBA.YOUNG_MODULUS]
        )
        self.element_properties.SetValue(
            POISSON_RATIO,
            pmesh.data[CUBA.POISSON_RATIO]
        )
        self.element_properties.SetValue(
            PARTICLE_FRICTION,
            pmesh.data[CUBA.PARTICLE_FRICTION]
        )
        self.element_properties.SetValue(
            PARTICLE_COHESION,
            pmesh.data[CUBA.PARTICLE_COHESION]
        )
        self.element_properties.SetValue(
            LN_OF_RESTITUTION_COEFF,
            pmesh.data[CUBA.LN_OF_RESTITUTION_COEFF]
        )
        self.element_properties.SetValue(
            PARTICLE_MATERIAL,
            pmesh.data[CUBA.PARTICLE_MATERIAL]
        )
        self.element_properties.SetValue(
            ROLLING_FRICTION,
            pmesh.data[CUBA.ROLLING_FRICTION]
        )
        self.element_properties.SetValue(
            DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME,
            pmesh.data[CUBA.DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME]
        )
        self.element_properties.SetValue(
            DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME,
            pmesh.data[CUBA.DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME]
        )

    def setConditionData(self):
        cmesh = self.get_mesh("Conditions")

        self.condition_properties.SetValue(
            WALL_FRICTION,
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
        self.creator_destructor = ParticleCreatorDestructor()
        self.demio = DEM_procedures.DEMIo()

        # Prepare ModelParts
        self.spheres_model_part = ModelPart("")
        self.rigid_face_model_part = ModelPart("")
        self.mixed_model_part = ModelPart("")
        self.cluster_model_part = ModelPart("")
        self.DEM_inlet_model_part = ModelPart("")
        self.mapping_model_part = ModelPart("")
        self.contact_model_part = ModelPart("")

        self.element_properties = Properties(0)
        self.condition_properties = Properties(1)

        self.elemNodeData = {CUBA.RADIUS: RADIUS}
        self.condNodeData = {}

        # Initialize GiD-IO
        self.demio.AddGlobalVariables()
        self.demio.AddSpheresVariables()
        self.demio.AddFEMBoundaryVariables()
        self.demio.AddClusterVariables()
        self.demio.AddContactVariables()
        #
        self.demio.AddMpiVariables()
        self.demio.EnableMpiVariables()

        self.demio.Configure(
            DEM_parameters.problem_name,
            DEM_parameters.OutputFileType,
            DEM_parameters.Multifile,
            DEM_parameters.ContactMeshOption
        )

        self.demio.SetOutputName(
            DEM_parameters.problem_name
        )

        self.demio.InitializeMesh(
            self.mixed_model_part,
            self.spheres_model_part,
            self.rigid_face_model_part,
            self.cluster_model_part,
            self.contact_model_part,
            self.mapping_model_part
        )

    def initializeTimeStep(self):

        SolverStrategy.AddDofs(self.spheres_model_part)

        # Create solver
        self.solver = SolverStrategy.ExplicitStrategy(
            self.spheres_model_part,
            self.rigid_face_model_part,
            self.cluster_model_part,
            self.DEM_inlet_model_part,
            self.creator_destructor,
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

        newSphereMp = ModelPart("")
        newRigidbMp = ModelPart("")

        self.spheres_model_part = newSphereMp
        self.rigid_face_model_part = newRigidbMp

        newSphereMp.Properties.append(self.element_properties)
        newRigidbMp.Properties.append(self.condition_properties)

        self.__addNodalVariablesToModelpart()

        # Import the into Kratos
        self.__importKratosElements(
            self.get_mesh("Particles"),
            newSphereMp
        )

        self.__importKratosConditions(
            self.get_mesh("Particles"),
            newRigidbMp
        )

        self.updateBackwardDicc()
        self.setElementData()
        # self.setConditionData()

        self.initializeTimeStep()

        print("NON: {}".format(self.solver.model_part.NumberOfNodes(0)))

        # Not sure what to do here :S
        newSphereMp.ProcessInfo[TIME] = self.time
        newSphereMp.ProcessInfo[DELTA_TIME] = 0.5  # NYI
        newSphereMp.ProcessInfo[TIME_STEPS] = self.step

        newRigidbMp.ProcessInfo[TIME] = self.time
        newRigidbMp.ProcessInfo[DELTA_TIME] = 0.5  # NYI
        newRigidbMp.ProcessInfo[TIME_STEPS] = self.step

        substeps = 3

        # Solve
        for n in xrange(0, substeps):
            self.solver.Solve()

        new_mesh = SMesh(name="Particles")

        # Add the problem data
        self.setMeshData(new_mesh)

        # Export data back to SimPhoNy
        self.__exportKratosElements(
            self.spheres_model_part,
            new_mesh
        )

        self.__exportKratosConditions(
            self.rigid_face_model_part,
            new_mesh
        )

        if (self.step % 5) == 0:
            self.demio.PrintResults(
                self.mixed_model_part,
                self.spheres_model_part,
                self.rigid_face_model_part,
                self.cluster_model_part,
                self.contact_model_part,
                self.mapping_model_part, self.time
            )

        for p in self.spheres_model_part.Nodes:
            print(p.GetSolutionStepValue(DISPLACEMENT))

        self.add_mesh(new_mesh)
        self.updateForwardDicc()

        self.time += self.spheres_model_part.ProcessInfo.GetValue(DELTA_TIME)
        self.step += 1

    def finalizeTimeStep(self):
        pass

    def finalize(self):
        self.demio.FinalizeMesh()

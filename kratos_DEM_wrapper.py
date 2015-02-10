""" Template to convert modelparts from kratos to simphony 

This file show an example of how to utilze the kratos wrappers
in order to import or export models from KratosMultiphysics

"""

# Kratos works with Python3 by default.
from __future__ import print_function, absolute_import, division

# Simphony Imports
from simphony.core.cuba import CUBA
from simphony.core.data_container import DataContainer

from simphony.cuds.mesh import Mesh as SimphonyMesh
from simphony.cuds.mesh import Point as SimphonyPoint
from simphony.cuds.mesh import Edge as SimphonyEdge
from simphony.cuds.mesh import Face as SimphonyFace
from simphony.cuds.mesh import Cell as SimphonyCell

from uuid import *

from kratosWrapper import KratosWrapper

# Kratos Imports 
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

import DEM_explicit_solver_var as DEM_parameters
import sphere_strategy as SolverStrategy
import DEM_procedures

class KratosWrapper(KratosWrapper):

    def __init__(self):
        super(KratosWrapper,self).__init__()

    def __addNodalVariablesToModelpart(self,model_part):
        """ Adds the IncompressiveFluidApplication nodal variables

        Adds the IncompressiveFluidApplication nodal variables to the
        kratos modelpart provided in order to be usable later while
        importing the mesh.

        """

        self.procedures.AddCommonVariables(model_part, DEM_parameters)
        self.procedures.AddSpheresVariables(model_part, DEM_parameters)
        self.procedures.AddMpiVariables(model_part)
        # SolverStrategy.AddAdditionalVariables(spheres_model_part, DEM_parameters)
        # self.procedures.AddCommonVariables(DEM_inlet_model_part, DEM_parameters)
        # self.procedures.AddSpheresVariables(DEM_inlet_model_part, DEM_parameters)
        # SolverStrategy.AddAdditionalVariables(DEM_inlet_model_part, DEM_parameters)  
        # self.procedures.AddCommonVariables(rigid_face_model_part, DEM_parameters)
        # self.procedures.AddRigidFaceVariables(rigid_face_model_part, DEM_parameters)
        # self.procedures.AddMpiVariables(rigid_face_model_part)

    def __exportKratosNodes(self,src,dst):
        """ Parses all kratos nodes to simphony points

        Iterates over all nodes in the kratos mesh (src) and
        converts them to simphony points (dst). While doing this operation
        any point that has not currently mapped will have his uuid
        added in the 'id_map' of the weapper

        """

        for node in src.GetNodes():
            
            data = {
                CUBA.RADIUS:    node.GetSolutionStepValue(RADIUS),
                CUBA.VELOCITY:  [
                    node.GetSolutionStepValue(VELOCITY_X),
                    node.GetSolutionStepValue(VELOCITY_Y),
                    node.GetSolutionStepValue(VELOCITY_Z)
                ]
            }

            point = SimphonyPoint(
                coordinates=(node.X, node.Y, node.Z),
                data=DataContainer(data)
            )

            point_uuid = dst.add_point(point)

            if point_uuid not in self.id_to_uuid_node_map.keys():
                self.id_to_uuid_node_map.update(
                    {node.Id:point_uuid}
                )

    def __exportKratosElements(self,src,dst):
        """ Parses all kratos elements to simphony cells

        Iterates over all nodes in the kratos mesh (src) and
        converts them to simphony cells (dst). While doing this operation
        any point that has not currently mapped will have his uuid
        added in the 'id_map' of the weapper

        """

        for element in src.GetElements():

            data = {
                CUBA.RADIUS:    element.GetValue(RADIUS),
            }

            point_list = []

            for point in element.GetNodes():
                point_list.append(self.id_to_uuid_node_map[point.Id])

            point = SimphonyCell(
                points=point_list,
                data=DataContainer(data)
            )

            element_uuid = dst.add_cell(point)

            if element_uuid not in self.id_to_uuid_element_map.keys():
                self.id_to_uuid_element_map.update(
                    {element.Id:element_uuid}
                )

    def __exportKratosConditions(self,src,dst):
        """ Parses all kratos conditions to simphony faces

        Iterates over all nodes in the kratos mesh ( src ) and
        converts them to simphony faces (dst). While doing this operation
        any point that has not currently mapped will have his uuid
        added in the 'id_map' of the weapper

        """

        for condition in src.GetConditions():
            data = {
                CUBA.RADIUS:    condition.GetValue(RADIUS),
            }

            point_list = []

            for point in condition.GetNodes():
                point_list.append(self.id_to_uuid_node_map[point.Id])

            point = SimphonyFace(
                points=point_list,
                data=DataContainer(data)
            )

            condition_uuid = dst.add_face(point)

            if condition_uuid not in self.id_to_uuid_condition_map.keys():
                self.id_to_uuid_condition_map.update(
                    {condition.Id:condition_uuid}
                )

    def __importKratosNodes(self,src,dst):
        """ Parses all simphony points to kratos nodes

        Iterates over all points in the simphony mesh ( src ) and
        converts them to kratos nodes. While doing this operation
        any point/node pair that has not currently mapped will have 
        his uuid added in the 'id_map' of the wrapper

        """
        
        for point in src.iter_points():

            if point.uuid not in self.uuid_to_id_node_map.keys():
                self.uuid_to_id_node_map.update(
                    {point.uuid:self.free_id}
                )

                self.free_id += 1

            node_id = self.uuid_to_id_node_map[point.uuid]

            dst.CreateNewNode(
                node_id,
                point.coordinates[0],
                point.coordinates[1],
                point.coordinates[2])

    def __importKratosElements(self,src,dst):
        """ Parses all simphony cells to kratos elements

        Iterates over all cells in the simphony mesh (src) and
        converts them to kratos SphericContinuumParticle3D elements (dst). 
        While doing this operation any point/node pair that has not 
        currently mapped will have  his uuid added in the 'id_map' 
        of the wrapper

        """

        properties = Properties(0)
        
        for element in src.iter_cells():

            if element.uuid not in self.uuid_to_id_element_map.keys():
                self.uuid_to_id_element_map.update(
                    {element.uuid:self.free_id}
                )

                self.free_id += 1

            element_id = self.uuid_to_id_element_map[element.uuid]

            dst.CreateNewElement(
                "SphericContinuumParticle3D",
                element_id,
                [self.uuid_to_id_node_map[element.points[0]]],
                properties)

    def __importKratosConditions(self,src,dst):
        """ Parses all simphony faces to kratos conditions

        Iterates over all faces in the simphony mesh (src) and
        converts them to kratos XXXXXX conditions (dst). 
        While doing this operation any point/node pair that has not 
        currently mapped will have  his uuid added in the 'id_map' 
        of the wrapper

        """

        properties = Properties(0)

        for condition in src.iter_faces():

            if condition.uuid not in self.uuid_to_id_condition_map.keys():
                self.uuid_to_id_condition_map.update(
                    {condition.uuid:self.free_id}
                )

                self.free_id += 1

            condition_id = self.uuid_to_id_condition_map[condition.uuid]

            dst.CreateNewCondition(
                "XXXXXXXX",
                condition_id,
                [self.uuid_to_id_node_map[condition.points[0]]],
                properties)

    def initialize(self):
        # DEMPack SubClasses
        self.procedures = DEM_procedures.Procedures(DEM_parameters)
        # self.demio = DEM_procedures.DEMIo()
        # self.report = DEM_procedures.Report()
        self.parallelutils = DEM_procedures.ParallelUtils()
        self.materialTest = DEM_procedures.MaterialTest()
        self.creator_destructor = ParticleCreatorDestructor()

        # Prepare ModelParts
        self.spheres_model_part =       ModelPart("") # self.importMesh(self.fluid)
        self.rigid_face_model_part =    ModelPart("") # self.importMesh(self.rigid)
        self.mixed_model_part =         ModelPart("")
        self.cluster_model_part =       ModelPart("") ## DUDA
        self.DEM_inlet_model_part =     ModelPart("") ## DUDA
        self.mapping_model_part =       ModelPart("") ## DUDA
        self.contact_model_part =       ModelPart("") ## DUDA

        # Not sure where to put this
        # SolverStrategy.AddDofs(self.spheres_model_part)

    def run(self):

        time = 0
        step = 0

        # Import the data to Kratos
        self.spheres_model_part = self.importMesh(self.meshes["Fluid"],"DEMParticles")
        # self.rigid_face_model_part = self.importMesh(self.meshes["Structure"],"Structure")

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

        # Not sure what to do here :S
        self.spheres_model_part.ProcessInfo[TIME] = time
        self.spheres_model_part.ProcessInfo[DELTA_TIME] = self.spheres_model_part.ProcessInfo.GetValue(DELTA_TIME)
        self.spheres_model_part.ProcessInfo[TIME_STEPS] = step
        
        self.rigid_face_model_part.ProcessInfo[TIME] = time
        self.rigid_face_model_part.ProcessInfo[DELTA_TIME] = self.spheres_model_part.ProcessInfo.GetValue(DELTA_TIME)
        self.rigid_face_model_part.ProcessInfo[TIME_STEPS] = step

        # Is this necessary?
        # mesh_motion.MoveAllMeshes(rigid_face_model_part, time, dt)

        # Solve
        self.solver.Initialize()
        self.solver.Solve()

        # Export data back to SimPhoNy
        self.meshes["Fluid"] = self.exportMesh(self.spheres_model_part)
        # self.meshes["Structure"] = self.exportMesh(self.rigid_face_model_part)
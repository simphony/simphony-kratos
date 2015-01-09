""" Template to convert modelparts from kratos to simphony 

This file show an example of how to utilze the kratos wrappers
in order to import or export models from KratosMultiphysics

"""

# Kratos will works with Python3 by default.
# This line makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Simphony Imports
from simphony.core.cuba import CUBA
from simphony.core.data_container import DataContainer

from simphony.cuds.mesh import Mesh as SimphonyMesh
from simphony.cuds.mesh import Point as SimphonyPoint
from simphony.cuds.mesh import Edge as SimphonyEdge
from simphony.cuds.mesh import Face as SimphonyFace
from simphony.cuds.mesh import Cell as SimphonyCell

from kratosWrapper import KratosWrapper

# Kratos Imports 
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

# Uuid and other dependences
from uuid import *

class KratosWrapper(KratosWrapper):

    def __init__(self):
        super(KratosWrapper,self).__init__()

    def __addNodalVariablesToModelpart(self,model_part):
        """ Adds the IncompressiveFluidApplication nodal variables

        Adds the IncompressiveFluidApplication nodal variables to the
        kratos modelpart provided in order to be usable later while
        importing the mesh.

        """

        model_part.AddNodalSolutionStepVariable(VELOCITY)
        model_part.AddNodalSolutionStepVariable(PRESSURE)
        model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(DENSITY)

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
                self.id_to_uuid_node_map.update({node.Id:point_uuid})

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
                self.id_to_uuid_element_map.update({element.Id:element_uuid})

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
                self.id_to_uuid_condition_map.update({condition.Id:condition_uuid})

    def __importKratosNodes(self,src,dst):
        """ Parses all simphony points to kratos nodes

        Iterates over all points in the simphony mesh ( src ) and
        converts them to kratos nodes. While doing this operation
        any point/node pair that has not currently mapped will have 
        his uuid added in the 'id_map' of the wrapper

        """
        
        for point in src.iter_points():

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

        for condition in src.iter_faces():

            condition_id = self.uuid_to_id_condition_map[condition.uuid]

            dst.CreateNewCondition(
                "XXXXXXXX",
                condition_id,
                [self.uuid_to_id_node_map[condition.points[0]]],
                properties)

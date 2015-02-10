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
from KratosMultiphysics.IncompressiveFluidApplication import *

# Uuid and other dependences
from uuid import *

class KratosIncompressiveFluidWrapper(KratosWrapper):

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
        model_part.AddNodalSolutionStepVariable(VISCOSITY)
        model_part.AddNodalSolutionStepVariable(DENSITY)
        model_part.AddNodalSolutionStepVariable(BODY_FORCE)
        model_part.AddNodalSolutionStepVariable(POROSITY)

    def __exportKratosNodes(self,src,dst):
        """ Parses all kratos nodes to simphony points

        Iterates over all nodes in the kratos mesh (src) and
        converts them to simphony points (dst). While doing this operation
        any point that has not currently mapped will have his uuid
        added in the 'id_map' of the weapper

        """

        for node in src.GetNodes():

            data = {
                CUBA.VELOCITY:  [
                    node.GetSolutionStepValue(VELOCITY_X),
                    node.GetSolutionStepValue(VELOCITY_Y),
                    node.GetSolutionStepValue(VELOCITY_Z)
                ]
                CUBA.PRESSURE:  [
                    node.GetSolutionStepValue(PRESSURE_X),
                    node.GetSolutionStepValue(PRESSURE_Y),
                    node.GetSolutionStepValue(PRESSURE_Z)
                ]
                CUBA.DISPLACEMENT:  [
                    node.GetSolutionStepValue(DISPLACEMENT_X),
                    node.GetSolutionStepValue(DISPLACEMENT_Y),
                    node.GetSolutionStepValue(DISPLACEMENT_Z)
                ]
                CUBA.VISCOSITY: node.GetSolutionStepValue(VISCOSITY),
                CUBA.DENSITY: node.GetSolutionStepValue(DENSITY),
                CUBA.EXTERNAL_APPLIED_FORCE: node.GetSolutionStepValue(BODY_FORCE),
            }

            point = SimphonyPoint(
                coordinates=(node.X, node.Y, node.Z),
                data=DataContainer(data)
            )

            point_uuid = dst.add_point(point)

            if point_uuid not in self.id_map.keys():
                self.id_map.update({node.Id:point_uuid})

    def __exportKratosElements(self,src,dst):
        """ Parses all kratos elements to simphony cells

        Iterates over all nodes in the kratos mesh (src) and
        converts them to simphony cells (dst). While doing this operation
        any point that has not currently mapped will have his uuid
        added in the 'id_map' of the weapper

        """

        for element in src.GetElements():

            data = {
                CUBA.UNDEFINED1:    element.GetValue(IS_STRUCTURE),
                CUBA.UNDEFINED2:    element.GetValue(IS_SLIP),
                CUBA.UNDEFINED2:    element.GetValue(Y_WALL)            # Can be left outside alone
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
                CUBA.UNDEFINED1:    element.GetValue(IS_STRUCTURE),
                CUBA.UNDEFINED2:    element.GetValue(IS_SLIP),
                CUBA.UNDEFINED2:    element.GetValue(Y_WALL)            # Can be left outside alone
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

            kratos_node = dst.CreateNewNode(
                node_id,
                point.coordinates[0],
                point.coordinates[1],
                point.coordinates[2])

            kratos_node.SetSolutionStepValue(VELOCITY_X) =
                point.data[CUBA.VELOCITY][0]
            kratos_node.SetSolutionStepValue(VELOCITY_Y) =
                point.data[CUBA.VELOCITY][1]
            kratos_node.SetSolutionStepValue(VELOCITY_Z) =
                point.data[CUBA.VELOCITY][2]

            kratos_node.SetSolutionStepValue(PRESSURE_X) =
                point.data[CUBA.PRESSURE][0]
            kratos_node.SetSolutionStepValue(PRESSURE_Y) =
                point.data[CUBA.PRESSURE][1]
            kratos_node.SetSolutionStepValue(PRESSURE_Z) =
                point.data[CUBA.PRESSURE][2]

            kratos_node.SetSolutionStepValue(DISPLACEMENT_X) =
                point.data[CUBA.DISPLACEMENT][0]
            kratos_node.SetSolutionStepValue(DISPLACEMENT_Y) =
                point.data[CUBA.DISPLACEMENT][1]
            kratos_node.SetSolutionStepValue(DISPLACEMENT_Z) =
                point.data[CUBA.DISPLACEMENT][2]

            kratos_node.SetSolutionStepValue(VISCOSITY) =
                point.data[CUBA.VISCOSITY]
            kratos_node.SetSolutionStepValue(DENSITY) =
                point.data[CUBA.DENSITY]
            kratos_node.SetSolutionStepValue(BODY_FORCE) =
                point.data[CUBA.EXTERNAL_APPLIED_FORCE]
            kratos_node.SetSolutionStepValue(POROSITY) =
                1

    def __importKratosElements(self,src,dst):
        """ Parses all simphony cells to kratos elements

        Iterates over all cells in the simphony mesh (src) and
        converts them to kratos SphericContinuumParticle3D elements (dst). 
        While doing this operation any point/node pair that has not 
        currently mapped will have  his uuid added in the 'id_map' 
        of the wrapper

        """

        for element in src.iter_cells():

            element_id = self.uuid_to_id_element_map[element.uuid]

            kratos_element = dst.CreateNewElement(
                element_id,
                uuid_to_id_node_map[element.points[0]],
                uuid_to_id_node_map[element.points[1]],
                uuid_to_id_node_map[element.points[2]])

            kratos_element.SetValue(IS_STRUCTURE) = UNDEFINED1
            kratos_element.SetValue(IS_SLIP) = UNDEFINED1
            kratos_element.SetValue(Y_WALL) = 0

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

            condition_id.SetValue(IS_STRUCTURE) = UNDEFINED1
            condition_id.SetValue(IS_SLIP) = UNDEFINED1
            condition_id.SetValue(Y_WALL) = 0
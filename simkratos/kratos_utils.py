from simphony.core.cuba import CUBA
from simphony.core.data_container import DataContainer

from simphony.cuds.mesh import Point as SPoint
from simphony.cuds.mesh import Mesh as SMesh
from simphony.cuds.mesh import Face as SFace
from simphony.cuds.mesh import Cell as SCell

from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MeshingApplication import *


class CFD_Utils(object):

    def __init__(self):
        self.id_to_uuid_node_map = {}
        self.uuid_to_id_node_map = {}
        self.id_to_uuid_element_map = {}
        self.uuid_to_id_element_map = {}
        self.id_to_uuid_condition_map = {}
        self.uuid_to_id_condition_map = {}

        self.variables_dictionary = {
            "PRESSURE": [
                CUBA.PRESSURE,
                PRESSURE
            ],
            "VELOCITY": [
                CUBA.VELOCITY,
                VELOCITY,
                VELOCITY_X,
                VELOCITY_Y,
                VELOCITY_Z
            ],
            "VISCOSITY": [
                None,
                VISCOSITY
            ],
            "DENSITY": [
                CUBA.DENSITY,
                DENSITY
            ]
        }

    def _getSolutionStepVariable1D(self, data, entity, variable):
        pair = self.variables_dictionary[variable]
        if(pair[0] is not None):
            data.update({
                pair[0]: entity.GetSolutionStepValue(pair[1])
            })

    def _getSolutionStepVariable3D(self, data, entity, variable):
        pair = self.variables_dictionary[variable]
        if(pair[0] is not None):
            data.update({
                pair[0]: [
                    entity.GetSolutionStepValue(pair[2]),
                    entity.GetSolutionStepValue(pair[3]),
                    entity.GetSolutionStepValue(pair[4])
                ]
            })

    def _getNodalData(self, data, node):
        """ Extracts the node data

        Extracts the node data and puts in ina format readable
        by the Simphony DataContainer

        """

        self._getSolutionStepVariable1D(data, node, "PRESSURE")
        self._getSolutionStepVariable3D(data, node, "VELOCITY")
        self._getSolutionStepVariable1D(data, node, "VISCOSITY")
        self._getSolutionStepVariable1D(data, node, "DENSITY")

    def _exportKratosNodes(self, src, dst, group):
        """ Parses all kratos nodes to simphony points

        Iterates over all nodes in the kratos mesh (src) and
        converts them to simphony points (dst). While doing this operation
        any node/point that has not currently been mapped will have his uuid
        added in the 'id_map' of the wrapper

        """

        for node in src.GetNodes(group):

            data = {}

            self._getNodalData(data, node)

            point_uid = None

            if node.Id in self.id_to_uuid_node_map:
                point_uid = self.id_to_uuid_node_map[node.Id]

            point = SPoint(
                coordinates=(node.X, node.Y, node.Z),
                data=DataContainer(data),
                uid=point_uid
            )

            pid = dst.add_points([point])

            self.id_to_uuid_node_map[node.Id] = pid[0]

    def _exportKratosElements(self, src, dst, group):
        """ Parses all kratos elements to simphony cells

        Iterates over all elements in the kratos mesh (src) and
        converts them to simphony cells (dst). While doing this operation
        any element/cell that has not currently been mapped will have his uuid
        added in the 'id_map' of the wrapper

        """

        for element in src.GetElements(group):

            element_uid = None

            if element.Id in self.id_to_uuid_element_map:
                element_uid = self.id_to_uuid_element_map[element.Id]

            point_list = [
                self.id_to_uuid_node_map[pointl.Id]
                for pointl in element.GetNodes()
            ]

            cell = SCell(
                points=point_list,
                uid=element_uid
            )

            cid = dst.add_cells([cell])

            self.id_to_uuid_element_map[element.Id] = cid[0]

    def _exportKratosConditions(self, src, dst, group):
        """ Parses all kratos conditions to simphony faces

        Iterates over all conditions in the kratos mesh (src) and
        converts them to simphony faces (dst). While doing this operation
        any condition/face that has not currently been mapped will have
        his uuid added in the 'id_map' of the wrapper

        """

        for condition in src.GetConditions(group):

            condition_uid = None

            if condition.Id in self.id_to_uuid_condition_map:
                condition_uid = self.id_to_uuid_condition_map[condition.Id]

            point_list = [
                self.id_to_uuid_node_map[point.Id]
                for point in condition.GetNodes()
            ]

            face = SFace(
                points=point_list,
                uid=condition_uid
            )

            fid = dst.add_faces([face])

            self.id_to_uuid_condition_map[condition.Id] = fid[0]

    def read_modelpart(self, filename):
        """ Reads a Kratos formated modelpart fro CFD

        """

        model_part = ModelPart("FluidPart")

        model_part.AddNodalSolutionStepVariable(PRESSURE)
        model_part.AddNodalSolutionStepVariable(VELOCITY)
        model_part.AddNodalSolutionStepVariable(VISCOSITY)
        model_part.AddNodalSolutionStepVariable(DENSITY)

        model_part_io_fluid = ModelPartIO(filename)
        model_part_io_fluid.ReadModelPart(model_part)

        smp_meshes = []
        smp_bcs = []

        for i in xrange(0, model_part.NumberOfMeshes()):

            mesh_name = 'fluid_' + str(i)

            smp_mesh = SMesh(name=mesh_name)

            # Export data back to SimPhoNy
            self._exportKratosNodes(
                model_part,
                smp_mesh,
                i
            )
            self._exportKratosElements(
                model_part,
                smp_mesh,
                i
            )
            self._exportKratosConditions(
                model_part,
                smp_mesh,
                i
            )

            data = DataContainer()
            data[CUBA.MATERIAL_ID] = i
            smp_mesh.data = data

            pressure = 'empty'
            velocity = 'empty'

            properties = model_part.GetProperties(0)[i]

            if properties.GetValue(IMPOSED_PRESSURE) == 1:
                pressure = model_part.GetProperties(0)[i].GetValue(PRESSURE)
            if properties.GetValue(IMPOSED_VELOCITY_X) == 1:
                velocity = (
                    properties.GetValue(IMPOSED_VELOCITY_X_VALUE),
                    properties.GetValue(IMPOSED_VELOCITY_Y_VALUE),
                    properties.GetValue(IMPOSED_VELOCITY_Z_VALUE)
                )

            smp_bc = {
                'name': mesh_name,
                'pressure': pressure,
                'velocity': velocity
            }

            smp_bcs.append(smp_bc)
            smp_meshes.append(smp_mesh)

        return {'meshes': smp_meshes, 'bcs': smp_bcs}


class DEM_Utils(object):

    def __init__(self):
        self.id_to_uuid_node_map = {}
        self.uuid_to_id_node_map = {}
        self.id_to_uuid_element_map = {}
        self.uuid_to_id_element_map = {}
        self.id_to_uuid_condition_map = {}
        self.uuid_to_id_condition_map = {}

        self.variables_dictionary = {
            "RADIUS": [
                CUBA.RADIUS,
                RADIUS
            ],
            "DENSITY": [
                CUBA.DENSITY,
                DENSITY
            ],
            "VELOCITY": [
                CUBA.VELOCITY,
                VELOCITY,
                VELOCITY_X,
                VELOCITY_Y,
                VELOCITY_Z
            ],
        }

    def _getSolutionStepVariable1D(self, data, entity, variable):
        pair = self.variables_dictionary[variable]
        if(pair[0] is not None):
            data.update({
                pair[0]: entity.GetSolutionStepValue(pair[1])
            })

    def _getSolutionStepVariable3D(self, data, entity, variable):
        pair = self.variables_dictionary[variable]
        if(pair[0] is not None):
            data.update({
                pair[0]: [
                    entity.GetSolutionStepValue(pair[2]),
                    entity.GetSolutionStepValue(pair[3]),
                    entity.GetSolutionStepValue(pair[4])
                ]
            })

    def _getNodalData(self, data, node):
        """ Extracts the node data

        Extracts the node data and puts in ina format readable
        by the Simphony DataContainer

        """

        self._getSolutionStepVariable1D(data, node, "RADIUS")
        self._getSolutionStepVariable1D(data, node, "DENSITY")
        self._getSolutionStepVariable3D(data, node, "VELOCITY")

    def _exportKratosNodes(self, src, dst, group):
        """ Parses all kratos nodes to simphony points

        Iterates over all nodes in the kratos mesh (src) and
        converts them to simphony points (dst). While doing this operation
        any node/point that has not currently been mapped will have his uuid
        added in the 'id_map' of the wrapper

        """

        for node in src.GetNodes(group):

            data = {}

            self._getNodalData(data, node)

            point_uid = None

            if node.Id in self.id_to_uuid_node_map:
                point_uid = self.id_to_uuid_node_map[node.Id]

            point = SPoint(
                coordinates=(node.X, node.Y, node.Z),
                data=DataContainer(data),
                uid=point_uid
            )

            pid = dst.add_points([point])

            self.id_to_uuid_node_map[node.Id] = pid[0]

    def _exportKratosElements(self, src, dst, group):
        """ Parses all kratos elements to simphony cells

        Iterates over all elements in the kratos mesh (src) and
        converts them to simphony cells (dst). While doing this operation
        any element/cell that has not currently been mapped will have his uuid
        added in the 'id_map' of the wrapper

        """

        for element in src.GetElements(group):

            element_uid = None

            if element.Id in self.id_to_uuid_element_map:
                element_uid = self.id_to_uuid_element_map[element.Id]

            point_list = [
                self.id_to_uuid_node_map[pointl.Id]
                for pointl in element.GetNodes()
            ]

            cell = SCell(
                points=point_list,
                uid=element_uid
            )

            cid = dst.add_cells([cell])

            self.id_to_uuid_element_map[element.Id] = cid[0]

    def _exportKratosConditions(self, src, dst, group):
        """ Parses all kratos conditions to simphony faces

        Iterates over all conditions in the kratos mesh (src) and
        converts them to simphony faces (dst). While doing this operation
        any condition/face that has not currently been mapped will have
        his uuid added in the 'id_map' of the wrapper

        """

        for condition in src.GetConditions(group):

            condition_uid = None

            if condition.Id in self.id_to_uuid_condition_map:
                condition_uid = self.id_to_uuid_condition_map[condition.Id]

            point_list = [
                self.id_to_uuid_node_map[point.Id]
                for point in condition.GetNodes()
            ]

            face = SFace(
                points=point_list,
                uid=condition_uid
            )

            fid = dst.add_faces([face])

            self.id_to_uuid_condition_map[condition.Id] = fid[0]

    def read_modelpart(self, filename, basename):
        """ Reads a Kratos formated modelpart fro CFD

        """

        model_part = ModelPart("FluidPart")

        model_part.AddNodalSolutionStepVariable(RADIUS)
        model_part.AddNodalSolutionStepVariable(DENSITY)
        model_part.AddNodalSolutionStepVariable(VELOCITY)

        model_part_io_fluid = ModelPartIO(filename)
        model_part_io_fluid.ReadModelPart(model_part)

        print(model_part)

        smp_meshes = []
        smp_bcs = []

        for i in xrange(0, model_part.NumberOfMeshes()):

            mesh_name = basename + '_' + str(i)

            smp_mesh = SMesh(name=mesh_name)

            # Export data back to SimPhoNy
            self._exportKratosNodes(
                model_part,
                smp_mesh,
                i
            )
            self._exportKratosElements(
                model_part,
                smp_mesh,
                i
            )
            self._exportKratosConditions(
                model_part,
                smp_mesh,
                i
            )

            data = DataContainer()
            data[CUBA.MATERIAL_ID] = i
            smp_mesh.data = data

            pressure = 'empty'
            velocity = 'empty'

            properties = model_part.GetProperties(0)[i]

            if properties.GetValue(IMPOSED_PRESSURE) == 1:
                pressure = model_part.GetProperties(0)[i].GetValue(PRESSURE)
            if properties.GetValue(IMPOSED_VELOCITY_X) == 1:
                velocity = (
                    properties.GetValue(IMPOSED_VELOCITY_X_VALUE),
                    properties.GetValue(IMPOSED_VELOCITY_Y_VALUE),
                    properties.GetValue(IMPOSED_VELOCITY_Z_VALUE)
                )

            smp_bc = {
                'name': mesh_name,
                'pressure': pressure,
                'velocity': velocity
            }

            smp_bcs.append(smp_bc)
            smp_meshes.append(smp_mesh)

        return {'meshes': smp_meshes, 'bcs': smp_bcs}

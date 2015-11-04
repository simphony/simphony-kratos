""" Template to convert modelparts from kratos to simphony

This file show an example of how to utilze the kratos wrappers
in order to import or export models from KratosMultiphysics

"""

# Kratos will works with Python3 by default.
# This line makes KratosMultiphysics backward
# compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Simphony Imports
from simphony.cuds.abc_modeling_engine import ABCModelingEngine
from simphony.cuds.abc_mesh import ABCMesh
from simphony.cuds.mesh import Mesh
from simphony.core.data_container import DataContainer

# Kratos Imports
# from KratosMultiphysics import *

# Uuid and other dependences


class KratosWrapper(ABCModelingEngine):

    def __init__(self):
        # Uid <-> KratosId mapping
        self.id_to_uuid_node_map = {}
        self.uuid_to_id_node_map = {}
        self.id_to_uuid_element_map = {}
        self.uuid_to_id_element_map = {}
        self.id_to_uuid_condition_map = {}
        self.uuid_to_id_condition_map = {}

        self.id_to_ref_node = {}
        self.id_to_ref_element = {}
        self.id_to_ref_condition = {}

        # Containers
        self.pcs = []
        self.meshes = []
        self.lattices = []

        # Id's stuff
        self.free_id = 1

        # Simphony
        self.CM = DataContainer()
        self.SP = DataContainer()
        self.BC = DataContainer()

        # Extended data containers
        self.SPE = {}

    # ABCModelingEngine Implementation

    def add_dataset(self, container):
        if not isinstance(container, ABCMesh):
            raise TypeError(
                "The type of the dataset container is not supported")
        self._add_mesh(container)

    def _add_mesh(self, src):
        c_mesh = Mesh(name=src.name)

        c_mesh.add_points(src.iter_points())
        c_mesh.add_edges(src.iter_edges())
        c_mesh.add_faces(src.iter_faces())
        c_mesh.add_cells(src.iter_cells())

        c_mesh.data = src.data

        self.meshes.append(c_mesh)

    def remove_dataset(self, name):
        for mesh in self.meshes:
            if mesh.name == name:
                self.meshes.remove(mesh)
                return
        else:
            message = 'Mesh not found'
            raise KeyError(message)

    def get_dataset(self, name):
        for mesh in self.meshes:
            if mesh.name == name:
                return mesh
        else:
            message = 'Mesh not found'
            raise KeyError(message)

    def iter_datasets(self, names=None):
        if names is None:
            for mesh in self.meshes:
                yield mesh
        else:
            mesh_names = [m.name for m in self.meshes]
            for name in names:
                if name not in mesh_names:
                    message = 'Mesh not found'
                    raise KeyError(message)
            for mesh in self.meshes:
                if mesh.name in names:
                    yield mesh

    def get_dataset_names(self):  # pragma: no cover
        """ Returns the names of the all the datasets in the engine workspace.

        """
        return [mesh.name for mesh in self.meshes]

    # KratosWrapper Internal

    def getSolutionStepVariable1D(self, data, entity, variable):
        pair = self.variables_dictionary[variable]
        if(pair[0] is not None):
            data.update({
                pair[0]: entity.GetSolutionStepValue(pair[1])
            })

    def getSolutionStepVariable3D(self, data, entity, variable):
        pair = self.variables_dictionary[variable]
        if(pair[0] is not None):
            data.update({
                pair[0]: [
                    entity.GetSolutionStepValue(pair[2]),
                    entity.GetSolutionStepValue(pair[3]),
                    entity.GetSolutionStepValue(pair[4])
                ]
            })

    def setSolutionStepVariable1D(self, data, entity, variable):
        pair = self.variables_dictionary[variable]
        if(pair[0] is not None):
            entity.SetSolutionStepValue(
                pair[1],
                data[pair[0]]
            )

    def setSolutionStepVariable3D(self, data, entity, variable):
        pair = self.variables_dictionary[variable]
        if(pair[0] is not None):
            for i in xrange(0, 3):
                entity.SetSolutionStepValue(
                    pair[2 + i],
                    data[pair[0]][0 + i]
                )

    def updateForwardDicc(self):
        if not self.uuid_to_id_node_map:
            self.uuid_to_id_node_map = {
                v: k for k, v in self.id_to_uuid_node_map.items()
            }
        if not self.uuid_to_id_element_map:
            self.uuid_to_id_element_map = {
                v: k for k, v in self.id_to_uuid_element_map.items()
            }
        if not self.uuid_to_id_condition_map:
            self.uuid_to_id_condition_map = {
                v: k for k, v in self.id_to_uuid_condition_map.items()
            }

    def updateBackwardDicc(self):
        if not self.id_to_uuid_node_map:
            self.id_to_uuid_node_map = {
                v: k for k, v in self.uuid_to_id_node_map.items()
            }
        if not self.id_to_uuid_element_map:
            self.id_to_uuid_element_map = {
                v: k for k, v in self.uuid_to_id_element_map.items()
            }
        if not self.id_to_uuid_condition_map:
            self.id_to_uuid_condition_map = {
                v: k for k, v in self.uuid_to_id_condition_map.items()
            }

    def run(self, mesh):
        pass

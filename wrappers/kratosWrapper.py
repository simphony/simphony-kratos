""" Template to convert modelparts from kratos to simphony

This file show an example of how to utilze the kratos wrappers
in order to import or export models from KratosMultiphysics

"""

# Kratos will works with Python3 by default.
# This line makes KratosMultiphysics backward
# compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Simphony Imports
from simphony.core.cuba import CUBA
from simphony.core.data_container import DataContainer

from simphony.cuds.mesh import Mesh as SimphonyMesh
from simphony.cuds.mesh import Point as SimphonyPoint
from simphony.cuds.mesh import Edge as SimphonyEdge
from simphony.cuds.mesh import Face as SimphonyFace
from simphony.cuds.mesh import Cell as SimphonyCell

# Kratos Imports
# from KratosMultiphysics import *

# Uuid and other dependences
from uuid import *


class KratosWrapper(object):

    def __init__(self):
        # Uid <-> KratosId mapping
        self.id_to_uuid_node_map = {}
        self.uuid_to_id_node_map = {}
        self.id_to_uuid_element_map = {}
        self.uuid_to_id_element_map = {}
        self.id_to_uuid_condition_map = {}
        self.uuid_to_id_condition_map = {}

        # Containers
        self.pcs = {}
        self.meshes = {}
        self.lattices = {}

        # Id's stuff
        self.free_id = 1

        # Simphony
        self.CM = DataContainer()
        self.SP = DataContainer()
        self.BC = DataContainer()

        # Initialization
        self.initialize()
        
    # ABCModelingEngine Implementation

    def add_particles(self, src):
        self.pcs[src.name] = src

    def add_mesh(self, src):
        self.meshes[src.name] = src

    def add_lattice(self, src):
        self.lattices[src.name] = src

    def delete_particles(self, name):
        self.pcs.remove(name)

    def delete_mesh(self, name):
        self.meshes.remove(name)

    def delete_lattice(self, name):
        self.lattices.remove(name)

    def get_particles(self, name):
        return self.pcs[name]

    def get_mesh(self, name):
        return self.meshes[name]

    def get_lattice(self, name):
        return self.lattices[name]

    def iter_particles(self):
        for pc in self.pcs:
            yield pc

    def iter_meshes(self):
        for mesh in self.meshes:
            yield mesh

    def iter_lattices(self):
        for latt in self.lattices:
            yield latt

    # KratosWrapper Internal

    def addNodalVariablesToModelpart(self):
        pass

    # Small kernels to get ( kratos to simp) and set ( simp to kratos)
    # entity data

    def getSolutionStepVariable1D(self, data, entity, variable):
        pair = self.variables_dictionary[variable]
        data.update({
            pair[0]: entity.GetSolutionStepValue(pair[1])
        })

    def getSolutionStepVariable3D(self, data, entity, variable):
        pair = self.variables_dictionary[variable]
        data.update({
            pair[0]: [
                entity.GetSolutionStepValue(pair[2]),
                entity.GetSolutionStepValue(pair[3]),
                entity.GetSolutionStepValue(pair[4])
            ]
        })

    def setSolutionStepVariable1D(self, data, entity, variable):
        pair = self.variables_dictionary[variable]
        entity.SetSolutionStepValue(
            pair[1],
            data[pair[0]]
        )

    def setSolutionStepVariable3D(self, data, entity, variable):
        pair = self.variables_dictionary[variable]
        for i in xrange(0, 3):
            entity.SetSolutionStepValue(
                pair[2 + i],
                data[pair[0]][0 + i]
            )

    def __getNodalData(self, data, node):
        pass

    def __setNodalData(self, data, node):
        pass

    def __exportKratosNodes(self, src, dst, entitylist=None, data=None):
        pass

    def __exportKratosElements(self, src, dst, entitylist=None, data=None):
        pass

    def __exportKratosConditions(self, src, dst, entitylist=None, data=None):
        pass

    def __importKratosNodes(self, src, dst, entitylist=None, data=None):
        pass

    def __importKratosElements(self, src, dst, entitylist=None, data=None):
        pass

    def __importKratosConditions(self, src, dst, entitylist=None, data=None):
        pass

    def read_modelpart(self, filename):
        pass

    def write_modelpart(self, filename):
        pass

    def add_pc(self, src):
        self.pcs[src.name] = src

    def add_mesh(self, src):
        self.meshes[src.name] = src

    def add_lattice(self, src):
        self.lattices[src.name] = src

    def get_pc(self, name):
        return self.pcs[name]

    def get_mesh(self, name):
        return self.meshes[name]

    def get_lattice(self, name):
        return self.lattices[name]

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

    def setMeshData(self, mesh):
        pass

    def initialize(self):
        pass

    def run(self, mesh):
        pass

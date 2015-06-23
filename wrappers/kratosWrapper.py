""" Template to convert modelparts from kratos to simphony

This file show an example of how to utilze the kratos wrappers
in order to import or export models from KratosMultiphysics

"""

# Kratos will works with Python3 by default.
# This line makes KratosMultiphysics backward
# compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Simphony Imports
from simphony.cuds.mesh import Mesh
from simphony.core.data_container import DataContainer

# Kratos Imports
# from KratosMultiphysics import *

# Uuid and other dependences


class KratosWrapper(object):

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

        # Initialization
        self.initialize()

    # ABCModelingEngine Implementation

    def add_particles(self, particle_container):
        message = 'KratosWrapper does not handle particle container'
        raise NotImplementedError(message)

    def add_mesh(self, src):
        c_mesh = Mesh(name=src.name)

        for p in src.iter_points():
            c_mesh.add_point(p)

        for e in src.iter_edges():
            c_mesh.add_edge(e)

        for f in src.iter_faces():
            c_mesh.add_face(f)

        for c in src.iter_cells():
            c_mesh.add_cell(c)

        c_mesh.data = src.data

        self.meshes.append(c_mesh)

        return c_mesh

    def add_lattice(self, lattice):
        message = 'KratosWrapper does not handle lattice'
        raise NotImplementedError(message)

    def delete_particles(self, name):
        message = 'KratosWrapper does not handle particle container'
        raise NotImplementedError(message)

    def delete_mesh(self, name):
        for mesh in self.meshes:
            if mesh.name == name:
                self.meshes.remove(mesh)
                return
        else:
            message = 'Mesh not found'
            raise KeyError(message)

    def delete_lattice(self, name):
        message = 'KratosWrapper does not handle lattice'
        raise NotImplementedError(message)

    def get_particles(self, name):
        message = 'KratosWrapper does not handle particle container'
        raise NotImplementedError(message)

    def get_mesh(self, name):
        for mesh in self.meshes:
            if mesh.name == name:
                return mesh
        else:
            message = 'Mesh not found'
            raise KeyError(message)

    def get_lattice(self, name):
        message = 'KratosWrapper does not handle lattice'
        raise NotImplementedError(message)

    def iter_particles(self, names=None):
        message = 'KratosWrapper does not handle particle container'
        raise NotImplementedError(message)

    def iter_meshes(self, names=None):
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

    def iter_lattices(self, names=None):
        message = 'KratosWrapper does not handle lattice'
        raise NotImplementedError(message)

    # KratosWrapper Internal

    def addNodalVariablesToModelpart(self):
        pass

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

    def getNodalData(self, data, node):
        pass

    def setNodalData(self, data, node):
        pass

    def exportKratosNodes(self, src, dst):
        pass

    def exportKratosElements(self, src, dst):
        pass

    def exportKratosConditions(self, src, dst):
        pass

    def exportKratosDof(self, src, dst):
        pass

    def importKratosNodes(self, src, dst):
        pass

    def importKratosElements(self, src, dst):
        pass

    def importKratosConditions(self, src, dst):
        pass

    def importKratosDof(self, src, dst):
        pass

    def read_modelpart(self, filename):
        pass

    def write_modelpart(self, filename):
        pass

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

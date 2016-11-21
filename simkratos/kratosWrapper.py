""" KRATOS Common Simphony Wrapper

This module provides the common functionality for Kratos based
wrappers

"""

# Kratos will works with Python3 by default.
# This line makes KratosMultiphysics backward
# compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

import itertools

# Simphony Imports
from simphony.cuds.abc_modeling_engine import ABCModelingEngine
from simphony.cuds.abc_mesh import ABCMesh
from simphony.cuds.abc_particles import ABCParticles
from simphony.cuds.mesh import Mesh
from simphony.cuds.particles import Particles
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

        if any(container.name == ps.name for ps in self.pcs):
            raise ValueError(
                'Particles \'{n}\` already exists'.format(container.name))

        if any(container.name == ms.name for ms in self.meshes):
            raise ValueError(
                'Mesh \'{n}\` already exists'.format(container.name))

        if isinstance(container, ABCMesh):
            self._add_mesh(container)
        elif isinstance(container, ABCParticles):
            self._add_particles(container)
        else:
            raise TypeError(
                "The type of the dataset container is not supported")

    def _add_mesh(self, src):
        c_mesh = Mesh(name=src.name)

        c_mesh.add_points(src.iter_points())
        c_mesh.add_edges(src.iter_edges())
        c_mesh.add_faces(src.iter_faces())
        c_mesh.add_cells(src.iter_cells())

        c_mesh.data = src.data

        self.meshes.append(c_mesh)

    def _add_particles(self, src):
        c_particles = Particles(name=src.name)

        c_particles.add_particles(src.iter_particles())
        c_particles.add_bonds(src.iter_bonds())

        c_particles.data = src.data

        self.pcs.append(c_particles)

    def remove_dataset(self, name):
        if name in [m.name for m in self.meshes]:
            for mesh in self.meshes:
                if mesh.name == name:
                    self.meshes.remove(mesh)
                    return
        elif name in [p.name for p in self.pcs]:
            for pc in self.pcs:
                if pc.name == name:
                    self.pcs.remove(pc)
                    return
        else:
            message = 'Dataset not found'
            raise KeyError(message)

    def get_dataset(self, name):
        if name in [m.name for m in self.meshes]:
            for mesh in self.meshes:
                if mesh.name == name:
                    return mesh
        elif name in [p.name for p in self.pcs]:
            for pc in self.pcs:
                if pc.name == name:
                    return pc
        else:
            message = 'Dataset not found'
            raise KeyError(message)

    def iter_datasets(self, names=None):

        ip = self._iter_particles(names)
        im = self._iter_meshes(names)

        iter_list = [i for i in ip] + [i for i in im]

        if names is not None:
            for name in names:
                if name not in [i.name for i in iter_list]:
                    raise ValueError(
                        'Container \'{n}\' does not exists'.format(n=name))
        for i in iter_list:
            yield i

    def get_dataset_names(self):  # pragma: no cover
        """ Returns the names of the all the datasets in the engine workspace.

        """

        ip = self._iter_particles()
        im = self._iter_meshes()

        iter_list = itertools.chain(ip, im)

        return [i.name for i in iter_list]

    # private

    def _iter_meshes(self, names=None):
        if names is None:
            for mesh in self.meshes:
                yield mesh
        else:
            for mesh in [mesh for mesh in self.meshes]:
                if mesh.name in names:
                    yield mesh

    def _iter_particles(self, names=None):
        if names is None:
            for pc in self.pcs:
                yield pc
        else:
            for pc in [pc.name for pc in self.pcs]:
                if pc.name in names:
                    yield pc

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

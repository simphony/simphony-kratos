""" KRATOS Common Simphony Wrapper

This module provides the common functionality for Kratos based
wrappers

"""

# Kratos works with Python3 by default.
# This line makes KratosMultiphysics backward
#  compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

import itertools

# Simphony Imports
from simphony.cuds.abc_modeling_engine import ABCModelingEngine
from simphony.cuds.abc_mesh import ABCMesh
from simphony.cuds.abc_particles import ABCParticles
from simphony.cuds.mesh import Mesh
from simphony.cuds.particles import Particles
from simphony.core.data_container import DataContainer


class KratosWrapper(ABCModelingEngine):
    """ BaseClass for Kratos based wrappers

    """

    def __init__(self, use_internal_interface=True, **kwargs):
        """ Constructor

        Parameters
        ----------
        use_internal_interface: bool, optional
            Must be always true. It indicates that the internal
            interface is used for comminication with Kratos.
            An Exception is raised otherwise

        """

        # Uid <-> KratosId mapping
        if use_internal_interface is False:
            raise Exception('Kratos wrappers cannot be File-io based')

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
        self.BC = DataContainer()
        self.CM = DataContainer()
        self.SP = DataContainer()
        self.CM_extension = {}
        self.SP_extension = {}
        self.BC_extension = {}

        # Call the base class in order to load CUDS
        super(KratosWrapper, self).__init__(**kwargs)

    def add_dataset(self, container):
        """Add a CUDS container

        Parameters
        ----------
        container : {ABCMesh, ABCParticles}
            The CUDS container to add to the engine.

        Raises
        ------
        TypeError:
            If the container type is not supported (i.e. ABCLattice, ABCMesh).
        ValueError:
            If there is already a dataset with the given name.
        """

        if any(container.name == ps.name for ps in self.pcs):
            raise ValueError(
                'Particle container \'{n}\' already exists'.format(
                    container.name
                )
            )

        if any(container.name == ms.name for ms in self.meshes):
            raise ValueError(
                'Mesh \'{n}\` already exists'.format(
                    container.name
                )
            )

        if isinstance(container, ABCMesh):
            self._add_mesh(container)
        elif isinstance(container, ABCParticles):
            self._add_particles(container)
        else:
            raise TypeError(
                "The type of the dataset container is not supported")

    def get_dataset(self, name):
        """ Get the dataset

        The returned dataset can be used to query
        and change the related data inside Kratos.

        Parameters
        ----------
        name: str
            name of CUDS container to be retrieved.

        Returns
        -------
        container :
            A proxy of the dataset named ``name`` that is stored
            internally in the Engine.

        Raises
        ------
        ValueError:
            If there is no dataset with the given name

        """

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

    def get_dataset_names(self):  # pragma: no cover
        """ Returns the names of the all the datasets in the engine workspace.

        """

        ip = self._iter_particles()
        im = self._iter_meshes()

        iter_list = itertools.chain(ip, im)

        return [i.name for i in iter_list]

    def remove_dataset(self, name):
        """ Remove a dataset

        Parameters
        ----------
        name: str
            name of CUDS container to be deleted

        Raises
        ------
        ValueError:
            If there is no dataset with the given name

        """

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

    def iter_datasets(self, names=None):
        """ Returns an iterator over a subset or all of the containers.

        Parameters
        ----------
        names : sequence of str, optional
            names of specific containers to be iterated over. If names is not
            given, then all containers will be iterated over.

        """

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

    # KratosWrapper Internal

    def getSolutionStepVariable1D(self, data, entity, variable):
        """ Updates simphony data 'data' with the scalar 'variable'
        of a Kratos entity

        Parameters
        ----------
        data: DataContainer
            The data container to be updated

        entity: Kratos node, element, condition
            Source of the data

        variable: Enum of kratos variables
            CUBA to be updated
        """

        pair = self.variables_dictionary[variable]
        if(pair[0] is not None):
            data.update({
                pair[0]: entity.GetSolutionStepValue(pair[1])
            })

    def getSolutionStepVariable3D(self, data, entity, variable):
        """ Updates simphony data 'data' with the vectorial 'variable'
        of a Kratos entity

        Parameters
        ----------
        data: DataContainer
            The data container to be updated

        entity: Kratos node, element, condition
            Source of the data

        variable: Enum of kratos variables
            CUBA to be updated
        """

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
        """ Updates a kratos entity scalar variable with a simphony
        variable

        Parameters
        ----------
        data: DataContainer
            The data container where the source is located

        entity: Kratos node, element, condition
            Destination of the data

        variable: Enum of kratos variables
            Identifier of the variable to be updated
        """

        pair = self.variables_dictionary[variable]
        if(pair[0] is not None):
            entity.SetSolutionStepValue(
                pair[1],
                data[pair[0]]
            )

    def setSolutionStepVariable3D(self, data, entity, variable):
        """ Updates a kratos entity scalar variable with a simphony
        variable

        Parameters
        ----------
        data: DataContainer
            The data container where the source is located

        entity: Kratos node, element, condition
            Destination of the data

        variable: Enum of kratos variables
            Identifier of the variable to be updated
        """

        pair = self.variables_dictionary[variable]
        if(pair[0] is not None):
            for i in xrange(0, 3):
                entity.SetSolutionStepValue(
                    pair[2 + i],
                    data[pair[0]][0 + i]
                )

    def updateForwardDicc(self):
        """ Updates the internal mapping of the wrapper

        Updates the internal mapping of the wrapper using the Simphony to
        Kratos mapping as a source
        """

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
        """ Updates the internal mapping of the wrapper

        Updates the internal mapping of the wrapper using the Kratos to
        Simphony mapping as a source
        """

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

    def _add_mesh(self, src):
        """Adds a mesh

        Parameters
        ----------
        container : {ABCMesh}
            The mesh to add to the engine.

        """

        c_mesh = Mesh(name=src.name)

        c_mesh.add_points(src.iter_points())
        c_mesh.add_edges(src.iter_edges())
        c_mesh.add_faces(src.iter_faces())
        c_mesh.add_cells(src.iter_cells())

        c_mesh.data = src.data

        self.meshes.append(c_mesh)

    def _add_particles(self, src):
        """Adds a particle container

        Parameters
        ----------
        container : {ABCParticles}
            The particle container to add to the engine.

        """

        c_particles = Particles(name=src.name)

        c_particles.add_particles(src.iter_particles())
        c_particles.add_bonds(src.iter_bonds())

        c_particles.data = src.data

        self.pcs.append(c_particles)

    def _iter_meshes(self, names=None):
        """ Returns an iterator over a subset or all of the meshes.

        Parameters
        ----------
        names : sequence of str, optional
            names of specific meshes to be iterated over. If names is not
            given, then all meshes will be iterated over.

        """

        if names is None:
            for mesh in self.meshes:
                yield mesh
        else:
            for mesh in [mesh for mesh in self.meshes]:
                if mesh.name in names:
                    yield mesh

    def _iter_particles(self, names=None):
        """ Returns an iterator over a subset or all of the particle containers.

        Parameters
        ----------
        names : sequence of str, optional
            names of specific particle containers to be iterated over.
            If names is not given, then all particle containers will be
            iterated over.

        """

        if names is None:
            for pc in self.pcs:
                yield pc
        else:
            for pc in [pc.name for pc in self.pcs]:
                if pc.name in names:
                    yield pc

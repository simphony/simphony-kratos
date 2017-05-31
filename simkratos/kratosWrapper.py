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
from simphony.core.cuba import CUBA
from simphony.cuds.abc_modeling_engine import ABCModelingEngine
from simphony.cuds.abc_mesh import ABCMesh
from simphony.cuds.abc_particles import ABCParticles
from simphony.cuds.mesh import Mesh
from simphony.cuds.particles import Particles
from simphony.core.data_container import DataContainer
from simphony.cuds.mesh import Point as SPoint
from simphony.cuds.mesh import Face as SFace
from simphony.cuds.mesh import Cell as SCell
from simphony.cuds.particles import Particle as SParticle

# Kratos Imports
import KratosMultiphysics as KRTS


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

        if use_internal_interface is False:
            raise Exception('Kratos wrappers cannot be File-io based')

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

        # Entity count
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

    def _load_cuds(self):
        """Load CUDS data into lammps engine."""
        pass

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
            if pair[0] not in data.keys():
                data[pair[0]] = 0
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
            if pair[0] not in data.keys():
                data[pair[0]] = (0, 0, 0)
            for i in xrange(0, 3):
                entity.SetSolutionStepValue(
                    pair[2 + i],
                    data[pair[0]][0 + i]
                )

    def setProperty(self, data, propertyContainer, variable):
        """ Updates a kratos property with a simphony material
        variable

        Parameters
        ----------
        data: DataContainer
            The data container where the source is located

        propertyContainer: Kratos property container
            Destination of the data

        variable: Enum of kratos variables
            Identifier of the variable to be updated
        """

        pair = self.properties_dictionary[variable]
        if(pair[0] is not None):
            if pair[0] not in data.keys():
                data[pair[0]] = 0
            propertyContainer.SetValue(
                pair[1],
                data[pair[0]]
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

    # Common

    def exportKratosNodes(self, src, dst, group):
        """ Parses all kratos nodes to simphony points

        Iterates over all nodes in the kratos mesh (src) and
        converts them to simphony points (dst). While doing this operation
        any node/point that has not currently been mapped will have his uuid
        added in the 'id_map' of the wrapper

        """

        for node in src.GetNodes(group):

            point_uid = None

            if node.Id in self.id_to_uuid_node_map.keys():
                point_uid = self.id_to_uuid_node_map[node.Id]

            hasNode = point_uid is not None and dst.has(point_uid)

            if not hasNode:

                data = {}
                self.getNodalData(data, node, src.Name)
                
                point = SPoint(
                    coordinates=(node.X, node.Y, node.Z),
                    data=DataContainer(data),
                    uid=point_uid
                )

                pid = dst.add([point])

                self.id_to_uuid_node_map[node.Id] = pid[0]

            else:

                point = dst.get(
                    uid=self.id_to_uuid_node_map[node.Id]
                )

                data = point.data
                self.getNodalData(data, node, src.Name)

                point.data = DataContainer(data)
                point.coordinates = (node.X, node.Y, node.Z)

                dst.update([point])

    def exportKratosElements(self, src, dst, group):
        """ Parses all kratos elements to simphony cells

        Iterates over all nodes in the kratos mesh (src) and
        converts them to simphony cells (dst). While doing this operation
        any point that has not currently mapped will have his uuid
        added in the 'id_map' of the weapper

        """

        for element in src.GetElements(group):

            element_uid = None

            data = {}

            if element.Id not in self.id_to_uuid_element_map:

                point_list = [
                    self.id_to_uuid_node_map[pointl.Id]
                    for pointl in element.GetNodes()
                ]

                cell = SCell(
                    points=point_list,
                    data=DataContainer(data),
                    uid=element_uid
                )

                cid = dst.add([cell])

                self.id_to_uuid_element_map[element.Id] = cid[0]

            else:

                # No data is stored in the element yet

                pass

    def exportKratosConditions(self, src, dst, group):
        """ Parses all kratos conditions to simphony faces

        Iterates over all nodes in the kratos mesh ( src ) and
        converts them to simphony faces (dst). While doing this operation
        any point that has not currently mapped will have his uuid
        added in the 'id_map' of the weapper

        """

        for condition in src.GetConditions(group):

            condition_uid = None

            data = {}

            if condition.Id not in self.id_to_uuid_condition_map:

                point_list = [
                    self.id_to_uuid_node_map[point.Id]
                    for point in condition.GetNodes()
                ]

                face = SFace(
                    points=point_list,
                    data=DataContainer(data),
                    uid=condition_uid
                )

                fid = dst.add([face])

                self.id_to_uuid_condition_map[condition.Id] = fid[0]

            else:

                # No data is stored in the condition yet

                pass

    def exportKratosParticles(self, src, dst, group):
        """ Parses all kratos nodes to simphony Particles

        Iterates over all nodes in the kratos mesh (src) and
        converts them to simphony Particles (dst). While doing this operation
        any node/point that has not currently been mapped will have his uuid
        added in the 'id_map' of the wrapper. Notice that Kratos Element
        information will be not updated as Particles in Simphony are stand
        alone objects

        """

        for node in src.GetNodes(group):

            particle_uid = None

            if node.Id not in self.id_to_uuid_node_map:

                data = {}
                self.getNodalData(data, node, src.Name)

                particle = SParticle(
                    coordinates=(node.X, node.Y, node.Z),
                    data=DataContainer(data),
                    uid=particle_uid
                )

                pid = dst.add([particle])

                self.id_to_uuid_node_map[node.Id] = pid[0]

            else:

                particle = dst.get(
                    uid=self.id_to_uuid_node_map[node.Id]
                )

                data = particle.data
                self.getNodalData(data, node, src.Name)

                particle.data = DataContainer(data)
                particle.coordinates = (node.X, node.Y, node.Z)

                dst.update([particle])

    def importKratosNodes(self, src, dst, group):
        """ Parses all simphony points to kratos nodes

        Iterates over all points in the simphony mesh (src) and
        converts them to kratos nodes (dst).
        While doing this operation any point/node pair that has not
        currently mapped will have  his uuid added in the 'id_map'
        of the wrapper

        """

        # Add the points in case they don't exist and update their value in
        # case they do.
        for point in src.iter(item_type=CUBA.POINT):

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

                self.setNodalData(data, node, dst.Name)

                self.id_to_ref_node[node_id] = node

            else:

                node = self.id_to_ref_node[self.uuid_to_id_node_map[point.uid]]
                data = point.data
                self.setNodalData(data, node, dst.Name)

        # If they belong to a different group, add them
        if group != 0:
            nodes = KRTS.NodesArray()
            for point in src.iter(item_type=CUBA.POINT):
                nodes.append(
                    dst.Nodes[self.uuid_to_id_node_map[point.uid]]
                )
            dst.GetMesh(group)
            dst.SetNodes(nodes, group)

    def importKratosElements(self, src, dst, group, element_type):
        """ Parses all simphony cells to kratos elements

        Iterates over all cells in the simphony mesh (src) and
        converts them to kratos SphericContinuumParticle3D elements (dst).
        While doing this operation any point/node pair that has not
        currently mapped will have  his uuid added in the 'id_map'
        of the wrapper

        """

        for element in src.iter(item_type=CUBA.CELL):

            if element.uid not in self.uuid_to_id_element_map.keys():
                self.uuid_to_id_element_map.update(
                    {element.uid: self.free_id}
                )

                self.free_id += 1

                element_id = self.uuid_to_id_element_map[element.uid]

                property_range = 0
                if group in self.kratos_props.keys():
                    property_range = group

                dst.CreateNewElement(
                    element_type,
                    element_id,
                    [self.uuid_to_id_node_map[p] for p in element.points],
                    self.kratos_props[property_range]
                )

        # If they belong to a different group, add them
        if group != 0:
            elements = KRTS.ElementsArray()
            for elem in src.iter(item_type=CUBA.CELL):
                elements.append(
                    dst.Elements[self.uuid_to_id_element_map[elem.uid]]
                )
            dst.GetMesh(group)
            dst.SetElements(elements, group)

    def importKratosConditions(self, src, dst, group, condition_type):
        """ Parses all simphony faces to kratos conditions

        Iterates over all faces in the simphony mesh (src) and
        converts them to kratos RigidFace3D3N conditions (dst).
        While doing this operation any point/node pair that has not
        currently mapped will have  his uuid added in the 'id_map'
        of the wrapper

        """

        for condition in src.iter(item_type=CUBA.FACE):

            if condition.uid not in self.uuid_to_id_condition_map.keys():
                self.uuid_to_id_condition_map.update(
                    {condition.uid: self.free_id}
                )

                self.free_id += 1

                condition_id = self.uuid_to_id_condition_map[condition.uid]

                property_range = 0
                if group in self.kratos_props.keys():
                    property_range = group

                dst.CreateNewCondition(
                    condition_type,
                    condition_id,
                    [self.uuid_to_id_node_map[p] for p in condition.points],
                    self.kratos_props[property_range]
                )

        # If they belong to a different group, add them
        if group != 0:
            conditions = KRTS.ConditionsArray()
            for cnd in src.iter(item_type=CUBA.FACE):
                conditions.append(
                    dst.Conditions[self.uuid_to_id_condition_map[cnd.uid]]
                )
            dst.GetMesh(group)
            dst.SetConditions(conditions, group)

    def importKratosParticles(self, src, dst, group, particle_type):
        """ Parses all simphony particles to kratos particle elements

        Iterates over all particles in the simphony mesh (src) and
        converts them to kratos particle elements (dst).
        While doing this operation any elements that has not
        currently mapped will have  his uuid added in the 'id_map'
        of the wrapper

        """

        for particle in src.iter(item_type=CUBA.PARTICLE):

            data = particle.data

            if particle.uid not in self.uuid_to_id_node_map.keys():
                self.uuid_to_id_node_map.update(
                    {particle.uid: self.free_id}
                )

                self.free_id += 1

                node_id = self.uuid_to_id_node_map[particle.uid]

                node = dst.CreateNewNode(
                    node_id,
                    particle.coordinates[0],
                    particle.coordinates[1],
                    particle.coordinates[2])

                self.setNodalData(data, node, dst.Name)

                self.id_to_ref_node[node_id] = node

            else:

                node = self.id_to_ref_node[
                    self.uuid_to_id_node_map[particle.uid]
                ]

                data = particle.data

                self.setNodalData(data, node, dst.Name)

            if particle.uid not in self.uuid_to_id_element_map.keys():
                self.uuid_to_id_element_map.update(
                    {particle.uid: self.free_id}
                )

                self.free_id += 1

                element_id = self.uuid_to_id_element_map[particle.uid]

                property_range = 0
                if group in self.kratos_props.keys():
                    property_range = group

                dst.CreateNewElement(
                    particle_type,
                    element_id,
                    [self.uuid_to_id_node_map[particle.uid]],
                    self.kratos_props[property_range]
                )

        # If they belong to a different group, add them
        if group != 0:
            nodes = KRTS.NodesArray()
            elements = KRTS.ElementsArray()
            for particle in src.iter(item_type=CUBA.PARTICLE):
                nodes.append(
                    dst.Nodes[self.uuid_to_id_node_map[particle.uid]]
                )
                elements.append(
                    dst.Elements[self.uuid_to_id_element_map[particle.uid]]
                )
            dst.GetMesh(group)
            dst.SetNodes(nodes, group)
            dst.SetElements(elements, group)

    # Run

    def run(self, mesh):
        pass

    # Undefined

    def _add_mesh(self, src):
        """Adds a mesh

        Parameters
        ----------
        container : {ABCMesh}
            The mesh to add to the engine.

        """

        c_mesh = Mesh(name=src.name)

        c_mesh.add(src.iter(item_type=CUBA.POINT))
        c_mesh.add(src.iter(item_type=CUBA.EDGE))
        c_mesh.add(src.iter(item_type=CUBA.FACE))
        c_mesh.add(src.iter(item_type=CUBA.CELL))

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

        c_particles.add(src.iter(item_type=CUBA.PARTICLE))
        c_particles.add(src.iter(item_type=CUBA.BOND))

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

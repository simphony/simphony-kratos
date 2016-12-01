from simphony.core.cuba import CUBA
from simphony.cuds.meta import api
from simphony.core.data_container import DataContainer

from simphony.cuds.mesh import Point as SPoint
from simphony.cuds.mesh import Mesh as SMesh
from simphony.cuds.mesh import Face as SFace
from simphony.cuds.mesh import Cell as SCell
from simphony.cuds.particles import Particle as SParticle
from simphony.cuds.particles import Particles as SParticles

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

        self.supportedMaterialProp = {
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

            pid = dst.add([point])

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

            cid = dst.add([cell])

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

            fid = dst.add([face])

            self.id_to_uuid_condition_map[condition.Id] = fid[0]

    def _convertBc(self, properties, mesh_name):
        condition = api.Condition(name='condition_' + mesh_name)
        conditionData = condition.data

        velocityCondition = [None, None, None]

        if properties.GetValue(IMPOSED_PRESSURE) == 1:
            conditionData[CUBA.PRESSURE] = properties.GetValue(
                PRESSURE
            )

        if properties.GetValue(IMPOSED_VELOCITY_X) == 1:
            velocityCondition[0] = properties.GetValue(
                IMPOSED_VELOCITY_X_VALUE
            )
        if properties.GetValue(IMPOSED_VELOCITY_Y) == 1:
            velocityCondition[1] = properties.GetValue(
                IMPOSED_VELOCITY_Y_VALUE
            )
        if properties.GetValue(IMPOSED_VELOCITY_Z) == 1:
            velocityCondition[2] = properties.GetValue(
                IMPOSED_VELOCITY_Z_VALUE
            )

        if not (None in velocityCondition):
            conditionData[CUBA.VELOCITY] = velocityCondition

        condition.data = conditionData
        return condition

    def _convertMaterial(self, properties, mesh_name):
        material = api.Material(name='material' + mesh_name)
        materialData = material.data

        for key, val in self.supportedMaterialProp:
            if properties.GetValue(val):
                materialData[key] = properties.GetValue(val)

        material.data = materialData
        return material

    def read_modelpart(self, filename):
        """ Reads a Kratos formated modelpart for KratosCFD wrapper

        This functions translates each of the modelpart meshes to
        one simphony mesh.

        Properties in the modelpart related to boundary conditions are
        stored in simphony conditions and referenced in CUBA.CONDITION
        inside the mesh.

        Properties in the modelpart not relared to boundary conditions
        are stored in simphony as materials and referenced in
        CUBA.MATERIAL inside the mesh.

        """

        model_part = ModelPart("FluidPart")

        model_part.AddNodalSolutionStepVariable(PRESSURE)
        model_part.AddNodalSolutionStepVariable(VELOCITY)
        model_part.AddNodalSolutionStepVariable(VISCOSITY)
        model_part.AddNodalSolutionStepVariable(DENSITY)

        model_part_io_fluid = ModelPartIO(filename)
        model_part_io_fluid.ReadModelPart(model_part)

        smp_meshes = []
        smp_conditions = []
        smp_materials = []

        for i in xrange(1, model_part.NumberOfMeshes()):

            mesh_name = 'fluid_' + str(i)

            mesh = SMesh(name=mesh_name)

            # Export mesh data to Simphony
            self._exportKratosNodes(model_part, mesh, i)
            self._exportKratosElements(model_part, mesh, i)
            self._exportKratosConditions(model_part, mesh, i)

            properties = model_part.GetProperties(0)[i]

            # Fill the boundary condition for the mesh
            condition = self._convertBc(properties, mesh_name)

            # Fill the material for the mesh
            material = self._convertMaterial(properties, mesh_name)

            # Save the relations
            meshData = mesh.data
            meshData[CUBA.CONDITION] = condition.name
            meshData[CUBA.MATERIAL] = material.name
            mesh.data = meshData

            # Pack the return objects
            smp_meshes.append(mesh)
            smp_conditions.append(condition)
            smp_materials.append(material)

        return {
            'datasets': smp_meshes,
            'conditions': smp_conditions,
            'materials': smp_materials,
        }


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

        self.supportedMaterialProp = {
            CUBA.DENSITY: PARTICLE_DENSITY,
            CUBA.YOUNG_MODULUS: YOUNG_MODULUS,
            CUBA.POISSON_RATIO: POISSON_RATIO,
            CUBA.FRICTION_COEFFICIENT: PARTICLE_FRICTION,
            CUBA.ROLLING_FRICTION: ROLLING_FRICTION
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

            pid = dst.add([point])

            self.id_to_uuid_node_map[node.Id] = pid[0]

    def _exportKratosParticles(self, src, dst, group):
        """ Parses all kratos nodes to simphony Particles

        Iterates over all nodes in the kratos mesh (src) and
        converts them to simphony Particles (dst). While doing this operation
        any node/point that has not currently been mapped will have his uuid
        added in the 'id_map' of the wrapper. Notice that Kratos Element
        information will be not updated as Particles in Simphony are stand
        alone objects

        """

        for particle in src.GetNodes(group):

            data = {}

            self._getNodalData(data, particle)

            particle_uid = None

            if particle.Id in self.id_to_uuid_node_map:
                particle_uid = self.id_to_uuid_node_map[particle.Id]

            sparticle = SParticle(
                coordinates=(particle.X, particle.Y, particle.Z),
                data=DataContainer(data),
                uid=particle_uid
            )

            pid = dst.add([sparticle])

            self.id_to_uuid_node_map[particle.Id] = pid[0]

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

            cid = dst.add([cell])

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

            fid = dst.add([face])

            self.id_to_uuid_condition_map[condition.Id] = fid[0]

    def _convertBc(self, properties, mesh_name):
        condition = api.Condition(name='condition_' + mesh_name)
        conditionData = condition.data

        velocityCondition = [None, None, None]

        if properties.GetValue(IMPOSED_PRESSURE) == 1:
            conditionData[CUBA.PRESSURE] = properties.GetValue(
                PRESSURE
            )

        if properties.GetValue(IMPOSED_VELOCITY_X) == 1:
            velocityCondition[0] = properties.GetValue(
                IMPOSED_VELOCITY_X_VALUE
            )
        if properties.GetValue(IMPOSED_VELOCITY_Y) == 1:
            velocityCondition[1] = properties.GetValue(
                IMPOSED_VELOCITY_Y_VALUE
            )
        if properties.GetValue(IMPOSED_VELOCITY_Z) == 1:
            velocityCondition[2] = properties.GetValue(
                IMPOSED_VELOCITY_Z_VALUE
            )

        if not (None in velocityCondition):
            conditionData[CUBA.VELOCITY] = velocityCondition

        condition.data = conditionData
        return condition

    def _convertMaterial(self, properties, mesh_name):
        material = api.Material(name='material' + mesh_name)
        materialData = material.data

        for key, val in self.supportedMaterialProp.items():
            if properties.GetValue(val):
                materialData[key] = properties.GetValue(val)

        material.data = materialData
        return material

    def read_modelpart_as_mesh(self, filename):
        """ Reads a Kratos formated modelpart from DEM

        """

        model_part = ModelPart("FluidPart")

        model_part.AddNodalSolutionStepVariable(RADIUS)
        model_part.AddNodalSolutionStepVariable(DENSITY)
        model_part.AddNodalSolutionStepVariable(VELOCITY)

        model_part_io_fluid = ModelPartIO(filename)
        model_part_io_fluid.ReadModelPart(model_part)

        smp_meshes = []
        smp_conditions = []
        smp_materials = []

        for i in xrange(1, model_part.NumberOfMeshes()):

            mesh_name = 'solid_' + str(i)

            mesh = SMesh(name=mesh_name)

            # Export mesh data to Simphony
            self._exportKratosNodes(model_part, mesh, i)
            self._exportKratosElements(model_part, mesh, i)
            self._exportKratosConditions(model_part, mesh, i)

            properties = model_part.GetProperties(0)[i]

            # Fill the boundary condition for the mesh
            condition = self._convertBc(properties, mesh_name)

            # Fill the material for the mesh
            material = self._convertMaterial(properties, mesh_name)

            # Save the relations
            meshData = mesh.data
            meshData[CUBA.CONDITION] = condition.name
            meshData[CUBA.MATERIAL] = material.name
            mesh.data = meshData

            # Pack the return objects
            smp_meshes.append(mesh)
            smp_conditions.append(condition)
            smp_materials.append(material)

        return {
            'datasets': smp_meshes,
            'conditions': smp_conditions,
            'materials': smp_materials,
        }

    def read_modelpart_as_particles(self, filename):
        """ Reads a Kratos formated modelpart from DEM

        """

        model_part = ModelPart("FluidPart")

        model_part.AddNodalSolutionStepVariable(RADIUS)
        model_part.AddNodalSolutionStepVariable(DENSITY)
        model_part.AddNodalSolutionStepVariable(VELOCITY)

        model_part_io_fluid = ModelPartIO(filename)
        model_part_io_fluid.ReadModelPart(model_part)

        smp_particles = []
        smp_conditions = []
        smp_materials = []

        for i in xrange(1, model_part.NumberOfMeshes()):

            particles_name = 'particles_' + str(i)

            particles = SParticles(name=particles_name)

            # Export particle data to Simphony
            self._exportKratosParticles(model_part, particles, i)

            properties = model_part.GetProperties(0)[i]

            # Fill the boundary condition for the mesh
            condition = self._convertBc(properties, particles_name)

            # Fill the material for the mesh
            material = self._convertMaterial(properties, particles_name)

            # Save the relations
            particlesData = particles.data
            particlesData[CUBA.CONDITION] = condition.name
            particlesData[CUBA.MATERIAL] = material.name
            particles.data = particlesData

            # Pack the return objects
            smp_particles.append(particles)
            smp_conditions.append(condition)
            smp_materials.append(material)

        return {
            'datasets': smp_particles,
            'conditions': smp_conditions,
            'materials': smp_materials,
        }

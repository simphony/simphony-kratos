""" TODO: Description was out of date

"""

# Kratos works with Python3 by default.
from __future__ import print_function, absolute_import, division

# Simphony Imports
from simphony.core.cuba import CUBA
from simphony.core.data_container import DataContainer

from simphony.cuds.mesh import Point as SPoint
from simphony.cuds.mesh import Face as SFace
from simphony.cuds.mesh import Cell as SCell
from simphony.cuds.particles import Particle as SParticle

# Wrapper Imports
from simkratos.kratosWrapper import KratosWrapper
from simkratos.DEM import DEM_explicit_solver_var as DEM_parameters

# Kratos Imports
from KratosMultiphysics import *                                                # noqa: F403
from KratosMultiphysics.DEMApplication import *                                 # noqa: F403



import sphere_strategy as SolverStrategy
import DEM_procedures


class DEMWrapper(KratosWrapper):

    def __init__(self, use_internal_interface=True, **kwargs):
        super(DEMWrapper, self).__init__(use_internal_interface, **kwargs)

        self.time = 0
        self.step = 0
        self.substeps = 0

        # The dictionary defines the relation between CUBA and
        # kratos variables

        self.variables_dictionary = {
            "RADIUS": [
                CUBA.RADIUS,
                RADIUS
            ],
            "NODAL_MASS": [
                None,
                NODAL_MASS
            ],
            "VELOCITY": [
                CUBA.VELOCITY,
                VELOCITY,
                VELOCITY_X,
                VELOCITY_Y,
                VELOCITY_Z
            ],
            "DISPLACEMENT": [
                None,
                DISPLACEMENT,
                DISPLACEMENT_X,
                DISPLACEMENT_Y,
                DISPLACEMENT_Z
            ],
            "TOTAL_FORCES": [
                None,
                TOTAL_FORCES,
                TOTAL_FORCES_X,
                TOTAL_FORCES_Y,
                TOTAL_FORCES_Z
            ]
        }

        self.initialize()

    def _load_cuds(self):
        """Load CUDS data into lammps engine."""
        cuds = self.get_cuds()
        if not cuds:
            return

        for component in cuds.iter(item_type=CUBA.MESH):
            self.add_dataset(component)

    def getNodalData(self, data, node, model):
        """ Extracts the node data

        Extracts the node data and puts in ina format readeable
        by the Simphony DataContainer

        """

        if model == "Particles":
            self.getSolutionStepVariable1D(data, node, "RADIUS")
            self.getSolutionStepVariable1D(data, node, "NODAL_MASS")
            self.getSolutionStepVariable3D(data, node, "VELOCITY")
            self.getSolutionStepVariable3D(data, node, "DISPLACEMENT")
            self.getSolutionStepVariable3D(data, node, "TOTAL_FORCES")

    def setNodalData(self, data, node, model):
        """ Assembles the point data

        Assembles the point data and puts in ina format readeable
        by the Kratos ModelPart

        """

        if model == "Particles":
            self.setSolutionStepVariable1D(data, node, "RADIUS")
            self.setSolutionStepVariable1D(data, node, "NODAL_MASS")
            self.setSolutionStepVariable3D(data, node, "VELOCITY")
            self.setSolutionStepVariable3D(data, node, "DISPLACEMENT")
            self.setSolutionStepVariable3D(data, node, "TOTAL_FORCES")

    def exportKratosNodes(self, src, dst, group):
        """ Parses all kratos nodes to simphony points

        Iterates over all nodes in the kratos mesh (src) and
        converts them to simphony points (dst). While doing this operation
        any node/point that has not currently been mapped will have his uuid
        added in the 'id_map' of the wrapper

        """

        for node in src.GetNodes(group):

            data = {}

            self.getNodalData(data, node, src.Name)

            point_uid = None

            if node.Id not in self.id_to_uuid_node_map:

                point = SPoint(
                    coordinates=(node.X, node.Y, node.Z),
                    data=DataContainer(data),
                    uid=point_uid
                )

                pid = dst.add([point])

                self.id_to_uuid_node_map[node.Id] = pid[0]

            else:

                point = dst.get(uid=self.id_to_uuid_node_map[node.Id])

                point.data = DataContainer(data)

                dst.update_points([point])

    def exportKratosParticles(self, src, dst, group):
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

            self.getNodalData(data, particle, src.Name)

            particle_uid = None

            if particle.Id not in self.id_to_uuid_node_map:

                particle = SParticle(
                    coordinates=(particle.X, particle.Y, particle.Z),
                    data=DataContainer(data),
                    uid=particle_uid
                )

                pid = dst.add([particle])

                self.id_to_uuid_node_map[particle.Id] = pid[0]

            else:

                particle = dst.get(
                    uid=self.id_to_uuid_node_map[particle.Id]
                )

                particle.data = DataContainer(data)

                dst.update_particles([particle])

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

                fid = dst.add(face)

                self.id_to_uuid_condition_map[condition.Id] = fid[0]

            else:

                # No data is stored in the condition yet

                pass

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
            nodes = NodesArray()
            for point in src.iter(item_type=CUBA.POINT):
                nodes.append(
                    dst.Nodes[self.uuid_to_id_node_map[point.uid]]
                )
            dst.SetNodes(nodes, group)

    def importKratosParticles(self, src, dst, group):
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

                dst.CreateNewElement(
                    "SphericParticle3D",
                    element_id,
                    [self.uuid_to_id_node_map[particle.uid]],
                    self.element_properties)

        # If they belong to a different group, add them
        if group != 0:
            nodes = NodesArray()
            elements = ElementsArray()
            for particle in src.iter(item_type=CUBA.PARTICLE):
                nodes.append(
                    dst.Nodes[self.uuid_to_id_node_map[particle.uid]]
                )
                elements.append(
                    dst.Elements[self.uuid_to_id_element_map[particle.uid]]
                )
            dst.SetNodes(nodes, group)
            dst.SetElements(elements, group)

    def importKratosElements(self, src, dst, group):
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

                dst.CreateNewElement(
                    "SphericParticle3D",
                    element_id,
                    [self.uuid_to_id_node_map[p] for p in element.points],
                    self.element_properties)

        # If they belong to a different group, add them
        if group != 0:
            elements = ElementsArray()
            for elem in src.iter(item_type=CUBA.CELL):
                elements.append(
                    dst.Elements[self.uuid_to_id_element_map[elem.uid]]
                )
            dst.SetElements(elements, group)

    def importKratosConditions(self, src, dst, group):
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

                dst.CreateNewCondition(
                    "RigidFace3D3N",
                    condition_id,
                    [self.uuid_to_id_node_map[p] for p in condition.points],
                    self.condition_properties)

        # If they belong to a different group, add them
        if group != 0:
            conditions = ConditionsArray()
            for cnd in src.iter(item_type=CUBA.FACE):
                conditions.append(
                    dst.Conditions[self.uuid_to_id_condition_map[cnd.uid]]
                )
            dst.SetConditions(conditions, group)

    def _setMeshData(self):
        " This probably needs to be done throug configuration"

        cLawString = "DEMContinuumConstitutiveLaw"
        dLawString = "DEMDiscontinuumConstitutiveLaw"

        self.SP[CUBA.DENSITY] = 2500.0
        self.SP[CUBA.YOUNG_MODULUS] = 1.0e5
        self.SP[CUBA.POISSON_RATIO] = 0.20
        self.SP[CUBA.ROLLING_FRICTION] = 0.01

        self.PARTICLE_FRICTION = 0.99
        self.PARTICLE_COHESION = 0.0
        self.LN_OF_RESTITUTION_COEFF = -1.6094379124341003
        self.PARTICLE_MATERIAL = 1
        self.WALL_FRICTION = 0.3
        self.DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME = cLawString
        self.DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME = dLawString

    def setElementData(self):
        self.element_properties.SetValue(
            PARTICLE_DENSITY,
            self.SP[CUBA.DENSITY]
        )
        self.element_properties.SetValue(
            YOUNG_MODULUS,
            self.SP[CUBA.YOUNG_MODULUS]
        )
        self.element_properties.SetValue(
            POISSON_RATIO,
            self.SP[CUBA.POISSON_RATIO]
        )
        self.element_properties.SetValue(
            PARTICLE_FRICTION,
            self.PARTICLE_FRICTION
        )
        self.element_properties.SetValue(
            PARTICLE_COHESION,
            self.PARTICLE_COHESION
        )
        self.element_properties.SetValue(
            LN_OF_RESTITUTION_COEFF,
            self.LN_OF_RESTITUTION_COEFF
        )
        self.element_properties.SetValue(
            PARTICLE_MATERIAL,
            self.PARTICLE_MATERIAL
        )
        self.element_properties.SetValue(
            ROLLING_FRICTION,
            self.SP[CUBA.ROLLING_FRICTION]
        )
        self.element_properties.SetValue(
            DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME,
            self.DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME
        )
        self.element_properties.SetValue(
            DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME,
            self.DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME
        )

    def setConditionData(self):
        cmesh = self.get_dataset("Conditions")

        self.condition_properties.SetValue(
            WALL_FRICTION,
            cmesh.data[CUBA.WALL_FRICTION]
        )

    def initialize(self):
        """ Initalizes the necessary kratos commponents

            Initalizes the necessary kratos commponents to
            execute Kratos' DEMPack solver
        """

        # DEMPack SubClasses
        self.procedures = DEM_procedures.Procedures(DEM_parameters)
        self.parallelutils = DEM_procedures.ParallelUtils()
        self.materialTest = DEM_procedures.MaterialTest()
        self.creator_destructor = ParticleCreatorDestructor()

        # Prepare ModelParts
        self.spheres_model_part = ModelPart("Particles")
        self.rigid_face_model_part = ModelPart("Conditions")
        self.mixed_model_part = ModelPart("")
        self.cluster_model_part = ModelPart("")
        self.DEM_inlet_model_part = ModelPart("")
        self.mapping_model_part = ModelPart("")
        self.contact_model_part = ModelPart("")

        # Create solver
        self.solver = SolverStrategy.ExplicitStrategy(
            self.spheres_model_part,
            self.rigid_face_model_part,
            self.cluster_model_part,
            self.DEM_inlet_model_part,
            self.creator_destructor,
            DEM_parameters
        )

        # Prepare properties
        self.element_properties = Properties(0)
        self.condition_properties = Properties(1)

        self._setMeshData()
        self.setElementData()

        # Prepare variables
        self.solver.AddAdditionalVariables(
            self.spheres_model_part,
            DEM_parameters
        )
        self.procedures.AddCommonVariables(
            self.spheres_model_part,
            DEM_parameters
        )
        self.procedures.AddSpheresVariables(
            self.spheres_model_part,
            DEM_parameters
        )
        self.procedures.AddCommonVariables(
            self.rigid_face_model_part,
            DEM_parameters
        )
        self.procedures.AddRigidFaceVariables(
            self.rigid_face_model_part,
            DEM_parameters
        )

        # Set a search strategy
        self.solver.search_strategy = self.parallelutils.GetSearchStrategy(
            self.solver,
            self.spheres_model_part
        )

        self.solver.Initialize()

    def run(self):
        """ Run a step of the wrapper """

        fluid_particles = self.pcs
        solid_meshes = self.meshes

        cuds = self.get_cuds()

        self.spheres_model_part.GetMesh(len(fluid_particles))
        self.rigid_face_model_part.GetMesh(len(solid_meshes))

        fluid_properties = PropertiesArray()
        meshNumber = 1
        meshDict = {}

        for particles in cuds.iter(item_type=CUBA.PARTICLE):

            group = meshNumber

            self.importKratosParticles(
                particles,
                self.spheres_model_part,
                group
            )

            meshDict[mesh.name] = meshNumber
            meshNumber += 1

        self.updateBackwardDicc()

        fluid_properties.append(self.element_properties)
        # solid_properties.append(self.condition_properties)

        self.spheres_model_part.SetProperties(fluid_properties)
        # self.rigid_face_model_part.SetProperties(solid_properties)

        SolverStrategy.AddDofs(self.spheres_model_part)

        self.solver.Initialize()

        self.dt = cuds.get_by_name('dem_integration_time').step

        # Start the simulation itself
        self.time = cuds.get_by_name('dem_integration_time').time
        self.final = cuds.get_by_name('dem_integration_time').final

        # Solve
        while self.time < self.final:

            self.dt = self.spheres_model_part.ProcessInfo.GetValue(DELTA_TIME)
            cuds.get_by_name('dem_integration_time').step = self.dt

            self.spheres_model_part.ProcessInfo[TIME] = self.time
            self.spheres_model_part.ProcessInfo[DELTA_TIME] = self.dt
            self.spheres_model_part.ProcessInfo[TIME_STEPS] = self.step

            self.rigid_face_model_part.ProcessInfo[TIME] = self.time
            self.rigid_face_model_part.ProcessInfo[DELTA_TIME] = self.dt
            self.rigid_face_model_part.ProcessInfo[TIME_STEPS] = self.step

            self.solver.Solve()

            self.step += 1
            self.time = self.time + self.dt

        cuds.get_by_name('dem_integration_time').time = self.time
        cuds.get_by_name('dem_integration_time').final = self.final

        for particles in cuds.iter(item_type=CUBA.PARTICLE):

            group = meshDict[mesh.name]

            self.exportKratosParticles(
                self.spheres_model_part,
                particles,
                group
            )

        self.updateForwardDicc()

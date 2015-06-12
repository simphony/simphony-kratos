""" TODO: Description was out of date

"""

# Kratos works with Python3 by default.
from __future__ import print_function, absolute_import, division

# Simphony Imports
from simphony.core.cuba import CUBA
from simphony.core.data_container import DataContainer

from simphony.cuds.mesh import Point as SPoint
from simphony.cuds.mesh import Mesh as SMesh
from simphony.cuds.mesh import Face as SFace
from simphony.cuds.mesh import Cell as SCell

# Wrapper Imports
from wrappers.kratosWrapper import KratosWrapper
from wrappers.cuba_extension import CUBAExt

# Kratos Imports
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

import DEM_explicit_solver_var as DEM_parameters
import sphere_strategy as SolverStrategy
import DEM_procedures


class DEMPackWrapper(KratosWrapper):

    def __init__(self):
        KratosWrapper.__init__(self)

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

    def exportKratosNodes(self, src, dst):
        """ Parses all kratos nodes to simphony points

        Iterates over all nodes in the kratos mesh (src) and
        converts them to simphony points (dst). While doing this operation
        any node/point that has not currently been mapped will have his uuid
        added in the 'id_map' of the wrapper

        """

        for node in src.GetNodes():

            data = {}

            # this is provisional
            data[CUBA.LABEL] = src.Name

            self.getNodalData(data, node, src.Name)

            point_uid = None

            if node.Id not in self.id_to_uuid_node_map:

                point = SPoint(
                    coordinates=(node.X, node.Y, node.Z),
                    data=DataContainer(data),
                    uid=point_uid
                )

                pid = dst.add_point(point)

                self.id_to_uuid_node_map[node.Id] = pid

            else:

                point = dst.get_point(uid=self.id_to_uuid_node_map[node.Id])

                # iterate over the correct data
                point.data = DataContainer(data)

                dst.update_point(point)

    def exportKratosElements(self, src, dst):
        """ Parses all kratos elements to simphony cells

        Iterates over all nodes in the kratos mesh (src) and
        converts them to simphony cells (dst). While doing this operation
        any point that has not currently mapped will have his uuid
        added in the 'id_map' of the weapper

        """

        for element in src.GetElements():

            element_uid = None

            data = {}

            # self.getElementalData(data, node)

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

                cid = dst.add_cell(cell)

                self.id_to_uuid_element_map[element.Id] = cid

            else:

                # No data is stored in the element yet

                pass

    def exportKratosConditions(self, src, dst):
        """ Parses all kratos conditions to simphony faces

        Iterates over all nodes in the kratos mesh ( src ) and
        converts them to simphony faces (dst). While doing this operation
        any point that has not currently mapped will have his uuid
        added in the 'id_map' of the weapper

        """

        for condition in src.GetConditions():

            condition_uid = None

            data = {}

            # self.getElementalData(data, node)

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

                fid = dst.add_face(face)

                self.id_to_uuid_condition_map[condition.Id] = fid

            else:

                # No data is stored in the condition yet

                pass

    def importKratosNodes(self, src, dst):
        """ Parses all simphony points to kratos nodes

        Iterates over all points in the simphony mesh (src) and
        converts them to kratos nodes (dst).
        While doing this operation any point/node pair that has not
        currently mapped will have  his uuid added in the 'id_map'
        of the wrapper

        """

        for point in src.iter_points():

            data = point.data

            # this condition is provisional
            if data[CUBA.LABEL] == dst:

                if point.uid not in self.uuid_to_id_node_map.keys():
                    self.uuid_to_id_node_map.update(
                        {point.uid: self.free_id}
                    )

                    self.free_id += 1

                    node_id = self.uuid_to_id_node_map[point.uid]

                    node = dst.CreateNewNode(
                        node_id,
                        point.coordinates[0],
                        point.coordinates[1],
                        point.coordinates[2])

                    self.setNodalData(data, node, dst.Name)

                    self.id_to_ref_node[node_id] = node

                else:

                    node = self.id_to_ref_node[
                        self.uuid_to_id_node_map[point.uid]
                    ]

                    data = point.data

                    self.setNodalData(data, node)

    def importKratosElements(self, src, dst):
        """ Parses all simphony cells to kratos elements

        Iterates over all cells in the simphony mesh (src) and
        converts them to kratos SphericContinuumParticle3D elements (dst).
        While doing this operation any point/node pair that has not
        currently mapped will have  his uuid added in the 'id_map'
        of the wrapper

        """

        for element in src.iter_cells():

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

    def importKratosConditions(self, src, dst):
        """ Parses all simphony faces to kratos conditions

        Iterates over all faces in the simphony mesh (src) and
        converts them to kratos RigidFace3D3N conditions (dst).
        While doing this operation any point/node pair that has not
        currently mapped will have  his uuid added in the 'id_map'
        of the wrapper

        """

        for condition in src.iter_faces():

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

    def read_modelpart(self, fluid_filename, rigid_face_filename):
        """ Reads a Kratos formated modelpart

        This adds partial support for the future FileIO
        """

        new_mesh = SMesh(name="Model")

        model_part_io_fluid = ModelPartIO(fluid_filename)
        model_part_io_fluid.ReadModelPart(self.spheres_model_part)

        model_part_io_cnd = ModelPartIO(rigid_face_filename)
        model_part_io_cnd.ReadModelPart(self.rigid_face_model_part)

        # Setting up the buffer size
        self.spheres_model_part.SetBufferSize(1)
        self.rigid_face_model_part.SetBufferSize(1)

        # Set dofs
        SolverStrategy.AddDofs(self.spheres_model_part)

        self.element_properties = Properties(0)
        self.condition_properties = Properties(1)

        self.elemNodeData = {CUBA.RADIUS: RADIUS}
        self.condNodeData = {}

        for n in self.spheres_model_part.Nodes:
            self.id_to_ref_node[n.Id] = n

        # Add the problem data
        self.setMeshData(new_mesh)

        # Export data back to SimPhoNy
        self.exportKratosNodes(
            self.spheres_model_part,
            new_mesh
        )
        self.exportKratosElements(
            self.spheres_model_part,
            new_mesh
        )

        self.exportKratosNodes(
            self.rigid_face_model_part,
            new_mesh
        )
        self.exportKratosConditions(
            self.rigid_face_model_part,
            new_mesh
        )

        self.solver.Initialize()

        self.updateForwardDicc()

        return new_mesh

    def write_modelpart(self, filename):
        """ Writes a Kratos formated modelpart WIP

        This adds partial support for the future FileIO
        """

        pmesh = self.get_mesh("Particles")

        f = open(filename, 'w')

        # Variable information
        f.write('Begin ModelPartData\n')
        f.write('//  VARIABLE_NAME value\n')
        f.write('End ModelPartData\n\n')

        # Properties ( out shared values )
        f.write('Begin Properties 1\n')
        f.write('PARTICLE_DENSITY {}\n'.format(pmesh.data[PARTICLE_DENSITY]))
        f.write('YOUNG_MODULUS {}\n'.format(pmesh.data[YOUNG_MODULUS]))
        f.write('POISSON_RATIO {}\n'.format(pmesh.data[POISSON_RATIO]))
        f.write('PARTICLE_FRICTION {}\n'.format(pmesh.data[PARTICLE_FRICTION]))
        f.write('PARTICLE_COHESION {}\n'.format(pmesh.data[PARTICLE_COHESION]))
        f.write('LN_OF_RESTITUTION_COEFF {}\n'.format(
            pmesh.data[LN_OF_RESTITUTION_COEFF]))
        f.write('PARTICLE_MATERIAL {}\n'.format(pmesh.data[PARTICLE_MATERIAL]))
        f.write('ROLLING_FRICTION {}\n'.format(pmesh.data[ROLLING_FRICTION]))
        f.write('DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME {}\n'.format(
            pmesh.data[DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME]))
        f.write('DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME {}\n'.format(
            pmesh.data[DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME]))
        f.write('End Properties\n\n')

        # Nodes
        f.write('Begin Nodes\n')
        for point in src.iter_points():
            f.write('{} {} {} {}\n').format(
                self.uuid_to_id_point_map[point.id],
                point.coordinates[0],
                point.coordinates[1],
                point.coordinates[2]
            )
        f.write('End Nodes\n\n')

        f.write('Begin Elements SphericParticle3D')
        for element in src.iter_cells():
            f.write('{} {} {}\n').format(
                self.uuid_to_id_point_map[elemet.id],
                self.uuid_to_id_point_map[element.points[0]],
                self.uuid_to_id_point_map[element.points[1]]
            )
        f.write('End Elements\n\n')

        f.write('Begin NodalData RADIUS')
        for point in src.iter_points():
            f.write('{} {} {} {}\n').format(
                self.uuid_to_id_point_map[point.id],
                0,
                point.data[CUBA.RADIUS],
            )
        f.write('End NodalData\n\n')

    def setMeshData(self, mesh):
        " This probably needs to be done throug configuration"

        cLawString = "DEMContinuumConstitutiveLaw"
        dLawString = "DEMDiscontinuumConstitutiveLaw"

        self.SP[CUBA.DENSITY] = 2500.0
        self.SP[CUBA.YOUNG_MODULUS] = 1.0e5
        self.SP[CUBA.POISSON_RATIO] = 0.20
        self.SPE[CUBAExt.PARTICLE_FRICTION] = 0.99
        self.SPE[CUBAExt.PARTICLE_COHESION] = 0.0
        self.SPE[CUBAExt.LN_OF_RESTITUTION_COEFF] = -1.6094379124341003
        self.SPE[CUBAExt.PARTICLE_MATERIAL] = 1
        self.SP[CUBA.ROLLING_FRICTION] = 0.01
        self.SPE[CUBAExt.WALL_FRICTION] = 0.3
        self.SPE[CUBAExt.DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME] = cLawString
        self.SPE[CUBAExt.DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME] = dLawString

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
            self.SPE[CUBAExt.PARTICLE_FRICTION]
        )
        self.element_properties.SetValue(
            PARTICLE_COHESION,
            self.SPE[CUBAExt.PARTICLE_COHESION]
        )
        self.element_properties.SetValue(
            LN_OF_RESTITUTION_COEFF,
            self.SPE[CUBAExt.LN_OF_RESTITUTION_COEFF]
        )
        self.element_properties.SetValue(
            PARTICLE_MATERIAL,
            self.SPE[CUBAExt.PARTICLE_MATERIAL]
        )
        self.element_properties.SetValue(
            ROLLING_FRICTION,
            self.SP[CUBA.ROLLING_FRICTION]
        )
        self.element_properties.SetValue(
            DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME,
            self.SPE[CUBAExt.DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME]
        )
        self.element_properties.SetValue(
            DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME,
            self.SPE[CUBAExt.DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME]
        )

    def setConditionData(self):
        cmesh = self.get_mesh("Conditions")

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
        self.demio = DEM_procedures.DEMIo()

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

        # Import the into Kratos
        self.importKratosNodes(
            self.get_mesh("Model"),
            self.spheres_model_part
        )
        self.importKratosElements(
            self.get_mesh("Model"),
            self.spheres_model_part
        )

        self.importKratosNodes(
            self.get_mesh("Model"),
            self.spheres_model_part
        )
        self.importKratosConditions(
            self.get_mesh("Model"),
            self.rigid_face_model_part
        )

        self.updateBackwardDicc()
        self.setElementData()
        # self.setConditionData()

        step = 0
        time = 0.0

        # Solve
        for n in xrange(0, self.CM[CUBA.NUMBER_OF_TIME_STEPS]):

            dt = self.spheres_model_part.ProcessInfo.GetValue(DELTA_TIME)
            time = time + dt
            step += 1

            self.spheres_model_part.ProcessInfo[TIME] = time
            self.spheres_model_part.ProcessInfo[DELTA_TIME] = dt
            self.spheres_model_part.ProcessInfo[TIME_STEPS] = step

            self.rigid_face_model_part.ProcessInfo[TIME] = time
            self.rigid_face_model_part.ProcessInfo[DELTA_TIME] = dt
            self.rigid_face_model_part.ProcessInfo[TIME_STEPS] = step

            self.solver.Solve()

            time += dt

        self.exportKratosNodes(
            self.spheres_model_part,
            self.get_mesh("Model")
        )
        self.exportKratosElements(
            self.spheres_model_part,
            self.get_mesh("Model")
        )

        self.exportKratosNodes(
            self.rigid_face_model_part,
            self.get_mesh("Model")
        )
        self.exportKratosConditions(
            self.rigid_face_model_part,
            self.get_mesh("Model")
        )

        self.updateForwardDicc()

    def finalize(self):
        self.demio.FinalizeMesh()

""" TODO: Description was out of date

"""

# Kratos works with Python3 by default.
from __future__ import print_function, absolute_import, division

# Simphony Imports
from simphony.core.cuba import CUBA
from simphony.core.data_container import DataContainer

from simphony.cuds.mesh import Mesh as SMesh
from simphony.cuds.mesh import Point as SPoint
from simphony.cuds.mesh import Edge as SEdge
from simphony.cuds.mesh import Face as SFace
from simphony.cuds.mesh import Cell as SCell

from uuid import *

# Wrapper Imports
from kratosWrapper import KratosWrapper

# Kratos Imports 
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *

import fractional_step_solver
import decoupled_solver_eulerian
import monolithic_solver_eulerian
import monolithic_solver_eulerian_compressible

class KratosWrapper(KratosWrapper):

    def __init__(self):
        super(KratosWrapper,self).__init__()

    def __addNodalVariablesToModelpart(self):
        """ Adds the DEMPack nodal variables

        Adds the DEMPack nodal variables to the ratos modelpart provided 
        in order to be usable later while importing the mesh.

        """

        self.procedures.AddCommonVariables(self.spheres_model_part, DEM_parameters)
        self.procedures.AddSpheresVariables(self.spheres_model_part, DEM_parameters)
        
        SolverStrategy.AddAdditionalVariables(self.spheres_model_part, DEM_parameters)

        self.procedures.AddCommonVariables(self.rigid_face_model_part, DEM_parameters)
        self.procedures.AddRigidFaceVariables(self.rigid_face_model_part, DEM_parameters)

    def __exportKratosElements(self,src,dst,entitylist=None):
        """ Parses all kratos elements to simphony cells

        Iterates over all nodes in the kratos mesh (src) and
        converts them to simphony cells (dst). While doing this operation
        any point that has not currently mapped will have his uuid
        added in the 'id_map' of the weapper

        """

        for element in src.GetElements():

            point_list = [self.id_to_uuid_node_map[point.Id] for point in element.GetNodes()]

            for node in element.GetNodes():
        
                data = {
                    CUBA.RADIUS: node.GetSolutionStepValue(RADIUS),
                    CUBA.NODAL_MASS: node.GetSolutionStepValue(NODAL_MASS),
                    CUBA.VELOCITY:  [
                        node.GetSolutionStepValue(VELOCITY_X),
                        node.GetSolutionStepValue(VELOCITY_Y),
                        node.GetSolutionStepValue(VELOCITY_Z)
                    ],
                    CUBA.DISPLACEMENT:  [
                        node.GetSolutionStepValue(DISPLACEMENT_X),
                        node.GetSolutionStepValue(DISPLACEMENT_Y),
                        node.GetSolutionStepValue(DISPLACEMENT_Z)
                    ],
                    CUBA.DELTA_DISPLACEMENT:  [
                        node.GetSolutionStepValue(DELTA_DISPLACEMENT_X),
                        node.GetSolutionStepValue(DELTA_DISPLACEMENT_Y),
                        node.GetSolutionStepValue(DELTA_DISPLACEMENT_Z)
                    ],
                    CUBA.TOTAL_FORCES:  [
                        node.GetSolutionStepValue(TOTAL_FORCES_X),
                        node.GetSolutionStepValue(TOTAL_FORCES_Y),
                        node.GetSolutionStepValue(TOTAL_FORCES_Z)
                    ]
                }

                point = SPoint(
                    coordinates=(node.X, node.Y, node.Z),
                    data=DataContainer(data)
                )

                point_uuid = dst.add_point(point)

                if point_uuid not in self.id_to_uuid_node_map.keys():
                    self.id_to_uuid_node_map.update(
                        {node.Id:point_uuid}
                    )

            cell = SCell(
                points=point_list,
                data=DataContainer(data)
            )

            element_uuid = dst.add_cell(cell)

            if element_uuid not in self.id_to_uuid_element_map.keys():
                self.id_to_uuid_element_map.update(
                    {element.Id:element_uuid}
                )

    def __exportKratosConditions(self,src,dst,entitylist=None):
        """ Parses all kratos conditions to simphony faces

        Iterates over all nodes in the kratos mesh ( src ) and
        converts them to simphony faces (dst). While doing this operation
        any point that has not currently mapped will have his uuid
        added in the 'id_map' of the weapper

        """

        for condition in src.GetConditions():

            point_list = [self.id_to_uuid_node_map[point.Id] for point in condition.GetNodes()]

            for node in src.GetNodes():
        
                data = {}

                point = SPoint(
                    coordinates=(node.X, node.Y, node.Z),
                    data=DataContainer(data)
                )

                point_uuid = dst.add_point(point)

                if point_uuid not in self.id_to_uuid_node_map.keys():
                    self.id_to_uuid_node_map.update(
                        {node.Id:point_uuid}
                    )

            face = SFace(
                points=point_list,
                data=DataContainer(data)
            )

            condition_uuid = dst.add_face(face)

            if condition_uuid not in self.id_to_uuid_condition_map.keys():
                self.id_to_uuid_condition_map.update(
                    {condition.Id:condition_uuid}
                )

    def __importKratosElements(self,src,dst,entitylist=None):
        """ Parses all simphony cells to kratos elements

        Iterates over all cells in the simphony mesh (src) and
        converts them to kratos SphericContinuumParticle3D elements (dst). 
        While doing this operation any point/node pair that has not 
        currently mapped will have  his uuid added in the 'id_map' 
        of the wrapper

        """
        
        for element in src.iter_cells(entitylist):

            if element.uuid not in self.uuid_to_id_element_map.keys():
                self.uuid_to_id_element_map.update(
                    {element.uuid:self.free_id}
                )

                self.free_id += 1

            element_id = self.uuid_to_id_element_map[element.uuid]

            for point in src.iter_points(element.points):

                if point.uuid not in self.uuid_to_id_node_map.keys():
                    self.uuid_to_id_node_map.update(
                        {point.uuid:self.free_id}
                    )

                    self.free_id += 1

                node_id = self.uuid_to_id_node_map[point.uuid]

                node = dst.CreateNewNode(
                    node_id,
                    point.coordinates[0],
                    point.coordinates[1],
                    point.coordinates[2])

                node.SetSolutionStepValue(RADIUS,point.data[CUBA.RADIUS])
                node.SetSolutionStepValue(NODAL_MASS,point.data[CUBA.NODAL_MASS])
                node.SetSolutionStepValue(VELOCITY_X,point.data[CUBA.VELOCITY][0])
                node.SetSolutionStepValue(VELOCITY_Y,point.data[CUBA.VELOCITY][1])
                node.SetSolutionStepValue(VELOCITY_Z,point.data[CUBA.VELOCITY][2])
                node.SetSolutionStepValue(DISPLACEMENT_X,point.data[CUBA.DISPLACEMENT][0])
                node.SetSolutionStepValue(DISPLACEMENT_Y,point.data[CUBA.DISPLACEMENT][1])
                node.SetSolutionStepValue(DISPLACEMENT_Z,point.data[CUBA.DISPLACEMENT][2])
                node.SetSolutionStepValue(DELTA_DISPLACEMENT_X,point.data[CUBA.DELTA_DISPLACEMENT][0])
                node.SetSolutionStepValue(DELTA_DISPLACEMENT_Y,point.data[CUBA.DELTA_DISPLACEMENT][1])
                node.SetSolutionStepValue(DELTA_DISPLACEMENT_Z,point.data[CUBA.DELTA_DISPLACEMENT][2])
                node.SetSolutionStepValue(TOTAL_FORCES_X,point.data[CUBA.TOTAL_FORCES][0])
                node.SetSolutionStepValue(TOTAL_FORCES_Y,point.data[CUBA.TOTAL_FORCES][1])
                node.SetSolutionStepValue(TOTAL_FORCES_Z,point.data[CUBA.TOTAL_FORCES][2])

            dst.CreateNewElement(
                "Fluid2D",
                element_id,
                [self.uuid_to_id_node_map[p] for p in element.points],
                self.element_properties)

    def __importKratosConditions(self,src,dst,entitylist=None):
        """ Parses all simphony faces to kratos conditions

        Iterates over all faces in the simphony mesh (src) and
        converts them to kratos RigidFace3D3N conditions (dst). 
        While doing this operation any point/node pair that has not 
        currently mapped will have  his uuid added in the 'id_map' 
        of the wrapper

        """
        pass 

    #### AAAA INTERAL WRAPPER AAAA 
    #### ||||                 ||||   ( Can be splited, I guess... )
    #### VVVV GENERAL WRAPPER VVVV

    def parseSp(self):
        self.element_properties.SetValue(PARTICLE_DENSITY,self.SP[PARTICLE_DENSITY])
        self.element_properties.SetValue(YOUNG_MODULUS,self.SP[YOUNG_MODULUS])
        self.element_properties.SetValue(POISSON_RATIO,self.SP[POISSON_RATIO])
        self.element_properties.SetValue(PARTICLE_FRICTION,self.SP[PARTICLE_FRICTION])
        self.element_properties.SetValue(PARTICLE_COHESION,self.SP[PARTICLE_COHESION])
        self.element_properties.SetValue(LN_OF_RESTITUTION_COEFF,self.SP[LN_OF_RESTITUTION_COEFF])
        self.element_properties.SetValue(PARTICLE_MATERIAL,self.SP[PARTICLE_MATERIAL])
        self.element_properties.SetValue(ROLLING_FRICTION,self.SP[ROLLING_FRICTION])
        self.element_properties.SetValue(DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME,self.SP[DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME])
        self.element_properties.SetValue(DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME,self.SP[DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME])

    def initialize(self):
        """ Initalizes the necessary kratos commponents

            Initalizes the necessary kratos commponents to 
            execute Kratos - DEMPack solver
        """

        # DEMPack SubClasses
        self.procedures = DEM_procedures.Procedures(DEM_parameters)
        self.parallelutils = DEM_procedures.ParallelUtils()
        self.materialTest = DEM_procedures.MaterialTest()
        self.creator_destructor = ParticleCreatorDestructor()

        # Prepare ModelParts
        self.spheres_model_part = ModelPart("")
        self.rigid_face_model_part = ModelPart("")
        self.mixed_model_part = ModelPart("")
        self.cluster_model_part = ModelPart("")
        self.DEM_inlet_model_part = ModelPart("")
        self.mapping_model_part = ModelPart("")
        self.contact_model_part = ModelPart("")

        self.element_properties = Properties(0)
        self.condition_properties = Properties(1)

        self.condition_properties.SetValue(WALL_FRICTION,self.SP[WALL_FRICTION])

        self.elemNodeData = {CUBA.RADIUS:RADIUS}
        self.condNodeData = {}

    def run(self):
        """ Run a step of the wrapper """ 

        self.spheres_model_part = ModelPart("")
        self.rigid_face_model_part = ModelPart("")

        self.spheres_model_part.Properties.append(self.element_properties)
        self.rigid_face_model_part.Properties.append(self.condition_properties)

        self.__addNodalVariablesToModelpart()

        # Import the data to Kratos
        self.__importKratosElements(
            self.get_mesh("CDF"), self.spheres_model_part
        )

        self.__importKratosConditions(
            self.get_mesh("CDF"), self.rigid_face_model_part
        )

        self.updateBackwardDicc()
        self.parseSp()

        SolverStrategy.AddDofs(self.spheres_model_part)

        # Create solver
        if(self.CM[CUDS.NAME] == "fractional_step"):
            fluid_solver = fractional_step_solver.IncompressibleFluidSolver(
                fluid_model_part, domain_size)
            fluid_solver.laplacian_form = laplacian_form
            # standard laplacian form
            fluid_solver.predictor_corrector = False
            fluid_solver.max_press_its = fluid_only_var.max_press_its
            fluid_solver.velocity_linear_solver = SkylineLUFactorizationSolver()
            fluid_solver.pressure_linear_solver = SkylineLUFactorizationSolver()

            fluid_solver.Initialize()
        elif(self.CM[CUDS.NAME] == "pressure_splitting"):
            fluid_solver = decoupled_solver_eulerian.\
                DecoupledSolver(fluid_model_part, domain_size)
            oss_switch = fluid_only_var.use_oss
            dynamic_tau = fluid_only_var.dynamic_tau
        # pPrecond = ILU0Preconditioner()
            pPrecond = DiagonalPreconditioner()
            fluid_solver.pressure_linear_solver = BICGSTABSolver(1e-3, 5000, pPrecond)
        # fluid_solver.linear_solver =  SuperLUSolver()
            fluid_solver.rel_vel_tol = 1e-4
            fluid_solver.abs_vel_tol = 1e-6
            fluid_solver.rel_pres_tol = 1e-4
            fluid_solver.abs_pres_tol = 1e-6
            fluid_solver.use_inexact_newton = False
            fluid_model_part.ProcessInfo.SetValue(OSS_SWITCH, oss_switch)
            fluid_model_part.ProcessInfo.SetValue(DYNAMIC_TAU, dynamic_tau)
            fluid_solver.Initialize()
        elif(self.CM[CUDS.NAME] == "monolithic_solver_eulerian"):
            fluid_solver = monolithic_solver_eulerian.MonolithicSolver(
                fluid_model_part, domain_size)
            oss_switch = fluid_only_var.use_oss
            dynamic_tau = fluid_only_var.dynamic_tau
            fluid_model_part.ProcessInfo.SetValue(OSS_SWITCH, oss_switch)
            fluid_model_part.ProcessInfo.SetValue(DYNAMIC_TAU, dynamic_tau)
            fluid_solver.Initialize()
        elif(self.CM[CUDS.NAME] == "monolithic_solver_eulerian_compressible"):
            fluid_solver = monolithic_solver_eulerian_compressible.MonolithicSolver(
                fluid_model_part, domain_size)
            oss_switch = fluid_only_var.use_oss
            dynamic_tau = fluid_only_var.dynamic_tau
            fluid_model_part.ProcessInfo.SetValue(OSS_SWITCH, oss_switch)
            fluid_model_part.ProcessInfo.SetValue(DYNAMIC_TAU, dynamic_tau)
            fluid_solver.Initialize()

        # Not sure what to do here :S
        self.spheres_model_part.ProcessInfo[TIME] = self.time
        self.spheres_model_part.ProcessInfo[DELTA_TIME] = self.spheres_model_part.ProcessInfo.GetValue(DELTA_TIME)
        self.spheres_model_part.ProcessInfo[TIME_STEPS] = self.step
        
        self.rigid_face_model_part.ProcessInfo[TIME] = self.time
        self.rigid_face_model_part.ProcessInfo[DELTA_TIME] = self.spheres_model_part.ProcessInfo.GetValue(DELTA_TIME)
        self.rigid_face_model_part.ProcessInfo[TIME_STEPS] = self.step

        # Solve
        self.solver.Initialize()
        self.solver.Solve()

        new_mesh = SMesh(name="CDF")

        # Export data back to SimPhoNy
        self.__exportKratosElements(
            self.spheres_model_part,
            new_mesh
        )

        self.__exportKratosConditions(
            self.rigid_face_model_part,
            new_mesh
        )

        self.add_mesh(SMesh(name="CDF"))

        # debug = [[n.X, n.Y, n.Z] for n in self.spheres_model_part.Nodes]

        # print(debug[0])
        # print(self.meshes["Simphony"].get_point(self.id_to_uuid_node_map[1]).coordinates)

        self.updateForwardDicc()

        self.time += 0
        self.step += 1
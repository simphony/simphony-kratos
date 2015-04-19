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

# Kratos Imports
import ProjectParameters

from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MeshingApplication import *

# TBD: Which solver?
# import vms_fractional_step_solver as solver

class CFDWrapper(KratosWrapper):

    def __init__(self):
        KratosWrapper.__init__(self)
        self.time = 0
        self.step = 0

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
            "REACTION": [
                CUBA.REACTION,
                REACTION
            ],
            "DISTANCE": [
                CUBA.DISTANCE,
                DISTANCE
            ],
            "DISPLACEMENT": [
                CUBA.DISPLACEMENT,
                DISPLACEMENT,
                DISPLACEMENT_X,
                DISPLACEMENT_Y,
                DISPLACEMENT_Z
            ],
            "VISCOSITY" : [
                CUBA.VISCOSITY,
                VISCOSITY
            ],
            "DENSITY" : [
                CUBA.DENSITY,
                DENSITY
            ],
            "BODY_FORCE" : [
                CUBA.BODY_FORCE,
                BODY_FORCE
            ],
            "FLAG_VARIABLE" : [
                CUBA.FLAG_VARIABLE,
                FLAG_VARIABLE
            ],
            "IS_STRUCTURE" : [
                CUBA.IS_STRUCTURE, 
                IS_STRUCTURE
            ],
            "IS_SLIP" : [
                CUBA.IS_SLIP, 
                IS_SLIP
            ]
        }

        # Falta:
        # - Y_WALL -- elemental? = 0
        # - POROSITY = 1

    def __addNodalVariablesToModelpart(self):
        """ Adds the Kratos CFD nodal variables

        Adds the Kratos CFD nodal variables to the particle and
        solid Kratos modelparts in order to be usable later
        while importing the mesh.

        """

        # if "REACTION" in ProjectParameters.nodal_results:
        self.fluid_model_part.AddNodalSolutionStepVariable(REACTION)
        # if "DISTANCE" in ProjectParameters.nodal_results:
        self.fluid_model_part.AddNodalSolutionStepVariable(DISTANCE)

        self.solver_module.AddVariables(self.fluid_model_part, self.SolverSettings)

    # Small kernels to get ( kratos to simp) and set ( simp to kratos)
    # entity data 

    def __getSolutionStepVariable1D(self, data, entity, variable):
        pair = self.variables_dictionary[variable]
        data.update({
            pair[0]: entity.GetSolutionStepValue(pair[1])
        })

    def __getSolutionStepVariable3D(self, data, entity, variable):
        pair = self.variables_dictionary[variable]
        data.update({
            pair[0]: [
                entity.GetSolutionStepValue(pair[2]), 
                entity.GetSolutionStepValue(pair[3]), 
                entity.GetSolutionStepValue(pair[4])
            ]
        })

    def __setSolutionStepVariable1D(self, data, entity, variable):
        pair = self.variables_dictionary[variable]
        entity.SetSolutionStepValue(
            pair[1],
            data[pair[0]]
        )

    def __setSolutionStepVariable3D(self, data, entity, variable):
        pair = self.variables_dictionary[variable]
        for i in xrange(0,3):
            entity.SetSolutionStepValue(
                pair[2+i],
                data[pair[0]][0+i]
            )

    # gets data for the nodes

    def __getNodalData(self, data, node):
        """ Extracts the node data

        Extracts the node data and puts in ina format readeable
        by the Simphony DataContainer

        """

        self.__getSolutionStepVariable1D(data,node,"PRESSURE")
        self.__getSolutionStepVariable3D(data,node,"VELOCITY")
        self.__getSolutionStepVariable1D(data,node,"REACTION")
        self.__getSolutionStepVariable1D(data,node,"DISTANCE")
        self.__getSolutionStepVariable3D(data,node,"DISPLACEMENT")
        self.__getSolutionStepVariable1D(data,node,"VISCOSITY")
        self.__getSolutionStepVariable1D(data,node,"DENSITY")
        self.__getSolutionStepVariable1D(data,node,"BODY_FORCE")
        self.__getSolutionStepVariable1D(data,node,"FLAG_VARIABLE")
        self.__getSolutionStepVariable1D(data,node,"IS_STRUCTURE")
        # self.__getSolutionStepVariable1D(data,node,"IS_SLIP")

    def __setNodalData(self, data, node):
        """ Assembles the point data

        Assembles the point data and puts in ina format readeable
        by the Kratos ModelPart

        """

        self.__setSolutionStepVariable1D(data,node,"PRESSURE")
        self.__setSolutionStepVariable3D(data,node,"VELOCITY")
        self.__setSolutionStepVariable1D(data,node,"REACTION")
        self.__setSolutionStepVariable1D(data,node,"DISTANCE")
        self.__setSolutionStepVariable3D(data,node,"DISPLACEMENT")
        self.__setSolutionStepVariable1D(data,node,"VISCOSITY")
        self.__setSolutionStepVariable1D(data,node,"DENSITY")
        self.__setSolutionStepVariable1D(data,node,"BODY_FORCE")
        self.__setSolutionStepVariable1D(data,node,"FLAG_VARIABLE")
        self.__setSolutionStepVariable1D(data,node,"IS_STRUCTURE")
        # self.__setSolutionStepVariable1D(data,node,"IS_SLIP")

    # export

    def __exportKratosElements2(self, src, dst, entitylist=None):
        """ Parses all kratos elements to simphony cells

        Iterates over all nodes in the kratos mesh (src) and
        converts them to simphony cells (dst). While doing this operation
        any point that has not currently mapped will have his uuid
        added in the 'id_map' of the weapper

        """

        for node in src.GetNodes():

            data = {}

            self.__getNodalData(data, node)

            point_uid = None
            if node.Id in self.id_to_uuid_node_map:
                point_uid = self.id_to_uuid_node_map[node.Id]

            point = SPoint(
                coordinates=(node.X, node.Y, node.Z),
                data=DataContainer(data),
                uid=point_uid
            )

            pid = dst.add_point(point)

            self.id_to_uuid_node_map[node.Id] = pid

        for element in src.GetElements():

            element_uid = None
            if element.Id in self.id_to_uuid_element_map:
                element_uid = self.id_to_uuid_element_map[element.Id]

            point_list = [
                self.id_to_uuid_node_map[point.Id]
                for point in element.GetNodes()
            ]

            cell = SCell(
                points=point_list,
                data=DataContainer(data),
                uid=element_uid
            )

            cid = dst.add_cell(cell)

            self.id_to_uuid_element_map[element.Id] = cid

        for condition in src.GetConditions():

            condition_uid = None
            if condition.Id in self.id_to_uuid_condition_map:
                condition_uid = self.id_to_uuid_condition_map[condition.Id]

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

    def __exportKratosElements(self, src, dst, entitylist=None):
        """ Parses all kratos elements to simphony cells

        Iterates over all nodes in the kratos mesh (src) and
        converts them to simphony cells (dst). While doing this operation
        any point that has not currently mapped will have his uuid
        added in the 'id_map' of the weapper

        """

        for element in src.GetElements():

            point_list = [
                self.id_to_uuid_node_map[point.Id]
                for point in element.GetNodes()
            ]

            for node in element.GetNodes():

                data = {}

                self.__getNodalData(data, node)

                point = SPoint(
                    coordinates=(node.X, node.Y, node.Z),
                    data=DataContainer(data),
                    uid=self.id_to_uuid_node_map[node.Id]
                )

                dst.add_point(point)

            cell = SCell(
                points=point_list,
                data=DataContainer(data),
                uid=self.id_to_uuid_element_map[element.Id]
            )

            dst.add_cell(cell)

    def __exportKratosConditions(self, src, dst, entitylist=None):
        pass

    # import

    def __importKratosElements(self, src, dst, entitylist=None):
        """ Parses all simphony cells to kratos elements

        Iterates over all cells in the simphony mesh (src) and
        converts them to kratos SphericContinuumParticle3D elements (dst).
        While doing this operation any point/node pair that has not
        currently mapped will have  his uuid added in the 'id_map'
        of the wrapper

        """

        for point in src.iter_points():

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

            self.__setNodalData(data, node)

        for element in src.iter_cells(entitylist):

            if element.uid not in self.uuid_to_id_element_map.keys():
                self.uuid_to_id_element_map.update(
                    {element.uid: self.free_id}
                )

                self.free_id += 1

            element_id = self.uuid_to_id_element_map[element.uid]

            dst.CreateNewElement(
                "FractionalStep3D",
                element_id,
                [self.uuid_to_id_node_map[p] for p in element.points],
                self.element_properties)

    def __importKratosConditions(self, src, dst, entitylist=None):
        for condition in src.iter_faces(entitylist):

            if condition.uid not in self.uuid_to_id_condition_map.keys():
                self.uuid_to_id_condition_map.update(
                    {condition.uid: self.free_id}
                )

                self.free_id += 1

            condition_id = self.uuid_to_id_condition_map[condition.uid]

            dst.CreateNewCondition(
                "WallCondition3D",
                condition_id,
                [self.uuid_to_id_node_map[p] for p in condition.points],
                self.element_properties)

    # FileIO

    def read_modelpart(self, filename):
        """ Reads a Kratos formated modelpart NYI

        This adds partial support for the future FileIO
        """

        new_mesh = SMesh(name="Model")
        self.fluid_model_part = ModelPart("FluidPart")

        self.__addNodalVariablesToModelpart()

        model_part_io_fluid = ModelPartIO(filename)
        model_part_io_fluid.ReadModelPart(self.fluid_model_part)

        print(self.fluid_model_part)

        # Add the problem data
        self.setMeshData(new_mesh)

        # Export data back to SimPhoNy
        self.__exportKratosElements2(
            self.fluid_model_part,
            new_mesh
        )

        self.uuid_to_id_node_map = {}
        self.uuid_to_id_element_map = {}
        self.id_to_uuid_node_map = {}
        self.id_to_uuid_element_map = {}

        return new_mesh


    def write_modelpart(self, filename):
        """ Writes a Kratos formated modelpart

        This adds partial support for the future FileIO
        """

        mesh = self.get_mesh("Model")

        f = open(filename, 'w')

        # Variable information
        f.write('Begin ModelPartData\n')
        f.write('//  VARIABLE_NAME value\n')
        f.write('End ModelPartData\n\n')

        # Properties ( out shared values )
        f.write('Begin Properties 1\n')
        f.write('PARTICLE_DENSITY {}\n'.format(mesh.data[PARTICLE_DENSITY]))
        f.write('YOUNG_MODULUS {}\n'.format(mesh.data[YOUNG_MODULUS]))
        f.write('POISSON_RATIO {}\n'.format(mesh.data[POISSON_RATIO]))
        f.write('PARTICLE_FRICTION {}\n'.format(mesh.data[PARTICLE_FRICTION]))
        f.write('PARTICLE_COHESION {}\n'.format(mesh.data[PARTICLE_COHESION]))
        f.write('LN_OF_RESTITUTION_COEFF {}\n'.format(
            mesh.data[LN_OF_RESTITUTION_COEFF]))
        f.write('PARTICLE_MATERIAL {}\n'.format(mesh.data[PARTICLE_MATERIAL]))
        f.write('ROLLING_FRICTION {}\n'.format(mesh.data[ROLLING_FRICTION]))
        f.write('DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME {}\n'.format(
            mesh.data[DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME]))
        f.write('DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME {}\n'.format(
            mesh.data[DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME]))
        f.write('End Properties\n\n')

        # Nodes
        f.write('Begin Nodes\n')
        for point in src.iter_points(entitylist=None):
            f.write('{} {} {} {}\n').format(
                self.uuid_to_id_point_map[point.id],
                point.coordinates[0],
                point.coordinates[1],
                point.coordinates[2]
            )
        f.write('End Nodes\n\n')

        f.write('Begin Elements FractionalStep2D')
        for element in src.iter_cells(entitylist=None):
            f.write('{} {} {}\n').format(
                self.uuid_to_id_point_map[elemet.id],
                self.uuid_to_id_point_map[element.points[0]],
                self.uuid_to_id_point_map[element.points[1]]
            )
        f.write('End Elements\n\n')

        f.write('Begin NodalData RADIUS')
        for point in src.iter_points(entitylist=None):
            f.write('{} {} {} {}\n').format(
                self.uuid_to_id_point_map[point.id],
                0,
                point.data[CUBA.RADIUS],
            )
        f.write('End NodalData\n\n')

    # Properties

    def setMeshData(self, mesh):
        pass

    def setElementData(self):
        pass

    def setConditionData(self):
        pass

    # Solver

    def initialize(self):
        """ Initalizes the necessary kratos commponents

            Initalizes the necessary kratos commponents to
            execute Kratos' CFD solver
        """

        self.fluid_model_part = ModelPart("")
        self.fluid_model_part.SetBufferSize(3)

        self.element_properties = Properties(0)

        self.SolverSettings = ProjectParameters.FluidSolverConfiguration
        self.solver_module = import_solver(self.SolverSettings)
        self.solver = self.solver_module.CreateSolver(
            self.fluid_model_part,
            self.SolverSettings)

        pass

    def initializeTimeStep(self):
        pass

    def run(self):
        """ Run a step of the wrapper """

        newFluidMp = ModelPart("")

        self.fluid_model_part = newFluidMp

        self.__addNodalVariablesToModelpart()

        # model_part_io_fluid = ModelPartIO("CDF-exampleFluid")
        # model_part_io_fluid.ReadModelPart(self.fluid_model_part)

        # Import the into Kratos
        self.__importKratosElements(
            self.get_mesh("Model"),
            newFluidMp
        )

        newFluidMp.Properties.append(self.element_properties)

        self.SolverSettings = ProjectParameters.FluidSolverConfiguration
        self.solver_module = import_solver(self.SolverSettings)

        self.solver_module.AddVariables(self.fluid_model_part, self.SolverSettings)

        self.updateBackwardDicc()
        self.setElementData()
        # self.setConditionData()

        self.fluid_model_part.SetBufferSize(3)

        self.solver_module.AddDofs(self.fluid_model_part, self.SolverSettings)

        # copy Y_WALL
        for node in self.fluid_model_part.Nodes:
            y = node.GetSolutionStepValue(Y_WALL, 0)
            node.SetValue(Y_WALL, y)

        self.solver = self.solver_module.CreateSolver(
            self.fluid_model_part,
            self.SolverSettings)

        self.solver.Initialize()

        substeps = 2

        Dt = ProjectParameters.Dt
        Nsteps = ProjectParameters.nsteps
        final_time = ProjectParameters.max_time
        output_time = ProjectParameters.output_time

        self.fluid_model_part.CloneTimeStep(self.time)
        self.time = self.time + Dt
        self.fluid_model_part.CloneTimeStep(self.time)
        self.time = self.time + Dt
        self.fluid_model_part.CloneTimeStep(self.time)
        self.time = self.time + Dt

        for n in xrange(0, substeps):
            self.fluid_model_part.CloneTimeStep(self.time)
            self.solver.Solve()
            self.time = self.time + Dt

        new_mesh = SMesh(name="Model")

        # Add the problem data
        self.setMeshData(new_mesh)

        # Export data back to SimPhoNy
        self.__exportKratosElements(
            self.fluid_model_part,
            new_mesh
        )

        self.add_mesh(new_mesh)
        self.updateForwardDicc()

        self.time += self.fluid_model_part.ProcessInfo.GetValue(DELTA_TIME)
        self.step += 1

    def finalizeTimeStep(self):
        pass

    def finalize(self):
        pass

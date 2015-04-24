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


class CFDWrapper(KratosWrapper):

    def __init__(self):
        KratosWrapper.__init__(self)

        self.time = 0
        self.step = 0
        self.substeps = 1

        # The dictionary defines the relation between CUBA and
        # kratos variables

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
            "VISCOSITY": [
                CUBA.VISCOSITY,
                VISCOSITY
            ],
            "DENSITY": [
                CUBA.DENSITY,
                DENSITY
            ],
            "BODY_FORCE": [
                CUBA.BODY_FORCE,
                BODY_FORCE
            ],
            "FLAG_VARIABLE": [
                CUBA.FLAG_VARIABLE,
                FLAG_VARIABLE
            ],
            "IS_STRUCTURE": [
                CUBA.IS_STRUCTURE,
                IS_STRUCTURE
            ],
            "IS_SLIP": [
                CUBA.IS_SLIP,
                IS_SLIP
            ]
        }

    def addNodalVariablesToModelpart(self, modelPart):
        """ Adds the Kratos CFD nodal variables

        Adds the Kratos CFD nodal variables to the particle and
        solid Kratos modelparts in order to be usable later
        while importing the mesh.

        """

        if "REACTION" in ProjectParameters.nodal_results:
            modelPart.AddNodalSolutionStepVariable(REACTION)
        if "DISTANCE" in ProjectParameters.nodal_results:
            modelPart.AddNodalSolutionStepVariable(DISTANCE)

        self.solver_module.AddVariables(modelPart, self.SolverSettings)

    # gets data for the nodes

    def getNodalData(self, data, node):
        """ Extracts the node data

        Extracts the node data and puts in ina format readable
        by the Simphony DataContainer

        """

        self.getSolutionStepVariable1D(data, node, "PRESSURE")
        self.getSolutionStepVariable3D(data, node, "VELOCITY")
        if "REACTION" in ProjectParameters.nodal_results:
            self.getSolutionStepVariable1D(data, node, "REACTION")
        if "DISTANCE" in ProjectParameters.nodal_results:
            self.getSolutionStepVariable1D(data, node, "DISTANCE")
        self.getSolutionStepVariable3D(data, node, "DISPLACEMENT")
        self.getSolutionStepVariable1D(data, node, "VISCOSITY")
        self.getSolutionStepVariable1D(data, node, "DENSITY")
        self.getSolutionStepVariable1D(data, node, "BODY_FORCE")
        self.getSolutionStepVariable1D(data, node, "FLAG_VARIABLE")
        self.getSolutionStepVariable1D(data, node, "IS_STRUCTURE")

        data.update({CUBA.IMPOSED_PRES: node.IsFixed(PRESSURE)})
        data.update({CUBA.IMPOSED_VEL: node.IsFixed(VELOCITY_X)})

    def setNodalData(self, data, node):
        """ Assembles the point data

        Assembles the point data and puts in ina format readable
        by the Kratos ModelPart

        """

        self.setSolutionStepVariable1D(data, node, "PRESSURE")
        self.setSolutionStepVariable3D(data, node, "VELOCITY")
        if "REACTION" in ProjectParameters.nodal_results:
            self.setSolutionStepVariable1D(data, node, "REACTION")
        if "DISTANCE" in ProjectParameters.nodal_results:
            self.setSolutionStepVariable1D(data, node, "DISTANCE")
        self.setSolutionStepVariable3D(data, node, "DISPLACEMENT")
        self.setSolutionStepVariable1D(data, node, "VISCOSITY")
        self.setSolutionStepVariable1D(data, node, "DENSITY")
        self.setSolutionStepVariable1D(data, node, "BODY_FORCE")
        self.setSolutionStepVariable1D(data, node, "FLAG_VARIABLE")
        self.setSolutionStepVariable1D(data, node, "IS_STRUCTURE")

    # export

    def exportKratosNodes(self, src, dst):
        """ Parses all kratos nodes to simphony points

        Iterates over all nodes in the kratos mesh (src) and
        converts them to simphony points (dst). While doing this operation
        any node/point that has not currently been mapped will have his uuid
        added in the 'id_map' of the wrapper

        """

        for node in src.GetNodes():

            data = {}

            self.getNodalData(data, node)

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

    def exportKratosElements(self, src, dst):
        """ Parses all kratos elements to simphony cells

        Iterates over all elements in the kratos mesh (src) and
        converts them to simphony cells (dst). While doing this operation
        any element/cell that has not currently been mapped will have his uuid
        added in the 'id_map' of the wrapper

        """

        for element in src.GetElements():

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

            cid = dst.add_cell(cell)

            self.id_to_uuid_element_map[element.Id] = cid

    def exportKratosConditions(self, src, dst):
        """ Parses all kratos conditions to simphony faces

        Iterates over all conditions in the kratos mesh (src) and
        converts them to simphony faces (dst). While doing this operation
        any condition/face that has not currently been mapped will have
        his uuid added in the 'id_map' of the wrapper

        """

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
                uid=condition_uid
            )

            fid = dst.add_face(face)

            self.id_to_uuid_condition_map[condition.Id] = fid

    def exportKratosDof(self, src, dst):
        """ Sets the Dof information for the appropiate points

        Iterates over all points in the simphony mesh and adds the Kratos
        DoF information in order to be able to recover it at the begining
        of the iteration.

        """

        for i in xrange(1, src.NumberOfMeshes()):

            mesh = src.GetMesh(i)
            prop = src.GetProperties()[i]

            nodelist = [self.id_to_uuid_node_map[n.Id] for n in mesh.Nodes]

            for p in dst.iter_points(nodelist):

                data = p.data

                if prop.GetValue(IMPOSED_PRESSURE) == 1:
                    data[CUBA.IMPOSED_PRES] = 1

                if prop.GetValue(IMPOSED_VELOCITY_X) == 1:
                    data[CUBA.IMPOSED_VEL] = 1

                p.data = data
                dst.update_point(p)

    # import

    def importKratosNodes(self, src, dst):
        """ Parses all simphony points to kratos nodes

        Iterates over all points in the simphony mesh (src) and
        converts them to kratos nodes (dst).
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

            self.setNodalData(data, node)

    def importKratosElements(self, src, dst):
        """ Parses all simphony cells to kratos elements

        Iterates over all cells in the simphony mesh (src) and
        converts them to kratos FractionalStep3D elements (dst).
        While doing this operation any cell/element pair that has not
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
                "FractionalStep3D",
                element_id,
                [self.uuid_to_id_node_map[p] for p in element.points],
                self.element_properties)

    def importKratosConditions(self, src, dst):
        """ Parses all simphony faces to kratos conditions

        Iterates over all faces in the simphony mesh (src) and
        converts them to kratos WallCondition3D conditions (dst).
        While doing this operation any face/condition pair that has not
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
                "WallCondition3D",
                condition_id,
                [self.uuid_to_id_node_map[p] for p in condition.points],
                self.element_properties)

    def importKratosDof(self, src, dst):
        """ Fix the Dof for the appropiate nodes

        Iterates over all nodes in the kratos mesh and fix its DoF
        according with the information of the SimPhoNy mesh.

        """

        for point in src.iter_points():

            node_id = self.uuid_to_id_node_map[point.uid]

            node = dst.GetNodes()[node_id]

            if point.data[CUBA.IMPOSED_VEL] == 1:
                node.Fix(VELOCITY_X)
                node.Fix(VELOCITY_Y)
                node.Fix(VELOCITY_Z)

            if point.data[CUBA.IMPOSED_PRES] == 1:
                node.Fix(PRESSURE)

    # FileIO

    def read_modelpart(self, filename):
        """ Reads a Kratos formated modelpart NYI

        This adds partial support for the future FileIO
        """

        new_mesh = SMesh(name="Model")
        self.fluid_model_part = ModelPart("FluidPart")

        self.addNodalVariablesToModelpart(self.fluid_model_part)

        model_part_io_fluid = ModelPartIO(filename)
        model_part_io_fluid.ReadModelPart(self.fluid_model_part)

        print(self.fluid_model_part)

        # Add the problem data
        self.setMeshData(new_mesh)

        # Export data back to SimPhoNy
        self.exportKratosNodes(
            self.fluid_model_part,
            new_mesh
        )
        self.exportKratosElements(
            self.fluid_model_part,
            new_mesh
        )
        self.exportKratosConditions(
            self.fluid_model_part,
            new_mesh
        )
        self.exportKratosDof(
            self.fluid_model_part,
            new_mesh
        )

        # Reset the dictionaries
        self.uuid_to_id_node_map = {}
        self.uuid_to_id_element_map = {}
        self.uuid_to_id_condition_map = {}
        self.id_to_uuid_node_map = {}
        self.id_to_uuid_element_map = {}
        self.id_to_uuid_condition_map = {}

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
        for point in src.iter_points():
            f.write('{} {} {} {}\n').format(
                self.uuid_to_id_point_map[point.id],
                point.coordinates[0],
                point.coordinates[1],
                point.coordinates[2]
            )
        f.write('End Nodes\n\n')

        f.write('Begin Elements FractionalStep2D')
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

        self.addNodalVariablesToModelpart(self.fluid_model_part)

        # Import the into Kratos
        self.importKratosNodes(
            self.get_mesh("Model"),
            self.fluid_model_part
        )
        self.importKratosElements(
            self.get_mesh("Model"),
            self.fluid_model_part
        )
        self.importKratosConditions(
            self.get_mesh("Model"),
            self.fluid_model_part
        )

        self.fluid_model_part.Properties.append(self.element_properties)

        self.SolverSettings = ProjectParameters.FluidSolverConfiguration
        self.solver_module = import_solver(self.SolverSettings)

        self.solver_module.AddVariables(
            self.fluid_model_part,
            self.SolverSettings
        )

        self.updateBackwardDicc()
        self.setElementData()
        # self.setConditionData()

        self.fluid_model_part.SetBufferSize(3)

        self.solver_module.AddDofs(self.fluid_model_part, self.SolverSettings)

        self.importKratosDof(
            self.get_mesh("Model"),
            self.fluid_model_part
        )

        # copy Y_WALL
        for node in self.fluid_model_part.Nodes:
            y = node.GetSolutionStepValue(Y_WALL, 0)
            node.SetValue(Y_WALL, y)

        self.solver = self.solver_module.CreateSolver(
            self.fluid_model_part,
            self.SolverSettings)

        self.solver.Initialize()

        Dt = self.CM[CUBA.TIME_STEP]
        # Nsteps = ProjectParameters.nsteps
        # final_time = ProjectParameters.max_time
        # output_time = ProjectParameters.output_time

        self.fluid_model_part.ProcessInfo.SetValue(DELTA_TIME, Dt)

        self.fluid_model_part.CloneTimeStep(self.time)
        self.time = self.time + Dt
        self.fluid_model_part.CloneTimeStep(self.time)
        self.time = self.time + Dt
        self.fluid_model_part.CloneTimeStep(self.time)
        self.time = self.time + Dt

        for n in xrange(0, self.substeps):
            self.fluid_model_part.CloneTimeStep(self.time)
            self.solver.Solve()
            self.time = self.time + Dt

        new_mesh = SMesh(name="Model")

        # Add the problem data
        self.setMeshData(new_mesh)

        # Export data back to SimPhoNy
        self.exportKratosNodes(
            self.fluid_model_part,
            new_mesh
        )
        self.exportKratosElements(
            self.fluid_model_part,
            new_mesh
        )
        self.exportKratosConditions(
            self.fluid_model_part,
            new_mesh
        )
        self.exportKratosDof(
            self.fluid_model_part,
            new_mesh
        )

        self.add_mesh(new_mesh)
        self.updateForwardDicc()

    def finalizeTimeStep(self):
        pass

    def finalize(self):
        pass

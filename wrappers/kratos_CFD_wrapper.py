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

        variables_dictionary = {
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
                DISPLACEMENT
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

        if "REACTION" in ProjectParameters.nodal_results:
            self.fluid_model_part.AddNodalSolutionStepVariable(REACTION)
        if "DISTANCE" in ProjectParameters.nodal_results:
            self.fluid_model_part.AddNodalSolutionStepVariable(DISTANCE)

    # Small kernels to get ( kratos to simp) and set ( simp to kratos)
    # entity data 

    def __getSolutionStepVariable1D(self, data, entity, variable):
        pair = variables_dictionary[variable]
        data.update({
            pair[0]: entity.GetSolutionStepValue(pair[1])
        })

    def __getSolutionStepVariable3D(self, data, entity, variable):
        pair = variables_dictionary[variable]
        data.update({
            pair[0]: [
                entity.GetSolutionStepValue(pair[2]), 
                entity.GetSolutionStepValue(pair[3]), 
                entity.GetSolutionStepValue(pair[4])
            ]
        })

    def __setSolutionStepVariable1D(self, data, entity, variable):
        pair = variables_dictionary[variable]
        entity.SetSolutionStepValue(
            pair[1],
            data[pair[0]]
        )

    def __setSolutionStepVariable3D(self, data, entity, variable):
        pair = variables_dictionary[variable]
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

        self.__getVariable1D(data,node,"PRESSURE")
        self.__getVariable3D(data,node,"VELOCITY")
        self.__getVariable1D(data,node,"REACTION")
        self.__getVariable1D(data,node,"DISTANCE")
        self.__getVariable3D(data,node,"DISPLACEMENT")
        self.__getVariable1D(data,node,"VISCOSITY")
        self.__getVariable1D(data,node,"DENSITY")
        self.__getVariable1D(data,node,"BODY_FORCE")
        self.__getVariable1D(data,node,"FLAG_VARIABLE")
        self.__getVariable1D(data,node,"IS_STRUCTURE")
        self.__getVariable1D(data,node,"IS_SLIP")

    def __setNodalData(self, data, node):
        """ Assembles the point data

        Assembles the point data and puts in ina format readeable
        by the Kratos ModelPart

        """

        self.__setVariable1D(data,node,"PRESSURE")
        self.__setVariable3D(data,node,"VELOCITY")
        self.__setVariable1D(data,node,"REACTION")
        self.__setVariable1D(data,node,"DISTANCE")
        self.__setVariable3D(data,node,"DISPLACEMENT")
        self.__setVariable1D(data,node,"VISCOSITY")
        self.__setVariable1D(data,node,"DENSITY")
        self.__setVariable1D(data,node,"BODY_FORCE")
        self.__setVariable1D(data,node,"FLAG_VARIABLE")
        self.__setVariable1D(data,node,"IS_STRUCTURE")
        self.__setVariable1D(data,node,"IS_SLIP")

    # export

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

        for element in src.iter_cells(entitylist):

            if element.uid not in self.uuid_to_id_element_map.keys():
                self.uuid_to_id_element_map.update(
                    {element.uid: self.free_id}
                )

                self.free_id += 1

            element_id = self.uuid_to_id_element_map[element.uid]

            for point in src.iter_points(element.points):

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

                self__.setNodalData(data, node)

            dst.CreateNewElement(
                "FractionalStep2D",
                element_id,
                [self.uuid_to_id_node_map[p] for p in element.points],
                self.element_properties)

    def __importKratosConditions(self, src, dst, entitylist=None):
        pass

    # FileIO

    def read_modelpart(self, filename):
        """ Reads a Kratos formated modelpart NYI

        This adds partial support for the future FileIO
        """

        mesh = self.get_mesh("Model")
        tmp = ModelPart("TemporalModelpart")

        f = open(filename, 'r')

        for l in f:
            if "Begin Nodes" in l:
                read_nodes = 1
                continue
        
            if "End Nodes" in l:
                read_nodes = 0
                continue

            if "Begin Elements" in l:
                read_elements = 1

            if "End Elements" in l:
                read_elements = 0

            if "Begin Conditions" in l:
                read_conditions = 1

            if "End Conditions" in l:
                read_conditions = 0
                  
            if read_nodes == 1:
                node_info = l.split(" ")

                node = tmp.CreateNewNode(
                    node_info[0],
                    node_info[1],
                    node_info[2],
                    node_info[3])

            if read_elements == 1:
                elem_info = l.split(" ")

                elem = tmp.CreateNewElement(
                    "FractionalStep3D",
                    elem_info[0],
                    [elem_info[1], elem_info[2], elem_info[3]],
                    self.element_properties
                )

            if read_conditions == 1:
                cndt_info = l.split(" ")

                cndt = tmp.CreateNewCondition(
                    "WallCondition3D",
                    cndt_info[0],
                    [cndt_info[1], cndt_info[2], cndt_info[3]],
                    self.element_properties
                )

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

        self.SolverSettings = ProjectParameters.FluidSolverConfiguration
        self.solver_module = import_solver(self.SolverSettings)

        pass

    def initializeTimeStep(self):
        pass

    def run(self):
        """ Run a step of the wrapper """

        newFluidMp = ModelPart("")

        self.fluid_model_part = newSphereMp

        newSphereMp.Properties.append(self.element_properties)
        newRigidbMp.Properties.append(self.condition_properties)

        self.__addNodalVariablesToModelpart()

        # Import the into Kratos
        self.__importKratosElements(
            self.get_mesh("Model"),
            newSphereMp
        )

        self.updateBackwardDicc()
        self.setElementData()
        # self.setConditionData()

        self.initializeTimeStep()

        # Not sure what to do here :S
        newSphereMp.ProcessInfo[TIME] = self.time
        newSphereMp.ProcessInfo[DELTA_TIME] = 0.5  # NYI
        newSphereMp.ProcessInfo[TIME_STEPS] = self.step

        newRigidbMp.ProcessInfo[TIME] = self.time
        newRigidbMp.ProcessInfo[DELTA_TIME] = 0.5  # NYI
        newRigidbMp.ProcessInfo[TIME_STEPS] = self.step

        substeps = 1

        # Solve
        for n in xrange(0, substeps):
            self.solver.Solve()

        new_mesh = SMesh(name="Model")

        # Add the problem data
        self.setMeshData(new_mesh)

        # Export data back to SimPhoNy
        self.__exportKratosElements(
            self.spheres_model_part,
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

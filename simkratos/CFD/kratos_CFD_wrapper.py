""" TODO: Description was out of date

"""

# Kratos works with Python3 by default.
from __future__ import print_function, absolute_import, division

# Simphony Imports
from simphony.core.cuba import CUBA
from simphony.core.data_container import DataContainer

from simphony.cuds.abc_mesh import ABCMesh

from simphony.cuds.mesh import Point as SPoint
from simphony.cuds.mesh import Face as SFace
from simphony.cuds.mesh import Cell as SCell

# Wrapper Imports
from simkratos.kratosWrapper import KratosWrapper
from simkratos.cuba_extension import CUBAExt

# Kratos Imports
from simkratos.CFD import ProjectParameters

from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MeshingApplication import *


class CFDWrapper(KratosWrapper):

    def __init__(self, use_internal_interface=True, **kwargs):
        super(CFDWrapper, self).__init__(use_internal_interface, **kwargs)

        self.time = 0
        self.step = 0

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
                None,
                REACTION
            ],
            "DISTANCE": [
                None,
                DISTANCE
            ],
            "DISPLACEMENT": [
                None,
                DISPLACEMENT,
                DISPLACEMENT_X,
                DISPLACEMENT_Y,
                DISPLACEMENT_Z
            ],
            "VISCOSITY": [
                None,
                VISCOSITY
            ],
            "DENSITY": [
                CUBA.DENSITY,
                DENSITY
            ],
            "BODY_FORCE": [
                None,
                BODY_FORCE
            ],
            "FLAG_VARIABLE": [
                None,
                FLAG_VARIABLE
            ],
            "IS_STRUCTURE": [
                None,
                IS_STRUCTURE
            ],
            "IS_SLIP": [
                None,
                IS_SLIP
            ]
        }

        self.initialize()

    def _load_cuds(self):
        """Load CUDS data into lammps engine."""
        cuds = self.get_cuds()
        if not cuds:
            return

        for component in cuds.iter(ABCMesh):
            self.add_dataset(component)

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

        # data.update({CUBAExtension.IMPOSED_PRES: node.IsFixed(PRESSURE)})
        # data.update({CUBAExtension.IMPOSED_VEL: node.IsFixed(VELOCITY_X)})

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

    def exportKratosNodes(self, src, dst, group):
        """ Parses all kratos nodes to simphony points

        Iterates over all nodes in the kratos mesh (src) and
        converts them to simphony points (dst). While doing this operation
        any node/point that has not currently been mapped will have his uuid
        added in the 'id_map' of the wrapper

        """

        for node in src.GetNodes(group):

            data = {}

            self.getNodalData(data, node)

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

                # iterate over the correct data
                point.data = DataContainer(data)

                dst.update([point])

    def exportKratosElements(self, src, dst, group):
        """ Parses all kratos elements to simphony cells

        Iterates over all elements in the kratos mesh (src) and
        converts them to simphony cells (dst). While doing this operation
        any element/cell that has not currently been mapped will have his uuid
        added in the 'id_map' of the wrapper

        """

        for element in src.GetElements(group):

            element_uid = None

            if element.Id not in self.id_to_uuid_element_map:

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

            else:

                # No data is stored in the element yet

                pass

    def exportKratosConditions(self, src, dst, group):
        """ Parses all kratos conditions to simphony faces

        Iterates over all conditions in the kratos mesh (src) and
        converts them to simphony faces (dst). While doing this operation
        any condition/face that has not currently been mapped will have
        his uuid added in the 'id_map' of the wrapper

        """

        for condition in src.GetConditions(group):

            condition_uid = None

            if condition.Id not in self.id_to_uuid_condition_map:

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

            else:

                # No data is stored in the condition yet

                pass

    def exportKratosDof(self, src, dst, group):
        """ Sets the Dof information for the appropiate points

        Iterates over all points in the simphony mesh and adds the Kratos
        DoF information in order to be able to recover it at the begining
        of the iteration.

        """

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

                self.setNodalData(data, node)

                self.id_to_ref_node[node_id] = node

            else:

                node = self.id_to_ref_node[self.uuid_to_id_node_map[point.uid]]
                data = point.data
                self.setNodalData(data, node)

        # If they belong to a different group, add them
        if group != 0:
            nodes = NodesArray()
            for point in src.iter(item_type=CUBA.POINT):
                nodes.append(
                    dst.Nodes[self.uuid_to_id_node_map[point.uid]]
                )
            dst.SetNodes(nodes, group)

    def importKratosElements(self, src, dst, group):
        """ Parses all simphony cells to kratos elements

        Iterates over all cells in the simphony mesh (src) and
        converts them to kratos FractionalStep3D elements (dst).
        While doing this operation any cell/element pair that has not
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
                    "FractionalStep3D",
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
        converts them to kratos WallCondition3D conditions (dst).
        While doing this operation any face/condition pair that has not
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
                    "WallCondition3D",
                    condition_id,
                    [self.uuid_to_id_node_map[p] for p in condition.points],
                    self.element_properties)

        # If they belong to a different group, add them
        if group != 0:
            conditions = ConditionsArray()
            for cnd in src.iter(item_type=CUBA.FACE):
                conditions.append(
                    dst.Conditions[self.uuid_to_id_condition_map[cnd.uid]]
                )
            dst.SetConditions(conditions, group)

    def importKratosDof(self, src, dst, group):

        bc = self.get_cuds().get(src.data[CUBA.CONDITION])

        mesh_prop = Properties(group)
        mesh_prop.SetValue(IS_SLIP, 0)

        print("Applying DOF to mesh", group, bc, src.data[CUBA.CONDITION])

        if CUBA.PRESSURE not in bc.data:
            mesh_prop.SetValue(IMPOSED_PRESSURE, 0)
        else:
            print(bc.data[CUBA.PRESSURE])
            mesh_prop.SetValue(IMPOSED_PRESSURE, 1)
            mesh_prop.SetValue(PRESSURE, bc.data[CUBA.PRESSURE])

            for node in self.fluid_model_part.GetNodes(group):
                node.Fix(PRESSURE)
                node.SetValue(PRESSURE, bc.data[CUBA.PRESSURE])

        if CUBA.VELOCITY not in bc.data:
            mesh_prop.SetValue(IMPOSED_VELOCITY_X, 0)
            mesh_prop.SetValue(IMPOSED_VELOCITY_Y, 0)
            mesh_prop.SetValue(IMPOSED_VELOCITY_Z, 0)
        else:
            print(bc.data[CUBA.VELOCITY])
            imposedVel = bc.data[CUBA.VELOCITY]
            if imposedVel[0] is not None:
                mesh_prop.SetValue(IMPOSED_VELOCITY_X, 1)
                mesh_prop.SetValue(
                    IMPOSED_VELOCITY_X_VALUE,
                    bc.data[CUBA.VELOCITY][0]
                )
            if imposedVel[1] is not None:
                mesh_prop.SetValue(IMPOSED_VELOCITY_Y, 1)
                mesh_prop.SetValue(
                    IMPOSED_VELOCITY_Y_VALUE,
                    bc.data[CUBA.VELOCITY][1]
                )
            if imposedVel[2] is not None:
                mesh_prop.SetValue(IMPOSED_VELOCITY_Z, 1)
                mesh_prop.SetValue(
                    IMPOSED_VELOCITY_Z_VALUE,
                    bc.data[CUBA.VELOCITY][2]
                )

            for node in self.fluid_model_part.GetNodes(group):
                if imposedVel[0] is not None:
                    node.Fix(VELOCITY_X)
                    node.SetValue(VELOCITY_X, bc.data[CUBA.VELOCITY][0])
                if imposedVel[1] is not None:
                    node.Fix(VELOCITY_Y)
                    node.SetValue(VELOCITY_Y, bc.data[CUBA.VELOCITY][1])
                if imposedVel[2] is not None:
                    node.Fix(VELOCITY_Z)
                    node.SetValue(VELOCITY_Z, bc.data[CUBA.VELOCITY][2])

        return mesh_prop

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

        self.addNodalVariablesToModelpart(self.fluid_model_part)

        self.SolverSettings = ProjectParameters.FluidSolverConfiguration
        self.solver_module = import_solver(self.SolverSettings)

        self.solver_module.AddVariables(
            self.fluid_model_part,
            self.SolverSettings
        )

        self.element_properties = Properties(0)

    def run(self):
        """ Run a step of the wrapper """

        fluid_meshes = self.meshes

        cuds = self.get_cuds()

        self.fluid_model_part.GetMesh(len(fluid_meshes))

        properties = PropertiesArray()
        meshNumber = 1
        meshDict = {}

        for mesh in cuds.iter(ABCMesh):

            group = meshNumber

            self.importKratosNodes(mesh, self.fluid_model_part, group)
            self.importKratosElements(mesh, self.fluid_model_part, group)
            self.importKratosConditions(mesh, self.fluid_model_part, group)

            mesh_prop = self.importKratosDof(
                mesh, self.fluid_model_part, group
            )

            meshDict[mesh.name] = meshNumber

            properties.append(mesh_prop)
            meshNumber += 1

        self.fluid_model_part.SetProperties(properties)
        print(self.fluid_model_part)

        self.updateBackwardDicc()
        self.setElementData()
        self.setConditionData()

        self.fluid_model_part.SetBufferSize(3)

        self.solver_module.AddDofs(self.fluid_model_part, self.SolverSettings)

        for node in self.fluid_model_part.Nodes:
            y = node.GetSolutionStepValue(Y_WALL, 0)
            node.SetValue(Y_WALL, y)

        self.solver = self.solver_module.CreateSolver(
            self.fluid_model_part,
            self.SolverSettings)

        self.solver.Initialize()

        Dt = cuds.get('cfd_integration_time').step

        self.fluid_model_part.ProcessInfo.SetValue(DELTA_TIME, Dt)

        # Start the simulation itself
        self.time = cuds.get('cfd_integration_time').time
        self.final = cuds.get('cfd_integration_time').final

        # Init the temporal db without starting the simulation since we
        # cannot make sure this is the first execution of kratos or not.
        # NOTE: Temporal db from previous kratos steps is lost.
        for i in xrange(0, 3):
            self.fluid_model_part.CloneTimeStep(self.time)
            self.time = self.time + Dt

        while self.time < self.final:
            self.fluid_model_part.CloneTimeStep(self.time)
            self.solver.Solve()
            self.time = self.time + Dt

        cuds.get('cfd_integration_time').time = self.time
        cuds.get('cfd_integration_time').final = self.final

        # Resotre the information to SimPhoNy
        for mesh in cuds.iter(ABCMesh):

            group = meshDict[mesh.name]

            self.exportKratosNodes(self.fluid_model_part, mesh, group)
            self.exportKratosElements(self.fluid_model_part, mesh, group)
            self.exportKratosConditions(self.fluid_model_part, mesh, group)

        self.updateForwardDicc()

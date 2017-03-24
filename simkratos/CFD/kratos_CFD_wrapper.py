""" TODO: Description was out of date

"""

# Kratos works with Python3 by default.
from __future__ import print_function, absolute_import, division

# Simphony Imports
from simphony.core.cuba import CUBA
from simphony.core.data_container import DataContainer


# Wrapper Imports
from simkratos.kratosWrapper import KratosWrapper

import KratosMultiphysics as KRTS                # noqa                         # noqa

# Kratos Imports
from simkratos.CFD import ProjectParameters


class CFDWrapper(KratosWrapper):

    def __init__(self, use_internal_interface=True, **kwargs):

        super(CFDWrapper, self).__init__(use_internal_interface, **kwargs)

        self.time = 0
        self.step = 0
        self.element_type = "FractionalStep3D"
        self.condition_type = "WallCondition3D"

        # The dictionary defines the relation between CUBA and
        # kratos variables

        self.variables_dictionary = {
            "PRESSURE": [
                CUBA.PRESSURE,
                KRTS.PRESSURE
            ],
            "VELOCITY": [
                CUBA.VELOCITY,
                KRTS.VELOCITY,
                KRTS.VELOCITY_X,
                KRTS.VELOCITY_Y,
                KRTS.VELOCITY_Z
            ],
            "REACTION": [
                None,
                KRTS.REACTION
            ],
            "DISTANCE": [
                None,
                KRTS.DISTANCE
            ],
            "DISPLACEMENT": [
                None,
                KRTS.DISPLACEMENT,
                KRTS.DISPLACEMENT_X,
                KRTS.DISPLACEMENT_Y,
                KRTS.DISPLACEMENT_Z
            ],
            "VISCOSITY": [
                None,
                KRTS.VISCOSITY
            ],
            "DENSITY": [
                CUBA.DENSITY,
                KRTS.DENSITY
            ],
            "BODY_FORCE": [
                None,
                KRTS.BODY_FORCE
            ],
            "FLAG_VARIABLE": [
                None,
                KRTS.FLAG_VARIABLE
            ],
            "IS_STRUCTURE": [
                None,
                KRTS.IS_STRUCTURE
            ],
            "IS_SLIP": [
                None,
                KRTS.IS_SLIP
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

    def addNodalVariablesToModelpart(self, modelPart):
        """ Adds the Kratos CFD nodal variables

        Adds the Kratos CFD nodal variables to the particle and
        solid Kratos modelparts in order to be usable later
        while importing the mesh.

        """

        if "REACTION" in ProjectParameters.nodal_results:
            modelPart.AddNodalSolutionStepVariable(KRTS.REACTION)
        if "DISTANCE" in ProjectParameters.nodal_results:
            modelPart.AddNodalSolutionStepVariable(KRTS.DISTANCE)

    # gets data for the nodes

    def getNodalData(self, data, node, model):
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

    def setNodalData(self, data, node, model):
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

    def exportKratosDof(self, src, dst, group):
        """ Sets the Dof information for the appropiate points

        Iterates over all points in the simphony mesh and adds the Kratos
        DoF information in order to be able to recover it at the begining
        of the iteration.

        """

        pass

    def importKratosDof(self, src, dst, group):

        print("=================================================")
        print("=================================================")

        bc = self.get_cuds().get_by_name(src.data[CUBA.CONDITION])

        mesh_prop = KRTS.Properties(group)
        mesh_prop.SetValue(KRTS.IS_SLIP, 0)

        if CUBA.PRESSURE not in bc.data:
            mesh_prop.SetValue(KRTS.IMPOSED_PRESSURE, 0)
        else:
            mesh_prop.SetValue(KRTS.IMPOSED_PRESSURE, 1)
            mesh_prop.SetValue(KRTS.PRESSURE, bc.data[CUBA.PRESSURE])

            for node in self.fluid_model_part.GetNodes(group):
                node.Fix(KRTS.PRESSURE)
                node.SetValue(KRTS.PRESSURE, bc.data[CUBA.PRESSURE])

        if CUBA.VELOCITY not in bc.data:
            mesh_prop.SetValue(KRTS.IMPOSED_VELOCITY_X, 0)
            mesh_prop.SetValue(KRTS.IMPOSED_VELOCITY_Y, 0)
            mesh_prop.SetValue(KRTS.IMPOSED_VELOCITY_Z, 0)
        else:
            imposedVel = bc.data[CUBA.VELOCITY]
            if imposedVel[0] is not None:
                mesh_prop.SetValue(KRTS.IMPOSED_VELOCITY_X, 1)
                mesh_prop.SetValue(
                    KRTS.IMPOSED_VELOCITY_X_VALUE,
                    bc.data[CUBA.VELOCITY][0]
                )
            if imposedVel[1] is not None:
                mesh_prop.SetValue(KRTS.IMPOSED_VELOCITY_Y, 1)
                mesh_prop.SetValue(
                    KRTS.IMPOSED_VELOCITY_Y_VALUE,
                    bc.data[CUBA.VELOCITY][1]
                )
            if imposedVel[2] is not None:
                mesh_prop.SetValue(KRTS.IMPOSED_VELOCITY_Z, 1)
                mesh_prop.SetValue(
                    KRTS.IMPOSED_VELOCITY_Z_VALUE,
                    bc.data[CUBA.VELOCITY][2]
                )

            for node in self.fluid_model_part.GetNodes(group):
                if imposedVel[0] is not None:
                    node.Fix(KRTS.VELOCITY_X)
                    node.SetValue(KRTS.VELOCITY_X, bc.data[CUBA.VELOCITY][0])
                if imposedVel[1] is not None:
                    node.Fix(KRTS.VELOCITY_Y)
                    node.SetValue(KRTS.VELOCITY_Y, bc.data[CUBA.VELOCITY][1])
                if imposedVel[2] is not None:
                    node.Fix(KRTS.VELOCITY_Z)
                    node.SetValue(KRTS.VELOCITY_Z, bc.data[CUBA.VELOCITY][2])

        print("=================================================")
        print("=================================================")

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

        self.internal_step = 0
        self.bufferSize = 3
        self.fluid_model_part = KRTS.ModelPart("")
        self.fluid_model_part.SetBufferSize(self.bufferSize)

        self.addNodalVariablesToModelpart(self.fluid_model_part)

        self.SolverSettings = ProjectParameters.FluidSolverConfiguration
        self.solver_module = KRTS.import_solver(self.SolverSettings)

        self.solver_module.AddVariables(
            self.fluid_model_part,
            self.SolverSettings
        )

        self.element_properties = KRTS.Properties(0)
        self.condition_properties = KRTS.Properties(0)

    def run(self):
        """ Run a step of the wrapper """

        fluid_meshes = self.meshes

        cuds = self.get_cuds()

        self.fluid_model_part.GetMesh(len(fluid_meshes))

        properties = KRTS.PropertiesArray()
        meshNumber = 1
        meshDict = {}

        # for mesh in cuds.iter(item_type=CUBA.MESH):

        # Get the CFD pe
        if cuds.count_of(item_type=CUBA.CFD) != 1:
            raise "KratosCFD needs exactly one CFD pe."

        for cfd_pe in cuds.iter(item_type=CUBA.CFD):
            if len(cfd_pe.data[CUBA.DATA_SET]) < 1:
                raise "CFD PE does not have any associated dataset"

            for name in cfd_pe.data[CUBA.DATA_SET]:

                mesh = cuds.get_by_name(name)
                group = meshNumber

                self.importKratosNodes(
                    mesh, self.fluid_model_part,
                    group
                )
                self.importKratosElements(
                    mesh, self.fluid_model_part,
                    group, self.element_type
                )
                self.importKratosConditions(
                    mesh, self.fluid_model_part,
                    group, self.condition_type
                )

                mesh_prop = self.importKratosDof(
                    mesh, self.fluid_model_part, group
                )

                meshDict[mesh.name] = meshNumber

                properties.append(mesh_prop)
                meshNumber += 1

        self.fluid_model_part.SetProperties(properties)

        self.updateBackwardDicc()
        self.setElementData()
        self.setConditionData()

        self.fluid_model_part.SetBufferSize(self.bufferSize)

        self.solver_module.AddDofs(self.fluid_model_part, self.SolverSettings)

        for node in self.fluid_model_part.Nodes:
            y = node.GetSolutionStepValue(KRTS.Y_WALL, 0)
            node.SetValue(KRTS.Y_WALL, y)

        self.solver = self.solver_module.CreateSolver(
            self.fluid_model_part,
            self.SolverSettings)

        self.solver.Initialize()

        if cuds.count_of(item_type=CUBA.INTEGRATION_TIME) < 0:
            raise Exception("Error: No integran time")

        if cuds.count_of(item_type=CUBA.INTEGRATION_TIME) > 1:
            raise Exception("Error: More than one integration time")

        iTime = [it for it in cuds.iter(item_type=CUBA.INTEGRATION_TIME)][0]
        #
        # if(len(iTime) != 1):
        #     raise "Error: Cuds has more than one or zero integration times."
        #
        # iTime = iTime[0]

        Dt = iTime.step

        self.fluid_model_part.ProcessInfo.SetValue(KRTS.DELTA_TIME, Dt)

        # Start the simulation itself
        self.time = iTime.time
        self.final = iTime.final

        print(self.time, self.final)

        while self.time < self.final:
            self.fluid_model_part.CloneTimeStep(self.time)

            if self.internal_step >= self.bufferSize:
                self.solver.Solve()

            self.internal_step += 1
            self.time = self.time + Dt

        iTime.time = self.time
        iTime.final = self.final

        # Resotre the information to SimPhoNy
        for cfd_pe in cuds.iter(item_type=CUBA.CFD):
            if len(cfd_pe.data[CUBA.DATA_SET]) < 1:
                raise "CFD PE does not have any associated dataset"

            for name in cfd_pe.data[CUBA.DATA_SET]:

                mesh = cuds.get_by_name(name)
                group = meshDict[mesh.name]

                self.exportKratosNodes(self.fluid_model_part, mesh, group)
                self.exportKratosElements(self.fluid_model_part, mesh, group)
                self.exportKratosConditions(self.fluid_model_part, mesh, group)

        self.updateForwardDicc()

        for p in mesh.iter(item_type=CUBA.POINT):
            print(p.data[CUBA.VELOCITY])

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
import KratosMultiphysics as KRTS
import KratosMultiphysics.DEMApplication as KRTSDEM
import KratosMultiphysics.SwimmingDEMApplication as KRTSSWDEM

import sphere_strategy as SolverStrategy
import DEM_procedures

import simkratos.PROJECT.DEM_explicit_solver_var as DEM_parameters
import simkratos.PROJECT.ProjectParameters as pp

import swimming_DEM_algorithm
import math
import swimming_DEM_procedures as swim_proc
import CFD_DEM_coupling
import embedded
import swimming_DEM_algorithm


class PROJECTWrapper(KratosWrapper):

    def __init__(self, use_internal_interface=True, **kwargs):
        super(PROJECTWrapper, self).__init__(use_internal_interface, **kwargs)

        self.time = 0
        self.step = 0
        self.substeps = 0
        self.element_type = "VMS3D4N"
        self.particle_type = "SphericParticle3D"
        self.condition_type = "RigidFace3D3N"

        # The dictionary defines the relation between CUBA and
        # kratos variables

        self.variables_dictionary = {
            "RADIUS": [
                CUBA.RADIUS,
                KRTS.RADIUS
            ],
            "NODAL_MASS": [
                None,
                KRTS.NODAL_MASS
            ],
            "NODAL_AREA": [
                None,
                KRTS.NODAL_AREA
            ],
            "VELOCITY": [
                CUBA.VELOCITY,
                KRTS.VELOCITY,
                KRTS.VELOCITY_X,
                KRTS.VELOCITY_Y,
                KRTS.VELOCITY_Z
            ],
            "DISPLACEMENT": [
                None,
                KRTS.DISPLACEMENT,
                KRTS.DISPLACEMENT_X,
                KRTS.DISPLACEMENT_Y,
                KRTS.DISPLACEMENT_Z
            ],
            "FLUID_VEL_PROJECTED": [
                CUBA.RELATIVE_VELOCITY,
                KRTS.FLUID_VEL_PROJECTED,
                KRTS.FLUID_VEL_PROJECTED_X,
                KRTS.FLUID_VEL_PROJECTED_Y,
                KRTS.FLUID_VEL_PROJECTED_Z
            ]
        }

        # self.algorithm = swimming_DEM_algorithm
        # self.varying_parameters = dict()
        # self.pp = pp
        # self.pp.CFD_DEM = DEM_parameters
        # self.alg = self.algorithm.Algorithm(self.pp)
        # self.alg.SetCustomBetaParamters(self.varying_parameters)

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
            self.getSolutionStepVariable3D(data, node, "FLUID_VEL_PROJECTED")
        else:
            self.getSolutionStepVariable1D(data, node, "NODAL_AREA")

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
            self.setSolutionStepVariable3D(data, node, "FLUID_VEL_PROJECTED")
        else:
            self.setSolutionStepVariable1D(data, node, "NODAL_AREA")

    def _setMeshData(self):
        " This probably needs to be done throug configuration"

        cLawString = "DEM_KDEM"
        dLawString = "DEM_D_Hertz_viscous_Coulomb"

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
            KRTSDEM.PARTICLE_DENSITY,
            self.SP[CUBA.DENSITY]
        )
        self.element_properties.SetValue(
            KRTS.YOUNG_MODULUS,
            self.SP[CUBA.YOUNG_MODULUS]
        )
        self.element_properties.SetValue(
            KRTS.POISSON_RATIO,
            self.SP[CUBA.POISSON_RATIO]
        )
        self.element_properties.SetValue(
            KRTSDEM.PARTICLE_FRICTION,
            self.PARTICLE_FRICTION
        )
        self.element_properties.SetValue(
            KRTSDEM.PARTICLE_COHESION,
            self.PARTICLE_COHESION
        )
        self.element_properties.SetValue(
            KRTS.PARTICLE_MATERIAL,
            self.PARTICLE_MATERIAL
        )
        self.element_properties.SetValue(
            KRTSDEM.ROLLING_FRICTION,
            self.SP[CUBA.ROLLING_FRICTION]
        )
        self.element_properties.SetValue(
            KRTSDEM.DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME,
            self.DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME
        )
        self.element_properties.SetValue(
            KRTSDEM.DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME,
            self.DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME
        )

    def setConditionData(self):
        cmesh = self.get_dataset("Conditions")

        self.condition_properties.SetValue(
            KRTS.WALL_FRICTION,
            cmesh.data[CUBA.WALL_FRICTION]
        )

    def exportKratosDof(self, src, dst, group):
        """ Sets the Dof information for the appropiate points

        Iterates over all points in the simphony mesh and adds the Kratos
        DoF information in order to be able to recover it at the begining
        of the iteration.

        """

        pass

    def importKratosDof(self, src, dst, group):

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

        return mesh_prop

    def initialize(self):
        # Prepare properties
        self.element_properties = KRTS.Properties(0)
        self.condition_properties = KRTS.Properties(1)

        self.fluid_model_part = KRTS.ModelPart("Fluid")
        self.spheres_model_part = KRTS.ModelPart("Particles")
        self.rigid_face_model_part = KRTS.ModelPart("P_Conditions")
        self.mixed_model_part = KRTS.ModelPart("Mixed")

        # self.alg.AddExtraVariables()

        # self.bounding_box_time_limits = []
        # if self.pp.CFD_DEM.BoundingBoxOption == "ON":
        #     self.alg.procedures.SetBoundingBox(self.alg.spheres_model_part, self.alg.cluster_model_part, self.alg.rigid_face_model_part, self.alg.creator_destructor)
        #     self.bounding_box_time_limits = [self.alg.solver.bounding_box_start_time, self.alg.solver.bounding_box_stop_time]

        # Creating the fluid solver
        # self.SolverSettings = self.pp.FluidSolverConfiguration

        # Reading the timestep
        # self.Dt_DEM = self.pp.CFD_DEM.MaxTimeStep

        # reading the fluid part
        # self.alg.Initialize()
        # self.alg.SetFluidBufferSizeAndAddAdditionalDofs()

    def run(self):
        """ Run a step of the wrapper """

        fluid_particles = self.pcs
        solid_meshes = self.meshes

        cuds = self.get_cuds()

        self.fluid_model_part = KRTS.ModelPart("Fluid")
        self.spheres_model_part = KRTS.ModelPart("Particles")
        self.rigid_face_model_part = KRTS.ModelPart("Bc")

        self.fluid_model_part.GetMesh(len(solid_meshes))
        self.spheres_model_part.GetMesh(len(fluid_particles))
        self.rigid_face_model_part.GetMesh(len(solid_meshes))

        self.spheres_model_part.AddNodalSolutionStepVariable(KRTS.RADIUS)
        self.spheres_model_part.AddNodalSolutionStepVariable(KRTS.VELOCITY)
        self.spheres_model_part.AddNodalSolutionStepVariable(KRTS.FLUID_VEL_PROJECTED)

        meshNumber = 1
        meshDict = {}

        # Get the CFD pe
        if cuds.count_of(item_type=CUBA.CFD) != 1:
            raise "KratosCFD needs exactly one CFD pe."

        for cfd_pe in cuds.iter(item_type=CUBA.CFD):
            if len(cfd_pe.data[CUBA.DATA_SET]) < 1:
                raise Exception("CFD PE does not have any associated dataset")

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

                # mesh_prop = self.importKratosDof(
                #     mesh, self.fluid_model_part, group
                # )

                meshDict[mesh.name] = meshNumber

                # properties.append(mesh_prop)
                meshNumber += 1

        # Get the DEM pe
        if cuds.count_of(item_type=CUBA.GRANULAR_DYNAMICS) != 1:
            raise Exception("KratosDEM needs exactly one GRANULAR_DYNAMICS pe.")

        # Reset the mesh number while importing the DEM Modelpart
        meshNumber = 1

        for gd_pe in cuds.iter(item_type=CUBA.GRANULAR_DYNAMICS):
            if len(gd_pe.data[CUBA.DATA_SET]) < 1:
                raise Exception("GD PE does not have any associated dataset")

            for name in gd_pe.data[CUBA.DATA_SET]:

                particles = cuds.get_by_name(name)
                group = meshNumber

                self.importKratosParticles(
                    particles, self.spheres_model_part,
                    group, self.particle_type
                )

                meshDict[particles.name] = meshNumber
                meshNumber += 1

        self.updateBackwardDicc()

        print("Before proj")

        print(self.fluid_model_part)

        # self.projector = KRTSSWDEM.BinBasedDEMFluidCoupledMapping3D(0.6, 1, 0, 1)
        #
        # self.bin_of_objects_fluid = KRTS.BinBasedFastPointLocator3D(self.fluid_model_part)
        # self.bin_of_objects_fluid.UpdateSearchDatabase()
        #
        # self.projector.InterpolateFromNewestFluidMesh(self.fluid_model_part, self.spheres_model_part, self.bin_of_objects_fluid)

        print("After proj")

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

        for gd_pe in cuds.iter(item_type=CUBA.GRANULAR_DYNAMICS):
            if len(gd_pe.data[CUBA.DATA_SET]) < 1:
                raise Exception("GD PE does not have any associated dataset")

            for name in gd_pe.data[CUBA.DATA_SET]:

                particles = cuds.get_by_name(name)
                group = meshDict[particles.name]

                self.exportKratosParticles(
                    self.spheres_model_part,
                    particles,
                    group
                )

        self.updateForwardDicc()

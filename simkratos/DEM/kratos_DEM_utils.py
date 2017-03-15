from simphony.core.cuba import CUBA
from simphony.cuds.meta import api

from simphony.cuds.mesh import Mesh as SMesh
from simphony.cuds.particles import Particles as SParticles

from simkratos.DEM.kratos_DEM_wrapper import DEMWrapper

# Kratos Imports
import KratosMultiphysics as KRTS
import KratosMultiphysics.DEMApplication as KRTSDEM


class DEM_Utils(DEMWrapper):

    def __init__(self, use_internal_interface=True, **kwargs):
        super(DEMWrapper, self).__init__(use_internal_interface, **kwargs)

        self.supportedMaterialProp = {
            CUBA.DENSITY: KRTSDEM.PARTICLE_DENSITY,
            CUBA.YOUNG_MODULUS: KRTS.YOUNG_MODULUS,
            CUBA.POISSON_RATIO: KRTS.POISSON_RATIO,
            CUBA.FRICTION_COEFFICIENT: KRTSDEM.PARTICLE_FRICTION,
            CUBA.ROLLING_FRICTION: KRTSDEM.ROLLING_FRICTION
        }

    def _getNodalData(self, data, node):
        """ Extracts the node data

        Extracts the node data and puts in ina format readable
        by the Simphony DataContainer

        """

        self._getSolutionStepVariable1D(data, node, "RADIUS")
        self._getSolutionStepVariable1D(data, node, "DENSITY")
        self._getSolutionStepVariable3D(data, node, "VELOCITY")

    def _convertBc(self, properties, mesh_name):
        condition = api.Condition(name='condition_' + mesh_name)
        conditionData = condition.data

        velocityCondition = [None, None, None]

        if properties.GetValue(KRTS.IMPOSED_PRESSURE) == 1:
            conditionData[CUBA.PRESSURE] = properties.GetValue(
                KRTS.PRESSURE
            )

        if properties.GetValue(KRTS.IMPOSED_VELOCITY_X) == 1:
            velocityCondition[0] = properties.GetValue(
                KRTS.IMPOSED_VELOCITY_X_VALUE
            )
        if properties.GetValue(KRTS.IMPOSED_VELOCITY_Y) == 1:
            velocityCondition[1] = properties.GetValue(
                KRTS.IMPOSED_VELOCITY_Y_VALUE
            )
        if properties.GetValue(KRTS.IMPOSED_VELOCITY_Z) == 1:
            velocityCondition[2] = properties.GetValue(
                KRTS.IMPOSED_VELOCITY_Z_VALUE
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

        model_part = KRTS.ModelPart("FluidPart")

        model_part.AddNodalSolutionStepVariable(KRTS.RADIUS)
        model_part.AddNodalSolutionStepVariable(KRTS.DENSITY)
        model_part.AddNodalSolutionStepVariable(KRTS.VELOCITY)

        model_part_io_fluid = KRTS.ModelPartIO(filename)
        model_part_io_fluid.ReadModelPart(model_part)

        smp_meshes = []
        smp_conditions = []
        smp_materials = []

        for i in xrange(1, model_part.NumberOfMeshes()):

            mesh_name = 'solid_' + str(i)

            mesh = SMesh(name=mesh_name)

            # Export mesh data to Simphony
            self.exportKratosNodes(model_part, mesh, i)
            self.exportKratosElements(model_part, mesh, i)
            self.exportKratosConditions(model_part, mesh, i)

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

        model_part = KRTS.ModelPart("FluidPart")

        model_part.AddNodalSolutionStepVariable(KRTS.RADIUS)
        model_part.AddNodalSolutionStepVariable(KRTS.DENSITY)
        model_part.AddNodalSolutionStepVariable(KRTS.VELOCITY)

        model_part_io_fluid = KRTS.ModelPartIO(filename)
        model_part_io_fluid.ReadModelPart(model_part)

        smp_particles = []
        smp_conditions = []
        smp_materials = []

        for i in xrange(1, model_part.NumberOfMeshes()):

            particles_name = 'particles_' + str(i)

            particles = SParticles(name=particles_name)

            # Export particle data to Simphony
            self.exportKratosParticles(model_part, particles, i)

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

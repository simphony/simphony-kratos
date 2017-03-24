from simphony.core.cuba import CUBA
from simphony.cuds.meta import api

from simphony.cuds.mesh import Mesh as SMesh

from simkratos.CFD.kratos_CFD_wrapper import CFDWrapper

# Kratos Imports
import KratosMultiphysics as KRTS


class CFD_Utils(CFDWrapper):

    def __init__(self, use_internal_interface=True, **kwargs):
        super(CFD_Utils, self).__init__(use_internal_interface, **kwargs)

        self.supportedMaterialProp = {
        }

    def _getNodalData(self, data, node):
        """ Extracts the node data

        Extracts the node data and puts in ina format readable
        by the Simphony DataContainer

        """

        self._getSolutionStepVariable1D(data, node, "PRESSURE")
        self._getSolutionStepVariable3D(data, node, "VELOCITY")
        self._getSolutionStepVariable1D(data, node, "VISCOSITY")
        self._getSolutionStepVariable1D(data, node, "DENSITY")

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

        for key, val in self.supportedMaterialProp:
            if properties.GetValue(val):
                materialData[key] = properties.GetValue(val)

        material.data = materialData
        return material

    def read_modelpart(self, filename):
        """ Reads a Kratos formated modelpart for KratosCFD wrapper

        This functions translates each of the modelpart meshes to
        one simphony mesh.

        Properties in the modelpart related to boundary conditions are
        stored in simphony conditions and referenced in CUBA.CONDITION
        inside the mesh.

        Properties in the modelpart not relared to boundary conditions
        are stored in simphony as materials and referenced in
        CUBA.MATERIAL inside the mesh.

        """

        model_part = KRTS.ModelPart("FluidPart")

        model_part.AddNodalSolutionStepVariable(KRTS.PRESSURE)
        model_part.AddNodalSolutionStepVariable(KRTS.VELOCITY)
        model_part.AddNodalSolutionStepVariable(KRTS.VISCOSITY)
        model_part.AddNodalSolutionStepVariable(KRTS.DENSITY)

        model_part_io_fluid = KRTS.ModelPartIO(filename)
        model_part_io_fluid.ReadModelPart(model_part)

        smp_meshes = []
        smp_conditions = []
        smp_materials = []

        cfd_pe = api.Cfd()
        cfd_pe.data[CUBA.DATA_SET] = []

        for i in xrange(0, model_part.NumberOfMeshes()):

            mesh_name = 'fluid_' + str(i)

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

            # Add to the pe?
            cfd_pe.data[CUBA.DATA_SET].append(mesh.name)

        return {
            'datasets': smp_meshes,
            'conditions': smp_conditions,
            'materials': smp_materials,
            'pe': cfd_pe,
        }

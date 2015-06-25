from simphony.core.cuba import CUBA
from simphony.core.data_container import DataContainer

from simphony.cuds.mesh import Point as SPoint
from simphony.cuds.mesh import Mesh as SMesh
from simphony.cuds.mesh import Face as SFace
from simphony.cuds.mesh import Cell as SCell

from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MeshingApplication import *

from wrappers.tests.cfd import ProjectParameters

id_to_uuid_node_map = {}
uuid_to_id_node_map = {}
id_to_uuid_element_map = {}
uuid_to_id_element_map = {}
id_to_uuid_condition_map = {}
uuid_to_id_condition_map = {}

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
            "VISCOSITY": [
                None,
                VISCOSITY
            ],
            "DENSITY": [
                CUBA.DENSITY,
                DENSITY
            ]
        }

def getSolutionStepVariable1D(data, entity, variable):
    pair = variables_dictionary[variable]
    if(pair[0] is not None):
        data.update({
            pair[0]: entity.GetSolutionStepValue(pair[1])
        })

def getSolutionStepVariable3D(data, entity, variable):
    pair = variables_dictionary[variable]
    if(pair[0] is not None):
        data.update({
            pair[0]: [
                entity.GetSolutionStepValue(pair[2]),
                entity.GetSolutionStepValue(pair[3]),
                entity.GetSolutionStepValue(pair[4])
            ]
        })

def getNodalData(data, node):
    """ Extracts the node data

    Extracts the node data and puts in ina format readable
    by the Simphony DataContainer

    """

    getSolutionStepVariable1D(data, node, "PRESSURE")
    getSolutionStepVariable3D(data, node, "VELOCITY")
    getSolutionStepVariable1D(data, node, "VISCOSITY")
    getSolutionStepVariable1D(data, node, "DENSITY")

def read_modelpart(filename):
    """ Reads a Kratos formated modelpart

    This adds partial support for the future FileIO
    """

    model_part = ModelPart("FluidPart")

    SolverSettings = ProjectParameters.FluidSolverConfiguration
    solver_module = import_solver(SolverSettings)

    model_part.AddNodalSolutionStepVariable(PRESSURE)
    model_part.AddNodalSolutionStepVariable(VELOCITY)
    model_part.AddNodalSolutionStepVariable(VISCOSITY)
    model_part.AddNodalSolutionStepVariable(DENSITY)

    model_part_io_fluid = ModelPartIO(filename)
    model_part_io_fluid.ReadModelPart(model_part)

    smp_meshes = []
    smp_bcs = []

    print(model_part)

    for i in xrange(0,model_part.NumberOfMeshes()):

        mesh_name = 'fluid_'+str(i)

        smp_mesh = SMesh(name = mesh_name)

        # Export data back to SimPhoNy
        exportKratosNodes(
            model_part,
            smp_mesh,
            i
        )
        exportKratosElements(
            model_part,
            smp_mesh,
            i
        )
        exportKratosConditions(
            model_part,
            smp_mesh,
            i
        )

        data = DataContainer()
        data[CUBA.MATERIAL_ID] = i
        smp_mesh.data = data

        pressure = 'empty' 
        velocity = 'empty'

        if model_part.GetProperties(0)[i].GetValue(IMPOSED_PRESSURE) == 1:
            pressure = model_part.GetProperties(0)[i].GetValue(PRESSURE)
        if model_part.GetProperties(0)[i].GetValue(IMPOSED_VELOCITY_X) == 1:
            velocity = (
                model_part.GetProperties(0)[i].GetValue(IMPOSED_VELOCITY_X_VALUE),
                model_part.GetProperties(0)[i].GetValue(IMPOSED_VELOCITY_Y_VALUE),
                model_part.GetProperties(0)[i].GetValue(IMPOSED_VELOCITY_Z_VALUE)
            )

        smp_bc = {
            'name':mesh_name,
            'pressure':pressure,
            'velocity':velocity
        }

        smp_bcs.append(smp_bc)
        smp_meshes.append(smp_mesh)

    return {'meshes':smp_meshes,'bcs':smp_bcs}

def exportKratosNodes(src, dst, num_mesh):
    """ Parses all kratos nodes to simphony points

    Iterates over all nodes in the kratos mesh (src) and
    converts them to simphony points (dst). While doing this operation
    any node/point that has not currently been mapped will have his uuid
    added in the 'id_map' of the wrapper

    """

    for node in src.GetNodes(num_mesh):

        data = {}

        getNodalData(data, node)

        point_uid = None

        if node.Id in id_to_uuid_node_map:
            point_uid = id_to_uuid_node_map[node.Id]

        point = SPoint(
            coordinates=(node.X, node.Y, node.Z),
            data=DataContainer(data),
            uid=point_uid
        )

        pid = dst.add_point(point)

        id_to_uuid_node_map[node.Id] = pid

def exportKratosElements(src, dst, num_mesh):
    """ Parses all kratos elements to simphony cells

    Iterates over all elements in the kratos mesh (src) and
    converts them to simphony cells (dst). While doing this operation
    any element/cell that has not currently been mapped will have his uuid
    added in the 'id_map' of the wrapper

    """

    for element in src.GetElements(num_mesh):

        element_uid = None

        if element.Id in id_to_uuid_element_map:
            element_uid = id_to_uuid_element_map[element.Id]

        point_list = [
            id_to_uuid_node_map[pointl.Id]
            for pointl in element.GetNodes()
        ]

        cell = SCell(
            points=point_list,
            uid=element_uid
        )

        cid = dst.add_cell(cell)

        id_to_uuid_element_map[element.Id] = cid

def exportKratosConditions(src, dst, num_mesh):
    """ Parses all kratos conditions to simphony faces

    Iterates over all conditions in the kratos mesh (src) and
    converts them to simphony faces (dst). While doing this operation
    any condition/face that has not currently been mapped will have
    his uuid added in the 'id_map' of the wrapper

    """

    for condition in src.GetConditions(num_mesh):

        condition_uid = None

        if condition.Id in id_to_uuid_condition_map:
            condition_uid = id_to_uuid_condition_map[condition.Id]

        point_list = [
            id_to_uuid_node_map[point.Id]
            for point in condition.GetNodes()
        ]

        face = SFace(
            points=point_list,
            uid=condition_uid
        )

        fid = dst.add_face(face)

        id_to_uuid_condition_map[condition.Id] = fid
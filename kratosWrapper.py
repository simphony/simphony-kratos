""" Template to convert modelparts from kratos to simphony 

This file show an example of how to utilze the kratos wrappers
in order to import or export models from KratosMultiphysics

"""

# Kratos will works with Python3 by default.
# This line makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Simphony Imports
from simphony.core.cuba import CUBA
from simphony.core.data_container import DataContainer

from simphony.cuds.mesh import Mesh as SimphonyMesh
from simphony.cuds.mesh import Point as SimphonyPoint
from simphony.cuds.mesh import Edge as SimphonyEdge
from simphony.cuds.mesh import Face as SimphonyFace
from simphony.cuds.mesh import Cell as SimphonyCell

# Kratos Imports 
from KratosMultiphysics import *

# Uuid and other dependences
from uuid import *

class KratosWrapper(object):

    def __init__(self):
        # Uid <-> KratosId mapping
        self.id_to_uuid_node_map = {}
        self.uuid_to_id_node_map = {}
        self.id_to_uuid_element_map = {}
        self.uuid_to_id_element_map = {}
        self.id_to_uuid_condition_map = {}
        self.uuid_to_id_condition_map = {}

        # Containers
        self.particleContainers = {}
        self.meshes = {}
        self.lattices = {}

        # Id's stuff
        self.free_id = 0

        # Initialization
        self.initialize()

    def __addNodalVariablesToModelpart(self,model_part):
        pass

    def __exportKratosNodes(self,src,dst):
        pass

    def __exportKratosElements(self,src,dst):
        pass

    def __exportKratosConditions(self,src,dst):
        pass

    def __importKratosNodes(self,src,dst):
        pass

    def __importKratosElements(self,src,dst):
        pass

    def __importKratosConditions(self,src,dst):
        pass

    def exportMesh(self,mesh):
        simphony_mesh = SimphonyMesh()

        self.__exportKratosNodes(mesh,simphony_mesh)
        self.__exportKratosElements(mesh,simphony_mesh)
        self.__exportKratosConditions(mesh,simphony_mesh)

        if not self.uuid_to_id_node_map:
            self.uuid_to_id_node_map = {v: k for k, v in self.id_to_uuid_node_map.items()}
        if not self.uuid_to_id_element_map:
            self.uuid_to_id_element_map = {v: k for k, v in self.id_to_uuid_element_map.items()}
        if not self.uuid_to_id_condition_map:
            self.uuid_to_id_condition_map = {v: k for k, v in self.id_to_uuid_condition_map.items()}

        return simphony_mesh

    def importMesh(self,mesh,name):
        model_part = ModelPart(name)

        self.__importKratosNodes(mesh,model_part)
        self.__importKratosElements(mesh,model_part)
        self.__importKratosConditions(mesh,model_part)

        self.__addNodalVariablesToModelpart(model_part)

        if not self.id_to_uuid_node_map:
            self.id_to_uuid_node_map = {v: k for k, v in self.uuid_to_id_node_map.items()}
        if not self.id_to_uuid_element_map:
            self.id_to_uuid_element_map = {v: k for k, v in self.uuid_to_id_element_map.items()}
        if not self.id_to_uuid_condition_map:
            self.id_to_uuid_condition_map = {v: k for k, v in self.uuid_to_id_condition_map.items()}

        return model_part

    def initialize(self):
        pass

    def run(self,mesh):
        pass




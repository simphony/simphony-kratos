""" test_kratosWrapper module

This module contains the unitary tests for the
kratosWrapper class.

"""
import os

from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MeshingApplication import *

import unittest

from simphony.core.cuba import CUBA
from simphony.engine import kratos


class TestKratosCFDWrapper(unittest.TestCase):

    def setUp(self):
        """ Creates the necessary information tu run the tests

        """

        self.fluid_path = os.path.join(
            os.path.dirname(__file__),
            '3balls'
        )
        self.rigid_path = os.path.join(
            os.path.dirname(__file__),
            '3ballsDEM_FEM_boundary'
        )

        self.time_step = 0.001
        self.num_steps = 600

    def test_run(self):
        """ Test if cfd can run

        """

        utils = kratos.DEM_Utils()
        wrapper = kratos.DEMWrapper()

        wrapper.CM[CUBA.TIME_STEP] = self.time_step
        wrapper.CM[CUBA.NUMBER_OF_TIME_STEPS] = self.num_steps

        # Set the meshes that are part of the fluid
        wrapper.SPE[kratos.CUBAExt.FLUID_PARTICLES] = [
            "fluid_0"
        ]

        # Set the meshes that are part of the fluid
        wrapper.SPE[kratos.CUBAExt.STRUCTURE_MESHES] = [
            "solid_0"
        ]

        kratos_model_f = utils.read_modelpart_as_particles(
            self.fluid_path, "fluid"
        )

        kratos_model_s = utils.read_modelpart_as_mesh(
            self.rigid_path, "solid"
        )

        for mesh in kratos_model_f['meshes']:
            wrapper.add_dataset(mesh)

        for mesh in kratos_model_s['meshes']:
            wrapper.add_dataset(mesh)

        for i in xrange(0, 1):
            wrapper.run()

    def tear_down(self):
        pass


if __name__ == '__main__':
    unittest.main()

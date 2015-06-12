""" test_kratosWrapper module

This module contains the unitary tests for the
kratosWrapper class.

"""
import sys
import os

sys.path.append(os.path.abspath('wrappers/tests/dem'))

from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MeshingApplication import *

import unittest

from simphony.core.cuba import CUBA

from wrappers import kratos_DEM_wrapper as DEMPackEngine


class TestKratosCFDWrapper(unittest.TestCase):

    def setUp(self):
        """ Creates the necessary information tu run the tests

        """

        self.fluid_path = "wrappers/tests/dem/3balls"
        self.rigid_path = "wrappers/tests/dem/3ballsDEM_FEM_boundary"

        self.time_step = 0.001
        self.num_steps = 600

    def test_run(self):
        """ Test if cfd can run

        """

        wrapper = DEMPackEngine.DEMPackWrapper()

        wrapper.CM[CUBA.TIME_STEP] = self.time_step
        wrapper.CM[CUBA.NUMBER_OF_TIME_STEPS] = self.num_steps

        mesh = wrapper.read_modelpart(
            self.fluid_path,
            self.rigid_path
        )

        wrapper.setMeshData(mesh)
        wrapper.add_mesh(mesh)

        for i in xrange(0, 1):
            wrapper.run()

        wrapper.finalize()

    def tear_down(self):
        pass


if __name__ == '__main__':
    unittest.main()

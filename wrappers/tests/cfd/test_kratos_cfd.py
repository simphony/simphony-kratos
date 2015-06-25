""" test_kratosWrapper module

This module contains the unitary tests for the
kratosWrapper class.

"""
import sys
import os

from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MeshingApplication import *

import unittest

from simphony.core.cuba import CUBA

from wrappers import kratos_CFD_wrapper as CFDengine
from wrappers.tests.cfd import ProjectParameters

class TestKratosCFDWrapper(unittest.TestCase):

    def setUp(self):
        """ Creates the necessary information tu run the tests

        """

        self.path = "wrappers/tests/cfd/CFD_exampleFluid"

        self.bc_vel = {'inlet': (1.0, 0.0, 0.0),
                       'outlet': 'zeroGradient',
                       'wall': (0.0, 0.0, 0.0),
                       'frontAndBack': 'empty'}

        self.bc_pre = {'inlet': 'zeroGradient',
                       'outlet': 0.0,
                       'wall': 'zeroGradient',
                       'frontAndBack': 'empty'}

        self.time_step = 0.001
        self.num_steps = 5

    def test_run(self):
        """ Test if cfd can run

        """

        wrapper = CFDengine.CFDWrapper()

        wrapper.CM[CUBA.TIME_STEP] = self.time_step
        wrapper.CM[CUBA.NUMBER_OF_TIME_STEPS] = self.num_steps
        wrapper.BC[CUBA.VELOCITY] = self.bc_vel
        wrapper.BC[CUBA.PRESSURE] = self.bc_pre

        mesh = wrapper.read_modelpart(self.path)

        wrapper.setMeshData(mesh)
        wrapper.add_mesh(mesh)

        for i in xrange(0, wrapper.CM[CUBA.NUMBER_OF_TIME_STEPS]):
            wrapper.run()

        wrapper.finalize()

    def tear_down(self):
        pass


if __name__ == '__main__':
    unittest.main()

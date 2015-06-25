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
from wrappers.kratos_utils import read_modelpart
from wrappers.cuba_extension import CUBAExt
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

        wrapper.SPE[CUBAExt.FLUID_MESHES] = "fluid_0,fluid_1,fluid_2,fluid_3,fluid_4"

        # reads kratos data so its interpretable by simphony
        kratos_model = read_modelpart(self.path)
        # mesh.name = "fluid"

        wrapper.BC[CUBA.VELOCITY] = {}
        wrapper.BC[CUBA.PRESSURE] = {}

        for mesh in kratos_model['meshes']:
            wrapper.add_mesh(mesh)

        for bc in kratos_model['bcs']:
            wrapper.BC[CUBA.VELOCITY][bc['name']] = bc['velocity']
            wrapper.BC[CUBA.PRESSURE][bc['name']] = bc['pressure']

        print(wrapper.BC[CUBA.VELOCITY])
        print(wrapper.BC[CUBA.PRESSURE])

        for i in xrange(0, wrapper.CM[CUBA.NUMBER_OF_TIME_STEPS]):
           wrapper.run()

    def tear_down(self):
        pass


if __name__ == '__main__':
    unittest.main()

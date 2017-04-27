""" test_kratosWrapper module

This module contains the unitary tests for the
kratosWrapper class.

"""
import os
import unittest

from simphony.api import CUDS, Simulation
from simphony.cuds.meta import api

# TODO: Utils now belong to probably another package
from simkratos.CFD.kratos_CFD_utils import CFD_Utils


class TestKratosCFDWrapper(unittest.TestCase):

    def setUp(self):
        """ Creates the necessary information tu run the tests

        """

        self.path = os.path.join(
            os.path.dirname(__file__),
            "CFD_exampleFluid"
        )

    def test_run(self):
        """ Test if cfd can run

        """

        path = os.path.join(os.path.dirname(__file__), "CFD_exampleFluid")

        cuds = CUDS(name='example_kratos_cfd_simulatiob')

        itime = api.IntegrationTime(name="cfd_integration_time")
        itime.time = 0.0001
        itime.step = 0.0025
        itime.final = 0.0075
        cuds.add([itime])

        cfd_utils = CFD_Utils()

        model_fluid = cfd_utils.read_modelpart(path)

        for model in [model_fluid]:
            cuds.add(list(model['datasets']))
            cuds.add(list(model['conditions']))
            cuds.add(list(model['materials']))
            cuds.add([model['pe']])

        # Create the simulation and run the problem
        sim = Simulation(cuds, "KRATOS_CFD", engine_interface=True)
        sim.run()

    def tear_down(self):
        pass


if __name__ == '__main__':
    unittest.main()

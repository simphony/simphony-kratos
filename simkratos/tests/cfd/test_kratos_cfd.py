""" test_kratosWrapper module

This module contains the unitary tests for the
kratosWrapper class.

"""
import os
import unittest

from simphony.core.cuba import CUBA

from simphony.api import CUDS, Simulation
from simphony.cuds.meta import api

# TODO: Utils now belong to probably another package
from simphony.engine import EngineInterface
from simphony.engine import kratos


class TestKratosCFDWrapper(unittest.TestCase):

    def setUp(self):
        """ Creates the necessary information tu run the tests

        """

        self.path = os.path.join(
            os.path.dirname(__file__),
            "CFD_exampleFluid"
        )

        self.time_step = 0.001
        self.num_steps = 5

    def test_run(self):
        """ Test if cfd can run

        """

        # Add the problem path to the script
        path = os.path.join(os.path.dirname(__file__), "CFD_exampleFluid")

        cuds = CUDS(name='example_kratos_cfd_simulatiob')

        # Integration time:
        itime = api.IntegrationTime(name="cfd_integration_time")
        itime.time = 0.0001
        itime.step = 0.0025
        itime.final = 0.0075  # 5 Kratos Timesteps
        cuds.add(itime)

        # Utils are used to read an existing Kratos model as raw data so we can
        # initialize the correct simphony datasets
        utils = kratos.CFD_Utils()

        # Reads Kratos mpda as a simphony data.
        model = utils.read_modelpart(path)

        # Add the datasets readed from the conversor.
        for dataset in model['datasets']:
            cuds.add(dataset)

        # Add the boundary contitions from the conversor
        for condition in model['conditions']:
            cuds.add(condition)

        # Add the materials contitions from the conversor
        for material in model['materials']:
            cuds.add(material)

        # Create the simulation and run the problem
        sim = Simulation(cuds, "KRATOS_CFD", engine_interface=True)
        sim.run()

    def tear_down(self):
        pass


if __name__ == '__main__':
    unittest.main()

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


def abs_path(relPath):
    return os.path.join(os.path.dirname(__file__), relPath)


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
        pathParticles = abs_path("3balls")
        pathSolid = abs_path("3ballsDEM_FEM_boundary")

        cuds = CUDS(name='example_kratos_dem_somulation')

        # Integration time:
        itime = api.IntegrationTime(name="dem_integration_time")
        itime.time = 0.0001
        itime.step = 0.001
        itime.final = 60 * itime.step
        cuds.add(itime)

        # Utils are used to read an existing Kratos model as raw data so we can
        # initialize the correct simphony datasets
        utils = kratos.DEM_Utils()

        # Reads Kratos mpda as a simphony data.
        model_particles = utils.read_modelpart_as_particles(pathParticles)
        model_solid = utils.read_modelpart_as_mesh(pathSolid)

        # Add all models
        for model in [model_particles, model_solid]:
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
        sim = Simulation(cuds, "KRATOS_DEM", engine_interface=True)
        sim.run()

    def tear_down(self):
        pass


if __name__ == '__main__':
    unittest.main()

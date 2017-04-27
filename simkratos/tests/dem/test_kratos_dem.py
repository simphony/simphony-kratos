""" test_kratosWrapper module

This module contains the unitary tests for the
kratosWrapper class.

"""
import os
import unittest

from simphony.api import CUDS, Simulation
from simphony.cuds.meta import api

# TODO: Utils now belong to probably another package
from simkratos.DEM.kratos_DEM_utils import DEM_Utils


class TestKratosDEMWrapper(unittest.TestCase):

    def setUp(self):
        """ Creates the necessary information tu run the tests

        """

        self.path = os.path.join(
            os.path.dirname(__file__),
            "DEM_exampleFluid"
        )

    def test_run(self):
        """ Test if cfd can run

        """

        pathParticles = os.path.join(
            "/home/travis/build/simphony/simphony-kratos/simkratos/tests/dem",
            "3balls"
        )
        pathSolid = os.path.join(
            "/home/travis/build/simphony/simphony-kratos/simkratos/tests/dem",
            "3ballsDEM_FEM_boundary"
        )

        cuds = CUDS(name='example_kratos_dem_somulation')

        itime = api.IntegrationTime(name="dem_integration_time")
        itime.time = 0.0001
        itime.step = 0.001
        itime.final = 60 * itime.step
        cuds.add([itime])

        utils = DEM_Utils()

        print(pathParticles)

        model_particles = utils.read_modelpart_as_particles(pathParticles)
        model_solid = utils.read_modelpart_as_mesh(pathSolid)

        for model in [model_particles, model_solid]:
            cuds.add(list(model['datasets']))
            cuds.add(list(model['conditions']))
            cuds.add(list(model['materials']))
            cuds.add([model['pe']])

        sim = Simulation(cuds, "KRATOS_DEM", engine_interface=True)
        sim.run()

    def tear_down(self):
        pass


if __name__ == '__main__':
    unittest.main()

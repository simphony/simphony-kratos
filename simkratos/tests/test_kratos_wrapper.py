""" test_kratosWrapper module

This module contains the unitary tests for the
kratosWrapper class.

"""

import unittest

from simphony.cuds.mesh import Mesh, Point, Edge, Face, Cell
from simphony.cuds.particles import Particles

from simkratos.kratosWrapper import KratosWrapper


class TestKratosWrapper(unittest.TestCase):

    def setUp(self):
        """ Creates a couple of empty meshes to perform the tests

        Creates a couple of empty meshes and and a set of points
        to tests all the mesh methods

        """

        pass

    def test_add_mesh(self):
        """ Test if a mesh can be added to the wrapper

        """

        pass

    def test_add_existing_mesh(self):
        pass

    def test_delete_mesh(self):
        """ Test if a mesh can be deleted from the wrapper

        """

        pass

    def test_get_mesh(self):
        """ Test if a mesh can be get from the wrapper

        """

        pass

    def test_iter_meshes(self):
        """ Test if meshes can be iterated.

        """

        pass

    def test_iter_meshes_subset(self):
        """ Test if meshes can be iterated.

        """

        pass

    def test_add_particles(self):
        """ Test if a particle container can be added to the wrapper

        """

        pass

    def test_change_mesh_name(self):
        """ Test if the wrapper can correctly handle mesh name
        changes

        """

        pass

    def test_change_mesh_name_iter(self):
        """ Test if the wrapper can correctly handle mesh name
        changes in iterators

        """

        pass


if __name__ == '__main__':
    unittest.main()

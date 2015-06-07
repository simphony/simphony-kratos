""" test_kratosWrapper module

This module contains the unitary tests for the
kratosWrapper class.

"""

import unittest

from simphony.cuds.mesh import Mesh, Point, Edge, Face, Cell
from simphony.core.data_container import DataContainer

from wrappers.kratosWrapper import KratosWrapper


class TestKratosWrapper(unittest.TestCase):

    def setUp(self):
        """ Creates a couple of empty meshes to perform the tests

        Creates a couple of empty meshes and and a set of points
        to tests all the mesh methods

        """

        self.mesh1 = Mesh(name="foo1")
        self.mesh2 = Mesh(name="foo2")

        self.points = [
            Point((0.0, 0.0, 0.0)),
            Point((1.0, 0.0, 0.0)),
            Point((0.0, 1.0, 0.0)),
            Point((0.0, 0.0, 1.0)),
            Point((1.0, 0.0, 1.0)),
            Point((0.0, 1.0, 1.0))
        ]

        puuid1 = []
        puuid2 = []
        euuids = []
        fuuids = []
        cuuids = []

        for point in self.points:
            puuid = self.mesh1.add_point(point)
            puuid1.append(puuid)

        for point in self.points:
            puuid = self.mesh2.add_point(point)
            puuid2.append(puuid)

        edges = [
            Edge(puuid1[0:2]),
            Edge(puuid1[1:3])
        ]

        faces = [
            Face(puuid1[0:3]),
            Face(puuid1[1:4])
        ]

        cells = [
            Cell(puuid2[0:4]),
            Cell(puuid2[1:5])
        ]

        for edge in edges:
            euuid = self.mesh1.add_edge(edge)
            euuids.append(euuid)

        for face in faces:
            fuuid = self.mesh1.add_face(face)
            fuuids.append(fuuid)

        for cell in cells:
            cuuid = self.mesh2.add_cell(cell)
            cuuids.append(cuuid)

    def test_add_mesh(self):
        """ Test if a mesh can be added to the wrapper

        """

        wrapper = KratosWrapper()
        wrapper.add_mesh(self.mesh1)

        num_meshes = 0
        for mesh in wrapper.iter_meshes():
            num_meshes = num_meshes + 1

        self.assertEqual(num_meshes, 1)

        wrapper.add_mesh(self.mesh2)

        num_meshes = 0
        for mesh in wrapper.iter_meshes():
            num_meshes = num_meshes + 1

        self.assertEqual(num_meshes, 2)

    def test_delete_mesh(self):
        """ Test if a mesh can be deleted from the wrapper

        """

        wrapper = KratosWrapper()
        wrapper.add_mesh(self.mesh1)

        num_meshes = 0
        for mesh in wrapper.iter_meshes():
            num_meshes = num_meshes + 1

        self.assertEqual(num_meshes, 1)

        wrapper.delete_mesh("foo1")

        num_meshes = 0
        for mesh in wrapper.iter_meshes():
            num_meshes = num_meshes + 1

        self.assertEqual(num_meshes, 0)

    def test_get_mesh(self):
        """ Test if a mesh can be get from the wrapper

        """

        wrapper = KratosWrapper()
        wrapper.add_mesh(self.mesh1)
        wrapper.add_mesh(self.mesh2)

        w_mesh1 = wrapper.get_mesh("foo1")

        for point in w_mesh1.iter_points():
            wp = self.mesh1.get_point(point.uid)
            self.assertEqual(point.uid, wp.uid)

        for edge in w_mesh1.iter_edges():
            we = self.mesh1.get_edge(edge.uid)
            self.assertEqual(edge.uid, we.uid)

        for face in w_mesh1.iter_faces():
            wf = self.mesh1.get_face(face.uid)
            self.assertEqual(face.uid, wf.uid)

        for cell in w_mesh1.iter_cells():
            wc = self.mesh1.get_cell(cell.uid)
            self.assertEqual(cell.uid, wc.uid)

    def test_iter_meshes(self):
        """ Test if meshes can be iterated.

        """

        wrapper = KratosWrapper()
        wrapper.add_mesh(self.mesh1)
        wrapper.add_mesh(self.mesh2)

        meshes_n = ["foo1", "foo2"]
        meshes_w = [m for m in wrapper.iter_meshes()]

        self.assertItemsEqual(meshes_n, meshes_w)

    def test_add_particles(self):
        """ Test if a particle container can be added to the wrapper

        """

        wrapper = KratosWrapper()
        with self.assertRaises(NotImplementedError):
            wrapper.add_particles(DataContainer())

    def test_get_particles(self):
        """ Test if a particle container can be get from the wrapper

        """

        wrapper = KratosWrapper()
        with self.assertRaises(NotImplementedError):
            wrapper.get_particles('')

    def test_delete_particles(self):
        """ Test if a particle container can be deleted from the wrapper

        """

        wrapper = KratosWrapper()
        with self.assertRaises(NotImplementedError):
            wrapper.delete_particles('')

    def test_iter_particles(self):
        """ Test if particle containers can be iterated.

        """

        wrapper = KratosWrapper()
        with self.assertRaises(NotImplementedError):
            wrapper.iter_particles()

    def test_add_lattice(self):
        """ Test if a lattice can be added to the wrapper

        """

        wrapper = KratosWrapper()
        with self.assertRaises(NotImplementedError):
            wrapper.add_lattice(DataContainer())

    def test_get_lattice(self):
        """ Test if a lattice can be get from the wrapper

        """

        wrapper = KratosWrapper()
        with self.assertRaises(NotImplementedError):
            wrapper.get_lattice('')

    def test_delete_lattice(self):
        """ Test if a lattice can be deleted from the wrapper

        """

        wrapper = KratosWrapper()
        with self.assertRaises(NotImplementedError):
            wrapper.delete_lattice('')

    def test_iter_lattices(self):
        """ Test if lattices can be iterated.

        """

        wrapper = KratosWrapper()
        with self.assertRaises(NotImplementedError):
            wrapper.iter_lattices()

    def test_run_cfd(self):
        """ Test the execution of one step of the simulation

        """

        pass

    def test_run_dem(self):
        """ Test the execution of one step of the simulation

        """

        pass


if __name__ == '__main__':
    unittest.main()

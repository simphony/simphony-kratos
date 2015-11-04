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
            puuid = self.mesh1.add_points([point])
            puuid1.append(puuid)

        for point in self.points:
            puuid = self.mesh2.add_points([point])
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
            euuid = self.mesh1.add_edges([edge])
            euuids.append(euuid)

        for face in faces:
            fuuid = self.mesh1.add_faces([face])
            fuuids.append(fuuid)

        for cell in cells:
            cuuid = self.mesh2.add_cells([cell])
            cuuids.append(cuuid[0])

    def test_add_mesh(self):
        """ Test if a mesh can be added to the wrapper

        """

        wrapper = KratosWrapper()
        wrapper.add_dataset(self.mesh1)

        num_meshes = 0
        for mesh in wrapper.iter_datasets():
            num_meshes = num_meshes + 1

        self.assertEqual(num_meshes, 1)

        wrapper.add_dataset(self.mesh2)

        num_meshes = 0
        for _ in wrapper.iter_datasets():
            num_meshes = num_meshes + 1

        self.assertEqual(num_meshes, 2)

    def test_delete_mesh(self):
        """ Test if a mesh can be deleted from the wrapper

        """

        wrapper = KratosWrapper()
        wrapper.add_dataset(self.mesh1)

        num_meshes = 0
        for _ in wrapper.iter_datasets():
            num_meshes = num_meshes + 1

        self.assertEqual(num_meshes, 1)

        wrapper.remove_dataset("foo1")

        num_meshes = 0
        for _ in wrapper.iter_datasets():
            num_meshes = num_meshes + 1

        self.assertEqual(num_meshes, 0)

    def test_get_mesh(self):
        """ Test if a mesh can be get from the wrapper

        """

        wrapper = KratosWrapper()
        wrapper.add_dataset(self.mesh1)
        wrapper.add_dataset(self.mesh2)

        w_mesh1 = wrapper.get_dataset("foo1")

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
        wrapper.add_dataset(self.mesh1)
        wrapper.add_dataset(self.mesh2)

        meshes_n = ["foo1", "foo2"]
        meshes_w = [m.name for m in wrapper.iter_datasets()]

        self.assertItemsEqual(meshes_n, meshes_w)

    def test_iter_meshes_subset(self):
        """ Test if meshes can be iterated.

        """

        wrapper = KratosWrapper()
        wrapper.add_dataset(self.mesh1)
        wrapper.add_dataset(self.mesh2)

        meshes_n = ["foo1"]
        meshes_w = [m.name for m in wrapper.iter_datasets(meshes_n)]

        self.assertItemsEqual(meshes_n, meshes_w)

    def test_add_particles(self):
        """ Test if a particle container can be added to the wrapper

        """

        wrapper = KratosWrapper()
        with self.assertRaises(TypeError):
            wrapper.add_dataset(Particles(name="test"))

    def test_change_mesh_name(self):
        """ Test if the wrapper can correctly handle mesh name
        changes

        """

        wrapper = KratosWrapper()
        wrapper.add_dataset(self.mesh1)
        proxy_mesh = wrapper.get_dataset(self.mesh1.name)

        proxy_mesh.name = "fooRenamed"

        w_mesh1 = wrapper.get_dataset("fooRenamed")

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

    def test_change_mesh_name_iter(self):
        """ Test if the wrapper can correctly handle mesh name
        changes in iterators

        """

        wrapper = KratosWrapper()
        wrapper.add_dataset(self.mesh1)
        wrapper.add_dataset(self.mesh2)
        proxy_mesh = wrapper.get_dataset(self.mesh2.name)

        proxy_mesh.name = "fooRenamed"

        meshes_n = ["foo1", "fooRenamed"]
        meshes_w = [m.name for m in wrapper.iter_datasets()]

        self.assertItemsEqual(meshes_n, meshes_w)


if __name__ == '__main__':
    unittest.main()

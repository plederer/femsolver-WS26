# Tests for your class Mesh and MakeStructuredMesh (from femsolver/mesh.py).
# They will only pass once you have implemented both.  Run with:  python3 test_mesh.py
import numpy as np
from femsolver import *


def signed_areas(mesh):
    p = mesh.points[mesh.elements]          # shape (ne, 3, 2)
    a, b, c = p[:, 0], p[:, 1], p[:, 2]
    return 0.5 * ((b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1])
                  - (c[:, 0] - a[:, 0]) * (b[:, 1] - a[:, 1]))


def check_mesh(mesh, area):
    assert mesh.points.shape == (mesh.nv, 2)
    assert mesh.elements.shape == (mesh.ne, 3)
    assert mesh.bnd_edges.shape == (mesh.nbnd, 2)
    assert len(mesh.bnd_labels) == mesh.nbnd
    areas = signed_areas(mesh)
    assert np.all(areas > 0), "elements must be ordered counterclockwise"
    assert np.isclose(areas.sum(), area), "areas do not add up"


def test_structured_mesh(nx, ny, xmin=0.0, xmax=1.0, ymin=0.0, ymax=1.0):
    mesh = MakeStructuredMesh(nx, ny, xmin, xmax, ymin, ymax)
    assert mesh.nv == (nx + 1) * (ny + 1), "wrong number of vertices"
    assert mesh.ne == 2 * nx * ny, "wrong number of elements"
    assert mesh.nbnd == 2 * (nx + ny), "wrong number of boundary edges"
    check_mesh(mesh, (xmax - xmin) * (ymax - ymin))

    for label, n, coord, value in [("bottom", nx, 1, ymin), ("right", ny, 0, xmax),
                                   ("top", nx, 1, ymax), ("left", ny, 0, xmin)]:
        edges = mesh.GetBoundaryEdges(label)
        assert len(edges) == n, f"wrong number of edges on '{label}'"
        assert np.allclose(mesh.points[edges][:, :, coord], value), \
            f"an edge labeled '{label}' does not lie on this boundary"


def test_ngsolve_mesh():
    try:
        from netgen.geom2d import unit_square
    except ImportError:
        print("NGSolve not installed, skipping test_ngsolve_mesh")
        return
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))
    check_mesh(mesh, 1.0)
    assert set(mesh.bnd_labels) == {"bottom", "right", "top", "left"}


if __name__ == "__main__":
    assert Mesh().nv == 0
    test_structured_mesh(1, 1)
    test_structured_mesh(3, 2)
    test_structured_mesh(4, 5, -1.0, 2.0, 0.0, 0.5)
    test_ngsolve_mesh()
    print("all tests passed")

# Tests for MakeStructuredMesh. Run with:  python3 test_mesh.py
import numpy as np
from myfem import *


def signed_areas(mesh):
    p = mesh.points[mesh.elements]          # shape (ne, 3, 2)
    a, b, c = p[:, 0], p[:, 1], p[:, 2]
    return 0.5 * ((b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1])
                  - (c[:, 0] - a[:, 0]) * (b[:, 1] - a[:, 1]))


def test_structured_mesh(nx, ny, xmin=0.0, xmax=1.0, ymin=0.0, ymax=1.0):
    mesh = MakeStructuredMesh(nx, ny, xmin, xmax, ymin, ymax)

    assert mesh.nv == (nx + 1) * (ny + 1), "wrong number of vertices"
    assert mesh.ne == 2 * nx * ny, "wrong number of elements"
    assert mesh.nbnd == 2 * (nx + ny), "wrong number of boundary edges"

    areas = signed_areas(mesh)
    assert np.all(areas > 0), "elements must be ordered counterclockwise"
    assert np.isclose(areas.sum(), (xmax - xmin) * (ymax - ymin)), "areas do not add up"

    for label, n, coord, value in [("bottom", nx, 1, ymin), ("right", ny, 0, xmax),
                                   ("top", nx, 1, ymax), ("left", ny, 0, xmin)]:
        edges = mesh.GetBoundaryEdges(label)
        assert len(edges) == n, f"wrong number of edges on '{label}'"
        assert np.allclose(mesh.points[edges][:, :, coord], value), \
            f"an edge labeled '{label}' does not lie on this boundary"


if __name__ == "__main__":
    test_structured_mesh(1, 1)
    test_structured_mesh(3, 2)
    test_structured_mesh(4, 5, -1.0, 2.0, 0.0, 0.5)
    print("all tests passed")

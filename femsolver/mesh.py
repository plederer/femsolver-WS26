import numpy as np

__all__ = ["Mesh", "MakeStructuredMesh"]


class Mesh:
    """Triangular mesh of a two-dimensional domain.

    The mesh stores plain numpy arrays:

    points     : float array of shape (nv, 2)
                 coordinates of the vertices, points[i] = (x_i, y_i)
    elements   : int array of shape (ne, 3)
                 vertex numbers of every triangle, ordered counterclockwise
    bnd_edges  : int array of shape (nb, 2)
                 vertex numbers of every boundary edge
    bnd_labels : list of nb strings
                 name of the boundary part an edge belongs to, e.g. "left"

    Mesh()        creates an empty mesh (filled e.g. by MakeStructuredMesh),
    Mesh(ngmesh)  converts a Netgen/NGSolve mesh.
    """

    def __init__(self, ngmesh=None):
        self.points = np.zeros((0, 2))
        self.elements = np.zeros((0, 3), dtype=int)
        self.bnd_edges = np.zeros((0, 2), dtype=int)
        self.bnd_labels = []

        if ngmesh is not None:
            from .ngsinterface import ConvertNGSolveMesh
            (self.points, self.elements,
             self.bnd_edges, self.bnd_labels) = ConvertNGSolveMesh(ngmesh)

    @property
    def nv(self):
        """number of vertices"""
        return len(self.points)

    @property
    def ne(self):
        """number of elements (triangles)"""
        return len(self.elements)

    @property
    def nbnd(self):
        """number of boundary edges"""
        return len(self.bnd_edges)

    def GetBoundaryEdges(self, label):
        """Returns the boundary edges with the given label."""
        edges = [e for e, l in zip(self.bnd_edges, self.bnd_labels)
                 if l == label]
        return np.array(edges, dtype=int).reshape(-1, 2)

    def __str__(self):
        return (f"Mesh with {self.nv} vertices, {self.ne} elements "
                f"and {self.nbnd} boundary edges")


def MakeStructuredMesh(nx, ny, xmin=0.0, xmax=1.0, ymin=0.0, ymax=1.0):
    """Structured triangular mesh of [xmin, xmax] x [ymin, ymax].

    The rectangle is divided into nx * ny equal cells, and every cell is
    split into two triangles. The boundary edges are labeled
    "bottom", "right", "top" and "left".
    """
    xs = np.linspace(xmin, xmax, nx + 1)
    ys = np.linspace(ymin, ymax, ny + 1)

    def v(i, j):
        """number of the vertex at (xs[i], ys[j])"""
        return j * (nx + 1) + i

    # create points in lexicographic order and save physical coordinates in a list
    points = []
    for j in range(ny+1):
        for i in range(nx+1):
            points.append([xs[i%(nx+1)], ys[j%(ny+1)]])

    # every cell is split along the diagonal from (i, j) to (i+1, j+1)
    elements = []
    for j in range(ny):
        for i in range(nx):
            elements.append((v(i, j), v(i + 1, j), v(i + 1, j + 1)))
            elements.append((v(i, j), v(i + 1, j + 1), v(i, j + 1)))

    # boundary edges, counterclockwise around the rectangle
    bnd_edges, bnd_labels = [], []
    for i in range(nx):
        bnd_edges.append((v(i, 0), v(i + 1, 0)))
        bnd_labels.append("bottom")
    for j in range(ny):
        bnd_edges.append((v(nx, j), v(nx, j + 1)))
        bnd_labels.append("right")
    for i in range(nx, 0, -1):
        bnd_edges.append((v(i, ny), v(i - 1, ny)))
        bnd_labels.append("top")
    for j in range(ny, 0, -1):
        bnd_edges.append((v(0, j), v(0, j - 1)))
        bnd_labels.append("left")

    mesh = Mesh()
    mesh.points = np.array(points)
    mesh.elements = np.array(elements, dtype=int)
    mesh.bnd_edges = np.array(bnd_edges, dtype=int)
    mesh.bnd_labels = bnd_labels
    return mesh

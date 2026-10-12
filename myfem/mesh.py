import numpy as np

__all__ = ["Mesh", "MakeStructuredMesh"]


class Mesh:
    """Triangular mesh of a two-dimensional domain.

    The mesh only stores plain numpy arrays:

    points     : float array of shape (nv, 2)
                 coordinates of the vertices, points[i] = (x_i, y_i)
    elements   : int array of shape (ne, 3)
                 vertex numbers of every triangle, ordered counterclockwise
    bnd_edges  : int array of shape (nb, 2)
                 vertex numbers of every boundary edge
    bnd_labels : list of nb strings
                 name of the boundary part an edge belongs to, e.g. "left"
    """

    def __init__(self, points, elements, bnd_edges, bnd_labels):
        self.points = np.asarray(points, dtype=float)
        self.elements = np.asarray(elements, dtype=int)
        self.bnd_edges = np.asarray(bnd_edges, dtype=int).reshape(-1, 2)
        self.bnd_labels = list(bnd_labels)

        if len(self.bnd_labels) != len(self.bnd_edges):
            raise ValueError("every boundary edge needs exactly one label")

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
        """Returns the boundary edges (pairs of vertex numbers) with the given label."""
        return np.array([e for e, l in zip(self.bnd_edges, self.bnd_labels) if l == label],
                        dtype=int).reshape(-1, 2)

    def __str__(self):
        return f"Mesh with {self.nv} vertices, {self.ne} elements and {self.nbnd} boundary edges"


def MakeStructuredMesh(nx, ny, xmin=0.0, xmax=1.0, ymin=0.0, ymax=1.0):
    """Structured triangular mesh of the rectangle [xmin, xmax] x [ymin, ymax].

    The rectangle is divided into nx * ny equal cells, and every cell is
    split into two triangles. The boundary edges are labeled
    "bottom", "right", "top" and "left".

    Returns a Mesh.
    """
    # TODO: implement (exercise sheet)
    raise NotImplementedError("MakeStructuredMesh is not implemented yet")

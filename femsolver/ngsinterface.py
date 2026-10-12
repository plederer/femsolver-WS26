"""Interface to meshes generated with Netgen/NGSolve.

This file is handed out; you do not need to change it.
NGSolve is only required if you actually convert an NGSolve mesh.
"""
import numpy as np

__all__ = ["ConvertNGSolveMesh"]


def ConvertNGSolveMesh(mesh):
    """Converts a two-dimensional Netgen or NGSolve mesh into plain arrays.

    mesh : ngsolve.Mesh or netgen.meshing.Mesh (2D, triangles only)

    Returns (points, elements, bnd_edges, bnd_labels):

    points     : float array (nv, 2), vertex coordinates
    elements   : int array (ne, 3), vertex numbers of the triangles, counterclockwise
    bnd_edges  : int array (nb, 2), vertex numbers of the boundary edges
    bnd_labels : list of nb strings, boundary names given in the geometry
    """
    import ngsolve

    if not isinstance(mesh, ngsolve.Mesh):
        mesh = ngsolve.Mesh(mesh)
    if mesh.dim != 2:
        raise ValueError("only two-dimensional meshes are supported")

    points = np.array([v.point for v in mesh.vertices], dtype=float)

    elements = []
    for el in mesh.Elements(ngsolve.VOL):
        if el.type != ngsolve.ET.TRIG:
            raise ValueError("only triangular meshes are supported")
        elements.append([v.nr for v in el.vertices])
    elements = np.array(elements, dtype=int)

    # make sure that all triangles are ordered counterclockwise
    p0, p1, p2 = (points[elements[:, k]] for k in range(3))
    det = (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1]) \
        - (p2[:, 0] - p0[:, 0]) * (p1[:, 1] - p0[:, 1])
    elements[det < 0] = elements[det < 0][:, [0, 2, 1]]

    bnd_edges, bnd_labels = [], []
    for el in mesh.Elements(ngsolve.BND):
        bnd_edges.append([v.nr for v in el.vertices])
        bnd_labels.append(el.mat)
    bnd_edges = np.array(bnd_edges, dtype=int).reshape(-1, 2)

    return points, elements, bnd_edges, bnd_labels

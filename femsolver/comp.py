"""Finite element spaces and grid functions."""
import re
import numpy as np
from .fem import H1FiniteElement
from .mesh import MeshTrafo

__all__ = ["H1", "GridFunction", "FindElement"]


class H1:
    """Continuous, piecewise linear finite element space on a mesh.

    The degrees of freedom (dofs) are the values at the vertices:
    dof number i belongs to vertex number i.

    dirichlet : regular expression for the boundary labels with Dirichlet
                conditions, e.g. "left|right" or ".*" (all boundaries)
    """

    def __init__(self, mesh, order=1, dirichlet=""):
        self.mesh = mesh
        self.order = order
        self.fe = H1FiniteElement(order)
        self.ndof = mesh.nv
        self.dirichlet = dirichlet

        # FreeDofs[i] is False for dofs on a Dirichlet boundary
        self.FreeDofs = np.ones(self.ndof, dtype=bool)
        for (a, b), label in zip(mesh.bnd_edges, mesh.bnd_labels):
            if dirichlet and re.fullmatch(dirichlet, label):
                self.FreeDofs[[a, b]] = False

    def GetDofs(self, nr):
        """Global dof numbers of element nr, in the order of the local shape functions."""
        return self.mesh.elements[nr]


class GridFunction:
    """A finite element function u_h = sum_i vec[i] * phi_i."""

    def __init__(self, space):
        self.space = space
        self.vec = np.zeros(space.ndof)

    def Set(self, f):
        """Interpolation: sets the dofs to the values f(x, y) at the vertices."""
        for i, (x, y) in enumerate(self.space.mesh.points):
            self.vec[i] = f(x, y)

    def Evaluate(self, x, y, nr=None):
        """Value of u_h at the physical point(s) (x, y) in element nr.
        If nr is None, the element containing the point (x, y) is searched."""
        mesh = self.space.mesh
        if nr is None:
            nr = FindElement(mesh, x, y)
        xhat, yhat = MeshTrafo(nr, mesh).mapinv(x, y)
        dofs = self.space.GetDofs(nr)
        shapes = self.space.fe.Evaluate(xhat, yhat)
        return self.vec[dofs] @ shapes


def FindElement(mesh, x, y, tol=1e-12):
    """Number of an element containing the point (x, y)."""
    for nr in range(mesh.ne):
        xhat, yhat = MeshTrafo(nr, mesh).mapinv(x, y)
        if xhat >= -tol and yhat >= -tol and xhat + yhat <= 1 + tol:
            return nr
    raise ValueError(f"point ({x}, {y}) is not in the mesh")

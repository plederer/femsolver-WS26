"""Finite element spaces (grid functions follow on the next sheet)."""
import re
import numpy as np
from .fem import H1FiniteElement

__all__ = ["H1"]


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

"""Finite elements on the reference triangle {(x,y): x, y >= 0, x + y <= 1}."""
import numpy as np

__all__ = ["H1FiniteElement"]


class H1FiniteElement:
    """Lagrange (H1-conforming) finite element on the reference triangle.

    Only order 1 for now: the three shape functions belong to the vertices
    (0,0), (1,0), (0,1), in this order:

        phi_0 = 1 - x - y,   phi_1 = x,   phi_2 = y
    """

    def __init__(self, order=1):
        if order != 1:
            raise NotImplementedError("only order 1 is implemented")
        self.order = order
        self.ndof = 3

    def Evaluate(self, xhat, yhat):
        """Shape functions at the point(s) (xhat, yhat).
        Returns an array of shape (ndof,) for a point, (ndof, npoints) for arrays."""
        xhat, yhat = np.asarray(xhat, dtype=float), np.asarray(yhat, dtype=float)
        return np.array([1 - xhat - yhat, xhat, yhat])

    def DEvaluate(self, xhat, yhat):
        """Gradients of the shape functions, array of shape (ndof, 2).
        For order 1 they are constant."""
        return np.array([[-1.0, -1.0],
                         [ 1.0,  0.0],
                         [ 0.0,  1.0]])

# Finite element space, grid functions and drawing.  Run with:  python example_gridfunction.py
import numpy as np
from femsolver import *

mesh = MakeStructuredMesh(8, 8)
V = H1(mesh, order=1, dirichlet="left|bottom")
print("ndof =", V.ndof, " free dofs:", np.sum(V.FreeDofs))
print("dofs of element 0:", V.GetDofs(0))

# interpolate a function and draw it
gfu = GridFunction(V)
gfu.Set(lambda x, y: np.sin(np.pi * x) * np.sin(np.pi * y))
Draw(gfu)
Draw(gfu, surface=True)

# the free dofs as a grid function
free = GridFunction(V)
free.vec[:] = V.FreeDofs
Draw(free)

# a linear function is reproduced exactly: evaluate at the (physical)
# integration points of every element
gfu.Set(lambda x, y: 1 + 2 * x - 3 * y)
ir = IntegrationRule(2)
for nr in range(mesh.ne):
    x, y = MeshTrafo(nr, mesh).map(ir.points[:, 0], ir.points[:, 1])
    assert np.allclose(gfu.Evaluate(x, y, nr), 1 + 2 * x - 3 * y)
print("linear function reproduced exactly")

# evaluation at a point without knowing the element
print("u_h(0.3, 0.4) =", gfu.Evaluate(0.3, 0.4), " exact:", 1 + 2 * 0.3 - 3 * 0.4)

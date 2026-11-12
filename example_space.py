# Finite element space H1 (femsolver IV).  Run with:  python example_space.py
import numpy as np
from femsolver import *

mesh = MakeStructuredMesh(8, 8)
for dirichlet in ["", "left|bottom", ".*"]:
    V = H1(mesh, order=1, dirichlet=dirichlet)
    print(f"dirichlet={dirichlet!r:15} ndof = {V.ndof}, free dofs = {np.sum(V.FreeDofs)}")
print("dofs of element 0:", V.GetDofs(0))

# shape functions on the reference triangle
fe = V.fe
print("shape functions at (1/3, 1/3):", fe.Evaluate(1/3, 1/3))
print("gradients:\n", fe.DEvaluate(1/3, 1/3))

# a mesh generated with NGSolve: unit square, all boundaries Dirichlet
from netgen.occ import unit_square
V = H1(Mesh(unit_square.GenerateMesh(maxh=0.3)), dirichlet=".*")
print("NGSolve mesh: ndof =", V.ndof, ", free dofs =", np.sum(V.FreeDofs))

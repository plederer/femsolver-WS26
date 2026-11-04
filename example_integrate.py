# Integration over meshes.  Run with:  python example_integrate.py
import numpy as np
from femsolver import *

mesh = MakeStructuredMesh(4, 4)
print("area of the unit square:     ", Integrate(lambda x, y: 1 + 0 * x, mesh))
print("int x*y over the unit square:", Integrate(lambda x, y: x * y, mesh, order=2), "(exact 0.25)")

# a smooth function: exact value 4 / pi^2
f = lambda x, y: np.sin(np.pi * x) * np.sin(np.pi * y)
exact = 4 / np.pi**2
print("\n   n   order 1     order 4")
for n in [2, 4, 8, 16]:
    mesh = MakeStructuredMesh(n, n)
    print(f"{n:4d}  {abs(Integrate(f, mesh, 1) - exact):.2e}  {abs(Integrate(f, mesh, 4) - exact):.2e}")

# a mesh from NGSolve: unit disk, area pi
from netgen.occ import Circle, OCCGeometry
for h in [0.4, 0.2, 0.1]:
    mesh = Mesh(OCCGeometry(Circle((0, 0), 1).Face(), dim=2).GenerateMesh(maxh=h))
    A = Integrate(lambda x, y: 1 + 0 * x, mesh)
    print(f"disk, maxh={h}: area {A:.6f}, error {abs(A - np.pi):.2e}")

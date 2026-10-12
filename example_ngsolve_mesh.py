# Generates a mesh with Netgen/NGSolve, converts it into plain numpy arrays
# and draws it. This example only uses the provided files and runs right away.
# Run with:  python3 example_ngsolve_mesh.py

from netgen.occ import Rectangle, OCCGeometry, X, Y

from femsolver.ngsinterface import ConvertNGSolveMesh
from femsolver.draw import DrawMesh

# geometry: unit square with named boundary parts
shape = Rectangle(1, 1).Face()
shape.edges.Min(Y).name = "bottom"
shape.edges.Max(X).name = "right"
shape.edges.Max(Y).name = "top"
shape.edges.Min(X).name = "left"

# mesh with maximal element size 0.25
ngmesh = OCCGeometry(shape, dim=2).GenerateMesh(maxh=0.25)

# conversion into plain arrays
points, elements, bnd_edges, bnd_labels = ConvertNGSolveMesh(ngmesh)

print("points     ", points.shape, "\n", points[:5], "...")
print("elements   ", elements.shape, "\n", elements[:5], "...")
print("bnd_edges  ", bnd_edges.shape, "\n", bnd_edges[:5], "...")
print("bnd_labels ", len(bnd_labels), "\n", bnd_labels[:5], "...")

# DrawMesh only needs an object with these four attributes.
# Later, this will be your own class Mesh in femsolver/mesh.py.
class SimpleMesh:
    pass

mesh = SimpleMesh()
mesh.points = points
mesh.elements = elements
mesh.bnd_edges = bnd_edges
mesh.bnd_labels = bnd_labels
DrawMesh(mesh, numbering=True)

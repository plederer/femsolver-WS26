# A mesh of the unit square consisting of two triangles, entered by hand.
#
#   3 -------- 2
#   |  1    /  |
#   |    /     |
#   |  /    0  |
#   0 -------- 1

from myfem import *

points = [(0, 0), (1, 0), (1, 1), (0, 1)]
elements = [(0, 1, 2),      # counterclockwise
            (0, 2, 3)]
bnd_edges = [(0, 1), (1, 2), (2, 3), (3, 0)]
bnd_labels = ["bottom", "right", "top", "left"]

mesh = Mesh(points, elements, bnd_edges, bnd_labels)
print(mesh)
print("edges on the boundary 'left':", mesh.GetBoundaryEdges("left"))

DrawMesh(mesh, numbering=True)

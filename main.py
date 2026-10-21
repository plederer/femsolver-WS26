from netgen.occ import Rectangle, OCCGeometry, X, Y

from femsolver.ngsinterface import ConvertNGSolveMesh
from femsolver.draw import DrawMesh

# geometry: unit square with named boundary parts
shape = Rectangle(1, 1).Face()
shape.edges.Min(Y).name = "bottom"
shape.edges.Max(X).name = "right"
shape.edges.Max(Y).name = "top"
shape.edges.Min(X).name = "left"

from femsolver import *

mesh = MakeStructuredMesh(4,4)

DrawMesh(mesh, numbering=True)
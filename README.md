# femsolver

## Setup

Requires Python 3 with `numpy` and `matplotlib`:

```
python3 -m pip install -r requirements.txt
```

[NGSolve](https://ngsolve.org) is only needed to import meshes generated with
Netgen/NGSolve (`python3 -m pip install ngsolve`).

## Provided files

| file                       | content                                                 |
|----------------------------|---------------------------------------------------------|
| `femsolver/draw.py`        | `DrawMesh` for plotting a mesh                          |
| `femsolver/ngsinterface.py`| `ConvertNGSolveMesh`, converts a Netgen/NGSolve mesh    |
| `example_ngsolve_mesh.py`  | generates a mesh with NGSolve, converts and draws it    |
| `test_mesh.py`             | tests for your `Mesh` and `MakeStructuredMesh`          |

All other files of the package `femsolver` are written by you, step by step,
following the exercise sheets.

## Usage

The example only uses the provided files and runs right away:

```
python3 example_ngsolve_mesh.py
```

Once you have written `femsolver/mesh.py`, a mesh is created and drawn with

```python
from femsolver import *

mesh = MakeStructuredMesh(4, 4)       # your own structured mesh
DrawMesh(mesh, numbering=True)

from netgen.occ import unit_square
mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))   # mesh from NGSolve
DrawMesh(mesh)
```

and checked with `python3 test_mesh.py`.

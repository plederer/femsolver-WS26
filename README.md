# femsolver

A finite element solver in Python 3, developed step by step in the exercises of
the lecture *Numerics of Partial Differential Equations* (Universität Hamburg).

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

All other files of the package `femsolver` are written by you, step by step,
following the exercise sheets.

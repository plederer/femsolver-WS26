# fem-solver

A finite element solver in Python 3, developed step by step in the exercises of
the lecture *Numerics of Partial Differential Equations* (Universität Hamburg).

## Setup

Requires Python 3 with `numpy` and `matplotlib`:

```
python3 -m pip install -r requirements.txt
```

## Structure

| file                 | content                                              |
|----------------------|------------------------------------------------------|
| `myfem/mesh.py`      | class `Mesh`, function `MakeStructuredMesh` (to do)  |
| `myfem/draw.py`      | `DrawMesh` for plotting meshes                       |
| `example_mesh.py`    | a small mesh entered by hand                         |
| `test_mesh.py`       | tests for `MakeStructuredMesh`                       |

Run the example and the tests with

```
python3 example_mesh.py
python3 test_mesh.py
```

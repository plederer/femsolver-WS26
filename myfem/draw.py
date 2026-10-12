import numpy as np
import matplotlib.pyplot as plt

__all__ = ["DrawMesh"]


def DrawMesh(mesh, numbering=False, boundary=True, ax=None, show=True):
    """Draws a triangular mesh with matplotlib.

    mesh      : Mesh
    numbering : if True, the vertex numbers (black) and the
                element numbers (blue) are shown
    boundary  : if True, the boundary edges are drawn in color,
                one color per boundary label
    ax        : matplotlib axes to draw into (a new figure if None)
    show      : if True, the figure is shown (plt.show())
    """
    if ax is None:
        fig, ax = plt.subplots()

    x, y = mesh.points[:, 0], mesh.points[:, 1]
    ax.triplot(x, y, mesh.elements, color="gray", linewidth=0.8)

    if boundary:
        colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
        for i, label in enumerate(dict.fromkeys(mesh.bnd_labels)):
            edges = mesh.GetBoundaryEdges(label)
            for k, (a, b) in enumerate(edges):
                ax.plot(x[[a, b]], y[[a, b]], color=colors[i % len(colors)],
                        linewidth=2.5, label=label if k == 0 else None)
        if mesh.nbnd > 0:
            ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5))

    if numbering:
        for i, (px, py) in enumerate(mesh.points):
            ax.annotate(str(i), (px, py), textcoords="offset points", xytext=(3, 3),
                        fontsize=8, color="black")
        for i, el in enumerate(mesh.elements):
            cx, cy = mesh.points[el].mean(axis=0)
            ax.text(cx, cy, str(i), fontsize=8, color="tab:blue",
                    ha="center", va="center")

    ax.set_aspect("equal")
    ax.set_title(str(mesh))
    if show:
        plt.show()
    return ax

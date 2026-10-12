import numpy as np
import matplotlib.pyplot as plt

__all__ = ["DrawMesh"]


def DrawMesh(mesh, numbering=False, boundary=True, ax=None, show=True):
    """Draws a triangular mesh with matplotlib.

    Only uses the arrays mesh.points, mesh.elements, mesh.bnd_edges
    and mesh.bnd_labels.

    mesh      : Mesh
    numbering : if True, the vertex numbers (black) and the
                element numbers (blue) are shown
    boundary  : if True, the boundary edges are drawn in color,
                one color per boundary label
    ax        : matplotlib axes to draw into (a new figure if None)
    show      : if True, the figure is shown (plt.show())
    """
    points = np.asarray(mesh.points, dtype=float)
    elements = np.asarray(mesh.elements, dtype=int)
    bnd_edges = np.asarray(mesh.bnd_edges, dtype=int).reshape(-1, 2)
    bnd_labels = list(mesh.bnd_labels)

    if ax is None:
        fig, ax = plt.subplots()

    x, y = points[:, 0], points[:, 1]
    ax.triplot(x, y, elements, color="gray", linewidth=0.8)

    if boundary and len(bnd_edges) > 0:
        colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
        for i, label in enumerate(dict.fromkeys(bnd_labels)):
            color = colors[i % len(colors)]
            edges = [e for e, l in zip(bnd_edges, bnd_labels) if l == label]
            for k, (a, b) in enumerate(edges):
                ax.plot(x[[a, b]], y[[a, b]], color=color, linewidth=2.5,
                        label=label if k == 0 else None)
        ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5))

    if numbering:
        for i, (px, py) in enumerate(points):
            ax.annotate(str(i), (px, py), textcoords="offset points", xytext=(3, 3),
                        fontsize=8, color="black")
        for i, el in enumerate(elements):
            cx, cy = points[el].mean(axis=0)
            ax.text(cx, cy, str(i), fontsize=8, color="tab:blue",
                    ha="center", va="center")

    ax.set_aspect("equal")
    ax.set_title(f"{len(points)} vertices, {len(elements)} elements, "
                 f"{len(bnd_edges)} boundary edges")
    if show:
        plt.show()
    return ax

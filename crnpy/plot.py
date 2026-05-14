"""Reaction graph visualization, phase portrait, and bifurcation diagrams."""
from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import matplotlib.axes
    import networkx as nx


def _complex_label(c) -> str:
    """Human-readable label for a Complex (frozenset of (name, coeff) pairs)."""
    parts = []
    for name, coeff in sorted(c, key=lambda x: x[0]):
        parts.append(f"{coeff}{name}" if coeff > 1 else name)
    return " + ".join(parts) if parts else "∅"


def _layout_components(
    G: nx.DiGraph,
    layout_fn,
    seed: int = 42,
    grid_spacing: float = 3.5,
) -> dict:
    """
    Lay out each weakly connected component independently, then arrange
    the components in a grid.  Prevents disconnected LCs from collapsing
    onto each other (spring layout has no inter-component repulsion).
    Components are sorted largest-first so they occupy the top-left cells.
    """
    import numpy as np
    import networkx as nx

    wccs = sorted(nx.weakly_connected_components(G), key=len, reverse=True)
    n = len(wccs)
    cols = max(1, int(np.ceil(np.sqrt(n))))

    pos: dict = {}
    for i, wcc in enumerate(wccs):
        row, col = divmod(i, cols)
        sub = G.subgraph(wcc)

        if len(wcc) == 1:
            sub_pos = {next(iter(wcc)): np.array([0.0, 0.0])}
        else:
            try:
                sub_pos = layout_fn(sub, seed=seed)
            except TypeError:
                sub_pos = layout_fn(sub)
            except Exception:
                sub_pos = nx.spring_layout(sub, seed=seed)

        # Normalise to [-1, 1] so each component occupies the same area
        coords = np.array(list(sub_pos.values()), dtype=float)
        if len(coords) > 1:
            coords -= coords.mean(axis=0)
            scale = np.max(np.abs(coords))
            if scale > 1e-9:
                coords /= scale

        center = np.array([col * grid_spacing, -row * grid_spacing])
        for node, coord in zip(sub_pos.keys(), coords):
            pos[node] = coord + center

    return pos


def draw_reaction_graph(
    G: nx.DiGraph,
    ax=None,
    layout: str = "spring",
    node_color: str = "#AED6F1",
    edge_color: str = "#2C3E50",
    node_size: int = 900,
    font_size: int = 8,
) -> matplotlib.axes.Axes:
    """
    Draw the complex reaction graph.

    Nodes are reaction complexes; directed edges are reactions.
    Disconnected linkage classes are arranged in a grid so they do not
    overlap.  Nodes, edges, and labels are drawn separately so that
    arrow heads are correctly placed at the node boundary and labels
    have a white background box for legibility.
    """
    import matplotlib.pyplot as plt
    import networkx as nx

    if ax is None:
        _, ax = plt.subplots(figsize=(10, 6))

    layouts = {
        "spring":   nx.spring_layout,
        "shell":    nx.shell_layout,
        "spectral": nx.spectral_layout,
        "planar":   nx.planar_layout,
        "kamada_kawai": nx.kamada_kawai_layout,
    }
    layout_fn = layouts.get(layout, nx.spring_layout)
    pos = _layout_components(G, layout_fn)

    # 1. Nodes
    nx.draw_networkx_nodes(
        G, pos, ax=ax,
        node_color=node_color,
        node_size=node_size,
        linewidths=0.8,
        edgecolors="#5D8AA8",
    )

    # 2. Edges — draw separately so node_size is used for correct arrowhead offset
    nx.draw_networkx_edges(
        G, pos, ax=ax,
        edge_color=edge_color,
        arrows=True,
        arrowsize=18,
        arrowstyle="-|>",
        node_size=node_size,
        connectionstyle="arc3,rad=0.18",
        width=1.4,
        min_source_margin=15,
        min_target_margin=15,
    )

    # 3. Labels — white bbox so they stay readable over edges/nodes
    labels = {node: _complex_label(node) for node in G.nodes()}
    nx.draw_networkx_labels(
        G, pos, labels, ax=ax,
        font_size=font_size,
        font_color="#1a1a2e",
        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.85),
    )

    ax.set_axis_off()
    ax.margins(0.15)
    return ax


def plot_bifurcation(
    result: BifurcationResult,
    species: str,
    ax=None,
    stable_color: str = "#2196F3",
    unstable_color: str = "#F44336",
) -> matplotlib.axes.Axes:
    """
    Bifurcation diagram: x = parameter value (log scale), y = steady-state [species].
    Solid lines = stable branches, dashed lines = unstable branches.

    Branches are identified by sorting steady states by concentration at each parameter
    value and connecting same-index entries across consecutive parameter points.
    Stability changes within a branch (e.g. Hopf bifurcation) are handled by splitting
    the branch into segments of consistent stability.
    """
    import matplotlib.pyplot as plt
    from .analysis import BifurcationResult

    if ax is None:
        _, ax = plt.subplots(figsize=(8, 5))

    # Determine the maximum number of simultaneous branches
    max_branches = max((len(ssl) for ssl in result.steady_states), default=0)
    if max_branches == 0:
        return ax

    # Collect (pval, conc, is_stable) per branch, sorted by concentration for
    # consistent branch ordering across parameter values.
    branch_data: list[list[tuple[float, float, bool]]] = [[] for _ in range(max_branches)]
    for pval, ss_list in zip(result.parameter_values, result.steady_states):
        sorted_ss = sorted(ss_list, key=lambda s: s.concentrations.get(species, 0.0))
        for i, ss in enumerate(sorted_ss):
            branch_data[i].append((
                float(pval),
                ss.concentrations.get(species, float("nan")),
                ss.is_stable,
            ))

    # Plot each branch as one or more line segments, splitting on stability changes.
    legend_added: set[str] = set()

    def _plot_segment(xs, ys, stable):
        color = stable_color if stable else unstable_color
        ls = "-" if stable else "--"
        label_str = "Stable" if stable else "Unstable"
        lbl = label_str if label_str not in legend_added else "_nolegend_"
        legend_added.add(label_str)
        ax.plot(xs, ys, ls, color=color, linewidth=2.0, label=lbl, zorder=3)

    for branch in branch_data:
        if not branch:
            continue
        cur_stable = branch[0][2]
        seg_x: list[float] = [branch[0][0]]
        seg_y: list[float] = [branch[0][1]]

        for pval, conc, is_stable in branch[1:]:
            if is_stable != cur_stable:
                # Close the current segment and start a new one at the transition
                seg_x.append(pval)
                seg_y.append(conc)
                _plot_segment(seg_x, seg_y, cur_stable)
                seg_x, seg_y = [pval], [conc]
                cur_stable = is_stable
            else:
                seg_x.append(pval)
                seg_y.append(conc)

        _plot_segment(seg_x, seg_y, cur_stable)

    ax.set_xscale("log")
    ax.set_xlabel(result.parameter_name)
    ax.set_ylabel(f"[{species}] (M)")
    ax.set_title(f"Bifurcation: {result.parameter_name}")
    ax.legend(framealpha=0.7)
    return ax


# Forward reference fix
try:
    from .analysis import BifurcationResult
except ImportError:
    pass

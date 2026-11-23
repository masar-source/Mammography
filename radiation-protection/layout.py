"""Compare virtual survey results at a fixed prospective workload."""

import csv
import math
from pathlib import Path
import numpy as np
from scipy.stats import t as student_t
from room_geometry import booth, reference, survey_points, screen_settings, validate_room_geometry


def csv_file(path, rows):
    if not rows:
        return
    with Path(path).open("w", encoding="utf-8", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def footprint(cfg, layout):
    _, centre, n, t = booth(cfg, layout)
    half = screen_settings(cfg, layout)["width_mm"] / 2
    return np.array(
        [
            centre - half * t,
            centre + half * t,
            centre + half * t + layout["return_mm"] * n,
        ]
    )


def fits(cfg, layout):
    h = np.array(cfg["room"]["inside_mm"][:2]) / 2
    if layout["screen_pb_mm"] > 0:
        vertices = footprint(cfg, layout)[:, :2]
        if np.any(abs(vertices) >= h) or screen_settings(cfg, layout)["height_mm"] > cfg["room"]["inside_mm"][2]:
            return False
    for p in survey_points(cfg, layout).values():
        if p["group"] == "booth" and np.any(
            np.abs(p["xyz_mm"][:2]) + np.array(p["size_mm"][:2]) / 2 >= h
        ):
            return False
    return True


def preview(cfg, layout_ids, out):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    for name in layout_ids:
        validate_room_geometry(cfg, cfg["layouts"][name])
    rm = cfg["room"]
    w, d, _ = np.array(rm["inside_mm"]) / 1000
    fig, ax = plt.subplots(figsize=(8, 6.2), layout="constrained")
    ax.add_patch(
        Rectangle((-w / 2, -d / 2), w, d, fill=False, edgecolor="#273d44", linewidth=2)
    )
    x0, wdoor = rm["door_centre_x_mm"] / 1000, rm["door_width_mm"] / 1000
    jamb = x0 - wdoor/2
    if rm.get("door_open"):
        ax.plot([jamb, jamb+wdoor], [-d/2]*2, color="white", lw=3)
        ax.plot([jamb+0.002]*2, [-d/2, -d/2-wdoor+0.004], color="#a15e2c", lw=3)
    else:
        ax.plot([jamb, jamb+wdoor], [-d/2]*2, color="#a15e2c", lw=4)
    ax.text(x0, -d/2+0.09, "Open door" if rm.get("door_open") else "Closed door",
            ha="center", fontsize=8)
    seen_rings = set()
    for i, id_ in enumerate(layout_ids):
        v = cfg["layouts"][id_]
        ref = reference(cfg, v)[:2] / 1000
        op = booth(cfg, v)[0] / 1000
        color = plt.get_cmap("tab10")(i % 10)
        sites = survey_points(cfg, v)
        if i == 0:
            door_xy = np.array([sites[k]["xyz_mm"][:2] for k in ("door_leaf", "door_jamb")])/1000
            ax.scatter(*door_xy.T, marker="x", color="#a15e2c", label="Door survey cells")
        trainee = np.asarray(sites["trainee_1700"]["xyz_mm"][:2])/1000
        ax.scatter(*trainee, marker="s", facecolor="none", edgecolor=color,
                   label="Second person" if i == 0 else None)
        ring = np.array([p["xyz_mm"][:2] for key, p in sites.items() if key.startswith("ring_")]) / 1000
        key = tuple(ring.round(8).ravel())
        if key not in seen_rings:
            first = not seen_rings
            ax.plot(*ref, marker="+", ms=12, mew=2, color="black",
                    label="Reference breast centre" if first else None)
            ax.scatter(*ring.T, facecolor="none", edgecolor="#34737c", s=23,
                       label="Scatter survey positions" if first else None)
            seen_rings.add(key)
        if v["screen_pb_mm"] > 0:
            pts = footprint(cfg, v) / 1000
            ax.plot(pts[:, 0], pts[:, 1], color=color, lw=2)
        ax.scatter(*op[:2], s=35, color=color, label=f"{id_}  {v['distance_m']:g} m")
    ax.set(
        xlabel="Room x (m)",
        ylabel="Room y (m)",
        aspect="equal",
        title="Proposed booth layouts and virtual survey positions",
    )
    ax.grid(alpha=0.15)
    ax.legend(loc="lower left", fontsize=8, framealpha=0.93, ncol=2)
    fig.savefig(out, dpi=155)
    plt.close(fig)

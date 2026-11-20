#!/usr/bin/env python3
"""
Show the scene that Geant4 itself exported (python revelation_gate.py --visu vrml_file_only -> geometry.wrl): the
geometry exactly as GATE/Geant4 built it, with the photon tracks. For builds without the Qt viewer.

Writes next to the .wrl file:
  geant4_scene.html   interactive 3D view in a browser (rotate, zoom; click legend entries to hide parts)
  geant4_scene.png    two still views: oblique, and side view (x-z) with the tracks from the focal spot
Usage: python view_vrml.py output/LE_phantom_PA12_MCT_visu_large/geometry.wrl [--no-html] [--no-png]
"""
import argparse, re
from pathlib import Path
import numpy as np

# display opacity per part (the VRML file carries Geant4's own transparency, which is set for the whole scene)
OPACITY = {"envelope": 0.10, "filling": 0.06, "plate": 0.35, "support": 0.5, "paddle": 0.25, "aperture": 0.4,
           "phsp_plane": 0.25}
HIDDEN = {"world", "absorber"}
LABEL = {"envelope": "container (PMMA)", "filling": "filling (liquid)", "plate": "plate (PMMA)", "disc": "iodine discs",
         "support": "breast support", "paddle": "paddle", "aperture": "field aperture (virtual)",
         "phsp_plane": "phase-space plane"}


def read_vrml(path):
    """shapes of a Geant4 VRML 2.0 file: solids (IndexedFaceSet) and tracks (IndexedLineSet), global coordinates (mm)"""
    txt = open(path).read()
    parts = re.split(r"#---------- (SOLID: |POLYLINE)", txt)
    shapes = []
    for sep, body in zip(parts[1::2], parts[2::2]):
        pm = re.search(r"point\s*\[(.*?)\]", body, re.S)
        if not pm:
            continue
        im = re.search(r"coordIndex\s*\[(.*?)\]", body, re.S)
        col = re.search(r"diffuseColor\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)", body)
        P = np.array(pm.group(1).replace(",", " ").split(), float).reshape(-1, 3)
        name = body.split("\n", 1)[0].strip().split(":")[0] if sep.startswith("SOLID") else "track"
        shapes.append(dict(name=name, kind="faces" if "IndexedFaceSet" in body else "line", points=P,
                           index=np.array(im.group(1).replace(",", " ").split(), int) if im else None,
                           color=tuple(float(c) for c in col.groups()) if col else (0.5, 0.5, 0.5)))
    return shapes


def group_of(name):
    name = re.sub(r"_rep_\d+$", "", name)              # repeated volumes: spheres_15p88_rep_12 -> spheres_15p88
    if name.startswith("disc"):
        return "disc"
    return name


def triangles(index):
    """fan triangulation of the polygons of an IndexedFaceSet (index lists separated by -1)"""
    tris, poly = [], []
    for k in index:
        if k == -1:
            tris += [(poly[0], poly[i], poly[i + 1]) for i in range(1, len(poly) - 1)]
            poly = []
        else:
            poly.append(int(k))
    return np.array(tris, int).reshape(-1, 3)


def merged(shapes):
    """{group: (vertices, triangles, colour)} with all solids of a group in one mesh"""
    groups = {}
    for s in shapes:
        if s["kind"] != "faces" or s["name"] in HIDDEN:
            continue
        g = group_of(s["name"])
        V, F, off = groups.get(g, ([], [], 0))[:3] if g in groups else ([], [], 0)
        V.append(s["points"]); F.append(triangles(s["index"]) + off)
        groups[g] = (V, F, off + len(s["points"]), s["color"])
    return {g: (np.concatenate(V), np.concatenate(F), c) for g, (V, F, _, c) in groups.items()}


def label(g, n=None):
    if g.startswith("spheres_"):
        return "spheres %s mm" % g.split("_", 1)[1].replace("p", ".") + ("" if n is None else " (%d)" % n)
    return LABEL.get(g, g)


def write_html(meshes, tracks, counts, out):
    import plotly.graph_objects as go
    data = []
    for g, (V, F, c) in meshes.items():
        rgb = "rgb(%d,%d,%d)" % tuple(int(255 * x) for x in c)
        data.append(go.Mesh3d(x=V[:, 0].round(3), y=V[:, 1].round(3), z=V[:, 2].round(3), i=F[:, 0], j=F[:, 1],
                              k=F[:, 2], color=rgb, opacity=OPACITY.get(g, 1.0), name=label(g, counts.get(g)),
                              showlegend=True, flatshading=True, hoverinfo="name"))
    lx, ly, lz = [], [], []
    for t in tracks:
        P = t["points"]
        lx += list(P[:, 0]) + [None]; ly += list(P[:, 1]) + [None]; lz += list(P[:, 2]) + [None]
    data.append(go.Scatter3d(x=lx, y=ly, z=lz, mode="lines", line=dict(color="#e6550d", width=2),
                             name="photon tracks (%d)" % len(tracks)))
    fig = go.Figure(data)
    fig.update_layout(title="Scene exported by Geant4 (GATE 10) - geometry as built, with photon tracks (mm)",
                      scene=dict(aspectmode="data", xaxis_title="x", yaxis_title="y", zaxis_title="z"),
                      legend=dict(itemsizing="constant"))
    fig.write_html(out, include_plotlyjs=True)


def write_png(meshes, tracks, out, sid=656.0):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection, LineCollection
    fig, ax = plt.subplots(1, 2, figsize=(17, 7.5), gridspec_kw=dict(width_ratios=[1.35, 1]))
    # (a) oblique view: rotate about z, then tilt; painter's algorithm (far triangles first)
    az, el = np.radians(-55), np.radians(28)
    Rz = np.array([[np.cos(az), -np.sin(az), 0], [np.sin(az), np.cos(az), 0], [0, 0, 1]])
    Rx = np.array([[1, 0, 0], [0, np.cos(el), -np.sin(el)], [0, np.sin(el), np.cos(el)]])
    M = Rx @ Rz
    polys, cols, depth = [], [], []
    for g, (V, F, c) in meshes.items():
        if g in ("aperture", "phsp_plane"):
            continue
        W = V @ M.T
        T = W[F]
        polys.append(T[:, :, [0, 2]]); depth.append(T[:, :, 1].mean(1))
        a = OPACITY.get(g, 1.0)
        shade = np.clip(0.55 + 0.45 * np.abs(np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0])[:, 1]) /
                        (np.linalg.norm(np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0]), axis=1) + 1e-12), 0, 1)
        cols.append(np.c_[np.outer(shade, c), np.full(len(F), a)])
    polys, cols, depth = np.concatenate(polys), np.concatenate(cols), np.concatenate(depth)
    o = np.argsort(-depth)
    ax[0].add_collection(PolyCollection(polys[o], facecolors=cols[o], edgecolors="none"))
    segs = [(t["points"] @ M.T)[:, [0, 2]] for t in tracks if len(t["points"]) > 1]
    ax[0].add_collection(LineCollection(segs, colors="#e6550d", linewidths=0.6))
    allp = np.concatenate([p.reshape(-1, 2) for p in [polys[o][::50]]])
    ax[0].set_xlim(allp[:, 0].min() - 20, allp[:, 0].max() + 20); ax[0].set_ylim(allp[:, 1].min() - 20, allp[:, 1].max() + 60)
    ax[0].set_aspect("equal"); ax[0].axis("off")
    ax[0].set_title("Geant4 scene: phantom as built (spheres >= 6.35 mm shown), tracks of %d photons" % len(segs))
    # (b) side view x-z with the tracks back to the focal spot
    for g, (V, F, c) in meshes.items():
        if g in ("support", "paddle", "envelope", "plate"):
            ax[1].fill_between([V[:, 0].min(), V[:, 0].max()], V[:, 2].min(), V[:, 2].max(), color=c, alpha=0.5, lw=0)
    ax[1].add_collection(LineCollection([t["points"][:, [0, 2]] for t in tracks if len(t["points"]) > 1],
                                        colors="#e6550d", linewidths=0.5))
    ax[1].plot(0, sid, "r*", ms=12)
    ax[1].set_xlim(-80, 280); ax[1].set_ylim(-110, sid + 30)
    ax[1].set_xlabel("x (mm)"); ax[1].set_ylabel("z (mm)"); ax[1].set_title("side view: tracks from the focal spot")
    fig.tight_layout(); fig.savefig(out, dpi=120); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("wrl")
    ap.add_argument("--no-html", action="store_true")
    ap.add_argument("--no-png", action="store_true")
    a = ap.parse_args()
    path = Path(a.wrl)
    shapes = read_vrml(path)
    tracks = [s for s in shapes if s["kind"] == "line"]
    counts = {}
    for s in shapes:
        if s["kind"] == "faces":
            counts[group_of(s["name"])] = counts.get(group_of(s["name"]), 0) + 1
    meshes = merged(shapes)
    print("read %s: %s; %d photon tracks" % (path.name, ", ".join("%s %d" % (label(g), n) for g, n in counts.items()
                                                                     if g not in HIDDEN), len(tracks)))
    if not a.no_html:
        write_html(meshes, tracks, counts, path.with_name("geant4_scene.html")); print("written", path.with_name("geant4_scene.html"))
    if not a.no_png:
        write_png(meshes, tracks, path.with_name("geant4_scene.png")); print("written", path.with_name("geant4_scene.png"))


if __name__ == "__main__":
    main()

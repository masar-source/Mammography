#!/usr/bin/env python3
"""
Look at the model before running GATE (no opengate needed). Built from the same files and constants that
revelation_gate.py gives to GATE: the STL meshes, the sphere packing and the disc layout of the phantom model
(--phantom L2v2 | L2, phantoms.py) and config.py.

Writes to output/preview/<case>/ (case = background for L2 v2, L2_<background> for the real L2):
  setup_side.png          side view at y = 0 (whole set-up to scale, and a zoom of phantom + detector stack)
  phantom_sections.png    phantom cut at y = 0, at x = 60 mm, a sphere layer seen from above, the plate plane
  projection.png          what the detector "sees": path length through spheres, PMMA and iodine per pixel
  phantom_3d.html         interactive 3D phantom: container, plate, discs, spheres (open in a browser)
  setup_3d.html           interactive 3D set-up: source, field, paddle, phantom, support, phase-space plane
To see the geometry as GATE actually builds it, use: python revelation_gate.py --visu vrml_file_only (or qt).

Usage: python preview.py [--phantom L2v2|L2] [--background PA12_MCT|PMMA_paraffin] [--no-3d] [--no-projection]
                         [--small-spheres markers|none] [--bin 4]
"""
import argparse, json
import numpy as np
from config import HERE, SYSTEM, SUPPORT, PADDLE, PHANTOM, OUT_DIR, BACKGROUNDS, add_case_args, apply_case_args
import phantoms
import sphere_packing as SP
import grid_detector as G
import revelation_gate as RG
from materials import iodine_rows

P = SP.P; DV = SP.D
H, R_FLAT, A_OUT = P["H"], P["R_flat"], P["a_out"]
SID = SYSTEM["sid_mm"]
Z = RG.geometry_z(); Z0 = Z["phantom"][0]
ROW_COLORS = ["#7f0000", "#b30000", "#d7301f", "#ef6548", "#fc8d59", "#fdbb84", "#fdd49e"]


# ----------------------------------------------------------------------------------------------- geometry helpers
def envelope_inside(x, y, z):
    if phantoms.is_l2():
        g = phantoms.l2()
        return (x >= 0) & (np.hypot(x, y) <= g["R_out"]) & (z >= 0) & (z <= g["H"])
    rho = np.where(x >= 0, np.hypot(x, y), np.abs(y))
    q, w = rho - R_FLAT, z - H / 2
    return (x >= DV["x_open"]) & (((q <= 0) & (np.abs(w) <= H / 2)) | ((q / A_OUT) ** 2 + (w / (H / 2)) ** 2 <= 1))


def plate_inside(x, y, z):
    if phantoms.is_l2():
        g = phantoms.l2()
        return (np.hypot(x, y) <= g["R_plate"]) & (x >= g["x_plate"]) & (np.abs(z - g["zc"]) <= g["t_plate"] / 2)
    disc = np.maximum(np.hypot(x, y) - SP.R_PLATE, SP.X_TAB_STOP - x)
    tab = np.maximum.reduce([np.abs(y) - SP.P_TAB, x - SP.X_TAB_STOP, SP.X_TAB_END - x])
    dh = np.maximum.reduce([np.minimum(disc, tab), SP.X_OPEN + 0.05 - x, np.abs(y) - SP.Y_CUT])
    return (dh <= 0) & (np.abs(z - SP.CI) <= SP.T_PLATE / 2)


def material_map(x, y, z):
    """0 outside, 1 container PMMA, 2 oil (free space), 3 plate (phantom coordinates; broadcastable arrays)"""
    x, y, z = np.broadcast_arrays(x, y, z)
    pts = np.c_[x.ravel(), y.ravel(), z.ravel()]
    env = envelope_inside(*pts.T)
    free = SP.free_distance(pts) > 0
    pl = plate_inside(*pts.T)
    m = np.zeros(len(pts), np.int8)
    m[env] = 1
    m[env & free] = 2
    m[pl] = 3
    return m.reshape(x.shape)


def load_spheres():
    f = phantoms.phantom_dir() / "spheres.npz"
    if not f.exists():
        print("note: %s not found (run sphere_packing.py --phantom %s); spheres are not shown" % (f, PHANTOM["model"]))
        return None
    S = np.load(f)
    return S["centres"], S["radii"], S["size_class"]


def load_discs():
    out, rows = [], iodine_rows()                                     # from data/L2v2_materials.db
    for d in phantoms.discs():
        io = rows[d["key"]]["mg_cm2"] if d["key"] in rows else 0.0
        out.append(dict(x=d["x"], y=d["y"], r=d["diam"] / 2, row=d["key"] if d["key"] != "large" else -1, iodine=io))
    return out


def section_circles(S, axis, value):
    """circles where the plane coord[axis] = value cuts the spheres: (centre in the other two coordinates, radius)"""
    if S is None:
        return np.zeros((0, 2)), np.zeros(0), np.zeros(0, int)
    C, R, cls = S
    d = np.abs(C[:, axis] - value)
    k = d < R
    other = [i for i in range(3) if i != axis]
    return C[k][:, other], np.sqrt(R[k] ** 2 - d[k] ** 2), cls[k]


# ----------------------------------------------------------------------------------------------- 2D figures
MAT_CMAP = None


def _mat_cmap():
    from matplotlib.colors import ListedColormap
    return ListedColormap(["white", "#9ecae1", "#fff2b2", "#bdbdbd"])


def draw_section(ax, h_axis, v_axis, h, v, fixed_axis, fixed_value, S, discs, z_world=False):
    """phantom section in the plane coord[fixed_axis] = fixed_value; h, v are 1D coordinate arrays (phantom
    coordinates; with z_world the section is drawn in world coordinates: z + z0 and x + the phantom's x offset)"""
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Circle, Rectangle
    dx = phantoms.x_offset() if (z_world and h_axis == 0) else 0.0
    zc, tp = phantoms.plate()
    h = np.asarray(h) - dx                                             # world x -> phantom x
    HH, VV = np.meshgrid(h, v)
    coords = [None] * 3
    coords[h_axis], coords[v_axis], coords[fixed_axis] = HH, VV, np.full(HH.shape, fixed_value)
    m = material_map(*coords)
    dz = Z0 if z_world else 0.0
    ext = [h[0] + dx, h[-1] + dx, v[0] + (dz if v_axis == 2 else 0), v[-1] + (dz if v_axis == 2 else 0)]
    ax.imshow(m, origin="lower", extent=ext, cmap=_mat_cmap(), vmin=0, vmax=3, interpolation="nearest")
    cen, rad, cls = section_circles(S, fixed_axis, fixed_value)
    if len(rad):
        hv = [h_axis, v_axis]
        other = [i for i in range(3) if i != fixed_axis]
        ch, cv = cen[:, other.index(hv[0])] + dx, cen[:, other.index(hv[1])] + (dz if v_axis == 2 else 0)
        ax.add_collection(PatchCollection([Circle((a, b), r) for a, b, r in zip(ch, cv, rad)],
                                          facecolor="#f08c2a", edgecolor="#8c4a0a", linewidth=0.2))
    # iodine discs cut by the plane (discs are vertical cylinders through the plate)
    if fixed_axis in (0, 1):
        for d in discs:
            off = abs((d["x"] if fixed_axis == 0 else d["y"]) - fixed_value)
            if off < d["r"] and d["iodine"] > 0:
                half = np.sqrt(d["r"] ** 2 - off ** 2)
                c = d["y"] if fixed_axis == 0 else d["x"] + dx
                ax.add_patch(Rectangle((c - half, zc - tp / 2 + dz), 2 * half, tp,
                                       facecolor=ROW_COLORS[d["row"]], edgecolor="k", linewidth=0.3))
    ax.set_aspect("equal")


def fig_setup_side(S, discs, out):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    fig = plt.figure(figsize=(15, 8.5))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 3.2])
    ax0, ax1 = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    th_ph, _, _ = RG.source_cone("phantom")
    th_f, _, _ = RG.source_cone("field")
    # (a) whole set-up to scale
    ax = ax0
    ax.plot(0, SID, "r*", ms=12, label="focal spot (0.3 mm)")
    for xe in RG.FIELD_X:
        ax.plot([0, xe], [SID, 0], "g-", lw=0.8)
    ax.plot([0, np.tan(np.radians(th_ph)) * SID], [SID, 0], "m--", lw=0.8, label="phantom-target cone")
    k = (SID - RG.APERTURE_Z) / SID
    ax.plot([RG.FIELD_X[0] * k, RG.FIELD_X[1] * k], [RG.APERTURE_Z] * 2, "g-", lw=3, label="field aperture (virtual)")
    ax.add_patch(Rectangle((RG.PADDLE_X[0], Z["paddle"][0]), np.ptp(RG.PADDLE_X), PADDLE["mm"], color="#6baed6"))
    x0p, x1p = (phantoms.x_offset(), phantoms.x_offset() + phantoms.l2()["R_out"]) if phantoms.is_l2() else (DV["x_open"], DV["R_out"])
    ax.add_patch(Rectangle((x0p, Z0), x1p - x0p, H, color="#fdae6b", label="phantom"))
    ax.add_patch(Rectangle((RG.SUPPORT_X[0], Z["support"][0]), np.ptp(RG.SUPPORT_X), SUPPORT["mm"], color="0.4"))
    ax.add_patch(Rectangle((SYSTEM["missed_tissue_mm"], -2), SYSTEM["detector_active_mm"][1], 2, color="k", label="detector"))
    ax.set_xlim(-60, 300); ax.set_ylim(-10, SID + 20); ax.set_aspect("equal")
    ax.set_xlabel("x (mm)"); ax.set_ylabel("z (mm)"); ax.set_title("(a) whole set-up, y = 0, to scale")
    ax.legend(fontsize=7, loc="upper right")
    # (b) zoom: phantom section + below-plane stack
    ax = ax1
    draw_section(ax, 0, 2, np.arange(-30, 140, 0.1), np.arange(-0.5, H + 0.5, 0.1), 1, 0.0, S, discs, z_world=True)
    ax.add_patch(Rectangle((RG.PADDLE_X[0], Z["paddle"][0]), np.ptp(RG.PADDLE_X), PADDLE["mm"], color="#6baed6", label="paddle %.0f mm %s [A]" % (PADDLE["mm"], PADDLE["material"])))
    ax.add_patch(Rectangle((RG.SUPPORT_X[0], Z["support"][0]), np.ptp(RG.SUPPORT_X), SUPPORT["mm"], color="0.45", label="support %.1f mm %s [A], 17 mm above detector [V]" % (SUPPORT["mm"], SUPPORT["material"])))
    ax.add_patch(Rectangle((RG.PLANE_X[0], RG.PLANE_Z[0]), np.ptp(RG.PLANE_X), np.ptp(RG.PLANE_Z), color="#c51b8a", alpha=0.6, label="phase-space plane (GATE stops here)"))
    gz0 = G.Z_GRID_BOT - G.c
    ax.add_patch(Rectangle((-30, gz0), 290, G.GRID_TOP - gz0, facecolor="none", edgecolor="k", hatch="|||", lw=0.5, label="grid (applied after GATE) [A]"))
    ax.add_patch(Rectangle((-30, G.DET_COVER[0]), 290, G.DET_COVER[1] - G.DET_COVER[0], color="#a1d99b", label="detector cover [A]"))
    ax.add_patch(Rectangle((SYSTEM["missed_tissue_mm"], -SYSTEM["a_se_um"] / 1000), SYSTEM["detector_active_mm"][1], SYSTEM["a_se_um"] / 1000, color="k", label="a-Se (image starts at x = %.0f mm [A])" % SYSTEM["missed_tissue_mm"]))
    for xd in (3, 60, 120, 140):                                        # primary rays to a few detector points
        zt = 80
        ax.plot([xd * (SID - zt) / SID, xd], [zt, 0], "r-", lw=0.4, alpha=0.6)
    ax.axvline(0, color="g", lw=0.8, ls=":"); ax.text(0.5, 78, "field edge / support face x = 0", fontsize=7, color="g")
    ax.set_xlim(-30, 160); ax.set_ylim(-1.5, 80)
    ax.set_xlabel("x (mm), chest wall -> nipple"); ax.set_ylabel("z (mm)")
    ax.set_title("(b) phantom and detector stack at y = 0: PMMA blue, %s yellow, plate grey, %s spheres orange, iodine red"
                 % (PHANTOM["liquid_material"], PHANTOM["sphere_material"]), fontsize=10)
    ax.legend(fontsize=7, loc="upper right", framealpha=0.9)
    fig.tight_layout(); fig.savefig(out / "setup_side.png", dpi=130); plt.close(fig)


def fig_sections(S, discs, out):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle
    fig, ax = plt.subplots(2, 2, figsize=(15, 13), gridspec_kw=dict(height_ratios=[1, 2.4]))
    draw_section(ax[0, 0], 0, 2, np.arange(-25, 132, 0.08), np.arange(-0.5, H + 0.5, 0.08), 1, 0.0, S, discs)
    ax[0, 0].set_title("cut at y = 0 (x-z)"); ax[0, 0].set_xlabel("x (mm)"); ax[0, 0].set_ylabel("z in phantom (mm)")
    x_cut = 60.4
    draw_section(ax[0, 1], 1, 2, np.arange(-132, 132, 0.1), np.arange(-0.5, H + 0.5, 0.1), 0, x_cut, S, discs)
    ax[0, 1].set_title("cut at x = %.1f mm (y-z), through disc column 4" % x_cut); ax[0, 1].set_xlabel("y (mm)")
    zc, tp = phantoms.plate()
    z_layer = (zc - tp / 2 + (phantoms.l2()["zb"] if phantoms.is_l2() else P["t_bot"])) / 2
    draw_section(ax[1, 0], 0, 1, np.arange(-25, 132, 0.15), np.arange(-132, 132, 0.15), 2, z_layer, S, discs)
    ax[1, 0].set_title("seen from above: cut at z = %.1f mm (middle of the lower sphere layer)" % z_layer)
    ax[1, 0].set_xlabel("x (mm)"); ax[1, 0].set_ylabel("y (mm)")
    a = ax[1, 1]
    draw_section(a, 0, 1, np.arange(-25, 132, 0.15), np.arange(-132, 132, 0.15), 2, zc, None, discs)
    for d in discs:
        a.add_patch(Circle((d["x"], d["y"]), d["r"], facecolor=ROW_COLORS[d["row"]] if d["iodine"] > 0 else "white",
                           edgecolor="k", lw=0.5))
    for row, io in ((k, v["mg_cm2"]) for k, v in iodine_rows().items() if k != "large"):
        ys = [d["y"] for d in discs if d["row"] == row]
        if ys:
            a.text(-3, np.mean(ys), "%.2f" % io, ha="right", va="center", fontsize=7)
    a.text(-3, 112, "mg I/cm$^2$ [A]", ha="right", fontsize=7)
    if not phantoms.is_l2():
        a.add_patch(Circle((0, 0), R_FLAT, fill=False, ls="--", lw=0.7, ec="b"))
        a.text(R_FLAT * 0.72, R_FLAT * 0.72, "end of flat region", fontsize=7, color="b")
    xi = SYSTEM["missed_tissue_mm"] - phantoms.x_offset()                # image start in phantom coordinates
    a.axvline(xi, color="k", lw=0.7, ls=":")
    a.text(xi + 1, -128, "image starts (missed tissue [A])", fontsize=7)
    a.set_title("plate plane z = %.2f mm: plate, iodine discs (9 mm discs white: no iodine [A])" % zc)
    a.set_xlabel("x (mm)")
    fig.suptitle("Phantom %s (%s), phantom coordinates; background %s: %s"
                 % (PHANTOM["model"], phantoms.model()["description"], PHANTOM["background"],
                    BACKGROUNDS[PHANTOM["background"]]["description"]))
    fig.tight_layout(); fig.savefig(out / "phantom_sections.png", dpi=120); plt.close(fig)


def fig_projection(binning, have_spheres, out):
    import matplotlib.pyplot as plt
    import projector as PJ
    maps, grid = PJ.path_maps(binning=binning, spheres=have_spheres, discs=True)
    x, y = grid["x"], grid["y"]
    ext = [x[0], x[-1], y[0], y[-1]]                                     # y up, as in phantom_sections.png
    panels = [("spheres", "path length through spheres (mm)", "viridis"), ("PMMA", "path length through PMMA (mm)", "Blues"),
              ("iodine_g_cm2", "iodine (mg/cm$^2$)", "Reds")]
    fig, ax = plt.subplots(1, 3, figsize=(16, 6.4))
    for a, (k, title, cmap) in zip(ax, panels):
        v = maps[k] * (1000 if k == "iodine_g_cm2" else 1)
        im = a.imshow(v, extent=ext, origin="lower", cmap=cmap, interpolation="nearest")
        fig.colorbar(im, ax=a, shrink=0.8); a.set_title(title, fontsize=10)
        a.set_xlabel("x on detector (mm)"); a.set_ylabel("y on detector (mm)")
    fig.suptitle("Projection onto the detector (image area %d x %d px at %.3f mm), geometry only - no spectrum needed"
                 % (len(y), len(x), grid["pixel_mm"]), fontsize=11)
    fig.tight_layout(); fig.savefig(out / "projection.png", dpi=110); plt.close(fig)


# ----------------------------------------------------------------------------------------------- 3D (plotly)
def mesh_from_triangles(T):
    V, inv = np.unique(np.round(T.reshape(-1, 3), 4), axis=0, return_inverse=True)
    F = inv.reshape(-1, 3)
    return V, F


def icosphere(subdiv=1):
    t = (1 + 5 ** 0.5) / 2
    V = np.array([[-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0], [0, -1, t], [0, 1, t], [0, -1, -t], [0, 1, -t],
                  [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1]], float)
    F = [[0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11], [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6],
         [7, 1, 8], [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9], [4, 9, 5], [2, 4, 11], [6, 2, 10],
         [8, 6, 7], [9, 8, 1]]
    V = list(V / np.linalg.norm(V, axis=1, keepdims=True))
    for _ in range(subdiv):
        cache, F2 = {}, []
        def mid(i, j):
            key = (min(i, j), max(i, j))
            if key not in cache:
                m = (V[i] + V[j]) / 2; V.append(m / np.linalg.norm(m)); cache[key] = len(V) - 1
            return cache[key]
        for a, b, c in F:
            ab, bc, ca = mid(a, b), mid(b, c), mid(c, a)
            F2 += [[a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]]
        F = F2
    return np.array(V), np.array(F)


def box_mesh(xr, yr, zr):
    V = np.array([[x, y, z] for z in zr for y in yr for x in xr], float)
    F = [[0, 1, 3], [0, 3, 2], [4, 6, 7], [4, 7, 5], [0, 4, 5], [0, 5, 1], [2, 3, 7], [2, 7, 6], [0, 2, 6], [0, 6, 4],
         [1, 5, 7], [1, 7, 3]]
    return V, np.array(F)


def cylinder_mesh(xc, yc, r, z0, z1, n=24):
    a = np.linspace(0, 2 * np.pi, n, endpoint=False)
    ring = np.c_[xc + r * np.cos(a), yc + r * np.sin(a)]
    V = np.r_[np.c_[ring, np.full(n, z0)], np.c_[ring, np.full(n, z1)], [[xc, yc, z0], [xc, yc, z1]]]
    F = []
    for i in range(n):
        j = (i + 1) % n
        F += [[i, j, n + j], [i, n + j, n + i], [2 * n, j, i], [2 * n + 1, n + i, n + j]]
    return V, np.array(F)


def merge(meshes):
    Vs, Fs, off = [], [], 0
    for V, F in meshes:
        Vs.append(V); Fs.append(F + off); off += len(V)
    return np.concatenate(Vs), np.concatenate(Fs)


def mesh3d(go, V, F, name, color, opacity, dz=0.0, dx=0.0, **kw):
    return go.Mesh3d(x=V[:, 0] + dx, y=V[:, 1], z=V[:, 2] + dz, i=F[:, 0], j=F[:, 1], k=F[:, 2], name=name, color=color,
                     opacity=opacity, showlegend=True, flatshading=True, hoverinfo="name", **kw)


def phantom_traces(go, S, discs, dz=0.0, small="markers", with_spheres=True, dx=0.0):
    tr = []
    zc, tp = phantoms.plate()
    for name, color, op in (("envelope", "#6baed6", 0.12), ("plate", "#969696", 0.35)):
        V, F = mesh_from_triangles(SP.load_stl(phantoms.phantom_dir() / f"{name}.stl"))
        tr.append(mesh3d(go, V, F, "container (PMMA)" if name == "envelope" else "plate (PMMA)", color, op, dz, dx))
    for row in sorted({d["row"] for d in discs}):
        ds = [d for d in discs if d["row"] == row and d["iodine"] > 0]
        if ds:
            V, F = merge([cylinder_mesh(d["x"], d["y"], d["r"], zc - tp / 2, zc + tp / 2) for d in ds])
            tr.append(mesh3d(go, V, F, "iodine row %d (%.2f mg/cm2)" % (row, ds[0]["iodine"]), ROW_COLORS[row], 1.0, dz, dx))
    if S is not None and with_spheres:
        C, R, cls = S
        Vi, Fi = icosphere(1)
        diam = PHANTOM["sphere_diameters_mm"]
        for k, dmm in enumerate(diam):
            c = C[cls == k]
            if len(c) == 0:
                continue
            if dmm >= 6.0:
                V = (Vi[None] * dmm / 2 + c[:, None, :]).reshape(-1, 3)
                F = (Fi[None] + (np.arange(len(c)) * len(Vi))[:, None, None]).reshape(-1, 3)
                tr.append(mesh3d(go, V, F, "%s spheres %.2f mm (%d)" % (PHANTOM["sphere_material"], dmm, len(c)), "#f08c2a", 1.0, dz, dx))
            elif small == "markers":
                tr.append(go.Scatter3d(x=c[:, 0] + dx, y=c[:, 1], z=c[:, 2] + dz, mode="markers", name="spheres %.2f mm (%d, markers)" % (dmm, len(c)),
                                       marker=dict(size=2.2 if dmm > 2 else 1.2, color="#d9711c"), hoverinfo="name"))
    return tr


def write_3d(S, discs, out, small):
    import plotly.graph_objects as go
    # (1) phantom
    fig = go.Figure(phantom_traces(go, S, discs, small=small))
    fig.update_layout(title="%s phantom as given to GATE (phantom coordinates, mm), background %s: %s. "
                            "Click legend entries to hide parts." % (PHANTOM["model"], PHANTOM["background"],
                                                                     BACKGROUNDS[PHANTOM["background"]]["description"]),
                      scene=dict(aspectmode="data", xaxis_title="x (chest wall -> nipple)", yaxis_title="y", zaxis_title="z"),
                      legend=dict(itemsizing="constant"))
    fig.write_html(out / "phantom_3d.html", include_plotlyjs=True)
    # (2) whole set-up (spheres left out; world coordinates)
    tr = phantom_traces(go, S, discs, dz=Z0, with_spheres=False, dx=phantoms.x_offset())
    for name, xr, yr, zr, col, op in (
            ("paddle", RG.PADDLE_X, RG.PADDLE_Y, Z["paddle"], "#9ecae1", 0.3),
            ("breast support", RG.SUPPORT_X, RG.SUPPORT_Y, Z["support"], "#636363", 0.6),
            ("phase-space plane", RG.PLANE_X, RG.PLANE_Y, RG.PLANE_Z, "#c51b8a", 0.25),
            ("detector (a-Se)", (SYSTEM["missed_tissue_mm"], SYSTEM["missed_tissue_mm"] + SYSTEM["detector_active_mm"][1]),
             (-SYSTEM["detector_active_mm"][0] / 2, SYSTEM["detector_active_mm"][0] / 2), (-0.2, 0.0), "#000000", 0.8)):
        V, F = box_mesh(xr, yr, zr)
        tr.append(mesh3d(go, V, F, name, col, op))
    k = (SID - RG.APERTURE_Z) / SID
    V, F = box_mesh((RG.FIELD_X[0] * k, RG.FIELD_X[1] * k), (RG.FIELD_Y[0] * k, RG.FIELD_Y[1] * k), (RG.APERTURE_Z - 0.05, RG.APERTURE_Z + 0.05))
    tr.append(mesh3d(go, V, F, "field aperture (virtual, z = %.0f)" % RG.APERTURE_Z, "#31a354", 0.5))
    lx, ly, lz = [], [], []
    for xe in RG.FIELD_X:
        for ye in RG.FIELD_Y:
            lx += [0, xe, None]; ly += [0, ye, None]; lz += [SID, 0, None]
    tr.append(go.Scatter3d(x=lx, y=ly, z=lz, mode="lines", line=dict(color="#31a354", width=2), name="collimated field edges"))
    th, _, _ = RG.source_cone("phantom")
    a = np.radians(np.linspace(90, 270, 61)); rr = np.tan(np.radians(th)) * (SID - Z0)
    tr.append(go.Scatter3d(x=list(-rr * np.cos(a)) + [0, -rr * np.cos(a[0])], y=list(-rr * np.sin(a)) + [0, -rr * np.sin(a[0])],
                           z=[Z0] * 61 + [SID, Z0], mode="lines", line=dict(color="#dd3497", width=2, dash="dash"),
                           name="emission cone (--target phantom), at the phantom bottom"))
    tr.append(go.Scatter3d(x=[0], y=[0], z=[SID], mode="markers", marker=dict(size=5, color="red"), name="focal spot (z = %.0f)" % SID))
    fig = go.Figure(tr)
    fig.update_layout(title="Revelation set-up in GATE world coordinates (mm). a-Se at z = 0, focal spot at z = %.0f." % SID,
                      scene=dict(aspectmode="data", xaxis_title="x", yaxis_title="y", zaxis_title="z"))
    fig.write_html(out / "setup_3d.html", include_plotlyjs=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-3d", action="store_true")
    ap.add_argument("--no-projection", action="store_true")
    ap.add_argument("--small-spheres", choices=["markers", "none"], default="markers",
                    help="3D view: show the 3.18 and 1.58 mm spheres as markers (default) or leave them out")
    ap.add_argument("--bin", type=int, default=4, help="pixel binning of the projection preview")
    add_case_args(ap, measured=False)
    a = ap.parse_args()
    import matplotlib; matplotlib.use("Agg")
    global H
    bg = apply_case_args(a)
    SP.use_model(PHANTOM["model"])
    H = phantoms.height()
    out = OUT_DIR / "preview" / bg; out.mkdir(parents=True, exist_ok=True)
    S, discs = load_spheres(), load_discs()
    if S is not None:
        C, R, cls = S
        print("spheres: %d (%s)" % (len(R), ", ".join("%g mm: %d" % (d, (cls == k).sum()) for k, d in enumerate(PHANTOM["sphere_diameters_mm"]))))
    fig_setup_side(S, discs, out); print("written", out / "setup_side.png")
    fig_sections(S, discs, out); print("written", out / "phantom_sections.png")
    if not a.no_projection:
        fig_projection(a.bin, S is not None, out); print("written", out / "projection.png")
    if not a.no_3d:
        write_3d(S, discs, out, a.small_spheres); print("written", out / "phantom_3d.html", out / "setup_3d.html")


if __name__ == "__main__":
    main()

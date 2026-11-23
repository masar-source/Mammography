#!/usr/bin/env python3
"""
Phantom models (config.PHANTOM_MODELS; chosen with --phantom in every script, or with REVELATION_PHANTOM):
  L2v2  the suggested redesign: revolved container with an elliptical rim, CadQuery meshes in phantom/; its geometry
        is in sphere_packing.py and projector.py (read from L2_container_v2/params.json)
  L2    the real L2 container with the new (July 2026) plate, defined here

Real L2, phantom coordinates (mm): x from the outer chest-wall (front) face towards the nipple, y along the chest
wall from the axis of the half-cylinder, z from the bottom outer surface. The axis lies in the front face (centre of
curvature of the outline on the front face in 0035 [M]); in the world the front face is at x = front_face_x_mm.
  envelope (PMMA)   rho <= R_out,               x >= 0,               0 <= z <= H
  filling (liquid)  rho <= R_in = R_out - t_wall, x >= t_front,       t_bot <= z <= H - t_top
  plate (PMMA)      rho <= R_in - clearance,     x >= t_front + 0.05,  |z - zc| <= t_plate / 2     (outline [A])
  discs             vertical cylinders through the plate at the positions measured in the July images
                    (L2_container_v2/disc_layout_measured.json), back-projected from the image to the plate plane
The three solids are prisms of a circular segment, so their GATE meshes are written here directly (inscribed polygons,
deviation from the exact surface <= 0.005 mm; no CadQuery needed):  python phantom_meshes.py --phantom L2
Usage (summary of the selected model):  python phantoms.py [--phantom L2]
"""
import argparse, csv, json, struct
import numpy as np
from config import HERE, SYSTEM, PHANTOM, PHANTOM_MODELS, set_phantom


def model():
    return PHANTOM_MODELS[PHANTOM["model"]]


def is_l2():
    return model()["geometry"] == "halfcylinder"


def phantom_dir():
    """meshes (envelope/filling/plate.stl, meshes.json) and sphere packing (spheres.npz) of the selected model"""
    return HERE / model()["dir"]


def l2():
    """derived dimensions of the real L2 (mm, phantom coordinates)"""
    m = PHANTOM_MODELS["L2"]
    R_in = m["R_out"] - m["t_wall"]
    zb, zt = m["t_bot"], m["H"] - m["t_top"]
    zc = m["plate_zc"] if m["plate_zc"] is not None else (zb + zt) / 2
    return dict(R_out=m["R_out"], H=m["H"], R_in=R_in, x_in=m["t_front"], zb=zb, zt=zt, zc=zc, t_plate=m["t_plate"],
                R_plate=R_in - m["plate_clearance"], x_plate=m["t_front"] + 0.05, x_f=m["front_face_x_mm"])


def _v2():
    return json.load(open(PHANTOM_MODELS["L2v2"]["params_json"]))


def height():
    return l2()["H"] if is_l2() else _v2()["H"]


def x_offset():
    """world x of the phantom origin (L2: outer front face; L2 v2: phantom coordinates are world coordinates)"""
    return l2()["x_f"] if is_l2() else 0.0


def plate():
    """(centre z, thickness) of the iodine plate in phantom coordinates (mm)"""
    if is_l2():
        g = l2()
        return g["zc"], g["t_plate"]
    P = _v2()
    return P["derived"]["ci"], P["t_plate"]


def reach_mm():
    """largest horizontal distance of the phantom from the focal-spot axis (x = y = 0) in the world (mm)"""
    if is_l2():
        g = l2()
        return g["x_f"] + g["R_out"]
    return _v2()["derived"]["R_out"]


def z0():
    """world z of the phantom bottom (on the breast support), as revelation_gate.geometry_z()"""
    return SYSTEM["support_to_detector_mm"]


# --------------------------------------------------------------------------------------------- real L2: free space
def quadrant_sdf(u, v):
    """signed distance to the region {u <= 0, v <= 0} (positive outside)"""
    return np.where((u > 0) | (v > 0), np.hypot(np.maximum(u, 0), np.maximum(v, 0)), np.maximum(u, v))


def l2_free_distance(Pt):
    """signed distance (mm) to the nearest wall or plate surface of the real L2, positive in the free space (exact
    inside the cavity; for the plate a lower bound near its rim, which only makes the packing more cautious)"""
    g = l2()
    x, y, z = Pt[:, 0], Pt[:, 1], Pt[:, 2]
    rho = np.hypot(x, y)
    d_cav = np.minimum.reduce([g["R_in"] - rho, x - g["x_in"], z - g["zb"], g["zt"] - z])
    dh = np.maximum(rho - g["R_plate"], g["x_plate"] - x)
    d_pl = quadrant_sdf(dh, np.abs(z - g["zc"]) - g["t_plate"] / 2)
    return np.minimum(d_cav, d_pl)


def l2_packing_box():
    g = l2()
    return np.array([g["x_in"], -g["R_in"], g["zb"]]), np.array([g["R_in"], g["R_in"], g["zt"]])


# ----------------------------------------------------------------------------------------- real L2: projector chords
def halfcyl_chord(xd, yd, Lray, z0w, xf, R, xa, zb, zt):
    """length (mm) of the ray from the focal spot (0, 0, SID) to the detector point (xd, yd, 0) inside
    {(x - xf)^2 + y^2 <= R^2, x - xf >= xa, zb <= z - z0w <= zt} (world x, y; phantom z). Ray point at fraction t:
    (t xd, t yd, SID (1 - t)); xd > 0 everywhere in the image."""
    sid = SYSTEM["sid_mm"]
    A = xd ** 2 + yd ** 2
    B = -2 * xd * xf
    C = xf ** 2 - R ** 2
    disc = B * B - 4 * A * C
    sq = np.sqrt(np.clip(disc, 0, None))
    ta, tb = (-B - sq) / (2 * A), (-B + sq) / (2 * A)
    lo = np.maximum.reduce([ta, (xf + xa) / xd, np.full_like(ta, 1 - (z0w + zt) / sid)])
    hi = np.minimum(tb, 1 - (z0w + zb) / sid)
    return np.where(disc > 0, np.clip(hi - lo, 0, None), 0.0) * Lray


def l2_chords(xd, yd, Lray, z0w):
    """path lengths (mm) through the envelope, the cavity and the plate of the real L2"""
    g = l2()
    env = halfcyl_chord(xd, yd, Lray, z0w, g["x_f"], g["R_out"], 0.0, 0.0, g["H"])
    cav = halfcyl_chord(xd, yd, Lray, z0w, g["x_f"], g["R_in"], g["x_in"], g["zb"], g["zt"])
    pl = halfcyl_chord(xd, yd, Lray, z0w, g["x_f"], g["R_plate"], g["x_plate"],
                       g["zc"] - g["t_plate"] / 2, g["zc"] + g["t_plate"] / 2)
    return env, cav, pl


# ---------------------------------------------------------------------------------------------------- iodine discs
def discs():
    """iodine discs of the selected model, phantom coordinates (mm): list of dict(name, kind, key, diam, x, y);
    key = row number (array) or 'large' (9 mm discs), as used by materials.iodine_rows()"""
    out = []
    if not is_l2():
        for r in csv.DictReader(open(model()["disc_layout"])):
            out.append(dict(name=f"disc_{r['kind']}_{r['row']}_{r['col']}", kind=r["kind"],
                            key=int(r["row"]) if r["kind"] == "array" else "large", diam=float(r["diameter_mm"]),
                            x=float(r["x_mm_from_support_chest_wall_face"]), y=float(r["y_mm"])))
        return out
    m, g = model(), l2()
    ref = m["layout_ref"]
    sid, p = SYSTEM["sid_mm"], SYSTEM["pixel_mm"]
    m_plate = sid / (sid - z0() - g["zc"])
    for k, d in enumerate(json.load(open(m["disc_layout"]))["discs"]):
        # measured image position (detector, world x from the focal-spot axis) -> plate plane
        x_det = SYSTEM["missed_tissue_mm"] + ref["col_px"] * p + d["x"] * ref["magnification"]
        y_det = d["y"] * ref["magnification"]
        arr = d["kind"] == "array"
        out.append(dict(name=f"disc_{d['kind']}_{d['row']}_{d['col']}" if arr else f"disc_{d['kind']}",
                        kind=d["kind"], key=int(d["row"]) if arr else "large", diam=float(d["diam"]),
                        x=x_det / m_plate - g["x_f"], y=y_det / m_plate))
    return out


# ------------------------------------------------------------------------------------------------ real L2: meshes
def segment_polygon(R, xa, tol=0.005):
    """vertices (counter-clockwise) of the polygon inscribed in the circular segment {rho <= R, x >= xa}, sagitta of
    the arc facets <= tol (mm)"""
    al = np.arccos(xa / R)
    n = int(np.ceil(2 * al / (2 * np.arccos(1 - tol / R))))
    ang = np.linspace(-al, al, n + 1)
    P = np.c_[R * np.cos(ang), R * np.sin(ang)]
    P[0, 0] = P[-1, 0] = xa                             # end points exactly on the flat face
    return P


def prism_triangles(poly, za, zb):
    """closed, outward-oriented triangle mesh of a convex polygon (CCW) extruded from za to zb: (n, 3, 3)"""
    c = poly.mean(0)
    T = []
    for i in range(len(poly)):
        a, b = poly[i], poly[(i + 1) % len(poly)]
        T.append([[c[0], c[1], zb], [a[0], a[1], zb], [b[0], b[1], zb]])          # top, normal +z
        T.append([[c[0], c[1], za], [b[0], b[1], za], [a[0], a[1], za]])          # bottom, normal -z
        T.append([[a[0], a[1], za], [b[0], b[1], za], [b[0], b[1], zb]])          # side, normal outwards
        T.append([[a[0], a[1], za], [b[0], b[1], zb], [a[0], a[1], zb]])
    return np.array(T, float)


def write_stl(path, T, name="solid"):
    """binary STL with facet normals"""
    n = np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0])
    n /= np.linalg.norm(n, axis=1, keepdims=True)
    with open(path, "wb") as f:
        f.write(name.encode()[:80].ljust(80, b" "))
        f.write(struct.pack("<I", len(T)))
        rec = np.zeros(len(T), dtype=[("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")])
        rec["n"], rec["v"] = n, T
        f.write(rec.tobytes())


def write_l2_meshes(out, tol=0.005):
    """envelope, filling and plate of the real L2 as STL (phantom coordinates) + meshes.json"""
    from phantom_meshes import mesh_check
    g = l2()
    out.mkdir(parents=True, exist_ok=True)
    solids = dict(envelope=(g["R_out"], 0.0, 0.0, g["H"]),
                  filling=(g["R_in"], g["x_in"], g["zb"], g["zt"]),
                  plate=(g["R_plate"], g["x_plate"], g["zc"] - g["t_plate"] / 2, g["zc"] + g["t_plate"] / 2))
    info = {}
    for name, (R, xa, za, zb) in solids.items():
        poly = segment_polygon(R, xa, tol)
        write_stl(out / f"{name}.stl", prism_triangles(poly, za, zb), name)
        al = np.arccos(xa / R)
        exact = (R * R * al - xa * np.sqrt(R * R - xa * xa)) * (zb - za) / 1000    # circular segment x height (ml)
        info[name] = dict(R_mm=R, x_flat_mm=xa, z_mm=[za, zb], exact_volume_ml=exact, **mesh_check(out / f"{name}.stl"))
        print(name, info[name], flush=True)
    info.update(plate_centre_z_mm=g["zc"], plate_thickness_mm=g["t_plate"], front_face_x_world_mm=g["x_f"],
                tessellation=dict(max_sagitta_mm=tol), model="L2")
    json.dump(info, open(out / "meshes.json", "w"), indent=1)
    return info


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phantom", choices=list(PHANTOM_MODELS), default=None)
    a = ap.parse_args()
    name = set_phantom(a.phantom)
    print("%s: %s" % (name, model()["description"]))
    if is_l2():
        print(json.dumps({k: round(v, 3) for k, v in l2().items()}))
    zc, t = plate()
    print("height %.1f mm, plate z %.2f +- %.2f mm, world x offset %.2f mm, reach %.1f mm" % (height(), zc, t / 2, x_offset(),
                                                                                           reach_mm()))
    for d in discs():
        print("  %-24s %-6s %4.1f mm  x %7.2f  y %7.2f" % (d["name"], d["key"], d["diam"], d["x"], d["y"]))


if __name__ == "__main__":
    main()

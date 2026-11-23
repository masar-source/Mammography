#!/usr/bin/env python3
"""
Meshes of the L2 v2 phantom for GATE (exact CadQuery solids from L2_container_v2, tessellated finely; the
deviation of the facets from the exact surfaces is measured by sphere_packing.py and sets its wall clearance).

Volume hierarchy used in GATE (all in phantom coordinates; z from the bottom outer surface):
  envelope  (PMMA)  = outer shape of the container incl. the 20 mm extension      -> placed in the world
  filling   (oil)   = envelope minus container material (cavity, slots)            -> daughter of envelope
  plate     (PMMA)  = new iodine plate, trimmed at the open end of the extension      -> daughter of filling
  spheres, discs    = analytic solids (sphere_packing.py, disc_layout.csv)          -> daughters of filling / plate
Lid, gasket and positioning bracket lie outside the beam (x < -20 mm or below the support) and are not modelled.
Simplifications outside the field (x < 0): the extension keeps the rim profile (its 10 mm end walls and the plate slots
are left out, because CadQuery's union with the thicker ends gives a non-watertight mesh), and the plate tab is
trimmed to the cavity width. Neither is in the primary beam.

Usage (needs CadQuery): python phantom_meshes.py   -> phantom/{envelope,filling,plate}.stl, phantom/meshes.json
Real L2 (half-cylinder prisms, written by phantoms.py, no CadQuery):
       python phantom_meshes.py --phantom L2   -> phantom_L2/{envelope,filling,plate}.stl, phantom_L2/meshes.json
"""
import argparse, json, sys
from pathlib import Path
from config import HERE, V2_DIR, PHANTOM_MODELS, set_phantom


def mesh_check(path):
    """closed (every edge shared by exactly two facets) and consistently oriented (each directed edge used once)"""
    import numpy as np
    raw = open(path, 'rb').read()
    if raw[:5] == b'solid':                                   # ASCII STL
        import re
        v = np.array(re.findall(rb'vertex\s+(\S+)\s+(\S+)\s+(\S+)', raw), float).reshape(-1, 3, 3)
    else:
        n = int(np.frombuffer(raw[80:84], '<u4')[0])
        v = np.frombuffer(raw[84:84 + 50 * n], dtype=[('n', '<f4', 3), ('v', '<f4', (3, 3)), ('a', '<u2')])['v'].astype(float)
    P = v.reshape(-1, 3); keys, inv = np.unique(np.round(P, 5), axis=0, return_inverse=True); F = inv.reshape(-1, 3)
    E = np.r_[F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]]
    _, cnt = np.unique(np.sort(E, 1), axis=0, return_counts=True)
    _, dcnt = np.unique(E, axis=0, return_counts=True)
    vol = float(np.einsum('ij,ij->i', keys[F[:, 0]], np.cross(keys[F[:, 1]], keys[F[:, 2]])).sum() / 6)
    return dict(triangles=int(len(F)), closed=bool(np.all(cnt == 2)), oriented=bool(np.all(dcnt == 1)), mesh_volume_ml=vol / 1000)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tol", type=float, default=0.005, help="linear tessellation tolerance, absolute (mm)")
    ap.add_argument("--ang", type=float, default=0.1, help="angular tessellation tolerance (rad)")
    ap.add_argument("--phantom", choices=list(PHANTOM_MODELS), default=None)
    ap.add_argument("--out", default=None, help="default: the model's folder (phantom/ or phantom_L2/)")
    a = ap.parse_args()
    set_phantom(a.phantom)
    import phantoms
    out = Path(a.out) if a.out else phantoms.phantom_dir()
    if phantoms.is_l2():
        phantoms.write_l2_meshes(out, a.tol)
        return
    sys.path.insert(0, str(V2_DIR))
    import l2v2_step as S          # exact solids of the v2 container (needs CadQuery)
    out.mkdir(parents=True, exist_ok=True)
    envelope = S.half_rev(S.a_out, S.b_out).union(S.stadium_prism(S.a_out, S.b_out, S.x_open, 0))
    # cavity (centred at ci) incl. the extension, minus the thickened walls of the bolting zone
    cav = S.half_rev(S.a_in, S.b_in, S.ci).union(S.stadium_prism(S.a_in, S.b_in, S.x_open, 0, S.ci))
    frame = S.stadium_prism(S.a_in, S.b_in, S.x_open - 1, S.x_open + S.P['d_frame'], S.ci).intersect(
        S.box(S.x_open - 2, S.x_open + S.P['d_frame'] + 1, -S.R_ext, S.R_ext, 0, S.z_bot_open).union(
            S.box(S.x_open - 2, S.x_open + S.P['d_frame'] + 1, -S.R_ext, S.R_ext, S.z_top_open, S.H)))
    filling = cav.cut(frame)
    y_cut = S.D['R_cav'] - 0.5
    plate = S.plate().intersect(S.box(S.x_open + 0.05, S.R_plate + 1, -y_cut, y_cut, S.ci - 5, S.ci + 5))
    info = {}
    for name, wp in (("envelope", envelope), ("filling", filling), ("plate", plate)):
        s = wp.val()
        # absolute deflection (CadQuery's default STL export scales the tolerance by the edge size)
        s.exportStl(str(out / f"{name}.stl"), tolerance=a.tol, angularTolerance=a.ang, relative=False)
        info[name] = dict(valid=bool(s.isValid()), volume_ml=s.Volume() / 1000, bbox=[round(v, 3) for v in S.tight_bbox(s)],
                          **mesh_check(out / f"{name}.stl"))
        print(name, info[name], flush=True)
    info["plate_centre_z_mm"] = S.ci
    info["plate_thickness_mm"] = S.P["t_plate"]
    info["tessellation"] = dict(tolerance_mm=a.tol, angular_tolerance_rad=a.ang)
    json.dump(info, open(out / "meshes.json", "w"), indent=1)


if __name__ == "__main__":
    main()

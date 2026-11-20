#!/usr/bin/env python3
"""
Random packing of the six sphere sizes (equal volume per size) in the free space of the phantom (--phantom: L2 v2 =
cavity minus the thickened walls of the bolting zone and minus the iodine plate; L2 = real half-cylinder cavity minus
the plate, geometry in phantoms.py), by overlap relaxation with gradual growth of the radii.

Spheres may touch the walls, the plate and each other, as in the real phantom. The only clearances kept are numerical:
WALL = 0.02 mm (the GATE mesh facets lie up to 0.007 mm inside the exact walls at 0.005 mm absolute tessellation,
and a sphere protruding through a facet is a geometry overlap) and GAP = 0.001 mm between spheres. The final check is
against the facets themselves (verify_meshes).

Distances to the walls are analytic (L2 v2: same geometry as phantom_meshes.py / l2v2_step.py, read from params.json;
L2: phantoms.l2_free_distance).

Output (phantom coordinates, mm; z from the bottom outer surface of the container), in the model's folder
(phantom/ for L2 v2, phantom_L2/ for L2):
  spheres.npz  centres (N,3), radii (N,), size_class (N,)
  spheres_summary.json
Usage: python sphere_packing.py [--phantom L2v2|L2] [--fraction 0.61] [--seed 1] [--resume]
"""
import argparse, json, time
import numpy as np
from scipy.spatial import cKDTree
from config import HERE, PHANTOM, PHANTOM_MODELS, set_phantom
import phantoms

WALL = 0.02
GAP = 0.001
TOL_GROW = 0.05     # mm, overlap allowed before the radii grow again (intermediate steps only)
RELIEF_MM = 0.02    # overlap (mm) that counts as jammed
MAX_REMOVED = 0.015 # at most 1.5 % of the sphere volume removed by jam relief

P = json.load(open(PHANTOM["params_json"]))
D = P["derived"]
R_FLAT, CI, A_IN, B_IN = P["R_flat"], D["ci"], D["a_in"], D["b_in"]
X_OPEN, D_FRAME, Z_BOT_OPEN, Z_TOP_OPEN = D["x_open"], P["d_frame"], D["z_bot_open"], D["z_top_open"]
T_PLATE, R_PLATE = P["t_plate"], D["R_plate"]
Y_CUT = D["R_cav"] - 0.5                     # plate trimmed to the cavity width in the GATE model (phantom_meshes.py)
X_TAB_STOP, X_TAB_END, P_TAB = D["x_tab_stop"], D["x_tab_end"], D["P_tab"]


def ellipse_distance(q, w, a, b, iters=64):
    """unsigned distance from points (q, w) to an ellipse with semi-axes a >= b (centred at 0), inside or outside;
    robust bisection on the foot-point equation (D. Eberly, Distance from a point to an ellipse)"""
    y0, y1 = np.abs(np.asarray(q, float)), np.abs(np.asarray(w, float))
    r0, r1 = a * y0, b * y1
    lo = -b * b + r1                                       # F(lo) >= 0
    hi = -b * b + np.hypot(r0, r1)                         # F(hi) <= 0
    with np.errstate(divide='ignore', invalid='ignore'):
        for _ in range(iters):
            t = 0.5 * (lo + hi)
            F = (r0 / (t + a * a)) ** 2 + (r1 / (t + b * b)) ** 2 - 1
            lo, hi = np.where(F > 0, t, lo), np.where(F > 0, hi, t)
        t = 0.5 * (lo + hi)
        x0, x1 = a * a * y0 / (t + a * a), b * b * y1 / (t + b * b)
    # points on the major axis (y1 = 0): the bisection bracket degenerates; closed form
    ax = y1 < 1e-12
    xa = np.where(y0 < (a * a - b * b) / a, a * a * y0 / (a * a - b * b), a)
    x0 = np.where(ax, xa, x0)
    x1 = np.where(ax, b * np.sqrt(np.clip(1 - (xa / a) ** 2, 0, None)), x1)
    return np.hypot(x0 - y0, x1 - y1)


def quadrant_sdf(u, v):
    """signed distance to the region {u <= 0, v <= 0} (positive outside)"""
    return np.where((u > 0) | (v > 0), np.hypot(np.maximum(u, 0), np.maximum(v, 0)), np.maximum(u, v))


def free_distance(Pt):
    """signed distance (mm) to the nearest wall or plate surface, positive in the free space (exact inside it,
    a lower bound outside)"""
    x, y, z = Pt[:, 0], Pt[:, 1], Pt[:, 2]
    rho = np.where(x >= 0, np.hypot(x, y), np.abs(y))       # revolved rim in the field, stadium prism in the extension
    q, w = rho - R_FLAT, z - CI
    inside_ell = (q / A_IN) ** 2 + (w / B_IN) ** 2 <= 1
    d_ell = ellipse_distance(q, w, A_IN, B_IN)
    d_cav = np.where(q > 0, np.where(inside_ell, d_ell, -d_ell), B_IN - np.abs(w))
    d_lid = x - X_OPEN                                        # inner face of the lid (open end of the extension)
    xf = X_OPEN + D_FRAME                                     # thickened walls of the bolting zone
    d_fbot = quadrant_sdf(x - xf, z - Z_BOT_OPEN)
    d_ftop = quadrant_sdf(x - xf, Z_TOP_OPEN - z)
    # plate (GATE mesh): half-disc for x >= x_tab_stop plus the tab, cut at x_open + 0.05 and |y| <= Y_CUT
    disc = np.maximum(np.hypot(x, y) - R_PLATE, X_TAB_STOP - x)
    tab = np.maximum.reduce([np.abs(y) - P_TAB, x - X_TAB_STOP, X_TAB_END - x])
    dh = np.maximum.reduce([np.minimum(disc, tab), X_OPEN + 0.05 - x, np.abs(y) - Y_CUT])
    d_pl = quadrant_sdf(dh, np.abs(z - CI) - T_PLATE / 2)
    return np.minimum.reduce([d_cav, d_lid, d_fbot, d_ftop, d_pl])


def free_gradient(Pt, eps=0.01):
    g = np.zeros_like(Pt)
    for k in range(3):
        e = np.zeros(3); e[k] = eps
        g[:, k] = (free_distance(Pt + e) - free_distance(Pt - e)) / (2 * eps)
    return g / (np.linalg.norm(g, axis=1, keepdims=True) + 1e-12)


def use_model(name=None):
    """select the phantom model; the real L2 replaces free_distance (used by every function below) by its own"""
    global free_distance
    set_phantom(name)
    if phantoms.is_l2():
        free_distance = phantoms.l2_free_distance
    return phantoms.phantom_dir()


def free_volume(seed=0):
    """Monte Carlo estimate of the free volume (ml)"""
    rng = np.random.default_rng(seed)
    if phantoms.is_l2():
        lo, hi = phantoms.l2_packing_box()
    else:
        lo = np.array([X_OPEN, -(R_FLAT + A_IN), CI - B_IN]); hi = np.array([R_FLAT + A_IN, R_FLAT + A_IN, CI + B_IN])
    n = 4_000_000
    Pt = rng.uniform(lo, hi, size=(n, 3))
    return float((free_distance(Pt) > 0).mean() * np.prod(hi - lo) / 1000), (lo, hi)


def sample_inside(n, rmin, rng, box):
    out = []
    while sum(len(o) for o in out) < n:
        Pt = rng.uniform(box[0], box[1], size=(4 * n + 2000, 3))
        out.append(Pt[free_distance(Pt) >= rmin])
    return np.concatenate(out)[:n]


def relax(C, R, cls, s, ncls, iters, tol):
    """overlap removal at radii s*R (Jacobi steps). A pair's overlap is shared in inverse proportion to an effective
    mass = radius, x20 for a sphere in contact with a wall, so that a small sphere caught between a large sphere and a
    wall pushes the large one away instead of bouncing between the two."""
    order = [np.nonzero(cls == k)[0] for k in range(ncls)]
    worst = 0.0
    for it in range(iters):
        r = s * R
        fd = free_distance(C)
        need = r + WALL - fd
        m = R * np.where(need > -0.02, 20.0, 1.0)
        trees = [cKDTree(C[o]) for o in order]
        disp = np.zeros_like(C); worst = 0.0
        for a in range(ncls):
            for b in range(a, ncls):
                cut = s * (R[order[a][0]] + R[order[b][0]]) + GAP
                M = trees[a].sparse_distance_matrix(trees[b], cut, output_type='ndarray')
                if len(M) == 0:
                    continue
                i, j, d = order[a][M['i']], order[b][M['j']], M['v']
                if a == b:
                    k = i < j; i, j, d = i[k], j[k], d[k]
                ov = r[i] + r[j] + GAP - d
                k = ov > 0
                if not k.any():
                    continue
                i, j, d, ov = i[k], j[k], np.maximum(d[k], 1e-6), ov[k]
                u = (C[j] - C[i]) / d[:, None]
                wi = m[j] / (m[i] + m[j])
                np.add.at(disp, i, -(ov * wi)[:, None] * u)
                np.add.at(disp, j, (ov * (1 - wi))[:, None] * u)
                worst = max(worst, float(ov.max()))
        wv = need > 0
        if wv.any():
            disp[wv] += need[wv][:, None] * free_gradient(C[wv])
            worst = max(worst, float(need[wv].max()))
        step = np.linalg.norm(disp, axis=1, keepdims=True)
        disp *= np.minimum(1.0, 0.5 * r[:, None] / np.maximum(step, 1e-12))   # limit the move per iteration
        C += disp
        if worst < tol:
            return C, it + 1, worst
    return C, iters, worst


def jammed(C, R, thresh):
    """spheres to remove: the smaller of each pair overlapping by more than thresh, and spheres pressed into a wall"""
    drop = (R + WALL - free_distance(C)) > thresh
    pairs = cKDTree(C).query_pairs(2 * R.max() + GAP, output_type='ndarray')
    d = np.linalg.norm(C[pairs[:, 0]] - C[pairs[:, 1]], axis=1)
    ov = R[pairs[:, 0]] + R[pairs[:, 1]] + GAP - d
    for i, j in pairs[ov > thresh]:
        if not (drop[i] or drop[j]):
            drop[i if R[i] < R[j] or (R[i] == R[j] and i > j) else j] = True
    return drop


def verify(C, R):
    bad_wall = free_distance(C) < R + WALL * 0.5
    tree = cKDTree(C)
    pairs = tree.query_pairs(2 * R.max(), output_type='ndarray')
    d = np.linalg.norm(C[pairs[:, 0]] - C[pairs[:, 1]], axis=1)
    ov = R[pairs[:, 0]] + R[pairs[:, 1]] - d
    return bad_wall, pairs[ov > 0], float(max(ov.max(), 0)) if len(ov) else 0.0


def load_stl(path):
    raw = open(path, 'rb').read()
    if raw[:5] == b'solid':
        import re
        return np.array(re.findall(rb'vertex\s+(\S+)\s+(\S+)\s+(\S+)', raw), float).reshape(-1, 3, 3)
    n = int(np.frombuffer(raw[80:84], '<u4')[0])
    return np.frombuffer(raw[84:84 + 50 * n], dtype=[('n', '<f4', 3), ('v', '<f4', (3, 3)), ('a', '<u2')])['v'].astype(float)


def point_triangle_distance(p, A, B, C):
    """exact distance from points p to triangles ABC (all (n,3)); Ericson, Real-Time Collision Detection 5.1.5"""
    ab, ac, ap = B - A, C - A, p - A
    d1, d2 = (ab * ap).sum(1), (ac * ap).sum(1)
    bp = p - B; d3, d4 = (ab * bp).sum(1), (ac * bp).sum(1)
    cp = p - C; d5, d6 = (ab * cp).sum(1), (ac * cp).sum(1)
    va, vb, vc = d3 * d6 - d5 * d4, d5 * d2 - d1 * d6, d1 * d4 - d3 * d2
    with np.errstate(divide='ignore', invalid='ignore'):
        den = va + vb + vc
        v, w = vb / den, vc / den
        X = A + ab * v[:, None] + ac * w[:, None]                                   # interior
        ce = (vb <= 0) & (d2 >= 0) & (d6 <= 0); t = d2 / (d2 - d6)                 # edge AC
        X = np.where(ce[:, None], A + ac * t[:, None], X)
        be = (va <= 0) & (d4 - d3 >= 0) & (d5 - d6 >= 0); t = (d4 - d3) / ((d4 - d3) + (d5 - d6))   # edge BC
        X = np.where(be[:, None], B + (C - B) * t[:, None], X)
        ae = (vc <= 0) & (d1 >= 0) & (d3 <= 0); t = d1 / (d1 - d3)                 # edge AB
        X = np.where(ae[:, None], A + ab * t[:, None], X)
    X = np.where(((d6 >= 0) & (d5 <= d6))[:, None], C, X)
    X = np.where(((d3 >= 0) & (d4 <= d3))[:, None], B, X)
    X = np.where(((d1 <= 0) & (d2 <= 0))[:, None], A, X)
    return np.linalg.norm(p - X, axis=1)


def verify_meshes(C, R, near=0.5):
    """clearance of every sphere near a wall from the actual GATE meshes (filling.stl, plate.stl): exact
    point-triangle distances minus the radius; must be > 0 (no sphere crosses a mesh facet)"""
    near_idx = np.nonzero(free_distance(C) - R < near)[0]
    clear = np.full(len(C), np.inf)
    tree = cKDTree(C[near_idx])
    for name in ("filling", "plate"):
        T = load_stl(phantoms.phantom_dir() / f"{name}.stl")
        cen = T.mean(1); rad = np.linalg.norm(T - cen[:, None], axis=2).max(1) + R.max() + near
        for s0 in range(0, len(T), 20000):
            lists = tree.query_ball_point(cen[s0:s0 + 20000], rad[s0:s0 + 20000])
            ti = np.repeat(np.arange(s0, s0 + len(lists)), [len(l) for l in lists])
            if len(ti) == 0:
                continue
            si = near_idx[np.concatenate([np.asarray(l, int) for l in lists])]
            for c0 in range(0, len(ti), 4_000_000):
                t_, s_ = ti[c0:c0 + 4_000_000], si[c0:c0 + 4_000_000]
                d = point_triangle_distance(C[s_], T[t_, 0], T[t_, 1], T[t_, 2]) - R[s_]
                np.minimum.at(clear, s_, d)
    return clear, near_idx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phantom", choices=list(PHANTOM_MODELS), default=None, help="default: REVELATION_PHANTOM or L2v2")
    ap.add_argument("--fraction", type=float, default=PHANTOM["sphere_fraction"])
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--max-iter", type=int, default=60000)
    ap.add_argument("--resume", action="store_true", help="continue from phantom/packing_checkpoint.npz (saved at every growth step)")
    a = ap.parse_args()
    out = use_model(a.phantom)
    print("phantom %s -> %s" % (PHANTOM["model"], out), flush=True)
    t0 = time.time(); rng = np.random.default_rng(a.seed)
    Vfree, box = free_volume()
    print('free space %.1f ml (Monte Carlo)' % Vfree, flush=True)
    diam = np.array(PHANTOM["sphere_diameters_mm"]); radii = diam / 2
    counts = np.round(a.fraction * Vfree * 1000 / len(diam) / (np.pi / 6 * diam ** 3)).astype(int)
    print('target counts per size:', dict(zip(diam.tolist(), counts.tolist())), flush=True)
    s0 = 0.55
    C, R, cls = [], [], []
    for k, (r, n) in enumerate(zip(radii, counts)):
        C.append(sample_inside(n, s0 * r + WALL, rng, box)); R.append(np.full(n, r)); cls.append(np.full(n, k))
    C, R, cls = np.concatenate(C), np.concatenate(R), np.concatenate(cls)
    # grow the radii in steps; overlaps need only be small (< TOL_GROW) before the next step, and fully removed at the end.
    # At full size a random packing at this density can jam locally: when the largest overlap has stopped shrinking
    # (< 10 % over the last 3 chunks of 200 iterations), the smaller sphere of each pair still overlapping by more than
    # max(RELIEF_MM, half the largest overlap) is removed (total removed volume capped at MAX_REMOVED).
    s, total, since_step, removed_vol, n_removed = s0, 0, 0, 0.0, 0
    v_target = float((4 / 3 * np.pi * R ** 3).sum())
    ckpt = out / "packing_checkpoint.npz"
    if a.resume and ckpt.exists():
        z = np.load(ckpt); C, R, cls, s, total = z["C"], z["R"], z["cls"], float(z["s"]), int(z["total"])
        print("  resumed from checkpoint at s = %.3f (%d spheres)" % (s, len(R)), flush=True)
    history = []
    while True:
        tol = 1e-3 if s >= 1.0 else TOL_GROW
        C, it, worst = relax(C, R, cls, s, len(diam), iters=200, tol=tol)
        total += it; since_step += it
        print('  s %.3f: %3d iterations, worst overlap %.4f mm, %d spheres (total %d, %.0f s)'
              % (s, it, worst, len(R), total, time.time() - t0), flush=True)
        if s >= 1.0 and worst < 1e-3:
            break
        if total > a.max_iter or (s < 1.0 and since_step > 4000):
            print('  stopped: iteration limit or jammed before full size (s = %.3f)' % s, flush=True)
            break
        if s >= 1.0:
            history.append(worst)
            if len(history) >= 4 and history[-1] > 0.9 * history[-4]:
                drop = jammed(C, R, max(RELIEF_MM, 0.5 * worst))
                dv = float((4 / 3 * np.pi * R[drop] ** 3).sum())
                if removed_vol + dv > MAX_REMOVED * v_target:
                    print('  removal cap reached; stopping relaxation', flush=True)
                    break
                removed_vol += dv; n_removed += int(drop.sum())
                print('  jam relief: removed %d spheres (%s), %.3f %% of the sphere volume so far'
                      % (drop.sum(), ", ".join("%g mm: %d" % (diam[k], (cls[drop] == k).sum()) for k in range(len(diam))
                                               if (cls[drop] == k).any()), 100 * removed_vol / v_target), flush=True)
                C, R, cls = C[~drop], R[~drop], cls[~drop]
                history = []
        if worst < tol and s < 1.0:
            s = min(1.0, s + (0.02 if s < 0.9 else 0.01)); since_step = 0
            np.savez_compressed(ckpt, C=C, R=R, cls=cls, s=s, total=total)   # resume point (--resume)
    bad_wall, bad_pair, _ = verify(C, R)
    drop = bad_wall.copy()
    for i, j in bad_pair:
        drop[i if R[i] <= R[j] else j] = True
    C, R, cls = C[~drop], R[~drop], cls[~drop]
    bad_wall2, bad_pair2, maxov2 = verify(C, R)
    out.mkdir(exist_ok=True)
    np.savez_compressed(out / "spheres.npz", centres=C, radii=R, size_class=cls)
    vol = float((4 / 3 * np.pi * R ** 3).sum())
    fd = free_distance(C) - R
    summ = dict(free_volume_ml=Vfree, target_fraction=a.fraction, achieved_fraction=vol / (Vfree * 1000),
                counts={str(d): int((cls == k).sum()) for k, d in enumerate(diam)}, removed_after_check=int(drop.sum()), final_scale=s,
                removed_jam_relief=dict(spheres=n_removed, volume_fraction=removed_vol / v_target),
                remaining_wall_violations=int(bad_wall2.sum()), remaining_overlaps=int(len(bad_pair2)), max_overlap_mm=maxov2,
                wall_clearance_mm=dict(min=float(fd.min()), median=float(np.median(fd)),
                                       touching_within_0p05mm=int((fd < 0.05).sum())),
                phantom=PHANTOM["model"], wall_mm=WALL, gap_mm=GAP, seed=a.seed, iterations=total, seconds=time.time() - t0)
    clear, near_idx = verify_meshes(C, R)
    summ["mesh_check"] = dict(spheres_checked=int(len(near_idx)), min_clearance_mm=float(clear[near_idx].min()),
                              crossing_a_facet=int((clear[near_idx] <= 0).sum()))
    json.dump(summ, open(out / "spheres_summary.json", "w"), indent=1)
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()

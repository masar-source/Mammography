#!/usr/bin/env python3
"""Background spheres for the L2 v2 phantom: one binary STL per diameter (mm), plus a check report.
Geodesic (icosphere) meshes, subdivided to < 2 µm chordal error, then scaled so the mesh volume equals the
ideal sphere volume (the vertices then lie a few µm outside the nominal radius).
Usage: python make_spheres.py  -> spheres/L2v2_sphere_D*.stl, spheres/sphere_checks.json"""
import json, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'spheres'); os.makedirs(OUT, exist_ok=True)
DIAMETERS = [15.88, 12.70, 9.52, 6.35, 3.18, 1.58]

def icosphere(level):
    t = (1 + 5 ** 0.5) / 2
    V = [[-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0], [0, -1, t], [0, 1, t], [0, -1, -t], [0, 1, -t], [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1]]
    F = [[0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11], [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8],
         [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9], [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1]]
    V = [np.array(v, float) / np.linalg.norm(v) for v in V]
    for _ in range(level):
        cache = {}; F2 = []
        def mid(i, j):
            k = (min(i, j), max(i, j))
            if k not in cache:
                m = V[i] + V[j]; V.append(m / np.linalg.norm(m)); cache[k] = len(V) - 1
            return cache[k]
        for a, b, c in F:
            ab, bc, ca = mid(a, b), mid(b, c), mid(c, a); F2 += [[a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]]
        F = F2
    return np.array(V), np.array(F)
def volume(V, F):
    return float(np.einsum('ij,ij->i', V[F[:, 0]], np.cross(V[F[:, 1]], V[F[:, 2]])).sum() / 6)
def write_stl(path, V, F, title):
    tri = V[F].astype(np.float32); n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]); n /= np.linalg.norm(n, axis=1, keepdims=True)
    rec = np.zeros(len(F), dtype=[('n', '<f4', 3), ('v', '<f4', (3, 3)), ('a', '<u2')]); rec['n'] = n; rec['v'] = tri
    with open(path, 'wb') as fh:
        fh.write(title.encode()[:80].ljust(80, b' ')); fh.write(np.uint32(len(F)).tobytes()); fh.write(rec.tobytes())
def read_stl(path):
    raw = open(path, 'rb').read(); n = int(np.frombuffer(raw[80:84], '<u4')[0])
    rec = np.frombuffer(raw[84:84 + 50 * n], dtype=[('n', '<f4', 3), ('v', '<f4', (3, 3)), ('a', '<u2')]); return rec['v'].astype(float), rec['n']

U, F = icosphere(5)                                   # 20 480 faces; edge angle ≈ 2.0° -> sagitta 1.5e-4 r
checks = {}
for D in DIAMETERS:
    r = D / 2; V = U * r
    V *= (4 / 3 * np.pi * r ** 3 / volume(V, F)) ** (1 / 3)      # volume-exact
    path = os.path.join(OUT, 'L2v2_sphere_D%.2fmm.stl' % D)
    write_stl(path, V, F, 'L2 v2 background sphere, diameter %.2f mm, units mm' % D)
    # --- checks on the file as written
    tri, nrm = read_stl(path); P = tri.reshape(-1, 3)
    keys, inv = np.unique(np.round(P, 6), axis=0, return_inverse=True); Fi = inv.reshape(-1, 3)
    E = np.sort(np.r_[Fi[:, [0, 1]], Fi[:, [1, 2]], Fi[:, [2, 0]]], axis=1); _, cnt = np.unique(E, axis=0, return_counts=True)
    cen = tri.mean(1); outward = float(np.mean(np.sum(nrm * cen, 1) > 0))
    rad = np.linalg.norm(P, axis=1); cr = np.linalg.norm(cen, axis=1)
    vol = volume(keys, Fi)
    checks['%.2f' % D] = dict(file=os.path.basename(path), triangles=len(tri), watertight=bool(np.all(cnt == 2)), outward_normals=outward,
                              volume_mm3=vol, volume_error_pct=100 * (vol / (np.pi / 6 * D ** 3) - 1),
                              vertex_diameter_mm=[2 * rad.min(), 2 * rad.max()], min_face_centre_diameter_mm=2 * cr.min(),
                              size_kB=os.path.getsize(path) / 1024)
    c = checks['%.2f' % D]
    print('D %5.2f mm: %d tris, watertight %s, outward %.3f, volume %.4f mm3 (%+.4f %%), vertex D %.4f-%.4f, facet D min %.4f, %.0f kB'
          % (D, c['triangles'], c['watertight'], c['outward_normals'], c['volume_mm3'], c['volume_error_pct'], *c['vertex_diameter_mm'], c['min_face_centre_diameter_mm'], c['size_kB']))
json.dump(checks, open(os.path.join(OUT, 'sphere_checks.json'), 'w'), indent=1)

#!/usr/bin/env python3
"""Illustrative renders of the L2 v2 container, rev B (software rasteriser, clear acrylic shown translucent).
Usage: python l2v2_render.py   -> renders/*.png"""
import os, time, numpy as np
from l2v2_geom import p, c, ci, R_out, R_plate, x_open, x_lid_out, x_lid_in, a_ext, SCREWS, PORTS, f_container, f_lid, f_gasket, f_cavity, f_plate, f_bracket, mesh, plate_mesh, disc_layout, PARTS
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'renders'); os.makedirs(OUT, exist_ok=True)

# ------------------------------------------------------------------ simple meshes
def cylinder(cen, axis, r, h, n=40):
    axis = np.asarray(axis, float); axis /= np.linalg.norm(axis)
    a = np.cross(axis, [0, 0, 1.0]) if abs(axis[2]) < 0.9 else np.cross(axis, [1.0, 0, 0]); a /= np.linalg.norm(a); b = np.cross(axis, a)
    t = np.linspace(0, 2 * np.pi, n, endpoint=False); ring = np.cos(t)[:, None] * a + np.sin(t)[:, None] * b
    c0 = np.asarray(cen, float); c1 = c0 + axis * h
    V = np.vstack([c0 + r * ring, c1 + r * ring, c0, c1]); F = []; N = np.vstack([ring, ring, -axis, axis])
    for i in range(n):
        j = (i + 1) % n; F += [[i, j, n + j], [i, n + j, n + i], [2 * n, j, i], [2 * n + 1, n + i, n + j]]
    # separate cap normals would be nicer; flat-ish shading is fine at this size
    return V.astype(np.float32), np.array(F), N.astype(np.float32)
def icosphere(level):
    t = (1 + 5 ** 0.5) / 2
    V = [[-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0], [0, -1, t], [0, 1, t], [0, -1, -t], [0, 1, -t], [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1]]
    F = [[0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11], [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8], [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9], [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1]]
    V = [np.array(v, float) / np.linalg.norm(v) for v in V]
    for _ in range(level):
        cache = {}; F2 = []
        def mid(i, j):
            k = (min(i, j), max(i, j))
            if k not in cache:
                m = V[i] + V[j]; V.append(m / np.linalg.norm(m)); cache[k] = len(V) - 1
            return cache[k]
        for a, b, cc in F:
            ab, bc, ca = mid(a, b), mid(b, cc), mid(cc, a); F2 += [[a, ab, ca], [b, bc, ab], [cc, ca, bc], [ab, bc, ca]]
        F = F2
    return np.array(V, np.float32), np.array(F)
ICO = {k: icosphere(k) for k in (1, 2, 3)}
def boxmesh(x0, x1, y0, y1, z0, z1):
    """axis-aligned box with flat-shaded faces"""
    Vs, Fs, Ns = [], [], []
    for ax_, sgn in ((0, -1), (0, 1), (1, -1), (1, 1), (2, -1), (2, 1)):
        lo = np.array([x0, y0, z0], float); hi = np.array([x1, y1, z1], float); u, w = [k for k in range(3) if k != ax_]
        q = []
        for a_, b_ in ((0, 0), (1, 0), (1, 1), (0, 1)):
            pnt = lo.copy(); pnt[ax_] = hi[ax_] if sgn > 0 else lo[ax_]; pnt[u] = hi[u] if a_ else lo[u]; pnt[w] = hi[w] if b_ else lo[w]; q.append(pnt)
        n = np.zeros(3); n[ax_] = sgn; o = len(Vs) * 4
        Vs.append(np.array(q)); Ns.append(np.tile(n, (4, 1))); Fs += [[o, o + 1, o + 2], [o, o + 2, o + 3]]
    return np.vstack(Vs).astype(np.float32), np.array(Fs), np.vstack(Ns).astype(np.float32)
def merge(parts):
    Vs, Fs, Ns, o = [], [], [], 0
    for V, F, N in parts:
        Vs.append(V); Fs.append(F + o); Ns.append(N); o += len(V)
    return np.vstack(Vs).astype(np.float32), np.vstack(Fs), np.vstack(Ns).astype(np.float32)
def shifted(m, d): return (m[0] + np.asarray(d, np.float32), m[1], m[2])

# ------------------------------------------------------------------ rasteriser
class Cam:
    def __init__(s, eye, target, W, H, fov=30, up=(0, 0, 1)):
        s.eye = np.array(eye, float); fw = np.array(target, float) - s.eye; fw /= np.linalg.norm(fw)
        rt = np.cross(fw, up); rt /= np.linalg.norm(rt); u = np.cross(rt, fw)
        s.R = np.vstack([rt, u, fw]); s.W, s.H = W, H; s.f = 0.5 * W / np.tan(np.radians(fov) / 2); s.fw = fw
        s.ox = s.oy = 0.0
    def proj(s, P):
        q = (np.asarray(P, float) - s.eye) @ s.R.T
        return np.c_[s.W / 2 + s.ox + s.f * q[:, 0] / q[:, 2], s.H / 2 + s.oy - s.f * q[:, 1] / q[:, 2], q[:, 2]]
    def fit(s, pts, margin=(0.06, 0.06, 0.12, 0.05)):
        """scale and shift so that pts fill the frame; margin = (left, right, top, bottom) fractions"""
        s.ox = s.oy = 0.0; S = s.proj(pts); x0, x1, y0, y1 = S[:, 0].min(), S[:, 0].max(), S[:, 1].min(), S[:, 1].max()
        k = min(s.W * (1 - margin[0] - margin[1]) / (x1 - x0), s.H * (1 - margin[2] - margin[3]) / (y1 - y0)); s.f *= k
        S = s.proj(pts); x0, x1, y0, y1 = S[:, 0].min(), S[:, 0].max(), S[:, 1].min(), S[:, 1].max()
        s.ox = (s.W * margin[0] + s.W * (1 - margin[1])) / 2 - (x0 + x1) / 2; s.oy = (s.H * margin[2] + s.H * (1 - margin[3])) / 2 - (y0 + y1) / 2
LIGHTS = [(np.array([-0.35, -0.55, 0.76]), 0.62), (np.array([0.6, 0.3, 0.75]), 0.28), (np.array([0.2, 0.9, -0.2]), 0.12)]
def shade(col, n, view):
    n = n / (np.linalg.norm(n, axis=-1, keepdims=True) + 1e-12)
    l = 0.32 + sum(w * np.clip(n @ (d / np.linalg.norm(d)), 0, None) for d, w in LIGHTS)
    h = (LIGHTS[0][0] / np.linalg.norm(LIGHTS[0][0]) - view); h /= np.linalg.norm(h)
    spec = 0.18 * np.clip(n @ h, 0, None) ** 40
    return np.clip(np.asarray(col)[None, :] * l[:, None] + spec[:, None], 0, 1)
def render(scene, cam, bg=(1, 1, 1)):
    W, H = cam.W, cam.H
    Z = np.full((H, W), np.inf); C = np.ones((H, W, 3)) * np.array(bg)
    acc = np.zeros((H, W, 3)); accA = np.zeros((H, W)); rev = np.ones((H, W))
    def raster(V, F, N, col, glass, facecol=None):
        fn = np.cross(V[F[:, 1]] - V[F[:, 0]], V[F[:, 2]] - V[F[:, 0]]); flip = np.sum(fn * N[F].sum(1), 1) < 0
        F = F.copy(); F[flip] = F[flip][:, ::-1]                    # consistent outward winding from the vertex normals
        S = cam.proj(V); tri = S[F]; ok = (tri[:, :, 2] > 1).all(1)
        x0 = np.floor(tri[:, :, 0].min(1)).astype(int); x1 = np.ceil(tri[:, :, 0].max(1)).astype(int)
        y0 = np.floor(tri[:, :, 1].min(1)).astype(int); y1 = np.ceil(tri[:, :, 1].max(1)).astype(int)
        area = (tri[:, 1, 0] - tri[:, 0, 0]) * (tri[:, 2, 1] - tri[:, 0, 1]) - (tri[:, 2, 0] - tri[:, 0, 0]) * (tri[:, 1, 1] - tri[:, 0, 1])
        ok &= (x1 >= 0) & (x0 < W) & (y1 >= 0) & (y0 < H) & (np.abs(area) > 1e-9)
        if not glass: ok &= area < 0          # back-face culling (screen y is flipped)
        for k in np.nonzero(ok)[0]:
            a, b, cc = tri[k]; xa, xb = max(x0[k], 0), min(x1[k], W - 1); ya, yb = max(y0[k], 0), min(y1[k], H - 1)
            if xb < xa or yb < ya: continue
            gx, gy = np.meshgrid(np.arange(xa, xb + 1) + 0.5, np.arange(ya, yb + 1) + 0.5)
            w0 = ((b[0] - gx) * (cc[1] - gy) - (cc[0] - gx) * (b[1] - gy)) / area[k]
            w1 = ((cc[0] - gx) * (a[1] - gy) - (a[0] - gx) * (cc[1] - gy)) / area[k]
            w2 = 1 - w0 - w1; m = (w0 >= -1e-6) & (w1 >= -1e-6) & (w2 >= -1e-6)
            if not m.any(): continue
            z = w0 * a[2] + w1 * b[2] + w2 * cc[2]
            sub = (slice(ya, yb + 1), slice(xa, xb + 1)); zb = Z[sub]
            m &= z < zb
            if not m.any(): continue
            n = (w0[m, None] * N[F[k, 0]] + w1[m, None] * N[F[k, 1]] + w2[m, None] * N[F[k, 2]])
            P3 = w0[m, None] * V[F[k, 0]] + w1[m, None] * V[F[k, 1]] + w2[m, None] * V[F[k, 2]]
            view = P3 - cam.eye; view /= np.linalg.norm(view, axis=1, keepdims=True)
            cl = facecol[k] if facecol is not None else col
            if glass:
                nn = n / (np.linalg.norm(n, axis=1, keepdims=True) + 1e-12)
                al = 0.05 + 0.62 * (1 - np.abs(np.sum(nn * view, 1))) ** 3
                sh = shade(cl, nn * np.sign(-np.sum(nn * view, 1))[:, None], view.mean(0))
                acc[sub][m] += al[:, None] * sh; accA[sub][m] += al; rev[sub][m] *= (1 - al)
            else:
                zb[m] = z[m]; C[sub][m] = shade(cl, n, view.mean(0))
    t0 = time.time()
    for it in [s for s in scene if not s.get('glass')]: raster(it['V'], it['F'], it['N'], it['col'], False, it.get('facecol'))
    gl = [s for s in scene if s.get('glass')]
    if gl:
        Zo = Z.copy()
        for it in gl:
            Z[:] = Zo                      # glass fragments only in front of the opaque surfaces
            raster(it['V'], it['F'], it['N'], it['col'], True)
        Z[:] = Zo
        cg = acc / np.maximum(accA, 1e-9)[..., None]
        C = C * rev[..., None] + cg * (1 - rev[..., None])
    print('  render %.0f s' % (time.time() - t0), flush=True)
    return C
def down(img, s):
    H, W, _ = img.shape; return img[:H // s * s, :W // s * s].reshape(H // s, s, W // s, s, 3).mean((1, 3))

# ------------------------------------------------------------------ scene parts
ACRYL = (0.70, 0.84, 0.95); LIDC = (0.62, 0.78, 0.94); CUT = (0.55, 0.70, 0.86); PLATE = (0.93, 0.93, 0.90); IOD = (0.42, 0.22, 0.55)
GASK = (0.18, 0.18, 0.2); NYL = (0.96, 0.95, 0.90); SPH = (0.86, 0.80, 0.70); BRK = (0.55, 0.62, 0.70); SUP = (0.42, 0.44, 0.47)
def hardware(dx=0.0):
    hs = [cylinder((x_lid_out + dx, y, z), (-1, 0, 0), 2.75, 2.0, 24) for (y, z) in SCREWS]
    hs += [cylinder((x_lid_out + dx, y, z), (-1, 0, 0), 6.0, 3.0, 6) for (y, z) in PORTS]
    return merge(hs)
def discs(dx=0.0):
    ds = []
    for d in disc_layout():
        ds.append(cylinder((d['x'] + dx, d['y'], ci - p['t_plate'] / 2 - 0.06), (0, 0, 1), d['diam'] / 2, p['t_plate'] + 0.12, 32))
    return merge(ds)
def cut_colours(V, F, N, axis=1, val=0.0, col=CUT, base=ACRYL):
    fc = np.tile(np.array(base, float), (len(F), 1)); cen = V[F].mean(1); fnm = N[F].mean(1)
    cap = (np.abs(cen[:, axis] - val) < 0.6) & (np.abs(fnm[:, axis]) > 0.95); fc[cap] = col; return fc
def spheres_rsa(yside=1, frac=0.45, seed=3):
    """random sequential packing (illustration only, < real 61 %), 6 sizes with equal volume shares; one side (y>0)"""
    rng = np.random.default_rng(seed); sizes = [16, 12.7, 9.5, 6.4, 4.8, 3.2]
    dvol = 1092e3 / 2; target = frac * dvol / len(sizes); placed = []; grid = {}
    def free(q, r):
        if f_cavity(*q) > -r - 0.05 or f_container(*q) < r + 0.05 or abs(q[2] - c) < p['t_plate'] / 2 + r + 0.05: return False
        g = tuple((np.array(q) // 16).astype(int))
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for (q2, r2) in grid.get((g[0] + dx, g[1] + dy, g[2] + dz), ()):
                        if (q[0] - q2[0]) ** 2 + (q[1] - q2[1]) ** 2 + (q[2] - q2[2]) ** 2 < (r + r2) ** 2: return False
        return True
    for D in sizes:
        r = D / 2; n = int(target / (np.pi / 6 * D ** 3)); got = 0; tries = 0
        while got < n and tries < 60 * n:
            tries += 1
            q = (rng.uniform(x_open, R_out), yside * rng.uniform(r, R_out), rng.uniform(0, p['H']))
            if free(q, r):
                placed.append((q, r)); grid.setdefault(tuple((np.array(q) // 16).astype(int)), []).append((q, r)); got += 1
    return placed
def sphere_mesh(placed):
    ms = []
    for q, r in placed:
        V, F = ICO[3 if r > 5 else (2 if r > 2.5 else 1)]; ms.append((V * r + np.array(q, np.float32), F, V))
    return merge(ms)

def annotate(img, cam, notes, path, title, ss, lines=()):
    """notes: (anchor_3d, text, 'top'|'bottom'); labels spread along a band, ordered by anchor position to avoid crossings.
    lines: (polyline_3d, colour, label) drawn as dashed overlays."""
    H, W, _ = img.shape
    fig = plt.figure(figsize=(W / 100, H / 100), dpi=100); ax = fig.add_axes([0, 0, 1, 1]); ax.imshow(img); ax.set_axis_off()
    for poly, col, lab in lines:
        S = cam.proj(np.asarray(poly, float)) / np.array([ss, ss, 1])
        ax.plot(S[:, 0], S[:, 1], color=col, lw=1.6, ls=(0, (6, 3)))
        if lab: ax.text(S[-1, 0] + 6, S[-1, 1], lab, color=col, fontsize=10.5, ha='left', va='center', weight='bold',
                        bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='none', alpha=0.85))
    for band, yb in (('top', 92), ('bottom', H - 52)):
        items = sorted([(cam.proj(np.array([P3]))[0] / np.array([ss, ss, 1]), txt) for P3, txt, bd in notes if bd == band], key=lambda t: t[0][0])
        k = len(items)
        for i, ((u, v, _), txt) in enumerate(items):
            lx = W * (i + 0.5) / k
            ax.annotate(txt, xy=(u, v), xytext=(lx, yb), fontsize=10.5, color='#1a1a1a', ha='center', va='center',
                        arrowprops=dict(arrowstyle='-', color='#444444', lw=0.9, shrinkA=2, shrinkB=0),
                        bbox=dict(boxstyle='round,pad=0.25', fc='white', ec='#bbbbbb', lw=0.6, alpha=0.95))
    ax.text(18, 26, title, fontsize=13, color='#111111', ha='left', va='top', weight='bold')
    fig.savefig(path, dpi=100); plt.close(fig)

if __name__ == '__main__':
    SS = 2; W, H = 1500 * SS, 1000 * SS
    def box(x0, x1, y0, y1, z0, z1): return np.array([[x, y, z] for x in (x0, x1) for y in (y0, y1) for z in (z0, z1)])
    tt = np.linspace(-np.pi / 2, np.pi / 2, 60)
    ring = lambda r, z: np.c_[r * np.cos(tt), r * np.sin(tt), np.full(60, z)]
    YE = p['R_flat'] + a_ext
    PTS1 = np.r_[ring(R_out, c), ring(p['R_flat'], p['H']), ring(p['R_flat'], 0), box(x_lid_out, 0, -YE, YE, -p['lip_h'], p['H'])]
    EDGE = [np.array([[0, -YE - 8, p['H']], [0, YE + 8, p['H']]]), '#c0392b', 'x = 0: chest-wall face of the support\n(the image starts a few mm inside)']
    print('meshing', flush=True)
    cont = mesh(f_container, PARTS['container'][1], 0.9)
    lid = mesh(f_lid, PARTS['lid'][1], 0.6); gas = mesh(f_gasket, PARTS['gasket'][1], 0.6); pl = plate_mesh(); dsk = discs(); hw = hardware()
    brk = mesh(f_bracket, PARTS['bracket'][1], 0.5); sup = boxmesh(0.0, 215.0, -160.0, 160.0, -42.0, -0.05)
    # 1) assembled, empty, clear acrylic
    cam = Cam((-235, -285, 215), (62, 6, 12), W, H, fov=33); cam.fit(PTS1, margin=(0.05, 0.14, 0.16, 0.12))
    sc = [dict(V=pl[0], F=pl[1], N=pl[2], col=PLATE), dict(V=dsk[0], F=dsk[1], N=dsk[2], col=IOD), dict(V=hw[0], F=hw[1], N=hw[2], col=NYL),
          dict(V=gas[0], F=gas[1], N=gas[2], col=GASK), dict(V=brk[0], F=brk[1], N=brk[2], col=BRK), dict(V=sup[0], F=sup[1], N=sup[2], col=SUP),
          dict(V=cont[0], F=cont[1], N=cont[2], col=ACRYL, glass=True), dict(V=lid[0], F=lid[1], N=lid[2], col=LIDC, glass=True)]
    img = down(render(sc, cam), SS)
    annotate(img, cam, [((60, -60, p['H']), 'flat top 6 mm, bottom 4 mm (as the current L2)\n(uniform thickness out to R %.0f mm)' % p['R_flat'], 'top'),
                        ((p['R_flat'] + 18, -75, p['H'] - 6), 'rim curved with a constant radius (%.2f mm);\nwall thins to %.1f mm (skin-equivalent)' % (p['a_out'], p['t_edge']), 'top'),
                        ((70, 10, ci + 1.5), 'iodine plate at mid-height\n(measured July layout)', 'top'),
                        ((-10, -YE + 3, c), '%.0f mm extension beyond the support face:\nbolting zone and plate slots' % p['L_ext'], 'bottom'),
                        ((x_lid_out, 20, 1), 'plain chest-wall lid, %.0f mm, 18 × M3 nylon/PEEK,\n1 mm FKM gasket, 4 fill/vent ports' % p['t_lid'], 'bottom'),
                        ((-p['lip_t'], 40, -p['lip_h'] + 2), 'positioning bracket: lip hooks against the\nsupport face (removable, ±%.0f mm adjustable)' % p['brk_adjust'], 'bottom'),
                        ((120, -150, -20), 'breast support (illustrative)', 'bottom')],
             os.path.join(OUT, 'L2v2_render_assembled.png'), 'L2 phantom v2 (rev B) on the breast support — container and lid shown as clear acrylic', SS, lines=[EDGE])
    # 2) exploded: lid and gasket pulled off, plate half out
    dxl, dxg, dxp = -95.0, -55.0, -75.0
    cam2 = Cam((-330, -330, 250), (25, 4, 10), W, H, fov=36)
    cam2.fit(np.r_[PTS1, box(x_lid_out + dxl, x_lid_in + dxl, -YE, YE, 0, p['H']), box(x_open + dxp - 1, R_plate + dxp, -R_plate, R_plate, ci - 2, ci + 2), box(x_open, 0, -85, 85, -p['lip_h'], 0)], margin=(0.04, 0.04, 0.15, 0.13))
    hw2 = hardware(dxl)
    sc = [dict(V=pl[0] + np.float32([dxp, 0, 0]), F=pl[1], N=pl[2], col=PLATE), dict(V=dsk[0] + np.float32([dxp, 0, 0]), F=dsk[1], N=dsk[2], col=IOD),
          dict(V=hw2[0], F=hw2[1], N=hw2[2], col=NYL), dict(V=gas[0] + np.float32([dxg, 0, 0]), F=gas[1], N=gas[2], col=GASK),
          dict(V=cont[0], F=cont[1], N=cont[2], col=ACRYL, glass=True), dict(V=lid[0] + np.float32([dxl, 0, 0]), F=lid[1], N=lid[2], col=LIDC, glass=True)]
    img = down(render(sc, cam2), SS)
    annotate(img, cam2, [((dxl + x_lid_out, -60, c), 'plain lid (no rib)', 'bottom'),
                         ((dxg - 0.5, 90, p['H'] - 3), 'FKM gasket 1 mm', 'top'),
                         ((dxp + x_open + 2, -125, ci), 'new plate: half-disc + full-width tab;\nthe tab runs in the slots', 'bottom'),
                         ((-12, 126, ci), 'slot in the thickened end wall\n(in the extension, outside the field)', 'top'),
                         ((x_open + 3, -60, p['H'] - 2), 'bolting zone, %.0f mm beyond the support face' % (-(x_open + p['d_frame'])), 'top'),
                         ((p['R_flat'] + 20, 40, ci), 'curved plate edge sits\nin the narrowing rim', 'top')],
             os.path.join(OUT, 'L2v2_render_exploded.png'), 'L2 phantom v2 (rev B) — exploded: lid, gasket and plate removed from the chest-wall side', SS)
    # 3) cut-away at y = 0 with spheres (illustrative packing) on the far side
    print('packing spheres', flush=True)
    sp = spheres_rsa(); print('  %d spheres' % len(sp), flush=True); smesh = sphere_mesh(sp)
    clip = lambda x, y, z: -y
    contc = mesh(f_container, PARTS['container'][1], 0.6, clip=clip); lidc = mesh(f_lid, PARTS['lid'][1], 0.6, clip=clip)
    brkc = mesh(f_bracket, PARTS['bracket'][1], 0.4, clip=clip); supc = boxmesh(0.0, R_out + 25, 0.0, 160.0, -42.0, -0.05)
    plc = mesh(f_plate, PARTS['plate'][1], 0.5, clip=clip)
    dl = [d for d in disc_layout() if d['y'] > 0]
    dskc = merge([cylinder((d['x'], d['y'], ci - p['t_plate'] / 2 - 0.06), (0, 0, 1), d['diam'] / 2, p['t_plate'] + 0.12, 32) for d in dl])
    cam3 = Cam((40, -300, 120), (50, 40, 18), W, H, fov=34); cam3.fit(box(x_lid_out, R_out, 0, R_out, -p['lip_h'] - 8, p['H']), margin=(0.06, 0.06, 0.2, 0.15))
    sc = [dict(V=contc[0], F=contc[1], N=contc[2], col=ACRYL, facecol=cut_colours(*contc)),
          dict(V=lidc[0], F=lidc[1], N=lidc[2], col=LIDC, facecol=cut_colours(*lidc, base=LIDC)),
          dict(V=plc[0], F=plc[1], N=plc[2], col=PLATE, facecol=cut_colours(*plc, col=(0.75, 0.75, 0.72), base=PLATE)),
          dict(V=dskc[0], F=dskc[1], N=dskc[2], col=IOD), dict(V=smesh[0], F=smesh[1], N=smesh[2], col=SPH),
          dict(V=brkc[0], F=brkc[1], N=brkc[2], col=BRK, facecol=cut_colours(*brkc, col=(0.45, 0.5, 0.56), base=BRK)), dict(V=supc[0], F=supc[1], N=supc[2], col=SUP)]
    img = down(render(sc, cam3), SS)
    EDGE3 = [np.array([[0, 0, -p['lip_h'] - 4], [0, 0, p['H'] + 6]]), '#c0392b', 'x = 0: support chest-wall face']
    annotate(img, cam3, [((p['R_flat'] + 22, 0, c + 8), 'curved rim: thickness falls to zero,\nlike a compressed breast edge', 'bottom'),
                         ((60, 0, ci + p['t_plate'] / 2), 'iodine plate at mid-height', 'top'),
                         ((x_lid_out + 2, 0, c - 10), 'plain lid', 'bottom'),
                         ((-p['lip_t'] / 2, 0, -p['lip_h'] + 3), 'positioning lip against the support face', 'bottom'),
                         ((-10, 0, ci + p['t_plate'] / 2), 'plate tab (outside the field)', 'top'),
                         ((p['R_flat'] - 40, 0, p['H'] - 3), '6 mm top wall (4 mm bottom)', 'top'),
                         ((x_open + 3, 0, p['H'] - 6), 'bolting zone', 'top'),
                         ((40, 30, 12), 'spheres and oil continue to the field edge and beyond\n(illustrative packing, not 61 %)', 'bottom')],
             os.path.join(OUT, 'L2v2_render_cutaway.png'), 'L2 phantom v2 (rev B) — section at the centre line (spheres drawn for illustration)', SS, lines=[EDGE3])
    print('done')

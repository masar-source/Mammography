#!/usr/bin/env python3
"""L2 phantom container v2 (rev B): parametric geometry (implicit solids) and design checks.
Walls (22 Oct 2026): top 6 mm and bottom 4 mm as measured on the current L2; rim wall 1.5 mm at the cavity equator.

Rev B: the container continues as a straight extension beyond the chest-wall edge. The bolting zone, the lid and
the plate slots sit in that extension, outside the field of view; the lid is a plain plate. A removable L-shaped
positioning bracket under the extension hooks against the chest-wall face of the breast support.

Coordinates (mm): x = 0 at the chest-wall face of the breast support (= contact face of the positioning lip); the image
starts at x = d_miss >= 0 (missed tissue, to be measured); x > 0 into the phantom, x < 0 = extension, lid and bracket
(outside the field); y = along the chest wall (0 = centre line); z = up from the bottom of the phantom (support surface).
Every solid is an implicit function f(x,y,z) <= 0 inside, built to be approximately a signed distance.
Usage: python l2v2_geom.py [resolution_mm] [--mesh]   (writes params.json, checks.json, disc_layout.csv; --mesh: preview STLs)
"""
import json, os, sys, csv, numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

# ------------------------------------------------------------------ parameters
P = dict(
    H=53.5,           # outer height (confirmed: outer)
    t_top=5.0,        # top wall: equal walls (user, 27 Oct 2026); same 10 mm total as the current L2 (lid 6, base 4 measured)
    t_bot=5.0,        # bottom wall
    t_edge=1.5,       # rim wall at the cavity equator: PMMA equivalent of 1.45 mm breast skin (Huang 2008), LE/HE/recombined compromise
    a_out=26.75,      # radial semi-axis of the outer rim profile; = H/2 gives a constant-radius rim
    R_flat=100.0,     # radius of the flat (uniform-thickness) region; outer radius = R_flat + a_out
    L_ext=20.0,       # straight extension beyond the detector chest-wall edge (outside the field of view)
    t_side_ext=10.0,   # extra wall at the rounded ends of the extension (carries the plate slots and end screws)
    d_frame=6.0,      # depth of the bolting zone at the open end of the extension
    ft_top=2.0,       # inward thickening of the top wall in the bolting zone (7 mm there)
    ft_bot=3.0,       # inward thickening of the bottom wall in the bolting zone (8 mm there; carries the bracket taps)
    t_lid=6.0,        # chest-wall lid (plain plate)
    t_gasket=1.0,     # flat FKM gasket
    t_plate=2.0,      # iodine plate thickness: 2 mm PMMA, as the July plate (user, 27 Oct 2026; was 3.0 assumed)
    plate_gap=0.4,    # radial gap between the plate's curved edge and the cavity equator
    slot_depth=3.0,   # plate slots cut into the thick end walls of the extension, beyond the cavity
    slot_clear=0.2,   # vertical clearance per side in the slots
    tab_protrusion=0.8,  # plate tab end protrudes into the gasket gap so the plain lid presses it
    slot_cap=5.0,     # outer lug closing each slot end: the slots run to x = 0 and reach beyond the imaged part's outer
                      # surface there, so without it they open through the step face (leak); lug length beyond x = 0
    slot_cap_margin=3.0,  # lug wall around the slot mouth (beyond the slot's outer face, above and below it)
    screw_tap=2.5,    # M3 tap drill
    screw_clear=3.4,  # M3 clearance hole in lid and gasket
    screw_depth=5.0,  # blind thread depth (must stay < d_frame)
    port_d=6.8,       # M8 tap drill, fill/vent ports in the lid
    brk_w=170.0,      # positioning bracket: width along the chest wall
    brk_t=4.0,        # thickness of its horizontal leg (under the extension)
    lip_t=5.0,        # thickness of the lip (x)
    lip_h=16.0,       # depth of the lip below the phantom bottom
    brk_x=-17.0,      # x of the bracket screws (in the bolting zone, bottom wall 7 mm thick)
    brk_depth=5.0,    # blind thread depth of the bracket taps (leaves 2 mm below the cavity)
    brk_adjust=2.0,   # +/- adjustment of the bracket in x (slotted holes)
)
p = P
c = p['H'] / 2                                   # centre of the outer profile
b_out = p['H'] / 2; a_out = p['a_out']
ci = (p['t_bot'] + p['H'] - p['t_top']) / 2      # centre of the cavity (mid-height of the sphere bed; plate plane)
b_in = (p['H'] - p['t_top'] - p['t_bot']) / 2
def _min_wall(a_i, a_o=None, c_o=None, b_o=None, n=1200):
    """shortest distance between the cavity profile (a_i, b_in, centre ci) and an outer profile (default: the imaged rim)"""
    a_o = a_out if a_o is None else a_o; c_o = c if c_o is None else c_o; b_o = b_out if b_o is None else b_o
    to = np.linspace(-np.pi / 2, np.pi / 2, 3 * n); O = np.c_[a_o * np.cos(to), c_o + b_o * np.sin(to)]
    ti = np.linspace(-np.pi / 2, np.pi / 2, n); I = np.c_[a_i * np.cos(ti), ci + b_in * np.sin(ti)]
    return float(min(np.min(np.hypot(*(O - q).T)) for q in I[::4]))
a_in = a_out * np.sqrt(1 - ((ci - c) / b_out) ** 2) - p['t_edge']
for _ in range(30):                               # solve: thinnest rim wall = t_edge exactly
    a_in -= p['t_edge'] - _min_wall(a_in)
a_ext = a_out + p['t_side_ext']                  # outer semi-axis of the extension's rounded ends
R_out = p['R_flat'] + a_out                      # outer radius of the imaged part
R_cav = p['R_flat'] + a_in                       # cavity radius at mid-height
R_plate = R_cav - p['plate_gap']
x_open = -p['L_ext']                             # open end of the extension (frame face)
x_lid_in = x_open - p['t_gasket']; x_lid_out = x_lid_in - p['t_lid']
y_slot_in, y_slot_out = R_cav - 1.0, R_cav + p['slot_depth']
P_tab = y_slot_out - 0.3                         # half-width of the plate tab that runs in the slots
z_top_open = p['H'] - p['t_top'] - p['ft_top']; z_bot_open = p['t_bot'] + p['ft_bot']   # opening of the bolting zone
c_open = (z_top_open + z_bot_open) / 2; b_open = (z_top_open - z_bot_open) / 2
a_gin = a_in + p['slot_depth'] + 0.2             # gasket inner boundary (encloses the slot mouths)
x_tab_end = x_open - p['tab_protrusion']; x_tab_stop = -0.3
x_cap_end = p['slot_cap']                        # outer lugs: x 0..x_cap_end, |y| y_cap_in..y_cap_out, |z - ci| <= z_cap_half
y_cap_in = R_cav + 0.3                           # inside the rim wall (cavity <= R_cav, outer surface >= 126.1 mm here)
y_cap_out = y_slot_out + p['slot_cap_margin']
z_cap_half = p['t_plate'] / 2 + p['slot_clear'] + p['slot_cap_margin']

# screws: 7 per straight edge (centre of the thickened top/bottom wall) + 2 per rounded end
z_s_bot = (p['t_bot'] + p['ft_bot']) / 2; z_s_top = p['H'] - (p['t_top'] + p['ft_top']) / 2
a_s, b_s, c_s = (a_in + a_ext) / 2, (b_in + b_out) / 2, (c + ci) / 2
SCREWS = [(y, z) for y in (-90, -60, -30, 0, 30, 60, 90) for z in (z_s_bot, z_s_top)]
def _end_screw(lo, hi):
    """position in the extension's end wall (between cavity and outer profile) with the most material around it"""
    to = np.linspace(-np.pi / 2, np.pi / 2, 1500); O = np.c_[a_ext * np.cos(to), c + b_out * np.sin(to)]
    ti = np.linspace(-np.pi / 2, np.pi / 2, 1500); I = np.c_[a_in * np.cos(ti), ci + b_in * np.sin(ti)]
    GO = np.c_[a_gin * np.cos(ti), c_open + b_open * np.sin(ti)]                                  # gasket opening
    best = None
    for th in np.radians(np.linspace(lo, hi, 36)):
        pi_ = np.array([a_in * np.cos(th), ci + b_in * np.sin(th)]); po = np.array([a_ext * np.cos(th), c + b_out * np.sin(th)])
        for f in np.linspace(0.2, 0.8, 25):
            q = pi_ + f * (po - pi_); d = min(np.min(np.hypot(*(O - q).T)), np.min(np.hypot(*(I - q).T)))
            if abs(q[1] - ci) < p['t_plate'] / 2 + p['slot_clear'] + p['screw_tap']: continue     # keep clear of the slot
            if np.min(np.hypot(*(GO - q).T)) < p['screw_clear'] / 2 + 2.0: continue              # keep a 2 mm gasket web
            if best is None or d > best[0]: best = (d, q)
    return best[1]
for s in (-1, 1):
    for lo, hi in ((-65, -20), (20, 65)):
        qy, qz = _end_screw(lo, hi); SCREWS.append((s * (p['R_flat'] + qy), qz))
PORTS = [(y, z) for y in (-50, 50) for z in (ci - 10, ci + 10)]
BRK_SCREWS = [(p['brk_x'], y) for y in (-75, -45, 45, 75)]   # vertical M3, between the lid screws

def disc_layout():
    """measured July 2026 layout: x from the outer chest-wall face of the July container, kept unchanged from the support face"""
    D = json.load(open(os.path.join(HERE, 'disc_layout_measured.json')))['discs']
    return [dict(d) for d in D]

# ------------------------------------------------------------------ implicit primitives
def ell(q, w, a, b):
    """approx. signed distance to an ellipse (semi-axes a radial, b vertical) for q = radial excess >= 0"""
    return (np.sqrt((q / a) ** 2 + (w / b) ** 2) - 1.0) * min(a, b)
def stadium(y, z, a, b, zc=None):
    """(y,z) cross-section: flat for |y| <= R_flat, elliptical ends with semi-axes (a, b), centred at height zc"""
    return ell(np.maximum(np.abs(y) - p['R_flat'], 0), z - (c if zc is None else zc), a, b)
def f_outer(x, y, z):
    rho = np.hypot(x, y)
    field = np.maximum(ell(np.maximum(rho - p['R_flat'], 0), z - c, a_out, b_out), -x - 0.5)   # 0.5 mm overlap with the extension
    ext = np.maximum.reduce([stadium(y, z, a_ext, b_out), x, x_open - x])
    return np.minimum(field, ext)
def f_cavity(x, y, z):
    rho = np.hypot(x, y)
    field = np.maximum(ell(np.maximum(rho - p['R_flat'], 0), z - ci, a_in, b_in), -x - 0.5)
    ext = np.maximum.reduce([stadium(y, z, a_in, b_in, ci), x, (x_open - 1.0) - x])
    return np.minimum(field, ext)
def f_slot(x, y, z):
    ay = np.abs(y)
    return np.maximum.reduce([ay - y_slot_out, y_slot_in - ay, np.abs(z - ci) - (p['t_plate'] / 2 + p['slot_clear']), x, (x_open - 1.0) - x])
def f_cap(x, y, z):
    """outer lugs on the step face at x = 0 that close the slot ends (from 0.5 mm inside the extension wall, as in the
    STEP model; the slots are cut after them)"""
    ay = np.abs(y)
    return np.maximum.reduce([ay - y_cap_out, y_cap_in - ay, np.abs(z - ci) - z_cap_half, -0.5 - x, x - x_cap_end])
def holes2d(y, z, pts, r):
    d = np.full(np.broadcast(y, z).shape, 1e9, dtype=float)
    for (yi, zi) in pts:
        d = np.minimum(d, np.hypot(y - yi, z - zi) - r)
    return d

def f_body(x, y, z):
    """container material without the tapped holes"""
    shell = np.maximum(f_outer(x, y, z), -f_cavity(x, y, z))
    band = np.maximum(z - z_top_open, z_bot_open - z)          # <= 0 inside the opening
    frame = np.maximum.reduce([stadium(y, z, a_in, b_in, ci), -band, x - (x_open + p['d_frame']), x_open - x])
    return np.maximum(np.minimum(np.minimum(shell, frame), f_cap(x, y, z)), -f_slot(x, y, z))
def holes_xy(x, y, pts, r):
    d = np.full(np.broadcast(x, y).shape, 1e9, dtype=float)
    for (xi, yi) in pts:
        d = np.minimum(d, np.hypot(x - xi, y - yi) - r)
    return d
def f_container(x, y, z):
    taps = np.maximum.reduce([holes2d(y, z, SCREWS, p['screw_tap'] / 2), x - (x_open + p['screw_depth']), (x_open - 1.0) - x])
    btaps = np.maximum.reduce([holes_xy(x, y, BRK_SCREWS, p['screw_tap'] / 2), z - p['brk_depth'], -1.0 - z])
    return np.maximum.reduce([f_body(x, y, z), -taps, -btaps])
def f_bracket(x, y, z):
    """L-shaped positioning bracket: horizontal leg under the extension, lip whose inner face (x = 0) hooks against the support"""
    leg = np.maximum.reduce([x - 0.0, x_open - x, np.abs(y) - p['brk_w'] / 2, z - 0.0, -p['brk_t'] - z])
    lip = np.maximum.reduce([x - 0.0, -p['lip_t'] - x, np.abs(y) - p['brk_w'] / 2, z - 0.0, -p['lip_h'] - z])
    a = p['brk_adjust']; r = p['screw_clear'] / 2; slot = np.full(np.broadcast(x, y).shape, 1e9)
    for (xi, yi) in BRK_SCREWS:      # slotted clearance holes along x (+/- adjustment)
        slot = np.minimum(slot, np.hypot(np.maximum(np.abs(x - xi) - a, 0), y - yi) - r)
    return np.maximum(np.minimum(leg, lip), -slot)
def f_lid(x, y, z):
    body = np.maximum.reduce([stadium(y, z, a_ext, b_out), x - x_lid_in, x_lid_out - x])
    return np.maximum.reduce([body, -holes2d(y, z, SCREWS, p['screw_clear'] / 2), -holes2d(y, z, PORTS, p['port_d'] / 2)])
def f_gasket(x, y, z):
    ring = np.maximum.reduce([stadium(y, z, a_ext, b_out), -stadium(y, z, a_gin, b_open, c_open), x - x_open, x_lid_in - x])
    return np.maximum(ring, -holes2d(y, z, SCREWS, p['screw_clear'] / 2))
def f_plate(x, y, z):
    disc = np.maximum(np.hypot(x, y) - R_plate, x_tab_stop - x)
    tab = np.maximum.reduce([np.abs(y) - P_tab, x - x_tab_stop, x_tab_end - x])
    return np.maximum(np.abs(z - ci) - p['t_plate'] / 2, np.minimum(disc, tab))

YB = p['R_flat'] + a_ext + 0.6
PARTS = {  # name: (function, bounding box)
    'container': (f_container, (x_open - 1.6, R_out + 0.6, -YB, YB, -0.6, p['H'] + 0.6)),
    'lid':       (f_lid, (x_lid_out - 0.6, x_lid_in + 0.6, -YB, YB, -0.6, p['H'] + 0.6)),
    'gasket':    (f_gasket, (x_lid_in - 0.6, x_open + 0.6, -YB, YB, -0.6, p['H'] + 0.6)),
    'plate':     (f_plate, (x_tab_end - 0.6, R_plate + 0.6, -P_tab - 0.6, P_tab + 0.6, ci - 2.5, ci + 2.5)),
    'bracket':   (f_bracket, (x_open - 0.6, 0.6, -p['brk_w'] / 2 - 0.6, p['brk_w'] / 2 + 0.6, -p['lip_h'] - 0.6, 0.6)),
}

# ------------------------------------------------------------------ meshing (preview only; exact solids: l2v2_step.py)
def sample(f, bb, h, clip=None):
    xs = np.arange(bb[0], bb[1] + h, h); ys = np.arange(bb[2], bb[3] + h, h); zs = np.arange(bb[4], bb[5] + h, h)
    vol = np.empty((len(xs), len(ys), len(zs)), np.float32)
    Y, Z = np.meshgrid(ys, zs, indexing='ij')
    for i, x in enumerate(xs):
        v = f(np.full_like(Y, x), Y, Z)
        if clip is not None: v = np.maximum(v, clip(np.full_like(Y, x), Y, Z))
        vol[i] = v
    vol[0] = vol[-1] = 1.0; vol[:, 0] = vol[:, -1] = 1.0; vol[:, :, 0] = vol[:, :, -1] = 1.0   # close the mesh
    return vol, (xs[0], ys[0], zs[0])
def mesh(f, bb, h, clip=None, shift=(0, 0, 0)):
    """marching cubes; winding made outward (checked by stepping along face normals); smooth vertex normals"""
    from skimage import measure
    vol, o = sample(f, bb, h, clip)
    V, F, _, _ = measure.marching_cubes(vol, 0.0, spacing=(h, h, h))
    V = V + np.array(o)
    fn = np.cross(V[F[:, 1]] - V[F[:, 0]], V[F[:, 2]] - V[F[:, 0]]); fn /= np.linalg.norm(fn, axis=1, keepdims=True) + 1e-12
    k = np.random.default_rng(0).choice(len(F), min(2000, len(F)), replace=False); cen = V[F[k]].mean(1); e = 0.3 * h
    fx = lambda P_: (f(P_[:, 0], P_[:, 1], P_[:, 2]) if clip is None else np.maximum(f(P_[:, 0], P_[:, 1], P_[:, 2]), clip(P_[:, 0], P_[:, 1], P_[:, 2])))
    if np.mean(fx(cen + e * fn[k]) > fx(cen - e * fn[k])) < 0.5:
        F = F[:, ::-1]; fn = -fn
    VN = np.zeros_like(V); np.add.at(VN, F.ravel(), np.repeat(fn, 3, axis=0)); VN /= np.linalg.norm(VN, axis=1, keepdims=True) + 1e-12
    return (V + np.array(shift)).astype(np.float32), F.astype(np.int64), VN.astype(np.float32)
def _prism(poly, z0, z1):
    """closed prism over a convex polygon (fan triangulation from its centroid)"""
    n = len(poly); cen = poly.mean(0)
    V = np.vstack([np.c_[poly, np.full(n, z0)], np.c_[poly, np.full(n, z1)], [[cen[0], cen[1], z0], [cen[0], cen[1], z1]]])
    F = []
    for i in range(n):
        j = (i + 1) % n
        F += [[2 * n, j, i], [2 * n + 1, n + i, n + j], [i, j, n + j], [i, n + j, n + i]]
    return V, np.array(F)
def plate_mesh(n_arc=400):
    """plate for rendering: half-disc (R_plate) and the full-width tab as two convex prisms"""
    z0, z1 = ci - p['t_plate'] / 2, ci + p['t_plate'] / 2
    t = np.linspace(-np.pi / 2, np.pi / 2, n_arc)
    disc = np.r_[[[x_tab_stop, -R_plate]], np.c_[np.maximum(R_plate * np.cos(t), x_tab_stop), R_plate * np.sin(t)][1:-1], [[x_tab_stop, R_plate]]]
    tab = np.array([[x_tab_end, -P_tab], [x_tab_stop, -P_tab], [x_tab_stop, P_tab], [x_tab_end, P_tab]])
    V1, F1 = _prism(disc, z0, z1); V2, F2 = _prism(tab, z0, z1)
    V = np.vstack([V1, V2]); F = np.vstack([F1, F2 + len(V1)])
    fn = np.cross(V[F[:, 1]] - V[F[:, 0]], V[F[:, 2]] - V[F[:, 0]]); fn /= np.linalg.norm(fn, axis=1, keepdims=True) + 1e-12
    VN = np.zeros_like(V); np.add.at(VN, F.ravel(), np.repeat(fn, 3, axis=0)); VN /= np.linalg.norm(VN, axis=1, keepdims=True) + 1e-12
    return V.astype(np.float32), F, VN.astype(np.float32)
def write_stl(path, V, F):
    tri = V[F]; n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]); n /= np.linalg.norm(n, axis=1, keepdims=True) + 1e-12
    rec = np.zeros(len(F), dtype=[('n', '<f4', 3), ('v', '<f4', (3, 3)), ('a', '<u2')]); rec['n'] = n; rec['v'] = tri
    with open(path, 'wb') as fh:
        fh.write(b'L2 phantom container v2 (mm)'.ljust(80, b' ')); fh.write(np.uint32(len(F)).tobytes()); fh.write(rec.tobytes())

# ------------------------------------------------------------------ checks
def slot_ends_sealed(h=0.1):
    """True if at both slot ends the cavity and the outside are separate void regions (local flood fill around x = 0,
    where the slots reach beyond the imaged part's outer surface)"""
    from scipy import ndimage as ndi
    o = 0.137 * h
    xs = np.arange(-4.0, x_cap_end + 3.0, h) + o; ys = np.arange(y_slot_in - 3.0, y_cap_out + 3.0, h) + o
    zs = np.arange(ci - z_cap_half - 3.0, ci + z_cap_half + 3.0, h) + o
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing='ij'); ok = []
    for s in (1, -1):
        lab, _ = ndi.label(f_container(X, s * Y, Z) > 0)
        cav = lab[np.argmin(abs(xs + 3.0)), np.argmin(abs(ys - (y_slot_in - 2.0))), np.argmin(abs(zs - ci))]   # in the cavity
        outside = lab[-1, -1, np.argmin(abs(zs - ci))]                                                       # beyond the lug
        ok.append(bool(cav and outside and cav != outside))
    return all(ok)

def checks():
    out = {}
    out['slot_ends_sealed'] = slot_ends_sealed()
    h = 0.5; o3 = 0.137 * h          # off-centre grid: cell centres never coincide with plate/wall faces
    xs = np.arange(x_open - 1, R_out, h) + o3; ys = np.arange(-YB, YB, h) + o3; zs = np.arange(0, p['H'], h) + o3
    Y, Z = np.meshgrid(ys, zs, indexing='ij'); vin = vext = 0
    for x in xs:
        X = np.full_like(Y, x)
        fill = (f_cavity(X, Y, Z) < 0) & (f_container(X, Y, Z) > 0) & (f_plate(X, Y, Z) > 0) & (X > x_open)
        vin += fill.sum(); vext += fill[:, :].sum() if x < 0 else 0
    out['filling_volume_ml'] = float(vin * h ** 3 / 1000); out['filling_in_extension_ml'] = float(vext * h ** 3 / 1000)
    out['outer_radius_imaged_part_mm'] = R_out; out['outer_width_imaged_part_mm'] = 2 * R_out
    out['extension_width_mm'] = 2 * (p['R_flat'] + a_ext); out['overall_depth_mm'] = R_out - x_lid_out
    out['overhang_beyond_support_face_mm'] = -x_lid_out
    out['bolting_zone_distance_from_support_face_mm'] = -(x_open + p['d_frame'])
    out['cavity_radius_mid_mm'] = R_cav; out['plate_radius_mm'] = R_plate; out['plate_tab_halfwidth_mm'] = P_tab
    zz = np.linspace(0, p['H'], 400); rr = p['R_flat'] + a_out * np.sqrt(np.clip(1 - ((zz - c) / b_out) ** 2, 0, 1))
    M = 656 / (656 - 17 - zz); out['projected_width_mm'] = float(2 * np.max(rr * M)); out['field_of_view_mm'] = 291
    D = disc_layout()
    out['outermost_disc_edge_radius_mm'] = float(max(np.hypot(d['x'], d['y']) + d['diam'] / 2 for d in D))
    out['flat_region_margin_mm'] = p['R_flat'] - out['outermost_disc_edge_radius_mm']
    out['closest_disc_edge_to_support_face_mm'] = float(min(d['x'] - d['diam'] / 2 for d in D))
    # plate insertion along +x: every plate point, shifted back by s, must keep clear of the container walls
    gx = np.linspace(x_tab_end, R_plate, 420); gy = np.linspace(-P_tab, P_tab, 520)
    XX, YY = np.meshgrid(gx, gy, indexing='ij'); m = f_plate(XX, YY, np.full_like(XX, ci)) <= 0
    worst = 1e9
    for s in np.linspace(0, R_plate - x_open + 2, 160):
        X = XX[m] - s; Yp = YY[m]; inside = X > x_open - 1
        for zoff in (-p['t_plate'] / 2, 0.0, p['t_plate'] / 2):
            if inside.any(): worst = min(worst, float(f_container(X[inside], Yp[inside], np.full(inside.sum(), ci + zoff)).min()))
    out['plate_insertion_min_clearance_mm'] = worst
    out['plate_edge_cavity_halfheight_mm'] = float(b_in * np.sqrt(1 - ((R_plate - p['R_flat']) / a_in) ** 2))
    # wall left around the tapped holes: largest ring (beyond the hole radius) that is still entirely container material
    xs_h = np.linspace(x_open + 0.2, x_open + p['screw_depth'] - 0.2, 6); th = np.linspace(0, 2 * np.pi, 120, endpoint=False); walls = []
    for (yi, zi) in SCREWS:
        wall = 0.0
        for dr in np.arange(0.05, 6.0, 0.05):
            r = p['screw_tap'] / 2 + dr; Yc = yi + r * np.cos(th); Zc = zi + r * np.sin(th)
            if any((f_body(np.full_like(Yc, xh), Yc, Zc) > 0).any() for xh in xs_h): break
            wall = dr
        walls.append(wall)
    out['min_wall_around_tapped_holes_mm'] = float(min(walls)); out['wall_around_tapped_holes_mm'] = [round(w, 2) for w in walls]
    # bracket screws: material left above the blind holes and clearance to the lid-screw holes
    out['wall_above_bracket_taps_mm'] = float(p['t_bot'] + p['ft_bot'] - p['brk_depth'])
    out['bracket_tap_to_lid_tap_mm'] = float(min(np.hypot(yb - ys, p['brk_depth'] / 2 - zs) for (_, yb) in BRK_SCREWS for (ys, zs) in SCREWS
                                                if zs < c) - p['screw_tap'])
    out['bracket_taps_inside_bolting_zone'] = bool(all(x_open + p['screw_tap'] / 2 < xb < x_open + p['d_frame'] - p['screw_tap'] / 2 for (xb, _) in BRK_SCREWS))
    # rim wall: shortest distance from the inner to the outer profile (exact, 2D), sampled around the rim
    to = np.linspace(-np.pi / 2, np.pi / 2, 4001); OUTP = np.c_[a_out * np.cos(to), c + b_out * np.sin(to)]
    ti = np.linspace(-np.pi / 2, np.pi / 2, 181); INP = np.c_[a_in * np.cos(ti), ci + b_in * np.sin(ti)]
    wall = np.array([np.min(np.hypot(*(OUTP - q).T)) for q in INP])
    out['rim_wall_bottom_equator_top_mm'] = [float(wall[0]), float(wall[90]), float(wall[-1])]; out['rim_wall_min_mm'] = float(wall.min())
    ti = np.linspace(-np.pi / 2, np.pi / 2, 3000); GO = np.c_[p['R_flat'] + a_gin * np.cos(ti), c_open + b_open * np.sin(ti)]
    out['min_gasket_web_mm'] = float(min(np.min(np.hypot(GO[:, 0] - abs(ys), GO[:, 1] - zs)) for (ys, zs) in SCREWS if abs(ys) > p['R_flat']) - p['screw_clear'] / 2)
    out['n_screws'] = len(SCREWS); out['n_ports'] = len(PORTS)
    return out

if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if not a.startswith('--')]; h = float(args[0]) if args else 0.5
    json.dump(dict(P, derived=dict(R_out=R_out, R_cav=R_cav, R_plate=R_plate, a_in=a_in, b_in=b_in, ci=ci, a_ext=a_ext, x_open=x_open,
                                   x_lid_in=x_lid_in, x_lid_out=x_lid_out, y_slot_in=y_slot_in, y_slot_out=y_slot_out, P_tab=P_tab,
                                   z_top_open=z_top_open, z_bot_open=z_bot_open, c_open=c_open, b_open=b_open, a_gin=a_gin,
                                   x_tab_end=x_tab_end, x_tab_stop=x_tab_stop, x_cap_end=x_cap_end, y_cap_in=y_cap_in,
                                   y_cap_out=y_cap_out, z_cap_half=z_cap_half),
                   screws=SCREWS, ports=PORTS, bracket_screws=BRK_SCREWS), open(os.path.join(HERE, 'params.json'), 'w'), indent=1)
    with open(os.path.join(HERE, 'disc_layout.csv'), 'w', newline='') as fh:
        w = csv.writer(fh); w.writerow(['kind', 'row', 'col', 'diameter_mm', 'x_mm_from_support_chest_wall_face', 'y_mm'])
        for d in disc_layout(): w.writerow([d['kind'], d['row'], d['col'], d['diam'], '%.1f' % d['x'], '%.1f' % d['y']])
    ck = checks(); json.dump(ck, open(os.path.join(HERE, 'checks.json'), 'w'), indent=1)
    for k, v in ck.items(): print('%-44s %s' % (k, v if isinstance(v, bool) else np.round(v, 2)))
    if '--mesh' in sys.argv:          # optional preview meshes from the implicit model (exact STL/STEP: l2v2_step.py)
        os.makedirs(os.path.join(HERE, 'stl_implicit'), exist_ok=True)
        for name, (f, bb) in PARTS.items():
            V, F, N = plate_mesh() if name == 'plate' else mesh(f, bb, h)
            write_stl(os.path.join(HERE, 'stl_implicit', 'L2v2_%s.stl' % name), V, F)
            print('%-10s %7d triangles' % (name, len(F)), flush=True)

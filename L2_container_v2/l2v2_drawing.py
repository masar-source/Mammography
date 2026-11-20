#!/usr/bin/env python3
"""Dimensioned 2D drawing sheet and edge-profile comparison for the L2 v2 container (rev B, final state of 27 Oct 2026:
equal 5 / 5 mm walls, 2 mm plate, lugs closing the slot ends). All numbers come from l2v2_geom.py and checks.json.
Every label points at its feature with a leader that ends in an arrow head on the feature. After drawing, the sheet is
checked for collisions: no text may overlap another text, a line, a leader or an iodine disc; the result is printed.
Usage: python l2v2_drawing.py -> drawings/L2v2_drawing.pdf/.png, drawings/L2v2_edge_profile.png"""
import os, json, numpy as np
from l2v2_geom import (p, c, ci, z_top_open, z_bot_open, c_open, a_out, b_out, a_in, b_in, a_ext, R_out, R_cav, R_plate, x_open, x_lid_in, x_lid_out, y_slot_in, y_slot_out,
                       P_tab, b_open, a_gin, x_tab_end, x_tab_stop, SCREWS, PORTS, BRK_SCREWS, disc_layout, x_cap_end, y_cap_in, y_cap_out, z_cap_half)
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Circle, Rectangle
from matplotlib.text import Text
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'drawings'); os.makedirs(OUT, exist_ok=True)
INK, THIN, CL, HAT, RED = '#111111', '#555555', '#0b6bcb', '#c9dcef', '#c0392b'
FS = 7.2                                                  # label font size
plt.rcParams.update({'font.size': 8, 'font.family': 'DejaVu Sans'})
CK = json.load(open(os.path.join(HERE, 'checks.json')))   # written with the current geometry (the STEP checks are from 22 Oct)
VOL, VEXT = CK['filling_volume_ml'], CK['filling_in_extension_ml']
YE = p['R_flat'] + a_ext
_n = [0]


def owner(*artists, tag=None):
    """mark artists that belong together (a label and its leader, a dimension and its lines): the collision check
    does not compare them with each other"""
    if tag is None:
        _n[0] += 1; tag = 'o%d' % _n[0]
    for a in artists:
        a._owner = tag
    return tag


def dim(ax, p0, p1, txt, off=(0, 0), fs=FS, ext=True, t=0.5, gap=1.5):
    """dimension line between p0 and p1, moved by off, with extension lines that start 'gap' mm from the feature;
    the text sits above the line at the fraction t of its length"""
    p0, p1, o = np.array(p0, float), np.array(p1, float), np.array(off, float)
    q0, q1 = p0 + o, p1 + o; arts = []
    if ext and np.linalg.norm(o) > 0:
        u = o / np.linalg.norm(o)
        for a, b in ((p0, q0), (p1, q1)):
            arts += ax.plot([a[0] + gap * u[0], b[0] + 1.5 * u[0]], [a[1] + gap * u[1], b[1] + 1.5 * u[1]], color=THIN, lw=0.4)
    arts.append(ax.annotate('', q0, q1, arrowprops=dict(arrowstyle='<|-|>', lw=0.6, color=INK, mutation_scale=6, shrinkA=0, shrinkB=0)))
    m = q0 + t * (q1 - q0); ang = np.degrees(np.arctan2(q1[1] - q0[1], q1[0] - q0[0]))
    if ang > 90 or ang < -90:
        ang += 180
    if txt:
        arts.append(ax.text(m[0], m[1], txt, fontsize=fs, ha='center', va='bottom', rotation=ang, rotation_mode='anchor', color=INK,
                            bbox=dict(fc='white', ec='none', pad=0.3)))
    owner(*arts)


def radius(ax, ang_deg, r, txt, t=0.55, fs=FS):
    """radius dimension: line from the centre to the arc, arrow head on the arc, text along the line"""
    a = np.radians(ang_deg); q = np.array([r * np.cos(a), r * np.sin(a)]); arts = []
    arts.append(ax.annotate('', q, (0, 0), arrowprops=dict(arrowstyle='-|>', lw=0.6, color=INK, mutation_scale=7, shrinkA=0, shrinkB=0)))
    arts += ax.plot([0], [0], marker='+', color=INK, ms=5, mew=0.6)
    ang = ang_deg if -90 <= ang_deg <= 90 else ang_deg + 180
    arts.append(ax.text(*(t * q), txt, fontsize=fs, ha='center', va='bottom', rotation=ang, rotation_mode='anchor', color=INK, bbox=dict(fc='white', ec='none', pad=0.3)))
    owner(*arts)


def leader(ax, xy, txt, xytext, ha='left', fs=FS, col=INK):
    """label with a straight leader from the side of the text that faces the feature to an arrow head on the feature"""
    dx, dy = xy[0] - xytext[0], xy[1] - xytext[1]
    if abs(dy) > 1.2 * abs(dx):                           # leave from the top or bottom edge, near the anchored side
        rel = ({'left': 0.08, 'right': 0.92}.get(ha, 0.5), 1.0 if dy > 0 else 0.0)
    else:                                                 # leave from the left or right edge
        rel = (1.0 if dx > 0 else 0.0, 0.5)
    a = ax.annotate(txt, xy, xytext, fontsize=fs, color=col, ha=ha, va='center', zorder=6, bbox=dict(fc='white', ec='none', pad=0.4),
                    arrowprops=dict(arrowstyle='-|>', lw=0.55, color=THIN, mutation_scale=7, shrinkA=2, shrinkB=0, relpos=rel))
    owner(a); return a


def ellipse_pts(cx, cz, a, b, t0, t1, n=120):
    t = np.radians(np.linspace(t0, t1, n)); return np.c_[cx + a * np.cos(t), cz + b * np.sin(t)]


def stadium_yz(a, b, n=120, zc=None):
    zc = c if zc is None else zc
    return np.r_[ellipse_pts(p['R_flat'], zc, a, b, -90, 90, n), ellipse_pts(-p['R_flat'], zc, a, b, 90, 270, n)]


def arrow_paths(t):
    """the paths of an annotation's arrow in display coordinates (matplotlib returns them nested)"""
    out = []
    def walk(x):
        if hasattr(x, 'vertices'):
            out.append(x)
        elif isinstance(x, (list, tuple)):
            for y in x:
                walk(y)
    walk(t.arrow_patch._get_path_in_displaycoord()); return out


def collisions(fig):
    """every text against every other text, every line, leader and patch outline, and every iodine disc (display coordinates)"""
    fig.canvas.draw(); r = fig.canvas.get_renderer(); texts, segs, discs = [], [], []
    def add_path(path, tr, own):
        for poly in path.to_polygons(transform=tr, closed_only=False):
            for a, b in zip(poly[:-1], poly[1:]):
                segs.append((a, b, own))
    for ax in fig.axes:
        for t in ax.texts:
            if t.get_text().strip():
                texts.append((Text.get_window_extent(t, r), getattr(t, '_owner', id(t)), t.get_text().split('\n')[0][:40]))
            if getattr(t, 'arrow_patch', None) is not None:
                for pth in arrow_paths(t):
                    add_path(pth, None, getattr(t, '_owner', id(t)))
        for ln in ax.lines:
            xy = ax.transData.transform(ln.get_xydata())
            for a, b in zip(xy[:-1], xy[1:]):
                segs.append((a, b, getattr(ln, '_owner', id(ln))))
        for pa in ax.patches:
            if isinstance(pa, Circle) and pa.get_edgecolor()[3] == 0:
                ctr = ax.transData.transform(pa.center); rad = abs(ax.transData.transform((pa.radius, 0))[0] - ax.transData.transform((0, 0))[0])
                discs.append((ctr, rad)); continue
            if pa.get_edgecolor()[3] > 0 and pa.get_linewidth() > 0:
                add_path(pa.get_path(), pa.get_transform(), getattr(pa, '_owner', id(pa)))
    for t in fig.texts:
        texts.append((t.get_window_extent(r), id(t), t.get_text()[:40]))
    def seg_hits_box(a, b, bb):
        x0, y0, x1, y1 = bb.x0 + 1, bb.y0 + 1, bb.x1 - 1, bb.y1 - 1
        if x1 <= x0 or y1 <= y0:
            return False
        t0, t1 = 0.0, 1.0; d = b - a
        for pp, q in ((-d[0], a[0] - x0), (d[0], x1 - a[0]), (-d[1], a[1] - y0), (d[1], y1 - a[1])):
            if pp == 0:
                if q < 0: return False
            else:
                tt = q / pp
                if pp < 0: t0 = max(t0, tt)
                else: t1 = min(t1, tt)
                if t0 > t1: return False
        return True
    out = []
    W, H_ = fig.bbox.width, fig.bbox.height
    for bb, own, name in texts:                           # every text inside the sheet, 5 px from the edge
        if bb.x0 < 5 or bb.y0 < 5 or bb.x1 > W - 5 or bb.y1 > H_ - 5:
            out.append('text "%s" reaches the edge of the sheet' % name)
    for i, (bb, own, name) in enumerate(texts):
        for bb2, own2, name2 in texts[i + 1:]:
            if own != own2 and bb.x0 < bb2.x1 - 1 and bb2.x0 < bb.x1 - 1 and bb.y0 < bb2.y1 - 1 and bb2.y0 < bb.y1 - 1:
                out.append('text "%s" overlaps text "%s"' % (name, name2))
        hit = sum(1 for a, b, o in segs if o != own and seg_hits_box(np.asarray(a), np.asarray(b), bb))
        if hit:
            out.append('text "%s" is crossed by %d line segment(s)' % (name, hit))
        for ctr, rad in discs:
            cx, cy = min(max(ctr[0], bb.x0), bb.x1), min(max(ctr[1], bb.y0), bb.y1)
            if np.hypot(cx - ctr[0], cy - ctr[1]) < rad:
                out.append('text "%s" covers an iodine disc' % name); break
    for ax in fig.axes:                                   # leaders must not run over a disc they do not point at
        for t in ax.texts:
            if getattr(t, 'arrow_patch', None) is None or not t.get_text().strip():
                continue
            a = np.asarray(ax.transData.transform(t.xy))
            pieces = [(q0, q1) for pp in arrow_paths(t) for poly in pp.to_polygons(closed_only=False) for q0, q1 in zip(poly[:-1], poly[1:])]
            for ctr, rad in discs:
                if np.hypot(*(a - ctr)) <= rad + 2:
                    continue
                for q0, q1 in pieces:
                    d = q1 - q0; L = np.dot(d, d)
                    s = 0 if L == 0 else np.clip(np.dot(ctr - q0, d) / L, 0, 1)
                    if np.hypot(*(q0 + s * d - ctr)) < rad:
                        out.append('leader of "%s" runs over an iodine disc' % t.get_text().split('\n')[0][:40]); break
    return out


fig = plt.figure(figsize=(16.5, 11.7))
# ======================= PLAN VIEW
ax = fig.add_axes([0.01, 0.40, 0.47, 0.575]); ax.set_aspect('equal'); ax.set_axis_off()
ax.add_patch(Rectangle((x_lid_out - 12, -YE - 14), -x_lid_out + 12, 2 * YE + 28, fc='#ececec', ec='none', zorder=0))
t = np.linspace(-np.pi / 2, np.pi / 2, 400)
ax.plot(R_out * np.cos(t), R_out * np.sin(t), color=INK, lw=1.2)
ax.plot([x_open, 0, 0, x_open, x_open], [-YE, -YE, YE, YE, -YE], color=INK, lw=1.2)
ax.plot([0, 0], [-R_out, -YE], color=INK, lw=1.2); ax.plot([0, 0], [R_out, YE], color=INK, lw=1.2)
ax.plot(p['R_flat'] * np.cos(t), p['R_flat'] * np.sin(t), color=INK, lw=0.7, ls=(0, (5, 3)))
ax.plot(R_cav * np.cos(t), R_cav * np.sin(t), color=THIN, lw=0.5, ls=':')
for s_ in (-1, 1): ax.plot([x_open, 0], [s_ * R_cav] * 2, color=THIN, lw=0.5, ls=':')
ax.add_patch(Rectangle((x_lid_out, -YE), p['t_lid'], 2 * YE, fc=HAT, ec=INK, lw=0.8))
ax.add_patch(Rectangle((x_lid_in, -YE), p['t_gasket'], 2 * YE, fc='#333333', ec='none'))
ax.add_patch(Rectangle((x_open, -YE), p['d_frame'], 2 * YE, fc='none', ec=INK, lw=0.5, hatch='////'))
for s_ in (-1, 1):
    ax.add_patch(Rectangle((x_open, min(s_ * R_cav, s_ * y_slot_out)), -x_open, abs(y_slot_out - R_cav), fc='none', ec=INK, lw=0.6, hatch='xxxx'))
    ax.add_patch(Rectangle((0, min(s_ * y_cap_in, s_ * y_cap_out)), x_cap_end, y_cap_out - y_cap_in, fc=HAT, ec=INK, lw=0.6))     # lug closing the slot end
tt = np.linspace(-np.pi / 2, np.pi / 2, 300)
pl = np.r_[[[x_tab_end, -P_tab], [x_tab_stop, -P_tab], [x_tab_stop, -R_plate]], np.c_[np.maximum(R_plate * np.cos(tt), x_tab_stop), R_plate * np.sin(tt)],
           [[x_tab_stop, R_plate], [x_tab_stop, P_tab], [x_tab_end, P_tab], [x_tab_end, -P_tab]]]
ax.plot(pl[:, 0], pl[:, 1], color=CL, lw=1.0, ls='-.')
DISCS = disc_layout()
for d in DISCS: ax.add_patch(Circle((d['x'], d['y']), d['diam'] / 2, fc='#5b3a78', ec='none'))
ax.add_patch(Rectangle((x_open, -p['brk_w'] / 2), -x_open, p['brk_w'], fc='none', ec='#6b7785', lw=0.8, ls=(0, (3, 2))))
for (xb, yb) in BRK_SCREWS: ax.add_patch(Rectangle((xb - p['brk_adjust'] - p['screw_clear'] / 2, yb - p['screw_clear'] / 2), 2 * p['brk_adjust'] + p['screw_clear'], p['screw_clear'], fc='none', ec='#6b7785', lw=0.6))
rl = ax.plot([0, 0], [-YE - 6, YE + 12], color=RED, lw=1.4, ls=(0, (6, 3)))
owner(*rl, ax.text(2, YE + 16, 'x = 0: chest-wall face of the breast support (lip contact)', color=RED, fontsize=7.8, weight='bold', va='bottom'))
owner(*ax.plot([x_lid_out - 13, x_lid_out - 3], [0, 0], color=INK, lw=1.0), *ax.plot([R_out + 3, R_out + 13], [0, 0], color=INK, lw=1.0),
      ax.text(R_out + 15, 0, 'A', fontsize=10, weight='bold', va='center'))
owner(*ax.plot([-10, -10], [YE + 2, YE + 12], color=INK, lw=1.0), ax.text(-10, YE + 14, 'C', fontsize=10, weight='bold', ha='center', va='bottom'))
dim(ax, (7, -R_out), (7, R_out), '%.1f (imaged part, at mid-height)' % (2 * R_out), off=(R_out + 25, 0), gap=0)
dim(ax, (x_lid_out, -YE), (x_lid_out, YE), '%.1f (extension and lid)' % (2 * YE), off=(-20, 0), t=0.22)
dim(ax, (x_open, -YE), (0, -YE), '%.0f' % p['L_ext'], off=(0, -12))
radius(ax, 80, R_out, 'R %.2f outer' % R_out)
radius(ax, -80, p['R_flat'], 'R %.0f flat region' % p['R_flat'])
leader(ax, (-7, p['brk_w'] / 2), 'positioning bracket under the\nextension (hidden; 4 × M3\nin slots, ±%.0f mm)' % p['brk_adjust'], (-66, 112), ha='right')
leader(ax, (x_open + p['d_frame'] / 2, 66), 'bolting zone, %.0f mm deep' % p['d_frame'], (-66, 74), ha='right')
leader(ax, (x_lid_out, 40), 'plain lid %.0f mm + gasket %.0f mm' % (p['t_lid'], p['t_gasket']), (-66, 40), ha='right')
leader(ax, (-10, -P_tab), 'plate tab %.1f wide, runs in slots\nin the thickened end walls' % (2 * P_tab), (-66, -124), ha='right')
leader(ax, (x_cap_end, y_cap_out - 1), 'lugs close the slot ends (×2)', (36, YE + 6))
leader(ax, (R_plate * np.cos(np.radians(35)), R_plate * np.sin(np.radians(35))), 'plate edge R %.2f\n(sits in the rim)' % R_plate, (R_out + 45, 96))
d_br = max((d for d in DISCS if d['kind'] == 'array'), key=lambda d: d['x'] - d['y'])
leader(ax, (d_br['x'] + d_br['diam'] / 2, d_br['y']), 'iodine discs: measured July\nlayout, same distance from\nthe chest-wall edge', (R_out + 45, -84))
ax.text(-165, YE + 44, 'PLAN VIEW (top), rev B', fontsize=11, weight='bold')
ax.text(-165, YE + 31, 'dashed: flat region; dotted: cavity at mid-height; dash-dot: plate\ncross-hatched: plate slots; grey dashed: bracket; grey band: outside the field',
        fontsize=7.2, color=THIN, va='bottom')
ax.set_xlim(-168, R_out + 112); ax.set_ylim(-YE - 30, YE + 56)

# ======================= CHEST-WALL VIEW (lid)
ax = fig.add_axes([0.50, 0.71, 0.49, 0.265]); ax.set_aspect('equal'); ax.set_axis_off()
ax.add_patch(Polygon(stadium_yz(a_ext, b_out), fc=HAT, ec=INK, lw=1.1))
ax.add_patch(Polygon(stadium_yz(a_gin, b_open, zc=c_open), fc='none', ec=INK, lw=0.6, ls=(0, (4, 2))))
for (y, z) in SCREWS: ax.add_patch(Circle((y, z), p['screw_clear'] / 2, fc='white', ec=INK, lw=0.6))
for (y, z) in PORTS:
    ax.add_patch(Circle((y, z), p['port_d'] / 2, fc='white', ec=INK, lw=0.8)); ax.plot([y - 5, y + 5], [z, z], color=THIN, lw=0.3); ax.plot([y, y], [z - 5, z + 5], color=THIN, lw=0.3)
hs = p['t_plate'] / 2 + p['slot_clear']
for s_ in (-1, 1): ax.add_patch(Rectangle((min(s_ * R_cav, s_ * y_slot_out), ci - hs), y_slot_out - R_cav, 2 * hs, fc='none', ec=INK, lw=0.5, ls=':'))
ax.plot([-P_tab, P_tab], [ci, ci], color=CL, lw=0.8, ls='-.')
br = Rectangle((-p['brk_w'] / 2, -p['lip_h']), p['brk_w'], p['lip_h'], fc='#b8c2cc', ec=INK, lw=0.7); ax.add_patch(br)
owner(br, ax.text(0, -p['lip_h'] / 2, 'positioning bracket (behind it: the lip at x = 0)', fontsize=FS, ha='center', va='center'))
dim(ax, (-YE, -p['lip_h']), (YE, -p['lip_h']), '%.1f' % (2 * YE), off=(0, -9))
dim(ax, (YE + 2, 0), (YE + 2, p['H']), '%.1f' % p['H'], off=(10, 0), gap=0)
dim(ax, (-90, p['H'] + 2), (-60, p['H'] + 2), '30 (×6)', off=(0, 5), gap=0)
dim(ax, (p['brk_w'] / 2, -p['lip_h']), (p['brk_w'] / 2, 0), '%.0f' % p['lip_h'], off=(8, 0))
leader(ax, (-50, 36.75 + p['port_d'] / 2), '4 × M8 fill/vent ports (2 above, 2 below the plate)', (-42, p['H'] + 14))
leader(ax, (90, 50 + p['screw_clear'] / 2), '18 × M3, nylon or PEEK', (98, p['H'] + 14))
leader(ax, (p['R_flat'] + a_ext * np.cos(np.radians(30)), c + b_out * np.sin(np.radians(30))), 'ends: semi-ellipse\n%.2f × %.2f, thickened\nin the extension' % (a_ext, b_out), (YE + 22, 46))
leader(ax, (120, ci), 'plate tab (dash-dot), pressed\nby the plain lid', (YE + 22, 9))
ax.text(-YE - 8, p['H'] + 33, 'CHEST-WALL VIEW (plain lid, seen from the patient side)', fontsize=11, weight='bold')
ax.text(-YE - 8, p['H'] + 25, 'dashed: gasket inner edge; dotted: slot mouths; dash-dot: plate tab', fontsize=7.2, color=THIN)
ax.set_xlim(-YE - 10, YE + 92); ax.set_ylim(-p['lip_h'] - 24, p['H'] + 40)

# ======================= SECTION A-A (centre line)
ax = fig.add_axes([0.50, 0.40, 0.49, 0.29]); ax.set_aspect('equal'); ax.set_axis_off()
o = np.r_[[[x_open, p['H']]], ellipse_pts(p['R_flat'], c, a_out, b_out, 90, -90), [[x_open, 0]]]
i_ = np.r_[[[x_open, ci + b_in]], ellipse_pts(p['R_flat'], ci, a_in, b_in, 90, -90), [[x_open, ci - b_in]]]
ax.add_patch(Polygon(o, fc=HAT, ec=INK, lw=1.1)); ax.add_patch(Polygon(i_, fc='white', ec=INK, lw=0.8))
for z0, z1 in ((p['t_bot'], z_bot_open), (z_top_open, p['H'] - p['t_top'])): ax.add_patch(Rectangle((x_open, z0), p['d_frame'], z1 - z0, fc=HAT, ec=INK, lw=0.6))
ax.add_patch(Rectangle((x_lid_out, 0), p['t_lid'], p['H'], fc=HAT, ec=INK, lw=0.9))
for z0, z1 in ((0, z_bot_open), (z_top_open, p['H'])): ax.add_patch(Rectangle((x_lid_in, z0), p['t_gasket'], z1 - z0, fc='#333333', ec='none'))
ax.add_patch(Rectangle((x_tab_end, ci - p['t_plate'] / 2), R_plate - x_tab_end, p['t_plate'], fc='#e9e6dc', ec=CL, lw=0.8))
ax.plot([x_lid_out - 4, R_out + 6], [c, c], color=THIN, lw=0.4, ls=(0, (10, 3, 2, 3)))
rl = ax.plot([0, 0], [-p['lip_h'] - 6, p['H'] + 7], color=RED, lw=1.3, ls=(0, (6, 3)))
owner(*rl, ax.text(1, p['H'] + 9, 'x = 0 (support face)', color=RED, fontsize=7.5, weight='bold', va='bottom'))
ax.add_patch(Polygon([[x_open, 0], [0, 0], [0, -p['lip_h']], [-p['lip_t'], -p['lip_h']], [-p['lip_t'], -p['brk_t']], [x_open, -p['brk_t']]], fc='#b8c2cc', ec=INK, lw=0.8))
sl = ax.plot([R_out + 4, 0.4, 0.4], [-0.4, -0.4, -p['lip_h'] - 6], color='#6b6b6b', lw=0.9, ls=(0, (4, 2)))
owner(*sl, ax.text(64, -7.5, 'breast support', fontsize=FS, color='#6b6b6b', ha='center', va='center'))
dim(ax, (R_out + 2, 0), (R_out + 2, p['H']), '%.1f' % p['H'], off=(12, 0), gap=0)
dim(ax, (60, p['H'] - p['t_top']), (60, p['H']), '%.0f' % p['t_top'], ext=False)
dim(ax, (60, 0), (60, p['t_bot']), '%.0f' % p['t_bot'], ext=False)
dim(ax, (0, -2), (p['R_flat'], -2), '%.0f flat' % p['R_flat'], off=(0, -15), gap=0)
dim(ax, (x_open, -2), (0, -2), '%.0f' % p['L_ext'], off=(0, -23), gap=0)
leader(ax, (25, ci + p['t_plate'] / 2), 'plate %.0f mm, centred in the cavity (z = %.2f)' % (p['t_plate'], ci), (8, ci + 12))
leader(ax, (p['R_flat'] + (a_in + a_out) / 2, ci), 'wall %.1f at the equator\n(= 1.45 mm skin)' % p['t_edge'], (R_out + 24, 8))
leader(ax, (p['R_flat'] + a_out * np.cos(np.radians(55)), c + b_out * np.sin(np.radians(55))), 'outer: R %.2f (constant radius)' % a_out, (96, p['H'] + 13))
leader(ax, (p['R_flat'] + a_in * np.cos(np.radians(-60)), ci + b_in * np.sin(np.radians(-60))),
       'inner: semi-ellipse %.2f × %.2f, centre z = %.2f\n(wall %.0f / %.0f → %.1f)' % (a_in, b_in, ci, p['t_top'], p['t_bot'], p['t_edge']), (30, -39))
leader(ax, (-p['lip_t'], -10), 'positioning lip, %.0f deep' % p['lip_h'], (x_lid_out - 8, -10), ha='right')
leader(ax, (x_open + p['d_frame'] / 2, (z_top_open + p['H'] - p['t_top']) / 2), 'bolting zone: top +%.0f / bottom +%.0f mm\nover %.0f mm' % (p['ft_top'], p['ft_bot'], p['d_frame']),
       (x_lid_out - 8, p['H'] + 8), ha='right')
leader(ax, (x_lid_out, 28), 'plain lid %.0f + gasket %.0f' % (p['t_lid'], p['t_gasket']), (x_lid_out - 8, 28), ha='right')
ax.text(x_lid_out - 95, p['H'] + 32, 'SECTION A–A (centre line, chest wall left)', fontsize=11, weight='bold')
ax.set_xlim(x_lid_out - 98, R_out + 84); ax.set_ylim(-p['lip_h'] - 34, p['H'] + 40)

# ======================= SECTION C-C (x = -10, in the extension): slot and plate tab
ax = fig.add_axes([0.02, 0.04, 0.29, 0.33]); ax.set_aspect('equal'); ax.set_axis_off()
y0 = 92.0
o = np.r_[[[y0, p['H']]], ellipse_pts(p['R_flat'], c, a_ext, b_out, 90, -90), [[y0, 0]]]
i_ = np.r_[[[y0, ci + b_in]], ellipse_pts(p['R_flat'], ci, a_in, b_in, 90, -90), [[y0, ci - b_in]]]
ax.add_patch(Polygon(o, fc=HAT, ec=INK, lw=1.1)); ax.add_patch(Polygon(i_, fc='white', ec=INK, lw=0.8))
ax.add_patch(Rectangle((R_cav - 0.3, ci - hs), y_slot_out - R_cav + 0.3, 2 * hs, fc='white', ec=INK, lw=0.6))
ax.add_patch(Rectangle((y0, ci - p['t_plate'] / 2), P_tab - y0, p['t_plate'], fc='#e9e6dc', ec=CL, lw=0.8))
leader(ax, (y_slot_out, ci + hs), 'slot %.1f high, %.1f deep into the end wall' % (2 * hs, y_slot_out - R_cav), (y0 + 2, p['H'] + 9))
leader(ax, (P_tab, ci - p['t_plate'] / 2), 'plate tab edge, y = ±%.2f' % P_tab, (y0 + 2, -9))
leader(ax, (p['R_flat'] + a_ext * np.cos(np.radians(-20)), c + b_out * np.sin(np.radians(-20))), 'end wall %.1f mm\nat the equator' % (a_ext - a_in), (YE + 10, 6))
ax.text(y0, p['H'] + 18, 'SECTION C–C, right end\n(x = −10 mm, in the extension)', fontsize=11, weight='bold', va='bottom')
ax.set_xlim(y0 - 2, YE + 42); ax.set_ylim(-20, p['H'] + 40)

# ======================= notes and parameters
ax = fig.add_axes([0.33, 0.03, 0.66, 0.33]); ax.set_axis_off()
rows = [('Outer height H', '%.1f mm (outer, confirmed)' % p['H']),
        ('Imaged part', 'half-disc: outer R %.2f, flat region R %.0f; top %.0f / bottom %.0f mm walls (the current L2 has 6 / 4); rim radius %.2f, wall → %.1f mm' % (R_out, p['R_flat'], p['t_top'], p['t_bot'], a_out, p['t_edge'])),
        ('Extension', '%.0f mm beyond the support face, %.1f mm wide; end walls %.1f mm at the equator; filled like the rest' % (p['L_ext'], 2 * YE, a_ext - a_in)),
        ('Bolting zone', 'last %.0f mm of the extension, top +%.0f / bottom +%.0f mm; 18 × M3 blind tapped %.0f mm' % (p['d_frame'], p['ft_top'], p['ft_bot'], p['screw_depth'])),
        ('Lid', 'plain %.0f mm cast PMMA (no rib); 4 × M8 ports; 1 mm FKM gasket' % p['t_lid']),
        ('Plate (new)', 'half-disc R %.2f + tab %.1f × %.1f mm, %.0f mm thick, centred in the cavity z = %.2f' % (R_plate, -x_tab_end, 2 * P_tab, p['t_plate'], ci)),
        ('Plate support', 'tab in two slots (%.1f × %.1f mm) in the extension, slot ends closed by lugs; curved edge in the rim; tab end pressed by the lid' % (2 * hs, y_slot_out - R_cav)),
        ('Positioning', 'removable L-bracket (PMMA): lip %.0f × %.0f × %.0f mm, inner face at x = 0 against the support; 4 × M3 in slots (±%.0f mm)' % (p['brk_w'], p['lip_t'], p['lip_h'], p['brk_adjust'])),
        ('Overall', 'depth %.2f mm (lid to rim); overhang beyond the support face %.0f mm' % (R_out - x_lid_out, -x_lid_out)),
        ('Filling volume', '≈ %.0f ml, of which ≈ %.0f ml in the extension' % (VOL, VEXT))]
y = 0.98
ax.text(0.0, y, 'MAIN PARAMETERS (all mm; parametric in l2v2_geom.py; current STL files in stl/)', fontsize=10, weight='bold', va='top', transform=ax.transAxes); y -= 0.07
for k, v in rows:
    ax.text(0.0, y, k, fontsize=8.5, va='top', transform=ax.transAxes, weight='bold'); ax.text(0.15, y, v, fontsize=8.5, va='top', transform=ax.transAxes); y -= 0.058
ax.text(0.0, y - 0.02, 'Not to be machined from this sheet without checking the open points in DESIGN_NOTES.md (machine clearance of the overhang, oil expansion, material, printer size).',
        fontsize=8, color='#a33', va='top', transform=ax.transAxes)
fig.text(0.985, 0.012, 'L2 phantom container v2, rev B with positioning bracket (equal walls, 2 mm plate, slot-end lugs) — drawing of 20 Nov 2026 — scale as shown, dimensions in mm',
         ha='right', fontsize=8, color=THIN)
problems = collisions(fig)
print('drawing sheet: %d collision(s)' % len(problems)); [print('  ' + m) for m in problems]
fig.savefig(os.path.join(OUT, 'L2v2_drawing.pdf')); fig.savefig(os.path.join(OUT, 'L2v2_drawing.png'), dpi=130); plt.close(fig)

# ======================= edge profile: current vs v2 vs compressed breast
fig, axs = plt.subplots(2, 3, figsize=(15, 7.2), gridspec_kw=dict(height_ratios=[1, 1.05]))
H = p['H']; d = np.linspace(-45, 1, 2000)          # distance from the outermost edge (0 = outer edge)
def cur(dd, tw=8.0, tb=10.0):  # current L2 (measured): vertical wall 8 mm, lid 6 + base 4 mm
    fill = np.where(dd < -tw, H - tb, 0.0); pm = np.where(dd < -tw, tb, np.where(dd <= 0, H, 0.0)); return fill, pm
def v2(dd):
    q = dd + a_out                                   # radial position relative to the flat-region edge
    out_ = np.where(q <= 0, H, 2 * b_out * np.sqrt(np.clip(1 - (q / a_out) ** 2, 0, 1)))
    inn = np.where(q <= 0, H - p['t_top'] - p['t_bot'], 2 * b_in * np.sqrt(np.clip(1 - (np.maximum(q, 0) / a_in) ** 2, 0, 1)) * (q < a_in))
    return inn, out_ - inn
def breast(dd, skin=1.5):
    q = dd + H / 2
    out_ = np.where(q <= 0, H, H * np.sqrt(np.clip(1 - (q / (H / 2)) ** 2, 0, 1)))
    ri = H / 2 - skin; inn = np.where(q <= 0, H - 2 * skin, 2 * ri * np.sqrt(np.clip(1 - (np.maximum(q, 0) / ri) ** 2, 0, 1)) * (q < ri))
    return inn, out_ - inn
cases = [('Current L2 (measured: 8 mm wall, 6 + 4 mm lid/base)', cur, 'PMMA wall'), ('L2 v2 (constant-radius rim, wall → 1.5 mm)', v2, 'PMMA wall'),
         ('Compressed breast (idealised, 1.5 mm skin)', breast, 'skin')]
for j, (ttl, fn, wname) in enumerate(cases):
    ax = axs[0, j]; ax.set_aspect('equal')
    # draw profile cross-section over the outer 45 mm
    if fn is cur:
        ax.add_patch(Rectangle((-45, 0), 45, H, fc=HAT, ec=INK, lw=0.8)); ax.add_patch(Rectangle((-45, 4), 37, H - 10, fc='#fbf3dc', ec=INK, lw=0.6))
    else:
        tt = np.linspace(-np.pi / 2, np.pi / 2, 200)
        a_o, b_o = (a_out, b_out) if fn is v2 else (H / 2, H / 2)
        a_i, b_i = (a_in, b_in) if fn is v2 else (H / 2 - 1.5, H / 2 - 1.5)
        x0 = -a_o
        ax.add_patch(Polygon(np.r_[[[-45, H]], np.c_[x0 + a_o * np.cos(tt[::-1]), c + b_o * np.sin(tt[::-1])], [[-45, 0]]], fc=HAT if fn is v2 else '#f1c9b8', ec=INK, lw=0.8))
        zi = ci if fn is v2 else c
        ax.add_patch(Polygon(np.r_[[[-45, zi + b_i]], np.c_[x0 + a_i * np.cos(tt[::-1]), zi + b_i * np.sin(tt[::-1])], [[-45, zi - b_i]]], fc='#fbf3dc', ec=INK, lw=0.6))
    ax.set_xlim(-45, 3); ax.set_ylim(-3, H + 3); ax.set_title(ttl, fontsize=9.5, loc='left'); ax.set_yticks([0, H]); ax.set_xticks([-40, -30, -20, -10, 0])
    ax.tick_params(labelsize=7.5)
    ax = axs[1, j]; fill, wall_ = fn(d)
    ax.plot(d, fill, color='#c68a12', lw=1.6, label='filling (spheres + oil) / breast tissue')
    ax.plot(d, wall_, color=CL, lw=1.6, label=wname + ' path')
    ax.set_xlim(-45, 3); ax.set_ylim(0, 60); ax.grid(alpha=0.25); ax.set_xlabel('distance from the outer edge (mm)', fontsize=8.5)
    if j == 0: ax.set_ylabel('vertical path length (mm)', fontsize=8.5)
    ax.legend(fontsize=7.5, loc='upper left', frameon=False); ax.tick_params(labelsize=7.5)
fig.suptitle('Edge profile and vertical path lengths: the curved rim gives a breast-like thickness fall-off and a skin-equivalent wall at the edge',
             fontsize=11, x=0.01, ha='left')
fig.text(0.01, 0.005, 'Path lengths along vertical rays, ignoring beam divergence. Current L2 as measured (wall 8 mm, lid 6 mm, base 4 mm). Breast: idealised compressed-breast edge (semicircle) with uniform 1.5 mm skin.',
         fontsize=8, color=THIN)
plt.tight_layout(rect=(0, 0.02, 1, 0.95)); fig.savefig(os.path.join(OUT, 'L2v2_edge_profile.png'), dpi=130); plt.close(fig)
print('done')

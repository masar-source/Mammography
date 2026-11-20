#!/usr/bin/env python3
"""L2 phantom container v2 (rev B + positioning bracket): exact B-rep solids with CadQuery, exported as STEP and STL.
Reads params.json (written by l2v2_geom.py), so both models use the same parameters.
Usage (needs CadQuery: pip install cadquery):  python l2v2_step.py [--check]
  -> step/*.step, stl/*.stl, step/brep_checks.json; --check also compares every part with the implicit model."""
import json, os, sys, cadquery as cq

HERE = os.path.dirname(os.path.abspath(__file__))
P = json.load(open(os.path.join(HERE, 'params.json'))); D = P['derived']
H, R_flat = P['H'], P['R_flat']
c = H / 2; b_out = H / 2; a_out = P['a_out']; a_in, b_in, a_ext, ci = D['a_in'], D['b_in'], D['a_ext'], D['ci']
R_plate, x_open, x_lid_in, x_lid_out = D['R_plate'], D['x_open'], D['x_lid_in'], D['x_lid_out']
y_slot_in, y_slot_out, P_tab, b_open, a_gin = D['y_slot_in'], D['y_slot_out'], D['P_tab'], D['b_open'], D['a_gin']
z_top_open, z_bot_open, c_open = D['z_top_open'], D['z_bot_open'], D['c_open']
x_tab_end, x_tab_stop = D['x_tab_end'], D['x_tab_stop']
SCREWS, PORTS = [tuple(s) for s in P['screws']], [tuple(s) for s in P['ports']]
BRK_SCREWS = [tuple(s) for s in P['bracket_screws']]
R_ext = R_flat + a_ext

def half_rev(a, b, zc=None):
    """solid of revolution: flat to R_flat, elliptical rim (semi-axes a radial, b vertical, centre height zc); x >= 0 half"""
    zc = c if zc is None else zc
    prof = (cq.Workplane('XZ').moveTo(0, zc - b).lineTo(R_flat, zc - b)
            .ellipseArc(a, b, -90, 90, startAtCurrent=True).lineTo(0, zc + b).close())
    s = prof.revolve(180, (0, 0, 0), (0, 1, 0)).rotate((0, 0, 0), (0, 0, 1), -90)
    if s.val().BoundingBox().xmin < -1:          # revolve sense: make sure it is the x >= 0 half
        s = s.rotate((0, 0, 0), (0, 0, 1), 180)
    return s
def stadium_prism(a, b, x0, x1, zc=None):
    """(y,z) stadium: flat for |y| <= R_flat, elliptical ends, centre height zc; extruded along x from x0 to x1"""
    zc = c if zc is None else zc
    prof = (cq.Workplane('YZ', origin=(x0, 0, 0)).moveTo(-R_flat, zc - b).lineTo(R_flat, zc - b)
            .ellipseArc(a, b, -90, 90, startAtCurrent=True).lineTo(-R_flat, zc + b)
            .ellipseArc(a, b, 90, 270, startAtCurrent=True).close())
    return prof.extrude(x1 - x0)
def box(x0, x1, y0, y1, z0, z1):
    return cq.Workplane('XY').box(x1 - x0, y1 - y0, z1 - z0, centered=False).translate((x0, y0, z0))
def holes(pts, d, x0, x1):
    out = None
    for (y, z) in pts:
        cyl = cq.Workplane('YZ', origin=(x0, 0, 0)).center(y, z).circle(d / 2).extrude(x1 - x0)
        out = cyl if out is None else out.union(cyl)
    return out

def container():
    outer = half_rev(a_out, b_out).union(stadium_prism(a_ext, b_out, x_open, 0))
    cav = half_rev(a_in, b_in, ci).union(stadium_prism(a_in, b_in, x_open - 1, 0.01, ci))
    body = outer.cut(cav)
    frame = stadium_prism(a_in, b_in, x_open, x_open + P['d_frame'], ci).intersect(box(x_open - 1, x_open + P['d_frame'] + 1, -R_ext, R_ext, 0, z_bot_open).union(
        box(x_open - 1, x_open + P['d_frame'] + 1, -R_ext, R_ext, z_top_open, H)))
    body = body.union(frame)
    # outer lugs closing the slot ends; they start 0.5 mm inside the extension wall (no coplanar faces with the step face
    # at x = 0) and are added before the slots are cut, so the slots still run to x = 0
    for ya, yb in ((D['y_cap_in'], D['y_cap_out']), (-D['y_cap_out'], -D['y_cap_in'])):
        body = body.union(box(-0.5, D['x_cap_end'], ya, yb, ci - D['z_cap_half'], ci + D['z_cap_half']))
    hs = P['t_plate'] / 2 + P['slot_clear']
    for ya, yb in ((y_slot_in, y_slot_out), (-y_slot_out, -y_slot_in)):
        body = body.cut(box(x_open - 1, 0, ya, yb, ci - hs, ci + hs))
    body = body.cut(holes(SCREWS, P['screw_tap'], x_open - 1, x_open + P['screw_depth']))
    for (xb, yb) in BRK_SCREWS:      # blind M3 taps from below for the positioning bracket
        body = body.cut(cq.Workplane('XY', origin=(0, 0, -1)).center(xb, yb).circle(P['screw_tap'] / 2).extrude(1 + P['brk_depth']))
    return body
def bracket():
    """removable L-shaped positioning bracket; the lip's inner face (x = 0) hooks against the support's chest-wall face"""
    w = P['brk_w']
    b = box(x_open, 0, -w / 2, w / 2, -P['brk_t'], 0).union(box(-P['lip_t'], 0, -w / 2, w / 2, -P['lip_h'], 0))
    for (xb, yb) in BRK_SCREWS:      # slotted clearance holes, +/- brk_adjust in x
        b = b.cut(cq.Workplane('XY', origin=(0, 0, -P['brk_t'] - 1)).center(xb, yb)
                  .slot2D(2 * P['brk_adjust'] + P['screw_clear'], P['screw_clear']).extrude(P['brk_t'] + 2))
    return b
def lid():
    return stadium_prism(a_ext, b_out, x_lid_out, x_lid_in).cut(holes(SCREWS, P['screw_clear'], x_lid_out - 1, x_lid_in + 1)).cut(
        holes(PORTS, P['port_d'], x_lid_out - 1, x_lid_in + 1))
def gasket():
    return stadium_prism(a_ext, b_out, x_lid_in, x_open).cut(stadium_prism(a_gin, b_open, x_lid_in - 1, x_open + 1, c_open)).cut(
        holes(SCREWS, P['screw_clear'], x_lid_in - 1, x_open + 1))
def plate():
    z0, z1 = ci - P['t_plate'] / 2, ci + P['t_plate'] / 2
    disc = box(x_tab_stop, R_plate, -R_plate, R_plate, z0, z1).intersect(cq.Workplane('XY').circle(R_plate).extrude(H))
    return disc.union(box(x_tab_end, x_tab_stop, -P_tab, P_tab, z0, z1))

def write_stl(wp, path, settings=((0.02, 0.1), (0.025, 0.1), (0.015, 0.1), (0.01, 0.1), (0.03, 0.15), (0.05, 0.2))):
    """binary STL with absolute chord deflection (mm) and angular deflection (rad); OCC's mesher can skip a face
    (it skipped the -y step face of the container at 0.02 mm, leaving a hole in the STL), so the first setting that
    triangulates every face is used"""
    from OCP.BRepMesh import BRepMesh_IncrementalMesh
    from OCP.BRep import BRep_Tool
    from OCP.BRepTools import BRepTools
    from OCP.TopLoc import TopLoc_Location
    from OCP.StlAPI import StlAPI_Writer
    s = wp.val()
    for defl, ang in settings:
        BRepTools.Clean_s(s.wrapped); BRepMesh_IncrementalMesh(s.wrapped, defl, False, ang, True)
        if all(BRep_Tool.Triangulation_s(f.wrapped, TopLoc_Location()) is not None for f in s.Faces()):
            w = StlAPI_Writer(); w.ASCIIMode = False; w.Write(s.wrapped, path)
            return dict(deflection_mm=defl, angular_rad=ang)
    raise RuntimeError('no mesh setting triangulates every face of %s' % path)
def tight_bbox(shape):
    """bounding box (xmin, xmax, ymin, ymax, zmin, zmax) from a fine tessellation; the OCC box can be loose on the
    trimmed elliptical/revolved faces of the boolean results"""
    import numpy as np
    V = np.array([(v.x, v.y, v.z) for v in shape.tessellate(0.05, 0.1)[0]])
    return tuple(float(v) for v in (V[:, 0].min(), V[:, 0].max(), V[:, 1].min(), V[:, 1].max(), V[:, 2].min(), V[:, 2].max()))
def crosscheck(parts, n=6000, seed=0):
    """inside/outside agreement with the implicit model of l2v2_geom.py (uniform points plus points near the surfaces)"""
    import numpy as np, l2v2_geom as G
    rng = np.random.default_rng(seed); fun = {'container': G.f_container, 'lid': G.f_lid, 'gasket': G.f_gasket, 'plate': G.f_plate, 'bracket': G.f_bracket}; res = {}
    for name, wp in parts.items():
        sol = wp.val(); f = fun[name]; bb = G.PARTS[name][1]
        pts = np.c_[rng.uniform(bb[0], bb[1], n), rng.uniform(bb[2], bb[3], n), rng.uniform(bb[4], bb[5], n)]
        fv = f(pts[:, 0], pts[:, 1], pts[:, 2]); pts = np.r_[pts[:n // 2], pts[np.abs(fv) < 1.5][:n // 2]]
        fv = f(pts[:, 0], pts[:, 1], pts[:, 2]); tested = mism = 0
        for q, v in zip(pts, fv):
            if abs(v) < 0.02: continue
            tested += 1; mism += int(sol.isInside(cq.Vector(*q), tolerance=1e-4) != (v < 0))
        res[name] = dict(points=tested, mismatches=mism); print('  crosscheck', name, res[name], flush=True)
    return res

if __name__ == '__main__':
    os.makedirs(os.path.join(HERE, 'step'), exist_ok=True); os.makedirs(os.path.join(HERE, 'stl'), exist_ok=True)
    parts = {'container': container(), 'lid': lid(), 'gasket': gasket(), 'plate': plate(), 'bracket': bracket()}
    checks = {}
    for name, wp in parts.items():
        s = wp.val(); bb = tight_bbox(s)
        checks[name] = dict(valid=bool(s.isValid()), n_solids=len(wp.solids().vals()), volume_ml=s.Volume() / 1000,
                            bbox=[round(v, 3) for v in bb])
        cq.exporters.export(wp, os.path.join(HERE, 'step', 'L2v2_%s.step' % name))
        checks[name]['stl_mesh'] = write_stl(wp, os.path.join(HERE, 'stl', 'L2v2_%s.stl' % name))
        print(name, checks[name], flush=True)
    cav = half_rev(a_in, b_in, ci).union(stadium_prism(a_in, b_in, x_open, 0.01, ci))
    inside = parts['container'].intersect(cav).val().Volume() + parts['plate'].intersect(cav).val().Volume()
    checks['filling_volume_ml'] = (cav.val().Volume() - inside) / 1000
    cav_ext = stadium_prism(a_in, b_in, x_open, 0, ci)
    checks['filling_in_extension_ml'] = (cav_ext.val().Volume() - parts['container'].intersect(cav_ext).val().Volume()
                                         - parts['plate'].intersect(cav_ext).val().Volume()) / 1000
    print('filling volume %.1f ml (of which %.1f ml in the extension)' % (checks['filling_volume_ml'], checks['filling_in_extension_ml']))
    assy = cq.Assembly(name='L2_phantom_v2_revB')
    cols = {'container': cq.Color(0.70, 0.84, 0.95, 0.35), 'lid': cq.Color(0.62, 0.78, 0.94, 0.45), 'gasket': cq.Color(0.15, 0.15, 0.15), 'plate': cq.Color(0.93, 0.93, 0.90), 'bracket': cq.Color(0.55, 0.62, 0.70)}
    for name, wp in parts.items(): assy.add(wp, name=name, color=cols[name])
    assy.export(os.path.join(HERE, 'step', 'L2v2_assembly.step'))
    if '--check' in sys.argv: checks['crosscheck_implicit_vs_step'] = crosscheck(parts)
    json.dump(checks, open(os.path.join(HERE, 'step', 'brep_checks.json'), 'w'), indent=1)
    print('done')

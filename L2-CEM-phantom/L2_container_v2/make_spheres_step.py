#!/usr/bin/env python3
"""Exact STEP spheres (CadQuery) for the six background diameters, plus the order sheet for the L2 v2 container.
Usage (needs CadQuery): python make_spheres_step.py -> spheres/L2v2_sphere_D*.step, spheres/sphere_order_sheet.csv"""
import os, json, csv, math, cadquery as cq
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'spheres')
D_LIST = [15.88, 12.70, 9.52, 6.35, 3.18, 1.58]
V_fill = json.load(open(os.path.join(HERE, 'step', 'brep_checks.json')))['filling_volume_ml']   # v2 container
V_size = 0.61 * V_fill / len(D_LIST)                                                           # ml per size (equal volumes)
rows = []
for D in D_LIST:
    s = cq.Workplane('XY').sphere(D / 2)
    cq.exporters.export(s, os.path.join(OUT, 'L2v2_sphere_D%.2fmm.step' % D))
    v1 = math.pi / 6 * D ** 3 / 1000                       # ml per sphere
    n = V_size / v1
    rows.append(dict(diameter_mm='%.2f' % D, stl='L2v2_sphere_D%.2fmm.stl' % D, step='L2v2_sphere_D%.2fmm.step' % D,
                     volume_per_size_ml='%.0f' % V_size, mass_per_size_g_at_1p046='%.0f' % (V_size * 1.046),
                     quantity_nominal=int(math.ceil(n)), quantity_with_10pct_spare=int(math.ceil(1.1 * n)),
                     step_volume_check_pct='%+.4f' % (100 * (s.val().Volume() / 1000 / v1 - 1))))
with open(os.path.join(OUT, 'sphere_order_sheet.csv'), 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
for r in rows: print(r)

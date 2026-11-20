#!/usr/bin/env python3
"""Underside view of the removable positioning bracket (dropped 30 mm, with its four M3 screws). Usage: python l2v2_render_bracket.py"""
import os, numpy as np
from l2v2_render import Cam, render, down, annotate, mesh, cylinder, merge, boxmesh, OUT, ACRYL, BRK, NYL, SUP
from l2v2_geom import p, c, R_out, x_open, x_lid_out, a_ext, BRK_SCREWS, f_container, f_bracket, PARTS

SS = 2; W, H = 1500 * SS, 1000 * SS; DZ = -30.0
cont = mesh(f_container, PARTS['container'][1], 0.9)
brk = mesh(f_bracket, PARTS['bracket'][1], 0.4)
scr = merge([cylinder((xb, yb, DZ - p['brk_t'] - 2.0), (0, 0, 1), 1.5, 12.0 - DZ - 6, 20) for (xb, yb) in BRK_SCREWS]
            + [cylinder((xb, yb, DZ - p['brk_t'] - 2.0), (0, 0, 1), 2.75, 2.0, 24) for (xb, yb) in BRK_SCREWS])
YE = p['R_flat'] + a_ext
pts = np.array([[x, y, z] for x in (x_lid_out, 60) for y in (-YE, YE) for z in (DZ - p['lip_h'] - 4, p['H'])])
cam = Cam((-260, -230, -150), (-10, 0, 0), W, H, fov=34); cam.fit(pts, margin=(0.05, 0.05, 0.14, 0.14))
sc = [dict(V=brk[0] + np.float32([0, 0, DZ]), F=brk[1], N=brk[2], col=BRK), dict(V=scr[0], F=scr[1], N=scr[2], col=NYL),
      dict(V=cont[0], F=cont[1], N=cont[2], col=ACRYL, glass=True)]
img = down(render(sc, cam), SS)
annotate(img, cam, [((-p['lip_t'] / 2, -40, DZ - p['lip_h'] + 2), 'lip: inner face rests against the\nchest-wall face of the support (x = 0)', 'bottom'),
                    ((BRK_SCREWS[0][0], BRK_SCREWS[0][1], DZ - p['brk_t']), 'slotted holes: ±%.0f mm adjustment' % p['brk_adjust'], 'bottom'),
                    ((BRK_SCREWS[3][0], BRK_SCREWS[3][1], 2.0), '4 × M3 nylon/PEEK into blind taps\n(7 mm bottom wall, 2 mm left above)', 'top'),
                    ((40, 0, 0), 'container (underside)', 'top')],
         os.path.join(OUT, 'L2v2_render_bracket.png'), 'L2 phantom v2 (rev B) — removable positioning bracket, seen from below (dropped 30 mm)', SS)
print('done')

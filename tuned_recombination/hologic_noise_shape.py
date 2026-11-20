#!/usr/bin/env python3
"""Noise shape g(b) for the Hologic profile of l2_tuned_recombination_multisystem.py, measured once from the open field
of the Hologic CEM images (no published MTF of this detector was at hand).

For every LE/HE pair: ROIs of the method's size (10.9 mm) lying completely in the open field (raw value above 85 % of
the open-field level, at least 8 mm from the object and 5 mm from the image border, no saturated pixel); in each ROI the
log difference D = Lam_H - w Lam_L (w = 0.41, the cancelling weight of the phantom on this system) is detrended and
Hann-windowed as in the method; its power per band is divided by the power in the last band. The table is the median
over the exposures. The open field also holds large-scale structure that is not quantum noise (calibration residue), so
the lowest bands are an upper estimate. The effect on the result is small: flat g = 1 instead of this table changes the
left-in share by about 0.4 percentage points; a Siemens-like MTF shape by 0.03 points (two exposures tested).
Usage: python3 hologic_noise_shape.py <folder with the Hologic DICOM files>"""
import os, sys, glob, numpy as np, pydicom
from scipy import ndimage as ndi

PX = 0.07
EDGES = [0.06, 0.1, 0.15, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0, 4.0, round(0.5 / PX, 2)]
N = int(2 * round(128 * 0.085 / PX / 2)); W = 0.41


def shape(L, H):
    sat = 2 ** 14 - 1
    lvl = np.percentile(H[H < sat], 99.5); obj = ndi.binary_opening(H < 0.6 * lvl, iterations=3)
    m = (H > 0.85 * lvl) & (H < sat) & (L < sat) & ~ndi.binary_dilation(obj, iterations=int(8 / PX))
    b = int(5 / PX); m[:b] = m[-b:] = False; m[:, :b] = m[:, -b:] = False
    win = np.outer(np.hanning(N), np.hanning(N)); yv, xv = np.mgrid[0:N, 0:N]
    Ad = np.c_[np.ones(N * N), xv.ravel(), yv.ravel(), (xv * xv).ravel(), (yv * yv).ravel(), (xv * yv).ravel()]; pinv = np.linalg.pinv(Ad)
    f = np.fft.fftfreq(N, d=PX); FR = np.hypot(f[:, None], f[None, :])
    bands = [np.flatnonzero(((FR >= a) & (FR < c)).ravel()) for a, c in zip(EDGES[:-1], EDGES[1:])]
    D = -np.log(H) + W * np.log(L); P = np.zeros(len(bands)); n = 0
    for y0 in range(0, H.shape[0] - N, N // 2):
        for x0 in range(0, H.shape[1] - N, N // 2):
            if m[y0:y0 + N, x0:x0 + N].all():
                a = D[y0:y0 + N, x0:x0 + N]; F = np.abs(np.fft.fft2(win * (a - (Ad @ (pinv @ a.ravel())).reshape(N, N)))).ravel() ** 2
                P += [F[k].mean() for k in bands]; n += 1
    return P / P[-1], n


if __name__ == '__main__':
    files = sorted(f for f in glob.glob(os.path.join(sys.argv[1], '**', 'MG*'), recursive=True) if os.path.isfile(f) and not os.path.basename(f).startswith('.'))
    hd = {f: pydicom.dcmread(f, stop_before_pixels=True) for f in files}
    raw = lambda f, k: k in list(hd[f].ImageType) and 'FOR PROCESSING' in str(hd[f].get('PresentationIntentType', ''))
    G = []
    for f in [f for f in files if raw(f, 'LOW_ENERGY')]:
        h = [g for g in files if raw(g, 'HIGH_ENERGY') and hd[g].AcquisitionTime == hd[f].AcquisitionTime]
        if len(h) != 1:
            continue
        g, n = shape(pydicom.dcmread(f).pixel_array.astype(float), pydicom.dcmread(h[0]).pixel_array.astype(float)); G.append(g)
        print('%s: %d open-field ROIs  %s' % (os.path.basename(f), n, ' '.join('%.2f' % v for v in g)))
    G = np.array(G); print('median of %d exposures: [%s]' % (len(G), ', '.join('%.2f' % v for v in np.median(G, 0))))
    print('lowest / highest exposure per band:', ' '.join('%.2f-%.2f' % (a, b) for a, b in zip(G.min(0), G.max(0))))

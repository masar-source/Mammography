#!/usr/bin/env python3
"""
L2 phantom: tuned recombination of Siemens TiCEM (Insight CEM) images.

What it does, for every exposure found (raw LE, raw HE and Siemens' recombined image R):
  Siemens' recombined image is kept exactly as it is (walls, edges, mask, grey scale) and only the weight on the
  background is changed: the sphere background is cancelled except for a soft remainder at coarse scales.
      R_tuned = R - IFFT{ c(f) FFT(Lam_L) },   c(f) = alpha(f) [ w_t(f) - w_eff(f) ],   Lam = -ln(stored value)
  alpha(f), w_eff(f)  Siemens' gain and effective weight per spatial-frequency band, measured from R itself by
                      regressing R on Lam_H and Lam_L in interior regions (R = alpha (Lam_H - w_eff Lam_L) + rest);
  w*(f)               the weight that cancels the spheres, S_HL/S_LL, from the same regions;
  w_t(f)              = w*(f) [1 - rho(f)]: the tuned weight. rho is the share of the sphere contrast left in:
                      one value rho0 at coarse scales (to 0.3 cy/mm), faded to zero at 1 cy/mm.
  Strength 1.0        rho0 is set so that the residual structure at 0.06-0.3 cy/mm has the power that Siemens'
                      recombination leaves in real breasts (BREAST table below: aggregated numbers from 9 views of
                      2 patients; no patient data). Other strengths scale the left-in share: 0 = full cancellation,
                      0.5 = half the share, 2 = double (about the grey-level texture of Siemens' breast images,
                      which Siemens displays with 1.6-2x more gain than the phantom).
Everything is measured from each image on its own; nothing is averaged over images. Regions that contain iodine discs,
bubbles or other content that LE does not predict are found and left out automatically (no disc template needed).

Output, in the folder where the script is run, one folder per strength and main folder of the data:
  <main>_strength_<s>/<session>/<id>_tuned.dcm        the tuned image (copy of Siemens' DICOM with new pixel data)
                                 <id>_comparison.png   three images: Siemens' R, ours, and the difference (ours minus
                                                       Siemens), in Siemens' default window and a narrow one
                                 <id>_curves.png/.csv  weights, left-in share and residual spectra of this image
  <main>_strength_<s>/summary.csv                      one line per image
Check of the method (option --siemens-weighting): the raw LE and HE images are recombined with our method at Siemens'
own weighting, R_ours = offset + IFFT{ alpha(f) [FFT(Lam_H) - w_eff(f) FFT(Lam_L)] }, with the gain alpha(f) and weight
w_eff(f) per band as measured in interior regions, applied to the whole images. If the method is valid, this reproduces
Siemens' recombined image. Not part of the method, so not reproduced: Siemens' processing of the outer ~15 mm.
  <main>_siemens_weighting/<session>/<id>_siemens_weighting.dcm, <id>_comparison.png (Siemens, ours, difference),
                                     <id>_curves.png/.csv (weight in both images, share reproduced per band, spectra)
  <main>_siemens_weighting/summary.csv   correlation and rms difference in the interior, share reproduced, edge zone
All numbers of all analyses together: L2_tuned_recombination_numbers.xlsx in the output folder (sheets: tuned images,
check at Siemens' weighting, the curves per frequency band, notes). It is rebuilt at the end of every run from the CSV
files and so also contains earlier runs. Needs the Python package openpyxl; without it only the CSV files are written.
  <main> = the folder named like "28kV" on the path of the image; <session> = the session folder;
  <id> = session + series number of Siemens' recombined image.

Usage (needs Python 3 with numpy, scipy, pydicom, matplotlib; openpyxl for the Excel file):
  python3 l2_tuned_recombination.py                              all exposures below the current folder, strength 1.0
  python3 l2_tuned_recombination.py 28kV 30kV --strengths 1.0 0.5 2.0
  python3 l2_tuned_recombination.py --siemens-weighting          strength 1.0 and the check at Siemens' weighting
  python3 l2_tuned_recombination.py 28kV --only HIGHMAS --limit 2          a first look: two exposures per session
About 16-20 s per exposure, plus a few seconds for each further strength; about 22 MB of output per exposure and strength.
Folders named *4AFC* (copies of the same images) and the output folders of earlier runs are not searched; an image that
exists in several copies is processed once. Results that already exist are kept (a stopped run can be continued by
starting it again; --overwrite makes everything anew).
Research tool: the tuned images are not Siemens output and not for diagnostic use.
"""
import argparse, csv, os, re, sys, time, traceback
import numpy as np
import pydicom
from pydicom.uid import generate_uid
from scipy import ndimage as ndi, fft as sfft

# ----------------------------------------------------------------------------------------------- fixed numbers
EDGES = [0.06, 0.1, 0.15, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0, 4.0, 5.88]      # spatial-frequency bands (cy/mm)
FC = np.sqrt(np.array(EDGES[:-1]) * np.array(EDGES[1:]))
# Real breasts under Siemens' recombination, per band (mm2, log-signal units): structured residual (total minus quantum
# noise), median of 9 views; and the total residual (median, lowest, highest view) for the curves. 28 Oct / 17 Nov 2026.
BREAST = dict(
    T=[7.961e-06, 1.099e-05, 1.147e-05, 8.156e-06, 3.854e-06, 1.797e-06, 9.085e-07, 3.779e-07, 0, 0, 0, 0, 0],
    PD=[1.027e-05, 1.328e-05, 1.373e-05, 1.037e-05, 6.005e-06, 3.854e-06, 2.844e-06, 2.163e-06, 1.692e-06, 1.327e-06, 1.12e-06, 9.037e-07, 6.037e-07],
    PD_LO=[4.633e-06, 5.46e-06, 5.439e-06, 4.235e-06, 3.274e-06, 2.787e-06, 2.501e-06, 2.033e-06, 1.462e-06, 1.153e-06, 1.001e-06, 8.347e-07, 5.647e-07],
    PD_HI=[2.457e-05, 3.372e-05, 2.704e-05, 1.631e-05, 1.007e-05, 6.22e-06, 3.516e-06, 2.574e-06, 1.82e-06, 1.387e-06, 1.167e-06, 9.7e-07, 6.738e-07])
MTF = (np.arange(15.0), np.array([1.00, 0.91, 0.86, 0.78, 0.67, 0.55, 0.44, 0.34, 0.25, 0.18, 0.13, 0.09, 0.06, 0.03, 0.0]))   # detector MTF (NHSBSP 2019)
F0 = 0.05            # cy/mm: Siemens' large-scale flattening, gain ~ f^2/(f^2 + F0^2) below the first band
N = 128              # ROI size (px) for the spectra, half overlap
MARGIN_MM = 25.5     # ROIs lie at least this far inside the phantom outline
FADE = (0.3, 1.0)    # cy/mm: the left-in share is constant to FADE[0] and zero from FADE[1]
NARROW = 800.0       # width of the narrow grey window in the comparison image


def noise_shape(px):
    """quantum-noise power per band relative to the last band: white noise through the detector blur and sampling
    (sum over aliases of MTF^2 per axis), averaged over angle"""
    fs = 1 / px; s1 = lambda f: sum(np.interp(np.abs(f + k * fs), *MTF) ** 2 for k in (-2, -1, 0, 1, 2))
    th = (np.arange(64) + 0.5) * np.pi / 128
    rad = lambda fr: float(np.mean(s1(fr * np.cos(th)) * s1(fr * np.sin(th))))
    g = np.array([np.mean([rad(f) for f in np.linspace(a, b, 9)]) for a, b in zip(EDGES[:-1], EDGES[1:])])
    return g / g[-1]


# ----------------------------------------------------------------------------------------------- finding the images
def scan(roots, ext, skip=()):
    """LE / HE / recombined triplets below the given folders, paired by the irradiation-event UIDs in the headers;
    folders whose name contains one of the 'skip' texts are not entered, and copies of the same image are used once"""
    recs = []; seen = set()
    for root in roots:
        root = os.path.abspath(root)
        for d, dirs, files in os.walk(root):
            dirs[:] = sorted(x for x in dirs if not x.startswith('.') and not re.search(r'_strength_[0-9.]+$', x) and not any(k.lower() in x.lower() for k in skip))
            for f in sorted(files):
                if f.startswith('.') or not f.lower().endswith(ext):
                    continue
                p = os.path.join(d, f)
                try:
                    ds = pydicom.dcmread(p, stop_before_pixels=True)
                except Exception:
                    continue
                it = [str(x) for x in ds.get('ImageType', [])]; intent = str(ds.get('PresentationIntentType', ''))
                kind = ('R' if 'RECOMBINED' in it else 'L' if 'LOW_ENERGY' in it and 'FOR PROCESSING' in intent
                        else 'H' if 'HIGH_ENERGY' in it and 'FOR PROCESSING' in intent else None)
                uid = str(ds.get('SOPInstanceUID', p))
                if kind is None or uid in seen:
                    continue
                seen.add(uid)
                iu = ds.get((0x0008, 0x3010)); ev = [str(u) for u in (iu.value if iu is not None and iu.VM > 1 else ([iu.value] if iu is not None else []))]
                pv = ds.get((0x0019, 0x1010)); m = re.search(r'_W\(([0-9.]+)\)', str(pv.value)) if pv is not None else None
                recs.append(dict(path=p, root=root, kind=kind, ev=ev, study=str(ds.get('StudyInstanceUID', '')), kv=float(ds.get('KVP') or 0),
                                 mAs=float(ds.get('ExposureInuAs') or 0) / 1000, W=float(m.group(1)) if m else None, series=int(ds.get('SeriesNumber') or 0)))
    by = {(r['study'], e): r for r in recs if r['kind'] in ('L', 'H') for e in r['ev']}; out = []
    for r in recs:
        if r['kind'] != 'R':
            continue
        L = next((by[(r['study'], e)] for e in r['ev'] if (r['study'], e) in by and by[(r['study'], e)]['kind'] == 'L'), None)
        H = next((by[(r['study'], e)] for e in r['ev'] if (r['study'], e) in by and by[(r['study'], e)]['kind'] == 'H'), None)
        if not L or not H:
            continue
        parts = os.path.relpath(r['path'], r['root']).split(os.sep); chain = [os.path.basename(r['root'])] + parts[:-1]
        main = next((c for c in chain if re.match(r'^\d+\s*kV', c, re.I)), parts[0] if len(parts) > 1 else chain[0])
        session = parts[-4] if len(parts) >= 4 else chain[0]
        out.append(dict(L=L['path'], H=H['path'], R=r['path'], kv=L['kv'], mAsL=L['mAs'], mAsH=H['mAs'], W=r['W'], main=re.sub(r'\s+', '', main),
                        session=session, id='%s_%04d' % (session, r['series'])))
    return sorted(out, key=lambda t: (t['main'], t['session'], t['id']))


# ----------------------------------------------------------------------------------------------- spectra and weights
def footprint(H):
    """the phantom: HE below 60 % of the open-field level, largest connected region, holes filled"""
    m = ndi.binary_opening(H < 0.6 * np.percentile(H, 99.5), iterations=3); lab, n = ndi.label(m)
    if n == 0:
        raise RuntimeError('no object found in the HE image')
    return ndi.binary_fill_holes(lab == (np.argmax(ndi.sum(m, lab, range(1, n + 1))) + 1))


class Rois:
    """Hann-windowed, quadratically detrended N px ROIs inside a mask; per ROI and band the sums of the auto and cross
    power of Lam_L, Lam_H and R"""
    def __init__(self, px):
        self.win = np.outer(np.hanning(N), np.hanning(N)); yv, xv = np.mgrid[0:N, 0:N]
        Ad = np.c_[np.ones(N * N), xv.ravel(), yv.ravel(), (xv * xv).ravel(), (yv * yv).ravel(), (xv * yv).ravel()]
        self.Ad, self.pinv = Ad, np.linalg.pinv(Ad)
        f = np.fft.fftfreq(N, d=px); FR = np.hypot(f[:, None], f[None, :])
        self.bands = [np.flatnonzero(((FR >= a) & (FR < b)).ravel()) for a, b in zip(EDGES[:-1], EDGES[1:])]
        self.norm = px * px / (self.win ** 2).sum() / np.array([len(b) for b in self.bands])     # band sum -> mean PSD (mm2)

    def ft(self, a):
        return np.fft.fft2(self.win * (a - (self.Ad @ (self.pinv @ a.ravel())).reshape(N, N))).ravel()

    def positions(self, ok):
        return [(y0, x0) for y0 in range(0, ok.shape[0] - N, N // 2) for x0 in range(0, ok.shape[1] - N, N // 2) if ok[y0:y0 + N, x0:x0 + N].all()]

    def stats(self, lL, lH, R, pos):
        keys = ('LL', 'HH', 'HL', 'RH', 'RL', 'RR'); S = {k: np.zeros((len(pos), len(self.bands))) for k in keys}
        for i, (y0, x0) in enumerate(pos):
            sl = (slice(y0, y0 + N), slice(x0, x0 + N)); l, h, r = self.ft(lL[sl]), self.ft(lH[sl]), self.ft(R[sl])
            prod = dict(LL=np.abs(l) ** 2, HH=np.abs(h) ** 2, HL=(h * np.conj(l)).real, RH=(r * np.conj(h)).real, RL=(r * np.conj(l)).real, RR=np.abs(r) ** 2)
            for k in keys:
                S[k][i] = [prod[k][b].sum() for b in self.bands]
        return S


def select(S):
    """ROIs kept for the weights: those whose residual at the cancelling weight (0.1-1 cy/mm) is not an outlier
    (iodine discs, bubbles, markers give content that LE does not predict); robust, upper tail only"""
    keep = np.ones(len(S['LL']), bool); mid = (FC >= 0.1) & (FC < 1.0)
    for _ in range(8):
        ws = S['HL'][keep].sum(0) / S['LL'][keep].sum(0)
        le = np.log(np.clip((S['HH'] - 2 * ws * S['HL'] + ws ** 2 * S['LL'])[:, mid].sum(1), 1e-30, None))
        med = np.median(le[keep]); mad = 1.4826 * np.median(np.abs(le[keep] - med)) + 1e-9
        new = le < med + 3.5 * mad
        if (new == keep).all():
            break
        keep = new
    return keep


def measure(S, keep, norm, nshape):
    """per band, from the kept ROIs: PSDs, cancelling weight w*, coherence, Siemens' gain alpha and effective weight"""
    s = {k: v[keep].sum(0) for k, v in S.items()}; n = int(keep.sum())
    LL, HH, HL = (s[k] * norm / n for k in ('LL', 'HH', 'HL'))
    det = s['HH'] * s['LL'] - s['HL'] ** 2
    alpha = (s['RH'] * s['LL'] - s['RL'] * s['HL']) / det; beta = (s['RL'] * s['HH'] - s['RH'] * s['HL']) / det
    R2 = (alpha * s['RH'] + beta * s['RL']) / s['RR']
    ws = HL / LL; C = HL ** 2 / LL
    E = np.clip((HH - C) - (HH[-1] + ws ** 2 * LL[-1]) * nshape, 0, None)      # non-noise residual already there at w*
    return dict(n=n, LL=LL, HH=HH, HL=HL, ws=ws, coh=C / HH, C=C, E=E, alpha=alpha, beta=beta, R2=R2, w_eff=-beta / alpha)


def design(m, strength):
    """left-in share rho(f) and tuned weight w_t(f) for one strength (see the header)"""
    T = np.array(BREAST['T']); lo = FC < FADE[0]
    rho0 = strength * float(np.sqrt(max(np.sum(T[lo] - m['E'][lo]), 0.0) / np.sum(m['C'][lo])))     # share at strength 1, scaled
    fade = np.clip((np.log10(FC) - np.log10(FADE[0])) / (np.log10(FADE[1]) - np.log10(FADE[0])), 0, 1)
    rho = min(rho0, 0.6) * 0.5 * (1 + np.cos(np.pi * fade)); wt = m['ws'] * (1 - rho)
    c = np.where(m['R2'] >= 0.2, m['alpha'] * wt + m['beta'], 0.0)             # change of Siemens' image per unit Lam_L
    pd = lambda w: m['HH'] - 2 * w * m['HL'] + w ** 2 * m['LL']
    return dict(rho0=min(rho0, 0.6), rho=rho, wt=wt, c=c, PD=pd(wt), PD_full=pd(m['ws']), PD_siemens=pd(np.where(m['R2'] >= 0.2, m['w_eff'], 0.0)))


class Padded:
    """FFT of a log image (smoothly continued outside the phantom, padded), to apply radial filters given per band"""
    def __init__(self, img, fp, px):
        m = ndi.binary_dilation(fp, iterations=10)
        num = ndi.gaussian_filter(np.where(m, img, 0.0), 60); den = ndi.gaussian_filter(m.astype(float), 60)
        f = np.where(m, img, np.where(den > 1e-3, num / np.maximum(den, 1e-3), img[m].mean()))
        self.shape = img.shape; self.pad = 256
        ny, nx = (sfft.next_fast_len(n + 2 * self.pad, real=True) for n in img.shape); self.n = (ny, nx)
        fpad = np.pad(f, ((self.pad, ny - img.shape[0] - self.pad), (self.pad, nx - img.shape[1] - self.pad)), mode='reflect')
        self.F = sfft.rfft2(fpad, workers=-1)
        self.FR = np.hypot(sfft.fftfreq(ny, d=px)[:, None], sfft.rfftfreq(nx, d=px)[None, :])

    def filtered(self, cb):
        g = lambda f: f ** 2 / (f ** 2 + F0 ** 2)
        filt = np.interp(np.log10(np.clip(self.FR, FC[0], FC[-1])), np.log10(FC), cb)
        low = self.FR < FC[0]; filt[low] = cb[0] * g(self.FR[low]) / g(FC[0])
        p = self.pad; return sfft.irfft2(self.F * filt, s=self.n, workers=-1)[p:p + self.shape[0], p:p + self.shape[1]]


# ----------------------------------------------------------------------------------------------- output
def binned(a, k=3):
    h, w = (a.shape[0] // k) * k, (a.shape[1] // k) * k
    return a[:h, :w].reshape(h // k, k, w // k, k).mean((1, 3))


def save_dicom(ds, img, label, path):
    ds.PixelData = np.ascontiguousarray(img, dtype='<u2').tobytes()
    ds.SOPInstanceUID = generate_uid(); ds.file_meta.MediaStorageSOPInstanceUID = ds.SOPInstanceUID; ds.SeriesInstanceUID = generate_uid()
    ds.SeriesDescription = ('%s, L2 %s' % (str(ds.get('SeriesDescription', 'Insight CEM')), label))[:64]
    ds.ImageComments = 'L2 recombination, %s (research; not Siemens output)' % label
    ds.save_as(path)


def fig_comparison(R, X, win, med, name, head, path):
    """Siemens' recombined image, ours and the difference (ours minus Siemens); top row in Siemens' default window,
    bottom row in a narrow one; the difference on the same grey scale (same window width), centred on zero"""
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    rows = np.flatnonzero((R > 0).any(1)); cols = np.flatnonzero((R > 0).any(0)); k = 30
    sl = (slice(max(rows.min() - k, 0), rows.max() + k), slice(max(cols.min() - k, 0), cols.max() + k))
    D = np.where(R > 0, X - R, 0.0)
    fig, ax = plt.subplots(2, 3, figsize=(18, 16.5))
    for i, (c, w) in enumerate(((win[0], win[1]), (med, NARROW))):
        for j, (a, c0) in enumerate(((R, c), (X, c), (D, 0.0))):
            ax[i, j].imshow(binned(a[sl]), cmap='gray', vmin=c0 - w / 2, vmax=c0 + w / 2, interpolation='antialiased'); ax[i, j].set_xticks([]); ax[i, j].set_yticks([])
    for j, n in enumerate(('Siemens recombined image', name, 'difference: ours minus Siemens (same grey scale, centred on zero)')):
        ax[0, j].set_title(n, fontsize=12, loc='left')
    ax[0, 0].set_ylabel("Siemens' default window (centre %.0f, width %.0f)" % win, fontsize=10)
    ax[1, 0].set_ylabel('narrow window (width %.0f, centred on the interior)' % NARROW, fontsize=10)
    fig.suptitle(head, fontsize=11, x=0.01, ha='left')
    fig.tight_layout(rect=(0, 0, 1, 0.965)); fig.savefig(path, dpi=80); plt.close(fig)


def fig_curves(t, m, d, strength, achieved, path):
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    ok = m['R2'] >= 0.2
    fig, ax = plt.subplots(1, 3, figsize=(18, 5.2))
    ax[0].plot(FC[ok], m['w_eff'][ok], 'o-', color='#d62728', label="Siemens' effective weight (measured from R)")
    ax[0].plot(FC, m['ws'], 'o-', color='k', label='w*(f): cancels the spheres'); ax[0].plot(FC, d['wt'], 'o-', color='#2ca02c', label='tuned w_t(f), strength %g' % strength)
    if achieved is not None:
        ax[0].plot(FC[ok], achieved[ok], 'x', color='#1f77b4', ms=7, label='measured in the tuned image (check)')
    ax[0].set_title('weight on the LE image per spatial frequency', fontsize=10); ax[0].set_ylabel('weight')
    ax[1].plot(FC, 100 * d['rho'], 'o-', color='#2ca02c', label='tuned, strength %g' % strength)
    ax[1].plot(FC[ok], 100 * (1 - m['w_eff'][ok] / m['ws'][ok]), 'o-', color='#d62728', label='Siemens')
    ax[1].set_title('share of the sphere contrast left in (%)', fontsize=10); ax[1].set_ylabel('%'); ax[1].set_ylim(-5, 105)
    ax[2].fill_between(FC, BREAST['PD_LO'], BREAST['PD_HI'], color='#1f77b4', alpha=0.2, lw=0, label='real breasts under Siemens (9 views, lowest-highest)')
    ax[2].plot(FC, BREAST['PD'], color='#1f77b4', label='real breasts, median')
    for y, c, lab in ((d['PD_siemens'], '#d62728', 'L2 under Siemens'), (d['PD_full'], 'k', 'L2 full cancellation'), (d['PD'], '#2ca02c', 'L2 tuned, strength %g' % strength)):
        ax[2].plot(FC, y, 'o-', ms=3, color=c, label=lab)
    ax[2].set_yscale('log'); ax[2].set_title('spectrum of the log-signal residual Lam_H - w(f) Lam_L (mm$^2$), interior', fontsize=10)
    for a in ax:
        a.set_xscale('log'); a.set_xlabel('spatial frequency (cy/mm)'); a.legend(fontsize=8); a.grid(alpha=0.25)
    fig.suptitle('%s   %.0f kV, %.0f / %.0f mAs   (%d regions used)' % (t['id'], t['kv'], t['mAsL'], t['mAsH'], m['n']), fontsize=11, x=0.01, ha='left')
    fig.tight_layout(rect=(0, 0, 1, 0.94)); fig.savefig(path, dpi=85); plt.close(fig)


def write_curves(m, d, achieved, path):
    with open(path, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['band_from_cy_mm', 'band_to_cy_mm', 'band_centre_cy_mm', 'siemens_effective_weight', 'siemens_gain_alpha', 'siemens_fit_R2', 'w_star', 'coherence',
                    'left_in_share_rho', 'tuned_weight', 'tuned_weight_measured', 'residual_siemens_mm2', 'residual_full_cancellation_mm2', 'residual_tuned_mm2',
                    'breast_median_mm2', 'breast_lowest_mm2', 'breast_highest_mm2'])
        for i in range(len(FC)):
            w.writerow([EDGES[i], EDGES[i + 1], '%.4f' % FC[i], '%.5f' % m['w_eff'][i], '%.1f' % m['alpha'][i], '%.4f' % m['R2'][i], '%.5f' % m['ws'][i], '%.4f' % m['coh'][i],
                        '%.5f' % d['rho'][i], '%.5f' % d['wt'][i], '' if achieved is None else '%.5f' % achieved[i], '%.4e' % d['PD_siemens'][i], '%.4e' % d['PD_full'][i],
                        '%.4e' % d['PD'][i], '%.4e' % BREAST['PD'][i], '%.4e' % BREAST['PD_LO'][i], '%.4e' % BREAST['PD_HI'][i]])


def fig_check_curves(t, m, mo, psd, path):
    """curves of the Siemens-weighting check: effective weight in Siemens' image and in ours, the share of Siemens'
    image that ours reproduces per band, and the power spectra of both images and of their difference"""
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    ok = m['R2'] >= 0.2
    fig, ax = plt.subplots(1, 3, figsize=(18, 5.2))
    ax[0].plot(FC[ok], m['w_eff'][ok], 'o-', color='#d62728', label="Siemens' recombined image")
    ax[0].plot(FC[ok], mo['w_eff'][ok], 'x', color='#1f77b4', ms=8, label="our recombination at Siemens' weighting")
    ax[0].set_title('effective weight on the LE image per spatial frequency', fontsize=10); ax[0].set_ylabel('weight')
    ax[1].plot(FC, 100 * (1 - psd['D'] / psd['R']), 'o-', color='#2ca02c', label="share of Siemens' image that ours reproduces")
    ax[1].set_title("agreement per spatial frequency: 100 x (1 - power of the difference / power of Siemens' image)", fontsize=10); ax[1].set_ylabel('%'); ax[1].set_ylim(-5, 105)
    for y, c, mk, lab in ((psd['R'], '#d62728', 'o-', "Siemens' recombined image"), (psd['O'], '#1f77b4', 'x-', 'our recombination'), (psd['D'], 'k', 'o-', 'difference')):
        ax[2].plot(FC, y, mk, ms=4, color=c, label=lab)
    ax[2].set_yscale('log'); ax[2].set_title('power spectrum of the images (grey levels$^2$ mm$^2$), interior', fontsize=10)
    for a in ax:
        a.set_xscale('log'); a.set_xlabel('spatial frequency (cy/mm)'); a.legend(fontsize=8); a.grid(alpha=0.25)
    fig.suptitle("%s   %.0f kV, %.0f / %.0f mAs   check: our recombination at Siemens' weighting (%d regions used)" % (t['id'], t['kv'], t['mAsL'], t['mAsH'], m['n']), fontsize=11, x=0.01, ha='left')
    fig.tight_layout(rect=(0, 0, 1, 0.94)); fig.savefig(path, dpi=85); plt.close(fig)


def write_check_curves(m, mo, psd, path):
    with open(path, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['band_from_cy_mm', 'band_to_cy_mm', 'band_centre_cy_mm', 'siemens_effective_weight', 'ours_effective_weight', 'siemens_gain_alpha', 'ours_gain_alpha',
                    'fit_R2_in_regions', 'power_siemens', 'power_ours', 'power_difference', 'share_reproduced_percent'])
        for i in range(len(FC)):
            w.writerow([EDGES[i], EDGES[i + 1], '%.4f' % FC[i], '%.5f' % m['w_eff'][i], '%.5f' % mo['w_eff'][i], '%.1f' % m['alpha'][i], '%.1f' % mo['alpha'][i], '%.4f' % m['R2'][i],
                        '%.4e' % psd['R'][i], '%.4e' % psd['O'][i], '%.4e' % psd['D'][i], '%.1f' % (100 * (1 - psd['D'][i] / psd['R'][i]))])


# ----------------------------------------------------------------------------------------------- one exposure
CHECK = 'siemens_weighting'          # folder and file label of the check image
EXT_T = ('_tuned.dcm', '_comparison.png', '_curves.png', '_curves.csv')
EXT_C = ('_%s.dcm' % CHECK, '_comparison.png', '_curves.png', '_curves.csv')


def out_base(t, s, out_root):
    name = CHECK if s == CHECK else 'strength_' + (('%.2f' % s).rstrip('0').rstrip('.') if s != int(s) else '%.1f' % s)
    folder = os.path.join(out_root, '%s_%s' % (t['main'], name), t['session'])
    return folder, os.path.join(folder, t['id'])


def process(job):
    t, strengths, out_root, overwrite, check = job; t0 = time.time(); rows = []
    if not overwrite:            # what is already there is not made again
        strengths = [s for s in strengths if not all(os.path.exists(out_base(t, s, out_root)[1] + e) for e in EXT_T)]
        check = check and not all(os.path.exists(out_base(t, CHECK, out_root)[1] + e) for e in EXT_C)
        if not strengths and not check:
            return rows
    try:
        dsL, dsH, dsR = (pydicom.dcmread(t[k]) for k in ('L', 'H', 'R'))
        L, H, R = (ds.pixel_array.astype(np.float64) for ds in (dsL, dsH, dsR))
        px = float(dsL.ImagerPixelSpacing[0]); vmax = 2 ** int(dsR.BitsStored) - 1
        lL, lH = -np.log(np.clip(L, 1, None)), -np.log(np.clip(H, 1, None))
        fp = footprint(H); dist = ndi.distance_transform_edt(fp) * px; inner = dist > MARGIN_MM; valid = (R > 0) & (R < vmax)
        rois = Rois(px); pos = rois.positions(inner & valid)
        if len(pos) < 30:
            raise RuntimeError('only %d interior regions found' % len(pos))
        nshape = noise_shape(px); S = rois.stats(lL, lH, R, pos); keep = select(S); m = measure(S, keep, rois.norm, nshape)
        le = Padded(lL, fp, px); med = float(np.median(R[inner & valid]))
        wc, ww = dsR.get('WindowCenter', 2600), dsR.get('WindowWidth', 2100)
        win = (float(wc[0] if isinstance(wc, pydicom.multival.MultiValue) else wc), float(ww[0] if isinstance(ww, pydicom.multival.MultiValue) else ww))
        kept = [p for p, k in zip(pos, keep) if k]; allk = np.ones(len(kept), bool); lo = slice(0, 2)
        common = dict(id=t['id'], session=t['session'], kV=t['kv'], mAs_LE=t['mAsL'], mAs_HE=t['mAsH'], regions_used=m['n'], regions_left_out=int((~keep).sum()),
                      siemens_header_W=t['W'], siemens_effective_weight=round(float(m['w_eff'][lo].mean()), 4), w_star=round(float(m['ws'][lo].mean()), 4), error='')
        head = '%s   %.0f kV, %.0f / %.0f mAs\n' % (t['id'], t['kv'], t['mAsL'], t['mAsH'])
        for s in strengths:
            d = design(m, s)
            Rt = np.where(R > 0, np.clip(np.rint(R - le.filtered(d['c'])), 0, vmax), 0)
            a = measure(rois.stats(lL, lH, Rt, kept), allk, rois.norm, nshape)['w_eff']                    # check on the result
            folder, base = out_base(t, s, out_root); os.makedirs(folder, exist_ok=True)
            save_dicom(pydicom.dcmread(t['R']), Rt, 'tuned s=%g' % s, base + '_tuned.dcm')
            fig_comparison(R, Rt, win, med, 'ours: tuned, strength %g' % s, head + 'weight on LE at coarse scales: Siemens %.3f, cancelling w* %.3f, tuned %.3f '
                           '(%.0f %% of the sphere contrast left; Siemens leaves %.0f %%)' % (m['w_eff'][lo].mean(), m['ws'][lo].mean(), d['wt'][lo].mean(), 100 * d['rho0'],
                                                                                         100 * (1 - m['w_eff'][lo].mean() / m['ws'][lo].mean())), base + '_comparison.png')
            fig_curves(t, m, d, s, a, base + '_curves.png'); write_curves(m, d, a, base + '_curves.csv')
            rows.append(dict(common, folder=os.path.dirname(folder), strength=s, tuned_weight=round(float(d['wt'][lo].mean()), 4), tuned_weight_measured=round(float(a[lo].mean()), 4),
                             left_in_share_percent=round(100 * d['rho0'], 2), siemens_left_in_share_percent=round(100 * float(1 - m['w_eff'][lo].mean() / m['ws'][lo].mean()), 2),
                             residual_0p2_siemens_mm2='%.3e' % d['PD_siemens'][3], residual_0p2_tuned_mm2='%.3e' % d['PD'][3], breast_0p2_mm2='%.3e' % BREAST['PD'][3],
                             seconds=round(time.time() - t0, 1)))
        if check:
            # our recombination of the raw LE and HE images at Siemens' weighting: gain alpha(f) and weight w_eff(f) per band as
            # measured (R = alpha Lam_H + beta Lam_L, beta = -alpha w_eff), applied to the whole images; offset from the interior
            ok = m['R2'] >= 0.2; rec = Padded(lH, fp, px).filtered(np.where(ok, m['alpha'], 0.0)) + le.filtered(np.where(ok, m['beta'], 0.0))
            ins = valid & (dist > 15); edge = valid & fp & (dist <= 15)
            Ro = np.where(R > 0, np.clip(np.rint(rec + np.median(R[ins] - rec[ins])), 0, vmax), 0)
            mo = measure(rois.stats(lL, lH, Ro, kept), allk, rois.norm, nshape); n = len(kept)
            psd = dict(R=S['RR'][keep].sum(0) * rois.norm / n, O=rois.stats(lL, lH, Ro, kept)['RR'].sum(0) * rois.norm / n,
                       D=rois.stats(lL, lH, Ro - R, kept)['RR'].sum(0) * rois.norm / n)
            st = dict(corr=float(np.corrcoef(R[ins], Ro[ins])[0, 1]), rmse=float(np.sqrt(np.mean((Ro[ins] - R[ins]) ** 2))), sd=float(R[ins].std()),
                      edge=float(np.corrcoef(R[edge], Ro[edge])[0, 1]), share=float(100 * (1 - psd['D'][:8].sum() / psd['R'][:8].sum())))
            folder, base = out_base(t, CHECK, out_root); os.makedirs(folder, exist_ok=True)
            save_dicom(pydicom.dcmread(t['R']), Ro, "our recombination at Siemens' weighting", base + EXT_C[0])
            fig_comparison(R, Ro, win, med, "ours: recombined from raw LE and HE at Siemens' weighting", head + "interior (more than 15 mm inside the outline): correlation %.3f, "
                           "difference %.0f grey levels rms (image SD %.0f); outer 15 mm: correlation %.2f; %.0f %% of Siemens' image reproduced at 0.06-1 cy/mm in the regions used"
                           % (st['corr'], st['rmse'], st['sd'], st['edge'], st['share']), base + '_comparison.png')
            fig_check_curves(t, m, mo, psd, base + '_curves.png'); write_check_curves(m, mo, psd, base + '_curves.csv')
            rows.append(dict(common, folder=os.path.dirname(folder), strength="Siemens' weighting (check)", check_interior_correlation=round(st['corr'], 4),
                             check_interior_difference_rms=round(st['rmse'], 1), check_interior_image_sd=round(st['sd'], 1), check_edge_15mm_correlation=round(st['edge'], 3),
                             check_share_reproduced_percent=round(st['share'], 1), check_ours_effective_weight=round(float(mo['w_eff'][lo].mean()), 4), seconds=round(time.time() - t0, 1)))
    except Exception as ex:
        rows.append(dict(folder=os.path.join(out_root, '%s_strength_error' % t['main']), id=t['id'], session=t['session'], kV=t['kv'], strength='', error='%s: %s' % (type(ex).__name__, ex)))
        traceback.print_exc()
    return rows


EXCEL = 'L2_tuned_recombination_numbers.xlsx'


def write_excel(out_root):
    """all numbers of all analyses in one Excel file in the output folder, rebuilt from the CSV files of this and earlier
    runs: one line per image (tuned images; check at Siemens' weighting) and one line per image and frequency band (curves)"""
    try:
        from openpyxl import Workbook
        from openpyxl.styles import Font
        from openpyxl.utils import get_column_letter
    except ImportError:
        print('Excel file not written: the Python package openpyxl is missing (pip install openpyxl). The numbers are in the CSV files.'); return

    def num(v):
        if v is None or v == '':
            return None
        try:
            return float(v)
        except ValueError:
            return v

    def read(f, extra):
        with open(f, newline='') as fh:
            return [dict(extra, **r) for r in csv.DictReader(fh)]

    sheets = {'tuned_images': [], 'siemens_weighting_check': [], 'tuned_curves': [], 'check_curves': [], 'not_processed': []}
    for name in sorted(os.listdir(out_root)):
        folder = os.path.join(out_root, name)
        if not os.path.isdir(folder) or not os.path.exists(os.path.join(folder, 'summary.csv')):
            continue
        if name.endswith('_strength_error'):
            sheets['not_processed'] += read(os.path.join(folder, 'summary.csv'), dict(folder=name)); continue
        check = name.endswith('_' + CHECK); m = re.match(r'^(.*)_strength_([0-9.]+)$', name)
        if not check and not m:
            continue
        main = name[:-len(CHECK) - 1] if check else m.group(1)
        sheets['siemens_weighting_check' if check else 'tuned_images'] += read(os.path.join(folder, 'summary.csv'), dict(main_folder=main, output_folder=name))
        for session in sorted(os.listdir(folder)):
            d = os.path.join(folder, session)
            if os.path.isdir(d):
                for f in sorted(os.listdir(d)):
                    if f.endswith('_curves.csv'):
                        extra = dict(main_folder=main, output_folder=name, session=session, id=f[:-len('_curves.csv')])
                        if not check:
                            extra['strength'] = float(m.group(2))
                        sheets['check_curves' if check else 'tuned_curves'] += read(os.path.join(d, f), extra)
    notes = [('sheet', 'content'),
             ('tuned_images', 'one line per image and strength: weights at coarse scales (mean of 0.06-0.15 cy/mm), share of the sphere contrast left in, residual power at 0.2-0.3 cy/mm'),
             ('siemens_weighting_check', "one line per image: our recombination of raw LE and HE at Siemens' weighting against Siemens' recombined image (option --siemens-weighting)"),
             ('tuned_curves', 'one line per image, strength and spatial-frequency band: the curves of <id>_curves.png'),
             ('check_curves', 'one line per image and band: the curves of the check'),
             ('not_processed', 'images that could not be processed, with the error'),
             ('', ''), ('quantity', 'meaning'),
             ('siemens_effective_weight', "weight on the LE image that Siemens' recombined image contains, measured from it (R = alpha (Lam_H - w_eff Lam_L))"),
             ('w_star', 'weight that cancels the sphere background'), ('tuned_weight / tuned_weight_measured', 'weight of the tuned image as designed / as measured in the tuned image'),
             ('left_in_share', 'share of the sphere contrast left in the image, 1 - weight / w_star'),
             ('residual ... mm2', 'power spectral density of the log-signal residual Lam_H - w Lam_L in the interior (Lam = -ln of the stored value)'),
             ('breast ... mm2', "the same quantity in real breasts under Siemens' recombination (9 views of 2 patients; aggregated values)"),
             ('check_interior_correlation / difference_rms / image_sd', "pixel-wise agreement more than 15 mm inside the outline, in grey levels of the recombined image"),
             ('check_share_reproduced_percent', "100 x (1 - power of the difference / power of Siemens' image) at 0.06-1 cy/mm in the regions used"),
             ('check_edge_15mm_correlation', "agreement in the outer 15 mm, where Siemens applies processing that is not part of our method"),
             ('power_siemens / power_ours / power_difference', 'power spectral density of the recombined images (grey levels squared times mm2), interior')]
    wb = Workbook(); wb.remove(wb.active)
    for title, rows in sheets.items():
        if not rows:
            continue
        ws = wb.create_sheet(title); keys = list(dict.fromkeys(k for r in rows for k in r))
        ws.append(keys)
        for r in rows:
            ws.append([num(r.get(k)) for k in keys])
        for c in ws[1]:
            c.font = Font(bold=True)
        ws.freeze_panes = 'A2'; ws.auto_filter.ref = ws.dimensions
        for j, k in enumerate(keys, 1):
            ws.column_dimensions[get_column_letter(j)].width = min(max(len(k) + 2, 12), 46)
    ws = wb.create_sheet('notes')
    for r in notes:
        ws.append(list(r))
    ws.column_dimensions['A'].width = 52; ws.column_dimensions['B'].width = 150
    for row in (1, 8):
        for c in ws[row]:
            c.font = Font(bold=True)
    path = os.path.join(out_root, EXCEL)
    try:
        wb.save(path); print('written', path)
    except OSError as ex:
        print('Excel file not written (%s). Close it if it is open and run the script again; the numbers are in the CSV files.' % ex)


def main():
    ap = argparse.ArgumentParser(description='L2 phantom: tuned recombination of Siemens TiCEM images (see the header of this file).')
    ap.add_argument('inputs', nargs='*', default=['.'], help='folders with the DICOM files (default: the current folder)')
    ap.add_argument('--strengths', type=float, nargs='+', default=[1.0], help='strengths to produce (default 1.0); 0 = full cancellation')
    ap.add_argument('--siemens-weighting', action='store_true', help="also recombine the raw LE and HE images with our method at Siemens' weighting and compare "
                    "the result with Siemens' image (check of the method) -> <main>_siemens_weighting/")
    ap.add_argument('--out', default='.', help='where the <main>_strength_<s> folders are written (default: the current folder)')
    ap.add_argument('--ext', default='.ima', help='file ending of the DICOM files (default .IMA)')
    ap.add_argument('--only', default=None, help='process only exposures whose path contains this text')
    ap.add_argument('--skip', nargs='*', default=['4AFC'], help='folders whose name contains one of these texts are not searched (default: 4AFC)')
    ap.add_argument('--limit', type=int, default=None, help='process at most this many exposures per session (for a first look)')
    ap.add_argument('--workers', type=int, default=1, help='exposures processed in parallel (default 1; usually not faster from a USB drive, the FFTs already use all cores)')
    ap.add_argument('--overwrite', action='store_true', help='make the files again even if they exist (default: existing results are kept, so a run can be continued)')
    a = ap.parse_args()
    T = scan(a.inputs, a.ext.lower(), a.skip)
    if a.only:
        T = [t for t in T if a.only.lower() in t['R'].lower()]
    if a.limit:
        seen = {}; T = [t for t in T if seen.setdefault(t['session'], []).append(1) or len(seen[t['session']]) <= a.limit]
    if not T:
        sys.exit('no LE / HE / recombined triplets found below: %s' % ', '.join(a.inputs))
    print('%d exposures in %d sessions; strengths %s' % (len(T), len(set(t['session'] for t in T)), a.strengths), flush=True)
    jobs = [(t, a.strengths, os.path.abspath(a.out), a.overwrite, a.siemens_weighting) for t in T]; rows = []; t_start = time.time()
    if a.workers > 1:
        from multiprocessing import Pool
        with Pool(a.workers) as pool:
            it = pool.imap(process, jobs)
            for k, r in enumerate(it):
                rows += r; report(k, len(jobs), r)
    else:
        for k, j in enumerate(jobs):
            r = process(j); rows += r; report(k, len(jobs), r)
    keys = ['id', 'session', 'kV', 'mAs_LE', 'mAs_HE', 'strength', 'regions_used', 'regions_left_out', 'siemens_header_W', 'siemens_effective_weight', 'w_star', 'tuned_weight',
            'tuned_weight_measured', 'left_in_share_percent', 'siemens_left_in_share_percent', 'residual_0p2_siemens_mm2', 'residual_0p2_tuned_mm2', 'breast_0p2_mm2', 'seconds', 'error']
    keys_check = keys[:5] + ['strength', 'regions_used', 'regions_left_out', 'siemens_header_W', 'siemens_effective_weight', 'check_ours_effective_weight', 'check_interior_correlation',
                             'check_interior_difference_rms', 'check_interior_image_sd', 'check_share_reproduced_percent', 'check_edge_15mm_correlation', 'seconds', 'error']
    for folder in sorted(set(r['folder'] for r in rows)):
        os.makedirs(folder, exist_ok=True)
        f = os.path.join(folder, 'summary.csv'); new = {r['id']: r for r in rows if r['folder'] == folder}
        old = {r['id']: r for r in csv.DictReader(open(f))} if os.path.exists(f) else {}          # lines of earlier runs are kept
        with open(f, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=keys_check if folder.endswith(CHECK) else keys, extrasaction='ignore'); w.writeheader(); w.writerows(v for _, v in sorted({**old, **new}.items()))
        print('written', folder)
    write_excel(os.path.abspath(a.out))
    print('%d exposures, %.0f min' % (len(jobs), (time.time() - t_start) / 60))
    bad = [r for r in rows if r['error']]
    if bad:
        print('%d exposures failed, see the summary.csv in *_strength_error' % len(bad))


def report(k, n, r):
    if not r:
        print('[%d/%d] already done' % (k + 1, n), flush=True)
    for x in r:
        if x['error']:
            msg = 'ERROR ' + x['error']
        elif 'check_interior_correlation' in x:
            msg = "check at Siemens' weighting: interior correlation %.3f, difference %.0f rms (image SD %.0f), %.0f %% reproduced, %.0f s" % (
                x['check_interior_correlation'], x['check_interior_difference_rms'], x['check_interior_image_sd'], x['check_share_reproduced_percent'], x['seconds'])
        else:
            msg = 'strength %g: w* %.3f, tuned %.3f (measured %.3f), %.1f %% left in, %s regions, %.0f s' % (
                x['strength'], x['w_star'], x['tuned_weight'], x['tuned_weight_measured'], x['left_in_share_percent'], x['regions_used'], x['seconds'])
        print('[%d/%d] %s  %s' % (k + 1, n, x['id'], msg), flush=True)


if __name__ == '__main__':
    main()

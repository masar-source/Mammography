#!/usr/bin/env python3
"""
L2 phantom: tuned recombination of contrast-enhanced mammography images, Siemens TiCEM and Hologic CEDM.
Multi-system version of l2_tuned_recombination.py (Siemens only). The method is the same, step by step; on Siemens images
it gives the same tuned images as that script. What depends on the system is listed in SYSTEMS below.

What it does, for every exposure found (raw LE, raw HE and the system's own recombined image R):
  the system's recombined image is kept exactly as it is (walls, edges, mask, grey scale) and only the weight on the
  background is changed: the sphere background is cancelled except for a soft remainder at coarse scales.
      R_tuned = R - IFFT{ c(f) FFT(Lam_L) },   c(f) = alpha(f) [ w_t(f) - w_eff(f) ],   Lam = -ln(stored value)
  alpha(f), w_eff(f)  the system's gain and effective weight per spatial-frequency band, measured from R itself by
                      regressing R on Lam_H and Lam_L in interior regions (R = alpha (Lam_H - w_eff Lam_L) + rest);
  w*(f)               the weight that cancels the spheres, S_HL/S_LL, from the same regions;
  w_t(f)              = w*(f) [1 - rho(f)]: the tuned weight. rho is the share of the sphere contrast left in:
                      one value rho0 at coarse scales (to 0.3 cy/mm), faded to zero at 1 cy/mm.
  Strength 1.0        rho0 is set so that the residual structure at 0.06-0.3 cy/mm has the power that Siemens'
                      recombination leaves in real breasts (BREAST table: aggregated numbers from 9 views of 2 patients;
                      no patient data). Other strengths scale the left-in share: 0 = full cancellation, 2 = double.
Everything is measured from each image on its own; nothing is averaged over images. Regions that contain iodine discs,
bubbles or other content that LE does not predict are found and left out automatically (no disc template needed).

Differences between the systems (all other numbers are the same; lengths are kept in mm, so pixel counts scale):
                      Siemens MAMMOMAT Revelation           Hologic Selenia Dimensions
  recombined image    ImageType RECOMBINED                  ImageType SUBTRACTION (DES)
  pairing             irradiation-event UIDs in R           R carries the HE event only; LE = same acquisition time
  pixel, last band    0.085 mm: 4-5.88 cy/mm                0.070 mm: 4-7.14 cy/mm (Nyquist); bands below 4 cy/mm equal
  ROI size            128 px (10.9 mm)                      156 px (10.9 mm)
  noise shape g(b)    white noise through the published     measured once in the open field of the 15 Hologic
                      MTF of this detector type (NHSBSP)    exposures (hologic_noise_shape.py); no MTF at hand
  weight in header    yes (private tag, _W(...))            no (only the measured effective weight is reported)
  breast target       measured on Siemens breast views      the same Siemens table, used unchanged [assumption:
                                                            no raw Hologic breast images; see below]
  object outline      HE < 0.6 x 99.5th percentile          the same, saturated pixels (open field at the top edge)
                                                            left out of the percentile
Hologic and the breast target: there are no raw Hologic breast images, so the residual that Hologic's recombination
leaves in breasts cannot be measured. The Siemens target is used as it is, in log-signal units. Both systems set their
weight to cancel breast tissue, and in the Siemens breasts most of the residual is content that LE does not predict.
On the L2 the result is a left-in share of 6-7 % on Hologic, inside the 5.4-7.8 % found on Siemens, so using the same
share instead would give nearly the same images. This is an assumption, not a measurement.

Output, in the folder given with --out (default: where the script is run), one folder per strength and main folder:
  <main>_strength_<s>/<session>/<id>_tuned.dcm        the tuned image (copy of the system's DICOM, new pixel data)
                                 <id>_comparison.png   the system's R, ours, and the difference (ours minus R)
                                 <id>_curves.png/.csv  weights, left-in share and residual spectra of this image
  <main>_strength_<s>/summary.csv                      one line per image, with the acquisition data from the headers
  <main>_<system>_weighting/...                        check of the method (option --system-weighting): the raw LE
                                                       and HE images recombined by our method at the system's own
                                                       weighting, against the system's R
  L2_tuned_recombination_multisystem_numbers.xlsx      all numbers of all runs in this output folder
  <main> = a folder named like "28kV" on the path (Siemens), or "<System>_<kV>kV" (other systems);
  <id>   = Siemens: session + series number of R; other systems: file name of the LE image + kV.

Usage (needs Python 3 with numpy, scipy, pydicom, matplotlib; openpyxl for the Excel file):
  python3 l2_tuned_recombination_multisystem.py <folders>                       strength 1.0
  python3 l2_tuned_recombination_multisystem.py hologic --out hologic --system-weighting
  python3 l2_tuned_recombination_multisystem.py 28kV --strengths 1.0 0.5 2.0 --out results_siemens
Files: DICOM files with the ending .IMA or .dcm or without ending (option --ext). Folders named *4AFC* and the output
folders of earlier runs are not searched; an image that exists in several copies is processed once; results that exist
are kept (--overwrite makes everything anew).
Research tool: the tuned images are not the manufacturer's output and not for diagnostic use.
"""
import argparse, csv, os, re, sys, time, traceback
import numpy as np
import pydicom
from pydicom.uid import generate_uid, ExplicitVRLittleEndian
from scipy import ndimage as ndi, fft as sfft

# ----------------------------------------------------------------------------------------------- fixed numbers
BASE_EDGES = [0.06, 0.1, 0.15, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0, 4.0]     # band edges (cy/mm); the last band ends at Nyquist
# Real breasts under Siemens' recombination, per band (mm2, log-signal units): structured residual T (total minus quantum
# noise; median of 9 views, lowest and highest view) and the total residual PD (Siemens' noise included; Siemens only).
# Bands below 4 cy/mm are the same for every pixel size; T is zero from 1 cy/mm. 28 Oct / 17 Nov 2026.
BREAST = dict(
    T=[7.961e-06, 1.099e-05, 1.147e-05, 8.156e-06, 3.854e-06, 1.797e-06, 9.085e-07, 3.779e-07, 0, 0, 0, 0, 0],
    T_LO=[2.25e-06, 3.075e-06, 3.144e-06, 1.966e-06, 1.055e-06, 6.663e-07, 5.071e-07, 2.11e-07, 0, 0, 0, 0, 0],
    T_HI=[2.218e-05, 3.122e-05, 2.48e-05, 1.412e-05, 8.05e-06, 4.287e-06, 1.481e-06, 5.732e-07, 9.322e-08, 0, 0, 0, 0],
    PD=[1.027e-05, 1.328e-05, 1.373e-05, 1.037e-05, 6.005e-06, 3.854e-06, 2.844e-06, 2.163e-06, 1.692e-06, 1.327e-06, 1.12e-06, 9.037e-07, 6.037e-07],
    PD_LO=[4.633e-06, 5.46e-06, 5.439e-06, 4.235e-06, 3.274e-06, 2.787e-06, 2.501e-06, 2.033e-06, 1.462e-06, 1.153e-06, 1.001e-06, 8.347e-07, 5.647e-07],
    PD_HI=[2.457e-05, 3.372e-05, 2.704e-05, 1.631e-05, 1.007e-05, 6.22e-06, 3.516e-06, 2.574e-06, 1.82e-06, 1.387e-06, 1.167e-06, 9.7e-07, 6.738e-07])
F0 = 0.05            # cy/mm: large-scale roll-off of the correction below the first band, gain ~ f^2/(f^2 + F0^2) (fitted on both systems)
ROI_MM = 128 * 0.085 # ROI size for the spectra (10.88 mm), half overlap
MARGIN_MM = 25.5     # ROIs lie at least this far inside the phantom outline
FADE = (0.3, 1.0)    # cy/mm: the left-in share is constant to FADE[0] and zero from FADE[1]
OPEN_MM, DIL_MM, SIGMA_MM, PAD_MM = 3 * 0.085, 10 * 0.085, 60 * 0.085, 256 * 0.085    # outline opening; continuation outside the object; padding

SYSTEMS = {
    'siemens': dict(
        label='Siemens', poss="Siemens'", match='SIEMENS', recombined=('RECOMBINED',), own="Siemens' recombined image",
        mtf=(np.arange(15.0), np.array([1.00, 0.91, 0.86, 0.78, 0.67, 0.55, 0.44, 0.34, 0.25, 0.18, 0.13, 0.09, 0.06, 0.03, 0.0])),  # NHSBSP 2019
        noise_table=None,
        narrow=800.0, narrow_label='narrow window (width %.0f, centred on the interior)'),
    'hologic': dict(
        label='Hologic', poss="Hologic's", match='HOLOGIC', recombined=('SUBTRACTION',), own="Hologic's recombined image (DES)",
        mtf=None,
        noise_table=(0.07, [2.11, 2.64, 3.11, 2.63, 2.50, 2.44, 2.30, 2.20, 2.01, 1.78, 1.49, 1.23, 1.00]),   # open field, median of 15 exposures
        narrow=380.0, narrow_label='second window (width %.0f, about the log-signal width of the Siemens narrow window; centred on the interior)'),
}


def system_of(ds):
    m = str(ds.get('Manufacturer', '')).upper()
    return next((k for k, s in SYSTEMS.items() if s['match'] in m), None)


def noise_shape(px, edges, mtf):
    """quantum-noise power per band relative to the last band: white noise through the detector blur and sampling
    (sum over aliases of MTF^2 per axis), averaged over angle"""
    fs = 1 / px; s1 = lambda f: sum(np.interp(np.abs(f + k * fs), *mtf) ** 2 for k in (-2, -1, 0, 1, 2))
    th = (np.arange(64) + 0.5) * np.pi / 128
    rad = lambda fr: float(np.mean(s1(fr * np.cos(th)) * s1(fr * np.sin(th))))
    g = np.array([np.mean([rad(f) for f in np.linspace(a, b, 9)]) for a, b in zip(edges[:-1], edges[1:])])
    return g / g[-1]


class Bands:
    """spatial-frequency bands, ROI size and noise shape for one pixel size and system"""
    def __init__(self, px, system):
        s = SYSTEMS[system]; self.px = px
        self.edges = BASE_EDGES + [round(0.5 / px, 2)]; self.fc = np.sqrt(np.array(self.edges[:-1]) * np.array(self.edges[1:]))
        self.N = int(2 * round(ROI_MM / px / 2))
        if s['noise_table'] is not None:
            tpx, tab = s['noise_table']
            if abs(tpx - px) > 1e-6:
                raise RuntimeError('the noise shape of %s was measured at %g mm pixels; this image has %g mm' % (s['label'], tpx, px))
            self.g = np.array(tab, float)
        else:
            self.g = noise_shape(px, self.edges, s['mtf'])


# ----------------------------------------------------------------------------------------------- finding the images
def scan(roots, exts, skip=()):
    """LE / HE / recombined triplets below the given folders. Paired by the irradiation-event UIDs in the headers; where
    the recombined image names only the HE event (Hologic), the LE image of the same acquisition time is taken. Folders
    whose name contains one of the 'skip' texts and output folders are not entered; copies of an image are used once."""
    recs = []; seen = set(); outdir = re.compile(r'(_strength_[0-9.]+|_(%s)_weighting)$' % '|'.join(SYSTEMS))
    for root in roots:
        root = os.path.abspath(root)
        for d, dirs, files in os.walk(root):
            dirs[:] = sorted(x for x in dirs if not x.startswith('.') and not outdir.search(x) and not any(k.lower() in x.lower() for k in skip))
            for f in sorted(files):
                e = os.path.splitext(f)[1].lower()
                if f.startswith('.') or f.upper() == 'DICOMDIR' or e not in exts:
                    continue
                p = os.path.join(d, f)
                try:
                    ds = pydicom.dcmread(p, stop_before_pixels=True)
                except Exception:
                    continue
                sysk = system_of(ds)
                if sysk is None:
                    continue
                it = [str(x).upper() for x in ds.get('ImageType', [])]; intent = str(ds.get('PresentationIntentType', '')).upper()
                kind = ('R' if any(k in it for k in SYSTEMS[sysk]['recombined']) else 'L' if 'LOW_ENERGY' in it and 'FOR PROCESSING' in intent
                        else 'H' if 'HIGH_ENERGY' in it and 'FOR PROCESSING' in intent else None)
                uid = str(ds.get('SOPInstanceUID', p))
                if kind is None or uid in seen:
                    continue
                seen.add(uid)
                iu = ds.get((0x0008, 0x3010)); ev = [str(u) for u in (iu.value if iu is not None and iu.VM > 1 else ([iu.value] if iu is not None else []))]
                pv = ds.get((0x0019, 0x1010)) if sysk == 'siemens' else None; m = re.search(r'_W\(([0-9.]+)\)', str(pv.value)) if pv is not None else None
                recs.append(dict(path=p, root=root, kind=kind, system=sysk, ev=ev, study=str(ds.get('StudyInstanceUID', '')), time=str(ds.get('AcquisitionTime', '')),
                                 kv=float(ds.get('KVP') or 0), mAs=float(ds.get('ExposureInuAs') or 0) / 1000, W=float(m.group(1)) if m else None, series=int(ds.get('SeriesNumber') or 0)))
    by = {(r['study'], e): r for r in recs if r['kind'] in ('L', 'H') for e in r['ev']}
    by_time = {(r['study'], r['time'], r['kind']): r for r in recs if r['kind'] in ('L', 'H') and r['time']}
    out = []
    for r in recs:
        if r['kind'] != 'R':
            continue
        L = next((by[(r['study'], e)] for e in r['ev'] if (r['study'], e) in by and by[(r['study'], e)]['kind'] == 'L'), None)
        H = next((by[(r['study'], e)] for e in r['ev'] if (r['study'], e) in by and by[(r['study'], e)]['kind'] == 'H'), None)
        if H and not L:
            L = by_time.get((r['study'], H['time'], 'L'))
        if L and not H:
            H = by_time.get((r['study'], L['time'], 'H'))
        if not L or not H:
            continue
        s = SYSTEMS[r['system']]
        parts = os.path.relpath(r['path'], r['root']).split(os.sep); chain = [os.path.basename(r['root'])] + parts[:-1]
        kvdir = next((c for c in chain if re.match(r'^\d+\s*kV', c, re.I)), None)
        if r['system'] == 'siemens':
            main = kvdir or (parts[0] if len(parts) > 1 else chain[0])
        else:
            main = '%s_%s' % (s['label'], kvdir or '%dkV' % round(L['kv']))
        session = re.sub(r'\s+', '_', parts[-4] if len(parts) >= 4 else chain[0])
        rid = '%s_%04d' % (session, r['series']) if r['system'] == 'siemens' else '%s_%dkV' % (os.path.basename(L['path']), round(L['kv']))
        out.append(dict(L=L['path'], H=H['path'], R=r['path'], system=r['system'], kv=L['kv'], mAsL=L['mAs'], mAsH=H['mAs'], W=r['W'], main=re.sub(r'\s+', '', main),
                        session=session, id=rid))
    return sorted(out, key=lambda t: (t['main'], t['session'], t['id']))


def acquisition(dsL, dsH, dsR):
    """acquisition data from the headers, for the comparison between systems"""
    def num(ds, k, scale=1.0, nd=3):
        v = ds.get(k)
        try:
            return round(float(v) * scale, nd)
        except (TypeError, ValueError):
            return ''
    sw = dsL.get('SoftwareVersions', ''); sw = str(sw[0] if isinstance(sw, pydicom.multival.MultiValue) else sw)
    aec = str(dsL.get('ExposureControlModeDescription', '')); lvl = re.search(r'Dose Level=(\d+)%', aec)
    grid = dsL.get('Grid', ''); grid = '/'.join(grid) if isinstance(grid, pydicom.multival.MultiValue) else str(grid)
    row = dict(model=str(dsL.get('ManufacturerModelName', '')), software=sw, detector=str(dsL.get('DetectorType', '')), pixel_mm=round(float(dsL.ImagerPixelSpacing[0]), 4),
               SID_mm=num(dsL, 'DistanceSourceToDetector', nd=1), magnification=num(dsL, 'EstimatedRadiographicMagnificationFactor', nd=4), grid=grid,
               AEC=('Dose Level=%s%%' % lvl.group(1)) if lvl else aec[:40], thickness_mm=num(dsL, 'BodyPartThickness', nd=1), force_N=num(dsL, 'CompressionForce', nd=1))
    for tag, ds in (('LE', dsL), ('HE', dsH)):
        row.update({'kV_' + tag: num(ds, 'KVP', nd=1), 'target_' + tag: str(ds.get('AnodeTargetMaterial', '')), 'filter_' + tag: str(ds.get('FilterMaterial', '')),
                    'filter_mm_' + tag: num(ds, 'FilterThicknessMinimum'), 'mAs_' + tag: num(ds, 'ExposureInuAs', 1e-3, 1), 'exposure_ms_' + tag: num(ds, 'ExposureTime', nd=0),
                    'tube_mA_' + tag: num(ds, 'XRayTubeCurrent', nd=0), 'AGD_mGy_' + tag: num(ds, 'OrganDose', 100.0, 3), 'ESD_mGy_' + tag: num(ds, 'EntranceDoseInmGy', nd=3),
                    'HVL_mmAl_' + tag: num(ds, 'HalfValueLayer', nd=3)})
    row['AGD_mGy_total'] = round(row['AGD_mGy_LE'] + row['AGD_mGy_HE'], 3) if row['AGD_mGy_LE'] != '' and row['AGD_mGy_HE'] != '' else ''
    row['bits_raw'], row['bits_recombined'] = int(dsL.BitsStored), int(dsR.BitsStored)
    return row


ACQ_KEYS = ['model', 'software', 'kV_LE', 'target_LE', 'filter_LE', 'filter_mm_LE', 'HVL_mmAl_LE', 'mAs_LE', 'exposure_ms_LE', 'tube_mA_LE', 'kV_HE', 'target_HE', 'filter_HE',
            'filter_mm_HE', 'HVL_mmAl_HE', 'mAs_HE', 'exposure_ms_HE', 'tube_mA_HE', 'AGD_mGy_LE', 'AGD_mGy_HE', 'AGD_mGy_total', 'ESD_mGy_LE', 'ESD_mGy_HE', 'AEC',
            'thickness_mm', 'force_N', 'detector', 'pixel_mm', 'SID_mm', 'magnification', 'grid', 'bits_raw', 'bits_recombined']


# ----------------------------------------------------------------------------------------------- spectra and weights
def footprint(H, sat, px):
    """the phantom: HE below 60 % of the open-field level (99.5th percentile of the HE image, saturated pixels left
    out), largest connected region, holes filled"""
    hv = H[H < sat] if (H >= sat).any() else H
    m = ndi.binary_opening(H < 0.6 * np.percentile(hv, 99.5), iterations=max(1, int(round(OPEN_MM / px)))); lab, n = ndi.label(m)
    if n == 0:
        raise RuntimeError('no object found in the HE image')
    return ndi.binary_fill_holes(lab == (np.argmax(ndi.sum(m, lab, range(1, n + 1))) + 1))


class Rois:
    """Hann-windowed, quadratically detrended N px ROIs inside a mask; per ROI and band the sums of the auto and cross
    power of Lam_L, Lam_H and R"""
    def __init__(self, B):
        N = self.N = B.N; self.win = np.outer(np.hanning(N), np.hanning(N)); yv, xv = np.mgrid[0:N, 0:N]
        Ad = np.c_[np.ones(N * N), xv.ravel(), yv.ravel(), (xv * xv).ravel(), (yv * yv).ravel(), (xv * yv).ravel()]
        self.Ad, self.pinv = Ad, np.linalg.pinv(Ad)
        f = np.fft.fftfreq(N, d=B.px); FR = np.hypot(f[:, None], f[None, :])
        self.bands = [np.flatnonzero(((FR >= a) & (FR < b)).ravel()) for a, b in zip(B.edges[:-1], B.edges[1:])]
        self.norm = B.px * B.px / (self.win ** 2).sum() / np.array([len(b) for b in self.bands])     # band sum -> mean PSD (mm2)

    def ft(self, a):
        N = self.N
        return np.fft.fft2(self.win * (a - (self.Ad @ (self.pinv @ a.ravel())).reshape(N, N))).ravel()

    def positions(self, ok):
        N = self.N
        return [(y0, x0) for y0 in range(0, ok.shape[0] - N, N // 2) for x0 in range(0, ok.shape[1] - N, N // 2) if ok[y0:y0 + N, x0:x0 + N].all()]

    def stats(self, lL, lH, R, pos):
        N = self.N; keys = ('LL', 'HH', 'HL', 'RH', 'RL', 'RR'); S = {k: np.zeros((len(pos), len(self.bands))) for k in keys}
        for i, (y0, x0) in enumerate(pos):
            sl = (slice(y0, y0 + N), slice(x0, x0 + N)); l, h, r = self.ft(lL[sl]), self.ft(lH[sl]), self.ft(R[sl])
            prod = dict(LL=np.abs(l) ** 2, HH=np.abs(h) ** 2, HL=(h * np.conj(l)).real, RH=(r * np.conj(h)).real, RL=(r * np.conj(l)).real, RR=np.abs(r) ** 2)
            for k in keys:
                S[k][i] = [prod[k][b].sum() for b in self.bands]
        return S


def select(S, fc):
    """ROIs kept for the weights: those whose residual at the cancelling weight (0.1-1 cy/mm) is not an outlier
    (iodine discs, bubbles, markers give content that LE does not predict); robust, upper tail only"""
    keep = np.ones(len(S['LL']), bool); mid = (fc >= 0.1) & (fc < 1.0)
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
    """per band, from the kept ROIs: PSDs, cancelling weight w*, coherence, the system's gain alpha and effective weight"""
    s = {k: v[keep].sum(0) for k, v in S.items()}; n = int(keep.sum())
    LL, HH, HL = (s[k] * norm / n for k in ('LL', 'HH', 'HL'))
    det = s['HH'] * s['LL'] - s['HL'] ** 2
    alpha = (s['RH'] * s['LL'] - s['RL'] * s['HL']) / det; beta = (s['RL'] * s['HH'] - s['RH'] * s['HL']) / det
    R2 = (alpha * s['RH'] + beta * s['RL']) / s['RR']
    ws = HL / LL; C = HL ** 2 / LL
    noise = lambda w: (HH[-1] + w ** 2 * LL[-1]) * nshape                      # quantum noise in Lam_H - w Lam_L
    E = np.clip((HH - C) - noise(ws), 0, None)                                 # non-noise residual already there at w*
    return dict(n=n, LL=LL, HH=HH, HL=HL, ws=ws, coh=C / HH, C=C, E=E, alpha=alpha, beta=beta, R2=R2, w_eff=-beta / alpha, noise=noise)


def design(m, strength, fc, share=None):
    """left-in share rho(f) and tuned weight w_t(f) for one strength (see the header); share: rho0 given directly"""
    T = np.array(BREAST['T']); lo = fc < FADE[0]
    rho0 = strength * float(np.sqrt(max(np.sum(T[lo] - m['E'][lo]), 0.0) / np.sum(m['C'][lo]))) if share is None else strength * share
    fade = np.clip((np.log10(fc) - np.log10(FADE[0])) / (np.log10(FADE[1]) - np.log10(FADE[0])), 0, 1)
    rho = min(rho0, 0.6) * 0.5 * (1 + np.cos(np.pi * fade)); wt = m['ws'] * (1 - rho)
    c = np.where(m['R2'] >= 0.2, m['alpha'] * wt + m['beta'], 0.0)             # change of the system's image per unit Lam_L
    pd = lambda w: m['HH'] - 2 * w * m['HL'] + w ** 2 * m['LL']
    w_sys = np.where(m['R2'] >= 0.2, m['w_eff'], 0.0)
    return dict(rho0=min(rho0, 0.6), rho=rho, wt=wt, c=c, PD=pd(wt), PD_full=pd(m['ws']), PD_system=pd(w_sys),
                N_tuned=m['noise'](wt), N_full=m['noise'](m['ws']), N_system=m['noise'](w_sys))


class Padded:
    """FFT of a log image (smoothly continued outside the phantom, padded), to apply radial filters given per band"""
    def __init__(self, img, fp, B):
        px = B.px; self.fc = B.fc
        m = ndi.binary_dilation(fp, iterations=int(round(DIL_MM / px))); sig = int(round(SIGMA_MM / px))
        num = ndi.gaussian_filter(np.where(m, img, 0.0), sig); den = ndi.gaussian_filter(m.astype(float), sig)
        f = np.where(m, img, np.where(den > 1e-3, num / np.maximum(den, 1e-3), img[m].mean()))
        self.shape = img.shape; self.pad = int(round(PAD_MM / px))
        ny, nx = (sfft.next_fast_len(n + 2 * self.pad, real=True) for n in img.shape); self.n = (ny, nx)
        fpad = np.pad(f, ((self.pad, ny - img.shape[0] - self.pad), (self.pad, nx - img.shape[1] - self.pad)), mode='reflect')
        self.F = sfft.rfft2(fpad, workers=-1)
        self.FR = np.hypot(sfft.fftfreq(ny, d=px)[:, None], sfft.rfftfreq(nx, d=px)[None, :])

    def filtered(self, cb):
        FC = self.fc; g = lambda f: f ** 2 / (f ** 2 + F0 ** 2)
        filt = np.interp(np.log10(np.clip(self.FR, FC[0], FC[-1])), np.log10(FC), cb)
        low = self.FR < FC[0]; filt[low] = cb[0] * g(self.FR[low]) / g(FC[0])
        p = self.pad; return sfft.irfft2(self.F * filt, s=self.n, workers=-1)[p:p + self.shape[0], p:p + self.shape[1]]


# ----------------------------------------------------------------------------------------------- output
def binned(a, k=3):
    h, w = (a.shape[0] // k) * k, (a.shape[1] // k) * k
    return a[:h, :w].reshape(h // k, k, w // k, k).mean((1, 3))


def save_dicom(ds, img, label, system, path):
    if getattr(ds.file_meta, 'TransferSyntaxUID', None) is not None and ds.file_meta.TransferSyntaxUID.is_compressed:
        ds.file_meta.TransferSyntaxUID = ExplicitVRLittleEndian
    ds.PixelData = np.ascontiguousarray(img, dtype='<u2').tobytes()
    ds.SOPInstanceUID = generate_uid(); ds.file_meta.MediaStorageSOPInstanceUID = ds.SOPInstanceUID; ds.SeriesInstanceUID = generate_uid()
    ds.SeriesDescription = ('%s, L2 %s' % (str(ds.get('SeriesDescription', 'CEM')), label))[:64]
    ds.ImageComments = 'L2 recombination, %s (research; not %s output)' % (label, SYSTEMS[system]['label'])
    ds.save_as(path)


def fig_comparison(R, X, win, med, name, head, path, system):
    """the system's recombined image, ours and the difference (ours minus the system's); top row in the system's default
    window, bottom row in a second fixed window; the difference on the same grey scale, centred on zero"""
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    s = SYSTEMS[system]; narrow = s['narrow']
    rows = np.flatnonzero((R > 0).any(1)); cols = np.flatnonzero((R > 0).any(0)); k = 30
    sl = (slice(max(rows.min() - k, 0), rows.max() + k), slice(max(cols.min() - k, 0), cols.max() + k))
    D = np.where(R > 0, X - R, 0.0)
    fig, ax = plt.subplots(2, 3, figsize=(18, 16.5))
    for i, (c, w) in enumerate(((win[0], win[1]), (med, narrow))):
        for j, (a, c0) in enumerate(((R, c), (X, c), (D, 0.0))):
            ax[i, j].imshow(binned(a[sl]), cmap='gray', vmin=c0 - w / 2, vmax=c0 + w / 2, interpolation='antialiased'); ax[i, j].set_xticks([]); ax[i, j].set_yticks([])
    for j, n in enumerate((s['own'], name, 'difference, ours minus %s (same scale, 0 = mid grey)' % s['label'])):
        ax[0, j].set_title(n, fontsize=12, loc='left')
    ax[0, 0].set_ylabel("%s default window (centre %.0f, width %.0f)" % ((s['poss'],) + tuple(win)), fontsize=10)
    ax[1, 0].set_ylabel(s['narrow_label'] % narrow, fontsize=10)
    fig.suptitle(head, fontsize=11, x=0.01, ha='left')
    fig.tight_layout(rect=(0, 0, 1, 0.965)); fig.savefig(path, dpi=80); plt.close(fig)


def fig_curves(t, m, d, strength, achieved, path, B):
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    FC = B.fc; lab = SYSTEMS[t['system']]['label']; poss = SYSTEMS[t['system']]['poss']; ok = m['R2'] >= 0.2
    fig, ax = plt.subplots(1, 3, figsize=(18, 5.2))
    ax[0].plot(FC[ok], m['w_eff'][ok], 'o-', color='#d62728', label='%s effective weight (measured from R)' % poss)
    ax[0].plot(FC, m['ws'], 'o-', color='k', label='w*(f): cancels the spheres'); ax[0].plot(FC, d['wt'], 'o-', color='#2ca02c', label='tuned w_t(f), strength %g' % strength)
    if achieved is not None:
        ax[0].plot(FC[ok], achieved[ok], 'x', color='#1f77b4', ms=7, label='measured in the tuned image (check)')
    ax[0].set_title('weight on the LE image per spatial frequency', fontsize=10); ax[0].set_ylabel('weight')
    ax[1].plot(FC, 100 * d['rho'], 'o-', color='#2ca02c', label='tuned, strength %g' % strength)
    ax[1].plot(FC[ok], 100 * (1 - m['w_eff'][ok] / m['ws'][ok]), 'o-', color='#d62728', label=lab)
    ax[1].set_title('share of the sphere contrast left in (%)', fontsize=10); ax[1].set_ylabel('%'); ax[1].set_ylim(-5, 105)
    if t['system'] == 'siemens':          # same system as the breasts: total residuals (noise included) compared directly
        ax[2].fill_between(FC, BREAST['PD_LO'], BREAST['PD_HI'], color='#1f77b4', alpha=0.2, lw=0, label='real breasts under Siemens (9 views, lowest-highest)')
        ax[2].plot(FC, BREAST['PD'], color='#1f77b4', label='real breasts, median')
        curves = ((d['PD_system'], '#d62728', 'L2 under Siemens'), (d['PD_full'], 'k', 'L2 full cancellation'), (d['PD'], '#2ca02c', 'L2 tuned, strength %g' % strength))
        ax[2].set_title('spectrum of the log-signal residual Lam_H - w(f) Lam_L (mm$^2$), interior', fontsize=10)
    else:                                 # other system: noise differs, so the structured parts (noise taken off) are compared
        sel = FC < FADE[1]
        ax[2].fill_between(FC[sel], np.array(BREAST['T_LO'])[sel], np.array(BREAST['T_HI'])[sel], color='#1f77b4', alpha=0.2, lw=0, label='real breasts under Siemens, structured (9 views, lowest-highest)')
        ax[2].plot(FC[sel], np.array(BREAST['T'])[sel], color='#1f77b4', label='breast target (Siemens, median; used unchanged)')
        pos = lambda x: np.where(x > 0, x, np.nan)          # where the residual is below the noise estimate nothing is drawn
        curves = ((pos(d['PD_system'] - d['N_system']), '#d62728', 'L2 under %s' % lab), (pos(d['PD_full'] - d['N_full']), 'k', 'L2 full cancellation'),
                  (pos(d['PD'] - d['N_tuned']), '#2ca02c', 'L2 tuned, strength %g' % strength))
        ax[2].set_title('structured log-signal residual (noise estimate taken off, mm$^2$), interior', fontsize=10)
    for y, c, lb in curves:
        ax[2].plot(FC, y, 'o-', ms=3, color=c, label=lb)
    ax[2].set_yscale('log')
    for a in ax:
        a.set_xscale('log'); a.set_xlabel('spatial frequency (cy/mm)'); a.legend(fontsize=8); a.grid(alpha=0.25)
    fig.suptitle('%s   %s   %.0f kV, %.0f / %.0f mAs   (%d regions used)' % (t['id'], lab, t['kv'], t['mAsL'], t['mAsH'], m['n']), fontsize=11, x=0.01, ha='left')
    fig.tight_layout(rect=(0, 0, 1, 0.94)); fig.savefig(path, dpi=85); plt.close(fig)


def write_curves(m, d, achieved, path, B):
    with open(path, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['band_from_cy_mm', 'band_to_cy_mm', 'band_centre_cy_mm', 'system_effective_weight', 'system_gain_alpha', 'system_fit_R2', 'w_star', 'coherence',
                    'left_in_share_rho', 'tuned_weight', 'tuned_weight_measured', 'residual_system_mm2', 'residual_full_cancellation_mm2', 'residual_tuned_mm2',
                    'noise_estimate_tuned_mm2', 'noise_shape_g', 'non_noise_residual_at_w_star_E_mm2', 'cancellable_power_C_mm2', 'breast_target_T_mm2', 'breast_target_lowest_mm2',
                    'breast_target_highest_mm2'])
        for i in range(len(B.fc)):
            w.writerow([B.edges[i], B.edges[i + 1], '%.4f' % B.fc[i], '%.5f' % m['w_eff'][i], '%.1f' % m['alpha'][i], '%.4f' % m['R2'][i], '%.5f' % m['ws'][i], '%.4f' % m['coh'][i],
                        '%.5f' % d['rho'][i], '%.5f' % d['wt'][i], '' if achieved is None else '%.5f' % achieved[i], '%.4e' % d['PD_system'][i], '%.4e' % d['PD_full'][i],
                        '%.4e' % d['PD'][i], '%.4e' % d['N_tuned'][i], '%.3f' % B.g[i], '%.4e' % m['E'][i], '%.4e' % m['C'][i], '%.4e' % BREAST['T'][i], '%.4e' % BREAST['T_LO'][i],
                        '%.4e' % BREAST['T_HI'][i]])


def fig_check_curves(t, m, mo, psd, path, B):
    """curves of the check: effective weight in the system's image and in ours, the share of the system's image that
    ours reproduces per band, and the power spectra of both images and of their difference"""
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    FC = B.fc; lab = SYSTEMS[t['system']]['label']; poss = SYSTEMS[t['system']]['poss']; ok = m['R2'] >= 0.2
    fig, ax = plt.subplots(1, 3, figsize=(18, 5.2))
    ax[0].plot(FC[ok], m['w_eff'][ok], 'o-', color='#d62728', label='%s recombined image' % poss)
    ax[0].plot(FC[ok], mo['w_eff'][ok], 'x', color='#1f77b4', ms=8, label='our recombination at %s weighting' % poss)
    ax[0].set_title('effective weight on the LE image per spatial frequency', fontsize=10); ax[0].set_ylabel('weight')
    ax[1].plot(FC, 100 * (1 - psd['D'] / psd['R']), 'o-', color='#2ca02c', label='share of %s image that ours reproduces' % poss)
    ax[1].set_title('agreement per spatial frequency: 100 x (1 - power of the difference / power of %s image)' % poss, fontsize=10); ax[1].set_ylabel('%'); ax[1].set_ylim(-5, 105)
    for y, c, mk, lb in ((psd['R'], '#d62728', 'o-', '%s recombined image' % poss), (psd['O'], '#1f77b4', 'x-', 'our recombination'), (psd['D'], 'k', 'o-', 'difference')):
        ax[2].plot(FC, y, mk, ms=4, color=c, label=lb)
    ax[2].set_yscale('log'); ax[2].set_title('power spectrum of the images (grey levels$^2$ mm$^2$), interior', fontsize=10)
    for a in ax:
        a.set_xscale('log'); a.set_xlabel('spatial frequency (cy/mm)'); a.legend(fontsize=8); a.grid(alpha=0.25)
    fig.suptitle('%s   %.0f kV, %.0f / %.0f mAs   check: our recombination at %s weighting (%d regions used)' % (t['id'], t['kv'], t['mAsL'], t['mAsH'], poss, m['n']), fontsize=11, x=0.01, ha='left')
    fig.tight_layout(rect=(0, 0, 1, 0.94)); fig.savefig(path, dpi=85); plt.close(fig)


def write_check_curves(m, mo, psd, path, B):
    with open(path, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['band_from_cy_mm', 'band_to_cy_mm', 'band_centre_cy_mm', 'system_effective_weight', 'ours_effective_weight', 'system_gain_alpha', 'ours_gain_alpha',
                    'fit_R2_in_regions', 'power_system', 'power_ours', 'power_difference', 'share_reproduced_percent'])
        for i in range(len(B.fc)):
            w.writerow([B.edges[i], B.edges[i + 1], '%.4f' % B.fc[i], '%.5f' % m['w_eff'][i], '%.5f' % mo['w_eff'][i], '%.1f' % m['alpha'][i], '%.1f' % mo['alpha'][i], '%.4f' % m['R2'][i],
                        '%.4e' % psd['R'][i], '%.4e' % psd['O'][i], '%.4e' % psd['D'][i], '%.1f' % (100 * (1 - psd['D'][i] / psd['R'][i]))])


# ----------------------------------------------------------------------------------------------- one exposure
def check_name(t):
    return '%s_weighting' % t['system']


def out_base(t, s, out_root):
    name = check_name(t) if s == 'check' else 'strength_' + (('%.2f' % s).rstrip('0').rstrip('.') if s != int(s) else '%.1f' % s)
    folder = os.path.join(out_root, '%s_%s' % (t['main'], name), t['session'])
    return folder, os.path.join(folder, t['id'])


def exts_of(t, s):
    return ('_%s.dcm' % check_name(t), '_comparison.png', '_curves.png', '_curves.csv') if s == 'check' else ('_tuned.dcm', '_comparison.png', '_curves.png', '_curves.csv')


def process(job):
    t, strengths, out_root, overwrite, check, share = job; t0 = time.time(); rows = []; lab = SYSTEMS[t['system']]['label']; poss = SYSTEMS[t['system']]['poss']
    if not overwrite:            # what is already there is not made again
        strengths = [s for s in strengths if not all(os.path.exists(out_base(t, s, out_root)[1] + e) for e in exts_of(t, s))]
        check = check and not all(os.path.exists(out_base(t, 'check', out_root)[1] + e) for e in exts_of(t, 'check'))
        if not strengths and not check:
            return rows
    try:
        dsL, dsH, dsR = (pydicom.dcmread(t[k]) for k in ('L', 'H', 'R'))
        L, H, R = (ds.pixel_array.astype(np.float64) for ds in (dsL, dsH, dsR))
        px = float(dsL.ImagerPixelSpacing[0]); vmax = 2 ** int(dsR.BitsStored) - 1; B = Bands(px, t['system'])
        lL, lH = -np.log(np.clip(L, 1, None)), -np.log(np.clip(H, 1, None))
        fp = footprint(H, 2 ** int(dsH.BitsStored) - 1, px); dist = ndi.distance_transform_edt(fp) * px; inner = dist > MARGIN_MM; valid = (R > 0) & (R < vmax)
        rois = Rois(B); pos = rois.positions(inner & valid)
        if len(pos) < 30:
            raise RuntimeError('only %d interior regions found' % len(pos))
        S = rois.stats(lL, lH, R, pos); keep = select(S, B.fc); m = measure(S, keep, rois.norm, B.g)
        le = Padded(lL, fp, B); med = float(np.median(R[inner & valid]))
        wc, ww = dsR.get('WindowCenter', 2600), dsR.get('WindowWidth', 2100)
        win = (float(wc[0] if isinstance(wc, pydicom.multival.MultiValue) else wc), float(ww[0] if isinstance(ww, pydicom.multival.MultiValue) else ww))
        kept = [p for p, k in zip(pos, keep) if k]; allk = np.ones(len(kept), bool); lo = slice(0, 2)
        common = dict(id=t['id'], session=t['session'], system=lab, kV=t['kv'], regions_used=m['n'], regions_left_out=int((~keep).sum()),
                      header_W=t['W'] if t['W'] is not None else '', system_effective_weight=round(float(m['w_eff'][lo].mean()), 4), w_star=round(float(m['ws'][lo].mean()), 4),
                      system_left_in_share_percent=round(100 * float(1 - m['w_eff'][lo].mean() / m['ws'][lo].mean()), 2), interior_median_R=round(med, 1),
                      interior_sd_R=round(float(R[inner & valid].std()), 1), LE_raw_median_interior=round(float(np.median(L[inner & valid])), 1),
                      HE_raw_median_interior=round(float(np.median(H[inner & valid])), 1), error='', **acquisition(dsL, dsH, dsR))
        head = '%s   %s   %.0f kV, %.0f / %.0f mAs\n' % (t['id'], lab, t['kv'], t['mAsL'], t['mAsH'])
        for s in strengths:
            d = design(m, s, B.fc, share)
            Rt = np.where(R > 0, np.clip(np.rint(R - le.filtered(d['c'])), 0, vmax), 0)
            a = measure(rois.stats(lL, lH, Rt, kept), allk, rois.norm, B.g)['w_eff']                    # check on the result
            folder, base = out_base(t, s, out_root); os.makedirs(folder, exist_ok=True)
            save_dicom(pydicom.dcmread(t['R']), Rt, 'tuned s=%g' % s, t['system'], base + '_tuned.dcm')
            fig_comparison(R, Rt, win, med, 'ours: tuned, strength %g' % s, head + 'weight on LE at coarse scales: %s %.3f, cancelling w* %.3f, tuned %.3f '
                           '(%.0f %% of the sphere contrast left; %s leaves %.0f %%)' % (lab, m['w_eff'][lo].mean(), m['ws'][lo].mean(), d['wt'][lo].mean(), 100 * d['rho0'], lab,
                                                                                       100 * (1 - m['w_eff'][lo].mean() / m['ws'][lo].mean())), base + '_comparison.png', t['system'])
            fig_curves(t, m, d, s, a, base + '_curves.png', B); write_curves(m, d, a, base + '_curves.csv', B)
            rows.append(dict(common, folder=os.path.dirname(folder), strength=s, target='left-in share given' if share is not None else 'breast target (Siemens)',
                             tuned_weight=round(float(d['wt'][lo].mean()), 4), tuned_weight_measured=round(float(a[lo].mean()), 4), left_in_share_percent=round(100 * d['rho0'], 2),
                             E_sum_0p06_0p3_mm2='%.3e' % m['E'][B.fc < FADE[0]].sum(), C_sum_0p06_0p3_mm2='%.3e' % m['C'][B.fc < FADE[0]].sum(),
                             residual_0p2_system_mm2='%.3e' % d['PD_system'][3], residual_0p2_tuned_mm2='%.3e' % d['PD'][3], noise_0p2_tuned_mm2='%.3e' % d['N_tuned'][3],
                             breast_target_0p2_mm2='%.3e' % BREAST['T'][3], seconds=round(time.time() - t0, 1)))
        if check:
            # our recombination of the raw LE and HE images at the system's weighting: gain alpha(f) and weight w_eff(f) per band as
            # measured (R = alpha Lam_H + beta Lam_L, beta = -alpha w_eff), applied to the whole images; offset from the interior
            ok = m['R2'] >= 0.2; rec = Padded(lH, fp, B).filtered(np.where(ok, m['alpha'], 0.0)) + le.filtered(np.where(ok, m['beta'], 0.0))
            ins = valid & (dist > 15); edge = valid & fp & (dist <= 15)
            Ro = np.where(R > 0, np.clip(np.rint(rec + np.median(R[ins] - rec[ins])), 0, vmax), 0)
            mo = measure(rois.stats(lL, lH, Ro, kept), allk, rois.norm, B.g); n = len(kept)
            psd = dict(R=S['RR'][keep].sum(0) * rois.norm / n, O=rois.stats(lL, lH, Ro, kept)['RR'].sum(0) * rois.norm / n,
                       D=rois.stats(lL, lH, Ro - R, kept)['RR'].sum(0) * rois.norm / n)
            st = dict(corr=float(np.corrcoef(R[ins], Ro[ins])[0, 1]), rmse=float(np.sqrt(np.mean((Ro[ins] - R[ins]) ** 2))), sd=float(R[ins].std()),
                      edge=float(np.corrcoef(R[edge], Ro[edge])[0, 1]), share=float(100 * (1 - psd['D'][:8].sum() / psd['R'][:8].sum())))
            folder, base = out_base(t, 'check', out_root); os.makedirs(folder, exist_ok=True)
            save_dicom(pydicom.dcmread(t['R']), Ro, 'our recombination at %s weighting' % poss, t['system'], base + exts_of(t, 'check')[0])
            fig_comparison(R, Ro, win, med, 'ours: recombined from raw LE and HE at %s weighting' % poss, head + 'interior (more than 15 mm inside the outline): correlation %.3f, '
                           'difference %.0f grey levels rms (image SD %.0f); outer 15 mm: correlation %.2f; %.0f %% of %s image reproduced at 0.06-1 cy/mm in the regions used'
                           % (st['corr'], st['rmse'], st['sd'], st['edge'], st['share'], poss), base + '_comparison.png', t['system'])
            fig_check_curves(t, m, mo, psd, base + '_curves.png', B); write_check_curves(m, mo, psd, base + '_curves.csv', B)
            rows.append(dict(common, folder=os.path.dirname(folder), strength='%s weighting (check)' % poss, check_interior_correlation=round(st['corr'], 4),
                             check_interior_difference_rms=round(st['rmse'], 1), check_interior_image_sd=round(st['sd'], 1), check_edge_15mm_correlation=round(st['edge'], 3),
                             check_share_reproduced_percent=round(st['share'], 1), check_ours_effective_weight=round(float(mo['w_eff'][lo].mean()), 4), seconds=round(time.time() - t0, 1)))
    except Exception as ex:
        rows.append(dict(folder=os.path.join(out_root, '%s_strength_error' % t['main']), id=t['id'], session=t['session'], system=lab, kV=t['kv'], strength='', error='%s: %s' % (type(ex).__name__, ex)))
        traceback.print_exc()
    return rows


EXCEL = 'L2_tuned_recombination_multisystem_numbers.xlsx'


def write_excel(out_root):
    """all numbers of all analyses in one Excel file in the output folder, rebuilt from the CSV files of this and earlier
    runs: one line per image (tuned images; check at the system's weighting) and one line per image and band (curves)"""
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

    sheets = {'tuned_images': [], 'own_weighting_check': [], 'tuned_curves': [], 'check_curves': [], 'not_processed': []}
    chk = re.compile(r'^(.*)_(%s)_weighting$' % '|'.join(SYSTEMS))
    for name in sorted(os.listdir(out_root)):
        folder = os.path.join(out_root, name)
        if not os.path.isdir(folder) or not os.path.exists(os.path.join(folder, 'summary.csv')):
            continue
        if name.endswith('_strength_error'):
            sheets['not_processed'] += read(os.path.join(folder, 'summary.csv'), dict(folder=name)); continue
        c = chk.match(name); m = re.match(r'^(.*)_strength_([0-9.]+)$', name)
        if not c and not m:
            continue
        main = c.group(1) if c else m.group(1)
        sheets['own_weighting_check' if c else 'tuned_images'] += read(os.path.join(folder, 'summary.csv'), dict(main_folder=main, output_folder=name))
        for session in sorted(os.listdir(folder)):
            d = os.path.join(folder, session)
            if os.path.isdir(d):
                for f in sorted(os.listdir(d)):
                    if f.endswith('_curves.csv') and not f.startswith('.'):
                        extra = dict(main_folder=main, output_folder=name, session=session, id=f[:-len('_curves.csv')])
                        if not c:
                            extra['strength'] = float(m.group(2))
                        sheets['check_curves' if c else 'tuned_curves'] += read(os.path.join(d, f), extra)
    notes = [('sheet', 'content'),
             ('tuned_images', 'one line per image and strength: acquisition data from the headers; weights at coarse scales (mean of 0.06-0.15 cy/mm); share of the sphere contrast left in; residual power at 0.2-0.3 cy/mm'),
             ('own_weighting_check', "one line per image: our recombination of raw LE and HE at the system's own weighting against the system's recombined image (option --system-weighting)"),
             ('tuned_curves', 'one line per image, strength and spatial-frequency band: the curves of <id>_curves.png'),
             ('check_curves', 'one line per image and band: the curves of the check'),
             ('not_processed', 'images that could not be processed, with the error'),
             ('', ''), ('quantity', 'meaning'),
             ('system_effective_weight', "weight on the LE image that the system's recombined image contains, measured from it (R = alpha (Lam_H - w_eff Lam_L))"),
             ('header_W', 'weight written in the header of the recombined image (Siemens only; Hologic writes none)'),
             ('w_star', 'weight that cancels the sphere background'), ('tuned_weight / tuned_weight_measured', 'weight of the tuned image as designed / as measured in the tuned image'),
             ('left_in_share / system_left_in_share', 'share of the sphere contrast left in the image, 1 - weight / w_star (tuned / the system itself)'),
             ('target', "strength 1 is set by the breast target measured on Siemens; for Hologic it is used unchanged (assumption, no raw Hologic breast images)"),
             ('residual ... mm2', 'power spectral density of the log-signal residual Lam_H - w Lam_L in the interior (Lam = -ln of the stored value); noise included'),
             ('noise_0p2_tuned_mm2', 'estimated quantum-noise part of residual_0p2_tuned_mm2 (last band times the noise shape g)'),
             ('E_sum / C_sum', 'non-noise residual at w* and cancellable power, summed over 0.06-0.3 cy/mm (they set the left-in share)'),
             ('AGD_mGy', 'Organ Dose (0040,0316) of the header, converted from dGy; the manufacturer\'s estimate'),
             ('ESD_mGy, HVL_mmAl', 'Entrance Dose in mGy (0040,8302) and Half Value Layer (0040,0314) as in the header (Siemens writes no HVL)'),
             ('check_interior_correlation / difference_rms / image_sd', "pixel-wise agreement more than 15 mm inside the outline, in grey levels of the recombined image"),
             ('check_share_reproduced_percent', "100 x (1 - power of the difference / power of the system's image) at 0.06-1 cy/mm in the regions used"),
             ('check_edge_15mm_correlation', "agreement in the outer 15 mm, where the systems apply processing that is not part of our method")]
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


KEYS = ['id', 'session', 'system', 'kV', 'mAs_LE', 'mAs_HE', 'strength', 'target', 'regions_used', 'regions_left_out', 'header_W', 'system_effective_weight', 'w_star', 'tuned_weight',
        'tuned_weight_measured', 'left_in_share_percent', 'system_left_in_share_percent', 'E_sum_0p06_0p3_mm2', 'C_sum_0p06_0p3_mm2', 'residual_0p2_system_mm2', 'residual_0p2_tuned_mm2',
        'noise_0p2_tuned_mm2', 'breast_target_0p2_mm2', 'interior_median_R', 'interior_sd_R', 'LE_raw_median_interior', 'HE_raw_median_interior', 'seconds', 'error'] + ACQ_KEYS
KEYS_CHECK = ['id', 'session', 'system', 'kV', 'mAs_LE', 'mAs_HE', 'strength', 'regions_used', 'regions_left_out', 'header_W', 'system_effective_weight', 'check_ours_effective_weight', 'w_star',
              'system_left_in_share_percent', 'check_interior_correlation', 'check_interior_difference_rms', 'check_interior_image_sd', 'check_share_reproduced_percent',
              'check_edge_15mm_correlation', 'interior_median_R', 'interior_sd_R', 'seconds', 'error'] + ACQ_KEYS


def main():
    ap = argparse.ArgumentParser(description='L2 phantom: tuned recombination of Siemens and Hologic CEM images (see the header of this file).')
    ap.add_argument('inputs', nargs='*', default=['.'], help='folders with the DICOM files (default: the current folder)')
    ap.add_argument('--strengths', type=float, nargs='+', default=[1.0], help='strengths to produce (default 1.0); 0 = full cancellation')
    ap.add_argument('--share', type=float, default=None, help='set the left-in share at strength 1 directly (e.g. 0.068) instead of matching the breast target')
    ap.add_argument('--system-weighting', '--siemens-weighting', dest='check', action='store_true',
                    help="also recombine the raw LE and HE images with our method at the system's own weighting and compare with the system's image (check of the method)")
    ap.add_argument('--out', default='.', help='where the output folders are written (default: the current folder)')
    ap.add_argument('--ext', nargs='+', default=['.ima', '.dcm', ''], help="file endings of the DICOM files (default: .IMA, .dcm and none; '' = no ending)")
    ap.add_argument('--only', default=None, help='process only exposures whose path contains this text')
    ap.add_argument('--skip', nargs='*', default=['4AFC'], help='folders whose name contains one of these texts are not searched (default: 4AFC)')
    ap.add_argument('--limit', type=int, default=None, help='process at most this many exposures per session (for a first look)')
    ap.add_argument('--overwrite', action='store_true', help='make the files again even if they exist (default: existing results are kept, so a run can be continued)')
    a = ap.parse_args()
    T = scan(a.inputs, tuple(e.lower() for e in a.ext), a.skip)
    if a.only:
        T = [t for t in T if a.only.lower() in t['R'].lower()]
    if a.limit:
        seen = {}; T = [t for t in T if seen.setdefault(t['session'], []).append(1) or len(seen[t['session']]) <= a.limit]
    if not T:
        sys.exit('no LE / HE / recombined triplets found below: %s' % ', '.join(a.inputs))
    print('%d exposures (%s) in %d sessions; strengths %s' % (len(T), ', '.join('%s %d' % (SYSTEMS[k]['label'], sum(t['system'] == k for t in T)) for k in SYSTEMS if any(t['system'] == k for t in T)),
                                                          len(set(t['session'] for t in T)), a.strengths), flush=True)
    jobs = [(t, a.strengths, os.path.abspath(a.out), a.overwrite, a.check, a.share) for t in T]; rows = []; t_start = time.time()
    for k, j in enumerate(jobs):
        r = process(j); rows += r; report(k, len(jobs), r)
    for folder in sorted(set(r['folder'] for r in rows)):
        os.makedirs(folder, exist_ok=True)
        f = os.path.join(folder, 'summary.csv'); new = {r['id']: r for r in rows if r['folder'] == folder}
        old = {r['id']: r for r in csv.DictReader(open(f))} if os.path.exists(f) else {}          # lines of earlier runs are kept
        with open(f, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=KEYS_CHECK if folder.endswith('_weighting') else KEYS, extrasaction='ignore'); w.writeheader(); w.writerows(v for _, v in sorted({**old, **new}.items()))
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
            msg = "check at the system's weighting: interior correlation %.3f, difference %.0f rms (image SD %.0f), %.0f %% reproduced, %.0f s" % (
                x['check_interior_correlation'], x['check_interior_difference_rms'], x['check_interior_image_sd'], x['check_share_reproduced_percent'], x['seconds'])
        else:
            msg = 'strength %g: w* %.3f, %s %.3f, tuned %.3f (measured %.3f), %.1f %% left in (%s %.1f %%), %s regions, %.0f s' % (
                x['strength'], x['w_star'], x['system'], x['system_effective_weight'], x['tuned_weight'], x['tuned_weight_measured'], x['left_in_share_percent'], x['system'],
                x['system_left_in_share_percent'], x['regions_used'], x['seconds'])
        print('[%d/%d] %s  %s' % (k + 1, n, x['id'], msg), flush=True)


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""
Measurement-like LE and HE images: the simulated detected energy turned into what the Revelation stores in its raw
(TICEM_RAW) images, so that simulation and measurement can be compared in the same domain (compare_0035.py).

Per beam:
  1. detected energy per pixel  E = N_sr (P + S)   keV
       P  primary from projector.py (keV/sr per pixel), S scatter from postprocess.py; N_sr = photons emitted per
       steradian for the measured exposure (SpekPy fluence per mAs at 1 m x 1e4 cm2/sr x mAs, the exact mAs of that
       pair's LE and HE image from their DICOM headers)
  2. detector blur (DETECTOR_RESPONSE["mtf"]): the presampling MTF published for the Revelation's a-Se detector
       (NHSBSP 2019, [R1] in config.py; it includes the pixel aperture, as the projector samples pixel centres), or
       the ideal pixel aperture only ("aperture"); optional extra Gaussian blur [A]
  3. flat-field normalisation, as the raw images show a flat open field: D = D0 E / E0, E0 = open-field energy
       D0 (HE) = HE open-field stored value of that measured pair; D0 (LE) from the same detector gain
       (gain_mode "constant") [A]
  4. quantum noise: variance N_sr (P2 + S P2/P), P2 = second moment from projector.py (energy-integrating detector;
       scatter assumed to have the primary's mean detected energy [A]), times a factor k. It is generated before the
       detector blur, so it is shaped by the aliased MTF, NPS(f) ~ sum_k MTF(f + k/a)^2 (white for the ideal
       aperture). Plus fixed-pattern (gain) noise of relative rms f and electronic noise of rms e, both white, and a
       multiplicative low-frequency structured noise with the spectrum measured in the open field (shared by LE and
       HE). All are calibrated on the measured HE open fields at three exposure levels (--calibrate-noise).
  4b. offset: the LE stored-value offset measured in REPORT §13.4 (40 units at 28 kV); none at HE
  5. stored values: rounded, clipped at 0 and at the measured ceiling (7500)

One simulated pair is made for ONE measured pair of the reference set (--set, --pair; July default HIGHMAS_28KV2_00),
with that pair's exact LE and HE mAs and HE open level (measured_set.exposures); the spectra are those fitted for that
X-ray unit. Another pair is a separate run with its own --pair (fixed noise seed per pair); nothing is averaged.
Output: output/measlike/<case>/<set>/<pair>/LE.npy, HE.npy (uint16, rows = y along the chest wall, columns = x from
the chest wall, as the projector grid), meta.json, and preview.png in <set>/; <case> = config.case_name()
Usage: python measlike.py [--phantom L2v2|L2] [--background PMMA_paraffin] [--set 0035|july28high]
                          [--scatter-le LE_PMMA_paraffin --scatter-he HE_PMMA_paraffin]
                          [--bin 1] [--seed 1] [--pair ID] [--no-noise] [--le-open STORED]
       python measlike.py --phantom L2 --scatter-le ... --scatter-he ... --calibrate-noise
                          noise model k, f, e from the measured HE open fields of the set's exposures at low / normal /
                          high mAs (0.6-5.9 cy/mm) -> data/detector_noise_<set>.json (no images written)
"""
import argparse, json
import numpy as np
from config import (HERE, OUT_DIR, MEASURED, MEASURED_SETS, DETECTOR_RESPONSE, PHANTOM, BEAMS, spectra_dir,
                    add_case_args, apply_case_args)
from projector import pixel_grid, out_dir as projector_dir
from measured_set import exposures


def n_per_sr(beam, mAs):
    spec = np.load(spectra_dir() / f"{beam}.npz")
    return float(spec["fluence_per_mAs_1m"]) * 1e4 * mAs


def mtf_1d(f, p):
    """detector MTF along one axis. 'nhsbsp2019': the published table up to 10 cy/mm; beyond it a fitted
    c sinc(f p) exp(-2 pi^2 s^2 f^2), scaled to join the table at 10 cy/mm (needed only for the aliasing sum).
    'aperture': sinc(f p)"""
    f = np.abs(f)
    if DETECTOR_RESPONSE["mtf"] == "aperture":
        return np.sinc(f * p)
    t = DETECTOR_RESPONSE["mtf_table"]; ft, mt = np.array(t["f_cy_mm"], float), np.array(t["mtf"], float)
    k = ft > 0
    q = np.polyfit(ft[k] ** 2, np.log(mt[k] / np.sinc(ft[k] * p)), 1)          # ln c - 2 pi^2 s^2 f^2
    model = lambda v: np.exp(np.polyval(q, v ** 2)) * np.sinc(v * p)
    return np.where(f <= ft[-1], np.interp(f, ft, mt), model(f) * mt[-1] / model(ft[-1]))


def detector_filters(shape, p, extra_sigma_mm=0.0):
    """signal filter (MTF, separable in x and y [A]) and noise filter sqrt(A(fx) A(fy)), A(f) = sum_k MTF(f + k/p)^2:
    white noise generated before the detector blur and sampled at pixel pitch p (aliasing). For the ideal aperture
    A = 1 (exactly), i.e. white noise as before."""
    fy = np.fft.fftfreq(shape[0], p)[:, None]
    fx = np.fft.rfftfreq(shape[1], p)[None, :]
    Hs = mtf_1d(fx, p) * mtf_1d(fy, p)
    if extra_sigma_mm > 0:
        Hs = Hs * np.exp(-2 * np.pi ** 2 * extra_sigma_mm ** 2 * (fx ** 2 + fy ** 2))
    if DETECTOR_RESPONSE["mtf"] == "aperture":
        return Hs, np.ones_like(Hs)                    # sum_k sinc^2(f p + k) = 1: white noise
    alias = lambda f: sum(mtf_1d(f + k / p, p) ** 2 for k in range(-4, 5))
    return Hs, np.sqrt(alias(fx) * alias(fy))


def blur(img, p, extra_sigma_mm):
    """detector MTF (detector_filters), in the Fourier domain (reflective padding against wrap-around at the borders)"""
    pad = 16
    a = np.pad(img, pad, mode="reflect")
    Hs, _ = detector_filters(a.shape, p, extra_sigma_mm)
    return np.fft.irfft2(np.fft.rfft2(a) * Hs, s=a.shape)[pad:-pad, pad:-pad]


def shaped_noise(shape, p, rng):
    """unit-variance-per-white-pixel Gaussian noise, shaped by the detector's aliased MTF (detector_filters)"""
    _, Hn = detector_filters(shape, p)
    return np.fft.irfft2(np.fft.rfft2(rng.standard_normal(shape)) * Hn, s=shape)


def nnps(img, mask, p, N=256, bands=(0.3, 0.6, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 5.9)):
    """normalised noise power spectrum (mm2) of a flat region: N px ROIs inside mask, quadratic detrend, divided by the
    ROI mean squared; radial averages in the frequency bands (cy/mm)"""
    yv, xv = np.mgrid[0:N, 0:N]
    Ad = np.c_[np.ones(N * N), xv.ravel(), yv.ravel(), (xv * xv).ravel(), (yv * yv).ravel(), (xv * yv).ravel()]
    pinv = np.linalg.pinv(Ad); P, n = 0, 0
    for y0 in range(0, img.shape[0] - N, N):
        for x0 in range(0, img.shape[1] - N, N):
            sl = (slice(y0, y0 + N), slice(x0, x0 + N))
            if not mask[sl].all():
                continue
            a = img[sl]; d = a - (Ad @ (pinv @ a.ravel())).reshape(N, N)
            P = P + np.abs(np.fft.fft2(d)) ** 2 / a.mean() ** 2; n += 1
    if n == 0:
        return 0, None
    P = P / n * p * p / (N * N)
    f = np.fft.fftfreq(N, p); FR = np.hypot(*np.meshgrid(f, f))
    return n, [float(P[(FR >= lo) & (FR < hi)].mean()) for lo, hi in zip(bands[:-1], bands[1:])]


def open_field_mask(img, p, ceiling):
    """open field beyond the phantom: below the ceiling, within 10 % of the open level, x > 160 mm"""
    cols = np.arange(img.shape[1])[None, :] * p
    return (img < ceiling - 100) & (img > 0.9 * np.median(img[:, -int(30 / p):])) & (cols > 160)


def noise_params():
    """k, f, e of the noise model: fixed in config, or the calibration of the selected measured set"""
    m = DETECTOR_RESPONSE["noise_model"]
    if m != "calibrated":
        return dict(k=float(m["k"]), f=float(m["f"]), e=float(m["e"]), source="config")
    path = HERE / "data" / ("detector_noise_%s.json" % MEASURED["set"])
    c = json.load(open(path)) if path.exists() else {}
    if "k" not in c:
        print("WARNING: no noise-model calibration in %s (run measlike.py --calibrate-noise): k 1, f 0, e 0" % path.name)
        return dict(k=1.0, f=0.0, e=0.0, source="not calibrated")
    if c.get("mtf") != DETECTOR_RESPONSE["mtf"]:
        print("WARNING: the noise calibration was made with mtf '%s', now '%s': recalibrate" % (c.get("mtf"), DETECTOR_RESPONSE["mtf"]))
    return dict(k=float(c["k"]), f=float(c["f"]), e=float(c["e"]), source="calibrated on %s (%d open fields)" % (c["set"], len(c["exposures"])))


def structured_filter(shape, p):
    """amplitude filter for the structured noise: sqrt(excess NNPS) / p (unit-variance white noise has NPS p^2),
    log-interpolated between the calibration bands, constant below the lowest, zero above the highest band"""
    if not DETECTOR_RESPONSE.get("structured_noise"):
        return None
    path = HERE / "data" / ("detector_noise_%s.json" % MEASURED["set"])
    c = json.load(open(path)) if path.exists() else {}
    if "structured" not in c:
        print("WARNING: no structured-noise calibration in %s (run measlike.py --calibrate-noise): left out" % path.name)
        return None
    b = np.array(c["structured"]["bands_cy_mm"]); v = np.clip(np.array(c["structured"]["excess_nnps_mm2"]), 1e-14, None)
    fc = np.sqrt(b[:-1] * b[1:])
    FR = np.hypot(np.fft.fftfreq(shape[0], p)[:, None], np.fft.rfftfreq(shape[1], p)[None, :])
    lv = np.interp(np.log(np.clip(FR, fc[0], fc[-1])), np.log(fc), np.log(v))
    return np.where(FR <= b[-1], np.sqrt(np.exp(lv)) / p, 0.0)


def offset(beam):
    o = DETECTOR_RESPONSE["offset_stored"][beam]
    return float(o[int(round(BEAMS[beam]["kvp"]))] if isinstance(o, dict) else o)


def load_beam(beam, tag, scatter_name, image_bin):
    pr = projector_dir()
    P = np.load(pr / f"primary_{beam}{tag}.npy").astype(np.float64)
    P2 = np.load(pr / f"primary2_{beam}{tag}.npy").astype(np.float64)
    P0 = np.load(pr / f"primary_{beam}_open{'_bin%d' % image_bin if image_bin > 1 else ''}.npy").astype(np.float64)
    if scatter_name:
        S = np.load(OUT_DIR / "mc" / f"scatter_{scatter_name}{'_bin%d' % image_bin if image_bin > 1 else ''}.npy")
        S = np.clip(S.astype(np.float64), 0, None)
    else:
        print("WARNING: no scatter image for %s: primary only" % beam)
        S = np.zeros_like(P)
    return P, P2, P0, S


def main():
    ap = argparse.ArgumentParser()
    add_case_args(ap)
    ap.add_argument("--scatter-le", default=None); ap.add_argument("--scatter-he", default=None)
    ap.add_argument("--bin", type=int, default=1, help="pixel binning of the projector images used")
    ap.add_argument("--seed", type=int, default=1, help="noise seed of the set's first pair (pair n: seed + n)")
    ap.add_argument("--no-noise", action="store_true")
    ap.add_argument("--le-open", type=float, default=None,
                    help="LE open-field stored value; default: same detector gain as HE (gain_mode 'constant')")
    ap.add_argument("--calibrate-noise", action="store_true",
                    help="fit the noise-variance factor to the measured HE open field of the pair (0.6-2 cy/mm) and stop")
    a = ap.parse_args()
    case = apply_case_args(a)                               # case name: background (L2 v2) or L2_<background>
    exps = exposures()                                      # the selected pair: exact LE / HE mAs from its headers
    seed_of = {e["id"]: a.seed + n for n, e in enumerate(exposures(all_pairs=True))}   # fixed noise seed per pair
    tag = f"_bin{a.bin}" if a.bin > 1 else ""
    x, y, p = pixel_grid(a.bin)
    det = DETECTOR_RESPONSE
    base = {b: load_beam(b, tag, sc, a.bin) for b, sc in (("LE", a.scatter_le), ("HE", a.scatter_he))}
    # open-field stored levels: HE measured (per pair); LE from the same gain (energy ratio at the image centre)
    ic, jc = len(y) // 2, len(x) // 4
    e0c = {b: float(np.median(base[b][2][ic - 20:ic + 20, jc - 20:jc + 20])) for b in base}     # keV/sr per pixel
    root = OUT_DIR / "measlike" / case / MEASURED["set"]
    if det["mtf"] != "aperture" and a.bin > 1:
        print("NOTE: --bin %d: the published MTF is for 85 um pixels; using the ideal aperture of the binned pixel" % a.bin)
        det = dict(det, mtf="aperture")
        DETECTOR_RESPONSE.update(mtf="aperture")
    npar = dict(k=1.0, f=0.0, e=0.0, source="calibration run (quantum noise only)") if a.calibrate_noise else noise_params()
    for k, ex in enumerate(exps):
        rng = np.random.default_rng(seed_of[ex["id"]])
        data = {}
        meta = dict(case=case, phantom=PHANTOM["model"], background=PHANTOM["background"], measured_set=MEASURED["set"],
                    pair=ex["id"], mAs=ex["mAs"], spectra=str(spectra_dir()), sphere_material=PHANTOM["sphere_material"],
                    liquid_material=PHANTOM["liquid_material"], pixel_mm=p, detector=det, beams={},
                    noise_model=npar, offset_stored={b: offset(b) for b in ("LE", "HE")})
        for beam, (P, P2, P0, S) in base.items():
            nsr = n_per_sr(beam, ex["mAs"][beam])
            data[beam] = dict(E=nsr * (P + S), V=nsr * P2 * (1 + S / np.maximum(P, 1e-30)), E0=nsr * P0, nsr=nsr)
        D0 = {"HE": ex["he_open_stored"]}
        D0["LE"] = a.le_open if a.le_open else D0["HE"] * data["LE"]["nsr"] * e0c["LE"] / (data["HE"]["nsr"] * e0c["HE"])
        fp_shared = rng.standard_normal(base["HE"][0].shape) if det["fixed_pattern_shared"] else None
        Hst = None if (a.calibrate_noise or a.no_noise) else structured_filter(base["HE"][0].shape, p)
        st_field = lambda: np.fft.irfft2(np.fft.rfft2(rng.standard_normal(base["HE"][0].shape)) * Hst, s=base["HE"][0].shape)
        st_shared = st_field() if (Hst is not None and det.get("structured_shared", True)) else None
        for beam, d in data.items():
            if a.calibrate_noise and beam != "HE":
                continue
            Eb = blur(d["E"], p, det["extra_blur_sigma_mm"])
            D = D0[beam] * Eb / d["E0"]                    # signal part of the stored value
            rel_noise = np.sqrt(npar["k"] * d["V"]) / np.maximum(d["E"], 1e-30)
            if not a.no_noise:
                Dn = D * rel_noise * shaped_noise(D.shape, p, rng)                       # quantum, MTF-shaped
                if npar["f"] > 0:                                                        # fixed pattern, white
                    Dn = Dn + npar["f"] * D * (fp_shared if fp_shared is not None else rng.standard_normal(D.shape))
                if npar["e"] > 0:                                                        # electronic, white
                    Dn = Dn + npar["e"] * rng.standard_normal(D.shape)
                if Hst is not None:                                                      # structured, low frequency
                    Dn = Dn + D * (st_shared if st_shared is not None else st_field())
                D = D + Dn
            D = D + (0.0 if a.calibrate_noise else offset(beam))
            stored = np.clip(np.rint(D), 0, MEASURED["ceiling"]).astype(np.uint16)
            data[beam]["stored"] = stored
            if a.calibrate_noise:
                calibrate(stored, ex, p)
                return
            inside = d["E"] > 0.05 * np.median(d["E0"])
            meta["beams"][beam] = dict(
                mAs=ex["mAs"][beam], photons_per_sr=d["nsr"], open_field_stored=D0[beam],
                open_field_keV_per_pixel_centre=e0c[beam] * d["nsr"],
                effective_detected_photons_per_pixel=float(np.median((d["E"] ** 2 / np.maximum(d["V"], 1e-30))[inside])),
                median_relative_noise=float(np.median(rel_noise[inside])),
                scatter=a.scatter_le if beam == "LE" else a.scatter_he,
                clipped_fraction=float((stored >= MEASURED["ceiling"]).mean()))
        out = root / ex["id"]; out.mkdir(parents=True, exist_ok=True)
        for beam in data:
            np.save(out / f"{beam}{tag}.npy", data[beam]["stored"])
        meta.update(x_mm=[float(x[0]), float(x[-1])], y_mm=[float(y[0]), float(y[-1])], noise=not a.no_noise, seed=seed_of[ex["id"]],
                    d0_note="HE open level = measured pair %s; LE " % ex["id"] +
                            ("given" if a.le_open else "from constant detector gain (open-field energy ratio)"))
        json.dump(meta, open(out / f"meta{tag}.json", "w"), indent=1)
        if k == 0:
            preview(data, x, y, root / f"preview{tag}.png", case)
        print("%s: " % ex["id"] + "; ".join(
            "%s %.2f mAs, open %.0f, %.0f photons/pixel, noise %.4f, %.1f %% at the ceiling"
            % (b, m_["mAs"], m_["open_field_stored"], m_["effective_detected_photons_per_pixel"], m_["median_relative_noise"],
               100 * m_["clipped_fraction"]) for b, m_ in meta["beams"].items()), flush=True)
    print("written", root)


def calibrate(H_sim, ex, p, bands=(0.6, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 5.9)):
    """fit the noise model on the measured HE open fields of the set's exposures at low / normal / high mAs (same kV and
    unit): NNPS_rel(f, S) = a Q(f) / S + B + C / S^2, Q = band mean of the aliased-MTF noise shape (quantum), B white
    fixed-pattern, C white electronic (S = open-field stored value; non-negative least squares on relative residuals).
    Then k = a Q / S over the simulated quantum NNPS (k = 1) at the pair's open level (0.6-2 cy/mm), f = sqrt(B) / p,
    e = sqrt(C) / p. Written to data/detector_noise_<set>.json"""
    import pydicom
    from scipy.optimize import nnls
    import measured_set as MS
    src = MEASURED_SETS[MEASURED["set"]]["source"]
    if src.get("kind") != "list":
        raise SystemExit("the noise-model fit needs a set with an exposure list at several mAs levels")
    E = [e for e in json.load(open(src["exposures"])) if e["kv"] == src["kvp"] and e["unit"] == MEASURED_SETS[MEASURED["set"]]["unit"]
         and e["level"] in ("low", "normal", "high")]
    rows = []
    for e in E:
        path = MS.locate(e["H"])
        try:
            H = pydicom.dcmread(path).pixel_array.astype(float)
        except FileNotFoundError:
            print("  missing %s (%s): skipped" % (e["id"], path)); continue
        if MEASURED_SETS[MEASURED["set"]]["mirror"]:
            H = H[:, ::-1]
        m = open_field_mask(H, p, MEASURED["ceiling"])
        n, v = nnps(H, m, p, bands=bands)
        if not n:
            print("  %s: no open-field ROIs" % e["id"]); continue
        rows.append(dict(id=e["id"], level=e["level"], mAsH=e["mAsH"], S=float(np.median(H[m])), rois=n, nnps_mm2=v))
        print("  %-18s %-6s %5.1f mAs  S %6.0f  NNPS(1e-7 mm2) %s" % (e["id"], e["level"], e["mAsH"], rows[-1]["S"], " ".join("%.2f" % (x * 1e7) for x in v)), flush=True)
    levels = sorted({r["level"] for r in rows})
    if len(levels) < 2:
        raise SystemExit("open fields at only %s mAs: connect the drive with the low / normal-mAs exposures" % levels)
    # quantum noise shape per band: radial band means of A(fx) A(fy) (detector_filters), normalised at f -> 0
    N = 256; _, Hn = detector_filters((N, N), p)
    FR = np.hypot(np.fft.fftfreq(N, p)[:, None], np.fft.rfftfreq(N, p)[None, :])
    Q = np.array([(Hn ** 2)[(FR >= lo) & (FR < hi)].mean() for lo, hi in zip(bands[:-1], bands[1:])])
    terms = [t for t in ("quantum", "fixed_pattern", "electronic") if t in DETECTOR_RESPONSE["noise_fit_terms"]]
    A, y = [], []
    for r in rows:
        for b, v in enumerate(r["nnps_mm2"]):
            col = dict(quantum=Q[b] / r["S"] / v, fixed_pattern=1.0 / v, electronic=1.0 / r["S"] ** 2 / v)
            A.append([col[t] for t in terms]); y.append(1.0)                                      # relative residuals
    sol, _ = nnls(np.array(A), np.array(y)); sol = dict(zip(terms, sol))
    a_, B, C = sol.get("quantum", 0.0), sol.get("fixed_pattern", 0.0), sol.get("electronic", 0.0)
    pred = lambda S: a_ * Q / S + B + C / S ** 2
    res = np.array([[v / pv - 1 for v, pv in zip(r["nnps_mm2"], pred(r["S"]))] for r in rows])
    # simulated quantum noise (k = 1) at the pair's open level
    pr = next(MS.pairs(MEASURED["set"], ids=[ex["id"]]), None)
    S_pair = float(np.median(pr["HE"][open_field_mask(pr["HE"], p, MEASURED["ceiling"])]))
    n, vs = nnps(H_sim.astype(float), open_field_mask(H_sim.astype(float), p, MEASURED["ceiling"]), p, bands=bands)
    sel = [b for b, (lo, hi) in enumerate(zip(bands[:-1], bands[1:])) if hi <= 2.0]
    k = float(np.mean([a_ * Q[b] / S_pair / vs[b] for b in sel]))
    f, e = float(np.sqrt(B) / p), float(np.sqrt(C) / p)
    share = lambda S: dict(quantum=float(a_ * Q[0] / S / pred(S)[0]), fixed_pattern=float(B / pred(S)[0]), electronic=float(C / S ** 2 / pred(S)[0]))
    # low-frequency structured noise: excess of the measured NNPS over the quantum + fixed-pattern model (median)
    sb = (0.05, 0.1, 0.15, 0.2, 0.3, 0.45, 0.6, 1.0)
    Qs = np.array([(Hn ** 2)[(FR >= lo) & (FR < hi)].mean() for lo, hi in zip(sb[:-1], sb[1:])])
    exs = []
    for e_ in [e for e in E if e["id"] in {r["id"] for r in rows}]:
        H = pydicom.dcmread(MS.locate(e_["H"])).pixel_array.astype(float)
        H = H[:, ::-1] if MEASURED_SETS[MEASURED["set"]]["mirror"] else H
        m = open_field_mask(H, p, MEASURED["ceiling"]); n, v = nnps(H, m, p, bands=sb); S = float(np.median(H[m]))
        exs.append([x - (a_ * q / S + B + C / S ** 2) for x, q in zip(v, Qs)])
    structured = dict(bands_cy_mm=list(sb), excess_nnps_mm2=np.median(exs, 0).tolist(), min=np.min(exs, 0).tolist(), max=np.max(exs, 0).tolist())
    print("  structured noise (excess, median of %d open fields, 1e-7 mm2): %s" % (len(exs), " ".join("%.2f" % (x * 1e7) for x in structured["excess_nnps_mm2"])))
    out = dict(set=MEASURED["set"], pair=ex["id"], mtf=DETECTOR_RESPONSE["mtf"], k=k, f=f, e=e, structured=structured,
               fit=dict(terms=terms, a=float(a_), B=float(B), C=float(C), bands_cy_mm=list(bands), Q=Q.tolist(),
                        rms_relative_residual=float(np.sqrt(np.mean(res ** 2))), residual_by_band=res.mean(0).tolist()),
               shares_0p6_1_cy_mm=dict(open_field=share(S_pair), phantom_interior_1650=share(1650.0)),
               simulated_quantum_nnps_mm2=vs, S_pair=S_pair, exposures=rows,
               note="fitted on HE open fields; the same k, f, e are used at LE [A]")
    path = HERE / "data" / ("detector_noise_%s.json" % MEASURED["set"])
    json.dump(out, open(path, "w"), indent=1)
    print("fit (%s) over %d open fields (%s): a %.3g, B %.3g mm2, C %.3g mm2; rms relative residual %.3f" % (" + ".join(terms), len(rows), "/".join(levels), a_, B, C, out["fit"]["rms_relative_residual"]))
    print("  residual by band (mean, measured/fit - 1): %s" % " ".join("%+.3f" % v for v in res.mean(0)))
    print("  -> k %.3f (quantum factor), f %.4f (fixed-pattern rms, relative), e %.1f (electronic rms, stored units)" % (k, f, e))
    for nm, sh in out["shares_0p6_1_cy_mm"].items():
        print("  noise shares at 0.6-1 cy/mm, %s: quantum %.2f, fixed pattern %.2f, electronic %.2f" % (nm, sh["quantum"], sh["fixed_pattern"], sh["electronic"]))
    print("written", path)


def preview(data, x, y, path, bg):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 6.5))
    ext = [x[0], x[-1], y[-1], y[0]]
    for a_, beam in zip(ax, ("LE", "HE")):
        im = data[beam]["stored"].astype(float)
        lo, hi = np.percentile(im, [1, 99.5])
        a_.imshow(im, cmap="gray_r", vmin=lo, vmax=hi, extent=ext)   # inverse display, as MONOCHROME1
        a_.set_title("%s simulated stored values (%s)" % (beam, bg)); a_.set_xlabel("x (mm)")
    fig.tight_layout(); fig.savefig(path, dpi=90); plt.close(fig)


if __name__ == "__main__":
    main()

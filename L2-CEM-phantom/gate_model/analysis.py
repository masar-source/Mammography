#!/usr/bin/env python3
"""
Dual-energy analysis of the simulated phantom images (L2 v2 or real L2, --phantom), with the same quantities as the
measured data
(investigation_0035_0037/REPORT.md):
  Lambda_k = -ln(I_k / I0_k), I = primary (projector.py) + scatter (postprocess.py), I0 = open field
  w*     structure-cancelling weight: regression of Lambda_H on Lambda_L in the iodine-free background, band 2-35 px
         (difference of Gaussians, as REPORT §5)
  r      residual of the Siemens weight W: r = 1 - W / w*
  DES    Lambda_H - W Lambda_L (Siemens _W for the LE kV) and Lambda_H - w* Lambda_L
  texture  power spectrum of the log images, 21.8 mm ROIs, quadratic detrend, Hann window, power-law fit 0.2-0.7 cy/mm
           (REPORT §11.1): beta and the power at 0.2 cy/mm. ROIs lie in the flat region of the container, away from the
           iodine discs, on a quarter-ROI lattice (overlapping)
  iodine   mean DES signal of each disc against a local background ring
  SPR    scatter-to-primary ratio in the phantom interior
Optional quantum noise (--noise): Gaussian with the variance of the detected energy, from the spectrum fluence and the
mAs of the measured exposure (approximation: mean detected energy per photon = detected-spectrum mean [A]).

Usage: python analysis.py [--phantom L2v2|L2] [--background PA12_MCT|PMMA_paraffin] [--set 0035|july28high]
                          [--le LE] [--he HE] [--bin 1] [--scatter-le LE_PA12_MCT --scatter-he HE_PA12_MCT] [--noise]
                          [--tag name]
Output: output/analysis/<tag>/summary.json and figures (tag default: <case>[_noise], case = config.case_name())
"""
import argparse, json
import numpy as np
from scipy import ndimage as ndi
from config import (HERE, SYSTEM, BEAMS, SIEMENS_W, PHANTOM, OUT_DIR, MEASURED, spectra_dir, add_case_args,
                    apply_case_args)
from projector import pixel_grid, P as PARAMS, SID, out_dir as projector_dir
import phantoms
from revelation_gate import geometry_z

ROI_MM = 256 * 0.085          # texture ROI 21.8 mm (256 px at full resolution, REPORT §11.1)


def load(a):
    pr = projector_dir(); mc = OUT_DIR / "mc"
    tag = f"_bin{a.bin}" if a.bin > 1 else ""
    d = {}
    for k, beam, sc in (("L", a.le, a.scatter_le), ("H", a.he, a.scatter_he)):
        d["P" + k] = np.load(pr / f"primary_{beam}{a.variant}{tag}.npy").astype(np.float64)
        d["P0" + k] = np.load(pr / f"primary_{beam}_open{tag}.npy").astype(np.float64)
        d["S" + k] = np.load(mc / f"scatter_{sc}{tag}.npy").astype(np.float64) if sc else np.zeros_like(d["P" + k])
        d["S0" + k] = np.load(mc / f"scatter_{a.scatter_open}_{beam}{tag}.npy").astype(np.float64) \
            if a.scatter_open else np.zeros_like(d["P" + k])
    d["maps"] = np.load(pr / f"maps{a.variant}{tag}.npz")
    return d


def add_noise(I, beam, rng):
    """detected energy per pixel per photon/sr -> add Poisson-like noise for the measured exposure"""
    spec = np.load(spectra_dir() / f"{beam}.npz")
    n_sr = spec["fluence_per_mAs_1m"] * 1e4 * BEAMS[beam]["mAs"]           # photons per sr (cm2 at 1 m = 1e-4 sr)
    e_mean = (spec["centres_keV"] * spec["weights"]).sum() * 1.15           # [A] hardening by the phantom
    E = I * n_sr                                                            # keV per pixel
    return (E + rng.normal(size=I.shape) * np.sqrt(E * e_mean)) / n_sr


def masks(d, x, y, p):
    """interior: rays that stay in the uniform-thickness region of the container over the whole phantom height and
    > 5 mm from the image edge; iodine discs excluded (+3 mm).
    L2 v2: inside the flat region (the ray radius is largest at the phantom bottom). Real L2: rays that cross the whole
    cavity height and the whole plate thickness, i.e. stay clear of the front wall, the curved wall and the plate rim."""
    z0 = geometry_z()["phantom"][0]
    xx, yy = np.meshgrid(x, y)
    if phantoms.is_l2():
        g = phantoms.l2()
        uz = d["maps"]["uz"]
        interior = ((d["maps"]["cavity"] >= (g["zt"] - g["zb"]) / uz - 0.01) &
                    (d["maps"]["plate"] >= g["t_plate"] / uz - 0.01) & (xx > SYSTEM["missed_tissue_mm"] + 5.0))
    else:
        rho_bot = np.hypot(xx, yy) * (SID - z0) / SID
        interior = (rho_bot < PARAMS["R_flat"]) & (xx > SYSTEM["missed_tissue_mm"] + 5.0)
    io = d["maps"]["iodine_g_cm2"] > 0
    io = ndi.binary_dilation(io, iterations=int(round(3.0 / p))) if io.any() else io
    return interior, interior & ~io, io


def band(img, p, lo_px=2, hi_px=35):
    return ndi.gaussian_filter(img, lo_px) - ndi.gaussian_filter(img, hi_px)


def w_star(lL, lH, m, p):
    bL, bH = band(lL, p)[m], band(lH, p)[m]
    return float((bL * bH).sum() / (bL * bL).sum())


# ---------------------------------------------------------------- texture (REPORT §11.1)
def _texture_tools(p):
    N_ROI = int(round(ROI_MM / p))
    win = np.outer(np.hanning(N_ROI), np.hanning(N_ROI))
    yv, xv = np.mgrid[0:N_ROI, 0:N_ROI]
    Ad = np.c_[np.ones(N_ROI * N_ROI), xv.ravel(), yv.ravel(), (xv * xv).ravel(), (yv * yv).ravel(), (xv * yv).ravel()]
    f = np.fft.fftfreq(N_ROI, d=p); FX, FY = np.meshgrid(f, f); FR = np.sqrt(FX ** 2 + FY ** 2)
    edges = np.logspace(np.log10(0.05), np.log10(5.0), 31)
    return N_ROI, win, Ad, FR, edges, np.sqrt(edges[:-1] * edges[1:])


def spectrum(img, mask, p):
    N_ROI, win, Ad, FR, edges, fc = _texture_tools(p)
    P, n = 0, 0
    for y0 in range(0, img.shape[0] - N_ROI, N_ROI // 4):           # quarter-ROI lattice: the iodine array leaves
        for x0 in range(0, img.shape[1] - N_ROI, N_ROI // 4):       # little room (ROIs overlap, not independent)
            sl = (slice(y0, y0 + N_ROI), slice(x0, x0 + N_ROI))
            if not mask[sl].all():
                continue
            a = img[sl].astype(float)
            c = np.linalg.lstsq(Ad, a.ravel(), rcond=None)[0]
            F = np.fft.fft2(win * (a - (Ad @ c).reshape(N_ROI, N_ROI)))
            P = P + np.abs(F) ** 2; n += 1
    if n == 0:
        return fc, np.full(len(fc), np.nan), 0
    P = P / n * p * p / (win ** 2).sum()
    rad = np.array([P[(FR >= lo) & (FR < hi)].mean() if ((FR >= lo) & (FR < hi)).any() else np.nan
                    for lo, hi in zip(edges[:-1], edges[1:])])
    return fc, rad, n


def fit_beta(fc, rad, lo=0.2, hi=0.7):
    s = (fc >= lo) & (fc <= hi) & np.isfinite(rad) & (rad > 0)
    if s.sum() < 3:
        return np.nan, np.nan
    q = np.polyfit(np.log10(fc[s]), np.log10(rad[s]), 1)
    return float(-q[0]), float(10 ** np.polyval(q, np.log10(0.2)))


def iodine_signals(img, x, y, p):
    z0 = geometry_z()["phantom"][0]
    m = SID / (SID - z0 - phantoms.plate()[0])
    out = []
    xx, yy = np.meshgrid(x, y)
    for r in phantoms.discs():
        xc, yc, rd = (r["x"] + phantoms.x_offset()) * m, r["y"] * m, r["diam"] / 2 * m
        i0, j0 = int(np.argmin(np.abs(y - yc))), int(np.argmin(np.abs(x - xc)))
        k = int(np.ceil((rd + 4.0) / p))
        sl = (slice(max(i0 - k, 0), i0 + k + 1), slice(max(j0 - k, 0), j0 + k + 1))
        rr = np.hypot(xx[sl] - xc, yy[sl] - yc)
        inner, ring = rr < 0.6 * rd, (rr > rd + 1.0) & (rr < rd + 3.0)
        if inner.sum() == 0 or ring.sum() == 0:
            continue
        sig = img[sl][inner].mean() - img[sl][ring].mean()
        out.append(dict(name=r["name"], kind=r["kind"], row=r["key"], diameter_mm=r["diam"],
                        signal=float(sig), ring_sd=float(img[sl][ring].std())))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--le", default="LE"); ap.add_argument("--he", default="HE")
    ap.add_argument("--bin", type=int, default=1)
    ap.add_argument("--variant", default="", help="projector variant suffix, e.g. _nospheres")
    ap.add_argument("--scatter-le", default=None); ap.add_argument("--scatter-he", default=None)
    ap.add_argument("--scatter-open", default=None, help="prefix of open-field scatter files: scatter_<prefix>_<beam>")
    ap.add_argument("--noise", action="store_true"); ap.add_argument("--seed", type=int, default=1)
    add_case_args(ap)                       # selects the projector images of output/projector/<case>/
    ap.add_argument("--tag", default=None, help="output folder name (default: <case>[_noise])")
    a = ap.parse_args()
    case = apply_case_args(a)
    bg_name = PHANTOM["background"]
    a.tag = a.tag or case + ("_noise" if a.noise else "")
    out = OUT_DIR / "analysis" / a.tag; out.mkdir(parents=True, exist_ok=True)
    d = load(a)
    x, y, p = pixel_grid(a.bin)
    rng = np.random.default_rng(a.seed)
    lam, spr = {}, {}
    interior, bg, io = masks(d, x, y, p)
    for k, beam in (("L", a.le), ("H", a.he)):
        I, I0 = d["P" + k] + d["S" + k], d["P0" + k] + d["S0" + k]
        if a.noise:
            I = add_noise(I, beam, rng)
        lam[k] = -np.log(np.clip(I, 1e-30, None) / I0)
        spr[k] = float(np.median((d["S" + k] / d["P" + k])[interior])) if (a.scatter_le or a.scatter_he) else None
    kv = int(round(BEAMS[a.le]["kvp"])) if a.le in BEAMS else int(a.le[2:])
    W = SIEMENS_W.get(kv)
    ws = w_star(lam["L"], lam["H"], bg, p)
    res = dict(phantom=PHANTOM["model"], case=case, measured_set=MEASURED["set"], background=bg_name, sphere_material=PHANTOM["sphere_material"], liquid_material=PHANTOM["liquid_material"],
               le=a.le, he=a.he, noise=a.noise, siemens_W=W, w_star=ws, residual_r=(1 - W / ws) if W else None,
               spr_interior=spr, mean_lambda_interior={k: float(lam[k][interior].mean()) for k in lam})
    imgs = {"LE": lam["L"], "HE": lam["H"], "DES_wstar": lam["H"] - ws * lam["L"]}
    if W:
        imgs["DES_siemensW"] = lam["H"] - W * lam["L"]
    tex = {}
    for name, im in imgs.items():
        fc, rad, n = spectrum(im, bg, p)
        beta, p02 = fit_beta(fc, rad)
        tex[name] = dict(beta=beta, power_at_0p2_mm2=p02, n_roi=n, f=fc.tolist(), nps=rad.tolist())
    res["texture"] = {k: {kk: v for kk, v in t.items() if kk not in ("f", "nps")} for k, t in tex.items()}
    res["iodine"] = {name: iodine_signals(imgs[name], x, y, p) for name in imgs if name.startswith("DES")}
    json.dump(res, open(out / "summary.json", "w"), indent=1)
    np.savez_compressed(out / "images.npz", **{k: v.astype(np.float32) for k, v in imgs.items()}, x=x, y=y)
    figures(imgs, tex, bg, x, y, out)
    print(json.dumps({k: v for k, v in res.items() if k != "iodine"}, indent=1))


def figures(imgs, tex, bg, x, y, out):
    """images.png: top row at one common contrast width (that of the LE image), so the subtracted images can be
    compared with LE; bottom row with each window stretched to its own image, which magnifies small residuals
    (a stretched DES shows ~1 % fluctuations, e.g. Monte Carlo scatter noise, as if they were texture)"""
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, len(imgs), figsize=(4.2 * len(imgs), 10.4), squeeze=False)
    ext = [x[0], x[-1], y[-1], y[0]]
    lo, hi = np.percentile(imgs["LE"][bg], [1, 99])
    width = 1.6 * (hi - lo)
    for j, (k, im) in enumerate(imgs.items()):
        v = im[bg]
        c = np.median(v)
        ax[0, j].imshow(im, cmap="gray", vmin=c - width / 2, vmax=c + width / 2, extent=ext)
        ax[0, j].set_title("%s\ncommon window (LE width %.3f)" % (k, width), fontsize=9)
        a, b = np.percentile(v, [1, 99])
        ax[1, j].imshow(im, cmap="gray", vmin=a - 0.3 * (b - a), vmax=b + 0.3 * (b - a), extent=ext)
        ax[1, j].set_title("%s\nown window, stretched (width %.3f)" % (k, 1.6 * (b - a)), fontsize=9)
        for i in (0, 1):
            ax[i, j].set_xlabel("x from chest wall (mm)")
    ax[0, 0].set_ylabel("y (mm)"); ax[1, 0].set_ylabel("y (mm)")
    fig.tight_layout(); fig.savefig(out / "images.png", dpi=100); plt.close(fig)
    fig, a_ = plt.subplots(figsize=(6, 4.5))
    for k, t in tex.items():
        a_.loglog(t["f"], t["nps"], label="%s  beta %.2f" % (k, t["beta"]))
    a_.axvspan(0.2, 0.7, color="0.9")
    a_.set_xlabel("spatial frequency (cy/mm)"); a_.set_ylabel("power (mm$^2$)"); a_.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(out / "texture_spectra.png", dpi=110); plt.close(fig)


if __name__ == "__main__":
    main()

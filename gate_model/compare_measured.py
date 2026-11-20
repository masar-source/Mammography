#!/usr/bin/env python3
"""
One-to-one comparison of the measurement-like simulated images (measlike.py) with a measured reference set
(config.MEASURED_SETS, read by measured_set.py). Both sides are analysed by the same code, in the same domain
(Lambda = -ln(stored value), as in REPORT.md), with masks built the same way:
  phantom   HE < 0.64 x HE open field (largest component, holes filled)
  interior  >= 15 mm from the phantom edge, iodine discs excluded, as REPORT §11.1
  core      >= 30 mm from the edge, discs excluded (weight estimate; REPORT §5: w* rises within 30 mm of the edge)
Disc exclusion: the 7 x 5 array as a box (+12.75 mm), the two 9 mm discs (if known) as circles (+8 mm).
  measured 0035:      7 x 5 lattice of 16 mm pitch fitted to the image (old plate, layout unknown; REPORT §1.3)
  measured July:      the July template (disc positions in HIGHMAS_28KV2_00), shifted to each pair by maximising the
                      iodine response at the disc centres
  simulated:          the model's disc layout (phantoms.discs), projected to the detector

Cases:
  --phantom L2  --set july28high   the real L2 with the new plate against one July pair (28 kV, high mAs, unit 10159,
                                   --pair, default HIGHMAS_28KV2_00): the same phantom on both sides (defaults for L2)
  --phantom L2v2 --set 0035        the v2 design against the old-L2 pair 0035 (different phantoms; interior background
                                   comparisons only)
One measured pair is compared with its simulated twin from measlike.py (the exact LE and HE mAs of that pair's images
and its HE open level). Each pair is treated on its own: another pair is a separate run with its own --pair; nothing
is averaged over pairs.

Metrics: open-field and interior stored values, HE interior transmission, high-frequency noise power (2.5-4.5 cy/mm),
texture beta and power at 0.2 cy/mm (REPORT §11.1 method), structure-cancelling weight w* (band 2-35 px), residual at
the Siemens weight r = 1 - W/w*, the scatter "bowl" of the recombined image versus distance from the edge, and the
iodine signal of each disc (recombined image at W, disc centre minus a ring).

Output: output/compare/<set>_<case>/<pair>/metrics.json, comparison.png, spectra_profiles.png
Usage: python compare_measured.py [--phantom L2] [--background PMMA_paraffin] [--set july28high] [--pair ID] [--bin 1]
                                  [--root "/Volumes/SanDisk/Masar/New L2 Measurements"]
"""
import argparse, json
import numpy as np
from scipy import ndimage as ndi
from config import OUT_DIR, MEASURED, MEASURED_SETS, PHANTOM, SIEMENS_W, BEAMS, SYSTEM, add_case_args, apply_case_args
import analysis as A
import measured_set as MS
import phantoms
from measured_set import phantom_mask
from projector import pixel_grid
from revelation_gate import geometry_z


def w_siemens():
    return SIEMENS_W[int(round(BEAMS["LE"]["kvp"]))]


def response(lL, lH, m, p):
    """iodine response: recombined image (w 0.355) band-passed (1-6 mm), clipped, zero outside a 10 mm eroded mask"""
    S0 = lH - 0.355 * lL
    resp = ndi.gaussian_filter(S0, 1.0 / p) - ndi.gaussian_filter(S0, 6.0 / p)
    inner = ndi.binary_erosion(m, iterations=int(10 / p))
    v = resp[inner]; med = np.median(v); sd = 1.4826 * np.median(np.abs(v - med))
    return np.where(inner, np.clip(resp - med, -5 * sd, 5 * sd), 0.0)


def discs_lattice(resp, p):
    """0035: 7 x 5 lattice (16 mm pitch) fitted to the response; disc list (row, col, radius px, kind)"""
    step = 16.0 / p; ry, rx = np.arange(7) * step, np.arange(5) * step
    def score(y0, x0):
        return resp[(y0 + ry[:, None]).astype(int), (x0 + rx[None, :]).astype(int)].sum()
    best = max(((score(y0, x0), y0, x0) for y0 in np.arange(0, resp.shape[0] - ry[-1] - 1, 3)
                for x0 in np.arange(0, resp.shape[1] - rx[-1] - 1, 3)))
    _, y0, x0 = max((score(best[1] + dy, best[2] + dx), best[1] + dy, best[2] + dx)
                    for dy in range(-3, 4) for dx in range(-3, 4))
    return [(y0 + a, x0 + b, 2.0 / p, "array") for a in ry for b in rx]


def discs_template(resp, template, p, search_mm=15.0):
    """July: template disc positions (pixels of the reference image) shifted to this image; returns the disc list and
    the shift (rows, cols). Only the iodine array is scored (the 9 mm discs follow the same shift)."""
    T = [d for d in template["discs"]]
    arr = np.array([(d["y"], d["x"]) for d in T if d["kind"] == "array"])
    def score(dy, dx):
        r, c = np.rint(arr[:, 0] + dy).astype(int), np.rint(arr[:, 1] + dx).astype(int)
        ok = (r >= 0) & (r < resp.shape[0]) & (c >= 0) & (c < resp.shape[1])
        return resp[r[ok], c[ok]].sum() if ok.all() else -np.inf
    k = int(search_mm / p)
    best = max((score(dy, dx), dy, dx) for dy in range(-k, k + 1, 2) for dx in range(-k, k + 1, 2))
    _, dy, dx = max((score(best[1] + a, best[2] + b), best[1] + a, best[2] + b) for a in (-1, 0, 1) for b in (-1, 0, 1))
    if max(abs(dy), abs(dx)) >= k - 1:
        print("WARNING: disc array found at the edge of the search range (shift %d, %d px): check the orientation of "
              "the images and the template" % (dy, dx))
    return [(d["y"] + dy, d["x"] + dx, d["radius_px"], d["kind"]) for d in T], (int(dy), int(dx))


def discs_simulated(x, y, p):
    """model layout projected to the detector (pixel indices of the simulated grid)"""
    sid = SYSTEM["sid_mm"]
    m = sid / (sid - geometry_z()["phantom"][0] - phantoms.plate()[0])
    out = []
    for d in phantoms.discs():
        xd, yd = (d["x"] + phantoms.x_offset()) * m, d["y"] * m
        out.append(((yd - y[0]) / p, (xd - x[0]) / p, d["diam"] / 2 * m / p, "array" if d["kind"] == "array" else "large"))
    return out


def exclusion(shape, discs, p):
    """array box + 12.75 mm (REPORT §11.1) and the 9 mm discs + 8 mm"""
    rows, cols = np.mgrid[0:shape[0], 0:shape[1]]
    arr = np.array([(r, c) for r, c, _, k in discs if k == "array"])
    marg = 12.75 / p
    ex = ((rows > arr[:, 0].min() - marg) & (rows < arr[:, 0].max() + marg) &
          (cols > arr[:, 1].min() - marg) & (cols < arr[:, 1].max() + marg))
    for r, c, rad, k in discs:
        if k != "array":
            ex |= np.hypot(rows - r, cols - c) < rad + 8.0 / p
    return ex


def iodine(des, discs, p):
    """recombined-image signal of each array disc: mean over 0.6 radius minus the ring 1-3 mm outside the edge"""
    out = []
    for r, c, rad, k in discs:
        if k != "array":
            continue
        n = int(rad + 3.5 / p) + 1
        i0, j0 = int(round(r)), int(round(c))
        sl = (slice(max(i0 - n, 0), i0 + n + 1), slice(max(j0 - n, 0), j0 + n + 1))
        rr = np.hypot(*np.mgrid[sl[0], sl[1]] - np.array([r, c])[:, None, None])
        inner, ring = rr < 0.6 * rad, (rr > rad + 1.0 / p) & (rr < rad + 3.0 / p)
        out.append(float(des[sl][inner].mean() - des[sl][ring].mean()) if inner.any() and ring.any() else np.nan)
    return out


def analyse(L, H, p, discs):
    W = w_siemens()
    m, open_he = phantom_mask(H)
    lL, lH = -np.log(np.clip(L, 1, None)), -np.log(np.clip(H, 1, None))
    dist = ndi.distance_transform_edt(m) * p
    ex = exclusion(L.shape, discs, p)
    interior = m & (dist >= 15) & ~ex
    core = m & (dist >= 30) & ~ex
    outside = (~ndi.binary_dilation(m, iterations=int(10 / p))) & (H < MEASURED["ceiling"] - 1) & (H > 0.9 * open_he)
    res = dict(open_field_HE=open_he,
               open_field_LE=float(np.median(L[outside])) if outside.any() else None,
               open_field_LE_clipped_fraction=float((L[outside] >= MEASURED["ceiling"] - 1).mean()) if outside.any() else None,
               interior_median_LE=float(np.median(L[interior])), interior_median_HE=float(np.median(H[interior])),
               pixels_interior=int(interior.sum()), pixels_core=int(core.sum()))
    res["HE_interior_transmission"] = res["interior_median_HE"] / open_he
    ws = A.w_star(lL, lH, core, p)
    res["w_star"] = ws
    res["residual_r_siemens"] = 1 - W / ws
    spec = {}
    des = lH - W * lL
    for name, im in (("LE", lL), ("HE", lH), ("DES_siemensW", des)):
        fc, rad, n = A.spectrum(im, interior, p)
        beta, p02 = A.fit_beta(fc, rad)
        hf = (fc >= 2.5) & (fc <= 4.5) & np.isfinite(rad)
        res[name] = dict(beta=beta, power_at_0p2_mm2=p02, n_roi=n, hf_noise_power_mm2=float(np.mean(rad[hf])) if hf.any() else None)
        spec[name] = (fc, rad)
    mk = m & ~ex
    ref = np.median(des[mk & (dist >= 45)])
    edges = np.arange(0, 45, 2.5)
    prof = [float(np.median(des[mk & (dist >= q) & (dist < q + 2.5)]) - ref) if (mk & (dist >= q) & (dist < q + 2.5)).any()
            else np.nan for q in edges]
    res["bowl_DES_vs_edge_distance_mm"] = dict(centres=(edges + 1.25).tolist(), values=prof)
    res["iodine_DES_signal"] = iodine(des, discs, p)
    res["spectra"] = dict(f=spec["LE"][0].tolist(), **{k: v[1].tolist() for k, v in spec.items()})
    return res, dict(lL=lL, lH=lH, des=des, mask=m, interior=interior, spec=spec)


def main():
    ap = argparse.ArgumentParser()
    add_case_args(ap)
    ap.add_argument("--bin", type=int, default=1)
    ap.add_argument("--root", default=None, help="folder of 'New L2 Measurements' if not at %s" % MS.JULY_ROOT)
    a = ap.parse_args()
    case = apply_case_args(a)
    st, pid = MEASURED["set"], MEASURED["pair"]
    tag = f"_bin{a.bin}" if a.bin > 1 else ""
    x, y, p_sim = pixel_grid(a.bin)
    template = json.load(open(MEASURED_SETS[st]["template"])) if MEASURED_SETS[st].get("template") else None
    if template is None and PHANTOM["model"] != "L2v2":
        print("NOTE: no disc template for set %s: the disc array is found as a 16 mm lattice" % st)
    # its simulated twin (same LE / HE mAs and open level, measlike.py)
    sd = OUT_DIR / "measlike" / case / st / pid
    if not (sd / f"LE{tag}.npy").exists():
        raise SystemExit("no simulated pair %s: run measlike.py with the same --phantom/--background/--set/--pair first" % sd)
    sim = {b: np.load(sd / f"{b}{tag}.npy").astype(np.float64) for b in ("LE", "HE")}
    meta = json.load(open(sd / f"meta{tag}.json"))
    # the measured pair
    pr = next(MS.pairs(st, a.root, ids=[pid]), None)
    if pr is None:
        raise SystemExit("measured pair %s not found in set %s" % (pid, st))
    p_meas = pr["pixel_mm"]
    lL, lH = -np.log(np.clip(pr["LE"], 1, None)), -np.log(np.clip(pr["HE"], 1, None))
    m, _ = phantom_mask(pr["HE"])
    resp = response(lL, lH, m, p_meas)
    if template:
        dm_discs, sh = discs_template(resp, template, p_meas)
    else:
        dm_discs, sh = discs_lattice(resp, p_meas), None
    rm, dmm = analyse(pr["LE"], pr["HE"], p_meas, dm_discs)
    rs, dsim = analyse(sim["LE"], sim["HE"], p_sim, discs_simulated(x, y, p_sim))
    print("  %s: mAs LE %.2f HE %.2f (simulated with %.2f / %.2f) | disc shift %s"
          % (pid, pr["header_LE"]["mAs"], pr["header_HE"]["mAs"], meta["mAs"]["LE"], meta["mAs"]["HE"], sh), flush=True)
    out = OUT_DIR / "compare" / f"{st}_{case}" / pid; out.mkdir(parents=True, exist_ok=True)
    json.dump(dict(case=case, phantom=PHANTOM["model"], background=PHANTOM["background"], measured_set=st, pair=pid,
                   mAs_measured=dict(LE=pr["header_LE"]["mAs"], HE=pr["header_HE"]["mAs"]), mAs_simulated=meta["mAs"],
                   disc_shift_px=sh, measured=rm, simulated=rs), open(out / "metrics.json", "w"), indent=1, default=float)
    print_table(rm, rs, case, st, pid)
    figures(pr, dmm, rm, sim, dsim, rs, p_sim, out, case, st)
    print("written", out)


def _fmt(v, f):
    return f % v if v is not None and np.isfinite(v) else "-"


def print_table(rm, rs, case, st, pid):
    rows = [("HE open field (stored)", "open_field_HE", "%.0f"), ("LE open field (stored)", "open_field_LE", "%.0f"),
            ("LE interior median (stored)", "interior_median_LE", "%.0f"), ("HE interior median (stored)", "interior_median_HE", "%.0f"),
            ("HE interior transmission", "HE_interior_transmission", "%.3f"), ("w* (>= 30 mm from edge)", "w_star", "%.3f"),
            ("residual r at W %.3f" % w_siemens(), "residual_r_siemens", "%+.3f")]
    print("\n%-34s %14s %14s   (%s pair %s and its simulated twin, case %s)" % ("", "measured", "simulated", st, pid, case))
    for lab, k, f in rows:
        print("%-34s %14s %14s" % (lab, _fmt(rm[k], f), _fmt(rs[k], f)))
    for img in ("LE", "HE", "DES_siemensW"):
        for k, lab, f in (("beta", "beta (0.2-0.7 cy/mm)", "%.2f"), ("power_at_0p2_mm2", "power at 0.2 cy/mm (mm2)", "%.2e"),
                          ("hf_noise_power_mm2", "noise power 2.5-4.5 cy/mm (mm2)", "%.2e")):
            print("%-34s %14s %14s" % ("%s %s" % (img, lab), _fmt(rm[img][k], f), _fmt(rs[img][k], f)))
    print("%-34s %14s %14s" % ("texture ROIs (LE)", rm["LE"]["n_roi"], rs["LE"]["n_roi"]))
    im, iv = np.array(rm["iodine_DES_signal"], float), np.array(rs["iodine_DES_signal"], float)
    if MEASURED_SETS[st].get("template") and len(im) == len(iv) == 35:
        q = (iv / im).reshape(7, 5)
        print("iodine DES signal, array discs: simulated / measured median ratio %.2f (rows 0-2, 4-7 mm discs: %.2f)"
              % (np.nanmedian(q), np.nanmedian(q[:3, 2:])))


def figures(pr, dm, rm, sim, ds, rs, ps, out, case, st):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    pm = pr["pixel_mm"]
    W = w_siemens()
    fig, ax = plt.subplots(3, 2, figsize=(12, 17))
    lab_m, lab_s = "measured %s (%s)" % (st, pr["id"]), "simulated twin (%s)" % case
    for i, beam in enumerate(("LE", "HE")):
        mm = np.median(pr[beam][dm["interior"]])
        lo, hi = np.percentile(pr[beam][dm["interior"]] / mm, [0.5, 99.5]); half = 0.65 * (hi - lo)
        for j, (im, d, p, lab) in enumerate(((pr[beam], dm, pm, lab_m), (sim[beam], ds, ps, lab_s))):
            c = np.median(im[d["interior"]])
            ax[i, j].imshow(im, cmap="gray_r", vmin=c * (1 - half), vmax=c * (1 + half), extent=[0, im.shape[1] * p, im.shape[0] * p, 0])
            ax[i, j].set_title("%s %s, stored values; interior median %.0f, window +-%.0f %%" % (beam, lab, c, 100 * half), fontsize=9)
    v = dm["des"][dm["interior"]]; w = 1.6 * (np.percentile(v, 99) - np.percentile(v, 1))
    for j, (d, p, lab) in enumerate(((dm, pm, lab_m), (ds, ps, lab_s))):
        c = np.median(d["des"][d["interior"]])
        ax[2, j].imshow(np.where(d["mask"], d["des"], np.nan), cmap="gray", vmin=c - w / 2, vmax=c + w / 2,
                        extent=[0, d["des"].shape[1] * p, d["des"].shape[0] * p, 0])
        ax[2, j].set_title("-ln(HE) - %.3f(-ln(LE)) %s, common window width %.3f" % (W, lab, w), fontsize=9)
    for a_ in ax.ravel():
        a_.set_xlabel("x from the chest-wall edge (mm)"); a_.set_ylabel("y (mm)")
    fig.tight_layout(); fig.savefig(out / "comparison.png", dpi=80); plt.close(fig)
    fig, ax = plt.subplots(1, 4, figsize=(20, 4.6))
    for k, name in enumerate(("LE", "HE", "DES_siemensW")):
        for r_, lab, col, st_ in ((rm, "measured %s" % pr["id"], "C0", "-"), (rs, "simulated twin", "C1", "--")):
            ax[k].loglog(r_["spectra"]["f"], r_["spectra"][name], st_, color=col, label=lab)
        ax[k].axvspan(0.2, 0.7, color="0.92"); ax[k].set_title(name); ax[k].set_xlabel("cy/mm"); ax[k].legend(fontsize=8)
    ax[0].set_ylabel("power spectrum of -ln(stored) (mm$^2$)")
    for r_, lab, fmt in ((rm, "measured %s" % pr["id"], "o-"), (rs, "simulated twin", "s--")):
        b = r_["bowl_DES_vs_edge_distance_mm"]
        ax[3].plot(b["centres"], b["values"], fmt, label=lab)
    ax[3].set_xlabel("distance from the phantom edge (mm)"); ax[3].set_ylabel("recombined - interior (log units)")
    ax[3].set_title("scatter bowl at W = %.3f" % W); ax[3].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(out / "spectra_profiles.png", dpi=90); plt.close(fig)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Predicted measurement of the redesigned background (REPORT §10.4: PA12 + 5.5 wt% CaCO3 spheres in MCT oil) from the
measured current L2 pair and two simulations of it that differ only in the sphere and liquid materials.

The simulation reproduces the current L2 inside the phantom only up to known, unexplained errors (texture power
1.2-1.8x, HE level -7 %, PROJECT_NOTES §6). They are taken as common to both backgrounds (same container, plate,
sphere positions, spectra, detector model, mAs and noise seed), so the redesign is read as a change from the
current phantom, never in absolute terms:
  positive quantities (weights, levels, texture and noise power, disc contrast and noise):
      predicted = measured current x (simulated redesign / simulated current)
  signed quantities (beta, rim heights):
      predicted = measured current + (simulated redesign - simulated current)
  r = 1 - W/w* from the predicted w*, with the measured pair's header W (Siemens' weight does not depend on the
  object). Errors that depend on the materials do not cancel this way [A].
Disc SDNR: per image, the measured SDNRs are scaled by (median contrast ratio of the 4-7 mm discs of rows 0-2) /
(median noise ratio of all 35 array discs), redesign vs current simulation.
R (recombined) is the Siemens emulation at the header W on both sides (base: measured_emuR); everything else comes
from the measured pair as it is.

Inputs (pipeline_compare.py of both cases; nothing is recomputed from images):
  output/pipeline/<set>_<phantom>_<current>/{measured,measured_emuR,simulated}/results/<pair>.json, report/<pair>/report.json
  output/pipeline/<set>_<phantom>_<redesign>/simulated/results/<pair>.json, report/<pair>/report.json
Checks: the two twins must share pair, mAs, spectra, detector and noise model, seed and grid (meta of the twins) and
give the same disc registration; a mismatch is printed and written to the summary.
Targets (REPORT §10, model [E]): r = r_tissue (+0.07 / +0.09 / +0.08 / +0.12 at 28 / 30 / 32 / 34 kV, glandular
tissue in adipose tissue under Siemens' W); LE texture power 0.94 x the current L2; beta unchanged (same spheres).
Tissue reference (--reference, default gland_adipose): when that twin has been simulated and analysed, it is predicted
from the measured current L2 in exactly the same way and listed next to the redesign (column, summary lines, figure), so
the redesign is judged against tissue in the same geometry, chain and method instead of against the REPORT model. It
also gives the LE attenuation of a ≈ 64 % glandular (by mass) bed for the density comparison.
Not covered: AEC (the twins use the measured pair's mAs).

One pair, no averaging. Output: output/redesign/<set>_<pair>_<current>_vs_<redesign>/ summary.md, prediction.json,
prediction.png
Usage: python redesign_compare.py [--phantom L2] [--set july28high] [--pair ID] [--current PMMA_paraffin] [--redesign PA12_MCT]
"""
import argparse, json, os
import numpy as np
import config
from config import OUT_DIR, MEASURED, BACKGROUNDS, add_case_args, apply_case_args
import pipeline_compare as PC

R_TISSUE = {28: 0.07, 30: 0.09, 32: 0.08, 34: 0.12}      # REPORT §10.1 / §10.3 [E], L2 geometry, Hammerstein compositions
LE_POWER_TARGET = 0.94                                   # REPORT §10.3 [E]: (dmu_LE new / dmu_LE current)^2 for PA12 + 5.5 % CaCO3 / MCT
TILE_BINS = [(10, 20), (20, 30), (30, 45), (45, 140)]
RIM_SHOW = (6.75, 8.75, 11.25, 13.75, 17.5, 22.5, 27.5)  # ring centres (mm from the curved edge) listed in the table
SAME = ("phantom", "measured_set", "pair", "mAs", "detector", "noise_model", "offset_stored", "seed", "x_mm", "y_mm", "noise")


def case_root(bg):
    config.set_background(bg)
    return PC.root(), config.case_name()


def load_case(root, kinds, pid):
    out = {}
    for k in kinds:
        f = root / k / "results" / ("%s.json" % pid)
        if not f.exists():
            raise SystemExit("missing %s\n  run on that case first:  python pipeline_compare.py all --pair %s  (with its --background)" % (f, pid))
        out[k] = json.load(open(f))
    rep = root / "report" / pid / "report.json"
    if not rep.exists():
        raise SystemExit("missing %s\n  run:  python pipeline_compare.py report --pair %s  (with its --background)" % (rep, pid))
    rim = json.load(open(rep))["rim"]
    return out, {k: rim[k][0] for k in kinds}


# ------------------------------------------------------------------------------------------------------ checks
def checks(cur, new):
    """the two twins must differ only in the materials"""
    msgs = []
    a, b = cur["meta"].get("simulation") or {}, new["meta"].get("simulation") or {}
    if not a or not b:
        msgs.append("twin metadata missing (meta.json of measlike.py): consistency not checked")
    else:
        for k in SAME:
            if json.dumps(a.get(k), sort_keys=True) != json.dumps(b.get(k), sort_keys=True):
                msgs.append("twins differ in '%s': %s vs %s" % (k, json.dumps(a.get(k))[:80], json.dumps(b.get(k))[:80]))
        if os.path.basename(str(a.get("spectra", "")).rstrip("/")) != os.path.basename(str(b.get("spectra", "")).rstrip("/")):
            msgs.append("twins use different spectra: %s vs %s" % (a.get("spectra"), b.get("spectra")))
        if (a.get("sphere_material"), a.get("liquid_material")) == (b.get("sphere_material"), b.get("liquid_material")):
            msgs.append("twins have the same materials (%s / %s): nothing to compare" % (a.get("sphere_material"), a.get("liquid_material")))
    pa, pb = ({(d["kind"], d.get("row"), d.get("col")): (d["y"], d["x"]) for d in r["discs"]} for r in (cur, new))
    shift = max(np.hypot(pb[k][0] - pa[k][0], pb[k][1] - pa[k][1]) for k in pa if k in pb)
    if shift > 2:
        msgs.append("disc registration differs by up to %.1f px between the twins (same geometry: expected 0)" % shift)
    return msgs, float(shift)


# ------------------------------------------------------------------------------------------------------ prediction
def transfer(m, c, n, mode):
    if any(v is None or not np.isfinite(v) for v in (m, c, n)):
        return None
    if mode == "ratio":
        return m * n / c if c else None
    return m + (n - c)


def band_mean(r, img, lo, hi):
    f = np.array(r["fc"]); v = np.array([np.nan if x is None else x for x in r["texture"][img]["rad"]], float)
    i = (f >= lo) & (f < hi)
    return float(np.nanmean(v[i])) if i.any() else None


def tile_median(r, lo, hi):
    T = np.array(r["tiles"]); s = (T[:, 0] >= lo) & (T[:, 0] < hi)
    return float(np.median(T[s, 1])) if s.any() else None


def disc_table(r, img):
    """{(row, col): (contrast, noise)} of the 35 array discs in one image"""
    return {(d["row"], d["col"]): tuple(d[img]) for d in r["discs"] if d["kind"] == "array" and d.get(img) and d[img][1] > 0}


def disc_scaling(cur, new, img):
    """contrast ratio (4-7 mm discs, rows 0-2) and noise ratio (all array discs), redesign / current simulation"""
    c, n = disc_table(cur, img), disc_table(new, img)
    keys = [k for k in c if k in n]
    qc = [n[k][0] / c[k][0] for k in keys if k[0] <= 2 and k[1] >= 2 and c[k][0] > 0]
    qn = [n[k][1] / c[k][1] for k in keys]
    return float(np.median(qc)), float(np.median(qn)), [float(np.percentile(qn, 25)), float(np.percentile(qn, 75))]


def predict(meas, emu, cur, new, rim_meas, rim_emu, rim_cur, rim_new):
    """rows: (group, quantity, measured, sim current, sim redesign, mode, predicted, target, fmt)"""
    W, kv = meas["meta"]["W"], meas["meta"].get("kv")
    rows = []

    def add(group, label, get, mode, fmt, base=meas, target=""):
        m, c, n = get(base), get(cur), get(new)
        rows.append(dict(group=group, quantity=label, measured=m, sim_current=c, sim_redesign=n, mode=mode,
                         predicted=transfer(m, c, n, mode), target=target, fmt=fmt))
        return rows[-1]

    rt = R_TISSUE.get(kv)
    ws = add("Weights", "w* (0.06-0.15 cy/mm)", lambda r: r["w_low"], "ratio", "%.4f",
             target="%.3f (= W/(1 - r_tissue))" % (W / (1 - rt)) if rt is not None else "")
    rows.append(dict(group="Weights", quantity="r = 1 - W/w* (W = %.3f, header)" % W, measured=1 - W / ws["measured"],
                     sim_current=1 - W / ws["sim_current"], sim_redesign=1 - W / ws["sim_redesign"], mode="from w*",
                     predicted=None if ws["predicted"] is None else 1 - W / ws["predicted"],
                     target="%+.2f (r_tissue, %s kV)" % (rt, kv) if rt is not None else "", fmt="%+.3f"))
    add("Weights", "Siemens effective weight (emulated R)", lambda r: float(np.mean([b["w_eff"] for b in r["cross"][:2]])), "ratio",
        "%.4f", base=emu, target="unchanged (Siemens processing)")
    for i, b in enumerate(meas["cross"]):
        if b["f2"] <= 2.0:
            add("w*(f)", "%.2f-%.2f cy/mm" % (b["f1"], b["f2"]), (lambda i: lambda r: r["cross"][i]["w_fwd"])(i), "ratio", "%.3f")
    for lo, hi in TILE_BINS:
        # in the L2 13 of the 14 tiles >= 45 mm lie in the strip 5-16 mm from the chest-wall edge (inner face of the PMMA
        # front wall); the disc and marker exclusions leave almost no interior tiles there
        add("Local w*", "%d-%d mm from the curved edge%s" % (lo, hi, " (mostly the strip next to the chest wall)" if lo >= 45 else ""),
            (lambda lo, hi: lambda r: tile_median(r, lo, hi))(lo, hi), "ratio", "%.3f")
    add("Levels", "interior median LE (stored, same mAs)", lambda r: r["medians"]["L"], "ratio", "%.0f")
    add("Levels", "interior median HE (stored, same mAs)", lambda r: r["medians"]["H"], "ratio", "%.0f")
    for img, lab in (("LE", "LE"), ("HE", "HE"), ("S_W", "H - W L"), ("S_wstar", "H - w* L"), ("R", "R (emulated)")):
        base = emu if img == "R" else meas
        tb = {"LE": "unchanged (same spheres); patients ~3.4, indicative", "S_W": "patients ~1.1 after subtraction, indicative"}.get(img, "")
        b = add("Texture", "beta %s" % lab, (lambda img: lambda r: r["texture"][img]["beta"])(img), "diff", "%.2f", base=base, target=tb)
        p = add("Texture", "P(0.2) %s (mm2)" % lab, (lambda img: lambda r: r["texture"][img]["P02"])(img), "ratio", "%.3g", base=base,
                target="x%.2f of current (REPORT §10)" % LE_POWER_TARGET if img == "LE" else "")
        q = p["sim_redesign"] / p["sim_current"]
        if not 0.5 <= q <= 2:      # the spectrum changes composition (texture vs noise): an additive beta offset does not carry over
            b.update(predicted=None, target="not transferable (texture x%.2f); simulated redesign %.2f%s" % (q, b["sim_redesign"], "; " + tb if tb else ""))
    for img in ("LE", "HE"):
        add("Noise", "power 2.5-4.5 cy/mm, %s (same mAs)" % img, (lambda img: lambda r: band_mean(r, img, 2.5, 4.5))(img), "ratio", "%.3g")
    # rim: ring median - interior vs distance from the curved edge (pipeline_compare.rim_profile)
    c_e = 0.5 * (np.array(PC.RIM_BINS[:-1]) + PC.RIM_BINS[1:])
    for nm, lab, rm, fmt in (("S_W", "H - W L", rim_meas, "%+.4f"), ("S_wstar", "H - w* L", rim_meas, "%+.4f"), ("R", "R (emulated)", rim_emu, "%+.0f")):
        P = {k: np.array([np.nan if v is None else v for v in p[nm]], float) for k, p in (("m", rm), ("c", rim_cur), ("n", rim_new))}
        pred = P["m"] + (P["n"] - P["c"])
        pk ={k: (float(np.nanmax(v)), float(c_e[int(np.nanargmax(v))])) for k, v in
              (("measured", P["m"]), ("sim_current", P["c"]), ("sim_redesign", P["n"]), ("predicted", pred))}
        rows.append(dict(group="Rim", quantity="%s peak height" % lab, measured=pk["measured"][0], sim_current=pk["sim_current"][0],
                         sim_redesign=pk["sim_redesign"][0], mode="diff (profile)", predicted=pk["predicted"][0], target="", fmt=fmt))
        rows.append(dict(group="Rim", quantity="%s peak position (mm)" % lab, measured=pk["measured"][1], sim_current=pk["sim_current"][1],
                         sim_redesign=pk["sim_redesign"][1], mode="diff (profile)", predicted=pk["predicted"][1], target="", fmt="%.2f"))
        for x in RIM_SHOW:
            j = int(np.argmin(np.abs(c_e - x)))
            rows.append(dict(group="Rim", quantity="%s at %.1f mm" % (lab, c_e[j]), measured=P["m"][j], sim_current=P["c"][j],
                             sim_redesign=P["n"][j], mode="diff", predicted=pred[j], target="", fmt=fmt))
    # discs
    discs = {}
    for img, lab in (("R", "R (emulated)"), ("S_W", "H - W L"), ("S_wstar", "H - w* L"), ("split", "split weight")):
        base = emu if img == "R" else meas
        qc, qn, qn_iqr = disc_scaling(cur, new, img)
        sd = {k: v[0] / v[1] for k, v in disc_table(base, img).items()}
        pred = {k: v * qc / qn for k, v in sd.items()}
        cnt = lambda t, th: sum(1 for v in t.values() if v >= th)
        simsd = {s: {k: v[0] / v[1] for k, v in disc_table(r, img).items()} for s, r in (("c", cur), ("n", new))}
        for th in (3, 5):
            rows.append(dict(group="Discs (of 35)", quantity="SDNR >= %d, %s" % (th, lab), measured=cnt(sd, th), sim_current=cnt(simsd["c"], th),
                             sim_redesign=cnt(simsd["n"], th), mode="scaled SDNR", predicted=cnt(pred, th), target="", fmt="%d"))
        rows.append(dict(group="Discs (of 35)", quantity="contrast / noise ratio, %s" % lab, measured=None, sim_current=None, sim_redesign=None,
                         mode="redesign/current", predicted=None, target="", fmt="%.3f", note="contrast x%.3f, noise x%.3f (IQR %.3f-%.3f)" % (qc, qn, *qn_iqr)))
        discs[img] = dict(contrast_ratio=qc, noise_ratio=qn, noise_ratio_iqr=qn_iqr, measured_sdnr={"%d,%d" % k: v for k, v in sd.items()},
                          predicted_sdnr={"%d,%d" % k: v for k, v in pred.items()})
    return rows, discs


# ------------------------------------------------------------------------------------------------------ output
def cell(v, fmt):
    return "-" if v is None or (isinstance(v, float) and not np.isfinite(v)) else fmt % v


def change(r):
    c, n = r["sim_current"], r["sim_redesign"]
    if r.get("note"):
        return r["note"]
    if c is None or n is None or not np.isfinite(c) or not np.isfinite(n):
        return "-"
    if r["mode"] == "ratio":
        return "x%.3f (%+.1f %%)" % (n / c, 100 * (n / c - 1)) if c else "-"
    return ("%+" + r["fmt"].lstrip("%+")) % (n - c)                     # signed difference in the row's format


def figure(meas, emu, cur, new, rim_meas, rim_emu, rim_cur, rim_new, discs, out, title, ref=None):
    """ref: dict(sim, rim, discs) of the tissue reference, drawn as its prediction (same method) where given"""
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    col = dict(m="#1f77b4", c="#d62728", n="#2ca02c", p="#000000", t="#9467bd")
    RL = "predicted tissue reference"
    fig, ax = plt.subplots(2, 3, figsize=(18, 9.5))
    W, kv = meas["meta"]["W"], meas["meta"].get("kv")
    # weights vs frequency
    fc = np.sqrt([b["f1"] * b["f2"] for b in meas["cross"]])
    y = {k: np.array([b["w_fwd"] for b in r["cross"]]) for k, r in (("m", meas), ("c", cur), ("n", new))}
    for k, lab in (("m", "measured, current"), ("c", "simulated, current"), ("n", "simulated, redesign")):
        ax[0, 0].plot(fc, y[k], marker="o", ms=3, color=col[k], label=lab, ls="-" if k == "m" else "--")
    ax[0, 0].plot(fc, y["m"] * y["n"] / y["c"], "-o", ms=3, color=col["p"], label="predicted redesign")
    if ref:
        ax[0, 0].plot(fc, y["m"] * np.array([b["w_fwd"] for b in ref["sim"]["cross"]]) / y["c"], "-.o", ms=3, color=col["t"], label=RL)
    ax[0, 0].axhline(W, color="k", ls=":", lw=0.8)
    if kv in R_TISSUE:
        ax[0, 0].axhline(W / (1 - R_TISSUE[kv]), color="0.5", ls="-.", lw=0.8, label="target w* (r_tissue %+.2f)" % R_TISSUE[kv])
    ax[0, 0].set_xscale("log"); ax[0, 0].set_xlabel("cy/mm"); ax[0, 0].set_title("structure-cancelling weight w*(f); dotted: header W", fontsize=10)
    # texture spectra
    f = np.array(meas["fc"])
    for a_, img, base in ((ax[0, 1], "LE", meas), (ax[0, 2], "R", emu)):
        sp = {k: np.array([np.nan if v is None else v for v in r["texture"][img]["rad"]], float) for k, r in (("m", base), ("c", cur), ("n", new))}
        for k, lab in (("m", "measured, current"), ("c", "simulated, current"), ("n", "simulated, redesign")):
            a_.loglog(f, sp[k], color=col[k], label=lab, ls="-" if k == "m" else "--")
        a_.loglog(f, sp["m"] * sp["n"] / sp["c"], color=col["p"], label="predicted redesign")
        if ref:
            st = np.array([np.nan if v is None else v for v in ref["sim"]["texture"][img]["rad"]], float)
            a_.loglog(f, sp["m"] * st / sp["c"], "-.", color=col["t"], label=RL)
        a_.axvspan(0.2, 0.7, color="0.92"); a_.set_xlabel("cy/mm")
        a_.set_title("texture spectrum %s (pipe.py texture())" % ("LE" if img == "LE" else "R, emulated at W"), fontsize=10)
    # rim profiles
    c_e = 0.5 * (np.array(PC.RIM_BINS[:-1]) + PC.RIM_BINS[1:])
    for a_, nm, rm, lab in ((ax[1, 0], "S_W", rim_meas, "H - W L"), (ax[1, 1], "R", rim_emu, "R (emulated)")):
        P = {k: np.array([np.nan if v is None else v for v in p[nm]], float) for k, p in (("m", rm), ("c", rim_cur), ("n", rim_new))}
        for k, l in (("m", "measured, current"), ("c", "simulated, current"), ("n", "simulated, redesign")):
            a_.plot(c_e, P[k], marker="o", ms=3, color=col[k], label=l, ls="-" if k == "m" else "--")
        a_.plot(c_e, P["m"] + P["n"] - P["c"], "-o", ms=3, color=col["p"], label="predicted redesign")
        if ref:
            Pt = np.array([np.nan if v is None else v for v in ref["rim"][nm]], float)
            a_.plot(c_e, P["m"] + Pt - P["c"], "-.o", ms=3, color=col["t"], label=RL)
        a_.axhline(0, color="k", lw=0.6); a_.set_xlim(-4.5, 50); a_.set_xlabel("mm from the curved edge")
        a_.set_title("rim: %s, ring median - interior" % lab, fontsize=10)
    # discs: median SDNR by diameter, measured and predicted
    for img, ls in (("R", "-"), ("S_wstar", "--")):
        sets = [(discs, "measured_sdnr", col["m"], "measured"), (discs, "predicted_sdnr", col["p"], "predicted")]
        sets += [(ref["discs"], "predicted_sdnr", col["t"], "tissue reference")] if ref else []
        for dd, key, c_, lab in sets:
            A = np.full((7, 5), np.nan)
            for k, v in dd[img][key].items():
                r_, c2 = map(int, k.split(",")); A[r_, c2] = v
            ax[1, 2].plot([2, 3, 4, 5, 7], np.nanmedian(A, 0), marker="o", ms=3, ls=ls, color=c_, label="%s, %s" % (lab, "R" if img == "R" else "H - w* L"))
    ax[1, 2].axhline(3, color="k", ls=":", lw=0.7); ax[1, 2].set_xlabel("disc diameter (mm)")
    ax[1, 2].set_title("iodine discs: median SDNR over the 7 rows", fontsize=10)
    for a_ in ax.ravel():
        a_.legend(fontsize=7)
    fig.suptitle(title, fontsize=11)
    fig.tight_layout(); fig.savefig(out / "prediction.png", dpi=100); plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    add_case_args(ap, background=False)
    ap.add_argument("--current", choices=list(BACKGROUNDS), default="PMMA_paraffin", help="background of the measured phantom")
    ap.add_argument("--redesign", choices=list(BACKGROUNDS), default="PA12_MCT", help="background to predict")
    ap.add_argument("--reference", choices=list(BACKGROUNDS), default="gland_adipose",
                    help="tissue reference background: when its twin exists, its prediction (same method) is the target")
    a = ap.parse_args()
    apply_case_args(a)
    pid = MEASURED["pair"]
    root_c, case_c = case_root(a.current)
    root_n, case_n = case_root(a.redesign)
    cur_all, rim_c = load_case(root_c, ("measured", "measured_emuR", "simulated"), pid)
    new_all, rim_n = load_case(root_n, ("simulated",), pid)
    meas, emu, cur, new = cur_all["measured"], cur_all["measured_emuR"], cur_all["simulated"], new_all["simulated"]
    rims = (rim_c["measured"], rim_c["measured_emuR"], rim_c["simulated"])
    msgs, shift = checks(cur, new)
    rows, discs = predict(meas, emu, cur, new, *rims, rim_n["simulated"])
    # tissue reference (glandular spheres in adipose tissue), predicted from the measured current L2 in the same way
    ref = None
    if a.reference not in (a.current, a.redesign):
        root_r, case_r = case_root(a.reference)
        if (root_r / "simulated" / "results" / ("%s.json" % pid)).exists() and (root_r / "report" / pid / "report.json").exists():
            ref_all, rim_r = load_case(root_r, ("simulated",), pid)
            ref = dict(sim=ref_all["simulated"], rim=rim_r["simulated"], case=case_r)
            msgs += ["reference: " + m for m in checks(cur, ref["sim"])[0]]
            rows_r, ref["discs"] = predict(meas, emu, cur, ref["sim"], *rims, ref["rim"])
            for r, rr in zip(rows, rows_r):
                r["reference"] = rr["predicted"]
        else:
            print("tissue reference %s not simulated yet: targets from the REPORT §10 model only" % case_r)
    for m in msgs:
        print("CHECK:", m)
    out = OUT_DIR / "redesign" / ("%s_%s_%s_vs_%s" % (MEASURED["set"], pid, a.current, a.redesign))
    out.mkdir(parents=True, exist_ok=True)
    rc = ["Predicted tissue reference (%s)" % a.reference] if ref else []
    head = ["| " + " | ".join(["Group", "Quantity", "Measured (current)", "Simulated current", "Simulated redesign", "Simulated change",
                               "Predicted redesign"] + rc + ["Target (REPORT §10 model)"]) + " |", "|---" * (8 + len(rc)) + "|"]
    lines = head + ["| " + " | ".join([r["group"], r["quantity"], cell(r["measured"], r["fmt"]), cell(r["sim_current"], r["fmt"]),
                                       cell(r["sim_redesign"], r["fmt"]), change(r), cell(r["predicted"], r["fmt"])]
                                      + ([cell(r.get("reference"), r["fmt"])] if ref else []) + [r["target"]]) + " |" for r in rows]
    verdict = ""
    if ref:
        # the redesign against the tissue reference, both predicted from the measured current L2 in the same way
        key = ("w* (0.06-0.15", "r = 1 - W/w*", "interior median LE", "interior median HE", "beta LE", "P(0.2) LE", "P(0.2) HE",
               "P(0.2) H - W L", "P(0.2) R", "H - W L peak height", "R (emulated) peak height", "SDNR >= 3, R", "SDNR >= 3, H - w* L")
        vl = []
        for r in rows:
            p, t = r["predicted"], r.get("reference")
            if r["quantity"].startswith(key) and p is not None and t is not None:
                d = ("x%.2f" % (p / t)) if r["mode"] == "ratio" and t else (("%+" + r["fmt"].lstrip("%+")) % (p - t))
                vl.append("- %s: redesign %s, tissue reference %s (redesign vs reference %s)" % (r["quantity"], cell(p, r["fmt"]), cell(t, r["fmt"]), d))
        verdict = ("Redesign against the tissue reference (glandular spheres in adipose tissue, Hammerstein; 61 %% of the volume, "
                   "≈ 64 %% glandular by mass), both predicted from the measured current L2 in the same way:\n%s\n\n" % "\n".join(vl))
    intro = ("Measured pair %s (set %s, current background %s) and its simulated twins %s and %s%s; single exposure, no averaging.\n\n"
             "Predicted redesign = measured x (simulated redesign / simulated current) for positive quantities, measured + "
             "(simulated redesign - simulated current) for beta and rim heights; r from the predicted w* with the header W. "
             "R = Siemens emulation at the header W (measured: measured_emuR). Disc SDNR: measured SDNR x contrast ratio / noise ratio "
             "(redesign / current simulation). Model errors common to both simulations cancel; errors that depend on the materials "
             "do not [A]. Same mAs as the measured pair (no AEC).\n\n"
             "Checks: %s\n\n" % (pid, MEASURED["set"], a.current, case_c, case_n, " (tissue reference %s)" % ref["case"] if ref else "",
                                 "; ".join(msgs) if msgs else "twins identical except the materials; disc registration shift %.1f px" % shift))
    (out / "summary.md").write_text(intro + verdict + "\n".join(lines) + "\n")
    json.dump(dict(pair=pid, set=MEASURED["set"], current=case_c, redesign=case_n, reference=ref["case"] if ref else None, checks=msgs,
                   registration_shift_px=shift, rows=rows, discs=discs, reference_discs=ref["discs"] if ref else None),
              open(out / "prediction.json", "w"), indent=1, default=float)
    figure(meas, emu, cur, new, *rims, rim_n["simulated"], discs, out,
           "Pair %s: measured current L2 (%s) and predicted redesign (%s)%s; single exposure"
           % (pid, a.current, a.redesign, " with the tissue reference" if ref else ""), ref=ref)
    print(verdict + "\n".join(lines))
    print("written", out)


if __name__ == "__main__":
    main()

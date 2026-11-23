#!/usr/bin/env python3
"""
The investigation's own analysis and recombination, applied unchanged to the simulated twins of the measured pairs.

The July analysis (investigation_0035_0037/scripts/july2026/pipe.py, REPORT §13) reads DICOM pairs from an exposure
list and computes, per exposure: weights vs spatial frequency (w*(f), coherence) and Siemens' effective weight from
regressing the recombined image R on H and L; the local weight vs distance from the curved edge (128 px tiles);
texture spectra, beta and P(0.2) of LE, HE, H - W L, H - w* L and R; disc contrast and SDNR in six images (LE, HE,
S_W, S_w*, split weight, R); the Siemens emulation (REPORT §12). agg.py and discagg.py summarise them.

This script makes the simulated twins look exactly like the measured input, so that the same files run on both:
  measured/        the measured pair as it is (LE, HE and Siemens' recombined image R)
  measured_emuR/   the same pair with R replaced by the emulation (pipe.emulate at the header weight, plus the
                   offset and the non-zero mask of the real R). Separates what the emulation cannot reproduce
                   (Siemens' border processing in the outer ~3 mm) from what differs between phantom and simulation
  simulated/       the twins from measlike.py (each with its pair's exact LE / HE mAs and HE open level), written as
                   DICOM in the raw orientation (R-CC: pipe.py mirrors them back), with R = the same emulation, using
                   the offset and mask of the measured twin
Disc registration of the twins: pipe.py places its disc template by the phantom outline plus the summed response of
H - 0.37 L in 5 x 5 px at the disc centres. On the twins this estimate is off (6 px in y for PMMA/paraffin; 14 px for
PA12/MCT, whose spheres do not cancel at 0.37), while their true disc positions are known from the simulation (centroids
of the projector's iodine map, agreeing with the template layout within 0.2 px). The twins are therefore analysed with
that known registration: pipe.py unchanged, only its register() replaced in its namespace (run_twins). The measured
pair keeps pipe.py's own estimate (an independent disc-shaped locator agrees with it within 1 px). twins.json records
both values.
Each workspace (output/pipeline/<set>_<case>/<kind>/) holds links to pipe.py, agg.py, discagg.py, offset.py and
template.json of july2026/ and its own exposures.json; pipe.py writes results/<id>.json there, agg.py / discagg.py
write agg.json, discagg.json (and their printed tables: agg.txt, discagg.txt).

The report adds one analysis that pipe.py does not do, applied identically to all three: the rim profile, i.e. the
median of S_W = H - W L, S_w* = H - w* L and R in rings of distance from the fitted curved edge (and from the chest
wall), relative to the interior. This is the ring-median method of scripts/shading.py (REPORT §4, Fig. 6), with the
0035 pixel coordinates replaced by pipe.py's own geometry fit and disc registration.

One measured pair per run (--pair, default HIGHMAS_28KV2_00), compared with its own simulated twin; each pair is
treated on its own and nothing is averaged over pairs (another pair = another run with its --pair).
Usage:
  python pipeline_compare.py prepare [--phantom L2] [--pair ID]   write the workspaces (needs the measlike twin)
  python pipeline_compare.py run [--workers 3]                     pipe.py, agg.py, discagg.py in each workspace
  python pipeline_compare.py report [--pair ID]                    tables and figures -> .../report/<pair>/
  python pipeline_compare.py all [--pair ID]                       the three in order
"""
import argparse, json, os, subprocess, sys, time
from pathlib import Path
import numpy as np
from scipy import ndimage as ndi
from config import OUT_DIR, PROJECT, MEASURED, MEASURED_SETS, PHANTOM, add_case_args, apply_case_args
import measured_set as MS

INV = PROJECT / "investigation_0035_0037" / "scripts" / "july2026"
LINKED = ("pipe.py", "agg.py", "discagg.py", "offset.py", "template.json")
KINDS = ("measured", "measured_emuR", "simulated")
LABEL = {"measured": "measured (Siemens R)", "measured_emuR": "measured (emulated R)", "simulated": "simulated (emulated R)"}
COLOR = {"measured": "#1f77b4", "measured_emuR": "#7f7f7f", "simulated": "#d62728"}
RIM_BINS = [-4, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7.5, 10, 12.5, 15, 20, 25, 30, 45, 70]
CW_BINS = [-2, 0, 1, 2, 3, 5, 8, 12, 16, 20, 30]


def root():
    return OUT_DIR / "pipeline" / ("%s_%s" % (MEASURED["set"], __import__("config").case_name()))


def pipe_ns(ws):
    """the definitions of pipe.py (everything before process()), executed as if pipe.py were in workspace ws"""
    ns = {"__file__": str(ws / "pipe.py"), "__name__": "pipe_definitions"}
    exec((INV / "pipe.py").read_text().split("def process(e):")[0], ns)
    return ns


def measured_entries(ids=None):
    src = MEASURED_SETS[MEASURED["set"]]["source"]
    if src.get("kind") != "list":
        raise SystemExit("set %s has no exposure list: the July pipeline needs july2026-style exposures" % MEASURED["set"])
    ex = json.load(open(src["exposures"]))
    sel = [dict(e) for e in ex if e["kv"] == src["kvp"] and e["level"] == src["level"]
           and e["unit"] == MEASURED_SETS[MEASURED["set"]]["unit"] and (not ids or e["id"] in ids)]
    for e in sel:
        for k in ("L", "H", "R"):
            e[k] = MS.locate(e[k])
    return sel


def make_workspace(ws, entries):
    ws.mkdir(parents=True, exist_ok=True)
    for f in LINKED:
        link = ws / f
        if link.is_symlink() or link.exists():
            link.unlink()
        link.symlink_to(INV / f)
    json.dump(entries, open(ws / "exposures.json", "w"), indent=1)


def write_dicom(path, img, desc, e, beam):
    """minimal mammography DICOM (uint16, stored-value image) that pydicom-based analysis reads like the real files"""
    import pydicom
    from pydicom.dataset import Dataset, FileMetaDataset
    from pydicom.uid import ExplicitVRLittleEndian, generate_uid
    meta = FileMetaDataset()
    meta.MediaStorageSOPClassUID = "1.2.840.10008.5.1.4.1.1.1.2.1"          # Digital Mammography X-Ray, for processing
    meta.MediaStorageSOPInstanceUID = generate_uid()
    meta.TransferSyntaxUID = ExplicitVRLittleEndian
    ds = Dataset(); ds.file_meta = meta
    ds.SOPClassUID, ds.SOPInstanceUID = meta.MediaStorageSOPClassUID, meta.MediaStorageSOPInstanceUID
    ds.Modality = "MG"; ds.SeriesDescription = desc; ds.PatientName = "SIMULATION^GATE10"; ds.PatientID = "SIMULATION"
    ds.ImageComments = "not a patient image: %s for measured pair %s" % (desc, e["id"])
    ds.Rows, ds.Columns = img.shape
    ds.SamplesPerPixel = 1; ds.PhotometricInterpretation = "MONOCHROME2"
    ds.BitsAllocated = 16; ds.BitsStored = 16; ds.HighBit = 15; ds.PixelRepresentation = 0
    ds.ImagerPixelSpacing = [0.085, 0.085]
    if beam in ("LE", "HE"):
        ds.KVP = 28.0 if beam == "LE" else 49.0
        ds.ExposureInuAs = int(round(1000 * (e["mAsL"] if beam == "LE" else e["mAsH"])))
    ds.PixelData = np.ascontiguousarray(img, dtype="<u2").tobytes()
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        ds.save_as(str(path), enforce_file_format=True)                       # pydicom 3
    except TypeError:
        ds.is_little_endian, ds.is_implicit_VR = True, False                   # pydicom 2
        ds.save_as(str(path), write_like_original=False)


def grown(fp, k, cache):
    """footprint grown (k > 0) or shrunk (k < 0) by k pixels"""
    if "out" not in cache:
        cache["out"], cache["in"] = ndi.distance_transform_edt(~fp), ndi.distance_transform_edt(fp)
    return (cache["out"] <= k) if k >= 0 else (cache["in"] > -k)


def true_registration(ns, fp, cw, cy):
    """disc registration of a simulated twin from the simulation itself: centroids of the projector's iodine map
    (output/projector/<case>/maps.npz, same grid and orientation as the twin) matched to the template's array discs;
    ([dy, dx], largest deviation of a single disc in px), or None if the map is missing or does not fit the template"""
    f = OUT_DIR / "projector" / __import__("config").case_name() / "maps.npz"
    if not f.exists():
        return None
    I = np.load(f)["iodine_g_cm2"]
    if I.shape != fp.shape:
        return None
    lab, n = ndi.label(I > 0)
    cen = np.array(ndi.center_of_mass(I > 0, lab, range(1, n + 1)))
    T = ns["T"]; dy0, dx0 = cy - T["ref"][1], cw - T["ref"][0]
    off = np.array([cen[int(np.argmin(np.hypot(*(cen - (d["y"] + dy0, d["x"] + dx0)).T)))] - (d["y"], d["x"])
                    for d in T["discs"] if d["kind"] == "array"])
    med = np.median(off, 0); dev = float(np.abs(off - med).max())
    return ([float(med[0]), float(med[1])], dev) if dev <= 2 else None


def emulated_R(ns, lL, lH, fp, W, c, k):
    """R as Siemens stores it, from the emulation: offset c, zero outside the footprint grown by k px, 12 bit"""
    E = ns["emulate"](lL, lH, fp, W)
    return np.where(grown(fp, k, {}), np.clip(np.rint(E + c), 0, 4095), 0).astype(np.uint16)


def prepare(a):
    meas = measured_entries([MEASURED["pair"]])
    if not meas:
        raise SystemExit("no measured pairs selected")
    R0 = root(); R0.mkdir(parents=True, exist_ok=True)
    make_workspace(R0 / "measured", meas)
    ns = pipe_ns(R0 / "measured")
    sim_dir = OUT_DIR / "measlike" / __import__("config").case_name() / MEASURED["set"]
    twins, emu_entries, sim_entries = {}, [], []
    for e in meas:
        t0 = time.time()
        L, H, R = ns["load"](e)                                                 # analysis frame (mirrored)
        lL, lH = -np.log(np.clip(L, 1, None)), -np.log(np.clip(H, 1, None))
        fp = ns["footprint"](H)
        E = ns["emulate"](lL, lH, fp, e["W"])
        dist = ndi.distance_transform_edt(fp) * ns["PX"]
        inner = (R > 0) & (R < 4095) & fp & (dist > 15)
        c = float(np.median(R[inner] - E[inner]))                               # as pipe.py / fig19_20.py
        cache, mR = {}, R > 0
        iou = [(float((grown(fp, k, cache) & mR).sum() / (grown(fp, k, cache) | mR).sum()), k) for k in range(-20, 81, 2)]
        best_iou, k = max(iou)
        twins[e["id"]] = dict(offset_c=c, mask_growth_px=k, mask_iou=best_iou, saturated_px=int((R >= 4095).sum()))
        # measured pair with the emulated R
        Rm = np.where(grown(fp, k, cache), np.clip(np.rint(E + c), 0, 4095), 0).astype(np.uint16)
        pR = R0 / "measured_emuR" / "dicom" / e["id"] / "R_emulated.dcm"
        write_dicom(pR, Rm[:, ::-1], "emulated Siemens recombination of the measured pair", e, "R")
        emu_entries.append(dict(e, R=str(pR), R_source="emulated"))
        # simulated twin
        d = sim_dir / e["id"]
        if (d / "LE.npy").exists():
            Ls, Hs = np.load(d / "LE.npy"), np.load(d / "HE.npy")               # simulation frame = analysis frame
            lLs, lHs = -np.log(np.clip(Ls.astype(float), 1, None)), -np.log(np.clip(Hs.astype(float), 1, None))
            fps = ns["footprint"](Hs.astype(float))
            Rs = emulated_R(ns, lLs, lHs, fps, e["W"], c, k)
            cws, _, cys, _ = ns["geom"](fps)
            S0 = lHs - 0.37 * lLs                                               # pipe.py's own estimate, for the record
            est = [float(v) for v in ns["register"](ndi.gaussian_filter(S0, 2.5) - ndi.gaussian_filter(S0, 40), cws, cys)]
            reg = true_registration(ns, fps, cws, cys)
            twins[e["id"]].update(registration_pipe=est, registration_true=reg[0] if reg else None,
                                  registration_true_max_dev_px=reg[1] if reg else None)
            print("  twin disc registration (dy, dx): pipe.py estimate (%.2f, %.2f), from the simulation %s"
                  % (*est, "(%.2f, %.2f)" % tuple(reg[0]) if reg else "not available: pipe.py estimate used"))
            out = R0 / "simulated" / "dicom" / e["id"]
            for nm, img, beam in (("LE.dcm", Ls, "LE"), ("HE.dcm", Hs, "HE"), ("R_emulated.dcm", Rs, "R")):
                write_dicom(out / nm, img[:, ::-1], "GATE 10 simulation (%s), twin of %s" % (nm[:-4], e["id"]), e, beam)
            sim_entries.append(dict(e, L=str(out / "LE.dcm"), H=str(out / "HE.dcm"), R=str(out / "R_emulated.dcm"),
                                    R_source="emulated", simulated=True, registration_from_simulation=reg[0] if reg else None,
                                    simulation=json.load(open(d / "meta.json")) if (d / "meta.json").exists() else None))
        else:
            print("  no simulated twin for %s in %s (run measlike.py first)" % (e["id"], sim_dir))
        print("  %s: R offset %.0f, mask = footprint + %d px (IoU %.3f), %d saturated px; %.0f s"
              % (e["id"], c, k, best_iou, twins[e["id"]]["saturated_px"], time.time() - t0), flush=True)
    make_workspace(R0 / "measured_emuR", emu_entries)
    if sim_entries:
        make_workspace(R0 / "simulated", sim_entries)
    json.dump(dict(set=MEASURED["set"], case=__import__("config").case_name(), phantom=PHANTOM["model"], twins=twins,
                   simulated_ids=[s["id"] for s in sim_entries]), open(R0 / "twins.json", "w"), indent=1)
    print("workspaces in", R0)


def run(a):
    R0 = root()
    for kind in KINDS:
        ws = R0 / kind
        if not (ws / "exposures.json").exists():
            print("skipping %s (not prepared)" % kind); continue
        t0 = time.time()
        if kind == "simulated":
            run_twins(ws)
        else:
            procs = [subprocess.Popen([sys.executable, "pipe.py", str(k), str(a.workers)], cwd=ws) for k in range(a.workers)]
            for p in procs:
                p.wait()
        n = len(list((ws / "results").glob("*.json")))
        print("%s: %d exposures analysed (%.0f s)" % (kind, n, time.time() - t0), flush=True)
        for script in ("agg.py", "discagg.py"):
            r = subprocess.run([sys.executable, script], cwd=ws, capture_output=True, text=True)
            (ws / script.replace(".py", ".txt")).write_text(r.stdout + r.stderr)
            if r.returncode:
                print("  %s failed in %s:\n%s" % (script, kind, r.stderr[-2000:]))


def run_twins(ws):
    """pipe.py on the simulated twins, unchanged except that register() returns each twin's disc registration known from
    the simulation (exposure entry 'registration_from_simulation'; pipe.py's own estimate where it is missing)"""
    ns = {"__file__": str(ws / "pipe.py"), "__name__": "pipe_twin"}
    exec((INV / "pipe.py").read_text(), ns)
    estimate = ns["register"]
    for e in ns["E"]:
        reg = e.get("registration_from_simulation")
        ns["register"] = (lambda r: lambda Sb, cw, cy: (r[0], r[1]))(reg) if reg else estimate
        s = ns["process"](e)
        print("  %s done %.0f s (disc registration: %s)" % (e["id"], s, "from the simulation" if reg else "pipe.py estimate"), flush=True)


# ------------------------------------------------------------------------------------------------------ report
def results(kind):
    ws = root() / kind
    rs = [json.load(open(f)) for f in sorted((ws / "results").glob("*.json"))] if (ws / "results").exists() else []
    return {r["id"]: r for r in rs}


def rim_profile(ns, e, w_star):
    """ring medians of S_W, S_w* and R vs distance from the curved edge (and from the chest wall), relative to the
    interior; the shading.py method with pipe.py's geometry and disc exclusion"""
    PX = ns["PX"]
    L, H, R = ns["load"](e)
    lL, lH = -np.log(np.clip(L, 1, None)), -np.log(np.clip(H, 1, None))
    fp = ns["footprint"](H); cw, cx, cy, rad = ns["geom"](fp)
    yy, xx = np.mgrid[0:L.shape[0], 0:L.shape[1]]
    de = (rad - np.hypot(xx - cx, yy - cy)) * PX; dcw = (xx - cw) * PX
    S0 = lH - 0.37 * lL; Sb = ndi.gaussian_filter(S0, 2.5) - ndi.gaussian_filter(S0, 40)       # as pipe.process
    dy, dx = ns["register"](Sb, cw, cy)
    excl = np.zeros(L.shape, bool)
    for d in ns["T"]["discs"]:
        excl |= ((yy - d["y"] - dy) ** 2 + (xx - d["x"] - dx) ** 2) < (d["radius_px"] + (60 if d["kind"] == "array" else 90)) ** 2
    excl |= ((yy - ns["T"]["marker"][0] - dy) ** 2 + (xx - ns["T"]["marker"][1] - dx) ** 2) < 120 ** 2
    interior = fp & (de > 25) & (dcw > 25) & ~excl
    W = e["W"]
    imgs = {"S_W": lH - W * lL, "S_wstar": lH - w_star * lL, "R": R}
    out = dict(edge_bins=RIM_BINS, cw_bins=CW_BINS)
    for nm, im in imgs.items():
        ref = float(np.median(im[interior]))
        ok = (R > 0) if nm == "R" else np.ones(im.shape, bool)
        out[nm] = [float(np.median(im[s]) - ref) if s.any() else None for s in
                   (ok & (de >= lo) & (de < hi) & (dcw > 20) & ~excl for lo, hi in zip(RIM_BINS[:-1], RIM_BINS[1:]))]
        out[nm + "_cw"] = [float(np.median(im[s]) - ref) if s.any() else None for s in
                           (ok & (dcw >= lo) & (dcw < hi) & (np.abs(yy - cy) < 600) & (de > 30) & ~excl
                            for lo, hi in zip(CW_BINS[:-1], CW_BINS[1:]))]
    out["R_saturated"] = [float((R[s] >= 4095).mean()) if s.any() else None for s in
                          ((R > 0) & (de >= lo) & (de < hi) & (dcw > 20) for lo, hi in zip(RIM_BINS[:-1], RIM_BINS[1:]))]
    return out


def report(a):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    R0 = root(); pid = MEASURED["pair"]
    res = {k: {pid: r} for k in KINDS for r in [results(k).get(pid)] if r is not None}
    kinds = [k for k in KINDS if k in res]
    if "measured" not in kinds:
        raise SystemExit("no results for pair %s yet: run  python pipeline_compare.py prepare/run --pair %s" % (pid, pid))
    ids = [pid]
    out = R0 / "report" / pid; out.mkdir(parents=True, exist_ok=True)
    print("pair %s analysed in: %s" % (pid, ", ".join(kinds)))
    summ = {}
    # ---- the values of this pair (same quantities as agg.py / discagg.py, for one exposure)
    rows = [("w* (0.06-0.15 cy/mm)", lambda r: r["w_low"], "%.4f"),
            ("Siemens effective w (R on H, L)", lambda r: np.mean([b["w_eff"] for b in r["cross"][:2]]), "%.4f"),
            ("r = 1 - W/w*", lambda r: 1 - r["meta"]["W"] / r["w_low"], "%+.3f"),
            ("disc registration dy (px)", lambda r: r["geom"]["reg_dy"], "%.2f"),
            ("disc registration dx (px)", lambda r: r["geom"]["reg_dx"], "%.2f"),
            ("interior median LE (stored)", lambda r: r["medians"]["L"], "%.0f"),
            ("interior median HE (stored)", lambda r: r["medians"]["H"], "%.0f")]
    rows += [("beta %s" % nm, (lambda nm: lambda r: r["texture"][nm]["beta"])(nm), "%.2f") for nm in ("LE", "HE", "S_W", "S_wstar", "R")]
    rows += [("P(0.2) %s (mm2)" % nm, (lambda nm: lambda r: r["texture"][nm]["P02"])(nm), "%.2e") for nm in ("LE", "HE", "S_W", "S_wstar", "R")]
    # emulation vs the real Siemens image: only meaningful where R is Siemens' own
    rows += [("emulation corr, interior (real R)", lambda r: r["emu"]["interior"]["corr"] if r["meta"].get("R_source") != "emulated" else None, "%.3f"),
             ("emulation corr, outer 15 mm (real R)", lambda r: r["emu"]["edge"]["corr"] if r["meta"].get("R_source") != "emulated" else None, "%.3f")]
    lines = ["%-38s" % "" + "".join("%24s" % LABEL[k] for k in kinds)]
    for lab, fn, f in rows:
        cells = []
        for k in kinds:
            v = fn(res[k][pid])
            v = None if v is None or not np.isfinite(v) else float(v)
            summ.setdefault(k, {})[lab] = v
            cells.append("%24s" % ("-" if v is None else f % v))
        lines.append("%-38s" % lab + "".join(cells))
    # discs: number of array discs with SDNR >= 3 / >= 5 (discagg.py definition), per image
    for img in ("R", "S_W", "S_wstar", "split"):
        for th in (3, 5):
            cells = []
            for k in kinds:
                n = sum(1 for d in res[k][pid]["discs"] if d["kind"] == "array" and d.get(img) and d[img][1] > 0
                        and d[img][0] / d[img][1] >= th)
                summ.setdefault(k, {})["discs SDNR>=%d %s" % (th, img)] = n
                cells.append("%24d" % n)
            lines.append("%-38s" % ("discs with SDNR >= %d in %s (of 35)" % (th, img)) + "".join(cells))
    txt = "\n".join(lines)
    print(txt)
    (out / "summary.txt").write_text("Set %s, case %s: measured pair %s and its simulated twin (single exposure, no averaging).\n\n%s\n"
                                     % (MEASURED["set"], __import__("config").case_name(), pid, txt))
    # ---- rim profiles (images needed)
    rim = {}
    for k in kinds:
        ns = pipe_ns(root() / k)
        E = {e["id"]: e for e in json.load(open(root() / k / "exposures.json"))}
        rim[k] = [rim_profile(ns, E[i], res[k][i]["w_low"]) for i in ids]
        print("  rim profiles %s done" % k, flush=True)
    json.dump(dict(summary=summ, rim=rim, ids=ids), open(out / "report.json", "w"), indent=1, default=float)
    figures(res, rim, kinds, ids, out, plt)
    print("written", out)


def _band(ax, x, Y, k, **kw):
    Y = np.array(Y, float)
    m, s = np.nanmean(Y, 0), np.nanstd(Y, 0)
    ax.plot(x, m, "-o", ms=3, color=COLOR[k], label=LABEL[k], **kw)
    if len(Y) > 1:
        ax.fill_between(x, m - s, m + s, color=COLOR[k], alpha=0.18, lw=0)


def figures(res, rim, kinds, ids, out, plt):
    import warnings
    warnings.filterwarnings("ignore", category=RuntimeWarning, message="(Mean of empty slice|Degrees of freedom|All-NaN)")
    # 1. weights vs frequency
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.6))
    r0 = res["measured"][ids[0]]
    fc = np.sqrt([b["f1"] * b["f2"] for b in r0["cross"]])
    for k in kinds:
        _band(ax[0], fc, [[b["w_fwd"] for b in res[k][i]["cross"]] for i in ids], k)
        _band(ax[1], fc, [[b["w_eff"] / res[k][i]["meta"]["W"] for b in res[k][i]["cross"]] for i in ids], k)
        _band(ax[2], fc, [[b["coh"] for b in res[k][i]["cross"]] for i in ids], k)
    for a_, t in zip(ax, ("structure-cancelling weight w*(f)", "Siemens effective weight / W  (K(f), from R)", "coherence of H and L")):
        a_.set_xscale("log"); a_.set_xlabel("spatial frequency (cy/mm)"); a_.set_title(t, fontsize=10)
    ax[0].axhline(r0["meta"]["W"], color="k", ls=":", lw=0.8); ax[0].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(out / "1_weights_vs_frequency.png", dpi=100); plt.close(fig)
    # 2. rim: local w* vs edge distance (tiles) and ring profiles
    fig, ax = plt.subplots(2, 3, figsize=(17, 9))
    tb = [10, 15, 20, 25, 30, 40, 50, 70, 100, 140]
    for k in kinds:
        Y = []
        for i in ids:
            T = np.array(res[k][i]["tiles"])
            Y.append([np.median(T[(T[:, 0] >= lo) & (T[:, 0] < hi), 1]) if ((T[:, 0] >= lo) & (T[:, 0] < hi)).any() else np.nan
                      for lo, hi in zip(tb[:-1], tb[1:])])
        _band(ax[0, 0], 0.5 * (np.array(tb[:-1]) + tb[1:]), Y, k)
    ax[0, 0].set_title("local w* (128 px tiles, pipe.py) vs distance from the curved edge", fontsize=10)
    ax[0, 0].set_xlabel("mm from the curved edge"); ax[0, 0].legend(fontsize=8)
    eb = 0.5 * (np.array(RIM_BINS[:-1]) + RIM_BINS[1:]); cb = 0.5 * (np.array(CW_BINS[:-1]) + CW_BINS[1:])
    for j, (nm, t) in enumerate((("S_W", "H - W L"), ("S_wstar", "H - w* L"))):
        for k in kinds:
            if k == "measured_emuR":
                continue                                       # raw images identical to 'measured'
            _band(ax[0, j + 1], eb, [[np.nan if v is None else v for v in p[nm]] for p in rim[k]], k)
            _band(ax[1, j], cb, [[np.nan if v is None else v for v in p[nm + "_cw"]] for p in rim[k]], k)
        ax[0, j + 1].set_title("rim: %s, ring median - interior (curved edge)" % t, fontsize=10)
        ax[1, j].set_title("chest-wall side: %s, median - interior" % t, fontsize=10)
        ax[0, j + 1].set_xlabel("mm from the curved edge"); ax[1, j].set_xlabel("mm from the chest-wall edge")
    for k in kinds:
        _band(ax[1, 2], eb, [[np.nan if v is None else v for v in p["R"]] for p in rim[k]], k)
    ax[1, 2].set_title("rim in the recombined image R (stored units), ring median - interior", fontsize=10)
    ax[1, 2].set_xlabel("mm from the curved edge")
    for a_ in ax.ravel():
        a_.axhline(0, color="k", lw=0.6)
    for a_ in (ax[0, 1], ax[0, 2], ax[1, 2]):
        a_.set_xlim(-4.5, 50); a_.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(out / "2_rim_and_edges.png", dpi=100); plt.close(fig)
    # 3. texture spectra
    fig, ax = plt.subplots(1, 5, figsize=(24, 4.6))
    f = np.array(r0["fc"])
    for j, nm in enumerate(("LE", "HE", "S_W", "S_wstar", "R")):
        for k in kinds:
            if nm != "R" and k == "measured_emuR":
                continue
            Y = np.array([[np.nan if v is None else v for v in res[k][i]["texture"][nm]["rad"]] for i in ids])
            m, s = np.nanmean(Y, 0), np.nanstd(Y, 0)
            ax[j].loglog(f, m, color=COLOR[k], label=LABEL[k])
            ax[j].fill_between(f, np.clip(m - s, m * 0.05, None), m + s, color=COLOR[k], alpha=0.18, lw=0)
        ax[j].axvspan(0.2, 0.7, color="0.92"); ax[j].set_title(nm); ax[j].set_xlabel("cy/mm"); ax[j].legend(fontsize=7)
    ax[0].set_ylabel("power (mm$^2$), pipe.py texture()")
    fig.tight_layout(); fig.savefig(out / "3_texture_spectra.png", dpi=100); plt.close(fig)
    # 4. discs: SDNR by size and by row, per image
    fig, ax = plt.subplots(2, 4, figsize=(20, 8))
    for j, img in enumerate(("R", "S_W", "S_wstar", "split")):
        for k in kinds:
            if img != "R" and k == "measured_emuR":
                continue
            A = np.full((len(ids), 7, 5), np.nan)
            for n, i in enumerate(ids):
                for d in res[k][i]["discs"]:
                    if d["kind"] == "array" and d.get(img) and d[img][1] > 0:
                        A[n, d["row"], d["col"]] = d[img][0] / d[img][1]
            ax[0, j].plot([2, 3, 4, 5, 7], np.nanmedian(A, axis=(0, 1)), "-o", color=COLOR[k], label=LABEL[k])
            ax[1, j].plot(range(7), np.nanmedian(A, axis=(0, 2)), "-o", color=COLOR[k], label=LABEL[k])
        ax[0, j].set_title("%s: median SDNR by disc diameter" % img, fontsize=10); ax[0, j].set_xlabel("mm")
        ax[1, j].set_title("%s: median SDNR by row (0 = most iodine)" % img, fontsize=10); ax[1, j].set_xlabel("row")
        for a_ in ax[:, j]:
            a_.axhline(3, color="k", ls=":", lw=0.7); a_.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(out / "4_iodine_discs.png", dpi=100); plt.close(fig)
    # 5. images of the first pair: Siemens / emulated R, H - W L, H - w* L, one window per column
    i = ids[0]
    show = [k for k in ("measured", "measured_emuR", "simulated") if k in kinds]
    fig, ax = plt.subplots(len(show), 3, figsize=(15, 6.2 * len(show)), squeeze=False)
    imgs = {}
    for k in show:
        ns = pipe_ns(root() / k)
        e = {x["id"]: x for x in json.load(open(root() / k / "exposures.json"))}[i]
        L, H, R = ns["load"](e)
        lL, lH = -np.log(np.clip(L, 1, None)), -np.log(np.clip(H, 1, None))
        imgs[k] = dict(R=R, SW=ndi.gaussian_filter(lH - e["W"] * lL, 1.0), Sw=ndi.gaussian_filter(lH - res[k][i]["w_low"] * lL, 1.0),
                       fp=ns["footprint"](H), W=e["W"], ws=res[k][i]["w_low"])
    crop = (slice(150, 3380), slice(0, 1800))
    for col, (key, title) in enumerate((("R", "recombined R"), ("SW", "H - W L"), ("Sw", "H - w* L"))):
        ref = imgs["measured"][key][imgs["measured"]["fp"]]
        width = np.subtract(*np.percentile(ref, [99, 1]))
        for row, k in enumerate(show):
            im = imgs[k][key]; c = np.median(im[imgs[k]["fp"]])
            ax[row, col].imshow(im[crop][::3, ::3], cmap="gray", vmin=c - 0.55 * width, vmax=c + 0.55 * width)
            ax[row, col].set_xticks([]); ax[row, col].set_yticks([])
            ax[row, col].set_title("%s - %s%s" % (LABEL[k], title, " (w* %.3f)" % imgs[k]["ws"] if key == "Sw" else
                                                 (" (W %.3f)" % imgs[k]["W"] if key == "SW" else "")), fontsize=9)
    fig.suptitle("Pair %s and its simulated twin; one grey-window width per column (from the measured image)" % i, fontsize=11)
    fig.tight_layout(); fig.savefig(out / "5_images.png", dpi=70); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["prepare", "run", "report", "all"])
    add_case_args(ap)
    ap.add_argument("--workers", type=int, default=3, help="parallel pipe.py processes per workspace")
    a = ap.parse_args()
    apply_case_args(a)
    if a.step in ("prepare", "all"):
        prepare(a)
    if a.step in ("run", "all"):
        run(a)
    if a.step in ("report", "all"):
        report(a)


if __name__ == "__main__":
    main()

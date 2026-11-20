#!/usr/bin/env python3
"""
Phase-space files of one or more GATE runs (same beam and mode, different seeds) -> primary and scatter images.

Each photon at the plane is carried through the grid/detector model (grid_detector.Response) to z = 0 and binned:
  - coarse bins (default 1 mm): Monte Carlo primary and scatter (scatter = everything with UnscatteredPrimaryFlag 0,
    i.e. Compton, Rayleigh and fluorescence photons);
  - the scatter image is smoothed (Gaussian, default sigma 5 mm: scatter is smooth on cm scales, and weaker smoothing
    lets the Monte Carlo noise dominate the texture of the subtracted image) and interpolated onto the projector grid.
Units as in projector.py: detected energy per pixel per emitted photon per steradian (keV/sr). The runs emit
`events` photons uniformly into the cone of solid angle omega_cone_sr, i.e. events/omega photons per steradian.

Usage: python postprocess.py output/LE_phantom_PA12_MCT_s1 output/LE_phantom_PA12_MCT_s2 ... --name LE_PA12_MCT
       [--delete-phsp]   (all runs must have the same beam, background and phantom model)
       [--bin-mm 1.0] [--sigma-mm 5.0] [--image-bin 1] [--accumulate]
Output: output/mc/<name>.npz (coarse images, counts, normalisation) and output/mc/scatter_<name>[_binK].npy
"""
import argparse, json
import numpy as np
from scipy import ndimage as ndi
from scipy.interpolate import RegularGridInterpolator
from config import SYSTEM, OUT_DIR
from grid_detector import Response, E_GRID
from projector import pixel_grid

BRANCHES = ["KineticEnergy", "PrePosition_X", "PrePosition_Y", "PrePosition_Z",
            "PreDirection_X", "PreDirection_Y", "PreDirection_Z", "UnscatteredPrimaryFlag", "PDGCode"]


def coarse_edges(bin_mm):
    X, Y, p = pixel_grid(1)
    x0, x1 = X[0] - p / 2, X[-1] + p / 2
    y0, y1 = Y[0] - p / 2, Y[-1] + p / 2
    return np.arange(x0, x1 + bin_mm, bin_mm), np.arange(y0, y1 + bin_mm, bin_mm)


def new_acc(ex, ey):
    shape = (len(ey) - 1, len(ex) - 1)
    return {k: np.zeros(shape) for k in ("primary", "scatter", "n_primary", "n_scatter", "w2_scatter")}


def add_batch(b, resp, ex, ey, acc):
    """b: dict of phase-space branches (numpy arrays, GATE units: MeV, mm) -> adds to the accumulators"""
    e = b["KineticEnergy"] * 1000.0                                      # MeV -> keV
    uz = b["PreDirection_Z"]
    ok = (b["PDGCode"] == 22) & (uz < 0) & (e > E_GRID[0])
    pos = np.c_[b["PrePosition_X"], b["PrePosition_Y"], b["PrePosition_Z"]][ok]
    u = np.c_[b["PreDirection_X"], b["PreDirection_Y"], uz][ok]
    w, x, y = resp.weight(e[ok], pos, u)
    prim = b["UnscatteredPrimaryFlag"][ok] == 1
    for key, m in (("primary", prim), ("scatter", ~prim)):
        acc[key] += np.histogram2d(y[m], x[m], bins=(ey, ex), weights=w[m])[0]
        acc["n_" + key] += np.histogram2d(y[m], x[m], bins=(ey, ex))[0]
    acc["w2_scatter"] += np.histogram2d(y[~prim], x[~prim], bins=(ey, ex), weights=w[~prim] ** 2)[0]
    return int(ok.sum())


def process_run(run_dir, resp, ex, ey, step="200 MB"):
    import uproot
    info = json.load(open(run_dir / "run_info.json"))
    acc = new_acc(ex, ey)
    f = uproot.open(run_dir / "phsp.root")
    tree = f[[k for k in f.keys() if k.split(";")[0] == "phsp"][0]]
    n_in = sum(add_batch(b, resp, ex, ey, acc) for b in tree.iterate(BRANCHES, step_size=step, library="np"))
    n_sr = info["events"] / info["omega_cone_sr"]                         # emitted photons per steradian
    return acc, n_sr, info, n_in


def to_pixels(coarse, ex, ey, sigma_mm, bin_mm, image_bin):
    """smooth the coarse image and interpolate it onto the projector grid (values per detector pixel)"""
    X, Y, p = pixel_grid(image_bin)
    sm = ndi.gaussian_filter(coarse, sigma_mm / bin_mm, mode="nearest")
    xc, yc = 0.5 * (ex[1:] + ex[:-1]), 0.5 * (ey[1:] + ey[:-1])
    f = RegularGridInterpolator((yc, xc), sm, bounds_error=False, fill_value=None)
    yy, xx = np.meshgrid(Y, X, indexing="ij")
    return (f(np.c_[yy.ravel(), xx.ravel()]).reshape(yy.shape) * (p * p) / (bin_mm * bin_mm)).astype(np.float32)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--name", required=True)
    ap.add_argument("--bin-mm", type=float, default=1.0)
    ap.add_argument("--sigma-mm", type=float, default=5.0,
                    help="Gaussian smoothing of the scatter image (mm). Scatter has no structure above the Monte Carlo "
                         "noise at 0.02-0.5 cy/mm (quick test, 22 Oct 2026); 5 mm keeps the scatter noise ~6000x below "
                         "the subtracted-image texture at 0.2 cy/mm")
    ap.add_argument("--image-bin", type=int, default=1)
    ap.add_argument("--accumulate", action="store_true", help="add to an existing output/mc/<name>.npz")
    ap.add_argument("--delete-phsp", action="store_true")
    a = ap.parse_args()
    from pathlib import Path
    out = OUT_DIR / "mc"; out.mkdir(parents=True, exist_ok=True)
    ex, ey = coarse_edges(a.bin_mm)
    resp = Response(E_GRID)
    tot, n_sr_tot, runs = None, 0.0, []
    fn = out / f"{a.name}.npz"
    beam, background, model = None, None, None
    if a.accumulate and fn.exists():
        z = np.load(fn, allow_pickle=True)
        tot = {k: z[k] for k in ("primary", "scatter", "n_primary", "n_scatter", "w2_scatter")}
        n_sr_tot, runs = float(z["n_sr"]), list(z["runs"])
        beam = str(z["beam"])
        background = str(z["background"]) if "background" in z and str(z["background"]) != "None" else None
        model = str(z["phantom"]) if "phantom" in z and str(z["phantom"]) != "None" else "L2v2"
    for r in a.runs:
        r = Path(r)
        acc, n_sr, info, n_in = process_run(r, resp, ex, ey)
        m_r = info.get("phantom_model") or ("L2v2" if info.get("mode", "phantom") == "phantom" else None)
        if beam is None:
            beam, background, model = info["beam"], info.get("background"), m_r
        elif info["beam"] != beam or info.get("background") != background or m_r != model:
            raise SystemExit(f"{r}: beam/background/phantom {info['beam']}/{info.get('background')}/{m_r} differ from "
                             f"{beam}/{background}/{model}")
        tot = acc if tot is None else {k: tot[k] + acc[k] for k in tot}
        n_sr_tot += n_sr; runs.append(str(r))
        print("%s: %d photons at the plane, %.3g photons/sr emitted" % (r, n_in, n_sr), flush=True)
        if a.delete_phsp:
            (r / "phsp.root").unlink()
    np.savez_compressed(fn, **tot, n_sr=n_sr_tot, runs=np.array(runs), x_edges=ex, y_edges=ey, beam=beam,
                        background=str(background), phantom=str(model))
    # normalised images (keV/sr per coarse bin)
    S = tot["scatter"] / n_sr_tot
    Pm = tot["primary"] / n_sr_tot
    rel = np.sqrt(tot["w2_scatter"]) / np.maximum(tot["scatter"], 1e-30)
    tag = f"_bin{a.image_bin}" if a.image_bin > 1 else ""
    np.save(out / f"scatter_{a.name}{tag}.npy", to_pixels(S, ex, ey, a.sigma_mm, a.bin_mm, a.image_bin))
    np.save(out / f"mcprimary_{a.name}{tag}.npy", to_pixels(Pm, ex, ey, 0.0, a.bin_mm, a.image_bin))
    inner = tot["n_scatter"] > 0
    print("scatter: median relative uncertainty per %.1f mm bin %.3f (before smoothing)" % (a.bin_mm, np.median(rel[inner])))
    print("written:", fn)


if __name__ == "__main__":
    main()

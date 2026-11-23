#!/usr/bin/env python3
"""
Deterministic (noise-free) primary image of the phantom (--phantom L2v2 or L2) on the Revelation geometry,
polychromatic, with the same below-plane response (grid, covers, a-Se) as the Monte Carlo photons. Complements the GATE
run, which is used for the scatter image: total image = primary (this) + scatter (postprocess.py).

Geometry is exact and analytic, identical to the GATE model:
  - L2 v2: the focal spot lies on the axis of the revolved container (x = y = 0), so every ray stays in one meridional
    plane and its chords through the envelope, cavity and plate follow from a line-ellipse intersection;
  - L2 (real container): chords through the half-cylinder prisms from a line-circle intersection (phantoms.py);
  - sphere chords from the ray-to-centre distance; iodine discs as vertical cylinders in the plate.
Iodine discs = PMMA + iodine (mixture rule), so the plate is PMMA everywhere plus an iodine mass-thickness map.
Focal-spot blur is neglected (0.3 mm x (M - 1) < 0.04 mm at the detector).

Output: output/projector/<case>/maps.npz (path lengths, mm; iodine g/cm2), primary_<beam>[_open].npy
(detected energy per pixel per emitted photon per steradian, keV/sr) and primary2_<beam>.npy (its second moment,
keV^2/sr, for the quantum noise of the measurement-like images). <case> = config.case_name(): the background for
L2 v2 (e.g. PA12_MCT), L2_<background> for the real L2.
Usage: python projector.py [--phantom L2v2|L2] [--background PA12_MCT|PMMA_paraffin] [--set 0035|july28high]
                           [--beams LE HE] [--bin 1] [--no-spheres] [--no-discs]
"""
import argparse, json, time
import numpy as np
from config import (HERE, SYSTEM, SUPPORT, PADDLE, PHANTOM, OUT_DIR, BACKGROUNDS, spectra_dir, add_case_args,
                    apply_case_args, case_name)
from materials import mu, mu_rho_el, iodine_rows
import phantoms
from grid_detector import Response, PLANE_Z
from revelation_gate import geometry_z

SID = SYSTEM["sid_mm"]
P = json.load(open(PHANTOM["params_json"])); DV = P["derived"]
H, R_FLAT = P["H"], P["R_flat"]
Y_CUT = DV["R_cav"] - 0.5


def pixel_grid(binning=1):
    p = SYSTEM["pixel_mm"] * binning
    ny, nx = SYSTEM["image_rows"] // binning, SYSTEM["image_cols"] // binning
    x = SYSTEM["missed_tissue_mm"] + (np.arange(nx) + 0.5) * p          # columns: chest wall -> nipple
    y = (np.arange(ny) - ny / 2 + 0.5) * p                              # rows: along the chest wall
    return x, y, p


def interval_len(a0, a1, b0, b1):
    return np.clip(np.minimum(a1, b1) - np.maximum(a0, b0), 0, None)


def revolved_chord(rho_d, z0, R, a, b, c):
    """fraction t of each ray (source -> detector point at radius rho_d) inside the solid
    {|z - c| <= b} and {rho <= R or ((rho - R)/a)^2 + ((z - c)/b)^2 <= 1}; z in phantom coordinates (bottom = z0)"""
    t1, t2 = 1 - (z0 + c + b) / SID, 1 - (z0 + c - b) / SID
    with np.errstate(divide="ignore"):
        tR = np.where(rho_d > 0, R / rho_d, np.inf)
    q0 = SID - z0 - c
    A = rho_d ** 2 / a ** 2 + SID ** 2 / b ** 2
    B = -2 * R * rho_d / a ** 2 - 2 * q0 * SID / b ** 2
    C = R ** 2 / a ** 2 + q0 ** 2 / b ** 2 - 1
    disc = B * B - 4 * A * C
    sq = np.sqrt(np.clip(disc, 0, None))
    ta, tb = np.where(disc > 0, (-B - sq) / (2 * A), 0), np.where(disc > 0, (-B + sq) / (2 * A), 0)
    flat = interval_len(t1, t2, -np.inf, tR)
    ell = interval_len(t1, t2, ta, tb)
    both = interval_len(np.maximum(t1, ta), np.minimum(t2, tb), -np.inf, tR)
    return flat + ell - both


def sphere_chords(X, Y, p, z0, centres, r, out, max_pairs=2e7, x_off=0.0):
    """add the chord length (mm) through spheres of radius r to the image out (ny, nx); centres in phantom
    coordinates, placed at world x + x_off"""
    ny, nx = out.shape
    Cw = centres + np.array([x_off, 0, z0])
    m = SID / (SID - Cw[:, 2])
    xc, yc = Cw[:, 0] * m, Cw[:, 1] * m
    k = int(np.ceil(r * m.max() * 1.25 / p)) + 1
    di, dj = np.mgrid[-k:k + 1, -k:k + 1]; di, dj = di.ravel(), dj.ravel()
    nchunk = max(1, int(max_pairs // len(di)))
    acc = np.zeros(ny * nx)
    for s in range(0, len(Cw), nchunk):
        i = np.rint((yc[s:s + nchunk] - Y[0]) / p).astype(int)[:, None] + di[None]
        j = np.rint((xc[s:s + nchunk] - X[0]) / p).astype(int)[:, None] + dj[None]
        ok = (i >= 0) & (i < ny) & (j >= 0) & (j < nx)
        cs = np.broadcast_to(Cw[s:s + nchunk, None, :], i.shape + (3,))[ok]
        i, j = i[ok], j[ok]
        d = np.c_[X[j], Y[i], np.full(len(i), -SID)]                    # pixel - source
        d /= np.linalg.norm(d, axis=1, keepdims=True)
        v = cs - np.array([0, 0, SID])
        perp2 = (v * v).sum(1) - (v * d).sum(1) ** 2
        ch = 2 * np.sqrt(np.clip(r * r - perp2, 0, None))
        acc += np.bincount(i * nx + j, weights=ch, minlength=ny * nx)
    out += acc.reshape(ny, nx)


def disc_chord(X, Y, p, z0, xc, yc, rd, zc, t, out_add, weight):
    """iodine discs: vertical cylinders (radius rd, centre xc, yc, height t at zc); adds weight x chord (mm)"""
    ny, nx = out_add.shape
    m = SID / (SID - z0 - zc)
    k = int(np.ceil((rd * m + 1.0) / p)) + 2
    i0, j0 = int(round((yc * m - Y[0]) / p)), int(round((xc * m - X[0]) / p))
    si, sj = slice(max(i0 - k, 0), min(i0 + k + 1, ny)), slice(max(j0 - k, 0), min(j0 + k + 1, nx))
    xd, yd = np.meshgrid(X[sj], Y[si])
    A = xd ** 2 + yd ** 2; B = -2 * (xd * xc + yd * yc); C = xc ** 2 + yc ** 2 - rd ** 2
    disc = B * B - 4 * A * C
    sq = np.sqrt(np.clip(disc, 0, None))
    ta, tb = (-B - sq) / (2 * A), (-B + sq) / (2 * A)
    t1, t2 = 1 - (z0 + zc + t / 2) / SID, 1 - (z0 + zc - t / 2) / SID
    L = np.where(disc > 0, interval_len(ta, tb, t1, t2), 0) * np.sqrt(A + SID ** 2)
    out_add[si, sj] += weight * L


def path_maps(binning=1, spheres=True, discs=True):
    """path lengths (mm) per material along the ray to every pixel, down to the phase-space plane"""
    t0 = time.time()
    X, Y, p = pixel_grid(binning)
    xd, yd = np.meshgrid(X, Y)
    rho = np.hypot(xd, yd)
    Lray = np.sqrt(rho ** 2 + SID ** 2)                                 # source -> pixel (z = 0)
    uz = SID / Lray
    z = geometry_z(); z0 = z["phantom"][0]
    maps = {}
    maps["paddle"] = np.full(rho.shape, PADDLE["mm"]) / uz
    maps["support"] = np.full(rho.shape, SUPPORT["mm"]) / uz
    x_off = phantoms.x_offset()
    zc, t_plate = phantoms.plate()
    if phantoms.is_l2():
        env, cav, plate = phantoms.l2_chords(xd, yd, Lray, z0)
    else:
        env = revolved_chord(rho, z0, R_FLAT, P["a_out"], H / 2, H / 2) * Lray
        cav = revolved_chord(rho, z0, R_FLAT, DV["a_in"], DV["b_in"], DV["ci"]) * Lray
        # plate: slab |z - ci| <= t/2 inside rho <= R_plate and |y| <= Y_CUT (the GATE plate is trimmed there)
        sinp = np.abs(yd) / np.maximum(rho, 1e-9)
        with np.errstate(divide="ignore"):
            r_eff = np.minimum(DV["R_plate"], np.where(sinp > 0, Y_CUT / sinp, np.inf))
            t_r = np.where(rho > 0, r_eff / rho, np.inf)
        t1, t2 = 1 - (z0 + DV["ci"] + P["t_plate"] / 2) / SID, 1 - (z0 + DV["ci"] - P["t_plate"] / 2) / SID
        plate = interval_len(t1, t2, -np.inf, t_r) * Lray
    sph = np.zeros(rho.shape)
    if spheres:
        S = np.load(phantoms.phantom_dir() / "spheres.npz")
        for kcls, dmm in enumerate(PHANTOM["sphere_diameters_mm"]):
            sphere_chords(X, Y, p, z0, S["centres"][S["size_class"] == kcls], dmm / 2, sph, x_off=x_off)
            print("  spheres %.2f mm done (%.0f s)" % (dmm, time.time() - t0), flush=True)
    iod = np.zeros(rho.shape)                                            # g/cm2
    if discs:
        rows = iodine_rows()                                             # from data/L2v2_materials.db
        for d in phantoms.discs():
            if d["key"] not in rows:
                continue                                                 # no iodine material: plate PMMA
            disc_chord(X, Y, p, z0, d["x"] + x_off, d["y"], d["diam"] / 2, zc, t_plate - 0.002, iod,
                       rows[d["key"]]["conc_g_cm3"] / 10.0)
    maps["PMMA"] = env - cav + plate
    maps["oil"] = cav - plate - sph
    maps["spheres"] = sph
    maps["iodine_g_cm2"] = iod
    maps["cavity"], maps["plate"] = cav, plate                          # for the analysis masks
    maps["air"] = (SID - PLANE_Z) / uz - maps["paddle"] - maps["support"] - env
    maps["uz"] = uz
    maps["domega_sr"] = p * p * SID / Lray ** 3
    if (maps["oil"] < -1e-3).any():
        print("  warning: negative oil path in %d pixels (min %.4f mm)" % ((maps["oil"] < -1e-3).sum(), maps["oil"].min()))
    print("path maps done (%.0f s)" % (time.time() - t0), flush=True)
    return {k: v.astype(np.float32) for k, v in maps.items()}, dict(x=X, y=Y, pixel_mm=p)


def material_of(part):
    """material of each path-length map; spheres and liquid follow the selected background (config.set_background)"""
    return {"paddle": PADDLE["material"], "support": SUPPORT["material"], "PMMA": PHANTOM["container_material"],
            "oil": PHANTOM["liquid_material"], "spheres": PHANTOM["sphere_material"], "air": "G4_AIR"}[part]


def out_dir():
    """projector outputs per case: output/projector/<case>/ (config.case_name())"""
    return OUT_DIR / "projector" / case_name()


def primary_image(maps, beam, open_field=False, min_weight=1e-6, second_moment=False):
    """detected energy per pixel per emitted photon per steradian (keV/sr). With second_moment=True also the second
    moment of the energy deposited per photon (keV^2/sr): N_sr times it is the quantum-noise variance of the pixel
    for an energy-integrating detector (each detected photon deposits its energy; K-escape neglected)"""
    spec = np.load(spectra_dir() / f"{beam}.npz")
    e, w = spec["centres_keV"], spec["weights"]
    keep = w > min_weight * w.max()
    e, w = e[keep], w[keep] / w[keep].sum()
    resp = Response(e)
    inv = 1.0 / maps["uz"].astype(np.float64)
    names = ["paddle", "support", "air"] + ([] if open_field else ["PMMA", "oil", "spheres"])
    air = maps["air"] + (maps["PMMA"] + maps["oil"] + maps["spheres"] if open_field else 0)  # open: phantom -> air
    img = np.zeros(inv.shape)
    img2 = np.zeros(inv.shape) if second_moment else None
    for k, (ek, wk) in enumerate(zip(e, w)):
        A = np.zeros(inv.shape)
        for nme in names:
            L = air if nme == "air" else maps[nme]
            A += mu(material_of(nme), ek) * L
        if not open_field:
            A += mu_rho_el("I", ek) * maps["iodine_g_cm2"]
        term = wk * np.exp(-A) * resp.primary_factor_e(k, inv)
        img += term
        if second_moment:
            img2 += term * ek
    img = (img * maps["domega_sr"]).astype(np.float32)
    if second_moment:
        return img, (img2 * maps["domega_sr"]).astype(np.float32)
    return img


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--beams", nargs="+", default=["LE", "HE"])
    ap.add_argument("--bin", type=int, default=1)
    ap.add_argument("--no-spheres", action="store_true")
    ap.add_argument("--no-discs", action="store_true")
    add_case_args(ap)
    a = ap.parse_args()
    case = apply_case_args(a)
    bg = PHANTOM["background"]
    print("case %s: phantom %s, background %s (%s), spectra %s" % (case, PHANTOM["model"], bg, BACKGROUNDS[bg]["description"],
                                                                  spectra_dir()), flush=True)
    out = out_dir(); out.mkdir(parents=True, exist_ok=True)
    tag = ("_nospheres" if a.no_spheres else "") + ("_nodiscs" if a.no_discs else "") + (f"_bin{a.bin}" if a.bin > 1 else "")
    maps, grid = path_maps(a.bin, spheres=not a.no_spheres, discs=not a.no_discs)
    np.savez_compressed(out / f"maps{tag}.npz", **maps, **grid)
    for b in a.beams:
        img, img2 = primary_image(maps, b, second_moment=True)
        np.save(out / f"primary_{b}{tag}.npy", img)
        np.save(out / f"primary2_{b}{tag}.npy", img2)          # second moment: quantum noise (measlike.py)
        np.save(out / f"primary_{b}_open{'_bin%d' % a.bin if a.bin > 1 else ''}.npy", primary_image(maps, b, open_field=True))
        print("primary images %s done" % b, flush=True)


if __name__ == "__main__":
    main()

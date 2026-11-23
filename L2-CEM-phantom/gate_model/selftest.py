#!/usr/bin/env python3
"""
Checks of the parts that run without GATE (no opengate needed):
  0. the material database data/L2v2_materials.db;
  1. grid table against a brute-force average over lamella phase and path;
  2. projector path lengths at the phantom centre against the wall/plate thicknesses;
  3. normalisation: synthetic "Monte Carlo" photons, emitted isotropically into the source cone like the GATE source
     and attenuated through paddle, support and air, are binned by postprocess.add_batch; the result must equal the
     projector's open-field primary image (same units, keV/sr per pixel) within the statistical error;
  4. the sphere packing (if phantom/spheres.npz exists): no overlaps, no sphere closer to a wall than the margin or
     crossing a mesh facet, volume fraction within 0.02 of the target;
  5. the real L2 model (if phantom_L2/ exists): closed meshes, path lengths at an interior pixel against the wall and
     plate thicknesses, and its sphere packing as in 4.
Uses spectra from the spectra folder if present, otherwise a rough Kramers spectrum (only for these checks).
Usage: python selftest.py
"""
import os, sys, tempfile, json
from pathlib import Path
import numpy as np

TMP = Path(tempfile.mkdtemp(prefix="revelation_selftest_"))
os.environ.setdefault("REVELATION_OUT_DIR", str(TMP / "out"))
from config import spectra_dir, HERE                              # noqa: E402
if not (spectra_dir() / "LE.npz").exists():
    os.environ["REVELATION_SPECTRA_DIR"] = str(TMP / "spectra")
    import config; config.SPECTRA_DIR = TMP / "spectra"
    from materials import mu_rho_el
    (TMP / "spectra").mkdir()
    for name, kvp, el, rho, t in (("LE", 28, "Rh", 12.41, 0.05), ("HE", 49, "Ti", 4.506, 1.0)):
        edges = np.arange(5.0, kvp + 0.01, 0.5); k = 0.5 * (edges[1:] + edges[:-1])
        f = np.clip(kvp - k, 0, None) / k * np.exp(-mu_rho_el("Be", k) * 1.85 * 0.05 - mu_rho_el(el, k) * rho * t / 10)
        np.savez(TMP / "spectra" / f"{name}.npz", edges_keV=edges, centres_keV=k, weights=f / f.sum(),
                 fluence_per_mAs_1m=1e7, kerma_uGy_per_mAs_1m=0.0, hvl_mm_al=0.0)
    print("(no SpekPy spectra found: using a rough Kramers spectrum for the checks)")

import grid_detector as G                                          # noqa: E402
import projector as PJ                                             # noqa: E402
import postprocess as PP                                           # noqa: E402
from materials import mu                                           # noqa: E402
from config import SYSTEM, SUPPORT, PADDLE                         # noqa: E402
import revelation_gate as RG                                       # noqa: E402

ok_all = True


def check(name, ok, msg):
    global ok_all
    ok_all &= bool(ok)
    print("%-4s %-28s %s" % ("OK" if ok else "FAIL", name, msg))


# 0. material database (data/L2v2_materials.db): parses, fractions sum to 1, every material used is defined, and the
#    sphere and oil compositions match their chemical formulas
import materials as MT                                             # noqa: E402
problems = MT.check_db()
AW = dict(C=12.011, H=1.008, N=14.007, O=15.999, Ca=40.078)
def formula(f):
    m = sum(AW[e] * n for e, n in f.items()); return {e: AW[e] * n / m for e, n in f.items()}
pa, cc = formula(dict(C=12, H=23, N=1, O=1)), formula(dict(Ca=1, C=1, O=3))
ref = {"PA12_CaCO3_5p5": {e: 0.945 * pa.get(e, 0) + 0.055 * cc.get(e, 0) for e in AW},
       "MCT_oil": formula(dict(C=29, H=54, O=6))}
dev = max(abs(MT.DB[m]["weights"].get(e, 0) - w) for m, r in ref.items() for e, w in r.items())
check("material database", not problems and dev < 1e-4,
      "%d materials, %s; max deviation from formulas %.1e" % (len(MT.DB), "no problems" if not problems else "; ".join(problems), dev))

# 1. grid
e, ang = 20.0, np.radians(8)
s = G.h * np.tan(ang)
p = (np.arange(5000) + 0.5) / 5000 * G.D; zz = (np.arange(2000) + 0.5) / 2000       # midpoint rule
lead = (((p[:, None] + s * zz[None]) % G.D) < G.d).mean(1)
T_bf = np.exp(-G.h * (mu("G4_Pb", e) * lead + mu("CarbonFibre", e) * (1 - lead))).mean()
T_tab = G.GridTable(np.array([19.0, 20.0, 21.0]))(np.array([e]), np.array([s]), np.array([1.0]))[0]
check("grid transmission (8 deg)", abs(T_tab / T_bf - 1) < 0.002, "table %.4f, brute force %.4f" % (T_tab, T_bf))

# 2. path lengths
import config; config.set_phantom("L2v2")                          # sections 2-4 test L2 v2, whatever REVELATION_PHANTOM says
maps, grid = PJ.path_maps(binning=8, spheres=False, discs=False)
x, y = grid["x"], grid["y"]
i, j = np.argmin(np.abs(y)), np.argmin(np.abs(x - 60))
uz = maps["uz"][i, j]
P = PJ.P
exp_pmma = (P["t_top"] + P["t_bot"] + P["t_plate"]) / uz
exp_oil = (P["H"] - P["t_top"] - P["t_bot"] - P["t_plate"]) / uz
check("path lengths at centre", abs(maps["PMMA"][i, j] - exp_pmma) < 1e-3 and abs(maps["oil"][i, j] - exp_oil) < 1e-3,
      "PMMA %.3f/%.3f, oil %.3f/%.3f mm" % (maps["PMMA"][i, j], exp_pmma, maps["oil"][i, j], exp_oil))

# 3. normalisation: synthetic open-field photons vs projector open field
beam = "LE"
P0 = PJ.primary_image(maps, beam, open_field=True)                # keV/sr per (binned) pixel
spec = np.load(spectra_dir() / f"{beam}.npz")
rng = np.random.default_rng(1)
theta_max, phi, omega = RG.source_cone()
n = 8_000_000
ct = rng.uniform(np.cos(np.radians(theta_max)), 1, n); st = np.sqrt(1 - ct * ct)
ph = np.radians(rng.uniform(phi[0], phi[1], n))
u = -np.c_[st * np.cos(ph), st * np.sin(ph), ct]                    # Geant4 GPS convention
src = np.c_[np.zeros((n, 2)), np.full(n, SYSTEM["sid_mm"])]
eph = rng.choice(spec["centres_keV"], n, p=spec["weights"])
z = RG.geometry_z()
tpl = (G.PLANE_Z - SYSTEM["sid_mm"]) / u[:, 2]
pos = src + tpl[:, None] * u
k = RG.APERTURE_Z; kk = (SYSTEM["sid_mm"] - k) / SYSTEM["sid_mm"]
xap = (k - SYSTEM["sid_mm"]) / u[:, 2] * u[:, 0]; yap = (k - SYSTEM["sid_mm"]) / u[:, 2] * u[:, 1]
inside = (xap >= RG.FIELD_X[0] * kk) & (xap <= RG.FIELD_X[1] * kk) & (yap >= RG.FIELD_Y[0] * kk) & (yap <= RG.FIELD_Y[1] * kk)
L = 1.0 / -u[:, 2]
att = (mu(PADDLE["material"], eph) * PADDLE["mm"] + mu(SUPPORT["material"], eph) * SUPPORT["mm"]) * L \
    + mu("G4_AIR", eph) * ((SYSTEM["sid_mm"] - G.PLANE_Z) * L - (PADDLE["mm"] + SUPPORT["mm"]) * L)
keep = inside & (rng.uniform(size=n) < np.exp(-att))
b = dict(KineticEnergy=eph[keep] / 1000.0, PrePosition_X=pos[keep, 0], PrePosition_Y=pos[keep, 1], PrePosition_Z=pos[keep, 2],
         PreDirection_X=u[keep, 0], PreDirection_Y=u[keep, 1], PreDirection_Z=u[keep, 2],
         UnscatteredPrimaryFlag=np.ones(keep.sum(), int), PDGCode=np.full(keep.sum(), 22))
ex, ey = PP.coarse_edges(10.0)
acc = PP.new_acc(ex, ey)
PP.add_batch(b, G.Response(G.E_GRID), ex, ey, acc)
mc = acc["primary"] / (n / omega) / 10.0 ** 2                      # keV/sr per mm2
xx, yy = np.meshgrid(x, y)
cnt = np.histogram2d(yy.ravel(), xx.ravel(), bins=(ey, ex))[0]
pj = np.histogram2d(yy.ravel(), xx.ravel(), bins=(ey, ex), weights=P0.ravel())[0] / np.maximum(cnt, 1) / grid["pixel_mm"] ** 2
sel = (cnt >= 0.95 * (10.0 / grid["pixel_mm"]) ** 2) & (acc["n_primary"] > 1000)
ratio = mc[sel] / pj[sel]
xc, yc = np.meshgrid(0.5 * (ex[1:] + ex[:-1]), 0.5 * (ey[1:] + ey[:-1]))
inner = (np.hypot(xc, yc) < 120)[sel]                               # under the phantom footprint vs outside
r_in, r_out = mc[sel][inner].sum() / pj[sel][inner].sum(), mc[sel][~inner].sum() / pj[sel][~inner].sum()
check("MC vs projector normalisation", abs(np.median(ratio) - 1) < 0.003 and abs(r_in - r_out) < 0.004,
      "median ratio %.4f (footprint %.4f, outside %.4f; stat. %.4f per bin, %d bins)"
      % (np.median(ratio), r_in, r_out, 1 / np.sqrt(np.median(acc["n_primary"][sel])), sel.sum()))

# 4. sphere packing
def packing_check(label, folder):
    f = folder / "spheres.npz"
    if not f.exists():
        print("(%s: no sphere packing yet, %s)" % (label, f))
        return
    import sphere_packing as SP
    S = np.load(f)
    bad_wall, bad_pair, maxov = SP.verify(S["centres"], S["radii"])
    summ = json.load(open(folder / "spheres_summary.json"))
    frac = summ["achieved_fraction"]
    mesh = summ.get("mesh_check", {})
    ok = bad_wall.sum() == 0 and len(bad_pair) == 0 and mesh.get("crossing_a_facet", 1) == 0 and abs(frac - summ["target_fraction"]) < 0.02
    check(label, ok, "%d spheres, %d wall violations, %d overlapping pairs, %s spheres crossing a mesh facet, "
          "volume fraction %.3f (target %.2f)" % (len(S["radii"]), bad_wall.sum(), len(bad_pair), mesh.get("crossing_a_facet", "?"),
                                                  frac, summ["target_fraction"]))


packing_check("sphere packing (L2 v2)", HERE / "phantom")

# 5. real L2 (half-cylinder container, new plate)
if (HERE / "phantom_L2" / "meshes.json").exists():
    import config, phantoms
    import sphere_packing as SP
    config.set_phantom("L2"); SP.use_model("L2")
    mj = json.load(open(HERE / "phantom_L2" / "meshes.json"))
    check("L2 meshes", all(mj[k]["closed"] and mj[k]["oriented"] and abs(mj[k]["mesh_volume_ml"] / mj[k]["exact_volume_ml"] - 1) < 1e-3
                          for k in ("envelope", "filling", "plate")),
          ", ".join("%s %.1f/%.1f ml" % (k, mj[k]["mesh_volume_ml"], mj[k]["exact_volume_ml"]) for k in ("envelope", "filling", "plate")))
    g = phantoms.l2()
    maps, grid = PJ.path_maps(binning=8, spheres=False, discs=False)
    m_mid = SYSTEM["sid_mm"] / (SYSTEM["sid_mm"] - RG.geometry_z()["phantom"][0] - g["zc"])
    j = int(np.argmin(np.abs(grid["x"] - (g["x_f"] + 60.0) * m_mid))); i = int(np.argmin(np.abs(grid["y"] - 20.0)))
    uz = maps["uz"][i, j]
    exp_pmma = (g["t_plate"] + g["H"] - (g["zt"] - g["zb"])) / uz
    exp_oil = (g["zt"] - g["zb"] - g["t_plate"]) / uz
    check("L2 path lengths at an interior pixel", abs(maps["PMMA"][i, j] - exp_pmma) < 1e-3 and abs(maps["oil"][i, j] - exp_oil) < 1e-3,
          "PMMA %.3f/%.3f, oil %.3f/%.3f mm" % (maps["PMMA"][i, j], exp_pmma, maps["oil"][i, j], exp_oil))
    packing_check("sphere packing (L2)", HERE / "phantom_L2")
    config.set_phantom("L2v2")
print("\nALL CHECKS PASSED" if ok_all else "\nSOME CHECKS FAILED")
sys.exit(0 if ok_all else 1)

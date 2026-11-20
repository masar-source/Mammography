#!/usr/bin/env python3
"""
X-ray spectra of the Revelation beams with SpekPy (installed with opengate): W anode, Be window, Rh or Ti filter,
optional extra Al-equivalent filtration. Air, paddle and phantom are NOT included in the saved spectra; they are
simulated in GATE (world = air). The dose check adds them (air to the entrance plane and the paddle), because the
DICOM entrance dose is taken to be the air kerma at the breast entrance, below the paddle [A].

Writes <spectra folder>/<beam>.npz (spectra/ for the measured set 0035, spectra_july28high/ for July; --set) with: edges_keV (N+1), weights (N, photons per bin, normalised to 1),
fluence_per_mAs_1m (photons / cm2 / mAs at 1 m, all energies) and kerma_uGy_per_mAs_1m, plus a printed check of the
entrance air kerma against the DICOM value (EntranceDoseInmGy).

Usage: python spectra.py [--set 0035|july28high] [--kvp 28 30 32 34] [--extra-al MM | --fit-extra-al]
       (always writes LE and HE)
       --fit-extra-al: extra Al that makes the LE entrance dose match DICOM; the HE ratio is then an independent check.
       --set: the measured reference whose mAs and entrance doses are used (each X-ray unit gets its own fit); for
       july28high run  python measured_set.py --set july28high  first (reads them from the DICOM headers).
Also writes spectra_info.json (extra Al, HVLs, dose ratios) in the same folder.
"""
import argparse, json
import numpy as np
from config import HERE, SYSTEM, BEAMS, PADDLE, MEASURED, MEASURED_SETS, set_measured, spectra_dir

# SpekPy is not compatible with numpy 2.x in some versions; the opengate test suite applies this patch.
_orig_arange = np.arange


def _arange(*args, **kw):
    def scalar(v):
        return np.asarray(v).item(0) if isinstance(v, (np.ndarray, list)) and np.size(v) > 0 else v
    return _orig_arange(*[scalar(a) for a in args], **kw)


np.arange = _arange


SPEKPY_NAMES = {"G4_POLYCARBONATE": "Polycarbonate (Makrolon or Lexan)", "G4_AIR": "Air"}


def _spek(kvp, filt, extra_al, dk=0.5, to_entrance=False):
    """SpekPy spectrum at the tube exit (W anode, Be window, beam filter, extra Al); with to_entrance=True also the
    air from the focal spot to the breast entrance and the compression paddle (only for the dose check)"""
    import spekpy as sp
    s = sp.Spek(kvp=kvp, th=SYSTEM["anode_angle_deg"], targ=SYSTEM["anode"], dk=dk, mu_data_source="nist")
    s.filter("Be", SYSTEM["be_window_mm"])
    s.filter(filt[0], filt[1])
    if extra_al > 0:
        s.filter("Al", extra_al)
    if to_entrance:
        s.filter(SPEKPY_NAMES["G4_AIR"], SYSTEM["source_to_entrance_mm"] - PADDLE["mm"])
        s.filter(SPEKPY_NAMES[PADDLE["material"]], PADDLE["mm"])
    return s


def make_spectrum(kvp, filt, extra_al=0.0, dk=0.5):
    s = _spek(kvp, filt, extra_al, dk)
    k, f = s.get_spectrum()                 # bin centres (keV), fluence per keV per cm2 per mAs at 1 m
    phi = f * dk                            # photons per bin
    edges = np.r_[k - dk / 2, k[-1] + dk / 2]
    keep = phi > 0
    first, last = np.argmax(keep), len(keep) - np.argmax(keep[::-1])
    edges, phi, k = edges[first:last + 1], phi[first:last], k[first:last]
    return dict(edges_keV=edges, centres_keV=k, weights=phi / phi.sum(), fluence_per_mAs_1m=phi.sum(),
                kerma_uGy_per_mAs_1m=float(s.get_kerma()), hvl_mm_al=float(s.get_hvl1()), extra_al_mm=extra_al)


def esak_model(b, extra_al, to_entrance=True):
    """entrance air kerma (mGy) for the DICOM exposure: kerma at 1 m (SpekPy) x inverse square to the entrance plane;
    with to_entrance the air path and the paddle are included (the DICOM value is assumed to be below the paddle)"""
    k1m = _spek(b["kvp"], b["filter"], extra_al, to_entrance=to_entrance).get_kerma()   # uGy/mAs at 1 m, in vacuo
    return k1m * (1000.0 / SYSTEM["source_to_entrance_mm"]) ** 2 * b["mAs"] / 1000.0


def fit_extra_al(lo=0.0, hi=2.0, tol=0.002):
    """Al thickness that makes the LE model entrance dose equal to the DICOM value (bisection)"""
    b = BEAMS["LE"]
    f = lambda t: esak_model(b, t) / b["esak_mGy"] - 1.0
    if f(lo) <= 0:
        return lo
    if f(hi) > 0:
        raise SystemExit("even %.1f mm extra Al does not bring the LE dose down to the DICOM value" % hi)
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if f(mid) > 0 else (lo, mid)
    return 0.5 * (lo + hi)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kvp", type=float, nargs="*", default=[], help="extra LE kVp values (Rh filter)")
    ap.add_argument("--extra-al", type=float, default=None,
                    help="extra Al-equivalent tube filtration (mm), overriding config SYSTEM['extra_filtration_al_mm']")
    ap.add_argument("--fit-extra-al", action="store_true",
                    help="choose the extra Al so that the LE entrance dose (with air and paddle) matches DICOM; "
                         "HE is then an independent check (this is the default when the config value is 'fit')")
    ap.add_argument("--set", choices=list(MEASURED_SETS), default=None,
                    help="measured reference set (default: REVELATION_SET or 0035)")
    a = ap.parse_args()
    name = set_measured(a.set)
    if not MEASURED["set_ready"]:
        raise SystemExit("set %s: mAs and entrance doses not known yet. Connect the drive with the images and run\n"
                         "  python measured_set.py --set %s" % (name, name))
    out = spectra_dir()
    print("measured set %s (%s) -> %s" % (name, MEASURED_SETS[name]["description"], out))
    out.mkdir(parents=True, exist_ok=True)
    setting = SYSTEM["extra_filtration_al_mm"]
    fitted = a.extra_al is None and (a.fit_extra_al or setting == "fit")
    if fitted:
        extra_al = fit_extra_al()
        print("fitted extra filtration: %.3f mm Al (LE entrance dose = DICOM)" % extra_al)
    else:
        extra_al = float(a.extra_al if a.extra_al is not None else setting)
    beams = dict(BEAMS)
    for kv in a.kvp:
        beams[f"LE{int(kv)}"] = dict(kvp=kv, filter=BEAMS["LE"]["filter"], mAs=None, esak_mGy=None)
    info = dict(measured_set=name, extra_filtration_al_mm=extra_al, extra_filtration_fitted=fitted, paddle=PADDLE,
                source_to_entrance_mm=SYSTEM["source_to_entrance_mm"], beams={})
    print("extra Al filtration %.3f mm; paddle %.1f mm %s; air %.0f mm" % (extra_al, PADDLE["mm"], PADDLE["material"],
                                                                         SYSTEM["source_to_entrance_mm"] - PADDLE["mm"]))
    for name, b in beams.items():
        s = make_spectrum(b["kvp"], b["filter"], extra_al)
        np.savez(out / f"{name}.npz", **s)
        msg = "%-5s %4.0f kV %s %.3f mm: HVL %.2f mm Al, mean %.1f keV" % (
            name, b["kvp"], b["filter"][0], b["filter"][1], s["hvl_mm_al"], (s["centres_keV"] * s["weights"]).sum())
        info["beams"][name] = dict(kvp=b["kvp"], hvl_mm_al=s["hvl_mm_al"])
        if b.get("esak_mGy"):
            e0 = esak_model(b, extra_al, to_entrance=False)
            e1 = esak_model(b, extra_al, to_entrance=True)
            msg += " | entrance dose vs DICOM %.2f mGy: tube only %.2f (x%.2f), with air + paddle %.2f (x%.2f)" % (
                b["esak_mGy"], e0, e0 / b["esak_mGy"], e1, e1 / b["esak_mGy"])
            info["beams"][name].update(esak_dicom_mGy=b["esak_mGy"], esak_model_mGy=e1, ratio=e1 / b["esak_mGy"])
        print(msg)
    he = info["beams"].get("HE", {}).get("ratio")
    if fitted and he is not None:
        info["he_check"] = "ok" if abs(he - 1) <= 0.10 else "not confirmed"
        if abs(he - 1) <= 0.10:
            print("HE check: model/DICOM = %.2f, within 10 %%: extra filtration consistent with both beams" % he)
        else:
            print("WARNING: HE model/DICOM = %.2f after the fit. Extra Al filtration alone does not explain both "
                  "doses (tube output, anode angle or the DICOM dose definition may differ); report this." % he)
    json.dump(info, open(out / "spectra_info.json", "w"), indent=1)


if __name__ == "__main__":
    main()

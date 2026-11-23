#!/usr/bin/env python3
"""
Measured reference sets (config.MEASURED_SETS): reads the raw LE/HE pairs and their DICOM headers.
  0035        the single pair TICEM_RAW_L-CC_DIAGNOSIS_0035 (unit 1024; old L2 container and old plate)
  july28high  the 12 pairs of the real L2 with the new plate, 28 kV high mAs (unit 10159, R-CC), listed in
              investigation_0035_0037/scripts/july2026/exposures.json; the images are on the SanDisk drive

For a set whose exposure values are not in config.py (july28high), this script reads them once from the headers and
the images and writes data/measured_<set>.json: for every pair the exact mAs of its LE and its HE image
(ExposureInuAs), the entrance doses (EntranceDoseInmGy) and the HE open-field stored value (median of the open field
next to the phantom). Nothing is averaged over pairs: config.set_measured() takes the values of the ONE selected pair
(--pair, default HIGHMAS_28KV2_00), which spectra.py (extra Al of that unit) and measlike.py then use. Geometry tags (SID, source-entrance distance, pixel spacing,
kVp and filters) are compared with config.py, and differences are printed.

Images are returned in the simulation's orientation: columns from the chest wall (x), rows along it (y). R-CC images
are mirrored left-right, as in the July analysis (july2026/pipe.py).

The July images are looked up on the SanDisk, under --root, or as session folders copied into the project folder
(CEM_NEWL2_MASAR_HIGHMAS_28KV_02 and ..._28KV2_02 are there).
Usage: python measured_set.py --set july28high [--root "/Volumes/SanDisk/Masar/New L2 Measurements"]
       python measured_set.py --set 0035          (reads the headers and prints them next to config.py; no file written)
"""
import argparse, glob, json
from pathlib import Path
import numpy as np
from config import MEASURED_SETS, MEASURED, BEAMS, SYSTEM, PROJECT

JULY_ROOT = "/Volumes/SanDisk/Masar/New L2 Measurements"     # prefix of the paths in exposures.json
ELEMENT = {"RHODIUM": "Rh", "TITANIUM": "Ti", "MOLYBDENUM": "Mo", "SILVER": "Ag", "ALUMINUM": "Al", "ALUMINIUM": "Al",
           "COPPER": "Cu"}


def phantom_mask(H, ceiling=None):
    """phantom = HE < 0.64 x HE open field (largest component, holes filled); returns (mask, HE open-field level)"""
    from scipy import ndimage as ndi
    ceiling = ceiling or MEASURED["ceiling"]
    hv = H[H < ceiling - 1]
    open_he = float(np.median(hv[hv > np.percentile(hv, 75)]))
    m = H < 0.64 * open_he
    m = ndi.binary_opening(m, iterations=3)
    lab, n = ndi.label(m)
    m = lab == (np.argmax(ndi.sum(m, lab, range(1, n + 1))) + 1)
    return ndi.binary_fill_holes(m), open_he


def _num(ds, *keys):
    for k in keys:
        v = ds.get(k)
        if v is not None and str(v) != "":
            try:
                return float(v)
            except (TypeError, ValueError):
                return str(v)
    return None


def header(ds):
    """the header values used by the model"""
    uas = _num(ds, "ExposureInuAs")
    mas = uas / 1000.0 if isinstance(uas, float) else _num(ds, "ExposureInmAs", "Exposure")
    return dict(kvp=_num(ds, "KVP"), mAs=mas, esak_mGy=_num(ds, "EntranceDoseInmGy"),
                filter=[str(ds.get("FilterMaterial", "")), _num(ds, "FilterThicknessMinimum"), _num(ds, "FilterThicknessMaximum")],
                anode=str(ds.get("AnodeTargetMaterial", "")), sid_mm=_num(ds, "DistanceSourceToDetector"),
                source_to_entrance_mm=_num(ds, "DistanceSourceToEntrance"), thickness_mm=_num(ds, "BodyPartThickness"),
                force_N=_num(ds, "CompressionForce"), pixel_mm=float(ds.ImagerPixelSpacing[0]),
                rows=int(ds.Rows), cols=int(ds.Columns), detector_id=str(ds.get("DetectorID", "")),
                device_serial=str(ds.get("DeviceSerialNumber", "")), laterality=str(ds.get("ImageLaterality", "")),
                view=str(ds.get("ViewPosition", "")), magnification=_num(ds, "EstimatedRadiographicMagnificationFactor"),
                date=str(ds.get("AcquisitionDate", ds.get("StudyDate", ""))))


def _files(name, root=None):
    """[(id, LE path, HE path)] of a set"""
    src = MEASURED_SETS[name]["source"]
    if src["kind"] == "dir":
        import pydicom
        L = H = None
        for f in sorted(glob.glob(str(Path(src["path"]) / "*.IMA"))):
            kv = float(pydicom.dcmread(f, stop_before_pixels=True).KVP)
            L, H = (f, H) if kv < 40 else (L, f)
        if L is None or H is None:
            raise SystemExit("could not find both LE and HE images in %s" % src["path"])
        return [(Path(src["path"]).name, L, H)]
    ex = json.load(open(src["exposures"]))
    sel = [e for e in ex if e["kv"] == src["kvp"] and e["level"] == src["level"] and e["unit"] == MEASURED_SETS[name]["unit"]]
    return [(e["id"], locate(e["L"], root), locate(e["H"], root)) for e in sel]


def locate(path, root=None):
    """a July image path of exposures.json (on the SanDisk) or its copy: under --root (same layout as 'New L2
    Measurements'), or the session folder copied directly into the project folder (as for the 28 kV high-mAs set)"""
    rest = path[len(JULY_ROOT):].lstrip("/") if path.startswith(JULY_ROOT) else None
    cands = [path]
    if rest is not None:
        parts = rest.split("/")
        for base in ([Path(root)] if root else []) + [PROJECT]:
            cands += [str(base / rest), str(base.joinpath(*parts[1:]))]       # with and without the '28kV' level
    for c in cands:
        if Path(c).exists():
            return c
    return cands[-1] if root is None else cands[2]


def pairs(name, root=None, images=True, ids=None):
    """yields dict(id, LE, HE, pixel_mm, header_LE, header_HE), images in the simulation orientation; only the pairs
    in ids if given"""
    import pydicom
    files = [f for f in _files(name, root) if not ids or f[0] in ids]
    missing = [f for _, l, h in files for f in (l, h) if not Path(f).exists()]
    if missing:
        raise SystemExit("%d of %d images of set %s not found, e.g. %s\nConnect the drive (or give --root)."
                         % (len(missing), 2 * len(files), name, missing[0]))
    for i, l, h in files:
        dl, dh = pydicom.dcmread(l, stop_before_pixels=not images), pydicom.dcmread(h, stop_before_pixels=not images)
        out = dict(id=i, header_LE=header(dl), header_HE=header(dh), pixel_mm=float(dl.ImagerPixelSpacing[0]))
        if images:
            flip = MEASURED_SETS[name]["mirror"] or out["header_LE"]["laterality"] == "R"
            out["LE"] = dl.pixel_array.astype(np.float64)[:, ::-1] if flip else dl.pixel_array.astype(np.float64)
            out["HE"] = dh.pixel_array.astype(np.float64)[:, ::-1] if flip else dh.pixel_array.astype(np.float64)
            out["mirrored"] = bool(flip)
        yield out


def summarise(name, root=None):
    rows = []
    for d in pairs(name, root):
        _, open_he = phantom_mask(d["HE"])
        rows.append(dict(id=d["id"], LE=d["header_LE"], HE=d["header_HE"], he_open_stored=open_he,
                         max_stored=dict(LE=float(d["LE"].max()), HE=float(d["HE"].max())), mirrored=d["mirrored"]))
        print("  %-18s LE %5.1f kV %6.1f mAs %6.3f mGy | HE %5.1f kV %6.2f mAs %6.3f mGy | HE open %6.0f | %s %s"
              % (d["id"], d["header_LE"]["kvp"], d["header_LE"]["mAs"], d["header_LE"]["esak_mGy"] or np.nan,
                 d["header_HE"]["kvp"], d["header_HE"]["mAs"], d["header_HE"]["esak_mGy"] or np.nan, open_he,
                 d["header_LE"]["laterality"], d["header_LE"]["view"]), flush=True)
    # geometry and beam quality against config.py
    warn = []
    h0l, h0h = rows[0]["LE"], rows[0]["HE"]
    for k, ref in (("sid_mm", SYSTEM["sid_mm"]), ("source_to_entrance_mm", SYSTEM["source_to_entrance_mm"]),
                   ("pixel_mm", SYSTEM["pixel_mm"])):
        if isinstance(h0l[k], float) and abs(h0l[k] - ref) > 1e-3 * max(ref, 1):
            warn.append("%s %.3f in the headers, %.3f in config.SYSTEM" % (k, h0l[k], ref))
    for b, h in (("LE", h0l), ("HE", h0h)):
        if h["kvp"] != BEAMS[b]["kvp"]:
            warn.append("%s kVp %s in the headers, %s in config.BEAMS" % (b, h["kvp"], BEAMS[b]["kvp"]))
        if h["filter"][0] and ELEMENT.get(h["filter"][0].upper(), h["filter"][0]).upper() != BEAMS[b]["filter"][0].upper():
            warn.append("%s filter %s in the headers, %s in config.BEAMS" % (b, h["filter"], BEAMS[b]["filter"]))
    if (h0l["rows"], h0l["cols"]) != (SYSTEM["image_rows"], SYSTEM["image_cols"]):
        warn.append("image %d x %d, config %d x %d" % (h0l["rows"], h0l["cols"], SYSTEM["image_rows"], SYSTEM["image_cols"]))
    mx = max(max(r["max_stored"].values()) for r in rows)
    if mx > MEASURED["ceiling"]:
        warn.append("stored values up to %.0f, above the ceiling %d of config.MEASURED" % (mx, MEASURED["ceiling"]))
    return dict(set=name, n_pairs=len(rows), warnings=warn, unit=h0l["device_serial"], detector_id=h0l["detector_id"],
                note="values per pair; each pair is treated on its own (config.set_measured selects one)", exposures=rows)


def exposures(name=None, all_pairs=False):
    """the measured pair to reproduce (the selected one, MEASURED["pair"]; all with all_pairs=True):
    [dict(id, mAs=dict(LE, HE), he_open_stored)] with the exact mAs of its LE and HE image (headers). Sets with values
    in config.py (0035) give their single pair from there."""
    name = name or MEASURED["set"]
    m = MEASURED_SETS[name]
    if m.get("values_json") and Path(m["values_json"]).exists():
        rows = json.load(open(m["values_json"]))["exposures"]
        return [dict(id=r["id"], mAs=dict(LE=r["LE"]["mAs"], HE=r["HE"]["mAs"]), he_open_stored=r["he_open_stored"])
                for r in rows if all_pairs or r["id"] == MEASURED.get("pair")]
    if m.get("mAs") is None:
        raise SystemExit("set %s: exposure values not read yet; run  python measured_set.py --set %s" % (name, name))
    pid = Path(m["source"]["path"]).name if m["source"]["kind"] == "dir" else name
    return [dict(id=pid, mAs=dict(m["mAs"]), he_open_stored=m["he_open_stored"])]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", choices=list(MEASURED_SETS), default="july28high")
    ap.add_argument("--root", default=None, help="where the 'New L2 Measurements' folder is, if not at %s" % JULY_ROOT)
    a = ap.parse_args()
    m = MEASURED_SETS[a.set]
    print("%s: %s" % (a.set, m["description"]), flush=True)
    s = summarise(a.set, a.root)
    if m.get("mAs"):
        print("config.py:  mAs LE %.2f / HE %.2f, entrance dose LE %.3f / HE %.3f mGy, HE open field %.0f"
              % (m["mAs"]["LE"], m["mAs"]["HE"], m["esak_mGy"]["LE"], m["esak_mGy"]["HE"], m["he_open_stored"]))
    for w in s["warnings"]:
        print("NOTE:", w)
    if m.get("values_json"):
        json.dump(s, open(m["values_json"], "w"), indent=1)
        print("written", m["values_json"], "(used by config.set_measured)")


if __name__ == "__main__":
    main()

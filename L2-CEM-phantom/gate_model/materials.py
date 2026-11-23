#!/usr/bin/env python3
"""
Materials of the model. data/L2v2_materials.db (GATE material-database format) is the single source of material data:
GATE reads it (register_gate_materials), and the analytic parts (projector, grid, detector) take their compositions
from it too (composition, mu). Geant4 NIST materials (G4_PLEXIGLASS, ...) are used by name; their compositions for the
analytic parts are in config.NIST_COMPOSITIONS. Attenuation for the analytic parts: NIST XCOM (data/xcom_elements.json).

Usage:
  python materials.py                                         list and check the materials in the .db
  python materials.py --iodine 1.2 0.9 0.7 0.4 0.27 0.17 0.10  (re)write the iodine-disc rows (mg I/cm2, rows 0..6)
                      [--large 0]                             iodine of the two 9 mm discs (0 = plain PMMA)
                      [--phantom L2v2|L2]                     for which phantom model (own plate thickness and
                                                              material names: PMMA_I_row* for L2 v2, PMMA_I_L2_row* for L2)
"""
import argparse, json, re
import numpy as np
from config import HERE, NIST_COMPOSITIONS, PHANTOM

DB_PATH = HERE / "data" / "L2v2_materials.db"
IODINE_MARKER = "# ==== iodine discs:"       # L2 v2 section; other models: "# ==== iodine discs (<model>):"
_X = json.load(open(HERE / "data" / "xcom_elements.json"))
PMMA_W = NIST_COMPOSITIONS["G4_PLEXIGLASS"]["weights"]
PMMA_RHO = NIST_COMPOSITIONS["G4_PLEXIGLASS"]["density"]
_DENSITY_UNITS = {"g/cm3": 1.0, "mg/cm3": 1e-3, "kg/m3": 1e-3}


# ----------------------------------------------------------------------------------------------- the .db file
def read_db(path=DB_PATH):
    """materials of a GATE .db file -> {name: dict(density (g/cm3), state, weights {element: mass fraction})}.
    Follows the GATE 10 parser (opengate.geometry.materials.MaterialBuilder); +mat components are expanded into
    elements with the compositions of the .db or config.NIST_COMPOSITIONS."""
    raw, section = {}, None
    lines = [l.strip().replace("\t", " ") for l in open(path)]
    i = 0
    while i < len(lines):
        line = lines[i]; i += 1
        if not line or line.startswith("#"):
            continue
        if line.split()[0] in ("[Elements]", "[Materials]", "[Isotopes]"):
            section = line.split()[0]
            continue
        if section != "[Materials]":
            continue
        name, rest = (s.strip() for s in line.split(":"))
        parts = [p.strip() for p in rest.split(";")]
        dval, dunit = parts[0].split("=")[1].split()
        n = int(parts[1].split("=")[1])
        state = parts[2].split("=")[1].strip().lower() if len(parts) > 2 else None
        comps = []
        for _ in range(n):
            c = lines[i]; i += 1
            kind = "el" if c.startswith("+el:") else "mat" if c.startswith("+mat:") else None
            if kind is None:
                raise ValueError(f"{path}: material {name}: expected {n} component lines, got '{c}'")
            a, b = re.split("[;,]", c.split(":", 1)[1])
            cname = a.split("=")[1].strip()
            key, val = (s.strip() for s in b.split("="))
            if key != "f":
                raise ValueError(f"{path}: material {name}: use mass fractions (f=), not '{key}='")
            comps.append((kind, cname, float(val)))
        raw[name] = dict(density=float(dval) * _DENSITY_UNITS[dunit], state=state, components=comps)
    out = {}

    def expand(name):
        if name in out:
            return out[name]["weights"]
        if name in NIST_COMPOSITIONS:
            return NIST_COMPOSITIONS[name]["weights"]
        if name not in raw:
            raise KeyError(f"material {name} is neither in {path.name} nor in config.NIST_COMPOSITIONS")
        w = {}
        for kind, cname, f in raw[name]["components"]:
            sub = {cname: 1.0} if kind == "el" else expand(cname)
            for el, fe in sub.items():
                w[el] = w.get(el, 0.0) + f * fe
        out[name] = dict(density=raw[name]["density"], state=raw[name]["state"], weights=w)
        return w

    for name in raw:
        expand(name)
    return out


DB = read_db()


def check_db(db=None):
    """problems in the .db: fractions not summing to 1, elements without attenuation data, missing materials"""
    db = DB if db is None else db
    problems = []
    for name, m in db.items():
        s = sum(m["weights"].values())
        if abs(s - 1) > 1e-4:
            problems.append(f"{name}: mass fractions sum to {s:.6f}")
        for el in m["weights"]:
            if el not in _X:
                problems.append(f"{name}: no attenuation data for element {el} in data/xcom_elements.json")
    from config import SUPPORT, PADDLE, GRID, SYSTEM, BACKGROUNDS
    used = [PHANTOM["container_material"], PHANTOM["plate_material"], SUPPORT["material"], PADDLE["material"],
            GRID["interspace_material"], GRID["cover_material"], GRID["lamella_material"],
            SYSTEM["detector_cover"]["material"], "aSe"]
    used += [m for b in BACKGROUNDS.values() for m in (b["sphere_material"], b["liquid_material"])]
    for name in used:
        if name not in db and name not in NIST_COMPOSITIONS:
            problems.append(f"material {name} (used in config.py) is not defined")
    if not iodine_rows(db):
        problems.append("no iodine-disc rows for phantom %s in the .db: run python materials.py --phantom %s --iodine ..."
                        % (PHANTOM["model"], PHANTOM["model"]))
    return problems


# ----------------------------------------------------------------------------------------------- iodine discs
def iodine_rows(db=None, plate_mm=None):
    """{row (int) or 'large': dict(material, conc_g_cm3, mg_cm2)} from the iodine materials of the selected phantom
    model in the .db (PMMA_I_row*/PMMA_I_large for L2 v2, PMMA_I_L2_row*/PMMA_I_L2_large for L2)"""
    import phantoms
    db = DB if db is None else db
    if plate_mm is None:
        plate_mm = phantoms.plate()[1]
    prefix = phantoms.model()["iodine_prefix"]
    out = {}
    for name, m in db.items():
        k = re.fullmatch(re.escape(prefix) + r"(row(\d+)|large)", name)
        if not k:
            continue
        c = m["weights"].get("I", 0.0) * m["density"]
        out[int(k.group(2)) if k.group(2) else "large"] = dict(material=name, conc_g_cm3=c, mg_cm2=c * plate_mm / 10 * 1000)
    return dict(sorted(out.items(), key=lambda kv: (isinstance(kv[0], str), kv[0])))


def iodine_material(mg_cm2, plate_mm):
    """PMMA disc loaded with iodine: areal density (mg/cm2) over the plate thickness -> density and PMMA/I fractions"""
    c = mg_cm2 / (plate_mm / 10.0) / 1000.0           # g iodine per cm3
    rho = PMMA_RHO + c
    return dict(density=rho, f_pmma=PMMA_RHO / rho, f_iodine=c / rho, conc_g_cm3=c)


def iodine_marker(model=None):
    return IODINE_MARKER if (model or PHANTOM["model"]) == "L2v2" else "# ==== iodine discs (%s):" % (model or PHANTOM["model"])


def write_iodine_rows(mg_cm2_rows, large=0.0, path=DB_PATH):
    """replace the iodine section of the selected phantom model in the .db (from its marker to the next '# ====' section
    or the end); the sections of the other models are kept"""
    import phantoms
    plate_mm = phantoms.plate()[1]
    prefix, marker = phantoms.model()["iodine_prefix"], iodine_marker()
    text = open(path).read()
    if marker in text:
        i0 = text.index(marker)
        i1 = text.find("\n# ==== ", i0 + 1)
        head, tail = text[:i0], ("" if i1 < 0 else text[i1 + 1:])
    else:
        head, tail = text.rstrip() + "\n\n", ""
    cmd = ("" if PHANTOM["model"] == "L2v2" else " --phantom " + PHANTOM["model"])
    out = [head + marker + " PMMA + iodine, one material per row. Generated by:  python materials.py%s --iodine " % cmd
           + " ".join("%g" % a for a in mg_cm2_rows) + (" --large %g" % large if large else "") + " ===="]
    out.append("# Each row: the areal density (mg I/cm2) spread over the %.1f mm plate thickness, mixed into PMMA "
               "(mixture rule). Values [A]:" % plate_mm)
    out.append("# lower-bound estimates from the July 2026 plate images (REPORT §13.2); replace with the nominal values.")
    rows = [(f"{prefix}row{k}", a, f"row {k}") for k, a in enumerate(mg_cm2_rows)]
    if large > 0:
        rows.append((f"{prefix}large", large, "the two 9 mm discs"))
    for name, a, label in rows:
        m = iodine_material(a, plate_mm)
        fp = round(m["f_pmma"], 7); fi = round(1 - fp, 7)
        out.append(f"# {label}: {a:g} mg I/cm2 = {m['conc_g_cm3'] * 1000:.3f} mg I/cm3")
        out.append(f"{name}: d={m['density']:.6f} g/cm3 ; n=2 ; state=solid")
        out.append(f"+mat: name=G4_PLEXIGLASS ; f={fp:.7f}")
        out.append(f"+el: name=I ; f={fi:.7f}")
    open(path, "w").write("\n".join(out) + "\n" + ("\n" + tail if tail else ""))


# ----------------------------------------------------------------------------------------------- attenuation
def mu_rho_el(el, e_kev):
    """total mass attenuation incl. coherent scattering (cm2/g); log-log interpolation, edge-aware"""
    a = np.array(_X[el])
    e = a[:, 0] * 1000.0
    m = a[:, 6]
    e = e.copy()
    for i in range(1, len(e)):
        if e[i] <= e[i - 1]:
            e[i - 1] -= 1e-6        # put the low side of an absorption edge just below it
    o = np.argsort(e)
    return np.exp(np.interp(np.log(e_kev), np.log(e[o]), np.log(m[o])))


def composition(name):
    if name in DB:
        return DB[name]
    if name in NIST_COMPOSITIONS:
        return NIST_COMPOSITIONS[name]
    raise KeyError(f"material {name} is neither in {DB_PATH.name} nor in config.NIST_COMPOSITIONS")


def mu(name, e_kev):
    """linear attenuation coefficient (1/mm) of a named material at energies e_kev"""
    c = composition(name)
    return c["density"] * sum(w * mu_rho_el(el, e_kev) for el, w in c["weights"].items()) / 10.0


def register_gate_materials(sim):
    """make the materials of the .db available to a GATE 10 simulation (Geant4 NIST materials need nothing)"""
    problems = check_db()
    if problems:
        raise SystemExit("material database %s:\n  " % DB_PATH + "\n  ".join(problems))
    sim.volume_manager.add_material_database(str(DB_PATH))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iodine", type=float, nargs=7, metavar="MG_CM2", help="iodine per disc row 0..6 (mg I/cm2)")
    ap.add_argument("--large", type=float, default=0.0, help="iodine of the two 9 mm discs (mg I/cm2; 0 = plain PMMA)")
    from config import PHANTOM_MODELS, set_phantom
    import phantoms
    ap.add_argument("--phantom", choices=list(PHANTOM_MODELS), default=None, help="default: REVELATION_PHANTOM or L2v2")
    a = ap.parse_args()
    set_phantom(a.phantom)
    global DB
    if a.iodine:
        write_iodine_rows(a.iodine, a.large)
        DB = read_db()
        print("iodine rows written to", DB_PATH)
    e = np.array([20.0, 30.0, 40.0])
    print("%-16s %8s %-6s  %-48s  mu (1/cm) at 20 / 30 / 40 keV" % ("material", "g/cm3", "state", "mass fractions"))
    for name in list(DB) + [k for k in NIST_COMPOSITIONS if k != "G4_Pb"]:
        c = composition(name)
        wt = " ".join("%s %.4f" % (el, w) for el, w in sorted(c["weights"].items(), key=lambda kv: -kv[1]))
        print("%-16s %8.4f %-6s  %-48s  %s" % (name, c["density"], (c.get("state") or "NIST")[:6], wt[:48],
                                               " / ".join("%.3f" % (v * 10) for v in mu(name, e))))
    from config import BACKGROUNDS, DEFAULT_BACKGROUND
    print("\nbackgrounds (choose with --background; default %s):" % DEFAULT_BACKGROUND)
    for k, b in BACKGROUNDS.items():
        print("  %-14s spheres %-15s in %-13s %s" % (k, b["sphere_material"], b["liquid_material"], b["description"]))
    print("\niodine discs of phantom %s (plate %.1f mm):" % (PHANTOM["model"], phantoms.plate()[1]))
    for k, r in iodine_rows().items():
        print("  %-6s %-13s %.3f mg I/cm3 = %.3f mg I/cm2" % (k, r["material"], r["conc_g_cm3"] * 1000, r["mg_cm2"]))
    problems = check_db()
    print("\ncheck:", "OK" if not problems else "\n  " + "\n  ".join(problems))


if __name__ == "__main__":
    main()

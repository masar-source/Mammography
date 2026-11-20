#!/usr/bin/env python3
"""
GATE 10 (opengate) simulation of one exposure of the Siemens MAMMOMAT Revelation with an L2 phantom (iodine plate
in the middle). Phantom model with --phantom: L2v2 (suggested redesign, default) or L2 (the real container with the
new July 2026 plate; phantoms.py). Background with --background: PA12_MCT (PA12 + 5.5 wt% CaCO3 spheres in MCT oil,
default for L2v2) or PMMA_paraffin (PMMA spheres in liquid paraffin, the real L2 materials, default for L2). Spectra
of the measured reference set (--set: 0035, or july28high for L2). Materials in data/L2v2_materials.db.

Simulated in Geant4: focal spot, air, compression paddle, phantom, breast support. Below the support a phase-space
plane records every photon (energy, position, direction, scattered or not); a kill volume below it stops them.
The anti-scatter grid, the detector cover and the a-Se response are applied afterwards (grid_detector.py,
postprocess.py), because they act on each photon independently and can then be varied without re-running GATE.

Source: isotropic emission from a 0.3 mm disc at the focal spot into a cone (x >= 0 only: chest-wall edge of the
field), trimmed by GATE's angular acceptance to directions that intersect a target volume: the phantom envelope
(default, efficient) or a virtual aperture = the collimated field (detector area [A]). Rejected directions are emitted
with zero energy, so every event is counted and the normalisation (events per steradian) is exact. Spectrum from
spectra.py (SpekPy).

World coordinates (mm, see config.py): a-Se surface at z = 0, focal spot at (0, 0, SID), x from the chest-wall face of
the support into the phantom.

Usage:
  python revelation_gate.py --beam LE --n 2e8 --threads 8                 phantom (default background PA12_MCT)
  python revelation_gate.py --beam LE --n 2e8 --background PMMA_paraffin  phantom with PMMA spheres in paraffin oil
  python revelation_gate.py --beam LE --n 2e8 --phantom L2                 real L2 + new plate (PMMA/paraffin, July spectra)
  python revelation_gate.py --beam LE --n 2e8 --mode open                  no phantom (flat field; paddle at the same height)
  python revelation_gate.py --visu vrml_file_only                          look at the geometry first (see --visu-*)
  options: --no-spheres (oil only), --no-discs, --seed, --tag, --check-overlaps, --flip-iso
Output: output/<tag>/phsp.root, stats.txt, run_info.json (with --visu also geometry.wrl / .gdml);
default tag <beam>_<mode>[_<case>]_s<seed>, case = background (L2v2) or L2_<background> (config.case_name)
To look at the model without GATE, use preview.py.
Large photon numbers are run as several runs with different seeds; postprocess.py accumulates them into images and can
delete each phase-space file after use (an LE phantom run of 2e8 photons gives about 1 GB).
"""
import argparse, csv, json
import numpy as np
from config import (HERE, SYSTEM, SUPPORT, PADDLE, BEAMS, PHANTOM, SIM, OUT_DIR, MEASURED, BACKGROUNDS, spectra_dir,
                    add_case_args, apply_case_args, case_name)
from materials import register_gate_materials, iodine_rows
import phantoms

# ------------------------------------------------------------------------------------------------- fixed layout [A]
# collimated field at z = 0 (mm) [A]. The images show no collimator edge [M]: the open field reaches every edge of the
# stored image (238.0 x 299.0 mm), so the field is at least that large; here it slightly exceeds the detector.
FIELD_X = (0.0, 245.0)          # chest-wall edge to beyond the detector (239.36 + 3)
FIELD_Y = (-155.0, 155.0)       # detector 304.64 mm wide
APERTURE_Z = 200.0              # plane of the virtual collimator aperture (above the paddle)
PADDLE_X, PADDLE_Y = (-25.0, 240.0), (-155.0, 155.0)
SUPPORT_X, SUPPORT_Y = (0.0, 265.0), (-165.0, 165.0)
PLANE_Z = (14.0, 14.4)          # phase-space plane below the support (support bottom at 15.49)
PLANE_X, PLANE_Y = (-60.0, 300.0), (-210.0, 210.0)
GAP = 0.01                      # air gap between touching volumes (avoids coincident surfaces)

# colours for the Geant4 visualisation: RGB + alpha (alpha 0 = hidden; 0 < alpha < 1 = transparent)
COLORS = dict(support=[0.35, 0.35, 0.35, 0.8], paddle=[0.6, 0.85, 1.0, 0.25], envelope=[0.55, 0.75, 1.0, 0.15],
              filling=[1.0, 0.9, 0.5, 0.12], plate=[0.9, 0.9, 0.9, 0.5], disc=[1.0, 0.1, 0.1, 1.0],
              spheres=[1.0, 0.55, 0.1, 0.9], aperture=[0.1, 0.8, 0.1, 0.6], phsp_plane=[0.8, 0.2, 0.8, 0.3],
              absorber=[0.2, 0.2, 0.2, 0.0])
SPHERES_IN_VISU = dict(all=0.0, large=6.0, none=None)    # smallest diameter shown with --visu (mm)


def geometry_z():
    """z positions (mm) of the stacked parts"""
    H = phantoms.height()
    z_sup_top = SYSTEM["support_to_detector_mm"]
    return dict(support=(z_sup_top - SUPPORT["mm"] - GAP, z_sup_top - GAP), phantom=(z_sup_top, z_sup_top + H),
                paddle=(z_sup_top + H + GAP, z_sup_top + H + GAP + PADDLE["mm"]), H=H)


def source_cone(target="field"):
    """theta/phi range (deg) of the emission cone and its solid angle. Geant4 GPS convention (G4SPSAngDistribution::
    GenerateIsotropicFlux): direction = -(sin t cos p, sin t sin p, cos t), so theta = 0 points to -z and
    phi in [90, 270] deg gives x >= 0 (the collimated field ends at the chest-wall edge, x = 0).
    target 'field': cone over the whole collimated field; 'phantom': cone over the phantom envelope only."""
    sid = SYSTEM["sid_mm"]
    if target == "field":
        far = max(np.hypot(x, y) for x in FIELD_X for y in FIELD_Y)
        theta_max = np.degrees(np.arctan((far + 1.0) / sid)) + 0.5
    else:
        z = geometry_z()
        theta_max = np.degrees(np.arctan((phantoms.reach_mm() + 1.0) / (sid - z["phantom"][1]))) + 0.5
    phi = (90.0, 270.0)
    omega = np.radians(phi[1] - phi[0]) * (1 - np.cos(np.radians(theta_max)))
    return theta_max, phi, omega


def check_qt():
    """GATE silently skips --visu qt when opengate_core was built without Qt: stop with a clear message instead"""
    import opengate_core as g4
    info = g4.GateInfo                     # static functions: called on the class, as opengate itself does
    if not info.get_G4VIS_USE_OPENGLQT():
        raise SystemExit("This opengate build has no Geant4 Qt/OpenGL visualisation. Use --visu vrml (window via "
                         "pyvista) or --visu vrml_file_only (open the .wrl in ParaView).")
    import os
    if os.name != "nt" and not os.environ.get("DISPLAY") and not os.environ.get("WAYLAND_DISPLAY"):
        raise SystemExit("No display found (DISPLAY is empty). Run on the workstation's own desktop, or connect "
                         "with 'ssh -Y' from a machine with an X server.")
    print("Geant4 Qt visualisation available (Qt %s)" % info.get_QT_VERSION(), flush=True)


def read_g4_commands(path):
    """Geant4 commands of a .mac file, without comments and empty lines"""
    return [l.strip() for l in open(path) if l.strip() and not l.strip().startswith("#")]


def add_box(sim, name, xr, yr, zr, material, mother=None):
    mm = __import__("opengate").g4_units.mm
    b = sim.add_volume("Box", name)
    if mother:
        b.mother = mother
    b.size = [(xr[1] - xr[0]) * mm, (yr[1] - yr[0]) * mm, (zr[1] - zr[0]) * mm]
    b.translation = [(xr[0] + xr[1]) / 2 * mm, (yr[0] + yr[1]) / 2 * mm, (zr[0] + zr[1]) / 2 * mm]
    b.material = material
    return b


def add_phantom(sim, z0, spheres=True, discs=True, min_sphere_mm=0.0):
    """phantom of the selected model (phantoms.py); mesh and sphere coordinates are phantom coordinates (z from the
    bottom outer surface); the envelope is placed at world (x_offset, 0, z0).
    min_sphere_mm > 0 leaves out the smaller spheres (only for visualisation)."""
    import opengate as gate
    mm = gate.g4_units.mm
    ph = phantoms.phantom_dir()
    env = sim.add_volume("Tesselated", "envelope")
    env.file_name = str(ph / "envelope.stl"); env.origin_at_cog = False
    env.material = PHANTOM["container_material"]; env.translation = [phantoms.x_offset() * mm, 0, z0 * mm]
    fil = sim.add_volume("Tesselated", "filling")
    fil.mother = env.name; fil.file_name = str(ph / "filling.stl"); fil.origin_at_cog = False
    fil.material = PHANTOM["liquid_material"]
    pl = sim.add_volume("Tesselated", "plate")
    pl.mother = fil.name; pl.file_name = str(ph / "plate.stl"); pl.origin_at_cog = False
    pl.material = PHANTOM["plate_material"]
    env.color, fil.color, pl.color = COLORS["envelope"], COLORS["filling"], COLORS["plate"]
    info = dict(spheres=0, discs=0)
    if spheres:
        S = np.load(ph / "spheres.npz")
        for k, d in enumerate(PHANTOM["sphere_diameters_mm"]):
            c = S["centres"][S["size_class"] == k]
            if len(c) == 0 or d < min_sphere_mm:
                continue
            sp = sim.add_volume("Sphere", f"spheres_{str(d).replace('.', 'p')}")
            sp.mother = fil.name; sp.rmin = 0; sp.rmax = d / 2 * mm
            sp.material = PHANTOM["sphere_material"]; sp.color = COLORS["spheres"]
            sp.translation = [[float(x) * mm, float(y) * mm, float(z) * mm] for x, y, z in c]
            info["spheres"] += len(c)
    if discs:
        zc, t = phantoms.plate()
        rows = iodine_rows()                              # iodine materials of data/L2v2_materials.db
        for d in phantoms.discs():
            if d["key"] not in rows:
                continue                                  # no iodine material: the disc is plate material (PMMA)
            dv = sim.add_volume("Tubs", d["name"])
            dv.mother = pl.name; dv.rmin = 0; dv.rmax = d["diam"] / 2 * mm
            dv.dz = (t / 2 - 0.001) * mm; dv.material = rows[d["key"]]["material"]; dv.color = COLORS["disc"]
            dv.translation = [d["x"] * mm, d["y"] * mm, zc * mm]
            info["discs"] += 1
    return info


def build(a):
    import opengate as gate
    u = gate.g4_units
    mm, keV, deg = u.mm, u.keV, u.deg
    sim = gate.Simulation()
    out = OUT_DIR / a.tag
    sim.output_dir = str(out)
    sim.random_seed = a.seed
    sim.number_of_threads = a.threads
    sim.check_volumes_overlap = a.check_overlaps      # off by default: 60 000 spheres (checked by sphere_packing.py)
    sim.progress_bar = True
    visu = getattr(a, "visu", "off") != "off"
    if visu:
        # Geant4 visualisation of the geometry as built, with a few photon tracks; single thread, few photons
        sim.visu = True
        sim.visu_type = a.visu
        if a.visu in ("vrml", "vrml_file_only", "gdml", "gdml_file_only"):
            out.mkdir(parents=True, exist_ok=True)
            # a str, not a Path: opengate 10.1.1 passes it to os.environ for the VRML driver
            sim.visu_filename = str(out / ("geometry.wrl" if a.visu.startswith("vrml") else "geometry.gdml"))
        if a.visu == "qt":
            check_qt()
            sim.visu_commands = read_g4_commands(HERE / "data" / "visu_qt.mac")
        sim.number_of_threads = 1
    sim.world.size = [1.2 * u.m, 1.2 * u.m, 1.6 * u.m]
    sim.world.material = "G4_AIR"
    register_gate_materials(sim)                      # data/L2v2_materials.db

    z = geometry_z()
    add_box(sim, "support", SUPPORT_X, SUPPORT_Y, z["support"], SUPPORT["material"])
    add_box(sim, "paddle", PADDLE_X, PADDLE_Y, z["paddle"], PADDLE["material"])
    ph = dict(spheres=0, discs=0)
    if a.mode == "phantom":
        smin = SPHERES_IN_VISU[a.visu_spheres] if visu else 0.0
        ph = add_phantom(sim, z["phantom"][0], spheres=not a.no_spheres and smin is not None, discs=not a.no_discs,
                         min_sphere_mm=smin or 0.0)
    sid = SYSTEM["sid_mm"]
    k = (sid - APERTURE_Z) / sid                     # aperture = field scaled to the aperture plane
    aperture = add_box(sim, "aperture", (FIELD_X[0] * k, FIELD_X[1] * k), (FIELD_Y[0] * k, FIELD_Y[1] * k),
                       (APERTURE_Z - 0.05, APERTURE_Z + 0.05), "G4_AIR")
    plane = add_box(sim, "phsp_plane", PLANE_X, PLANE_Y, PLANE_Z, "G4_AIR")
    absorber = add_box(sim, "absorber", (-500, 500), (-500, 500), (-100, PLANE_Z[0] - GAP), "G4_AIR")
    for v in ("support", "paddle", "aperture", "phsp_plane", "absorber"):
        sim.volume_manager.get_volume(v).color = COLORS[v]

    # physics: Livermore photon models (Rayleigh, photoelectric, Compton with binding), fluorescence on
    sim.physics_manager.physics_list_name = SIM["physics_list"]
    sim.physics_manager.em_parameters.fluo = SIM["fluorescence"]
    sim.physics_manager.set_production_cut(sim.world.name, "all", SIM["production_cut_mm"] * mm)
    if a.kill_electrons:
        # electrons (< 49 keV, range < 0.05 mm in plastic) deposit locally and do not reach the image; killing them at
        # creation saves time. Bremsstrahlung from them (yield ~1e-4) is neglected.
        sim.physics_manager.user_limits_particles = ["electron"]
        for v in ("support", "paddle") + (("envelope",) if a.mode == "phantom" else ()):
            sim.physics_manager.set_min_ekine(v, 1 * u.MeV)

    # source
    spec = np.load(spectra_dir() / f"{a.beam}.npz")
    target = "field" if a.mode == "open" else a.target
    theta_max, phi, omega = source_cone(target)
    src = sim.add_source("GenericSource", "tube")
    src.particle = "gamma"
    # GATE 10.1.1 generates source.n photons in EACH thread (verified on the workstation: --n 1e6 with 8 threads gave
    # 8e6 events), so --n (total) is split over the threads. The normalisation uses the counted events anyway.
    n_total = int(min(a.n, a.visu_n)) if visu else int(a.n)
    src.n = int(np.ceil(n_total / sim.number_of_threads))
    src.position.type = "disc"
    src.position.radius = SYSTEM["focal_spot_mm"] / 2 * mm
    src.position.translation = [0, 0, sid * mm]
    src.energy.type = "spectrum_histogram"
    src.energy.spectrum_energy_bin_edges = [float(e) * keV for e in spec["edges_keV"]]
    src.energy.spectrum_weights = [float(w) for w in spec["weights"]]
    src.direction.type = "iso"
    th = (180.0 - theta_max, 180.0) if a.flip_iso else (0.0, theta_max)
    src.direction.theta = [th[0] * deg, th[1] * deg]
    src.direction.phi = [phi[0] * deg, phi[1] * deg]
    aa = src.direction.angular_acceptance
    aa.policy = "Rejection"
    aa.skip_policy = "ZeroEnergy"                    # rejected directions -> zero-energy events (counted, not tracked)
    aa.target_volumes = [aperture.name if target == "field" else "envelope"]
    aa.enable_intersection_check = True

    # actors
    stats = sim.add_actor("SimulationStatisticsActor", "stats")
    stats.output_filename = "stats.txt"
    phsp = sim.add_actor("PhaseSpaceActor", "phsp")
    phsp.attached_to = plane.name
    phsp.attributes = SIM["phsp_attributes"]
    phsp.steps_to_store = "entering"
    phsp.output_filename = "phsp.root"
    kill = sim.add_actor("KillActor", "kill")
    kill.attached_to = absorber.name

    b = BEAMS[a.beam] if a.beam in BEAMS else dict(kvp=float(a.beam[2:]), filter=BEAMS["LE"]["filter"], mAs=None)
    info = dict(beam=a.beam, kvp=b["kvp"], filter=b["filter"], mode=a.mode, target=target,
                background=PHANTOM["background"] if a.mode == "phantom" else None,
                phantom_model=PHANTOM["model"] if a.mode == "phantom" else None, case=case_name(),
                measured_set=MEASURED["set"], spectra_dir=str(spectra_dir()),
                sphere_material=PHANTOM["sphere_material"], liquid_material=PHANTOM["liquid_material"],
                spheres=not a.no_spheres, discs=not a.no_discs, n=int(src.n) * sim.number_of_threads,
                n_per_thread=int(src.n), seed=a.seed,
                threads=sim.number_of_threads, visu=getattr(a, "visu", "off"),
                theta_deg=list(th), phi_deg=list(phi), omega_cone_sr=omega, field_x_mm=FIELD_X, field_y_mm=FIELD_Y,
                plane_z_mm=PLANE_Z[1], z_mm=z, phantom=ph, physics=SIM["physics_list"], fluo=SIM["fluorescence"],
                spectrum=dict(fluence_per_mAs_1m=float(spec["fluence_per_mAs_1m"]), mAs=b["mAs"],
                              extra_al_mm=float(spec["extra_al_mm"]) if "extra_al_mm" in spec else 0.0))
    return sim, stats, info, out


def sanity_check(path, n=200000, min_primaries=5):
    """Unscattered primary photons must come straight from the focal spot: traced back along their direction they must
    pass within 1 mm of it, and they must arrive at x >= 0 (chest-wall edge of the field). Scattered photons may go
    anywhere and are not tested. Catches a wrong source position or direction convention, which would otherwise give
    a silently empty or misplaced phase space. Too few unscattered photons: no verdict (ok = None)."""
    keys = ["PrePosition_X", "PrePosition_Y", "PrePosition_Z", "PreDirection_X", "PreDirection_Y", "PreDirection_Z",
            "UnscatteredPrimaryFlag"]
    try:
        import uproot
        f = uproot.open(path)
        t = f[[k for k in f.keys() if k.split(";")[0] == "phsp"][0]]
        b = t.arrays(keys, entry_stop=n, library="np")
    except Exception as e:                                         # noqa: BLE001
        print("WARNING: sanity check could not read the phase space:", e)
        return dict(ok=False, error=str(e))
    nin = len(b["PrePosition_X"])
    prim = b["UnscatteredPrimaryFlag"] == 1
    npr = int(prim.sum())
    out = dict(entries_checked=int(nin), unscattered_checked=npr)
    if npr < min_primaries:
        print("sanity check: only %d unscattered photons at the plane, too few to judge" % npr)
        return dict(ok=None if nin else False, **out)
    x, y, z = b["PrePosition_X"][prim], b["PrePosition_Y"][prim], b["PrePosition_Z"][prim]
    ux, uy, uz = b["PreDirection_X"][prim], b["PreDirection_Y"][prim], b["PreDirection_Z"][prim]
    back = (SYSTEM["sid_mm"] - z) / np.maximum(-uz, 1e-9)              # path length back to the focal-spot plane
    r_spot = np.hypot(x - ux * back, y - uy * back)                    # distance from the focal-spot centre there
    out.update(fraction_from_focal_spot=float((r_spot < 1.0).mean()), fraction_in_field=float((x > -1.0).mean()),
               fraction_downward=float((b["PreDirection_Z"] < 0).mean()), max_backprojected_offset_mm=float(r_spot.max()))
    ok = out["fraction_from_focal_spot"] > 0.99 and out["fraction_in_field"] > 0.99 and out["fraction_downward"] > 0.99
    print("sanity check %s: %d unscattered photons, %.1f %% trace back to the focal spot (max offset %.2f mm), "
          "%.1f %% in the field" % ("OK" if ok else "FAILED", npr, 100 * out["fraction_from_focal_spot"],
                                     out["max_backprojected_offset_mm"], 100 * out["fraction_in_field"]))
    if not ok:
        print("WARNING: unscattered photons do not come from the focal spot into the field. Check the source "
              "direction convention (--flip-iso) before using this run.")
    return dict(ok=bool(ok), **out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--beam", choices=list(BEAMS) + [f"LE{k}" for k in (30, 32, 34)], default="LE")
    ap.add_argument("--n", type=float, default=1e7, help="photons emitted into the cone, in total (split over threads)")
    ap.add_argument("--mode", choices=["phantom", "open"], default="phantom")
    add_case_args(ap)
    ap.add_argument("--no-spheres", action="store_true")
    ap.add_argument("--no-discs", action="store_true")
    ap.add_argument("--threads", type=int, default=SIM["threads"])
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--tag", default=None)
    ap.add_argument("--target", choices=["phantom", "field"], default="phantom",
                    help="phantom: emit only towards the phantom envelope (default; photons that miss it are not "
                         "simulated, which drops paddle/support scatter from the open field); field: whole field")
    ap.add_argument("--check-overlaps", action="store_true")
    ap.add_argument("--kill-electrons", action=argparse.BooleanOptionalAction, default=True)
    ap.add_argument("--visu", choices=["off", "qt", "vrml", "vrml_file_only", "gdml_file_only"], default="off",
                    help="show the geometry as GATE builds it, with a few photon tracks, instead of a production run: "
                         "qt = Geant4 Qt/OpenGL viewer (needs a display; commands in data/visu_qt.mac); vrml = window via pyvista "
                         "(pip install pyvista); *_file_only = write output/<tag>/geometry.wrl or .gdml only")
    ap.add_argument("--visu-spheres", choices=list(SPHERES_IN_VISU), default="large",
                    help="spheres shown with --visu: large (>= 6.35 mm, ~1 200 spheres; default), all (~60 000, slow), none")
    ap.add_argument("--visu-n", type=float, default=200, help="photons tracked with --visu")
    ap.add_argument("--flip-iso", action="store_true",
                    help="use if GATE stops with no accepted direction (opposite iso sign convention)")
    a = ap.parse_args()
    case = apply_case_args(a)
    bg = PHANTOM["background"]
    print("phantom %s, background %s: %s; spectra %s (measured set %s)" % (PHANTOM["model"], bg, BACKGROUNDS[bg]["description"],
                                                                          spectra_dir(), MEASURED["set"]), flush=True)
    if a.tag is None:
        a.tag = f"{a.beam}_{a.mode}" + (f"_{case}" if a.mode == "phantom" else "") + \
            ("_nospheres" if a.no_spheres else "") + ("_nodiscs" if a.no_discs else "") + \
            (f"_visu_{a.visu_spheres}" if a.visu != "off" else f"_s{a.seed}")
    sim, stats, info, out = build(a)
    sim.run()
    info["events"] = int(stats.counts.events)
    info["sanity"] = sanity_check(out / "phsp.root")
    out.mkdir(parents=True, exist_ok=True)
    json.dump(info, open(out / "run_info.json", "w"), indent=1)
    print(stats)
    print("written:", out)


if __name__ == "__main__":
    main()

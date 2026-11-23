"""A room-transport extension for the existing Revelation GATE 10 project."""

from contextlib import nullcontext
import hashlib
import math
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np

from acquisition import dbt_source
from air_kerma import mu_en, read_cell
from patient_cc import PHANTOM_MODES, make_patient, projection_for_scatter
from room_geometry import (
    box,
    materials,
    reference,
    rot,
    scanner_installation,
    screen_panels,
    survey_points,
    track_segments,
    walls,
)


def check_survey_clearance(sim, sites, intrinsic, rotation, shift):
    """Reject cells touching a transport volume's conservative room-frame bounding box."""
    bounds = []
    for name, volume in sim.volume_manager.volumes.items():
        if name == "world" or volume.mother != "world":
            continue
        if name not in intrinsic and name != "rp_torso":
            continue
        if name == "rp_patient_breast":
            half = volume.size_pix * volume.spacing / 2
            lo, hi = -half, half
        elif type(volume).__name__ == "BoxVolume":
            half = np.asarray(volume.size) / 2
            lo, hi = -half, half
        else:
            lo, hi = (np.array([v[i] for i in range(3)]) for v in volume.bounding_limits)
        frame = rotation if name in intrinsic else np.eye(3)
        axes = frame @ volume.rotation
        centre = frame @ np.asarray(volume.translation) + axes @ ((lo + hi) / 2)
        if name in intrinsic:
            centre += shift
        bounds.append((name, centre, np.abs(axes) @ ((hi - lo) / 2)))
    for label, site in sites.items():
        point, half = np.array(site["xyz_mm"]), np.array(site["size_mm"]) / 2
        for name, centre, extent in bounds:
            if np.all(np.abs(point - centre) <= half + extent):
                raise ValueError(f"Survey cell {label} touches the bounds of {name}. "
                                 "Review placement for this phantom and view before transport.")


def receptor(sim, native, geometry, options, focus, grid_enabled):
    """Transport counterpart of the grid_detector.py contact detector response."""
    from grid_detector import DET_COVER, GRID_TOP, Z_GRID_BOT, Z_GRID_MID, Z_GRID_TOP

    sysconf, grid = native.SYSTEM, native.GRID
    width_y, length_x = sysconf["detector_active_mm"]

    def layer(name, material, z0, z1):
        return box(
            sim,
            name,
            [length_x, width_y, z1 - z0],
            [length_x / 2, 0, (z0 + z1) / 2],
            material,
        )

    if grid["enabled"] and grid_enabled:
        c = grid["cover_mm_each"]
        layer("rp_grid_top", grid["cover_material"], Z_GRID_TOP, GRID_TOP)
        layer("rp_grid_bottom", grid["cover_material"], Z_GRID_BOT - c, Z_GRID_BOT)
        core = layer("rp_grid", grid["interspace_material"], Z_GRID_BOT, Z_GRID_TOP)
        pitch, lamella = 10 / grid["lines_per_cm"], grid["lamella_um"] / 1000
        h = Z_GRID_TOP - Z_GRID_BOT
        focal_distance = grid["focus_mm"] - Z_GRID_MID
        # Omit partial strips at the two lateral boundaries of the grid.
        outermost = int(
            (width_y / 2 - lamella / 2) / (pitch * (1 + h / (2 * focal_distance)))
        )
        centres = np.arange(-outermost, outermost + 1) * pitch
        for i, group in enumerate(
            np.array_split(centres, int(options["grid_slope_groups"]))
        ):
            if not len(group):
                continue
            y = float(np.mean(group))
            strip = sim.add_volume("Trap", f"rp_lamellae_{i:02d}")
            strip.mother = core.name
            strip.material = grid["lamella_material"]
            strip.dz, strip.dy1, strip.dy2 = h / 2 - 0.00001, lamella / 2, lamella / 2
            strip.dx1 = strip.dx2 = strip.dx3 = strip.dx4 = length_x / 2 - 0.00001
            strip.alp1 = strip.alp2 = 0
            strip.theta = math.atan(abs(y) / focal_distance)
            strip.phi = -math.copysign(math.pi / 2, y)
            strip.translation = [[0, float(v), 0] for v in group]
    layer("rp_detector_cover", sysconf["detector_cover"]["material"], *DET_COVER)
    t = sysconf["a_se_um"] / 1000
    layer("rp_selenium", "G4_Se", -t, 0)
    cx = options["housing_centre_x_mm"]
    wx, wy = options["backstop_size_xy_mm"]
    frame = options.get("image_receptor")
    if frame is not None:
        front = -frame["chest_wall_setback_mm"]
        rear = front + frame["depth_mm"]
        top = sysconf["support_to_detector_mm"]
        bottom = top - frame["height_mm"]
        bottom_layer = -t - options["backing_al_mm"] - options["backstop_pb_mm"]
        if not (wx > 0 and wy > 0 and front < cx - wx / 2
                and cx + wx / 2 < rear and wy < frame["width_mm"]
                and bottom < bottom_layer <= -t < top):
            raise ValueError("The modelled detector backing and backstop must fit inside the declared image receptor.")
    depth = -t
    for name, mat, thick in (
        ("backing", "G4_Al", options["backing_al_mm"]),
        ("backstop", "G4_Pb", options["backstop_pb_mm"]),
    ):
        if thick:
            box(sim, "rp_" + name, [wx, wy, thick], [cx, 0, depth - thick / 2], mat)
        depth -= thick
    # Geometric interception of the full nominal rectangle, projected to the back of the assembly.
    k = (focus[2] - depth) / focus[2]
    margin = min(
        min(wx / 2 - abs(k * x - cx), wy / 2 - abs(focus[1] + k * (y - focus[1])))
        for x in geometry.FIELD_X
        for y in geometry.FIELD_Y
    )
    return float(margin)


def alternate_phantom(sim, geometry, choice):
    import opengate as gate

    if choice == "open":
        return
    thickness = {"pmma20": 20.0, "pmma45": 45.0, "pmma70": 70.0, "breast45": 45.0}[
        choice
    ]
    base = geometry.geometry_z()["phantom"][0]
    obj = sim.add_volume("Tubs", "rp_compressed_phantom")
    obj.material = "BreastICRU44" if choice == "breast45" else "G4_PLEXIGLASS"
    obj.rmin, obj.rmax, obj.dz = 0, 120, thickness / 2
    obj.sphi, obj.dphi = -90 * gate.g4_units.deg, 180 * gate.g4_units.deg
    obj.translation = [0, 0, base + thickness / 2]
    paddle = sim.volume_manager.get_volume("paddle")
    centre_z = base + thickness + geometry.GAP + geometry.PADDLE["mm"] / 2
    paddle.translation = [*paddle.translation[:2], centre_z]


def source_scale(fluence, solid_angle, primaries):
    # Native analysis.py uses the fluence at 1 m in photons / (cm2 mAs).
    # With ZeroEnergy rejection the denominator includes attempts outside the rectangular aperture.
    return fluence * 10000.0 * solid_angle / primaries


def leakage_source(sim, native, cfg, beam, n):
    """Equivalent isotropic emerging housing field, supplied independently for each beam."""
    import opengate as gate
    from scipy.integrate import quad

    data = cfg["planning"]["leakage_model"][beam]
    edges = np.asarray(data["edges_keV"], float) / 1000
    weights = np.asarray(data["photon_bin_fractions"], float)
    if (
        len(edges) != len(weights) + 1
        or len(weights) == 0
        or np.any(np.diff(edges) <= 0)
        or edges[0] < 0.001
        or edges[-1] > native.BEAMS[beam]["kvp"] / 1000 + 1e-10
        or not np.all(np.isfinite(weights))
        or np.any(weights < 0)
        or not np.isclose(weights.sum(), 1)
        or not math.isfinite(data["unattenuated_collision_kerma_uGy_per_mAs_at_1m"])
        or data["unattenuated_collision_kerma_uGy_per_mAs_at_1m"] <= 0
    ):
        raise ValueError(
            f"Invalid emerging leakage spectrum or normalization for {beam}"
        )
    response = 0.0
    for lo, hi, frac in zip(edges[:-1], edges[1:], weights):
        points = [0.00320290] if lo < 0.00320290 < hi else []
        response += (
            frac
            * quad(lambda e: e * float(mu_en(e)), lo, hi, points=points)[0]
            / (hi - lo)
        )
    # At 1 m: N/(4*pi m2); cm2/g to m2/kg and MeV to J, then to microGy.
    n_per_mAs = (
        data["unattenuated_collision_kerma_uGy_per_mAs_at_1m"]
        * (4 * math.pi)
        / (response * 0.1 * 1.602176634e-13 * 1e6)
    )
    # Removing the native source applies to this leakage-only replicate.
    sim.source_manager.remove_source("tube")
    q = sim.add_source("GenericSource", "rp_housing_leakage")
    q.particle, q.n = "gamma", int(n)
    q.position.type = "point"
    q.position.translation = [0, 0, native.SYSTEM["sid_mm"]]
    q.direction.type = "iso"
    q.energy.type = "spectrum_histogram"
    q.energy.spectrum_energy_bin_edges = list(edges * gate.g4_units.MeV)
    q.energy.spectrum_weights = list(weights)
    return n_per_mAs / n


def run_case(
    root,
    cfg,
    out,
    layout_name,
    beam,
    view,
    phantom,
    component,
    n,
    seed,
    overlaps=False,
    keep_tracks=False,
):
    import uproot
    from importlib.metadata import version

    root = Path(root).resolve()
    if not (root / "revelation_gate.py").is_file():
        raise FileNotFoundError(
            f"Expected the existing model at {root / 'revelation_gate.py'}"
        )
    sys.path.insert(0, str(root))
    import config as native
    import revelation_gate as geometry

    native.set_phantom("L2")
    native.set_background("PMMA_paraffin")
    native.set_measured("july28high")
    out = Path(out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    args = SimpleNamespace(
        mode="phantom" if phantom == "L2" else "open",
        beam="LE" if beam == "DBT" else beam,
        seed=seed,
        n=n,
        threads=1,
        check_overlaps=overlaps,
        tag=str(out),
        target="field",
        flip_iso=False,
        no_spheres=False,
        no_discs=False,
        kill_electrons=False,
        visu="off",
    )
    sim, stats, info, _ = geometry.build(args)
    sim.progress_bar = False
    for actor in ("phsp", "kill"):
        sim.actor_manager.remove_actor(actor)
    for volume in ("absorber", "phsp_plane"):
        sim.volume_manager.remove_volume(volume)
    materials(sim, cfg)
    shot = cfg.get("acquisition", {}).get("shots", {}).get(view, {})
    focus = np.array([0., 0., native.SYSTEM["sid_mm"]])
    dbt_spec = cfg.get("acquisition", {}).get("spectrum_file")
    if beam == "DBT":
        focus = dbt_source(sim, geometry, info, focus[2], native.SYSTEM["support_to_detector_mm"] + 30, shot["tube_offset_deg"], dbt_spec)
    patient_anatomy = None
    patient_projection = None
    if phantom in PHANTOM_MODES:
        patient_projection = projection_for_scatter(shot.get("base_view", view), cfg["views_deg"][view], custom=bool(shot))
        patient_anatomy = make_patient(sim, geometry, phantom)
    elif phantom != "L2":
        alternate_phantom(sim, geometry, phantom)
    primary_margin = receptor(sim, native, geometry, cfg["detector_extension"], focus, beam != "DBT")

    clearance = cfg["room"]["breast_reference_local_mm"]
    layout = cfg["layouts"][layout_name]
    yaw, displacement = scanner_installation(layout)
    # The floor placement carries the scanner and patient. Only the arm tilts for MLO.
    # Transport stays in the native scanner frame, so inverse-transform room volumes.
    R = yaw @ rot("x", cfg["views_deg"][view])
    shift = reference(cfg, layout) - R @ np.asarray(clearance)
    monitor = box(sim, "rp_beam_monitor", [70, 80, 1], [60, focus[1] * 135 / focus[2], 135], "G4_AIR")
    intrinsic = set(sim.volume_manager.volumes)
    walls(sim, cfg, layout)
    screen_panels(sim, cfg, layout)
    if cfg["room"]["torso"] and phantom != "open":
        box(
            sim,
            "rp_torso",
            [150, 300, 300],
            yaw @ (np.asarray(cfg["room"]["detector_chestwall_mm"]) + [-115, 0, 130]) + displacement,
            "G4_TISSUE_SOFT_ICRP",
            rotation=yaw,
        )
    sites = survey_points(cfg, layout)
    if "survey_geometry" in layout:
        check_survey_clearance(sim, sites, intrinsic, R, shift)
    scorers = []
    for site, p in sites.items():
        p["breast_reference_to_cell_m"] = float(np.linalg.norm(np.asarray(p["xyz_mm"]) - reference(cfg, layout)) / 1000)
        volume = box(sim, "rp_cell_" + site, p["size_mm"], p["xyz_mm"], "G4_AIR")
        scorers.append(track_segments(sim, site, volume))
    # Keep the original world for preset rooms. Include the full rotated extent
    # of every room box and survey cell when a larger geometry needs more space.
    world_half_mm = np.full(3, 8000.0)
    for name, v in sim.volume_manager.volumes.items():
        if name not in intrinsic and v.mother == "world":
            centre = R.T @ (np.asarray(v.translation) - shift)
            rotation = R.T @ v.rotation
            v.translation = list(centre)
            v.rotation = rotation
            farthest = np.abs(centre) + np.abs(rotation) @ (np.asarray(v.size) / 2)
            world_half_mm = np.maximum(world_half_mm, farthest + 1.0)
    sim.world.size = (2 * world_half_mm).tolist()
    sites["beam_monitor"] = {
        "group": "beam",
        "xyz_mm": list(shift + R @ np.asarray(monitor.translation)),
        "size_mm": list(monitor.size),
    }
    scorers.append(track_segments(sim, "beam_monitor", monitor))
    spec = Path(dbt_spec) if dbt_spec else native.spectra_dir() / f"{args.beam}.npz"
    if component == "field":
        weight = source_scale(
            info["spectrum"]["fluence_per_mAs_1m"], info["omega_cone_sr"], n
        )
    else:
        weight = leakage_source(sim, native, cfg, beam, n)
    sim.run()
    count = int(stats.counts.events)
    if count != n:
        raise RuntimeError(
            f"Native single-thread run counted {count} histories, expected {n}."
        )
    tracks = out / "tracks.root"
    if not tracks.exists() and any(a.total_number_of_entries for a in scorers):
        raise RuntimeError("GATE recorded tracks but the output is missing.")
    with uproot.open(tracks) if tracks.exists() else nullcontext({}) as fp:
        for name, p in sites.items():
            key = "tracks_" + name
            p.update(
                read_cell(
                    fp[key] if key in fp else None, float(np.prod(p["size_mm"])), weight
                )
            )
    if component == "field" and sites["beam_monitor"]["photon_steps"] == 0:
        raise RuntimeError(
            "No photons reached the forward beam monitor. Check the native source and aperture or increase histories."
        )
    if tracks.exists() and not keep_tracks:
        tracks.unlink()
    return dict(
        model="Siemens MAMMOMAT Revelation",
        gate_version=version("opengate"),
        layout=layout_name,
        beam=beam,
        view=view,
        view_deg=cfg["views_deg"][view],
        tube_offset_deg=shot.get("tube_offset_deg", 0), projection_number=shot.get("projection_number"),
        focal_spot_native_mm=focus.tolist(), focal_spot_room_mm=(shift + R @ focus).tolist(),
        grid_present=beam != "DBT" and native.GRID["enabled"],
        torso_present=phantom != "open" and cfg["room"]["torso"], door_open=cfg["room"].get("door_open", False),
        spectrum_model=cfg.get("acquisition", {}).get("spectrum", beam),
        scanner_to_room=dict(rotation=R.tolist(), translation_mm=shift.tolist()),
        room_inside_mm=list(cfg["room"]["inside_mm"]),
        world_size_mm=list(sim.world.size),
        phantom=phantom,
        patient_anatomy=patient_anatomy,
        patient_projection=patient_projection,
        component=component,
        seed=seed,
        events=count,
        shape_sha256=hashlib.sha256(spec.read_bytes()).hexdigest(),
        field=info["target"],
        source_to_detector_mm=native.SYSTEM["sid_mm"],  # nominal perpendicular SID at zero
        focal_spot_to_detector_plane_mm=float(focus[2]),
        dbt_pivot_z_mm=native.SYSTEM["support_to_detector_mm"] + 30 if beam == "DBT" else None,
        detector_extension=cfg["detector_extension"],
        grid_slope_groups=cfg["detector_extension"]["grid_slope_groups"],
        nominal_primary_margin_mm=primary_margin,
        cone_sr=info["omega_cone_sr"],
        native_fluence_cm2_mAs_1m=info["spectrum"]["fluence_per_mAs_1m"],
        per_history_photons_per_mAs=weight,
        native_phantom_counts=info["phantom"],
        sites=sites,
    )

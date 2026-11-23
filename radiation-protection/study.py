#!/usr/bin/env python3
"""Run and compare room layouts around an existing Siemens MAMMOMAT Revelation GATE model."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import uuid
from importlib.metadata import version

import numpy as np
from scipy.stats import t
from acquisition import configure, read_spectrum
from layout import csv_file, fits, footprint, preview
from room_geometry import booth, reference, scanner_installation, survey_points, screen_settings, validate_room_geometry
from patient_cc import MANIFEST, PHANTOM_MODES, image_path, projection_for_scatter

HERE = Path(__file__).resolve().parent
PHANTOMS = ("L2", "pmma20", "pmma45", "pmma70", "breast45", "open", *PHANTOM_MODES)


def read(path):
    return json.loads(Path(path).read_text())


def write(path, data):
    path = Path(path)
    temp = path.with_name(path.name + "." + uuid.uuid4().hex + ".tmp")
    temp.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    os.replace(temp, path)


def digest(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True).encode()).hexdigest()


def inputs(scanner, phantom, cfg):
    paths = list(HERE.glob("*.py"))
    paths += [
        scanner / p
        for p in (
            "revelation_gate.py",
            "config.py",
            "materials.py",
            "phantoms.py",
            "grid_detector.py",
            "data/L2v2_materials.db",
            "data/xcom_elements.json",
            "data/measured_july28high.json",
        )
    ]
    spectra = Path(
        os.environ.get("REVELATION_SPECTRA_DIR", scanner / "spectra_july28high")
    )
    paths += [spectra / (b + ".npz") for b in ("LE", "HE")]
    if cfg.get("acquisition", {}).get("spectrum_file"):
        paths.append(Path(cfg["acquisition"]["spectrum_file"]))
    if phantom in PHANTOM_MODES:
        paths += [MANIFEST, image_path(phantom)]
    if phantom == "L2":
        paths += [
            scanner / "phantom_L2" / p
            for p in ("envelope.stl", "filling.stl", "plate.stl", "spheres.npz")
        ]
        paths += [scanner.parent / "L2_container_v2" / "disc_layout_measured.json"]
    hashes = {}
    for f in paths:
        if not f.is_file():
            raise FileNotFoundError(f"Expected model input {f}")
        hashes[str(f)] = hashlib.sha256(f.read_bytes()).hexdigest()
    return dict(files=hashes, configuration=cfg, gate_version=version("opengate"))


def validate(cfg, work, layouts):
    p = cfg["planning"]
    beam_keys = set(work)
    if not work or not beam_keys.issubset({"LE", "HE", "DBT"}):
        raise ValueError("Supply an LE or HE workload in mAs/week.")
    for b, views in work.items():
        if (
            not views
            or not set(views).issubset(cfg["views_deg"])
            or any(not math.isfinite(v) or v <= 0 for v in views.values())
        ):
            raise ValueError(
                "Each selected view needs a positive workload in mAs/week."
            )
        if b != "DBT" and (not math.isfinite(p["source_scale"][b]) or p["source_scale"][b] <= 0):
            raise ValueError("Beam fluence scale must be positive.")
    if any(not 0 <= o <= 1 for o in p["occupancy"].values()) or any(
        not math.isfinite(g) or g <= 0 for g in p["goal_uGy_week"].values()
    ):
        raise ValueError(
            "Check area occupancies (0–1) and positive weekly air-kerma goals."
        )
    if any(k not in p["occupancy"] for k in p["goal_uGy_week"]):
        raise ValueError("Each area with a goal must have an occupancy.")
    leak = p["leakage_model"]
    if "acquisition" not in cfg and leak is not None and not beam_keys.issubset(leak):
        raise ValueError(
            "Supply a leakage model for every beam in this workload or set leakage_model to null."
        )
    if any(name not in cfg["layouts"] for name in layouts):
        raise ValueError("A selected layout is missing from study.json.")
    for name in layouts:
        validate_room_geometry(cfg, cfg["layouts"][name])
        screen_settings(cfg, cfg["layouts"][name])
        scanner_installation(cfg["layouts"][name])
    if not 0 < p["max_relative_se"] <= 1 or p["min_independent_batches"] < 3:
        raise ValueError("Check the statistical criteria in study.json.")
    det = cfg["detector_extension"]
    if (
        det["grid_slope_groups"] < 1
        or min(det["backing_al_mm"], det["backstop_pb_mm"]) < 0
    ):
        raise ValueError("Check the receptor extension parameters.")


def execute_one(task):
    out = Path(task["out"])
    out.mkdir(parents=True, exist_ok=True)
    spec = out / "task.json"
    result = out / "case.json"
    if (
        result.exists()
        and spec.exists()
        and read(spec) == task
        and read(result).get("task_sha256") == digest(task)
    ):
        return read(result)
    write(spec, task)
    with (out / "gate.log").open("w") as log:
        process = subprocess.run(
            [sys.executable, str(HERE / "study.py"), "_worker", str(spec)],
            stdout=log,
            stderr=subprocess.STDOUT,
        )
    if process.returncode:
        raise RuntimeError(f"GATE failed. See {out / 'gate.log'}")
    return read(result)


def worker(spec):
    from transport import run_case

    v = read(spec)
    cfg = read(v["config_path"])
    if digest(cfg) != v["config_digest"]:
        raise ValueError("Study configuration changed while jobs were running.")
    result = run_case(
        v["scanner"],
        cfg,
        v["out"],
        v["layout"],
        v["beam"],
        v["view"],
        v["phantom"],
        v["component"],
        v["n"],
        v["seed"],
        v["check_overlaps"],
        v["keep_tracks"],
    )
    result.update(task_sha256=digest(v), batch=v["batch"])
    write(Path(v["out"]) / "case.json", result)


def sightlines(cfg, layout):
    if not layout["screen_pb_mm"]:
        return True
    _, sc, n, tan = booth(cfg, layout)
    ref = reference(cfg, layout)
    sites = survey_points(cfg, layout)
    scr = screen_settings(cfg, layout)
    for name in (
        "operator_1700",
        "trainee_1700",
        "left_reach_1700",
        "outer_reach_1700",
    ):
        ray = np.asarray(sites[name]["xyz_mm"]) - ref
        den = np.dot(ray, n)
        if abs(den) < 1e-10:
            return False
        f = np.dot(sc - ref, n) / den
        point = ref + f * ray
        if not (
            0 < f < 1
            and abs((point - sc) @ tan) < scr["width_mm"] / 2
            and scr["sill_mm"] + scr["horizontal_joint_gap_mm"] / 2
            < point[2]
            < scr["height_mm"]
        ):
            return False
    return True


def write_projection_survey(results, out):
    """Keep the individual beam and view estimates before any weekly weighting."""
    rows = []
    for case in results:
        if case["component"] != "field":
            continue
        projection = case.get("patient_projection") or {}
        for name, site in case["sites"].items():
            if site["group"] == "beam":
                continue
            rows.append(dict(
                layout=case["layout"], phantom=case["phantom"], beam=case["beam"],
                view=case["view"], arm_angle_deg=case["view_deg"],
                tube_offset_deg=case["tube_offset_deg"], projection_number=case["projection_number"],
                torso_present=case["torso_present"], door_open=case["door_open"],
                grid_present=case["grid_present"], spectrum_model=case["spectrum_model"],
                source_compression_view=projection.get("source_compression_view"),
                anatomy_model=projection.get("anatomy_model"),
                batch=case["batch"], histories=case["events"],
                position=name, group=site["group"],
                breast_reference_to_cell_m=site.get("breast_reference_to_cell_m"),
                x_mm=site["xyz_mm"][0], y_mm=site["xyz_mm"][1], z_mm=site["xyz_mm"][2],
                collision_air_kerma_uGy_mAs=site["collision_air_kerma_uGy_mAs"],
                secondary_field_uGy_mAs=site["secondary_field_uGy_mAs"],
                photon_steps=site["photon_steps"], sub_keV_steps=site["sub_keV_steps"],
            ))
    csv_file(out / "projection_survey.csv", rows)


def exposure_report(results, work, out, scan_mas):
    """Combine projections inside each independent batch before estimating the SE."""
    totals = {}
    for case in results:
        weight = work[case["beam"]][case["view"]]
        for name, site in case["sites"].items():
            if site["group"] == "beam":
                continue
            batches = totals.setdefault((case["layout"], name), {})
            estimate = batches.setdefault(case["batch"], np.zeros(3))
            estimate += [weight*site["collision_air_kerma_uGy_mAs"],
                         weight*site["secondary_field_uGy_mAs"], site["photon_steps"]]
    rows = []
    for (layout, position), batches in totals.items():
        values = np.array(list(batches.values()))
        mean, secondary, _ = values.mean(axis=0)
        se = float(values[:, 0].std(ddof=1)/np.sqrt(len(values))) if len(values)>1 else None
        rows.append(dict(layout=layout, position=position, independent_batches=len(values),
                         field_uGy_per_total_mAs=mean, secondary_uGy_per_total_mAs=secondary,
                         batch_SE_uGy_per_total_mAs=se, relative_SE=se/mean if mean and se is not None else None,
                         batches_with_photon_steps=int(np.count_nonzero(values[:, 2])),
                         photon_steps=int(values[:, 2].sum()),
                         field_uGy_per_sweep=mean*scan_mas if scan_mas is not None else None))
    csv_file(out / "exposure_survey.csv", rows)


def compare(cfg, work, layouts, results, out):
    p = cfg["planning"]
    batches = sorted({r["batch"] for r in results})
    index = {
        (r["layout"], r["beam"], r["view"], r["component"], r["batch"]): r
        for r in results
    }
    sources = ["field"] if p["leakage_model"] is None else ["field", "leakage"]
    scenarios = [(b, v, c) for b, views in work.items() for v in views for c in sources]
    family = max(
        1,
        sum(
            sum(
                s["group"] in p["goal_uGy_week"]
                for s in survey_points(cfg, cfg["layouts"][l]).values()
            )
            for l in layouts
            if not cfg["layouts"][l].get("reference_only")
        ),
    )
    rows, decisions = [], []
    for l in layouts:
        v = cfg["layouts"][l]
        base = index[(l, *scenarios[0], batches[0])]
        for name, site in base["sites"].items():
            group = site["group"]
            if group == "beam":
                continue
            occ = p["occupancy"].get(group, 1.0)
            goal = p["goal_uGy_week"].get(group)
            vals = np.zeros((len(batches), 3))
            positives = {key: 0 for key in scenarios}
            low = 0
            hits = 0
            for j, rep in enumerate(batches):
                for b, view, src in scenarios:
                    item = index[(l, b, view, src, rep)]["sites"][name]
                    w = (
                        work[b][view]
                        * occ
                        * (p["source_scale"][b] if src == "field" else 1.0)
                    )
                    if src == "field":
                        vals[j, 0] += w * item["collision_air_kerma_uGy_mAs"]
                        vals[j, 1] += w * item["secondary_field_uGy_mAs"]
                    else:
                        vals[j, 2] += w * item["collision_air_kerma_uGy_mAs"]
                    positives[b, view, src] += int(item["photon_steps"] > 0)
                    low += item["sub_keV_steps"]
                    hits += item["photon_steps"]
            combined = vals[:, 0] + vals[:, 2]
            mean = float(np.mean(combined))
            se = (
                float(np.std(combined, ddof=1) / np.sqrt(len(batches)))
                if len(batches) > 1
                else None
            )
            rse = se / mean if mean > 0 and se is not None else None
            precise = (
                len(batches) >= p["min_independent_batches"]
                and min(positives.values()) >= max(3, math.ceil(len(batches) / 2))
                and rse is not None
                and 0 < rse <= p["max_relative_se"]
                and low == 0
            )
            upper = (
                float(mean + t.ppf(1 - 0.05 / family, len(batches) - 1) * se)
                if precise and goal
                else None
            )
            status = (
                "map"
                if goal is None
                else (
                    "unresolved"
                    if upper is None
                    else ("within_goal" if upper <= goal else "upper_exceeds_goal")
                )
            )
            rows.append(
                dict(
                    layout=l,
                    position=name,
                    group=group,
                    x_mm=site["xyz_mm"][0],
                    y_mm=site["xyz_mm"][1],
                    z_mm=site["xyz_mm"][2],
                    occupancy=occ,
                    field_uGy_week=float(vals[:, 0].mean()),
                    secondary_field_uGy_week=float(vals[:, 1].mean()),
                    housing_leakage_uGy_week=float(vals[:, 2].mean())
                    if "leakage" in sources
                    else None,
                    modelled_sum_uGy_week=mean,
                    batch_SE=se,
                    relative_SE=rse,
                    upper95_simultaneous_uGy_week=upper,
                    goal_uGy_week=goal,
                    photon_steps=hits,
                    sub_keV_steps=low,
                    batches=len(batches),
                    status=status,
                )
            )
        scored = [
            row
            for row in rows
            if row["layout"] == l and row["goal_uGy_week"] is not None
        ]
        unresolved = sum(row["status"] == "unresolved" for row in scored)
        exceed = sum(row["status"] == "upper_exceeds_goal" for row in scored)
        margin = base["nominal_primary_margin_mm"]
        fit = bool(fits(cfg, v))
        sight = sightlines(cfg, v)
        missing = []
        if not p["basis_documented"]:
            missing.append("design_inputs")
        if not p["coverage_reviewed"]:
            missing.append("projections_and_access")
        if "leakage" not in sources:
            missing.append("housing_leakage")
        if not fit or not sight or margin <= 0:
            missing.append("geometry")
        if unresolved:
            missing.append("statistical_precision")
        if exceed:
            missing.append("upper_exceeds_goal")
        status = (
            "reference_only"
            if v.get("reference_only")
            else ("conditional_candidate" if not missing else "unresolved")
        )
        decisions.append(
            dict(
                layout=l,
                distance_m=v["distance_m"],
                reference_to_screen_plane_m=float(
                    abs(np.dot(booth(cfg, v)[1] - reference(cfg, v), booth(cfg, v)[2]))
                    / 1000
                ),
                booth_fits=fit,
                nominal_sightlines=sight,
                primary_backstop_margin_mm=margin,
                positions_with_goals=len(scored),
                unresolved_positions=unresolved,
                upper_exceeds_goal=exceed,
                status=status,
                outstanding=" | ".join(missing),
            )
        )
    csv_file(out / "weekly_survey.csv", rows)
    csv_file(out / "layout_comparison.csv", decisions)
    candidates = [r for r in decisions if r["status"] == "conditional_candidate"]
    best = min((r["distance_m"] for r in candidates), default=None)
    write(
        out / "decision.json",
        dict(
            shortest_qualifying_tested_distance_m=best,
            layouts_at_that_distance=[
                r["layout"] for r in candidates if r["distance_m"] == best
            ],
            quantities="Occupancy-weighted collision air kerma in dry air. Missing leakage is unassessed.",
            comparisons=decisions,
        ),
    )
    return rows


def scatter_map(cfg, rows, layout, out):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    items = [
        r
        for r in rows
        if r["layout"] == layout
        and (r["position"].startswith("ring_") or r["position"].endswith("_1700"))
    ]
    vals = [
        r["secondary_field_uGy_week"]
        for r in items
        if r["secondary_field_uGy_week"] > 0
    ]
    fig, ax = plt.subplots(figsize=(7, 5), layout="constrained")
    w, d, _ = np.asarray(cfg["room"]["inside_mm"]) / 1000
    ax.plot(
        np.array([-1, 1, 1, -1, -1]) * w / 2,
        np.array([-1, -1, 1, 1, -1]) * d / 2,
        color="#34555c",
    )
    if cfg["layouts"][layout]["screen_pb_mm"]:
        f = footprint(cfg, cfg["layouts"][layout]) / 1000
        ax.plot(f[:, 0], f[:, 1], color="#29444c")
    if vals:
        norm = LogNorm(min(vals), max(max(vals), min(vals) * 1.01))
        fig.colorbar(
            plt.cm.ScalarMappable(norm=norm, cmap="viridis"),
            ax=ax,
            label="Secondary-field collision air kerma (µGy/week)",
        )
    for r in items:
        value = r["secondary_field_uGy_week"]
        pt = (r["x_mm"] / 1000, r["y_mm"] / 1000)
        if value:
            ax.scatter(
                *pt,
                c=[value],
                cmap="viridis",
                norm=norm,
                marker="o" if r["position"].startswith("ring_") else "s",
                s=45,
            )
        else:
            ax.scatter(*pt, marker="x", color="gray")
    ax.set(
        aspect="equal",
        xlabel="Room x (m)",
        ylabel="Room y (m)",
        title=f"{layout}   |   Breast-height ring and staff eye positions",
    )
    ax.text(
        0.01,
        0.01,
        "Cross = no sampled contribution",
        transform=ax.transAxes,
        fontsize=8,
    )
    fig.savefig(out, dpi=145)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("_worker", help=argparse.SUPPRESS).add_argument("spec")
    for cmd in ("preview", "run"):
        q = sub.add_parser(cmd)
        q.add_argument("--config", type=Path, required=(cmd == "run"), default=HERE / "study.json")
        q.add_argument(
            "--layouts", nargs="+"
        )
        q.add_argument("--out", type=Path, default=HERE / "room_results")
        q.add_argument("--programme", help="viewer, dbt_sweep, or a named weekly workload")
        q.add_argument("--phantom", choices=PHANTOMS)
        q.add_argument("--beam", choices=("LE", "HE", "DBT"), help="override the single viewer exposure")
        q.add_argument("--view", choices=("CC", "MLO_R", "MLO_L"), help="base view; uses its configured arm angle")
        q.add_argument("--door", choices=("open", "closed"))
        q.add_argument("--no-shield", action="store_true")
        q.add_argument("--shield-width-mm", type=float)
        if cmd == "run":
            q.add_argument(
                "--scanner",
                type=Path,
                default=HERE / "Revelation_Scanner_Inputs" / "gate_mammography",
                help="scanner directory (defaults to the bundled MAMMOMAT Revelation model)",
            )
            q.add_argument("--saved-export", action="store_true",
                           help="run the exact viewer selection and group outputs by export under --out")
            q.add_argument("--dbt-spectrum", type=Path, help="DBT spectrum and fluence in native NPZ format")
            q.add_argument("--scan-mas", type=float, help="total mAs of the 25 DBT pulses combined")
            q.add_argument(
                "--histories",
                type=int,
                default=1000000,
                help="primary events per beam/view/component/batch",
            )
            q.add_argument("--batches", type=int, default=2)
            q.add_argument(
                "--jobs",
                type=int,
                default=1,
                help="independent single-thread processes",
            )
            q.add_argument(
                "--check-overlaps",
                action="store_true",
                help="expensive with all L2 spheres",
            )
            q.add_argument("--keep-tracks", action="store_true")
    args = parser.parse_args()
    if args.command == "_worker":
        worker(args.spec)
        return
    try:
        source = args.config.read_bytes()
        cfg = json.loads(source)
    except (OSError, ValueError) as exc:
        parser.error(f"Cannot read settings file {args.config}: {exc}")
    if not isinstance(cfg, dict):
        parser.error("The settings file must contain a JSON object.")
    source_hash = hashlib.sha256(source).hexdigest()
    viewer = cfg.get("revelation_viewer", {})
    selected = viewer.get("selected_layout") if isinstance(viewer, dict) else None
    out = args.out.resolve()
    saved_export = getattr(args, "saved_export", False)
    if saved_export:
        if (not isinstance(selected, str) or selected not in cfg.get("layouts", {})
                or viewer.get("phantom") not in PHANTOMS or viewer.get("projection") not in cfg.get("views_deg", {})):
            parser.error("Choose a geometry JSON saved from the Revelation viewer with a selected layout, phantom and projection.")
        if (args.layouts or args.phantom or args.programme or args.beam or args.view or args.door
                or args.no_shield or args.shield_width_mm is not None):
            parser.error("A saved-export run takes its layout and exposure selections directly from the file.")
        out = out / args.config.stem / source_hash[:12] / f"n{args.histories}_b{args.batches}"
    out.mkdir(parents=True, exist_ok=True)
    defaults = [selected] if selected else ["siemens_2m", "siemens_2p5m", "siemens_3m"]
    layouts = list(dict.fromkeys(args.layouts or defaults))
    args.programme = args.programme or ("viewer" if selected else "illustrative_CEM_CC")
    args.phantom = args.phantom or cfg.get("revelation_viewer", {}).get("phantom", "L2")
    if args.phantom not in PHANTOMS:
        parser.error("Unknown saved phantom.")
    if args.programme not in ("viewer", "dbt_sweep") and args.programme not in cfg["workloads_mAs_per_week"]:
        parser.error("Unknown acquisition programme.")
    if args.beam and args.programme != "viewer":
        parser.error("--beam applies to a single viewer exposure.")
    if args.view and args.programme not in ("viewer", "dbt_sweep"):
        parser.error("--view applies to a single exposure or DBT sweep.")
    viewer = cfg.setdefault("revelation_viewer", {})
    options = viewer.setdefault("study_preview", {})
    if args.beam:
        options["acquisition"] = args.beam
    if args.view:
        viewer["projection"] = args.view
        options["base_angle_override_deg"] = None
    if args.door:
        options["door_open"] = args.door == "open"
    if args.no_shield:
        options["no_shield"] = True
    if args.shield_width_mm is not None:
        if not math.isfinite(args.shield_width_mm) or args.shield_width_mm <= 0:
            parser.error("Shield width must be positive and finite.")
        for name in layouts:
            cfg["layouts"][name].setdefault("screen", {})["width_mm"] = args.shield_width_mm
    work = configure(cfg, args.programme, layouts, args.phantom)
    if saved_export:
        profile = cfg["layouts"][selected]
        screen = screen_settings(cfg, profile)
        shield = f'{screen["width_mm"]:g} x {screen["height_mm"]:g} mm, {profile["screen_pb_mm"]:g} mm Pb' if profile["screen_pb_mm"] else "absent"
        beam = cfg["acquisition"]["beam"]
        shot = next(iter(cfg["acquisition"]["shots"].values()))
        projection = f', DBT projection {shot["projection_number"]}/25' if beam == "DBT" else ""
        print(f"Settings: {args.config.resolve()}\nLayout: {selected} | room {cfg['room']['inside_mm']} mm | "
              f"breast to operator {profile['distance_m']:.3f} m | shield {shield}\n"
              f"Exposure: {args.phantom}, {beam}, {shot['base_view']} at {shot['base_angle_deg']:g} degrees{projection} | "
              f"door {'open' if cfg['room']['door_open'] else 'closed'}\nResults: {out}", flush=True)
    for name in layouts:
        profile = cfg["layouts"][name]
        if profile["screen_pb_mm"] == 0 or cfg["room"]["door_open"] or args.phantom == "open":
            profile["reference_only"] = True
    research = args.programme in ("viewer", "dbt_sweep")
    if args.command == "run":
        if args.scan_mas is not None and (args.programme != "dbt_sweep" or not math.isfinite(args.scan_mas) or args.scan_mas <= 0):
            parser.error("--scan-mas requires a full DBT sweep and a positive finite tube load.")
        if args.dbt_spectrum:
            if "DBT" not in work:
                parser.error("--dbt-spectrum requires a DBT acquisition.")
            read_spectrum(args.dbt_spectrum)
            cfg["acquisition"].update(spectrum="user-supplied DBT spectrum", spectrum_file=str(args.dbt_spectrum.resolve()))
    audit = {}
    for name in layouts:
        layout = cfg["layouts"][name]
        if "survey_geometry" in layout:
            if not fits(cfg, layout):
                parser.error(f"The shield or staff cells extend beyond the room in {name}. Review the placement.")
            excluded = {}
            sites = survey_points(cfg, layout, excluded)
            audit[name] = {"active_cells": sites, "unavailable_cells": excluded}
    if audit:
        write(out / "survey_geometry.json", audit)
    if args.command == "preview":
        preview(cfg, layouts, out / "floorplan.png")
        return
    if min(args.batches, args.jobs, args.histories) < 1:
        parser.error("Use positive integer histories, batches, and jobs.")
    validate(cfg, work, layouts)
    if args.phantom in PHANTOM_MODES:
        for view in dict.fromkeys(view for views in work.values() for view in views):
            try:
                shot = cfg.get("acquisition", {}).get("shots", {}).get(view, {})
                projection_for_scatter(shot.get("base_view", view), cfg["views_deg"][view], custom=bool(shot))
            except ValueError as error:
                parser.error(str(error))
    if args.phantom == "open" and not research and cfg["planning"]["basis_documented"]:
        parser.error(
            "The empty-field run is a scatter sensitivity, separate from workload qualification."
        )
    scanner = args.scanner.resolve()
    archived = inputs(scanner, args.phantom, cfg)
    fingerprint = digest(archived)
    root = out / "cases" / fingerprint[:12]
    root.mkdir(parents=True, exist_ok=True)
    # Retain the exact source bytes as well as the effective transport configuration.
    source_copy = root / f"source_{source_hash}.json"
    if not source_copy.is_file() or source_copy.read_bytes() != source:
        temp = source_copy.with_suffix(f".{uuid.uuid4().hex}.tmp")
        temp.write_bytes(source)
        os.replace(temp, source_copy)
    write(root / "inputs.json", archived)
    config_path = root / "study.json"
    write(config_path, cfg)
    sources = (
        ["field"] if research or cfg["planning"]["leakage_model"] is None else ["field", "leakage"]
    )
    tasks = []
    for l in layouts:
        for b, views in work.items():
            for view in views:
                for component in sources:
                    for r in range(1, args.batches + 1):
                        key = [l, b, view, component, args.phantom, r]
                        seed = int(digest(key)[:8], 16) % 2000000000 + 1
                        case_dir = root / (
                            f"{l}_{b}_{view}_{component}_{args.phantom}_n{args.histories}_r{r:03d}"
                        )
                        tasks.append(
                            dict(
                                scanner=str(scanner),
                                config_path=str(config_path),
                                config_digest=digest(cfg),
                                input_sha256=fingerprint,
                                out=str(case_dir),
                                layout=l,
                                beam=b,
                                view=view,
                                component=component,
                                phantom=args.phantom,
                                n=args.histories,
                                seed=seed,
                                batch=r,
                                check_overlaps=args.check_overlaps,
                                keep_tracks=args.keep_tracks,
                            )
                        )
    print(f"{args.phantom}: {len(tasks)} GATE cases, {args.histories:,} histories per case", flush=True)
    results = []
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        for done, result in enumerate(pool.map(execute_one, tasks), 1):
            results.append(result)
            print(f"{done}/{len(tasks)}  {result['layout']}  {result['beam']}  {result['view']}  batch {result['batch']}", flush=True)
    if research:
        exposure_report(results, work, out, args.scan_mas)
        write(out / "decision.json", dict(status="research_only", programme=args.programme,
              quantity="Field collision air kerma per mAs. No weekly workload or housing leakage applied."))
    else:
        rows = compare(cfg, work, layouts, results, out)
        for l in layouts:
            scatter_map(cfg, rows, l, out / f"scatter_{l}.png")
    write_projection_survey(results, out)
    preview(cfg, layouts, out / "floorplan.png")
    write(
        out / "run_manifest.json",
        dict(
            model=str(scanner),
            source_file=str(args.config.resolve()), source_sha256=source_hash,
            source_archive=str(source_copy.relative_to(out)),
            input_sha256=fingerprint,
            programme=args.programme,
            workload_mAs_week=None if research else work,
            exposure_weights=work if research else None, scan_mAs=args.scan_mas,
            acquisition=cfg.get("acquisition"), door_open=cfg["room"]["door_open"],
            saved_case=cfg.get("revelation_viewer", {}).get("study_preview", {}).get("case_id"),
            placements={l: dict(breast_reference_room_mm=reference(cfg, cfg["layouts"][l]).tolist(),
                                operator_floor_room_mm=booth(cfg, cfg["layouts"][l])[0].tolist(),
                                breast_to_operator_horizontal_m=cfg["layouts"][l]["distance_m"],
                                screen_pb_mm=cfg["layouts"][l]["screen_pb_mm"],
                                screen=screen_settings(cfg, cfg["layouts"][l])) for l in layouts},
            phantom=args.phantom,
            patient_anatomy=results[0].get("patient_anatomy"),
            patient_projections={
                r["view"]: r["patient_projection"]
                for r in results if r.get("patient_projection")
            },
            jobs=len(results),
            total_events=sum(r["events"] for r in results),
            config_snapshot=str(config_path),
        ),
    )
    print((out / "decision.json").read_text())


if __name__ == "__main__":
    main()

"""Resolve viewer settings and the fixed-detector DBT research acquisition."""
import math
from pathlib import Path

import numpy as np


def configure(cfg, programme, layouts, phantom):
    viewer = cfg.get("revelation_viewer", {})
    opts = viewer.get("study_preview", {})
    for name in layouts:
        if opts.get("no_shield"):
            cfg["layouts"][name]["screen_pb_mm"] = 0
    cfg["room"]["door_open"] = bool(opts.get("door_open", False))
    cfg["room"]["torso"] = phantom != "open"
    if programme not in ("viewer", "dbt_sweep"):
        return cfg["workloads_mAs_per_week"][programme]
    base = viewer.get("projection", "CC")
    angle = opts.get("base_angle_override_deg")
    angle = cfg["views_deg"][base] if angle is None else angle
    if base not in ("CC", "MLO_R", "MLO_L") or not math.isfinite(angle) or abs(angle) > 90:
        raise ValueError("Unsupported base view or arm angle.")
    beam = "DBT" if programme == "dbt_sweep" else opts.get("acquisition", "LE")
    if beam not in ("LE", "HE", "DBT"):
        raise ValueError("Unsupported viewer acquisition.")
    indices = range(25) if programme == "dbt_sweep" else [opts.get("dbt_index", 12)]
    shots = {}
    for i in indices:
        if type(i) is not int or not 0 <= i <= 24:
            raise ValueError("DBT index must be an integer from 0 to 24.")
        key = f"{base}_DBT_{i+1:02d}" if beam == "DBT" else base
        cfg["views_deg"][key] = angle
        shots[key] = dict(base_view=base, base_angle_deg=angle,
                          tube_offset_deg=25 - i * 50 / 24 if beam == "DBT" else 0,
                          projection_number=i+1 if beam == "DBT" else None)
    cfg["acquisition"] = dict(programme=programme, beam=beam, shots=shots,
                              spectrum="bundled LE proxy" if beam == "DBT" else beam)
    return {beam: {key: 1/len(shots) for key in shots}}


def read_spectrum(path):
    with np.load(Path(path), allow_pickle=False) as data:
        edges = np.asarray(data["edges_keV"], float)
        weights = np.asarray(data["weights"], float)
        fluence = float(data["fluence_per_mAs_1m"].item())
    if (edges.ndim != 1 or weights.ndim != 1 or len(edges) != len(weights)+1
            or not len(weights) or not np.all(np.isfinite(edges))
            or edges[0] < 1 or edges[-1] > 100 or not np.all(np.diff(edges) > 0)
            or not np.all(np.isfinite(weights)) or np.min(weights) < 0 or weights.sum() <= 0
            or not math.isfinite(fluence) or fluence <= 0):
        raise ValueError("Invalid DBT spectrum bins, photon weights or fluence.")
    return edges, weights, fluence


def dbt_source(sim, geometry, info, sid, pivot_z, angle, spectrum):
    """Move focus on the viewer arc. Reproject the fixed detector rectangle to the aperture."""
    source = sim.source_manager.get_source("tube")
    aperture = sim.volume_manager.get_volume("aperture")
    a = math.radians(angle)
    radius = sid - pivot_z
    focus = np.array([0., -radius*math.sin(a), pivot_z + radius*math.cos(a)])
    target = np.array([np.mean(geometry.FIELD_X), np.mean(geometry.FIELD_Y), 0.])
    k = (sid - geometry.APERTURE_Z)/sid
    aperture.translation = (focus + k*(target-focus)).tolist()
    source.position.translation = focus.tolist()
    # Keep the disc normal to the nominal ray, and sample directions in the native frame.
    from room_geometry import rot
    source.position.rotation = rot("x", angle)
    far = max(math.hypot(x, y-focus[1]) for x in geometry.FIELD_X for y in geometry.FIELD_Y)
    theta = math.atan2(far+1, focus[2])+math.radians(0.5)
    # Full azimuth covers the finite spot on either side of the chest-wall boundary.
    source.direction.theta = [0., theta]
    source.direction.phi = [0., 2*math.pi]
    info["omega_cone_sr"] = 2*math.pi*(1-math.cos(theta))
    if spectrum:
        edges, weights, fluence = read_spectrum(spectrum)
        source.energy.spectrum_energy_bin_edges = (edges*0.001).tolist()
        source.energy.spectrum_weights = weights.tolist()
        info["spectrum"]["fluence_per_mAs_1m"] = fluence
    return focus

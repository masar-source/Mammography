"""OpenGATE 10 transport. All length arguments are mm, positions use a right-handed frame."""

from pathlib import Path
import math
import numpy as np
from scipy.spatial.transform import Rotation


def rot(axis, angle):
    return Rotation.from_euler(axis, angle, degrees=True).as_matrix()


def simulation(out, seed, overlaps=False):
    import opengate as gate

    s = gate.Simulation()
    s.output_dir = str(Path(out).resolve())
    s.number_of_threads = 1
    s.random_seed = int(seed)
    s.check_volumes_overlap = bool(overlaps)
    s.physics_manager.physics_list_name = "G4EmLivermorePhysics"
    s.physics_manager.global_production_cuts.all = 0.1 * gate.g4_units.mm
    em = s.physics_manager.em_parameters
    em.fluo = True
    em.auger = True
    em.deexcitation_ignore_cut = True
    s.g4_commands_before_init.append("/cuts/setLowEdge 100 eV")
    return s


def box(sim, name, size, centre, material, rotation=None, mother="world"):
    v = sim.add_volume("Box", name)
    v.mother = mother
    v.size = list(map(float, size))
    v.translation = list(map(float, centre))
    v.material = material
    if rotation is not None:
        v.rotation = rotation
    return v


def materials(sim, cfg):
    import opengate as gate

    db = sim.volume_manager.material_database
    density = gate.g4_units.g / gate.g4_units.cm3
    # CaSO4.2H2O at board bulk density; chemistry and density are separate inputs.
    db.add_material_nb_atoms(
        "RoomBoard",
        ["Ca", "S", "O", "H"],
        [1, 1, 6, 4],
        cfg["room"]["board_density_g_cm3"] * density,
    )
    db.add_material_nb_atoms(
        "DoorWood",
        ["C", "H", "O"],
        [6, 10, 5],
        cfg["room"]["timber_density_g_cm3"] * density,
    )
    # ICRU-44 breast composition from NIST X-Ray Mass Attenuation Coefficients, Table 2.
    db.add_material_weights(
        "BreastICRU44",
        ["H", "C", "N", "O", "Na", "P", "S", "Cl"],
        [0.106, 0.332, 0.030, 0.527, 0.001, 0.001, 0.002, 0.001],
        1.02 * density,
    )


def scanner_installation(layout):
    """Rigid floor placement in GATE room coordinates. Old layouts use the identity."""
    if "scanner_pose" not in layout:
        return np.eye(3), np.zeros(3)
    pose = layout["scanner_pose"]
    try:
        pivot = np.asarray(pose["pivot_xy_mm"], dtype=float)
        move = np.asarray(pose["translation_xy_mm"], dtype=float)
        angle = float(pose["yaw_deg"])
        if (pivot.shape != (2,) or move.shape != (2,) or not np.isfinite(pivot).all()
                or not np.isfinite(move).all() or not math.isfinite(angle)):
            raise ValueError()
    except (KeyError, TypeError, ValueError, OverflowError) as exc:
        raise ValueError("scanner_pose requires finite pivot_xy_mm, translation_xy_mm and yaw_deg.") from exc
    pivot = np.r_[pivot, 0.0]
    yaw = rot("z", angle)
    return yaw, pivot + np.r_[move, 0.0] - yaw @ pivot


def reference(cfg, layout=None):
    neutral = np.asarray(cfg["room"]["detector_chestwall_mm"], float) + np.asarray(
        cfg["room"]["breast_reference_local_mm"], float
    )
    if layout is None:
        return neutral
    yaw, displacement = scanner_installation(layout)
    return yaw @ neutral + displacement


def screen_settings(cfg, layout):
    screen = {**cfg["screen"], **layout.get("screen", {})}
    pb = layout["screen_pb_mm"]
    if not math.isfinite(pb) or pb < 0:
        raise ValueError("Screen lead equivalence must be finite and non-negative.")
    if pb == 0:
        return screen
    width, height, sill = screen["width_mm"], screen["height_mm"], screen["sill_mm"]
    floor, gap = screen["floor_gap_mm"], screen["horizontal_joint_gap_mm"]
    side, factor = layout["return_mm"], screen["window_thickness_factor"]
    if not all(math.isfinite(v) for v in (width, height, sill, floor, gap, side, factor)):
        raise ValueError("Screen dimensions and shielding values must be finite.")
    if (width <= 0 or gap < 0 or floor < 0 or side < 0 or factor <= 0
            or sill - gap / 2 <= floor or sill + gap / 2 >= height):
        raise ValueError("Check screen width, height, window sill, gaps, side return and shielding values.")
    return screen


def equipment_headings(layout):
    """Floor yaw in degrees, CCW from room +Y along the panel and worktop edges."""
    legacy = float(layout["azimuth_deg"]) + float(layout["screen_yaw_deg"])
    screen = float(layout.get("screen_heading_deg", legacy))
    console = float(layout.get("console_heading_deg", screen))
    if not all(math.isfinite(x) for x in (legacy, screen, console)):
        raise ValueError("The operator bearing, screen heading and console heading must be finite.")
    return screen, console


def booth(cfg, layout):
    a = math.radians(layout["azimuth_deg"])
    b = math.radians(equipment_headings(layout)[0])
    n = np.array([math.cos(b), math.sin(b), 0.0])
    tangent = np.array([-math.sin(b), math.cos(b), 0.0])
    operator = reference(cfg, layout) + 1000 * layout["distance_m"] * np.array(
        [math.cos(a), math.sin(a), 0]
    )
    operator[2] = 0
    scr = screen_settings(cfg, layout)
    screen = operator - scr["offset_mm"] * n + scr.get("lateral_offset_mm", 0) * tangent
    return operator, screen, n, tangent


def screen_intersects_cell(cfg, layout, point, half):
    """Check a room-aligned air cube against the oriented front panel and return."""
    pb = layout["screen_pb_mm"]
    if pb <= 0:
        return False
    scr = screen_settings(cfg, layout)
    _, centre, normal, tangent = booth(cfg, layout)
    if point[2] + half < scr["floor_gap_mm"] or point[2] - half > scr["height_mm"]:
        return False
    # A solid envelope conservatively includes the front panel's horizontal joint.
    thick = pb * max(1, scr["window_thickness_factor"])
    panels = [(centre, thick, scr["width_mm"])]
    length = layout["return_mm"]
    if length:
        corner = centre + tangent * (scr["width_mm"] + pb) / 2 + normal * length / 2
        panels.append((corner, length, pb))
    for origin, depth, width in panels:
        for axis in ([1, 0, 0], [0, 1, 0], normal, tangent):
            extent = (abs(np.dot(normal, axis)) * depth + abs(np.dot(tangent, axis)) * width) / 2
            if abs(np.dot(point - origin, axis)) > extent + half * np.abs(axis).sum():
                break
        else:
            return True
    return False


def survey_points(cfg, layout, exclusions=None):
    pts = {}

    def add(name, g, p, size):
        pts[name] = dict(
            group=g,
            xyz_mm=[float(v) for v in p],
            size_mm=list(np.broadcast_to(size, (3,)).astype(float)),
        )

    ref = reference(cfg, layout)
    yaw, _ = scanner_installation(layout)
    settings = layout.get("survey_geometry")
    radii = [1000] if settings is None else [
        settings["common_ring_radius_mm"], settings["supplementary_ring_radius_mm"]
    ]
    if any(not np.isfinite(r) or r <= 0 for r in radii) or len(set(radii)) != len(radii):
        raise ValueError("Supply distinct positive scatter-ring radii.")
    width, depth, height = cfg["room"]["inside_mm"]
    lower = np.array([-width / 2, -depth / 2, 0])
    upper = np.array([width / 2, depth / 2, height])
    for index, radius in enumerate(radii):
        for degrees in range(0, 360, 45):
            angle = math.radians(degrees)
            point = ref + yaw @ [-radius * math.cos(angle), radius * math.sin(angle), 0]
            name = f"ring_{degrees:03d}" if settings is None else f"ring_{radius:g}mm_{degrees:03d}"
            if settings is not None:
                reason = None
                if np.any(point - 50 <= lower) or np.any(point + 50 >= upper):
                    reason = "room boundary"
                elif screen_intersects_cell(cfg, layout, point, 50):
                    reason = "shield"
                if reason:
                    if index == 0:
                        raise ValueError(f"Required {name} meets {reason}.")
                    if exclusions is not None:
                        exclusions[name] = dict(reason=reason, xyz_mm=point.tolist())
                    continue
            add(name, "map", point, 100)
    op, sc, n, t = booth(cfg, layout)
    console_angle = math.radians(equipment_headings(layout)[1])
    console_tangent = np.array([-math.sin(console_angle), math.cos(console_angle), 0.0])
    for role, side in (
        ("operator", 0),
        ("trainee", 440),
        ("left_reach", -220),
        ("outer_reach", 715 if settings is None else settings["outer_reach_mm"]),
    ):
        for z in (500, 1200, 1700):
            add(f"{role}_{z}", "booth", op + side * console_tangent + [0, 0, z], 65)
    for tag, point in (
        ("source_side", sc - 160 * n),
        ("window_protected", sc + 140 * n),
        (
            "open_end",
            sc + 140 * n - 0.57 * screen_settings(cfg, layout)["width_mm"] * t,
        ),
    ):
        add(tag, "map", point + [0, 0, 1500], 60)
    width, depth, height = cfg["room"]["inside_mm"]
    # Fixed physical evaluation coordinates across wall thickness variants.
    hx, hy = width / 2, depth / 2
    add("west_office", "office", [-hx - 350, 0, 1200], 120)
    add("east_office", "office", [hx + 350, 0, 1200], 120)
    if layout.get("open_service_duct"):
        add(
            "east_duct", "office", [hx + 175, *cfg["room"]["service_opening_yz_mm"]], 65
        )
    add("north_office", "office", [0, hy + 350, 1200], 120)
    dx = cfg["room"]["door_centre_x_mm"]
    # Keep the survey station on the solid section between the west corner and the door.
    # In the preset it is 2200 mm along that 4150 mm section, at room x = -1000 mm.
    left_wall_length = dx - cfg["room"]["door_width_mm"] / 2 + hx
    south_wall_x = -hx + left_wall_length * (2200.0 / 4150.0)
    add("south_wall", "corridor", [south_wall_x, -hy - 350, 1200], 120)
    add("door_leaf", "corridor", [dx, -hy - 350, 1200], 100)
    add(
        "door_jamb",
        "corridor",
        [dx + cfg["room"]["door_width_mm"] / 2, -hy - 160, 1200],
        70,
    )
    add("below", "floor", [ref[0], ref[1], -cfg["room"]["concrete_slab_mm"] - 300], 160)
    add(
        "above",
        "ceiling",
        [ref[0], ref[1], height + cfg["room"]["concrete_slab_mm"] + 300],
        160,
    )
    if settings is not None:
        checked = []
        for name, site in pts.items():
            if site["group"] not in ("booth", "map"):
                continue
            centre = np.asarray(site["xyz_mm"])
            half = np.asarray(site["size_mm"]) / 2
            if not np.all(np.isfinite(centre)) or np.any(centre - half <= lower) or np.any(centre + half >= upper):
                raise ValueError(f"Required survey cell {name} crosses the inner room boundary.")
            if screen_intersects_cell(cfg, layout, centre, half[0]):
                raise ValueError(f"Survey cell {name} intersects the shield.")
            for other, other_centre, other_half in checked:
                if np.all(np.abs(centre - other_centre) < half + other_half):
                    raise ValueError(f"Survey cells {name} and {other} overlap.")
            checked.append((name, centre, half))
    return pts


def validate_room_geometry(cfg, layout):
    equipment_headings(layout)
    room = cfg["room"]
    inside = np.asarray(room["inside_mm"], dtype=float)
    if inside.shape != (3,) or not np.all(np.isfinite(inside)) or np.any(inside <= 0):
        raise ValueError("room.inside_mm must contain positive, finite width, depth and height in mm.")
    w, d, h = inside
    door_x, door_w, door_h = room["door_centre_x_mm"], room["door_width_mm"], room["door_height_mm"]
    if not (door_w > 4 and abs(door_x) + door_w / 2 < w / 2 and 4 < door_h < h):
        raise ValueError("The south wall needs room for the door, both jambs and the header.")
    if layout.get("open_service_duct"):
        dy, dz = room["service_opening_size_mm"]
        y, z = room["service_opening_yz_mm"]
        if not (dy > 0 and dz > 0 and abs(y) + dy / 2 < d / 2 and dz / 2 < z < h - dz / 2):
            raise ValueError("The service opening must fit inside the east wall.")


def walls(sim, cfg, layout):
    validate_room_geometry(cfg, layout)
    rm = cfg["room"]
    w, d, h = rm["inside_mm"]
    dx, dw, dh = rm["door_centre_x_mm"], rm["door_width_mm"], rm["door_height_mm"]
    outer = 0.0
    for mat, thick, tag in (
        ("G4_Pb", layout["wall_pb_mm"], "lining"),
        ("RoomBoard", layout["wall_board_mm"], "board"),
    ):
        if thick <= 0:
            continue
        x = w / 2 + outer + thick / 2
        y = d / 2 + outer + thick / 2
        for side in (-1, 1):
            sidewall = box(
                sim,
                f"{tag}_eastwest_{side}",
                [thick, d + 2 * (outer + thick), h],
                [side * x, 0, h / 2],
                mat,
            )
            if side == 1 and layout.get("open_service_duct"):
                dy, dz = rm["service_opening_size_mm"]
                py, pz = rm["service_opening_yz_mm"]
                box(
                    sim,
                    "duct_" + tag,
                    [thick, dy, dz],
                    [0, py, pz - h / 2],
                    "G4_AIR",
                    mother=sidewall.name,
                )
        box(sim, f"{tag}_north", [w + 2 * outer, thick, h], [0, y, h / 2], mat)
        for name, lo, hi in (
            ("left", -w / 2 - outer, dx - dw / 2),
            ("right", dx + dw / 2, w / 2 + outer),
        ):
            box(
                sim,
                f"{tag}_south_{name}",
                [hi - lo, thick, h],
                [(hi + lo) / 2, -y, h / 2],
                mat,
            )
        box(sim, f"{tag}_header", [dw, thick, h - dh], [dx, -y, (h + dh) / 2], mat)
        outer += thick
    for name, position in (
        ("floor", -rm["concrete_slab_mm"] / 2),
        ("ceiling", h + rm["concrete_slab_mm"] / 2),
    ):
        box(
            sim,
            name,
            [w + 2 * outer, d + 2 * outer, rm["concrete_slab_mm"]],
            [0, 0, position],
            "G4_CONCRETE",
        )
    # Closed door. Explicit nominal 2 mm clearance at jambs and sill/head.
    for name, material, start, end in (
        ("leaf", "DoorWood", 0, rm["door_wood_mm"]),
        (
            "door_lead",
            "G4_Pb",
            rm["door_wood_mm"],
            rm["door_wood_mm"] + rm["door_pb_mm"],
        ),
    ):
        if end <= start:
            continue
        if rm.get("door_open"):
            # Same 90-degree outward pose as the viewer. The jamb occupies the
            # intersecting hinge strip, so trim just that strip from the leaf.
            hinge_x = dx - (dw - 4)/2
            left, right = hinge_x - end, hinge_x - start
            length = dw - 4
            near = min(outer, length)
            for part, x0, y0, y1 in (
                ("hinge", max(left, dx-dw/2), 0, near),
                ("open", left, near, length),
            ):
                if right > x0 and y1 > y0:
                    box(sim, f"{name}_{part}", [right-x0, y1-y0, dh-4],
                        [(right+x0)/2, -d/2-(y0+y1)/2, dh/2], material)
        else:
            box(
                sim,
                name,
                [dw - 4, end - start, dh - 4],
                [dx, -d / 2 - (start + end) / 2, dh / 2],
                material,
            )


def screen_panels(sim, cfg, layout):
    pb = layout["screen_pb_mm"]
    if pb <= 0:
        return
    scr = screen_settings(cfg, layout)
    _, sc, n, t = booth(cfg, layout)
    rot_z = np.column_stack((n, t, [0, 0, 1]))
    width, height, sill = scr["width_mm"], scr["height_mm"], scr["sill_mm"]
    gap = scr["horizontal_joint_gap_mm"]
    top_thick = pb * scr["window_thickness_factor"]
    for mat, name, z0, z1, thick in (
        ("G4_Pb", "opaque", scr["floor_gap_mm"], sill - gap / 2, pb),
        (scr["window_material"], "window", sill + gap / 2, height, top_thick),
    ):
        box(sim, name, [thick, width, z1 - z0], sc + [0, 0, (z0 + z1) / 2], mat, rot_z)
    length = layout["return_mm"]
    if length > 0:
        corner = sc + t * (width / 2 + pb / 2) + n * (length / 2)
        box(
            sim,
            "booth_return",
            [length, pb, height],
            corner + [0, 0, height / 2],
            "G4_Pb",
            rot_z,
        )


def track_segments(sim, name, volume):
    from opengate.actors.filters import GateFilterBuilder

    a = sim.add_actor("PhaseSpaceActor", "tracks_" + name)
    a.attached_to = volume.name
    a.steps_to_store = "all"
    a.filter = GateFilterBuilder().ParticleName == "gamma"
    a.attributes = [
        "PDGCode",
        "ParentID",
        "PrePosition",
        "PostPosition",
        "PreDirection",
        "EventDirection",
        "PreKineticEnergy",
        "EventKineticEnergy",
        "Weight",
    ]
    a.output_filename = "tracks.root"
    return a

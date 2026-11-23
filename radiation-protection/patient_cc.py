"""Place the CC-compressed breast in the native frame of the Revelation arm."""

import hashlib
import json
import math
from pathlib import Path

import numpy as np


DATA = Path(__file__).resolve().parent / "anatomy"
MANIFEST = DATA / "cases.json"
PHANTOM_MODES = {
    "rad_fatty": ("rad_fatty", "segmented"),
    "rad_intermediate": ("rad_intermediate", "segmented"),
    "rad_dense": ("rad_dense", "segmented"),
    "rad_highest": ("rad_highest", "segmented"),
    "rad_cohort_min": ("rad_cohort_min", "segmented"),
    "rad_cohort_max": ("rad_cohort_max", "segmented"),
    "rad_uniform": ("rad_intermediate", "uniform_mixture"),
    "rad_fat_only": ("rad_fatty", "adipose_only"),
}

# Ivanov et al., Phys. Med. Biol. 63 (2018) 175020, Table 1.
# Mass fractions H, C, N, O, P and density in g/cm3.
ELEMENTS = ["H", "C", "N", "O", "P"]
TISSUES = {
    "RadAdipose": (0.93, [0.112, 0.619, 0.017, 0.251, 0.001]),
    "RadGlandular": (1.04, [0.102, 0.184, 0.032, 0.677, 0.005]),
    "RadSkin": (1.09, [0.098, 0.178, 0.050, 0.667, 0.007]),
}


def patient_case(choice):
    key, mode = PHANTOM_MODES[choice]
    case = json.loads(MANIFEST.read_text())["cases"][key]
    return case, mode


def image_path(choice):
    case, _ = patient_case(choice)
    return DATA / case["gate_image"]


def projection_for_scatter(view, angle_deg, custom=False):
    """Describe how the source breast is used at the requested arm angle."""
    if not math.isfinite(angle_deg) or not (
        (custom and view in ("CC", "MLO_R", "MLO_L") and abs(angle_deg) <= 90)
        or (view == "CC" and angle_deg == 0)
        or (view == "MLO_R" and 0 < angle_deg < 90)
        or (view == "MLO_L" and -90 < angle_deg < 0)
    ):
        raise ValueError("Use CC at 0°, MLO_R between 0° and +90°, or MLO_L between −90° and 0° in the study coordinate system.")
    return {
        "source_compression_view": "CC",
        "simulated_view": view,
        "arm_angle_deg": angle_deg,
        "anatomy_model": "native_cc_compression" if angle_deg == 0 else "rigid_cc_surrogate_for_mlo_scatter",
    }


def make_patient(sim, geometry, choice):
    import opengate as gate

    case, mode = patient_case(choice)
    path = image_path(choice)
    if hashlib.sha256(path.read_bytes()).hexdigest() != case["gate_sha256"]:
        raise ValueError(f"The supplied phantom image has changed: {path}")

    organ = sim.add_volume("Image", "rp_patient_breast")
    organ.mother = "world"
    organ.image = str(path)
    organ.material = "G4_AIR"
    organ.load_input_image()
    if list(organ.size_pix) != case["gate_dimensions_xyz"] or not np.allclose(
        organ.spacing, case["gate_spacing_mm"], atol=1e-8, rtol=0
    ):
        raise ValueError(f"Image geometry disagrees with the case manifest for {choice}.")

    width_x, width_y, thickness = map(float, organ.size_pix * organ.spacing)
    support_x, support_y = geometry.SUPPORT_X, geometry.SUPPORT_Y
    paddle_x, paddle_y = geometry.PADDLE_X, geometry.PADDLE_Y
    # Each file is cropped to the entire non-air bounding box; +X points to nipple.
    x0 = 0.0
    if (
        max(support_x[0], paddle_x[0]) > x0
        or min(support_x[1], paddle_x[1]) < x0 + width_x
        or max(support_y[0], paddle_y[0]) > -width_y / 2
        or min(support_y[1], paddle_y[1]) < width_y / 2
    ):
        raise ValueError("The full patient footprint exceeds this support or paddle.")
    base = float(geometry.geometry_z()["phantom"][0])
    organ.translation = [x0 + width_x / 2, 0, base + thickness / 2]
    paddle = sim.volume_manager.get_volume("paddle")
    paddle_base = base + thickness + geometry.GAP
    paddle.translation = [*paddle.translation[:2], paddle_base + geometry.PADDLE["mm"] / 2]

    density_unit = gate.g4_units.g / gate.g4_units.cm3
    db = sim.volume_manager.material_database
    for name, (rho, fractions) in TISSUES.items():
        db.add_material_weights(name, ELEMENTS, fractions, rho * density_unit)
    gland_fraction = case["glandular_volume_fraction"]
    materials = [[0, 1, "G4_AIR"], [1, 2, "RadAdipose"], [2, 3, "RadGlandular"], [3, 4, "RadSkin"]]
    mixture = None
    if mode == "adipose_only":
        materials[2][2] = "RadAdipose"
        gland_fraction = 0.0
    elif mode == "uniform_mixture":
        rho_a, wa = TISSUES["RadAdipose"]
        rho_g, wg = TISSUES["RadGlandular"]
        adipose_mass = (1 - gland_fraction) * rho_a
        gland_mass = gland_fraction * rho_g
        rho = adipose_mass + gland_mass
        weights = (adipose_mass * np.array(wa) + gland_mass * np.array(wg)) / rho
        db.add_material_weights("RadUniform", ELEMENTS, list(weights), rho * density_unit)
        materials[1][2] = materials[2][2] = "RadUniform"
        mixture = dict(density_g_cm3=rho, mass_fractions=dict(zip(ELEMENTS, weights)))
    organ.voxel_materials = materials
    return dict(
        patient_id=case["source_id"], material_mode=mode,
        source_compression_view="CC",
        original_glandular_volume_fraction=case["glandular_volume_fraction"],
        modelled_glandular_volume_fraction=gland_fraction,
        voxel_materials=materials, uniform_mixture=mixture,
        gate_dimensions_xyz=case["gate_dimensions_xyz"], spacing_mm=case["gate_spacing_mm"],
        tissue_extent_xyz_mm=[width_x, width_y, thickness],
        centre_in_scanner_mm=list(organ.translation),
        breast_base_mm=base, paddle_underside_mm=paddle_base,
        gate_image_sha256=case["gate_sha256"],
    )

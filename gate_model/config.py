"""
Parameters of the GATE 10 model of the Siemens MAMMOMAT Revelation (TiCEM) with the L2 v2 phantom. Background
materials selectable (BACKGROUNDS): PA12 + 5.5 wt% CaCO3 spheres in MCT oil (default) or PMMA spheres in paraffin oil.

Every value carries its source:
  [V] read from the DICOM headers of the L2 phantom images on this unit (all 15 raw/recombined images of 18 and 22 Aug
      2025 agree; checked 22 Oct 2026, including the Siemens private and CSA headers)
  [M] measured / derived in this project (REPORT.md, L2_container_v2)
  [A] assumption - not documented; replace with the real value when known

World coordinates (mm): x from the chest-wall face of the breast support into the phantom (x = 0 at the support face,
the phantom's positioning lip rests there), y along the chest wall (0 = centre of the phantom), z up.
The a-Se surface is the plane z = 0, the focal spot is at z = SID above the chest-wall edge (x = 0, y = 0).
"""
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent                                   # .../L2_background_analysis
V2_DIR = PROJECT / "L2_container_v2"

# ----------------------------------------------------------------------------------------------------------- system
SYSTEM = dict(
    sid_mm=656.0,                # [V] DistanceSourceToDetector
    source_to_entrance_mm=585.0,  # [V] DistanceSourceToEntrance (with 54 mm body-part thickness)
    support_to_detector_mm=17.0,  # [V] 656 - 585 - 54 (and 656 - 586 - 53 where 53 mm is read); confirmed independently by
                                  #     EstimatedRadiographicMagnificationFactor 1.0266 = 656 / 639 (magnification at the support)
    focal_spot_mm=0.3,            # [V] FocalSpots (nominal); modelled as a uniform disc of this diameter
    anode="W",                    # [V] AnodeTargetMaterial
    anode_angle_deg=16.0,         # [A] not in DICOM (value used in the earlier spectral model)
    be_window_mm=0.5,             # [A]
    extra_filtration_al_mm="fit",  # [M] tube filtration not otherwise modelled (mirror, window, ...), Al-equivalent.
                                  #     "fit": spectra.py chooses it so that the LE entrance dose (air + paddle included)
                                  #     equals the DICOM value; HE is the independent check. User decision 22 Oct 2026
                                  #     (no QC report with HVLs available). A number (mm) fixes it instead.
    pixel_mm=0.085,               # [V] ImagerPixelSpacing
    image_rows=3518,              # [V] along y (chest wall)
    image_cols=2800,              # [V] along x (chest wall -> nipple)
    detector_active_mm=(304.64, 239.36),  # [V] DetectorActiveDimensions (y, x) = 3584 x 2816 elements of 85 um
    detector_type="direct a-Se, ANRAD, contact mode",   # [V] DetectorType DIRECT, DetectorDescription "Contact;ANRAD",
                                                        #     DetectorID LV2-27960
    image_origin_on_detector_px=(33, 0),  # [V] FieldOfViewOrigin: the stored 3518 x 2800 image is centred along the chest
                                          #     wall (33-element margins) and starts at detector column 0 (side: see REPORT
                                          #     §1.3 - the rotation tag is not consistent with the pixels)
    missed_tissue_mm=3.0,         # [A] the image starts this far inside the support face; measure it (DESIGN_NOTES §8)
    a_se_um=200.0,                # [A] a-Se thickness: not in the headers (ANRAD direct a-Se [V])
    detector_cover=dict(material="CarbonFibre", mm=1.0),  # [A] detector entrance cover: not in the headers
)

# breast support and compression paddle. From the images [V]: support plane 17 mm above the detector (see above),
# TableID 1, PaddleID 14 (PaddleDescription "14"), compression force 49-74 N, thickness reading 53-54 mm. Thickness and
# material of both are NOT in the headers or measurable in the images (no paddle edge inside the field) -> [A]
SUPPORT = dict(material="CarbonFibre", mm=1.5, table_id=1)          # [A] material and thickness
PADDLE = dict(material="G4_POLYCARBONATE", mm=2.0, paddle_id=14)    # [A] material and thickness; ID 14 [V]

# anti-scatter grid: [V] present, focused, linear ("FOCUSED", "PARALLEL" in the Grid tag). No grid dimensions in any
# header. [M] The open field of all 6 unclipped raw images shows no residual grid-line peak in either direction (moving
# grid), so line density and orientation cannot be measured from the images either. All dimensions [A]
GRID = dict(
    enabled=True,
    lines_per_cm=31.0,            # [A] typical mammography value
    ratio=5.0,                    # [A] typical mammography value
    lamella_um=16.0,              # [A] lead lamella thickness
    lamella_material="G4_Pb",
    interspace_material="CarbonFibre",   # [A] (paper/cellulose is the other common choice)
    cover_material="CarbonFibre",
    cover_mm_each=0.5,            # [A] top and bottom cover plates
    focus_mm=656.0,               # [A] focal distance = SID
    lamellae_along="x",           # [A] lamellae parallel to the chest-wall -> nipple axis; selectivity in the y-z plane
)

# measurement-like images (measlike.py) and the comparison with the measured raw pair 0035 (compare_0035.py)
MEASURED = dict(
    dir_0035=PROJECT / "TICEM_RAW_L-CC_DIAGNOSIS_0035",   # raw LE (28 kV) and HE (49 kV) of 18 Aug 2025 [V]
    he_open_stored=5158.0,        # [M] HE open-field stored value in 0035 (REPORT §1.3: 5136-5180)
    ceiling=7500,                 # [M] hard ceiling of the stored values (LE open field is clipped)
)
# measured reference sets: exposure parameters, open-field level and where the images are. The July values come from
# the DICOMs (python measured_set.py --set july28high, needs the SanDisk drive) and are stored in data/.
MEASURED_SETS = {
    "0035": dict(description="old L2 container, old plate; unit 1024, 18 Aug 2025 (one LE/HE pair)",
                 unit="1024", kv=28, mAs=dict(LE=198.39, HE=38.76), esak_mGy=dict(LE=7.513, HE=0.722),
                 he_open_stored=5158.0, mirror=False, spectra_dir="spectra", default_pair="TICEM_RAW_L-CC_DIAGNOSIS_0035",
                 source=dict(kind="dir", path=PROJECT / "TICEM_RAW_L-CC_DIAGNOSIS_0035")),
    "july28high": dict(description="real L2 with the new plate; unit 10159, 28 kV high mAs, 12 pairs (21-30 Jul 2026)",
                       unit="10159", kv=28, mAs=None, esak_mGy=None, he_open_stored=None, mirror=True,
                       spectra_dir="spectra_july28high",
                       default_pair="HIGHMAS_28KV2_00",     # the image the disc template was made from
                       source=dict(kind="list", kvp=28, level="high",
                                   exposures=PROJECT / "investigation_0035_0037" / "scripts" / "july2026" / "exposures.json"),
                       values_json=HERE / "data" / "measured_july28high.json",
                       # disc positions in the (mirrored) image HIGHMAS_28KV2_00, used to find the disc array in each pair
                       template=PROJECT / "investigation_0035_0037" / "scripts" / "july2026" / "template.json"),
}
DEFAULT_MEASURED_SET = "0035"


def set_measured(name=None, pair=None):
    """select the measured reference set and ONE measured pair of it (each pair is treated on its own, never averaged):
    updates BEAMS (that pair's exact LE / HE mAs and entrance dose), MEASURED (its HE open-field level) and the spectra
    folder (spectra_dir()); returns the set name. Pair: --pair, REVELATION_PAIR or the set's default_pair."""
    import os, json as _json
    name = name or os.environ.get("REVELATION_SET") or PHANTOM_MODELS[PHANTOM["model"]].get("default_set",
                                                                                           DEFAULT_MEASURED_SET)
    if name not in MEASURED_SETS:
        raise SystemExit(f"unknown measured set '{name}'; choose one of {list(MEASURED_SETS)}")
    m = MEASURED_SETS[name]
    pair = pair or os.environ.get("REVELATION_PAIR") or m["default_pair"]
    vals = dict(mAs=m["mAs"], esak_mGy=m["esak_mGy"], he_open_stored=m["he_open_stored"])
    if m.get("values_json") and Path(m["values_json"]).exists():
        rows = {r["id"]: r for r in _json.load(open(m["values_json"]))["exposures"]}
        if pair not in rows:
            raise SystemExit("pair %s is not in set %s; pairs: %s" % (pair, name, ", ".join(rows)))
        r = rows[pair]
        vals = dict(mAs={b: r[b]["mAs"] for b in ("LE", "HE")}, esak_mGy={b: r[b]["esak_mGy"] for b in ("LE", "HE")},
                    he_open_stored=r["he_open_stored"])
    for b in ("LE", "HE"):
        if vals["mAs"]:
            BEAMS[b]["mAs"] = vals["mAs"][b]
        if vals["esak_mGy"]:
            BEAMS[b]["esak_mGy"] = vals["esak_mGy"][b]
    if vals["he_open_stored"]:
        MEASURED["he_open_stored"] = vals["he_open_stored"]
    MEASURED.update(set=name, pair=pair, set_ready=bool(vals["mAs"] and vals["esak_mGy"] and vals["he_open_stored"]))
    return name


def spectra_dir():
    """folder of the X-ray spectra for the selected measured set (env REVELATION_SPECTRA_DIR overrides)"""
    import os
    if os.environ.get("REVELATION_SPECTRA_DIR"):
        return Path(os.environ["REVELATION_SPECTRA_DIR"])
    return HERE / MEASURED_SETS[MEASURED.get("set", DEFAULT_MEASURED_SET)]["spectra_dir"]


DETECTOR_RESPONSE = dict(
    gain_mode="constant",         # [A] same stored units per deposited keV at LE and HE (LE scale follows from physics)
    # detector blur, applied to the signal and to the quantum noise (measlike.py):
    #   "nhsbsp2019": presampling MTF published for the Revelation's a-Se detector [R1] (it includes the 85 um pixel
    #                 aperture); separable in x and y [A]; the noise is shaped by the aliased MTF, sum_k MTF(f + k/a)^2
    #   "aperture":   ideal 85 um pixel aperture only, white noise (the model used before 23 Oct 2026)
    mtf="nhsbsp2019",
    # [R1] Table 10: edge on the breast support (0 mm), u direction (no tube motion), 29 kV W/Rh + 2 mm Al, tomosynthesis
    # projections without grid, detector LV2-00007, software VC10B. Assumed for our unit (LV2-33695, 2D, grid) [A].
    mtf_table=dict(f_cy_mm=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                   mtf=[1.00, 0.91, 0.86, 0.78, 0.67, 0.55, 0.44, 0.34, 0.25, 0.18, 0.13]),
    # noise model, variance of a stored pixel value: k x the quantum variance from the photon numbers (MTF-shaped)
    #   + (f S)^2 fixed-pattern (gain) noise + e^2 electronic noise, S = signal part of the stored value (no offset).
    # "calibrated": k, f, e fitted on the HE open fields of the set's exposures at low / normal / high mAs (July
    # 28 kV, unit 10159), python measlike.py --calibrate-noise -> data/detector_noise_<set>.json [M]. A dict
    # dict(k=..., f=..., e=...) fixes them. Same k, f, e at LE [A] (the LE open field is at the ceiling).
    noise_model="calibrated",
    # terms fitted on the open fields. Electronic noise is fixed at 0 [M]: a free fit gave e = 13 stored units, but the
    # measured LE phantom interior (stored ~620) bounds it to <= 9 units even if all its high-frequency noise were
    # electronic, and e = 13 tripled the simulated LE noise there (validation of 23 Oct 2026)
    noise_fit_terms=("quantum", "fixed_pattern"),
    fixed_pattern_shared=False,   # [A] True: the same fixed-pattern realisation in LE and HE (one gain map)
    # low-frequency structured noise (< 1 cy/mm), multiplicative: the excess of the measured HE open-field NNPS over the
    # quantum + fixed-pattern model, median over the calibration exposures [M] (stored by --calibrate-noise). Shared
    # by LE and HE [A]: supported by the test of 23 Oct 2026 (shared: residual at w* inside the measured range of the
    # 12 July arrangements; independent: about twice too high). None switches it off.
    structured_noise="calibrated",
    structured_shared=True,
    extra_blur_sigma_mm=0.0,      # [A] optional extra Gaussian blur on top of the MTF model
    # stored-value offset per beam [M]: REPORT §13.4, the LE offset that makes the July w* independent of mAs
    # (c_LE = 40 / 35 / 36 at 28 / 30 / 32 kV); no HE offset (the HE intercept of about -150 is not a detector offset)
    offset_stored=dict(LE={28: 40.0, 30: 35.0, 32: 36.0}, HE=0.0),
)
# References for the detector model
# [R1] NHS Breast Screening Programme. Technical evaluation of Siemens Revelation digital breast tomosynthesis system.
#      Public Health England / NCCPM, March 2019 (Table 10, MTF of the central projections). https://assets.publishing.
#      service.gov.uk/government/uploads/system/uploads/attachment_data/file/785503/Technical_evaluation_of_Siemens_
#      Revelation_digital_breast_tomosynthesis_system_report.pdf
# [R2] Zhao W, Ji WG, Debrie A, Rowlands JA. Imaging performance of amorphous selenium based flat-panel detectors for
#      digital mammography: characterization of a small area prototype detector. Med Phys 2003;30(2):254-263.
#      (85 um a-Se: K-fluorescence lowers the MTF at Nyquist by ~15 %; extra blur from charge trapping at the blocking
#      layer lowers the high-frequency MTF and NPS)
# [R3] Zhao W, Ji WG, Rowlands JA. Effects of characteristic x rays on the noise power spectra and detective quantum
#      efficiency of photoconductive x-ray detectors. Med Phys 28:2039-2049 (2001) [author order and year as listed by
#      the journal index; to confirm]. (K-fluorescence reabsorption lowers the high-frequency NPS of a-Se by ~12-18 %)
# Checked, no detector data: NHSBSP Practical Equipment Report "Practical evaluation of the Siemens Mammomat
# Revelation", Public Health England, August 2020 (reports no MTF, NNPS or DQE).

# ----------------------------------------------------------------------------------------------------------- beams
BEAMS = {
    # name: kVp, filter material, filter thickness (mm), mAs and entrance air kerma from the headers of 0035 [V]
    "LE": dict(kvp=28.0, filter=("Rh", 0.050), mAs=198.39, esak_mGy=7.513),
    "HE": dict(kvp=49.0, filter=("Ti", 1.000), mAs=38.76, esak_mGy=0.722),
}
# Siemens recombination weight _W per LE kV (private header (0019,1010)) [V/M]; HE is 49 kV in all series
SIEMENS_W = {28: 0.255, 30: 0.265, 32: 0.292, 34: 0.308}

# ----------------------------------------------------------------------------------------------------------- phantom
# geometry of the L2 v2 container (rev B, walls as measured) for the GATE model: a frozen copy of the
# L2_container_v2/params.json from which phantom/*.stl and the packing were made (22 Oct 2026: 3 mm plate, walls 6/4).
# The design has moved on since (2 mm plate, equal 5/5 walls, slot lugs; 27 Oct 2026); the GATE L2v2 model follows only
# when its meshes and packing are regenerated (phantom_meshes.py, sphere_packing.py) and this path is switched back.
V2_PARAMS = HERE / "data" / "L2v2_params_rev_B.json"
PHANTOM = dict(
    params_json=V2_PARAMS,
    disc_layout_csv=V2_DIR / "disc_layout.csv",
    sphere_fraction=0.61,         # [M] 61 % spheres / 39 % liquid by volume (as the current L2)
    sphere_diameters_mm=(15.88, 12.70, 9.52, 6.35, 3.18, 1.58),  # [A] published L1 set; L2: "6 sizes 1.6-16 mm"
    sphere_material=None,         # set by set_background() below
    liquid_material=None,
    background=None,
    container_material="G4_PLEXIGLASS",
    plate_material="G4_PLEXIGLASS",   # [A] new plate in PMMA
    # materials are defined in data/L2v2_materials.db (the names above refer to it or to Geant4 NIST materials).
    # Iodine discs: materials PMMA_I_row0..6 there (1.2 / 0.9 / 0.7 / 0.4 / 0.27 / 0.17 / 0.10 mg I/cm2, lower-bound
    # estimates [A], REPORT §13.2); change them with  python materials.py --iodine ...  The two 9 mm discs are plain
    # PMMA unless PMMA_I_large is defined (--large).
)

# phantom model = container, plate and sphere positions. Chosen with --phantom in every script, or with the environment
# variable REVELATION_PHANTOM. Geometry code: phantoms.py. [U] user statement, [M] measured from images, [A] assumed.
PHANTOM_MODELS = {
    "L2v2": dict(description="suggested redesign L2 v2 (rev B; 3 mm plate [A] - the real plate is 2 mm [U])",
                 geometry="v2", dir="phantom", params_json=V2_PARAMS, disc_layout=V2_DIR / "disc_layout.csv",
                 iodine_prefix="PMMA_I_", default_background="PA12_MCT", default_set="0035"),
    "L2": dict(description="real L2 container with the new (July 2026) plate",
               geometry="halfcylinder", dir="phantom_L2",
               R_out=120.0,          # [U] half-cylinder, outer diameter 24 cm; [M] outline in 0035
               H=53.5,               # [U] outer height
               t_top=6.0, t_bot=4.0,  # [U] lid and base, measured by the user
               t_wall=8.0,           # [U] curved wall, measured by the user
               t_front=8.0,          # [M] chest-wall wall: 8.2 mm band on the detector in 0035 (user 23 Oct: keep)
               t_plate=2.0,          # [U] "2 mm PMMA plate"; the July plate has the same thickness and outline
               plate_zc=None,        # [U] exactly in the middle of the cavity (user 23 Oct 2026): (t_bot + H - t_top) / 2
               plate_clearance=0.5,  # [A] plate = half-disc 0.5 mm inside the curved wall, from the front wall
                                     #     (user 23 Oct: use the assumptions; sphere sizes and densities likewise)
               # position in the world and disc positions from the July images [M]: the outer front face is at
               # column 56 of the 28 kV template (july2026/template.json), back-projected to mid-height (M = 1.072),
               # as for the disc layout (x from that face, y from the centre of curvature, both at M = 1.072).
               # The discs are placed by back-projecting their measured image positions to the plate plane.
               layout_ref=dict(col_px=56.0, magnification=1.072),
               front_face_x_mm=(SYSTEM["missed_tissue_mm"] + 56.0 * SYSTEM["pixel_mm"]) / 1.072,   # ~7.24 mm
               disc_layout=V2_DIR / "disc_layout_measured.json",
               iodine_prefix="PMMA_I_L2_", default_background="PMMA_paraffin", default_set="july28high"),
}
DEFAULT_PHANTOM = "L2v2"


def set_phantom(name=None):
    """select the phantom model (updates PHANTOM["model"]); returns its name"""
    import os
    name = name or os.environ.get("REVELATION_PHANTOM", DEFAULT_PHANTOM)
    if name not in PHANTOM_MODELS:
        raise SystemExit(f"unknown phantom '{name}'; choose one of {list(PHANTOM_MODELS)}")
    PHANTOM["model"] = name
    return name


set_phantom()

# background = sphere and liquid materials (both defined in data/L2v2_materials.db). Chosen with --background in every
# script, or with the environment variable REVELATION_BACKGROUND; same container, plate and sphere positions for all.
BACKGROUNDS = {
    "PA12_MCT": dict(sphere_material="PA12_CaCO3_5p5", liquid_material="MCT_oil",
                     description="suggested v2 background: PA12 + 5.5 wt% CaCO3 spheres in MCT oil"),
    "PMMA_paraffin": dict(sphere_material="PMMA_sphere", liquid_material="paraffin_oil",
                          description="current L2 background: PMMA spheres in liquid paraffin"),
    "gland_adipose": dict(sphere_material="glandular_tissue", liquid_material="adipose_tissue",
                          description="tissue reference (not buildable): glandular-tissue spheres in adipose tissue, Hammerstein compositions; the REPORT §10 target"),
}
DEFAULT_BACKGROUND = "PA12_MCT"


def set_background(name=None):
    """select the sphere/liquid materials (updates PHANTOM); returns the background name"""
    import os
    name = name or os.environ.get("REVELATION_BACKGROUND") or PHANTOM_MODELS[PHANTOM["model"]]["default_background"]
    if name not in BACKGROUNDS:
        raise SystemExit(f"unknown background '{name}'; choose one of {list(BACKGROUNDS)}")
    PHANTOM.update(sphere_material=BACKGROUNDS[name]["sphere_material"],
                   liquid_material=BACKGROUNDS[name]["liquid_material"], background=name)
    return name


set_background()

# ----------------------------------------------------------------------------------------------------------- materials
# All non-NIST materials (spheres, oil, carbon fibre, a-Se, iodine discs) are in data/L2v2_materials.db (GATE format),
# read by GATE and by materials.py. Compositions of the Geant4 NIST materials, needed by the analytic parts:
NIST_COMPOSITIONS = {
    "G4_PLEXIGLASS":    dict(density=1.19,  weights={"H": 0.080538, "C": 0.599848, "O": 0.319614}),
    "G4_POLYCARBONATE": dict(density=1.20,  weights={"H": 0.055491, "C": 0.755751, "O": 0.188758}),
    "G4_AIR":           dict(density=0.00120479, weights={"C": 0.000124, "N": 0.755268, "O": 0.231781, "Ar": 0.012827}),
    "G4_Pb":            dict(density=11.35, weights={"Pb": 1.0}),
}

# ----------------------------------------------------------------------------------------------------------- simulation
SIM = dict(
    physics_list="G4EmLivermorePhysics",   # Rayleigh + photoelectric + Compton with Livermore data
    fluorescence=True,
    production_cut_mm=0.1,
    threads=8,
    # stored per photon at the phase-space plane (no event-level attributes: keeps the ROOT files small)
    phsp_attributes=["KineticEnergy", "PrePosition", "PreDirection", "UnscatteredPrimaryFlag", "PDGCode"],
)
import os as _os
OUT_DIR = Path(_os.environ.get("REVELATION_OUT_DIR", HERE / "output"))
# X-ray spectra written by spectra.py; the environment variables let a test use other spectra/outputs
SPECTRA_DIR = Path(_os.environ.get("REVELATION_SPECTRA_DIR", HERE / "spectra"))

set_measured()


def add_case_args(ap, background=True, measured=True):
    """command-line choice of the simulated case: --phantom, --background, --set (each also by environment variable)"""
    ap.add_argument("--phantom", choices=list(PHANTOM_MODELS), default=None,
                    help="phantom model: " + "; ".join("%s = %s" % (k, v["description"].replace("%", "%%")) for k, v in PHANTOM_MODELS.items())
                         + " (default: REVELATION_PHANTOM or %s)" % DEFAULT_PHANTOM)
    if background:
        ap.add_argument("--background", choices=list(BACKGROUNDS), default=None,
                        help="sphere/liquid materials: " + "; ".join("%s = %s" % (k, v["description"].replace("%", "%%")) for k, v in BACKGROUNDS.items())
                             + " (default: REVELATION_BACKGROUND or the phantom model's default)")
    if measured:
        ap.add_argument("--set", choices=list(MEASURED_SETS), default=None,
                        help="measured reference (spectra, mAs, open-field level): " +
                             "; ".join("%s = %s" % (k, v["description"].replace("%", "%%")) for k, v in MEASURED_SETS.items()) +
                             " (default: REVELATION_SET or the phantom model's default)")
        ap.add_argument("--pair", default=None,
                        help="the one measured pair to reproduce (id as in measured_set.py; default: REVELATION_PAIR or "
                             "the set's default, July: HIGHMAS_28KV2_00). Each pair is treated on its own.")


def apply_case_args(a):
    """apply --phantom / --background / --set / --pair (phantom first: the others default per phantom); returns
    case_name()"""
    set_phantom(getattr(a, "phantom", None))
    set_background(getattr(a, "background", None))
    set_measured(getattr(a, "set", None), getattr(a, "pair", None))
    return case_name()


def case_name():
    """name of the simulated case in output folders: <background> for the L2 v2 (as before), L2_<background> for the
    real L2"""
    return PHANTOM["background"] if PHANTOM["model"] == "L2v2" else "%s_%s" % (PHANTOM["model"], PHANTOM["background"])

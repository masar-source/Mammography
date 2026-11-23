# Parameters of the GATE 10 model: Siemens MAMMOMAT Revelation (TiCEM) with the real L2 phantom

Status 23 Oct 2026. This lists every parameter of the current simulation case and where it came from:
- phantom `L2` (real container with the new July 2026 plate);
- background `PMMA_paraffin`;
- measured reference set `july28high`, pair `HIGHMAS_28KV2_00`.

"Code" names the place where each value is set: `config.py` unless stated otherwise.

**Source codes:**
- **[V]** read from the DICOM headers (unit and date given);
- **[P]** published: Salomon et al. 2020 (SPIE) and 2023 (Med Phys), PDFs in the project folder;
- **[U]** stated or measured by the user in this project (date given);
- **[M]** measured or derived in this project (REPORT.md, project notes, image analysis);
- **[R1–R3]** literature, see References;
- **[F]** fitted to measured data, method given;
- **[A]** assumption: not documented, to be replaced when known;
- **[N]** numerical or model choice.

---

## 1. Coordinates

| Item | Value | Source |
|---|---|---|
| Origin and axes | a-Se surface at z = 0; x from the chest-wall face of the breast support towards the nipple; y along the chest wall (0 = image centre); z up | [N] |
| Focal spot | (0, 0, SID), above the chest-wall edge | [N] from the geometry below |
| Image grid | columns along x from x = 3 mm (pixel centres at 3 + (i + 0.5)·0.085 mm); rows along y, centred | [N]; 3 mm = missed tissue [A], see §2 |
| Orientation of the July images | R-CC; mirrored left–right so that the chest wall is at column 0 | [V] ImageLaterality R; [M] as in july2026/pipe.py |

## 2. X-ray system geometry and detector format

| Parameter | Value | Source | Code |
|---|---|---|---|
| Source–detector distance (SID) | 656 mm | [V] DistanceSourceToDetector (unit 1024, Aug 2025; identical on unit 10159, July 2026) | `SYSTEM["sid_mm"]` |
| Source–entrance distance | 585 mm (at 54 mm thickness) | [V] DistanceSourceToEntrance (both units) | `SYSTEM["source_to_entrance_mm"]` |
| Support surface above detector | 17 mm | [V] 656 − 585 − 54; confirmed by EstimatedRadiographicMagnificationFactor 1.0266 = 656/639 | `SYSTEM["support_to_detector_mm"]` |
| Focal spot | 0.3 mm, modelled as a uniform disc | [V] FocalSpots (nominal); disc shape [A] | `SYSTEM["focal_spot_mm"]` |
| Pixel spacing | 0.085 mm | [V] ImagerPixelSpacing | `SYSTEM["pixel_mm"]` |
| Stored image | 3518 rows (y) × 2800 columns (x) | [V] Rows, Columns | `SYSTEM["image_rows"/"image_cols"]` |
| Detector active area | 304.64 × 239.36 mm (3584 × 2816 elements) | [V] DetectorActiveDimensions | `SYSTEM["detector_active_mm"]` |
| Image origin on the detector | element (33, 0) | [V] FieldOfViewOrigin | `SYSTEM["image_origin_on_detector_px"]` |
| Detector type | direct a-Se, ANRAD, contact mode | [V] DetectorType, DetectorDescription | `SYSTEM["detector_type"]` |
| Detector IDs | LV2-27960 (unit 1024); LV2-33695 (unit 10159, July) | [V] DetectorID | data/measured_july28high.json |
| Missed tissue at the chest wall (image starts inside the support face) | 3 mm | [A] | `SYSTEM["missed_tissue_mm"]` |
| Collimated field at the detector | x 0–245 mm, y ±155 mm | [A]; [M] no collimator edge in the images | revelation_gate.py `FIELD_X/Y` |
| Virtual field aperture | plane z = 200 mm, the field scaled to that plane | [N] | revelation_gate.py `APERTURE_Z` |

## 3. Breast support, compression paddle, below-plane stack

| Parameter | Value | Source | Code |
|---|---|---|---|
| Breast support | carbon fibre, 1.5 mm, z 15.49–16.99 mm, x 0–265, y ±165 mm | thickness, material and extent [A]; TableID 1 [V] | `SUPPORT`, revelation_gate.py |
| Compression paddle | polycarbonate (G4_POLYCARBONATE), 2.0 mm, z 70.51–72.51 mm, x −25 to 240, y ±155 mm | thickness, material and extent [A]; PaddleID 14 [V]; placed on the phantom top | `PADDLE`, revelation_gate.py |
| Gap between touching volumes | 0.01 mm | [N] avoids coincident surfaces | revelation_gate.py `GAP` |
| Phase-space plane | z 14.0–14.4 mm (x −60 to 300, y ±210) | [N] | revelation_gate.py `PLANE_Z` |
| Anti-scatter grid | focused, linear, moving | [V] Grid tag "FOCUSED", "PARALLEL"; [M] no grid lines in the images (moving) | `GRID` |
| Grid: line density / ratio / lamellae | 31 lines/cm, ratio 5, 16 µm lead; lamella height 1.533 mm | [A] typical values; height [N] = ratio × (period − lamella) | `GRID`, grid_detector.py |
| Grid: interspace, covers, focus | carbon fibre; 2 × 0.5 mm carbon-fibre covers; focused at 656 mm; lamellae along x | [A] | `GRID` |
| Grid position | top of assembly at z = 6.533 mm | [A] | grid_detector.py `GRID_TOP` |
| Detector cover | carbon fibre 1.0 mm at z 1–2 mm | [A] | `SYSTEM["detector_cover"]`, grid_detector.py |
| a-Se layer | 200 µm, energy-integrating; Se K-fluorescence escape neglected | thickness [A]; type [V] | `SYSTEM["a_se_um"]` |

## 4. Real L2 phantom: container and plate

Phantom coordinates: x from the outer chest-wall (front) face, y from the axis of the half-cylinder, z from the bottom outer surface.

| Parameter | Value | Source | Code |
|---|---|---|---|
| Shape | half-cylinder; its axis lies in the flat chest-wall face | [P] "half cylinder"; [M] centre of curvature on the front face in image 0035 | `PHANTOM_MODELS["L2"]`, phantoms.py |
| Outer diameter / radius | 240 mm / 120 mm | [P]; [U]; [M] outline in 0035 falls at 134.5 mm = 120 × 1.121 | `R_out` |
| Outer height | 53.5 mm | [P]; [U] 21 Oct 2026 | `H` |
| Lid (top) | 6 mm PMMA | [U] measured, 22 Oct 2026 | `t_top` |
| Base (bottom) | 4 mm PMMA | [U] measured, 22 Oct 2026 | `t_bot` |
| Curved wall | 8 mm PMMA | [U] measured, 22 Oct 2026 | `t_wall` |
| Chest-wall (front) wall | 8 mm PMMA | [M] 8.2 mm band in image 0035; kept by the user on 23 Oct 2026 | `t_front` |
| Container material | PMMA, G4_PLEXIGLASS, 1.19 g/cm³ | material [U] (this project; the papers [P] do not state it); density = NIST value [A] | `PHANTOM["container_material"]` |
| Cavity (filling) | radius 112 mm, x ≥ 8 mm, z 4–47.5 mm; 779.2 ml | [N] derived from the above | phantoms.py `l2()` |
| Plate | PMMA, 2 mm thick | [U] "2 mm PMMA plate"; the July plate has the same outline and thickness as the old one (21 Oct 2026) | `t_plate` |
| Plate height | centre at z = 25.75 mm, exactly in the middle of the cavity | [U] 23 Oct 2026 | `plate_zc` (derived) |
| Plate outline | half-disc of radius 111.5 mm (0.5 mm inside the curved wall), from x = 8.05 mm (the front wall) | [A] (the drawing is inaccurate, per the user) | `plate_clearance` |
| How the plate is held | fixed at mid-height; holes and supports not modelled | [A] | — |
| Position in the world | front face at x = 7.24 mm; axis at y = 0 | [M] column 56 of the July template, (3 + 56 × 0.085)/1.072 | `front_face_x_mm` |
| GATE meshes | prisms of a circular segment, sagitta ≤ 0.005 mm; 696 / 640 / 640 facets (envelope / filling / plate); closed; volumes within 0.006 % of the exact solids | [N] | phantoms.py, phantom_L2/*.stl |

## 5. Real L2 phantom: spheres and liquid

| Parameter | Value | Source | Code |
|---|---|---|---|
| Fill | 61 % PMMA spheres / 39 % liquid paraffin oil by volume | [P] Salomon 2020 and 2023; [U] | `PHANTOM["sphere_fraction"]` |
| Sphere sizes | 15.88 / 12.70 / 9.52 / 6.35 / 3.18 / 1.58 mm (standard inch balls) | [P] "six diameters 1.6–16 mm"; exact sizes [A] | `PHANTOM["sphere_diameters_mm"]` |
| Size distribution | equal volume per size | [P] "equal volumes" | sphere_packing.py |
| Sphere counts (seed 1) | 36 / 71 / 167 / 564 / 4490 / 36615 = 41 943 spheres; 61.0 % of the 743.8 ml free volume | [N] packing result | phantom_L2/spheres_summary.json |
| Packing method | random packing by overlap relaxation with growing radii; no gravity settling | [N]; homogeneity [A] | sphere_packing.py |
| Packing clearances | 0.02 mm to walls and plate (above the 0.005 mm mesh deviation), 0.001 mm between spheres | [N] | sphere_packing.py `WALL`, `GAP` |
| Sphere material | PMMA (`PMMA_sphere` = G4_PLEXIGLASS), 1.19 g/cm³ | material [U]/[P]; density [A] (NIST) | data/L2v2_materials.db |
| Liquid | liquid paraffin (`paraffin_oil`): 0.86 g/cm³, C 0.852 / H 0.148 by mass | material [P] (paraffin oil, WHC GmbH, Salomon 2023); density and composition [A] (Ph. Eur. range 0.827–0.890 g/cm³, REPORT §6) | data/L2v2_materials.db |
| Cross-check: attenuation at 20 keV | model: paraffin 0.371, PMMA 0.680 cm⁻¹. Published: PMMA 0.68 (literature), paraffin 0.34 measured (Salomon 2023) or 0.29 estimated (Salomon 2020) | [P] |
| | With our composition, 0.34 cm⁻¹ would need 0.79 g/cm³, below the Ph. Eur. range. It would raise the sphere contrast, i.e. the texture excess. Not adopted; recorded as an open discrepancy | [M] |
| Rearrangement | spheres rearranged between exposures; one packing per simulation | [U]; seed 1 [N] | — |

## 6. Iodine discs (new plate)

| Parameter | Value | Source | Code |
|---|---|---|---|
| Array | 7 rows × 5 diameters 2 / 3 / 4 / 5 / 7 mm, plus two 9 mm discs | [M] July 2026 images (REPORT §13) | L2_container_v2/disc_layout_measured.json |
| Positions | measured image positions in the July 28 kV template (HIGHMAS_28KV2_00), back-projected to the plate plane (M = 656/(656 − 17 − 25.75)); array x 21.4–77.7 mm from the front face, y −46.1 to +45.2 mm | [M] | phantoms.py `discs()`, `layout_ref` |
| Disc thickness | through the 2 mm plate (vertical cylinders, 0.002 mm short of the plate faces) | [A] | revelation_gate.py, projector.py |
| Iodine per row (row 0 to 6) | 1.2 / 0.9 / 0.7 / 0.4 / 0.27 / 0.17 / 0.10 mg I/cm² | [M] lower-bound estimates from the July plate images (REPORT §13.2) | data/L2v2_materials.db (`PMMA_I_L2_row0..6`) |
| Disc material | PMMA + iodine by mass fraction over 2 mm (6.0 to 0.5 mg I/cm³; densities 1.1960 to 1.1905 g/cm³) | [N] mixture rule | materials.py |
| 9 mm discs | plain PMMA (no iodine) | [A] content unknown | — |
| Marker | not modelled (only excluded in the analysis) | [N] | — |

## 7. X-ray spectra (SpekPy)

| Parameter | Value | Source | Code |
|---|---|---|---|
| Anode | tungsten (W) | [V] AnodeTargetMaterial TUNGSTEN (both units) | `SYSTEM["anode"]` |
| Anode angle | 16° | [A] not in DICOM | `SYSTEM["anode_angle_deg"]` |
| Be window | 0.5 mm | [A] | `SYSTEM["be_window_mm"]` |
| LE beam | 28 kV, rhodium 0.050 mm | [V] KVP, FilterMaterial RHODIUM, FilterThickness 0.05 | `BEAMS["LE"]` |
| HE beam | 49 kV, titanium 1.000 mm | [V] KVP, FilterMaterial TITANIUM, FilterThickness 1.0 | `BEAMS["HE"]` |
| Extra Al-equivalent filtration | 0.136 mm, both beams | [F] chosen so that the modelled LE entrance dose (air + paddle included) equals the DICOM value of the pair; user decision 22 Oct 2026 (no QC report) | spectra.py, spectra_july28high/spectra_info.json |
| Entrance dose, model / DICOM | LE 7.597 / 7.598 mGy (1.00, by the fit); HE 0.704 / 0.771 mGy (0.91, independent check) | [M] | spectra_info.json |
| Half-value layers (model) | LE 0.50 mm Al; HE 3.04 mm Al | [N] from the spectra | spectra_info.json |
| Energy bins | 0.5 keV | [N] | spectra.py |
| Siemens recombination weight W | 0.255 at 28 kV (0.265 / 0.292 / 0.308 at 30 / 32 / 34 kV) | [V] private header (0019,1010) `_W(...)` | `SIEMENS_W` |

## 8. Exposure of the measured pair reproduced (HIGHMAS_28KV2_00)

| Parameter | Value | Source |
|---|---|---|
| X-ray unit / date | serial 10159, 28 Jul 2026, R-CC | [V] |
| LE / HE tube load | 200.64 / 41.39 mAs | [V] ExposureInuAs of the LE and of the HE image |
| LE / HE entrance dose | 7.598 / 0.771 mGy | [V] EntranceDoseInmGy |
| Compressed thickness / force | 54 mm / 46.9 N | [V] BodyPartThickness, CompressionForce |
| HE open-field stored value | 5784 | [M] median of the open field next to the phantom (measured_set.py) |
| Stored-value ceiling | 7500 | [M] (the LE open field is clipped) |
| Other pairs | the other 11 high-mAs pairs are available with `--pair` (196–221 mAs LE, 41–47 mAs HE, open level 5779–6572) | [V]/[M] data/measured_july28high.json |

## 9. Detector response (measurement-like images, measlike.py)

| Parameter | Value | Source | Code |
|---|---|---|---|
| Signal | stored value ∝ deposited energy (flat-fielded), scaled to the measured HE open level; same gain per keV for LE | [A] constant gain | `gain_mode` |
| Presampling MTF | 1.00 / 0.91 / 0.86 / 0.78 / 0.67 / 0.55 / 0.44 / 0.34 / 0.25 / 0.18 / 0.13 at 0–10 cy/mm; separable in x and y | [R1] Table 10 (edge on the breast support, u-direction); separability [A]; measured on another Revelation (LV2-00007) [A] | `mtf`, `mtf_table` |
| Photon-noise factor k | 1.78 × the photon noise from the photon numbers; shaped by the aliased MTF | [F] fitted on the HE open fields of all 36 July 28 kV exposures (low / normal / high mAs), 0.6–5.9 cy/mm, rms misfit 8.5 % | data/detector_noise_july28high.json |
| Fixed-pattern noise f | 0.28 % of the signal (rms), white, independent in LE and HE | [F] same fit; independence [A] | same |
| Electronic noise e | 0 | [M] a free fit gave 13 units; the measured LE phantom interior allows at most ≈ 9 | `noise_fit_terms` |
| Low-frequency structured noise | measured excess spectrum: 1.40·10⁻⁵ / 7.0·10⁻⁶ / 3.3·10⁻⁶ / 1.8·10⁻⁶ / 7.9·10⁻⁷ / 3.0·10⁻⁷ / 4.6·10⁻⁸ mm² in the bands 0.05–0.1 / 0.1–0.15 / 0.15–0.2 / 0.2–0.3 / 0.3–0.45 / 0.45–0.6 / 0.6–1.0 cy/mm; multiplicative; the same field in LE and HE | [M] median over the 36 open fields; sharing [A], supported by the w* residual test | `structured_noise`, `structured_shared` |
| Noise at LE | same k, f and structured spectrum as HE | [A] the LE open field is at the ceiling | — |
| Offset | LE + 40 stored units at 28 kV (35 / 36 at 30 / 32 kV); HE 0 | [M] REPORT §13.4 (offset that makes w* independent of mAs; the HE intercept is not a detector offset) | `offset_stored` |
| Scatter noise | scatter assumed to have the primary's mean detected energy | [A] | measlike.py |
| Noise seed | 1 + pair index | [N] | measlike.py |

## 10. Monte Carlo (GATE 10.1.1)

| Parameter | Value | Source | Code |
|---|---|---|---|
| Software | opengate 10.1.1 (Geant4 11.4 data), Python 3.10, novis build | [M] workstation | — |
| Physics | G4EmLivermorePhysics, fluorescence on, production cut 0.1 mm | [N] | `SIM` |
| Electrons | killed at creation in support, paddle and phantom (min. kinetic energy 1 MeV) | [N] range < 0.05 mm; bremsstrahlung neglected | revelation_gate.py |
| World | 1.2 × 1.2 × 1.6 m, air | [N] | revelation_gate.py |
| Source | isotropic from a 0.3 mm disc at the focal spot; cone θ 0–12.85°, φ 90–270° (x ≥ 0); acceptance aimed at the phantom envelope, rejected directions emitted with zero energy | [N] | revelation_gate.py `source_cone` |
| Photons | 5 runs × 2·10⁷ (LE) and 5 runs × 4·10⁶ (HE), seeds 1–5, 8 threads | [N] | run_all.sh |
| Recorded | every photon crossing the phase-space plane: energy, position, direction, unscattered flag, particle type | [N] | `SIM["phsp_attributes"]` |
| Sanity check | unscattered photons trace back to the focal spot (all 10 runs: 100 %, max. offset 0.15 mm) | [M] | revelation_gate.py |

## 11. Image formation

| Parameter | Value | Source | Code |
|---|---|---|---|
| Primary image | analytic, noise-free, polychromatic; exact chords through container, plate, spheres and discs; focal-spot blur neglected | [N] | projector.py |
| Scatter image | Monte Carlo photons (scattered or fluorescence) weighted by grid, cover and a-Se response; 1 mm bins; Gaussian smoothing σ = 5 mm | [N] | postprocess.py |
| Open field | the same chain with the phantom replaced by air | [N] | projector.py |
| Normalisation | per emitted photon and steradian; detected photons per sr = SpekPy fluence per mAs at 1 m × 10⁴ × mAs | [N] | measlike.py |

## 12. Analysis applied to measured and simulated images

| Parameter | Value | Source |
|---|---|---|
| Recombination | −ln(stored value); H − W·L (W = 0.255), H − w*·L, split weight (LE smoothed 2 px) | [M] REPORT §5, §13 |
| w* | pipe.py: cross-spectral regression of H on L in the bands 0.06–0.15 cy/mm, 128 px tiles ≥ 300 px from the edge. compare_measured.py / analysis.py: regression of band-passed images (difference of Gaussians 2–35 px), ≥ 30 mm from the edge | [M] REPORT §5, §13, july2026/pipe.py |
| Phantom footprint | HE < 3300 stored units (pipe.py) or < 0.64 × HE open field (compare_measured.py) | [M] REPORT |
| Texture | 256 px (21.8 mm) ROIs, quadratic detrend, Hann window; β fitted at 0.2–0.7 cy/mm; power at 0.2 cy/mm | [M] REPORT §11 |
| Rim profile | ring medians vs distance from the fitted curved edge and from the chest wall | [M] scripts/shading.py method |
| Disc registration | measured: pipe.py estimate (outline + summed H − 0.37·L response; (3, 3) px for HIGHMAS_28KV2_00, independent locator within 1 px). Simulated twins: true positions from the projector iodine map (pipe.py's estimate was 6–14 px off) | [M] pipeline_compare.py, 24 Oct 2026 |
| Discs | contrast = mean within 0.6 r minus annulus r+15 to r+45 px; noise = SD of 12 background ROIs; SDNR thresholds 3 and 5 | [M] july2026/pipe.py |
| Siemens emulation (for the simulated recombined image) | fitted LE filter K(f), band gain A(f), gain and offset | [M] REPORT §12, siemens_cem_emulation_params.json |
| Redesign prediction (redesign_compare.py) | measured current × (sim redesign / sim current) for positive quantities; + (sim redesign − sim current) for β and rim; r from the predicted w* with the header W; disc SDNR scaled by contrast / noise ratios | [A] errors common to both simulations cancel |
| Redesign targets | r_tissue +0.07 / +0.09 / +0.08 / +0.12 at 28 / 30 / 32 / 34 kV; LE texture power 0.94 × current; β unchanged | [E] REPORT §10.1, §10.3 |

## 13. Options in the code that the L2 case does not use

| Item | Value | Source |
|---|---|---|
| L2 v2 redesign | revolved container, rev B as of 22 Oct 2026 (frozen copy `data/L2v2_params_rev_B.json`: 3 mm plate, walls 6/4, no slot lugs), 61 799 spheres. The design has since changed (2 mm plate, 5/5 walls, slot lugs); not yet in the GATE model | [M] L2_container_v2 |
| Tissue reference (option gland_adipose, not buildable) | glandular 1.04 g/cm³ (H .102 C .184 N .032 O .677 P .005); adipose 0.93 g/cm³ (H .112 C .619 N .017 O .251 P .001) | [M] Hammerstein 1979 as used by the REPORT model (scripts/model.py; REPORT [R6], secondary table, primary not accessed; trace elements as P) |
| PA12 + CaCO₃ spheres in MCT oil | 1.046 g/cm³ (C .697 H .111 O .103 N .067 Ca .022); MCT 0.945 g/cm³ (C .698 H .109 O .192) | [M] design values (MATERIALS.md) |
| Measured set 0035 | unit 1024, 18 Aug 2025: 198.39 / 38.76 mAs, 7.513 / 0.722 mGy, HE open 5158 | [V]/[M] |
| Carbon fibre (support, grid, cover) | 1.55 g/cm³, C 0.85 / H 0.05 / O 0.10 | [A] |
| a-Se | 4.27 g/cm³ | [N] standard value |

## References

- **[P] Salomon et al. 2020 (SPIE):** "Equivalent breast thickness and dose sensitivity of a next iteration 3D structured breast phantom with lesion models" (PDF in the project folder).
- **[P] Salomon et al. 2023 (Med Phys):** "Technical note: Realization and uncertainty analysis for an adjustable 3D structured breast phantom…" (PDF in the project folder).
- **[R1] NHS Breast Screening Programme:** *Technical evaluation of Siemens Revelation digital breast tomosynthesis system.* Public Health England / NCCPM, March 2019 (Table 10). https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/785503/Technical_evaluation_of_Siemens_Revelation_digital_breast_tomosynthesis_system_report.pdf
- **[R2] Zhao W, Ji WG, Debrie A, Rowlands JA:** Imaging performance of amorphous selenium based flat-panel detectors for digital mammography: characterization of a small area prototype detector. *Med Phys* 2003;30(2):254–263. Abstract level; background only.
- **[R3] Zhao W, Ji WG, Rowlands JA:** Effects of characteristic x rays on the noise power spectra and detective quantum efficiency of photoconductive x-ray detectors. *Med Phys* 28:2039–2049 (2001). Abstract level; author order and year to confirm; background only.
- **REPORT.md** (investigation_0035_0037) and **PROJECT_NOTES.md**: this project's own measurements and decisions.

## Known limitations (inside the phantom)

- **HE interior level:** 6–8 % low. Not explained by realistic material, a-Se or filtration ranges. The grid (all [A]), anode angle and Be window are untested.
- **Sphere texture at 0.2 cy/mm:** ≈ 1.2× the highest and ≈ 1.4× the median of the 12 measured arrangements. It is the same for two simulated packings and not explained by realistic densities.
- **w\*:** 0.351 simulated against 0.354–0.360 measured.
- **Paraffin attenuation:** the model has 0.371 cm⁻¹ at 20 keV against 0.34 measured by Salomon 2023 (see §5).
- **Outside the phantom (accepted):** the open-field noise shape, the outer edge zone and Siemens' border processing.

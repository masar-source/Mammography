# GATE 10 model: Siemens MAMMOMAT Revelation (TiCEM) with the L2 phantoms

A Monte Carlo model of one contrast-enhanced exposure pair (LE 28 kV Rh, HE 49 kV Ti) on the Revelation, with one of two phantom models (`--phantom`, see "Phantom models"):
- **`L2`: the real L2 container with the new (July 2026) plate.** This is the phantom of the July measurements, so the simulation can be compared with them one to one (`compare_measured.py`).
- **`L2v2`: the suggested redesign** (rev B, walls as measured; default).

Both have spheres in a liquid (61 : 39 by volume, six sizes of equal volume) in one of two backgrounds (see "Background options"):
- **PMMA spheres in liquid paraffin** (the real L2 materials; default for `L2`);
- **PA12 + 5.5 wt % CaCO₃ spheres in MCT oil** (suggested design; default for `L2v2`).

Each has its iodine plate in the middle.

**Status:** code only. GATE was not installed or run here, by request. The non-GATE parts were run and checked locally:
- geometry and sphere packing;
- grid/detector model;
- primary projector;
- normalisation between Monte Carlo and projector (`selftest.py`, all checks pass).

The GATE script was checked call by call against the source of the opengate **10.1.1** release (the version installed on the workstation, in a Python 3.10 virtual environment) and with a mock of its API, not with GATE itself. The material file was read with 10.1.1's own parser, and all scripts parse under Python 3.10 rules. The first real run must still be checked (see "First run").

## How the model works

| Part | Where | How |
|---|---|---|
| Spectra | `spectra.py` | SpekPy: W anode, Be window, Rh 50 µm / Ti 1 mm. Air and paddle are simulated in GATE. |
| Source, air, paddle, phantom, support | `revelation_gate.py` (GATE 10) | Isotropic emission from a 0.3 mm disc at the focal spot, aimed at the phantom by GATE's angular acceptance. Livermore photon physics with fluorescence. |
| Phase space | `revelation_gate.py` | Every photon crossing a plane 1 mm below the support: energy, position, direction, scattered or not. A kill volume stops photons below it. |
| Anti-scatter grid, covers, a-Se | `grid_detector.py` | Applied per photon after the run. Focused linear grid, averaged over lamella phase (moving grid). Energy-integrating a-Se absorption. |
| Primary image | `projector.py` | Exact analytic, noise-free primary. Same geometry and same below-plane response as the Monte Carlo photons. |
| Scatter image | `postprocess.py` | Monte Carlo scatter (Compton, Rayleigh, fluorescence) in 1 mm bins, smoothed (σ = 5 mm), on the image grid. In the quick test the scatter showed no structure above the Monte Carlo noise at 0.02–0.5 cy/mm. With σ = 2 mm the Monte Carlo noise dominated the texture of the subtracted image; with 5 mm it is ≈ 6000× below it at 0.2 cy/mm. |
| Analysis | `analysis.py` | −ln(I/I₀), w*, residual under the Siemens weight, texture β and power (REPORT §11 method), iodine signal, SPR. Optional quantum noise. |

Why the image is split into an analytic primary and a Monte Carlo scatter:
- A Monte Carlo primary image at 85 µm pixels with mammography-level noise would need about 10¹¹ photons.
- Scatter is smooth, so ~10⁹ photons suffice for it.
- The two parts use the same geometry, materials and detector model. `selftest.py` shows they share the same normalisation: Monte Carlo/projector = 1.000 ± 0.001.

**Geometry** (world mm):
- a-Se surface at z = 0; focal spot at (0, 0, 656);
- x from the support's chest-wall face into the phantom, y along the chest wall;
- the phantom stands on the support at z = 17;
- image columns run along x from 3 mm (missed tissue [A]), rows along y.

**Phantom in GATE** (both models):
- `envelope` (PMMA) → `filling` (liquid) → `plate` (PMMA) → iodine discs (PMMA + iodine, one material per row).
- The spheres are daughters of `filling`, as six `Sphere` volumes with translation lists.
- L2 v2: 61 799 spheres (53 / 104 / 247 / 831 / 6 616 / 53 948 of 15.88 / 12.70 / 9.52 / 6.35 / 3.18 / 1.58 mm). They fill 61.0 % of the 1 096 ml free volume with no overlaps; 12 770 touch a wall (within 0.05 mm).
- Real L2: see "Phantom models".

**Meshes:**
- L2 v2: from the exact CadQuery solids, tessellated with an absolute tolerance of 0.005 mm (`phantom_meshes.py`). About 120 000 facets each; the facets deviate at most 0.007 mm from the exact surfaces (measured). The lid, gasket and positioning bracket are outside the field and not modelled.
- Real L2: prisms of a circular segment, written directly by `phantoms.py` (`python phantom_meshes.py --phantom L2`, no CadQuery). Inscribed polygons with a sagitta of at most 0.005 mm, about 650 facets each; closed and consistently oriented, volumes within 0.006 % of the exact solids.

**Spheres** (`sphere_packing.py`):
- They touch the walls, the plate and each other, as in the real phantom.
- The only clearances are numerical: 0.02 mm to walls (about 3× the mesh deviation) and 0.001 mm between spheres.
- Every sphere is verified against the actual mesh facets.

**Materials** (`data/L2v2_materials.db`):
- One file in GATE's material-database format holds every non-NIST material:
  - spheres: PA12 + 5.5 wt % CaCO₃ (1.046 g/cm³) and PMMA (1.19 g/cm³);
  - liquids: MCT oil (0.945 g/cm³) and liquid paraffin (0.86 g/cm³, assumed);
  - carbon fibre;
  - a-Se;
  - the seven iodine-disc rows (PMMA + iodine, 0.33–4.0 mg I/cm³).
- Each material is listed with its density, state, mass fractions and source.
- GATE loads it (`sim.volume_manager.add_material_database`). The analytic parts (projector, grid, detector) read the same file, so Geant4 and the analytic model use identical compositions.
- PMMA (container, plate), polycarbonate (paddle) and air are Geant4 NIST materials, used by name.
- Check the file after editing:

  ```bash
  python materials.py
  ```

  This lists densities, compositions and attenuation, and checks that fractions sum to 1 and that every material used is defined.
- Change the iodine rows with:

  ```bash
  python materials.py --iodine 1.2 0.9 0.7 0.4 0.27 0.17 0.10
  ```

  To give the 9 mm discs iodine, add `--large <mg/cm²>`.
- To add a third background, put its materials in the `.db` and add an entry to `BACKGROUNDS` in `config.py`.

## Phantom models

`--phantom` in every script (or `export REVELATION_PHANTOM=L2`); geometry in `phantoms.py` (real L2) and `sphere_packing.py`/`projector.py` (L2 v2); parameters in `config.PHANTOM_MODELS`. [U] = stated or measured by the user, [P] = Salomon et al. 2020/2023, [M] = measured in the images, [A] = assumed (user, 23 Oct 2026: "use the assumptions").

| Real L2 (`--phantom L2`) | Value | Source |
|---|---|---|
| Shape | half-cylinder, outer Ø 240 mm, axis in the flat chest-wall face | [P], [U]; centre of curvature on the front face in 0035 [M] |
| Outer height | 53.5 mm | [P], [U] |
| Lid / base / curved wall | 6 / 4 / 8 mm PMMA | [U] measured |
| Chest-wall (front) wall | 8 mm | [M] 8.2 mm band in 0035 |
| Cavity | R 112 mm, x ≥ 8 mm, z 4–47.5 mm (779 ml) | from the above |
| Plate | 2 mm PMMA, exactly in the middle of the cavity (centre z = 25.75 mm) | [U] |
| Plate outline | half-disc R 111.5 mm (0.5 mm inside the curved wall), from the front wall | [A] |
| Fill | 61 % PMMA spheres / 39 % liquid paraffin; 15.88 / 12.70 / 9.52 / 6.35 / 3.18 / 1.58 mm, equal volume per size | [P], [U]; sizes [A] |
| Densities | PMMA 1.19, paraffin 0.86 g/cm³ | [A] |
| Discs | 7 rows × 2/3/4/5/7 mm + two 9 mm, at the positions measured in the July images, back-projected to the plate plane | [M] `disc_layout_measured.json` |
| Iodine | 1.2 / 0.9 / 0.7 / 0.4 / 0.27 / 0.17 / 0.10 mg/cm² per row (`PMMA_I_L2_row*`, 2 mm plate); 9 mm discs plain PMMA | lower bounds [M]; 9 mm [A] |
| Position | front face at world x = 7.24 mm (column 56 of the July template, back-projected at M = 1.072), axis at y = 0 | [M] |

**Sphere packing of the real L2** (`python sphere_packing.py --phantom L2`, 23 min on the Mac; `phantom_L2/spheres.npz`, `spheres_summary.json`):
- 41 943 spheres: 36 / 71 / 167 / 564 / 4 490 / 36 615 of 15.88 / 12.70 / 9.52 / 6.35 / 3.18 / 1.58 mm.
- They fill 61.0 % of the 743.8 ml free volume (cavity 779.2 ml minus the plate) with no overlaps; 8 798 touch a wall or the plate (within 0.05 mm).
- Every sphere near a wall was checked against the mesh facets: none crosses one; the smallest clearance is 0.015 mm.

## Background options

| `--background` | Spheres | Liquid | Use |
|---|---|---|---|
| `PA12_MCT` (default) | PA12 + 5.5 wt % CaCO₃ (`PA12_CaCO3_5p5`) | MCT oil (`MCT_oil`) | the suggested new background: glandular-equivalent spheres, tissue-like residual |
| `PMMA_paraffin` | PMMA (`PMMA_sphere`) | liquid paraffin (`paraffin_oil`) | the current L2 materials: a reference that can be compared with the measured images |
| `gland_adipose` | glandular tissue (`glandular_tissue`, 1.04 g/cm³) | adipose tissue (`adipose_tissue`, 0.93 g/cm³) | tissue reference, not buildable: the REPORT §10 target pair (Hammerstein compositions as in the REPORT model), simulated so that `redesign_compare.py` judges the redesign against tissue in the same geometry, chain and method |

- **How to choose:** pass `--background` to `revelation_gate.py`, `projector.py`, `analysis.py`, `measlike.py`, `compare_measured.py` and `preview.py`, or set it once:

  ```bash
  export REVELATION_BACKGROUND=PMMA_paraffin
  ```

  In the full chain: `BACKGROUND=PMMA_paraffin ./run_all.sh`.
- **Outputs are kept apart:**
  - the GATE run tag contains the background (`LE_phantom_PMMA_paraffin_s1`);
  - `postprocess.py` refuses to mix runs of different backgrounds;
  - projector, analysis and preview outputs go to per-background folders.
- **Same geometry for both backgrounds:** within one phantom model, both backgrounds use the same container, plate and sphere positions. With `--phantom L2v2` the PMMA/paraffin option therefore isolates the effect of the materials; the real container is `--phantom L2`.
- **Plausibility check of the option** (stand-in spectrum, no scatter, so not results). The PMMA/paraffin option gives:
  - w* = 0.350, against 0.354–0.358 measured on the real L2;
  - a residual under Siemens' weight of 0.27, against ≈ 0.28 measured.

  PA12/MCT gives w* = 0.272 and a residual of 0.06, against 0.07 from the design calculation.

## Measurement-like images and the comparison with the measurements

`measlike.py` turns the simulated detected energy into what the Revelation stores in its raw (`TICEM_RAW`) images:
- detector blur: the presampling MTF published for the Revelation's a-Se detector [R1] (see "Detector model");
- flat-field normalisation, as the raw open field is flat;
- stored-value scale anchored to the measured HE open field of the pair, with LE at the same detector gain;
- quantum noise from the calibrated photon numbers (second moment of the detected energy from `projector.py`), times a noise-variance factor calibrated on the measured HE open field, and shaped by the same MTF;
- the 7500 ceiling.

### Detector model (23 Oct 2026)

`config.DETECTOR_RESPONSE`, applied in `measlike.py`.
- **Signal blur:** the presampling MTF of the Revelation's 85 µm a-Se detector measured by the NHS Breast Screening Programme [R1, Table 10]. It was measured with an edge on the breast support, in the u-direction (no tube motion): 0.91 / 0.86 / 0.78 / 0.67 / 0.55 / 0.44 at 1 / 2 / 3 / 4 / 5 / 6 cy/mm.
  - It includes the pixel aperture, which alone would give 0.99 / 0.95 / 0.90 / 0.82 / 0.73 / 0.62. The extra loss is consistent with K-fluorescence reabsorption and blocking-layer charge trapping in a-Se [R2, R3].
  - It is applied separably in x and y [A]. It was measured on another Revelation (detector LV2-00007, software VC10B, tomosynthesis projections without grid), and is assumed to hold for our unit (LV2-33695, 2D, grid) [A].
  - `mtf="aperture"` restores the previous ideal-aperture model.
- **Noise shape:** the quantum noise is generated before the detector blur and sampled at the pixel pitch. Its spectrum is therefore Σ_k MTF(f + k/a)², which is white for the ideal aperture.
  - With the published MTF it falls to 0.51 at Nyquist relative to 1 cy/mm.
  - In the measured HE open field it falls to 0.57.
- **Offset:** the stored LE values carry the offset measured in REPORT §13.4 from the exposure dependence of w*: 40 units at 28 kV (35 / 36 at 30 / 32 kV). There is no HE offset; §13.4 found the HE intercept of about −150 not to be a detector offset.
- **Noise level:** the variance of a stored pixel value is k × (photon noise from the photon numbers, MTF-shaped) + (f·S)², the fixed-pattern noise, white. Electronic noise is fixed at 0.
  - k and f are fitted on the measured HE open fields of all 36 July 28 kV exposures (12 each at low, normal and high mAs; levels 3700–6500 stored units; 0.6–5.9 cy/mm). The command is `python measlike.py --phantom L2 --scatter-le ... --scatter-he ... --calibrate-noise`, which writes `data/detector_noise_july28high.json`.
  - Result: **k = 1.78, f = 0.28 %**, 8.5 % rms misfit. In the open field 13 % of the noise is fixed pattern; behind the phantom about 4 %.
  - **Low-frequency structured noise** (below 1 cy/mm): the measured open-field NNPS minus the photon + fixed-pattern model, median of the 36 open fields. It is about 1.8·10⁻⁶ mm² at 0.2–0.3 cy/mm, steeply rising towards lower frequencies, and similar in all images.
    - It is applied as a multiplicative field shared by LE and HE (`structured_noise`, `structured_shared`) [A].
    - Sharing is supported by a test: shared, the residual at w* lies inside the measured range of the 12 July arrangements; independent, it is about twice too high.
  - A free fit with an electronic term gave e = 13 stored units. The measured LE phantom interior allows at most about 9, even if all of its high-frequency noise were electronic, and e = 13 tripled the simulated LE noise there. So electronic noise is fixed at 0 (`noise_fit_terms`); that variant is kept in `output/noise_model_with_electronic/`.
  - LE uses the same k and f [A], because the LE open field is at the 7500 ceiling.
- **Validation on HIGHMAS_28KV2_00** (pipe.py analysis; nothing fitted on the phantom):

  | Model | w* | LE interior | Noise 2.5–4.5 cy/mm, simulated / measured: LE / HE / H − W·L / H − w*·L |
  |---|---|---|---|
  | measured | 0.359 | 663 | — |
  | ideal aperture, white noise (before 23 Oct) | 0.326 | 622 | 1.29 / 1.12 / 1.13 / 1.11 |
  | published MTF, one factor 2.07 | 0.326 | 624 | 1.27 / 1.25 / 1.30 / 1.28 |
  | published MTF + LE offset + k 1.78, f 0.28 % | 0.347 | 664 | 1.04 / 1.14 / 1.19 / 1.18 |
  | **+ measured structured noise, shared (current)** | **0.351** | **664** | see below |

  With the structured noise, the residual at w* matches: power at 0.2 cy/mm 3.76·10⁻⁶ simulated against 3.70·10⁻⁶ measured. Discs with SDNR ≥ 3, simulated / measured: 4 / 3 in the recombined image R, 1 / 1 in H − W·L, 8 / 8 in H − w*·L, 9 / 8 with the split weight.

  - **The LE offset** closes most of the w* gap and matches the LE level.
  - **The high-frequency noise** behind the phantom is now within 4–19 % of the measured noise.
  - **At 1–2 cy/mm** the simulated LE and HE power is 13–21 % lower than measured. That band also contains texture, so the noise shape there is not separately validated.
  - **Tests of realistic parameters (23 Oct 2026), nothing fitted on the phantom:**
    - **HE level** (simulated ≈ 6–8 % low): no realistic change explains it with LE kept right. Paraffin density 0.827–0.890, H fraction 0.143–0.152 and PMMA 1.18 move HE by at most ±1 %. a-Se thickness and extra HE filtration move HE by at most +3 %, but LE by 11–24 % (through the LE/HE detector-gain relation). The anti-scatter grid (all parameters assumed) is untested: it needs a GATE run with the phase space kept. It is a limitation for now.
    - **Texture vs sphere arrangement:** the measured texture varies strongly between the 12 July arrangements, each image on its own. LE power at 0.2 cy/mm is 5.5–8.9·10⁻³ and LE β 1.81–2.69. The comparison pair HIGHMAS_28KV2_00 lies at the low end (5.9·10⁻³).
      - Two simulated arrangements (packing seeds 1 and 2) give 10.5 and 10.8·10⁻³, with β 2.54 and 2.37.
      - The simulated texture is therefore systematically ≈ 1.2× the highest and ≈ 1.4× the median measured value. It varies only ≈ 3 % between simulated arrangements, against up to ×1.6 between measured ones.
      - Not explained by realistic parameters; a limitation. Possibly the real sphere bed is not a homogeneous random packing (size segregation or clustering, voids).
      - β, the residual at w* and the disc SDNR are within the measured ranges.
    - **Texture contrast:** realistic densities change the texture power by at most ×0.88 (LE) / ×0.81 (HE). They do not explain the excess.
    - **Residual at w*:** a sphere-to-sphere density spread would need about 1.3 % (plausible). But its spectrum would be steep (β ≈ 3.5), while the measured residual is nearly white (β 0.2–0.9). The detector's low-frequency structured noise explains it instead (see above).
  - **Earlier results** are kept in `output/before_detector_model/` and `output/before_noise_model/`.

**References (detector model):**
- [R1] NHS Breast Screening Programme. *Technical evaluation of Siemens Revelation digital breast tomosynthesis system.* Public Health England / NCCPM, March 2019. Table 10: MTF of the central projections. https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/785503/Technical_evaluation_of_Siemens_Revelation_digital_breast_tomosynthesis_system_report.pdf
- [R2] Zhao W, Ji WG, Debrie A, Rowlands JA. Imaging performance of amorphous selenium based flat-panel detectors for digital mammography: characterization of a small area prototype detector. *Med Phys* 2003;30(2):254–263. Abstract: K-fluorescence lowers the MTF at Nyquist by ≈ 15 %; charge trapping at the blocking layer adds high-frequency blur in MTF and NPS.
- [R3] Zhao W, Ji WG, Rowlands JA. Effects of characteristic x rays on the noise power spectra and detective quantum efficiency of photoconductive x-ray detectors. *Med Phys* 28:2039–2049 (2001). Abstract: K-fluorescence reabsorption lowers the high-frequency NPS of a-Se by ≈ 12–18 %. Author order and year are as listed by a journal index; to confirm.
- Checked, without detector data: NHSBSP Practical Equipment Report *Practical evaluation of the Siemens Mammomat Revelation*, Public Health England, August 2020 (no MTF, NNPS or DQE).
- R2 and R3 were read at abstract level only. The model uses only the numbers of R1. The LE offset comes from this project's REPORT §13.4. The noise parameters come from the July measurements.

**One measured pair at a time, with its exact mAs.** Every script works on one measured pair, selected with `--pair` (or `REVELATION_PAIR`). The July default is HIGHMAS_28KV2_00, the image the disc template was made from.
- `measlike.py` makes that pair's simulated twin: the exact mAs of its LE image and of its HE image (`ExposureInuAs` in the DICOM headers), its HE open-field level, and a noise seed fixed for that pair.
- Nothing is averaged over pairs. Another pair is a separate run with its own `--pair`.
- The GATE runs do not depend on the mAs (they are normalised per emitted photon), so the mAs enters only here, in the spectra fit (that pair's entrance dose) and in `analysis.py --noise`.

**Measured reference sets** (`--set`, `config.MEASURED_SETS`, read by `measured_set.py`):

| `--set` | Images | Phantom | Unit | Default for |
|---|---|---|---|---|
| `july28high` | 12 pairs, 28 kV high mAs, 21–30 Jul 2026 (R-CC, mirrored like the July analysis); one is used per run (`--pair`) | real L2 with the new plate | 10159 | `--phantom L2` |
| `0035` | 1 pair, 18 Aug 2025 | real L2 container with the old plate | 1024 | `--phantom L2v2` |

The July exposure values are not typed in. `measured_set.py --set july28high` reads them once from the headers and images of the 12 pairs and writes `data/measured_july28high.json`.
- It stores per pair (no set averages) the LE and HE mAs, the entrance doses and the HE open-field level.
- It compares the header geometry (SID, source–entrance distance, pixel size, kVp, filters) with `config.py`.

`config.set_measured()` takes the values of the selected pair. `spectra.py --set july28high` then fits the extra Al filtration of unit 10159 to that pair's LE entrance dose (0.136 mm Al for HIGHMAS_28KV2_00, the same as for unit 1024). The spectra go to `spectra_july28high/`.

`compare_measured.py` analyses the selected measured pair and its simulated twin with identical code:
- **Masks:**
  - HE below 0.64 × open field;
  - interior ≥ 15 mm from the edge, iodine discs excluded (array box + 12.75 mm; 9 mm discs + 8 mm);
  - core ≥ 30 mm from the edge.
- **Disc locations:**
  - July: the July template (`july2026/template.json`), shifted to the pair;
  - 0035: a 7 × 5, 16 mm lattice fit;
  - simulated: the model's layout.
- **Metrics** (values of that pair and of its twin):
  - open-field and interior levels, HE transmission;
  - noise power at 2.5–4.5 cy/mm;
  - texture β and power at 0.2 cy/mm;
  - w* at ≥ 30 mm from the edge, and r at the Siemens weight;
  - the scatter "bowl" of the recombined image near the edge;
  - the recombined-image signal of each iodine disc (July: simulated/measured ratio).

Outputs: `output/measlike/<case>/<set>/<pair>/` and `output/compare/<set>_<case>/<pair>/` (`metrics.json`, `comparison.png`, `spectra_profiles.png`). `<case>` is the background for L2 v2 (e.g. `PMMA_paraffin`) and `L2_<background>` for the real L2.

Real L2 against a July pair, after the GATE runs (`run_all.sh` does all of it with `PHANTOM=L2`; add `--pair <id>` for a pair other than HIGHMAS_28KV2_00):

```bash
python measlike.py --phantom L2 --scatter-le LE_L2_PMMA_paraffin --scatter-he HE_L2_PMMA_paraffin
```

```bash
python compare_measured.py --phantom L2
```

The July images are found in the project folder (the session folders `CEM_NEWL2_MASAR_HIGHMAS_28KV_02` and `…28KV2_02` are there) or on the SanDisk (paths in `july2026/exposures.json`); elsewhere, add `--root "<folder 'New L2 Measurements'>"`. `compare_0035.py` still works and is the same as `compare_measured.py --set 0035`; it needs `../TICEM_RAW_L-CC_DIAGNOSIS_0035`.

**Checked on the measured side (0035):**
- The script reproduces the report: w* 0.353 (report 0.351–0.357), interior LE/HE 683/1475 (682/1492), HE transmission 0.285 (0.288–0.290).
- `measured_set.py --set 0035` reads exactly the mAs (198.39 / 38.76) and entrance doses (7.513 / 0.722 mGy) that are in `config.py`.

**What the comparison can and cannot show:**
- **Real L2 vs July:** the same container, plate and materials on both sides. What still differs is one sphere arrangement against twelve, the assumed sphere sizes, plate outline and front wall (see "Phantom models"), and the iodine contents, which are lower bounds.
- **L2 v2 vs 0035:** different phantoms (v2 container against the old L2 with the old plate). Only interior background quantities (levels, noise, texture, w*) are comparable.
- The measured spectra keep falling above 1 cy/mm. Since 23 Oct 2026 the simulation includes the published detector MTF, the measured LE offset and a noise model fitted on the measured open fields (see "Detector model"). Behind the phantom the high-frequency noise is now within 4–19 % of the measurement.

## The investigation's own analysis on the simulation (July pipeline)

`pipeline_compare.py` runs the analysis and recombination used for the measured July images, **unchanged**, on a measured pair and its simulated twin. It uses `investigation_0035_0037/scripts/july2026/pipe.py`, `agg.py` and `discagg.py` (REPORT §13).
- The twin is written as DICOM files that look like the measured input (raw orientation, uint16 stored values; LE, HE and a recombined image R), so the same files read both.
- One pair per run (`--pair`, default HIGHMAS_28KV2_00); each image is analysed on its own and nothing is averaged.
- Three workspaces, each with links to those scripts and its own `exposures.json`:

| Workspace | LE, HE | Recombined image R | What it is for |
|---|---|---|---|
| `measured` | the measured pair | Siemens' own | the reference |
| `measured_emuR` | the measured pair | the investigation's emulation of Siemens' recombination (REPORT §12), with the offset and mask of the real R | shows what the emulation itself cannot reproduce |
| `simulated` | its twin from `measlike.py` (that pair's exact LE / HE mAs and HE open level) | the same emulation, with the offset and mask of the measured pair | the simulation |

What `pipe.py` computes for each image, and the report compares (values of the one pair, side by side):
- weights vs spatial frequency: w*(f), coherence, and Siemens' effective weight from regressing R on H and L;
- **rim:** local w* in 128-px tiles vs distance from the curved edge;
- texture spectra, β and P(0.2) of LE, HE, H − W·L, H − w*·L and R;
- contrast and SDNR of every disc in LE, HE, H − W·L, H − w*·L, the split weight and R.

The report also adds the **rim profile**: ring medians of H − W·L, H − w*·L and R vs distance from the fitted curved edge and from the chest wall, relative to the interior. This is the ring-median method of `scripts/shading.py` (REPORT Fig. 6), with the 0035 pixel coordinates replaced by `pipe.py`'s geometry fit and disc registration. It is applied identically to all three workspaces.

**What the measured rim already shows** (pair HIGHMAS_28KV2_00):
- In the raw H − W·L there is a bright band 5–10 mm inside the curved edge (+0.065 above the interior), fading by ≈ 30 mm. Being in the raw image, it comes from the physics (scatter, edge geometry), so the simulation can reproduce it or not.
- In Siemens' R the band is weaker than in the emulation (+270 against +500 stored units at 5–10 mm).
- In the outermost ≈ 3 mm, Siemens keeps R near the interior level, while the emulation drops by ≈ 1000. This is Siemens' border processing, which the emulation does not contain. In that zone the simulation can only be compared with `measured_emuR`, not with Siemens' image.

```bash
python pipeline_compare.py all --phantom L2
```

- **Steps:**
  - `prepare`: workspaces and the DICOM twin;
  - `run`: `pipe.py` in each workspace, then `agg.py` and `discagg.py`, which with one pair just list that exposure;
  - `report`.

  `run_all.sh` runs it at the end for the July set. For another pair, add `--pair <id>` to `measlike.py` and to this command.
- **Runtime:** ≈ 3.5 min per pair on the Mac (≈ 50 s `prepare`, ≈ 40 s `pipe.py` per workspace, ≈ 1 min report).
- **Needs:** pydicom; the July session folders in the project folder (as now) or on the SanDisk.
- **Outputs:** `output/pipeline/july28high_L2_PMMA_paraffin/`:
  - `measured/`, `measured_emuR/`, `simulated/`: per workspace `results/<id>.json`, `agg.json` / `agg.txt`, `discagg.json` / `discagg.txt`;
  - `report/<pair>/`: `summary.txt`, `1_weights_vs_frequency.png`, `2_rim_and_edges.png`, `3_texture_spectra.png`, `4_iodine_discs.png`, `5_images.png`, `report.json`.
- **Disc registration of the twins:** `pipe.py`'s estimate (outline + summed H − 0.37·L response) was 6 px off on the PMMA twin and 14 px on the PA12 twin. The twins are now analysed with the registration known from the simulation (projector iodine map; `run_twins`, `pipe.py` unchanged except its `register()`). The measured pair keeps `pipe.py`'s estimate (an independent disc-shaped locator agrees within 1 px). `twins.json` records both.
- **One change to `pipe.py`:** it loaded the emulation parameters from a temporary folder of an earlier session. It now falls back to `scripts/siemens_cem_emulation_params.json`, which is byte-for-byte the same file, so it runs on the workstation too. The computation is unchanged.
- **Tested** on HIGHMAS_28KV2_00 with a stand-in twin (rough spectra, no scatter). All steps run, and the twin goes through every analysis. Without scatter the twin shows almost none of the measured rim band, which is what the real scatter run has to test. The numbers from that test are not results.

## Predicting the redesigned background from the current L2 (`redesign_compare.py`)

The simulation reproduces the current L2 inside the phantom only up to known errors (texture power 1.2–1.8×, HE level −7 %). `redesign_compare.py` therefore reads another background, by default the REPORT §10.4 design (PA12 + 5.5 wt % CaCO₃ spheres in MCT oil), as a **change from the current phantom**, never in absolute terms:
- positive quantities (weights, levels, texture and noise power, disc contrast and noise): predicted = measured current × (simulated redesign / simulated current);
- signed quantities (β, rim heights): predicted = measured current + (simulated redesign − simulated current);
- r = 1 − W/w* from the predicted w*, with the pair's header W (Siemens' weight does not depend on the object);
- disc SDNR: the measured SDNRs scaled, per image, by the median contrast ratio (4–7 mm discs, rows 0–2) over the median noise ratio (all 35 discs).

The two twins share container, plate, sphere positions, spectra, detector model, mAs and noise seed, so errors common to both cancel. Errors that depend on the materials do not [A]. The script checks that the twins differ only in the materials and give the same disc registration, and prints and records any mismatch.

Targets in the table (REPORT §10, model [E]): r = r_tissue (+0.07 at 28 kV, glandular in adipose tissue under Siemens' W, i.e. w* ≈ 0.274); LE texture power 0.94 × the current L2; β unchanged (same spheres). Patient β values are listed as indicative only (other method). Not covered: glandular equivalence (63 % in the model) and AEC (the twins use the measured pair's mAs).

```bash
python redesign_compare.py --phantom L2 --redesign PA12_MCT
```

- **Needs:** `pipeline_compare.py` results (`results/<pair>.json` and `report/<pair>/report.json`) of both cases for the same pair: `output/pipeline/july28high_L2_PMMA_paraffin/` and `output/pipeline/july28high_L2_PA12_MCT/`. `run_all.sh` calls it at the end of `PHANTOM=L2 BACKGROUND=PA12_MCT ./run_all.sh`; if the PMMA results there come from an older detector model, the check reports it. In that case copy `output/pipeline/july28high_L2_PA12_MCT/` to the Mac and run it there.
- **Outputs:** `output/redesign/july28high_<pair>_PMMA_paraffin_vs_PA12_MCT/`: `summary.md` (measured, both simulations, simulated change, predicted redesign, target), `prediction.json`, `prediction.png` (w*(f), LE and R texture spectra, rim of H − W·L and R, disc SDNR by size).
- **Tested** (24 Oct 2026, on the Mac): with the current background as a stand-in redesign, every predicted value equals the measured one and the material check fires; with a synthetic redesign of known changes (w* × 0.78, LE level × 1.1, LE texture × 0.9, β + 0.1, R disc contrast × 1.05 and noise × 0.8, rim + 0.01), every prediction follows the rules above. A β is not transferred when the image's texture power changes by more than 2× (its spectrum then changes composition); the simulated value is shown instead.
- **Tissue reference:** `--reference gland_adipose` (default). Once `PHANTOM=L2 BACKGROUND=gland_adipose` has been simulated and run through `measlike.py` and `pipeline_compare.py`, the reference is predicted from the measured current L2 in the same way and added as a column, as summary lines (redesign vs reference) and in the figure. Until then the REPORT §10 model values are the only targets. Tested with a stand-in reference (every reference value equal to the redesign, as it must be).
- **First result** (PA12_MCT, 24 Oct 2026): see `PROJECT_NOTES.md` §6. In the L2, 13 of the 14 tiles ≥ 45 mm from the curved edge lie in the strip next to the chest wall, so that local-w* bin measures the PMMA front wall against the liquid, not the interior.

## Install and run

```bash
pip install opengate==10.1.1 numpy scipy matplotlib
```

opengate brings Geant4 (data downloaded on first use), SpekPy and uproot. CadQuery is only needed to regenerate the meshes.

```bash
python selftest.py
```

```bash
python revelation_gate.py --beam LE --n 1e6 --tag speedtest
```

```bash
PHANTOM=L2 ./run_all.sh
```

**Real L2 against the July measurements, step by step:**
1. `data/measured_july28high.json` is already written: per pair, the exact mAs of the LE and HE image, the doses and the HE open level (`python measured_set.py --set july28high` rewrites it from the images).
2. On the workstation the run needs:
   - this folder (without `output/`);
   - `../L2_container_v2/disc_layout_measured.json`;
   - `../investigation_0035_0037/scripts/` (with `july2026/` and the emulation parameters);
   - for the comparison, the July session folders `CEM_NEWL2_MASAR_HIGHMAS_28KV_02` and `…28KV2_02` in the main folder (as on the Mac).
3. `PHANTOM=L2 ./run_all.sh` (one pair, HIGHMAS_28KV2_00; another pair with `PAIR=<id>`).

`run_all.sh` runs one case (phantom model + background + measured set) in order:
1. the exposure values of the measured set, if not read yet (`measured_set.py`; July: needs the images);
2. spectra for that X-ray unit (fitted to the selected pair's entrance dose);
3. self-test;
4. projector;
5. `RUNS` LE and `RUNS` HE GATE runs, each post-processed and its phase space deleted;
6. analysis, without and with noise;
7. the measurement-like twin of the selected pair (its exact mAs), and its comparison with that pair;
8. for the July set: the investigation's own analysis on the pair and its twin (`pipeline_compare.py`);
9. for another background in the real L2: its measurement predicted from the current L2 (`redesign_compare.py`).

Variables:
- `PHANTOM` (`L2v2` default, or `L2`), `BACKGROUND`, `SET` (defaults per phantom), `PAIR` (the measured pair, default per set), `ROOT` (July images folder if neither in the project folder nor on the SanDisk);
- `THREADS`, `RUNS` (default 5);
- `N_LE` (default 2·10⁷ per run), `N_HE` (4·10⁶).

Outputs go to `output/`, kept apart per case (`<case>` = background for L2 v2, `L2_<background>` for the real L2):
- `projector/<case>/`;
- `mc/<beam>_<case>.npz`, `scatter_…`;
- `analysis/<case>[_noise]/summary.json`, `images.png` (top row at the LE contrast width, bottom row stretched per image; a stretched subtracted image magnifies ~1 % residuals), `texture_spectra.png`;
- `measlike/<case>/<set>/<pair>/`, `compare/<set>_<case>/`;
- `preview/<case>/`.

## Look at the model before running

There are two options. The first works now, without GATE.

**1. `preview.py` (no GATE needed).**

```bash
python preview.py --phantom L2
```

It reads the same meshes, sphere positions, disc layout and constants that the GATE script uses, and writes to `output/preview/<case>/` (without `--phantom L2`: the L2 v2; `--background` changes only the labels):

| File | What it shows |
|---|---|
| `setup_side.png` | The whole set-up to scale (focal spot, field, aperture, paddle, phantom, support, detector), and a zoom of the phantom and the detector stack at y = 0, including the grid and covers applied after GATE |
| `phantom_sections.png` | The phantom cut at y = 0 and at x = 60 mm, a sphere layer seen from above, and the plate plane with the iodine discs coloured by row |
| `projection.png` | What the detector sees geometrically: path length through spheres and PMMA, and iodine per pixel. No spectrum needed. |
| `phantom_3d.html` | Interactive 3D phantom in a browser: container, plate, discs, spheres (large ones as spheres, small ones as markers). Parts can be hidden by clicking the legend. |
| `setup_3d.html` | Interactive 3D set-up: focal spot, field edges, emission cone, aperture, paddle, phantom, support, phase-space plane, detector |

Options:
- `--no-3d`, `--no-projection`;
- `--small-spheres none` (lighter 3D file);
- `--bin` (projection pixel binning, default 4).

**2. GATE's own visualisation (needs opengate).** This shows the geometry as Geant4 actually builds it, with a few photon tracks:

```bash
python revelation_gate.py --phantom L2 --visu vrml_file_only
```

This writes `output/LE_phantom_L2_PMMA_paraffin_visu_large/geometry.wrl` (L2 v2: `LE_phantom_PA12_MCT_visu_large`); open it in any VRML viewer, e.g. ParaView.

To look at that Geant4 scene without a VRML viewer, convert it into an interactive 3D page and a picture. This runs anywhere with plotly and matplotlib, including the Mac:

```bash
python view_vrml.py output/LE_phantom_L2_PMMA_paraffin_visu_large/geometry.wrl
```

It writes `geant4_scene.html` and `geant4_scene.png` next to the `.wrl`. The workstation's opengate is the `novis` build without Qt (checked 22 Oct 2026), so this is the Geant4 view available there.

Other `--visu` values:
- `qt`: the interactive Geant4 viewer (Qt/OpenGL), opened after the photons are tracked.
  - It needs a display: the workstation's own desktop, or `ssh -Y` from a machine with an X server. OpenGL over SSH may be slow.
  - The script stops with a message if the opengate build has no Qt or no display is found; GATE itself would silently show nothing.
  - The view is set in `data/visu_qt.mac`: world hidden, no facet edges, z up, phantom centred, photon tracks in yellow. Edit it freely.
  - In the window: drag to rotate, scroll to zoom; Geant4 commands can be typed at the bottom.
- `vrml`: a window via pyvista (`pip install pyvista`);
- `gdml_file_only`: the geometry as GDML.

Defaults with `--visu`:
- single thread, 200 photons (`--visu-n`);
- only the spheres ≥ 6.35 mm (`--visu-spheres large`, about 1 200), because all 60 000 make the scene very slow. `--visu-spheres all` or `none` change this.

A `--visu` run is for looking only; do not post-process it.

**Runtime and disk** (measured on the workstation, 8 threads, opengate 10.1.1, 22 Oct 2026):
- **Speed:** LE 36 000 photons/s, HE 12 000 photons/s. HE photons cross the phantom and many sphere boundaries, so they need about twice the steps. Add ≈ 1 min of set-up per run (placing 61 799 spheres).
- **Phase space:** ≈ 3 MB (LE) and ≈ 19 MB (HE) per 10⁶ photons; each file is deleted after post-processing.
- **`run_all.sh` defaults:** 5 runs × 2·10⁷ LE and 5 × 4·10⁶ HE photons, giving 1·10⁸ LE and 2·10⁷ HE in total. That takes ≈ 50 + 35 min per background, with scatter noise ≈ 1.5 % (LE) and < 1 % (HE) after smoothing.
- **Photon count:** `--n` is the total number of photons. GATE 10.1.1 counts `source.n` per thread, so the script divides by the number of threads.
- `sphere_packing.py` took ≈ 45 min here (≈ 20 min to grow the spheres, ≈ 25 min to remove the last overlaps at full size). It saves a resume point at every step (`--resume`).
- `projector.py` at full resolution: 97 s and 3.7 GB memory here. `preview.py`: ≈ 10 s.

## First run: what to check

1. **`run_info.json` → `sanity.ok` is true.** Every unscattered photon at the plane must trace back along its direction to the focal spot (within 1 mm) and arrive in the field (x ≥ 0). Scattered photons are not tested, since they may go anywhere. With fewer than 5 unscattered photons there is no verdict (`ok: null`).
   - The direction convention was verified in the Geant4 source (`G4SPSAngDistribution::GenerateIsotropicFlux`).
   - A failure would mean the installed version differs.
2. **Spectra and entrance dose:** `spectra.py` prints the modelled entrance air kerma against the DICOM value (7.51 mGy LE, 0.72 mGy HE), for the tube alone and with the air and paddle in front of the breast.
   - **First result on the workstation (tube only):** LE ×1.64, HE ×1.03.
   - **Interpretation:** an LE-only excess points to missing low-energy filtration (paddle, air, mirror or window). That makes the LE spectrum too soft, which affects w*, not only the noise level.
   - **Calibration (default since 22 Oct 2026, user decision, no QC report available):** `config.SYSTEM["extra_filtration_al_mm"] = "fit"`. `python spectra.py`, and therefore `run_all.sh`, chooses the extra Al-equivalent tube filtration that brings the LE entrance dose (air and paddle included) to the DICOM value.
   - **HE check:** HE is the independent check. `spectra.py` warns if HE then differs from DICOM by more than 10 %, because extra filtration alone would not explain both beams.
   - **Fixing a value:** `--extra-al <mm>` uses a given value instead.
   - **Best constraint:** measured HVLs from the system's physics QC report would pin the spectra down better than the DICOM dose, which the system calculates itself.
   - **Traceability:** the extra filtration is stored in `spectra/spectra_info.json` and in every run's `run_info.json`.
3. **Projector vs Monte Carlo primary:** keep one phase space (no `--delete-phsp`) and compare `mc/mcprimary_*.npy` with `projector/primary_*.npy` on 1 mm bins. They should agree within statistics everywhere inside the phantom. This tests the tessellated walls and sphere placement in GATE against the analytic geometry.
4. **Overlap check:** one run with `--check-overlaps`. This is slow with 60 000 spheres; the packing is already verified against the meshes.

## Assumptions

[V] = from the DICOM headers of this unit; [M] = measured or derived in this project; [A] = assumed, to replace when known. All are in `config.py` and `grid_detector.py`.

| Item | Value | Source |
|---|---|---|
| SID, source–entrance, support–detector | 656, 585, 17 mm | [V]; support–detector confirmed by the magnification factor 1.0266 = 656/639 |
| Pixel, matrix | 85 µm, 3518 × 2800 | [V] |
| LE / HE beams | 28 kV Rh 50 µm / 49 kV Ti 1 mm; mAs = exact header values of the selected measured pair (0035: 198.39 / 38.76 mAs; HIGHMAS_28KV2_00: 200.64 / 41.39 mAs; others in `data/measured_july28high.json`) | [V] |
| Siemens weight `_W` | 0.255 at 28 kV (0.265 / 0.292 / 0.308 at 30 / 32 / 34 kV) | [V] |
| Focal spot | 0.3 mm, uniform disc | [V] nominal size, [A] shape |
| Anode angle, Be window | 16°, 0.5 mm | [A] |
| Extra tube filtration | Al-equivalent, fitted so that the LE entrance dose matches DICOM (air + 2 mm paddle included; DICOM dose assumed to be below the paddle) | [M] calibrated to the DICOM dose, user decision 22 Oct 2026; value in `spectra/spectra_info.json` |
| Heel effect | not modelled (isotropic emission) | [A]; cancels in −ln(I/I₀) except for small spectral changes |
| Collimated field | x 0–245 mm, y ±155 mm at the detector | [A]; [M] no collimator edge in the images: the field covers at least the whole stored image |
| Compression paddle | 2 mm polycarbonate | [A] thickness and material; [V] paddle ID 14, force 49–74 N |
| Breast support | 1.5 mm carbon fibre (CFRP approximation) | [A] thickness and material; [V] table ID 1, 17 mm above the detector |
| Missed tissue at the chest wall | 3 mm | [A]; measure |
| Grid | focused linear, 31 lines/cm, ratio 5, 16 µm Pb, carbon-fibre interspace, 2 × 0.5 mm covers, lamellae along x | [V] focused and linear; all dimensions [A]. [M] No residual grid lines in the images (moving grid), so line density and orientation are not measurable there |
| Detector | direct a-Se (ANRAD), contact mode, 3584 × 2816 elements of 85 µm; image 3518 × 2800 placed at FOV origin (33, 0) | [V] |
| Detector layers | a-Se 200 µm, 1 mm carbon-fibre cover, energy-integrating, Se K-escape neglected | [A] |
| Container (L2 v2) | rev B: walls 6 / 4 mm (top / bottom), 1.5 mm rim, plate 3 mm at mid-height | [M]; the real plate is 2 mm [U] |
| Container (real L2) | half-cylinder Ø 240 × 53.5 mm, walls 6 / 4 / 8 mm, front wall 8 mm, 2 mm plate in the middle of the cavity | see "Phantom models" |
| Spheres | PA12 + 5.5 wt % CaCO₃, 1.046 g/cm³ (no voids); 15.88 / 12.70 / 9.52 / 6.35 / 3.18 / 1.58 mm, equal volume per size, 61 % of the free volume | [M] design values (MATERIALS.md, `data/L2v2_materials.db`); sizes [A]; porosity / oil uptake not modelled |
| Liquid | MCT oil, C₂₉H₅₄O₆ approximation, 0.945 g/cm³ | [M] design value |
| Current-L2 option | PMMA spheres 1.19 g/cm³ (NIST PMMA); liquid paraffin H 0.148 / C 0.852, 0.86 g/cm³ | [M] materials of the current L2; PMMA sphere density = NIST value [A]; paraffin composition and density [A] (Ph. Eur. range 0.827–0.890 g/cm³, REPORT §6) |
| Iodine | 1.2 / 0.9 / 0.7 / 0.4 / 0.27 / 0.17 / 0.10 mg/cm² per row; the two 9 mm discs without iodine | [A] lower-bound estimates (REPORT §13.2) |
| Detector response (measurement-like images) | flat-fielded stored values ∝ deposited energy; HE open field of the selected measured pair (0035: 5158; HIGHMAS_28KV2_00: 5784), LE at the same gain; presampling MTF of the Revelation detector, separable; photon noise shaped by the aliased MTF × k 1.78, fixed-pattern noise 0.28 %, no electronic noise, low-frequency structured noise measured in the open field, shared by LE and HE (all fitted or measured on 36 July HE open fields; same at LE); LE offset 40 units at 28 kV, none at HE; scatter has the primary's mean detected energy | MTF [R1] (other unit) [A]; k, f, offset, HE level [M]; LE noise parameters and the rest [A] |
| Physics | G4EmLivermorePhysics, fluorescence on, 0.1 mm cut. Electrons killed at creation (they cannot reach the image; bremsstrahlung neglected). | model choice |

## Limitations

- **Simplifications:**
  - The sphere arrangement is one random packing, with no gravity settling. In the real phantom the spheres sink (PA12 1.046 vs oil 0.945).
  - At 61 % packing the bed is nearly jammed, so the difference should be small. A layer of lower density may form near the top wall.
  - Spheres are rearranged between real acquisitions. Rerun `sphere_packing.py --seed N` for a new arrangement.
- **Source targeting:** with `--target phantom` (default), photons that miss the phantom are not simulated. This drops paddle and support scatter from the open field, which only reaches the outer ~15 mm of the phantom image. The analysis region is further inside. One `--target field` run quantifies it.
- **What the model does not include:**
  - Siemens' own processing (LE blur K(f), scaling, offsets; REPORT §12);
  - detector blur (MTF) and electronic noise;
  - ghosting and lag.

  `analysis.py` uses the plain weighted log subtraction. Apply the Siemens emulation from the investigation scripts on top if needed.
- **Material check:** PMMA compatibility of MCT is untested (MATERIALS.md). The model assumes the design materials. To model the low-cost nylon/glycerol-PG option, add those materials to `data/L2v2_materials.db` and select them in `PHANTOM` (`config.py`), ideally with measured densities.

## Open questions (to confirm with the team or Siemens)

1. Grid specifications: lines/cm, ratio, lamella and interspace material, position. Not in the image headers and not visible in the images.
2. a-Se thickness, detector cover, support and paddle thickness and material. Not in the image headers. Paddle ID 14 and table ID 1 are in the headers, so Siemens' paddle/table list would identify them.
3. Missed tissue at the chest wall (measure with a ruler phantom).
4. Nominal iodine concentrations of the new plate.
5. Anode angle and inherent filtration (used by SpekPy).
6. Actual sphere sizes and the measured density after oil saturation, once the spheres are made.
7. Real L2: the six sphere diameters, the chest-wall wall thickness (8 mm from the images), the plate outline and how it is held, the nominal iodine contents of the July plate and what the two 9 mm discs contain.

## Files

| File | Purpose |
|---|---|
| `config.py` | all parameters with their source |
| `data/L2v2_materials.db` | all non-NIST materials, GATE format; read by GATE and by the analytic parts |
| `materials.py` | reads, checks and lists the `.db`, writes the iodine rows, loads the `.db` into GATE; attenuation coefficients from NIST XCOM (`data/xcom_elements.json`) |
| `spectra.py` | SpekPy spectra, kerma check |
| `phantoms.py` | phantom models: real L2 geometry (free space, ray chords, disc placement, meshes) and the common interface |
| `phantom_meshes.py` | STL meshes of envelope, filling and plate (L2 v2: CadQuery; real L2: `phantoms.py`) |
| `sphere_packing.py` | sphere positions (`phantom/` or `phantom_L2/spheres.npz`), with verification |
| `revelation_gate.py` | GATE 10 simulation |
| `grid_detector.py` | grid, covers, a-Se response |
| `projector.py` | analytic primary images |
| `postprocess.py` | phase space → scatter and Monte Carlo primary images |
| `analysis.py` | dual-energy analysis and figures |
| `selftest.py` | checks without GATE |
| `preview.py` | figures and interactive 3D views of the model, without GATE |
| `measlike.py` | measurement-like stored-value images (aperture, flat field, quantum noise, ceiling) |
| `measured_set.py` | measured reference sets: DICOM headers (exact mAs, doses, geometry), images, HE open level |
| `compare_measured.py` | pair-by-pair comparison of simulated and measured raw images (July set or 0035) |
| `compare_0035.py` | same as `compare_measured.py --set 0035` (kept for the earlier commands) |
| `pipeline_compare.py` | the investigation's July analysis (`pipe.py`, `agg.py`, `discagg.py`) unchanged on one measured pair, the same pair with emulated R, and its simulated twin; rim profiles; side-by-side report |
| `redesign_compare.py` | another background in the real L2 (default PA12_MCT) predicted from the measured current L2 and the two twins (ratio or difference of the simulations), against the REPORT §10 targets |
| `view_vrml.py` | interactive 3D page and picture of the scene Geant4 exported with `--visu vrml_file_only` |
| `run_all.sh` | full chain |

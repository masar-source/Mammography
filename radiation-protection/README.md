# Radiation protection: a GATE 10 room model for the Siemens MAMMOMAT Revelation

Room shielding and staff protection for a mammography room, simulated before installation: the GATE 10 model of the Revelation, patient-derived compressed breasts, walls, door and console shield, and air kerma scored at staff and neighbouring positions.

This package adds six segmented breasts from Radboud University Medical Center and two material controls to the GATE 10 room study. All eight options support CC, right MLO and left MLO for room-scatter comparisons. It includes the six original compressed NRRDs, six registered GATE label images, the room simulation and an offline 3D viewer with independent scanner and console placement.

## Run a saved room design

1. Open **Revelation_Patient_CC_3D.html** and change the room, scanner, console, shield, phantom, projection and acquisition to the setup you want to simulate.
2. Choose **File → Save geometry settings (.json)**. Each exported file represents one setup. You can keep it in Downloads or save it inside the project.
3. Open Terminal in the extracted project folder. Type `bash run_cases.sh` followed by a space, drag the saved JSON from Finder into Terminal and press Return. The example below assumes you named the file `my_room.json` in Downloads.

~~~bash
cd Mammography/radiation-protection   # the folder of this README
bash run_cases.sh "$HOME/Downloads/my_room.json"
~~~

The launcher reads the file contents, so renaming the file is fine. It prints the input path, selected layout, room dimensions, breast-to-operator distance, shield, phantom, exposure, door state and output directory. It starts one batch of 10,000 histories for the saved exposure. The positional-file command takes its physical settings directly from the file. It ignores preset-related environment variables. With no input file or explicit CONFIG it prints usage and exits.

4. Open the printed results directory. It is under `results/<export name>/<file hash>/n<histories>_b<batches>/`. `projection_survey.csv` contains the station results. `run_manifest.json` identifies the input and its exact archived copy. Editing and resaving a design under the same name creates a different results folder.
5. For additional setups, repeat the export. To process several in sequence, drag all selected JSONs after `bash run_cases.sh `. You can also save only your own exports in a folder and run:

~~~bash
bash run_cases.sh "$HOME/Downloads/my_room_exports/"*.json
~~~

For a DBT export, the launcher transports the projection selected in the viewer. For a complete 25-projection sweep, use the explicit command below with your own export. The 14 bundled examples are available as starting geometries in the viewer. Their current viewer exports are in `cases/`, named `case_01_...` through `case_14_...` with the selected geometry and exposure in each filename. For the student experiment, pass the JSON saved after making your edits.

The named weekly programme uses the angles and mAs in its workload table. For a specific exported custom angle, select the `viewer` programme. The three MLO/CC choices refer to arm tilt. Floor rotation of the installed scanner and the workstation heading are separate settings. The breast reference is fixed for a given layout. The manifest records its horizontal distance to the operator. The projection CSV also gives its three-dimensional distance to each cell centre. A breast exposure always includes the upright torso. Selecting `open` as the phantom removes both the breast and torso. The paddle, support, detector and room still contribute to this background run.

### Compare shield widths, door states and exposures

Use your existing GATE Python. On a Mac with the virtual environment in Downloads:

~~~bash
PY="$HOME/Downloads/.venv-gate/bin/python"
CASE="$HOME/Downloads/my_room.json"
"$PY" study.py run --config "$CASE" --beam HE --view MLO_L --shield-width-mm 1000 --histories 10000 --batches 1 --out results/07_HE_left_w1000
"$PY" study.py run --config "$CASE" --beam HE --view MLO_L --no-shield --histories 10000 --batches 1 --out results/07_HE_left_no_shield
"$PY" study.py run --config "$CASE" --door open --histories 10000 --batches 1 --out results/07_LE_door_open
"$PY" study.py run --config "$CASE" --phantom open --histories 10000 --batches 1 --out results/07_LE_empty
~~~

The 1000 mm width is an example study variable. Changing it preserves both personnel positions and the centre and heading of the shield. Use a distinct output folder for each comparison. Staff cells cover the operator and second person (named `trainee`), plus left and outer reach, at 500, 1200 and 1700 mm above the floor. Compare the largest kerma across both people and their reach positions for each design. Treat a width as a candidate only after checking uncertainty and the planning assumptions. For an open-door comparison use `door_leaf` at the centreline and `door_jamb` near the far jamb. These are fixed air cells in the corridor. The leaf and its lead layer turn 90° outward about the left jamb. In the small region where the viewer's open leaf intersects the wall at its hinge, GATE assigns the wall material once.

### Run a complete oblique DBT sweep

~~~bash
"$PY" study.py run --config "$CASE" --programme dbt_sweep --histories 10000 --batches 1 --jobs 1 --out results/07_DBT_right
~~~

This command uses the saved base projection and custom tilt. Use `--programme viewer` on a DBT export to transport just the selected position. In the viewer, select DBT and move the projection slider to turn the visible tube housing and swivel arm through the selected positions. The focal spot and source-to-phantom line are hidden by default. Select **Show focal spot and beam line** in the display settings to show both together in LE, HE or DBT. Switching back to LE or HE returns the arm to the base view. Each sweep uses 25 independent source positions, with the detector and compressed breast fixed at the base arm angle. The nominal central source-to-detector distance is the existing system model's DICOM-derived 656 mm. The rotation centre is 30 mm above the support, or 47 mm above the native detector reference. The grid is withdrawn. The model samples 25 equally spaced angles from +25° to −25° in viewer order and assigns equal mAs shares. The collimation projects the native rectangular field onto the stationary detector plane. The individual pulse coordinates, static rectangular field, and constant-output spectrum are research assumptions. The UK technical evaluation gives the sweep range, projection count, pivot, grid state and equal division of tube loading. It also describes an initial 5 mAs AEC exposure. The reported sweep excludes that preliminary exposure. [UK technical evaluation, section 2.1 and Table 1](https://assets.publishing.service.gov.uk/media/60dafd13d3bf7f7c2c6ba253/Technical_evaluation_of_Siemens_Revelation_digital_breast_tomosynthesis_system_report.pdf) · [Siemens Revelation operator manual, tomosynthesis option](https://doclib.siemens-healthineers.com/rest/v1/view?document-id=878411)

By default DBT uses the bundled 28 kVp W/Rh LE spectrum and fluence as an explicit proxy. To supply a DBT-specific beam, append `--dbt-spectrum /path/to/DBT.npz`. The archive must contain one-dimensional `edges_keV`, non-negative photon-bin `weights` and scalar `fluence_per_mAs_1m` (photons/cm²/mAs at 1 m). To calculate µGy per sweep, also pass `--scan-mas` with the total mAs of the 25 pulses. Geometry alone does not establish clinical DBT output. MLO uses the native CC-compressed anatomy rotated with the receptor as a scatter surrogate.

### Read the results

| File | Interpretation |
|---|---|
| `survey_geometry.json` | Active air cells and excluded supplementary stations, with room coordinates and reasons. The 500 mm common ring and 1000 mm supplementary ring have eight scanner-relative azimuths each. An excluded position has no estimate. |
| `projection_survey.csv` | Field collision air kerma and secondary-field kerma in µGy/mAs per exposure, station and independent batch. It includes beam, base angle, DBT offset, coordinates, spectrum model and photon steps. The secondary field includes all modeled scattering structures. |
| `exposure_survey.csv` | Single-exposure or sweep means by station and standard errors across independent batches. For a full sweep the program averages the 25 estimates using mAs fractions inside each batch. One batch gives no standard error. Zero sampled tracks cannot establish zero exposure. |
| `run_manifest.json` | Geometry, horizontal breast-to-operator distance, effective exposure weights, source inputs and archived configuration path. Individual `case.json` files include the focal spot in native and room coordinates. |
| `weekly_survey.csv`, `layout_comparison.csv` | Produced for explicitly named weekly workload programmes. Single viewer exposures and DBT sweeps give field-only research estimates. A clinical room decision also needs the workload, housing leakage and planning basis. |

Use `study.py preview --config "$CASE" --out results/07_preview` with the same Python interpreter to inspect the floorplan and the list of active stations. Before transport the runner checks complete cell volumes against the room, shield and other survey cells, and checks equipment and patient bounds for the selected view on corrected layouts. An invalid required station stops the run with its name. The supplementary ring omits stations that intersect the wall or the shield. `--check-overlaps` additionally enables the Geant4 geometry check. Since field fluence can be low behind shielding, use smoke runs to verify geometry and output; increase histories and independent batches before comparing designs.

## Cases

Glandular volume fraction is the number of fibroglandular voxels divided by the sum of adipose and fibroglandular voxels. All source voxels have equal volume. The fraction excludes skin and air. Thickness in the table is the occupied voxel extent along the compression axis. Total tissue volume includes skin.

| GATE phantom | Source ID | Glandular volume | Nominal thickness | Voxel thickness | Total tissue |
|---|---|---:|---:|---:|---:|
| `rad_fatty` | KBCT410079R | 7.45% | 52 mm | 52.14 mm | 817 cm³ |
| `rad_intermediate` | KBCT110136R | 20.04% | 53 mm | 52.96 mm | 504 cm³ |
| `rad_dense` | KBCT110288R | 40.79% | 54 mm | 54.05 mm | 479 cm³ |
| `rad_highest` | KBCT310033L | 58.90% | 53 mm | 52.96 mm | 331 cm³ |
| `rad_cohort_min` | KBCT310116R | 6.32% | 63 mm | 63.06 mm | 1,235 cm³ |
| `rad_cohort_max` | KBCT010152L | 82.24% | 40 mm | 40.13 mm | 229 cm³ |

We counted the four labels in all 821 compressed volumes. The first four examples come from the 79 breasts with nominal thickness 52–54 mm, close to the 53.5 mm L2 phantom. Within that group we selected the nearest fractions to 0%, 20%, 40% and 100%. The last two examples are the minimum and maximum across all 821. These are quantitative research selections. BI-RADS composition is assigned from mammographic appearance. The release provides tissue labels and does not supply individual BI-RADS assessments. Even in the thickness-matched series, the breasts have different projected areas, tissue volumes and spatial distributions.

### Same-shape material controls

| GATE phantom | Parent anatomy | Material assignment |
|---|---|---|
| `rad_fat_only` | `rad_fatty` | Both interior labels use adipose tissue. Interior glandular fraction is exactly 0%. |
| `rad_uniform` | `rad_intermediate` | Both interior labels use one uniform adipose–glandular mixture with 20.04% glandular tissue by volume. |

Both controls retain their parent skin, external shape, compression thickness and placement. Comparing `rad_intermediate` with `rad_uniform` isolates the redistribution of the same tissue inventory. Comparing `rad_fatty` with `rad_fat_only` isolates the replacement of its glandular tissue with fat.

## Run on the existing Mac installation

Extract the ZIP and open the project folder in Terminal. The scanner modules, spectra, materials and L2 geometry are in `Revelation_Scanner_Inputs/` inside the project. The launcher uses a project-local `.venv-gate` if present and also recognises your existing `Downloads/.venv-gate` when the project is in Downloads. For a GATE installation elsewhere, set `GATE_PYTHON` to its Python executable.

Use the saved-file command at the top of this guide. To change the computational budget for that same design:

~~~bash
HISTORIES=1000000 BATCHES=8 JOBS=1 bash run_cases.sh "$HOME/Downloads/my_room.json"
~~~

Choose histories and batches using the uncertainty at the protected positions. The launcher also accepts `OUTDIR` to change the results root.

## Oblique exposures for the scatter study

The release contains 821 patient breasts mechanically compressed in CC. For each selected case, CC is the native compression at 0°. In the two MLO simulations the same voxel volume rotates with the X-ray tube, paddle, breast support and detector. The preset angles are +45° for `MLO_R` and −45° for `MLO_L` in this project's coordinates. Both values come from `study.json → views_deg`. They are study choices. The room, shield, operator and survey points keep their coordinates. The optional torso block stays upright. [Radboud dataset](https://doi.org/10.5281/zenodo.18617629)

Compare each phantom with itself at all three angles, using the same spectrum, voxel map, breast thickness and room layout. The resulting change estimates the effect of beam and detector orientation on the room scatter field. The method preserves the same breast for both sides and applies no mirroring or recompression. It omits the change in chest-wall tissue coverage and pectoral anatomy associated with clinical MLO positioning. That limits absolute clinical estimates. A CEM run transports the LE or HE spectrum through the supplied unenhanced tissue labels. Housing leakage, if assessed, uses a separate source term. [IAEA-TECDOC-1447, §7.3.11](https://www-pub.iaea.org/MTCD/Publications/PDF/TE_1447_web.pdf)

`run_manifest.json` records the patient, tissue mapping, placement and the actual angle and anatomy assumption for each view. Each `case.json` and each row in `projection_survey.csv` also identifies the CC-derived oblique surrogate.

## Inspect the anatomy and room

Open `Revelation_Patient_CC_3D.html` in a WebGL 2 browser. Select a **Breast phantom** and choose **CC**, **R-MLO** or **L-MLO**. You can switch breasts while viewing either oblique angle. Below the selector, the note gives the anatomy and projection assumption. **Phantom** zooms to the anatomy. **Phantom and grid internals** reveals the fibroglandular surface or the assigned control interior. **Room** returns to the whole room. Drag the desk, monitor or shield to move the console and staff. **GATE torso surrogate** displays the block used by the transport model. The file exports the selected scene to GLB and the room settings to JSON.

### Gantry and detector geometry

The black detector table uses the overall image-receptor dimensions given by Siemens for the Revelation VC20: **359.2 mm wide × 346.5 mm deep × 45.0 mm high**. Width is parallel to the chest-wall edge. Depth runs from the patient towards the upright. The moulded corners approximate the product photographs. The manual gives a **maximum** of 5 mm for chest-wall dead space with the object table included. The viewer uses that maximum to position the patient-facing edge at native x = −5 mm. The surface supporting the breast is at native z = 17 mm in the supplied scanner model. [Siemens MAMMOMAT Revelation Operator Manual, VC20, Print No. XPW7-340G.620.60.02.24, pp. 172–173](https://doclib.siemens-healthineers.com/rest/v1/view?document-id=878404#page=172)

The exterior view shows a continuous black tabletop. Directly beneath it, the original GATE breast support is a 330 × 265 × 1.5 mm carbon-fibre layer. Its top is at z = 16.99 mm, the original numerical clearance below the breast. The active selenium area remains 304.64 mm wide × 239.36 mm deep. Switch off **Equipment exterior and console** to see the transport volumes. Switch off **Phantom and grid internals** as well for an opaque view of the support.

The internal backing and beam stopper use a **modelled** footprint of **330 mm wide × 295 mm deep**. The width matches the native breast support. The retained 295 mm depth is a study input. Their patient-facing edge is 3 mm inside the housing, at x = −2 mm, with the rear edge at x = 293 mm. The project assigns 3.0 mm aluminium and 0.5 mm lead to these layers. Siemens identifies the detector and object table as the beam stopper. Its published external dimensions do not establish the proprietary layer sizes, thicknesses or materials. The value for an operator shield is a separate specification. [Siemens Operator Manual, pp. 28, 172–173](https://doclib.siemens-healthineers.com/rest/v1/view?document-id=878404)

Siemens drawing 18072, sheet A-101, details 1–2, gives a 779.5 × 704.9 mm floor base and 1,135.1 mm from the stand orientation mark to the rear datum. The reconstruction registers the mark to the detector chest-wall plane, placing the front of the base 430.2 mm behind that plane. A white support carries the detector back to the rotating arm. The room dimensions and scanner placement determine the room-wall clearance.

The GATE focal spot is at (0, 0, 656) mm and the detector reference plane is z = 0 mm. The source, collimation, breast, paddle and receptor share the same rigid transform in CC and both MLO views. The external cover is used for visualisation and mechanical contact checks. The existing GATE support supplies the physical material at the breast surface.

**Earlier configurations and results.** Earlier viewer files used 400 mm wide internal plates. To run a saved room with the corrected plates, open its JSON in this viewer and save it again. This preserves the saved scanner, console and shield placement. An unmodified older GATE JSON retains its original plate settings for reproduction. The console pilot report in `references/` and the original protocol describe earlier work. Use fresh transport results with the revised backing when comparing shielding or dose. Each new `case.json` records the detector parameters actually used.

### PDFs included with the project

The `references/` folder contains the earlier project protocol (`Mammography_Room_Radiation_Protection_Project.pdf`) and the L2 console-shielding pilot report (`Revelation_L2_Console_Pilot_Review.pdf`). The Siemens site-planning drawing (MAMMOMAT Revelation, drawing 18072, sheet A-101) is a Siemens document and is not included here; the scanner footprint above was taken from it. The research plan for the current student project is not included here.

### Set the internal room dimensions

Under **Room dimensions**, enter **Width (X)**, **Depth (Y)** and **Height (Z)** in metres. Press Enter or leave the field to apply the change. Width and depth move opposite walls about the floor centre. Height moves the ceiling above the finished floor. The equipment keeps its size and room coordinates. **Reset size** restores the 6.40 × 5.20 × 3.00 m preset while preserving equipment edits.

The fixed-size doorway stays on the south wall. Its centre retains the same proportional position across the available wall width. The service opening retains its size and proportional location on the east wall in presets that include it. Wall layers, slabs and door shielding retain their specified thicknesses. Resizing changes the boundary survey positions in GATE. A room size applies to every shielding preset in the saved study.

Choose **File → Save geometry settings (.json)** to run the design in the included GATE model, or **Export current 3D model (.glb)** to inspect it in another 3D application. JSON stores the internal dimensions under `room.inside_mm` in millimetres. Each GATE case records `room_inside_mm` and `world_size_mm`.

### Resize the shield and compare layouts from a fixed workstation

Use **Shield size → Width** or **Height** to resize the shield in millimetres. The front panel grows equally to either side of its horizontal centre. Height is measured from the floor. A height edit scales the window sill in proportion. Open **Window and side return** to set the sill height or the return length explicitly. Edited dimensions carry across shielded layouts. **Use preset size** restores the selected layout’s shield dimensions at the current position. Lead equivalence comes from the selected preset.

Changing layouts preserves the position and heading of the console and shield, including any separation you set with the independent drag controls. **Keep scanner fixed** also preserves the scanner position and rotation for every change. The default **Follow distance presets** moves only the scanner when you select a distance-labelled preset whose nominal distance differs from the previous layout. It travels along the existing breast-reference-to-operator line and retains its floor rotation. Equal-distance and shielding-only changes preserve the scanner placement. **Apply** sets the displayed preset distance explicitly, even when the scanner was moved manually. In follow mode, editing **Breast to operator** also moves the scanner. In fixed mode, editing that field explicitly moves the workstation. Choosing the open reference hides the shield and retains its placement and custom size for the next shielded comparison.

### Place the scanner

Select **Plan** for an overhead view. Under **Equipment placement → Drag in 3D**, choose **Move scanner** and drag the gantry or pedestal. Choose **Rotate scanner** and drag left or right to turn it. The X and Y offsets and floor rotation fields accept exact values. Positive rotation is counterclockwise when viewed from above. **Reset scanner** returns the scanner to its preset while preserving the workstation position. **Reset layout** restores both.

Under **Room angles**, 0° follows room +Y, parallel to the east and west walls. Positive angles turn counterclockwise in Plan view. At +90° the selected edge points along −X, parallel to the north and south walls. The three angles control the detector chest-wall edge, front edge of the console worktop and width of the front shield panel. The scanner turns about the breast reference, the console about the operator floor point, and the shield about the centre of its front face. Every angle edit preserves the horizontal **Breast to operator** distance. The X and Y fields locate the pedestal after rotation. CC and MLO tilt the arm independently.

**Room dimensions** displays simultaneous contacts. Amber marks touching and red marks overlaps or extension beyond a named boundary. Equipment checks include the scanner exterior, console, shield panels, frame, feet and side return at the current arm angle. The breast and torso are also checked against the console, shield and room limits. Personnel messages identify scoring centres at their GATE heights. The mechanical casing and workstation use the viewer's reconstructed geometry, with a 0.01 mm tolerance for contact. The colored notices describe placement. GATE calculates the radiation field.

Saved geometry stores the absolute floor angles in `scanner_pose.yaw_deg`, `console_heading_deg` and `screen_heading_deg`. The two explicit heading fields take precedence over the legacy combined orientation. Staff and reach positions follow the console tangent. Panel and shield-edge scoring positions follow the shield tangent. The GATE `azimuth_deg` field locates the operator from the breast using +X as its native zero. The viewer converts that value for the **Operator bearing** field. Earlier saved layouts open at their original angles and positions.

Use **File → Save geometry settings (.json)** to keep a layout. **Export current 3D model (.glb)** saves its visible geometry at the same placement. Use the copy of `study.py` included in this release to run a saved scanner placement. The selected layout’s `scanner_pose` expresses the installation as yaw about the preset plinth point plus the stored translation, including the compensation that holds the breast reference fixed during an angle edit. The one-metre ring is centred and oriented with the scanner. The wall and door survey points remain attached to the room.

For a saved room edit, pass the exported JSON with `--config`. Optional `--layouts` and `--phantom` arguments override the saved selections. The 3D viewer uses reduced surfaces from a 1.09 mm grid for the outer contour and 1.37 mm for glandular detail. The colour slice overview in `Radboud_CC_Anatomy.png` shows full-resolution central slices. The GATE simulation always uses the full label image.

## Registration, materials and verification

The original NRRDs are in `anatomy/original/`. The files in `anatomy/gate/` retain all tissue voxels and the native 0.273 mm isotropic spacing. We removed exterior air padding and permuted the voxel axes with a right-handed rotation. Source axes map to scanner coordinates as `(X, Y, Z) = (+i, −k, +j)`. Minimum native i is the chest-wall side and native j is the compression axis. The registered chest-wall face lies at scanner X = 0, the bounding box is centred laterally, and the inferior face lies at the native phantom base, 17.00 mm above the detector reference. The paddle follows the actual voxel thickness with the scanner’s 0.01 mm gap. The room’s reference point for comparing console distances stays fixed across breasts.

| Source label | Tissue | GATE density | Mass fractions H, C, N, O, P |
|---|---|---:|---|
| 0 | Air | Native `G4_AIR` | Geant4 air composition |
| 1 | Adipose | 0.93 g/cm³ | 0.112, 0.619, 0.017, 0.251, 0.001 |
| 2 | Fibroglandular | 1.04 g/cm³ | 0.102, 0.184, 0.032, 0.677, 0.005 |
| 3 | Skin | 1.09 g/cm³ | 0.098, 0.178, 0.050, 0.667, 0.007 |

The uniform control uses volume-weighted density and mass-weighted elemental fractions. Each original and registered image has a SHA-256 digest in `anatomy/cases.json`, alongside its counts, dimensions, orientation and crop indices. `anatomy/cohort_density.csv` gives the ranked inventory of all 821 cases. `anatomy/reorient.py` can rebuild a registered image from its original.

All six registered arrays were compared voxel by voxel with their source arrays after reorientation and cropping. The CC transport checks with the revised backing use 160,000 histories across eight breast options and two beams. Their results are recorded in `transport_check.csv`.

The three-view verification used OpenGATE 10.1.0 on Linux. Six segmented breasts and two material controls each completed LE and HE at 0°, +45° and −45°, giving 48 runs of 10,000 histories each. All 1,488 per-view and per-position entries agreed with the transport output. With the higher-density patient, the exported 3D poses matched GATE for the breast and paddle at all three angles and for the upright torso and fixed workstation. These short runs verify transport and reporting. Evaluation of shielded-room kerma requires independent batches and adequate statistical precision.

For the viewer-linked runner, the geometry build reached the transport boundary for all 14 saved placements at CC and both ±45° MLO angles using the intermediate breast. All layouts kept the same operator floor coordinates. The eight common ring cells fit in all 14 cases. Fifteen of the 112 supplementary stations met the wall or shield. Actual OpenGATE 10.1.0 smoke runs included LE and HE, an empty control, both oblique 25-projection DBT sweeps with the door open and shield removed, and a DBT case with the Geant4 overlap check enabled. The two sweep tests used 1,600 histories per position. All 50 source positions registered forward photons. The combined estimates at all 39 active stations in each sweep were independently recalculated from the per-projection files.

## Data attribution and references

- Pacheco G, Michielsen K, Sechopoulos I. *Patient-derived compressed digital breast phantom dataset for mammography and digital breast tomosynthesis simulations*. Version 1.0, 2026. Original NRRDs redistributed under CC BY 4.0. The registrations, material controls, visual surfaces, selection tables and overview are adaptations prepared for this project. [Dataset](https://doi.org/10.5281/zenodo.18617629) · [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
- Pacheco et al. *Optimization and validation of a mechanical compression model for digital breast phantoms in mammography and tomosynthesis simulations*. See Figure 2 and the description of chest-wall fixation and compression. [Medical Physics](https://doi.org/10.1002/mp.17980)
- Ivanov et al. *Suitability of low density materials for 3D printing of physical breast phantoms*. Phys. Med. Biol. 63 (2018) 175020. The material compositions and densities follow Table 1, which attributes the tissue data to Hammerstein et al. [Article](https://doi.org/10.1088/1361-6560/aad315)
- American College of Radiology. *BI-RADS Contrast Enhanced Mammography*. [Breast composition categories](https://cs.acr.org/-/media/ACR/Files/RADS/BI-RADS/BIRADS_CEM_2022.pdf)
- OpenGATE 10.1.0. [Image volume and voxel-to-material mapping](https://opengate-python.readthedocs.io/en/10.1.0/user_guide/user_guide_reference_volumes.html#image-volume)

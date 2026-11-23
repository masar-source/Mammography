# L2 phantom container v2, rev B: design notes (concept, 21 October 2026)

**Status:** concept for discussion, not yet checked against the physical L2 container or the Revelation's chest-wall side. The open points are in §8.

**Revision history:**
- **22 October 2026:** walls as measured on the current L2 (top 6 mm, bottom 4 mm). Rim wall 1.5 mm, skin-equivalent. Extension end walls 10 mm.
- **Rev A:** rounded rim, lid and plate supports at the chest wall.
- **Rev B:** the user's proposal of 21 October 2026. A straight extension beyond the chest-wall edge carries the lid, the bolting zone and the plate supports. The lid is plain. A new, longer plate is used.
- **Rev B + bracket:** a removable positioning bracket whose lip hooks against the support's chest-wall face (user's choice, 21 October 2026).
- **27 October 2026, slot leak fixed:** the plate slots reached beyond the imaged part's outer surface at x = 0 and opened through the step face there (found by the user in the STL). Each slot end is now closed by an outer lug (user's choice instead of a shortened slot, which would leave a fragile 2 mm web). Only the container STL has been re-exported so far; STEP, renders and drawings still show the old geometry.
- **27 October 2026, plate 2 mm:** the new plate is 2 mm PMMA like the July plate (user: "it should be 2 mm even for v2"; was 3 mm assumed). Slots 2.4 mm high, lugs 8.4 mm high; screws, ports, lid and gasket unchanged. Container and plate STLs re-exported (3 mm versions in `stl/superseded/`); STEP, renders and drawings still show the 3 mm plate and no lugs.
- **27 October 2026, equal walls:** top and bottom walls 5 mm each (user; was 6/4 as on the current L2), same 10 mm total. Cavity and plate now at the outer mid-height (z = 26.75 mm); rim symmetric (5 / 1.5 / 5 mm). Screws, ports and gasket moved with it. Container, lid, gasket and plate STLs re-exported and verified (6/4 versions in `stl/superseded/`); STEP and renders not updated.
- **20 November 2026, drawing sheet redrawn:** `drawings/L2v2_drawing.pdf/.png` and `L2v2_edge_profile.png` now show the current geometry (equal walls, 2 mm plate, slot-end lugs, filling volume from `checks.json`). Every label has a leader ending in an arrow head on its feature; `l2v2_drawing.py` checks the sheet for text that overlaps text, lines, leaders or discs and reports 0 collisions.

**Coordinates:**
- **x = 0** is the chest-wall face of the breast support, where the inner face of the positioning lip sits. The image starts a few mm inside, at x = d_miss (missed tissue, to be measured; §8).
- **x > 0** runs into the phantom. **x < 0** holds the extension, lid and bracket, all outside the field.
- **z = 0** is the support surface.

## 1. What changes and why

**Rounded rim (since rev A).**
- The curved side has a constant-radius rim, like the periphery of a compressed breast. The current L2's vertical wall projects as ≈ 53.5 mm of PMMA over its outer ~5 mm.
- Top and bottom are seamless.

**Extension beyond the chest-wall edge (rev B).**
- The container continues **20 mm beyond the support's chest-wall face** as a straight section with the same cross-section.
- The spheres and oil continue through it. The image therefore shows "tissue" right up to the chest-wall edge, with no container wall, as with a real breast.
- In the extension, outside the field, sit:
  - the bolting zone;
  - the plain lid with its gasket and screws;
  - the plate slots.
- **The only fixed structure left in the image is the rounded rim.**

**New, longer plate (rev B).**
- A half-disc plus a full-width tab that reaches into the extension.
- The disc layout is the measured July layout, measured from the container's chest-wall face.

**Positioning bracket.**
- A removable L-shaped PMMA bracket is screwed under the extension. Its 16 mm-deep lip hangs in front of the support. The phantom is pushed towards the chest wall until the lip touches the support's chest-wall face; this sets x = 0 reproducibly.
- **Why removable:**
  - the phantom can still stand flat on a bench;
  - its position can be fine-tuned by ±2 mm (slotted holes);
  - the bracket can be re-made if the Revelation's support face is shaped differently.
- It lies entirely at x ≤ 0, so it is never in the image.

**Expected effect on the edge [E, geometry only]** (`drawings/L2v2_edge_profile.png`):
- The filling thickness falls smoothly to zero over the outer 27 mm.
- The PMMA path through the rim peaks at ≈ 18 mm (1.5 mm wall), about 1.5 mm inside the edge. An idealised breast with 1.45 mm skin gives ≈ 17 mm of skin there (§4a).
- Not predicted: whether Siemens' recombined image then shows a skin-line rim like the clinical ones.

## 2. Geometry (defaults; all parametric)

| Parameter | Value | Note |
|---|---|---|
| Outer height H | 53.5 mm | outer (confirmed 21 Oct 2026) |
| Top / bottom wall | **5 mm / 5 mm** | equal walls (decided 27 Oct 2026); same 10 mm total as the current L2 (lid 6 / base 4 measured) |
| Imaged part | half-disc, outer R 126.75 mm, flat region R 100 mm | footprint 253.5 mm (decided 21 Oct 2026) |
| Rim | outer radius 26.75 mm (= H/2), constant; inner semi-ellipse 25.25 × 21.75 mm centred at z = 26.75 mm (concentric) | wall 5 mm (top and bottom) → **1.5 mm** at the equator (skin-equivalent, §4a), smooth |
| Extension | 20 mm beyond the support face, 273.5 mm wide | end walls thickened to 11.6 mm at the equator (outside the field) |
| Bolting zone | last 6 mm of the extension | top wall +2 mm (7 mm), bottom wall +3 mm (8 mm); 18 × M3, blind tapped 5 mm |
| Lid | plain 6 mm cast PMMA | 4 × M8 fill/vent ports, 2 per half; 1 mm FKM gasket |
| Plate (new) | half-disc R 124.85 mm + tab 20.8 × 255.9 mm, **2 mm thick** (as the July plate) | centred in the cavity at z = 26.75 mm (= the outer mid-height, with equal walls) |
| Plate slots | 2.4 mm high (0.2 mm clearance per side), 3.0 mm deep into the end walls | from the open end to x = 0 |
| Slot lugs | x 0–5 mm, \|y\| 125.55–131.25 mm, z 22.55–30.95 mm | close the slot ends on the outside of the step face at x = 0 (5 mm beyond the slot end, 3 mm around it); at the chest-wall corners, just outside the rim |
| Positioning bracket | lip 170 × 5 × 16 mm, horizontal leg 20 × 170 × 4 mm | 4 × M3 from below into the 7 mm bottom wall of the bolting zone, slotted ±2 mm |
| Overall | depth 153.75 mm (lid to rim) | overhangs the support face by 27 mm |
| Filling volume | ≈ 1128 ml (2 mm plate, 5/5 walls; 1096 ml in the original rev B) | ≈ 192 ml of it in the extension |

## 3. Plate at mid-height

**Three supports:**
1. **Tab in slots, outside the field.** The tab runs in the two slots in the extension's end walls, which hold it to ±0.2 mm.
2. **Rim seat.** The curved edge sits in the narrowing rim (cavity half-height 3.9 mm at the plate edge). It self-centres there, with at most ≈ ±2.9 mm vertical play at the far edge (2 mm plate; ±2.4 mm with 3 mm).
3. **Plain lid.** The tab end protrudes 0.8 mm into the gasket gap, so the plain lid presses the plate home.

**Behaviour:**
- **The two halves are sealed from each other.**
- **Insertion checked numerically:** the plate slides in with ≥ 0.2 mm clearance.
- **Sphere filling:** below and above the plate the opening is 17.75 / 18.75 mm high (2 mm plate), so 16 mm spheres pass.

## 4. Imaging-relevant consequences

- **No fixed structure inside the field** except the rounded rim.
- **Uniform-thickness region:** R 100 mm.
  - The 7 mm discs end at R 93.0 mm, 7.0 mm inside it. The filling is only 0.04 mm thinner 1 mm into the rim and 0.3 mm thinner 3 mm in.
  - Keep ROIs inside R ≈ 95 mm.
- **Disc positions relative to the image edge.** In the July images the container's chest-wall face sat ≈ 5 mm inside the image (an air gap is visible at the image edge). With the bracket, the phantom touches the support face. The discs therefore appear about 5 mm + d_miss closer to the image's chest-wall edge than in July. The 2 mm column will be ≈ 16–19 mm from the image edge for d_miss = 2–5 mm, against ≈ 26 mm in July.
  - To reproduce the July image positions instead, shift the layout away from the chest wall (fixed on the plate) or use the bracket adjustment. Both reduce the margin to the rim.
- **Oil-rich periphery:** only small spheres fit near the rim, so the outer ~15 mm is oil-rich.

### 4a. Wall thickness and skin (added 22 October 2026)

**Measured on the current L2:** lid 6 mm, base 4 mm, side wall 8 mm.

**Is the lid/base skin-equivalent?** The L2 papers do not state a reason for these thicknesses.
- Expressed as skin (ICRU 46), 6 mm PMMA ≈ 4.9 mm and 4 mm PMMA ≈ 3.3 mm of skin at the LE mean energy.
- Measured breast skin is 1.45 ± 0.30 mm ([Huang et al. 2008](https://pmc.ncbi.nlm.nih.gov/articles/PMC2811552/)).
- Older dose models assumed 4–5 mm "skin". That may be where the lid/base came from **[?]**.

In the flat region the lid and base add only uniform attenuation. v2 first kept them unchanged (6/4 mm); on 27 Oct 2026 they were made equal (5/5 mm), keeping the same 10 mm total, so the transmission in the flat region is unchanged while the cavity, plate and rim become symmetric about the outer mid-height.

**The rim is where the wall shows up in the image** (as the skin line). Its PMMA thickness was therefore matched to real skin (`investigation_0035_0037/scripts/skin_equiv.py`, mono estimates **[E]**):

| | LE | HE | Recombined (H − W·L) |
|---|---|---|---|
| PMMA thickness equivalent to 1.45 mm skin | 1.76 mm | 1.50 mm | 1.21 mm |
| Rim peak with a **1.5 mm** wall, × skin | ×0.84 | — | ×1.22 |
| with a 2.0 mm wall | ×0.96 | — | ×1.40 |
| with a 5 mm wall (= lid/base average) | ×1.48 | — | ×2.15 |
| with the current 8 mm wall | ×1.81 | — | ×2.63 |

**1.5 mm is the best compromise across LE and the recombined image.** PMMA cannot match skin in all three at once, because its energy dependence differs (μ_H/μ_L 0.48 vs 0.41).

## 5. Materials and manufacturing

**Container: cast PMMA, CNC-machined as two halves split at mid-height.**
- The inner rim is an undercut that can only be machined from the split plane.
- The slots straddle the split plane, so each half carries a 1.7 mm recess.
- Bond the halves with a polymerising PMMA adhesive, then anneal.

**Prototype alternative: a one-piece SLA or PolyJet print.** Resin attenuation, leak-tightness and oil compatibility must be tested.

**Other parts:**
- **Lid:** plain 6 mm cast PMMA, with M8 tapped ports, nylon plugs and FKM O-rings.
- **Gasket:** 1 mm FKM (Viton); check it for MCT if the oil changes.
- **Screws:** M3 × 12 pan head, nylon (PA66) or PEEK, for the lid; M3 × 8 for the bracket.
- **Bracket:** PMMA (or POM), CNC-cut, with 4 slotted holes. It is outside the beam, so its material does not matter for imaging. If it will be removed often, threaded inserts in the bottom wall are more durable than tapped PMMA; metal is acceptable there (outside the field). Check with the workshop whether the 7 mm wall can take them.
- **Plate:** a new 2 mm PMMA plate, as the July plate (decided 27 Oct 2026). Make the iodine discs as for the July plate, at the positions in `disc_layout.csv`.

**If the filling changes to MCT oil** (report §10), an environmental-stress-cracking test of the bonded PMMA container is mandatory.

## 6. Assembly and filling (proposal)

1. Remove the bracket. Stand the container lid-side up in a cradle, with the rim resting in it. The cradle is not designed yet.
2. Slide the plate in with the tab edges in the slots until the curved edge sits in the rim.
3. Split the sphere set into two equal portions by weight. Pour them alternately into both halves, tapping.
4. Add oil slowly to each half. Tap until no more bubbles rise.
5. Fit the gasket and the plain lid. Tighten the screws crosswise, hand-tight.
6. Top up each half through its ports until oil overflows. Close the plugs.
7. Refit the bracket. Lay the phantom on the support and push it towards the chest wall until the lip rests against the support face.

To rearrange the spheres, rotate or shake the phantom as now. Stand it lid-up and vent before imaging.

## 7. Checks done (`checks.json`, `step/brep_checks.json`)

| Check | Result |
|---|---|
| Slot ends sealed (local flood fill, `slot_ends_sealed()`; STL ray test) | yes since 27 Oct 2026 (before: open at both ends) |
| Exact solids (STEP) | 5 valid single solids (container, lid, gasket, plate, bracket). 0 mismatches against the implicit model in ≈ 27,000 test points. |
| Plate insertion, minimum clearance | 0.2 mm |
| Wall around the lid taps / above the bracket taps | ≥ 2.25 mm / 3.0 mm |
| Gasket web next to the end screws | ≥ 2.02 mm |
| Rim wall: bottom / equator / top / thinnest | 5.0 / 1.50 / 5.0 / 1.50 mm |
| Bracket taps to the nearest lid tap | 12.5 mm |
| Outermost disc edge / flat-region margin | R 93.0 mm / 7.0 mm |
| Bolting zone to the support face | 14 mm |
| Projected width of the imaged part | 273 mm (field of view 291 mm) |
| Filling volume (voxel check, implicit model) | 1128.1 ml (2 mm plate, 5/5 walls; original rev B: 1095.9 STEP / 1097.6 voxel) |

## 8. Open points (to confirm with the team)

1. ~~Outer or inner height; lid and base~~ **53.5 mm outer; lid 6 mm, base 4 mm, side wall 8 mm (measured 22 Oct 2026).** Adopted for v2 (top 6, bottom 4); **changed 27 Oct 2026 to equal 5 / 5 mm** (same total).
2. The July plate has the same outline and thickness as the old plate. Please send the thickness and material so the new plate matches.
   **Conflict found 27 Oct 2026:** the user has since stated "2 mm PMMA plate" (used in the GATE model, `gate10_revelation/list.md` §4), but v2 is still sized for a 3 mm plate (`t_plate` 3.0: slots 3.4 mm high, plate STL 3 mm). Decide the new plate's thickness before printing the container. **Decided 27 Oct 2026: 2 mm.**
3. Wall thickness of the current container on the chest-wall side: **to be checked by the team.**
4. ~~Footprint~~ **Decided: 253.5 mm.**
5. ~~Rim wall~~ **Decided 22 October 2026: 1.5 mm at the equator, to mimic real skin** (§4a). Confirm with the workshop that 1.5 mm can be machined reliably; 2.0 mm is the fallback.
6. **Oil expansion:** about 0.8 ml per K for ≈ 1.09 L (at a typical 7 × 10⁻⁴ K⁻¹). How does the current L2 handle it? An expansion bellows on a port is an option.
7. ~~STEP export~~ **Done.**
8. **Paddle load:** 46.9 N in July. A 5 mm top wall supported by the sphere bed is assumed to be adequate.
9. **Chest-wall side of the Revelation:**
   - Is the support's chest-wall face flat and vertical over the 16 mm lip depth?
   - Is the bracket clear of anything on that face?
   - Do the 27 mm overhang and the lid clear the paddle's chest-wall flange and the face shield?
10. ~~Positioning reference~~ **Decided: positioning bracket with lip.**
11. **Disc layout:** keep it as measured from the container face (default; discs appear ≈ 5 mm + d_miss closer to the image edge than in July), or shift it to reproduce the July image positions?
12. **Missed tissue d_miss** at the chest wall: measure it, e.g. with a marker held against the support face. It sets where the image starts. To my knowledge the European guidelines limit missed tissue to 5 mm; to confirm.

## 9. Files

| File | Content |
|---|---|
| `step/L2v2_*.step` | **Exact solids for the workshop:** container, lid, gasket, plate, bracket, and the coloured assembly. |
| `stl/L2v2_*.stl` | STL tessellated from the exact solids (0.02 mm tolerance). |
| `l2v2_geom.py` | Parameters, implicit geometry and design checks; writes `params.json`, `checks.json`, `disc_layout.csv`. |
| `l2v2_step.py` | CadQuery script: exact solids, STEP/STL export; `--check` runs the cross-check. Needs `pip install cadquery`. |
| `drawings/L2v2_drawing.pdf` / `.png` | Dimensioned concept drawing: plan, chest-wall view, sections A–A and C–C, parameters. |
| `drawings/L2v2_edge_profile.png` | Edge profile and path lengths: current L2 vs v2 vs an idealised breast. |
| `renders/L2v2_render_*.png` | Illustrations: on the support, exploded, cut-away with spheres, bracket seen from below. |
| `spheres/` | Background spheres: STL and STEP per diameter, order sheet and print specification (`spheres/README.md`); made by `make_spheres.py` / `make_spheres_step.py`. |
| `MATERIALS.md` | Sphere and liquid materials, properties, amounts, PMMA-safe oil options. |
| `disc_layout_measured.json`, `disc_layout.csv` | Measured July disc positions (x from the container's chest-wall face). |
| `l2v2_render.py`, `l2v2_render_bracket.py`, `l2v2_drawing.py` | Scripts for the renders and drawings. |

The rev A OpenSCAD file was removed. It was never test-rendered and no longer matches the design. The CadQuery script is the parametric CAD source.

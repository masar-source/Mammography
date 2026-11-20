# L2 v2 background spheres: print files and order sheet (22 October 2026)

**Files:**
- `L2v2_sphere_D*.stl` and `.step`, one per diameter: 15.88 / 12.70 / 9.52 / 6.35 / 3.18 / 1.58 mm. Units are mm and each sphere is centred on the origin.
- `sphere_order_sheet.csv`: quantities for the L2 v2 container.
- `sphere_checks.json`: STL checks.

**STL checks:**
- 20,480 triangles each.
- Watertight, with outward normals.
- Mesh volume equal to the ideal sphere volume.
- Vertex diameters within +0.003 mm of nominal.

The STEP files are exact spheres from CadQuery.

## Order quantities (L2 v2 container, 61 % spheres by volume, equal volume per size)

| Diameter (mm) | 15.88 | 12.70 | 9.52 | 6.35 | 3.18 | 1.58 |
|---|---|---|---|---|---|---|
| Quantity, nominal | 54 | 104 | 247 | 832 | 6,618 | 53,951 |
| Quantity with 10 % spare | 59 | 115 | 272 | 915 | 7,280 | 59,346 |
| Volume / mass per size | 111 ml / 117 g | ← same for all sizes | | | | |

The diameters are those published for the L1 phantom; the L2 papers state only "6 sizes, 1.6–16 mm, equal volumes". **Confirm them against the current sphere set.**

## Specification to send with the files

| Item | Requirement |
|---|---|
| Process | Laser sintering (SLS) of PA12 powder, dry-blended with CaCO₃, all sizes in the same build and powder lot |
| Material | PA12, natural (unpigmented), **5.5 ± 0.5 wt % CaCO₃** (5.8 % if high-oleic sunflower oil is used instead of MCT), particle size ≤ 10 µm |
| Density | **1.046 ± 0.01 g/cm³ after vacuum saturation in the filling oil** (Archimedes, ≥ 20 pieces per size) |
| Filler check | ash content (TGA, ISO 3451) |
| Geometry | solid (no internal voids), no supports; diameter tolerance to be stated by the supplier, especially for 1.58 and 3.18 mm |
| Finish | powder fully removed (blasting or tumbling); no coating |

## Points to clear with the print service

1. **Accuracy at 1.58 and 3.18 mm.** The typical SLS tolerances I have seen quoted are in the ±0.2–0.3 mm range (unverified for any specific service). At 1.58 mm that is ±15–20 % in diameter. Ask for measured diameters on a sample, or consider sourcing the two smallest sizes differently. Dropping the 1.58 mm size is also an option under consideration (REPORT §11.4).
2. **Handling ~60,000 small spheres.** Depowdering and counting are laborious. Agree on delivery by weight per size (117 g each) rather than by count.
3. **Porosity.** SLS PA12 takes up oil. Ask for the part density, and measure it again after oil saturation.

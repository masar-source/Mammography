# L2 v2: sphere and liquid materials (22 October 2026)

**Recommended background:** PA12 + 5.5 wt % CaCO₃ spheres (laser-sintered) in MCT oil, **61 : 39 by volume**, in 6 sphere sizes of equal volume. The alternatives and references below are for comparison.

**Sources and status:**
- Compositions are from REPORT §6/§10; μ is calculated from NIST XCOM data.
- μ is mono-energetic at the mean detected energies of the TiCEM beams: LE 28 kV Rh, **21.1 keV**; HE 49 kV Ti, **37.2 keV**. These values are indicative only.
- The residual, contrast power and glandular equivalence come from the calibrated spectral model (REPORT §10.3) **[E]**. The calibration-tank test still has to confirm them.
- Script: `investigation_0035_0037/scripts/materials_table.py`.

## 1. Sphere materials

| Material | Use | Density (g/cm³) | Composition (wt %) | e⁻ density (rel. water) | Z_eff | μ at LE (cm⁻¹) | μ at HE (cm⁻¹) |
|---|---|---|---|---|---|---|---|
| **PA12 + 5.5 % CaCO₃** | **recommended** (SLS, all sizes) | **1.046** | C 69.7, H 11.1, O 10.3, N 6.7, Ca 2.2 | 1.045 | 7.08 | 0.715 | 0.297 |
| PA12 + 3.0 % TiO₂ | alternative filler | 1.034 | C 70.9, H 11.4, O 9.1, N 6.9, Ti 1.8 | 1.033 | 7.05 | 0.702 | 0.294 |
| HIPS + 6.5 % CaCO₃ | FDM, spheres ≥ ~10 mm only | 1.094 | C 87.0, H 7.2, O 3.1, Ca 2.6 | 1.054 | 7.21 | 0.757 | 0.307 |
| ABS + 6.5 % CaCO₃ | FDM, spheres ≥ ~10 mm only | 1.094 | C 81.8, H 7.5, N 4.9, O 3.1, Ca 2.6 | 1.057 | 7.24 | 0.765 | 0.309 |
| CIRS glandular-equivalent epoxy | if custom spheres can be bought (unconfirmed) | 1.05 | C 69.1, O 17.7, H 9.4, N 1.8, Ca 1.75, Cl 0.14 | 1.032 | 7.06 | 0.696 | 0.291 |
| *Glandular tissue* | *reference (Hammerstein)* | *1.04* | *O 67.7, C 18.4, H 10.2, N 3.2, P 0.5* | *1.031* | *7.21* | *0.695* | *0.289* |
| *PMMA* | *current* | *1.19* | *C 60.0, O 32.0, H 8.0* | *1.156* | *6.47* | *0.613* | *0.295* |

Densities of the doped materials follow the mixing rule and assume no voids. **Specify them after oil saturation:** SLS PA12 is porous.

## 2. Liquids

| Liquid | Use | Density (g/cm³) | Composition (wt %) | e⁻ density (rel. water) | Z_eff | μ at LE (cm⁻¹) | μ at HE (cm⁻¹) | Notes |
|---|---|---|---|---|---|---|---|---|
| **MCT oil** (caprylic/capric triglyceride, Ph. Eur.) | **recommended** | **0.945** | C 69.8, O 19.3, H 10.9 | 0.942 | 6.07 | 0.442 | 0.230 | 25–33 mPa·s; oxidation-stable; an ester, so a PMMA stress-cracking test is mandatory |
| Olive / high-oleic oil | alternative | 0.915 | C 77.3, H 11.8, O 10.8 | 0.920 | 5.83 | 0.401 | 0.219 | texture too strong with these spheres |
| *Liquid paraffin* | *current* | *0.86 (assumed)* | *C 85.2, H 14.8 (assumed)* | *0.887* | *5.43* | *0.342* | *0.203* | |
| *Adipose tissue* | *reference (Hammerstein)* | *0.93* | *C 61.9, O 25.1, H 11.2, N 1.7, P 0.1* | *0.930* | *6.24* | *0.461* | *0.231* | |

The MCT composition is approximated as C₂₉H₅₄O₆.

## 3. Sphere/liquid pairs at 61 : 39 by volume

| Pair | Δμ LE / HE (cm⁻¹) | Δμ_HE / Δμ_LE | Residual under Siemens' weight, 28 / 30 / 32 / 34 kV | LE texture power vs current L2 | Glandular-equivalent (LE) | Δρ sphere − liquid (g/cm³) |
|---|---|---|---|---|---|---|
| **PA12 + 5.5 % CaCO₃ / MCT** | 0.273 / 0.067 | **0.246** | **+0.07 / +0.09 / +0.08 / +0.12** | 0.94 | **63 %** | +0.10 (spheres sink) |
| PA12 + 3.0 % TiO₂ / MCT | 0.259 / 0.064 | 0.246 | +0.07 / +0.09 / +0.07 / +0.11 | 0.86 | 61 % | +0.09 |
| HIPS or ABS + 6.5 % CaCO₃ / MCT | 0.315 / 0.077 (HIPS) | 0.245 | +0.07 / +0.09 / +0.08 / +0.12 | 1.24 | 72–74 % | +0.15 |
| CIRS glandular epoxy / MCT | 0.254 / 0.061 | 0.242 | +0.06 / +0.07 / +0.06 / +0.10 | 0.82 | 60 % | +0.11 |
| PA12 + 5.75 % CaCO₃ / olive oil | — | — | +0.07 / +0.09 / +0.08 / +0.12 | 1.37 | 61 % | — |
| *Glandular / adipose (target)* | *0.234 / 0.058* | ***0.246*** | ***+0.07 / +0.09 / +0.08 / +0.12*** | *0.69* | *62 %* | *+0.11* |
| *PMMA / paraffin (current)* | *0.271 / 0.091* | *0.338* | *+0.29 / +0.29 / +0.26 / +0.28* | *1.00* | *29 %* | *+0.33* |

The Δμ_HE/Δμ_LE ratio is what the dual-energy weight has to match. The recommended pair matches the tissue target (0.246); the current pair does not (0.338).

## 4. Amounts for the L2 v2 container

The filling volume is **1096 ml**, from the STEP solids (walls as measured, 22 Oct 2026). It includes 187 ml in the extension beyond the chest wall.

| Item | Volume | Mass |
|---|---|---|
| Spheres, PA12 + 5.5 % CaCO₃ (61 %) | 669 ml | ≈ 699 g (≈ 661 g PA12 + 38 g CaCO₃) |
| MCT oil (39 %) | 427 ml | ≈ 404 g |
| Per half (above / below the plate) | 334 ml spheres + 214 ml oil | — |

**Per sphere size** (equal volumes, 111 ml ≈ 117 g each):

| Diameter (mm) | 15.88 | 12.70 | 9.52 | 6.35 | 3.18 | 1.58 |
|---|---|---|---|---|---|---|
| Number of spheres | 54 | 104 | 247 | 832 | 6,618 | 53,951 |

**Notes:**
- **Sizes are an assumption.** The L2 papers state only "6 diameters, 1.6–16 mm, equal volumes". The listed set (5/8″ … 1/16″) is the one published for the L1 phantom (Cockmartin 2017). Confirm against the actual sphere set.
- **61 % is the bulk packing of the current L2.** Near the rounded rim and in the extension corners only small spheres fit, so somewhat less than 669 ml will fit. Weigh equal portions per size, split them equally between the two halves, and fill the oil to the brim rather than to a fixed volume.
- **Buy extra oil:** ≈ 0.6 L of MCT. Porous SLS spheres take up oil, and topping up is needed.
- **1.58 mm spheres:** whether SLS can make them accurately is open (test plan). Dropping them was suggested to bring β closer to patients (REPORT §11.4).

## 5. Oils that do not attack PMMA (added 22 October 2026)

**Why this matters.** MCT oil is a medium-chain ester. Acrylic chemical-resistance tables rate short esters (acetates, butyrate, benzoate) "not recommended". They rate mineral oil, silicone oil, glycols, glycerine, fatty acids and a long-chain oleate ester "excellent" at 20 and 50 °C. Cooking oils are listed as resistant. No data were found for MCT.

The tables are secondary copies of a PLEXIGLAS list, not Röhm's own document: [eplastics](https://www.eplastics.com/blog/chemical-resistance-acrylic-plexiglass), [Darwin Microfluidics](https://blog.darwin-microfluidics.com/pmma-polymethylmethacrylate-chemical-resistance-chart/). Whichever oil is chosen still needs the stress-cracking test on the bonded container.

**Matching the spheres to each PMMA-safe oil.** The sphere filler was re-tuned so that the contrast ratio Δμ_HE/Δμ_LE equals glandular/adipose (0.246), which gives a tissue-like residual under Siemens' weight. The table compares what results.

- **Model:** mono-energetic estimates at 21.1 / 37.2 keV, 61 : 39 by volume (`scripts/oil_alternatives.py`, output in `oil_alternatives.txt`).
- **Agreement with the spectral model:** where REPORT §10.3 has a spectral result for the same pair, they agree. For the olive-oil pair the spectral model gives 5.75 % CaCO₃, ×1.37 texture power and 61 % glandular.

| Oil (PMMA rating) | Spheres | Density (g/cm³); spheres sink by | LE contrast (× tissue) | Glandular-equivalent | Verdict |
|---|---|---|---|---|---|
| **High-oleic sunflower oil** (long-chain triglyceride; cooking oils "resistant") | **PA12 + 5.8 % CaCO₃** | 1.048; +0.13 | 0.33 /cm (×1.4) | **61 %** | **best PMMA-safe option**: dense and tissue-like residual; texture power ≈ 1.4–1.5× that of the MCT design, at the top of the patient range |
| High-oleic sunflower oil | PA12 + 3.1 % TiO₂ | 1.035; +0.12 | 0.31 /cm (×1.3) | 56 % | slightly weaker texture, slightly less dense |
| High-oleic sunflower oil | HDPE + 5.7 % CaCO₃ | 0.991; +0.08 | 0.25 /cm (×1.05) | 38 % | tissue-like texture, but not a dense breast |
| Liquid paraffin, 0.86 (current oil; years of use with PMMA; mineral oil "excellent") | PA12 + 6.0–6.4 % CaCO₃ | 1.05; +0.19 | 0.42 /cm (×1.8) | 56–59 % | proven compatibility, but texture ≈ 2× too strong |
| Heavy liquid paraffin, 0.88 | HDPE + 5.1 % CaCO₃ | 0.987; +0.11 | 0.27 /cm (×1.15) | 22 % | tissue-like texture, but fatty, not dense |
| *MCT oil (reference; PMMA compatibility untested)* | *PA12 + 5.5 % CaCO₃* | *1.046; +0.10* | *0.27 /cm (×1.17)* | *63 %* | *radiologically best* |
| Silicone oil (PDMS) ("excellent" for PMMA) | — | — | — | — | **unsuitable:** silicon makes it attenuate ≈ 2.5× glandular tissue (μ 1.71 cm⁻¹ at 21 keV) |
| Glycols, glycerine ("excellent") | — | — | — | — | **unsuitable:** water-like attenuation, not adipose; as dense as or denser than the spheres |

**The trade-off.** All PMMA-safe oils attenuate less than adipose tissue. With 61 % spheres, they can match the tissue residual, but not both the tissue texture strength and dense-breast attenuation at once.

**Recommendation:**
1. **Keep MCT if the stress-cracking test passes;** it is radiologically the best match.
2. **If PMMA safety must be guaranteed without testing, use high-oleic sunflower oil with PA12 + 5.8 % CaCO₃ spheres.** The CaCO₃ level lies within the rib-plate range already planned (4.5 / 5.5 / 6.5 %), so the calibration tank only needs a second oil.
3. **Handling of sunflower oil:**
   - use a refined, food- or pharma-grade high-oleic oil with tocopherol;
   - keep the phantom sealed and free of air, and away from light;
   - re-check the density and colour every year.

   It oxidises more easily than MCT or paraffin. Its viscosity is roughly 2–3× that of MCT (typical values; confirm with the supplier's data sheet).
4. **Amounts for the v2 container with this option:** 669 ml of spheres (≈ 701 g, including ≈ 41 g CaCO₃) and 427 ml of oil (≈ 391 g).

### 5a. Keeping the current PMMA spheres: no oil can do it

This was checked with the same mono-energetic estimate (`scripts/pmma_liquids.py`). Glandular-equivalent percentages from this estimate run ≈ 9 points below the spectral model: the current L2 comes out at 20 % here against 29 % in the spectral model.

**Cancellation.** A tissue-like residual with PMMA spheres needs a liquid with μ_HE/μ_LE ≈ 0.63–0.67. That is a hydrogen-rich hydrocarbon at a density of ≈ 1.04–1.05 g/cm³, and no such liquid exists: liquid paraffins are 0.83–0.89. The candidates give:

| Liquid | Contrast ratio with PMMA (target 0.246) |
|---|---|
| Heavy paraffin (0.89) | 0.327, the best available |
| Liquid paraffin (0.86, current) | 0.338 |
| Squalane | 0.356 |
| High-oleic sunflower | 0.357 |
| Castor oil | 0.357 |
| MCT | 0.380 |

Oxygen-containing oils are denser but move the ratio the wrong way.

**Density.** PMMA alone is only ≈ 65 % glandular-equivalent at LE. With 61 % spheres, the mixture is PMMA minus 0.39 × the contrast. A tissue-like contrast therefore caps the mixture at ≈ 20–26 % glandular (≈ 30–35 % spectral). Reaching ≈ 60 % needs a liquid nearly as attenuating as PMMA, which leaves almost no texture: propylene glycol gives 58 % glandular but 0.18 × tissue contrast.

**K-edge dopants.** Iodine or tin dissolved in the oil could force the ratio to 0.246, at about 0.15 wt % iodine. They are rejected:
- the result is extremely sensitive: 0.2 % already overshoots;
- it depends strongly on kV;
- the liquid would carry an iodine-like signal into the recombined image.

**Conclusion.** With PMMA spheres the best oil (heavy paraffin) only lowers the residual from ≈ 0.28 to ≈ 0.26 **[E]**. Tissue-like cancellation plus dense-breast attenuation require new spheres (§3, §5).

### 5b. Low-cost option: stock nylon balls in glycerol / propylene glycol (added 22 October 2026)

**Why.** Custom-blended PA12 + CaCO₃ spheres are expensive. This option needs no printing.

**How it works.** Texture statistics (power spectrum, β, κ) do not depend on whether the spheres look brighter or darker than the liquid. A liquid more attenuating than the spheres is therefore allowed. Glycerol is the only common liquid whose energy dependence lies as far above the tissue line as nylon does, so a glycerol-rich blend is needed.

The table is a mono-energetic estimate (`scripts/cheap_pairs.py`, `scripts/cheap_pairs_check.py`). It supersedes the "nylon in glycerol" entry rejected in REPORT §10.6, which was tuned for zero residual.

| Spheres | Liquid (by mass) | Liquid density (g/cm³) | Contrast ratio (tissue 0.246) | Contrast | Glandular-equivalent |
|---|---|---|---|---|---|
| **PA6 nylon balls** (≈ 1.13) | **glycerol / PG 82 / 18** | 1.214 | 0.246 | 0.75 × tissue | ≈ 58 % |
| **PA66 nylon balls** (≈ 1.14) | **glycerol / PG 88 / 12** | 1.229 | 0.247 | 0.78 × tissue | ≈ 62 % |
| PC balls (1.20) | glycerol / PG 90 / 10 | 1.234 | 0.247 | 0.67 × tissue | ≈ 70 % |
| PMMA balls (current) | glycerol-rich (≈ 92 % glycerol) | 1.25 | 0.24 | 0.73 × tissue | ≈ 93 %, too dense |

- **Texture strength:** ≈ 0.4 × the current L2's power, within the published patient range [R28].
- **Appearance:** the balls look *darker* than the liquid, like fat lobules in dense parenchyma.
- **Balls float** (1.13–1.14 in ≈ 1.21–1.23) and pack against the top of each half. The v2 plate is held mechanically. In the current L2 container a PMMA plate (1.19) would also float.
- **PMMA-safe:** glycerol and propylene glycol are rated "excellent" for acrylic.
- **Tunable after purchase.** Measure the actual ball density first (one listing gives 1.19 g/cm³ for nylon 6/6), then set the glycerol/PG ratio: ±3 % glycerol shifts the ratio by ≈ ±0.005.

**Handling:**
- **Viscous:** ≈ 0.5–0.8 Pa·s at 20 °C (estimate; glycerol alone is ≈ 1.4 Pa·s). Fill and degas warm (≈ 40 °C) under vacuum. Spheres rearrange more slowly when shaken.
- **Hygroscopic:** keep the phantom sealed. 3 % absorbed water shifts the ratio by ≈ −0.007.
- **Nylon takes up moisture:** 4 % water shifts the ratio by ≈ +0.01. Pre-soak the balls in the liquid for several weeks, then measure density and diameter.
- **Shelf life:** the blend is too dry for microbial growth (glycerol-rich, with PG).

**Amounts for the v2 container:** 669 ml of balls (≈ 126 g per size for PA6, same counts as `spheres/sphere_order_sheet.csv`) and 427 ml of liquid (≈ 520 g).

**Sources (sizes to confirm):** nylon 6/6 balls are catalogue items in 1/16″–3/16″ ([McMaster-Carr](https://www.mcmaster.com/beads/comparable-to~ertalon/), [Polysciences 1/16″](https://polysciences.com/products/polyballs-polyamide-66-nylon-116-diameter), [Miniature Bearings](https://minibearings.com.au/collections/balls-plastic-nylon)). Confirm 1/4″–5/8″ with a ball supplier (e.g. those listed in REPORT [R26]).

**Cheapest validation.** Refill the existing L2 (or a small tank) with nylon balls in the glycerol/PG blend. Image at 28–32 kV and run the existing analysis (w*, residual, texture, Siemens image) before committing to the v2 build. Test plan: `investigation_0035_0037/NYLON_REFILL_TEST_PLAN.md`.

### 5c. Resin (SLA/DLP) printed spheres: not recommended

- **In oil:** methacrylate resins behave like PMMA. The contrast ratio is ≈ 0.38–0.40 (tissue 0.246), and CaCO₃ cannot correct this for oxygen-rich materials (REPORT §10.3).
- **In the glycerol liquid:** the ratio can be matched, but the background becomes ≈ 80–83 % glandular-equivalent with only ≈ 0.5 × tissue contrast.
- **Composition unknown:** resin compositions are proprietary and would have to be measured.
- **Impractical to print:** ≈ 60,000 spheres down to 1.58 mm each need supports, washing and post-curing. Resin printing is useful at most for a few large test spheres.

## 6. Specification checks before production

| Check | Requirement |
|---|---|
| CaCO₃ content | 5.5 ± 0.5 wt %, by ash content (TGA, ISO 3451) |
| Filler particle size | CaCO₃ ≤ 10 µm, evenly dispersed; natural (unpigmented) PA12 |
| Sphere density | 1.046 ± 0.01 g/cm³ **after vacuum saturation in MCT** (Archimedes) |
| Compatibility | Soak tests of PA12 + CaCO₃ in MCT; stress-cracking test of the bonded PMMA container with MCT; FKM gasket with MCT |
| Tolerance effect [E] | ± 0.5 wt % CaCO₃ → r ± 0.01–0.02; ± 0.01 g/cm³ → r ± 0.015 |
| Radiological confirmation | Calibration-tank test plan: rib plates at 4.5 / 5.5 / 6.5 % CaCO₃, predicted densities 1.039 / 1.046 / 1.053 |

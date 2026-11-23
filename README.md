# Mammography

This repository holds two of my projects on the Siemens MAMMOMAT Revelation. Both use the same GATE 10 Monte Carlo
model of the system. Each folder has its own README and runs on its own.

| Folder | What it is |
|---|---|
| [`L2-CEM-phantom/`](L2-CEM-phantom/) | The L2 phantom for contrast-enhanced mammography: a tuned recombination of the dual-energy images (Siemens and Hologic), a redesigned container (L2 v2) and the GATE 10 model of the Revelation with the phantom. |
| [`radiation-protection/`](radiation-protection/) | A GATE 10 room model for planning the shielding of a Revelation room before installation: the scanner model, patient-derived compressed breasts (Radboud dataset, CC BY 4.0), walls, door and console shield, air kerma scored at staff and neighbouring positions, and an offline 3D viewer to place the scanner and the console. |

## Status

**L2-CEM-phantom.** The tuned recombination is finished and was applied to 112 Siemens and 15 Hologic exposures. The
GATE model agrees with a measured L2 exposure within 3 % inside the phantom. The redesigned container has not been
built yet.

**radiation-protection.** The model runs, and its parts were checked: the air-kerma scoring, the source
normalisation, the room geometry and the registration of the breast images. In a first pilot with the L2 phantom, the
scatter 1 m from the unit came out at the level that NCRP Report No. 147 gives for mammography rooms. The positions
behind the console shield and outside the room still need many more histories, or a two-step method, before they
give usable numbers.

## Requirements

Python 3 with numpy, scipy and the packages listed in each folder. The simulations need opengate (GATE 10.1).

## Licence

The code, models and documents in this repository are under the MIT licence (`LICENSE`). Two parts keep their own
licences: the six Radboud breast images in `radiation-protection/anatomy/` and the label images made from them
(CC BY 4.0; Pacheco, Michielsen and Sechopoulos, see the attribution in `radiation-protection/README.md`), and the
three.js library and the Inter font in the 3D viewer (MIT and SIL Open Font License, see
`radiation-protection/Third_party_notices.txt`).

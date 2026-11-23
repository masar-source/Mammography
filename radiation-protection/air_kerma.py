"""Collision air kerma from transported photon track segments in dry air."""

import numpy as np

# NIST X-Ray Mass Attenuation Coefficients, dry air near sea level, table 4.
# https://physics.nist.gov/PhysRefData/XrayMassCoef/ComTab/air.html
# Columns are MeV and mass energy-absorption coefficient in cm2/g.
TABLE = np.array(
    [
        [0.001, 3599],
        [0.0015, 1188],
        [0.002, 526.2],
        [0.003, 161.4],
        [0.00320290, 133.0],
        [0.00320290, 146.0],
        [0.004, 76.36],
        [0.005, 39.31],
        [0.006, 22.70],
        [0.008, 9.446],
        [0.010, 4.742],
        [0.015, 1.334],
        [0.020, 0.5389],
        [0.030, 0.1537],
        [0.040, 0.06833],
        [0.050, 0.04098],
        [0.060, 0.03041],
        [0.080, 0.02407],
        [0.100, 0.02325],
    ],
    dtype=float,
)


def mu_en(energies):
    energy = np.asarray(energies, dtype=float)
    if np.any((energy < TABLE[0, 0]) | (energy > TABLE[-1, 0])):
        raise ValueError("Photon energy outside the dry-air coefficient table")
    # The second value at the duplicated argon edge applies just above the edge.
    left = np.clip(
        np.searchsorted(TABLE[:, 0], energy, side="right") - 1, 0, len(TABLE) - 2
    )
    e0, c0 = TABLE[left].T
    e1, c1 = TABLE[left + 1].T
    return c0 * np.exp(np.log(energy / e0) / np.log(e1 / e0) * np.log(c1 / c0))


def read_cell(tree, volume_mm3, per_mAs):
    total = 0.0
    scattered_total = 0.0
    steps = 0
    below = 0
    for a in [] if tree is None else tree.iterate(step_size="16 MB", library="np"):
        photons = a["PDGCode"] == 22
        length = np.linalg.norm(
            np.column_stack(
                [a["PostPosition_" + k] - a["PrePosition_" + k] for k in "XYZ"]
            ),
            axis=1,
        )
        energy = a["PreKineticEnergy"]
        below += int(np.count_nonzero(photons & (length > 0) & (energy < TABLE[0, 0])))
        use = photons & (length > 0) & (energy >= TABLE[0, 0])
        e = energy[use]
        direction_change = np.column_stack(
            [a["PreDirection_" + k][use] - a["EventDirection_" + k][use] for k in "XYZ"]
        )
        scattered = (
            (a["ParentID"][use] != 0)
            | (np.linalg.norm(direction_change, axis=1) > 1e-6)
            | (abs(e - a["EventKineticEnergy"][use]) > 1e-8)
        )
        contribution = a["Weight"][use] * e * length[use] * mu_en(e)
        total += float(contribution.sum(dtype=float))
        scattered_total += float(contribution[scattered].sum(dtype=float))
        steps += len(e)
    # MeV->J, mm/mm3->1/m2, cm2/g->m2/kg, Gy->microGy.
    factor = 1.602176634e-13 * 1e6 * 0.1 * 1e6 * per_mAs / volume_mm3
    return {
        "collision_air_kerma_uGy_mAs": total * factor,
        "secondary_field_uGy_mAs": scattered_total * factor,
        "photon_steps": steps,
        "sub_keV_steps": below,
    }

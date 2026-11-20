"""
Everything below the phase-space plane, applied per photon: air, focused linear anti-scatter grid (moving, so averaged
over the lamella phase), detector cover and a-Se absorption. Used identically for the Monte Carlo photons
(postprocess.py) and for the deterministic primary image (projector.py).

Stack below the plane [A] (z in mm, a-Se surface at z = 0):
  plane z = 14.4 | air | grid top cover CF 0.5 | grid (Pb lamellae + CF interspace, height h) | bottom cover CF 0.5 |
  air | detector cover CF 1.0 at z = 1..2 | air | a-Se 200 um below z = 0
Signal = E x (fraction absorbed in a-Se): energy-integrating detector; Se K-fluorescence escape neglected [A].

Grid: lamellae parallel to x, focused at the focal spot in the y-z plane. A photon whose slope in the y-z plane differs
from the lamella slope at its position by d crosses the grid height h with a lateral shift s = h |d|; the fraction of its
path in lead, averaged over where it enters (the grid moves), gives the transmission T(E, s, L), L = path length.
"""
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from config import SYSTEM, GRID
from materials import mu

PLANE_Z = 14.4                    # must match revelation_gate.PLANE_Z[1]
GRID_TOP = 6.533                  # [A] top of the grid assembly (top cover)
DET_COVER = (1.0, 2.0)            # [A]

D = 10.0 / GRID["lines_per_cm"]                     # period (mm)
d = GRID["lamella_um"] / 1000.0                     # lamella thickness (mm)
h = GRID["ratio"] * (D - d)                          # lamella height (mm)
c = GRID["cover_mm_each"]
Z_GRID_TOP, Z_GRID_BOT = GRID_TOP - c, GRID_TOP - c - h
Z_GRID_MID = 0.5 * (Z_GRID_TOP + Z_GRID_BOT)
F_FOCUS = GRID["focus_mm"] - Z_GRID_MID             # lamellae point at the focal spot


def _lead_fraction(p, s):
    """mean fraction of the lateral interval [p, p+s] lying in lead (lead occupies [kD, kD+d)); s > 0"""
    G = lambda x: np.floor(x / D) * d + np.minimum(np.mod(x, D), d)
    return (G(p + s) - G(p)) / s


class GridTable:
    """T(E, s, L/h) on a grid; trilinear interpolation"""

    def __init__(self, e_kev, n_phase=256):
        self.e = np.asarray(e_kev, float)
        self.s = np.r_[0.0, np.geomspace(1e-4, 6 * D, 120)]
        self.l = np.linspace(1.0, 1.25, 11)                       # L/h for |uz| >= 0.8 (theta <= 37 deg)
        mpb, mi = mu(GRID["lamella_material"], self.e), mu(GRID["interspace_material"], self.e)
        p = (np.arange(n_phase) + 0.5) / n_phase * D
        f = _lead_fraction(p[None, :], self.s[1:, None])                 # (s>0, phase)
        att = mpb[:, None, None] * f[None] + mi[:, None, None] * (1 - f[None])   # (E, s>0, phase)
        T = np.empty((len(self.e), len(self.s), len(self.l)))
        for k, L in enumerate(self.l * h):
            T[:, 1:, k] = np.exp(-att * L).mean(-1)
            T[:, 0, k] = (1 - d / D) * np.exp(-mi * L) + d / D * np.exp(-mpb * L)   # s = 0 exactly
        self.T = T
        self._interp = RegularGridInterpolator((self.e, self.s, self.l), T, bounds_error=False, fill_value=None)

    def __call__(self, e_kev, s, l_over_h):
        pts = np.c_[np.broadcast_to(e_kev, np.shape(s)).ravel(), np.clip(np.ravel(s), 0, self.s[-1]),
                    np.clip(np.ravel(l_over_h), 1.0, self.l[-1])]
        return self._interp(pts).reshape(np.shape(s))


class Response:
    """detected energy per photon below the plane: R = E * T_air * T_covers * T_grid * A_Se"""

    def __init__(self, e_kev, grid=True):
        self.e = np.asarray(e_kev, float)
        self.grid = GridTable(self.e) if (grid and GRID["enabled"]) else None
        self.mu_air = mu("G4_AIR", self.e)
        self.mu_cf = mu(GRID["cover_material"], self.e)
        self.mu_dc = mu(SYSTEM["detector_cover"]["material"], self.e)
        self.mu_se = mu("aSe", self.e)
        self.mu_pb, self.mu_i = mu(GRID["lamella_material"], self.e), mu(GRID["interspace_material"], self.e)
        self.t_se = SYSTEM["a_se_um"] / 1000.0
        self.dc = SYSTEM["detector_cover"]["mm"]
        self.air = PLANE_Z - (2 * c + h if self.grid else 0) - self.dc

    def weight(self, e_kev, pos, u):
        """detected energy for Monte Carlo photons at the plane: energies (n,) in keV, positions (n,3), unit
        directions (n,3) with uz < 0. Attenuation coefficients at the exact energies; grid table interpolated.
        Returns (weight in keV, x, y) with x, y where the photon reaches z = 0."""
        e = np.asarray(e_kev, float)
        uz = -u[:, 2]
        tx, ty = u[:, 0] / uz, u[:, 1] / uz
        x0 = pos[:, 0] + tx * pos[:, 2]; y0 = pos[:, 1] + ty * pos[:, 2]
        inv = 1.0 / uz
        w = e * np.exp(-(mu("G4_AIR", e) * self.air + mu(SYSTEM["detector_cover"]["material"], e) * self.dc) * inv)
        if self.grid is not None:
            yg = pos[:, 1] + ty * (pos[:, 2] - Z_GRID_MID)
            s = h * np.abs(ty - yg / F_FOCUS)
            w = w * np.exp(-mu(GRID["cover_material"], e) * 2 * c * inv) * self.grid(e, s, inv)
        w = w * (1 - np.exp(-mu("aSe", e) * self.t_se * inv))
        return w, x0, y0

    def primary_factor_e(self, k, inv):
        """detected energy per photon (keV) at energy self.e[k] for primary rays from the focal spot (focused grid:
        no lateral shift, closed form), inv = 1/|uz| (any shape)"""
        out = self.e[k] * np.exp(-(self.mu_air[k] * self.air + self.mu_dc[k] * self.dc) * inv)
        if self.grid is not None:
            L = h * inv
            g = (1 - d / D) * np.exp(-self.mu_i[k] * L) + d / D * np.exp(-self.mu_pb[k] * L)
            out = out * np.exp(-self.mu_cf[k] * 2 * c * inv) * g
        return out * (1 - np.exp(-self.mu_se[k] * self.t_se * inv))


E_GRID = np.round(np.arange(1.0, 50.05, 0.1), 1)       # keV; energy axis of the grid table for Monte Carlo photons


def describe():
    e = np.array([20.0, 30.0, 40.0])
    r = Response(e)
    tp = r.grid(e, np.zeros(3), np.ones(3))
    ts = [r.grid(e, np.full(3, h * np.tan(np.radians(a))), np.ones(3)) for a in (5, 10, 20)]
    print("grid: period %.4f mm, lamella %.3f mm, height %.3f mm, ratio %.1f" % (D, d, h, GRID["ratio"]))
    print("primary transmission at 20/30/40 keV:", np.round(tp, 3))
    for a, t in zip((5, 10, 20), ts):
        print("transmission at %2d deg off focus:" % a, np.round(t, 3))


if __name__ == "__main__":
    describe()

"""Reader for VMEC ``wout`` equilibrium files (NetCDF).

Conventions (VMEC):

* ``R(s, theta, phi) = sum rmnc * cos(m*theta - xn*phi)``
* ``Z(s, theta, phi) = sum zmns * sin(m*theta - xn*phi)``
* ``xn`` already contains the number of field periods, ``xn = n * nfp``.
* ``phi`` is the cylindrical toroidal angle.
* Radial index 0 is the magnetic axis, the last index is the boundary.
  Geometry coefficients live on the full radial grid, which is uniform in
  normalised toroidal flux ``s``.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from pathlib import Path

import numpy as np


class EquilibriumError(ValueError):
    """The equilibrium file is missing, unreadable or not usable."""


@dataclass(frozen=True)
class Equilibrium:
    """Flux-surface geometry and axis field of one equilibrium."""

    nfp: int
    lasym: bool
    xm: np.ndarray  # [mn] poloidal mode numbers
    xn: np.ndarray  # [mn] toroidal mode numbers, including nfp
    rmnc: np.ndarray  # [ns, mn]
    zmns: np.ndarray  # [ns, mn]
    rmns: np.ndarray | None  # [ns, mn], only if lasym
    zmnc: np.ndarray | None  # [ns, mn], only if lasym
    iota: np.ndarray  # [ns] rotational transform on the full grid
    #: Toroidal mode numbers (including nfp) of the field strength on axis.
    b_axis_xn: np.ndarray
    #: cos coefficients of the field strength on axis [T].
    b_axis_cos: np.ndarray
    #: sin coefficients of the field strength on axis [T]; zeros if not lasym.
    b_axis_sin: np.ndarray
    #: Scalars as reported in the file, kept for cross-checks and the report.
    reported: dict = field(default_factory=dict)
    source: str = ""

    # ------------------------------------------------------------------ grid
    @property
    def ns(self) -> int:
        return self.rmnc.shape[0]

    @property
    def s_grid(self) -> np.ndarray:
        """Normalised toroidal flux of each radial grid point."""
        return np.linspace(0.0, 1.0, self.ns)

    # ------------------------------------------------------------ evaluation
    def coefficients(self, s: float):
        """Fourier coefficients at normalised flux ``s``.

        Linear interpolation between the two neighbouring grid surfaces.
        Returns ``(rmnc, zmns, rmns, zmnc)``; the last two are ``None`` for
        stellarator-symmetric equilibria.
        """
        if not 0.0 <= s <= 1.0:
            raise EquilibriumError(f"s must be in [0, 1], got {s}")
        x = s * (self.ns - 1)
        i0 = min(int(np.floor(x)), self.ns - 2)
        w = x - i0

        def interp(arr):
            if arr is None:
                return None
            return (1.0 - w) * arr[i0] + w * arr[i0 + 1]

        return interp(self.rmnc), interp(self.zmns), interp(self.rmns), interp(self.zmnc)

    def surface(self, s: float, theta: np.ndarray, phi: np.ndarray) -> dict:
        """Evaluate a flux surface on the tensor grid ``theta`` x ``phi``.

        Returns a dict of arrays with shape ``[len(theta), len(phi)]``:
        ``R``, ``Z`` and their derivatives ``R_t``, ``R_p``, ``Z_t``, ``Z_p``
        with respect to ``theta`` and ``phi``.
        """
        theta = np.atleast_1d(np.asarray(theta, dtype=float))
        phi = np.atleast_1d(np.asarray(phi, dtype=float))
        rmnc, zmns, rmns, zmnc = self.coefficients(s)
        m = self.xm[:, None, None]
        n = self.xn[:, None, None]
        angle = m * theta[None, :, None] - n * phi[None, None, :]
        cos, sin = np.cos(angle), np.sin(angle)

        def contract(coeff, basis):
            return np.tensordot(coeff, basis, axes=(0, 0))

        out = {
            "R": contract(rmnc, cos),
            "Z": contract(zmns, sin),
            "R_t": contract(-rmnc * self.xm, sin),
            "R_p": contract(rmnc * self.xn, sin),
            "Z_t": contract(zmns * self.xm, cos),
            "Z_p": contract(-zmns * self.xn, cos),
        }
        if self.lasym:
            out["R"] += contract(rmns, sin)
            out["Z"] += contract(zmnc, cos)
            out["R_t"] += contract(rmns * self.xm, cos)
            out["R_p"] += contract(-rmns * self.xn, cos)
            out["Z_t"] += contract(-zmnc * self.xm, sin)
            out["Z_p"] += contract(zmnc * self.xn, sin)
        return out

    def surface_xyz(self, s: float, theta: np.ndarray, phi: np.ndarray) -> np.ndarray:
        """Cartesian points of a flux surface, shape ``[len(theta), len(phi), 3]``."""
        surf = self.surface(s, theta, phi)
        phi = np.atleast_1d(np.asarray(phi, dtype=float))
        return np.stack(
            [surf["R"] * np.cos(phi)[None, :], surf["R"] * np.sin(phi)[None, :], surf["Z"]],
            axis=-1,
        )

    def axis_xyz(self, phi: np.ndarray) -> np.ndarray:
        """Cartesian points of the magnetic axis, shape ``[len(phi), 3]``."""
        return self.surface_xyz(0.0, np.array([0.0]), phi)[0]

    def axis_field_strength(self, phi: np.ndarray) -> np.ndarray:
        """Field strength on the magnetic axis [T] at toroidal angles ``phi``."""
        phi = np.atleast_1d(np.asarray(phi, dtype=float))
        angle = -self.b_axis_xn[:, None] * phi[None, :]
        return self.b_axis_cos @ np.cos(angle) + self.b_axis_sin @ np.sin(angle)

    # --------------------------------------------------------------- scaling
    def scaled(self, length_factor: float = 1.0, field_factor: float = 1.0) -> Equilibrium:
        """Return a copy with all lengths and fields multiplied by the factors."""

        def mul(arr, f):
            return None if arr is None else arr * f

        rep = dict(self.reported)
        powers = {"Rmajor_p": (1, 0), "Aminor_p": (1, 0), "volume_p": (3, 0),
                  "b0": (0, 1), "volavgB": (0, 1)}
        for key, (pl, pb) in powers.items():
            if key in rep:
                rep[key] = rep[key] * length_factor**pl * field_factor**pb
        return replace(
            self,
            rmnc=self.rmnc * length_factor,
            zmns=self.zmns * length_factor,
            rmns=mul(self.rmns, length_factor),
            zmnc=mul(self.zmnc, length_factor),
            b_axis_cos=self.b_axis_cos * field_factor,
            b_axis_sin=self.b_axis_sin * field_factor,
            reported=rep,
        )


def _scalar(ds, name, default=None):
    if name not in ds.variables:
        if default is None:
            raise EquilibriumError(f"wout file has no variable '{name}'")
        return default
    return np.asarray(ds.variables[name][...]).item()


def _array(ds, name):
    if name not in ds.variables:
        raise EquilibriumError(f"wout file has no variable '{name}'")
    return np.ma.filled(ds.variables[name][...], np.nan).astype(float)


def read_wout(path: str | Path) -> Equilibrium:
    """Read a VMEC ``wout`` NetCDF file."""
    import netCDF4

    path = Path(path)
    if not path.is_file():
        raise EquilibriumError(f"equilibrium file not found: {path}")
    try:
        ds = netCDF4.Dataset(str(path))
    except OSError as exc:
        raise EquilibriumError(f"cannot open '{path}' as NetCDF: {exc}") from exc

    with ds:
        nfp = int(_scalar(ds, "nfp"))
        lasym = bool(_scalar(ds, "lasym__logical__", 0))
        xm = _array(ds, "xm")
        xn = _array(ds, "xn")
        rmnc = _array(ds, "rmnc")
        zmns = _array(ds, "zmns")
        rmns = _array(ds, "rmns") if lasym else None
        zmnc = _array(ds, "zmnc") if lasym else None
        iota = _array(ds, "iotaf")

        # Field strength lives on the half grid: index 0 is unused, index j is
        # centred between full-grid points j-1 and j. Only m = 0 modes are
        # non-zero on the axis, so extrapolate those linearly to s = 0.
        xm_nyq = _array(ds, "xm_nyq")
        xn_nyq = _array(ds, "xn_nyq")
        bmnc = _array(ds, "bmnc")
        m0 = xm_nyq == 0
        if bmnc.shape[0] < 3:
            raise EquilibriumError("wout file has fewer than 3 radial points")
        b_axis_cos = 1.5 * bmnc[1, m0] - 0.5 * bmnc[2, m0]
        if lasym and "bmns" in ds.variables:
            bmns = _array(ds, "bmns")
            b_axis_sin = 1.5 * bmns[1, m0] - 0.5 * bmns[2, m0]
        else:
            b_axis_sin = np.zeros_like(b_axis_cos)

        reported = {}
        for key in ("Rmajor_p", "Aminor_p", "aspect", "volume_p", "b0", "volavgB",
                    "betatotal", "rbtor", "ier_flag"):
            if key in ds.variables:
                reported[key] = float(np.asarray(ds.variables[key][...]).item())

    if rmnc.ndim != 2 or rmnc.shape != zmns.shape or rmnc.shape[1] != xm.size:
        raise EquilibriumError("inconsistent array shapes in wout file")
    if not np.all(np.isfinite(rmnc)) or not np.all(np.isfinite(zmns)):
        raise EquilibriumError("non-finite geometry coefficients in wout file")
    if nfp < 1:
        raise EquilibriumError(f"invalid number of field periods: {nfp}")

    return Equilibrium(
        nfp=nfp, lasym=lasym, xm=xm, xn=xn, rmnc=rmnc, zmns=zmns, rmns=rmns,
        zmnc=zmnc, iota=iota, b_axis_xn=xn_nyq[m0], b_axis_cos=b_axis_cos,
        b_axis_sin=b_axis_sin, reported=reported, source=str(path),
    )

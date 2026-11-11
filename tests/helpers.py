"""Builders for synthetic test inputs with known analytic properties."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from stellarator_configprocessor.constants import MU0
from stellarator_configprocessor.io.coils import Coil, CoilSet

DATA_LOCAL = Path(__file__).resolve().parents[1] / "data" / "local"
W7X_WOUT = DATA_LOCAL / "wout_w7x_standard.nc"
W7X_COILS = DATA_LOCAL / "coils.w7x_standard"


def write_torus_wout(path, r0=10.0, a=2.0, nfp=5, b0=5.0, ns=11):
    """Write a wout file for a circular torus: R = r0 + a*sqrt(s)*cos(theta)."""
    import netCDF4

    s = np.linspace(0.0, 1.0, ns)
    xm = np.array([0.0, 1.0])
    xn = np.array([0.0, 0.0])
    rmnc = np.column_stack([np.full(ns, r0), a * np.sqrt(s)])
    zmns = np.column_stack([np.zeros(ns), a * np.sqrt(s)])
    bmnc = np.zeros((ns, 1))
    bmnc[1:, 0] = b0
    with netCDF4.Dataset(str(path), "w") as ds:
        ds.createDimension("radius", ns)
        ds.createDimension("mn_mode", 2)
        ds.createDimension("mn_mode_nyq", 1)

        def scalar(name, value, kind="f8"):
            ds.createVariable(name, kind)[...] = value

        scalar("nfp", nfp, "i4")
        scalar("lasym__logical__", 0, "i4")
        scalar("Rmajor_p", r0)
        scalar("Aminor_p", a)
        scalar("volume_p", 2.0 * np.pi**2 * r0 * a**2)
        scalar("b0", b0)
        ds.createVariable("xm", "f8", ("mn_mode",))[:] = xm
        ds.createVariable("xn", "f8", ("mn_mode",))[:] = xn
        ds.createVariable("xm_nyq", "f8", ("mn_mode_nyq",))[:] = [0.0]
        ds.createVariable("xn_nyq", "f8", ("mn_mode_nyq",))[:] = [0.0]
        ds.createVariable("rmnc", "f8", ("radius", "mn_mode"))[:] = rmnc
        ds.createVariable("zmns", "f8", ("radius", "mn_mode"))[:] = zmns
        ds.createVariable("bmnc", "f8", ("radius", "mn_mode_nyq"))[:] = bmnc
        ds.createVariable("iotaf", "f8", ("radius",))[:] = np.full(ns, 0.9)
    return Path(path)


def circular_tf_coils(n_coils=20, r0=10.0, a_coil=4.0, current=1.0e6, n_points=120) -> CoilSet:
    """Planar circular coils of radius ``a_coil`` centred on the circle R = r0."""
    theta = np.linspace(0.0, 2.0 * np.pi, n_points, endpoint=False)
    coils = []
    for k in range(n_coils):
        phi = 2.0 * np.pi * (k + 0.5) / n_coils
        r = r0 + a_coil * np.cos(theta)
        pts = np.column_stack([r * np.cos(phi), r * np.sin(phi), a_coil * np.sin(theta)])
        coils.append(Coil(points=pts, current=current, group=1, name="TF"))
    return CoilSet(coils=tuple(coils), periods=None)


def ideal_toroid_field(n_coils, current, r):
    """Field of an ideal toroidal solenoid at major radius ``r`` [T]."""
    return MU0 * n_coils * current / (2.0 * np.pi * r)

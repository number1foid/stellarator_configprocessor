"""Download the public W7-X files used for the reference check.

Files are written to ``data/local/`` (not tracked by git):

* ``wout_w7x_standard.nc``: free-boundary vacuum VMEC equilibrium of the W7-X
  standard configuration, from J. Geiger's database via the repository
  ``landreman/vmec_equilibria`` (which states no licence, so the file is
  downloaded and not redistributed here).
* ``W7-X.dat``: Fourier coefficients of the W7-X coil filaments shipped with
  SIMSOPT (MIT licence).
* ``coils.w7x_standard``: MAKEGRID file built from ``W7-X.dat`` by this script:
  the 50 non-planar coils with 15 kA x 108 turns each. The planar coils carry
  no current in the standard configuration and are left out.
* ``epseff_w7x_standard.txt``: effective ripple profile of the same
  configuration, also from ``landreman/vmec_equilibria``.

Every download is pinned to a commit and verified against a SHA-256 checksum.

Usage::

    python scripts/fetch_w7x_data.py
"""

from __future__ import annotations

import hashlib
import sys
import urllib.request
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "local"

SIMSOPT = "https://raw.githubusercontent.com/hiddenSymmetries/simsopt/912bdcc729d95433518c71df82b8cc45287908b3"
VMEC_EQ = "https://raw.githubusercontent.com/landreman/vmec_equilibria/41fcf8beb9fb743f0ee6e9c1397981e7ac7b29aa"

FILES = {
    "W7-X.dat": (
        f"{SIMSOPT}/src/simsopt/configs/W7-X.dat",
        "7c17c0c589026915414676fe47c10cd0aefa581c23dccbb390b08f8b73580de6",
    ),
    "wout_w7x_standard.nc": (
        f"{VMEC_EQ}/W7-X/Standard/wout.nc",
        "681aa2a94e4e55d88b91d8130ae6181cdd432154396db11e2d215693908e4131",
    ),
    "epseff_w7x_standard.txt": (
        f"{VMEC_EQ}/W7-X/Standard/epseff_bc.txt",
        "d1084938fd48a5db689ce5f3ed25235cca3b1ffa2474fbd87622383736f8bc2b",
    ),
}

#: W7-X standard configuration: current per non-planar coil [A] and turns.
W7X_CURRENT_PER_TURN = 15000.0
W7X_TURNS = 108
W7X_NFP = 5
#: The first five curves in W7-X.dat are the non-planar coils of one half period.
W7X_N_NONPLANAR = 5


def download(name: str) -> Path:
    url, sha = FILES[name]
    target = OUT / name
    if target.is_file() and hashlib.sha256(target.read_bytes()).hexdigest() == sha:
        print(f"  {name}: already present")
        return target
    print(f"  {name}: downloading")
    with urllib.request.urlopen(url, timeout=120) as response:
        data = response.read()
    digest = hashlib.sha256(data).hexdigest()
    if digest != sha:
        raise SystemExit(f"checksum mismatch for {name}: got {digest}, expected {sha}")
    target.write_bytes(data)
    return target


def fourier_curves(path: Path, n_points: int) -> list[np.ndarray]:
    """Evaluate the curves stored in SIMSOPT's Fourier table.

    Each coil has six columns: sin and cos coefficients of x, then y, then z.
    Row ``m`` holds the coefficients of mode ``m``.
    """
    table = np.loadtxt(path, delimiter=",")
    n_curves = table.shape[1] // 6
    modes = np.arange(table.shape[0])
    t = np.linspace(0.0, 2.0 * np.pi, n_points, endpoint=False)
    sin, cos = np.sin(np.outer(t, modes)), np.cos(np.outer(t, modes))
    curves = []
    for i in range(n_curves):
        cols = table[:, 6 * i:6 * i + 6]
        xyz = [sin @ cols[:, 2 * k] + cos @ cols[:, 2 * k + 1] for k in range(3)]
        curves.append(np.column_stack(xyz))
    return curves


def build_coils_file(dat: Path, n_points: int = 200) -> Path:
    sys.path.insert(0, str(ROOT / "src"))
    from stellarator_configprocessor.io.coils import (
        Coil, CoilSet, expand_symmetry, write_makegrid)

    curves = fourier_curves(dat, n_points)[:W7X_N_NONPLANAR]
    current = W7X_CURRENT_PER_TURN * W7X_TURNS
    base = CoilSet(coils=tuple(
        Coil(points=pts, current=current, group=1, name=f"NPC{i + 1}")
        for i, pts in enumerate(curves)))
    full = expand_symmetry(base, W7X_NFP, stellsym=True)
    target = OUT / "coils.w7x_standard"
    write_makegrid(target, full, W7X_NFP)
    print(f"  coils.w7x_standard: {len(full)} coils, {n_points} points each")
    return target


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"Writing to {OUT}")
    for name in FILES:
        download(name)
    build_coils_file(OUT / "W7-X.dat")
    print("Done.")


if __name__ == "__main__":
    main()

"""Download the public reference files used for the real-machine checks.

Files are written to ``data/local/`` (not tracked by git).

W7-X, standard configuration (5 field periods, 50 coils):

* ``wout_w7x_standard.nc``: free-boundary vacuum VMEC equilibrium, from
  J. Geiger's database via the repository ``landreman/vmec_equilibria``.
* ``coils.w7x_standard``: MAKEGRID file built by this script from SIMSOPT's
  ``W7-X.dat``: the 50 non-planar coils with 15 kA x 108 turns each. The planar
  coils carry no current in the standard configuration and are left out.
* ``epseff_w7x_standard.txt``: effective ripple profile of that configuration.

HSX, quasi-helically symmetric configuration (4 field periods, 48 coils):

* ``wout_hsx_qhs.nc``: free-boundary vacuum VMEC equilibrium, from
  ``landreman/vmec_equilibria``.
* ``coils.hsx_qhs``: MAKEGRID file built from SIMSOPT's ``HSX.dat``: the 48
  main coils with 150.07 kA-turns each.

``landreman/vmec_equilibria`` states no licence, so its files are downloaded
and not redistributed here. SIMSOPT is MIT-licensed. Every download is pinned
to a commit and verified against a SHA-256 checksum.

Usage::

    python scripts/fetch_reference_data.py
"""

from __future__ import annotations

import hashlib
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
    "HSX.dat": (
        f"{SIMSOPT}/src/simsopt/configs/HSX.dat",
        "d60968bec9138d1f7fd5bead349e6d58ee0aaa36dfc0a52b573c04021f48a477",
    ),
    "wout_hsx_qhs.nc": (
        f"{VMEC_EQ}/HSX/free_boundary/wout_HSX.nc",
        "f809f5f276f990bf900f4b78ecb887f2231521c62dea9513b88b653191e6c38d",
    ),
}

#: Coil sets to build: output name -> (Fourier table, number of base coils used,
#: field periods, current per coil in ampere-turns, coil name prefix).
#: W7-X: the first five curves are the non-planar coils of one half period,
#: 15 kA x 108 turns. HSX: six main coils per half period, current as in SIMSOPT.
COIL_SETS = {
    "coils.w7x_standard": ("W7-X.dat", 5, 5, 15000.0 * 108, "NPC"),
    "coils.hsx_qhs": ("HSX.dat", 6, 4, -1.500725500000000e05, "MC"),
}


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


def build_coils_file(name: str, n_points: int = 200) -> Path:
    """Write the full coil set in MAKEGRID format.

    Each base coil is rotated into every field period. With stellarator
    symmetry it is also mirrored, ``(x, y, z) -> (x, -y, -z)``, with the
    current reversed, which keeps the field stellarator-symmetric.
    """
    table, n_base, nfp, current, prefix = COIL_SETS[name]
    curves = fourier_curves(OUT / table, n_points)[:n_base]
    lines = [f"periods {nfp}", "begin filament", "mirror NIL"]
    count = 0
    for k in range(nfp):
        ang = 2.0 * np.pi * k / nfp
        rot = np.array([[np.cos(ang), -np.sin(ang), 0.0],
                        [np.sin(ang), np.cos(ang), 0.0],
                        [0.0, 0.0, 1.0]])
        for flip in (False, True):
            for i, base in enumerate(curves):
                pts = base * np.array([1.0, -1.0, -1.0]) if flip else base
                cur = -current if flip else current
                pts = pts @ rot.T
                for p in pts:
                    lines.append(f"{p[0]: .15E} {p[1]: .15E} {p[2]: .15E} {cur: .15E}")
                p = pts[0]
                lines.append(f"{p[0]: .15E} {p[1]: .15E} {p[2]: .15E} {0.0: .15E} 1 {prefix}{i + 1}")
                count += 1
    lines.append("end")
    (OUT / name).write_text("\n".join(lines) + "\n")
    print(f"  {name}: {count} coils, {n_points} points each")
    return OUT / name


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"Writing to {OUT}")
    for name in FILES:
        download(name)
    for name in COIL_SETS:
        build_coils_file(name)
    print("Done.")


if __name__ == "__main__":
    main()

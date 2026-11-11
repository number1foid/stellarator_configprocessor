"""Run PROCESS on native Windows.

PROCESS (tested with 3.4.2) is developed for Linux. On native Windows three
things stop a run, none related to the stellarator configuration file:

1. PROCESS writes non-ASCII characters to its output files using the default
   text encoding, which fails under cp1252. Fix: Python UTF-8 mode.
2. numba's on-disk cache can exceed the 260-character path limit. Fix: a
   short cache directory.
3. PROCESS deletes its temporary "idempotence" output files while they are
   still open, which Windows forbids. Fix: close them first.

This script applies the three fixes at run time, without changing the
installed PROCESS, and then calls the normal PROCESS command line.

Usage (with the Python environment that has PROCESS installed)::

    python scripts/run_process_windows.py -i path/to/run.IN.DAT

On Linux, macOS or WSL use the ``process`` command directly.
"""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from pathlib import Path


def main() -> int:
    if os.environ.get("PYTHONUTF8") != "1" or "NUMBA_CACHE_DIR" not in os.environ:
        env = dict(os.environ)
        env["PYTHONUTF8"] = "1"
        env.setdefault("NUMBA_CACHE_DIR", str(Path(tempfile.gettempdir()) / "nbc"))
        Path(env["NUMBA_CACHE_DIR"]).mkdir(parents=True, exist_ok=True)
        return subprocess.call([sys.executable, __file__, *sys.argv[1:]], env=env)

    from process.core.process_output import OutputFileManager
    from process.main import process_cli

    def close_idempotence_files(cls, output_prefix: str):
        names = [cls._outfile.name, cls._mfile.name]
        cls._outfile.close()
        cls._mfile.close()
        for name in names:
            Path(name).unlink(missing_ok=True)
        cls.open_files(mode="a", output_prefix=output_prefix)

    OutputFileManager.close_idempotence_files = classmethod(close_idempotence_files)
    sys.argv = ["process", *sys.argv[1:]]
    return process_cli()


if __name__ == "__main__":
    raise SystemExit(main())

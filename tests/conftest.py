import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))

from helpers import HSX_COILS, HSX_WOUT, W7X_COILS, W7X_WOUT  # noqa: E402

REFERENCE_FILES = {"w7x": (W7X_WOUT, W7X_COILS), "hsx": (HSX_WOUT, HSX_COILS)}


def pytest_collection_modifyitems(config, items):
    for marker, files in REFERENCE_FILES.items():
        if all(f.is_file() for f in files):
            continue
        skip = pytest.mark.skip(
            reason=f"{marker} files missing; run scripts/fetch_reference_data.py")
        for item in items:
            if marker in item.keywords:
                item.add_marker(skip)

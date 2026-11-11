import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))

from helpers import W7X_COILS, W7X_WOUT  # noqa: E402


def pytest_collection_modifyitems(config, items):
    if W7X_WOUT.is_file() and W7X_COILS.is_file():
        return
    skip = pytest.mark.skip(reason="W7-X files missing; run scripts/fetch_w7x_data.py")
    for item in items:
        if "w7x" in item.keywords:
            item.add_marker(skip)

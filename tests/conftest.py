import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import pytest  # noqa: E402
from model.microglia_state_model import integrate  # noqa: E402


@pytest.fixture(scope="session")
def sol():
    """One integrated v0 solution, reused across bound checks (cheap)."""
    return integrate()

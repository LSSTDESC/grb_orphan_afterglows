'''Utility functions to load SkySurvey data files.

The original data lives in ``masson/skysurvey/data`` and consists of
pickled objects (typically pandas DataFrames or dictionaries).  These
helpers provide a thin wrapper that locates the files relative to the
installed package and returns the deserialized objects.
'''

import os
import pickle
from pathlib import Path
from typing import Any

# Resolve the directory containing the data files.  When the package is
# installed in editable mode the path is ``<repo>/orphans/skysurvey``.
_DATA_DIR = Path(__file__).resolve().parents[2] / "data" / "skysurvey"

def _load_pickle(filename: str) -> Any:
    """Load a pickle file from the package data directory.

    Args:
        filename: Name of the file inside ``orphans/skysurvey/data``.
    Returns:
        The Python object stored in the pickle.
    """
    path = _DATA_DIR / filename
    if not path.is_file():
        raise FileNotFoundError(f"SkySurvey data file not found: {path}")
    with open(path, "rb") as f:
        return pickle.load(f)

# Public loaders -----------------------------------------------------------

def load_orphan_configs() -> Any:
    """Return the ``orphan_configs_ztf.pkl`` object.
    """
    return _load_pickle("orphan_configs_ztf.pkl")

def load_orphan_pseudo_obs_features() -> Any:
    """Return the ``orphan_pseudo_obs_features_ztf.pkl`` object.
    """
    return _load_pickle("orphan_pseudo_obs_features_ztf.pkl")

def load_orphan_pseudo_obs() -> Any:
    """Return the ``orphan_pseudo_obs_ztf.pkl`` object.
    """
    return _load_pickle("orphan_pseudo_obs_ztf.pkl")

def load_ztf_alerts_lc_features() -> Any:
    """Return the ``ztf_alerts_lc_features.pkl`` object.
    """
    return _load_pickle("ztf_alerts_lc_features.pkl")

def load_ztf_alerts_lc() -> Any:
    """Return the ``ztf_alerts_lc.pkl`` object.
    """
    return _load_pickle("ztf_alerts_lc.pkl")

"""Small provenance helpers shared by the benchmark preparation and runner."""

import hashlib
import subprocess
from pathlib import Path


def revision(directory):
    """Identify a checkout and whether it has tracked or untracked changes.

    Arguments:
        directory = path inside a Git repository.
    Returns:
        Commit hash and porcelain status; a dirty checkout stays visible.
    Raises:
        CalledProcessError if the directory is not a readable repository.
    """
    directory = Path(directory)
    commit = subprocess.check_output(
        ["git", "-C", str(directory), "rev-parse", "HEAD"], text=True,
    ).strip()
    status = subprocess.check_output(
        ["git", "-C", str(directory), "status", "--porcelain"], text=True,
    ).strip()
    return {"commit": commit, "status": status}


def sha256(filename):
    """Return a file's SHA-256 digest to detect changed scientific inputs.

    Arguments:
        filename = input file path.
    Returns:
        Hexadecimal digest of its bytes; no file is modified.
    """
    return hashlib.sha256(Path(filename).read_bytes()).hexdigest()


def power_table_record(tables, settings):
    """Fingerprint the installed covariance power grids without refining again.

    initialize() returns the installed, already refined set_cosmology inputs.
    Their k axis is log10(k/[h/Mpc]); each lnP array advances redshift first.
    Hash float64 bytes as well as shapes, so changing either grid is visible.
    """
    import numpy as np

    fingerprints = {}
    for name in ("log10k_2D", "z_2D", "lnP_linear", "lnP_nonlinear",
                 "lnP_linear_cb"):
        values = np.ascontiguousarray(tables[name], dtype=np.float64)
        fingerprints[name] = {
            "shape": list(values.shape),
            "sha256_float64": hashlib.sha256(values.tobytes()).hexdigest(),
        }
    return {
        "power_refinement": int(settings["power_refinement"]),
        "nk": len(tables["log10k_2D"]), "nz": len(tables["z_2D"]),
        "preparation": "project initialize; natural-cubic global refinement",
        "arrays": fingerprints,
    }

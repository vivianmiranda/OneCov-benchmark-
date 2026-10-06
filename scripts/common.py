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

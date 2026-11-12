"""Shared atomic HTTP download helper."""
from collections.abc import Callable
import hashlib
import os
from pathlib import Path
import tempfile
from typing import Optional, Union
from urllib.request import Request, urlopen


def download_atomic(
    source: Union[str, Request],
    destination: Union[str, Path],
    validate: Optional[Callable[[Path], None]] = None,
    timeout: float = 180,
) -> str:
    """Stream `source` into `destination` and return the payload's SHA256 hex digest.

    Data goes to a uniquely named temporary sibling file that replaces the
    destination only after the transfer completes (and `validate`, if given,
    accepts it), so an interrupted or concurrent download never leaves a
    truncated file in place.
    """
    destination = Path(destination)
    digest = hashlib.sha256()
    fd, temp_name = tempfile.mkstemp(prefix=f".{destination.name}.", suffix=".part", dir=destination.parent)
    temp_path = Path(temp_name)
    try:
        with os.fdopen(fd, "wb") as out, urlopen(source, timeout=timeout) as response:
            while block := response.read(1024 * 1024):
                out.write(block)
                digest.update(block)
        if validate is not None:
            validate(temp_path)
        temp_path.replace(destination)
    finally:
        temp_path.unlink(missing_ok=True)
    return digest.hexdigest()

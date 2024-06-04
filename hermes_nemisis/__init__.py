# Licensed under Apache License v2 - see LICENSE.rst
from pathlib import Path

from hermes_core import log
from hermes_nemisis.io.file_tools import read_file

try:
    from ._version import version as __version__
    from ._version import version_tuple
except ImportError:
    __version__ = "unknown version"
    version_tuple = (0, 0, "unknown version")

__all__ = ["log", "read_file"]

INST_NAME = "nemisis"
INST_SHORTNAME = "nem"
INST_TARGETNAME = "NEM"
INST_TO_SHORTNAME = {INST_NAME: INST_SHORTNAME}
INST_TO_TARGETNAME = {INST_NAME: INST_TARGETNAME}

_package_directory = Path(__file__).parent
_data_directory = _package_directory / "data"

log.info(f"hermes_nemesis version: {__version__}")

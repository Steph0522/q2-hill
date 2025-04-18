# flake8: noqa
# ----------------------------------------------------------------------------
# Copyright (c) 2024, Stephanie Hereira-Pacheco.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from . import _version
from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions


__version__ = _version.get_versions()["version"]

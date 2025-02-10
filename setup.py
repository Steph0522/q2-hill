# ----------------------------------------------------------------------------
# Copyright (c) 2024, Stephanie Hereira-Pacheco.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import find_packages, setup

import versioneer

description = (
    "Hill numbers diversity"
)

setup(
    name="q2-hill",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license="BSD-3-Clause",
    packages=find_packages(),
    author="Stephanie Hereira-Pacheco",
    author_email="shereirap@gmail.com",
    description=description,
    url="https://steph0522.github.io/website/",
    entry_points={
        "qiime2.plugins": [
            "q2_hill="
            "q2_hill"
            ".plugin_setup:plugin"]
    },
    package_data={
        "q2_hill": ["citations.bib"],
        "q2_hill.tests": ["data/*"],
    },
    zip_safe=False,
)

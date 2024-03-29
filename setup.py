#!/usr/bin/env python

from collections import defaultdict
import os
from setuptools import setup, Extension

# Build the C / Cython code
extensions = []

# Get Gala path, Numpy path
import gala
gala_base_path = os.path.split(gala.__file__)[0]

import numpy
numpy_base_path = os.path.split(numpy.__file__)[0]

cfg = defaultdict(list)
cfg['include_dirs'].append(os.path.join(numpy_base_path, 'core', 'include'))
cfg['include_dirs'].append(os.path.join(gala_base_path, 'potential'))
cfg['extra_compile_args'].append('--std=gnu99')
cfg['sources'].append('truncatedhernquist/cpotential.pyx')
cfg['sources'].append('truncatedhernquist/src/potential.c')
extensions.append(Extension('truncatedhernquist.potential', **cfg))

pkg_data = dict()
pkg_data[""] = ["LICENSE", "AUTHORS"]
pkg_data["truncatedhernquist"] = ["src/*.h", "src/*.c"]

VERSION_TEMPLATE = """
# Note that we need to fall back to the hard-coded version if either
# setuptools_scm can't be imported or setuptools_scm can't determine the
# version, so we catch the generic 'Exception'.
try:
    from setuptools_scm import get_version
    version = get_version(root='..', relative_to=__file__)
except Exception:
    version = '{version}'
""".lstrip()

setup(use_scm_version={'write_to': os.path.join('truncatedhernquist', 'version.py'),
                       'write_to_template': VERSION_TEMPLATE},
      ext_modules=extensions)

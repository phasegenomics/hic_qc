import os
import subprocess
from distutils.core import setup
from setuptools import setup

setup(
    name='bam_to_mate_hist',
    use_scm_version={
        'write_to': 'version.py',
    },
    setup_requires=['setuptools_scm'],
)

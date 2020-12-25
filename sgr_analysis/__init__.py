"""
sgr_analysis
Tools and scripts for the 2D-IR collaboration
"""

# Add imports here
from .sgr_analysis import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions

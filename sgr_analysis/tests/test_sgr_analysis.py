"""
Unit and regression test for the sgr_analysis package.
"""

# Import package, test suite, and other packages as needed
import sgr_analysis
import pytest
import sys

def test_sgr_analysis_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "sgr_analysis" in sys.modules

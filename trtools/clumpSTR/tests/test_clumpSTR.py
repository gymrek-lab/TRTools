import argparse
import os

import pytest

from ..clumpSTR import *

# Note vcfdir gives access to
# testsupport/sample_vcfs
# Can add directories for other test file types
# See pytest.ini

# Set up base argparser
@pytest.fixture
def args(tmpdir):
    args = argparse.ArgumentParser()
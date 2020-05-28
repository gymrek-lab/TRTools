import os, sys
import pytest

import trtools.utils.common as common

def test_MSG():
    common.MSG("Writing a test message", debug=False)
    common.MSG("Writing a test message", debug=True)

def test_WARNING():
    common.WARNING("Writing a test warning")
    common.WARNING("Writing a test warning")

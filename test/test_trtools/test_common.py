import os, sys
import pytest
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..','..'))

import trtools.utils.common as common

def test_MSG():
    common.MSG("Writing a test message", debug=False)
    common.MSG("Writing a test message", debug=True)

def test_WARNING():
    common.WARNING("Writing a test warning")
    common.WARNING("Writing a test warning")

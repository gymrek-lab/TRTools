"""
Common util functions
"""

import sys

def ERROR(msg):
    sys.stderr.write(msg.strip()+"\n")
    sys.exit(1)

def MSG(msg, debug=True):
    if debug:
        sys.stderr.write(msg.strip()+"\n")

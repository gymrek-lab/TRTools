"""
Common util functions
"""

import sys

def WARNING(msg):
    sys.stderr.write(msg.strip()+"\n")

def ERROR(msg):
    sys.stderr.write(msg.strip()+"\n")
    sys.exit(1)

def MSG(msg, debug=False):
    if debug:
        sys.stderr.write(msg.strip()+"\n")

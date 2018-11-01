"""
Common util functions
"""

import sys

def ERROR(msg):
    sys.stderr.write(msg.strip()+"\n")
    sys.exit(1)


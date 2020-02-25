"""
Version number
"""
import os
version_file = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "_version"))
__version__ = version_file.read().strip()


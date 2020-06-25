"""
Common util functions
"""

import sys

def WARNING(msg):
    r"""Write a warning message to standard error

    Parameters
    ----------
    msg : str
          A descriptive warning message

    Examples
    --------
    >>> WARNING("Something unexpected happened")
    """
    sys.stderr.write(msg.strip()+"\n")

def MSG(msg, debug=False):
    r"""Write a status message to standard error

    Parameters
    ----------
    msg : str
          A descriptive message
    debug : bool, optional
          Only print the message if debug is True

    Examples
    --------
    >>> MSG("Something unexpected happened")
    """
    if debug:
        sys.stderr.write(msg.strip()+"\n")

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

def ERROR(msg):
    r"""Write an error message to standard error

    .. deprecated:: 1.0.0
       Functions should instead raise exceptions that will be caught by main

    Parameters
    ----------
    msg : str
          A descriptive error message

    Examples
    --------
    >>> ERROR("Something bad happened. Existing")
    """
    sys.stderr.write(msg.strip()+"\n")
    sys.exit(1)

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

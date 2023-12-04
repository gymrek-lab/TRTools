try:
    from importlib.metadata import version, PackageNotFoundError
except ImportError:
    # handles py3.7, since importlib.metadata was introduced in py3.8
    from importlib_metadata import version, PackageNotFoundError

try:
    __version__ = version(__name__)
except PackageNotFoundError:
    __version__ = "unknown"

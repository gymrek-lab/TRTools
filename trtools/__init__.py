try:
    from .version import __version__
except ModuleNotFoundError:
    __version__ = "unknown"

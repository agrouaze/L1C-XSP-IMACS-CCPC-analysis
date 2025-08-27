from importlib.metadata import version

try:
    __version__ = version("L1C-XSP-IMACS-CCPC-analysis")
except Exception:
    __version__ = "999"
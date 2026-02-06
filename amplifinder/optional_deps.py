"""Check availability of optional dependencies."""

try:
    import matplotlib
    matplotlib.use('Agg')  # Use non-GUI backend
    import matplotlib.pyplot as plt
except ImportError:
    matplotlib = None
    plt = None

__all__ = ["matplotlib", "plt"]

"""AmpliFinder: Detect IS-mediated gene amplifications from WGS data."""

__version__ = "1.0.0"

from amplifinder.config import Config, load_config
from amplifinder.logger import setup_logger

__all__ = ["Config", "load_config", "setup_logger", "__version__"]


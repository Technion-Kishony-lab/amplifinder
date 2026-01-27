#!/usr/bin/env python3
"""Examples of logger color usage in AmpliFinder.

Demonstrates three ways to use colors:
1. Auto-highlighting (default rich behavior)
2. Inline colors for specific text (colorize/c helper)
3. Entire message color (color parameter)
"""

from amplifinder import logger, colorize, c

# Set up logger with colors enabled
from amplifinder.logger import setup_logger
setup_logger(use_colors=True, verbose=True)

print("\n=== 1. Auto-highlighting (default) ===")
logger.info("Found 10 records in data.csv")
logger.info("Processing paths: /path/to/file and /another/path")
logger.info("Values: 42, 3.14, 'string', True")

print("\n=== 2. Inline colors (colorize/c helper) ===")
logger.info(f"Status: {colorize('OK', 'green')}")
logger.info(f"Found {c('10', 'cyan')} records in {c('data.csv', 'magenta')}")
logger.info(f"Result: {c('PASS', 'bold green')} - {c('100%', 'cyan')} success")
logger.warning(f"Using fallback: {c('default.csv', 'yellow')}")

print("\n=== 3. Entire message color (color parameter) ===")
logger.info("Important notice", color="bold cyan")
logger.info("Critical step completed", color="bold magenta")
logger.info("Performance: 1.23s", color="dim")

print("\n=== 4. Mixed approach ===")
# Note: color= parameter overrides inline colors (entire message gets one color)
logger.info(f"Found {c('5', 'cyan')} items: processing...", color="bold")
# Without color parameter, inline colors work:
logger.info(f"File {c('output.csv', 'magenta')} saved")
# Combine colorize with auto-highlighting:
logger.info(f"Processed {c('output.csv', 'magenta')} with 100 records")

print("\n=== 5. Standard levels (show auto-coloring) ===")
logger.debug("Debug message with 42 and 'test'")
logger.info("Info message with 3.14")
logger.warning("Warning message with /path/to/file")
logger.error("Error message with 123 items")

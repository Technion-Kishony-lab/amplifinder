"""Common utilities for finding and executing external tools."""

import subprocess
import shutil
from pathlib import Path
from typing import List, Optional, Union

from amplifinder.logger import warning


def find_tool(
    tool_name: str,
    config_path: Optional[Path] = None,
    check_executable: bool = True,
) -> Optional[Path]:
    """Find tool executable path.

    Checks in order:
    1. config_path (if provided and exists)
    2. System PATH (using shutil.which)

    Args:
        tool_name: Name of the tool executable
        config_path: Optional path from config (can be file or directory)
        check_executable: If True, verify the path is executable

    Returns:
        Path to tool executable, or None if not found
    """
    # Try config path first
    if config_path is not None:
        config_path = Path(config_path)

        # If it's a directory, append tool_name
        if config_path.is_dir():
            tool_path = config_path / tool_name
        else:
            tool_path = config_path

        # Check if it exists and is executable
        if tool_path.exists():
            if not check_executable or tool_path.is_file():
                return tool_path.resolve()

    # Fall back to PATH
    path = shutil.which(tool_name)
    if path:
        return Path(path)

    return None


def get_tool_path(
    tool_name: str,
    config_path: Optional[Path] = None,
    required: bool = True,
) -> Path:
    """Get tool executable path, raising error if not found.

    Args:
        tool_name: Name of the tool executable
        config_path: Optional path from config
        required: If True, raise FileNotFoundError if tool not found

    Returns:
        Path to tool executable

    Raises:
        FileNotFoundError: If tool not found and required=True
    """
    tool_path = find_tool(tool_name, config_path)

    if tool_path is None:
        if required:
            config_hint = f" (config path: {config_path})" if config_path else ""
            raise FileNotFoundError(
                f"{tool_name} not found. "
                f"Set path in amplifinder.yaml or ensure {tool_name} is in PATH.{config_hint}"
            )
        return Path(tool_name)  # Fallback to name only

    return tool_path


def run_command(
    cmd: List[Union[str, Path]],
    check: bool = True,
    capture_output: bool = False,
    text: bool = False,
    error_msg: Optional[str] = None,
    verbose: bool = False,
) -> subprocess.CompletedProcess:
    """Run subprocess command with consistent error handling.

    Args:
        cmd: Command and arguments as list
        check: If True, raise RuntimeError on non-zero exit
        capture_output: If True, capture stdout/stderr (ignored if verbose=True)
        text: If True, decode output as text
        error_msg: Custom error message (default: uses stderr)
        verbose: If True, show command output in real-time (overrides capture_output)

    Returns:
        CompletedProcess result

    Raises:
        RuntimeError: If check=True and command fails
    """
    # Convert Path objects to strings
    cmd_str = [str(c) for c in cmd]

    # If verbose, let output pass through; otherwise capture
    if verbose:
        result = subprocess.run(
            cmd_str,
            check=False,
        )
    else:
        result = subprocess.run(
            cmd_str,
            check=False,
            capture_output=capture_output,
            text=text,
        )

    if check and result.returncode != 0:
        if error_msg:
            msg = error_msg
        elif not verbose and capture_output and result.stderr:
            msg = result.stderr.strip()
        else:
            msg = f"Command failed with exit code {result.returncode}"

        if not verbose and capture_output and result.stderr:
            warning(f"Command stderr: {result.stderr}")

        raise RuntimeError(msg)

    return result

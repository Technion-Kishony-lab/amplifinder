"""Utilities."""

from amplifinder.utils.fasta import (
    read_fasta_lengths,
    read_fastq_lengths,
    get_read_length_stats,
    ReadLengthStats,
)
from amplifinder.utils.file_lock import (
    locked_resource,
    DEFAULT_LOCK_TIMEOUT,
)
from amplifinder.utils.run_utils import (
    find_tool,
    get_tool_path,
    run_command,
)
from amplifinder.utils.file_utils import (
    ensure_dir,
    ensure_parent_dir,
    remove_file_or_dir,
)
from amplifinder.utils.json_utils import (
    compact_short_lists,
)

__all__ = [
    # FASTA utilities
    "read_fasta_lengths",
    "read_fastq_lengths",
    "get_read_length_stats",
    "ReadLengthStats",

    # File locking utilities
    "locked_resource",
    "DEFAULT_LOCK_TIMEOUT",

    # Tool utilities
    "find_tool",
    "get_tool_path",
    "run_command",

    # Path utilities
    "ensure_dir",
    "ensure_parent_dir",
    "remove_file_or_dir",

    # JSON utilities
    "compact_short_lists",

    # Flag utilities
    "MutableFlag",
]

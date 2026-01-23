from amplifinder.steps.jct_coverage.alignment_data import AlignmentData


class HitsContainer:
    """Container for alignment hits, supporting both deduplicated and non-deduplicated modes."""

    def __init__(self, ignore_duplicates: bool):
        self._data: dict[str, AlignmentData] | list[AlignmentData] = {} if ignore_duplicates else []
        self._ignore_duplicates = ignore_duplicates

    def add_hit(self, alignment: AlignmentData, read_id: str) -> bool:
        """
        Add an alignment hit. In deduplicated mode, skips if read_id exists.
        Returns True if the hit was added, False if it was skipped.
        """
        if self._ignore_duplicates:
            if read_id not in self._data:
                self._data[read_id] = alignment
                return True
        else:
            self._data.append(alignment)
            return True
        return False

    def __len__(self) -> int:
        return len(self._data)

    def __iter__(self):
        """Iterate over alignments."""
        return iter(self._data.values() if self._ignore_duplicates else self._data)

    def keys(self):
        """Get read_ids (only valid in deduplicated mode)."""
        assert self._ignore_duplicates, "keys() only available in deduplicated mode"
        return self._data.keys()

    def pop(self, read_id: str) -> AlignmentData:
        """Pop alignment by read_id (only valid in deduplicated mode)."""
        assert self._ignore_duplicates, "pop() only available in deduplicated mode"
        return self._data.pop(read_id)

    def __setitem__(self, read_id: str, alignment: AlignmentData):
        """Set alignment by read_id (only valid in deduplicated mode)."""
        assert self._ignore_duplicates, "setitem only available in deduplicated mode"
        self._data[read_id] = alignment

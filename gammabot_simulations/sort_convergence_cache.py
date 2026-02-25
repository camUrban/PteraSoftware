"""Sort convergence cache JSON files for human readability.

Simulation entries are sorted by their key components
(wake, wake_length, panel_ar, num_chordwise), and the _optimized_dt entry is sorted
by (panel_ar, num_chordwise) and placed last. Uses the same file locking mechanism as
the convergence analysis to allow safe concurrent access.
"""

from __future__ import annotations

import contextlib
import json
import sys
from pathlib import Path
from typing import Generator


@contextlib.contextmanager
def _lock_cache_file(cache_path: Path) -> Generator[None, None, None]:
    """Acquire an exclusive file lock for safe concurrent cache access.

    Uses a .lock file adjacent to the cache file. On Windows uses msvcrt locking, on
    Unix uses fcntl file locking.

    :param cache_path: Path to the cache file to lock.
    """
    lock_path = cache_path.with_suffix(cache_path.suffix + ".lock")
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    lock_fh = open(lock_path, "w")
    try:
        if sys.platform == "win32":
            import msvcrt

            # Write a byte so there is content to lock.
            lock_fh.write(" ")
            lock_fh.flush()
            lock_fh.seek(0)
            # LK_LOCK retries for approximately 10 seconds before raising OSError.
            msvcrt.locking(lock_fh.fileno(), msvcrt.LK_LOCK, 1)
        else:
            import fcntl

            fcntl.flock(lock_fh, fcntl.LOCK_EX)
        yield
    finally:
        try:
            if sys.platform == "win32":
                import msvcrt

                lock_fh.seek(0)
                msvcrt.locking(lock_fh.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                import fcntl

                fcntl.flock(lock_fh, fcntl.LOCK_UN)
        except OSError:
            pass
        lock_fh.close()


def _simulation_sort_key(key: str) -> tuple:
    """Parse a simulation cache key into a sortable tuple.

    Keys have the format "wake,wake_length,panel_ar,num_chordwise", for example
    "True,1,4,5". Each component is parsed to its natural type (bool for wake, int for
    the rest) so that sorting is numerical rather than lexicographic.

    :param key: The simulation cache key string.
    :return: A tuple (wake_as_int, wake_length, panel_ar, num_chordwise) for sorting.
    """
    parts = key.split(",")
    # wake is a bool string: "True" sorts after "False" (1 > 0)
    wake = 0 if parts[0].strip() == "False" else 1
    return (wake,) + tuple(int(p.strip()) for p in parts[1:])


def _optimized_dt_sort_key(key: str) -> tuple[int, ...]:
    """Parse an _optimized_dt key into a sortable tuple.

    Keys have the format "panel_ar,num_chordwise", for example "4,5".

    :param key: The _optimized_dt key string.
    :return: A tuple of ints for sorting.
    """
    return tuple(int(p.strip()) for p in key.split(","))


def sort_cache_file(cache_path: Path) -> None:
    """Sort a single convergence cache JSON file in place.

    :param cache_path: Path to the JSON cache file to sort.
    """
    with _lock_cache_file(cache_path):
        with open(cache_path, "r") as f:
            data = json.load(f)

        # Separate simulation entries from the _optimized_dt entry.
        optimized_dt = data.pop("_optimized_dt", None)

        # Sort simulation entries by parsed key components.
        sorted_data: dict = dict(
            sorted(data.items(), key=lambda item: _simulation_sort_key(item[0]))
        )

        # Sort and append _optimized_dt at the end.
        if optimized_dt is not None:
            sorted_data["_optimized_dt"] = dict(
                sorted(
                    optimized_dt.items(),
                    key=lambda item: _optimized_dt_sort_key(item[0]),
                )
            )

        with open(cache_path, "w") as f:
            json.dump(sorted_data, f, indent=2)


def main() -> None:
    cache_dir = Path(__file__).parent / "convergence_cache"

    if not cache_dir.is_dir():
        print(f"Cache directory not found: {cache_dir}")
        return

    json_files = sorted(cache_dir.glob("*.json"))

    if not json_files:
        print("No JSON cache files found.")
        return

    for json_file in json_files:
        print(f"Sorting {json_file.name}...")
        sort_cache_file(json_file)
        print(f"  Done.")

    print(f"\nSorted {len(json_files)} cache file(s).")


if __name__ == "__main__":
    main()

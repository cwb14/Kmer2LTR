"""Write a derived T_BITS_SCHEDULE into src/ltrk2p/align.py.

Doing this with a script rather than by hand keeps the shipped constant provably
identical to what bench/derive_schedule.py produced from the sweep, which is the
whole point of deriving it under a pre-registered rule.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bench.derive_schedule import as_python  # noqa: E402

BLOCK = re.compile(
    r"T_BITS_SCHEDULE: tuple\[tuple\[float, float\], \.\.\.\] = \((?:[^)]|\)[^\n])*?\n\)",
    re.S)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--schedule-json", default="bench/out/improve/schedule.json")
    ap.add_argument("--target", default="src/ltrk2p/align.py")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)

    data = json.loads(Path(args.schedule_json).read_text())
    literal = as_python([(float(b), float(t)) for b, t in data["schedule"]])
    literal = literal.replace("T_BITS_SCHEDULE = (",
                              "T_BITS_SCHEDULE: tuple[tuple[float, float], ...] = (")
    src = Path(args.target).read_text()
    if not BLOCK.search(src):
        print("apply_schedule: could not locate T_BITS_SCHEDULE", file=sys.stderr)
        return 1
    new = BLOCK.sub(literal.rstrip(), src, count=1)
    print(literal)
    if args.dry_run:
        return 0
    Path(args.target).write_text(new)
    print(f"apply_schedule: wrote {args.target}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())

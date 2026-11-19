"""Command-line access to the fixed research archive."""

import argparse
import json
import subprocess
import sys
from pathlib import Path

from .research import BANNER, MAX_CSV_BYTES, Archive, UploadError, csv_text, workspace


def main() -> int:
    parser = argparse.ArgumentParser(prog="discern", description=BANNER)
    commands = parser.add_subparsers(dest="command", required=True)
    replay = commands.add_parser("replay", help="Open a saved held-out result")
    replay.add_argument("--sample-id", required=True)
    batch = commands.add_parser("batch", help="Rank a CSV of recorded isolate IDs")
    batch.add_argument("input", type=Path)
    batch.add_argument("output", type=Path)
    evaluation = commands.add_parser("eval", help="Run the repository evaluation pipeline")
    evaluation.add_argument("--all", action="store_true", required=True)
    commands.add_parser("verify", help="Verify registered research results")
    args = parser.parse_args()
    root = workspace()
    if args.command in ("eval", "verify"):
        script = "eval/run_all.py" if args.command == "eval" else "eval/verify.py"
        return subprocess.run([sys.executable, str(root / script)], cwd=root, check=False).returncode

    archive = Archive(root)
    try:
        if args.command == "replay":
            print(json.dumps(archive.replay(args.sample_id), ensure_ascii=True, indent=2))
        else:
            if args.input.resolve() == args.output.resolve():
                parser.error("Choose a separate output file.")
            with args.input.open("rb") as stream:
                ids = archive.imported_ids(stream.read(MAX_CSV_BYTES + 1))
            rows = archive.queue(ids)
            args.output.write_text(
                csv_text(["rank", "sample_id", "call", "call_label", "probability", "fold"], rows),
                encoding="utf-8",
            )
    except KeyError:
        parser.error("Choose sample IDs from the bundled replay collection.")
    except UploadError as error:
        parser.error(str(error))
    except OSError:
        parser.error("Check the CSV input and output paths.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

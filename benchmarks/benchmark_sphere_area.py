"""Compatibility wrapper for the sphere-area benchmark."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from benchmarks.cases import get_case, run_case
from benchmarks.run_benchmarks import write_csv


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--suite", action="store_true", help="run all registered sphere cases")
    args = parser.parse_args()

    case_names = ["sphere_n104_gl", "sphere_n124_pullback"] if args.suite else ["sphere_n104_gl"]
    write_csv([run_case(get_case(name)) for name in case_names], sys.stdout)


if __name__ == "__main__":
    main()

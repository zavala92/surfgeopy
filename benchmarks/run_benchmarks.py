"""Run surfgeopy benchmark cases and emit CSV or JSON."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional

ROOT_DIR = Path(__file__).resolve().parents[1]
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from benchmarks.cases import BenchmarkValue, ROOT, SUITES, list_cases, run_case, select_cases


FIELDNAMES = (
    "name",
    "surface",
    "quantity",
    "mesh",
    "n_vertices",
    "n_faces",
    "interpolation_degree",
    "lp_degree",
    "refinement_level",
    "integration_degree",
    "quadrature_rule",
    "value",
    "reference_value",
    "absolute_error",
    "relative_error",
    "seconds",
    "n_quadrature_points",
)


def write_csv(records: Iterable[Dict[str, BenchmarkValue]], stream) -> None:
    writer = csv.DictWriter(stream, fieldnames=FIELDNAMES)
    writer.writeheader()
    for record in records:
        writer.writerow(record)


def write_json(records: Iterable[Dict[str, BenchmarkValue]], stream) -> None:
    json.dump(list(records), stream, indent=2)
    stream.write("\n")


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--suite", choices=sorted(SUITES), default="quick")
    parser.add_argument("--case", action="append", default=[], help="case name; may be passed more than once")
    parser.add_argument("--format", choices=("csv", "json"), default="csv")
    parser.add_argument("--output", help="write results to this file instead of stdout")
    parser.add_argument("--list", action="store_true", help="list available benchmark cases")
    args = parser.parse_args(argv)

    if args.list:
        for case in list_cases():
            print(f"{case.name}\t{case.surface}\t{case.mesh_path.relative_to(ROOT)}")
        return 0

    records = [run_case(case) for case in select_cases(args.suite, args.case)]
    if args.output:
        with open(args.output, "w", newline="") as stream:
            write_json(records, stream) if args.format == "json" else write_csv(records, stream)
    else:
        write_json(records, sys.stdout) if args.format == "json" else write_csv(records, sys.stdout)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

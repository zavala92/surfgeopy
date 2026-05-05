import io
import json

from benchmarks.cases import clone_case, get_case, run_case, select_cases
from benchmarks.run_benchmarks import write_csv, write_json


def test_select_cases_from_suite_and_names():
    quick_cases = select_cases("quick")
    named_cases = select_cases(names=["torus_n260_gl"])

    assert [case.name for case in quick_cases] == ["sphere_n104_gl"]
    assert [case.name for case in named_cases] == ["torus_n260_gl"]


def test_benchmark_writers_emit_csv_and_json():
    records = [{"name": "example", "surface": "sphere", "value": 1.0}]

    csv_stream = io.StringIO()
    json_stream = io.StringIO()
    write_csv(records, csv_stream)
    write_json(records, json_stream)

    assert csv_stream.getvalue().splitlines()[0].startswith("name,surface")
    assert json.loads(json_stream.getvalue())[0]["name"] == "example"


def test_run_case_smoke():
    smoke_case = clone_case(
        get_case("sphere_n104_gl"),
        name="sphere_smoke",
        interpolation_degree=2,
        refinement_level=0,
        integration_degree=2,
    )

    record = run_case(smoke_case)

    assert record["name"] == "sphere_smoke"
    assert record["surface"] == "sphere"
    assert record["n_faces"] > 0
    assert record["n_quadrature_points"] > 0
    assert record["absolute_error"] is not None


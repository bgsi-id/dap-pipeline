from pathlib import Path

from conftest import REPO_ROOT
from check_process_containers import find_missing


def write(root: Path, rel: str, text: str):
    path = root / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def catalog(root: Path, *dirs: str):
    entries = "".join(
        f"  - id: {d}\n    path: pipelines/{d}\n    descriptor: pipelines/{d}/pipeline.yaml\n"
        for d in dirs
    )
    write(root, "catalog.yaml", f'schema_version: "1.0"\npipelines:\n{entries}')


def test_real_catalogued_pipelines_declare_a_container_for_every_process():
    # The k8s executor refuses a process without an image
    # ("Missing container image for process ..."), after the run was already billed.
    assert find_missing(REPO_ROOT) == []


def test_process_without_container_is_reported(tmp_path):
    catalog(tmp_path, "bare")
    write(
        tmp_path,
        "pipelines/bare/main.nf",
        "process A {\n    script:\n    \"\"\"echo hi\"\"\"\n}\n",
    )

    assert find_missing(tmp_path) == [("bare", "A", "pipelines/bare/main.nf")]


def test_container_directive_is_accepted(tmp_path):
    catalog(tmp_path, "ok")
    write(
        tmp_path,
        "pipelines/ok/main.nf",
        "process A {\n    container params.img\n    script:\n    \"\"\"echo hi\"\"\"\n}\n",
    )

    assert find_missing(tmp_path) == []


def test_config_wide_default_container_is_accepted(tmp_path):
    catalog(tmp_path, "cfg")
    write(tmp_path, "pipelines/cfg/main.nf", "process A {\n    script:\n    \"\"\"x\"\"\"\n}\n")
    write(tmp_path, "pipelines/cfg/nextflow.config", "process.container = 'python:3.12-slim'\n")

    assert find_missing(tmp_path) == []


def test_included_module_processes_are_checked(tmp_path):
    catalog(tmp_path, "inc")
    write(
        tmp_path,
        "pipelines/inc/main.nf",
        "include { M } from '../base/m'\nprocess A {\n    container 'x'\n    script:\n    \"\"\"x\"\"\"\n}\n",
    )
    write(tmp_path, "pipelines/base/m.nf", "process M {\n    script:\n    \"\"\"x\"\"\"\n}\n")

    assert find_missing(tmp_path) == [("inc", "M", "pipelines/base/m.nf")]


def test_uncatalogued_pipelines_are_ignored(tmp_path):
    catalog(tmp_path, "ok")
    write(tmp_path, "pipelines/ok/main.nf", "process A {\n    container 'x'\n    script:\n    \"\"\"x\"\"\"\n}\n")
    write(tmp_path, "pipelines/draft/main.nf", "process B {\n    script:\n    \"\"\"x\"\"\"\n}\n")

    assert find_missing(tmp_path) == []

from pathlib import Path

import pytest
import yaml

from conftest import REPO_ROOT
from validate_catalog import validate


def make_repo(tmp_path: Path, entries, descriptor=None, entrypoint="main.nf", write_entrypoint=True):
    (tmp_path / "pipelines" / "demo").mkdir(parents=True)
    (tmp_path / "catalog.yaml").write_text(yaml.safe_dump({"schema_version": "1.0", "pipelines": entries}))
    if descriptor is not False:
        (tmp_path / "pipelines" / "demo" / "pipeline.yaml").write_text(
            yaml.safe_dump(descriptor or {"id": "demo", "entrypoint": entrypoint})
        )
    if write_entrypoint:
        (tmp_path / "pipelines" / "demo" / "main.nf").write_text("workflow {}")
    return tmp_path


GOOD = {"id": "demo", "path": "pipelines/demo", "descriptor": "pipelines/demo/pipeline.yaml", "active": True}


def errors_for(tmp_path, entries=None, **kwargs):
    return validate(make_repo(tmp_path, entries or [dict(GOOD)], **kwargs))[0]


def test_real_repository_catalog_is_valid():
    errors, infos = validate(REPO_ROOT)
    assert errors == []


def test_valid_catalog_passes_and_unlisted_dirs_are_informational(tmp_path):
    root = make_repo(tmp_path, [dict(GOOD)])
    (root / "pipelines" / "other").mkdir()
    errors, infos = validate(root)
    assert errors == []
    assert any("other" in i for i in infos)


def test_missing_catalog(tmp_path):
    assert validate(tmp_path)[0] == ["catalog.yaml not found"]


@pytest.mark.parametrize("missing", ["id", "path", "descriptor"])
def test_missing_required_key(tmp_path, missing):
    entry = dict(GOOD)
    del entry[missing]
    assert any("missing or empty" in e for e in errors_for(tmp_path, [entry]))


def test_active_must_be_boolean(tmp_path):
    entry = dict(GOOD, active="yes")
    assert any("'active' must be a boolean" in e for e in errors_for(tmp_path, [entry]))


def test_duplicate_ids(tmp_path):
    assert any("duplicate id" in e for e in errors_for(tmp_path, [dict(GOOD), dict(GOOD)]))


@pytest.mark.parametrize(
    "bad_path",
    ["/etc", "pipelines/../x", "pipelines", "pipelines/demo/sub", "other/demo", "./pipelines/demo"],
)
def test_unsafe_pipeline_path(tmp_path, bad_path):
    entry = dict(GOOD, path=bad_path)
    assert errors_for(tmp_path, [entry])


def test_descriptor_outside_path(tmp_path):
    entry = dict(GOOD, descriptor="catalog.yaml")
    assert any("not inside" in e for e in errors_for(tmp_path, [entry]))


def test_descriptor_missing_file(tmp_path):
    assert any("does not exist" in e for e in errors_for(tmp_path, descriptor=False))


def test_descriptor_id_mismatch(tmp_path):
    errors = errors_for(tmp_path, descriptor={"id": "other", "entrypoint": "main.nf"})
    assert any("differs from catalog id" in e for e in errors)


def test_entrypoint_missing(tmp_path):
    assert any("entrypoint" in e for e in errors_for(tmp_path, write_entrypoint=False))


def test_entrypoint_traversal(tmp_path):
    assert any("entrypoint" in e for e in errors_for(tmp_path, entrypoint="../x.nf"))


def test_descriptor_invalid_yaml(tmp_path):
    root = make_repo(tmp_path, [dict(GOOD)])
    (root / "pipelines/demo/pipeline.yaml").write_text("id: [unclosed")
    assert any("not valid YAML" in e for e in validate(root)[0])

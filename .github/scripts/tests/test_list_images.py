import json
from pathlib import Path

from conftest import REPO_ROOT
from list_images import REF_RE, discover, image_name, main


def write(root: Path, rel: str, text: str):
    path = root / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


DIGEST = "@sha256:" + "a" * 64


def test_real_repository_images():
    images, skipped, errors = discover(REPO_ROOT)
    refs = {i["image_ref"] for i in images}
    assert errors == []
    assert "python:3.12-slim" in refs
    assert "ensemblorg/ensembl-vep:release_116.0" in refs
    assert "hailgenetics/hail:0.2.133" in refs
    assert "mambaorg/micromamba:1.5.10" in refs
    assert any(i["pinned"] for i in images)
    assert any("variant_etl" in s for s in skipped)  # variant_image = null


def test_image_name():
    assert image_name("python:3.12-slim") == "python"
    assert image_name("quay.io/biocontainers/plink2:2.00a5.12--h9948957_1") == "biocontainers/plink2"
    assert image_name(f"quay.io/biocontainers/echtvar{DIGEST}") == "biocontainers/echtvar"
    assert image_name("localhost:5000/team/app:1") == "team/app"
    assert image_name("ensemblorg/ensembl-vep:release_116.0") == "ensemblorg/ensembl-vep"


def test_ref_pattern_rejects_injection():
    for bad in ["a b", "x;rm -rf /", "$(id)", "img:tag`id`", "UP per", ""]:
        assert not REF_RE.match(bad), bad


def test_groovy_sources_dedupe_and_pinned(tmp_path):
    write(tmp_path, "pipelines/a/nextflow.config", f"params.x_image = 'quay.io/o/x:1'\nparams.y_image = \"quay.io/o/y{DIGEST}\"\n")
    write(tmp_path, "pipelines/a/main.nf", "process P {\n    container 'quay.io/o/x:1'\n}\n")
    write(tmp_path, "pipelines/a/params.example.json", json.dumps({"x_image": "quay.io/o/x:1", "other": "no"}))
    images, _, errors = discover(tmp_path)
    assert errors == []
    by_ref = {i["image_ref"]: i for i in images}
    assert len(images) == 2
    assert len(by_ref["quay.io/o/x:1"]["sources"]) == 3
    assert by_ref["quay.io/o/x:1"]["pinned"] is False
    assert by_ref[f"quay.io/o/y{DIGEST}"]["pinned"] is True


def test_commented_lines_are_ignored(tmp_path):
    write(tmp_path, "pipelines/a/main.nf", "// container 'quay.io/o/old:1'\n    container 'quay.io/o/new:1'\n")
    images, _, _ = discover(tmp_path)
    assert [i["image_ref"] for i in images] == ["quay.io/o/new:1"]


def test_null_and_dynamic_values_are_skipped_not_dropped_silently(tmp_path):
    write(tmp_path, "pipelines/a/nextflow.config", "params.v_image = null\nparams.w_image = '${registry}/w:1'\nparams.z_image = ''\n")
    write(tmp_path, "pipelines/a/params.example.json", json.dumps({"q_image": None}))
    images, skipped, errors = discover(tmp_path)
    assert images == [] and errors == []
    assert len(skipped) == 4


def test_dockerfile_from_skips_scratch_and_stage_aliases(tmp_path):
    write(
        tmp_path,
        "pipelines/a/Dockerfile",
        "FROM --platform=linux/amd64 golang:1.22 AS build\nFROM build AS test\nFROM scratch\nFROM alpine:3.20\n",
    )
    images, _, _ = discover(tmp_path)
    assert sorted(i["image_ref"] for i in images) == ["alpine:3.20", "golang:1.22"]


def test_malformed_reference_fails(tmp_path):
    write(tmp_path, "pipelines/a/nextflow.config", "params.x_image = 'bad ref;id'\n")
    _, _, errors = discover(tmp_path)
    assert errors and "malformed" in errors[0]


def test_main_exit_codes_and_enforce_flag(tmp_path, monkeypatch):
    write(tmp_path, "pipelines/a/nextflow.config", "params.x_image = 'quay.io/o/x:1'\n")
    monkeypatch.delenv("GITHUB_OUTPUT", raising=False)
    monkeypatch.delenv("GITHUB_STEP_SUMMARY", raising=False)
    assert main(["x", str(tmp_path)]) == 0
    monkeypatch.setenv("ENFORCE_DIGEST_PINS", "true")
    assert main(["x", str(tmp_path)]) == 1


def test_main_writes_github_output(tmp_path, monkeypatch):
    write(tmp_path, "pipelines/a/nextflow.config", "params.x_image = 'quay.io/o/x:1'\n")
    out = tmp_path / "out.txt"
    monkeypatch.setenv("GITHUB_OUTPUT", str(out))
    monkeypatch.delenv("GITHUB_STEP_SUMMARY", raising=False)
    monkeypatch.delenv("ENFORCE_DIGEST_PINS", raising=False)
    assert main(["x", str(tmp_path)]) == 0
    text = out.read_text()
    assert text.startswith("images=") and "count=1" in text

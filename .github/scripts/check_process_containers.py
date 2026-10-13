#!/usr/bin/env python3
"""Fail when a catalogued pipeline has a Nextflow process without a container.

DAP runs workflows on the Kubernetes executor, which needs an image for every
process. Without one the run is started (and its koin reserved) and then dies
with "Missing container image for process ...". A process is accepted when its
block has a `container` directive or the pipeline's nextflow.config sets a
process-wide `process.container`. Processes pulled in through `include` are
checked too. Pipelines that are not in catalog.yaml are ignored.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import yaml

PROCESS_RE = re.compile(r"^\s*process\s+(\w+)\s*\{", re.M)
INCLUDE_RE = re.compile(r"include\s*\{[^}]*\}\s*from\s*['\"]([^'\"]+)['\"]")
CONFIG_DEFAULT_RE = re.compile(
    r"process\s*\.\s*container\s*=|process\s*\{[^}]*?\bcontainer\s*=", re.S
)
DIRECTIVE_RE = re.compile(r"^\s*container\s+\S", re.M)


def _process_blocks(text: str):
    for match in PROCESS_RE.finditer(text):
        depth, index = 1, match.end()
        while index < len(text) and depth:
            depth += text[index] == "{"
            depth -= text[index] == "}"
            index += 1
        yield match.group(1), text[match.end():index]


def _included_files(source: Path):
    for target in INCLUDE_RE.findall(source.read_text()):
        base = (source.parent / target).resolve()
        for candidate in (base, base.with_suffix(".nf")):
            if candidate.is_file():
                yield candidate
                break


def find_missing(root: Path) -> list[tuple[str, str, str]]:
    """Return (pipeline directory, process, file) for each process lacking an image."""
    root = Path(root).resolve()
    catalog = yaml.safe_load((root / "catalog.yaml").read_text()) or {}
    missing: list[tuple[str, str, str]] = []
    for entry in catalog.get("pipelines") or []:
        pipeline_dir = root / entry["path"]
        config = pipeline_dir / "nextflow.config"
        config_default = config.is_file() and bool(
            CONFIG_DEFAULT_RE.search(config.read_text())
        )
        files = sorted(pipeline_dir.rglob("*.nf"))
        for source in list(files):
            for included in _included_files(source):
                if included not in files:
                    files.append(included)
        for source in files:
            for name, body in _process_blocks(source.read_text()):
                if not config_default and not DIRECTIVE_RE.search(body):
                    missing.append(
                        (pipeline_dir.name, name, str(source.relative_to(root)))
                    )
    return missing


def main(argv: list[str]) -> int:
    root = Path(argv[1]) if len(argv) > 1 else Path(".")
    missing = find_missing(root)
    for pipeline, process, source in missing:
        print(
            f"{source}: process {process} (pipeline {pipeline}) has no container "
            "directive; the Kubernetes executor will refuse it"
        )
    return 1 if missing else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
